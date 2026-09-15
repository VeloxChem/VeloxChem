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


#include "SimdElectronRepulsionGeom10VrrRecKI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_ki_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_0 = buffer.data(ii + 0);
    const auto *ii_1 = buffer.data(ii + 1);
    const auto *ii_2 = buffer.data(ii + 2);
    const auto *ii_3 = buffer.data(ii + 3);
    const auto *ii_4 = buffer.data(ii + 4);
    const auto *ii_5 = buffer.data(ii + 5);
    const auto *ii_6 = buffer.data(ii + 6);
    const auto *ii_7 = buffer.data(ii + 7);
    const auto *ii_8 = buffer.data(ii + 8);
    const auto *ii_9 = buffer.data(ii + 9);
    const auto *ii_10 = buffer.data(ii + 10);
    const auto *ii_11 = buffer.data(ii + 11);
    const auto *ii_12 = buffer.data(ii + 12);
    const auto *ii_13 = buffer.data(ii + 13);
    const auto *ii_14 = buffer.data(ii + 14);
    const auto *ii_15 = buffer.data(ii + 15);
    const auto *ii_16 = buffer.data(ii + 16);
    const auto *ii_17 = buffer.data(ii + 17);
    const auto *ii_18 = buffer.data(ii + 18);
    const auto *ii_19 = buffer.data(ii + 19);
    const auto *ii_20 = buffer.data(ii + 20);
    const auto *ii_21 = buffer.data(ii + 21);
    const auto *ii_22 = buffer.data(ii + 22);
    const auto *ii_23 = buffer.data(ii + 23);
    const auto *ii_24 = buffer.data(ii + 24);
    const auto *ii_25 = buffer.data(ii + 25);
    const auto *ii_26 = buffer.data(ii + 26);
    const auto *ii_27 = buffer.data(ii + 27);
    const auto *ii_28 = buffer.data(ii + 28);
    const auto *ii_29 = buffer.data(ii + 29);
    const auto *ii_30 = buffer.data(ii + 30);
    const auto *ii_31 = buffer.data(ii + 31);
    const auto *ii_32 = buffer.data(ii + 32);
    const auto *ii_33 = buffer.data(ii + 33);
    const auto *ii_34 = buffer.data(ii + 34);
    const auto *ii_35 = buffer.data(ii + 35);
    const auto *ii_36 = buffer.data(ii + 36);
    const auto *ii_37 = buffer.data(ii + 37);
    const auto *ii_38 = buffer.data(ii + 38);
    const auto *ii_39 = buffer.data(ii + 39);
    const auto *ii_40 = buffer.data(ii + 40);
    const auto *ii_41 = buffer.data(ii + 41);
    const auto *ii_42 = buffer.data(ii + 42);
    const auto *ii_43 = buffer.data(ii + 43);
    const auto *ii_44 = buffer.data(ii + 44);
    const auto *ii_45 = buffer.data(ii + 45);
    const auto *ii_46 = buffer.data(ii + 46);
    const auto *ii_47 = buffer.data(ii + 47);
    const auto *ii_48 = buffer.data(ii + 48);
    const auto *ii_49 = buffer.data(ii + 49);
    const auto *ii_50 = buffer.data(ii + 50);
    const auto *ii_51 = buffer.data(ii + 51);
    const auto *ii_52 = buffer.data(ii + 52);
    const auto *ii_53 = buffer.data(ii + 53);
    const auto *ii_54 = buffer.data(ii + 54);
    const auto *ii_55 = buffer.data(ii + 55);
    const auto *ii_56 = buffer.data(ii + 56);
    const auto *ii_57 = buffer.data(ii + 57);
    const auto *ii_58 = buffer.data(ii + 58);
    const auto *ii_59 = buffer.data(ii + 59);
    const auto *ii_60 = buffer.data(ii + 60);
    const auto *ii_61 = buffer.data(ii + 61);
    const auto *ii_62 = buffer.data(ii + 62);
    const auto *ii_63 = buffer.data(ii + 63);
    const auto *ii_64 = buffer.data(ii + 64);
    const auto *ii_65 = buffer.data(ii + 65);
    const auto *ii_66 = buffer.data(ii + 66);
    const auto *ii_67 = buffer.data(ii + 67);
    const auto *ii_68 = buffer.data(ii + 68);
    const auto *ii_69 = buffer.data(ii + 69);
    const auto *ii_70 = buffer.data(ii + 70);
    const auto *ii_71 = buffer.data(ii + 71);
    const auto *ii_72 = buffer.data(ii + 72);
    const auto *ii_73 = buffer.data(ii + 73);
    const auto *ii_74 = buffer.data(ii + 74);
    const auto *ii_75 = buffer.data(ii + 75);
    const auto *ii_76 = buffer.data(ii + 76);
    const auto *ii_77 = buffer.data(ii + 77);
    const auto *ii_78 = buffer.data(ii + 78);
    const auto *ii_79 = buffer.data(ii + 79);
    const auto *ii_80 = buffer.data(ii + 80);
    const auto *ii_81 = buffer.data(ii + 81);
    const auto *ii_82 = buffer.data(ii + 82);
    const auto *ii_83 = buffer.data(ii + 83);
    const auto *ii_84 = buffer.data(ii + 84);
    const auto *ii_85 = buffer.data(ii + 85);
    const auto *ii_86 = buffer.data(ii + 86);
    const auto *ii_87 = buffer.data(ii + 87);
    const auto *ii_88 = buffer.data(ii + 88);
    const auto *ii_89 = buffer.data(ii + 89);
    const auto *ii_90 = buffer.data(ii + 90);
    const auto *ii_91 = buffer.data(ii + 91);
    const auto *ii_92 = buffer.data(ii + 92);
    const auto *ii_93 = buffer.data(ii + 93);
    const auto *ii_94 = buffer.data(ii + 94);
    const auto *ii_95 = buffer.data(ii + 95);
    const auto *ii_96 = buffer.data(ii + 96);
    const auto *ii_97 = buffer.data(ii + 97);
    const auto *ii_98 = buffer.data(ii + 98);
    const auto *ii_99 = buffer.data(ii + 99);
    const auto *ii_100 = buffer.data(ii + 100);
    const auto *ii_101 = buffer.data(ii + 101);
    const auto *ii_102 = buffer.data(ii + 102);
    const auto *ii_103 = buffer.data(ii + 103);
    const auto *ii_104 = buffer.data(ii + 104);
    const auto *ii_105 = buffer.data(ii + 105);
    const auto *ii_106 = buffer.data(ii + 106);
    const auto *ii_107 = buffer.data(ii + 107);
    const auto *ii_108 = buffer.data(ii + 108);
    const auto *ii_109 = buffer.data(ii + 109);
    const auto *ii_110 = buffer.data(ii + 110);
    const auto *ii_111 = buffer.data(ii + 111);
    const auto *ii_112 = buffer.data(ii + 112);
    const auto *ii_113 = buffer.data(ii + 113);
    const auto *ii_114 = buffer.data(ii + 114);
    const auto *ii_115 = buffer.data(ii + 115);
    const auto *ii_116 = buffer.data(ii + 116);
    const auto *ii_117 = buffer.data(ii + 117);
    const auto *ii_118 = buffer.data(ii + 118);
    const auto *ii_119 = buffer.data(ii + 119);
    const auto *ii_120 = buffer.data(ii + 120);
    const auto *ii_121 = buffer.data(ii + 121);
    const auto *ii_122 = buffer.data(ii + 122);
    const auto *ii_123 = buffer.data(ii + 123);
    const auto *ii_124 = buffer.data(ii + 124);
    const auto *ii_125 = buffer.data(ii + 125);
    const auto *ii_126 = buffer.data(ii + 126);
    const auto *ii_127 = buffer.data(ii + 127);
    const auto *ii_128 = buffer.data(ii + 128);
    const auto *ii_129 = buffer.data(ii + 129);
    const auto *ii_130 = buffer.data(ii + 130);
    const auto *ii_131 = buffer.data(ii + 131);
    const auto *ii_132 = buffer.data(ii + 132);
    const auto *ii_133 = buffer.data(ii + 133);
    const auto *ii_134 = buffer.data(ii + 134);
    const auto *ii_135 = buffer.data(ii + 135);
    const auto *ii_136 = buffer.data(ii + 136);
    const auto *ii_137 = buffer.data(ii + 137);
    const auto *ii_138 = buffer.data(ii + 138);
    const auto *ii_139 = buffer.data(ii + 139);
    const auto *ii_140 = buffer.data(ii + 140);
    const auto *ii_141 = buffer.data(ii + 141);
    const auto *ii_142 = buffer.data(ii + 142);
    const auto *ii_143 = buffer.data(ii + 143);
    const auto *ii_144 = buffer.data(ii + 144);
    const auto *ii_145 = buffer.data(ii + 145);
    const auto *ii_146 = buffer.data(ii + 146);
    const auto *ii_147 = buffer.data(ii + 147);
    const auto *ii_148 = buffer.data(ii + 148);
    const auto *ii_149 = buffer.data(ii + 149);

    const auto *li_0 = buffer.data(li + 0);
    const auto *li_1 = buffer.data(li + 1);
    const auto *li_2 = buffer.data(li + 2);
    const auto *li_3 = buffer.data(li + 3);
    const auto *li_4 = buffer.data(li + 4);
    const auto *li_5 = buffer.data(li + 5);
    const auto *li_6 = buffer.data(li + 6);
    const auto *li_7 = buffer.data(li + 7);
    const auto *li_8 = buffer.data(li + 8);
    const auto *li_9 = buffer.data(li + 9);
    const auto *li_10 = buffer.data(li + 10);
    const auto *li_11 = buffer.data(li + 11);
    const auto *li_12 = buffer.data(li + 12);
    const auto *li_13 = buffer.data(li + 13);
    const auto *li_14 = buffer.data(li + 14);
    const auto *li_15 = buffer.data(li + 15);
    const auto *li_16 = buffer.data(li + 16);
    const auto *li_17 = buffer.data(li + 17);
    const auto *li_18 = buffer.data(li + 18);
    const auto *li_19 = buffer.data(li + 19);
    const auto *li_20 = buffer.data(li + 20);
    const auto *li_21 = buffer.data(li + 21);
    const auto *li_22 = buffer.data(li + 22);
    const auto *li_23 = buffer.data(li + 23);
    const auto *li_24 = buffer.data(li + 24);
    const auto *li_25 = buffer.data(li + 25);
    const auto *li_26 = buffer.data(li + 26);
    const auto *li_27 = buffer.data(li + 27);
    const auto *li_28 = buffer.data(li + 28);
    const auto *li_29 = buffer.data(li + 29);
    const auto *li_30 = buffer.data(li + 30);
    const auto *li_31 = buffer.data(li + 31);
    const auto *li_32 = buffer.data(li + 32);
    const auto *li_33 = buffer.data(li + 33);
    const auto *li_34 = buffer.data(li + 34);
    const auto *li_35 = buffer.data(li + 35);
    const auto *li_36 = buffer.data(li + 36);
    const auto *li_37 = buffer.data(li + 37);
    const auto *li_38 = buffer.data(li + 38);
    const auto *li_39 = buffer.data(li + 39);
    const auto *li_40 = buffer.data(li + 40);
    const auto *li_41 = buffer.data(li + 41);
    const auto *li_42 = buffer.data(li + 42);
    const auto *li_43 = buffer.data(li + 43);
    const auto *li_44 = buffer.data(li + 44);
    const auto *li_45 = buffer.data(li + 45);
    const auto *li_46 = buffer.data(li + 46);
    const auto *li_47 = buffer.data(li + 47);
    const auto *li_48 = buffer.data(li + 48);
    const auto *li_49 = buffer.data(li + 49);
    const auto *li_50 = buffer.data(li + 50);
    const auto *li_51 = buffer.data(li + 51);
    const auto *li_52 = buffer.data(li + 52);
    const auto *li_53 = buffer.data(li + 53);
    const auto *li_54 = buffer.data(li + 54);
    const auto *li_55 = buffer.data(li + 55);
    const auto *li_56 = buffer.data(li + 56);
    const auto *li_57 = buffer.data(li + 57);
    const auto *li_58 = buffer.data(li + 58);
    const auto *li_59 = buffer.data(li + 59);
    const auto *li_60 = buffer.data(li + 60);
    const auto *li_61 = buffer.data(li + 61);
    const auto *li_62 = buffer.data(li + 62);
    const auto *li_63 = buffer.data(li + 63);
    const auto *li_64 = buffer.data(li + 64);
    const auto *li_65 = buffer.data(li + 65);
    const auto *li_66 = buffer.data(li + 66);
    const auto *li_67 = buffer.data(li + 67);
    const auto *li_68 = buffer.data(li + 68);
    const auto *li_69 = buffer.data(li + 69);
    const auto *li_70 = buffer.data(li + 70);
    const auto *li_71 = buffer.data(li + 71);
    const auto *li_72 = buffer.data(li + 72);
    const auto *li_73 = buffer.data(li + 73);
    const auto *li_74 = buffer.data(li + 74);
    const auto *li_75 = buffer.data(li + 75);
    const auto *li_76 = buffer.data(li + 76);
    const auto *li_77 = buffer.data(li + 77);
    const auto *li_78 = buffer.data(li + 78);
    const auto *li_79 = buffer.data(li + 79);
    const auto *li_80 = buffer.data(li + 80);
    const auto *li_81 = buffer.data(li + 81);
    const auto *li_82 = buffer.data(li + 82);
    const auto *li_83 = buffer.data(li + 83);
    const auto *li_84 = buffer.data(li + 84);
    const auto *li_85 = buffer.data(li + 85);
    const auto *li_86 = buffer.data(li + 86);
    const auto *li_87 = buffer.data(li + 87);
    const auto *li_88 = buffer.data(li + 88);
    const auto *li_89 = buffer.data(li + 89);
    const auto *li_90 = buffer.data(li + 90);
    const auto *li_91 = buffer.data(li + 91);
    const auto *li_92 = buffer.data(li + 92);
    const auto *li_93 = buffer.data(li + 93);
    const auto *li_94 = buffer.data(li + 94);
    const auto *li_95 = buffer.data(li + 95);
    const auto *li_96 = buffer.data(li + 96);
    const auto *li_97 = buffer.data(li + 97);
    const auto *li_98 = buffer.data(li + 98);
    const auto *li_99 = buffer.data(li + 99);
    const auto *li_100 = buffer.data(li + 100);
    const auto *li_101 = buffer.data(li + 101);
    const auto *li_102 = buffer.data(li + 102);
    const auto *li_103 = buffer.data(li + 103);
    const auto *li_104 = buffer.data(li + 104);
    const auto *li_105 = buffer.data(li + 105);
    const auto *li_106 = buffer.data(li + 106);
    const auto *li_107 = buffer.data(li + 107);
    const auto *li_108 = buffer.data(li + 108);
    const auto *li_109 = buffer.data(li + 109);
    const auto *li_110 = buffer.data(li + 110);
    const auto *li_111 = buffer.data(li + 111);
    const auto *li_112 = buffer.data(li + 112);
    const auto *li_113 = buffer.data(li + 113);
    const auto *li_114 = buffer.data(li + 114);
    const auto *li_115 = buffer.data(li + 115);
    const auto *li_116 = buffer.data(li + 116);
    const auto *li_117 = buffer.data(li + 117);
    const auto *li_118 = buffer.data(li + 118);
    const auto *li_119 = buffer.data(li + 119);
    const auto *li_120 = buffer.data(li + 120);
    const auto *li_121 = buffer.data(li + 121);
    const auto *li_122 = buffer.data(li + 122);
    const auto *li_123 = buffer.data(li + 123);
    const auto *li_124 = buffer.data(li + 124);
    const auto *li_125 = buffer.data(li + 125);
    const auto *li_126 = buffer.data(li + 126);
    const auto *li_127 = buffer.data(li + 127);
    const auto *li_128 = buffer.data(li + 128);
    const auto *li_129 = buffer.data(li + 129);
    const auto *li_130 = buffer.data(li + 130);
    const auto *li_131 = buffer.data(li + 131);
    const auto *li_132 = buffer.data(li + 132);
    const auto *li_133 = buffer.data(li + 133);
    const auto *li_134 = buffer.data(li + 134);
    const auto *li_135 = buffer.data(li + 135);
    const auto *li_136 = buffer.data(li + 136);
    const auto *li_137 = buffer.data(li + 137);
    const auto *li_138 = buffer.data(li + 138);
    const auto *li_139 = buffer.data(li + 139);
    const auto *li_140 = buffer.data(li + 140);
    const auto *li_141 = buffer.data(li + 141);
    const auto *li_142 = buffer.data(li + 142);
    const auto *li_143 = buffer.data(li + 143);
    const auto *li_144 = buffer.data(li + 144);
    const auto *li_145 = buffer.data(li + 145);
    const auto *li_146 = buffer.data(li + 146);
    const auto *li_147 = buffer.data(li + 147);
    const auto *li_148 = buffer.data(li + 148);
    const auto *li_149 = buffer.data(li + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ii_0, ii_1, ii_2, ii_3, ii_4, li_0, li_1, \
                         li_2, li_3, li_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -7.0 * ii_0[k]
                 + f_0 * li_0[k];

        t_1[k] = -7.0 * ii_1[k]
                 + f_0 * li_1[k];

        t_2[k] = -7.0 * ii_2[k]
                 + f_0 * li_2[k];

        t_3[k] = -7.0 * ii_3[k]
                 + f_0 * li_3[k];

        t_4[k] = -7.0 * ii_4[k]
                 + f_0 * li_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ii_5, ii_6, ii_7, ii_8, ii_9, li_5, li_6, \
                         li_7, li_8, li_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -7.0 * ii_5[k]
                 + f_0 * li_5[k];

        t_6[k] = -7.0 * ii_6[k]
                 + f_0 * li_6[k];

        t_7[k] = -7.0 * ii_7[k]
                 + f_0 * li_7[k];

        t_8[k] = -7.0 * ii_8[k]
                 + f_0 * li_8[k];

        t_9[k] = -7.0 * ii_9[k]
                 + f_0 * li_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ii_10, ii_11, ii_12, ii_13, ii_14, \
                         li_10, li_11, li_12, li_13, li_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -7.0 * ii_10[k]
                  + f_0 * li_10[k];

        t_11[k] = -7.0 * ii_11[k]
                  + f_0 * li_11[k];

        t_12[k] = -7.0 * ii_12[k]
                  + f_0 * li_12[k];

        t_13[k] = -7.0 * ii_13[k]
                  + f_0 * li_13[k];

        t_14[k] = -7.0 * ii_14[k]
                  + f_0 * li_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ii_15, ii_16, ii_17, ii_18, ii_19, \
                         li_15, li_16, li_17, li_18, li_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -7.0 * ii_15[k]
                  + f_0 * li_15[k];

        t_16[k] = -7.0 * ii_16[k]
                  + f_0 * li_16[k];

        t_17[k] = -7.0 * ii_17[k]
                  + f_0 * li_17[k];

        t_18[k] = -7.0 * ii_18[k]
                  + f_0 * li_18[k];

        t_19[k] = -7.0 * ii_19[k]
                  + f_0 * li_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ii_20, ii_21, ii_22, ii_23, ii_24, \
                         li_20, li_21, li_22, li_23, li_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -7.0 * ii_20[k]
                  + f_0 * li_20[k];

        t_21[k] = -7.0 * ii_21[k]
                  + f_0 * li_21[k];

        t_22[k] = -7.0 * ii_22[k]
                  + f_0 * li_22[k];

        t_23[k] = -7.0 * ii_23[k]
                  + f_0 * li_23[k];

        t_24[k] = -7.0 * ii_24[k]
                  + f_0 * li_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ii_25, ii_26, ii_27, ii_28, ii_29, \
                         li_25, li_26, li_27, li_28, li_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -7.0 * ii_25[k]
                  + f_0 * li_25[k];

        t_26[k] = -7.0 * ii_26[k]
                  + f_0 * li_26[k];

        t_27[k] = -7.0 * ii_27[k]
                  + f_0 * li_27[k];

        t_28[k] = -6.0 * ii_28[k]
                  + f_0 * li_28[k];

        t_29[k] = -6.0 * ii_29[k]
                  + f_0 * li_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ii_30, ii_31, ii_32, ii_33, ii_34, \
                         li_30, li_31, li_32, li_33, li_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -6.0 * ii_30[k]
                  + f_0 * li_30[k];

        t_31[k] = -6.0 * ii_31[k]
                  + f_0 * li_31[k];

        t_32[k] = -6.0 * ii_32[k]
                  + f_0 * li_32[k];

        t_33[k] = -6.0 * ii_33[k]
                  + f_0 * li_33[k];

        t_34[k] = -6.0 * ii_34[k]
                  + f_0 * li_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ii_35, ii_36, ii_37, ii_38, ii_39, \
                         li_35, li_36, li_37, li_38, li_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -6.0 * ii_35[k]
                  + f_0 * li_35[k];

        t_36[k] = -6.0 * ii_36[k]
                  + f_0 * li_36[k];

        t_37[k] = -6.0 * ii_37[k]
                  + f_0 * li_37[k];

        t_38[k] = -6.0 * ii_38[k]
                  + f_0 * li_38[k];

        t_39[k] = -6.0 * ii_39[k]
                  + f_0 * li_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ii_40, ii_41, ii_42, ii_43, ii_44, \
                         li_40, li_41, li_42, li_43, li_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -6.0 * ii_40[k]
                  + f_0 * li_40[k];

        t_41[k] = -6.0 * ii_41[k]
                  + f_0 * li_41[k];

        t_42[k] = -6.0 * ii_42[k]
                  + f_0 * li_42[k];

        t_43[k] = -6.0 * ii_43[k]
                  + f_0 * li_43[k];

        t_44[k] = -6.0 * ii_44[k]
                  + f_0 * li_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ii_45, ii_46, ii_47, ii_48, ii_49, \
                         li_45, li_46, li_47, li_48, li_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -6.0 * ii_45[k]
                  + f_0 * li_45[k];

        t_46[k] = -6.0 * ii_46[k]
                  + f_0 * li_46[k];

        t_47[k] = -6.0 * ii_47[k]
                  + f_0 * li_47[k];

        t_48[k] = -6.0 * ii_48[k]
                  + f_0 * li_48[k];

        t_49[k] = -6.0 * ii_49[k]
                  + f_0 * li_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ii_50, ii_51, ii_52, ii_53, ii_54, \
                         li_50, li_51, li_52, li_53, li_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -6.0 * ii_50[k]
                  + f_0 * li_50[k];

        t_51[k] = -6.0 * ii_51[k]
                  + f_0 * li_51[k];

        t_52[k] = -6.0 * ii_52[k]
                  + f_0 * li_52[k];

        t_53[k] = -6.0 * ii_53[k]
                  + f_0 * li_53[k];

        t_54[k] = -6.0 * ii_54[k]
                  + f_0 * li_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ii_55, ii_56, ii_57, ii_58, ii_59, \
                         li_55, li_56, li_57, li_58, li_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -6.0 * ii_55[k]
                  + f_0 * li_55[k];

        t_56[k] = -6.0 * ii_56[k]
                  + f_0 * li_56[k];

        t_57[k] = -6.0 * ii_57[k]
                  + f_0 * li_57[k];

        t_58[k] = -6.0 * ii_58[k]
                  + f_0 * li_58[k];

        t_59[k] = -6.0 * ii_59[k]
                  + f_0 * li_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ii_60, ii_61, ii_62, ii_63, ii_64, \
                         li_60, li_61, li_62, li_63, li_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -6.0 * ii_60[k]
                  + f_0 * li_60[k];

        t_61[k] = -6.0 * ii_61[k]
                  + f_0 * li_61[k];

        t_62[k] = -6.0 * ii_62[k]
                  + f_0 * li_62[k];

        t_63[k] = -6.0 * ii_63[k]
                  + f_0 * li_63[k];

        t_64[k] = -6.0 * ii_64[k]
                  + f_0 * li_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ii_65, ii_66, ii_67, ii_68, ii_69, \
                         li_65, li_66, li_67, li_68, li_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -6.0 * ii_65[k]
                  + f_0 * li_65[k];

        t_66[k] = -6.0 * ii_66[k]
                  + f_0 * li_66[k];

        t_67[k] = -6.0 * ii_67[k]
                  + f_0 * li_67[k];

        t_68[k] = -6.0 * ii_68[k]
                  + f_0 * li_68[k];

        t_69[k] = -6.0 * ii_69[k]
                  + f_0 * li_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ii_70, ii_71, ii_72, ii_73, ii_74, \
                         li_70, li_71, li_72, li_73, li_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -6.0 * ii_70[k]
                  + f_0 * li_70[k];

        t_71[k] = -6.0 * ii_71[k]
                  + f_0 * li_71[k];

        t_72[k] = -6.0 * ii_72[k]
                  + f_0 * li_72[k];

        t_73[k] = -6.0 * ii_73[k]
                  + f_0 * li_73[k];

        t_74[k] = -6.0 * ii_74[k]
                  + f_0 * li_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ii_75, ii_76, ii_77, ii_78, ii_79, \
                         li_75, li_76, li_77, li_78, li_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -6.0 * ii_75[k]
                  + f_0 * li_75[k];

        t_76[k] = -6.0 * ii_76[k]
                  + f_0 * li_76[k];

        t_77[k] = -6.0 * ii_77[k]
                  + f_0 * li_77[k];

        t_78[k] = -6.0 * ii_78[k]
                  + f_0 * li_78[k];

        t_79[k] = -6.0 * ii_79[k]
                  + f_0 * li_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ii_80, ii_81, ii_82, ii_83, ii_84, \
                         li_80, li_81, li_82, li_83, li_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -6.0 * ii_80[k]
                  + f_0 * li_80[k];

        t_81[k] = -6.0 * ii_81[k]
                  + f_0 * li_81[k];

        t_82[k] = -6.0 * ii_82[k]
                  + f_0 * li_82[k];

        t_83[k] = -6.0 * ii_83[k]
                  + f_0 * li_83[k];

        t_84[k] = -5.0 * ii_84[k]
                  + f_0 * li_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ii_85, ii_86, ii_87, ii_88, ii_89, \
                         li_85, li_86, li_87, li_88, li_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -5.0 * ii_85[k]
                  + f_0 * li_85[k];

        t_86[k] = -5.0 * ii_86[k]
                  + f_0 * li_86[k];

        t_87[k] = -5.0 * ii_87[k]
                  + f_0 * li_87[k];

        t_88[k] = -5.0 * ii_88[k]
                  + f_0 * li_88[k];

        t_89[k] = -5.0 * ii_89[k]
                  + f_0 * li_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ii_90, ii_91, ii_92, ii_93, ii_94, \
                         li_90, li_91, li_92, li_93, li_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * ii_90[k]
                  + f_0 * li_90[k];

        t_91[k] = -5.0 * ii_91[k]
                  + f_0 * li_91[k];

        t_92[k] = -5.0 * ii_92[k]
                  + f_0 * li_92[k];

        t_93[k] = -5.0 * ii_93[k]
                  + f_0 * li_93[k];

        t_94[k] = -5.0 * ii_94[k]
                  + f_0 * li_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ii_95, ii_96, ii_97, ii_98, ii_99, \
                         li_95, li_96, li_97, li_98, li_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * ii_95[k]
                  + f_0 * li_95[k];

        t_96[k] = -5.0 * ii_96[k]
                  + f_0 * li_96[k];

        t_97[k] = -5.0 * ii_97[k]
                  + f_0 * li_97[k];

        t_98[k] = -5.0 * ii_98[k]
                  + f_0 * li_98[k];

        t_99[k] = -5.0 * ii_99[k]
                  + f_0 * li_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ii_100, ii_101, ii_102, ii_103, \
                         ii_104, li_100, li_101, li_102, li_103, \
                         li_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -5.0 * ii_100[k]
                   + f_0 * li_100[k];

        t_101[k] = -5.0 * ii_101[k]
                   + f_0 * li_101[k];

        t_102[k] = -5.0 * ii_102[k]
                   + f_0 * li_102[k];

        t_103[k] = -5.0 * ii_103[k]
                   + f_0 * li_103[k];

        t_104[k] = -5.0 * ii_104[k]
                   + f_0 * li_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ii_105, ii_106, ii_107, ii_108, \
                         ii_109, li_105, li_106, li_107, li_108, \
                         li_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -5.0 * ii_105[k]
                   + f_0 * li_105[k];

        t_106[k] = -5.0 * ii_106[k]
                   + f_0 * li_106[k];

        t_107[k] = -5.0 * ii_107[k]
                   + f_0 * li_107[k];

        t_108[k] = -5.0 * ii_108[k]
                   + f_0 * li_108[k];

        t_109[k] = -5.0 * ii_109[k]
                   + f_0 * li_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ii_110, ii_111, ii_112, ii_113, \
                         ii_114, li_110, li_111, li_112, li_113, \
                         li_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -5.0 * ii_110[k]
                   + f_0 * li_110[k];

        t_111[k] = -5.0 * ii_111[k]
                   + f_0 * li_111[k];

        t_112[k] = -5.0 * ii_112[k]
                   + f_0 * li_112[k];

        t_113[k] = -5.0 * ii_113[k]
                   + f_0 * li_113[k];

        t_114[k] = -5.0 * ii_114[k]
                   + f_0 * li_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ii_115, ii_116, ii_117, ii_118, \
                         ii_119, li_115, li_116, li_117, li_118, \
                         li_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -5.0 * ii_115[k]
                   + f_0 * li_115[k];

        t_116[k] = -5.0 * ii_116[k]
                   + f_0 * li_116[k];

        t_117[k] = -5.0 * ii_117[k]
                   + f_0 * li_117[k];

        t_118[k] = -5.0 * ii_118[k]
                   + f_0 * li_118[k];

        t_119[k] = -5.0 * ii_119[k]
                   + f_0 * li_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ii_120, ii_121, ii_122, ii_123, \
                         ii_124, li_120, li_121, li_122, li_123, \
                         li_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -5.0 * ii_120[k]
                   + f_0 * li_120[k];

        t_121[k] = -5.0 * ii_121[k]
                   + f_0 * li_121[k];

        t_122[k] = -5.0 * ii_122[k]
                   + f_0 * li_122[k];

        t_123[k] = -5.0 * ii_123[k]
                   + f_0 * li_123[k];

        t_124[k] = -5.0 * ii_124[k]
                   + f_0 * li_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ii_125, ii_126, ii_127, ii_128, \
                         ii_129, li_125, li_126, li_127, li_128, \
                         li_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -5.0 * ii_125[k]
                   + f_0 * li_125[k];

        t_126[k] = -5.0 * ii_126[k]
                   + f_0 * li_126[k];

        t_127[k] = -5.0 * ii_127[k]
                   + f_0 * li_127[k];

        t_128[k] = -5.0 * ii_128[k]
                   + f_0 * li_128[k];

        t_129[k] = -5.0 * ii_129[k]
                   + f_0 * li_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ii_130, ii_131, ii_132, ii_133, \
                         ii_134, li_130, li_131, li_132, li_133, \
                         li_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -5.0 * ii_130[k]
                   + f_0 * li_130[k];

        t_131[k] = -5.0 * ii_131[k]
                   + f_0 * li_131[k];

        t_132[k] = -5.0 * ii_132[k]
                   + f_0 * li_132[k];

        t_133[k] = -5.0 * ii_133[k]
                   + f_0 * li_133[k];

        t_134[k] = -5.0 * ii_134[k]
                   + f_0 * li_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ii_135, ii_136, ii_137, ii_138, \
                         ii_139, li_135, li_136, li_137, li_138, \
                         li_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -5.0 * ii_135[k]
                   + f_0 * li_135[k];

        t_136[k] = -5.0 * ii_136[k]
                   + f_0 * li_136[k];

        t_137[k] = -5.0 * ii_137[k]
                   + f_0 * li_137[k];

        t_138[k] = -5.0 * ii_138[k]
                   + f_0 * li_138[k];

        t_139[k] = -5.0 * ii_139[k]
                   + f_0 * li_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ii_140, ii_141, ii_142, ii_143, \
                         ii_144, li_140, li_141, li_142, li_143, \
                         li_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -5.0 * ii_140[k]
                   + f_0 * li_140[k];

        t_141[k] = -5.0 * ii_141[k]
                   + f_0 * li_141[k];

        t_142[k] = -5.0 * ii_142[k]
                   + f_0 * li_142[k];

        t_143[k] = -5.0 * ii_143[k]
                   + f_0 * li_143[k];

        t_144[k] = -5.0 * ii_144[k]
                   + f_0 * li_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ii_145, ii_146, ii_147, ii_148, \
                         ii_149, li_145, li_146, li_147, li_148, \
                         li_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -5.0 * ii_145[k]
                   + f_0 * li_145[k];

        t_146[k] = -5.0 * ii_146[k]
                   + f_0 * li_146[k];

        t_147[k] = -5.0 * ii_147[k]
                   + f_0 * li_147[k];

        t_148[k] = -5.0 * ii_148[k]
                   + f_0 * li_148[k];

        t_149[k] = -5.0 * ii_149[k]
                   + f_0 * li_149[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_150 = buffer.data(ii + 150);
    const auto *ii_151 = buffer.data(ii + 151);
    const auto *ii_152 = buffer.data(ii + 152);
    const auto *ii_153 = buffer.data(ii + 153);
    const auto *ii_154 = buffer.data(ii + 154);
    const auto *ii_155 = buffer.data(ii + 155);
    const auto *ii_156 = buffer.data(ii + 156);
    const auto *ii_157 = buffer.data(ii + 157);
    const auto *ii_158 = buffer.data(ii + 158);
    const auto *ii_159 = buffer.data(ii + 159);
    const auto *ii_160 = buffer.data(ii + 160);
    const auto *ii_161 = buffer.data(ii + 161);
    const auto *ii_162 = buffer.data(ii + 162);
    const auto *ii_163 = buffer.data(ii + 163);
    const auto *ii_164 = buffer.data(ii + 164);
    const auto *ii_165 = buffer.data(ii + 165);
    const auto *ii_166 = buffer.data(ii + 166);
    const auto *ii_167 = buffer.data(ii + 167);
    const auto *ii_168 = buffer.data(ii + 168);
    const auto *ii_169 = buffer.data(ii + 169);
    const auto *ii_170 = buffer.data(ii + 170);
    const auto *ii_171 = buffer.data(ii + 171);
    const auto *ii_172 = buffer.data(ii + 172);
    const auto *ii_173 = buffer.data(ii + 173);
    const auto *ii_174 = buffer.data(ii + 174);
    const auto *ii_175 = buffer.data(ii + 175);
    const auto *ii_176 = buffer.data(ii + 176);
    const auto *ii_177 = buffer.data(ii + 177);
    const auto *ii_178 = buffer.data(ii + 178);
    const auto *ii_179 = buffer.data(ii + 179);
    const auto *ii_180 = buffer.data(ii + 180);
    const auto *ii_181 = buffer.data(ii + 181);
    const auto *ii_182 = buffer.data(ii + 182);
    const auto *ii_183 = buffer.data(ii + 183);
    const auto *ii_184 = buffer.data(ii + 184);
    const auto *ii_185 = buffer.data(ii + 185);
    const auto *ii_186 = buffer.data(ii + 186);
    const auto *ii_187 = buffer.data(ii + 187);
    const auto *ii_188 = buffer.data(ii + 188);
    const auto *ii_189 = buffer.data(ii + 189);
    const auto *ii_190 = buffer.data(ii + 190);
    const auto *ii_191 = buffer.data(ii + 191);
    const auto *ii_192 = buffer.data(ii + 192);
    const auto *ii_193 = buffer.data(ii + 193);
    const auto *ii_194 = buffer.data(ii + 194);
    const auto *ii_195 = buffer.data(ii + 195);
    const auto *ii_196 = buffer.data(ii + 196);
    const auto *ii_197 = buffer.data(ii + 197);
    const auto *ii_198 = buffer.data(ii + 198);
    const auto *ii_199 = buffer.data(ii + 199);
    const auto *ii_200 = buffer.data(ii + 200);
    const auto *ii_201 = buffer.data(ii + 201);
    const auto *ii_202 = buffer.data(ii + 202);
    const auto *ii_203 = buffer.data(ii + 203);
    const auto *ii_204 = buffer.data(ii + 204);
    const auto *ii_205 = buffer.data(ii + 205);
    const auto *ii_206 = buffer.data(ii + 206);
    const auto *ii_207 = buffer.data(ii + 207);
    const auto *ii_208 = buffer.data(ii + 208);
    const auto *ii_209 = buffer.data(ii + 209);
    const auto *ii_210 = buffer.data(ii + 210);
    const auto *ii_211 = buffer.data(ii + 211);
    const auto *ii_212 = buffer.data(ii + 212);
    const auto *ii_213 = buffer.data(ii + 213);
    const auto *ii_214 = buffer.data(ii + 214);
    const auto *ii_215 = buffer.data(ii + 215);
    const auto *ii_216 = buffer.data(ii + 216);
    const auto *ii_217 = buffer.data(ii + 217);
    const auto *ii_218 = buffer.data(ii + 218);
    const auto *ii_219 = buffer.data(ii + 219);
    const auto *ii_220 = buffer.data(ii + 220);
    const auto *ii_221 = buffer.data(ii + 221);
    const auto *ii_222 = buffer.data(ii + 222);
    const auto *ii_223 = buffer.data(ii + 223);
    const auto *ii_224 = buffer.data(ii + 224);
    const auto *ii_225 = buffer.data(ii + 225);
    const auto *ii_226 = buffer.data(ii + 226);
    const auto *ii_227 = buffer.data(ii + 227);
    const auto *ii_228 = buffer.data(ii + 228);
    const auto *ii_229 = buffer.data(ii + 229);
    const auto *ii_230 = buffer.data(ii + 230);
    const auto *ii_231 = buffer.data(ii + 231);
    const auto *ii_232 = buffer.data(ii + 232);
    const auto *ii_233 = buffer.data(ii + 233);
    const auto *ii_234 = buffer.data(ii + 234);
    const auto *ii_235 = buffer.data(ii + 235);
    const auto *ii_236 = buffer.data(ii + 236);
    const auto *ii_237 = buffer.data(ii + 237);
    const auto *ii_238 = buffer.data(ii + 238);
    const auto *ii_239 = buffer.data(ii + 239);
    const auto *ii_240 = buffer.data(ii + 240);
    const auto *ii_241 = buffer.data(ii + 241);
    const auto *ii_242 = buffer.data(ii + 242);
    const auto *ii_243 = buffer.data(ii + 243);
    const auto *ii_244 = buffer.data(ii + 244);
    const auto *ii_245 = buffer.data(ii + 245);
    const auto *ii_246 = buffer.data(ii + 246);
    const auto *ii_247 = buffer.data(ii + 247);
    const auto *ii_248 = buffer.data(ii + 248);
    const auto *ii_249 = buffer.data(ii + 249);
    const auto *ii_250 = buffer.data(ii + 250);
    const auto *ii_251 = buffer.data(ii + 251);
    const auto *ii_252 = buffer.data(ii + 252);
    const auto *ii_253 = buffer.data(ii + 253);
    const auto *ii_254 = buffer.data(ii + 254);
    const auto *ii_255 = buffer.data(ii + 255);
    const auto *ii_256 = buffer.data(ii + 256);
    const auto *ii_257 = buffer.data(ii + 257);
    const auto *ii_258 = buffer.data(ii + 258);
    const auto *ii_259 = buffer.data(ii + 259);
    const auto *ii_260 = buffer.data(ii + 260);
    const auto *ii_261 = buffer.data(ii + 261);
    const auto *ii_262 = buffer.data(ii + 262);
    const auto *ii_263 = buffer.data(ii + 263);
    const auto *ii_264 = buffer.data(ii + 264);
    const auto *ii_265 = buffer.data(ii + 265);
    const auto *ii_266 = buffer.data(ii + 266);
    const auto *ii_267 = buffer.data(ii + 267);
    const auto *ii_268 = buffer.data(ii + 268);
    const auto *ii_269 = buffer.data(ii + 269);
    const auto *ii_270 = buffer.data(ii + 270);
    const auto *ii_271 = buffer.data(ii + 271);
    const auto *ii_272 = buffer.data(ii + 272);
    const auto *ii_273 = buffer.data(ii + 273);
    const auto *ii_274 = buffer.data(ii + 274);
    const auto *ii_275 = buffer.data(ii + 275);
    const auto *ii_276 = buffer.data(ii + 276);
    const auto *ii_277 = buffer.data(ii + 277);
    const auto *ii_278 = buffer.data(ii + 278);
    const auto *ii_279 = buffer.data(ii + 279);
    const auto *ii_280 = buffer.data(ii + 280);
    const auto *ii_281 = buffer.data(ii + 281);
    const auto *ii_282 = buffer.data(ii + 282);
    const auto *ii_283 = buffer.data(ii + 283);
    const auto *ii_284 = buffer.data(ii + 284);
    const auto *ii_285 = buffer.data(ii + 285);
    const auto *ii_286 = buffer.data(ii + 286);
    const auto *ii_287 = buffer.data(ii + 287);
    const auto *ii_288 = buffer.data(ii + 288);
    const auto *ii_289 = buffer.data(ii + 289);
    const auto *ii_290 = buffer.data(ii + 290);
    const auto *ii_291 = buffer.data(ii + 291);
    const auto *ii_292 = buffer.data(ii + 292);
    const auto *ii_293 = buffer.data(ii + 293);
    const auto *ii_294 = buffer.data(ii + 294);
    const auto *ii_295 = buffer.data(ii + 295);
    const auto *ii_296 = buffer.data(ii + 296);
    const auto *ii_297 = buffer.data(ii + 297);
    const auto *ii_298 = buffer.data(ii + 298);
    const auto *ii_299 = buffer.data(ii + 299);

    const auto *li_150 = buffer.data(li + 150);
    const auto *li_151 = buffer.data(li + 151);
    const auto *li_152 = buffer.data(li + 152);
    const auto *li_153 = buffer.data(li + 153);
    const auto *li_154 = buffer.data(li + 154);
    const auto *li_155 = buffer.data(li + 155);
    const auto *li_156 = buffer.data(li + 156);
    const auto *li_157 = buffer.data(li + 157);
    const auto *li_158 = buffer.data(li + 158);
    const auto *li_159 = buffer.data(li + 159);
    const auto *li_160 = buffer.data(li + 160);
    const auto *li_161 = buffer.data(li + 161);
    const auto *li_162 = buffer.data(li + 162);
    const auto *li_163 = buffer.data(li + 163);
    const auto *li_164 = buffer.data(li + 164);
    const auto *li_165 = buffer.data(li + 165);
    const auto *li_166 = buffer.data(li + 166);
    const auto *li_167 = buffer.data(li + 167);
    const auto *li_168 = buffer.data(li + 168);
    const auto *li_169 = buffer.data(li + 169);
    const auto *li_170 = buffer.data(li + 170);
    const auto *li_171 = buffer.data(li + 171);
    const auto *li_172 = buffer.data(li + 172);
    const auto *li_173 = buffer.data(li + 173);
    const auto *li_174 = buffer.data(li + 174);
    const auto *li_175 = buffer.data(li + 175);
    const auto *li_176 = buffer.data(li + 176);
    const auto *li_177 = buffer.data(li + 177);
    const auto *li_178 = buffer.data(li + 178);
    const auto *li_179 = buffer.data(li + 179);
    const auto *li_180 = buffer.data(li + 180);
    const auto *li_181 = buffer.data(li + 181);
    const auto *li_182 = buffer.data(li + 182);
    const auto *li_183 = buffer.data(li + 183);
    const auto *li_184 = buffer.data(li + 184);
    const auto *li_185 = buffer.data(li + 185);
    const auto *li_186 = buffer.data(li + 186);
    const auto *li_187 = buffer.data(li + 187);
    const auto *li_188 = buffer.data(li + 188);
    const auto *li_189 = buffer.data(li + 189);
    const auto *li_190 = buffer.data(li + 190);
    const auto *li_191 = buffer.data(li + 191);
    const auto *li_192 = buffer.data(li + 192);
    const auto *li_193 = buffer.data(li + 193);
    const auto *li_194 = buffer.data(li + 194);
    const auto *li_195 = buffer.data(li + 195);
    const auto *li_196 = buffer.data(li + 196);
    const auto *li_197 = buffer.data(li + 197);
    const auto *li_198 = buffer.data(li + 198);
    const auto *li_199 = buffer.data(li + 199);
    const auto *li_200 = buffer.data(li + 200);
    const auto *li_201 = buffer.data(li + 201);
    const auto *li_202 = buffer.data(li + 202);
    const auto *li_203 = buffer.data(li + 203);
    const auto *li_204 = buffer.data(li + 204);
    const auto *li_205 = buffer.data(li + 205);
    const auto *li_206 = buffer.data(li + 206);
    const auto *li_207 = buffer.data(li + 207);
    const auto *li_208 = buffer.data(li + 208);
    const auto *li_209 = buffer.data(li + 209);
    const auto *li_210 = buffer.data(li + 210);
    const auto *li_211 = buffer.data(li + 211);
    const auto *li_212 = buffer.data(li + 212);
    const auto *li_213 = buffer.data(li + 213);
    const auto *li_214 = buffer.data(li + 214);
    const auto *li_215 = buffer.data(li + 215);
    const auto *li_216 = buffer.data(li + 216);
    const auto *li_217 = buffer.data(li + 217);
    const auto *li_218 = buffer.data(li + 218);
    const auto *li_219 = buffer.data(li + 219);
    const auto *li_220 = buffer.data(li + 220);
    const auto *li_221 = buffer.data(li + 221);
    const auto *li_222 = buffer.data(li + 222);
    const auto *li_223 = buffer.data(li + 223);
    const auto *li_224 = buffer.data(li + 224);
    const auto *li_225 = buffer.data(li + 225);
    const auto *li_226 = buffer.data(li + 226);
    const auto *li_227 = buffer.data(li + 227);
    const auto *li_228 = buffer.data(li + 228);
    const auto *li_229 = buffer.data(li + 229);
    const auto *li_230 = buffer.data(li + 230);
    const auto *li_231 = buffer.data(li + 231);
    const auto *li_232 = buffer.data(li + 232);
    const auto *li_233 = buffer.data(li + 233);
    const auto *li_234 = buffer.data(li + 234);
    const auto *li_235 = buffer.data(li + 235);
    const auto *li_236 = buffer.data(li + 236);
    const auto *li_237 = buffer.data(li + 237);
    const auto *li_238 = buffer.data(li + 238);
    const auto *li_239 = buffer.data(li + 239);
    const auto *li_240 = buffer.data(li + 240);
    const auto *li_241 = buffer.data(li + 241);
    const auto *li_242 = buffer.data(li + 242);
    const auto *li_243 = buffer.data(li + 243);
    const auto *li_244 = buffer.data(li + 244);
    const auto *li_245 = buffer.data(li + 245);
    const auto *li_246 = buffer.data(li + 246);
    const auto *li_247 = buffer.data(li + 247);
    const auto *li_248 = buffer.data(li + 248);
    const auto *li_249 = buffer.data(li + 249);
    const auto *li_250 = buffer.data(li + 250);
    const auto *li_251 = buffer.data(li + 251);
    const auto *li_252 = buffer.data(li + 252);
    const auto *li_253 = buffer.data(li + 253);
    const auto *li_254 = buffer.data(li + 254);
    const auto *li_255 = buffer.data(li + 255);
    const auto *li_256 = buffer.data(li + 256);
    const auto *li_257 = buffer.data(li + 257);
    const auto *li_258 = buffer.data(li + 258);
    const auto *li_259 = buffer.data(li + 259);
    const auto *li_260 = buffer.data(li + 260);
    const auto *li_261 = buffer.data(li + 261);
    const auto *li_262 = buffer.data(li + 262);
    const auto *li_263 = buffer.data(li + 263);
    const auto *li_264 = buffer.data(li + 264);
    const auto *li_265 = buffer.data(li + 265);
    const auto *li_266 = buffer.data(li + 266);
    const auto *li_267 = buffer.data(li + 267);
    const auto *li_268 = buffer.data(li + 268);
    const auto *li_269 = buffer.data(li + 269);
    const auto *li_270 = buffer.data(li + 270);
    const auto *li_271 = buffer.data(li + 271);
    const auto *li_272 = buffer.data(li + 272);
    const auto *li_273 = buffer.data(li + 273);
    const auto *li_274 = buffer.data(li + 274);
    const auto *li_275 = buffer.data(li + 275);
    const auto *li_276 = buffer.data(li + 276);
    const auto *li_277 = buffer.data(li + 277);
    const auto *li_278 = buffer.data(li + 278);
    const auto *li_279 = buffer.data(li + 279);
    const auto *li_280 = buffer.data(li + 280);
    const auto *li_281 = buffer.data(li + 281);
    const auto *li_282 = buffer.data(li + 282);
    const auto *li_283 = buffer.data(li + 283);
    const auto *li_284 = buffer.data(li + 284);
    const auto *li_285 = buffer.data(li + 285);
    const auto *li_286 = buffer.data(li + 286);
    const auto *li_287 = buffer.data(li + 287);
    const auto *li_288 = buffer.data(li + 288);
    const auto *li_289 = buffer.data(li + 289);
    const auto *li_290 = buffer.data(li + 290);
    const auto *li_291 = buffer.data(li + 291);
    const auto *li_292 = buffer.data(li + 292);
    const auto *li_293 = buffer.data(li + 293);
    const auto *li_294 = buffer.data(li + 294);
    const auto *li_295 = buffer.data(li + 295);
    const auto *li_296 = buffer.data(li + 296);
    const auto *li_297 = buffer.data(li + 297);
    const auto *li_298 = buffer.data(li + 298);
    const auto *li_299 = buffer.data(li + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ii_150, ii_151, ii_152, ii_153, \
                         ii_154, li_150, li_151, li_152, li_153, \
                         li_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -5.0 * ii_150[k]
                   + f_0 * li_150[k];

        t_151[k] = -5.0 * ii_151[k]
                   + f_0 * li_151[k];

        t_152[k] = -5.0 * ii_152[k]
                   + f_0 * li_152[k];

        t_153[k] = -5.0 * ii_153[k]
                   + f_0 * li_153[k];

        t_154[k] = -5.0 * ii_154[k]
                   + f_0 * li_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ii_155, ii_156, ii_157, ii_158, \
                         ii_159, li_155, li_156, li_157, li_158, \
                         li_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -5.0 * ii_155[k]
                   + f_0 * li_155[k];

        t_156[k] = -5.0 * ii_156[k]
                   + f_0 * li_156[k];

        t_157[k] = -5.0 * ii_157[k]
                   + f_0 * li_157[k];

        t_158[k] = -5.0 * ii_158[k]
                   + f_0 * li_158[k];

        t_159[k] = -5.0 * ii_159[k]
                   + f_0 * li_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ii_160, ii_161, ii_162, ii_163, \
                         ii_164, li_160, li_161, li_162, li_163, \
                         li_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -5.0 * ii_160[k]
                   + f_0 * li_160[k];

        t_161[k] = -5.0 * ii_161[k]
                   + f_0 * li_161[k];

        t_162[k] = -5.0 * ii_162[k]
                   + f_0 * li_162[k];

        t_163[k] = -5.0 * ii_163[k]
                   + f_0 * li_163[k];

        t_164[k] = -5.0 * ii_164[k]
                   + f_0 * li_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ii_165, ii_166, ii_167, ii_168, \
                         ii_169, li_165, li_166, li_167, li_168, \
                         li_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -5.0 * ii_165[k]
                   + f_0 * li_165[k];

        t_166[k] = -5.0 * ii_166[k]
                   + f_0 * li_166[k];

        t_167[k] = -5.0 * ii_167[k]
                   + f_0 * li_167[k];

        t_168[k] = -4.0 * ii_168[k]
                   + f_0 * li_168[k];

        t_169[k] = -4.0 * ii_169[k]
                   + f_0 * li_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ii_170, ii_171, ii_172, ii_173, \
                         ii_174, li_170, li_171, li_172, li_173, \
                         li_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -4.0 * ii_170[k]
                   + f_0 * li_170[k];

        t_171[k] = -4.0 * ii_171[k]
                   + f_0 * li_171[k];

        t_172[k] = -4.0 * ii_172[k]
                   + f_0 * li_172[k];

        t_173[k] = -4.0 * ii_173[k]
                   + f_0 * li_173[k];

        t_174[k] = -4.0 * ii_174[k]
                   + f_0 * li_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ii_175, ii_176, ii_177, ii_178, \
                         ii_179, li_175, li_176, li_177, li_178, \
                         li_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -4.0 * ii_175[k]
                   + f_0 * li_175[k];

        t_176[k] = -4.0 * ii_176[k]
                   + f_0 * li_176[k];

        t_177[k] = -4.0 * ii_177[k]
                   + f_0 * li_177[k];

        t_178[k] = -4.0 * ii_178[k]
                   + f_0 * li_178[k];

        t_179[k] = -4.0 * ii_179[k]
                   + f_0 * li_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ii_180, ii_181, ii_182, ii_183, \
                         ii_184, li_180, li_181, li_182, li_183, \
                         li_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -4.0 * ii_180[k]
                   + f_0 * li_180[k];

        t_181[k] = -4.0 * ii_181[k]
                   + f_0 * li_181[k];

        t_182[k] = -4.0 * ii_182[k]
                   + f_0 * li_182[k];

        t_183[k] = -4.0 * ii_183[k]
                   + f_0 * li_183[k];

        t_184[k] = -4.0 * ii_184[k]
                   + f_0 * li_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ii_185, ii_186, ii_187, ii_188, \
                         ii_189, li_185, li_186, li_187, li_188, \
                         li_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -4.0 * ii_185[k]
                   + f_0 * li_185[k];

        t_186[k] = -4.0 * ii_186[k]
                   + f_0 * li_186[k];

        t_187[k] = -4.0 * ii_187[k]
                   + f_0 * li_187[k];

        t_188[k] = -4.0 * ii_188[k]
                   + f_0 * li_188[k];

        t_189[k] = -4.0 * ii_189[k]
                   + f_0 * li_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ii_190, ii_191, ii_192, ii_193, \
                         ii_194, li_190, li_191, li_192, li_193, \
                         li_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -4.0 * ii_190[k]
                   + f_0 * li_190[k];

        t_191[k] = -4.0 * ii_191[k]
                   + f_0 * li_191[k];

        t_192[k] = -4.0 * ii_192[k]
                   + f_0 * li_192[k];

        t_193[k] = -4.0 * ii_193[k]
                   + f_0 * li_193[k];

        t_194[k] = -4.0 * ii_194[k]
                   + f_0 * li_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ii_195, ii_196, ii_197, ii_198, \
                         ii_199, li_195, li_196, li_197, li_198, \
                         li_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -4.0 * ii_195[k]
                   + f_0 * li_195[k];

        t_196[k] = -4.0 * ii_196[k]
                   + f_0 * li_196[k];

        t_197[k] = -4.0 * ii_197[k]
                   + f_0 * li_197[k];

        t_198[k] = -4.0 * ii_198[k]
                   + f_0 * li_198[k];

        t_199[k] = -4.0 * ii_199[k]
                   + f_0 * li_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ii_200, ii_201, ii_202, ii_203, \
                         ii_204, li_200, li_201, li_202, li_203, \
                         li_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -4.0 * ii_200[k]
                   + f_0 * li_200[k];

        t_201[k] = -4.0 * ii_201[k]
                   + f_0 * li_201[k];

        t_202[k] = -4.0 * ii_202[k]
                   + f_0 * li_202[k];

        t_203[k] = -4.0 * ii_203[k]
                   + f_0 * li_203[k];

        t_204[k] = -4.0 * ii_204[k]
                   + f_0 * li_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ii_205, ii_206, ii_207, ii_208, \
                         ii_209, li_205, li_206, li_207, li_208, \
                         li_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -4.0 * ii_205[k]
                   + f_0 * li_205[k];

        t_206[k] = -4.0 * ii_206[k]
                   + f_0 * li_206[k];

        t_207[k] = -4.0 * ii_207[k]
                   + f_0 * li_207[k];

        t_208[k] = -4.0 * ii_208[k]
                   + f_0 * li_208[k];

        t_209[k] = -4.0 * ii_209[k]
                   + f_0 * li_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ii_210, ii_211, ii_212, ii_213, \
                         ii_214, li_210, li_211, li_212, li_213, \
                         li_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -4.0 * ii_210[k]
                   + f_0 * li_210[k];

        t_211[k] = -4.0 * ii_211[k]
                   + f_0 * li_211[k];

        t_212[k] = -4.0 * ii_212[k]
                   + f_0 * li_212[k];

        t_213[k] = -4.0 * ii_213[k]
                   + f_0 * li_213[k];

        t_214[k] = -4.0 * ii_214[k]
                   + f_0 * li_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ii_215, ii_216, ii_217, ii_218, \
                         ii_219, li_215, li_216, li_217, li_218, \
                         li_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -4.0 * ii_215[k]
                   + f_0 * li_215[k];

        t_216[k] = -4.0 * ii_216[k]
                   + f_0 * li_216[k];

        t_217[k] = -4.0 * ii_217[k]
                   + f_0 * li_217[k];

        t_218[k] = -4.0 * ii_218[k]
                   + f_0 * li_218[k];

        t_219[k] = -4.0 * ii_219[k]
                   + f_0 * li_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ii_220, ii_221, ii_222, ii_223, \
                         ii_224, li_220, li_221, li_222, li_223, \
                         li_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -4.0 * ii_220[k]
                   + f_0 * li_220[k];

        t_221[k] = -4.0 * ii_221[k]
                   + f_0 * li_221[k];

        t_222[k] = -4.0 * ii_222[k]
                   + f_0 * li_222[k];

        t_223[k] = -4.0 * ii_223[k]
                   + f_0 * li_223[k];

        t_224[k] = -4.0 * ii_224[k]
                   + f_0 * li_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ii_225, ii_226, ii_227, ii_228, \
                         ii_229, li_225, li_226, li_227, li_228, \
                         li_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -4.0 * ii_225[k]
                   + f_0 * li_225[k];

        t_226[k] = -4.0 * ii_226[k]
                   + f_0 * li_226[k];

        t_227[k] = -4.0 * ii_227[k]
                   + f_0 * li_227[k];

        t_228[k] = -4.0 * ii_228[k]
                   + f_0 * li_228[k];

        t_229[k] = -4.0 * ii_229[k]
                   + f_0 * li_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ii_230, ii_231, ii_232, ii_233, \
                         ii_234, li_230, li_231, li_232, li_233, \
                         li_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -4.0 * ii_230[k]
                   + f_0 * li_230[k];

        t_231[k] = -4.0 * ii_231[k]
                   + f_0 * li_231[k];

        t_232[k] = -4.0 * ii_232[k]
                   + f_0 * li_232[k];

        t_233[k] = -4.0 * ii_233[k]
                   + f_0 * li_233[k];

        t_234[k] = -4.0 * ii_234[k]
                   + f_0 * li_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ii_235, ii_236, ii_237, ii_238, \
                         ii_239, li_235, li_236, li_237, li_238, \
                         li_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -4.0 * ii_235[k]
                   + f_0 * li_235[k];

        t_236[k] = -4.0 * ii_236[k]
                   + f_0 * li_236[k];

        t_237[k] = -4.0 * ii_237[k]
                   + f_0 * li_237[k];

        t_238[k] = -4.0 * ii_238[k]
                   + f_0 * li_238[k];

        t_239[k] = -4.0 * ii_239[k]
                   + f_0 * li_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ii_240, ii_241, ii_242, ii_243, \
                         ii_244, li_240, li_241, li_242, li_243, \
                         li_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -4.0 * ii_240[k]
                   + f_0 * li_240[k];

        t_241[k] = -4.0 * ii_241[k]
                   + f_0 * li_241[k];

        t_242[k] = -4.0 * ii_242[k]
                   + f_0 * li_242[k];

        t_243[k] = -4.0 * ii_243[k]
                   + f_0 * li_243[k];

        t_244[k] = -4.0 * ii_244[k]
                   + f_0 * li_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ii_245, ii_246, ii_247, ii_248, \
                         ii_249, li_245, li_246, li_247, li_248, \
                         li_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -4.0 * ii_245[k]
                   + f_0 * li_245[k];

        t_246[k] = -4.0 * ii_246[k]
                   + f_0 * li_246[k];

        t_247[k] = -4.0 * ii_247[k]
                   + f_0 * li_247[k];

        t_248[k] = -4.0 * ii_248[k]
                   + f_0 * li_248[k];

        t_249[k] = -4.0 * ii_249[k]
                   + f_0 * li_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ii_250, ii_251, ii_252, ii_253, \
                         ii_254, li_250, li_251, li_252, li_253, \
                         li_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -4.0 * ii_250[k]
                   + f_0 * li_250[k];

        t_251[k] = -4.0 * ii_251[k]
                   + f_0 * li_251[k];

        t_252[k] = -4.0 * ii_252[k]
                   + f_0 * li_252[k];

        t_253[k] = -4.0 * ii_253[k]
                   + f_0 * li_253[k];

        t_254[k] = -4.0 * ii_254[k]
                   + f_0 * li_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ii_255, ii_256, ii_257, ii_258, \
                         ii_259, li_255, li_256, li_257, li_258, \
                         li_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -4.0 * ii_255[k]
                   + f_0 * li_255[k];

        t_256[k] = -4.0 * ii_256[k]
                   + f_0 * li_256[k];

        t_257[k] = -4.0 * ii_257[k]
                   + f_0 * li_257[k];

        t_258[k] = -4.0 * ii_258[k]
                   + f_0 * li_258[k];

        t_259[k] = -4.0 * ii_259[k]
                   + f_0 * li_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ii_260, ii_261, ii_262, ii_263, \
                         ii_264, li_260, li_261, li_262, li_263, \
                         li_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -4.0 * ii_260[k]
                   + f_0 * li_260[k];

        t_261[k] = -4.0 * ii_261[k]
                   + f_0 * li_261[k];

        t_262[k] = -4.0 * ii_262[k]
                   + f_0 * li_262[k];

        t_263[k] = -4.0 * ii_263[k]
                   + f_0 * li_263[k];

        t_264[k] = -4.0 * ii_264[k]
                   + f_0 * li_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ii_265, ii_266, ii_267, ii_268, \
                         ii_269, li_265, li_266, li_267, li_268, \
                         li_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -4.0 * ii_265[k]
                   + f_0 * li_265[k];

        t_266[k] = -4.0 * ii_266[k]
                   + f_0 * li_266[k];

        t_267[k] = -4.0 * ii_267[k]
                   + f_0 * li_267[k];

        t_268[k] = -4.0 * ii_268[k]
                   + f_0 * li_268[k];

        t_269[k] = -4.0 * ii_269[k]
                   + f_0 * li_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ii_270, ii_271, ii_272, ii_273, \
                         ii_274, li_270, li_271, li_272, li_273, \
                         li_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -4.0 * ii_270[k]
                   + f_0 * li_270[k];

        t_271[k] = -4.0 * ii_271[k]
                   + f_0 * li_271[k];

        t_272[k] = -4.0 * ii_272[k]
                   + f_0 * li_272[k];

        t_273[k] = -4.0 * ii_273[k]
                   + f_0 * li_273[k];

        t_274[k] = -4.0 * ii_274[k]
                   + f_0 * li_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ii_275, ii_276, ii_277, ii_278, \
                         ii_279, li_275, li_276, li_277, li_278, \
                         li_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -4.0 * ii_275[k]
                   + f_0 * li_275[k];

        t_276[k] = -4.0 * ii_276[k]
                   + f_0 * li_276[k];

        t_277[k] = -4.0 * ii_277[k]
                   + f_0 * li_277[k];

        t_278[k] = -4.0 * ii_278[k]
                   + f_0 * li_278[k];

        t_279[k] = -4.0 * ii_279[k]
                   + f_0 * li_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ii_280, ii_281, ii_282, ii_283, \
                         ii_284, li_280, li_281, li_282, li_283, \
                         li_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * ii_280[k]
                   + f_0 * li_280[k];

        t_281[k] = -3.0 * ii_281[k]
                   + f_0 * li_281[k];

        t_282[k] = -3.0 * ii_282[k]
                   + f_0 * li_282[k];

        t_283[k] = -3.0 * ii_283[k]
                   + f_0 * li_283[k];

        t_284[k] = -3.0 * ii_284[k]
                   + f_0 * li_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ii_285, ii_286, ii_287, ii_288, \
                         ii_289, li_285, li_286, li_287, li_288, \
                         li_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * ii_285[k]
                   + f_0 * li_285[k];

        t_286[k] = -3.0 * ii_286[k]
                   + f_0 * li_286[k];

        t_287[k] = -3.0 * ii_287[k]
                   + f_0 * li_287[k];

        t_288[k] = -3.0 * ii_288[k]
                   + f_0 * li_288[k];

        t_289[k] = -3.0 * ii_289[k]
                   + f_0 * li_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ii_290, ii_291, ii_292, ii_293, \
                         ii_294, li_290, li_291, li_292, li_293, \
                         li_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * ii_290[k]
                   + f_0 * li_290[k];

        t_291[k] = -3.0 * ii_291[k]
                   + f_0 * li_291[k];

        t_292[k] = -3.0 * ii_292[k]
                   + f_0 * li_292[k];

        t_293[k] = -3.0 * ii_293[k]
                   + f_0 * li_293[k];

        t_294[k] = -3.0 * ii_294[k]
                   + f_0 * li_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ii_295, ii_296, ii_297, ii_298, \
                         ii_299, li_295, li_296, li_297, li_298, \
                         li_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -3.0 * ii_295[k]
                   + f_0 * li_295[k];

        t_296[k] = -3.0 * ii_296[k]
                   + f_0 * li_296[k];

        t_297[k] = -3.0 * ii_297[k]
                   + f_0 * li_297[k];

        t_298[k] = -3.0 * ii_298[k]
                   + f_0 * li_298[k];

        t_299[k] = -3.0 * ii_299[k]
                   + f_0 * li_299[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_300 = buffer.data(ii + 300);
    const auto *ii_301 = buffer.data(ii + 301);
    const auto *ii_302 = buffer.data(ii + 302);
    const auto *ii_303 = buffer.data(ii + 303);
    const auto *ii_304 = buffer.data(ii + 304);
    const auto *ii_305 = buffer.data(ii + 305);
    const auto *ii_306 = buffer.data(ii + 306);
    const auto *ii_307 = buffer.data(ii + 307);
    const auto *ii_308 = buffer.data(ii + 308);
    const auto *ii_309 = buffer.data(ii + 309);
    const auto *ii_310 = buffer.data(ii + 310);
    const auto *ii_311 = buffer.data(ii + 311);
    const auto *ii_312 = buffer.data(ii + 312);
    const auto *ii_313 = buffer.data(ii + 313);
    const auto *ii_314 = buffer.data(ii + 314);
    const auto *ii_315 = buffer.data(ii + 315);
    const auto *ii_316 = buffer.data(ii + 316);
    const auto *ii_317 = buffer.data(ii + 317);
    const auto *ii_318 = buffer.data(ii + 318);
    const auto *ii_319 = buffer.data(ii + 319);
    const auto *ii_320 = buffer.data(ii + 320);
    const auto *ii_321 = buffer.data(ii + 321);
    const auto *ii_322 = buffer.data(ii + 322);
    const auto *ii_323 = buffer.data(ii + 323);
    const auto *ii_324 = buffer.data(ii + 324);
    const auto *ii_325 = buffer.data(ii + 325);
    const auto *ii_326 = buffer.data(ii + 326);
    const auto *ii_327 = buffer.data(ii + 327);
    const auto *ii_328 = buffer.data(ii + 328);
    const auto *ii_329 = buffer.data(ii + 329);
    const auto *ii_330 = buffer.data(ii + 330);
    const auto *ii_331 = buffer.data(ii + 331);
    const auto *ii_332 = buffer.data(ii + 332);
    const auto *ii_333 = buffer.data(ii + 333);
    const auto *ii_334 = buffer.data(ii + 334);
    const auto *ii_335 = buffer.data(ii + 335);
    const auto *ii_336 = buffer.data(ii + 336);
    const auto *ii_337 = buffer.data(ii + 337);
    const auto *ii_338 = buffer.data(ii + 338);
    const auto *ii_339 = buffer.data(ii + 339);
    const auto *ii_340 = buffer.data(ii + 340);
    const auto *ii_341 = buffer.data(ii + 341);
    const auto *ii_342 = buffer.data(ii + 342);
    const auto *ii_343 = buffer.data(ii + 343);
    const auto *ii_344 = buffer.data(ii + 344);
    const auto *ii_345 = buffer.data(ii + 345);
    const auto *ii_346 = buffer.data(ii + 346);
    const auto *ii_347 = buffer.data(ii + 347);
    const auto *ii_348 = buffer.data(ii + 348);
    const auto *ii_349 = buffer.data(ii + 349);
    const auto *ii_350 = buffer.data(ii + 350);
    const auto *ii_351 = buffer.data(ii + 351);
    const auto *ii_352 = buffer.data(ii + 352);
    const auto *ii_353 = buffer.data(ii + 353);
    const auto *ii_354 = buffer.data(ii + 354);
    const auto *ii_355 = buffer.data(ii + 355);
    const auto *ii_356 = buffer.data(ii + 356);
    const auto *ii_357 = buffer.data(ii + 357);
    const auto *ii_358 = buffer.data(ii + 358);
    const auto *ii_359 = buffer.data(ii + 359);
    const auto *ii_360 = buffer.data(ii + 360);
    const auto *ii_361 = buffer.data(ii + 361);
    const auto *ii_362 = buffer.data(ii + 362);
    const auto *ii_363 = buffer.data(ii + 363);
    const auto *ii_364 = buffer.data(ii + 364);
    const auto *ii_365 = buffer.data(ii + 365);
    const auto *ii_366 = buffer.data(ii + 366);
    const auto *ii_367 = buffer.data(ii + 367);
    const auto *ii_368 = buffer.data(ii + 368);
    const auto *ii_369 = buffer.data(ii + 369);
    const auto *ii_370 = buffer.data(ii + 370);
    const auto *ii_371 = buffer.data(ii + 371);
    const auto *ii_372 = buffer.data(ii + 372);
    const auto *ii_373 = buffer.data(ii + 373);
    const auto *ii_374 = buffer.data(ii + 374);
    const auto *ii_375 = buffer.data(ii + 375);
    const auto *ii_376 = buffer.data(ii + 376);
    const auto *ii_377 = buffer.data(ii + 377);
    const auto *ii_378 = buffer.data(ii + 378);
    const auto *ii_379 = buffer.data(ii + 379);
    const auto *ii_380 = buffer.data(ii + 380);
    const auto *ii_381 = buffer.data(ii + 381);
    const auto *ii_382 = buffer.data(ii + 382);
    const auto *ii_383 = buffer.data(ii + 383);
    const auto *ii_384 = buffer.data(ii + 384);
    const auto *ii_385 = buffer.data(ii + 385);
    const auto *ii_386 = buffer.data(ii + 386);
    const auto *ii_387 = buffer.data(ii + 387);
    const auto *ii_388 = buffer.data(ii + 388);
    const auto *ii_389 = buffer.data(ii + 389);
    const auto *ii_390 = buffer.data(ii + 390);
    const auto *ii_391 = buffer.data(ii + 391);
    const auto *ii_392 = buffer.data(ii + 392);
    const auto *ii_393 = buffer.data(ii + 393);
    const auto *ii_394 = buffer.data(ii + 394);
    const auto *ii_395 = buffer.data(ii + 395);
    const auto *ii_396 = buffer.data(ii + 396);
    const auto *ii_397 = buffer.data(ii + 397);
    const auto *ii_398 = buffer.data(ii + 398);
    const auto *ii_399 = buffer.data(ii + 399);
    const auto *ii_400 = buffer.data(ii + 400);
    const auto *ii_401 = buffer.data(ii + 401);
    const auto *ii_402 = buffer.data(ii + 402);
    const auto *ii_403 = buffer.data(ii + 403);
    const auto *ii_404 = buffer.data(ii + 404);
    const auto *ii_405 = buffer.data(ii + 405);
    const auto *ii_406 = buffer.data(ii + 406);
    const auto *ii_407 = buffer.data(ii + 407);
    const auto *ii_408 = buffer.data(ii + 408);
    const auto *ii_409 = buffer.data(ii + 409);
    const auto *ii_410 = buffer.data(ii + 410);
    const auto *ii_411 = buffer.data(ii + 411);
    const auto *ii_412 = buffer.data(ii + 412);
    const auto *ii_413 = buffer.data(ii + 413);
    const auto *ii_414 = buffer.data(ii + 414);
    const auto *ii_415 = buffer.data(ii + 415);
    const auto *ii_416 = buffer.data(ii + 416);
    const auto *ii_417 = buffer.data(ii + 417);
    const auto *ii_418 = buffer.data(ii + 418);
    const auto *ii_419 = buffer.data(ii + 419);
    const auto *ii_420 = buffer.data(ii + 420);
    const auto *ii_421 = buffer.data(ii + 421);
    const auto *ii_422 = buffer.data(ii + 422);
    const auto *ii_423 = buffer.data(ii + 423);
    const auto *ii_424 = buffer.data(ii + 424);
    const auto *ii_425 = buffer.data(ii + 425);
    const auto *ii_426 = buffer.data(ii + 426);
    const auto *ii_427 = buffer.data(ii + 427);
    const auto *ii_428 = buffer.data(ii + 428);
    const auto *ii_429 = buffer.data(ii + 429);
    const auto *ii_430 = buffer.data(ii + 430);
    const auto *ii_431 = buffer.data(ii + 431);
    const auto *ii_432 = buffer.data(ii + 432);
    const auto *ii_433 = buffer.data(ii + 433);
    const auto *ii_434 = buffer.data(ii + 434);
    const auto *ii_435 = buffer.data(ii + 435);
    const auto *ii_436 = buffer.data(ii + 436);
    const auto *ii_437 = buffer.data(ii + 437);
    const auto *ii_438 = buffer.data(ii + 438);
    const auto *ii_439 = buffer.data(ii + 439);
    const auto *ii_440 = buffer.data(ii + 440);
    const auto *ii_441 = buffer.data(ii + 441);
    const auto *ii_442 = buffer.data(ii + 442);
    const auto *ii_443 = buffer.data(ii + 443);
    const auto *ii_444 = buffer.data(ii + 444);
    const auto *ii_445 = buffer.data(ii + 445);
    const auto *ii_446 = buffer.data(ii + 446);
    const auto *ii_447 = buffer.data(ii + 447);
    const auto *ii_448 = buffer.data(ii + 448);
    const auto *ii_449 = buffer.data(ii + 449);

    const auto *li_300 = buffer.data(li + 300);
    const auto *li_301 = buffer.data(li + 301);
    const auto *li_302 = buffer.data(li + 302);
    const auto *li_303 = buffer.data(li + 303);
    const auto *li_304 = buffer.data(li + 304);
    const auto *li_305 = buffer.data(li + 305);
    const auto *li_306 = buffer.data(li + 306);
    const auto *li_307 = buffer.data(li + 307);
    const auto *li_308 = buffer.data(li + 308);
    const auto *li_309 = buffer.data(li + 309);
    const auto *li_310 = buffer.data(li + 310);
    const auto *li_311 = buffer.data(li + 311);
    const auto *li_312 = buffer.data(li + 312);
    const auto *li_313 = buffer.data(li + 313);
    const auto *li_314 = buffer.data(li + 314);
    const auto *li_315 = buffer.data(li + 315);
    const auto *li_316 = buffer.data(li + 316);
    const auto *li_317 = buffer.data(li + 317);
    const auto *li_318 = buffer.data(li + 318);
    const auto *li_319 = buffer.data(li + 319);
    const auto *li_320 = buffer.data(li + 320);
    const auto *li_321 = buffer.data(li + 321);
    const auto *li_322 = buffer.data(li + 322);
    const auto *li_323 = buffer.data(li + 323);
    const auto *li_324 = buffer.data(li + 324);
    const auto *li_325 = buffer.data(li + 325);
    const auto *li_326 = buffer.data(li + 326);
    const auto *li_327 = buffer.data(li + 327);
    const auto *li_328 = buffer.data(li + 328);
    const auto *li_329 = buffer.data(li + 329);
    const auto *li_330 = buffer.data(li + 330);
    const auto *li_331 = buffer.data(li + 331);
    const auto *li_332 = buffer.data(li + 332);
    const auto *li_333 = buffer.data(li + 333);
    const auto *li_334 = buffer.data(li + 334);
    const auto *li_335 = buffer.data(li + 335);
    const auto *li_336 = buffer.data(li + 336);
    const auto *li_337 = buffer.data(li + 337);
    const auto *li_338 = buffer.data(li + 338);
    const auto *li_339 = buffer.data(li + 339);
    const auto *li_340 = buffer.data(li + 340);
    const auto *li_341 = buffer.data(li + 341);
    const auto *li_342 = buffer.data(li + 342);
    const auto *li_343 = buffer.data(li + 343);
    const auto *li_344 = buffer.data(li + 344);
    const auto *li_345 = buffer.data(li + 345);
    const auto *li_346 = buffer.data(li + 346);
    const auto *li_347 = buffer.data(li + 347);
    const auto *li_348 = buffer.data(li + 348);
    const auto *li_349 = buffer.data(li + 349);
    const auto *li_350 = buffer.data(li + 350);
    const auto *li_351 = buffer.data(li + 351);
    const auto *li_352 = buffer.data(li + 352);
    const auto *li_353 = buffer.data(li + 353);
    const auto *li_354 = buffer.data(li + 354);
    const auto *li_355 = buffer.data(li + 355);
    const auto *li_356 = buffer.data(li + 356);
    const auto *li_357 = buffer.data(li + 357);
    const auto *li_358 = buffer.data(li + 358);
    const auto *li_359 = buffer.data(li + 359);
    const auto *li_360 = buffer.data(li + 360);
    const auto *li_361 = buffer.data(li + 361);
    const auto *li_362 = buffer.data(li + 362);
    const auto *li_363 = buffer.data(li + 363);
    const auto *li_364 = buffer.data(li + 364);
    const auto *li_365 = buffer.data(li + 365);
    const auto *li_366 = buffer.data(li + 366);
    const auto *li_367 = buffer.data(li + 367);
    const auto *li_368 = buffer.data(li + 368);
    const auto *li_369 = buffer.data(li + 369);
    const auto *li_370 = buffer.data(li + 370);
    const auto *li_371 = buffer.data(li + 371);
    const auto *li_372 = buffer.data(li + 372);
    const auto *li_373 = buffer.data(li + 373);
    const auto *li_374 = buffer.data(li + 374);
    const auto *li_375 = buffer.data(li + 375);
    const auto *li_376 = buffer.data(li + 376);
    const auto *li_377 = buffer.data(li + 377);
    const auto *li_378 = buffer.data(li + 378);
    const auto *li_379 = buffer.data(li + 379);
    const auto *li_380 = buffer.data(li + 380);
    const auto *li_381 = buffer.data(li + 381);
    const auto *li_382 = buffer.data(li + 382);
    const auto *li_383 = buffer.data(li + 383);
    const auto *li_384 = buffer.data(li + 384);
    const auto *li_385 = buffer.data(li + 385);
    const auto *li_386 = buffer.data(li + 386);
    const auto *li_387 = buffer.data(li + 387);
    const auto *li_388 = buffer.data(li + 388);
    const auto *li_389 = buffer.data(li + 389);
    const auto *li_390 = buffer.data(li + 390);
    const auto *li_391 = buffer.data(li + 391);
    const auto *li_392 = buffer.data(li + 392);
    const auto *li_393 = buffer.data(li + 393);
    const auto *li_394 = buffer.data(li + 394);
    const auto *li_395 = buffer.data(li + 395);
    const auto *li_396 = buffer.data(li + 396);
    const auto *li_397 = buffer.data(li + 397);
    const auto *li_398 = buffer.data(li + 398);
    const auto *li_399 = buffer.data(li + 399);
    const auto *li_400 = buffer.data(li + 400);
    const auto *li_401 = buffer.data(li + 401);
    const auto *li_402 = buffer.data(li + 402);
    const auto *li_403 = buffer.data(li + 403);
    const auto *li_404 = buffer.data(li + 404);
    const auto *li_405 = buffer.data(li + 405);
    const auto *li_406 = buffer.data(li + 406);
    const auto *li_407 = buffer.data(li + 407);
    const auto *li_408 = buffer.data(li + 408);
    const auto *li_409 = buffer.data(li + 409);
    const auto *li_410 = buffer.data(li + 410);
    const auto *li_411 = buffer.data(li + 411);
    const auto *li_412 = buffer.data(li + 412);
    const auto *li_413 = buffer.data(li + 413);
    const auto *li_414 = buffer.data(li + 414);
    const auto *li_415 = buffer.data(li + 415);
    const auto *li_416 = buffer.data(li + 416);
    const auto *li_417 = buffer.data(li + 417);
    const auto *li_418 = buffer.data(li + 418);
    const auto *li_419 = buffer.data(li + 419);
    const auto *li_420 = buffer.data(li + 420);
    const auto *li_421 = buffer.data(li + 421);
    const auto *li_422 = buffer.data(li + 422);
    const auto *li_423 = buffer.data(li + 423);
    const auto *li_424 = buffer.data(li + 424);
    const auto *li_425 = buffer.data(li + 425);
    const auto *li_426 = buffer.data(li + 426);
    const auto *li_427 = buffer.data(li + 427);
    const auto *li_428 = buffer.data(li + 428);
    const auto *li_429 = buffer.data(li + 429);
    const auto *li_430 = buffer.data(li + 430);
    const auto *li_431 = buffer.data(li + 431);
    const auto *li_432 = buffer.data(li + 432);
    const auto *li_433 = buffer.data(li + 433);
    const auto *li_434 = buffer.data(li + 434);
    const auto *li_435 = buffer.data(li + 435);
    const auto *li_436 = buffer.data(li + 436);
    const auto *li_437 = buffer.data(li + 437);
    const auto *li_438 = buffer.data(li + 438);
    const auto *li_439 = buffer.data(li + 439);
    const auto *li_440 = buffer.data(li + 440);
    const auto *li_441 = buffer.data(li + 441);
    const auto *li_442 = buffer.data(li + 442);
    const auto *li_443 = buffer.data(li + 443);
    const auto *li_444 = buffer.data(li + 444);
    const auto *li_445 = buffer.data(li + 445);
    const auto *li_446 = buffer.data(li + 446);
    const auto *li_447 = buffer.data(li + 447);
    const auto *li_448 = buffer.data(li + 448);
    const auto *li_449 = buffer.data(li + 449);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ii_300, ii_301, ii_302, ii_303, \
                         ii_304, li_300, li_301, li_302, li_303, \
                         li_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -3.0 * ii_300[k]
                   + f_0 * li_300[k];

        t_301[k] = -3.0 * ii_301[k]
                   + f_0 * li_301[k];

        t_302[k] = -3.0 * ii_302[k]
                   + f_0 * li_302[k];

        t_303[k] = -3.0 * ii_303[k]
                   + f_0 * li_303[k];

        t_304[k] = -3.0 * ii_304[k]
                   + f_0 * li_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ii_305, ii_306, ii_307, ii_308, \
                         ii_309, li_305, li_306, li_307, li_308, \
                         li_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -3.0 * ii_305[k]
                   + f_0 * li_305[k];

        t_306[k] = -3.0 * ii_306[k]
                   + f_0 * li_306[k];

        t_307[k] = -3.0 * ii_307[k]
                   + f_0 * li_307[k];

        t_308[k] = -3.0 * ii_308[k]
                   + f_0 * li_308[k];

        t_309[k] = -3.0 * ii_309[k]
                   + f_0 * li_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ii_310, ii_311, ii_312, ii_313, \
                         ii_314, li_310, li_311, li_312, li_313, \
                         li_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * ii_310[k]
                   + f_0 * li_310[k];

        t_311[k] = -3.0 * ii_311[k]
                   + f_0 * li_311[k];

        t_312[k] = -3.0 * ii_312[k]
                   + f_0 * li_312[k];

        t_313[k] = -3.0 * ii_313[k]
                   + f_0 * li_313[k];

        t_314[k] = -3.0 * ii_314[k]
                   + f_0 * li_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ii_315, ii_316, ii_317, ii_318, \
                         ii_319, li_315, li_316, li_317, li_318, \
                         li_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -3.0 * ii_315[k]
                   + f_0 * li_315[k];

        t_316[k] = -3.0 * ii_316[k]
                   + f_0 * li_316[k];

        t_317[k] = -3.0 * ii_317[k]
                   + f_0 * li_317[k];

        t_318[k] = -3.0 * ii_318[k]
                   + f_0 * li_318[k];

        t_319[k] = -3.0 * ii_319[k]
                   + f_0 * li_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ii_320, ii_321, ii_322, ii_323, \
                         ii_324, li_320, li_321, li_322, li_323, \
                         li_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -3.0 * ii_320[k]
                   + f_0 * li_320[k];

        t_321[k] = -3.0 * ii_321[k]
                   + f_0 * li_321[k];

        t_322[k] = -3.0 * ii_322[k]
                   + f_0 * li_322[k];

        t_323[k] = -3.0 * ii_323[k]
                   + f_0 * li_323[k];

        t_324[k] = -3.0 * ii_324[k]
                   + f_0 * li_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ii_325, ii_326, ii_327, ii_328, \
                         ii_329, li_325, li_326, li_327, li_328, \
                         li_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -3.0 * ii_325[k]
                   + f_0 * li_325[k];

        t_326[k] = -3.0 * ii_326[k]
                   + f_0 * li_326[k];

        t_327[k] = -3.0 * ii_327[k]
                   + f_0 * li_327[k];

        t_328[k] = -3.0 * ii_328[k]
                   + f_0 * li_328[k];

        t_329[k] = -3.0 * ii_329[k]
                   + f_0 * li_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ii_330, ii_331, ii_332, ii_333, \
                         ii_334, li_330, li_331, li_332, li_333, \
                         li_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -3.0 * ii_330[k]
                   + f_0 * li_330[k];

        t_331[k] = -3.0 * ii_331[k]
                   + f_0 * li_331[k];

        t_332[k] = -3.0 * ii_332[k]
                   + f_0 * li_332[k];

        t_333[k] = -3.0 * ii_333[k]
                   + f_0 * li_333[k];

        t_334[k] = -3.0 * ii_334[k]
                   + f_0 * li_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ii_335, ii_336, ii_337, ii_338, \
                         ii_339, li_335, li_336, li_337, li_338, \
                         li_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -3.0 * ii_335[k]
                   + f_0 * li_335[k];

        t_336[k] = -3.0 * ii_336[k]
                   + f_0 * li_336[k];

        t_337[k] = -3.0 * ii_337[k]
                   + f_0 * li_337[k];

        t_338[k] = -3.0 * ii_338[k]
                   + f_0 * li_338[k];

        t_339[k] = -3.0 * ii_339[k]
                   + f_0 * li_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ii_340, ii_341, ii_342, ii_343, \
                         ii_344, li_340, li_341, li_342, li_343, \
                         li_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -3.0 * ii_340[k]
                   + f_0 * li_340[k];

        t_341[k] = -3.0 * ii_341[k]
                   + f_0 * li_341[k];

        t_342[k] = -3.0 * ii_342[k]
                   + f_0 * li_342[k];

        t_343[k] = -3.0 * ii_343[k]
                   + f_0 * li_343[k];

        t_344[k] = -3.0 * ii_344[k]
                   + f_0 * li_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ii_345, ii_346, ii_347, ii_348, \
                         ii_349, li_345, li_346, li_347, li_348, \
                         li_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -3.0 * ii_345[k]
                   + f_0 * li_345[k];

        t_346[k] = -3.0 * ii_346[k]
                   + f_0 * li_346[k];

        t_347[k] = -3.0 * ii_347[k]
                   + f_0 * li_347[k];

        t_348[k] = -3.0 * ii_348[k]
                   + f_0 * li_348[k];

        t_349[k] = -3.0 * ii_349[k]
                   + f_0 * li_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ii_350, ii_351, ii_352, ii_353, \
                         ii_354, li_350, li_351, li_352, li_353, \
                         li_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -3.0 * ii_350[k]
                   + f_0 * li_350[k];

        t_351[k] = -3.0 * ii_351[k]
                   + f_0 * li_351[k];

        t_352[k] = -3.0 * ii_352[k]
                   + f_0 * li_352[k];

        t_353[k] = -3.0 * ii_353[k]
                   + f_0 * li_353[k];

        t_354[k] = -3.0 * ii_354[k]
                   + f_0 * li_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ii_355, ii_356, ii_357, ii_358, \
                         ii_359, li_355, li_356, li_357, li_358, \
                         li_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -3.0 * ii_355[k]
                   + f_0 * li_355[k];

        t_356[k] = -3.0 * ii_356[k]
                   + f_0 * li_356[k];

        t_357[k] = -3.0 * ii_357[k]
                   + f_0 * li_357[k];

        t_358[k] = -3.0 * ii_358[k]
                   + f_0 * li_358[k];

        t_359[k] = -3.0 * ii_359[k]
                   + f_0 * li_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ii_360, ii_361, ii_362, ii_363, \
                         ii_364, li_360, li_361, li_362, li_363, \
                         li_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * ii_360[k]
                   + f_0 * li_360[k];

        t_361[k] = -3.0 * ii_361[k]
                   + f_0 * li_361[k];

        t_362[k] = -3.0 * ii_362[k]
                   + f_0 * li_362[k];

        t_363[k] = -3.0 * ii_363[k]
                   + f_0 * li_363[k];

        t_364[k] = -3.0 * ii_364[k]
                   + f_0 * li_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ii_365, ii_366, ii_367, ii_368, \
                         ii_369, li_365, li_366, li_367, li_368, \
                         li_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * ii_365[k]
                   + f_0 * li_365[k];

        t_366[k] = -3.0 * ii_366[k]
                   + f_0 * li_366[k];

        t_367[k] = -3.0 * ii_367[k]
                   + f_0 * li_367[k];

        t_368[k] = -3.0 * ii_368[k]
                   + f_0 * li_368[k];

        t_369[k] = -3.0 * ii_369[k]
                   + f_0 * li_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ii_370, ii_371, ii_372, ii_373, \
                         ii_374, li_370, li_371, li_372, li_373, \
                         li_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * ii_370[k]
                   + f_0 * li_370[k];

        t_371[k] = -3.0 * ii_371[k]
                   + f_0 * li_371[k];

        t_372[k] = -3.0 * ii_372[k]
                   + f_0 * li_372[k];

        t_373[k] = -3.0 * ii_373[k]
                   + f_0 * li_373[k];

        t_374[k] = -3.0 * ii_374[k]
                   + f_0 * li_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ii_375, ii_376, ii_377, ii_378, \
                         ii_379, li_375, li_376, li_377, li_378, \
                         li_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -3.0 * ii_375[k]
                   + f_0 * li_375[k];

        t_376[k] = -3.0 * ii_376[k]
                   + f_0 * li_376[k];

        t_377[k] = -3.0 * ii_377[k]
                   + f_0 * li_377[k];

        t_378[k] = -3.0 * ii_378[k]
                   + f_0 * li_378[k];

        t_379[k] = -3.0 * ii_379[k]
                   + f_0 * li_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ii_380, ii_381, ii_382, ii_383, \
                         ii_384, li_380, li_381, li_382, li_383, \
                         li_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -3.0 * ii_380[k]
                   + f_0 * li_380[k];

        t_381[k] = -3.0 * ii_381[k]
                   + f_0 * li_381[k];

        t_382[k] = -3.0 * ii_382[k]
                   + f_0 * li_382[k];

        t_383[k] = -3.0 * ii_383[k]
                   + f_0 * li_383[k];

        t_384[k] = -3.0 * ii_384[k]
                   + f_0 * li_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ii_385, ii_386, ii_387, ii_388, \
                         ii_389, li_385, li_386, li_387, li_388, \
                         li_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -3.0 * ii_385[k]
                   + f_0 * li_385[k];

        t_386[k] = -3.0 * ii_386[k]
                   + f_0 * li_386[k];

        t_387[k] = -3.0 * ii_387[k]
                   + f_0 * li_387[k];

        t_388[k] = -3.0 * ii_388[k]
                   + f_0 * li_388[k];

        t_389[k] = -3.0 * ii_389[k]
                   + f_0 * li_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ii_390, ii_391, ii_392, ii_393, \
                         ii_394, li_390, li_391, li_392, li_393, \
                         li_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -3.0 * ii_390[k]
                   + f_0 * li_390[k];

        t_391[k] = -3.0 * ii_391[k]
                   + f_0 * li_391[k];

        t_392[k] = -3.0 * ii_392[k]
                   + f_0 * li_392[k];

        t_393[k] = -3.0 * ii_393[k]
                   + f_0 * li_393[k];

        t_394[k] = -3.0 * ii_394[k]
                   + f_0 * li_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ii_395, ii_396, ii_397, ii_398, \
                         ii_399, li_395, li_396, li_397, li_398, \
                         li_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -3.0 * ii_395[k]
                   + f_0 * li_395[k];

        t_396[k] = -3.0 * ii_396[k]
                   + f_0 * li_396[k];

        t_397[k] = -3.0 * ii_397[k]
                   + f_0 * li_397[k];

        t_398[k] = -3.0 * ii_398[k]
                   + f_0 * li_398[k];

        t_399[k] = -3.0 * ii_399[k]
                   + f_0 * li_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ii_400, ii_401, ii_402, ii_403, \
                         ii_404, li_400, li_401, li_402, li_403, \
                         li_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -3.0 * ii_400[k]
                   + f_0 * li_400[k];

        t_401[k] = -3.0 * ii_401[k]
                   + f_0 * li_401[k];

        t_402[k] = -3.0 * ii_402[k]
                   + f_0 * li_402[k];

        t_403[k] = -3.0 * ii_403[k]
                   + f_0 * li_403[k];

        t_404[k] = -3.0 * ii_404[k]
                   + f_0 * li_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ii_405, ii_406, ii_407, ii_408, \
                         ii_409, li_405, li_406, li_407, li_408, \
                         li_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -3.0 * ii_405[k]
                   + f_0 * li_405[k];

        t_406[k] = -3.0 * ii_406[k]
                   + f_0 * li_406[k];

        t_407[k] = -3.0 * ii_407[k]
                   + f_0 * li_407[k];

        t_408[k] = -3.0 * ii_408[k]
                   + f_0 * li_408[k];

        t_409[k] = -3.0 * ii_409[k]
                   + f_0 * li_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ii_410, ii_411, ii_412, ii_413, \
                         ii_414, li_410, li_411, li_412, li_413, \
                         li_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -3.0 * ii_410[k]
                   + f_0 * li_410[k];

        t_411[k] = -3.0 * ii_411[k]
                   + f_0 * li_411[k];

        t_412[k] = -3.0 * ii_412[k]
                   + f_0 * li_412[k];

        t_413[k] = -3.0 * ii_413[k]
                   + f_0 * li_413[k];

        t_414[k] = -3.0 * ii_414[k]
                   + f_0 * li_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ii_415, ii_416, ii_417, ii_418, \
                         ii_419, li_415, li_416, li_417, li_418, \
                         li_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -3.0 * ii_415[k]
                   + f_0 * li_415[k];

        t_416[k] = -3.0 * ii_416[k]
                   + f_0 * li_416[k];

        t_417[k] = -3.0 * ii_417[k]
                   + f_0 * li_417[k];

        t_418[k] = -3.0 * ii_418[k]
                   + f_0 * li_418[k];

        t_419[k] = -3.0 * ii_419[k]
                   + f_0 * li_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ii_420, ii_421, ii_422, ii_423, \
                         ii_424, li_420, li_421, li_422, li_423, \
                         li_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -2.0 * ii_420[k]
                   + f_0 * li_420[k];

        t_421[k] = -2.0 * ii_421[k]
                   + f_0 * li_421[k];

        t_422[k] = -2.0 * ii_422[k]
                   + f_0 * li_422[k];

        t_423[k] = -2.0 * ii_423[k]
                   + f_0 * li_423[k];

        t_424[k] = -2.0 * ii_424[k]
                   + f_0 * li_424[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ii_425, ii_426, ii_427, ii_428, \
                         ii_429, li_425, li_426, li_427, li_428, \
                         li_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -2.0 * ii_425[k]
                   + f_0 * li_425[k];

        t_426[k] = -2.0 * ii_426[k]
                   + f_0 * li_426[k];

        t_427[k] = -2.0 * ii_427[k]
                   + f_0 * li_427[k];

        t_428[k] = -2.0 * ii_428[k]
                   + f_0 * li_428[k];

        t_429[k] = -2.0 * ii_429[k]
                   + f_0 * li_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ii_430, ii_431, ii_432, ii_433, \
                         ii_434, li_430, li_431, li_432, li_433, \
                         li_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -2.0 * ii_430[k]
                   + f_0 * li_430[k];

        t_431[k] = -2.0 * ii_431[k]
                   + f_0 * li_431[k];

        t_432[k] = -2.0 * ii_432[k]
                   + f_0 * li_432[k];

        t_433[k] = -2.0 * ii_433[k]
                   + f_0 * li_433[k];

        t_434[k] = -2.0 * ii_434[k]
                   + f_0 * li_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ii_435, ii_436, ii_437, ii_438, \
                         ii_439, li_435, li_436, li_437, li_438, \
                         li_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -2.0 * ii_435[k]
                   + f_0 * li_435[k];

        t_436[k] = -2.0 * ii_436[k]
                   + f_0 * li_436[k];

        t_437[k] = -2.0 * ii_437[k]
                   + f_0 * li_437[k];

        t_438[k] = -2.0 * ii_438[k]
                   + f_0 * li_438[k];

        t_439[k] = -2.0 * ii_439[k]
                   + f_0 * li_439[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ii_440, ii_441, ii_442, ii_443, \
                         ii_444, li_440, li_441, li_442, li_443, \
                         li_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -2.0 * ii_440[k]
                   + f_0 * li_440[k];

        t_441[k] = -2.0 * ii_441[k]
                   + f_0 * li_441[k];

        t_442[k] = -2.0 * ii_442[k]
                   + f_0 * li_442[k];

        t_443[k] = -2.0 * ii_443[k]
                   + f_0 * li_443[k];

        t_444[k] = -2.0 * ii_444[k]
                   + f_0 * li_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ii_445, ii_446, ii_447, ii_448, \
                         ii_449, li_445, li_446, li_447, li_448, \
                         li_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -2.0 * ii_445[k]
                   + f_0 * li_445[k];

        t_446[k] = -2.0 * ii_446[k]
                   + f_0 * li_446[k];

        t_447[k] = -2.0 * ii_447[k]
                   + f_0 * li_447[k];

        t_448[k] = -2.0 * ii_448[k]
                   + f_0 * li_448[k];

        t_449[k] = -2.0 * ii_449[k]
                   + f_0 * li_449[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_450 = buffer.data(ii + 450);
    const auto *ii_451 = buffer.data(ii + 451);
    const auto *ii_452 = buffer.data(ii + 452);
    const auto *ii_453 = buffer.data(ii + 453);
    const auto *ii_454 = buffer.data(ii + 454);
    const auto *ii_455 = buffer.data(ii + 455);
    const auto *ii_456 = buffer.data(ii + 456);
    const auto *ii_457 = buffer.data(ii + 457);
    const auto *ii_458 = buffer.data(ii + 458);
    const auto *ii_459 = buffer.data(ii + 459);
    const auto *ii_460 = buffer.data(ii + 460);
    const auto *ii_461 = buffer.data(ii + 461);
    const auto *ii_462 = buffer.data(ii + 462);
    const auto *ii_463 = buffer.data(ii + 463);
    const auto *ii_464 = buffer.data(ii + 464);
    const auto *ii_465 = buffer.data(ii + 465);
    const auto *ii_466 = buffer.data(ii + 466);
    const auto *ii_467 = buffer.data(ii + 467);
    const auto *ii_468 = buffer.data(ii + 468);
    const auto *ii_469 = buffer.data(ii + 469);
    const auto *ii_470 = buffer.data(ii + 470);
    const auto *ii_471 = buffer.data(ii + 471);
    const auto *ii_472 = buffer.data(ii + 472);
    const auto *ii_473 = buffer.data(ii + 473);
    const auto *ii_474 = buffer.data(ii + 474);
    const auto *ii_475 = buffer.data(ii + 475);
    const auto *ii_476 = buffer.data(ii + 476);
    const auto *ii_477 = buffer.data(ii + 477);
    const auto *ii_478 = buffer.data(ii + 478);
    const auto *ii_479 = buffer.data(ii + 479);
    const auto *ii_480 = buffer.data(ii + 480);
    const auto *ii_481 = buffer.data(ii + 481);
    const auto *ii_482 = buffer.data(ii + 482);
    const auto *ii_483 = buffer.data(ii + 483);
    const auto *ii_484 = buffer.data(ii + 484);
    const auto *ii_485 = buffer.data(ii + 485);
    const auto *ii_486 = buffer.data(ii + 486);
    const auto *ii_487 = buffer.data(ii + 487);
    const auto *ii_488 = buffer.data(ii + 488);
    const auto *ii_489 = buffer.data(ii + 489);
    const auto *ii_490 = buffer.data(ii + 490);
    const auto *ii_491 = buffer.data(ii + 491);
    const auto *ii_492 = buffer.data(ii + 492);
    const auto *ii_493 = buffer.data(ii + 493);
    const auto *ii_494 = buffer.data(ii + 494);
    const auto *ii_495 = buffer.data(ii + 495);
    const auto *ii_496 = buffer.data(ii + 496);
    const auto *ii_497 = buffer.data(ii + 497);
    const auto *ii_498 = buffer.data(ii + 498);
    const auto *ii_499 = buffer.data(ii + 499);
    const auto *ii_500 = buffer.data(ii + 500);
    const auto *ii_501 = buffer.data(ii + 501);
    const auto *ii_502 = buffer.data(ii + 502);
    const auto *ii_503 = buffer.data(ii + 503);
    const auto *ii_504 = buffer.data(ii + 504);
    const auto *ii_505 = buffer.data(ii + 505);
    const auto *ii_506 = buffer.data(ii + 506);
    const auto *ii_507 = buffer.data(ii + 507);
    const auto *ii_508 = buffer.data(ii + 508);
    const auto *ii_509 = buffer.data(ii + 509);
    const auto *ii_510 = buffer.data(ii + 510);
    const auto *ii_511 = buffer.data(ii + 511);
    const auto *ii_512 = buffer.data(ii + 512);
    const auto *ii_513 = buffer.data(ii + 513);
    const auto *ii_514 = buffer.data(ii + 514);
    const auto *ii_515 = buffer.data(ii + 515);
    const auto *ii_516 = buffer.data(ii + 516);
    const auto *ii_517 = buffer.data(ii + 517);
    const auto *ii_518 = buffer.data(ii + 518);
    const auto *ii_519 = buffer.data(ii + 519);
    const auto *ii_520 = buffer.data(ii + 520);
    const auto *ii_521 = buffer.data(ii + 521);
    const auto *ii_522 = buffer.data(ii + 522);
    const auto *ii_523 = buffer.data(ii + 523);
    const auto *ii_524 = buffer.data(ii + 524);
    const auto *ii_525 = buffer.data(ii + 525);
    const auto *ii_526 = buffer.data(ii + 526);
    const auto *ii_527 = buffer.data(ii + 527);
    const auto *ii_528 = buffer.data(ii + 528);
    const auto *ii_529 = buffer.data(ii + 529);
    const auto *ii_530 = buffer.data(ii + 530);
    const auto *ii_531 = buffer.data(ii + 531);
    const auto *ii_532 = buffer.data(ii + 532);
    const auto *ii_533 = buffer.data(ii + 533);
    const auto *ii_534 = buffer.data(ii + 534);
    const auto *ii_535 = buffer.data(ii + 535);
    const auto *ii_536 = buffer.data(ii + 536);
    const auto *ii_537 = buffer.data(ii + 537);
    const auto *ii_538 = buffer.data(ii + 538);
    const auto *ii_539 = buffer.data(ii + 539);
    const auto *ii_540 = buffer.data(ii + 540);
    const auto *ii_541 = buffer.data(ii + 541);
    const auto *ii_542 = buffer.data(ii + 542);
    const auto *ii_543 = buffer.data(ii + 543);
    const auto *ii_544 = buffer.data(ii + 544);
    const auto *ii_545 = buffer.data(ii + 545);
    const auto *ii_546 = buffer.data(ii + 546);
    const auto *ii_547 = buffer.data(ii + 547);
    const auto *ii_548 = buffer.data(ii + 548);
    const auto *ii_549 = buffer.data(ii + 549);
    const auto *ii_550 = buffer.data(ii + 550);
    const auto *ii_551 = buffer.data(ii + 551);
    const auto *ii_552 = buffer.data(ii + 552);
    const auto *ii_553 = buffer.data(ii + 553);
    const auto *ii_554 = buffer.data(ii + 554);
    const auto *ii_555 = buffer.data(ii + 555);
    const auto *ii_556 = buffer.data(ii + 556);
    const auto *ii_557 = buffer.data(ii + 557);
    const auto *ii_558 = buffer.data(ii + 558);
    const auto *ii_559 = buffer.data(ii + 559);
    const auto *ii_560 = buffer.data(ii + 560);
    const auto *ii_561 = buffer.data(ii + 561);
    const auto *ii_562 = buffer.data(ii + 562);
    const auto *ii_563 = buffer.data(ii + 563);
    const auto *ii_564 = buffer.data(ii + 564);
    const auto *ii_565 = buffer.data(ii + 565);
    const auto *ii_566 = buffer.data(ii + 566);
    const auto *ii_567 = buffer.data(ii + 567);
    const auto *ii_568 = buffer.data(ii + 568);
    const auto *ii_569 = buffer.data(ii + 569);
    const auto *ii_570 = buffer.data(ii + 570);
    const auto *ii_571 = buffer.data(ii + 571);
    const auto *ii_572 = buffer.data(ii + 572);
    const auto *ii_573 = buffer.data(ii + 573);
    const auto *ii_574 = buffer.data(ii + 574);
    const auto *ii_575 = buffer.data(ii + 575);
    const auto *ii_576 = buffer.data(ii + 576);
    const auto *ii_577 = buffer.data(ii + 577);
    const auto *ii_578 = buffer.data(ii + 578);
    const auto *ii_579 = buffer.data(ii + 579);
    const auto *ii_580 = buffer.data(ii + 580);
    const auto *ii_581 = buffer.data(ii + 581);
    const auto *ii_582 = buffer.data(ii + 582);
    const auto *ii_583 = buffer.data(ii + 583);
    const auto *ii_584 = buffer.data(ii + 584);
    const auto *ii_585 = buffer.data(ii + 585);
    const auto *ii_586 = buffer.data(ii + 586);
    const auto *ii_587 = buffer.data(ii + 587);
    const auto *ii_588 = buffer.data(ii + 588);
    const auto *ii_589 = buffer.data(ii + 589);
    const auto *ii_590 = buffer.data(ii + 590);
    const auto *ii_591 = buffer.data(ii + 591);
    const auto *ii_592 = buffer.data(ii + 592);
    const auto *ii_593 = buffer.data(ii + 593);
    const auto *ii_594 = buffer.data(ii + 594);
    const auto *ii_595 = buffer.data(ii + 595);
    const auto *ii_596 = buffer.data(ii + 596);
    const auto *ii_597 = buffer.data(ii + 597);
    const auto *ii_598 = buffer.data(ii + 598);
    const auto *ii_599 = buffer.data(ii + 599);

    const auto *li_450 = buffer.data(li + 450);
    const auto *li_451 = buffer.data(li + 451);
    const auto *li_452 = buffer.data(li + 452);
    const auto *li_453 = buffer.data(li + 453);
    const auto *li_454 = buffer.data(li + 454);
    const auto *li_455 = buffer.data(li + 455);
    const auto *li_456 = buffer.data(li + 456);
    const auto *li_457 = buffer.data(li + 457);
    const auto *li_458 = buffer.data(li + 458);
    const auto *li_459 = buffer.data(li + 459);
    const auto *li_460 = buffer.data(li + 460);
    const auto *li_461 = buffer.data(li + 461);
    const auto *li_462 = buffer.data(li + 462);
    const auto *li_463 = buffer.data(li + 463);
    const auto *li_464 = buffer.data(li + 464);
    const auto *li_465 = buffer.data(li + 465);
    const auto *li_466 = buffer.data(li + 466);
    const auto *li_467 = buffer.data(li + 467);
    const auto *li_468 = buffer.data(li + 468);
    const auto *li_469 = buffer.data(li + 469);
    const auto *li_470 = buffer.data(li + 470);
    const auto *li_471 = buffer.data(li + 471);
    const auto *li_472 = buffer.data(li + 472);
    const auto *li_473 = buffer.data(li + 473);
    const auto *li_474 = buffer.data(li + 474);
    const auto *li_475 = buffer.data(li + 475);
    const auto *li_476 = buffer.data(li + 476);
    const auto *li_477 = buffer.data(li + 477);
    const auto *li_478 = buffer.data(li + 478);
    const auto *li_479 = buffer.data(li + 479);
    const auto *li_480 = buffer.data(li + 480);
    const auto *li_481 = buffer.data(li + 481);
    const auto *li_482 = buffer.data(li + 482);
    const auto *li_483 = buffer.data(li + 483);
    const auto *li_484 = buffer.data(li + 484);
    const auto *li_485 = buffer.data(li + 485);
    const auto *li_486 = buffer.data(li + 486);
    const auto *li_487 = buffer.data(li + 487);
    const auto *li_488 = buffer.data(li + 488);
    const auto *li_489 = buffer.data(li + 489);
    const auto *li_490 = buffer.data(li + 490);
    const auto *li_491 = buffer.data(li + 491);
    const auto *li_492 = buffer.data(li + 492);
    const auto *li_493 = buffer.data(li + 493);
    const auto *li_494 = buffer.data(li + 494);
    const auto *li_495 = buffer.data(li + 495);
    const auto *li_496 = buffer.data(li + 496);
    const auto *li_497 = buffer.data(li + 497);
    const auto *li_498 = buffer.data(li + 498);
    const auto *li_499 = buffer.data(li + 499);
    const auto *li_500 = buffer.data(li + 500);
    const auto *li_501 = buffer.data(li + 501);
    const auto *li_502 = buffer.data(li + 502);
    const auto *li_503 = buffer.data(li + 503);
    const auto *li_504 = buffer.data(li + 504);
    const auto *li_505 = buffer.data(li + 505);
    const auto *li_506 = buffer.data(li + 506);
    const auto *li_507 = buffer.data(li + 507);
    const auto *li_508 = buffer.data(li + 508);
    const auto *li_509 = buffer.data(li + 509);
    const auto *li_510 = buffer.data(li + 510);
    const auto *li_511 = buffer.data(li + 511);
    const auto *li_512 = buffer.data(li + 512);
    const auto *li_513 = buffer.data(li + 513);
    const auto *li_514 = buffer.data(li + 514);
    const auto *li_515 = buffer.data(li + 515);
    const auto *li_516 = buffer.data(li + 516);
    const auto *li_517 = buffer.data(li + 517);
    const auto *li_518 = buffer.data(li + 518);
    const auto *li_519 = buffer.data(li + 519);
    const auto *li_520 = buffer.data(li + 520);
    const auto *li_521 = buffer.data(li + 521);
    const auto *li_522 = buffer.data(li + 522);
    const auto *li_523 = buffer.data(li + 523);
    const auto *li_524 = buffer.data(li + 524);
    const auto *li_525 = buffer.data(li + 525);
    const auto *li_526 = buffer.data(li + 526);
    const auto *li_527 = buffer.data(li + 527);
    const auto *li_528 = buffer.data(li + 528);
    const auto *li_529 = buffer.data(li + 529);
    const auto *li_530 = buffer.data(li + 530);
    const auto *li_531 = buffer.data(li + 531);
    const auto *li_532 = buffer.data(li + 532);
    const auto *li_533 = buffer.data(li + 533);
    const auto *li_534 = buffer.data(li + 534);
    const auto *li_535 = buffer.data(li + 535);
    const auto *li_536 = buffer.data(li + 536);
    const auto *li_537 = buffer.data(li + 537);
    const auto *li_538 = buffer.data(li + 538);
    const auto *li_539 = buffer.data(li + 539);
    const auto *li_540 = buffer.data(li + 540);
    const auto *li_541 = buffer.data(li + 541);
    const auto *li_542 = buffer.data(li + 542);
    const auto *li_543 = buffer.data(li + 543);
    const auto *li_544 = buffer.data(li + 544);
    const auto *li_545 = buffer.data(li + 545);
    const auto *li_546 = buffer.data(li + 546);
    const auto *li_547 = buffer.data(li + 547);
    const auto *li_548 = buffer.data(li + 548);
    const auto *li_549 = buffer.data(li + 549);
    const auto *li_550 = buffer.data(li + 550);
    const auto *li_551 = buffer.data(li + 551);
    const auto *li_552 = buffer.data(li + 552);
    const auto *li_553 = buffer.data(li + 553);
    const auto *li_554 = buffer.data(li + 554);
    const auto *li_555 = buffer.data(li + 555);
    const auto *li_556 = buffer.data(li + 556);
    const auto *li_557 = buffer.data(li + 557);
    const auto *li_558 = buffer.data(li + 558);
    const auto *li_559 = buffer.data(li + 559);
    const auto *li_560 = buffer.data(li + 560);
    const auto *li_561 = buffer.data(li + 561);
    const auto *li_562 = buffer.data(li + 562);
    const auto *li_563 = buffer.data(li + 563);
    const auto *li_564 = buffer.data(li + 564);
    const auto *li_565 = buffer.data(li + 565);
    const auto *li_566 = buffer.data(li + 566);
    const auto *li_567 = buffer.data(li + 567);
    const auto *li_568 = buffer.data(li + 568);
    const auto *li_569 = buffer.data(li + 569);
    const auto *li_570 = buffer.data(li + 570);
    const auto *li_571 = buffer.data(li + 571);
    const auto *li_572 = buffer.data(li + 572);
    const auto *li_573 = buffer.data(li + 573);
    const auto *li_574 = buffer.data(li + 574);
    const auto *li_575 = buffer.data(li + 575);
    const auto *li_576 = buffer.data(li + 576);
    const auto *li_577 = buffer.data(li + 577);
    const auto *li_578 = buffer.data(li + 578);
    const auto *li_579 = buffer.data(li + 579);
    const auto *li_580 = buffer.data(li + 580);
    const auto *li_581 = buffer.data(li + 581);
    const auto *li_582 = buffer.data(li + 582);
    const auto *li_583 = buffer.data(li + 583);
    const auto *li_584 = buffer.data(li + 584);
    const auto *li_585 = buffer.data(li + 585);
    const auto *li_586 = buffer.data(li + 586);
    const auto *li_587 = buffer.data(li + 587);
    const auto *li_588 = buffer.data(li + 588);
    const auto *li_589 = buffer.data(li + 589);
    const auto *li_590 = buffer.data(li + 590);
    const auto *li_591 = buffer.data(li + 591);
    const auto *li_592 = buffer.data(li + 592);
    const auto *li_593 = buffer.data(li + 593);
    const auto *li_594 = buffer.data(li + 594);
    const auto *li_595 = buffer.data(li + 595);
    const auto *li_596 = buffer.data(li + 596);
    const auto *li_597 = buffer.data(li + 597);
    const auto *li_598 = buffer.data(li + 598);
    const auto *li_599 = buffer.data(li + 599);

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ii_450, ii_451, ii_452, ii_453, \
                         ii_454, li_450, li_451, li_452, li_453, \
                         li_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -2.0 * ii_450[k]
                   + f_0 * li_450[k];

        t_451[k] = -2.0 * ii_451[k]
                   + f_0 * li_451[k];

        t_452[k] = -2.0 * ii_452[k]
                   + f_0 * li_452[k];

        t_453[k] = -2.0 * ii_453[k]
                   + f_0 * li_453[k];

        t_454[k] = -2.0 * ii_454[k]
                   + f_0 * li_454[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ii_455, ii_456, ii_457, ii_458, \
                         ii_459, li_455, li_456, li_457, li_458, \
                         li_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -2.0 * ii_455[k]
                   + f_0 * li_455[k];

        t_456[k] = -2.0 * ii_456[k]
                   + f_0 * li_456[k];

        t_457[k] = -2.0 * ii_457[k]
                   + f_0 * li_457[k];

        t_458[k] = -2.0 * ii_458[k]
                   + f_0 * li_458[k];

        t_459[k] = -2.0 * ii_459[k]
                   + f_0 * li_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ii_460, ii_461, ii_462, ii_463, \
                         ii_464, li_460, li_461, li_462, li_463, \
                         li_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -2.0 * ii_460[k]
                   + f_0 * li_460[k];

        t_461[k] = -2.0 * ii_461[k]
                   + f_0 * li_461[k];

        t_462[k] = -2.0 * ii_462[k]
                   + f_0 * li_462[k];

        t_463[k] = -2.0 * ii_463[k]
                   + f_0 * li_463[k];

        t_464[k] = -2.0 * ii_464[k]
                   + f_0 * li_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ii_465, ii_466, ii_467, ii_468, \
                         ii_469, li_465, li_466, li_467, li_468, \
                         li_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -2.0 * ii_465[k]
                   + f_0 * li_465[k];

        t_466[k] = -2.0 * ii_466[k]
                   + f_0 * li_466[k];

        t_467[k] = -2.0 * ii_467[k]
                   + f_0 * li_467[k];

        t_468[k] = -2.0 * ii_468[k]
                   + f_0 * li_468[k];

        t_469[k] = -2.0 * ii_469[k]
                   + f_0 * li_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ii_470, ii_471, ii_472, ii_473, \
                         ii_474, li_470, li_471, li_472, li_473, \
                         li_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -2.0 * ii_470[k]
                   + f_0 * li_470[k];

        t_471[k] = -2.0 * ii_471[k]
                   + f_0 * li_471[k];

        t_472[k] = -2.0 * ii_472[k]
                   + f_0 * li_472[k];

        t_473[k] = -2.0 * ii_473[k]
                   + f_0 * li_473[k];

        t_474[k] = -2.0 * ii_474[k]
                   + f_0 * li_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ii_475, ii_476, ii_477, ii_478, \
                         ii_479, li_475, li_476, li_477, li_478, \
                         li_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -2.0 * ii_475[k]
                   + f_0 * li_475[k];

        t_476[k] = -2.0 * ii_476[k]
                   + f_0 * li_476[k];

        t_477[k] = -2.0 * ii_477[k]
                   + f_0 * li_477[k];

        t_478[k] = -2.0 * ii_478[k]
                   + f_0 * li_478[k];

        t_479[k] = -2.0 * ii_479[k]
                   + f_0 * li_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ii_480, ii_481, ii_482, ii_483, \
                         ii_484, li_480, li_481, li_482, li_483, \
                         li_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -2.0 * ii_480[k]
                   + f_0 * li_480[k];

        t_481[k] = -2.0 * ii_481[k]
                   + f_0 * li_481[k];

        t_482[k] = -2.0 * ii_482[k]
                   + f_0 * li_482[k];

        t_483[k] = -2.0 * ii_483[k]
                   + f_0 * li_483[k];

        t_484[k] = -2.0 * ii_484[k]
                   + f_0 * li_484[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ii_485, ii_486, ii_487, ii_488, \
                         ii_489, li_485, li_486, li_487, li_488, \
                         li_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -2.0 * ii_485[k]
                   + f_0 * li_485[k];

        t_486[k] = -2.0 * ii_486[k]
                   + f_0 * li_486[k];

        t_487[k] = -2.0 * ii_487[k]
                   + f_0 * li_487[k];

        t_488[k] = -2.0 * ii_488[k]
                   + f_0 * li_488[k];

        t_489[k] = -2.0 * ii_489[k]
                   + f_0 * li_489[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ii_490, ii_491, ii_492, ii_493, \
                         ii_494, li_490, li_491, li_492, li_493, \
                         li_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -2.0 * ii_490[k]
                   + f_0 * li_490[k];

        t_491[k] = -2.0 * ii_491[k]
                   + f_0 * li_491[k];

        t_492[k] = -2.0 * ii_492[k]
                   + f_0 * li_492[k];

        t_493[k] = -2.0 * ii_493[k]
                   + f_0 * li_493[k];

        t_494[k] = -2.0 * ii_494[k]
                   + f_0 * li_494[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ii_495, ii_496, ii_497, ii_498, \
                         ii_499, li_495, li_496, li_497, li_498, \
                         li_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -2.0 * ii_495[k]
                   + f_0 * li_495[k];

        t_496[k] = -2.0 * ii_496[k]
                   + f_0 * li_496[k];

        t_497[k] = -2.0 * ii_497[k]
                   + f_0 * li_497[k];

        t_498[k] = -2.0 * ii_498[k]
                   + f_0 * li_498[k];

        t_499[k] = -2.0 * ii_499[k]
                   + f_0 * li_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ii_500, ii_501, ii_502, ii_503, \
                         ii_504, li_500, li_501, li_502, li_503, \
                         li_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -2.0 * ii_500[k]
                   + f_0 * li_500[k];

        t_501[k] = -2.0 * ii_501[k]
                   + f_0 * li_501[k];

        t_502[k] = -2.0 * ii_502[k]
                   + f_0 * li_502[k];

        t_503[k] = -2.0 * ii_503[k]
                   + f_0 * li_503[k];

        t_504[k] = -2.0 * ii_504[k]
                   + f_0 * li_504[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ii_505, ii_506, ii_507, ii_508, \
                         ii_509, li_505, li_506, li_507, li_508, \
                         li_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * ii_505[k]
                   + f_0 * li_505[k];

        t_506[k] = -2.0 * ii_506[k]
                   + f_0 * li_506[k];

        t_507[k] = -2.0 * ii_507[k]
                   + f_0 * li_507[k];

        t_508[k] = -2.0 * ii_508[k]
                   + f_0 * li_508[k];

        t_509[k] = -2.0 * ii_509[k]
                   + f_0 * li_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ii_510, ii_511, ii_512, ii_513, \
                         ii_514, li_510, li_511, li_512, li_513, \
                         li_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -2.0 * ii_510[k]
                   + f_0 * li_510[k];

        t_511[k] = -2.0 * ii_511[k]
                   + f_0 * li_511[k];

        t_512[k] = -2.0 * ii_512[k]
                   + f_0 * li_512[k];

        t_513[k] = -2.0 * ii_513[k]
                   + f_0 * li_513[k];

        t_514[k] = -2.0 * ii_514[k]
                   + f_0 * li_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ii_515, ii_516, ii_517, ii_518, \
                         ii_519, li_515, li_516, li_517, li_518, \
                         li_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -2.0 * ii_515[k]
                   + f_0 * li_515[k];

        t_516[k] = -2.0 * ii_516[k]
                   + f_0 * li_516[k];

        t_517[k] = -2.0 * ii_517[k]
                   + f_0 * li_517[k];

        t_518[k] = -2.0 * ii_518[k]
                   + f_0 * li_518[k];

        t_519[k] = -2.0 * ii_519[k]
                   + f_0 * li_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ii_520, ii_521, ii_522, ii_523, \
                         ii_524, li_520, li_521, li_522, li_523, \
                         li_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -2.0 * ii_520[k]
                   + f_0 * li_520[k];

        t_521[k] = -2.0 * ii_521[k]
                   + f_0 * li_521[k];

        t_522[k] = -2.0 * ii_522[k]
                   + f_0 * li_522[k];

        t_523[k] = -2.0 * ii_523[k]
                   + f_0 * li_523[k];

        t_524[k] = -2.0 * ii_524[k]
                   + f_0 * li_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ii_525, ii_526, ii_527, ii_528, \
                         ii_529, li_525, li_526, li_527, li_528, \
                         li_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -2.0 * ii_525[k]
                   + f_0 * li_525[k];

        t_526[k] = -2.0 * ii_526[k]
                   + f_0 * li_526[k];

        t_527[k] = -2.0 * ii_527[k]
                   + f_0 * li_527[k];

        t_528[k] = -2.0 * ii_528[k]
                   + f_0 * li_528[k];

        t_529[k] = -2.0 * ii_529[k]
                   + f_0 * li_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ii_530, ii_531, ii_532, ii_533, \
                         ii_534, li_530, li_531, li_532, li_533, \
                         li_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -2.0 * ii_530[k]
                   + f_0 * li_530[k];

        t_531[k] = -2.0 * ii_531[k]
                   + f_0 * li_531[k];

        t_532[k] = -2.0 * ii_532[k]
                   + f_0 * li_532[k];

        t_533[k] = -2.0 * ii_533[k]
                   + f_0 * li_533[k];

        t_534[k] = -2.0 * ii_534[k]
                   + f_0 * li_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ii_535, ii_536, ii_537, ii_538, \
                         ii_539, li_535, li_536, li_537, li_538, \
                         li_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -2.0 * ii_535[k]
                   + f_0 * li_535[k];

        t_536[k] = -2.0 * ii_536[k]
                   + f_0 * li_536[k];

        t_537[k] = -2.0 * ii_537[k]
                   + f_0 * li_537[k];

        t_538[k] = -2.0 * ii_538[k]
                   + f_0 * li_538[k];

        t_539[k] = -2.0 * ii_539[k]
                   + f_0 * li_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ii_540, ii_541, ii_542, ii_543, \
                         ii_544, li_540, li_541, li_542, li_543, \
                         li_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -2.0 * ii_540[k]
                   + f_0 * li_540[k];

        t_541[k] = -2.0 * ii_541[k]
                   + f_0 * li_541[k];

        t_542[k] = -2.0 * ii_542[k]
                   + f_0 * li_542[k];

        t_543[k] = -2.0 * ii_543[k]
                   + f_0 * li_543[k];

        t_544[k] = -2.0 * ii_544[k]
                   + f_0 * li_544[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ii_545, ii_546, ii_547, ii_548, \
                         ii_549, li_545, li_546, li_547, li_548, \
                         li_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -2.0 * ii_545[k]
                   + f_0 * li_545[k];

        t_546[k] = -2.0 * ii_546[k]
                   + f_0 * li_546[k];

        t_547[k] = -2.0 * ii_547[k]
                   + f_0 * li_547[k];

        t_548[k] = -2.0 * ii_548[k]
                   + f_0 * li_548[k];

        t_549[k] = -2.0 * ii_549[k]
                   + f_0 * li_549[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ii_550, ii_551, ii_552, ii_553, \
                         ii_554, li_550, li_551, li_552, li_553, \
                         li_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -2.0 * ii_550[k]
                   + f_0 * li_550[k];

        t_551[k] = -2.0 * ii_551[k]
                   + f_0 * li_551[k];

        t_552[k] = -2.0 * ii_552[k]
                   + f_0 * li_552[k];

        t_553[k] = -2.0 * ii_553[k]
                   + f_0 * li_553[k];

        t_554[k] = -2.0 * ii_554[k]
                   + f_0 * li_554[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ii_555, ii_556, ii_557, ii_558, \
                         ii_559, li_555, li_556, li_557, li_558, \
                         li_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -2.0 * ii_555[k]
                   + f_0 * li_555[k];

        t_556[k] = -2.0 * ii_556[k]
                   + f_0 * li_556[k];

        t_557[k] = -2.0 * ii_557[k]
                   + f_0 * li_557[k];

        t_558[k] = -2.0 * ii_558[k]
                   + f_0 * li_558[k];

        t_559[k] = -2.0 * ii_559[k]
                   + f_0 * li_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ii_560, ii_561, ii_562, ii_563, \
                         ii_564, li_560, li_561, li_562, li_563, \
                         li_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = -2.0 * ii_560[k]
                   + f_0 * li_560[k];

        t_561[k] = -2.0 * ii_561[k]
                   + f_0 * li_561[k];

        t_562[k] = -2.0 * ii_562[k]
                   + f_0 * li_562[k];

        t_563[k] = -2.0 * ii_563[k]
                   + f_0 * li_563[k];

        t_564[k] = -2.0 * ii_564[k]
                   + f_0 * li_564[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ii_565, ii_566, ii_567, ii_568, \
                         ii_569, li_565, li_566, li_567, li_568, \
                         li_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = -2.0 * ii_565[k]
                   + f_0 * li_565[k];

        t_566[k] = -2.0 * ii_566[k]
                   + f_0 * li_566[k];

        t_567[k] = -2.0 * ii_567[k]
                   + f_0 * li_567[k];

        t_568[k] = -2.0 * ii_568[k]
                   + f_0 * li_568[k];

        t_569[k] = -2.0 * ii_569[k]
                   + f_0 * li_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ii_570, ii_571, ii_572, ii_573, \
                         ii_574, li_570, li_571, li_572, li_573, \
                         li_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -2.0 * ii_570[k]
                   + f_0 * li_570[k];

        t_571[k] = -2.0 * ii_571[k]
                   + f_0 * li_571[k];

        t_572[k] = -2.0 * ii_572[k]
                   + f_0 * li_572[k];

        t_573[k] = -2.0 * ii_573[k]
                   + f_0 * li_573[k];

        t_574[k] = -2.0 * ii_574[k]
                   + f_0 * li_574[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ii_575, ii_576, ii_577, ii_578, \
                         ii_579, li_575, li_576, li_577, li_578, \
                         li_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = -2.0 * ii_575[k]
                   + f_0 * li_575[k];

        t_576[k] = -2.0 * ii_576[k]
                   + f_0 * li_576[k];

        t_577[k] = -2.0 * ii_577[k]
                   + f_0 * li_577[k];

        t_578[k] = -2.0 * ii_578[k]
                   + f_0 * li_578[k];

        t_579[k] = -2.0 * ii_579[k]
                   + f_0 * li_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ii_580, ii_581, ii_582, ii_583, \
                         ii_584, li_580, li_581, li_582, li_583, \
                         li_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = -2.0 * ii_580[k]
                   + f_0 * li_580[k];

        t_581[k] = -2.0 * ii_581[k]
                   + f_0 * li_581[k];

        t_582[k] = -2.0 * ii_582[k]
                   + f_0 * li_582[k];

        t_583[k] = -2.0 * ii_583[k]
                   + f_0 * li_583[k];

        t_584[k] = -2.0 * ii_584[k]
                   + f_0 * li_584[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, t_589, ii_585, ii_586, ii_587, ii_588, \
                         ii_589, li_585, li_586, li_587, li_588, \
                         li_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = -2.0 * ii_585[k]
                   + f_0 * li_585[k];

        t_586[k] = -2.0 * ii_586[k]
                   + f_0 * li_586[k];

        t_587[k] = -2.0 * ii_587[k]
                   + f_0 * li_587[k];

        t_588[k] = -ii_588[k]
                   + f_0 * li_588[k];

        t_589[k] = -ii_589[k]
                   + f_0 * li_589[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, ii_590, ii_591, ii_592, ii_593, \
                         ii_594, li_590, li_591, li_592, li_593, \
                         li_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -ii_590[k]
                   + f_0 * li_590[k];

        t_591[k] = -ii_591[k]
                   + f_0 * li_591[k];

        t_592[k] = -ii_592[k]
                   + f_0 * li_592[k];

        t_593[k] = -ii_593[k]
                   + f_0 * li_593[k];

        t_594[k] = -ii_594[k]
                   + f_0 * li_594[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, ii_595, ii_596, ii_597, ii_598, \
                         ii_599, li_595, li_596, li_597, li_598, \
                         li_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -ii_595[k]
                   + f_0 * li_595[k];

        t_596[k] = -ii_596[k]
                   + f_0 * li_596[k];

        t_597[k] = -ii_597[k]
                   + f_0 * li_597[k];

        t_598[k] = -ii_598[k]
                   + f_0 * li_598[k];

        t_599[k] = -ii_599[k]
                   + f_0 * li_599[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_600 = buffer.data(ii + 600);
    const auto *ii_601 = buffer.data(ii + 601);
    const auto *ii_602 = buffer.data(ii + 602);
    const auto *ii_603 = buffer.data(ii + 603);
    const auto *ii_604 = buffer.data(ii + 604);
    const auto *ii_605 = buffer.data(ii + 605);
    const auto *ii_606 = buffer.data(ii + 606);
    const auto *ii_607 = buffer.data(ii + 607);
    const auto *ii_608 = buffer.data(ii + 608);
    const auto *ii_609 = buffer.data(ii + 609);
    const auto *ii_610 = buffer.data(ii + 610);
    const auto *ii_611 = buffer.data(ii + 611);
    const auto *ii_612 = buffer.data(ii + 612);
    const auto *ii_613 = buffer.data(ii + 613);
    const auto *ii_614 = buffer.data(ii + 614);
    const auto *ii_615 = buffer.data(ii + 615);
    const auto *ii_616 = buffer.data(ii + 616);
    const auto *ii_617 = buffer.data(ii + 617);
    const auto *ii_618 = buffer.data(ii + 618);
    const auto *ii_619 = buffer.data(ii + 619);
    const auto *ii_620 = buffer.data(ii + 620);
    const auto *ii_621 = buffer.data(ii + 621);
    const auto *ii_622 = buffer.data(ii + 622);
    const auto *ii_623 = buffer.data(ii + 623);
    const auto *ii_624 = buffer.data(ii + 624);
    const auto *ii_625 = buffer.data(ii + 625);
    const auto *ii_626 = buffer.data(ii + 626);
    const auto *ii_627 = buffer.data(ii + 627);
    const auto *ii_628 = buffer.data(ii + 628);
    const auto *ii_629 = buffer.data(ii + 629);
    const auto *ii_630 = buffer.data(ii + 630);
    const auto *ii_631 = buffer.data(ii + 631);
    const auto *ii_632 = buffer.data(ii + 632);
    const auto *ii_633 = buffer.data(ii + 633);
    const auto *ii_634 = buffer.data(ii + 634);
    const auto *ii_635 = buffer.data(ii + 635);
    const auto *ii_636 = buffer.data(ii + 636);
    const auto *ii_637 = buffer.data(ii + 637);
    const auto *ii_638 = buffer.data(ii + 638);
    const auto *ii_639 = buffer.data(ii + 639);
    const auto *ii_640 = buffer.data(ii + 640);
    const auto *ii_641 = buffer.data(ii + 641);
    const auto *ii_642 = buffer.data(ii + 642);
    const auto *ii_643 = buffer.data(ii + 643);
    const auto *ii_644 = buffer.data(ii + 644);
    const auto *ii_645 = buffer.data(ii + 645);
    const auto *ii_646 = buffer.data(ii + 646);
    const auto *ii_647 = buffer.data(ii + 647);
    const auto *ii_648 = buffer.data(ii + 648);
    const auto *ii_649 = buffer.data(ii + 649);
    const auto *ii_650 = buffer.data(ii + 650);
    const auto *ii_651 = buffer.data(ii + 651);
    const auto *ii_652 = buffer.data(ii + 652);
    const auto *ii_653 = buffer.data(ii + 653);
    const auto *ii_654 = buffer.data(ii + 654);
    const auto *ii_655 = buffer.data(ii + 655);
    const auto *ii_656 = buffer.data(ii + 656);
    const auto *ii_657 = buffer.data(ii + 657);
    const auto *ii_658 = buffer.data(ii + 658);
    const auto *ii_659 = buffer.data(ii + 659);
    const auto *ii_660 = buffer.data(ii + 660);
    const auto *ii_661 = buffer.data(ii + 661);
    const auto *ii_662 = buffer.data(ii + 662);
    const auto *ii_663 = buffer.data(ii + 663);
    const auto *ii_664 = buffer.data(ii + 664);
    const auto *ii_665 = buffer.data(ii + 665);
    const auto *ii_666 = buffer.data(ii + 666);
    const auto *ii_667 = buffer.data(ii + 667);
    const auto *ii_668 = buffer.data(ii + 668);
    const auto *ii_669 = buffer.data(ii + 669);
    const auto *ii_670 = buffer.data(ii + 670);
    const auto *ii_671 = buffer.data(ii + 671);
    const auto *ii_672 = buffer.data(ii + 672);
    const auto *ii_673 = buffer.data(ii + 673);
    const auto *ii_674 = buffer.data(ii + 674);
    const auto *ii_675 = buffer.data(ii + 675);
    const auto *ii_676 = buffer.data(ii + 676);
    const auto *ii_677 = buffer.data(ii + 677);
    const auto *ii_678 = buffer.data(ii + 678);
    const auto *ii_679 = buffer.data(ii + 679);
    const auto *ii_680 = buffer.data(ii + 680);
    const auto *ii_681 = buffer.data(ii + 681);
    const auto *ii_682 = buffer.data(ii + 682);
    const auto *ii_683 = buffer.data(ii + 683);
    const auto *ii_684 = buffer.data(ii + 684);
    const auto *ii_685 = buffer.data(ii + 685);
    const auto *ii_686 = buffer.data(ii + 686);
    const auto *ii_687 = buffer.data(ii + 687);
    const auto *ii_688 = buffer.data(ii + 688);
    const auto *ii_689 = buffer.data(ii + 689);
    const auto *ii_690 = buffer.data(ii + 690);
    const auto *ii_691 = buffer.data(ii + 691);
    const auto *ii_692 = buffer.data(ii + 692);
    const auto *ii_693 = buffer.data(ii + 693);
    const auto *ii_694 = buffer.data(ii + 694);
    const auto *ii_695 = buffer.data(ii + 695);
    const auto *ii_696 = buffer.data(ii + 696);
    const auto *ii_697 = buffer.data(ii + 697);
    const auto *ii_698 = buffer.data(ii + 698);
    const auto *ii_699 = buffer.data(ii + 699);
    const auto *ii_700 = buffer.data(ii + 700);
    const auto *ii_701 = buffer.data(ii + 701);
    const auto *ii_702 = buffer.data(ii + 702);
    const auto *ii_703 = buffer.data(ii + 703);
    const auto *ii_704 = buffer.data(ii + 704);
    const auto *ii_705 = buffer.data(ii + 705);
    const auto *ii_706 = buffer.data(ii + 706);
    const auto *ii_707 = buffer.data(ii + 707);
    const auto *ii_708 = buffer.data(ii + 708);
    const auto *ii_709 = buffer.data(ii + 709);
    const auto *ii_710 = buffer.data(ii + 710);
    const auto *ii_711 = buffer.data(ii + 711);
    const auto *ii_712 = buffer.data(ii + 712);
    const auto *ii_713 = buffer.data(ii + 713);
    const auto *ii_714 = buffer.data(ii + 714);
    const auto *ii_715 = buffer.data(ii + 715);
    const auto *ii_716 = buffer.data(ii + 716);
    const auto *ii_717 = buffer.data(ii + 717);
    const auto *ii_718 = buffer.data(ii + 718);
    const auto *ii_719 = buffer.data(ii + 719);
    const auto *ii_720 = buffer.data(ii + 720);
    const auto *ii_721 = buffer.data(ii + 721);
    const auto *ii_722 = buffer.data(ii + 722);
    const auto *ii_723 = buffer.data(ii + 723);
    const auto *ii_724 = buffer.data(ii + 724);
    const auto *ii_725 = buffer.data(ii + 725);
    const auto *ii_726 = buffer.data(ii + 726);
    const auto *ii_727 = buffer.data(ii + 727);
    const auto *ii_728 = buffer.data(ii + 728);
    const auto *ii_729 = buffer.data(ii + 729);
    const auto *ii_730 = buffer.data(ii + 730);
    const auto *ii_731 = buffer.data(ii + 731);
    const auto *ii_732 = buffer.data(ii + 732);
    const auto *ii_733 = buffer.data(ii + 733);
    const auto *ii_734 = buffer.data(ii + 734);
    const auto *ii_735 = buffer.data(ii + 735);
    const auto *ii_736 = buffer.data(ii + 736);
    const auto *ii_737 = buffer.data(ii + 737);
    const auto *ii_738 = buffer.data(ii + 738);
    const auto *ii_739 = buffer.data(ii + 739);
    const auto *ii_740 = buffer.data(ii + 740);
    const auto *ii_741 = buffer.data(ii + 741);
    const auto *ii_742 = buffer.data(ii + 742);
    const auto *ii_743 = buffer.data(ii + 743);
    const auto *ii_744 = buffer.data(ii + 744);
    const auto *ii_745 = buffer.data(ii + 745);
    const auto *ii_746 = buffer.data(ii + 746);
    const auto *ii_747 = buffer.data(ii + 747);
    const auto *ii_748 = buffer.data(ii + 748);
    const auto *ii_749 = buffer.data(ii + 749);

    const auto *li_600 = buffer.data(li + 600);
    const auto *li_601 = buffer.data(li + 601);
    const auto *li_602 = buffer.data(li + 602);
    const auto *li_603 = buffer.data(li + 603);
    const auto *li_604 = buffer.data(li + 604);
    const auto *li_605 = buffer.data(li + 605);
    const auto *li_606 = buffer.data(li + 606);
    const auto *li_607 = buffer.data(li + 607);
    const auto *li_608 = buffer.data(li + 608);
    const auto *li_609 = buffer.data(li + 609);
    const auto *li_610 = buffer.data(li + 610);
    const auto *li_611 = buffer.data(li + 611);
    const auto *li_612 = buffer.data(li + 612);
    const auto *li_613 = buffer.data(li + 613);
    const auto *li_614 = buffer.data(li + 614);
    const auto *li_615 = buffer.data(li + 615);
    const auto *li_616 = buffer.data(li + 616);
    const auto *li_617 = buffer.data(li + 617);
    const auto *li_618 = buffer.data(li + 618);
    const auto *li_619 = buffer.data(li + 619);
    const auto *li_620 = buffer.data(li + 620);
    const auto *li_621 = buffer.data(li + 621);
    const auto *li_622 = buffer.data(li + 622);
    const auto *li_623 = buffer.data(li + 623);
    const auto *li_624 = buffer.data(li + 624);
    const auto *li_625 = buffer.data(li + 625);
    const auto *li_626 = buffer.data(li + 626);
    const auto *li_627 = buffer.data(li + 627);
    const auto *li_628 = buffer.data(li + 628);
    const auto *li_629 = buffer.data(li + 629);
    const auto *li_630 = buffer.data(li + 630);
    const auto *li_631 = buffer.data(li + 631);
    const auto *li_632 = buffer.data(li + 632);
    const auto *li_633 = buffer.data(li + 633);
    const auto *li_634 = buffer.data(li + 634);
    const auto *li_635 = buffer.data(li + 635);
    const auto *li_636 = buffer.data(li + 636);
    const auto *li_637 = buffer.data(li + 637);
    const auto *li_638 = buffer.data(li + 638);
    const auto *li_639 = buffer.data(li + 639);
    const auto *li_640 = buffer.data(li + 640);
    const auto *li_641 = buffer.data(li + 641);
    const auto *li_642 = buffer.data(li + 642);
    const auto *li_643 = buffer.data(li + 643);
    const auto *li_644 = buffer.data(li + 644);
    const auto *li_645 = buffer.data(li + 645);
    const auto *li_646 = buffer.data(li + 646);
    const auto *li_647 = buffer.data(li + 647);
    const auto *li_648 = buffer.data(li + 648);
    const auto *li_649 = buffer.data(li + 649);
    const auto *li_650 = buffer.data(li + 650);
    const auto *li_651 = buffer.data(li + 651);
    const auto *li_652 = buffer.data(li + 652);
    const auto *li_653 = buffer.data(li + 653);
    const auto *li_654 = buffer.data(li + 654);
    const auto *li_655 = buffer.data(li + 655);
    const auto *li_656 = buffer.data(li + 656);
    const auto *li_657 = buffer.data(li + 657);
    const auto *li_658 = buffer.data(li + 658);
    const auto *li_659 = buffer.data(li + 659);
    const auto *li_660 = buffer.data(li + 660);
    const auto *li_661 = buffer.data(li + 661);
    const auto *li_662 = buffer.data(li + 662);
    const auto *li_663 = buffer.data(li + 663);
    const auto *li_664 = buffer.data(li + 664);
    const auto *li_665 = buffer.data(li + 665);
    const auto *li_666 = buffer.data(li + 666);
    const auto *li_667 = buffer.data(li + 667);
    const auto *li_668 = buffer.data(li + 668);
    const auto *li_669 = buffer.data(li + 669);
    const auto *li_670 = buffer.data(li + 670);
    const auto *li_671 = buffer.data(li + 671);
    const auto *li_672 = buffer.data(li + 672);
    const auto *li_673 = buffer.data(li + 673);
    const auto *li_674 = buffer.data(li + 674);
    const auto *li_675 = buffer.data(li + 675);
    const auto *li_676 = buffer.data(li + 676);
    const auto *li_677 = buffer.data(li + 677);
    const auto *li_678 = buffer.data(li + 678);
    const auto *li_679 = buffer.data(li + 679);
    const auto *li_680 = buffer.data(li + 680);
    const auto *li_681 = buffer.data(li + 681);
    const auto *li_682 = buffer.data(li + 682);
    const auto *li_683 = buffer.data(li + 683);
    const auto *li_684 = buffer.data(li + 684);
    const auto *li_685 = buffer.data(li + 685);
    const auto *li_686 = buffer.data(li + 686);
    const auto *li_687 = buffer.data(li + 687);
    const auto *li_688 = buffer.data(li + 688);
    const auto *li_689 = buffer.data(li + 689);
    const auto *li_690 = buffer.data(li + 690);
    const auto *li_691 = buffer.data(li + 691);
    const auto *li_692 = buffer.data(li + 692);
    const auto *li_693 = buffer.data(li + 693);
    const auto *li_694 = buffer.data(li + 694);
    const auto *li_695 = buffer.data(li + 695);
    const auto *li_696 = buffer.data(li + 696);
    const auto *li_697 = buffer.data(li + 697);
    const auto *li_698 = buffer.data(li + 698);
    const auto *li_699 = buffer.data(li + 699);
    const auto *li_700 = buffer.data(li + 700);
    const auto *li_701 = buffer.data(li + 701);
    const auto *li_702 = buffer.data(li + 702);
    const auto *li_703 = buffer.data(li + 703);
    const auto *li_704 = buffer.data(li + 704);
    const auto *li_705 = buffer.data(li + 705);
    const auto *li_706 = buffer.data(li + 706);
    const auto *li_707 = buffer.data(li + 707);
    const auto *li_708 = buffer.data(li + 708);
    const auto *li_709 = buffer.data(li + 709);
    const auto *li_710 = buffer.data(li + 710);
    const auto *li_711 = buffer.data(li + 711);
    const auto *li_712 = buffer.data(li + 712);
    const auto *li_713 = buffer.data(li + 713);
    const auto *li_714 = buffer.data(li + 714);
    const auto *li_715 = buffer.data(li + 715);
    const auto *li_716 = buffer.data(li + 716);
    const auto *li_717 = buffer.data(li + 717);
    const auto *li_718 = buffer.data(li + 718);
    const auto *li_719 = buffer.data(li + 719);
    const auto *li_720 = buffer.data(li + 720);
    const auto *li_721 = buffer.data(li + 721);
    const auto *li_722 = buffer.data(li + 722);
    const auto *li_723 = buffer.data(li + 723);
    const auto *li_724 = buffer.data(li + 724);
    const auto *li_725 = buffer.data(li + 725);
    const auto *li_726 = buffer.data(li + 726);
    const auto *li_727 = buffer.data(li + 727);
    const auto *li_728 = buffer.data(li + 728);
    const auto *li_729 = buffer.data(li + 729);
    const auto *li_730 = buffer.data(li + 730);
    const auto *li_731 = buffer.data(li + 731);
    const auto *li_732 = buffer.data(li + 732);
    const auto *li_733 = buffer.data(li + 733);
    const auto *li_734 = buffer.data(li + 734);
    const auto *li_735 = buffer.data(li + 735);
    const auto *li_736 = buffer.data(li + 736);
    const auto *li_737 = buffer.data(li + 737);
    const auto *li_738 = buffer.data(li + 738);
    const auto *li_739 = buffer.data(li + 739);
    const auto *li_740 = buffer.data(li + 740);
    const auto *li_741 = buffer.data(li + 741);
    const auto *li_742 = buffer.data(li + 742);
    const auto *li_743 = buffer.data(li + 743);
    const auto *li_744 = buffer.data(li + 744);
    const auto *li_745 = buffer.data(li + 745);
    const auto *li_746 = buffer.data(li + 746);
    const auto *li_747 = buffer.data(li + 747);
    const auto *li_748 = buffer.data(li + 748);
    const auto *li_749 = buffer.data(li + 749);

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, ii_600, ii_601, ii_602, ii_603, \
                         ii_604, li_600, li_601, li_602, li_603, \
                         li_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -ii_600[k]
                   + f_0 * li_600[k];

        t_601[k] = -ii_601[k]
                   + f_0 * li_601[k];

        t_602[k] = -ii_602[k]
                   + f_0 * li_602[k];

        t_603[k] = -ii_603[k]
                   + f_0 * li_603[k];

        t_604[k] = -ii_604[k]
                   + f_0 * li_604[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, ii_605, ii_606, ii_607, ii_608, \
                         ii_609, li_605, li_606, li_607, li_608, \
                         li_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -ii_605[k]
                   + f_0 * li_605[k];

        t_606[k] = -ii_606[k]
                   + f_0 * li_606[k];

        t_607[k] = -ii_607[k]
                   + f_0 * li_607[k];

        t_608[k] = -ii_608[k]
                   + f_0 * li_608[k];

        t_609[k] = -ii_609[k]
                   + f_0 * li_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, ii_610, ii_611, ii_612, ii_613, \
                         ii_614, li_610, li_611, li_612, li_613, \
                         li_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -ii_610[k]
                   + f_0 * li_610[k];

        t_611[k] = -ii_611[k]
                   + f_0 * li_611[k];

        t_612[k] = -ii_612[k]
                   + f_0 * li_612[k];

        t_613[k] = -ii_613[k]
                   + f_0 * li_613[k];

        t_614[k] = -ii_614[k]
                   + f_0 * li_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, ii_615, ii_616, ii_617, ii_618, \
                         ii_619, li_615, li_616, li_617, li_618, \
                         li_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -ii_615[k]
                   + f_0 * li_615[k];

        t_616[k] = -ii_616[k]
                   + f_0 * li_616[k];

        t_617[k] = -ii_617[k]
                   + f_0 * li_617[k];

        t_618[k] = -ii_618[k]
                   + f_0 * li_618[k];

        t_619[k] = -ii_619[k]
                   + f_0 * li_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, ii_620, ii_621, ii_622, ii_623, \
                         ii_624, li_620, li_621, li_622, li_623, \
                         li_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -ii_620[k]
                   + f_0 * li_620[k];

        t_621[k] = -ii_621[k]
                   + f_0 * li_621[k];

        t_622[k] = -ii_622[k]
                   + f_0 * li_622[k];

        t_623[k] = -ii_623[k]
                   + f_0 * li_623[k];

        t_624[k] = -ii_624[k]
                   + f_0 * li_624[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, ii_625, ii_626, ii_627, ii_628, \
                         ii_629, li_625, li_626, li_627, li_628, \
                         li_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -ii_625[k]
                   + f_0 * li_625[k];

        t_626[k] = -ii_626[k]
                   + f_0 * li_626[k];

        t_627[k] = -ii_627[k]
                   + f_0 * li_627[k];

        t_628[k] = -ii_628[k]
                   + f_0 * li_628[k];

        t_629[k] = -ii_629[k]
                   + f_0 * li_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ii_630, ii_631, ii_632, ii_633, \
                         ii_634, li_630, li_631, li_632, li_633, \
                         li_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -ii_630[k]
                   + f_0 * li_630[k];

        t_631[k] = -ii_631[k]
                   + f_0 * li_631[k];

        t_632[k] = -ii_632[k]
                   + f_0 * li_632[k];

        t_633[k] = -ii_633[k]
                   + f_0 * li_633[k];

        t_634[k] = -ii_634[k]
                   + f_0 * li_634[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ii_635, ii_636, ii_637, ii_638, \
                         ii_639, li_635, li_636, li_637, li_638, \
                         li_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -ii_635[k]
                   + f_0 * li_635[k];

        t_636[k] = -ii_636[k]
                   + f_0 * li_636[k];

        t_637[k] = -ii_637[k]
                   + f_0 * li_637[k];

        t_638[k] = -ii_638[k]
                   + f_0 * li_638[k];

        t_639[k] = -ii_639[k]
                   + f_0 * li_639[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ii_640, ii_641, ii_642, ii_643, \
                         ii_644, li_640, li_641, li_642, li_643, \
                         li_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -ii_640[k]
                   + f_0 * li_640[k];

        t_641[k] = -ii_641[k]
                   + f_0 * li_641[k];

        t_642[k] = -ii_642[k]
                   + f_0 * li_642[k];

        t_643[k] = -ii_643[k]
                   + f_0 * li_643[k];

        t_644[k] = -ii_644[k]
                   + f_0 * li_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ii_645, ii_646, ii_647, ii_648, \
                         ii_649, li_645, li_646, li_647, li_648, \
                         li_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -ii_645[k]
                   + f_0 * li_645[k];

        t_646[k] = -ii_646[k]
                   + f_0 * li_646[k];

        t_647[k] = -ii_647[k]
                   + f_0 * li_647[k];

        t_648[k] = -ii_648[k]
                   + f_0 * li_648[k];

        t_649[k] = -ii_649[k]
                   + f_0 * li_649[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, ii_650, ii_651, ii_652, ii_653, \
                         ii_654, li_650, li_651, li_652, li_653, \
                         li_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -ii_650[k]
                   + f_0 * li_650[k];

        t_651[k] = -ii_651[k]
                   + f_0 * li_651[k];

        t_652[k] = -ii_652[k]
                   + f_0 * li_652[k];

        t_653[k] = -ii_653[k]
                   + f_0 * li_653[k];

        t_654[k] = -ii_654[k]
                   + f_0 * li_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, ii_655, ii_656, ii_657, ii_658, \
                         ii_659, li_655, li_656, li_657, li_658, \
                         li_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -ii_655[k]
                   + f_0 * li_655[k];

        t_656[k] = -ii_656[k]
                   + f_0 * li_656[k];

        t_657[k] = -ii_657[k]
                   + f_0 * li_657[k];

        t_658[k] = -ii_658[k]
                   + f_0 * li_658[k];

        t_659[k] = -ii_659[k]
                   + f_0 * li_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, ii_660, ii_661, ii_662, ii_663, \
                         ii_664, li_660, li_661, li_662, li_663, \
                         li_664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -ii_660[k]
                   + f_0 * li_660[k];

        t_661[k] = -ii_661[k]
                   + f_0 * li_661[k];

        t_662[k] = -ii_662[k]
                   + f_0 * li_662[k];

        t_663[k] = -ii_663[k]
                   + f_0 * li_663[k];

        t_664[k] = -ii_664[k]
                   + f_0 * li_664[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, ii_665, ii_666, ii_667, ii_668, \
                         ii_669, li_665, li_666, li_667, li_668, \
                         li_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -ii_665[k]
                   + f_0 * li_665[k];

        t_666[k] = -ii_666[k]
                   + f_0 * li_666[k];

        t_667[k] = -ii_667[k]
                   + f_0 * li_667[k];

        t_668[k] = -ii_668[k]
                   + f_0 * li_668[k];

        t_669[k] = -ii_669[k]
                   + f_0 * li_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, ii_670, ii_671, ii_672, ii_673, \
                         ii_674, li_670, li_671, li_672, li_673, \
                         li_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -ii_670[k]
                   + f_0 * li_670[k];

        t_671[k] = -ii_671[k]
                   + f_0 * li_671[k];

        t_672[k] = -ii_672[k]
                   + f_0 * li_672[k];

        t_673[k] = -ii_673[k]
                   + f_0 * li_673[k];

        t_674[k] = -ii_674[k]
                   + f_0 * li_674[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, ii_675, ii_676, ii_677, ii_678, \
                         ii_679, li_675, li_676, li_677, li_678, \
                         li_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -ii_675[k]
                   + f_0 * li_675[k];

        t_676[k] = -ii_676[k]
                   + f_0 * li_676[k];

        t_677[k] = -ii_677[k]
                   + f_0 * li_677[k];

        t_678[k] = -ii_678[k]
                   + f_0 * li_678[k];

        t_679[k] = -ii_679[k]
                   + f_0 * li_679[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, ii_680, ii_681, ii_682, ii_683, \
                         ii_684, li_680, li_681, li_682, li_683, \
                         li_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -ii_680[k]
                   + f_0 * li_680[k];

        t_681[k] = -ii_681[k]
                   + f_0 * li_681[k];

        t_682[k] = -ii_682[k]
                   + f_0 * li_682[k];

        t_683[k] = -ii_683[k]
                   + f_0 * li_683[k];

        t_684[k] = -ii_684[k]
                   + f_0 * li_684[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, ii_685, ii_686, ii_687, ii_688, \
                         ii_689, li_685, li_686, li_687, li_688, \
                         li_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -ii_685[k]
                   + f_0 * li_685[k];

        t_686[k] = -ii_686[k]
                   + f_0 * li_686[k];

        t_687[k] = -ii_687[k]
                   + f_0 * li_687[k];

        t_688[k] = -ii_688[k]
                   + f_0 * li_688[k];

        t_689[k] = -ii_689[k]
                   + f_0 * li_689[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, ii_690, ii_691, ii_692, ii_693, \
                         ii_694, li_690, li_691, li_692, li_693, \
                         li_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -ii_690[k]
                   + f_0 * li_690[k];

        t_691[k] = -ii_691[k]
                   + f_0 * li_691[k];

        t_692[k] = -ii_692[k]
                   + f_0 * li_692[k];

        t_693[k] = -ii_693[k]
                   + f_0 * li_693[k];

        t_694[k] = -ii_694[k]
                   + f_0 * li_694[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, ii_695, ii_696, ii_697, ii_698, \
                         ii_699, li_695, li_696, li_697, li_698, \
                         li_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -ii_695[k]
                   + f_0 * li_695[k];

        t_696[k] = -ii_696[k]
                   + f_0 * li_696[k];

        t_697[k] = -ii_697[k]
                   + f_0 * li_697[k];

        t_698[k] = -ii_698[k]
                   + f_0 * li_698[k];

        t_699[k] = -ii_699[k]
                   + f_0 * li_699[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, ii_700, ii_701, ii_702, ii_703, \
                         ii_704, li_700, li_701, li_702, li_703, \
                         li_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -ii_700[k]
                   + f_0 * li_700[k];

        t_701[k] = -ii_701[k]
                   + f_0 * li_701[k];

        t_702[k] = -ii_702[k]
                   + f_0 * li_702[k];

        t_703[k] = -ii_703[k]
                   + f_0 * li_703[k];

        t_704[k] = -ii_704[k]
                   + f_0 * li_704[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, ii_705, ii_706, ii_707, ii_708, \
                         ii_709, li_705, li_706, li_707, li_708, \
                         li_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -ii_705[k]
                   + f_0 * li_705[k];

        t_706[k] = -ii_706[k]
                   + f_0 * li_706[k];

        t_707[k] = -ii_707[k]
                   + f_0 * li_707[k];

        t_708[k] = -ii_708[k]
                   + f_0 * li_708[k];

        t_709[k] = -ii_709[k]
                   + f_0 * li_709[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, ii_710, ii_711, ii_712, ii_713, \
                         ii_714, li_710, li_711, li_712, li_713, \
                         li_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -ii_710[k]
                   + f_0 * li_710[k];

        t_711[k] = -ii_711[k]
                   + f_0 * li_711[k];

        t_712[k] = -ii_712[k]
                   + f_0 * li_712[k];

        t_713[k] = -ii_713[k]
                   + f_0 * li_713[k];

        t_714[k] = -ii_714[k]
                   + f_0 * li_714[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, ii_715, ii_716, ii_717, ii_718, \
                         ii_719, li_715, li_716, li_717, li_718, \
                         li_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -ii_715[k]
                   + f_0 * li_715[k];

        t_716[k] = -ii_716[k]
                   + f_0 * li_716[k];

        t_717[k] = -ii_717[k]
                   + f_0 * li_717[k];

        t_718[k] = -ii_718[k]
                   + f_0 * li_718[k];

        t_719[k] = -ii_719[k]
                   + f_0 * li_719[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, ii_720, ii_721, ii_722, ii_723, \
                         ii_724, li_720, li_721, li_722, li_723, \
                         li_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -ii_720[k]
                   + f_0 * li_720[k];

        t_721[k] = -ii_721[k]
                   + f_0 * li_721[k];

        t_722[k] = -ii_722[k]
                   + f_0 * li_722[k];

        t_723[k] = -ii_723[k]
                   + f_0 * li_723[k];

        t_724[k] = -ii_724[k]
                   + f_0 * li_724[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, ii_725, ii_726, ii_727, ii_728, \
                         ii_729, li_725, li_726, li_727, li_728, \
                         li_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -ii_725[k]
                   + f_0 * li_725[k];

        t_726[k] = -ii_726[k]
                   + f_0 * li_726[k];

        t_727[k] = -ii_727[k]
                   + f_0 * li_727[k];

        t_728[k] = -ii_728[k]
                   + f_0 * li_728[k];

        t_729[k] = -ii_729[k]
                   + f_0 * li_729[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, ii_730, ii_731, ii_732, ii_733, \
                         ii_734, li_730, li_731, li_732, li_733, \
                         li_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -ii_730[k]
                   + f_0 * li_730[k];

        t_731[k] = -ii_731[k]
                   + f_0 * li_731[k];

        t_732[k] = -ii_732[k]
                   + f_0 * li_732[k];

        t_733[k] = -ii_733[k]
                   + f_0 * li_733[k];

        t_734[k] = -ii_734[k]
                   + f_0 * li_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, ii_735, ii_736, ii_737, ii_738, \
                         ii_739, li_735, li_736, li_737, li_738, \
                         li_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -ii_735[k]
                   + f_0 * li_735[k];

        t_736[k] = -ii_736[k]
                   + f_0 * li_736[k];

        t_737[k] = -ii_737[k]
                   + f_0 * li_737[k];

        t_738[k] = -ii_738[k]
                   + f_0 * li_738[k];

        t_739[k] = -ii_739[k]
                   + f_0 * li_739[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, ii_740, ii_741, ii_742, ii_743, \
                         ii_744, li_740, li_741, li_742, li_743, \
                         li_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -ii_740[k]
                   + f_0 * li_740[k];

        t_741[k] = -ii_741[k]
                   + f_0 * li_741[k];

        t_742[k] = -ii_742[k]
                   + f_0 * li_742[k];

        t_743[k] = -ii_743[k]
                   + f_0 * li_743[k];

        t_744[k] = -ii_744[k]
                   + f_0 * li_744[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, ii_745, ii_746, ii_747, ii_748, \
                         ii_749, li_745, li_746, li_747, li_748, \
                         li_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -ii_745[k]
                   + f_0 * li_745[k];

        t_746[k] = -ii_746[k]
                   + f_0 * li_746[k];

        t_747[k] = -ii_747[k]
                   + f_0 * li_747[k];

        t_748[k] = -ii_748[k]
                   + f_0 * li_748[k];

        t_749[k] = -ii_749[k]
                   + f_0 * li_749[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_750 = buffer.data(ii + 750);
    const auto *ii_751 = buffer.data(ii + 751);
    const auto *ii_752 = buffer.data(ii + 752);
    const auto *ii_753 = buffer.data(ii + 753);
    const auto *ii_754 = buffer.data(ii + 754);
    const auto *ii_755 = buffer.data(ii + 755);
    const auto *ii_756 = buffer.data(ii + 756);
    const auto *ii_757 = buffer.data(ii + 757);
    const auto *ii_758 = buffer.data(ii + 758);
    const auto *ii_759 = buffer.data(ii + 759);
    const auto *ii_760 = buffer.data(ii + 760);
    const auto *ii_761 = buffer.data(ii + 761);
    const auto *ii_762 = buffer.data(ii + 762);
    const auto *ii_763 = buffer.data(ii + 763);
    const auto *ii_764 = buffer.data(ii + 764);
    const auto *ii_765 = buffer.data(ii + 765);
    const auto *ii_766 = buffer.data(ii + 766);
    const auto *ii_767 = buffer.data(ii + 767);
    const auto *ii_768 = buffer.data(ii + 768);
    const auto *ii_769 = buffer.data(ii + 769);
    const auto *ii_770 = buffer.data(ii + 770);
    const auto *ii_771 = buffer.data(ii + 771);
    const auto *ii_772 = buffer.data(ii + 772);
    const auto *ii_773 = buffer.data(ii + 773);
    const auto *ii_774 = buffer.data(ii + 774);
    const auto *ii_775 = buffer.data(ii + 775);
    const auto *ii_776 = buffer.data(ii + 776);
    const auto *ii_777 = buffer.data(ii + 777);
    const auto *ii_778 = buffer.data(ii + 778);
    const auto *ii_779 = buffer.data(ii + 779);
    const auto *ii_780 = buffer.data(ii + 780);
    const auto *ii_781 = buffer.data(ii + 781);
    const auto *ii_782 = buffer.data(ii + 782);
    const auto *ii_783 = buffer.data(ii + 783);

    const auto *li_750 = buffer.data(li + 750);
    const auto *li_751 = buffer.data(li + 751);
    const auto *li_752 = buffer.data(li + 752);
    const auto *li_753 = buffer.data(li + 753);
    const auto *li_754 = buffer.data(li + 754);
    const auto *li_755 = buffer.data(li + 755);
    const auto *li_756 = buffer.data(li + 756);
    const auto *li_757 = buffer.data(li + 757);
    const auto *li_758 = buffer.data(li + 758);
    const auto *li_759 = buffer.data(li + 759);
    const auto *li_760 = buffer.data(li + 760);
    const auto *li_761 = buffer.data(li + 761);
    const auto *li_762 = buffer.data(li + 762);
    const auto *li_763 = buffer.data(li + 763);
    const auto *li_764 = buffer.data(li + 764);
    const auto *li_765 = buffer.data(li + 765);
    const auto *li_766 = buffer.data(li + 766);
    const auto *li_767 = buffer.data(li + 767);
    const auto *li_768 = buffer.data(li + 768);
    const auto *li_769 = buffer.data(li + 769);
    const auto *li_770 = buffer.data(li + 770);
    const auto *li_771 = buffer.data(li + 771);
    const auto *li_772 = buffer.data(li + 772);
    const auto *li_773 = buffer.data(li + 773);
    const auto *li_774 = buffer.data(li + 774);
    const auto *li_775 = buffer.data(li + 775);
    const auto *li_776 = buffer.data(li + 776);
    const auto *li_777 = buffer.data(li + 777);
    const auto *li_778 = buffer.data(li + 778);
    const auto *li_779 = buffer.data(li + 779);
    const auto *li_780 = buffer.data(li + 780);
    const auto *li_781 = buffer.data(li + 781);
    const auto *li_782 = buffer.data(li + 782);
    const auto *li_783 = buffer.data(li + 783);
    const auto *li_784 = buffer.data(li + 784);
    const auto *li_785 = buffer.data(li + 785);
    const auto *li_786 = buffer.data(li + 786);
    const auto *li_787 = buffer.data(li + 787);
    const auto *li_788 = buffer.data(li + 788);
    const auto *li_789 = buffer.data(li + 789);
    const auto *li_790 = buffer.data(li + 790);
    const auto *li_791 = buffer.data(li + 791);
    const auto *li_792 = buffer.data(li + 792);
    const auto *li_793 = buffer.data(li + 793);
    const auto *li_794 = buffer.data(li + 794);
    const auto *li_795 = buffer.data(li + 795);
    const auto *li_796 = buffer.data(li + 796);
    const auto *li_797 = buffer.data(li + 797);
    const auto *li_798 = buffer.data(li + 798);
    const auto *li_799 = buffer.data(li + 799);
    const auto *li_800 = buffer.data(li + 800);
    const auto *li_801 = buffer.data(li + 801);
    const auto *li_802 = buffer.data(li + 802);
    const auto *li_803 = buffer.data(li + 803);
    const auto *li_804 = buffer.data(li + 804);
    const auto *li_805 = buffer.data(li + 805);
    const auto *li_806 = buffer.data(li + 806);
    const auto *li_807 = buffer.data(li + 807);
    const auto *li_808 = buffer.data(li + 808);
    const auto *li_809 = buffer.data(li + 809);
    const auto *li_810 = buffer.data(li + 810);
    const auto *li_811 = buffer.data(li + 811);
    const auto *li_812 = buffer.data(li + 812);
    const auto *li_813 = buffer.data(li + 813);
    const auto *li_814 = buffer.data(li + 814);
    const auto *li_815 = buffer.data(li + 815);
    const auto *li_816 = buffer.data(li + 816);
    const auto *li_817 = buffer.data(li + 817);
    const auto *li_818 = buffer.data(li + 818);
    const auto *li_819 = buffer.data(li + 819);
    const auto *li_820 = buffer.data(li + 820);
    const auto *li_821 = buffer.data(li + 821);
    const auto *li_822 = buffer.data(li + 822);
    const auto *li_823 = buffer.data(li + 823);
    const auto *li_824 = buffer.data(li + 824);
    const auto *li_825 = buffer.data(li + 825);
    const auto *li_826 = buffer.data(li + 826);
    const auto *li_827 = buffer.data(li + 827);
    const auto *li_828 = buffer.data(li + 828);
    const auto *li_829 = buffer.data(li + 829);
    const auto *li_830 = buffer.data(li + 830);
    const auto *li_831 = buffer.data(li + 831);
    const auto *li_832 = buffer.data(li + 832);
    const auto *li_833 = buffer.data(li + 833);
    const auto *li_834 = buffer.data(li + 834);
    const auto *li_835 = buffer.data(li + 835);
    const auto *li_836 = buffer.data(li + 836);
    const auto *li_837 = buffer.data(li + 837);
    const auto *li_838 = buffer.data(li + 838);
    const auto *li_839 = buffer.data(li + 839);
    const auto *li_840 = buffer.data(li + 840);
    const auto *li_841 = buffer.data(li + 841);
    const auto *li_842 = buffer.data(li + 842);
    const auto *li_843 = buffer.data(li + 843);
    const auto *li_844 = buffer.data(li + 844);
    const auto *li_845 = buffer.data(li + 845);
    const auto *li_846 = buffer.data(li + 846);
    const auto *li_847 = buffer.data(li + 847);
    const auto *li_848 = buffer.data(li + 848);
    const auto *li_849 = buffer.data(li + 849);
    const auto *li_850 = buffer.data(li + 850);
    const auto *li_851 = buffer.data(li + 851);
    const auto *li_852 = buffer.data(li + 852);
    const auto *li_853 = buffer.data(li + 853);
    const auto *li_854 = buffer.data(li + 854);
    const auto *li_855 = buffer.data(li + 855);
    const auto *li_856 = buffer.data(li + 856);
    const auto *li_857 = buffer.data(li + 857);
    const auto *li_858 = buffer.data(li + 858);
    const auto *li_859 = buffer.data(li + 859);
    const auto *li_860 = buffer.data(li + 860);
    const auto *li_861 = buffer.data(li + 861);
    const auto *li_862 = buffer.data(li + 862);
    const auto *li_863 = buffer.data(li + 863);
    const auto *li_864 = buffer.data(li + 864);
    const auto *li_865 = buffer.data(li + 865);
    const auto *li_866 = buffer.data(li + 866);
    const auto *li_867 = buffer.data(li + 867);
    const auto *li_868 = buffer.data(li + 868);
    const auto *li_869 = buffer.data(li + 869);
    const auto *li_870 = buffer.data(li + 870);
    const auto *li_871 = buffer.data(li + 871);
    const auto *li_872 = buffer.data(li + 872);
    const auto *li_873 = buffer.data(li + 873);
    const auto *li_874 = buffer.data(li + 874);
    const auto *li_875 = buffer.data(li + 875);
    const auto *li_876 = buffer.data(li + 876);
    const auto *li_877 = buffer.data(li + 877);
    const auto *li_878 = buffer.data(li + 878);
    const auto *li_879 = buffer.data(li + 879);
    const auto *li_880 = buffer.data(li + 880);
    const auto *li_881 = buffer.data(li + 881);
    const auto *li_882 = buffer.data(li + 882);
    const auto *li_883 = buffer.data(li + 883);
    const auto *li_884 = buffer.data(li + 884);
    const auto *li_885 = buffer.data(li + 885);
    const auto *li_886 = buffer.data(li + 886);
    const auto *li_887 = buffer.data(li + 887);
    const auto *li_888 = buffer.data(li + 888);
    const auto *li_889 = buffer.data(li + 889);
    const auto *li_890 = buffer.data(li + 890);
    const auto *li_891 = buffer.data(li + 891);
    const auto *li_892 = buffer.data(li + 892);
    const auto *li_893 = buffer.data(li + 893);
    const auto *li_894 = buffer.data(li + 894);
    const auto *li_895 = buffer.data(li + 895);
    const auto *li_896 = buffer.data(li + 896);
    const auto *li_897 = buffer.data(li + 897);
    const auto *li_898 = buffer.data(li + 898);
    const auto *li_899 = buffer.data(li + 899);
    const auto *li_900 = buffer.data(li + 900);
    const auto *li_901 = buffer.data(li + 901);
    const auto *li_902 = buffer.data(li + 902);
    const auto *li_903 = buffer.data(li + 903);
    const auto *li_904 = buffer.data(li + 904);
    const auto *li_905 = buffer.data(li + 905);
    const auto *li_906 = buffer.data(li + 906);
    const auto *li_907 = buffer.data(li + 907);
    const auto *li_908 = buffer.data(li + 908);
    const auto *li_909 = buffer.data(li + 909);
    const auto *li_910 = buffer.data(li + 910);
    const auto *li_911 = buffer.data(li + 911);
    const auto *li_912 = buffer.data(li + 912);
    const auto *li_913 = buffer.data(li + 913);
    const auto *li_914 = buffer.data(li + 914);
    const auto *li_915 = buffer.data(li + 915);
    const auto *li_916 = buffer.data(li + 916);
    const auto *li_917 = buffer.data(li + 917);
    const auto *li_918 = buffer.data(li + 918);
    const auto *li_919 = buffer.data(li + 919);
    const auto *li_920 = buffer.data(li + 920);
    const auto *li_921 = buffer.data(li + 921);
    const auto *li_922 = buffer.data(li + 922);
    const auto *li_923 = buffer.data(li + 923);
    const auto *li_924 = buffer.data(li + 924);
    const auto *li_925 = buffer.data(li + 925);
    const auto *li_926 = buffer.data(li + 926);
    const auto *li_927 = buffer.data(li + 927);
    const auto *li_928 = buffer.data(li + 928);
    const auto *li_929 = buffer.data(li + 929);
    const auto *li_930 = buffer.data(li + 930);
    const auto *li_931 = buffer.data(li + 931);
    const auto *li_932 = buffer.data(li + 932);
    const auto *li_933 = buffer.data(li + 933);
    const auto *li_934 = buffer.data(li + 934);
    const auto *li_935 = buffer.data(li + 935);
    const auto *li_936 = buffer.data(li + 936);
    const auto *li_937 = buffer.data(li + 937);
    const auto *li_938 = buffer.data(li + 938);
    const auto *li_939 = buffer.data(li + 939);
    const auto *li_940 = buffer.data(li + 940);
    const auto *li_941 = buffer.data(li + 941);
    const auto *li_942 = buffer.data(li + 942);
    const auto *li_943 = buffer.data(li + 943);
    const auto *li_944 = buffer.data(li + 944);
    const auto *li_945 = buffer.data(li + 945);
    const auto *li_946 = buffer.data(li + 946);
    const auto *li_947 = buffer.data(li + 947);
    const auto *li_948 = buffer.data(li + 948);
    const auto *li_949 = buffer.data(li + 949);
    const auto *li_950 = buffer.data(li + 950);
    const auto *li_951 = buffer.data(li + 951);
    const auto *li_952 = buffer.data(li + 952);
    const auto *li_953 = buffer.data(li + 953);

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, ii_750, ii_751, ii_752, ii_753, \
                         ii_754, li_750, li_751, li_752, li_753, \
                         li_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -ii_750[k]
                   + f_0 * li_750[k];

        t_751[k] = -ii_751[k]
                   + f_0 * li_751[k];

        t_752[k] = -ii_752[k]
                   + f_0 * li_752[k];

        t_753[k] = -ii_753[k]
                   + f_0 * li_753[k];

        t_754[k] = -ii_754[k]
                   + f_0 * li_754[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, ii_755, ii_756, ii_757, ii_758, \
                         ii_759, li_755, li_756, li_757, li_758, \
                         li_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -ii_755[k]
                   + f_0 * li_755[k];

        t_756[k] = -ii_756[k]
                   + f_0 * li_756[k];

        t_757[k] = -ii_757[k]
                   + f_0 * li_757[k];

        t_758[k] = -ii_758[k]
                   + f_0 * li_758[k];

        t_759[k] = -ii_759[k]
                   + f_0 * li_759[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, t_764, ii_760, ii_761, ii_762, ii_763, \
                         ii_764, li_760, li_761, li_762, li_763, \
                         li_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = -ii_760[k]
                   + f_0 * li_760[k];

        t_761[k] = -ii_761[k]
                   + f_0 * li_761[k];

        t_762[k] = -ii_762[k]
                   + f_0 * li_762[k];

        t_763[k] = -ii_763[k]
                   + f_0 * li_763[k];

        t_764[k] = -ii_764[k]
                   + f_0 * li_764[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, t_769, ii_765, ii_766, ii_767, ii_768, \
                         ii_769, li_765, li_766, li_767, li_768, \
                         li_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = -ii_765[k]
                   + f_0 * li_765[k];

        t_766[k] = -ii_766[k]
                   + f_0 * li_766[k];

        t_767[k] = -ii_767[k]
                   + f_0 * li_767[k];

        t_768[k] = -ii_768[k]
                   + f_0 * li_768[k];

        t_769[k] = -ii_769[k]
                   + f_0 * li_769[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, ii_770, ii_771, ii_772, ii_773, \
                         ii_774, li_770, li_771, li_772, li_773, \
                         li_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = -ii_770[k]
                   + f_0 * li_770[k];

        t_771[k] = -ii_771[k]
                   + f_0 * li_771[k];

        t_772[k] = -ii_772[k]
                   + f_0 * li_772[k];

        t_773[k] = -ii_773[k]
                   + f_0 * li_773[k];

        t_774[k] = -ii_774[k]
                   + f_0 * li_774[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, ii_775, ii_776, ii_777, ii_778, \
                         ii_779, li_775, li_776, li_777, li_778, \
                         li_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = -ii_775[k]
                   + f_0 * li_775[k];

        t_776[k] = -ii_776[k]
                   + f_0 * li_776[k];

        t_777[k] = -ii_777[k]
                   + f_0 * li_777[k];

        t_778[k] = -ii_778[k]
                   + f_0 * li_778[k];

        t_779[k] = -ii_779[k]
                   + f_0 * li_779[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, t_785, ii_780, ii_781, ii_782, \
                         ii_783, li_780, li_781, li_782, li_783, li_784, \
                         li_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = -ii_780[k]
                   + f_0 * li_780[k];

        t_781[k] = -ii_781[k]
                   + f_0 * li_781[k];

        t_782[k] = -ii_782[k]
                   + f_0 * li_782[k];

        t_783[k] = -ii_783[k]
                   + f_0 * li_783[k];

        t_784[k] = f_0 * li_784[k];

        t_785[k] = f_0 * li_785[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, t_790, t_791, t_792, t_793, li_786, \
                         li_787, li_788, li_789, li_790, li_791, li_792, \
                         li_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_0 * li_786[k];

        t_787[k] = f_0 * li_787[k];

        t_788[k] = f_0 * li_788[k];

        t_789[k] = f_0 * li_789[k];

        t_790[k] = f_0 * li_790[k];

        t_791[k] = f_0 * li_791[k];

        t_792[k] = f_0 * li_792[k];

        t_793[k] = f_0 * li_793[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, t_798, t_799, t_800, t_801, li_794, \
                         li_795, li_796, li_797, li_798, li_799, li_800, \
                         li_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_0 * li_794[k];

        t_795[k] = f_0 * li_795[k];

        t_796[k] = f_0 * li_796[k];

        t_797[k] = f_0 * li_797[k];

        t_798[k] = f_0 * li_798[k];

        t_799[k] = f_0 * li_799[k];

        t_800[k] = f_0 * li_800[k];

        t_801[k] = f_0 * li_801[k];
    }

#pragma omp simd aligned(t_802, t_803, t_804, t_805, t_806, t_807, t_808, t_809, li_802, \
                         li_803, li_804, li_805, li_806, li_807, li_808, \
                         li_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_802[k] = f_0 * li_802[k];

        t_803[k] = f_0 * li_803[k];

        t_804[k] = f_0 * li_804[k];

        t_805[k] = f_0 * li_805[k];

        t_806[k] = f_0 * li_806[k];

        t_807[k] = f_0 * li_807[k];

        t_808[k] = f_0 * li_808[k];

        t_809[k] = f_0 * li_809[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, t_815, t_816, t_817, li_810, \
                         li_811, li_812, li_813, li_814, li_815, li_816, \
                         li_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_0 * li_810[k];

        t_811[k] = f_0 * li_811[k];

        t_812[k] = f_0 * li_812[k];

        t_813[k] = f_0 * li_813[k];

        t_814[k] = f_0 * li_814[k];

        t_815[k] = f_0 * li_815[k];

        t_816[k] = f_0 * li_816[k];

        t_817[k] = f_0 * li_817[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, t_822, t_823, t_824, t_825, li_818, \
                         li_819, li_820, li_821, li_822, li_823, li_824, \
                         li_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_0 * li_818[k];

        t_819[k] = f_0 * li_819[k];

        t_820[k] = f_0 * li_820[k];

        t_821[k] = f_0 * li_821[k];

        t_822[k] = f_0 * li_822[k];

        t_823[k] = f_0 * li_823[k];

        t_824[k] = f_0 * li_824[k];

        t_825[k] = f_0 * li_825[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, t_830, t_831, t_832, t_833, li_826, \
                         li_827, li_828, li_829, li_830, li_831, li_832, \
                         li_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_0 * li_826[k];

        t_827[k] = f_0 * li_827[k];

        t_828[k] = f_0 * li_828[k];

        t_829[k] = f_0 * li_829[k];

        t_830[k] = f_0 * li_830[k];

        t_831[k] = f_0 * li_831[k];

        t_832[k] = f_0 * li_832[k];

        t_833[k] = f_0 * li_833[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, t_837, t_838, t_839, t_840, t_841, li_834, \
                         li_835, li_836, li_837, li_838, li_839, li_840, \
                         li_841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_0 * li_834[k];

        t_835[k] = f_0 * li_835[k];

        t_836[k] = f_0 * li_836[k];

        t_837[k] = f_0 * li_837[k];

        t_838[k] = f_0 * li_838[k];

        t_839[k] = f_0 * li_839[k];

        t_840[k] = f_0 * li_840[k];

        t_841[k] = f_0 * li_841[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, t_845, t_846, t_847, t_848, t_849, li_842, \
                         li_843, li_844, li_845, li_846, li_847, li_848, \
                         li_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_0 * li_842[k];

        t_843[k] = f_0 * li_843[k];

        t_844[k] = f_0 * li_844[k];

        t_845[k] = f_0 * li_845[k];

        t_846[k] = f_0 * li_846[k];

        t_847[k] = f_0 * li_847[k];

        t_848[k] = f_0 * li_848[k];

        t_849[k] = f_0 * li_849[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, t_855, t_856, t_857, li_850, \
                         li_851, li_852, li_853, li_854, li_855, li_856, \
                         li_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = f_0 * li_850[k];

        t_851[k] = f_0 * li_851[k];

        t_852[k] = f_0 * li_852[k];

        t_853[k] = f_0 * li_853[k];

        t_854[k] = f_0 * li_854[k];

        t_855[k] = f_0 * li_855[k];

        t_856[k] = f_0 * li_856[k];

        t_857[k] = f_0 * li_857[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, t_862, t_863, t_864, t_865, li_858, \
                         li_859, li_860, li_861, li_862, li_863, li_864, \
                         li_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_0 * li_858[k];

        t_859[k] = f_0 * li_859[k];

        t_860[k] = f_0 * li_860[k];

        t_861[k] = f_0 * li_861[k];

        t_862[k] = f_0 * li_862[k];

        t_863[k] = f_0 * li_863[k];

        t_864[k] = f_0 * li_864[k];

        t_865[k] = f_0 * li_865[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, t_870, t_871, t_872, t_873, li_866, \
                         li_867, li_868, li_869, li_870, li_871, li_872, \
                         li_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_0 * li_866[k];

        t_867[k] = f_0 * li_867[k];

        t_868[k] = f_0 * li_868[k];

        t_869[k] = f_0 * li_869[k];

        t_870[k] = f_0 * li_870[k];

        t_871[k] = f_0 * li_871[k];

        t_872[k] = f_0 * li_872[k];

        t_873[k] = f_0 * li_873[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, t_878, t_879, t_880, t_881, li_874, \
                         li_875, li_876, li_877, li_878, li_879, li_880, \
                         li_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_0 * li_874[k];

        t_875[k] = f_0 * li_875[k];

        t_876[k] = f_0 * li_876[k];

        t_877[k] = f_0 * li_877[k];

        t_878[k] = f_0 * li_878[k];

        t_879[k] = f_0 * li_879[k];

        t_880[k] = f_0 * li_880[k];

        t_881[k] = f_0 * li_881[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, t_886, t_887, t_888, t_889, li_882, \
                         li_883, li_884, li_885, li_886, li_887, li_888, \
                         li_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_0 * li_882[k];

        t_883[k] = f_0 * li_883[k];

        t_884[k] = f_0 * li_884[k];

        t_885[k] = f_0 * li_885[k];

        t_886[k] = f_0 * li_886[k];

        t_887[k] = f_0 * li_887[k];

        t_888[k] = f_0 * li_888[k];

        t_889[k] = f_0 * li_889[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, t_895, t_896, t_897, li_890, \
                         li_891, li_892, li_893, li_894, li_895, li_896, \
                         li_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_0 * li_890[k];

        t_891[k] = f_0 * li_891[k];

        t_892[k] = f_0 * li_892[k];

        t_893[k] = f_0 * li_893[k];

        t_894[k] = f_0 * li_894[k];

        t_895[k] = f_0 * li_895[k];

        t_896[k] = f_0 * li_896[k];

        t_897[k] = f_0 * li_897[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, t_902, t_903, t_904, t_905, li_898, \
                         li_899, li_900, li_901, li_902, li_903, li_904, \
                         li_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_0 * li_898[k];

        t_899[k] = f_0 * li_899[k];

        t_900[k] = f_0 * li_900[k];

        t_901[k] = f_0 * li_901[k];

        t_902[k] = f_0 * li_902[k];

        t_903[k] = f_0 * li_903[k];

        t_904[k] = f_0 * li_904[k];

        t_905[k] = f_0 * li_905[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, t_909, t_910, t_911, t_912, t_913, li_906, \
                         li_907, li_908, li_909, li_910, li_911, li_912, \
                         li_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_0 * li_906[k];

        t_907[k] = f_0 * li_907[k];

        t_908[k] = f_0 * li_908[k];

        t_909[k] = f_0 * li_909[k];

        t_910[k] = f_0 * li_910[k];

        t_911[k] = f_0 * li_911[k];

        t_912[k] = f_0 * li_912[k];

        t_913[k] = f_0 * li_913[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, t_919, t_920, t_921, li_914, \
                         li_915, li_916, li_917, li_918, li_919, li_920, \
                         li_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_0 * li_914[k];

        t_915[k] = f_0 * li_915[k];

        t_916[k] = f_0 * li_916[k];

        t_917[k] = f_0 * li_917[k];

        t_918[k] = f_0 * li_918[k];

        t_919[k] = f_0 * li_919[k];

        t_920[k] = f_0 * li_920[k];

        t_921[k] = f_0 * li_921[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, t_927, t_928, t_929, li_922, \
                         li_923, li_924, li_925, li_926, li_927, li_928, \
                         li_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_0 * li_922[k];

        t_923[k] = f_0 * li_923[k];

        t_924[k] = f_0 * li_924[k];

        t_925[k] = f_0 * li_925[k];

        t_926[k] = f_0 * li_926[k];

        t_927[k] = f_0 * li_927[k];

        t_928[k] = f_0 * li_928[k];

        t_929[k] = f_0 * li_929[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, t_935, t_936, t_937, li_930, \
                         li_931, li_932, li_933, li_934, li_935, li_936, \
                         li_937 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_0 * li_930[k];

        t_931[k] = f_0 * li_931[k];

        t_932[k] = f_0 * li_932[k];

        t_933[k] = f_0 * li_933[k];

        t_934[k] = f_0 * li_934[k];

        t_935[k] = f_0 * li_935[k];

        t_936[k] = f_0 * li_936[k];

        t_937[k] = f_0 * li_937[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, t_942, t_943, t_944, t_945, li_938, \
                         li_939, li_940, li_941, li_942, li_943, li_944, \
                         li_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_0 * li_938[k];

        t_939[k] = f_0 * li_939[k];

        t_940[k] = f_0 * li_940[k];

        t_941[k] = f_0 * li_941[k];

        t_942[k] = f_0 * li_942[k];

        t_943[k] = f_0 * li_943[k];

        t_944[k] = f_0 * li_944[k];

        t_945[k] = f_0 * li_945[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, t_949, t_950, t_951, t_952, t_953, li_946, \
                         li_947, li_948, li_949, li_950, li_951, li_952, \
                         li_953 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_0 * li_946[k];

        t_947[k] = f_0 * li_947[k];

        t_948[k] = f_0 * li_948[k];

        t_949[k] = f_0 * li_949[k];

        t_950[k] = f_0 * li_950[k];

        t_951[k] = f_0 * li_951[k];

        t_952[k] = f_0 * li_952[k];

        t_953[k] = f_0 * li_953[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t li, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *li_954 = buffer.data(li + 954);
    const auto *li_955 = buffer.data(li + 955);
    const auto *li_956 = buffer.data(li + 956);
    const auto *li_957 = buffer.data(li + 957);
    const auto *li_958 = buffer.data(li + 958);
    const auto *li_959 = buffer.data(li + 959);
    const auto *li_960 = buffer.data(li + 960);
    const auto *li_961 = buffer.data(li + 961);
    const auto *li_962 = buffer.data(li + 962);
    const auto *li_963 = buffer.data(li + 963);
    const auto *li_964 = buffer.data(li + 964);
    const auto *li_965 = buffer.data(li + 965);
    const auto *li_966 = buffer.data(li + 966);
    const auto *li_967 = buffer.data(li + 967);
    const auto *li_968 = buffer.data(li + 968);
    const auto *li_969 = buffer.data(li + 969);
    const auto *li_970 = buffer.data(li + 970);
    const auto *li_971 = buffer.data(li + 971);
    const auto *li_972 = buffer.data(li + 972);
    const auto *li_973 = buffer.data(li + 973);
    const auto *li_974 = buffer.data(li + 974);
    const auto *li_975 = buffer.data(li + 975);
    const auto *li_976 = buffer.data(li + 976);
    const auto *li_977 = buffer.data(li + 977);
    const auto *li_978 = buffer.data(li + 978);
    const auto *li_979 = buffer.data(li + 979);
    const auto *li_980 = buffer.data(li + 980);
    const auto *li_981 = buffer.data(li + 981);
    const auto *li_982 = buffer.data(li + 982);
    const auto *li_983 = buffer.data(li + 983);
    const auto *li_984 = buffer.data(li + 984);
    const auto *li_985 = buffer.data(li + 985);
    const auto *li_986 = buffer.data(li + 986);
    const auto *li_987 = buffer.data(li + 987);
    const auto *li_988 = buffer.data(li + 988);
    const auto *li_989 = buffer.data(li + 989);
    const auto *li_990 = buffer.data(li + 990);
    const auto *li_991 = buffer.data(li + 991);
    const auto *li_992 = buffer.data(li + 992);
    const auto *li_993 = buffer.data(li + 993);
    const auto *li_994 = buffer.data(li + 994);
    const auto *li_995 = buffer.data(li + 995);
    const auto *li_996 = buffer.data(li + 996);
    const auto *li_997 = buffer.data(li + 997);
    const auto *li_998 = buffer.data(li + 998);
    const auto *li_999 = buffer.data(li + 999);
    const auto *li_1000 = buffer.data(li + 1000);
    const auto *li_1001 = buffer.data(li + 1001);
    const auto *li_1002 = buffer.data(li + 1002);
    const auto *li_1003 = buffer.data(li + 1003);
    const auto *li_1004 = buffer.data(li + 1004);
    const auto *li_1005 = buffer.data(li + 1005);
    const auto *li_1006 = buffer.data(li + 1006);
    const auto *li_1007 = buffer.data(li + 1007);

#pragma omp simd aligned(t_954, t_955, t_956, t_957, t_958, t_959, t_960, t_961, li_954, \
                         li_955, li_956, li_957, li_958, li_959, li_960, \
                         li_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_0 * li_954[k];

        t_955[k] = f_0 * li_955[k];

        t_956[k] = f_0 * li_956[k];

        t_957[k] = f_0 * li_957[k];

        t_958[k] = f_0 * li_958[k];

        t_959[k] = f_0 * li_959[k];

        t_960[k] = f_0 * li_960[k];

        t_961[k] = f_0 * li_961[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, t_966, t_967, t_968, t_969, li_962, \
                         li_963, li_964, li_965, li_966, li_967, li_968, \
                         li_969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_0 * li_962[k];

        t_963[k] = f_0 * li_963[k];

        t_964[k] = f_0 * li_964[k];

        t_965[k] = f_0 * li_965[k];

        t_966[k] = f_0 * li_966[k];

        t_967[k] = f_0 * li_967[k];

        t_968[k] = f_0 * li_968[k];

        t_969[k] = f_0 * li_969[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, t_975, t_976, t_977, li_970, \
                         li_971, li_972, li_973, li_974, li_975, li_976, \
                         li_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = f_0 * li_970[k];

        t_971[k] = f_0 * li_971[k];

        t_972[k] = f_0 * li_972[k];

        t_973[k] = f_0 * li_973[k];

        t_974[k] = f_0 * li_974[k];

        t_975[k] = f_0 * li_975[k];

        t_976[k] = f_0 * li_976[k];

        t_977[k] = f_0 * li_977[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, t_982, t_983, t_984, t_985, li_978, \
                         li_979, li_980, li_981, li_982, li_983, li_984, \
                         li_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_0 * li_978[k];

        t_979[k] = f_0 * li_979[k];

        t_980[k] = f_0 * li_980[k];

        t_981[k] = f_0 * li_981[k];

        t_982[k] = f_0 * li_982[k];

        t_983[k] = f_0 * li_983[k];

        t_984[k] = f_0 * li_984[k];

        t_985[k] = f_0 * li_985[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, t_989, t_990, t_991, t_992, t_993, li_986, \
                         li_987, li_988, li_989, li_990, li_991, li_992, \
                         li_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_0 * li_986[k];

        t_987[k] = f_0 * li_987[k];

        t_988[k] = f_0 * li_988[k];

        t_989[k] = f_0 * li_989[k];

        t_990[k] = f_0 * li_990[k];

        t_991[k] = f_0 * li_991[k];

        t_992[k] = f_0 * li_992[k];

        t_993[k] = f_0 * li_993[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, t_997, t_998, t_999, t_1000, t_1001, li_994, \
                         li_995, li_996, li_997, li_998, li_999, li_1000, \
                         li_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = f_0 * li_994[k];

        t_995[k] = f_0 * li_995[k];

        t_996[k] = f_0 * li_996[k];

        t_997[k] = f_0 * li_997[k];

        t_998[k] = f_0 * li_998[k];

        t_999[k] = f_0 * li_999[k];

        t_1000[k] = f_0 * li_1000[k];

        t_1001[k] = f_0 * li_1001[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, t_1005, t_1006, t_1007, li_1002, li_1003, \
                         li_1004, li_1005, li_1006, li_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_0 * li_1002[k];

        t_1003[k] = f_0 * li_1003[k];

        t_1004[k] = f_0 * li_1004[k];

        t_1005[k] = f_0 * li_1005[k];

        t_1006[k] = f_0 * li_1006[k];

        t_1007[k] = f_0 * li_1007[k];
    }
}

auto
compute_prim_geom_10_ki_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ii, const size_t li,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ki_electron_repulsion_0_piece0(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_0_piece1(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_0_piece2(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_0_piece3(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_0_piece4(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_0_piece5(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_0_piece6(buffer, target, li, ncols, alpha);
}

static auto
compute_prim_geom_10_ki_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_0 = buffer.data(ii + 0);
    const auto *ii_1 = buffer.data(ii + 1);
    const auto *ii_2 = buffer.data(ii + 2);
    const auto *ii_3 = buffer.data(ii + 3);
    const auto *ii_4 = buffer.data(ii + 4);
    const auto *ii_5 = buffer.data(ii + 5);
    const auto *ii_6 = buffer.data(ii + 6);
    const auto *ii_7 = buffer.data(ii + 7);
    const auto *ii_8 = buffer.data(ii + 8);
    const auto *ii_9 = buffer.data(ii + 9);
    const auto *ii_10 = buffer.data(ii + 10);
    const auto *ii_11 = buffer.data(ii + 11);
    const auto *ii_12 = buffer.data(ii + 12);
    const auto *ii_13 = buffer.data(ii + 13);
    const auto *ii_14 = buffer.data(ii + 14);
    const auto *ii_15 = buffer.data(ii + 15);
    const auto *ii_16 = buffer.data(ii + 16);
    const auto *ii_17 = buffer.data(ii + 17);
    const auto *ii_18 = buffer.data(ii + 18);
    const auto *ii_19 = buffer.data(ii + 19);
    const auto *ii_20 = buffer.data(ii + 20);
    const auto *ii_21 = buffer.data(ii + 21);
    const auto *ii_22 = buffer.data(ii + 22);
    const auto *ii_23 = buffer.data(ii + 23);
    const auto *ii_24 = buffer.data(ii + 24);
    const auto *ii_25 = buffer.data(ii + 25);
    const auto *ii_26 = buffer.data(ii + 26);
    const auto *ii_27 = buffer.data(ii + 27);
    const auto *ii_28 = buffer.data(ii + 28);
    const auto *ii_29 = buffer.data(ii + 29);
    const auto *ii_30 = buffer.data(ii + 30);
    const auto *ii_31 = buffer.data(ii + 31);
    const auto *ii_32 = buffer.data(ii + 32);
    const auto *ii_33 = buffer.data(ii + 33);
    const auto *ii_34 = buffer.data(ii + 34);
    const auto *ii_35 = buffer.data(ii + 35);
    const auto *ii_36 = buffer.data(ii + 36);
    const auto *ii_37 = buffer.data(ii + 37);
    const auto *ii_38 = buffer.data(ii + 38);
    const auto *ii_39 = buffer.data(ii + 39);
    const auto *ii_40 = buffer.data(ii + 40);
    const auto *ii_41 = buffer.data(ii + 41);
    const auto *ii_42 = buffer.data(ii + 42);
    const auto *ii_43 = buffer.data(ii + 43);
    const auto *ii_44 = buffer.data(ii + 44);
    const auto *ii_45 = buffer.data(ii + 45);
    const auto *ii_46 = buffer.data(ii + 46);
    const auto *ii_47 = buffer.data(ii + 47);
    const auto *ii_48 = buffer.data(ii + 48);
    const auto *ii_49 = buffer.data(ii + 49);
    const auto *ii_50 = buffer.data(ii + 50);
    const auto *ii_51 = buffer.data(ii + 51);
    const auto *ii_52 = buffer.data(ii + 52);
    const auto *ii_53 = buffer.data(ii + 53);
    const auto *ii_54 = buffer.data(ii + 54);
    const auto *ii_55 = buffer.data(ii + 55);
    const auto *ii_56 = buffer.data(ii + 56);
    const auto *ii_57 = buffer.data(ii + 57);
    const auto *ii_58 = buffer.data(ii + 58);
    const auto *ii_59 = buffer.data(ii + 59);
    const auto *ii_60 = buffer.data(ii + 60);
    const auto *ii_61 = buffer.data(ii + 61);
    const auto *ii_62 = buffer.data(ii + 62);
    const auto *ii_63 = buffer.data(ii + 63);
    const auto *ii_64 = buffer.data(ii + 64);
    const auto *ii_65 = buffer.data(ii + 65);
    const auto *ii_66 = buffer.data(ii + 66);
    const auto *ii_67 = buffer.data(ii + 67);
    const auto *ii_68 = buffer.data(ii + 68);
    const auto *ii_69 = buffer.data(ii + 69);
    const auto *ii_70 = buffer.data(ii + 70);
    const auto *ii_71 = buffer.data(ii + 71);
    const auto *ii_72 = buffer.data(ii + 72);
    const auto *ii_73 = buffer.data(ii + 73);
    const auto *ii_74 = buffer.data(ii + 74);
    const auto *ii_75 = buffer.data(ii + 75);
    const auto *ii_76 = buffer.data(ii + 76);
    const auto *ii_77 = buffer.data(ii + 77);
    const auto *ii_78 = buffer.data(ii + 78);
    const auto *ii_79 = buffer.data(ii + 79);
    const auto *ii_80 = buffer.data(ii + 80);
    const auto *ii_81 = buffer.data(ii + 81);
    const auto *ii_82 = buffer.data(ii + 82);
    const auto *ii_83 = buffer.data(ii + 83);
    const auto *ii_84 = buffer.data(ii + 84);
    const auto *ii_85 = buffer.data(ii + 85);
    const auto *ii_86 = buffer.data(ii + 86);
    const auto *ii_87 = buffer.data(ii + 87);
    const auto *ii_88 = buffer.data(ii + 88);
    const auto *ii_89 = buffer.data(ii + 89);
    const auto *ii_90 = buffer.data(ii + 90);

    const auto *li_28 = buffer.data(li + 28);
    const auto *li_29 = buffer.data(li + 29);
    const auto *li_30 = buffer.data(li + 30);
    const auto *li_31 = buffer.data(li + 31);
    const auto *li_32 = buffer.data(li + 32);
    const auto *li_33 = buffer.data(li + 33);
    const auto *li_34 = buffer.data(li + 34);
    const auto *li_35 = buffer.data(li + 35);
    const auto *li_36 = buffer.data(li + 36);
    const auto *li_37 = buffer.data(li + 37);
    const auto *li_38 = buffer.data(li + 38);
    const auto *li_39 = buffer.data(li + 39);
    const auto *li_40 = buffer.data(li + 40);
    const auto *li_41 = buffer.data(li + 41);
    const auto *li_42 = buffer.data(li + 42);
    const auto *li_43 = buffer.data(li + 43);
    const auto *li_44 = buffer.data(li + 44);
    const auto *li_45 = buffer.data(li + 45);
    const auto *li_46 = buffer.data(li + 46);
    const auto *li_47 = buffer.data(li + 47);
    const auto *li_48 = buffer.data(li + 48);
    const auto *li_49 = buffer.data(li + 49);
    const auto *li_50 = buffer.data(li + 50);
    const auto *li_51 = buffer.data(li + 51);
    const auto *li_52 = buffer.data(li + 52);
    const auto *li_53 = buffer.data(li + 53);
    const auto *li_54 = buffer.data(li + 54);
    const auto *li_55 = buffer.data(li + 55);
    const auto *li_84 = buffer.data(li + 84);
    const auto *li_85 = buffer.data(li + 85);
    const auto *li_86 = buffer.data(li + 86);
    const auto *li_87 = buffer.data(li + 87);
    const auto *li_88 = buffer.data(li + 88);
    const auto *li_89 = buffer.data(li + 89);
    const auto *li_90 = buffer.data(li + 90);
    const auto *li_91 = buffer.data(li + 91);
    const auto *li_92 = buffer.data(li + 92);
    const auto *li_93 = buffer.data(li + 93);
    const auto *li_94 = buffer.data(li + 94);
    const auto *li_95 = buffer.data(li + 95);
    const auto *li_96 = buffer.data(li + 96);
    const auto *li_97 = buffer.data(li + 97);
    const auto *li_98 = buffer.data(li + 98);
    const auto *li_99 = buffer.data(li + 99);
    const auto *li_100 = buffer.data(li + 100);
    const auto *li_101 = buffer.data(li + 101);
    const auto *li_102 = buffer.data(li + 102);
    const auto *li_103 = buffer.data(li + 103);
    const auto *li_104 = buffer.data(li + 104);
    const auto *li_105 = buffer.data(li + 105);
    const auto *li_106 = buffer.data(li + 106);
    const auto *li_107 = buffer.data(li + 107);
    const auto *li_108 = buffer.data(li + 108);
    const auto *li_109 = buffer.data(li + 109);
    const auto *li_110 = buffer.data(li + 110);
    const auto *li_111 = buffer.data(li + 111);
    const auto *li_112 = buffer.data(li + 112);
    const auto *li_113 = buffer.data(li + 113);
    const auto *li_114 = buffer.data(li + 114);
    const auto *li_115 = buffer.data(li + 115);
    const auto *li_116 = buffer.data(li + 116);
    const auto *li_117 = buffer.data(li + 117);
    const auto *li_118 = buffer.data(li + 118);
    const auto *li_119 = buffer.data(li + 119);
    const auto *li_120 = buffer.data(li + 120);
    const auto *li_121 = buffer.data(li + 121);
    const auto *li_122 = buffer.data(li + 122);
    const auto *li_123 = buffer.data(li + 123);
    const auto *li_124 = buffer.data(li + 124);
    const auto *li_125 = buffer.data(li + 125);
    const auto *li_126 = buffer.data(li + 126);
    const auto *li_127 = buffer.data(li + 127);
    const auto *li_128 = buffer.data(li + 128);
    const auto *li_129 = buffer.data(li + 129);
    const auto *li_130 = buffer.data(li + 130);
    const auto *li_131 = buffer.data(li + 131);
    const auto *li_132 = buffer.data(li + 132);
    const auto *li_133 = buffer.data(li + 133);
    const auto *li_134 = buffer.data(li + 134);
    const auto *li_135 = buffer.data(li + 135);
    const auto *li_136 = buffer.data(li + 136);
    const auto *li_137 = buffer.data(li + 137);
    const auto *li_138 = buffer.data(li + 138);
    const auto *li_139 = buffer.data(li + 139);
    const auto *li_168 = buffer.data(li + 168);
    const auto *li_169 = buffer.data(li + 169);
    const auto *li_170 = buffer.data(li + 170);
    const auto *li_171 = buffer.data(li + 171);
    const auto *li_172 = buffer.data(li + 172);
    const auto *li_173 = buffer.data(li + 173);
    const auto *li_174 = buffer.data(li + 174);
    const auto *li_175 = buffer.data(li + 175);
    const auto *li_176 = buffer.data(li + 176);
    const auto *li_177 = buffer.data(li + 177);
    const auto *li_178 = buffer.data(li + 178);
    const auto *li_179 = buffer.data(li + 179);
    const auto *li_180 = buffer.data(li + 180);
    const auto *li_181 = buffer.data(li + 181);
    const auto *li_182 = buffer.data(li + 182);
    const auto *li_183 = buffer.data(li + 183);
    const auto *li_184 = buffer.data(li + 184);
    const auto *li_185 = buffer.data(li + 185);
    const auto *li_186 = buffer.data(li + 186);
    const auto *li_187 = buffer.data(li + 187);
    const auto *li_188 = buffer.data(li + 188);
    const auto *li_189 = buffer.data(li + 189);
    const auto *li_190 = buffer.data(li + 190);
    const auto *li_191 = buffer.data(li + 191);
    const auto *li_192 = buffer.data(li + 192);
    const auto *li_193 = buffer.data(li + 193);
    const auto *li_194 = buffer.data(li + 194);
    const auto *li_195 = buffer.data(li + 195);
    const auto *li_196 = buffer.data(li + 196);
    const auto *li_197 = buffer.data(li + 197);
    const auto *li_198 = buffer.data(li + 198);
    const auto *li_199 = buffer.data(li + 199);
    const auto *li_200 = buffer.data(li + 200);
    const auto *li_201 = buffer.data(li + 201);
    const auto *li_202 = buffer.data(li + 202);
    const auto *li_203 = buffer.data(li + 203);
    const auto *li_204 = buffer.data(li + 204);
    const auto *li_205 = buffer.data(li + 205);
    const auto *li_206 = buffer.data(li + 206);
    const auto *li_207 = buffer.data(li + 207);
    const auto *li_208 = buffer.data(li + 208);
    const auto *li_209 = buffer.data(li + 209);
    const auto *li_210 = buffer.data(li + 210);
    const auto *li_211 = buffer.data(li + 211);
    const auto *li_212 = buffer.data(li + 212);
    const auto *li_213 = buffer.data(li + 213);
    const auto *li_214 = buffer.data(li + 214);
    const auto *li_215 = buffer.data(li + 215);
    const auto *li_216 = buffer.data(li + 216);
    const auto *li_217 = buffer.data(li + 217);
    const auto *li_218 = buffer.data(li + 218);
    const auto *li_219 = buffer.data(li + 219);
    const auto *li_220 = buffer.data(li + 220);
    const auto *li_221 = buffer.data(li + 221);
    const auto *li_222 = buffer.data(li + 222);
    const auto *li_223 = buffer.data(li + 223);
    const auto *li_224 = buffer.data(li + 224);
    const auto *li_225 = buffer.data(li + 225);
    const auto *li_226 = buffer.data(li + 226);
    const auto *li_227 = buffer.data(li + 227);
    const auto *li_228 = buffer.data(li + 228);
    const auto *li_229 = buffer.data(li + 229);
    const auto *li_230 = buffer.data(li + 230);
    const auto *li_231 = buffer.data(li + 231);
    const auto *li_232 = buffer.data(li + 232);
    const auto *li_233 = buffer.data(li + 233);
    const auto *li_234 = buffer.data(li + 234);
    const auto *li_235 = buffer.data(li + 235);
    const auto *li_236 = buffer.data(li + 236);
    const auto *li_237 = buffer.data(li + 237);
    const auto *li_238 = buffer.data(li + 238);
    const auto *li_239 = buffer.data(li + 239);
    const auto *li_240 = buffer.data(li + 240);
    const auto *li_241 = buffer.data(li + 241);
    const auto *li_242 = buffer.data(li + 242);
    const auto *li_243 = buffer.data(li + 243);
    const auto *li_244 = buffer.data(li + 244);
    const auto *li_245 = buffer.data(li + 245);
    const auto *li_246 = buffer.data(li + 246);
    const auto *li_247 = buffer.data(li + 247);
    const auto *li_248 = buffer.data(li + 248);
    const auto *li_249 = buffer.data(li + 249);
    const auto *li_250 = buffer.data(li + 250);
    const auto *li_251 = buffer.data(li + 251);
    const auto *li_280 = buffer.data(li + 280);
    const auto *li_281 = buffer.data(li + 281);
    const auto *li_282 = buffer.data(li + 282);
    const auto *li_283 = buffer.data(li + 283);
    const auto *li_284 = buffer.data(li + 284);
    const auto *li_285 = buffer.data(li + 285);
    const auto *li_286 = buffer.data(li + 286);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, li_28, li_29, li_30, li_31, \
                         li_32, li_33, li_34, li_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * li_28[k];

        t_1[k] = f_0 * li_29[k];

        t_2[k] = f_0 * li_30[k];

        t_3[k] = f_0 * li_31[k];

        t_4[k] = f_0 * li_32[k];

        t_5[k] = f_0 * li_33[k];

        t_6[k] = f_0 * li_34[k];

        t_7[k] = f_0 * li_35[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, li_36, li_37, li_38, \
                         li_39, li_40, li_41, li_42, li_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * li_36[k];

        t_9[k] = f_0 * li_37[k];

        t_10[k] = f_0 * li_38[k];

        t_11[k] = f_0 * li_39[k];

        t_12[k] = f_0 * li_40[k];

        t_13[k] = f_0 * li_41[k];

        t_14[k] = f_0 * li_42[k];

        t_15[k] = f_0 * li_43[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, li_44, li_45, li_46, \
                         li_47, li_48, li_49, li_50, li_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * li_44[k];

        t_17[k] = f_0 * li_45[k];

        t_18[k] = f_0 * li_46[k];

        t_19[k] = f_0 * li_47[k];

        t_20[k] = f_0 * li_48[k];

        t_21[k] = f_0 * li_49[k];

        t_22[k] = f_0 * li_50[k];

        t_23[k] = f_0 * li_51[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, ii_0, ii_1, li_52, li_53, li_54, \
                         li_55, li_84, li_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * li_52[k];

        t_25[k] = f_0 * li_53[k];

        t_26[k] = f_0 * li_54[k];

        t_27[k] = f_0 * li_55[k];

        t_28[k] = -ii_0[k]
                  + f_0 * li_84[k];

        t_29[k] = -ii_1[k]
                  + f_0 * li_85[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ii_2, ii_3, ii_4, ii_5, ii_6, li_86, \
                         li_87, li_88, li_89, li_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -ii_2[k]
                  + f_0 * li_86[k];

        t_31[k] = -ii_3[k]
                  + f_0 * li_87[k];

        t_32[k] = -ii_4[k]
                  + f_0 * li_88[k];

        t_33[k] = -ii_5[k]
                  + f_0 * li_89[k];

        t_34[k] = -ii_6[k]
                  + f_0 * li_90[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ii_7, ii_8, ii_9, ii_10, ii_11, li_91, \
                         li_92, li_93, li_94, li_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -ii_7[k]
                  + f_0 * li_91[k];

        t_36[k] = -ii_8[k]
                  + f_0 * li_92[k];

        t_37[k] = -ii_9[k]
                  + f_0 * li_93[k];

        t_38[k] = -ii_10[k]
                  + f_0 * li_94[k];

        t_39[k] = -ii_11[k]
                  + f_0 * li_95[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ii_12, ii_13, ii_14, ii_15, ii_16, \
                         li_96, li_97, li_98, li_99, li_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -ii_12[k]
                  + f_0 * li_96[k];

        t_41[k] = -ii_13[k]
                  + f_0 * li_97[k];

        t_42[k] = -ii_14[k]
                  + f_0 * li_98[k];

        t_43[k] = -ii_15[k]
                  + f_0 * li_99[k];

        t_44[k] = -ii_16[k]
                  + f_0 * li_100[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ii_17, ii_18, ii_19, ii_20, ii_21, \
                         li_101, li_102, li_103, li_104, li_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -ii_17[k]
                  + f_0 * li_101[k];

        t_46[k] = -ii_18[k]
                  + f_0 * li_102[k];

        t_47[k] = -ii_19[k]
                  + f_0 * li_103[k];

        t_48[k] = -ii_20[k]
                  + f_0 * li_104[k];

        t_49[k] = -ii_21[k]
                  + f_0 * li_105[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ii_22, ii_23, ii_24, ii_25, ii_26, \
                         li_106, li_107, li_108, li_109, li_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -ii_22[k]
                  + f_0 * li_106[k];

        t_51[k] = -ii_23[k]
                  + f_0 * li_107[k];

        t_52[k] = -ii_24[k]
                  + f_0 * li_108[k];

        t_53[k] = -ii_25[k]
                  + f_0 * li_109[k];

        t_54[k] = -ii_26[k]
                  + f_0 * li_110[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, ii_27, li_111, li_112, \
                         li_113, li_114, li_115, li_116, li_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -ii_27[k]
                  + f_0 * li_111[k];

        t_56[k] = f_0 * li_112[k];

        t_57[k] = f_0 * li_113[k];

        t_58[k] = f_0 * li_114[k];

        t_59[k] = f_0 * li_115[k];

        t_60[k] = f_0 * li_116[k];

        t_61[k] = f_0 * li_117[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, t_68, t_69, li_118, li_119, \
                         li_120, li_121, li_122, li_123, li_124, \
                         li_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * li_118[k];

        t_63[k] = f_0 * li_119[k];

        t_64[k] = f_0 * li_120[k];

        t_65[k] = f_0 * li_121[k];

        t_66[k] = f_0 * li_122[k];

        t_67[k] = f_0 * li_123[k];

        t_68[k] = f_0 * li_124[k];

        t_69[k] = f_0 * li_125[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, t_76, t_77, li_126, li_127, \
                         li_128, li_129, li_130, li_131, li_132, \
                         li_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * li_126[k];

        t_71[k] = f_0 * li_127[k];

        t_72[k] = f_0 * li_128[k];

        t_73[k] = f_0 * li_129[k];

        t_74[k] = f_0 * li_130[k];

        t_75[k] = f_0 * li_131[k];

        t_76[k] = f_0 * li_132[k];

        t_77[k] = f_0 * li_133[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, t_84, ii_28, li_134, li_135, \
                         li_136, li_137, li_138, li_139, li_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * li_134[k];

        t_79[k] = f_0 * li_135[k];

        t_80[k] = f_0 * li_136[k];

        t_81[k] = f_0 * li_137[k];

        t_82[k] = f_0 * li_138[k];

        t_83[k] = f_0 * li_139[k];

        t_84[k] = -2.0 * ii_28[k]
                  + f_0 * li_168[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ii_29, ii_30, ii_31, ii_32, ii_33, \
                         li_169, li_170, li_171, li_172, li_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * ii_29[k]
                  + f_0 * li_169[k];

        t_86[k] = -2.0 * ii_30[k]
                  + f_0 * li_170[k];

        t_87[k] = -2.0 * ii_31[k]
                  + f_0 * li_171[k];

        t_88[k] = -2.0 * ii_32[k]
                  + f_0 * li_172[k];

        t_89[k] = -2.0 * ii_33[k]
                  + f_0 * li_173[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ii_34, ii_35, ii_36, ii_37, ii_38, \
                         li_174, li_175, li_176, li_177, li_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * ii_34[k]
                  + f_0 * li_174[k];

        t_91[k] = -2.0 * ii_35[k]
                  + f_0 * li_175[k];

        t_92[k] = -2.0 * ii_36[k]
                  + f_0 * li_176[k];

        t_93[k] = -2.0 * ii_37[k]
                  + f_0 * li_177[k];

        t_94[k] = -2.0 * ii_38[k]
                  + f_0 * li_178[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ii_39, ii_40, ii_41, ii_42, ii_43, \
                         li_179, li_180, li_181, li_182, li_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * ii_39[k]
                  + f_0 * li_179[k];

        t_96[k] = -2.0 * ii_40[k]
                  + f_0 * li_180[k];

        t_97[k] = -2.0 * ii_41[k]
                  + f_0 * li_181[k];

        t_98[k] = -2.0 * ii_42[k]
                  + f_0 * li_182[k];

        t_99[k] = -2.0 * ii_43[k]
                  + f_0 * li_183[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ii_44, ii_45, ii_46, ii_47, ii_48, \
                         li_184, li_185, li_186, li_187, li_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * ii_44[k]
                   + f_0 * li_184[k];

        t_101[k] = -2.0 * ii_45[k]
                   + f_0 * li_185[k];

        t_102[k] = -2.0 * ii_46[k]
                   + f_0 * li_186[k];

        t_103[k] = -2.0 * ii_47[k]
                   + f_0 * li_187[k];

        t_104[k] = -2.0 * ii_48[k]
                   + f_0 * li_188[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ii_49, ii_50, ii_51, ii_52, ii_53, \
                         li_189, li_190, li_191, li_192, li_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * ii_49[k]
                   + f_0 * li_189[k];

        t_106[k] = -2.0 * ii_50[k]
                   + f_0 * li_190[k];

        t_107[k] = -2.0 * ii_51[k]
                   + f_0 * li_191[k];

        t_108[k] = -2.0 * ii_52[k]
                   + f_0 * li_192[k];

        t_109[k] = -2.0 * ii_53[k]
                   + f_0 * li_193[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ii_54, ii_55, ii_56, ii_57, ii_58, \
                         li_194, li_195, li_196, li_197, li_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * ii_54[k]
                   + f_0 * li_194[k];

        t_111[k] = -2.0 * ii_55[k]
                   + f_0 * li_195[k];

        t_112[k] = -ii_56[k]
                   + f_0 * li_196[k];

        t_113[k] = -ii_57[k]
                   + f_0 * li_197[k];

        t_114[k] = -ii_58[k]
                   + f_0 * li_198[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ii_59, ii_60, ii_61, ii_62, ii_63, \
                         li_199, li_200, li_201, li_202, li_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -ii_59[k]
                   + f_0 * li_199[k];

        t_116[k] = -ii_60[k]
                   + f_0 * li_200[k];

        t_117[k] = -ii_61[k]
                   + f_0 * li_201[k];

        t_118[k] = -ii_62[k]
                   + f_0 * li_202[k];

        t_119[k] = -ii_63[k]
                   + f_0 * li_203[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ii_64, ii_65, ii_66, ii_67, ii_68, \
                         li_204, li_205, li_206, li_207, li_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -ii_64[k]
                   + f_0 * li_204[k];

        t_121[k] = -ii_65[k]
                   + f_0 * li_205[k];

        t_122[k] = -ii_66[k]
                   + f_0 * li_206[k];

        t_123[k] = -ii_67[k]
                   + f_0 * li_207[k];

        t_124[k] = -ii_68[k]
                   + f_0 * li_208[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ii_69, ii_70, ii_71, ii_72, ii_73, \
                         li_209, li_210, li_211, li_212, li_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -ii_69[k]
                   + f_0 * li_209[k];

        t_126[k] = -ii_70[k]
                   + f_0 * li_210[k];

        t_127[k] = -ii_71[k]
                   + f_0 * li_211[k];

        t_128[k] = -ii_72[k]
                   + f_0 * li_212[k];

        t_129[k] = -ii_73[k]
                   + f_0 * li_213[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ii_74, ii_75, ii_76, ii_77, ii_78, \
                         li_214, li_215, li_216, li_217, li_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -ii_74[k]
                   + f_0 * li_214[k];

        t_131[k] = -ii_75[k]
                   + f_0 * li_215[k];

        t_132[k] = -ii_76[k]
                   + f_0 * li_216[k];

        t_133[k] = -ii_77[k]
                   + f_0 * li_217[k];

        t_134[k] = -ii_78[k]
                   + f_0 * li_218[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ii_79, ii_80, ii_81, ii_82, ii_83, \
                         li_219, li_220, li_221, li_222, li_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -ii_79[k]
                   + f_0 * li_219[k];

        t_136[k] = -ii_80[k]
                   + f_0 * li_220[k];

        t_137[k] = -ii_81[k]
                   + f_0 * li_221[k];

        t_138[k] = -ii_82[k]
                   + f_0 * li_222[k];

        t_139[k] = -ii_83[k]
                   + f_0 * li_223[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, li_224, \
                         li_225, li_226, li_227, li_228, li_229, li_230, \
                         li_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * li_224[k];

        t_141[k] = f_0 * li_225[k];

        t_142[k] = f_0 * li_226[k];

        t_143[k] = f_0 * li_227[k];

        t_144[k] = f_0 * li_228[k];

        t_145[k] = f_0 * li_229[k];

        t_146[k] = f_0 * li_230[k];

        t_147[k] = f_0 * li_231[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, t_154, t_155, li_232, \
                         li_233, li_234, li_235, li_236, li_237, li_238, \
                         li_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * li_232[k];

        t_149[k] = f_0 * li_233[k];

        t_150[k] = f_0 * li_234[k];

        t_151[k] = f_0 * li_235[k];

        t_152[k] = f_0 * li_236[k];

        t_153[k] = f_0 * li_237[k];

        t_154[k] = f_0 * li_238[k];

        t_155[k] = f_0 * li_239[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, t_162, t_163, li_240, \
                         li_241, li_242, li_243, li_244, li_245, li_246, \
                         li_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * li_240[k];

        t_157[k] = f_0 * li_241[k];

        t_158[k] = f_0 * li_242[k];

        t_159[k] = f_0 * li_243[k];

        t_160[k] = f_0 * li_244[k];

        t_161[k] = f_0 * li_245[k];

        t_162[k] = f_0 * li_246[k];

        t_163[k] = f_0 * li_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, ii_84, ii_85, li_248, \
                         li_249, li_250, li_251, li_280, li_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * li_248[k];

        t_165[k] = f_0 * li_249[k];

        t_166[k] = f_0 * li_250[k];

        t_167[k] = f_0 * li_251[k];

        t_168[k] = -3.0 * ii_84[k]
                   + f_0 * li_280[k];

        t_169[k] = -3.0 * ii_85[k]
                   + f_0 * li_281[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ii_86, ii_87, ii_88, ii_89, ii_90, \
                         li_282, li_283, li_284, li_285, li_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * ii_86[k]
                   + f_0 * li_282[k];

        t_171[k] = -3.0 * ii_87[k]
                   + f_0 * li_283[k];

        t_172[k] = -3.0 * ii_88[k]
                   + f_0 * li_284[k];

        t_173[k] = -3.0 * ii_89[k]
                   + f_0 * li_285[k];

        t_174[k] = -3.0 * ii_90[k]
                   + f_0 * li_286[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_91 = buffer.data(ii + 91);
    const auto *ii_92 = buffer.data(ii + 92);
    const auto *ii_93 = buffer.data(ii + 93);
    const auto *ii_94 = buffer.data(ii + 94);
    const auto *ii_95 = buffer.data(ii + 95);
    const auto *ii_96 = buffer.data(ii + 96);
    const auto *ii_97 = buffer.data(ii + 97);
    const auto *ii_98 = buffer.data(ii + 98);
    const auto *ii_99 = buffer.data(ii + 99);
    const auto *ii_100 = buffer.data(ii + 100);
    const auto *ii_101 = buffer.data(ii + 101);
    const auto *ii_102 = buffer.data(ii + 102);
    const auto *ii_103 = buffer.data(ii + 103);
    const auto *ii_104 = buffer.data(ii + 104);
    const auto *ii_105 = buffer.data(ii + 105);
    const auto *ii_106 = buffer.data(ii + 106);
    const auto *ii_107 = buffer.data(ii + 107);
    const auto *ii_108 = buffer.data(ii + 108);
    const auto *ii_109 = buffer.data(ii + 109);
    const auto *ii_110 = buffer.data(ii + 110);
    const auto *ii_111 = buffer.data(ii + 111);
    const auto *ii_112 = buffer.data(ii + 112);
    const auto *ii_113 = buffer.data(ii + 113);
    const auto *ii_114 = buffer.data(ii + 114);
    const auto *ii_115 = buffer.data(ii + 115);
    const auto *ii_116 = buffer.data(ii + 116);
    const auto *ii_117 = buffer.data(ii + 117);
    const auto *ii_118 = buffer.data(ii + 118);
    const auto *ii_119 = buffer.data(ii + 119);
    const auto *ii_120 = buffer.data(ii + 120);
    const auto *ii_121 = buffer.data(ii + 121);
    const auto *ii_122 = buffer.data(ii + 122);
    const auto *ii_123 = buffer.data(ii + 123);
    const auto *ii_124 = buffer.data(ii + 124);
    const auto *ii_125 = buffer.data(ii + 125);
    const auto *ii_126 = buffer.data(ii + 126);
    const auto *ii_127 = buffer.data(ii + 127);
    const auto *ii_128 = buffer.data(ii + 128);
    const auto *ii_129 = buffer.data(ii + 129);
    const auto *ii_130 = buffer.data(ii + 130);
    const auto *ii_131 = buffer.data(ii + 131);
    const auto *ii_132 = buffer.data(ii + 132);
    const auto *ii_133 = buffer.data(ii + 133);
    const auto *ii_134 = buffer.data(ii + 134);
    const auto *ii_135 = buffer.data(ii + 135);
    const auto *ii_136 = buffer.data(ii + 136);
    const auto *ii_137 = buffer.data(ii + 137);
    const auto *ii_138 = buffer.data(ii + 138);
    const auto *ii_139 = buffer.data(ii + 139);
    const auto *ii_140 = buffer.data(ii + 140);
    const auto *ii_141 = buffer.data(ii + 141);
    const auto *ii_142 = buffer.data(ii + 142);
    const auto *ii_143 = buffer.data(ii + 143);
    const auto *ii_144 = buffer.data(ii + 144);
    const auto *ii_145 = buffer.data(ii + 145);
    const auto *ii_146 = buffer.data(ii + 146);
    const auto *ii_147 = buffer.data(ii + 147);
    const auto *ii_148 = buffer.data(ii + 148);
    const auto *ii_149 = buffer.data(ii + 149);
    const auto *ii_150 = buffer.data(ii + 150);
    const auto *ii_151 = buffer.data(ii + 151);
    const auto *ii_152 = buffer.data(ii + 152);
    const auto *ii_153 = buffer.data(ii + 153);
    const auto *ii_154 = buffer.data(ii + 154);
    const auto *ii_155 = buffer.data(ii + 155);
    const auto *ii_156 = buffer.data(ii + 156);
    const auto *ii_157 = buffer.data(ii + 157);
    const auto *ii_158 = buffer.data(ii + 158);
    const auto *ii_159 = buffer.data(ii + 159);
    const auto *ii_160 = buffer.data(ii + 160);
    const auto *ii_161 = buffer.data(ii + 161);
    const auto *ii_162 = buffer.data(ii + 162);
    const auto *ii_163 = buffer.data(ii + 163);
    const auto *ii_164 = buffer.data(ii + 164);
    const auto *ii_165 = buffer.data(ii + 165);
    const auto *ii_166 = buffer.data(ii + 166);
    const auto *ii_167 = buffer.data(ii + 167);
    const auto *ii_168 = buffer.data(ii + 168);
    const auto *ii_169 = buffer.data(ii + 169);
    const auto *ii_170 = buffer.data(ii + 170);
    const auto *ii_171 = buffer.data(ii + 171);
    const auto *ii_172 = buffer.data(ii + 172);
    const auto *ii_173 = buffer.data(ii + 173);
    const auto *ii_174 = buffer.data(ii + 174);
    const auto *ii_175 = buffer.data(ii + 175);
    const auto *ii_176 = buffer.data(ii + 176);
    const auto *ii_177 = buffer.data(ii + 177);
    const auto *ii_178 = buffer.data(ii + 178);
    const auto *ii_179 = buffer.data(ii + 179);
    const auto *ii_180 = buffer.data(ii + 180);
    const auto *ii_181 = buffer.data(ii + 181);
    const auto *ii_182 = buffer.data(ii + 182);
    const auto *ii_183 = buffer.data(ii + 183);
    const auto *ii_184 = buffer.data(ii + 184);
    const auto *ii_185 = buffer.data(ii + 185);
    const auto *ii_186 = buffer.data(ii + 186);
    const auto *ii_187 = buffer.data(ii + 187);
    const auto *ii_188 = buffer.data(ii + 188);
    const auto *ii_189 = buffer.data(ii + 189);
    const auto *ii_190 = buffer.data(ii + 190);
    const auto *ii_191 = buffer.data(ii + 191);
    const auto *ii_192 = buffer.data(ii + 192);
    const auto *ii_193 = buffer.data(ii + 193);
    const auto *ii_194 = buffer.data(ii + 194);
    const auto *ii_195 = buffer.data(ii + 195);
    const auto *ii_196 = buffer.data(ii + 196);
    const auto *ii_197 = buffer.data(ii + 197);
    const auto *ii_198 = buffer.data(ii + 198);
    const auto *ii_199 = buffer.data(ii + 199);
    const auto *ii_200 = buffer.data(ii + 200);
    const auto *ii_201 = buffer.data(ii + 201);
    const auto *ii_202 = buffer.data(ii + 202);
    const auto *ii_203 = buffer.data(ii + 203);
    const auto *ii_204 = buffer.data(ii + 204);
    const auto *ii_205 = buffer.data(ii + 205);
    const auto *ii_206 = buffer.data(ii + 206);
    const auto *ii_207 = buffer.data(ii + 207);
    const auto *ii_208 = buffer.data(ii + 208);
    const auto *ii_209 = buffer.data(ii + 209);
    const auto *ii_210 = buffer.data(ii + 210);
    const auto *ii_211 = buffer.data(ii + 211);
    const auto *ii_212 = buffer.data(ii + 212);
    const auto *ii_213 = buffer.data(ii + 213);
    const auto *ii_214 = buffer.data(ii + 214);
    const auto *ii_215 = buffer.data(ii + 215);
    const auto *ii_216 = buffer.data(ii + 216);
    const auto *ii_217 = buffer.data(ii + 217);
    const auto *ii_218 = buffer.data(ii + 218);
    const auto *ii_219 = buffer.data(ii + 219);
    const auto *ii_220 = buffer.data(ii + 220);
    const auto *ii_221 = buffer.data(ii + 221);
    const auto *ii_222 = buffer.data(ii + 222);

    const auto *li_287 = buffer.data(li + 287);
    const auto *li_288 = buffer.data(li + 288);
    const auto *li_289 = buffer.data(li + 289);
    const auto *li_290 = buffer.data(li + 290);
    const auto *li_291 = buffer.data(li + 291);
    const auto *li_292 = buffer.data(li + 292);
    const auto *li_293 = buffer.data(li + 293);
    const auto *li_294 = buffer.data(li + 294);
    const auto *li_295 = buffer.data(li + 295);
    const auto *li_296 = buffer.data(li + 296);
    const auto *li_297 = buffer.data(li + 297);
    const auto *li_298 = buffer.data(li + 298);
    const auto *li_299 = buffer.data(li + 299);
    const auto *li_300 = buffer.data(li + 300);
    const auto *li_301 = buffer.data(li + 301);
    const auto *li_302 = buffer.data(li + 302);
    const auto *li_303 = buffer.data(li + 303);
    const auto *li_304 = buffer.data(li + 304);
    const auto *li_305 = buffer.data(li + 305);
    const auto *li_306 = buffer.data(li + 306);
    const auto *li_307 = buffer.data(li + 307);
    const auto *li_308 = buffer.data(li + 308);
    const auto *li_309 = buffer.data(li + 309);
    const auto *li_310 = buffer.data(li + 310);
    const auto *li_311 = buffer.data(li + 311);
    const auto *li_312 = buffer.data(li + 312);
    const auto *li_313 = buffer.data(li + 313);
    const auto *li_314 = buffer.data(li + 314);
    const auto *li_315 = buffer.data(li + 315);
    const auto *li_316 = buffer.data(li + 316);
    const auto *li_317 = buffer.data(li + 317);
    const auto *li_318 = buffer.data(li + 318);
    const auto *li_319 = buffer.data(li + 319);
    const auto *li_320 = buffer.data(li + 320);
    const auto *li_321 = buffer.data(li + 321);
    const auto *li_322 = buffer.data(li + 322);
    const auto *li_323 = buffer.data(li + 323);
    const auto *li_324 = buffer.data(li + 324);
    const auto *li_325 = buffer.data(li + 325);
    const auto *li_326 = buffer.data(li + 326);
    const auto *li_327 = buffer.data(li + 327);
    const auto *li_328 = buffer.data(li + 328);
    const auto *li_329 = buffer.data(li + 329);
    const auto *li_330 = buffer.data(li + 330);
    const auto *li_331 = buffer.data(li + 331);
    const auto *li_332 = buffer.data(li + 332);
    const auto *li_333 = buffer.data(li + 333);
    const auto *li_334 = buffer.data(li + 334);
    const auto *li_335 = buffer.data(li + 335);
    const auto *li_336 = buffer.data(li + 336);
    const auto *li_337 = buffer.data(li + 337);
    const auto *li_338 = buffer.data(li + 338);
    const auto *li_339 = buffer.data(li + 339);
    const auto *li_340 = buffer.data(li + 340);
    const auto *li_341 = buffer.data(li + 341);
    const auto *li_342 = buffer.data(li + 342);
    const auto *li_343 = buffer.data(li + 343);
    const auto *li_344 = buffer.data(li + 344);
    const auto *li_345 = buffer.data(li + 345);
    const auto *li_346 = buffer.data(li + 346);
    const auto *li_347 = buffer.data(li + 347);
    const auto *li_348 = buffer.data(li + 348);
    const auto *li_349 = buffer.data(li + 349);
    const auto *li_350 = buffer.data(li + 350);
    const auto *li_351 = buffer.data(li + 351);
    const auto *li_352 = buffer.data(li + 352);
    const auto *li_353 = buffer.data(li + 353);
    const auto *li_354 = buffer.data(li + 354);
    const auto *li_355 = buffer.data(li + 355);
    const auto *li_356 = buffer.data(li + 356);
    const auto *li_357 = buffer.data(li + 357);
    const auto *li_358 = buffer.data(li + 358);
    const auto *li_359 = buffer.data(li + 359);
    const auto *li_360 = buffer.data(li + 360);
    const auto *li_361 = buffer.data(li + 361);
    const auto *li_362 = buffer.data(li + 362);
    const auto *li_363 = buffer.data(li + 363);
    const auto *li_364 = buffer.data(li + 364);
    const auto *li_365 = buffer.data(li + 365);
    const auto *li_366 = buffer.data(li + 366);
    const auto *li_367 = buffer.data(li + 367);
    const auto *li_368 = buffer.data(li + 368);
    const auto *li_369 = buffer.data(li + 369);
    const auto *li_370 = buffer.data(li + 370);
    const auto *li_371 = buffer.data(li + 371);
    const auto *li_372 = buffer.data(li + 372);
    const auto *li_373 = buffer.data(li + 373);
    const auto *li_374 = buffer.data(li + 374);
    const auto *li_375 = buffer.data(li + 375);
    const auto *li_376 = buffer.data(li + 376);
    const auto *li_377 = buffer.data(li + 377);
    const auto *li_378 = buffer.data(li + 378);
    const auto *li_379 = buffer.data(li + 379);
    const auto *li_380 = buffer.data(li + 380);
    const auto *li_381 = buffer.data(li + 381);
    const auto *li_382 = buffer.data(li + 382);
    const auto *li_383 = buffer.data(li + 383);
    const auto *li_384 = buffer.data(li + 384);
    const auto *li_385 = buffer.data(li + 385);
    const auto *li_386 = buffer.data(li + 386);
    const auto *li_387 = buffer.data(li + 387);
    const auto *li_388 = buffer.data(li + 388);
    const auto *li_389 = buffer.data(li + 389);
    const auto *li_390 = buffer.data(li + 390);
    const auto *li_391 = buffer.data(li + 391);
    const auto *li_420 = buffer.data(li + 420);
    const auto *li_421 = buffer.data(li + 421);
    const auto *li_422 = buffer.data(li + 422);
    const auto *li_423 = buffer.data(li + 423);
    const auto *li_424 = buffer.data(li + 424);
    const auto *li_425 = buffer.data(li + 425);
    const auto *li_426 = buffer.data(li + 426);
    const auto *li_427 = buffer.data(li + 427);
    const auto *li_428 = buffer.data(li + 428);
    const auto *li_429 = buffer.data(li + 429);
    const auto *li_430 = buffer.data(li + 430);
    const auto *li_431 = buffer.data(li + 431);
    const auto *li_432 = buffer.data(li + 432);
    const auto *li_433 = buffer.data(li + 433);
    const auto *li_434 = buffer.data(li + 434);
    const auto *li_435 = buffer.data(li + 435);
    const auto *li_436 = buffer.data(li + 436);
    const auto *li_437 = buffer.data(li + 437);
    const auto *li_438 = buffer.data(li + 438);
    const auto *li_439 = buffer.data(li + 439);
    const auto *li_440 = buffer.data(li + 440);
    const auto *li_441 = buffer.data(li + 441);
    const auto *li_442 = buffer.data(li + 442);
    const auto *li_443 = buffer.data(li + 443);
    const auto *li_444 = buffer.data(li + 444);
    const auto *li_445 = buffer.data(li + 445);
    const auto *li_446 = buffer.data(li + 446);
    const auto *li_447 = buffer.data(li + 447);
    const auto *li_448 = buffer.data(li + 448);
    const auto *li_449 = buffer.data(li + 449);
    const auto *li_450 = buffer.data(li + 450);
    const auto *li_451 = buffer.data(li + 451);
    const auto *li_452 = buffer.data(li + 452);
    const auto *li_453 = buffer.data(li + 453);
    const auto *li_454 = buffer.data(li + 454);
    const auto *li_455 = buffer.data(li + 455);
    const auto *li_456 = buffer.data(li + 456);
    const auto *li_457 = buffer.data(li + 457);
    const auto *li_458 = buffer.data(li + 458);
    const auto *li_459 = buffer.data(li + 459);
    const auto *li_460 = buffer.data(li + 460);
    const auto *li_461 = buffer.data(li + 461);
    const auto *li_462 = buffer.data(li + 462);
    const auto *li_463 = buffer.data(li + 463);
    const auto *li_464 = buffer.data(li + 464);
    const auto *li_465 = buffer.data(li + 465);
    const auto *li_466 = buffer.data(li + 466);
    const auto *li_467 = buffer.data(li + 467);
    const auto *li_468 = buffer.data(li + 468);
    const auto *li_469 = buffer.data(li + 469);
    const auto *li_470 = buffer.data(li + 470);
    const auto *li_471 = buffer.data(li + 471);
    const auto *li_472 = buffer.data(li + 472);
    const auto *li_473 = buffer.data(li + 473);
    const auto *li_474 = buffer.data(li + 474);

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ii_91, ii_92, ii_93, ii_94, ii_95, \
                         li_287, li_288, li_289, li_290, li_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * ii_91[k]
                   + f_0 * li_287[k];

        t_176[k] = -3.0 * ii_92[k]
                   + f_0 * li_288[k];

        t_177[k] = -3.0 * ii_93[k]
                   + f_0 * li_289[k];

        t_178[k] = -3.0 * ii_94[k]
                   + f_0 * li_290[k];

        t_179[k] = -3.0 * ii_95[k]
                   + f_0 * li_291[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ii_96, ii_97, ii_98, ii_99, \
                         ii_100, li_292, li_293, li_294, li_295, \
                         li_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * ii_96[k]
                   + f_0 * li_292[k];

        t_181[k] = -3.0 * ii_97[k]
                   + f_0 * li_293[k];

        t_182[k] = -3.0 * ii_98[k]
                   + f_0 * li_294[k];

        t_183[k] = -3.0 * ii_99[k]
                   + f_0 * li_295[k];

        t_184[k] = -3.0 * ii_100[k]
                   + f_0 * li_296[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ii_101, ii_102, ii_103, ii_104, \
                         ii_105, li_297, li_298, li_299, li_300, \
                         li_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * ii_101[k]
                   + f_0 * li_297[k];

        t_186[k] = -3.0 * ii_102[k]
                   + f_0 * li_298[k];

        t_187[k] = -3.0 * ii_103[k]
                   + f_0 * li_299[k];

        t_188[k] = -3.0 * ii_104[k]
                   + f_0 * li_300[k];

        t_189[k] = -3.0 * ii_105[k]
                   + f_0 * li_301[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ii_106, ii_107, ii_108, ii_109, \
                         ii_110, li_302, li_303, li_304, li_305, \
                         li_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * ii_106[k]
                   + f_0 * li_302[k];

        t_191[k] = -3.0 * ii_107[k]
                   + f_0 * li_303[k];

        t_192[k] = -3.0 * ii_108[k]
                   + f_0 * li_304[k];

        t_193[k] = -3.0 * ii_109[k]
                   + f_0 * li_305[k];

        t_194[k] = -3.0 * ii_110[k]
                   + f_0 * li_306[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ii_111, ii_112, ii_113, ii_114, \
                         ii_115, li_307, li_308, li_309, li_310, \
                         li_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * ii_111[k]
                   + f_0 * li_307[k];

        t_196[k] = -2.0 * ii_112[k]
                   + f_0 * li_308[k];

        t_197[k] = -2.0 * ii_113[k]
                   + f_0 * li_309[k];

        t_198[k] = -2.0 * ii_114[k]
                   + f_0 * li_310[k];

        t_199[k] = -2.0 * ii_115[k]
                   + f_0 * li_311[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ii_116, ii_117, ii_118, ii_119, \
                         ii_120, li_312, li_313, li_314, li_315, \
                         li_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * ii_116[k]
                   + f_0 * li_312[k];

        t_201[k] = -2.0 * ii_117[k]
                   + f_0 * li_313[k];

        t_202[k] = -2.0 * ii_118[k]
                   + f_0 * li_314[k];

        t_203[k] = -2.0 * ii_119[k]
                   + f_0 * li_315[k];

        t_204[k] = -2.0 * ii_120[k]
                   + f_0 * li_316[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ii_121, ii_122, ii_123, ii_124, \
                         ii_125, li_317, li_318, li_319, li_320, \
                         li_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * ii_121[k]
                   + f_0 * li_317[k];

        t_206[k] = -2.0 * ii_122[k]
                   + f_0 * li_318[k];

        t_207[k] = -2.0 * ii_123[k]
                   + f_0 * li_319[k];

        t_208[k] = -2.0 * ii_124[k]
                   + f_0 * li_320[k];

        t_209[k] = -2.0 * ii_125[k]
                   + f_0 * li_321[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ii_126, ii_127, ii_128, ii_129, \
                         ii_130, li_322, li_323, li_324, li_325, \
                         li_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * ii_126[k]
                   + f_0 * li_322[k];

        t_211[k] = -2.0 * ii_127[k]
                   + f_0 * li_323[k];

        t_212[k] = -2.0 * ii_128[k]
                   + f_0 * li_324[k];

        t_213[k] = -2.0 * ii_129[k]
                   + f_0 * li_325[k];

        t_214[k] = -2.0 * ii_130[k]
                   + f_0 * li_326[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ii_131, ii_132, ii_133, ii_134, \
                         ii_135, li_327, li_328, li_329, li_330, \
                         li_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * ii_131[k]
                   + f_0 * li_327[k];

        t_216[k] = -2.0 * ii_132[k]
                   + f_0 * li_328[k];

        t_217[k] = -2.0 * ii_133[k]
                   + f_0 * li_329[k];

        t_218[k] = -2.0 * ii_134[k]
                   + f_0 * li_330[k];

        t_219[k] = -2.0 * ii_135[k]
                   + f_0 * li_331[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ii_136, ii_137, ii_138, ii_139, \
                         ii_140, li_332, li_333, li_334, li_335, \
                         li_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * ii_136[k]
                   + f_0 * li_332[k];

        t_221[k] = -2.0 * ii_137[k]
                   + f_0 * li_333[k];

        t_222[k] = -2.0 * ii_138[k]
                   + f_0 * li_334[k];

        t_223[k] = -2.0 * ii_139[k]
                   + f_0 * li_335[k];

        t_224[k] = -ii_140[k]
                   + f_0 * li_336[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ii_141, ii_142, ii_143, ii_144, \
                         ii_145, li_337, li_338, li_339, li_340, \
                         li_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -ii_141[k]
                   + f_0 * li_337[k];

        t_226[k] = -ii_142[k]
                   + f_0 * li_338[k];

        t_227[k] = -ii_143[k]
                   + f_0 * li_339[k];

        t_228[k] = -ii_144[k]
                   + f_0 * li_340[k];

        t_229[k] = -ii_145[k]
                   + f_0 * li_341[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ii_146, ii_147, ii_148, ii_149, \
                         ii_150, li_342, li_343, li_344, li_345, \
                         li_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -ii_146[k]
                   + f_0 * li_342[k];

        t_231[k] = -ii_147[k]
                   + f_0 * li_343[k];

        t_232[k] = -ii_148[k]
                   + f_0 * li_344[k];

        t_233[k] = -ii_149[k]
                   + f_0 * li_345[k];

        t_234[k] = -ii_150[k]
                   + f_0 * li_346[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ii_151, ii_152, ii_153, ii_154, \
                         ii_155, li_347, li_348, li_349, li_350, \
                         li_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -ii_151[k]
                   + f_0 * li_347[k];

        t_236[k] = -ii_152[k]
                   + f_0 * li_348[k];

        t_237[k] = -ii_153[k]
                   + f_0 * li_349[k];

        t_238[k] = -ii_154[k]
                   + f_0 * li_350[k];

        t_239[k] = -ii_155[k]
                   + f_0 * li_351[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ii_156, ii_157, ii_158, ii_159, \
                         ii_160, li_352, li_353, li_354, li_355, \
                         li_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -ii_156[k]
                   + f_0 * li_352[k];

        t_241[k] = -ii_157[k]
                   + f_0 * li_353[k];

        t_242[k] = -ii_158[k]
                   + f_0 * li_354[k];

        t_243[k] = -ii_159[k]
                   + f_0 * li_355[k];

        t_244[k] = -ii_160[k]
                   + f_0 * li_356[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ii_161, ii_162, ii_163, ii_164, \
                         ii_165, li_357, li_358, li_359, li_360, \
                         li_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -ii_161[k]
                   + f_0 * li_357[k];

        t_246[k] = -ii_162[k]
                   + f_0 * li_358[k];

        t_247[k] = -ii_163[k]
                   + f_0 * li_359[k];

        t_248[k] = -ii_164[k]
                   + f_0 * li_360[k];

        t_249[k] = -ii_165[k]
                   + f_0 * li_361[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, t_255, t_256, ii_166, ii_167, \
                         li_362, li_363, li_364, li_365, li_366, li_367, \
                         li_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -ii_166[k]
                   + f_0 * li_362[k];

        t_251[k] = -ii_167[k]
                   + f_0 * li_363[k];

        t_252[k] = f_0 * li_364[k];

        t_253[k] = f_0 * li_365[k];

        t_254[k] = f_0 * li_366[k];

        t_255[k] = f_0 * li_367[k];

        t_256[k] = f_0 * li_368[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, t_262, t_263, t_264, li_369, \
                         li_370, li_371, li_372, li_373, li_374, li_375, \
                         li_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_0 * li_369[k];

        t_258[k] = f_0 * li_370[k];

        t_259[k] = f_0 * li_371[k];

        t_260[k] = f_0 * li_372[k];

        t_261[k] = f_0 * li_373[k];

        t_262[k] = f_0 * li_374[k];

        t_263[k] = f_0 * li_375[k];

        t_264[k] = f_0 * li_376[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, t_271, t_272, li_377, \
                         li_378, li_379, li_380, li_381, li_382, li_383, \
                         li_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_0 * li_377[k];

        t_266[k] = f_0 * li_378[k];

        t_267[k] = f_0 * li_379[k];

        t_268[k] = f_0 * li_380[k];

        t_269[k] = f_0 * li_381[k];

        t_270[k] = f_0 * li_382[k];

        t_271[k] = f_0 * li_383[k];

        t_272[k] = f_0 * li_384[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, t_278, t_279, li_385, li_386, \
                         li_387, li_388, li_389, li_390, li_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_0 * li_385[k];

        t_274[k] = f_0 * li_386[k];

        t_275[k] = f_0 * li_387[k];

        t_276[k] = f_0 * li_388[k];

        t_277[k] = f_0 * li_389[k];

        t_278[k] = f_0 * li_390[k];

        t_279[k] = f_0 * li_391[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ii_168, ii_169, ii_170, ii_171, \
                         ii_172, li_420, li_421, li_422, li_423, \
                         li_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -4.0 * ii_168[k]
                   + f_0 * li_420[k];

        t_281[k] = -4.0 * ii_169[k]
                   + f_0 * li_421[k];

        t_282[k] = -4.0 * ii_170[k]
                   + f_0 * li_422[k];

        t_283[k] = -4.0 * ii_171[k]
                   + f_0 * li_423[k];

        t_284[k] = -4.0 * ii_172[k]
                   + f_0 * li_424[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ii_173, ii_174, ii_175, ii_176, \
                         ii_177, li_425, li_426, li_427, li_428, \
                         li_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -4.0 * ii_173[k]
                   + f_0 * li_425[k];

        t_286[k] = -4.0 * ii_174[k]
                   + f_0 * li_426[k];

        t_287[k] = -4.0 * ii_175[k]
                   + f_0 * li_427[k];

        t_288[k] = -4.0 * ii_176[k]
                   + f_0 * li_428[k];

        t_289[k] = -4.0 * ii_177[k]
                   + f_0 * li_429[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ii_178, ii_179, ii_180, ii_181, \
                         ii_182, li_430, li_431, li_432, li_433, \
                         li_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -4.0 * ii_178[k]
                   + f_0 * li_430[k];

        t_291[k] = -4.0 * ii_179[k]
                   + f_0 * li_431[k];

        t_292[k] = -4.0 * ii_180[k]
                   + f_0 * li_432[k];

        t_293[k] = -4.0 * ii_181[k]
                   + f_0 * li_433[k];

        t_294[k] = -4.0 * ii_182[k]
                   + f_0 * li_434[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ii_183, ii_184, ii_185, ii_186, \
                         ii_187, li_435, li_436, li_437, li_438, \
                         li_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * ii_183[k]
                   + f_0 * li_435[k];

        t_296[k] = -4.0 * ii_184[k]
                   + f_0 * li_436[k];

        t_297[k] = -4.0 * ii_185[k]
                   + f_0 * li_437[k];

        t_298[k] = -4.0 * ii_186[k]
                   + f_0 * li_438[k];

        t_299[k] = -4.0 * ii_187[k]
                   + f_0 * li_439[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ii_188, ii_189, ii_190, ii_191, \
                         ii_192, li_440, li_441, li_442, li_443, \
                         li_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * ii_188[k]
                   + f_0 * li_440[k];

        t_301[k] = -4.0 * ii_189[k]
                   + f_0 * li_441[k];

        t_302[k] = -4.0 * ii_190[k]
                   + f_0 * li_442[k];

        t_303[k] = -4.0 * ii_191[k]
                   + f_0 * li_443[k];

        t_304[k] = -4.0 * ii_192[k]
                   + f_0 * li_444[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ii_193, ii_194, ii_195, ii_196, \
                         ii_197, li_445, li_446, li_447, li_448, \
                         li_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * ii_193[k]
                   + f_0 * li_445[k];

        t_306[k] = -4.0 * ii_194[k]
                   + f_0 * li_446[k];

        t_307[k] = -4.0 * ii_195[k]
                   + f_0 * li_447[k];

        t_308[k] = -3.0 * ii_196[k]
                   + f_0 * li_448[k];

        t_309[k] = -3.0 * ii_197[k]
                   + f_0 * li_449[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ii_198, ii_199, ii_200, ii_201, \
                         ii_202, li_450, li_451, li_452, li_453, \
                         li_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * ii_198[k]
                   + f_0 * li_450[k];

        t_311[k] = -3.0 * ii_199[k]
                   + f_0 * li_451[k];

        t_312[k] = -3.0 * ii_200[k]
                   + f_0 * li_452[k];

        t_313[k] = -3.0 * ii_201[k]
                   + f_0 * li_453[k];

        t_314[k] = -3.0 * ii_202[k]
                   + f_0 * li_454[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ii_203, ii_204, ii_205, ii_206, \
                         ii_207, li_455, li_456, li_457, li_458, \
                         li_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -3.0 * ii_203[k]
                   + f_0 * li_455[k];

        t_316[k] = -3.0 * ii_204[k]
                   + f_0 * li_456[k];

        t_317[k] = -3.0 * ii_205[k]
                   + f_0 * li_457[k];

        t_318[k] = -3.0 * ii_206[k]
                   + f_0 * li_458[k];

        t_319[k] = -3.0 * ii_207[k]
                   + f_0 * li_459[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ii_208, ii_209, ii_210, ii_211, \
                         ii_212, li_460, li_461, li_462, li_463, \
                         li_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -3.0 * ii_208[k]
                   + f_0 * li_460[k];

        t_321[k] = -3.0 * ii_209[k]
                   + f_0 * li_461[k];

        t_322[k] = -3.0 * ii_210[k]
                   + f_0 * li_462[k];

        t_323[k] = -3.0 * ii_211[k]
                   + f_0 * li_463[k];

        t_324[k] = -3.0 * ii_212[k]
                   + f_0 * li_464[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ii_213, ii_214, ii_215, ii_216, \
                         ii_217, li_465, li_466, li_467, li_468, \
                         li_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -3.0 * ii_213[k]
                   + f_0 * li_465[k];

        t_326[k] = -3.0 * ii_214[k]
                   + f_0 * li_466[k];

        t_327[k] = -3.0 * ii_215[k]
                   + f_0 * li_467[k];

        t_328[k] = -3.0 * ii_216[k]
                   + f_0 * li_468[k];

        t_329[k] = -3.0 * ii_217[k]
                   + f_0 * li_469[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ii_218, ii_219, ii_220, ii_221, \
                         ii_222, li_470, li_471, li_472, li_473, \
                         li_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -3.0 * ii_218[k]
                   + f_0 * li_470[k];

        t_331[k] = -3.0 * ii_219[k]
                   + f_0 * li_471[k];

        t_332[k] = -3.0 * ii_220[k]
                   + f_0 * li_472[k];

        t_333[k] = -3.0 * ii_221[k]
                   + f_0 * li_473[k];

        t_334[k] = -3.0 * ii_222[k]
                   + f_0 * li_474[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_223 = buffer.data(ii + 223);
    const auto *ii_224 = buffer.data(ii + 224);
    const auto *ii_225 = buffer.data(ii + 225);
    const auto *ii_226 = buffer.data(ii + 226);
    const auto *ii_227 = buffer.data(ii + 227);
    const auto *ii_228 = buffer.data(ii + 228);
    const auto *ii_229 = buffer.data(ii + 229);
    const auto *ii_230 = buffer.data(ii + 230);
    const auto *ii_231 = buffer.data(ii + 231);
    const auto *ii_232 = buffer.data(ii + 232);
    const auto *ii_233 = buffer.data(ii + 233);
    const auto *ii_234 = buffer.data(ii + 234);
    const auto *ii_235 = buffer.data(ii + 235);
    const auto *ii_236 = buffer.data(ii + 236);
    const auto *ii_237 = buffer.data(ii + 237);
    const auto *ii_238 = buffer.data(ii + 238);
    const auto *ii_239 = buffer.data(ii + 239);
    const auto *ii_240 = buffer.data(ii + 240);
    const auto *ii_241 = buffer.data(ii + 241);
    const auto *ii_242 = buffer.data(ii + 242);
    const auto *ii_243 = buffer.data(ii + 243);
    const auto *ii_244 = buffer.data(ii + 244);
    const auto *ii_245 = buffer.data(ii + 245);
    const auto *ii_246 = buffer.data(ii + 246);
    const auto *ii_247 = buffer.data(ii + 247);
    const auto *ii_248 = buffer.data(ii + 248);
    const auto *ii_249 = buffer.data(ii + 249);
    const auto *ii_250 = buffer.data(ii + 250);
    const auto *ii_251 = buffer.data(ii + 251);
    const auto *ii_252 = buffer.data(ii + 252);
    const auto *ii_253 = buffer.data(ii + 253);
    const auto *ii_254 = buffer.data(ii + 254);
    const auto *ii_255 = buffer.data(ii + 255);
    const auto *ii_256 = buffer.data(ii + 256);
    const auto *ii_257 = buffer.data(ii + 257);
    const auto *ii_258 = buffer.data(ii + 258);
    const auto *ii_259 = buffer.data(ii + 259);
    const auto *ii_260 = buffer.data(ii + 260);
    const auto *ii_261 = buffer.data(ii + 261);
    const auto *ii_262 = buffer.data(ii + 262);
    const auto *ii_263 = buffer.data(ii + 263);
    const auto *ii_264 = buffer.data(ii + 264);
    const auto *ii_265 = buffer.data(ii + 265);
    const auto *ii_266 = buffer.data(ii + 266);
    const auto *ii_267 = buffer.data(ii + 267);
    const auto *ii_268 = buffer.data(ii + 268);
    const auto *ii_269 = buffer.data(ii + 269);
    const auto *ii_270 = buffer.data(ii + 270);
    const auto *ii_271 = buffer.data(ii + 271);
    const auto *ii_272 = buffer.data(ii + 272);
    const auto *ii_273 = buffer.data(ii + 273);
    const auto *ii_274 = buffer.data(ii + 274);
    const auto *ii_275 = buffer.data(ii + 275);
    const auto *ii_276 = buffer.data(ii + 276);
    const auto *ii_277 = buffer.data(ii + 277);
    const auto *ii_278 = buffer.data(ii + 278);
    const auto *ii_279 = buffer.data(ii + 279);
    const auto *ii_280 = buffer.data(ii + 280);
    const auto *ii_281 = buffer.data(ii + 281);
    const auto *ii_282 = buffer.data(ii + 282);
    const auto *ii_283 = buffer.data(ii + 283);
    const auto *ii_284 = buffer.data(ii + 284);
    const auto *ii_285 = buffer.data(ii + 285);
    const auto *ii_286 = buffer.data(ii + 286);
    const auto *ii_287 = buffer.data(ii + 287);
    const auto *ii_288 = buffer.data(ii + 288);
    const auto *ii_289 = buffer.data(ii + 289);
    const auto *ii_290 = buffer.data(ii + 290);
    const auto *ii_291 = buffer.data(ii + 291);
    const auto *ii_292 = buffer.data(ii + 292);
    const auto *ii_293 = buffer.data(ii + 293);
    const auto *ii_294 = buffer.data(ii + 294);
    const auto *ii_295 = buffer.data(ii + 295);
    const auto *ii_296 = buffer.data(ii + 296);
    const auto *ii_297 = buffer.data(ii + 297);
    const auto *ii_298 = buffer.data(ii + 298);
    const auto *ii_299 = buffer.data(ii + 299);
    const auto *ii_300 = buffer.data(ii + 300);
    const auto *ii_301 = buffer.data(ii + 301);
    const auto *ii_302 = buffer.data(ii + 302);
    const auto *ii_303 = buffer.data(ii + 303);
    const auto *ii_304 = buffer.data(ii + 304);
    const auto *ii_305 = buffer.data(ii + 305);
    const auto *ii_306 = buffer.data(ii + 306);
    const auto *ii_307 = buffer.data(ii + 307);
    const auto *ii_308 = buffer.data(ii + 308);
    const auto *ii_309 = buffer.data(ii + 309);
    const auto *ii_310 = buffer.data(ii + 310);
    const auto *ii_311 = buffer.data(ii + 311);
    const auto *ii_312 = buffer.data(ii + 312);
    const auto *ii_313 = buffer.data(ii + 313);
    const auto *ii_314 = buffer.data(ii + 314);
    const auto *ii_315 = buffer.data(ii + 315);
    const auto *ii_316 = buffer.data(ii + 316);
    const auto *ii_317 = buffer.data(ii + 317);
    const auto *ii_318 = buffer.data(ii + 318);
    const auto *ii_319 = buffer.data(ii + 319);
    const auto *ii_320 = buffer.data(ii + 320);
    const auto *ii_321 = buffer.data(ii + 321);
    const auto *ii_322 = buffer.data(ii + 322);
    const auto *ii_323 = buffer.data(ii + 323);
    const auto *ii_324 = buffer.data(ii + 324);
    const auto *ii_325 = buffer.data(ii + 325);
    const auto *ii_326 = buffer.data(ii + 326);
    const auto *ii_327 = buffer.data(ii + 327);
    const auto *ii_328 = buffer.data(ii + 328);
    const auto *ii_329 = buffer.data(ii + 329);
    const auto *ii_330 = buffer.data(ii + 330);
    const auto *ii_331 = buffer.data(ii + 331);
    const auto *ii_332 = buffer.data(ii + 332);
    const auto *ii_333 = buffer.data(ii + 333);
    const auto *ii_334 = buffer.data(ii + 334);
    const auto *ii_335 = buffer.data(ii + 335);
    const auto *ii_336 = buffer.data(ii + 336);
    const auto *ii_337 = buffer.data(ii + 337);
    const auto *ii_338 = buffer.data(ii + 338);
    const auto *ii_339 = buffer.data(ii + 339);
    const auto *ii_340 = buffer.data(ii + 340);
    const auto *ii_341 = buffer.data(ii + 341);
    const auto *ii_342 = buffer.data(ii + 342);
    const auto *ii_343 = buffer.data(ii + 343);
    const auto *ii_344 = buffer.data(ii + 344);
    const auto *ii_345 = buffer.data(ii + 345);
    const auto *ii_346 = buffer.data(ii + 346);
    const auto *ii_347 = buffer.data(ii + 347);
    const auto *ii_348 = buffer.data(ii + 348);
    const auto *ii_349 = buffer.data(ii + 349);
    const auto *ii_350 = buffer.data(ii + 350);
    const auto *ii_351 = buffer.data(ii + 351);
    const auto *ii_352 = buffer.data(ii + 352);
    const auto *ii_353 = buffer.data(ii + 353);
    const auto *ii_354 = buffer.data(ii + 354);

    const auto *li_475 = buffer.data(li + 475);
    const auto *li_476 = buffer.data(li + 476);
    const auto *li_477 = buffer.data(li + 477);
    const auto *li_478 = buffer.data(li + 478);
    const auto *li_479 = buffer.data(li + 479);
    const auto *li_480 = buffer.data(li + 480);
    const auto *li_481 = buffer.data(li + 481);
    const auto *li_482 = buffer.data(li + 482);
    const auto *li_483 = buffer.data(li + 483);
    const auto *li_484 = buffer.data(li + 484);
    const auto *li_485 = buffer.data(li + 485);
    const auto *li_486 = buffer.data(li + 486);
    const auto *li_487 = buffer.data(li + 487);
    const auto *li_488 = buffer.data(li + 488);
    const auto *li_489 = buffer.data(li + 489);
    const auto *li_490 = buffer.data(li + 490);
    const auto *li_491 = buffer.data(li + 491);
    const auto *li_492 = buffer.data(li + 492);
    const auto *li_493 = buffer.data(li + 493);
    const auto *li_494 = buffer.data(li + 494);
    const auto *li_495 = buffer.data(li + 495);
    const auto *li_496 = buffer.data(li + 496);
    const auto *li_497 = buffer.data(li + 497);
    const auto *li_498 = buffer.data(li + 498);
    const auto *li_499 = buffer.data(li + 499);
    const auto *li_500 = buffer.data(li + 500);
    const auto *li_501 = buffer.data(li + 501);
    const auto *li_502 = buffer.data(li + 502);
    const auto *li_503 = buffer.data(li + 503);
    const auto *li_504 = buffer.data(li + 504);
    const auto *li_505 = buffer.data(li + 505);
    const auto *li_506 = buffer.data(li + 506);
    const auto *li_507 = buffer.data(li + 507);
    const auto *li_508 = buffer.data(li + 508);
    const auto *li_509 = buffer.data(li + 509);
    const auto *li_510 = buffer.data(li + 510);
    const auto *li_511 = buffer.data(li + 511);
    const auto *li_512 = buffer.data(li + 512);
    const auto *li_513 = buffer.data(li + 513);
    const auto *li_514 = buffer.data(li + 514);
    const auto *li_515 = buffer.data(li + 515);
    const auto *li_516 = buffer.data(li + 516);
    const auto *li_517 = buffer.data(li + 517);
    const auto *li_518 = buffer.data(li + 518);
    const auto *li_519 = buffer.data(li + 519);
    const auto *li_520 = buffer.data(li + 520);
    const auto *li_521 = buffer.data(li + 521);
    const auto *li_522 = buffer.data(li + 522);
    const auto *li_523 = buffer.data(li + 523);
    const auto *li_524 = buffer.data(li + 524);
    const auto *li_525 = buffer.data(li + 525);
    const auto *li_526 = buffer.data(li + 526);
    const auto *li_527 = buffer.data(li + 527);
    const auto *li_528 = buffer.data(li + 528);
    const auto *li_529 = buffer.data(li + 529);
    const auto *li_530 = buffer.data(li + 530);
    const auto *li_531 = buffer.data(li + 531);
    const auto *li_532 = buffer.data(li + 532);
    const auto *li_533 = buffer.data(li + 533);
    const auto *li_534 = buffer.data(li + 534);
    const auto *li_535 = buffer.data(li + 535);
    const auto *li_536 = buffer.data(li + 536);
    const auto *li_537 = buffer.data(li + 537);
    const auto *li_538 = buffer.data(li + 538);
    const auto *li_539 = buffer.data(li + 539);
    const auto *li_540 = buffer.data(li + 540);
    const auto *li_541 = buffer.data(li + 541);
    const auto *li_542 = buffer.data(li + 542);
    const auto *li_543 = buffer.data(li + 543);
    const auto *li_544 = buffer.data(li + 544);
    const auto *li_545 = buffer.data(li + 545);
    const auto *li_546 = buffer.data(li + 546);
    const auto *li_547 = buffer.data(li + 547);
    const auto *li_548 = buffer.data(li + 548);
    const auto *li_549 = buffer.data(li + 549);
    const auto *li_550 = buffer.data(li + 550);
    const auto *li_551 = buffer.data(li + 551);
    const auto *li_552 = buffer.data(li + 552);
    const auto *li_553 = buffer.data(li + 553);
    const auto *li_554 = buffer.data(li + 554);
    const auto *li_555 = buffer.data(li + 555);
    const auto *li_556 = buffer.data(li + 556);
    const auto *li_557 = buffer.data(li + 557);
    const auto *li_558 = buffer.data(li + 558);
    const auto *li_559 = buffer.data(li + 559);
    const auto *li_588 = buffer.data(li + 588);
    const auto *li_589 = buffer.data(li + 589);
    const auto *li_590 = buffer.data(li + 590);
    const auto *li_591 = buffer.data(li + 591);
    const auto *li_592 = buffer.data(li + 592);
    const auto *li_593 = buffer.data(li + 593);
    const auto *li_594 = buffer.data(li + 594);
    const auto *li_595 = buffer.data(li + 595);
    const auto *li_596 = buffer.data(li + 596);
    const auto *li_597 = buffer.data(li + 597);
    const auto *li_598 = buffer.data(li + 598);
    const auto *li_599 = buffer.data(li + 599);
    const auto *li_600 = buffer.data(li + 600);
    const auto *li_601 = buffer.data(li + 601);
    const auto *li_602 = buffer.data(li + 602);
    const auto *li_603 = buffer.data(li + 603);
    const auto *li_604 = buffer.data(li + 604);
    const auto *li_605 = buffer.data(li + 605);
    const auto *li_606 = buffer.data(li + 606);
    const auto *li_607 = buffer.data(li + 607);
    const auto *li_608 = buffer.data(li + 608);
    const auto *li_609 = buffer.data(li + 609);
    const auto *li_610 = buffer.data(li + 610);
    const auto *li_611 = buffer.data(li + 611);
    const auto *li_612 = buffer.data(li + 612);
    const auto *li_613 = buffer.data(li + 613);
    const auto *li_614 = buffer.data(li + 614);
    const auto *li_615 = buffer.data(li + 615);
    const auto *li_616 = buffer.data(li + 616);
    const auto *li_617 = buffer.data(li + 617);
    const auto *li_618 = buffer.data(li + 618);
    const auto *li_619 = buffer.data(li + 619);
    const auto *li_620 = buffer.data(li + 620);
    const auto *li_621 = buffer.data(li + 621);
    const auto *li_622 = buffer.data(li + 622);
    const auto *li_623 = buffer.data(li + 623);
    const auto *li_624 = buffer.data(li + 624);
    const auto *li_625 = buffer.data(li + 625);
    const auto *li_626 = buffer.data(li + 626);
    const auto *li_627 = buffer.data(li + 627);
    const auto *li_628 = buffer.data(li + 628);
    const auto *li_629 = buffer.data(li + 629);
    const auto *li_630 = buffer.data(li + 630);
    const auto *li_631 = buffer.data(li + 631);
    const auto *li_632 = buffer.data(li + 632);
    const auto *li_633 = buffer.data(li + 633);
    const auto *li_634 = buffer.data(li + 634);
    const auto *li_635 = buffer.data(li + 635);
    const auto *li_636 = buffer.data(li + 636);
    const auto *li_637 = buffer.data(li + 637);
    const auto *li_638 = buffer.data(li + 638);
    const auto *li_639 = buffer.data(li + 639);
    const auto *li_640 = buffer.data(li + 640);
    const auto *li_641 = buffer.data(li + 641);
    const auto *li_642 = buffer.data(li + 642);
    const auto *li_643 = buffer.data(li + 643);
    const auto *li_644 = buffer.data(li + 644);
    const auto *li_645 = buffer.data(li + 645);
    const auto *li_646 = buffer.data(li + 646);
    const auto *li_647 = buffer.data(li + 647);
    const auto *li_648 = buffer.data(li + 648);
    const auto *li_649 = buffer.data(li + 649);
    const auto *li_650 = buffer.data(li + 650);
    const auto *li_651 = buffer.data(li + 651);
    const auto *li_652 = buffer.data(li + 652);
    const auto *li_653 = buffer.data(li + 653);
    const auto *li_654 = buffer.data(li + 654);
    const auto *li_655 = buffer.data(li + 655);
    const auto *li_656 = buffer.data(li + 656);
    const auto *li_657 = buffer.data(li + 657);
    const auto *li_658 = buffer.data(li + 658);
    const auto *li_659 = buffer.data(li + 659);
    const auto *li_660 = buffer.data(li + 660);
    const auto *li_661 = buffer.data(li + 661);
    const auto *li_662 = buffer.data(li + 662);

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ii_223, ii_224, ii_225, ii_226, \
                         ii_227, li_475, li_476, li_477, li_478, \
                         li_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -3.0 * ii_223[k]
                   + f_0 * li_475[k];

        t_336[k] = -2.0 * ii_224[k]
                   + f_0 * li_476[k];

        t_337[k] = -2.0 * ii_225[k]
                   + f_0 * li_477[k];

        t_338[k] = -2.0 * ii_226[k]
                   + f_0 * li_478[k];

        t_339[k] = -2.0 * ii_227[k]
                   + f_0 * li_479[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ii_228, ii_229, ii_230, ii_231, \
                         ii_232, li_480, li_481, li_482, li_483, \
                         li_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * ii_228[k]
                   + f_0 * li_480[k];

        t_341[k] = -2.0 * ii_229[k]
                   + f_0 * li_481[k];

        t_342[k] = -2.0 * ii_230[k]
                   + f_0 * li_482[k];

        t_343[k] = -2.0 * ii_231[k]
                   + f_0 * li_483[k];

        t_344[k] = -2.0 * ii_232[k]
                   + f_0 * li_484[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ii_233, ii_234, ii_235, ii_236, \
                         ii_237, li_485, li_486, li_487, li_488, \
                         li_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * ii_233[k]
                   + f_0 * li_485[k];

        t_346[k] = -2.0 * ii_234[k]
                   + f_0 * li_486[k];

        t_347[k] = -2.0 * ii_235[k]
                   + f_0 * li_487[k];

        t_348[k] = -2.0 * ii_236[k]
                   + f_0 * li_488[k];

        t_349[k] = -2.0 * ii_237[k]
                   + f_0 * li_489[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ii_238, ii_239, ii_240, ii_241, \
                         ii_242, li_490, li_491, li_492, li_493, \
                         li_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * ii_238[k]
                   + f_0 * li_490[k];

        t_351[k] = -2.0 * ii_239[k]
                   + f_0 * li_491[k];

        t_352[k] = -2.0 * ii_240[k]
                   + f_0 * li_492[k];

        t_353[k] = -2.0 * ii_241[k]
                   + f_0 * li_493[k];

        t_354[k] = -2.0 * ii_242[k]
                   + f_0 * li_494[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ii_243, ii_244, ii_245, ii_246, \
                         ii_247, li_495, li_496, li_497, li_498, \
                         li_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * ii_243[k]
                   + f_0 * li_495[k];

        t_356[k] = -2.0 * ii_244[k]
                   + f_0 * li_496[k];

        t_357[k] = -2.0 * ii_245[k]
                   + f_0 * li_497[k];

        t_358[k] = -2.0 * ii_246[k]
                   + f_0 * li_498[k];

        t_359[k] = -2.0 * ii_247[k]
                   + f_0 * li_499[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ii_248, ii_249, ii_250, ii_251, \
                         ii_252, li_500, li_501, li_502, li_503, \
                         li_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * ii_248[k]
                   + f_0 * li_500[k];

        t_361[k] = -2.0 * ii_249[k]
                   + f_0 * li_501[k];

        t_362[k] = -2.0 * ii_250[k]
                   + f_0 * li_502[k];

        t_363[k] = -2.0 * ii_251[k]
                   + f_0 * li_503[k];

        t_364[k] = -ii_252[k]
                   + f_0 * li_504[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ii_253, ii_254, ii_255, ii_256, \
                         ii_257, li_505, li_506, li_507, li_508, \
                         li_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -ii_253[k]
                   + f_0 * li_505[k];

        t_366[k] = -ii_254[k]
                   + f_0 * li_506[k];

        t_367[k] = -ii_255[k]
                   + f_0 * li_507[k];

        t_368[k] = -ii_256[k]
                   + f_0 * li_508[k];

        t_369[k] = -ii_257[k]
                   + f_0 * li_509[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ii_258, ii_259, ii_260, ii_261, \
                         ii_262, li_510, li_511, li_512, li_513, \
                         li_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -ii_258[k]
                   + f_0 * li_510[k];

        t_371[k] = -ii_259[k]
                   + f_0 * li_511[k];

        t_372[k] = -ii_260[k]
                   + f_0 * li_512[k];

        t_373[k] = -ii_261[k]
                   + f_0 * li_513[k];

        t_374[k] = -ii_262[k]
                   + f_0 * li_514[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ii_263, ii_264, ii_265, ii_266, \
                         ii_267, li_515, li_516, li_517, li_518, \
                         li_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -ii_263[k]
                   + f_0 * li_515[k];

        t_376[k] = -ii_264[k]
                   + f_0 * li_516[k];

        t_377[k] = -ii_265[k]
                   + f_0 * li_517[k];

        t_378[k] = -ii_266[k]
                   + f_0 * li_518[k];

        t_379[k] = -ii_267[k]
                   + f_0 * li_519[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ii_268, ii_269, ii_270, ii_271, \
                         ii_272, li_520, li_521, li_522, li_523, \
                         li_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -ii_268[k]
                   + f_0 * li_520[k];

        t_381[k] = -ii_269[k]
                   + f_0 * li_521[k];

        t_382[k] = -ii_270[k]
                   + f_0 * li_522[k];

        t_383[k] = -ii_271[k]
                   + f_0 * li_523[k];

        t_384[k] = -ii_272[k]
                   + f_0 * li_524[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ii_273, ii_274, ii_275, ii_276, \
                         ii_277, li_525, li_526, li_527, li_528, \
                         li_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -ii_273[k]
                   + f_0 * li_525[k];

        t_386[k] = -ii_274[k]
                   + f_0 * li_526[k];

        t_387[k] = -ii_275[k]
                   + f_0 * li_527[k];

        t_388[k] = -ii_276[k]
                   + f_0 * li_528[k];

        t_389[k] = -ii_277[k]
                   + f_0 * li_529[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, t_395, t_396, ii_278, ii_279, \
                         li_530, li_531, li_532, li_533, li_534, li_535, \
                         li_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -ii_278[k]
                   + f_0 * li_530[k];

        t_391[k] = -ii_279[k]
                   + f_0 * li_531[k];

        t_392[k] = f_0 * li_532[k];

        t_393[k] = f_0 * li_533[k];

        t_394[k] = f_0 * li_534[k];

        t_395[k] = f_0 * li_535[k];

        t_396[k] = f_0 * li_536[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, t_402, t_403, t_404, li_537, \
                         li_538, li_539, li_540, li_541, li_542, li_543, \
                         li_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_0 * li_537[k];

        t_398[k] = f_0 * li_538[k];

        t_399[k] = f_0 * li_539[k];

        t_400[k] = f_0 * li_540[k];

        t_401[k] = f_0 * li_541[k];

        t_402[k] = f_0 * li_542[k];

        t_403[k] = f_0 * li_543[k];

        t_404[k] = f_0 * li_544[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, li_545, \
                         li_546, li_547, li_548, li_549, li_550, li_551, \
                         li_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * li_545[k];

        t_406[k] = f_0 * li_546[k];

        t_407[k] = f_0 * li_547[k];

        t_408[k] = f_0 * li_548[k];

        t_409[k] = f_0 * li_549[k];

        t_410[k] = f_0 * li_550[k];

        t_411[k] = f_0 * li_551[k];

        t_412[k] = f_0 * li_552[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, li_553, li_554, \
                         li_555, li_556, li_557, li_558, li_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * li_553[k];

        t_414[k] = f_0 * li_554[k];

        t_415[k] = f_0 * li_555[k];

        t_416[k] = f_0 * li_556[k];

        t_417[k] = f_0 * li_557[k];

        t_418[k] = f_0 * li_558[k];

        t_419[k] = f_0 * li_559[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ii_280, ii_281, ii_282, ii_283, \
                         ii_284, li_588, li_589, li_590, li_591, \
                         li_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -5.0 * ii_280[k]
                   + f_0 * li_588[k];

        t_421[k] = -5.0 * ii_281[k]
                   + f_0 * li_589[k];

        t_422[k] = -5.0 * ii_282[k]
                   + f_0 * li_590[k];

        t_423[k] = -5.0 * ii_283[k]
                   + f_0 * li_591[k];

        t_424[k] = -5.0 * ii_284[k]
                   + f_0 * li_592[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ii_285, ii_286, ii_287, ii_288, \
                         ii_289, li_593, li_594, li_595, li_596, \
                         li_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -5.0 * ii_285[k]
                   + f_0 * li_593[k];

        t_426[k] = -5.0 * ii_286[k]
                   + f_0 * li_594[k];

        t_427[k] = -5.0 * ii_287[k]
                   + f_0 * li_595[k];

        t_428[k] = -5.0 * ii_288[k]
                   + f_0 * li_596[k];

        t_429[k] = -5.0 * ii_289[k]
                   + f_0 * li_597[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ii_290, ii_291, ii_292, ii_293, \
                         ii_294, li_598, li_599, li_600, li_601, \
                         li_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -5.0 * ii_290[k]
                   + f_0 * li_598[k];

        t_431[k] = -5.0 * ii_291[k]
                   + f_0 * li_599[k];

        t_432[k] = -5.0 * ii_292[k]
                   + f_0 * li_600[k];

        t_433[k] = -5.0 * ii_293[k]
                   + f_0 * li_601[k];

        t_434[k] = -5.0 * ii_294[k]
                   + f_0 * li_602[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ii_295, ii_296, ii_297, ii_298, \
                         ii_299, li_603, li_604, li_605, li_606, \
                         li_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -5.0 * ii_295[k]
                   + f_0 * li_603[k];

        t_436[k] = -5.0 * ii_296[k]
                   + f_0 * li_604[k];

        t_437[k] = -5.0 * ii_297[k]
                   + f_0 * li_605[k];

        t_438[k] = -5.0 * ii_298[k]
                   + f_0 * li_606[k];

        t_439[k] = -5.0 * ii_299[k]
                   + f_0 * li_607[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ii_300, ii_301, ii_302, ii_303, \
                         ii_304, li_608, li_609, li_610, li_611, \
                         li_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -5.0 * ii_300[k]
                   + f_0 * li_608[k];

        t_441[k] = -5.0 * ii_301[k]
                   + f_0 * li_609[k];

        t_442[k] = -5.0 * ii_302[k]
                   + f_0 * li_610[k];

        t_443[k] = -5.0 * ii_303[k]
                   + f_0 * li_611[k];

        t_444[k] = -5.0 * ii_304[k]
                   + f_0 * li_612[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ii_305, ii_306, ii_307, ii_308, \
                         ii_309, li_613, li_614, li_615, li_616, \
                         li_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -5.0 * ii_305[k]
                   + f_0 * li_613[k];

        t_446[k] = -5.0 * ii_306[k]
                   + f_0 * li_614[k];

        t_447[k] = -5.0 * ii_307[k]
                   + f_0 * li_615[k];

        t_448[k] = -4.0 * ii_308[k]
                   + f_0 * li_616[k];

        t_449[k] = -4.0 * ii_309[k]
                   + f_0 * li_617[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ii_310, ii_311, ii_312, ii_313, \
                         ii_314, li_618, li_619, li_620, li_621, \
                         li_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -4.0 * ii_310[k]
                   + f_0 * li_618[k];

        t_451[k] = -4.0 * ii_311[k]
                   + f_0 * li_619[k];

        t_452[k] = -4.0 * ii_312[k]
                   + f_0 * li_620[k];

        t_453[k] = -4.0 * ii_313[k]
                   + f_0 * li_621[k];

        t_454[k] = -4.0 * ii_314[k]
                   + f_0 * li_622[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ii_315, ii_316, ii_317, ii_318, \
                         ii_319, li_623, li_624, li_625, li_626, \
                         li_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -4.0 * ii_315[k]
                   + f_0 * li_623[k];

        t_456[k] = -4.0 * ii_316[k]
                   + f_0 * li_624[k];

        t_457[k] = -4.0 * ii_317[k]
                   + f_0 * li_625[k];

        t_458[k] = -4.0 * ii_318[k]
                   + f_0 * li_626[k];

        t_459[k] = -4.0 * ii_319[k]
                   + f_0 * li_627[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ii_320, ii_321, ii_322, ii_323, \
                         ii_324, li_628, li_629, li_630, li_631, \
                         li_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -4.0 * ii_320[k]
                   + f_0 * li_628[k];

        t_461[k] = -4.0 * ii_321[k]
                   + f_0 * li_629[k];

        t_462[k] = -4.0 * ii_322[k]
                   + f_0 * li_630[k];

        t_463[k] = -4.0 * ii_323[k]
                   + f_0 * li_631[k];

        t_464[k] = -4.0 * ii_324[k]
                   + f_0 * li_632[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ii_325, ii_326, ii_327, ii_328, \
                         ii_329, li_633, li_634, li_635, li_636, \
                         li_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -4.0 * ii_325[k]
                   + f_0 * li_633[k];

        t_466[k] = -4.0 * ii_326[k]
                   + f_0 * li_634[k];

        t_467[k] = -4.0 * ii_327[k]
                   + f_0 * li_635[k];

        t_468[k] = -4.0 * ii_328[k]
                   + f_0 * li_636[k];

        t_469[k] = -4.0 * ii_329[k]
                   + f_0 * li_637[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ii_330, ii_331, ii_332, ii_333, \
                         ii_334, li_638, li_639, li_640, li_641, \
                         li_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -4.0 * ii_330[k]
                   + f_0 * li_638[k];

        t_471[k] = -4.0 * ii_331[k]
                   + f_0 * li_639[k];

        t_472[k] = -4.0 * ii_332[k]
                   + f_0 * li_640[k];

        t_473[k] = -4.0 * ii_333[k]
                   + f_0 * li_641[k];

        t_474[k] = -4.0 * ii_334[k]
                   + f_0 * li_642[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ii_335, ii_336, ii_337, ii_338, \
                         ii_339, li_643, li_644, li_645, li_646, \
                         li_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -4.0 * ii_335[k]
                   + f_0 * li_643[k];

        t_476[k] = -3.0 * ii_336[k]
                   + f_0 * li_644[k];

        t_477[k] = -3.0 * ii_337[k]
                   + f_0 * li_645[k];

        t_478[k] = -3.0 * ii_338[k]
                   + f_0 * li_646[k];

        t_479[k] = -3.0 * ii_339[k]
                   + f_0 * li_647[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ii_340, ii_341, ii_342, ii_343, \
                         ii_344, li_648, li_649, li_650, li_651, \
                         li_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -3.0 * ii_340[k]
                   + f_0 * li_648[k];

        t_481[k] = -3.0 * ii_341[k]
                   + f_0 * li_649[k];

        t_482[k] = -3.0 * ii_342[k]
                   + f_0 * li_650[k];

        t_483[k] = -3.0 * ii_343[k]
                   + f_0 * li_651[k];

        t_484[k] = -3.0 * ii_344[k]
                   + f_0 * li_652[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ii_345, ii_346, ii_347, ii_348, \
                         ii_349, li_653, li_654, li_655, li_656, \
                         li_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -3.0 * ii_345[k]
                   + f_0 * li_653[k];

        t_486[k] = -3.0 * ii_346[k]
                   + f_0 * li_654[k];

        t_487[k] = -3.0 * ii_347[k]
                   + f_0 * li_655[k];

        t_488[k] = -3.0 * ii_348[k]
                   + f_0 * li_656[k];

        t_489[k] = -3.0 * ii_349[k]
                   + f_0 * li_657[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ii_350, ii_351, ii_352, ii_353, \
                         ii_354, li_658, li_659, li_660, li_661, \
                         li_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -3.0 * ii_350[k]
                   + f_0 * li_658[k];

        t_491[k] = -3.0 * ii_351[k]
                   + f_0 * li_659[k];

        t_492[k] = -3.0 * ii_352[k]
                   + f_0 * li_660[k];

        t_493[k] = -3.0 * ii_353[k]
                   + f_0 * li_661[k];

        t_494[k] = -3.0 * ii_354[k]
                   + f_0 * li_662[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_355 = buffer.data(ii + 355);
    const auto *ii_356 = buffer.data(ii + 356);
    const auto *ii_357 = buffer.data(ii + 357);
    const auto *ii_358 = buffer.data(ii + 358);
    const auto *ii_359 = buffer.data(ii + 359);
    const auto *ii_360 = buffer.data(ii + 360);
    const auto *ii_361 = buffer.data(ii + 361);
    const auto *ii_362 = buffer.data(ii + 362);
    const auto *ii_363 = buffer.data(ii + 363);
    const auto *ii_364 = buffer.data(ii + 364);
    const auto *ii_365 = buffer.data(ii + 365);
    const auto *ii_366 = buffer.data(ii + 366);
    const auto *ii_367 = buffer.data(ii + 367);
    const auto *ii_368 = buffer.data(ii + 368);
    const auto *ii_369 = buffer.data(ii + 369);
    const auto *ii_370 = buffer.data(ii + 370);
    const auto *ii_371 = buffer.data(ii + 371);
    const auto *ii_372 = buffer.data(ii + 372);
    const auto *ii_373 = buffer.data(ii + 373);
    const auto *ii_374 = buffer.data(ii + 374);
    const auto *ii_375 = buffer.data(ii + 375);
    const auto *ii_376 = buffer.data(ii + 376);
    const auto *ii_377 = buffer.data(ii + 377);
    const auto *ii_378 = buffer.data(ii + 378);
    const auto *ii_379 = buffer.data(ii + 379);
    const auto *ii_380 = buffer.data(ii + 380);
    const auto *ii_381 = buffer.data(ii + 381);
    const auto *ii_382 = buffer.data(ii + 382);
    const auto *ii_383 = buffer.data(ii + 383);
    const auto *ii_384 = buffer.data(ii + 384);
    const auto *ii_385 = buffer.data(ii + 385);
    const auto *ii_386 = buffer.data(ii + 386);
    const auto *ii_387 = buffer.data(ii + 387);
    const auto *ii_388 = buffer.data(ii + 388);
    const auto *ii_389 = buffer.data(ii + 389);
    const auto *ii_390 = buffer.data(ii + 390);
    const auto *ii_391 = buffer.data(ii + 391);
    const auto *ii_392 = buffer.data(ii + 392);
    const auto *ii_393 = buffer.data(ii + 393);
    const auto *ii_394 = buffer.data(ii + 394);
    const auto *ii_395 = buffer.data(ii + 395);
    const auto *ii_396 = buffer.data(ii + 396);
    const auto *ii_397 = buffer.data(ii + 397);
    const auto *ii_398 = buffer.data(ii + 398);
    const auto *ii_399 = buffer.data(ii + 399);
    const auto *ii_400 = buffer.data(ii + 400);
    const auto *ii_401 = buffer.data(ii + 401);
    const auto *ii_402 = buffer.data(ii + 402);
    const auto *ii_403 = buffer.data(ii + 403);
    const auto *ii_404 = buffer.data(ii + 404);
    const auto *ii_405 = buffer.data(ii + 405);
    const auto *ii_406 = buffer.data(ii + 406);
    const auto *ii_407 = buffer.data(ii + 407);
    const auto *ii_408 = buffer.data(ii + 408);
    const auto *ii_409 = buffer.data(ii + 409);
    const auto *ii_410 = buffer.data(ii + 410);
    const auto *ii_411 = buffer.data(ii + 411);
    const auto *ii_412 = buffer.data(ii + 412);
    const auto *ii_413 = buffer.data(ii + 413);
    const auto *ii_414 = buffer.data(ii + 414);
    const auto *ii_415 = buffer.data(ii + 415);
    const auto *ii_416 = buffer.data(ii + 416);
    const auto *ii_417 = buffer.data(ii + 417);
    const auto *ii_418 = buffer.data(ii + 418);
    const auto *ii_419 = buffer.data(ii + 419);
    const auto *ii_420 = buffer.data(ii + 420);
    const auto *ii_421 = buffer.data(ii + 421);
    const auto *ii_422 = buffer.data(ii + 422);
    const auto *ii_423 = buffer.data(ii + 423);
    const auto *ii_424 = buffer.data(ii + 424);
    const auto *ii_425 = buffer.data(ii + 425);
    const auto *ii_426 = buffer.data(ii + 426);
    const auto *ii_427 = buffer.data(ii + 427);
    const auto *ii_428 = buffer.data(ii + 428);
    const auto *ii_429 = buffer.data(ii + 429);
    const auto *ii_430 = buffer.data(ii + 430);
    const auto *ii_431 = buffer.data(ii + 431);
    const auto *ii_432 = buffer.data(ii + 432);
    const auto *ii_433 = buffer.data(ii + 433);
    const auto *ii_434 = buffer.data(ii + 434);
    const auto *ii_435 = buffer.data(ii + 435);
    const auto *ii_436 = buffer.data(ii + 436);
    const auto *ii_437 = buffer.data(ii + 437);
    const auto *ii_438 = buffer.data(ii + 438);
    const auto *ii_439 = buffer.data(ii + 439);
    const auto *ii_440 = buffer.data(ii + 440);
    const auto *ii_441 = buffer.data(ii + 441);
    const auto *ii_442 = buffer.data(ii + 442);
    const auto *ii_443 = buffer.data(ii + 443);
    const auto *ii_444 = buffer.data(ii + 444);
    const auto *ii_445 = buffer.data(ii + 445);
    const auto *ii_446 = buffer.data(ii + 446);
    const auto *ii_447 = buffer.data(ii + 447);
    const auto *ii_448 = buffer.data(ii + 448);
    const auto *ii_449 = buffer.data(ii + 449);
    const auto *ii_450 = buffer.data(ii + 450);
    const auto *ii_451 = buffer.data(ii + 451);
    const auto *ii_452 = buffer.data(ii + 452);
    const auto *ii_453 = buffer.data(ii + 453);
    const auto *ii_454 = buffer.data(ii + 454);
    const auto *ii_455 = buffer.data(ii + 455);
    const auto *ii_456 = buffer.data(ii + 456);
    const auto *ii_457 = buffer.data(ii + 457);
    const auto *ii_458 = buffer.data(ii + 458);
    const auto *ii_459 = buffer.data(ii + 459);
    const auto *ii_460 = buffer.data(ii + 460);
    const auto *ii_461 = buffer.data(ii + 461);
    const auto *ii_462 = buffer.data(ii + 462);
    const auto *ii_463 = buffer.data(ii + 463);
    const auto *ii_464 = buffer.data(ii + 464);
    const auto *ii_465 = buffer.data(ii + 465);
    const auto *ii_466 = buffer.data(ii + 466);
    const auto *ii_467 = buffer.data(ii + 467);
    const auto *ii_468 = buffer.data(ii + 468);
    const auto *ii_469 = buffer.data(ii + 469);
    const auto *ii_470 = buffer.data(ii + 470);
    const auto *ii_471 = buffer.data(ii + 471);
    const auto *ii_472 = buffer.data(ii + 472);
    const auto *ii_473 = buffer.data(ii + 473);
    const auto *ii_474 = buffer.data(ii + 474);
    const auto *ii_475 = buffer.data(ii + 475);
    const auto *ii_476 = buffer.data(ii + 476);
    const auto *ii_477 = buffer.data(ii + 477);
    const auto *ii_478 = buffer.data(ii + 478);
    const auto *ii_479 = buffer.data(ii + 479);
    const auto *ii_480 = buffer.data(ii + 480);
    const auto *ii_481 = buffer.data(ii + 481);
    const auto *ii_482 = buffer.data(ii + 482);
    const auto *ii_483 = buffer.data(ii + 483);
    const auto *ii_484 = buffer.data(ii + 484);
    const auto *ii_485 = buffer.data(ii + 485);
    const auto *ii_486 = buffer.data(ii + 486);

    const auto *li_663 = buffer.data(li + 663);
    const auto *li_664 = buffer.data(li + 664);
    const auto *li_665 = buffer.data(li + 665);
    const auto *li_666 = buffer.data(li + 666);
    const auto *li_667 = buffer.data(li + 667);
    const auto *li_668 = buffer.data(li + 668);
    const auto *li_669 = buffer.data(li + 669);
    const auto *li_670 = buffer.data(li + 670);
    const auto *li_671 = buffer.data(li + 671);
    const auto *li_672 = buffer.data(li + 672);
    const auto *li_673 = buffer.data(li + 673);
    const auto *li_674 = buffer.data(li + 674);
    const auto *li_675 = buffer.data(li + 675);
    const auto *li_676 = buffer.data(li + 676);
    const auto *li_677 = buffer.data(li + 677);
    const auto *li_678 = buffer.data(li + 678);
    const auto *li_679 = buffer.data(li + 679);
    const auto *li_680 = buffer.data(li + 680);
    const auto *li_681 = buffer.data(li + 681);
    const auto *li_682 = buffer.data(li + 682);
    const auto *li_683 = buffer.data(li + 683);
    const auto *li_684 = buffer.data(li + 684);
    const auto *li_685 = buffer.data(li + 685);
    const auto *li_686 = buffer.data(li + 686);
    const auto *li_687 = buffer.data(li + 687);
    const auto *li_688 = buffer.data(li + 688);
    const auto *li_689 = buffer.data(li + 689);
    const auto *li_690 = buffer.data(li + 690);
    const auto *li_691 = buffer.data(li + 691);
    const auto *li_692 = buffer.data(li + 692);
    const auto *li_693 = buffer.data(li + 693);
    const auto *li_694 = buffer.data(li + 694);
    const auto *li_695 = buffer.data(li + 695);
    const auto *li_696 = buffer.data(li + 696);
    const auto *li_697 = buffer.data(li + 697);
    const auto *li_698 = buffer.data(li + 698);
    const auto *li_699 = buffer.data(li + 699);
    const auto *li_700 = buffer.data(li + 700);
    const auto *li_701 = buffer.data(li + 701);
    const auto *li_702 = buffer.data(li + 702);
    const auto *li_703 = buffer.data(li + 703);
    const auto *li_704 = buffer.data(li + 704);
    const auto *li_705 = buffer.data(li + 705);
    const auto *li_706 = buffer.data(li + 706);
    const auto *li_707 = buffer.data(li + 707);
    const auto *li_708 = buffer.data(li + 708);
    const auto *li_709 = buffer.data(li + 709);
    const auto *li_710 = buffer.data(li + 710);
    const auto *li_711 = buffer.data(li + 711);
    const auto *li_712 = buffer.data(li + 712);
    const auto *li_713 = buffer.data(li + 713);
    const auto *li_714 = buffer.data(li + 714);
    const auto *li_715 = buffer.data(li + 715);
    const auto *li_716 = buffer.data(li + 716);
    const auto *li_717 = buffer.data(li + 717);
    const auto *li_718 = buffer.data(li + 718);
    const auto *li_719 = buffer.data(li + 719);
    const auto *li_720 = buffer.data(li + 720);
    const auto *li_721 = buffer.data(li + 721);
    const auto *li_722 = buffer.data(li + 722);
    const auto *li_723 = buffer.data(li + 723);
    const auto *li_724 = buffer.data(li + 724);
    const auto *li_725 = buffer.data(li + 725);
    const auto *li_726 = buffer.data(li + 726);
    const auto *li_727 = buffer.data(li + 727);
    const auto *li_728 = buffer.data(li + 728);
    const auto *li_729 = buffer.data(li + 729);
    const auto *li_730 = buffer.data(li + 730);
    const auto *li_731 = buffer.data(li + 731);
    const auto *li_732 = buffer.data(li + 732);
    const auto *li_733 = buffer.data(li + 733);
    const auto *li_734 = buffer.data(li + 734);
    const auto *li_735 = buffer.data(li + 735);
    const auto *li_736 = buffer.data(li + 736);
    const auto *li_737 = buffer.data(li + 737);
    const auto *li_738 = buffer.data(li + 738);
    const auto *li_739 = buffer.data(li + 739);
    const auto *li_740 = buffer.data(li + 740);
    const auto *li_741 = buffer.data(li + 741);
    const auto *li_742 = buffer.data(li + 742);
    const auto *li_743 = buffer.data(li + 743);
    const auto *li_744 = buffer.data(li + 744);
    const auto *li_745 = buffer.data(li + 745);
    const auto *li_746 = buffer.data(li + 746);
    const auto *li_747 = buffer.data(li + 747);
    const auto *li_748 = buffer.data(li + 748);
    const auto *li_749 = buffer.data(li + 749);
    const auto *li_750 = buffer.data(li + 750);
    const auto *li_751 = buffer.data(li + 751);
    const auto *li_752 = buffer.data(li + 752);
    const auto *li_753 = buffer.data(li + 753);
    const auto *li_754 = buffer.data(li + 754);
    const auto *li_755 = buffer.data(li + 755);
    const auto *li_784 = buffer.data(li + 784);
    const auto *li_785 = buffer.data(li + 785);
    const auto *li_786 = buffer.data(li + 786);
    const auto *li_787 = buffer.data(li + 787);
    const auto *li_788 = buffer.data(li + 788);
    const auto *li_789 = buffer.data(li + 789);
    const auto *li_790 = buffer.data(li + 790);
    const auto *li_791 = buffer.data(li + 791);
    const auto *li_792 = buffer.data(li + 792);
    const auto *li_793 = buffer.data(li + 793);
    const auto *li_794 = buffer.data(li + 794);
    const auto *li_795 = buffer.data(li + 795);
    const auto *li_796 = buffer.data(li + 796);
    const auto *li_797 = buffer.data(li + 797);
    const auto *li_798 = buffer.data(li + 798);
    const auto *li_799 = buffer.data(li + 799);
    const auto *li_800 = buffer.data(li + 800);
    const auto *li_801 = buffer.data(li + 801);
    const auto *li_802 = buffer.data(li + 802);
    const auto *li_803 = buffer.data(li + 803);
    const auto *li_804 = buffer.data(li + 804);
    const auto *li_805 = buffer.data(li + 805);
    const auto *li_806 = buffer.data(li + 806);
    const auto *li_807 = buffer.data(li + 807);
    const auto *li_808 = buffer.data(li + 808);
    const auto *li_809 = buffer.data(li + 809);
    const auto *li_810 = buffer.data(li + 810);
    const auto *li_811 = buffer.data(li + 811);
    const auto *li_812 = buffer.data(li + 812);
    const auto *li_813 = buffer.data(li + 813);
    const auto *li_814 = buffer.data(li + 814);
    const auto *li_815 = buffer.data(li + 815);
    const auto *li_816 = buffer.data(li + 816);
    const auto *li_817 = buffer.data(li + 817);
    const auto *li_818 = buffer.data(li + 818);
    const auto *li_819 = buffer.data(li + 819);
    const auto *li_820 = buffer.data(li + 820);
    const auto *li_821 = buffer.data(li + 821);
    const auto *li_822 = buffer.data(li + 822);
    const auto *li_823 = buffer.data(li + 823);
    const auto *li_824 = buffer.data(li + 824);
    const auto *li_825 = buffer.data(li + 825);
    const auto *li_826 = buffer.data(li + 826);
    const auto *li_827 = buffer.data(li + 827);
    const auto *li_828 = buffer.data(li + 828);
    const auto *li_829 = buffer.data(li + 829);
    const auto *li_830 = buffer.data(li + 830);
    const auto *li_831 = buffer.data(li + 831);
    const auto *li_832 = buffer.data(li + 832);
    const auto *li_833 = buffer.data(li + 833);
    const auto *li_834 = buffer.data(li + 834);
    const auto *li_835 = buffer.data(li + 835);
    const auto *li_836 = buffer.data(li + 836);
    const auto *li_837 = buffer.data(li + 837);
    const auto *li_838 = buffer.data(li + 838);
    const auto *li_839 = buffer.data(li + 839);
    const auto *li_840 = buffer.data(li + 840);
    const auto *li_841 = buffer.data(li + 841);
    const auto *li_842 = buffer.data(li + 842);
    const auto *li_843 = buffer.data(li + 843);
    const auto *li_844 = buffer.data(li + 844);
    const auto *li_845 = buffer.data(li + 845);
    const auto *li_846 = buffer.data(li + 846);
    const auto *li_847 = buffer.data(li + 847);
    const auto *li_848 = buffer.data(li + 848);
    const auto *li_849 = buffer.data(li + 849);
    const auto *li_850 = buffer.data(li + 850);

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ii_355, ii_356, ii_357, ii_358, \
                         ii_359, li_663, li_664, li_665, li_666, \
                         li_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -3.0 * ii_355[k]
                   + f_0 * li_663[k];

        t_496[k] = -3.0 * ii_356[k]
                   + f_0 * li_664[k];

        t_497[k] = -3.0 * ii_357[k]
                   + f_0 * li_665[k];

        t_498[k] = -3.0 * ii_358[k]
                   + f_0 * li_666[k];

        t_499[k] = -3.0 * ii_359[k]
                   + f_0 * li_667[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ii_360, ii_361, ii_362, ii_363, \
                         ii_364, li_668, li_669, li_670, li_671, \
                         li_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -3.0 * ii_360[k]
                   + f_0 * li_668[k];

        t_501[k] = -3.0 * ii_361[k]
                   + f_0 * li_669[k];

        t_502[k] = -3.0 * ii_362[k]
                   + f_0 * li_670[k];

        t_503[k] = -3.0 * ii_363[k]
                   + f_0 * li_671[k];

        t_504[k] = -2.0 * ii_364[k]
                   + f_0 * li_672[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ii_365, ii_366, ii_367, ii_368, \
                         ii_369, li_673, li_674, li_675, li_676, \
                         li_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * ii_365[k]
                   + f_0 * li_673[k];

        t_506[k] = -2.0 * ii_366[k]
                   + f_0 * li_674[k];

        t_507[k] = -2.0 * ii_367[k]
                   + f_0 * li_675[k];

        t_508[k] = -2.0 * ii_368[k]
                   + f_0 * li_676[k];

        t_509[k] = -2.0 * ii_369[k]
                   + f_0 * li_677[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ii_370, ii_371, ii_372, ii_373, \
                         ii_374, li_678, li_679, li_680, li_681, \
                         li_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -2.0 * ii_370[k]
                   + f_0 * li_678[k];

        t_511[k] = -2.0 * ii_371[k]
                   + f_0 * li_679[k];

        t_512[k] = -2.0 * ii_372[k]
                   + f_0 * li_680[k];

        t_513[k] = -2.0 * ii_373[k]
                   + f_0 * li_681[k];

        t_514[k] = -2.0 * ii_374[k]
                   + f_0 * li_682[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ii_375, ii_376, ii_377, ii_378, \
                         ii_379, li_683, li_684, li_685, li_686, \
                         li_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -2.0 * ii_375[k]
                   + f_0 * li_683[k];

        t_516[k] = -2.0 * ii_376[k]
                   + f_0 * li_684[k];

        t_517[k] = -2.0 * ii_377[k]
                   + f_0 * li_685[k];

        t_518[k] = -2.0 * ii_378[k]
                   + f_0 * li_686[k];

        t_519[k] = -2.0 * ii_379[k]
                   + f_0 * li_687[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ii_380, ii_381, ii_382, ii_383, \
                         ii_384, li_688, li_689, li_690, li_691, \
                         li_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -2.0 * ii_380[k]
                   + f_0 * li_688[k];

        t_521[k] = -2.0 * ii_381[k]
                   + f_0 * li_689[k];

        t_522[k] = -2.0 * ii_382[k]
                   + f_0 * li_690[k];

        t_523[k] = -2.0 * ii_383[k]
                   + f_0 * li_691[k];

        t_524[k] = -2.0 * ii_384[k]
                   + f_0 * li_692[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ii_385, ii_386, ii_387, ii_388, \
                         ii_389, li_693, li_694, li_695, li_696, \
                         li_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -2.0 * ii_385[k]
                   + f_0 * li_693[k];

        t_526[k] = -2.0 * ii_386[k]
                   + f_0 * li_694[k];

        t_527[k] = -2.0 * ii_387[k]
                   + f_0 * li_695[k];

        t_528[k] = -2.0 * ii_388[k]
                   + f_0 * li_696[k];

        t_529[k] = -2.0 * ii_389[k]
                   + f_0 * li_697[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ii_390, ii_391, ii_392, ii_393, \
                         ii_394, li_698, li_699, li_700, li_701, \
                         li_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -2.0 * ii_390[k]
                   + f_0 * li_698[k];

        t_531[k] = -2.0 * ii_391[k]
                   + f_0 * li_699[k];

        t_532[k] = -ii_392[k]
                   + f_0 * li_700[k];

        t_533[k] = -ii_393[k]
                   + f_0 * li_701[k];

        t_534[k] = -ii_394[k]
                   + f_0 * li_702[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ii_395, ii_396, ii_397, ii_398, \
                         ii_399, li_703, li_704, li_705, li_706, \
                         li_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -ii_395[k]
                   + f_0 * li_703[k];

        t_536[k] = -ii_396[k]
                   + f_0 * li_704[k];

        t_537[k] = -ii_397[k]
                   + f_0 * li_705[k];

        t_538[k] = -ii_398[k]
                   + f_0 * li_706[k];

        t_539[k] = -ii_399[k]
                   + f_0 * li_707[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ii_400, ii_401, ii_402, ii_403, \
                         ii_404, li_708, li_709, li_710, li_711, \
                         li_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -ii_400[k]
                   + f_0 * li_708[k];

        t_541[k] = -ii_401[k]
                   + f_0 * li_709[k];

        t_542[k] = -ii_402[k]
                   + f_0 * li_710[k];

        t_543[k] = -ii_403[k]
                   + f_0 * li_711[k];

        t_544[k] = -ii_404[k]
                   + f_0 * li_712[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ii_405, ii_406, ii_407, ii_408, \
                         ii_409, li_713, li_714, li_715, li_716, \
                         li_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -ii_405[k]
                   + f_0 * li_713[k];

        t_546[k] = -ii_406[k]
                   + f_0 * li_714[k];

        t_547[k] = -ii_407[k]
                   + f_0 * li_715[k];

        t_548[k] = -ii_408[k]
                   + f_0 * li_716[k];

        t_549[k] = -ii_409[k]
                   + f_0 * li_717[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ii_410, ii_411, ii_412, ii_413, \
                         ii_414, li_718, li_719, li_720, li_721, \
                         li_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -ii_410[k]
                   + f_0 * li_718[k];

        t_551[k] = -ii_411[k]
                   + f_0 * li_719[k];

        t_552[k] = -ii_412[k]
                   + f_0 * li_720[k];

        t_553[k] = -ii_413[k]
                   + f_0 * li_721[k];

        t_554[k] = -ii_414[k]
                   + f_0 * li_722[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ii_415, ii_416, ii_417, ii_418, \
                         ii_419, li_723, li_724, li_725, li_726, \
                         li_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -ii_415[k]
                   + f_0 * li_723[k];

        t_556[k] = -ii_416[k]
                   + f_0 * li_724[k];

        t_557[k] = -ii_417[k]
                   + f_0 * li_725[k];

        t_558[k] = -ii_418[k]
                   + f_0 * li_726[k];

        t_559[k] = -ii_419[k]
                   + f_0 * li_727[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, t_565, t_566, t_567, li_728, \
                         li_729, li_730, li_731, li_732, li_733, li_734, \
                         li_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_0 * li_728[k];

        t_561[k] = f_0 * li_729[k];

        t_562[k] = f_0 * li_730[k];

        t_563[k] = f_0 * li_731[k];

        t_564[k] = f_0 * li_732[k];

        t_565[k] = f_0 * li_733[k];

        t_566[k] = f_0 * li_734[k];

        t_567[k] = f_0 * li_735[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, t_572, t_573, t_574, t_575, li_736, \
                         li_737, li_738, li_739, li_740, li_741, li_742, \
                         li_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_0 * li_736[k];

        t_569[k] = f_0 * li_737[k];

        t_570[k] = f_0 * li_738[k];

        t_571[k] = f_0 * li_739[k];

        t_572[k] = f_0 * li_740[k];

        t_573[k] = f_0 * li_741[k];

        t_574[k] = f_0 * li_742[k];

        t_575[k] = f_0 * li_743[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, t_581, t_582, t_583, li_744, \
                         li_745, li_746, li_747, li_748, li_749, li_750, \
                         li_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_0 * li_744[k];

        t_577[k] = f_0 * li_745[k];

        t_578[k] = f_0 * li_746[k];

        t_579[k] = f_0 * li_747[k];

        t_580[k] = f_0 * li_748[k];

        t_581[k] = f_0 * li_749[k];

        t_582[k] = f_0 * li_750[k];

        t_583[k] = f_0 * li_751[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, t_588, t_589, ii_420, ii_421, li_752, \
                         li_753, li_754, li_755, li_784, li_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_0 * li_752[k];

        t_585[k] = f_0 * li_753[k];

        t_586[k] = f_0 * li_754[k];

        t_587[k] = f_0 * li_755[k];

        t_588[k] = -6.0 * ii_420[k]
                   + f_0 * li_784[k];

        t_589[k] = -6.0 * ii_421[k]
                   + f_0 * li_785[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, t_594, ii_422, ii_423, ii_424, ii_425, \
                         ii_426, li_786, li_787, li_788, li_789, \
                         li_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = -6.0 * ii_422[k]
                   + f_0 * li_786[k];

        t_591[k] = -6.0 * ii_423[k]
                   + f_0 * li_787[k];

        t_592[k] = -6.0 * ii_424[k]
                   + f_0 * li_788[k];

        t_593[k] = -6.0 * ii_425[k]
                   + f_0 * li_789[k];

        t_594[k] = -6.0 * ii_426[k]
                   + f_0 * li_790[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, t_599, ii_427, ii_428, ii_429, ii_430, \
                         ii_431, li_791, li_792, li_793, li_794, \
                         li_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = -6.0 * ii_427[k]
                   + f_0 * li_791[k];

        t_596[k] = -6.0 * ii_428[k]
                   + f_0 * li_792[k];

        t_597[k] = -6.0 * ii_429[k]
                   + f_0 * li_793[k];

        t_598[k] = -6.0 * ii_430[k]
                   + f_0 * li_794[k];

        t_599[k] = -6.0 * ii_431[k]
                   + f_0 * li_795[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, ii_432, ii_433, ii_434, ii_435, \
                         ii_436, li_796, li_797, li_798, li_799, \
                         li_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = -6.0 * ii_432[k]
                   + f_0 * li_796[k];

        t_601[k] = -6.0 * ii_433[k]
                   + f_0 * li_797[k];

        t_602[k] = -6.0 * ii_434[k]
                   + f_0 * li_798[k];

        t_603[k] = -6.0 * ii_435[k]
                   + f_0 * li_799[k];

        t_604[k] = -6.0 * ii_436[k]
                   + f_0 * li_800[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, ii_437, ii_438, ii_439, ii_440, \
                         ii_441, li_801, li_802, li_803, li_804, \
                         li_805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = -6.0 * ii_437[k]
                   + f_0 * li_801[k];

        t_606[k] = -6.0 * ii_438[k]
                   + f_0 * li_802[k];

        t_607[k] = -6.0 * ii_439[k]
                   + f_0 * li_803[k];

        t_608[k] = -6.0 * ii_440[k]
                   + f_0 * li_804[k];

        t_609[k] = -6.0 * ii_441[k]
                   + f_0 * li_805[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, ii_442, ii_443, ii_444, ii_445, \
                         ii_446, li_806, li_807, li_808, li_809, \
                         li_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = -6.0 * ii_442[k]
                   + f_0 * li_806[k];

        t_611[k] = -6.0 * ii_443[k]
                   + f_0 * li_807[k];

        t_612[k] = -6.0 * ii_444[k]
                   + f_0 * li_808[k];

        t_613[k] = -6.0 * ii_445[k]
                   + f_0 * li_809[k];

        t_614[k] = -6.0 * ii_446[k]
                   + f_0 * li_810[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, t_619, ii_447, ii_448, ii_449, ii_450, \
                         ii_451, li_811, li_812, li_813, li_814, \
                         li_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = -6.0 * ii_447[k]
                   + f_0 * li_811[k];

        t_616[k] = -5.0 * ii_448[k]
                   + f_0 * li_812[k];

        t_617[k] = -5.0 * ii_449[k]
                   + f_0 * li_813[k];

        t_618[k] = -5.0 * ii_450[k]
                   + f_0 * li_814[k];

        t_619[k] = -5.0 * ii_451[k]
                   + f_0 * li_815[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, ii_452, ii_453, ii_454, ii_455, \
                         ii_456, li_816, li_817, li_818, li_819, \
                         li_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = -5.0 * ii_452[k]
                   + f_0 * li_816[k];

        t_621[k] = -5.0 * ii_453[k]
                   + f_0 * li_817[k];

        t_622[k] = -5.0 * ii_454[k]
                   + f_0 * li_818[k];

        t_623[k] = -5.0 * ii_455[k]
                   + f_0 * li_819[k];

        t_624[k] = -5.0 * ii_456[k]
                   + f_0 * li_820[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, t_629, ii_457, ii_458, ii_459, ii_460, \
                         ii_461, li_821, li_822, li_823, li_824, \
                         li_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = -5.0 * ii_457[k]
                   + f_0 * li_821[k];

        t_626[k] = -5.0 * ii_458[k]
                   + f_0 * li_822[k];

        t_627[k] = -5.0 * ii_459[k]
                   + f_0 * li_823[k];

        t_628[k] = -5.0 * ii_460[k]
                   + f_0 * li_824[k];

        t_629[k] = -5.0 * ii_461[k]
                   + f_0 * li_825[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, t_634, ii_462, ii_463, ii_464, ii_465, \
                         ii_466, li_826, li_827, li_828, li_829, \
                         li_830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = -5.0 * ii_462[k]
                   + f_0 * li_826[k];

        t_631[k] = -5.0 * ii_463[k]
                   + f_0 * li_827[k];

        t_632[k] = -5.0 * ii_464[k]
                   + f_0 * li_828[k];

        t_633[k] = -5.0 * ii_465[k]
                   + f_0 * li_829[k];

        t_634[k] = -5.0 * ii_466[k]
                   + f_0 * li_830[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, ii_467, ii_468, ii_469, ii_470, \
                         ii_471, li_831, li_832, li_833, li_834, \
                         li_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = -5.0 * ii_467[k]
                   + f_0 * li_831[k];

        t_636[k] = -5.0 * ii_468[k]
                   + f_0 * li_832[k];

        t_637[k] = -5.0 * ii_469[k]
                   + f_0 * li_833[k];

        t_638[k] = -5.0 * ii_470[k]
                   + f_0 * li_834[k];

        t_639[k] = -5.0 * ii_471[k]
                   + f_0 * li_835[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, t_644, ii_472, ii_473, ii_474, ii_475, \
                         ii_476, li_836, li_837, li_838, li_839, \
                         li_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = -5.0 * ii_472[k]
                   + f_0 * li_836[k];

        t_641[k] = -5.0 * ii_473[k]
                   + f_0 * li_837[k];

        t_642[k] = -5.0 * ii_474[k]
                   + f_0 * li_838[k];

        t_643[k] = -5.0 * ii_475[k]
                   + f_0 * li_839[k];

        t_644[k] = -4.0 * ii_476[k]
                   + f_0 * li_840[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, t_649, ii_477, ii_478, ii_479, ii_480, \
                         ii_481, li_841, li_842, li_843, li_844, \
                         li_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = -4.0 * ii_477[k]
                   + f_0 * li_841[k];

        t_646[k] = -4.0 * ii_478[k]
                   + f_0 * li_842[k];

        t_647[k] = -4.0 * ii_479[k]
                   + f_0 * li_843[k];

        t_648[k] = -4.0 * ii_480[k]
                   + f_0 * li_844[k];

        t_649[k] = -4.0 * ii_481[k]
                   + f_0 * li_845[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, ii_482, ii_483, ii_484, ii_485, \
                         ii_486, li_846, li_847, li_848, li_849, \
                         li_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = -4.0 * ii_482[k]
                   + f_0 * li_846[k];

        t_651[k] = -4.0 * ii_483[k]
                   + f_0 * li_847[k];

        t_652[k] = -4.0 * ii_484[k]
                   + f_0 * li_848[k];

        t_653[k] = -4.0 * ii_485[k]
                   + f_0 * li_849[k];

        t_654[k] = -4.0 * ii_486[k]
                   + f_0 * li_850[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_1_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_487 = buffer.data(ii + 487);
    const auto *ii_488 = buffer.data(ii + 488);
    const auto *ii_489 = buffer.data(ii + 489);
    const auto *ii_490 = buffer.data(ii + 490);
    const auto *ii_491 = buffer.data(ii + 491);
    const auto *ii_492 = buffer.data(ii + 492);
    const auto *ii_493 = buffer.data(ii + 493);
    const auto *ii_494 = buffer.data(ii + 494);
    const auto *ii_495 = buffer.data(ii + 495);
    const auto *ii_496 = buffer.data(ii + 496);
    const auto *ii_497 = buffer.data(ii + 497);
    const auto *ii_498 = buffer.data(ii + 498);
    const auto *ii_499 = buffer.data(ii + 499);
    const auto *ii_500 = buffer.data(ii + 500);
    const auto *ii_501 = buffer.data(ii + 501);
    const auto *ii_502 = buffer.data(ii + 502);
    const auto *ii_503 = buffer.data(ii + 503);
    const auto *ii_504 = buffer.data(ii + 504);
    const auto *ii_505 = buffer.data(ii + 505);
    const auto *ii_506 = buffer.data(ii + 506);
    const auto *ii_507 = buffer.data(ii + 507);
    const auto *ii_508 = buffer.data(ii + 508);
    const auto *ii_509 = buffer.data(ii + 509);
    const auto *ii_510 = buffer.data(ii + 510);
    const auto *ii_511 = buffer.data(ii + 511);
    const auto *ii_512 = buffer.data(ii + 512);
    const auto *ii_513 = buffer.data(ii + 513);
    const auto *ii_514 = buffer.data(ii + 514);
    const auto *ii_515 = buffer.data(ii + 515);
    const auto *ii_516 = buffer.data(ii + 516);
    const auto *ii_517 = buffer.data(ii + 517);
    const auto *ii_518 = buffer.data(ii + 518);
    const auto *ii_519 = buffer.data(ii + 519);
    const auto *ii_520 = buffer.data(ii + 520);
    const auto *ii_521 = buffer.data(ii + 521);
    const auto *ii_522 = buffer.data(ii + 522);
    const auto *ii_523 = buffer.data(ii + 523);
    const auto *ii_524 = buffer.data(ii + 524);
    const auto *ii_525 = buffer.data(ii + 525);
    const auto *ii_526 = buffer.data(ii + 526);
    const auto *ii_527 = buffer.data(ii + 527);
    const auto *ii_528 = buffer.data(ii + 528);
    const auto *ii_529 = buffer.data(ii + 529);
    const auto *ii_530 = buffer.data(ii + 530);
    const auto *ii_531 = buffer.data(ii + 531);
    const auto *ii_532 = buffer.data(ii + 532);
    const auto *ii_533 = buffer.data(ii + 533);
    const auto *ii_534 = buffer.data(ii + 534);
    const auto *ii_535 = buffer.data(ii + 535);
    const auto *ii_536 = buffer.data(ii + 536);
    const auto *ii_537 = buffer.data(ii + 537);
    const auto *ii_538 = buffer.data(ii + 538);
    const auto *ii_539 = buffer.data(ii + 539);
    const auto *ii_540 = buffer.data(ii + 540);
    const auto *ii_541 = buffer.data(ii + 541);
    const auto *ii_542 = buffer.data(ii + 542);
    const auto *ii_543 = buffer.data(ii + 543);
    const auto *ii_544 = buffer.data(ii + 544);
    const auto *ii_545 = buffer.data(ii + 545);
    const auto *ii_546 = buffer.data(ii + 546);
    const auto *ii_547 = buffer.data(ii + 547);
    const auto *ii_548 = buffer.data(ii + 548);
    const auto *ii_549 = buffer.data(ii + 549);
    const auto *ii_550 = buffer.data(ii + 550);
    const auto *ii_551 = buffer.data(ii + 551);
    const auto *ii_552 = buffer.data(ii + 552);
    const auto *ii_553 = buffer.data(ii + 553);
    const auto *ii_554 = buffer.data(ii + 554);
    const auto *ii_555 = buffer.data(ii + 555);
    const auto *ii_556 = buffer.data(ii + 556);
    const auto *ii_557 = buffer.data(ii + 557);
    const auto *ii_558 = buffer.data(ii + 558);
    const auto *ii_559 = buffer.data(ii + 559);
    const auto *ii_560 = buffer.data(ii + 560);
    const auto *ii_561 = buffer.data(ii + 561);
    const auto *ii_562 = buffer.data(ii + 562);
    const auto *ii_563 = buffer.data(ii + 563);
    const auto *ii_564 = buffer.data(ii + 564);
    const auto *ii_565 = buffer.data(ii + 565);
    const auto *ii_566 = buffer.data(ii + 566);
    const auto *ii_567 = buffer.data(ii + 567);
    const auto *ii_568 = buffer.data(ii + 568);
    const auto *ii_569 = buffer.data(ii + 569);
    const auto *ii_570 = buffer.data(ii + 570);
    const auto *ii_571 = buffer.data(ii + 571);
    const auto *ii_572 = buffer.data(ii + 572);
    const auto *ii_573 = buffer.data(ii + 573);
    const auto *ii_574 = buffer.data(ii + 574);
    const auto *ii_575 = buffer.data(ii + 575);
    const auto *ii_576 = buffer.data(ii + 576);
    const auto *ii_577 = buffer.data(ii + 577);
    const auto *ii_578 = buffer.data(ii + 578);
    const auto *ii_579 = buffer.data(ii + 579);
    const auto *ii_580 = buffer.data(ii + 580);
    const auto *ii_581 = buffer.data(ii + 581);
    const auto *ii_582 = buffer.data(ii + 582);
    const auto *ii_583 = buffer.data(ii + 583);
    const auto *ii_584 = buffer.data(ii + 584);
    const auto *ii_585 = buffer.data(ii + 585);
    const auto *ii_586 = buffer.data(ii + 586);
    const auto *ii_587 = buffer.data(ii + 587);
    const auto *ii_588 = buffer.data(ii + 588);
    const auto *ii_589 = buffer.data(ii + 589);
    const auto *ii_590 = buffer.data(ii + 590);
    const auto *ii_591 = buffer.data(ii + 591);
    const auto *ii_592 = buffer.data(ii + 592);
    const auto *ii_593 = buffer.data(ii + 593);
    const auto *ii_594 = buffer.data(ii + 594);
    const auto *ii_595 = buffer.data(ii + 595);
    const auto *ii_596 = buffer.data(ii + 596);
    const auto *ii_597 = buffer.data(ii + 597);
    const auto *ii_598 = buffer.data(ii + 598);
    const auto *ii_599 = buffer.data(ii + 599);
    const auto *ii_600 = buffer.data(ii + 600);
    const auto *ii_601 = buffer.data(ii + 601);
    const auto *ii_602 = buffer.data(ii + 602);
    const auto *ii_603 = buffer.data(ii + 603);
    const auto *ii_604 = buffer.data(ii + 604);
    const auto *ii_605 = buffer.data(ii + 605);
    const auto *ii_606 = buffer.data(ii + 606);
    const auto *ii_607 = buffer.data(ii + 607);
    const auto *ii_608 = buffer.data(ii + 608);
    const auto *ii_609 = buffer.data(ii + 609);
    const auto *ii_610 = buffer.data(ii + 610);
    const auto *ii_611 = buffer.data(ii + 611);
    const auto *ii_612 = buffer.data(ii + 612);
    const auto *ii_613 = buffer.data(ii + 613);
    const auto *ii_614 = buffer.data(ii + 614);
    const auto *ii_615 = buffer.data(ii + 615);
    const auto *ii_616 = buffer.data(ii + 616);
    const auto *ii_617 = buffer.data(ii + 617);
    const auto *ii_618 = buffer.data(ii + 618);

    const auto *li_851 = buffer.data(li + 851);
    const auto *li_852 = buffer.data(li + 852);
    const auto *li_853 = buffer.data(li + 853);
    const auto *li_854 = buffer.data(li + 854);
    const auto *li_855 = buffer.data(li + 855);
    const auto *li_856 = buffer.data(li + 856);
    const auto *li_857 = buffer.data(li + 857);
    const auto *li_858 = buffer.data(li + 858);
    const auto *li_859 = buffer.data(li + 859);
    const auto *li_860 = buffer.data(li + 860);
    const auto *li_861 = buffer.data(li + 861);
    const auto *li_862 = buffer.data(li + 862);
    const auto *li_863 = buffer.data(li + 863);
    const auto *li_864 = buffer.data(li + 864);
    const auto *li_865 = buffer.data(li + 865);
    const auto *li_866 = buffer.data(li + 866);
    const auto *li_867 = buffer.data(li + 867);
    const auto *li_868 = buffer.data(li + 868);
    const auto *li_869 = buffer.data(li + 869);
    const auto *li_870 = buffer.data(li + 870);
    const auto *li_871 = buffer.data(li + 871);
    const auto *li_872 = buffer.data(li + 872);
    const auto *li_873 = buffer.data(li + 873);
    const auto *li_874 = buffer.data(li + 874);
    const auto *li_875 = buffer.data(li + 875);
    const auto *li_876 = buffer.data(li + 876);
    const auto *li_877 = buffer.data(li + 877);
    const auto *li_878 = buffer.data(li + 878);
    const auto *li_879 = buffer.data(li + 879);
    const auto *li_880 = buffer.data(li + 880);
    const auto *li_881 = buffer.data(li + 881);
    const auto *li_882 = buffer.data(li + 882);
    const auto *li_883 = buffer.data(li + 883);
    const auto *li_884 = buffer.data(li + 884);
    const auto *li_885 = buffer.data(li + 885);
    const auto *li_886 = buffer.data(li + 886);
    const auto *li_887 = buffer.data(li + 887);
    const auto *li_888 = buffer.data(li + 888);
    const auto *li_889 = buffer.data(li + 889);
    const auto *li_890 = buffer.data(li + 890);
    const auto *li_891 = buffer.data(li + 891);
    const auto *li_892 = buffer.data(li + 892);
    const auto *li_893 = buffer.data(li + 893);
    const auto *li_894 = buffer.data(li + 894);
    const auto *li_895 = buffer.data(li + 895);
    const auto *li_896 = buffer.data(li + 896);
    const auto *li_897 = buffer.data(li + 897);
    const auto *li_898 = buffer.data(li + 898);
    const auto *li_899 = buffer.data(li + 899);
    const auto *li_900 = buffer.data(li + 900);
    const auto *li_901 = buffer.data(li + 901);
    const auto *li_902 = buffer.data(li + 902);
    const auto *li_903 = buffer.data(li + 903);
    const auto *li_904 = buffer.data(li + 904);
    const auto *li_905 = buffer.data(li + 905);
    const auto *li_906 = buffer.data(li + 906);
    const auto *li_907 = buffer.data(li + 907);
    const auto *li_908 = buffer.data(li + 908);
    const auto *li_909 = buffer.data(li + 909);
    const auto *li_910 = buffer.data(li + 910);
    const auto *li_911 = buffer.data(li + 911);
    const auto *li_912 = buffer.data(li + 912);
    const auto *li_913 = buffer.data(li + 913);
    const auto *li_914 = buffer.data(li + 914);
    const auto *li_915 = buffer.data(li + 915);
    const auto *li_916 = buffer.data(li + 916);
    const auto *li_917 = buffer.data(li + 917);
    const auto *li_918 = buffer.data(li + 918);
    const auto *li_919 = buffer.data(li + 919);
    const auto *li_920 = buffer.data(li + 920);
    const auto *li_921 = buffer.data(li + 921);
    const auto *li_922 = buffer.data(li + 922);
    const auto *li_923 = buffer.data(li + 923);
    const auto *li_924 = buffer.data(li + 924);
    const auto *li_925 = buffer.data(li + 925);
    const auto *li_926 = buffer.data(li + 926);
    const auto *li_927 = buffer.data(li + 927);
    const auto *li_928 = buffer.data(li + 928);
    const auto *li_929 = buffer.data(li + 929);
    const auto *li_930 = buffer.data(li + 930);
    const auto *li_931 = buffer.data(li + 931);
    const auto *li_932 = buffer.data(li + 932);
    const auto *li_933 = buffer.data(li + 933);
    const auto *li_934 = buffer.data(li + 934);
    const auto *li_935 = buffer.data(li + 935);
    const auto *li_936 = buffer.data(li + 936);
    const auto *li_937 = buffer.data(li + 937);
    const auto *li_938 = buffer.data(li + 938);
    const auto *li_939 = buffer.data(li + 939);
    const auto *li_940 = buffer.data(li + 940);
    const auto *li_941 = buffer.data(li + 941);
    const auto *li_942 = buffer.data(li + 942);
    const auto *li_943 = buffer.data(li + 943);
    const auto *li_944 = buffer.data(li + 944);
    const auto *li_945 = buffer.data(li + 945);
    const auto *li_946 = buffer.data(li + 946);
    const auto *li_947 = buffer.data(li + 947);
    const auto *li_948 = buffer.data(li + 948);
    const auto *li_949 = buffer.data(li + 949);
    const auto *li_950 = buffer.data(li + 950);
    const auto *li_951 = buffer.data(li + 951);
    const auto *li_952 = buffer.data(li + 952);
    const auto *li_953 = buffer.data(li + 953);
    const auto *li_954 = buffer.data(li + 954);
    const auto *li_955 = buffer.data(li + 955);
    const auto *li_956 = buffer.data(li + 956);
    const auto *li_957 = buffer.data(li + 957);
    const auto *li_958 = buffer.data(li + 958);
    const auto *li_959 = buffer.data(li + 959);
    const auto *li_960 = buffer.data(li + 960);
    const auto *li_961 = buffer.data(li + 961);
    const auto *li_962 = buffer.data(li + 962);
    const auto *li_963 = buffer.data(li + 963);
    const auto *li_964 = buffer.data(li + 964);
    const auto *li_965 = buffer.data(li + 965);
    const auto *li_966 = buffer.data(li + 966);
    const auto *li_967 = buffer.data(li + 967);
    const auto *li_968 = buffer.data(li + 968);
    const auto *li_969 = buffer.data(li + 969);
    const auto *li_970 = buffer.data(li + 970);
    const auto *li_971 = buffer.data(li + 971);
    const auto *li_972 = buffer.data(li + 972);
    const auto *li_973 = buffer.data(li + 973);
    const auto *li_974 = buffer.data(li + 974);
    const auto *li_975 = buffer.data(li + 975);
    const auto *li_976 = buffer.data(li + 976);
    const auto *li_977 = buffer.data(li + 977);
    const auto *li_978 = buffer.data(li + 978);
    const auto *li_979 = buffer.data(li + 979);
    const auto *li_1008 = buffer.data(li + 1008);
    const auto *li_1009 = buffer.data(li + 1009);
    const auto *li_1010 = buffer.data(li + 1010);
    const auto *li_1011 = buffer.data(li + 1011);
    const auto *li_1012 = buffer.data(li + 1012);
    const auto *li_1013 = buffer.data(li + 1013);
    const auto *li_1014 = buffer.data(li + 1014);
    const auto *li_1015 = buffer.data(li + 1015);
    const auto *li_1016 = buffer.data(li + 1016);
    const auto *li_1017 = buffer.data(li + 1017);
    const auto *li_1018 = buffer.data(li + 1018);
    const auto *li_1019 = buffer.data(li + 1019);
    const auto *li_1020 = buffer.data(li + 1020);
    const auto *li_1021 = buffer.data(li + 1021);
    const auto *li_1022 = buffer.data(li + 1022);
    const auto *li_1023 = buffer.data(li + 1023);
    const auto *li_1024 = buffer.data(li + 1024);
    const auto *li_1025 = buffer.data(li + 1025);
    const auto *li_1026 = buffer.data(li + 1026);
    const auto *li_1027 = buffer.data(li + 1027);
    const auto *li_1028 = buffer.data(li + 1028);
    const auto *li_1029 = buffer.data(li + 1029);
    const auto *li_1030 = buffer.data(li + 1030);
    const auto *li_1031 = buffer.data(li + 1031);
    const auto *li_1032 = buffer.data(li + 1032);
    const auto *li_1033 = buffer.data(li + 1033);
    const auto *li_1034 = buffer.data(li + 1034);
    const auto *li_1035 = buffer.data(li + 1035);
    const auto *li_1036 = buffer.data(li + 1036);
    const auto *li_1037 = buffer.data(li + 1037);
    const auto *li_1038 = buffer.data(li + 1038);

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, ii_487, ii_488, ii_489, ii_490, \
                         ii_491, li_851, li_852, li_853, li_854, \
                         li_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = -4.0 * ii_487[k]
                   + f_0 * li_851[k];

        t_656[k] = -4.0 * ii_488[k]
                   + f_0 * li_852[k];

        t_657[k] = -4.0 * ii_489[k]
                   + f_0 * li_853[k];

        t_658[k] = -4.0 * ii_490[k]
                   + f_0 * li_854[k];

        t_659[k] = -4.0 * ii_491[k]
                   + f_0 * li_855[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, ii_492, ii_493, ii_494, ii_495, \
                         ii_496, li_856, li_857, li_858, li_859, \
                         li_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = -4.0 * ii_492[k]
                   + f_0 * li_856[k];

        t_661[k] = -4.0 * ii_493[k]
                   + f_0 * li_857[k];

        t_662[k] = -4.0 * ii_494[k]
                   + f_0 * li_858[k];

        t_663[k] = -4.0 * ii_495[k]
                   + f_0 * li_859[k];

        t_664[k] = -4.0 * ii_496[k]
                   + f_0 * li_860[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, ii_497, ii_498, ii_499, ii_500, \
                         ii_501, li_861, li_862, li_863, li_864, \
                         li_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = -4.0 * ii_497[k]
                   + f_0 * li_861[k];

        t_666[k] = -4.0 * ii_498[k]
                   + f_0 * li_862[k];

        t_667[k] = -4.0 * ii_499[k]
                   + f_0 * li_863[k];

        t_668[k] = -4.0 * ii_500[k]
                   + f_0 * li_864[k];

        t_669[k] = -4.0 * ii_501[k]
                   + f_0 * li_865[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, ii_502, ii_503, ii_504, ii_505, \
                         ii_506, li_866, li_867, li_868, li_869, \
                         li_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = -4.0 * ii_502[k]
                   + f_0 * li_866[k];

        t_671[k] = -4.0 * ii_503[k]
                   + f_0 * li_867[k];

        t_672[k] = -3.0 * ii_504[k]
                   + f_0 * li_868[k];

        t_673[k] = -3.0 * ii_505[k]
                   + f_0 * li_869[k];

        t_674[k] = -3.0 * ii_506[k]
                   + f_0 * li_870[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, t_679, ii_507, ii_508, ii_509, ii_510, \
                         ii_511, li_871, li_872, li_873, li_874, \
                         li_875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = -3.0 * ii_507[k]
                   + f_0 * li_871[k];

        t_676[k] = -3.0 * ii_508[k]
                   + f_0 * li_872[k];

        t_677[k] = -3.0 * ii_509[k]
                   + f_0 * li_873[k];

        t_678[k] = -3.0 * ii_510[k]
                   + f_0 * li_874[k];

        t_679[k] = -3.0 * ii_511[k]
                   + f_0 * li_875[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, ii_512, ii_513, ii_514, ii_515, \
                         ii_516, li_876, li_877, li_878, li_879, \
                         li_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = -3.0 * ii_512[k]
                   + f_0 * li_876[k];

        t_681[k] = -3.0 * ii_513[k]
                   + f_0 * li_877[k];

        t_682[k] = -3.0 * ii_514[k]
                   + f_0 * li_878[k];

        t_683[k] = -3.0 * ii_515[k]
                   + f_0 * li_879[k];

        t_684[k] = -3.0 * ii_516[k]
                   + f_0 * li_880[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, ii_517, ii_518, ii_519, ii_520, \
                         ii_521, li_881, li_882, li_883, li_884, \
                         li_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = -3.0 * ii_517[k]
                   + f_0 * li_881[k];

        t_686[k] = -3.0 * ii_518[k]
                   + f_0 * li_882[k];

        t_687[k] = -3.0 * ii_519[k]
                   + f_0 * li_883[k];

        t_688[k] = -3.0 * ii_520[k]
                   + f_0 * li_884[k];

        t_689[k] = -3.0 * ii_521[k]
                   + f_0 * li_885[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, ii_522, ii_523, ii_524, ii_525, \
                         ii_526, li_886, li_887, li_888, li_889, \
                         li_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = -3.0 * ii_522[k]
                   + f_0 * li_886[k];

        t_691[k] = -3.0 * ii_523[k]
                   + f_0 * li_887[k];

        t_692[k] = -3.0 * ii_524[k]
                   + f_0 * li_888[k];

        t_693[k] = -3.0 * ii_525[k]
                   + f_0 * li_889[k];

        t_694[k] = -3.0 * ii_526[k]
                   + f_0 * li_890[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, ii_527, ii_528, ii_529, ii_530, \
                         ii_531, li_891, li_892, li_893, li_894, \
                         li_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = -3.0 * ii_527[k]
                   + f_0 * li_891[k];

        t_696[k] = -3.0 * ii_528[k]
                   + f_0 * li_892[k];

        t_697[k] = -3.0 * ii_529[k]
                   + f_0 * li_893[k];

        t_698[k] = -3.0 * ii_530[k]
                   + f_0 * li_894[k];

        t_699[k] = -3.0 * ii_531[k]
                   + f_0 * li_895[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, ii_532, ii_533, ii_534, ii_535, \
                         ii_536, li_896, li_897, li_898, li_899, \
                         li_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = -2.0 * ii_532[k]
                   + f_0 * li_896[k];

        t_701[k] = -2.0 * ii_533[k]
                   + f_0 * li_897[k];

        t_702[k] = -2.0 * ii_534[k]
                   + f_0 * li_898[k];

        t_703[k] = -2.0 * ii_535[k]
                   + f_0 * li_899[k];

        t_704[k] = -2.0 * ii_536[k]
                   + f_0 * li_900[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, ii_537, ii_538, ii_539, ii_540, \
                         ii_541, li_901, li_902, li_903, li_904, \
                         li_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = -2.0 * ii_537[k]
                   + f_0 * li_901[k];

        t_706[k] = -2.0 * ii_538[k]
                   + f_0 * li_902[k];

        t_707[k] = -2.0 * ii_539[k]
                   + f_0 * li_903[k];

        t_708[k] = -2.0 * ii_540[k]
                   + f_0 * li_904[k];

        t_709[k] = -2.0 * ii_541[k]
                   + f_0 * li_905[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, t_714, ii_542, ii_543, ii_544, ii_545, \
                         ii_546, li_906, li_907, li_908, li_909, \
                         li_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = -2.0 * ii_542[k]
                   + f_0 * li_906[k];

        t_711[k] = -2.0 * ii_543[k]
                   + f_0 * li_907[k];

        t_712[k] = -2.0 * ii_544[k]
                   + f_0 * li_908[k];

        t_713[k] = -2.0 * ii_545[k]
                   + f_0 * li_909[k];

        t_714[k] = -2.0 * ii_546[k]
                   + f_0 * li_910[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, ii_547, ii_548, ii_549, ii_550, \
                         ii_551, li_911, li_912, li_913, li_914, \
                         li_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = -2.0 * ii_547[k]
                   + f_0 * li_911[k];

        t_716[k] = -2.0 * ii_548[k]
                   + f_0 * li_912[k];

        t_717[k] = -2.0 * ii_549[k]
                   + f_0 * li_913[k];

        t_718[k] = -2.0 * ii_550[k]
                   + f_0 * li_914[k];

        t_719[k] = -2.0 * ii_551[k]
                   + f_0 * li_915[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, ii_552, ii_553, ii_554, ii_555, \
                         ii_556, li_916, li_917, li_918, li_919, \
                         li_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = -2.0 * ii_552[k]
                   + f_0 * li_916[k];

        t_721[k] = -2.0 * ii_553[k]
                   + f_0 * li_917[k];

        t_722[k] = -2.0 * ii_554[k]
                   + f_0 * li_918[k];

        t_723[k] = -2.0 * ii_555[k]
                   + f_0 * li_919[k];

        t_724[k] = -2.0 * ii_556[k]
                   + f_0 * li_920[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, ii_557, ii_558, ii_559, ii_560, \
                         ii_561, li_921, li_922, li_923, li_924, \
                         li_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = -2.0 * ii_557[k]
                   + f_0 * li_921[k];

        t_726[k] = -2.0 * ii_558[k]
                   + f_0 * li_922[k];

        t_727[k] = -2.0 * ii_559[k]
                   + f_0 * li_923[k];

        t_728[k] = -ii_560[k]
                   + f_0 * li_924[k];

        t_729[k] = -ii_561[k]
                   + f_0 * li_925[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, ii_562, ii_563, ii_564, ii_565, \
                         ii_566, li_926, li_927, li_928, li_929, \
                         li_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = -ii_562[k]
                   + f_0 * li_926[k];

        t_731[k] = -ii_563[k]
                   + f_0 * li_927[k];

        t_732[k] = -ii_564[k]
                   + f_0 * li_928[k];

        t_733[k] = -ii_565[k]
                   + f_0 * li_929[k];

        t_734[k] = -ii_566[k]
                   + f_0 * li_930[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, ii_567, ii_568, ii_569, ii_570, \
                         ii_571, li_931, li_932, li_933, li_934, \
                         li_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = -ii_567[k]
                   + f_0 * li_931[k];

        t_736[k] = -ii_568[k]
                   + f_0 * li_932[k];

        t_737[k] = -ii_569[k]
                   + f_0 * li_933[k];

        t_738[k] = -ii_570[k]
                   + f_0 * li_934[k];

        t_739[k] = -ii_571[k]
                   + f_0 * li_935[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, t_744, ii_572, ii_573, ii_574, ii_575, \
                         ii_576, li_936, li_937, li_938, li_939, \
                         li_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = -ii_572[k]
                   + f_0 * li_936[k];

        t_741[k] = -ii_573[k]
                   + f_0 * li_937[k];

        t_742[k] = -ii_574[k]
                   + f_0 * li_938[k];

        t_743[k] = -ii_575[k]
                   + f_0 * li_939[k];

        t_744[k] = -ii_576[k]
                   + f_0 * li_940[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, ii_577, ii_578, ii_579, ii_580, \
                         ii_581, li_941, li_942, li_943, li_944, \
                         li_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = -ii_577[k]
                   + f_0 * li_941[k];

        t_746[k] = -ii_578[k]
                   + f_0 * li_942[k];

        t_747[k] = -ii_579[k]
                   + f_0 * li_943[k];

        t_748[k] = -ii_580[k]
                   + f_0 * li_944[k];

        t_749[k] = -ii_581[k]
                   + f_0 * li_945[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, ii_582, ii_583, ii_584, ii_585, \
                         ii_586, li_946, li_947, li_948, li_949, \
                         li_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = -ii_582[k]
                   + f_0 * li_946[k];

        t_751[k] = -ii_583[k]
                   + f_0 * li_947[k];

        t_752[k] = -ii_584[k]
                   + f_0 * li_948[k];

        t_753[k] = -ii_585[k]
                   + f_0 * li_949[k];

        t_754[k] = -ii_586[k]
                   + f_0 * li_950[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, t_759, t_760, t_761, ii_587, li_951, \
                         li_952, li_953, li_954, li_955, li_956, \
                         li_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = -ii_587[k]
                   + f_0 * li_951[k];

        t_756[k] = f_0 * li_952[k];

        t_757[k] = f_0 * li_953[k];

        t_758[k] = f_0 * li_954[k];

        t_759[k] = f_0 * li_955[k];

        t_760[k] = f_0 * li_956[k];

        t_761[k] = f_0 * li_957[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, t_767, t_768, t_769, li_958, \
                         li_959, li_960, li_961, li_962, li_963, li_964, \
                         li_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_0 * li_958[k];

        t_763[k] = f_0 * li_959[k];

        t_764[k] = f_0 * li_960[k];

        t_765[k] = f_0 * li_961[k];

        t_766[k] = f_0 * li_962[k];

        t_767[k] = f_0 * li_963[k];

        t_768[k] = f_0 * li_964[k];

        t_769[k] = f_0 * li_965[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, t_775, t_776, t_777, li_966, \
                         li_967, li_968, li_969, li_970, li_971, li_972, \
                         li_973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_0 * li_966[k];

        t_771[k] = f_0 * li_967[k];

        t_772[k] = f_0 * li_968[k];

        t_773[k] = f_0 * li_969[k];

        t_774[k] = f_0 * li_970[k];

        t_775[k] = f_0 * li_971[k];

        t_776[k] = f_0 * li_972[k];

        t_777[k] = f_0 * li_973[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, t_782, t_783, t_784, ii_588, li_974, \
                         li_975, li_976, li_977, li_978, li_979, \
                         li_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_0 * li_974[k];

        t_779[k] = f_0 * li_975[k];

        t_780[k] = f_0 * li_976[k];

        t_781[k] = f_0 * li_977[k];

        t_782[k] = f_0 * li_978[k];

        t_783[k] = f_0 * li_979[k];

        t_784[k] = -7.0 * ii_588[k]
                   + f_0 * li_1008[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, ii_589, ii_590, ii_591, ii_592, \
                         ii_593, li_1009, li_1010, li_1011, li_1012, \
                         li_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = -7.0 * ii_589[k]
                   + f_0 * li_1009[k];

        t_786[k] = -7.0 * ii_590[k]
                   + f_0 * li_1010[k];

        t_787[k] = -7.0 * ii_591[k]
                   + f_0 * li_1011[k];

        t_788[k] = -7.0 * ii_592[k]
                   + f_0 * li_1012[k];

        t_789[k] = -7.0 * ii_593[k]
                   + f_0 * li_1013[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, ii_594, ii_595, ii_596, ii_597, \
                         ii_598, li_1014, li_1015, li_1016, li_1017, \
                         li_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = -7.0 * ii_594[k]
                   + f_0 * li_1014[k];

        t_791[k] = -7.0 * ii_595[k]
                   + f_0 * li_1015[k];

        t_792[k] = -7.0 * ii_596[k]
                   + f_0 * li_1016[k];

        t_793[k] = -7.0 * ii_597[k]
                   + f_0 * li_1017[k];

        t_794[k] = -7.0 * ii_598[k]
                   + f_0 * li_1018[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, ii_599, ii_600, ii_601, ii_602, \
                         ii_603, li_1019, li_1020, li_1021, li_1022, \
                         li_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = -7.0 * ii_599[k]
                   + f_0 * li_1019[k];

        t_796[k] = -7.0 * ii_600[k]
                   + f_0 * li_1020[k];

        t_797[k] = -7.0 * ii_601[k]
                   + f_0 * li_1021[k];

        t_798[k] = -7.0 * ii_602[k]
                   + f_0 * li_1022[k];

        t_799[k] = -7.0 * ii_603[k]
                   + f_0 * li_1023[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, ii_604, ii_605, ii_606, ii_607, \
                         ii_608, li_1024, li_1025, li_1026, li_1027, \
                         li_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = -7.0 * ii_604[k]
                   + f_0 * li_1024[k];

        t_801[k] = -7.0 * ii_605[k]
                   + f_0 * li_1025[k];

        t_802[k] = -7.0 * ii_606[k]
                   + f_0 * li_1026[k];

        t_803[k] = -7.0 * ii_607[k]
                   + f_0 * li_1027[k];

        t_804[k] = -7.0 * ii_608[k]
                   + f_0 * li_1028[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, ii_609, ii_610, ii_611, ii_612, \
                         ii_613, li_1029, li_1030, li_1031, li_1032, \
                         li_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = -7.0 * ii_609[k]
                   + f_0 * li_1029[k];

        t_806[k] = -7.0 * ii_610[k]
                   + f_0 * li_1030[k];

        t_807[k] = -7.0 * ii_611[k]
                   + f_0 * li_1031[k];

        t_808[k] = -7.0 * ii_612[k]
                   + f_0 * li_1032[k];

        t_809[k] = -7.0 * ii_613[k]
                   + f_0 * li_1033[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, t_814, ii_614, ii_615, ii_616, ii_617, \
                         ii_618, li_1034, li_1035, li_1036, li_1037, \
                         li_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = -7.0 * ii_614[k]
                   + f_0 * li_1034[k];

        t_811[k] = -7.0 * ii_615[k]
                   + f_0 * li_1035[k];

        t_812[k] = -6.0 * ii_616[k]
                   + f_0 * li_1036[k];

        t_813[k] = -6.0 * ii_617[k]
                   + f_0 * li_1037[k];

        t_814[k] = -6.0 * ii_618[k]
                   + f_0 * li_1038[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_1_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ii_619 = buffer.data(ii + 619);
    const auto *ii_620 = buffer.data(ii + 620);
    const auto *ii_621 = buffer.data(ii + 621);
    const auto *ii_622 = buffer.data(ii + 622);
    const auto *ii_623 = buffer.data(ii + 623);
    const auto *ii_624 = buffer.data(ii + 624);
    const auto *ii_625 = buffer.data(ii + 625);
    const auto *ii_626 = buffer.data(ii + 626);
    const auto *ii_627 = buffer.data(ii + 627);
    const auto *ii_628 = buffer.data(ii + 628);
    const auto *ii_629 = buffer.data(ii + 629);
    const auto *ii_630 = buffer.data(ii + 630);
    const auto *ii_631 = buffer.data(ii + 631);
    const auto *ii_632 = buffer.data(ii + 632);
    const auto *ii_633 = buffer.data(ii + 633);
    const auto *ii_634 = buffer.data(ii + 634);
    const auto *ii_635 = buffer.data(ii + 635);
    const auto *ii_636 = buffer.data(ii + 636);
    const auto *ii_637 = buffer.data(ii + 637);
    const auto *ii_638 = buffer.data(ii + 638);
    const auto *ii_639 = buffer.data(ii + 639);
    const auto *ii_640 = buffer.data(ii + 640);
    const auto *ii_641 = buffer.data(ii + 641);
    const auto *ii_642 = buffer.data(ii + 642);
    const auto *ii_643 = buffer.data(ii + 643);
    const auto *ii_644 = buffer.data(ii + 644);
    const auto *ii_645 = buffer.data(ii + 645);
    const auto *ii_646 = buffer.data(ii + 646);
    const auto *ii_647 = buffer.data(ii + 647);
    const auto *ii_648 = buffer.data(ii + 648);
    const auto *ii_649 = buffer.data(ii + 649);
    const auto *ii_650 = buffer.data(ii + 650);
    const auto *ii_651 = buffer.data(ii + 651);
    const auto *ii_652 = buffer.data(ii + 652);
    const auto *ii_653 = buffer.data(ii + 653);
    const auto *ii_654 = buffer.data(ii + 654);
    const auto *ii_655 = buffer.data(ii + 655);
    const auto *ii_656 = buffer.data(ii + 656);
    const auto *ii_657 = buffer.data(ii + 657);
    const auto *ii_658 = buffer.data(ii + 658);
    const auto *ii_659 = buffer.data(ii + 659);
    const auto *ii_660 = buffer.data(ii + 660);
    const auto *ii_661 = buffer.data(ii + 661);
    const auto *ii_662 = buffer.data(ii + 662);
    const auto *ii_663 = buffer.data(ii + 663);
    const auto *ii_664 = buffer.data(ii + 664);
    const auto *ii_665 = buffer.data(ii + 665);
    const auto *ii_666 = buffer.data(ii + 666);
    const auto *ii_667 = buffer.data(ii + 667);
    const auto *ii_668 = buffer.data(ii + 668);
    const auto *ii_669 = buffer.data(ii + 669);
    const auto *ii_670 = buffer.data(ii + 670);
    const auto *ii_671 = buffer.data(ii + 671);
    const auto *ii_672 = buffer.data(ii + 672);
    const auto *ii_673 = buffer.data(ii + 673);
    const auto *ii_674 = buffer.data(ii + 674);
    const auto *ii_675 = buffer.data(ii + 675);
    const auto *ii_676 = buffer.data(ii + 676);
    const auto *ii_677 = buffer.data(ii + 677);
    const auto *ii_678 = buffer.data(ii + 678);
    const auto *ii_679 = buffer.data(ii + 679);
    const auto *ii_680 = buffer.data(ii + 680);
    const auto *ii_681 = buffer.data(ii + 681);
    const auto *ii_682 = buffer.data(ii + 682);
    const auto *ii_683 = buffer.data(ii + 683);
    const auto *ii_684 = buffer.data(ii + 684);
    const auto *ii_685 = buffer.data(ii + 685);
    const auto *ii_686 = buffer.data(ii + 686);
    const auto *ii_687 = buffer.data(ii + 687);
    const auto *ii_688 = buffer.data(ii + 688);
    const auto *ii_689 = buffer.data(ii + 689);
    const auto *ii_690 = buffer.data(ii + 690);
    const auto *ii_691 = buffer.data(ii + 691);
    const auto *ii_692 = buffer.data(ii + 692);
    const auto *ii_693 = buffer.data(ii + 693);
    const auto *ii_694 = buffer.data(ii + 694);
    const auto *ii_695 = buffer.data(ii + 695);
    const auto *ii_696 = buffer.data(ii + 696);
    const auto *ii_697 = buffer.data(ii + 697);
    const auto *ii_698 = buffer.data(ii + 698);
    const auto *ii_699 = buffer.data(ii + 699);
    const auto *ii_700 = buffer.data(ii + 700);
    const auto *ii_701 = buffer.data(ii + 701);
    const auto *ii_702 = buffer.data(ii + 702);
    const auto *ii_703 = buffer.data(ii + 703);
    const auto *ii_704 = buffer.data(ii + 704);
    const auto *ii_705 = buffer.data(ii + 705);
    const auto *ii_706 = buffer.data(ii + 706);
    const auto *ii_707 = buffer.data(ii + 707);
    const auto *ii_708 = buffer.data(ii + 708);
    const auto *ii_709 = buffer.data(ii + 709);
    const auto *ii_710 = buffer.data(ii + 710);
    const auto *ii_711 = buffer.data(ii + 711);
    const auto *ii_712 = buffer.data(ii + 712);
    const auto *ii_713 = buffer.data(ii + 713);
    const auto *ii_714 = buffer.data(ii + 714);
    const auto *ii_715 = buffer.data(ii + 715);
    const auto *ii_716 = buffer.data(ii + 716);
    const auto *ii_717 = buffer.data(ii + 717);
    const auto *ii_718 = buffer.data(ii + 718);
    const auto *ii_719 = buffer.data(ii + 719);
    const auto *ii_720 = buffer.data(ii + 720);
    const auto *ii_721 = buffer.data(ii + 721);
    const auto *ii_722 = buffer.data(ii + 722);
    const auto *ii_723 = buffer.data(ii + 723);
    const auto *ii_724 = buffer.data(ii + 724);
    const auto *ii_725 = buffer.data(ii + 725);
    const auto *ii_726 = buffer.data(ii + 726);
    const auto *ii_727 = buffer.data(ii + 727);
    const auto *ii_728 = buffer.data(ii + 728);
    const auto *ii_729 = buffer.data(ii + 729);
    const auto *ii_730 = buffer.data(ii + 730);
    const auto *ii_731 = buffer.data(ii + 731);
    const auto *ii_732 = buffer.data(ii + 732);
    const auto *ii_733 = buffer.data(ii + 733);
    const auto *ii_734 = buffer.data(ii + 734);
    const auto *ii_735 = buffer.data(ii + 735);
    const auto *ii_736 = buffer.data(ii + 736);
    const auto *ii_737 = buffer.data(ii + 737);
    const auto *ii_738 = buffer.data(ii + 738);
    const auto *ii_739 = buffer.data(ii + 739);
    const auto *ii_740 = buffer.data(ii + 740);
    const auto *ii_741 = buffer.data(ii + 741);
    const auto *ii_742 = buffer.data(ii + 742);
    const auto *ii_743 = buffer.data(ii + 743);
    const auto *ii_744 = buffer.data(ii + 744);
    const auto *ii_745 = buffer.data(ii + 745);
    const auto *ii_746 = buffer.data(ii + 746);
    const auto *ii_747 = buffer.data(ii + 747);
    const auto *ii_748 = buffer.data(ii + 748);
    const auto *ii_749 = buffer.data(ii + 749);
    const auto *ii_750 = buffer.data(ii + 750);
    const auto *ii_751 = buffer.data(ii + 751);
    const auto *ii_752 = buffer.data(ii + 752);
    const auto *ii_753 = buffer.data(ii + 753);
    const auto *ii_754 = buffer.data(ii + 754);
    const auto *ii_755 = buffer.data(ii + 755);
    const auto *ii_756 = buffer.data(ii + 756);
    const auto *ii_757 = buffer.data(ii + 757);
    const auto *ii_758 = buffer.data(ii + 758);
    const auto *ii_759 = buffer.data(ii + 759);
    const auto *ii_760 = buffer.data(ii + 760);
    const auto *ii_761 = buffer.data(ii + 761);
    const auto *ii_762 = buffer.data(ii + 762);
    const auto *ii_763 = buffer.data(ii + 763);
    const auto *ii_764 = buffer.data(ii + 764);
    const auto *ii_765 = buffer.data(ii + 765);
    const auto *ii_766 = buffer.data(ii + 766);
    const auto *ii_767 = buffer.data(ii + 767);
    const auto *ii_768 = buffer.data(ii + 768);

    const auto *li_1039 = buffer.data(li + 1039);
    const auto *li_1040 = buffer.data(li + 1040);
    const auto *li_1041 = buffer.data(li + 1041);
    const auto *li_1042 = buffer.data(li + 1042);
    const auto *li_1043 = buffer.data(li + 1043);
    const auto *li_1044 = buffer.data(li + 1044);
    const auto *li_1045 = buffer.data(li + 1045);
    const auto *li_1046 = buffer.data(li + 1046);
    const auto *li_1047 = buffer.data(li + 1047);
    const auto *li_1048 = buffer.data(li + 1048);
    const auto *li_1049 = buffer.data(li + 1049);
    const auto *li_1050 = buffer.data(li + 1050);
    const auto *li_1051 = buffer.data(li + 1051);
    const auto *li_1052 = buffer.data(li + 1052);
    const auto *li_1053 = buffer.data(li + 1053);
    const auto *li_1054 = buffer.data(li + 1054);
    const auto *li_1055 = buffer.data(li + 1055);
    const auto *li_1056 = buffer.data(li + 1056);
    const auto *li_1057 = buffer.data(li + 1057);
    const auto *li_1058 = buffer.data(li + 1058);
    const auto *li_1059 = buffer.data(li + 1059);
    const auto *li_1060 = buffer.data(li + 1060);
    const auto *li_1061 = buffer.data(li + 1061);
    const auto *li_1062 = buffer.data(li + 1062);
    const auto *li_1063 = buffer.data(li + 1063);
    const auto *li_1064 = buffer.data(li + 1064);
    const auto *li_1065 = buffer.data(li + 1065);
    const auto *li_1066 = buffer.data(li + 1066);
    const auto *li_1067 = buffer.data(li + 1067);
    const auto *li_1068 = buffer.data(li + 1068);
    const auto *li_1069 = buffer.data(li + 1069);
    const auto *li_1070 = buffer.data(li + 1070);
    const auto *li_1071 = buffer.data(li + 1071);
    const auto *li_1072 = buffer.data(li + 1072);
    const auto *li_1073 = buffer.data(li + 1073);
    const auto *li_1074 = buffer.data(li + 1074);
    const auto *li_1075 = buffer.data(li + 1075);
    const auto *li_1076 = buffer.data(li + 1076);
    const auto *li_1077 = buffer.data(li + 1077);
    const auto *li_1078 = buffer.data(li + 1078);
    const auto *li_1079 = buffer.data(li + 1079);
    const auto *li_1080 = buffer.data(li + 1080);
    const auto *li_1081 = buffer.data(li + 1081);
    const auto *li_1082 = buffer.data(li + 1082);
    const auto *li_1083 = buffer.data(li + 1083);
    const auto *li_1084 = buffer.data(li + 1084);
    const auto *li_1085 = buffer.data(li + 1085);
    const auto *li_1086 = buffer.data(li + 1086);
    const auto *li_1087 = buffer.data(li + 1087);
    const auto *li_1088 = buffer.data(li + 1088);
    const auto *li_1089 = buffer.data(li + 1089);
    const auto *li_1090 = buffer.data(li + 1090);
    const auto *li_1091 = buffer.data(li + 1091);
    const auto *li_1092 = buffer.data(li + 1092);
    const auto *li_1093 = buffer.data(li + 1093);
    const auto *li_1094 = buffer.data(li + 1094);
    const auto *li_1095 = buffer.data(li + 1095);
    const auto *li_1096 = buffer.data(li + 1096);
    const auto *li_1097 = buffer.data(li + 1097);
    const auto *li_1098 = buffer.data(li + 1098);
    const auto *li_1099 = buffer.data(li + 1099);
    const auto *li_1100 = buffer.data(li + 1100);
    const auto *li_1101 = buffer.data(li + 1101);
    const auto *li_1102 = buffer.data(li + 1102);
    const auto *li_1103 = buffer.data(li + 1103);
    const auto *li_1104 = buffer.data(li + 1104);
    const auto *li_1105 = buffer.data(li + 1105);
    const auto *li_1106 = buffer.data(li + 1106);
    const auto *li_1107 = buffer.data(li + 1107);
    const auto *li_1108 = buffer.data(li + 1108);
    const auto *li_1109 = buffer.data(li + 1109);
    const auto *li_1110 = buffer.data(li + 1110);
    const auto *li_1111 = buffer.data(li + 1111);
    const auto *li_1112 = buffer.data(li + 1112);
    const auto *li_1113 = buffer.data(li + 1113);
    const auto *li_1114 = buffer.data(li + 1114);
    const auto *li_1115 = buffer.data(li + 1115);
    const auto *li_1116 = buffer.data(li + 1116);
    const auto *li_1117 = buffer.data(li + 1117);
    const auto *li_1118 = buffer.data(li + 1118);
    const auto *li_1119 = buffer.data(li + 1119);
    const auto *li_1120 = buffer.data(li + 1120);
    const auto *li_1121 = buffer.data(li + 1121);
    const auto *li_1122 = buffer.data(li + 1122);
    const auto *li_1123 = buffer.data(li + 1123);
    const auto *li_1124 = buffer.data(li + 1124);
    const auto *li_1125 = buffer.data(li + 1125);
    const auto *li_1126 = buffer.data(li + 1126);
    const auto *li_1127 = buffer.data(li + 1127);
    const auto *li_1128 = buffer.data(li + 1128);
    const auto *li_1129 = buffer.data(li + 1129);
    const auto *li_1130 = buffer.data(li + 1130);
    const auto *li_1131 = buffer.data(li + 1131);
    const auto *li_1132 = buffer.data(li + 1132);
    const auto *li_1133 = buffer.data(li + 1133);
    const auto *li_1134 = buffer.data(li + 1134);
    const auto *li_1135 = buffer.data(li + 1135);
    const auto *li_1136 = buffer.data(li + 1136);
    const auto *li_1137 = buffer.data(li + 1137);
    const auto *li_1138 = buffer.data(li + 1138);
    const auto *li_1139 = buffer.data(li + 1139);
    const auto *li_1140 = buffer.data(li + 1140);
    const auto *li_1141 = buffer.data(li + 1141);
    const auto *li_1142 = buffer.data(li + 1142);
    const auto *li_1143 = buffer.data(li + 1143);
    const auto *li_1144 = buffer.data(li + 1144);
    const auto *li_1145 = buffer.data(li + 1145);
    const auto *li_1146 = buffer.data(li + 1146);
    const auto *li_1147 = buffer.data(li + 1147);
    const auto *li_1148 = buffer.data(li + 1148);
    const auto *li_1149 = buffer.data(li + 1149);
    const auto *li_1150 = buffer.data(li + 1150);
    const auto *li_1151 = buffer.data(li + 1151);
    const auto *li_1152 = buffer.data(li + 1152);
    const auto *li_1153 = buffer.data(li + 1153);
    const auto *li_1154 = buffer.data(li + 1154);
    const auto *li_1155 = buffer.data(li + 1155);
    const auto *li_1156 = buffer.data(li + 1156);
    const auto *li_1157 = buffer.data(li + 1157);
    const auto *li_1158 = buffer.data(li + 1158);
    const auto *li_1159 = buffer.data(li + 1159);
    const auto *li_1160 = buffer.data(li + 1160);
    const auto *li_1161 = buffer.data(li + 1161);
    const auto *li_1162 = buffer.data(li + 1162);
    const auto *li_1163 = buffer.data(li + 1163);
    const auto *li_1164 = buffer.data(li + 1164);
    const auto *li_1165 = buffer.data(li + 1165);
    const auto *li_1166 = buffer.data(li + 1166);
    const auto *li_1167 = buffer.data(li + 1167);
    const auto *li_1168 = buffer.data(li + 1168);
    const auto *li_1169 = buffer.data(li + 1169);
    const auto *li_1170 = buffer.data(li + 1170);
    const auto *li_1171 = buffer.data(li + 1171);
    const auto *li_1172 = buffer.data(li + 1172);
    const auto *li_1173 = buffer.data(li + 1173);
    const auto *li_1174 = buffer.data(li + 1174);
    const auto *li_1175 = buffer.data(li + 1175);
    const auto *li_1176 = buffer.data(li + 1176);
    const auto *li_1177 = buffer.data(li + 1177);
    const auto *li_1178 = buffer.data(li + 1178);
    const auto *li_1179 = buffer.data(li + 1179);
    const auto *li_1180 = buffer.data(li + 1180);
    const auto *li_1181 = buffer.data(li + 1181);
    const auto *li_1182 = buffer.data(li + 1182);
    const auto *li_1183 = buffer.data(li + 1183);
    const auto *li_1184 = buffer.data(li + 1184);
    const auto *li_1185 = buffer.data(li + 1185);
    const auto *li_1186 = buffer.data(li + 1186);
    const auto *li_1187 = buffer.data(li + 1187);
    const auto *li_1188 = buffer.data(li + 1188);

#pragma omp simd aligned(t_815, t_816, t_817, t_818, t_819, ii_619, ii_620, ii_621, ii_622, \
                         ii_623, li_1039, li_1040, li_1041, li_1042, \
                         li_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = -6.0 * ii_619[k]
                   + f_0 * li_1039[k];

        t_816[k] = -6.0 * ii_620[k]
                   + f_0 * li_1040[k];

        t_817[k] = -6.0 * ii_621[k]
                   + f_0 * li_1041[k];

        t_818[k] = -6.0 * ii_622[k]
                   + f_0 * li_1042[k];

        t_819[k] = -6.0 * ii_623[k]
                   + f_0 * li_1043[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, t_824, ii_624, ii_625, ii_626, ii_627, \
                         ii_628, li_1044, li_1045, li_1046, li_1047, \
                         li_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = -6.0 * ii_624[k]
                   + f_0 * li_1044[k];

        t_821[k] = -6.0 * ii_625[k]
                   + f_0 * li_1045[k];

        t_822[k] = -6.0 * ii_626[k]
                   + f_0 * li_1046[k];

        t_823[k] = -6.0 * ii_627[k]
                   + f_0 * li_1047[k];

        t_824[k] = -6.0 * ii_628[k]
                   + f_0 * li_1048[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, ii_629, ii_630, ii_631, ii_632, \
                         ii_633, li_1049, li_1050, li_1051, li_1052, \
                         li_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = -6.0 * ii_629[k]
                   + f_0 * li_1049[k];

        t_826[k] = -6.0 * ii_630[k]
                   + f_0 * li_1050[k];

        t_827[k] = -6.0 * ii_631[k]
                   + f_0 * li_1051[k];

        t_828[k] = -6.0 * ii_632[k]
                   + f_0 * li_1052[k];

        t_829[k] = -6.0 * ii_633[k]
                   + f_0 * li_1053[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, ii_634, ii_635, ii_636, ii_637, \
                         ii_638, li_1054, li_1055, li_1056, li_1057, \
                         li_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = -6.0 * ii_634[k]
                   + f_0 * li_1054[k];

        t_831[k] = -6.0 * ii_635[k]
                   + f_0 * li_1055[k];

        t_832[k] = -6.0 * ii_636[k]
                   + f_0 * li_1056[k];

        t_833[k] = -6.0 * ii_637[k]
                   + f_0 * li_1057[k];

        t_834[k] = -6.0 * ii_638[k]
                   + f_0 * li_1058[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, t_839, ii_639, ii_640, ii_641, ii_642, \
                         ii_643, li_1059, li_1060, li_1061, li_1062, \
                         li_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = -6.0 * ii_639[k]
                   + f_0 * li_1059[k];

        t_836[k] = -6.0 * ii_640[k]
                   + f_0 * li_1060[k];

        t_837[k] = -6.0 * ii_641[k]
                   + f_0 * li_1061[k];

        t_838[k] = -6.0 * ii_642[k]
                   + f_0 * li_1062[k];

        t_839[k] = -6.0 * ii_643[k]
                   + f_0 * li_1063[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, ii_644, ii_645, ii_646, ii_647, \
                         ii_648, li_1064, li_1065, li_1066, li_1067, \
                         li_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = -5.0 * ii_644[k]
                   + f_0 * li_1064[k];

        t_841[k] = -5.0 * ii_645[k]
                   + f_0 * li_1065[k];

        t_842[k] = -5.0 * ii_646[k]
                   + f_0 * li_1066[k];

        t_843[k] = -5.0 * ii_647[k]
                   + f_0 * li_1067[k];

        t_844[k] = -5.0 * ii_648[k]
                   + f_0 * li_1068[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, ii_649, ii_650, ii_651, ii_652, \
                         ii_653, li_1069, li_1070, li_1071, li_1072, \
                         li_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = -5.0 * ii_649[k]
                   + f_0 * li_1069[k];

        t_846[k] = -5.0 * ii_650[k]
                   + f_0 * li_1070[k];

        t_847[k] = -5.0 * ii_651[k]
                   + f_0 * li_1071[k];

        t_848[k] = -5.0 * ii_652[k]
                   + f_0 * li_1072[k];

        t_849[k] = -5.0 * ii_653[k]
                   + f_0 * li_1073[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, ii_654, ii_655, ii_656, ii_657, \
                         ii_658, li_1074, li_1075, li_1076, li_1077, \
                         li_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = -5.0 * ii_654[k]
                   + f_0 * li_1074[k];

        t_851[k] = -5.0 * ii_655[k]
                   + f_0 * li_1075[k];

        t_852[k] = -5.0 * ii_656[k]
                   + f_0 * li_1076[k];

        t_853[k] = -5.0 * ii_657[k]
                   + f_0 * li_1077[k];

        t_854[k] = -5.0 * ii_658[k]
                   + f_0 * li_1078[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, ii_659, ii_660, ii_661, ii_662, \
                         ii_663, li_1079, li_1080, li_1081, li_1082, \
                         li_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = -5.0 * ii_659[k]
                   + f_0 * li_1079[k];

        t_856[k] = -5.0 * ii_660[k]
                   + f_0 * li_1080[k];

        t_857[k] = -5.0 * ii_661[k]
                   + f_0 * li_1081[k];

        t_858[k] = -5.0 * ii_662[k]
                   + f_0 * li_1082[k];

        t_859[k] = -5.0 * ii_663[k]
                   + f_0 * li_1083[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, ii_664, ii_665, ii_666, ii_667, \
                         ii_668, li_1084, li_1085, li_1086, li_1087, \
                         li_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = -5.0 * ii_664[k]
                   + f_0 * li_1084[k];

        t_861[k] = -5.0 * ii_665[k]
                   + f_0 * li_1085[k];

        t_862[k] = -5.0 * ii_666[k]
                   + f_0 * li_1086[k];

        t_863[k] = -5.0 * ii_667[k]
                   + f_0 * li_1087[k];

        t_864[k] = -5.0 * ii_668[k]
                   + f_0 * li_1088[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, t_869, ii_669, ii_670, ii_671, ii_672, \
                         ii_673, li_1089, li_1090, li_1091, li_1092, \
                         li_1093 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = -5.0 * ii_669[k]
                   + f_0 * li_1089[k];

        t_866[k] = -5.0 * ii_670[k]
                   + f_0 * li_1090[k];

        t_867[k] = -5.0 * ii_671[k]
                   + f_0 * li_1091[k];

        t_868[k] = -4.0 * ii_672[k]
                   + f_0 * li_1092[k];

        t_869[k] = -4.0 * ii_673[k]
                   + f_0 * li_1093[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, ii_674, ii_675, ii_676, ii_677, \
                         ii_678, li_1094, li_1095, li_1096, li_1097, \
                         li_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = -4.0 * ii_674[k]
                   + f_0 * li_1094[k];

        t_871[k] = -4.0 * ii_675[k]
                   + f_0 * li_1095[k];

        t_872[k] = -4.0 * ii_676[k]
                   + f_0 * li_1096[k];

        t_873[k] = -4.0 * ii_677[k]
                   + f_0 * li_1097[k];

        t_874[k] = -4.0 * ii_678[k]
                   + f_0 * li_1098[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, t_879, ii_679, ii_680, ii_681, ii_682, \
                         ii_683, li_1099, li_1100, li_1101, li_1102, \
                         li_1103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = -4.0 * ii_679[k]
                   + f_0 * li_1099[k];

        t_876[k] = -4.0 * ii_680[k]
                   + f_0 * li_1100[k];

        t_877[k] = -4.0 * ii_681[k]
                   + f_0 * li_1101[k];

        t_878[k] = -4.0 * ii_682[k]
                   + f_0 * li_1102[k];

        t_879[k] = -4.0 * ii_683[k]
                   + f_0 * li_1103[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, t_884, ii_684, ii_685, ii_686, ii_687, \
                         ii_688, li_1104, li_1105, li_1106, li_1107, \
                         li_1108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = -4.0 * ii_684[k]
                   + f_0 * li_1104[k];

        t_881[k] = -4.0 * ii_685[k]
                   + f_0 * li_1105[k];

        t_882[k] = -4.0 * ii_686[k]
                   + f_0 * li_1106[k];

        t_883[k] = -4.0 * ii_687[k]
                   + f_0 * li_1107[k];

        t_884[k] = -4.0 * ii_688[k]
                   + f_0 * li_1108[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, ii_689, ii_690, ii_691, ii_692, \
                         ii_693, li_1109, li_1110, li_1111, li_1112, \
                         li_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = -4.0 * ii_689[k]
                   + f_0 * li_1109[k];

        t_886[k] = -4.0 * ii_690[k]
                   + f_0 * li_1110[k];

        t_887[k] = -4.0 * ii_691[k]
                   + f_0 * li_1111[k];

        t_888[k] = -4.0 * ii_692[k]
                   + f_0 * li_1112[k];

        t_889[k] = -4.0 * ii_693[k]
                   + f_0 * li_1113[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, ii_694, ii_695, ii_696, ii_697, \
                         ii_698, li_1114, li_1115, li_1116, li_1117, \
                         li_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = -4.0 * ii_694[k]
                   + f_0 * li_1114[k];

        t_891[k] = -4.0 * ii_695[k]
                   + f_0 * li_1115[k];

        t_892[k] = -4.0 * ii_696[k]
                   + f_0 * li_1116[k];

        t_893[k] = -4.0 * ii_697[k]
                   + f_0 * li_1117[k];

        t_894[k] = -4.0 * ii_698[k]
                   + f_0 * li_1118[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, t_899, ii_699, ii_700, ii_701, ii_702, \
                         ii_703, li_1119, li_1120, li_1121, li_1122, \
                         li_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = -4.0 * ii_699[k]
                   + f_0 * li_1119[k];

        t_896[k] = -3.0 * ii_700[k]
                   + f_0 * li_1120[k];

        t_897[k] = -3.0 * ii_701[k]
                   + f_0 * li_1121[k];

        t_898[k] = -3.0 * ii_702[k]
                   + f_0 * li_1122[k];

        t_899[k] = -3.0 * ii_703[k]
                   + f_0 * li_1123[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, t_904, ii_704, ii_705, ii_706, ii_707, \
                         ii_708, li_1124, li_1125, li_1126, li_1127, \
                         li_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = -3.0 * ii_704[k]
                   + f_0 * li_1124[k];

        t_901[k] = -3.0 * ii_705[k]
                   + f_0 * li_1125[k];

        t_902[k] = -3.0 * ii_706[k]
                   + f_0 * li_1126[k];

        t_903[k] = -3.0 * ii_707[k]
                   + f_0 * li_1127[k];

        t_904[k] = -3.0 * ii_708[k]
                   + f_0 * li_1128[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, t_909, ii_709, ii_710, ii_711, ii_712, \
                         ii_713, li_1129, li_1130, li_1131, li_1132, \
                         li_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = -3.0 * ii_709[k]
                   + f_0 * li_1129[k];

        t_906[k] = -3.0 * ii_710[k]
                   + f_0 * li_1130[k];

        t_907[k] = -3.0 * ii_711[k]
                   + f_0 * li_1131[k];

        t_908[k] = -3.0 * ii_712[k]
                   + f_0 * li_1132[k];

        t_909[k] = -3.0 * ii_713[k]
                   + f_0 * li_1133[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, t_914, ii_714, ii_715, ii_716, ii_717, \
                         ii_718, li_1134, li_1135, li_1136, li_1137, \
                         li_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = -3.0 * ii_714[k]
                   + f_0 * li_1134[k];

        t_911[k] = -3.0 * ii_715[k]
                   + f_0 * li_1135[k];

        t_912[k] = -3.0 * ii_716[k]
                   + f_0 * li_1136[k];

        t_913[k] = -3.0 * ii_717[k]
                   + f_0 * li_1137[k];

        t_914[k] = -3.0 * ii_718[k]
                   + f_0 * li_1138[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, t_919, ii_719, ii_720, ii_721, ii_722, \
                         ii_723, li_1139, li_1140, li_1141, li_1142, \
                         li_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = -3.0 * ii_719[k]
                   + f_0 * li_1139[k];

        t_916[k] = -3.0 * ii_720[k]
                   + f_0 * li_1140[k];

        t_917[k] = -3.0 * ii_721[k]
                   + f_0 * li_1141[k];

        t_918[k] = -3.0 * ii_722[k]
                   + f_0 * li_1142[k];

        t_919[k] = -3.0 * ii_723[k]
                   + f_0 * li_1143[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, t_924, ii_724, ii_725, ii_726, ii_727, \
                         ii_728, li_1144, li_1145, li_1146, li_1147, \
                         li_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = -3.0 * ii_724[k]
                   + f_0 * li_1144[k];

        t_921[k] = -3.0 * ii_725[k]
                   + f_0 * li_1145[k];

        t_922[k] = -3.0 * ii_726[k]
                   + f_0 * li_1146[k];

        t_923[k] = -3.0 * ii_727[k]
                   + f_0 * li_1147[k];

        t_924[k] = -2.0 * ii_728[k]
                   + f_0 * li_1148[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, t_929, ii_729, ii_730, ii_731, ii_732, \
                         ii_733, li_1149, li_1150, li_1151, li_1152, \
                         li_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = -2.0 * ii_729[k]
                   + f_0 * li_1149[k];

        t_926[k] = -2.0 * ii_730[k]
                   + f_0 * li_1150[k];

        t_927[k] = -2.0 * ii_731[k]
                   + f_0 * li_1151[k];

        t_928[k] = -2.0 * ii_732[k]
                   + f_0 * li_1152[k];

        t_929[k] = -2.0 * ii_733[k]
                   + f_0 * li_1153[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, ii_734, ii_735, ii_736, ii_737, \
                         ii_738, li_1154, li_1155, li_1156, li_1157, \
                         li_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = -2.0 * ii_734[k]
                   + f_0 * li_1154[k];

        t_931[k] = -2.0 * ii_735[k]
                   + f_0 * li_1155[k];

        t_932[k] = -2.0 * ii_736[k]
                   + f_0 * li_1156[k];

        t_933[k] = -2.0 * ii_737[k]
                   + f_0 * li_1157[k];

        t_934[k] = -2.0 * ii_738[k]
                   + f_0 * li_1158[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, t_939, ii_739, ii_740, ii_741, ii_742, \
                         ii_743, li_1159, li_1160, li_1161, li_1162, \
                         li_1163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = -2.0 * ii_739[k]
                   + f_0 * li_1159[k];

        t_936[k] = -2.0 * ii_740[k]
                   + f_0 * li_1160[k];

        t_937[k] = -2.0 * ii_741[k]
                   + f_0 * li_1161[k];

        t_938[k] = -2.0 * ii_742[k]
                   + f_0 * li_1162[k];

        t_939[k] = -2.0 * ii_743[k]
                   + f_0 * li_1163[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, t_944, ii_744, ii_745, ii_746, ii_747, \
                         ii_748, li_1164, li_1165, li_1166, li_1167, \
                         li_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = -2.0 * ii_744[k]
                   + f_0 * li_1164[k];

        t_941[k] = -2.0 * ii_745[k]
                   + f_0 * li_1165[k];

        t_942[k] = -2.0 * ii_746[k]
                   + f_0 * li_1166[k];

        t_943[k] = -2.0 * ii_747[k]
                   + f_0 * li_1167[k];

        t_944[k] = -2.0 * ii_748[k]
                   + f_0 * li_1168[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, ii_749, ii_750, ii_751, ii_752, \
                         ii_753, li_1169, li_1170, li_1171, li_1172, \
                         li_1173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = -2.0 * ii_749[k]
                   + f_0 * li_1169[k];

        t_946[k] = -2.0 * ii_750[k]
                   + f_0 * li_1170[k];

        t_947[k] = -2.0 * ii_751[k]
                   + f_0 * li_1171[k];

        t_948[k] = -2.0 * ii_752[k]
                   + f_0 * li_1172[k];

        t_949[k] = -2.0 * ii_753[k]
                   + f_0 * li_1173[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, t_954, ii_754, ii_755, ii_756, ii_757, \
                         ii_758, li_1174, li_1175, li_1176, li_1177, \
                         li_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = -2.0 * ii_754[k]
                   + f_0 * li_1174[k];

        t_951[k] = -2.0 * ii_755[k]
                   + f_0 * li_1175[k];

        t_952[k] = -ii_756[k]
                   + f_0 * li_1176[k];

        t_953[k] = -ii_757[k]
                   + f_0 * li_1177[k];

        t_954[k] = -ii_758[k]
                   + f_0 * li_1178[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, t_959, ii_759, ii_760, ii_761, ii_762, \
                         ii_763, li_1179, li_1180, li_1181, li_1182, \
                         li_1183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = -ii_759[k]
                   + f_0 * li_1179[k];

        t_956[k] = -ii_760[k]
                   + f_0 * li_1180[k];

        t_957[k] = -ii_761[k]
                   + f_0 * li_1181[k];

        t_958[k] = -ii_762[k]
                   + f_0 * li_1182[k];

        t_959[k] = -ii_763[k]
                   + f_0 * li_1183[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, t_963, t_964, ii_764, ii_765, ii_766, ii_767, \
                         ii_768, li_1184, li_1185, li_1186, li_1187, \
                         li_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = -ii_764[k]
                   + f_0 * li_1184[k];

        t_961[k] = -ii_765[k]
                   + f_0 * li_1185[k];

        t_962[k] = -ii_766[k]
                   + f_0 * li_1186[k];

        t_963[k] = -ii_767[k]
                   + f_0 * li_1187[k];

        t_964[k] = -ii_768[k]
                   + f_0 * li_1188[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_1_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ii_769 = buffer.data(ii + 769);
    const auto *ii_770 = buffer.data(ii + 770);
    const auto *ii_771 = buffer.data(ii + 771);
    const auto *ii_772 = buffer.data(ii + 772);
    const auto *ii_773 = buffer.data(ii + 773);
    const auto *ii_774 = buffer.data(ii + 774);
    const auto *ii_775 = buffer.data(ii + 775);
    const auto *ii_776 = buffer.data(ii + 776);
    const auto *ii_777 = buffer.data(ii + 777);
    const auto *ii_778 = buffer.data(ii + 778);
    const auto *ii_779 = buffer.data(ii + 779);
    const auto *ii_780 = buffer.data(ii + 780);
    const auto *ii_781 = buffer.data(ii + 781);
    const auto *ii_782 = buffer.data(ii + 782);
    const auto *ii_783 = buffer.data(ii + 783);

    const auto *li_1189 = buffer.data(li + 1189);
    const auto *li_1190 = buffer.data(li + 1190);
    const auto *li_1191 = buffer.data(li + 1191);
    const auto *li_1192 = buffer.data(li + 1192);
    const auto *li_1193 = buffer.data(li + 1193);
    const auto *li_1194 = buffer.data(li + 1194);
    const auto *li_1195 = buffer.data(li + 1195);
    const auto *li_1196 = buffer.data(li + 1196);
    const auto *li_1197 = buffer.data(li + 1197);
    const auto *li_1198 = buffer.data(li + 1198);
    const auto *li_1199 = buffer.data(li + 1199);
    const auto *li_1200 = buffer.data(li + 1200);
    const auto *li_1201 = buffer.data(li + 1201);
    const auto *li_1202 = buffer.data(li + 1202);
    const auto *li_1203 = buffer.data(li + 1203);
    const auto *li_1204 = buffer.data(li + 1204);
    const auto *li_1205 = buffer.data(li + 1205);
    const auto *li_1206 = buffer.data(li + 1206);
    const auto *li_1207 = buffer.data(li + 1207);
    const auto *li_1208 = buffer.data(li + 1208);
    const auto *li_1209 = buffer.data(li + 1209);
    const auto *li_1210 = buffer.data(li + 1210);
    const auto *li_1211 = buffer.data(li + 1211);
    const auto *li_1212 = buffer.data(li + 1212);
    const auto *li_1213 = buffer.data(li + 1213);
    const auto *li_1214 = buffer.data(li + 1214);
    const auto *li_1215 = buffer.data(li + 1215);
    const auto *li_1216 = buffer.data(li + 1216);
    const auto *li_1217 = buffer.data(li + 1217);
    const auto *li_1218 = buffer.data(li + 1218);
    const auto *li_1219 = buffer.data(li + 1219);
    const auto *li_1220 = buffer.data(li + 1220);
    const auto *li_1221 = buffer.data(li + 1221);
    const auto *li_1222 = buffer.data(li + 1222);
    const auto *li_1223 = buffer.data(li + 1223);
    const auto *li_1224 = buffer.data(li + 1224);
    const auto *li_1225 = buffer.data(li + 1225);
    const auto *li_1226 = buffer.data(li + 1226);
    const auto *li_1227 = buffer.data(li + 1227);
    const auto *li_1228 = buffer.data(li + 1228);
    const auto *li_1229 = buffer.data(li + 1229);
    const auto *li_1230 = buffer.data(li + 1230);
    const auto *li_1231 = buffer.data(li + 1231);

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, ii_769, ii_770, ii_771, ii_772, \
                         ii_773, li_1189, li_1190, li_1191, li_1192, \
                         li_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = -ii_769[k]
                   + f_0 * li_1189[k];

        t_966[k] = -ii_770[k]
                   + f_0 * li_1190[k];

        t_967[k] = -ii_771[k]
                   + f_0 * li_1191[k];

        t_968[k] = -ii_772[k]
                   + f_0 * li_1192[k];

        t_969[k] = -ii_773[k]
                   + f_0 * li_1193[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, ii_774, ii_775, ii_776, ii_777, \
                         ii_778, li_1194, li_1195, li_1196, li_1197, \
                         li_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = -ii_774[k]
                   + f_0 * li_1194[k];

        t_971[k] = -ii_775[k]
                   + f_0 * li_1195[k];

        t_972[k] = -ii_776[k]
                   + f_0 * li_1196[k];

        t_973[k] = -ii_777[k]
                   + f_0 * li_1197[k];

        t_974[k] = -ii_778[k]
                   + f_0 * li_1198[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, ii_779, ii_780, ii_781, ii_782, \
                         ii_783, li_1199, li_1200, li_1201, li_1202, \
                         li_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = -ii_779[k]
                   + f_0 * li_1199[k];

        t_976[k] = -ii_780[k]
                   + f_0 * li_1200[k];

        t_977[k] = -ii_781[k]
                   + f_0 * li_1201[k];

        t_978[k] = -ii_782[k]
                   + f_0 * li_1202[k];

        t_979[k] = -ii_783[k]
                   + f_0 * li_1203[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, t_984, t_985, t_986, t_987, li_1204, \
                         li_1205, li_1206, li_1207, li_1208, li_1209, li_1210, \
                         li_1211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = f_0 * li_1204[k];

        t_981[k] = f_0 * li_1205[k];

        t_982[k] = f_0 * li_1206[k];

        t_983[k] = f_0 * li_1207[k];

        t_984[k] = f_0 * li_1208[k];

        t_985[k] = f_0 * li_1209[k];

        t_986[k] = f_0 * li_1210[k];

        t_987[k] = f_0 * li_1211[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, t_992, t_993, t_994, t_995, li_1212, \
                         li_1213, li_1214, li_1215, li_1216, li_1217, li_1218, \
                         li_1219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_0 * li_1212[k];

        t_989[k] = f_0 * li_1213[k];

        t_990[k] = f_0 * li_1214[k];

        t_991[k] = f_0 * li_1215[k];

        t_992[k] = f_0 * li_1216[k];

        t_993[k] = f_0 * li_1217[k];

        t_994[k] = f_0 * li_1218[k];

        t_995[k] = f_0 * li_1219[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, t_1000, t_1001, t_1002, t_1003, li_1220, \
                         li_1221, li_1222, li_1223, li_1224, li_1225, li_1226, \
                         li_1227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_0 * li_1220[k];

        t_997[k] = f_0 * li_1221[k];

        t_998[k] = f_0 * li_1222[k];

        t_999[k] = f_0 * li_1223[k];

        t_1000[k] = f_0 * li_1224[k];

        t_1001[k] = f_0 * li_1225[k];

        t_1002[k] = f_0 * li_1226[k];

        t_1003[k] = f_0 * li_1227[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, li_1228, li_1229, li_1230, \
                         li_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_0 * li_1228[k];

        t_1005[k] = f_0 * li_1229[k];

        t_1006[k] = f_0 * li_1230[k];

        t_1007[k] = f_0 * li_1231[k];
    }
}

auto
compute_prim_geom_10_ki_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ii, const size_t li,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ki_electron_repulsion_1_piece0(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_1_piece1(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_1_piece2(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_1_piece3(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_1_piece4(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_1_piece5(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_1_piece6(buffer, target, ii, li, ncols, alpha);
}

static auto
compute_prim_geom_10_ki_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_0 = buffer.data(ii + 0);
    const auto *ii_1 = buffer.data(ii + 1);
    const auto *ii_2 = buffer.data(ii + 2);
    const auto *ii_3 = buffer.data(ii + 3);
    const auto *ii_4 = buffer.data(ii + 4);
    const auto *ii_5 = buffer.data(ii + 5);
    const auto *ii_6 = buffer.data(ii + 6);
    const auto *ii_7 = buffer.data(ii + 7);
    const auto *ii_8 = buffer.data(ii + 8);
    const auto *ii_9 = buffer.data(ii + 9);
    const auto *ii_10 = buffer.data(ii + 10);
    const auto *ii_11 = buffer.data(ii + 11);
    const auto *ii_12 = buffer.data(ii + 12);
    const auto *ii_13 = buffer.data(ii + 13);
    const auto *ii_14 = buffer.data(ii + 14);
    const auto *ii_15 = buffer.data(ii + 15);
    const auto *ii_16 = buffer.data(ii + 16);
    const auto *ii_17 = buffer.data(ii + 17);
    const auto *ii_18 = buffer.data(ii + 18);
    const auto *ii_19 = buffer.data(ii + 19);
    const auto *ii_20 = buffer.data(ii + 20);
    const auto *ii_21 = buffer.data(ii + 21);
    const auto *ii_22 = buffer.data(ii + 22);
    const auto *ii_23 = buffer.data(ii + 23);
    const auto *ii_24 = buffer.data(ii + 24);
    const auto *ii_25 = buffer.data(ii + 25);
    const auto *ii_26 = buffer.data(ii + 26);
    const auto *ii_27 = buffer.data(ii + 27);
    const auto *ii_28 = buffer.data(ii + 28);
    const auto *ii_29 = buffer.data(ii + 29);
    const auto *ii_30 = buffer.data(ii + 30);
    const auto *ii_31 = buffer.data(ii + 31);
    const auto *ii_32 = buffer.data(ii + 32);
    const auto *ii_33 = buffer.data(ii + 33);
    const auto *ii_34 = buffer.data(ii + 34);
    const auto *ii_35 = buffer.data(ii + 35);
    const auto *ii_36 = buffer.data(ii + 36);
    const auto *ii_37 = buffer.data(ii + 37);
    const auto *ii_38 = buffer.data(ii + 38);
    const auto *ii_39 = buffer.data(ii + 39);
    const auto *ii_40 = buffer.data(ii + 40);
    const auto *ii_41 = buffer.data(ii + 41);
    const auto *ii_42 = buffer.data(ii + 42);
    const auto *ii_43 = buffer.data(ii + 43);
    const auto *ii_44 = buffer.data(ii + 44);
    const auto *ii_45 = buffer.data(ii + 45);
    const auto *ii_46 = buffer.data(ii + 46);
    const auto *ii_47 = buffer.data(ii + 47);
    const auto *ii_48 = buffer.data(ii + 48);
    const auto *ii_49 = buffer.data(ii + 49);
    const auto *ii_50 = buffer.data(ii + 50);
    const auto *ii_51 = buffer.data(ii + 51);
    const auto *ii_52 = buffer.data(ii + 52);
    const auto *ii_53 = buffer.data(ii + 53);
    const auto *ii_54 = buffer.data(ii + 54);
    const auto *ii_55 = buffer.data(ii + 55);
    const auto *ii_56 = buffer.data(ii + 56);
    const auto *ii_57 = buffer.data(ii + 57);
    const auto *ii_58 = buffer.data(ii + 58);
    const auto *ii_59 = buffer.data(ii + 59);
    const auto *ii_60 = buffer.data(ii + 60);
    const auto *ii_61 = buffer.data(ii + 61);
    const auto *ii_62 = buffer.data(ii + 62);
    const auto *ii_63 = buffer.data(ii + 63);
    const auto *ii_64 = buffer.data(ii + 64);
    const auto *ii_65 = buffer.data(ii + 65);
    const auto *ii_66 = buffer.data(ii + 66);
    const auto *ii_67 = buffer.data(ii + 67);
    const auto *ii_68 = buffer.data(ii + 68);
    const auto *ii_69 = buffer.data(ii + 69);
    const auto *ii_70 = buffer.data(ii + 70);
    const auto *ii_71 = buffer.data(ii + 71);
    const auto *ii_72 = buffer.data(ii + 72);
    const auto *ii_73 = buffer.data(ii + 73);
    const auto *ii_74 = buffer.data(ii + 74);
    const auto *ii_75 = buffer.data(ii + 75);
    const auto *ii_76 = buffer.data(ii + 76);
    const auto *ii_77 = buffer.data(ii + 77);
    const auto *ii_78 = buffer.data(ii + 78);
    const auto *ii_79 = buffer.data(ii + 79);
    const auto *ii_80 = buffer.data(ii + 80);
    const auto *ii_81 = buffer.data(ii + 81);
    const auto *ii_82 = buffer.data(ii + 82);
    const auto *ii_83 = buffer.data(ii + 83);

    const auto *li_56 = buffer.data(li + 56);
    const auto *li_57 = buffer.data(li + 57);
    const auto *li_58 = buffer.data(li + 58);
    const auto *li_59 = buffer.data(li + 59);
    const auto *li_60 = buffer.data(li + 60);
    const auto *li_61 = buffer.data(li + 61);
    const auto *li_62 = buffer.data(li + 62);
    const auto *li_63 = buffer.data(li + 63);
    const auto *li_64 = buffer.data(li + 64);
    const auto *li_65 = buffer.data(li + 65);
    const auto *li_66 = buffer.data(li + 66);
    const auto *li_67 = buffer.data(li + 67);
    const auto *li_68 = buffer.data(li + 68);
    const auto *li_69 = buffer.data(li + 69);
    const auto *li_70 = buffer.data(li + 70);
    const auto *li_71 = buffer.data(li + 71);
    const auto *li_72 = buffer.data(li + 72);
    const auto *li_73 = buffer.data(li + 73);
    const auto *li_74 = buffer.data(li + 74);
    const auto *li_75 = buffer.data(li + 75);
    const auto *li_76 = buffer.data(li + 76);
    const auto *li_77 = buffer.data(li + 77);
    const auto *li_78 = buffer.data(li + 78);
    const auto *li_79 = buffer.data(li + 79);
    const auto *li_80 = buffer.data(li + 80);
    const auto *li_81 = buffer.data(li + 81);
    const auto *li_82 = buffer.data(li + 82);
    const auto *li_83 = buffer.data(li + 83);
    const auto *li_112 = buffer.data(li + 112);
    const auto *li_113 = buffer.data(li + 113);
    const auto *li_114 = buffer.data(li + 114);
    const auto *li_115 = buffer.data(li + 115);
    const auto *li_116 = buffer.data(li + 116);
    const auto *li_117 = buffer.data(li + 117);
    const auto *li_118 = buffer.data(li + 118);
    const auto *li_119 = buffer.data(li + 119);
    const auto *li_120 = buffer.data(li + 120);
    const auto *li_121 = buffer.data(li + 121);
    const auto *li_122 = buffer.data(li + 122);
    const auto *li_123 = buffer.data(li + 123);
    const auto *li_124 = buffer.data(li + 124);
    const auto *li_125 = buffer.data(li + 125);
    const auto *li_126 = buffer.data(li + 126);
    const auto *li_127 = buffer.data(li + 127);
    const auto *li_128 = buffer.data(li + 128);
    const auto *li_129 = buffer.data(li + 129);
    const auto *li_130 = buffer.data(li + 130);
    const auto *li_131 = buffer.data(li + 131);
    const auto *li_132 = buffer.data(li + 132);
    const auto *li_133 = buffer.data(li + 133);
    const auto *li_134 = buffer.data(li + 134);
    const auto *li_135 = buffer.data(li + 135);
    const auto *li_136 = buffer.data(li + 136);
    const auto *li_137 = buffer.data(li + 137);
    const auto *li_138 = buffer.data(li + 138);
    const auto *li_139 = buffer.data(li + 139);
    const auto *li_140 = buffer.data(li + 140);
    const auto *li_141 = buffer.data(li + 141);
    const auto *li_142 = buffer.data(li + 142);
    const auto *li_143 = buffer.data(li + 143);
    const auto *li_144 = buffer.data(li + 144);
    const auto *li_145 = buffer.data(li + 145);
    const auto *li_146 = buffer.data(li + 146);
    const auto *li_147 = buffer.data(li + 147);
    const auto *li_148 = buffer.data(li + 148);
    const auto *li_149 = buffer.data(li + 149);
    const auto *li_150 = buffer.data(li + 150);
    const auto *li_151 = buffer.data(li + 151);
    const auto *li_152 = buffer.data(li + 152);
    const auto *li_153 = buffer.data(li + 153);
    const auto *li_154 = buffer.data(li + 154);
    const auto *li_155 = buffer.data(li + 155);
    const auto *li_156 = buffer.data(li + 156);
    const auto *li_157 = buffer.data(li + 157);
    const auto *li_158 = buffer.data(li + 158);
    const auto *li_159 = buffer.data(li + 159);
    const auto *li_160 = buffer.data(li + 160);
    const auto *li_161 = buffer.data(li + 161);
    const auto *li_162 = buffer.data(li + 162);
    const auto *li_163 = buffer.data(li + 163);
    const auto *li_164 = buffer.data(li + 164);
    const auto *li_165 = buffer.data(li + 165);
    const auto *li_166 = buffer.data(li + 166);
    const auto *li_167 = buffer.data(li + 167);
    const auto *li_196 = buffer.data(li + 196);
    const auto *li_197 = buffer.data(li + 197);
    const auto *li_198 = buffer.data(li + 198);
    const auto *li_199 = buffer.data(li + 199);
    const auto *li_200 = buffer.data(li + 200);
    const auto *li_201 = buffer.data(li + 201);
    const auto *li_202 = buffer.data(li + 202);
    const auto *li_203 = buffer.data(li + 203);
    const auto *li_204 = buffer.data(li + 204);
    const auto *li_205 = buffer.data(li + 205);
    const auto *li_206 = buffer.data(li + 206);
    const auto *li_207 = buffer.data(li + 207);
    const auto *li_208 = buffer.data(li + 208);
    const auto *li_209 = buffer.data(li + 209);
    const auto *li_210 = buffer.data(li + 210);
    const auto *li_211 = buffer.data(li + 211);
    const auto *li_212 = buffer.data(li + 212);
    const auto *li_213 = buffer.data(li + 213);
    const auto *li_214 = buffer.data(li + 214);
    const auto *li_215 = buffer.data(li + 215);
    const auto *li_216 = buffer.data(li + 216);
    const auto *li_217 = buffer.data(li + 217);
    const auto *li_218 = buffer.data(li + 218);
    const auto *li_219 = buffer.data(li + 219);
    const auto *li_220 = buffer.data(li + 220);
    const auto *li_221 = buffer.data(li + 221);
    const auto *li_222 = buffer.data(li + 222);
    const auto *li_223 = buffer.data(li + 223);
    const auto *li_224 = buffer.data(li + 224);
    const auto *li_225 = buffer.data(li + 225);
    const auto *li_226 = buffer.data(li + 226);
    const auto *li_227 = buffer.data(li + 227);
    const auto *li_228 = buffer.data(li + 228);
    const auto *li_229 = buffer.data(li + 229);
    const auto *li_230 = buffer.data(li + 230);
    const auto *li_231 = buffer.data(li + 231);
    const auto *li_232 = buffer.data(li + 232);
    const auto *li_233 = buffer.data(li + 233);
    const auto *li_234 = buffer.data(li + 234);
    const auto *li_235 = buffer.data(li + 235);
    const auto *li_236 = buffer.data(li + 236);
    const auto *li_237 = buffer.data(li + 237);
    const auto *li_238 = buffer.data(li + 238);
    const auto *li_239 = buffer.data(li + 239);
    const auto *li_240 = buffer.data(li + 240);
    const auto *li_241 = buffer.data(li + 241);
    const auto *li_242 = buffer.data(li + 242);
    const auto *li_243 = buffer.data(li + 243);
    const auto *li_244 = buffer.data(li + 244);
    const auto *li_245 = buffer.data(li + 245);
    const auto *li_246 = buffer.data(li + 246);
    const auto *li_247 = buffer.data(li + 247);
    const auto *li_248 = buffer.data(li + 248);
    const auto *li_249 = buffer.data(li + 249);
    const auto *li_250 = buffer.data(li + 250);
    const auto *li_251 = buffer.data(li + 251);
    const auto *li_252 = buffer.data(li + 252);
    const auto *li_253 = buffer.data(li + 253);
    const auto *li_254 = buffer.data(li + 254);
    const auto *li_255 = buffer.data(li + 255);
    const auto *li_256 = buffer.data(li + 256);
    const auto *li_257 = buffer.data(li + 257);
    const auto *li_258 = buffer.data(li + 258);
    const auto *li_259 = buffer.data(li + 259);
    const auto *li_260 = buffer.data(li + 260);
    const auto *li_261 = buffer.data(li + 261);
    const auto *li_262 = buffer.data(li + 262);
    const auto *li_263 = buffer.data(li + 263);
    const auto *li_264 = buffer.data(li + 264);
    const auto *li_265 = buffer.data(li + 265);
    const auto *li_266 = buffer.data(li + 266);
    const auto *li_267 = buffer.data(li + 267);
    const auto *li_268 = buffer.data(li + 268);
    const auto *li_269 = buffer.data(li + 269);
    const auto *li_270 = buffer.data(li + 270);
    const auto *li_271 = buffer.data(li + 271);
    const auto *li_272 = buffer.data(li + 272);
    const auto *li_273 = buffer.data(li + 273);
    const auto *li_274 = buffer.data(li + 274);
    const auto *li_275 = buffer.data(li + 275);
    const auto *li_276 = buffer.data(li + 276);
    const auto *li_277 = buffer.data(li + 277);
    const auto *li_278 = buffer.data(li + 278);
    const auto *li_279 = buffer.data(li + 279);
    const auto *li_308 = buffer.data(li + 308);
    const auto *li_309 = buffer.data(li + 309);
    const auto *li_310 = buffer.data(li + 310);
    const auto *li_311 = buffer.data(li + 311);
    const auto *li_312 = buffer.data(li + 312);
    const auto *li_313 = buffer.data(li + 313);
    const auto *li_314 = buffer.data(li + 314);
    const auto *li_315 = buffer.data(li + 315);
    const auto *li_316 = buffer.data(li + 316);
    const auto *li_317 = buffer.data(li + 317);
    const auto *li_318 = buffer.data(li + 318);
    const auto *li_319 = buffer.data(li + 319);
    const auto *li_320 = buffer.data(li + 320);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, li_56, li_57, li_58, li_59, \
                         li_60, li_61, li_62, li_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * li_56[k];

        t_1[k] = f_0 * li_57[k];

        t_2[k] = f_0 * li_58[k];

        t_3[k] = f_0 * li_59[k];

        t_4[k] = f_0 * li_60[k];

        t_5[k] = f_0 * li_61[k];

        t_6[k] = f_0 * li_62[k];

        t_7[k] = f_0 * li_63[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, li_64, li_65, li_66, \
                         li_67, li_68, li_69, li_70, li_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * li_64[k];

        t_9[k] = f_0 * li_65[k];

        t_10[k] = f_0 * li_66[k];

        t_11[k] = f_0 * li_67[k];

        t_12[k] = f_0 * li_68[k];

        t_13[k] = f_0 * li_69[k];

        t_14[k] = f_0 * li_70[k];

        t_15[k] = f_0 * li_71[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, li_72, li_73, li_74, \
                         li_75, li_76, li_77, li_78, li_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * li_72[k];

        t_17[k] = f_0 * li_73[k];

        t_18[k] = f_0 * li_74[k];

        t_19[k] = f_0 * li_75[k];

        t_20[k] = f_0 * li_76[k];

        t_21[k] = f_0 * li_77[k];

        t_22[k] = f_0 * li_78[k];

        t_23[k] = f_0 * li_79[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, li_80, li_81, li_82, \
                         li_83, li_112, li_113, li_114, li_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * li_80[k];

        t_25[k] = f_0 * li_81[k];

        t_26[k] = f_0 * li_82[k];

        t_27[k] = f_0 * li_83[k];

        t_28[k] = f_0 * li_112[k];

        t_29[k] = f_0 * li_113[k];

        t_30[k] = f_0 * li_114[k];

        t_31[k] = f_0 * li_115[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, li_116, li_117, \
                         li_118, li_119, li_120, li_121, li_122, \
                         li_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * li_116[k];

        t_33[k] = f_0 * li_117[k];

        t_34[k] = f_0 * li_118[k];

        t_35[k] = f_0 * li_119[k];

        t_36[k] = f_0 * li_120[k];

        t_37[k] = f_0 * li_121[k];

        t_38[k] = f_0 * li_122[k];

        t_39[k] = f_0 * li_123[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, li_124, li_125, \
                         li_126, li_127, li_128, li_129, li_130, \
                         li_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * li_124[k];

        t_41[k] = f_0 * li_125[k];

        t_42[k] = f_0 * li_126[k];

        t_43[k] = f_0 * li_127[k];

        t_44[k] = f_0 * li_128[k];

        t_45[k] = f_0 * li_129[k];

        t_46[k] = f_0 * li_130[k];

        t_47[k] = f_0 * li_131[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, li_132, li_133, \
                         li_134, li_135, li_136, li_137, li_138, \
                         li_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * li_132[k];

        t_49[k] = f_0 * li_133[k];

        t_50[k] = f_0 * li_134[k];

        t_51[k] = f_0 * li_135[k];

        t_52[k] = f_0 * li_136[k];

        t_53[k] = f_0 * li_137[k];

        t_54[k] = f_0 * li_138[k];

        t_55[k] = f_0 * li_139[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, ii_0, ii_1, ii_2, ii_3, ii_4, li_140, \
                         li_141, li_142, li_143, li_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -ii_0[k]
                  + f_0 * li_140[k];

        t_57[k] = -ii_1[k]
                  + f_0 * li_141[k];

        t_58[k] = -ii_2[k]
                  + f_0 * li_142[k];

        t_59[k] = -ii_3[k]
                  + f_0 * li_143[k];

        t_60[k] = -ii_4[k]
                  + f_0 * li_144[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, ii_5, ii_6, ii_7, ii_8, ii_9, li_145, \
                         li_146, li_147, li_148, li_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -ii_5[k]
                  + f_0 * li_145[k];

        t_62[k] = -ii_6[k]
                  + f_0 * li_146[k];

        t_63[k] = -ii_7[k]
                  + f_0 * li_147[k];

        t_64[k] = -ii_8[k]
                  + f_0 * li_148[k];

        t_65[k] = -ii_9[k]
                  + f_0 * li_149[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ii_10, ii_11, ii_12, ii_13, ii_14, \
                         li_150, li_151, li_152, li_153, li_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -ii_10[k]
                  + f_0 * li_150[k];

        t_67[k] = -ii_11[k]
                  + f_0 * li_151[k];

        t_68[k] = -ii_12[k]
                  + f_0 * li_152[k];

        t_69[k] = -ii_13[k]
                  + f_0 * li_153[k];

        t_70[k] = -ii_14[k]
                  + f_0 * li_154[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ii_15, ii_16, ii_17, ii_18, ii_19, \
                         li_155, li_156, li_157, li_158, li_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -ii_15[k]
                  + f_0 * li_155[k];

        t_72[k] = -ii_16[k]
                  + f_0 * li_156[k];

        t_73[k] = -ii_17[k]
                  + f_0 * li_157[k];

        t_74[k] = -ii_18[k]
                  + f_0 * li_158[k];

        t_75[k] = -ii_19[k]
                  + f_0 * li_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, ii_20, ii_21, ii_22, ii_23, ii_24, \
                         li_160, li_161, li_162, li_163, li_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -ii_20[k]
                  + f_0 * li_160[k];

        t_77[k] = -ii_21[k]
                  + f_0 * li_161[k];

        t_78[k] = -ii_22[k]
                  + f_0 * li_162[k];

        t_79[k] = -ii_23[k]
                  + f_0 * li_163[k];

        t_80[k] = -ii_24[k]
                  + f_0 * li_164[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, ii_25, ii_26, ii_27, li_165, \
                         li_166, li_167, li_196, li_197, li_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -ii_25[k]
                  + f_0 * li_165[k];

        t_82[k] = -ii_26[k]
                  + f_0 * li_166[k];

        t_83[k] = -ii_27[k]
                  + f_0 * li_167[k];

        t_84[k] = f_0 * li_196[k];

        t_85[k] = f_0 * li_197[k];

        t_86[k] = f_0 * li_198[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, t_93, t_94, li_199, li_200, \
                         li_201, li_202, li_203, li_204, li_205, \
                         li_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_0 * li_199[k];

        t_88[k] = f_0 * li_200[k];

        t_89[k] = f_0 * li_201[k];

        t_90[k] = f_0 * li_202[k];

        t_91[k] = f_0 * li_203[k];

        t_92[k] = f_0 * li_204[k];

        t_93[k] = f_0 * li_205[k];

        t_94[k] = f_0 * li_206[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, t_101, t_102, li_207, li_208, \
                         li_209, li_210, li_211, li_212, li_213, \
                         li_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * li_207[k];

        t_96[k] = f_0 * li_208[k];

        t_97[k] = f_0 * li_209[k];

        t_98[k] = f_0 * li_210[k];

        t_99[k] = f_0 * li_211[k];

        t_100[k] = f_0 * li_212[k];

        t_101[k] = f_0 * li_213[k];

        t_102[k] = f_0 * li_214[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, t_109, t_110, li_215, \
                         li_216, li_217, li_218, li_219, li_220, li_221, \
                         li_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * li_215[k];

        t_104[k] = f_0 * li_216[k];

        t_105[k] = f_0 * li_217[k];

        t_106[k] = f_0 * li_218[k];

        t_107[k] = f_0 * li_219[k];

        t_108[k] = f_0 * li_220[k];

        t_109[k] = f_0 * li_221[k];

        t_110[k] = f_0 * li_222[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, ii_28, ii_29, ii_30, ii_31, \
                         li_223, li_224, li_225, li_226, li_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * li_223[k];

        t_112[k] = -ii_28[k]
                   + f_0 * li_224[k];

        t_113[k] = -ii_29[k]
                   + f_0 * li_225[k];

        t_114[k] = -ii_30[k]
                   + f_0 * li_226[k];

        t_115[k] = -ii_31[k]
                   + f_0 * li_227[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, ii_32, ii_33, ii_34, ii_35, ii_36, \
                         li_228, li_229, li_230, li_231, li_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -ii_32[k]
                   + f_0 * li_228[k];

        t_117[k] = -ii_33[k]
                   + f_0 * li_229[k];

        t_118[k] = -ii_34[k]
                   + f_0 * li_230[k];

        t_119[k] = -ii_35[k]
                   + f_0 * li_231[k];

        t_120[k] = -ii_36[k]
                   + f_0 * li_232[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, ii_37, ii_38, ii_39, ii_40, ii_41, \
                         li_233, li_234, li_235, li_236, li_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -ii_37[k]
                   + f_0 * li_233[k];

        t_122[k] = -ii_38[k]
                   + f_0 * li_234[k];

        t_123[k] = -ii_39[k]
                   + f_0 * li_235[k];

        t_124[k] = -ii_40[k]
                   + f_0 * li_236[k];

        t_125[k] = -ii_41[k]
                   + f_0 * li_237[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, ii_42, ii_43, ii_44, ii_45, ii_46, \
                         li_238, li_239, li_240, li_241, li_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -ii_42[k]
                   + f_0 * li_238[k];

        t_127[k] = -ii_43[k]
                   + f_0 * li_239[k];

        t_128[k] = -ii_44[k]
                   + f_0 * li_240[k];

        t_129[k] = -ii_45[k]
                   + f_0 * li_241[k];

        t_130[k] = -ii_46[k]
                   + f_0 * li_242[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, ii_47, ii_48, ii_49, ii_50, ii_51, \
                         li_243, li_244, li_245, li_246, li_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -ii_47[k]
                   + f_0 * li_243[k];

        t_132[k] = -ii_48[k]
                   + f_0 * li_244[k];

        t_133[k] = -ii_49[k]
                   + f_0 * li_245[k];

        t_134[k] = -ii_50[k]
                   + f_0 * li_246[k];

        t_135[k] = -ii_51[k]
                   + f_0 * li_247[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, ii_52, ii_53, ii_54, ii_55, ii_56, \
                         li_248, li_249, li_250, li_251, li_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -ii_52[k]
                   + f_0 * li_248[k];

        t_137[k] = -ii_53[k]
                   + f_0 * li_249[k];

        t_138[k] = -ii_54[k]
                   + f_0 * li_250[k];

        t_139[k] = -ii_55[k]
                   + f_0 * li_251[k];

        t_140[k] = -2.0 * ii_56[k]
                   + f_0 * li_252[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, ii_57, ii_58, ii_59, ii_60, ii_61, \
                         li_253, li_254, li_255, li_256, li_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = -2.0 * ii_57[k]
                   + f_0 * li_253[k];

        t_142[k] = -2.0 * ii_58[k]
                   + f_0 * li_254[k];

        t_143[k] = -2.0 * ii_59[k]
                   + f_0 * li_255[k];

        t_144[k] = -2.0 * ii_60[k]
                   + f_0 * li_256[k];

        t_145[k] = -2.0 * ii_61[k]
                   + f_0 * li_257[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, ii_62, ii_63, ii_64, ii_65, ii_66, \
                         li_258, li_259, li_260, li_261, li_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -2.0 * ii_62[k]
                   + f_0 * li_258[k];

        t_147[k] = -2.0 * ii_63[k]
                   + f_0 * li_259[k];

        t_148[k] = -2.0 * ii_64[k]
                   + f_0 * li_260[k];

        t_149[k] = -2.0 * ii_65[k]
                   + f_0 * li_261[k];

        t_150[k] = -2.0 * ii_66[k]
                   + f_0 * li_262[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, ii_67, ii_68, ii_69, ii_70, ii_71, \
                         li_263, li_264, li_265, li_266, li_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = -2.0 * ii_67[k]
                   + f_0 * li_263[k];

        t_152[k] = -2.0 * ii_68[k]
                   + f_0 * li_264[k];

        t_153[k] = -2.0 * ii_69[k]
                   + f_0 * li_265[k];

        t_154[k] = -2.0 * ii_70[k]
                   + f_0 * li_266[k];

        t_155[k] = -2.0 * ii_71[k]
                   + f_0 * li_267[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, ii_72, ii_73, ii_74, ii_75, ii_76, \
                         li_268, li_269, li_270, li_271, li_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = -2.0 * ii_72[k]
                   + f_0 * li_268[k];

        t_157[k] = -2.0 * ii_73[k]
                   + f_0 * li_269[k];

        t_158[k] = -2.0 * ii_74[k]
                   + f_0 * li_270[k];

        t_159[k] = -2.0 * ii_75[k]
                   + f_0 * li_271[k];

        t_160[k] = -2.0 * ii_76[k]
                   + f_0 * li_272[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, ii_77, ii_78, ii_79, ii_80, ii_81, \
                         li_273, li_274, li_275, li_276, li_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = -2.0 * ii_77[k]
                   + f_0 * li_273[k];

        t_162[k] = -2.0 * ii_78[k]
                   + f_0 * li_274[k];

        t_163[k] = -2.0 * ii_79[k]
                   + f_0 * li_275[k];

        t_164[k] = -2.0 * ii_80[k]
                   + f_0 * li_276[k];

        t_165[k] = -2.0 * ii_81[k]
                   + f_0 * li_277[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, t_172, ii_82, ii_83, \
                         li_278, li_279, li_308, li_309, li_310, li_311, \
                         li_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = -2.0 * ii_82[k]
                   + f_0 * li_278[k];

        t_167[k] = -2.0 * ii_83[k]
                   + f_0 * li_279[k];

        t_168[k] = f_0 * li_308[k];

        t_169[k] = f_0 * li_309[k];

        t_170[k] = f_0 * li_310[k];

        t_171[k] = f_0 * li_311[k];

        t_172[k] = f_0 * li_312[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, t_180, li_313, \
                         li_314, li_315, li_316, li_317, li_318, li_319, \
                         li_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * li_313[k];

        t_174[k] = f_0 * li_314[k];

        t_175[k] = f_0 * li_315[k];

        t_176[k] = f_0 * li_316[k];

        t_177[k] = f_0 * li_317[k];

        t_178[k] = f_0 * li_318[k];

        t_179[k] = f_0 * li_319[k];

        t_180[k] = f_0 * li_320[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_84 = buffer.data(ii + 84);
    const auto *ii_85 = buffer.data(ii + 85);
    const auto *ii_86 = buffer.data(ii + 86);
    const auto *ii_87 = buffer.data(ii + 87);
    const auto *ii_88 = buffer.data(ii + 88);
    const auto *ii_89 = buffer.data(ii + 89);
    const auto *ii_90 = buffer.data(ii + 90);
    const auto *ii_91 = buffer.data(ii + 91);
    const auto *ii_92 = buffer.data(ii + 92);
    const auto *ii_93 = buffer.data(ii + 93);
    const auto *ii_94 = buffer.data(ii + 94);
    const auto *ii_95 = buffer.data(ii + 95);
    const auto *ii_96 = buffer.data(ii + 96);
    const auto *ii_97 = buffer.data(ii + 97);
    const auto *ii_98 = buffer.data(ii + 98);
    const auto *ii_99 = buffer.data(ii + 99);
    const auto *ii_100 = buffer.data(ii + 100);
    const auto *ii_101 = buffer.data(ii + 101);
    const auto *ii_102 = buffer.data(ii + 102);
    const auto *ii_103 = buffer.data(ii + 103);
    const auto *ii_104 = buffer.data(ii + 104);
    const auto *ii_105 = buffer.data(ii + 105);
    const auto *ii_106 = buffer.data(ii + 106);
    const auto *ii_107 = buffer.data(ii + 107);
    const auto *ii_108 = buffer.data(ii + 108);
    const auto *ii_109 = buffer.data(ii + 109);
    const auto *ii_110 = buffer.data(ii + 110);
    const auto *ii_111 = buffer.data(ii + 111);
    const auto *ii_112 = buffer.data(ii + 112);
    const auto *ii_113 = buffer.data(ii + 113);
    const auto *ii_114 = buffer.data(ii + 114);
    const auto *ii_115 = buffer.data(ii + 115);
    const auto *ii_116 = buffer.data(ii + 116);
    const auto *ii_117 = buffer.data(ii + 117);
    const auto *ii_118 = buffer.data(ii + 118);
    const auto *ii_119 = buffer.data(ii + 119);
    const auto *ii_120 = buffer.data(ii + 120);
    const auto *ii_121 = buffer.data(ii + 121);
    const auto *ii_122 = buffer.data(ii + 122);
    const auto *ii_123 = buffer.data(ii + 123);
    const auto *ii_124 = buffer.data(ii + 124);
    const auto *ii_125 = buffer.data(ii + 125);
    const auto *ii_126 = buffer.data(ii + 126);
    const auto *ii_127 = buffer.data(ii + 127);
    const auto *ii_128 = buffer.data(ii + 128);
    const auto *ii_129 = buffer.data(ii + 129);
    const auto *ii_130 = buffer.data(ii + 130);
    const auto *ii_131 = buffer.data(ii + 131);
    const auto *ii_132 = buffer.data(ii + 132);
    const auto *ii_133 = buffer.data(ii + 133);
    const auto *ii_134 = buffer.data(ii + 134);
    const auto *ii_135 = buffer.data(ii + 135);
    const auto *ii_136 = buffer.data(ii + 136);
    const auto *ii_137 = buffer.data(ii + 137);
    const auto *ii_138 = buffer.data(ii + 138);
    const auto *ii_139 = buffer.data(ii + 139);
    const auto *ii_140 = buffer.data(ii + 140);
    const auto *ii_141 = buffer.data(ii + 141);
    const auto *ii_142 = buffer.data(ii + 142);
    const auto *ii_143 = buffer.data(ii + 143);
    const auto *ii_144 = buffer.data(ii + 144);
    const auto *ii_145 = buffer.data(ii + 145);
    const auto *ii_146 = buffer.data(ii + 146);
    const auto *ii_147 = buffer.data(ii + 147);
    const auto *ii_148 = buffer.data(ii + 148);
    const auto *ii_149 = buffer.data(ii + 149);
    const auto *ii_150 = buffer.data(ii + 150);
    const auto *ii_151 = buffer.data(ii + 151);
    const auto *ii_152 = buffer.data(ii + 152);
    const auto *ii_153 = buffer.data(ii + 153);
    const auto *ii_154 = buffer.data(ii + 154);
    const auto *ii_155 = buffer.data(ii + 155);
    const auto *ii_156 = buffer.data(ii + 156);
    const auto *ii_157 = buffer.data(ii + 157);
    const auto *ii_158 = buffer.data(ii + 158);
    const auto *ii_159 = buffer.data(ii + 159);
    const auto *ii_160 = buffer.data(ii + 160);
    const auto *ii_161 = buffer.data(ii + 161);
    const auto *ii_162 = buffer.data(ii + 162);
    const auto *ii_163 = buffer.data(ii + 163);
    const auto *ii_164 = buffer.data(ii + 164);
    const auto *ii_165 = buffer.data(ii + 165);
    const auto *ii_166 = buffer.data(ii + 166);
    const auto *ii_167 = buffer.data(ii + 167);
    const auto *ii_168 = buffer.data(ii + 168);
    const auto *ii_169 = buffer.data(ii + 169);
    const auto *ii_170 = buffer.data(ii + 170);
    const auto *ii_171 = buffer.data(ii + 171);
    const auto *ii_172 = buffer.data(ii + 172);
    const auto *ii_173 = buffer.data(ii + 173);
    const auto *ii_174 = buffer.data(ii + 174);
    const auto *ii_175 = buffer.data(ii + 175);
    const auto *ii_176 = buffer.data(ii + 176);
    const auto *ii_177 = buffer.data(ii + 177);
    const auto *ii_178 = buffer.data(ii + 178);
    const auto *ii_179 = buffer.data(ii + 179);
    const auto *ii_180 = buffer.data(ii + 180);
    const auto *ii_181 = buffer.data(ii + 181);
    const auto *ii_182 = buffer.data(ii + 182);
    const auto *ii_183 = buffer.data(ii + 183);
    const auto *ii_184 = buffer.data(ii + 184);
    const auto *ii_185 = buffer.data(ii + 185);
    const auto *ii_186 = buffer.data(ii + 186);
    const auto *ii_187 = buffer.data(ii + 187);
    const auto *ii_188 = buffer.data(ii + 188);
    const auto *ii_189 = buffer.data(ii + 189);
    const auto *ii_190 = buffer.data(ii + 190);
    const auto *ii_191 = buffer.data(ii + 191);
    const auto *ii_192 = buffer.data(ii + 192);
    const auto *ii_193 = buffer.data(ii + 193);
    const auto *ii_194 = buffer.data(ii + 194);
    const auto *ii_195 = buffer.data(ii + 195);
    const auto *ii_196 = buffer.data(ii + 196);
    const auto *ii_197 = buffer.data(ii + 197);
    const auto *ii_198 = buffer.data(ii + 198);
    const auto *ii_199 = buffer.data(ii + 199);
    const auto *ii_200 = buffer.data(ii + 200);
    const auto *ii_201 = buffer.data(ii + 201);
    const auto *ii_202 = buffer.data(ii + 202);
    const auto *ii_203 = buffer.data(ii + 203);
    const auto *ii_204 = buffer.data(ii + 204);
    const auto *ii_205 = buffer.data(ii + 205);
    const auto *ii_206 = buffer.data(ii + 206);

    const auto *li_321 = buffer.data(li + 321);
    const auto *li_322 = buffer.data(li + 322);
    const auto *li_323 = buffer.data(li + 323);
    const auto *li_324 = buffer.data(li + 324);
    const auto *li_325 = buffer.data(li + 325);
    const auto *li_326 = buffer.data(li + 326);
    const auto *li_327 = buffer.data(li + 327);
    const auto *li_328 = buffer.data(li + 328);
    const auto *li_329 = buffer.data(li + 329);
    const auto *li_330 = buffer.data(li + 330);
    const auto *li_331 = buffer.data(li + 331);
    const auto *li_332 = buffer.data(li + 332);
    const auto *li_333 = buffer.data(li + 333);
    const auto *li_334 = buffer.data(li + 334);
    const auto *li_335 = buffer.data(li + 335);
    const auto *li_336 = buffer.data(li + 336);
    const auto *li_337 = buffer.data(li + 337);
    const auto *li_338 = buffer.data(li + 338);
    const auto *li_339 = buffer.data(li + 339);
    const auto *li_340 = buffer.data(li + 340);
    const auto *li_341 = buffer.data(li + 341);
    const auto *li_342 = buffer.data(li + 342);
    const auto *li_343 = buffer.data(li + 343);
    const auto *li_344 = buffer.data(li + 344);
    const auto *li_345 = buffer.data(li + 345);
    const auto *li_346 = buffer.data(li + 346);
    const auto *li_347 = buffer.data(li + 347);
    const auto *li_348 = buffer.data(li + 348);
    const auto *li_349 = buffer.data(li + 349);
    const auto *li_350 = buffer.data(li + 350);
    const auto *li_351 = buffer.data(li + 351);
    const auto *li_352 = buffer.data(li + 352);
    const auto *li_353 = buffer.data(li + 353);
    const auto *li_354 = buffer.data(li + 354);
    const auto *li_355 = buffer.data(li + 355);
    const auto *li_356 = buffer.data(li + 356);
    const auto *li_357 = buffer.data(li + 357);
    const auto *li_358 = buffer.data(li + 358);
    const auto *li_359 = buffer.data(li + 359);
    const auto *li_360 = buffer.data(li + 360);
    const auto *li_361 = buffer.data(li + 361);
    const auto *li_362 = buffer.data(li + 362);
    const auto *li_363 = buffer.data(li + 363);
    const auto *li_364 = buffer.data(li + 364);
    const auto *li_365 = buffer.data(li + 365);
    const auto *li_366 = buffer.data(li + 366);
    const auto *li_367 = buffer.data(li + 367);
    const auto *li_368 = buffer.data(li + 368);
    const auto *li_369 = buffer.data(li + 369);
    const auto *li_370 = buffer.data(li + 370);
    const auto *li_371 = buffer.data(li + 371);
    const auto *li_372 = buffer.data(li + 372);
    const auto *li_373 = buffer.data(li + 373);
    const auto *li_374 = buffer.data(li + 374);
    const auto *li_375 = buffer.data(li + 375);
    const auto *li_376 = buffer.data(li + 376);
    const auto *li_377 = buffer.data(li + 377);
    const auto *li_378 = buffer.data(li + 378);
    const auto *li_379 = buffer.data(li + 379);
    const auto *li_380 = buffer.data(li + 380);
    const auto *li_381 = buffer.data(li + 381);
    const auto *li_382 = buffer.data(li + 382);
    const auto *li_383 = buffer.data(li + 383);
    const auto *li_384 = buffer.data(li + 384);
    const auto *li_385 = buffer.data(li + 385);
    const auto *li_386 = buffer.data(li + 386);
    const auto *li_387 = buffer.data(li + 387);
    const auto *li_388 = buffer.data(li + 388);
    const auto *li_389 = buffer.data(li + 389);
    const auto *li_390 = buffer.data(li + 390);
    const auto *li_391 = buffer.data(li + 391);
    const auto *li_392 = buffer.data(li + 392);
    const auto *li_393 = buffer.data(li + 393);
    const auto *li_394 = buffer.data(li + 394);
    const auto *li_395 = buffer.data(li + 395);
    const auto *li_396 = buffer.data(li + 396);
    const auto *li_397 = buffer.data(li + 397);
    const auto *li_398 = buffer.data(li + 398);
    const auto *li_399 = buffer.data(li + 399);
    const auto *li_400 = buffer.data(li + 400);
    const auto *li_401 = buffer.data(li + 401);
    const auto *li_402 = buffer.data(li + 402);
    const auto *li_403 = buffer.data(li + 403);
    const auto *li_404 = buffer.data(li + 404);
    const auto *li_405 = buffer.data(li + 405);
    const auto *li_406 = buffer.data(li + 406);
    const auto *li_407 = buffer.data(li + 407);
    const auto *li_408 = buffer.data(li + 408);
    const auto *li_409 = buffer.data(li + 409);
    const auto *li_410 = buffer.data(li + 410);
    const auto *li_411 = buffer.data(li + 411);
    const auto *li_412 = buffer.data(li + 412);
    const auto *li_413 = buffer.data(li + 413);
    const auto *li_414 = buffer.data(li + 414);
    const auto *li_415 = buffer.data(li + 415);
    const auto *li_416 = buffer.data(li + 416);
    const auto *li_417 = buffer.data(li + 417);
    const auto *li_418 = buffer.data(li + 418);
    const auto *li_419 = buffer.data(li + 419);
    const auto *li_448 = buffer.data(li + 448);
    const auto *li_449 = buffer.data(li + 449);
    const auto *li_450 = buffer.data(li + 450);
    const auto *li_451 = buffer.data(li + 451);
    const auto *li_452 = buffer.data(li + 452);
    const auto *li_453 = buffer.data(li + 453);
    const auto *li_454 = buffer.data(li + 454);
    const auto *li_455 = buffer.data(li + 455);
    const auto *li_456 = buffer.data(li + 456);
    const auto *li_457 = buffer.data(li + 457);
    const auto *li_458 = buffer.data(li + 458);
    const auto *li_459 = buffer.data(li + 459);
    const auto *li_460 = buffer.data(li + 460);
    const auto *li_461 = buffer.data(li + 461);
    const auto *li_462 = buffer.data(li + 462);
    const auto *li_463 = buffer.data(li + 463);
    const auto *li_464 = buffer.data(li + 464);
    const auto *li_465 = buffer.data(li + 465);
    const auto *li_466 = buffer.data(li + 466);
    const auto *li_467 = buffer.data(li + 467);
    const auto *li_468 = buffer.data(li + 468);
    const auto *li_469 = buffer.data(li + 469);
    const auto *li_470 = buffer.data(li + 470);
    const auto *li_471 = buffer.data(li + 471);
    const auto *li_472 = buffer.data(li + 472);
    const auto *li_473 = buffer.data(li + 473);
    const auto *li_474 = buffer.data(li + 474);
    const auto *li_475 = buffer.data(li + 475);
    const auto *li_476 = buffer.data(li + 476);
    const auto *li_477 = buffer.data(li + 477);
    const auto *li_478 = buffer.data(li + 478);
    const auto *li_479 = buffer.data(li + 479);
    const auto *li_480 = buffer.data(li + 480);
    const auto *li_481 = buffer.data(li + 481);
    const auto *li_482 = buffer.data(li + 482);
    const auto *li_483 = buffer.data(li + 483);
    const auto *li_484 = buffer.data(li + 484);
    const auto *li_485 = buffer.data(li + 485);
    const auto *li_486 = buffer.data(li + 486);
    const auto *li_487 = buffer.data(li + 487);
    const auto *li_488 = buffer.data(li + 488);
    const auto *li_489 = buffer.data(li + 489);
    const auto *li_490 = buffer.data(li + 490);
    const auto *li_491 = buffer.data(li + 491);
    const auto *li_492 = buffer.data(li + 492);
    const auto *li_493 = buffer.data(li + 493);
    const auto *li_494 = buffer.data(li + 494);
    const auto *li_495 = buffer.data(li + 495);
    const auto *li_496 = buffer.data(li + 496);
    const auto *li_497 = buffer.data(li + 497);
    const auto *li_498 = buffer.data(li + 498);
    const auto *li_499 = buffer.data(li + 499);
    const auto *li_500 = buffer.data(li + 500);
    const auto *li_501 = buffer.data(li + 501);
    const auto *li_502 = buffer.data(li + 502);
    const auto *li_503 = buffer.data(li + 503);
    const auto *li_504 = buffer.data(li + 504);
    const auto *li_505 = buffer.data(li + 505);
    const auto *li_506 = buffer.data(li + 506);
    const auto *li_507 = buffer.data(li + 507);
    const auto *li_508 = buffer.data(li + 508);
    const auto *li_509 = buffer.data(li + 509);
    const auto *li_510 = buffer.data(li + 510);
    const auto *li_511 = buffer.data(li + 511);
    const auto *li_512 = buffer.data(li + 512);
    const auto *li_513 = buffer.data(li + 513);
    const auto *li_514 = buffer.data(li + 514);

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, t_186, t_187, t_188, li_321, \
                         li_322, li_323, li_324, li_325, li_326, li_327, \
                         li_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_0 * li_321[k];

        t_182[k] = f_0 * li_322[k];

        t_183[k] = f_0 * li_323[k];

        t_184[k] = f_0 * li_324[k];

        t_185[k] = f_0 * li_325[k];

        t_186[k] = f_0 * li_326[k];

        t_187[k] = f_0 * li_327[k];

        t_188[k] = f_0 * li_328[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, li_329, li_330, \
                         li_331, li_332, li_333, li_334, li_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * li_329[k];

        t_190[k] = f_0 * li_330[k];

        t_191[k] = f_0 * li_331[k];

        t_192[k] = f_0 * li_332[k];

        t_193[k] = f_0 * li_333[k];

        t_194[k] = f_0 * li_334[k];

        t_195[k] = f_0 * li_335[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, ii_84, ii_85, ii_86, ii_87, ii_88, \
                         li_336, li_337, li_338, li_339, li_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -ii_84[k]
                   + f_0 * li_336[k];

        t_197[k] = -ii_85[k]
                   + f_0 * li_337[k];

        t_198[k] = -ii_86[k]
                   + f_0 * li_338[k];

        t_199[k] = -ii_87[k]
                   + f_0 * li_339[k];

        t_200[k] = -ii_88[k]
                   + f_0 * li_340[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, ii_89, ii_90, ii_91, ii_92, ii_93, \
                         li_341, li_342, li_343, li_344, li_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -ii_89[k]
                   + f_0 * li_341[k];

        t_202[k] = -ii_90[k]
                   + f_0 * li_342[k];

        t_203[k] = -ii_91[k]
                   + f_0 * li_343[k];

        t_204[k] = -ii_92[k]
                   + f_0 * li_344[k];

        t_205[k] = -ii_93[k]
                   + f_0 * li_345[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, ii_94, ii_95, ii_96, ii_97, ii_98, \
                         li_346, li_347, li_348, li_349, li_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = -ii_94[k]
                   + f_0 * li_346[k];

        t_207[k] = -ii_95[k]
                   + f_0 * li_347[k];

        t_208[k] = -ii_96[k]
                   + f_0 * li_348[k];

        t_209[k] = -ii_97[k]
                   + f_0 * li_349[k];

        t_210[k] = -ii_98[k]
                   + f_0 * li_350[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, ii_99, ii_100, ii_101, ii_102, \
                         ii_103, li_351, li_352, li_353, li_354, \
                         li_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -ii_99[k]
                   + f_0 * li_351[k];

        t_212[k] = -ii_100[k]
                   + f_0 * li_352[k];

        t_213[k] = -ii_101[k]
                   + f_0 * li_353[k];

        t_214[k] = -ii_102[k]
                   + f_0 * li_354[k];

        t_215[k] = -ii_103[k]
                   + f_0 * li_355[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, ii_104, ii_105, ii_106, ii_107, \
                         ii_108, li_356, li_357, li_358, li_359, \
                         li_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = -ii_104[k]
                   + f_0 * li_356[k];

        t_217[k] = -ii_105[k]
                   + f_0 * li_357[k];

        t_218[k] = -ii_106[k]
                   + f_0 * li_358[k];

        t_219[k] = -ii_107[k]
                   + f_0 * li_359[k];

        t_220[k] = -ii_108[k]
                   + f_0 * li_360[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, ii_109, ii_110, ii_111, ii_112, \
                         ii_113, li_361, li_362, li_363, li_364, \
                         li_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = -ii_109[k]
                   + f_0 * li_361[k];

        t_222[k] = -ii_110[k]
                   + f_0 * li_362[k];

        t_223[k] = -ii_111[k]
                   + f_0 * li_363[k];

        t_224[k] = -2.0 * ii_112[k]
                   + f_0 * li_364[k];

        t_225[k] = -2.0 * ii_113[k]
                   + f_0 * li_365[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, t_230, ii_114, ii_115, ii_116, ii_117, \
                         ii_118, li_366, li_367, li_368, li_369, \
                         li_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = -2.0 * ii_114[k]
                   + f_0 * li_366[k];

        t_227[k] = -2.0 * ii_115[k]
                   + f_0 * li_367[k];

        t_228[k] = -2.0 * ii_116[k]
                   + f_0 * li_368[k];

        t_229[k] = -2.0 * ii_117[k]
                   + f_0 * li_369[k];

        t_230[k] = -2.0 * ii_118[k]
                   + f_0 * li_370[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, ii_119, ii_120, ii_121, ii_122, \
                         ii_123, li_371, li_372, li_373, li_374, \
                         li_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = -2.0 * ii_119[k]
                   + f_0 * li_371[k];

        t_232[k] = -2.0 * ii_120[k]
                   + f_0 * li_372[k];

        t_233[k] = -2.0 * ii_121[k]
                   + f_0 * li_373[k];

        t_234[k] = -2.0 * ii_122[k]
                   + f_0 * li_374[k];

        t_235[k] = -2.0 * ii_123[k]
                   + f_0 * li_375[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, ii_124, ii_125, ii_126, ii_127, \
                         ii_128, li_376, li_377, li_378, li_379, \
                         li_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = -2.0 * ii_124[k]
                   + f_0 * li_376[k];

        t_237[k] = -2.0 * ii_125[k]
                   + f_0 * li_377[k];

        t_238[k] = -2.0 * ii_126[k]
                   + f_0 * li_378[k];

        t_239[k] = -2.0 * ii_127[k]
                   + f_0 * li_379[k];

        t_240[k] = -2.0 * ii_128[k]
                   + f_0 * li_380[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, ii_129, ii_130, ii_131, ii_132, \
                         ii_133, li_381, li_382, li_383, li_384, \
                         li_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = -2.0 * ii_129[k]
                   + f_0 * li_381[k];

        t_242[k] = -2.0 * ii_130[k]
                   + f_0 * li_382[k];

        t_243[k] = -2.0 * ii_131[k]
                   + f_0 * li_383[k];

        t_244[k] = -2.0 * ii_132[k]
                   + f_0 * li_384[k];

        t_245[k] = -2.0 * ii_133[k]
                   + f_0 * li_385[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, ii_134, ii_135, ii_136, ii_137, \
                         ii_138, li_386, li_387, li_388, li_389, \
                         li_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = -2.0 * ii_134[k]
                   + f_0 * li_386[k];

        t_247[k] = -2.0 * ii_135[k]
                   + f_0 * li_387[k];

        t_248[k] = -2.0 * ii_136[k]
                   + f_0 * li_388[k];

        t_249[k] = -2.0 * ii_137[k]
                   + f_0 * li_389[k];

        t_250[k] = -2.0 * ii_138[k]
                   + f_0 * li_390[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, ii_139, ii_140, ii_141, ii_142, \
                         ii_143, li_391, li_392, li_393, li_394, \
                         li_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = -2.0 * ii_139[k]
                   + f_0 * li_391[k];

        t_252[k] = -3.0 * ii_140[k]
                   + f_0 * li_392[k];

        t_253[k] = -3.0 * ii_141[k]
                   + f_0 * li_393[k];

        t_254[k] = -3.0 * ii_142[k]
                   + f_0 * li_394[k];

        t_255[k] = -3.0 * ii_143[k]
                   + f_0 * li_395[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, ii_144, ii_145, ii_146, ii_147, \
                         ii_148, li_396, li_397, li_398, li_399, \
                         li_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = -3.0 * ii_144[k]
                   + f_0 * li_396[k];

        t_257[k] = -3.0 * ii_145[k]
                   + f_0 * li_397[k];

        t_258[k] = -3.0 * ii_146[k]
                   + f_0 * li_398[k];

        t_259[k] = -3.0 * ii_147[k]
                   + f_0 * li_399[k];

        t_260[k] = -3.0 * ii_148[k]
                   + f_0 * li_400[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, ii_149, ii_150, ii_151, ii_152, \
                         ii_153, li_401, li_402, li_403, li_404, \
                         li_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = -3.0 * ii_149[k]
                   + f_0 * li_401[k];

        t_262[k] = -3.0 * ii_150[k]
                   + f_0 * li_402[k];

        t_263[k] = -3.0 * ii_151[k]
                   + f_0 * li_403[k];

        t_264[k] = -3.0 * ii_152[k]
                   + f_0 * li_404[k];

        t_265[k] = -3.0 * ii_153[k]
                   + f_0 * li_405[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, ii_154, ii_155, ii_156, ii_157, \
                         ii_158, li_406, li_407, li_408, li_409, \
                         li_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = -3.0 * ii_154[k]
                   + f_0 * li_406[k];

        t_267[k] = -3.0 * ii_155[k]
                   + f_0 * li_407[k];

        t_268[k] = -3.0 * ii_156[k]
                   + f_0 * li_408[k];

        t_269[k] = -3.0 * ii_157[k]
                   + f_0 * li_409[k];

        t_270[k] = -3.0 * ii_158[k]
                   + f_0 * li_410[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, ii_159, ii_160, ii_161, ii_162, \
                         ii_163, li_411, li_412, li_413, li_414, \
                         li_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = -3.0 * ii_159[k]
                   + f_0 * li_411[k];

        t_272[k] = -3.0 * ii_160[k]
                   + f_0 * li_412[k];

        t_273[k] = -3.0 * ii_161[k]
                   + f_0 * li_413[k];

        t_274[k] = -3.0 * ii_162[k]
                   + f_0 * li_414[k];

        t_275[k] = -3.0 * ii_163[k]
                   + f_0 * li_415[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, t_281, ii_164, ii_165, ii_166, \
                         ii_167, li_416, li_417, li_418, li_419, li_448, \
                         li_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = -3.0 * ii_164[k]
                   + f_0 * li_416[k];

        t_277[k] = -3.0 * ii_165[k]
                   + f_0 * li_417[k];

        t_278[k] = -3.0 * ii_166[k]
                   + f_0 * li_418[k];

        t_279[k] = -3.0 * ii_167[k]
                   + f_0 * li_419[k];

        t_280[k] = f_0 * li_448[k];

        t_281[k] = f_0 * li_449[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, t_287, t_288, t_289, li_450, \
                         li_451, li_452, li_453, li_454, li_455, li_456, \
                         li_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_0 * li_450[k];

        t_283[k] = f_0 * li_451[k];

        t_284[k] = f_0 * li_452[k];

        t_285[k] = f_0 * li_453[k];

        t_286[k] = f_0 * li_454[k];

        t_287[k] = f_0 * li_455[k];

        t_288[k] = f_0 * li_456[k];

        t_289[k] = f_0 * li_457[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, t_295, t_296, t_297, li_458, \
                         li_459, li_460, li_461, li_462, li_463, li_464, \
                         li_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_0 * li_458[k];

        t_291[k] = f_0 * li_459[k];

        t_292[k] = f_0 * li_460[k];

        t_293[k] = f_0 * li_461[k];

        t_294[k] = f_0 * li_462[k];

        t_295[k] = f_0 * li_463[k];

        t_296[k] = f_0 * li_464[k];

        t_297[k] = f_0 * li_465[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, t_303, t_304, t_305, li_466, \
                         li_467, li_468, li_469, li_470, li_471, li_472, \
                         li_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_0 * li_466[k];

        t_299[k] = f_0 * li_467[k];

        t_300[k] = f_0 * li_468[k];

        t_301[k] = f_0 * li_469[k];

        t_302[k] = f_0 * li_470[k];

        t_303[k] = f_0 * li_471[k];

        t_304[k] = f_0 * li_472[k];

        t_305[k] = f_0 * li_473[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, ii_168, ii_169, ii_170, \
                         ii_171, li_474, li_475, li_476, li_477, li_478, \
                         li_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_0 * li_474[k];

        t_307[k] = f_0 * li_475[k];

        t_308[k] = -ii_168[k]
                   + f_0 * li_476[k];

        t_309[k] = -ii_169[k]
                   + f_0 * li_477[k];

        t_310[k] = -ii_170[k]
                   + f_0 * li_478[k];

        t_311[k] = -ii_171[k]
                   + f_0 * li_479[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, ii_172, ii_173, ii_174, ii_175, \
                         ii_176, li_480, li_481, li_482, li_483, \
                         li_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -ii_172[k]
                   + f_0 * li_480[k];

        t_313[k] = -ii_173[k]
                   + f_0 * li_481[k];

        t_314[k] = -ii_174[k]
                   + f_0 * li_482[k];

        t_315[k] = -ii_175[k]
                   + f_0 * li_483[k];

        t_316[k] = -ii_176[k]
                   + f_0 * li_484[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, ii_177, ii_178, ii_179, ii_180, \
                         ii_181, li_485, li_486, li_487, li_488, \
                         li_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -ii_177[k]
                   + f_0 * li_485[k];

        t_318[k] = -ii_178[k]
                   + f_0 * li_486[k];

        t_319[k] = -ii_179[k]
                   + f_0 * li_487[k];

        t_320[k] = -ii_180[k]
                   + f_0 * li_488[k];

        t_321[k] = -ii_181[k]
                   + f_0 * li_489[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, ii_182, ii_183, ii_184, ii_185, \
                         ii_186, li_490, li_491, li_492, li_493, \
                         li_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -ii_182[k]
                   + f_0 * li_490[k];

        t_323[k] = -ii_183[k]
                   + f_0 * li_491[k];

        t_324[k] = -ii_184[k]
                   + f_0 * li_492[k];

        t_325[k] = -ii_185[k]
                   + f_0 * li_493[k];

        t_326[k] = -ii_186[k]
                   + f_0 * li_494[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, ii_187, ii_188, ii_189, ii_190, \
                         ii_191, li_495, li_496, li_497, li_498, \
                         li_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -ii_187[k]
                   + f_0 * li_495[k];

        t_328[k] = -ii_188[k]
                   + f_0 * li_496[k];

        t_329[k] = -ii_189[k]
                   + f_0 * li_497[k];

        t_330[k] = -ii_190[k]
                   + f_0 * li_498[k];

        t_331[k] = -ii_191[k]
                   + f_0 * li_499[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, ii_192, ii_193, ii_194, ii_195, \
                         ii_196, li_500, li_501, li_502, li_503, \
                         li_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -ii_192[k]
                   + f_0 * li_500[k];

        t_333[k] = -ii_193[k]
                   + f_0 * li_501[k];

        t_334[k] = -ii_194[k]
                   + f_0 * li_502[k];

        t_335[k] = -ii_195[k]
                   + f_0 * li_503[k];

        t_336[k] = -2.0 * ii_196[k]
                   + f_0 * li_504[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, ii_197, ii_198, ii_199, ii_200, \
                         ii_201, li_505, li_506, li_507, li_508, \
                         li_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -2.0 * ii_197[k]
                   + f_0 * li_505[k];

        t_338[k] = -2.0 * ii_198[k]
                   + f_0 * li_506[k];

        t_339[k] = -2.0 * ii_199[k]
                   + f_0 * li_507[k];

        t_340[k] = -2.0 * ii_200[k]
                   + f_0 * li_508[k];

        t_341[k] = -2.0 * ii_201[k]
                   + f_0 * li_509[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, ii_202, ii_203, ii_204, ii_205, \
                         ii_206, li_510, li_511, li_512, li_513, \
                         li_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -2.0 * ii_202[k]
                   + f_0 * li_510[k];

        t_343[k] = -2.0 * ii_203[k]
                   + f_0 * li_511[k];

        t_344[k] = -2.0 * ii_204[k]
                   + f_0 * li_512[k];

        t_345[k] = -2.0 * ii_205[k]
                   + f_0 * li_513[k];

        t_346[k] = -2.0 * ii_206[k]
                   + f_0 * li_514[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_207 = buffer.data(ii + 207);
    const auto *ii_208 = buffer.data(ii + 208);
    const auto *ii_209 = buffer.data(ii + 209);
    const auto *ii_210 = buffer.data(ii + 210);
    const auto *ii_211 = buffer.data(ii + 211);
    const auto *ii_212 = buffer.data(ii + 212);
    const auto *ii_213 = buffer.data(ii + 213);
    const auto *ii_214 = buffer.data(ii + 214);
    const auto *ii_215 = buffer.data(ii + 215);
    const auto *ii_216 = buffer.data(ii + 216);
    const auto *ii_217 = buffer.data(ii + 217);
    const auto *ii_218 = buffer.data(ii + 218);
    const auto *ii_219 = buffer.data(ii + 219);
    const auto *ii_220 = buffer.data(ii + 220);
    const auto *ii_221 = buffer.data(ii + 221);
    const auto *ii_222 = buffer.data(ii + 222);
    const auto *ii_223 = buffer.data(ii + 223);
    const auto *ii_224 = buffer.data(ii + 224);
    const auto *ii_225 = buffer.data(ii + 225);
    const auto *ii_226 = buffer.data(ii + 226);
    const auto *ii_227 = buffer.data(ii + 227);
    const auto *ii_228 = buffer.data(ii + 228);
    const auto *ii_229 = buffer.data(ii + 229);
    const auto *ii_230 = buffer.data(ii + 230);
    const auto *ii_231 = buffer.data(ii + 231);
    const auto *ii_232 = buffer.data(ii + 232);
    const auto *ii_233 = buffer.data(ii + 233);
    const auto *ii_234 = buffer.data(ii + 234);
    const auto *ii_235 = buffer.data(ii + 235);
    const auto *ii_236 = buffer.data(ii + 236);
    const auto *ii_237 = buffer.data(ii + 237);
    const auto *ii_238 = buffer.data(ii + 238);
    const auto *ii_239 = buffer.data(ii + 239);
    const auto *ii_240 = buffer.data(ii + 240);
    const auto *ii_241 = buffer.data(ii + 241);
    const auto *ii_242 = buffer.data(ii + 242);
    const auto *ii_243 = buffer.data(ii + 243);
    const auto *ii_244 = buffer.data(ii + 244);
    const auto *ii_245 = buffer.data(ii + 245);
    const auto *ii_246 = buffer.data(ii + 246);
    const auto *ii_247 = buffer.data(ii + 247);
    const auto *ii_248 = buffer.data(ii + 248);
    const auto *ii_249 = buffer.data(ii + 249);
    const auto *ii_250 = buffer.data(ii + 250);
    const auto *ii_251 = buffer.data(ii + 251);
    const auto *ii_252 = buffer.data(ii + 252);
    const auto *ii_253 = buffer.data(ii + 253);
    const auto *ii_254 = buffer.data(ii + 254);
    const auto *ii_255 = buffer.data(ii + 255);
    const auto *ii_256 = buffer.data(ii + 256);
    const auto *ii_257 = buffer.data(ii + 257);
    const auto *ii_258 = buffer.data(ii + 258);
    const auto *ii_259 = buffer.data(ii + 259);
    const auto *ii_260 = buffer.data(ii + 260);
    const auto *ii_261 = buffer.data(ii + 261);
    const auto *ii_262 = buffer.data(ii + 262);
    const auto *ii_263 = buffer.data(ii + 263);
    const auto *ii_264 = buffer.data(ii + 264);
    const auto *ii_265 = buffer.data(ii + 265);
    const auto *ii_266 = buffer.data(ii + 266);
    const auto *ii_267 = buffer.data(ii + 267);
    const auto *ii_268 = buffer.data(ii + 268);
    const auto *ii_269 = buffer.data(ii + 269);
    const auto *ii_270 = buffer.data(ii + 270);
    const auto *ii_271 = buffer.data(ii + 271);
    const auto *ii_272 = buffer.data(ii + 272);
    const auto *ii_273 = buffer.data(ii + 273);
    const auto *ii_274 = buffer.data(ii + 274);
    const auto *ii_275 = buffer.data(ii + 275);
    const auto *ii_276 = buffer.data(ii + 276);
    const auto *ii_277 = buffer.data(ii + 277);
    const auto *ii_278 = buffer.data(ii + 278);
    const auto *ii_279 = buffer.data(ii + 279);
    const auto *ii_280 = buffer.data(ii + 280);
    const auto *ii_281 = buffer.data(ii + 281);
    const auto *ii_282 = buffer.data(ii + 282);
    const auto *ii_283 = buffer.data(ii + 283);
    const auto *ii_284 = buffer.data(ii + 284);
    const auto *ii_285 = buffer.data(ii + 285);
    const auto *ii_286 = buffer.data(ii + 286);
    const auto *ii_287 = buffer.data(ii + 287);
    const auto *ii_288 = buffer.data(ii + 288);
    const auto *ii_289 = buffer.data(ii + 289);
    const auto *ii_290 = buffer.data(ii + 290);
    const auto *ii_291 = buffer.data(ii + 291);
    const auto *ii_292 = buffer.data(ii + 292);
    const auto *ii_293 = buffer.data(ii + 293);
    const auto *ii_294 = buffer.data(ii + 294);
    const auto *ii_295 = buffer.data(ii + 295);
    const auto *ii_296 = buffer.data(ii + 296);
    const auto *ii_297 = buffer.data(ii + 297);
    const auto *ii_298 = buffer.data(ii + 298);
    const auto *ii_299 = buffer.data(ii + 299);
    const auto *ii_300 = buffer.data(ii + 300);
    const auto *ii_301 = buffer.data(ii + 301);
    const auto *ii_302 = buffer.data(ii + 302);
    const auto *ii_303 = buffer.data(ii + 303);
    const auto *ii_304 = buffer.data(ii + 304);
    const auto *ii_305 = buffer.data(ii + 305);
    const auto *ii_306 = buffer.data(ii + 306);
    const auto *ii_307 = buffer.data(ii + 307);
    const auto *ii_308 = buffer.data(ii + 308);
    const auto *ii_309 = buffer.data(ii + 309);
    const auto *ii_310 = buffer.data(ii + 310);
    const auto *ii_311 = buffer.data(ii + 311);
    const auto *ii_312 = buffer.data(ii + 312);
    const auto *ii_313 = buffer.data(ii + 313);
    const auto *ii_314 = buffer.data(ii + 314);
    const auto *ii_315 = buffer.data(ii + 315);
    const auto *ii_316 = buffer.data(ii + 316);
    const auto *ii_317 = buffer.data(ii + 317);
    const auto *ii_318 = buffer.data(ii + 318);
    const auto *ii_319 = buffer.data(ii + 319);
    const auto *ii_320 = buffer.data(ii + 320);
    const auto *ii_321 = buffer.data(ii + 321);
    const auto *ii_322 = buffer.data(ii + 322);
    const auto *ii_323 = buffer.data(ii + 323);
    const auto *ii_324 = buffer.data(ii + 324);
    const auto *ii_325 = buffer.data(ii + 325);
    const auto *ii_326 = buffer.data(ii + 326);
    const auto *ii_327 = buffer.data(ii + 327);
    const auto *ii_328 = buffer.data(ii + 328);
    const auto *ii_329 = buffer.data(ii + 329);
    const auto *ii_330 = buffer.data(ii + 330);
    const auto *ii_331 = buffer.data(ii + 331);
    const auto *ii_332 = buffer.data(ii + 332);
    const auto *ii_333 = buffer.data(ii + 333);
    const auto *ii_334 = buffer.data(ii + 334);
    const auto *ii_335 = buffer.data(ii + 335);
    const auto *ii_336 = buffer.data(ii + 336);
    const auto *ii_337 = buffer.data(ii + 337);
    const auto *ii_338 = buffer.data(ii + 338);

    const auto *li_515 = buffer.data(li + 515);
    const auto *li_516 = buffer.data(li + 516);
    const auto *li_517 = buffer.data(li + 517);
    const auto *li_518 = buffer.data(li + 518);
    const auto *li_519 = buffer.data(li + 519);
    const auto *li_520 = buffer.data(li + 520);
    const auto *li_521 = buffer.data(li + 521);
    const auto *li_522 = buffer.data(li + 522);
    const auto *li_523 = buffer.data(li + 523);
    const auto *li_524 = buffer.data(li + 524);
    const auto *li_525 = buffer.data(li + 525);
    const auto *li_526 = buffer.data(li + 526);
    const auto *li_527 = buffer.data(li + 527);
    const auto *li_528 = buffer.data(li + 528);
    const auto *li_529 = buffer.data(li + 529);
    const auto *li_530 = buffer.data(li + 530);
    const auto *li_531 = buffer.data(li + 531);
    const auto *li_532 = buffer.data(li + 532);
    const auto *li_533 = buffer.data(li + 533);
    const auto *li_534 = buffer.data(li + 534);
    const auto *li_535 = buffer.data(li + 535);
    const auto *li_536 = buffer.data(li + 536);
    const auto *li_537 = buffer.data(li + 537);
    const auto *li_538 = buffer.data(li + 538);
    const auto *li_539 = buffer.data(li + 539);
    const auto *li_540 = buffer.data(li + 540);
    const auto *li_541 = buffer.data(li + 541);
    const auto *li_542 = buffer.data(li + 542);
    const auto *li_543 = buffer.data(li + 543);
    const auto *li_544 = buffer.data(li + 544);
    const auto *li_545 = buffer.data(li + 545);
    const auto *li_546 = buffer.data(li + 546);
    const auto *li_547 = buffer.data(li + 547);
    const auto *li_548 = buffer.data(li + 548);
    const auto *li_549 = buffer.data(li + 549);
    const auto *li_550 = buffer.data(li + 550);
    const auto *li_551 = buffer.data(li + 551);
    const auto *li_552 = buffer.data(li + 552);
    const auto *li_553 = buffer.data(li + 553);
    const auto *li_554 = buffer.data(li + 554);
    const auto *li_555 = buffer.data(li + 555);
    const auto *li_556 = buffer.data(li + 556);
    const auto *li_557 = buffer.data(li + 557);
    const auto *li_558 = buffer.data(li + 558);
    const auto *li_559 = buffer.data(li + 559);
    const auto *li_560 = buffer.data(li + 560);
    const auto *li_561 = buffer.data(li + 561);
    const auto *li_562 = buffer.data(li + 562);
    const auto *li_563 = buffer.data(li + 563);
    const auto *li_564 = buffer.data(li + 564);
    const auto *li_565 = buffer.data(li + 565);
    const auto *li_566 = buffer.data(li + 566);
    const auto *li_567 = buffer.data(li + 567);
    const auto *li_568 = buffer.data(li + 568);
    const auto *li_569 = buffer.data(li + 569);
    const auto *li_570 = buffer.data(li + 570);
    const auto *li_571 = buffer.data(li + 571);
    const auto *li_572 = buffer.data(li + 572);
    const auto *li_573 = buffer.data(li + 573);
    const auto *li_574 = buffer.data(li + 574);
    const auto *li_575 = buffer.data(li + 575);
    const auto *li_576 = buffer.data(li + 576);
    const auto *li_577 = buffer.data(li + 577);
    const auto *li_578 = buffer.data(li + 578);
    const auto *li_579 = buffer.data(li + 579);
    const auto *li_580 = buffer.data(li + 580);
    const auto *li_581 = buffer.data(li + 581);
    const auto *li_582 = buffer.data(li + 582);
    const auto *li_583 = buffer.data(li + 583);
    const auto *li_584 = buffer.data(li + 584);
    const auto *li_585 = buffer.data(li + 585);
    const auto *li_586 = buffer.data(li + 586);
    const auto *li_587 = buffer.data(li + 587);
    const auto *li_616 = buffer.data(li + 616);
    const auto *li_617 = buffer.data(li + 617);
    const auto *li_618 = buffer.data(li + 618);
    const auto *li_619 = buffer.data(li + 619);
    const auto *li_620 = buffer.data(li + 620);
    const auto *li_621 = buffer.data(li + 621);
    const auto *li_622 = buffer.data(li + 622);
    const auto *li_623 = buffer.data(li + 623);
    const auto *li_624 = buffer.data(li + 624);
    const auto *li_625 = buffer.data(li + 625);
    const auto *li_626 = buffer.data(li + 626);
    const auto *li_627 = buffer.data(li + 627);
    const auto *li_628 = buffer.data(li + 628);
    const auto *li_629 = buffer.data(li + 629);
    const auto *li_630 = buffer.data(li + 630);
    const auto *li_631 = buffer.data(li + 631);
    const auto *li_632 = buffer.data(li + 632);
    const auto *li_633 = buffer.data(li + 633);
    const auto *li_634 = buffer.data(li + 634);
    const auto *li_635 = buffer.data(li + 635);
    const auto *li_636 = buffer.data(li + 636);
    const auto *li_637 = buffer.data(li + 637);
    const auto *li_638 = buffer.data(li + 638);
    const auto *li_639 = buffer.data(li + 639);
    const auto *li_640 = buffer.data(li + 640);
    const auto *li_641 = buffer.data(li + 641);
    const auto *li_642 = buffer.data(li + 642);
    const auto *li_643 = buffer.data(li + 643);
    const auto *li_644 = buffer.data(li + 644);
    const auto *li_645 = buffer.data(li + 645);
    const auto *li_646 = buffer.data(li + 646);
    const auto *li_647 = buffer.data(li + 647);
    const auto *li_648 = buffer.data(li + 648);
    const auto *li_649 = buffer.data(li + 649);
    const auto *li_650 = buffer.data(li + 650);
    const auto *li_651 = buffer.data(li + 651);
    const auto *li_652 = buffer.data(li + 652);
    const auto *li_653 = buffer.data(li + 653);
    const auto *li_654 = buffer.data(li + 654);
    const auto *li_655 = buffer.data(li + 655);
    const auto *li_656 = buffer.data(li + 656);
    const auto *li_657 = buffer.data(li + 657);
    const auto *li_658 = buffer.data(li + 658);
    const auto *li_659 = buffer.data(li + 659);
    const auto *li_660 = buffer.data(li + 660);
    const auto *li_661 = buffer.data(li + 661);
    const auto *li_662 = buffer.data(li + 662);
    const auto *li_663 = buffer.data(li + 663);
    const auto *li_664 = buffer.data(li + 664);
    const auto *li_665 = buffer.data(li + 665);
    const auto *li_666 = buffer.data(li + 666);
    const auto *li_667 = buffer.data(li + 667);
    const auto *li_668 = buffer.data(li + 668);
    const auto *li_669 = buffer.data(li + 669);
    const auto *li_670 = buffer.data(li + 670);
    const auto *li_671 = buffer.data(li + 671);
    const auto *li_672 = buffer.data(li + 672);
    const auto *li_673 = buffer.data(li + 673);
    const auto *li_674 = buffer.data(li + 674);
    const auto *li_675 = buffer.data(li + 675);
    const auto *li_676 = buffer.data(li + 676);
    const auto *li_677 = buffer.data(li + 677);
    const auto *li_678 = buffer.data(li + 678);
    const auto *li_679 = buffer.data(li + 679);
    const auto *li_680 = buffer.data(li + 680);
    const auto *li_681 = buffer.data(li + 681);
    const auto *li_682 = buffer.data(li + 682);
    const auto *li_683 = buffer.data(li + 683);
    const auto *li_684 = buffer.data(li + 684);
    const auto *li_685 = buffer.data(li + 685);
    const auto *li_686 = buffer.data(li + 686);
    const auto *li_687 = buffer.data(li + 687);
    const auto *li_688 = buffer.data(li + 688);
    const auto *li_689 = buffer.data(li + 689);
    const auto *li_690 = buffer.data(li + 690);
    const auto *li_691 = buffer.data(li + 691);
    const auto *li_692 = buffer.data(li + 692);
    const auto *li_693 = buffer.data(li + 693);
    const auto *li_694 = buffer.data(li + 694);
    const auto *li_695 = buffer.data(li + 695);
    const auto *li_696 = buffer.data(li + 696);
    const auto *li_697 = buffer.data(li + 697);
    const auto *li_698 = buffer.data(li + 698);
    const auto *li_699 = buffer.data(li + 699);
    const auto *li_700 = buffer.data(li + 700);
    const auto *li_701 = buffer.data(li + 701);
    const auto *li_702 = buffer.data(li + 702);

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, ii_207, ii_208, ii_209, ii_210, \
                         ii_211, li_515, li_516, li_517, li_518, \
                         li_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -2.0 * ii_207[k]
                   + f_0 * li_515[k];

        t_348[k] = -2.0 * ii_208[k]
                   + f_0 * li_516[k];

        t_349[k] = -2.0 * ii_209[k]
                   + f_0 * li_517[k];

        t_350[k] = -2.0 * ii_210[k]
                   + f_0 * li_518[k];

        t_351[k] = -2.0 * ii_211[k]
                   + f_0 * li_519[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, ii_212, ii_213, ii_214, ii_215, \
                         ii_216, li_520, li_521, li_522, li_523, \
                         li_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -2.0 * ii_212[k]
                   + f_0 * li_520[k];

        t_353[k] = -2.0 * ii_213[k]
                   + f_0 * li_521[k];

        t_354[k] = -2.0 * ii_214[k]
                   + f_0 * li_522[k];

        t_355[k] = -2.0 * ii_215[k]
                   + f_0 * li_523[k];

        t_356[k] = -2.0 * ii_216[k]
                   + f_0 * li_524[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, ii_217, ii_218, ii_219, ii_220, \
                         ii_221, li_525, li_526, li_527, li_528, \
                         li_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -2.0 * ii_217[k]
                   + f_0 * li_525[k];

        t_358[k] = -2.0 * ii_218[k]
                   + f_0 * li_526[k];

        t_359[k] = -2.0 * ii_219[k]
                   + f_0 * li_527[k];

        t_360[k] = -2.0 * ii_220[k]
                   + f_0 * li_528[k];

        t_361[k] = -2.0 * ii_221[k]
                   + f_0 * li_529[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, ii_222, ii_223, ii_224, ii_225, \
                         ii_226, li_530, li_531, li_532, li_533, \
                         li_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -2.0 * ii_222[k]
                   + f_0 * li_530[k];

        t_363[k] = -2.0 * ii_223[k]
                   + f_0 * li_531[k];

        t_364[k] = -3.0 * ii_224[k]
                   + f_0 * li_532[k];

        t_365[k] = -3.0 * ii_225[k]
                   + f_0 * li_533[k];

        t_366[k] = -3.0 * ii_226[k]
                   + f_0 * li_534[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, ii_227, ii_228, ii_229, ii_230, \
                         ii_231, li_535, li_536, li_537, li_538, \
                         li_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -3.0 * ii_227[k]
                   + f_0 * li_535[k];

        t_368[k] = -3.0 * ii_228[k]
                   + f_0 * li_536[k];

        t_369[k] = -3.0 * ii_229[k]
                   + f_0 * li_537[k];

        t_370[k] = -3.0 * ii_230[k]
                   + f_0 * li_538[k];

        t_371[k] = -3.0 * ii_231[k]
                   + f_0 * li_539[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, ii_232, ii_233, ii_234, ii_235, \
                         ii_236, li_540, li_541, li_542, li_543, \
                         li_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -3.0 * ii_232[k]
                   + f_0 * li_540[k];

        t_373[k] = -3.0 * ii_233[k]
                   + f_0 * li_541[k];

        t_374[k] = -3.0 * ii_234[k]
                   + f_0 * li_542[k];

        t_375[k] = -3.0 * ii_235[k]
                   + f_0 * li_543[k];

        t_376[k] = -3.0 * ii_236[k]
                   + f_0 * li_544[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, ii_237, ii_238, ii_239, ii_240, \
                         ii_241, li_545, li_546, li_547, li_548, \
                         li_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -3.0 * ii_237[k]
                   + f_0 * li_545[k];

        t_378[k] = -3.0 * ii_238[k]
                   + f_0 * li_546[k];

        t_379[k] = -3.0 * ii_239[k]
                   + f_0 * li_547[k];

        t_380[k] = -3.0 * ii_240[k]
                   + f_0 * li_548[k];

        t_381[k] = -3.0 * ii_241[k]
                   + f_0 * li_549[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, ii_242, ii_243, ii_244, ii_245, \
                         ii_246, li_550, li_551, li_552, li_553, \
                         li_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -3.0 * ii_242[k]
                   + f_0 * li_550[k];

        t_383[k] = -3.0 * ii_243[k]
                   + f_0 * li_551[k];

        t_384[k] = -3.0 * ii_244[k]
                   + f_0 * li_552[k];

        t_385[k] = -3.0 * ii_245[k]
                   + f_0 * li_553[k];

        t_386[k] = -3.0 * ii_246[k]
                   + f_0 * li_554[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, ii_247, ii_248, ii_249, ii_250, \
                         ii_251, li_555, li_556, li_557, li_558, \
                         li_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -3.0 * ii_247[k]
                   + f_0 * li_555[k];

        t_388[k] = -3.0 * ii_248[k]
                   + f_0 * li_556[k];

        t_389[k] = -3.0 * ii_249[k]
                   + f_0 * li_557[k];

        t_390[k] = -3.0 * ii_250[k]
                   + f_0 * li_558[k];

        t_391[k] = -3.0 * ii_251[k]
                   + f_0 * li_559[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, ii_252, ii_253, ii_254, ii_255, \
                         ii_256, li_560, li_561, li_562, li_563, \
                         li_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -4.0 * ii_252[k]
                   + f_0 * li_560[k];

        t_393[k] = -4.0 * ii_253[k]
                   + f_0 * li_561[k];

        t_394[k] = -4.0 * ii_254[k]
                   + f_0 * li_562[k];

        t_395[k] = -4.0 * ii_255[k]
                   + f_0 * li_563[k];

        t_396[k] = -4.0 * ii_256[k]
                   + f_0 * li_564[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, ii_257, ii_258, ii_259, ii_260, \
                         ii_261, li_565, li_566, li_567, li_568, \
                         li_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -4.0 * ii_257[k]
                   + f_0 * li_565[k];

        t_398[k] = -4.0 * ii_258[k]
                   + f_0 * li_566[k];

        t_399[k] = -4.0 * ii_259[k]
                   + f_0 * li_567[k];

        t_400[k] = -4.0 * ii_260[k]
                   + f_0 * li_568[k];

        t_401[k] = -4.0 * ii_261[k]
                   + f_0 * li_569[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, ii_262, ii_263, ii_264, ii_265, \
                         ii_266, li_570, li_571, li_572, li_573, \
                         li_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -4.0 * ii_262[k]
                   + f_0 * li_570[k];

        t_403[k] = -4.0 * ii_263[k]
                   + f_0 * li_571[k];

        t_404[k] = -4.0 * ii_264[k]
                   + f_0 * li_572[k];

        t_405[k] = -4.0 * ii_265[k]
                   + f_0 * li_573[k];

        t_406[k] = -4.0 * ii_266[k]
                   + f_0 * li_574[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, ii_267, ii_268, ii_269, ii_270, \
                         ii_271, li_575, li_576, li_577, li_578, \
                         li_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -4.0 * ii_267[k]
                   + f_0 * li_575[k];

        t_408[k] = -4.0 * ii_268[k]
                   + f_0 * li_576[k];

        t_409[k] = -4.0 * ii_269[k]
                   + f_0 * li_577[k];

        t_410[k] = -4.0 * ii_270[k]
                   + f_0 * li_578[k];

        t_411[k] = -4.0 * ii_271[k]
                   + f_0 * li_579[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, ii_272, ii_273, ii_274, ii_275, \
                         ii_276, li_580, li_581, li_582, li_583, \
                         li_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -4.0 * ii_272[k]
                   + f_0 * li_580[k];

        t_413[k] = -4.0 * ii_273[k]
                   + f_0 * li_581[k];

        t_414[k] = -4.0 * ii_274[k]
                   + f_0 * li_582[k];

        t_415[k] = -4.0 * ii_275[k]
                   + f_0 * li_583[k];

        t_416[k] = -4.0 * ii_276[k]
                   + f_0 * li_584[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, t_422, ii_277, ii_278, ii_279, \
                         li_585, li_586, li_587, li_616, li_617, \
                         li_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -4.0 * ii_277[k]
                   + f_0 * li_585[k];

        t_418[k] = -4.0 * ii_278[k]
                   + f_0 * li_586[k];

        t_419[k] = -4.0 * ii_279[k]
                   + f_0 * li_587[k];

        t_420[k] = f_0 * li_616[k];

        t_421[k] = f_0 * li_617[k];

        t_422[k] = f_0 * li_618[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, t_428, t_429, t_430, li_619, \
                         li_620, li_621, li_622, li_623, li_624, li_625, \
                         li_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_0 * li_619[k];

        t_424[k] = f_0 * li_620[k];

        t_425[k] = f_0 * li_621[k];

        t_426[k] = f_0 * li_622[k];

        t_427[k] = f_0 * li_623[k];

        t_428[k] = f_0 * li_624[k];

        t_429[k] = f_0 * li_625[k];

        t_430[k] = f_0 * li_626[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, t_436, t_437, t_438, li_627, \
                         li_628, li_629, li_630, li_631, li_632, li_633, \
                         li_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_0 * li_627[k];

        t_432[k] = f_0 * li_628[k];

        t_433[k] = f_0 * li_629[k];

        t_434[k] = f_0 * li_630[k];

        t_435[k] = f_0 * li_631[k];

        t_436[k] = f_0 * li_632[k];

        t_437[k] = f_0 * li_633[k];

        t_438[k] = f_0 * li_634[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, t_444, t_445, t_446, li_635, \
                         li_636, li_637, li_638, li_639, li_640, li_641, \
                         li_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_0 * li_635[k];

        t_440[k] = f_0 * li_636[k];

        t_441[k] = f_0 * li_637[k];

        t_442[k] = f_0 * li_638[k];

        t_443[k] = f_0 * li_639[k];

        t_444[k] = f_0 * li_640[k];

        t_445[k] = f_0 * li_641[k];

        t_446[k] = f_0 * li_642[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, ii_280, ii_281, ii_282, ii_283, \
                         li_643, li_644, li_645, li_646, li_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_0 * li_643[k];

        t_448[k] = -ii_280[k]
                   + f_0 * li_644[k];

        t_449[k] = -ii_281[k]
                   + f_0 * li_645[k];

        t_450[k] = -ii_282[k]
                   + f_0 * li_646[k];

        t_451[k] = -ii_283[k]
                   + f_0 * li_647[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, ii_284, ii_285, ii_286, ii_287, \
                         ii_288, li_648, li_649, li_650, li_651, \
                         li_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = -ii_284[k]
                   + f_0 * li_648[k];

        t_453[k] = -ii_285[k]
                   + f_0 * li_649[k];

        t_454[k] = -ii_286[k]
                   + f_0 * li_650[k];

        t_455[k] = -ii_287[k]
                   + f_0 * li_651[k];

        t_456[k] = -ii_288[k]
                   + f_0 * li_652[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, ii_289, ii_290, ii_291, ii_292, \
                         ii_293, li_653, li_654, li_655, li_656, \
                         li_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = -ii_289[k]
                   + f_0 * li_653[k];

        t_458[k] = -ii_290[k]
                   + f_0 * li_654[k];

        t_459[k] = -ii_291[k]
                   + f_0 * li_655[k];

        t_460[k] = -ii_292[k]
                   + f_0 * li_656[k];

        t_461[k] = -ii_293[k]
                   + f_0 * li_657[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, ii_294, ii_295, ii_296, ii_297, \
                         ii_298, li_658, li_659, li_660, li_661, \
                         li_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -ii_294[k]
                   + f_0 * li_658[k];

        t_463[k] = -ii_295[k]
                   + f_0 * li_659[k];

        t_464[k] = -ii_296[k]
                   + f_0 * li_660[k];

        t_465[k] = -ii_297[k]
                   + f_0 * li_661[k];

        t_466[k] = -ii_298[k]
                   + f_0 * li_662[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, ii_299, ii_300, ii_301, ii_302, \
                         ii_303, li_663, li_664, li_665, li_666, \
                         li_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = -ii_299[k]
                   + f_0 * li_663[k];

        t_468[k] = -ii_300[k]
                   + f_0 * li_664[k];

        t_469[k] = -ii_301[k]
                   + f_0 * li_665[k];

        t_470[k] = -ii_302[k]
                   + f_0 * li_666[k];

        t_471[k] = -ii_303[k]
                   + f_0 * li_667[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, ii_304, ii_305, ii_306, ii_307, \
                         ii_308, li_668, li_669, li_670, li_671, \
                         li_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = -ii_304[k]
                   + f_0 * li_668[k];

        t_473[k] = -ii_305[k]
                   + f_0 * li_669[k];

        t_474[k] = -ii_306[k]
                   + f_0 * li_670[k];

        t_475[k] = -ii_307[k]
                   + f_0 * li_671[k];

        t_476[k] = -2.0 * ii_308[k]
                   + f_0 * li_672[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, ii_309, ii_310, ii_311, ii_312, \
                         ii_313, li_673, li_674, li_675, li_676, \
                         li_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = -2.0 * ii_309[k]
                   + f_0 * li_673[k];

        t_478[k] = -2.0 * ii_310[k]
                   + f_0 * li_674[k];

        t_479[k] = -2.0 * ii_311[k]
                   + f_0 * li_675[k];

        t_480[k] = -2.0 * ii_312[k]
                   + f_0 * li_676[k];

        t_481[k] = -2.0 * ii_313[k]
                   + f_0 * li_677[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, ii_314, ii_315, ii_316, ii_317, \
                         ii_318, li_678, li_679, li_680, li_681, \
                         li_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -2.0 * ii_314[k]
                   + f_0 * li_678[k];

        t_483[k] = -2.0 * ii_315[k]
                   + f_0 * li_679[k];

        t_484[k] = -2.0 * ii_316[k]
                   + f_0 * li_680[k];

        t_485[k] = -2.0 * ii_317[k]
                   + f_0 * li_681[k];

        t_486[k] = -2.0 * ii_318[k]
                   + f_0 * li_682[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, ii_319, ii_320, ii_321, ii_322, \
                         ii_323, li_683, li_684, li_685, li_686, \
                         li_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = -2.0 * ii_319[k]
                   + f_0 * li_683[k];

        t_488[k] = -2.0 * ii_320[k]
                   + f_0 * li_684[k];

        t_489[k] = -2.0 * ii_321[k]
                   + f_0 * li_685[k];

        t_490[k] = -2.0 * ii_322[k]
                   + f_0 * li_686[k];

        t_491[k] = -2.0 * ii_323[k]
                   + f_0 * li_687[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, ii_324, ii_325, ii_326, ii_327, \
                         ii_328, li_688, li_689, li_690, li_691, \
                         li_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = -2.0 * ii_324[k]
                   + f_0 * li_688[k];

        t_493[k] = -2.0 * ii_325[k]
                   + f_0 * li_689[k];

        t_494[k] = -2.0 * ii_326[k]
                   + f_0 * li_690[k];

        t_495[k] = -2.0 * ii_327[k]
                   + f_0 * li_691[k];

        t_496[k] = -2.0 * ii_328[k]
                   + f_0 * li_692[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, t_501, ii_329, ii_330, ii_331, ii_332, \
                         ii_333, li_693, li_694, li_695, li_696, \
                         li_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = -2.0 * ii_329[k]
                   + f_0 * li_693[k];

        t_498[k] = -2.0 * ii_330[k]
                   + f_0 * li_694[k];

        t_499[k] = -2.0 * ii_331[k]
                   + f_0 * li_695[k];

        t_500[k] = -2.0 * ii_332[k]
                   + f_0 * li_696[k];

        t_501[k] = -2.0 * ii_333[k]
                   + f_0 * li_697[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, ii_334, ii_335, ii_336, ii_337, \
                         ii_338, li_698, li_699, li_700, li_701, \
                         li_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = -2.0 * ii_334[k]
                   + f_0 * li_698[k];

        t_503[k] = -2.0 * ii_335[k]
                   + f_0 * li_699[k];

        t_504[k] = -3.0 * ii_336[k]
                   + f_0 * li_700[k];

        t_505[k] = -3.0 * ii_337[k]
                   + f_0 * li_701[k];

        t_506[k] = -3.0 * ii_338[k]
                   + f_0 * li_702[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_339 = buffer.data(ii + 339);
    const auto *ii_340 = buffer.data(ii + 340);
    const auto *ii_341 = buffer.data(ii + 341);
    const auto *ii_342 = buffer.data(ii + 342);
    const auto *ii_343 = buffer.data(ii + 343);
    const auto *ii_344 = buffer.data(ii + 344);
    const auto *ii_345 = buffer.data(ii + 345);
    const auto *ii_346 = buffer.data(ii + 346);
    const auto *ii_347 = buffer.data(ii + 347);
    const auto *ii_348 = buffer.data(ii + 348);
    const auto *ii_349 = buffer.data(ii + 349);
    const auto *ii_350 = buffer.data(ii + 350);
    const auto *ii_351 = buffer.data(ii + 351);
    const auto *ii_352 = buffer.data(ii + 352);
    const auto *ii_353 = buffer.data(ii + 353);
    const auto *ii_354 = buffer.data(ii + 354);
    const auto *ii_355 = buffer.data(ii + 355);
    const auto *ii_356 = buffer.data(ii + 356);
    const auto *ii_357 = buffer.data(ii + 357);
    const auto *ii_358 = buffer.data(ii + 358);
    const auto *ii_359 = buffer.data(ii + 359);
    const auto *ii_360 = buffer.data(ii + 360);
    const auto *ii_361 = buffer.data(ii + 361);
    const auto *ii_362 = buffer.data(ii + 362);
    const auto *ii_363 = buffer.data(ii + 363);
    const auto *ii_364 = buffer.data(ii + 364);
    const auto *ii_365 = buffer.data(ii + 365);
    const auto *ii_366 = buffer.data(ii + 366);
    const auto *ii_367 = buffer.data(ii + 367);
    const auto *ii_368 = buffer.data(ii + 368);
    const auto *ii_369 = buffer.data(ii + 369);
    const auto *ii_370 = buffer.data(ii + 370);
    const auto *ii_371 = buffer.data(ii + 371);
    const auto *ii_372 = buffer.data(ii + 372);
    const auto *ii_373 = buffer.data(ii + 373);
    const auto *ii_374 = buffer.data(ii + 374);
    const auto *ii_375 = buffer.data(ii + 375);
    const auto *ii_376 = buffer.data(ii + 376);
    const auto *ii_377 = buffer.data(ii + 377);
    const auto *ii_378 = buffer.data(ii + 378);
    const auto *ii_379 = buffer.data(ii + 379);
    const auto *ii_380 = buffer.data(ii + 380);
    const auto *ii_381 = buffer.data(ii + 381);
    const auto *ii_382 = buffer.data(ii + 382);
    const auto *ii_383 = buffer.data(ii + 383);
    const auto *ii_384 = buffer.data(ii + 384);
    const auto *ii_385 = buffer.data(ii + 385);
    const auto *ii_386 = buffer.data(ii + 386);
    const auto *ii_387 = buffer.data(ii + 387);
    const auto *ii_388 = buffer.data(ii + 388);
    const auto *ii_389 = buffer.data(ii + 389);
    const auto *ii_390 = buffer.data(ii + 390);
    const auto *ii_391 = buffer.data(ii + 391);
    const auto *ii_392 = buffer.data(ii + 392);
    const auto *ii_393 = buffer.data(ii + 393);
    const auto *ii_394 = buffer.data(ii + 394);
    const auto *ii_395 = buffer.data(ii + 395);
    const auto *ii_396 = buffer.data(ii + 396);
    const auto *ii_397 = buffer.data(ii + 397);
    const auto *ii_398 = buffer.data(ii + 398);
    const auto *ii_399 = buffer.data(ii + 399);
    const auto *ii_400 = buffer.data(ii + 400);
    const auto *ii_401 = buffer.data(ii + 401);
    const auto *ii_402 = buffer.data(ii + 402);
    const auto *ii_403 = buffer.data(ii + 403);
    const auto *ii_404 = buffer.data(ii + 404);
    const auto *ii_405 = buffer.data(ii + 405);
    const auto *ii_406 = buffer.data(ii + 406);
    const auto *ii_407 = buffer.data(ii + 407);
    const auto *ii_408 = buffer.data(ii + 408);
    const auto *ii_409 = buffer.data(ii + 409);
    const auto *ii_410 = buffer.data(ii + 410);
    const auto *ii_411 = buffer.data(ii + 411);
    const auto *ii_412 = buffer.data(ii + 412);
    const auto *ii_413 = buffer.data(ii + 413);
    const auto *ii_414 = buffer.data(ii + 414);
    const auto *ii_415 = buffer.data(ii + 415);
    const auto *ii_416 = buffer.data(ii + 416);
    const auto *ii_417 = buffer.data(ii + 417);
    const auto *ii_418 = buffer.data(ii + 418);
    const auto *ii_419 = buffer.data(ii + 419);
    const auto *ii_420 = buffer.data(ii + 420);
    const auto *ii_421 = buffer.data(ii + 421);
    const auto *ii_422 = buffer.data(ii + 422);
    const auto *ii_423 = buffer.data(ii + 423);
    const auto *ii_424 = buffer.data(ii + 424);
    const auto *ii_425 = buffer.data(ii + 425);
    const auto *ii_426 = buffer.data(ii + 426);
    const auto *ii_427 = buffer.data(ii + 427);
    const auto *ii_428 = buffer.data(ii + 428);
    const auto *ii_429 = buffer.data(ii + 429);
    const auto *ii_430 = buffer.data(ii + 430);
    const auto *ii_431 = buffer.data(ii + 431);
    const auto *ii_432 = buffer.data(ii + 432);
    const auto *ii_433 = buffer.data(ii + 433);
    const auto *ii_434 = buffer.data(ii + 434);
    const auto *ii_435 = buffer.data(ii + 435);
    const auto *ii_436 = buffer.data(ii + 436);
    const auto *ii_437 = buffer.data(ii + 437);
    const auto *ii_438 = buffer.data(ii + 438);
    const auto *ii_439 = buffer.data(ii + 439);
    const auto *ii_440 = buffer.data(ii + 440);
    const auto *ii_441 = buffer.data(ii + 441);
    const auto *ii_442 = buffer.data(ii + 442);
    const auto *ii_443 = buffer.data(ii + 443);
    const auto *ii_444 = buffer.data(ii + 444);
    const auto *ii_445 = buffer.data(ii + 445);
    const auto *ii_446 = buffer.data(ii + 446);
    const auto *ii_447 = buffer.data(ii + 447);
    const auto *ii_448 = buffer.data(ii + 448);
    const auto *ii_449 = buffer.data(ii + 449);
    const auto *ii_450 = buffer.data(ii + 450);
    const auto *ii_451 = buffer.data(ii + 451);
    const auto *ii_452 = buffer.data(ii + 452);
    const auto *ii_453 = buffer.data(ii + 453);
    const auto *ii_454 = buffer.data(ii + 454);
    const auto *ii_455 = buffer.data(ii + 455);
    const auto *ii_456 = buffer.data(ii + 456);
    const auto *ii_457 = buffer.data(ii + 457);
    const auto *ii_458 = buffer.data(ii + 458);
    const auto *ii_459 = buffer.data(ii + 459);
    const auto *ii_460 = buffer.data(ii + 460);
    const auto *ii_461 = buffer.data(ii + 461);
    const auto *ii_462 = buffer.data(ii + 462);
    const auto *ii_463 = buffer.data(ii + 463);
    const auto *ii_464 = buffer.data(ii + 464);
    const auto *ii_465 = buffer.data(ii + 465);
    const auto *ii_466 = buffer.data(ii + 466);
    const auto *ii_467 = buffer.data(ii + 467);
    const auto *ii_468 = buffer.data(ii + 468);
    const auto *ii_469 = buffer.data(ii + 469);
    const auto *ii_470 = buffer.data(ii + 470);

    const auto *li_703 = buffer.data(li + 703);
    const auto *li_704 = buffer.data(li + 704);
    const auto *li_705 = buffer.data(li + 705);
    const auto *li_706 = buffer.data(li + 706);
    const auto *li_707 = buffer.data(li + 707);
    const auto *li_708 = buffer.data(li + 708);
    const auto *li_709 = buffer.data(li + 709);
    const auto *li_710 = buffer.data(li + 710);
    const auto *li_711 = buffer.data(li + 711);
    const auto *li_712 = buffer.data(li + 712);
    const auto *li_713 = buffer.data(li + 713);
    const auto *li_714 = buffer.data(li + 714);
    const auto *li_715 = buffer.data(li + 715);
    const auto *li_716 = buffer.data(li + 716);
    const auto *li_717 = buffer.data(li + 717);
    const auto *li_718 = buffer.data(li + 718);
    const auto *li_719 = buffer.data(li + 719);
    const auto *li_720 = buffer.data(li + 720);
    const auto *li_721 = buffer.data(li + 721);
    const auto *li_722 = buffer.data(li + 722);
    const auto *li_723 = buffer.data(li + 723);
    const auto *li_724 = buffer.data(li + 724);
    const auto *li_725 = buffer.data(li + 725);
    const auto *li_726 = buffer.data(li + 726);
    const auto *li_727 = buffer.data(li + 727);
    const auto *li_728 = buffer.data(li + 728);
    const auto *li_729 = buffer.data(li + 729);
    const auto *li_730 = buffer.data(li + 730);
    const auto *li_731 = buffer.data(li + 731);
    const auto *li_732 = buffer.data(li + 732);
    const auto *li_733 = buffer.data(li + 733);
    const auto *li_734 = buffer.data(li + 734);
    const auto *li_735 = buffer.data(li + 735);
    const auto *li_736 = buffer.data(li + 736);
    const auto *li_737 = buffer.data(li + 737);
    const auto *li_738 = buffer.data(li + 738);
    const auto *li_739 = buffer.data(li + 739);
    const auto *li_740 = buffer.data(li + 740);
    const auto *li_741 = buffer.data(li + 741);
    const auto *li_742 = buffer.data(li + 742);
    const auto *li_743 = buffer.data(li + 743);
    const auto *li_744 = buffer.data(li + 744);
    const auto *li_745 = buffer.data(li + 745);
    const auto *li_746 = buffer.data(li + 746);
    const auto *li_747 = buffer.data(li + 747);
    const auto *li_748 = buffer.data(li + 748);
    const auto *li_749 = buffer.data(li + 749);
    const auto *li_750 = buffer.data(li + 750);
    const auto *li_751 = buffer.data(li + 751);
    const auto *li_752 = buffer.data(li + 752);
    const auto *li_753 = buffer.data(li + 753);
    const auto *li_754 = buffer.data(li + 754);
    const auto *li_755 = buffer.data(li + 755);
    const auto *li_756 = buffer.data(li + 756);
    const auto *li_757 = buffer.data(li + 757);
    const auto *li_758 = buffer.data(li + 758);
    const auto *li_759 = buffer.data(li + 759);
    const auto *li_760 = buffer.data(li + 760);
    const auto *li_761 = buffer.data(li + 761);
    const auto *li_762 = buffer.data(li + 762);
    const auto *li_763 = buffer.data(li + 763);
    const auto *li_764 = buffer.data(li + 764);
    const auto *li_765 = buffer.data(li + 765);
    const auto *li_766 = buffer.data(li + 766);
    const auto *li_767 = buffer.data(li + 767);
    const auto *li_768 = buffer.data(li + 768);
    const auto *li_769 = buffer.data(li + 769);
    const auto *li_770 = buffer.data(li + 770);
    const auto *li_771 = buffer.data(li + 771);
    const auto *li_772 = buffer.data(li + 772);
    const auto *li_773 = buffer.data(li + 773);
    const auto *li_774 = buffer.data(li + 774);
    const auto *li_775 = buffer.data(li + 775);
    const auto *li_776 = buffer.data(li + 776);
    const auto *li_777 = buffer.data(li + 777);
    const auto *li_778 = buffer.data(li + 778);
    const auto *li_779 = buffer.data(li + 779);
    const auto *li_780 = buffer.data(li + 780);
    const auto *li_781 = buffer.data(li + 781);
    const auto *li_782 = buffer.data(li + 782);
    const auto *li_783 = buffer.data(li + 783);
    const auto *li_812 = buffer.data(li + 812);
    const auto *li_813 = buffer.data(li + 813);
    const auto *li_814 = buffer.data(li + 814);
    const auto *li_815 = buffer.data(li + 815);
    const auto *li_816 = buffer.data(li + 816);
    const auto *li_817 = buffer.data(li + 817);
    const auto *li_818 = buffer.data(li + 818);
    const auto *li_819 = buffer.data(li + 819);
    const auto *li_820 = buffer.data(li + 820);
    const auto *li_821 = buffer.data(li + 821);
    const auto *li_822 = buffer.data(li + 822);
    const auto *li_823 = buffer.data(li + 823);
    const auto *li_824 = buffer.data(li + 824);
    const auto *li_825 = buffer.data(li + 825);
    const auto *li_826 = buffer.data(li + 826);
    const auto *li_827 = buffer.data(li + 827);
    const auto *li_828 = buffer.data(li + 828);
    const auto *li_829 = buffer.data(li + 829);
    const auto *li_830 = buffer.data(li + 830);
    const auto *li_831 = buffer.data(li + 831);
    const auto *li_832 = buffer.data(li + 832);
    const auto *li_833 = buffer.data(li + 833);
    const auto *li_834 = buffer.data(li + 834);
    const auto *li_835 = buffer.data(li + 835);
    const auto *li_836 = buffer.data(li + 836);
    const auto *li_837 = buffer.data(li + 837);
    const auto *li_838 = buffer.data(li + 838);
    const auto *li_839 = buffer.data(li + 839);
    const auto *li_840 = buffer.data(li + 840);
    const auto *li_841 = buffer.data(li + 841);
    const auto *li_842 = buffer.data(li + 842);
    const auto *li_843 = buffer.data(li + 843);
    const auto *li_844 = buffer.data(li + 844);
    const auto *li_845 = buffer.data(li + 845);
    const auto *li_846 = buffer.data(li + 846);
    const auto *li_847 = buffer.data(li + 847);
    const auto *li_848 = buffer.data(li + 848);
    const auto *li_849 = buffer.data(li + 849);
    const auto *li_850 = buffer.data(li + 850);
    const auto *li_851 = buffer.data(li + 851);
    const auto *li_852 = buffer.data(li + 852);
    const auto *li_853 = buffer.data(li + 853);
    const auto *li_854 = buffer.data(li + 854);
    const auto *li_855 = buffer.data(li + 855);
    const auto *li_856 = buffer.data(li + 856);
    const auto *li_857 = buffer.data(li + 857);
    const auto *li_858 = buffer.data(li + 858);
    const auto *li_859 = buffer.data(li + 859);
    const auto *li_860 = buffer.data(li + 860);
    const auto *li_861 = buffer.data(li + 861);
    const auto *li_862 = buffer.data(li + 862);
    const auto *li_863 = buffer.data(li + 863);
    const auto *li_864 = buffer.data(li + 864);
    const auto *li_865 = buffer.data(li + 865);
    const auto *li_866 = buffer.data(li + 866);
    const auto *li_867 = buffer.data(li + 867);
    const auto *li_868 = buffer.data(li + 868);
    const auto *li_869 = buffer.data(li + 869);
    const auto *li_870 = buffer.data(li + 870);
    const auto *li_871 = buffer.data(li + 871);
    const auto *li_872 = buffer.data(li + 872);
    const auto *li_873 = buffer.data(li + 873);
    const auto *li_874 = buffer.data(li + 874);
    const auto *li_875 = buffer.data(li + 875);
    const auto *li_876 = buffer.data(li + 876);
    const auto *li_877 = buffer.data(li + 877);
    const auto *li_878 = buffer.data(li + 878);
    const auto *li_879 = buffer.data(li + 879);
    const auto *li_880 = buffer.data(li + 880);
    const auto *li_881 = buffer.data(li + 881);
    const auto *li_882 = buffer.data(li + 882);
    const auto *li_883 = buffer.data(li + 883);
    const auto *li_884 = buffer.data(li + 884);
    const auto *li_885 = buffer.data(li + 885);
    const auto *li_886 = buffer.data(li + 886);
    const auto *li_887 = buffer.data(li + 887);
    const auto *li_888 = buffer.data(li + 888);
    const auto *li_889 = buffer.data(li + 889);
    const auto *li_890 = buffer.data(li + 890);

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, ii_339, ii_340, ii_341, ii_342, \
                         ii_343, li_703, li_704, li_705, li_706, \
                         li_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = -3.0 * ii_339[k]
                   + f_0 * li_703[k];

        t_508[k] = -3.0 * ii_340[k]
                   + f_0 * li_704[k];

        t_509[k] = -3.0 * ii_341[k]
                   + f_0 * li_705[k];

        t_510[k] = -3.0 * ii_342[k]
                   + f_0 * li_706[k];

        t_511[k] = -3.0 * ii_343[k]
                   + f_0 * li_707[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, ii_344, ii_345, ii_346, ii_347, \
                         ii_348, li_708, li_709, li_710, li_711, \
                         li_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = -3.0 * ii_344[k]
                   + f_0 * li_708[k];

        t_513[k] = -3.0 * ii_345[k]
                   + f_0 * li_709[k];

        t_514[k] = -3.0 * ii_346[k]
                   + f_0 * li_710[k];

        t_515[k] = -3.0 * ii_347[k]
                   + f_0 * li_711[k];

        t_516[k] = -3.0 * ii_348[k]
                   + f_0 * li_712[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, ii_349, ii_350, ii_351, ii_352, \
                         ii_353, li_713, li_714, li_715, li_716, \
                         li_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -3.0 * ii_349[k]
                   + f_0 * li_713[k];

        t_518[k] = -3.0 * ii_350[k]
                   + f_0 * li_714[k];

        t_519[k] = -3.0 * ii_351[k]
                   + f_0 * li_715[k];

        t_520[k] = -3.0 * ii_352[k]
                   + f_0 * li_716[k];

        t_521[k] = -3.0 * ii_353[k]
                   + f_0 * li_717[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, ii_354, ii_355, ii_356, ii_357, \
                         ii_358, li_718, li_719, li_720, li_721, \
                         li_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = -3.0 * ii_354[k]
                   + f_0 * li_718[k];

        t_523[k] = -3.0 * ii_355[k]
                   + f_0 * li_719[k];

        t_524[k] = -3.0 * ii_356[k]
                   + f_0 * li_720[k];

        t_525[k] = -3.0 * ii_357[k]
                   + f_0 * li_721[k];

        t_526[k] = -3.0 * ii_358[k]
                   + f_0 * li_722[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, ii_359, ii_360, ii_361, ii_362, \
                         ii_363, li_723, li_724, li_725, li_726, \
                         li_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = -3.0 * ii_359[k]
                   + f_0 * li_723[k];

        t_528[k] = -3.0 * ii_360[k]
                   + f_0 * li_724[k];

        t_529[k] = -3.0 * ii_361[k]
                   + f_0 * li_725[k];

        t_530[k] = -3.0 * ii_362[k]
                   + f_0 * li_726[k];

        t_531[k] = -3.0 * ii_363[k]
                   + f_0 * li_727[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, ii_364, ii_365, ii_366, ii_367, \
                         ii_368, li_728, li_729, li_730, li_731, \
                         li_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = -4.0 * ii_364[k]
                   + f_0 * li_728[k];

        t_533[k] = -4.0 * ii_365[k]
                   + f_0 * li_729[k];

        t_534[k] = -4.0 * ii_366[k]
                   + f_0 * li_730[k];

        t_535[k] = -4.0 * ii_367[k]
                   + f_0 * li_731[k];

        t_536[k] = -4.0 * ii_368[k]
                   + f_0 * li_732[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, t_541, ii_369, ii_370, ii_371, ii_372, \
                         ii_373, li_733, li_734, li_735, li_736, \
                         li_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = -4.0 * ii_369[k]
                   + f_0 * li_733[k];

        t_538[k] = -4.0 * ii_370[k]
                   + f_0 * li_734[k];

        t_539[k] = -4.0 * ii_371[k]
                   + f_0 * li_735[k];

        t_540[k] = -4.0 * ii_372[k]
                   + f_0 * li_736[k];

        t_541[k] = -4.0 * ii_373[k]
                   + f_0 * li_737[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, ii_374, ii_375, ii_376, ii_377, \
                         ii_378, li_738, li_739, li_740, li_741, \
                         li_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -4.0 * ii_374[k]
                   + f_0 * li_738[k];

        t_543[k] = -4.0 * ii_375[k]
                   + f_0 * li_739[k];

        t_544[k] = -4.0 * ii_376[k]
                   + f_0 * li_740[k];

        t_545[k] = -4.0 * ii_377[k]
                   + f_0 * li_741[k];

        t_546[k] = -4.0 * ii_378[k]
                   + f_0 * li_742[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, ii_379, ii_380, ii_381, ii_382, \
                         ii_383, li_743, li_744, li_745, li_746, \
                         li_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -4.0 * ii_379[k]
                   + f_0 * li_743[k];

        t_548[k] = -4.0 * ii_380[k]
                   + f_0 * li_744[k];

        t_549[k] = -4.0 * ii_381[k]
                   + f_0 * li_745[k];

        t_550[k] = -4.0 * ii_382[k]
                   + f_0 * li_746[k];

        t_551[k] = -4.0 * ii_383[k]
                   + f_0 * li_747[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, ii_384, ii_385, ii_386, ii_387, \
                         ii_388, li_748, li_749, li_750, li_751, \
                         li_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -4.0 * ii_384[k]
                   + f_0 * li_748[k];

        t_553[k] = -4.0 * ii_385[k]
                   + f_0 * li_749[k];

        t_554[k] = -4.0 * ii_386[k]
                   + f_0 * li_750[k];

        t_555[k] = -4.0 * ii_387[k]
                   + f_0 * li_751[k];

        t_556[k] = -4.0 * ii_388[k]
                   + f_0 * li_752[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, ii_389, ii_390, ii_391, ii_392, \
                         ii_393, li_753, li_754, li_755, li_756, \
                         li_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -4.0 * ii_389[k]
                   + f_0 * li_753[k];

        t_558[k] = -4.0 * ii_390[k]
                   + f_0 * li_754[k];

        t_559[k] = -4.0 * ii_391[k]
                   + f_0 * li_755[k];

        t_560[k] = -5.0 * ii_392[k]
                   + f_0 * li_756[k];

        t_561[k] = -5.0 * ii_393[k]
                   + f_0 * li_757[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, ii_394, ii_395, ii_396, ii_397, \
                         ii_398, li_758, li_759, li_760, li_761, \
                         li_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -5.0 * ii_394[k]
                   + f_0 * li_758[k];

        t_563[k] = -5.0 * ii_395[k]
                   + f_0 * li_759[k];

        t_564[k] = -5.0 * ii_396[k]
                   + f_0 * li_760[k];

        t_565[k] = -5.0 * ii_397[k]
                   + f_0 * li_761[k];

        t_566[k] = -5.0 * ii_398[k]
                   + f_0 * li_762[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, ii_399, ii_400, ii_401, ii_402, \
                         ii_403, li_763, li_764, li_765, li_766, \
                         li_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -5.0 * ii_399[k]
                   + f_0 * li_763[k];

        t_568[k] = -5.0 * ii_400[k]
                   + f_0 * li_764[k];

        t_569[k] = -5.0 * ii_401[k]
                   + f_0 * li_765[k];

        t_570[k] = -5.0 * ii_402[k]
                   + f_0 * li_766[k];

        t_571[k] = -5.0 * ii_403[k]
                   + f_0 * li_767[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, ii_404, ii_405, ii_406, ii_407, \
                         ii_408, li_768, li_769, li_770, li_771, \
                         li_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -5.0 * ii_404[k]
                   + f_0 * li_768[k];

        t_573[k] = -5.0 * ii_405[k]
                   + f_0 * li_769[k];

        t_574[k] = -5.0 * ii_406[k]
                   + f_0 * li_770[k];

        t_575[k] = -5.0 * ii_407[k]
                   + f_0 * li_771[k];

        t_576[k] = -5.0 * ii_408[k]
                   + f_0 * li_772[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, ii_409, ii_410, ii_411, ii_412, \
                         ii_413, li_773, li_774, li_775, li_776, \
                         li_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -5.0 * ii_409[k]
                   + f_0 * li_773[k];

        t_578[k] = -5.0 * ii_410[k]
                   + f_0 * li_774[k];

        t_579[k] = -5.0 * ii_411[k]
                   + f_0 * li_775[k];

        t_580[k] = -5.0 * ii_412[k]
                   + f_0 * li_776[k];

        t_581[k] = -5.0 * ii_413[k]
                   + f_0 * li_777[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, ii_414, ii_415, ii_416, ii_417, \
                         ii_418, li_778, li_779, li_780, li_781, \
                         li_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -5.0 * ii_414[k]
                   + f_0 * li_778[k];

        t_583[k] = -5.0 * ii_415[k]
                   + f_0 * li_779[k];

        t_584[k] = -5.0 * ii_416[k]
                   + f_0 * li_780[k];

        t_585[k] = -5.0 * ii_417[k]
                   + f_0 * li_781[k];

        t_586[k] = -5.0 * ii_418[k]
                   + f_0 * li_782[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, t_591, t_592, t_593, ii_419, li_783, \
                         li_812, li_813, li_814, li_815, li_816, \
                         li_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -5.0 * ii_419[k]
                   + f_0 * li_783[k];

        t_588[k] = f_0 * li_812[k];

        t_589[k] = f_0 * li_813[k];

        t_590[k] = f_0 * li_814[k];

        t_591[k] = f_0 * li_815[k];

        t_592[k] = f_0 * li_816[k];

        t_593[k] = f_0 * li_817[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, t_599, t_600, t_601, li_818, \
                         li_819, li_820, li_821, li_822, li_823, li_824, \
                         li_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_0 * li_818[k];

        t_595[k] = f_0 * li_819[k];

        t_596[k] = f_0 * li_820[k];

        t_597[k] = f_0 * li_821[k];

        t_598[k] = f_0 * li_822[k];

        t_599[k] = f_0 * li_823[k];

        t_600[k] = f_0 * li_824[k];

        t_601[k] = f_0 * li_825[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, t_607, t_608, t_609, li_826, \
                         li_827, li_828, li_829, li_830, li_831, li_832, \
                         li_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_0 * li_826[k];

        t_603[k] = f_0 * li_827[k];

        t_604[k] = f_0 * li_828[k];

        t_605[k] = f_0 * li_829[k];

        t_606[k] = f_0 * li_830[k];

        t_607[k] = f_0 * li_831[k];

        t_608[k] = f_0 * li_832[k];

        t_609[k] = f_0 * li_833[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, t_615, t_616, ii_420, li_834, \
                         li_835, li_836, li_837, li_838, li_839, \
                         li_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_0 * li_834[k];

        t_611[k] = f_0 * li_835[k];

        t_612[k] = f_0 * li_836[k];

        t_613[k] = f_0 * li_837[k];

        t_614[k] = f_0 * li_838[k];

        t_615[k] = f_0 * li_839[k];

        t_616[k] = -ii_420[k]
                   + f_0 * li_840[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, ii_421, ii_422, ii_423, ii_424, \
                         ii_425, li_841, li_842, li_843, li_844, \
                         li_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = -ii_421[k]
                   + f_0 * li_841[k];

        t_618[k] = -ii_422[k]
                   + f_0 * li_842[k];

        t_619[k] = -ii_423[k]
                   + f_0 * li_843[k];

        t_620[k] = -ii_424[k]
                   + f_0 * li_844[k];

        t_621[k] = -ii_425[k]
                   + f_0 * li_845[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, ii_426, ii_427, ii_428, ii_429, \
                         ii_430, li_846, li_847, li_848, li_849, \
                         li_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = -ii_426[k]
                   + f_0 * li_846[k];

        t_623[k] = -ii_427[k]
                   + f_0 * li_847[k];

        t_624[k] = -ii_428[k]
                   + f_0 * li_848[k];

        t_625[k] = -ii_429[k]
                   + f_0 * li_849[k];

        t_626[k] = -ii_430[k]
                   + f_0 * li_850[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, t_630, t_631, ii_431, ii_432, ii_433, ii_434, \
                         ii_435, li_851, li_852, li_853, li_854, \
                         li_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = -ii_431[k]
                   + f_0 * li_851[k];

        t_628[k] = -ii_432[k]
                   + f_0 * li_852[k];

        t_629[k] = -ii_433[k]
                   + f_0 * li_853[k];

        t_630[k] = -ii_434[k]
                   + f_0 * li_854[k];

        t_631[k] = -ii_435[k]
                   + f_0 * li_855[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, t_636, ii_436, ii_437, ii_438, ii_439, \
                         ii_440, li_856, li_857, li_858, li_859, \
                         li_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = -ii_436[k]
                   + f_0 * li_856[k];

        t_633[k] = -ii_437[k]
                   + f_0 * li_857[k];

        t_634[k] = -ii_438[k]
                   + f_0 * li_858[k];

        t_635[k] = -ii_439[k]
                   + f_0 * li_859[k];

        t_636[k] = -ii_440[k]
                   + f_0 * li_860[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, t_641, ii_441, ii_442, ii_443, ii_444, \
                         ii_445, li_861, li_862, li_863, li_864, \
                         li_865 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = -ii_441[k]
                   + f_0 * li_861[k];

        t_638[k] = -ii_442[k]
                   + f_0 * li_862[k];

        t_639[k] = -ii_443[k]
                   + f_0 * li_863[k];

        t_640[k] = -ii_444[k]
                   + f_0 * li_864[k];

        t_641[k] = -ii_445[k]
                   + f_0 * li_865[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, t_646, ii_446, ii_447, ii_448, ii_449, \
                         ii_450, li_866, li_867, li_868, li_869, \
                         li_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = -ii_446[k]
                   + f_0 * li_866[k];

        t_643[k] = -ii_447[k]
                   + f_0 * li_867[k];

        t_644[k] = -2.0 * ii_448[k]
                   + f_0 * li_868[k];

        t_645[k] = -2.0 * ii_449[k]
                   + f_0 * li_869[k];

        t_646[k] = -2.0 * ii_450[k]
                   + f_0 * li_870[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, t_651, ii_451, ii_452, ii_453, ii_454, \
                         ii_455, li_871, li_872, li_873, li_874, \
                         li_875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = -2.0 * ii_451[k]
                   + f_0 * li_871[k];

        t_648[k] = -2.0 * ii_452[k]
                   + f_0 * li_872[k];

        t_649[k] = -2.0 * ii_453[k]
                   + f_0 * li_873[k];

        t_650[k] = -2.0 * ii_454[k]
                   + f_0 * li_874[k];

        t_651[k] = -2.0 * ii_455[k]
                   + f_0 * li_875[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, t_656, ii_456, ii_457, ii_458, ii_459, \
                         ii_460, li_876, li_877, li_878, li_879, \
                         li_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = -2.0 * ii_456[k]
                   + f_0 * li_876[k];

        t_653[k] = -2.0 * ii_457[k]
                   + f_0 * li_877[k];

        t_654[k] = -2.0 * ii_458[k]
                   + f_0 * li_878[k];

        t_655[k] = -2.0 * ii_459[k]
                   + f_0 * li_879[k];

        t_656[k] = -2.0 * ii_460[k]
                   + f_0 * li_880[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, ii_461, ii_462, ii_463, ii_464, \
                         ii_465, li_881, li_882, li_883, li_884, \
                         li_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = -2.0 * ii_461[k]
                   + f_0 * li_881[k];

        t_658[k] = -2.0 * ii_462[k]
                   + f_0 * li_882[k];

        t_659[k] = -2.0 * ii_463[k]
                   + f_0 * li_883[k];

        t_660[k] = -2.0 * ii_464[k]
                   + f_0 * li_884[k];

        t_661[k] = -2.0 * ii_465[k]
                   + f_0 * li_885[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, ii_466, ii_467, ii_468, ii_469, \
                         ii_470, li_886, li_887, li_888, li_889, \
                         li_890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = -2.0 * ii_466[k]
                   + f_0 * li_886[k];

        t_663[k] = -2.0 * ii_467[k]
                   + f_0 * li_887[k];

        t_664[k] = -2.0 * ii_468[k]
                   + f_0 * li_888[k];

        t_665[k] = -2.0 * ii_469[k]
                   + f_0 * li_889[k];

        t_666[k] = -2.0 * ii_470[k]
                   + f_0 * li_890[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_2_piece4(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
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

    const auto *ii_471 = buffer.data(ii + 471);
    const auto *ii_472 = buffer.data(ii + 472);
    const auto *ii_473 = buffer.data(ii + 473);
    const auto *ii_474 = buffer.data(ii + 474);
    const auto *ii_475 = buffer.data(ii + 475);
    const auto *ii_476 = buffer.data(ii + 476);
    const auto *ii_477 = buffer.data(ii + 477);
    const auto *ii_478 = buffer.data(ii + 478);
    const auto *ii_479 = buffer.data(ii + 479);
    const auto *ii_480 = buffer.data(ii + 480);
    const auto *ii_481 = buffer.data(ii + 481);
    const auto *ii_482 = buffer.data(ii + 482);
    const auto *ii_483 = buffer.data(ii + 483);
    const auto *ii_484 = buffer.data(ii + 484);
    const auto *ii_485 = buffer.data(ii + 485);
    const auto *ii_486 = buffer.data(ii + 486);
    const auto *ii_487 = buffer.data(ii + 487);
    const auto *ii_488 = buffer.data(ii + 488);
    const auto *ii_489 = buffer.data(ii + 489);
    const auto *ii_490 = buffer.data(ii + 490);
    const auto *ii_491 = buffer.data(ii + 491);
    const auto *ii_492 = buffer.data(ii + 492);
    const auto *ii_493 = buffer.data(ii + 493);
    const auto *ii_494 = buffer.data(ii + 494);
    const auto *ii_495 = buffer.data(ii + 495);
    const auto *ii_496 = buffer.data(ii + 496);
    const auto *ii_497 = buffer.data(ii + 497);
    const auto *ii_498 = buffer.data(ii + 498);
    const auto *ii_499 = buffer.data(ii + 499);
    const auto *ii_500 = buffer.data(ii + 500);
    const auto *ii_501 = buffer.data(ii + 501);
    const auto *ii_502 = buffer.data(ii + 502);
    const auto *ii_503 = buffer.data(ii + 503);
    const auto *ii_504 = buffer.data(ii + 504);
    const auto *ii_505 = buffer.data(ii + 505);
    const auto *ii_506 = buffer.data(ii + 506);
    const auto *ii_507 = buffer.data(ii + 507);
    const auto *ii_508 = buffer.data(ii + 508);
    const auto *ii_509 = buffer.data(ii + 509);
    const auto *ii_510 = buffer.data(ii + 510);
    const auto *ii_511 = buffer.data(ii + 511);
    const auto *ii_512 = buffer.data(ii + 512);
    const auto *ii_513 = buffer.data(ii + 513);
    const auto *ii_514 = buffer.data(ii + 514);
    const auto *ii_515 = buffer.data(ii + 515);
    const auto *ii_516 = buffer.data(ii + 516);
    const auto *ii_517 = buffer.data(ii + 517);
    const auto *ii_518 = buffer.data(ii + 518);
    const auto *ii_519 = buffer.data(ii + 519);
    const auto *ii_520 = buffer.data(ii + 520);
    const auto *ii_521 = buffer.data(ii + 521);
    const auto *ii_522 = buffer.data(ii + 522);
    const auto *ii_523 = buffer.data(ii + 523);
    const auto *ii_524 = buffer.data(ii + 524);
    const auto *ii_525 = buffer.data(ii + 525);
    const auto *ii_526 = buffer.data(ii + 526);
    const auto *ii_527 = buffer.data(ii + 527);
    const auto *ii_528 = buffer.data(ii + 528);
    const auto *ii_529 = buffer.data(ii + 529);
    const auto *ii_530 = buffer.data(ii + 530);
    const auto *ii_531 = buffer.data(ii + 531);
    const auto *ii_532 = buffer.data(ii + 532);
    const auto *ii_533 = buffer.data(ii + 533);
    const auto *ii_534 = buffer.data(ii + 534);
    const auto *ii_535 = buffer.data(ii + 535);
    const auto *ii_536 = buffer.data(ii + 536);
    const auto *ii_537 = buffer.data(ii + 537);
    const auto *ii_538 = buffer.data(ii + 538);
    const auto *ii_539 = buffer.data(ii + 539);
    const auto *ii_540 = buffer.data(ii + 540);
    const auto *ii_541 = buffer.data(ii + 541);
    const auto *ii_542 = buffer.data(ii + 542);
    const auto *ii_543 = buffer.data(ii + 543);
    const auto *ii_544 = buffer.data(ii + 544);
    const auto *ii_545 = buffer.data(ii + 545);
    const auto *ii_546 = buffer.data(ii + 546);
    const auto *ii_547 = buffer.data(ii + 547);
    const auto *ii_548 = buffer.data(ii + 548);
    const auto *ii_549 = buffer.data(ii + 549);
    const auto *ii_550 = buffer.data(ii + 550);
    const auto *ii_551 = buffer.data(ii + 551);
    const auto *ii_552 = buffer.data(ii + 552);
    const auto *ii_553 = buffer.data(ii + 553);
    const auto *ii_554 = buffer.data(ii + 554);
    const auto *ii_555 = buffer.data(ii + 555);
    const auto *ii_556 = buffer.data(ii + 556);
    const auto *ii_557 = buffer.data(ii + 557);
    const auto *ii_558 = buffer.data(ii + 558);
    const auto *ii_559 = buffer.data(ii + 559);
    const auto *ii_560 = buffer.data(ii + 560);
    const auto *ii_561 = buffer.data(ii + 561);
    const auto *ii_562 = buffer.data(ii + 562);
    const auto *ii_563 = buffer.data(ii + 563);
    const auto *ii_564 = buffer.data(ii + 564);
    const auto *ii_565 = buffer.data(ii + 565);
    const auto *ii_566 = buffer.data(ii + 566);
    const auto *ii_567 = buffer.data(ii + 567);
    const auto *ii_568 = buffer.data(ii + 568);
    const auto *ii_569 = buffer.data(ii + 569);
    const auto *ii_570 = buffer.data(ii + 570);
    const auto *ii_571 = buffer.data(ii + 571);
    const auto *ii_572 = buffer.data(ii + 572);
    const auto *ii_573 = buffer.data(ii + 573);
    const auto *ii_574 = buffer.data(ii + 574);
    const auto *ii_575 = buffer.data(ii + 575);
    const auto *ii_576 = buffer.data(ii + 576);
    const auto *ii_577 = buffer.data(ii + 577);
    const auto *ii_578 = buffer.data(ii + 578);
    const auto *ii_579 = buffer.data(ii + 579);
    const auto *ii_580 = buffer.data(ii + 580);
    const auto *ii_581 = buffer.data(ii + 581);
    const auto *ii_582 = buffer.data(ii + 582);
    const auto *ii_583 = buffer.data(ii + 583);
    const auto *ii_584 = buffer.data(ii + 584);
    const auto *ii_585 = buffer.data(ii + 585);
    const auto *ii_586 = buffer.data(ii + 586);
    const auto *ii_587 = buffer.data(ii + 587);
    const auto *ii_588 = buffer.data(ii + 588);
    const auto *ii_589 = buffer.data(ii + 589);
    const auto *ii_590 = buffer.data(ii + 590);
    const auto *ii_591 = buffer.data(ii + 591);
    const auto *ii_592 = buffer.data(ii + 592);
    const auto *ii_593 = buffer.data(ii + 593);
    const auto *ii_594 = buffer.data(ii + 594);
    const auto *ii_595 = buffer.data(ii + 595);
    const auto *ii_596 = buffer.data(ii + 596);
    const auto *ii_597 = buffer.data(ii + 597);
    const auto *ii_598 = buffer.data(ii + 598);
    const auto *ii_599 = buffer.data(ii + 599);
    const auto *ii_600 = buffer.data(ii + 600);
    const auto *ii_601 = buffer.data(ii + 601);
    const auto *ii_602 = buffer.data(ii + 602);

    const auto *li_891 = buffer.data(li + 891);
    const auto *li_892 = buffer.data(li + 892);
    const auto *li_893 = buffer.data(li + 893);
    const auto *li_894 = buffer.data(li + 894);
    const auto *li_895 = buffer.data(li + 895);
    const auto *li_896 = buffer.data(li + 896);
    const auto *li_897 = buffer.data(li + 897);
    const auto *li_898 = buffer.data(li + 898);
    const auto *li_899 = buffer.data(li + 899);
    const auto *li_900 = buffer.data(li + 900);
    const auto *li_901 = buffer.data(li + 901);
    const auto *li_902 = buffer.data(li + 902);
    const auto *li_903 = buffer.data(li + 903);
    const auto *li_904 = buffer.data(li + 904);
    const auto *li_905 = buffer.data(li + 905);
    const auto *li_906 = buffer.data(li + 906);
    const auto *li_907 = buffer.data(li + 907);
    const auto *li_908 = buffer.data(li + 908);
    const auto *li_909 = buffer.data(li + 909);
    const auto *li_910 = buffer.data(li + 910);
    const auto *li_911 = buffer.data(li + 911);
    const auto *li_912 = buffer.data(li + 912);
    const auto *li_913 = buffer.data(li + 913);
    const auto *li_914 = buffer.data(li + 914);
    const auto *li_915 = buffer.data(li + 915);
    const auto *li_916 = buffer.data(li + 916);
    const auto *li_917 = buffer.data(li + 917);
    const auto *li_918 = buffer.data(li + 918);
    const auto *li_919 = buffer.data(li + 919);
    const auto *li_920 = buffer.data(li + 920);
    const auto *li_921 = buffer.data(li + 921);
    const auto *li_922 = buffer.data(li + 922);
    const auto *li_923 = buffer.data(li + 923);
    const auto *li_924 = buffer.data(li + 924);
    const auto *li_925 = buffer.data(li + 925);
    const auto *li_926 = buffer.data(li + 926);
    const auto *li_927 = buffer.data(li + 927);
    const auto *li_928 = buffer.data(li + 928);
    const auto *li_929 = buffer.data(li + 929);
    const auto *li_930 = buffer.data(li + 930);
    const auto *li_931 = buffer.data(li + 931);
    const auto *li_932 = buffer.data(li + 932);
    const auto *li_933 = buffer.data(li + 933);
    const auto *li_934 = buffer.data(li + 934);
    const auto *li_935 = buffer.data(li + 935);
    const auto *li_936 = buffer.data(li + 936);
    const auto *li_937 = buffer.data(li + 937);
    const auto *li_938 = buffer.data(li + 938);
    const auto *li_939 = buffer.data(li + 939);
    const auto *li_940 = buffer.data(li + 940);
    const auto *li_941 = buffer.data(li + 941);
    const auto *li_942 = buffer.data(li + 942);
    const auto *li_943 = buffer.data(li + 943);
    const auto *li_944 = buffer.data(li + 944);
    const auto *li_945 = buffer.data(li + 945);
    const auto *li_946 = buffer.data(li + 946);
    const auto *li_947 = buffer.data(li + 947);
    const auto *li_948 = buffer.data(li + 948);
    const auto *li_949 = buffer.data(li + 949);
    const auto *li_950 = buffer.data(li + 950);
    const auto *li_951 = buffer.data(li + 951);
    const auto *li_952 = buffer.data(li + 952);
    const auto *li_953 = buffer.data(li + 953);
    const auto *li_954 = buffer.data(li + 954);
    const auto *li_955 = buffer.data(li + 955);
    const auto *li_956 = buffer.data(li + 956);
    const auto *li_957 = buffer.data(li + 957);
    const auto *li_958 = buffer.data(li + 958);
    const auto *li_959 = buffer.data(li + 959);
    const auto *li_960 = buffer.data(li + 960);
    const auto *li_961 = buffer.data(li + 961);
    const auto *li_962 = buffer.data(li + 962);
    const auto *li_963 = buffer.data(li + 963);
    const auto *li_964 = buffer.data(li + 964);
    const auto *li_965 = buffer.data(li + 965);
    const auto *li_966 = buffer.data(li + 966);
    const auto *li_967 = buffer.data(li + 967);
    const auto *li_968 = buffer.data(li + 968);
    const auto *li_969 = buffer.data(li + 969);
    const auto *li_970 = buffer.data(li + 970);
    const auto *li_971 = buffer.data(li + 971);
    const auto *li_972 = buffer.data(li + 972);
    const auto *li_973 = buffer.data(li + 973);
    const auto *li_974 = buffer.data(li + 974);
    const auto *li_975 = buffer.data(li + 975);
    const auto *li_976 = buffer.data(li + 976);
    const auto *li_977 = buffer.data(li + 977);
    const auto *li_978 = buffer.data(li + 978);
    const auto *li_979 = buffer.data(li + 979);
    const auto *li_980 = buffer.data(li + 980);
    const auto *li_981 = buffer.data(li + 981);
    const auto *li_982 = buffer.data(li + 982);
    const auto *li_983 = buffer.data(li + 983);
    const auto *li_984 = buffer.data(li + 984);
    const auto *li_985 = buffer.data(li + 985);
    const auto *li_986 = buffer.data(li + 986);
    const auto *li_987 = buffer.data(li + 987);
    const auto *li_988 = buffer.data(li + 988);
    const auto *li_989 = buffer.data(li + 989);
    const auto *li_990 = buffer.data(li + 990);
    const auto *li_991 = buffer.data(li + 991);
    const auto *li_992 = buffer.data(li + 992);
    const auto *li_993 = buffer.data(li + 993);
    const auto *li_994 = buffer.data(li + 994);
    const auto *li_995 = buffer.data(li + 995);
    const auto *li_996 = buffer.data(li + 996);
    const auto *li_997 = buffer.data(li + 997);
    const auto *li_998 = buffer.data(li + 998);
    const auto *li_999 = buffer.data(li + 999);
    const auto *li_1000 = buffer.data(li + 1000);
    const auto *li_1001 = buffer.data(li + 1001);
    const auto *li_1002 = buffer.data(li + 1002);
    const auto *li_1003 = buffer.data(li + 1003);
    const auto *li_1004 = buffer.data(li + 1004);
    const auto *li_1005 = buffer.data(li + 1005);
    const auto *li_1006 = buffer.data(li + 1006);
    const auto *li_1007 = buffer.data(li + 1007);
    const auto *li_1036 = buffer.data(li + 1036);
    const auto *li_1037 = buffer.data(li + 1037);
    const auto *li_1038 = buffer.data(li + 1038);
    const auto *li_1039 = buffer.data(li + 1039);
    const auto *li_1040 = buffer.data(li + 1040);
    const auto *li_1041 = buffer.data(li + 1041);
    const auto *li_1042 = buffer.data(li + 1042);
    const auto *li_1043 = buffer.data(li + 1043);
    const auto *li_1044 = buffer.data(li + 1044);
    const auto *li_1045 = buffer.data(li + 1045);
    const auto *li_1046 = buffer.data(li + 1046);
    const auto *li_1047 = buffer.data(li + 1047);
    const auto *li_1048 = buffer.data(li + 1048);
    const auto *li_1049 = buffer.data(li + 1049);
    const auto *li_1050 = buffer.data(li + 1050);
    const auto *li_1051 = buffer.data(li + 1051);
    const auto *li_1052 = buffer.data(li + 1052);
    const auto *li_1053 = buffer.data(li + 1053);
    const auto *li_1054 = buffer.data(li + 1054);
    const auto *li_1055 = buffer.data(li + 1055);
    const auto *li_1056 = buffer.data(li + 1056);
    const auto *li_1057 = buffer.data(li + 1057);
    const auto *li_1058 = buffer.data(li + 1058);
    const auto *li_1059 = buffer.data(li + 1059);
    const auto *li_1060 = buffer.data(li + 1060);
    const auto *li_1061 = buffer.data(li + 1061);
    const auto *li_1062 = buffer.data(li + 1062);
    const auto *li_1063 = buffer.data(li + 1063);
    const auto *li_1064 = buffer.data(li + 1064);
    const auto *li_1065 = buffer.data(li + 1065);
    const auto *li_1066 = buffer.data(li + 1066);
    const auto *li_1067 = buffer.data(li + 1067);
    const auto *li_1068 = buffer.data(li + 1068);
    const auto *li_1069 = buffer.data(li + 1069);
    const auto *li_1070 = buffer.data(li + 1070);
    const auto *li_1071 = buffer.data(li + 1071);
    const auto *li_1072 = buffer.data(li + 1072);
    const auto *li_1073 = buffer.data(li + 1073);
    const auto *li_1074 = buffer.data(li + 1074);
    const auto *li_1075 = buffer.data(li + 1075);
    const auto *li_1076 = buffer.data(li + 1076);
    const auto *li_1077 = buffer.data(li + 1077);
    const auto *li_1078 = buffer.data(li + 1078);

#pragma omp simd aligned(t_667, t_668, t_669, t_670, t_671, ii_471, ii_472, ii_473, ii_474, \
                         ii_475, li_891, li_892, li_893, li_894, \
                         li_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = -2.0 * ii_471[k]
                   + f_0 * li_891[k];

        t_668[k] = -2.0 * ii_472[k]
                   + f_0 * li_892[k];

        t_669[k] = -2.0 * ii_473[k]
                   + f_0 * li_893[k];

        t_670[k] = -2.0 * ii_474[k]
                   + f_0 * li_894[k];

        t_671[k] = -2.0 * ii_475[k]
                   + f_0 * li_895[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, ii_476, ii_477, ii_478, ii_479, \
                         ii_480, li_896, li_897, li_898, li_899, \
                         li_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = -3.0 * ii_476[k]
                   + f_0 * li_896[k];

        t_673[k] = -3.0 * ii_477[k]
                   + f_0 * li_897[k];

        t_674[k] = -3.0 * ii_478[k]
                   + f_0 * li_898[k];

        t_675[k] = -3.0 * ii_479[k]
                   + f_0 * li_899[k];

        t_676[k] = -3.0 * ii_480[k]
                   + f_0 * li_900[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, t_681, ii_481, ii_482, ii_483, ii_484, \
                         ii_485, li_901, li_902, li_903, li_904, \
                         li_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = -3.0 * ii_481[k]
                   + f_0 * li_901[k];

        t_678[k] = -3.0 * ii_482[k]
                   + f_0 * li_902[k];

        t_679[k] = -3.0 * ii_483[k]
                   + f_0 * li_903[k];

        t_680[k] = -3.0 * ii_484[k]
                   + f_0 * li_904[k];

        t_681[k] = -3.0 * ii_485[k]
                   + f_0 * li_905[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, ii_486, ii_487, ii_488, ii_489, \
                         ii_490, li_906, li_907, li_908, li_909, \
                         li_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = -3.0 * ii_486[k]
                   + f_0 * li_906[k];

        t_683[k] = -3.0 * ii_487[k]
                   + f_0 * li_907[k];

        t_684[k] = -3.0 * ii_488[k]
                   + f_0 * li_908[k];

        t_685[k] = -3.0 * ii_489[k]
                   + f_0 * li_909[k];

        t_686[k] = -3.0 * ii_490[k]
                   + f_0 * li_910[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, t_691, ii_491, ii_492, ii_493, ii_494, \
                         ii_495, li_911, li_912, li_913, li_914, \
                         li_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = -3.0 * ii_491[k]
                   + f_0 * li_911[k];

        t_688[k] = -3.0 * ii_492[k]
                   + f_0 * li_912[k];

        t_689[k] = -3.0 * ii_493[k]
                   + f_0 * li_913[k];

        t_690[k] = -3.0 * ii_494[k]
                   + f_0 * li_914[k];

        t_691[k] = -3.0 * ii_495[k]
                   + f_0 * li_915[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, t_696, ii_496, ii_497, ii_498, ii_499, \
                         ii_500, li_916, li_917, li_918, li_919, \
                         li_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = -3.0 * ii_496[k]
                   + f_0 * li_916[k];

        t_693[k] = -3.0 * ii_497[k]
                   + f_0 * li_917[k];

        t_694[k] = -3.0 * ii_498[k]
                   + f_0 * li_918[k];

        t_695[k] = -3.0 * ii_499[k]
                   + f_0 * li_919[k];

        t_696[k] = -3.0 * ii_500[k]
                   + f_0 * li_920[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, t_701, ii_501, ii_502, ii_503, ii_504, \
                         ii_505, li_921, li_922, li_923, li_924, \
                         li_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = -3.0 * ii_501[k]
                   + f_0 * li_921[k];

        t_698[k] = -3.0 * ii_502[k]
                   + f_0 * li_922[k];

        t_699[k] = -3.0 * ii_503[k]
                   + f_0 * li_923[k];

        t_700[k] = -4.0 * ii_504[k]
                   + f_0 * li_924[k];

        t_701[k] = -4.0 * ii_505[k]
                   + f_0 * li_925[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, t_706, ii_506, ii_507, ii_508, ii_509, \
                         ii_510, li_926, li_927, li_928, li_929, \
                         li_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = -4.0 * ii_506[k]
                   + f_0 * li_926[k];

        t_703[k] = -4.0 * ii_507[k]
                   + f_0 * li_927[k];

        t_704[k] = -4.0 * ii_508[k]
                   + f_0 * li_928[k];

        t_705[k] = -4.0 * ii_509[k]
                   + f_0 * li_929[k];

        t_706[k] = -4.0 * ii_510[k]
                   + f_0 * li_930[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, t_710, t_711, ii_511, ii_512, ii_513, ii_514, \
                         ii_515, li_931, li_932, li_933, li_934, \
                         li_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = -4.0 * ii_511[k]
                   + f_0 * li_931[k];

        t_708[k] = -4.0 * ii_512[k]
                   + f_0 * li_932[k];

        t_709[k] = -4.0 * ii_513[k]
                   + f_0 * li_933[k];

        t_710[k] = -4.0 * ii_514[k]
                   + f_0 * li_934[k];

        t_711[k] = -4.0 * ii_515[k]
                   + f_0 * li_935[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, t_716, ii_516, ii_517, ii_518, ii_519, \
                         ii_520, li_936, li_937, li_938, li_939, \
                         li_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = -4.0 * ii_516[k]
                   + f_0 * li_936[k];

        t_713[k] = -4.0 * ii_517[k]
                   + f_0 * li_937[k];

        t_714[k] = -4.0 * ii_518[k]
                   + f_0 * li_938[k];

        t_715[k] = -4.0 * ii_519[k]
                   + f_0 * li_939[k];

        t_716[k] = -4.0 * ii_520[k]
                   + f_0 * li_940[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, t_721, ii_521, ii_522, ii_523, ii_524, \
                         ii_525, li_941, li_942, li_943, li_944, \
                         li_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = -4.0 * ii_521[k]
                   + f_0 * li_941[k];

        t_718[k] = -4.0 * ii_522[k]
                   + f_0 * li_942[k];

        t_719[k] = -4.0 * ii_523[k]
                   + f_0 * li_943[k];

        t_720[k] = -4.0 * ii_524[k]
                   + f_0 * li_944[k];

        t_721[k] = -4.0 * ii_525[k]
                   + f_0 * li_945[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, t_726, ii_526, ii_527, ii_528, ii_529, \
                         ii_530, li_946, li_947, li_948, li_949, \
                         li_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = -4.0 * ii_526[k]
                   + f_0 * li_946[k];

        t_723[k] = -4.0 * ii_527[k]
                   + f_0 * li_947[k];

        t_724[k] = -4.0 * ii_528[k]
                   + f_0 * li_948[k];

        t_725[k] = -4.0 * ii_529[k]
                   + f_0 * li_949[k];

        t_726[k] = -4.0 * ii_530[k]
                   + f_0 * li_950[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, t_730, t_731, ii_531, ii_532, ii_533, ii_534, \
                         ii_535, li_951, li_952, li_953, li_954, \
                         li_955 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = -4.0 * ii_531[k]
                   + f_0 * li_951[k];

        t_728[k] = -5.0 * ii_532[k]
                   + f_0 * li_952[k];

        t_729[k] = -5.0 * ii_533[k]
                   + f_0 * li_953[k];

        t_730[k] = -5.0 * ii_534[k]
                   + f_0 * li_954[k];

        t_731[k] = -5.0 * ii_535[k]
                   + f_0 * li_955[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, t_735, t_736, ii_536, ii_537, ii_538, ii_539, \
                         ii_540, li_956, li_957, li_958, li_959, \
                         li_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = -5.0 * ii_536[k]
                   + f_0 * li_956[k];

        t_733[k] = -5.0 * ii_537[k]
                   + f_0 * li_957[k];

        t_734[k] = -5.0 * ii_538[k]
                   + f_0 * li_958[k];

        t_735[k] = -5.0 * ii_539[k]
                   + f_0 * li_959[k];

        t_736[k] = -5.0 * ii_540[k]
                   + f_0 * li_960[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, t_740, t_741, ii_541, ii_542, ii_543, ii_544, \
                         ii_545, li_961, li_962, li_963, li_964, \
                         li_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = -5.0 * ii_541[k]
                   + f_0 * li_961[k];

        t_738[k] = -5.0 * ii_542[k]
                   + f_0 * li_962[k];

        t_739[k] = -5.0 * ii_543[k]
                   + f_0 * li_963[k];

        t_740[k] = -5.0 * ii_544[k]
                   + f_0 * li_964[k];

        t_741[k] = -5.0 * ii_545[k]
                   + f_0 * li_965[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, t_746, ii_546, ii_547, ii_548, ii_549, \
                         ii_550, li_966, li_967, li_968, li_969, \
                         li_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = -5.0 * ii_546[k]
                   + f_0 * li_966[k];

        t_743[k] = -5.0 * ii_547[k]
                   + f_0 * li_967[k];

        t_744[k] = -5.0 * ii_548[k]
                   + f_0 * li_968[k];

        t_745[k] = -5.0 * ii_549[k]
                   + f_0 * li_969[k];

        t_746[k] = -5.0 * ii_550[k]
                   + f_0 * li_970[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, t_750, t_751, ii_551, ii_552, ii_553, ii_554, \
                         ii_555, li_971, li_972, li_973, li_974, \
                         li_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = -5.0 * ii_551[k]
                   + f_0 * li_971[k];

        t_748[k] = -5.0 * ii_552[k]
                   + f_0 * li_972[k];

        t_749[k] = -5.0 * ii_553[k]
                   + f_0 * li_973[k];

        t_750[k] = -5.0 * ii_554[k]
                   + f_0 * li_974[k];

        t_751[k] = -5.0 * ii_555[k]
                   + f_0 * li_975[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, t_756, ii_556, ii_557, ii_558, ii_559, \
                         ii_560, li_976, li_977, li_978, li_979, \
                         li_980 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = -5.0 * ii_556[k]
                   + f_0 * li_976[k];

        t_753[k] = -5.0 * ii_557[k]
                   + f_0 * li_977[k];

        t_754[k] = -5.0 * ii_558[k]
                   + f_0 * li_978[k];

        t_755[k] = -5.0 * ii_559[k]
                   + f_0 * li_979[k];

        t_756[k] = -6.0 * ii_560[k]
                   + f_0 * li_980[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, t_761, ii_561, ii_562, ii_563, ii_564, \
                         ii_565, li_981, li_982, li_983, li_984, \
                         li_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = -6.0 * ii_561[k]
                   + f_0 * li_981[k];

        t_758[k] = -6.0 * ii_562[k]
                   + f_0 * li_982[k];

        t_759[k] = -6.0 * ii_563[k]
                   + f_0 * li_983[k];

        t_760[k] = -6.0 * ii_564[k]
                   + f_0 * li_984[k];

        t_761[k] = -6.0 * ii_565[k]
                   + f_0 * li_985[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, t_766, ii_566, ii_567, ii_568, ii_569, \
                         ii_570, li_986, li_987, li_988, li_989, \
                         li_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = -6.0 * ii_566[k]
                   + f_0 * li_986[k];

        t_763[k] = -6.0 * ii_567[k]
                   + f_0 * li_987[k];

        t_764[k] = -6.0 * ii_568[k]
                   + f_0 * li_988[k];

        t_765[k] = -6.0 * ii_569[k]
                   + f_0 * li_989[k];

        t_766[k] = -6.0 * ii_570[k]
                   + f_0 * li_990[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, t_770, t_771, ii_571, ii_572, ii_573, ii_574, \
                         ii_575, li_991, li_992, li_993, li_994, \
                         li_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = -6.0 * ii_571[k]
                   + f_0 * li_991[k];

        t_768[k] = -6.0 * ii_572[k]
                   + f_0 * li_992[k];

        t_769[k] = -6.0 * ii_573[k]
                   + f_0 * li_993[k];

        t_770[k] = -6.0 * ii_574[k]
                   + f_0 * li_994[k];

        t_771[k] = -6.0 * ii_575[k]
                   + f_0 * li_995[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, ii_576, ii_577, ii_578, ii_579, \
                         ii_580, li_996, li_997, li_998, li_999, \
                         li_1000 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = -6.0 * ii_576[k]
                   + f_0 * li_996[k];

        t_773[k] = -6.0 * ii_577[k]
                   + f_0 * li_997[k];

        t_774[k] = -6.0 * ii_578[k]
                   + f_0 * li_998[k];

        t_775[k] = -6.0 * ii_579[k]
                   + f_0 * li_999[k];

        t_776[k] = -6.0 * ii_580[k]
                   + f_0 * li_1000[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, t_780, t_781, ii_581, ii_582, ii_583, ii_584, \
                         ii_585, li_1001, li_1002, li_1003, li_1004, \
                         li_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = -6.0 * ii_581[k]
                   + f_0 * li_1001[k];

        t_778[k] = -6.0 * ii_582[k]
                   + f_0 * li_1002[k];

        t_779[k] = -6.0 * ii_583[k]
                   + f_0 * li_1003[k];

        t_780[k] = -6.0 * ii_584[k]
                   + f_0 * li_1004[k];

        t_781[k] = -6.0 * ii_585[k]
                   + f_0 * li_1005[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, t_786, t_787, t_788, ii_586, ii_587, \
                         li_1006, li_1007, li_1036, li_1037, li_1038, li_1039, \
                         li_1040 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = -6.0 * ii_586[k]
                   + f_0 * li_1006[k];

        t_783[k] = -6.0 * ii_587[k]
                   + f_0 * li_1007[k];

        t_784[k] = f_0 * li_1036[k];

        t_785[k] = f_0 * li_1037[k];

        t_786[k] = f_0 * li_1038[k];

        t_787[k] = f_0 * li_1039[k];

        t_788[k] = f_0 * li_1040[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, t_793, t_794, t_795, t_796, li_1041, \
                         li_1042, li_1043, li_1044, li_1045, li_1046, li_1047, \
                         li_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_0 * li_1041[k];

        t_790[k] = f_0 * li_1042[k];

        t_791[k] = f_0 * li_1043[k];

        t_792[k] = f_0 * li_1044[k];

        t_793[k] = f_0 * li_1045[k];

        t_794[k] = f_0 * li_1046[k];

        t_795[k] = f_0 * li_1047[k];

        t_796[k] = f_0 * li_1048[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, t_800, t_801, t_802, t_803, t_804, li_1049, \
                         li_1050, li_1051, li_1052, li_1053, li_1054, li_1055, \
                         li_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_0 * li_1049[k];

        t_798[k] = f_0 * li_1050[k];

        t_799[k] = f_0 * li_1051[k];

        t_800[k] = f_0 * li_1052[k];

        t_801[k] = f_0 * li_1053[k];

        t_802[k] = f_0 * li_1054[k];

        t_803[k] = f_0 * li_1055[k];

        t_804[k] = f_0 * li_1056[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, t_810, t_811, li_1057, li_1058, \
                         li_1059, li_1060, li_1061, li_1062, li_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_0 * li_1057[k];

        t_806[k] = f_0 * li_1058[k];

        t_807[k] = f_0 * li_1059[k];

        t_808[k] = f_0 * li_1060[k];

        t_809[k] = f_0 * li_1061[k];

        t_810[k] = f_0 * li_1062[k];

        t_811[k] = f_0 * li_1063[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, t_816, ii_588, ii_589, ii_590, ii_591, \
                         ii_592, li_1064, li_1065, li_1066, li_1067, \
                         li_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = -ii_588[k]
                   + f_0 * li_1064[k];

        t_813[k] = -ii_589[k]
                   + f_0 * li_1065[k];

        t_814[k] = -ii_590[k]
                   + f_0 * li_1066[k];

        t_815[k] = -ii_591[k]
                   + f_0 * li_1067[k];

        t_816[k] = -ii_592[k]
                   + f_0 * li_1068[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, t_820, t_821, ii_593, ii_594, ii_595, ii_596, \
                         ii_597, li_1069, li_1070, li_1071, li_1072, \
                         li_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = -ii_593[k]
                   + f_0 * li_1069[k];

        t_818[k] = -ii_594[k]
                   + f_0 * li_1070[k];

        t_819[k] = -ii_595[k]
                   + f_0 * li_1071[k];

        t_820[k] = -ii_596[k]
                   + f_0 * li_1072[k];

        t_821[k] = -ii_597[k]
                   + f_0 * li_1073[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, t_826, ii_598, ii_599, ii_600, ii_601, \
                         ii_602, li_1074, li_1075, li_1076, li_1077, \
                         li_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = -ii_598[k]
                   + f_0 * li_1074[k];

        t_823[k] = -ii_599[k]
                   + f_0 * li_1075[k];

        t_824[k] = -ii_600[k]
                   + f_0 * li_1076[k];

        t_825[k] = -ii_601[k]
                   + f_0 * li_1077[k];

        t_826[k] = -ii_602[k]
                   + f_0 * li_1078[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_2_piece5(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ii_603 = buffer.data(ii + 603);
    const auto *ii_604 = buffer.data(ii + 604);
    const auto *ii_605 = buffer.data(ii + 605);
    const auto *ii_606 = buffer.data(ii + 606);
    const auto *ii_607 = buffer.data(ii + 607);
    const auto *ii_608 = buffer.data(ii + 608);
    const auto *ii_609 = buffer.data(ii + 609);
    const auto *ii_610 = buffer.data(ii + 610);
    const auto *ii_611 = buffer.data(ii + 611);
    const auto *ii_612 = buffer.data(ii + 612);
    const auto *ii_613 = buffer.data(ii + 613);
    const auto *ii_614 = buffer.data(ii + 614);
    const auto *ii_615 = buffer.data(ii + 615);
    const auto *ii_616 = buffer.data(ii + 616);
    const auto *ii_617 = buffer.data(ii + 617);
    const auto *ii_618 = buffer.data(ii + 618);
    const auto *ii_619 = buffer.data(ii + 619);
    const auto *ii_620 = buffer.data(ii + 620);
    const auto *ii_621 = buffer.data(ii + 621);
    const auto *ii_622 = buffer.data(ii + 622);
    const auto *ii_623 = buffer.data(ii + 623);
    const auto *ii_624 = buffer.data(ii + 624);
    const auto *ii_625 = buffer.data(ii + 625);
    const auto *ii_626 = buffer.data(ii + 626);
    const auto *ii_627 = buffer.data(ii + 627);
    const auto *ii_628 = buffer.data(ii + 628);
    const auto *ii_629 = buffer.data(ii + 629);
    const auto *ii_630 = buffer.data(ii + 630);
    const auto *ii_631 = buffer.data(ii + 631);
    const auto *ii_632 = buffer.data(ii + 632);
    const auto *ii_633 = buffer.data(ii + 633);
    const auto *ii_634 = buffer.data(ii + 634);
    const auto *ii_635 = buffer.data(ii + 635);
    const auto *ii_636 = buffer.data(ii + 636);
    const auto *ii_637 = buffer.data(ii + 637);
    const auto *ii_638 = buffer.data(ii + 638);
    const auto *ii_639 = buffer.data(ii + 639);
    const auto *ii_640 = buffer.data(ii + 640);
    const auto *ii_641 = buffer.data(ii + 641);
    const auto *ii_642 = buffer.data(ii + 642);
    const auto *ii_643 = buffer.data(ii + 643);
    const auto *ii_644 = buffer.data(ii + 644);
    const auto *ii_645 = buffer.data(ii + 645);
    const auto *ii_646 = buffer.data(ii + 646);
    const auto *ii_647 = buffer.data(ii + 647);
    const auto *ii_648 = buffer.data(ii + 648);
    const auto *ii_649 = buffer.data(ii + 649);
    const auto *ii_650 = buffer.data(ii + 650);
    const auto *ii_651 = buffer.data(ii + 651);
    const auto *ii_652 = buffer.data(ii + 652);
    const auto *ii_653 = buffer.data(ii + 653);
    const auto *ii_654 = buffer.data(ii + 654);
    const auto *ii_655 = buffer.data(ii + 655);
    const auto *ii_656 = buffer.data(ii + 656);
    const auto *ii_657 = buffer.data(ii + 657);
    const auto *ii_658 = buffer.data(ii + 658);
    const auto *ii_659 = buffer.data(ii + 659);
    const auto *ii_660 = buffer.data(ii + 660);
    const auto *ii_661 = buffer.data(ii + 661);
    const auto *ii_662 = buffer.data(ii + 662);
    const auto *ii_663 = buffer.data(ii + 663);
    const auto *ii_664 = buffer.data(ii + 664);
    const auto *ii_665 = buffer.data(ii + 665);
    const auto *ii_666 = buffer.data(ii + 666);
    const auto *ii_667 = buffer.data(ii + 667);
    const auto *ii_668 = buffer.data(ii + 668);
    const auto *ii_669 = buffer.data(ii + 669);
    const auto *ii_670 = buffer.data(ii + 670);
    const auto *ii_671 = buffer.data(ii + 671);
    const auto *ii_672 = buffer.data(ii + 672);
    const auto *ii_673 = buffer.data(ii + 673);
    const auto *ii_674 = buffer.data(ii + 674);
    const auto *ii_675 = buffer.data(ii + 675);
    const auto *ii_676 = buffer.data(ii + 676);
    const auto *ii_677 = buffer.data(ii + 677);
    const auto *ii_678 = buffer.data(ii + 678);
    const auto *ii_679 = buffer.data(ii + 679);
    const auto *ii_680 = buffer.data(ii + 680);
    const auto *ii_681 = buffer.data(ii + 681);
    const auto *ii_682 = buffer.data(ii + 682);
    const auto *ii_683 = buffer.data(ii + 683);
    const auto *ii_684 = buffer.data(ii + 684);
    const auto *ii_685 = buffer.data(ii + 685);
    const auto *ii_686 = buffer.data(ii + 686);
    const auto *ii_687 = buffer.data(ii + 687);
    const auto *ii_688 = buffer.data(ii + 688);
    const auto *ii_689 = buffer.data(ii + 689);
    const auto *ii_690 = buffer.data(ii + 690);
    const auto *ii_691 = buffer.data(ii + 691);
    const auto *ii_692 = buffer.data(ii + 692);
    const auto *ii_693 = buffer.data(ii + 693);
    const auto *ii_694 = buffer.data(ii + 694);
    const auto *ii_695 = buffer.data(ii + 695);
    const auto *ii_696 = buffer.data(ii + 696);
    const auto *ii_697 = buffer.data(ii + 697);
    const auto *ii_698 = buffer.data(ii + 698);
    const auto *ii_699 = buffer.data(ii + 699);
    const auto *ii_700 = buffer.data(ii + 700);
    const auto *ii_701 = buffer.data(ii + 701);
    const auto *ii_702 = buffer.data(ii + 702);
    const auto *ii_703 = buffer.data(ii + 703);
    const auto *ii_704 = buffer.data(ii + 704);
    const auto *ii_705 = buffer.data(ii + 705);
    const auto *ii_706 = buffer.data(ii + 706);
    const auto *ii_707 = buffer.data(ii + 707);
    const auto *ii_708 = buffer.data(ii + 708);
    const auto *ii_709 = buffer.data(ii + 709);
    const auto *ii_710 = buffer.data(ii + 710);
    const auto *ii_711 = buffer.data(ii + 711);
    const auto *ii_712 = buffer.data(ii + 712);
    const auto *ii_713 = buffer.data(ii + 713);
    const auto *ii_714 = buffer.data(ii + 714);
    const auto *ii_715 = buffer.data(ii + 715);
    const auto *ii_716 = buffer.data(ii + 716);
    const auto *ii_717 = buffer.data(ii + 717);
    const auto *ii_718 = buffer.data(ii + 718);
    const auto *ii_719 = buffer.data(ii + 719);
    const auto *ii_720 = buffer.data(ii + 720);
    const auto *ii_721 = buffer.data(ii + 721);
    const auto *ii_722 = buffer.data(ii + 722);
    const auto *ii_723 = buffer.data(ii + 723);
    const auto *ii_724 = buffer.data(ii + 724);
    const auto *ii_725 = buffer.data(ii + 725);
    const auto *ii_726 = buffer.data(ii + 726);
    const auto *ii_727 = buffer.data(ii + 727);
    const auto *ii_728 = buffer.data(ii + 728);
    const auto *ii_729 = buffer.data(ii + 729);
    const auto *ii_730 = buffer.data(ii + 730);
    const auto *ii_731 = buffer.data(ii + 731);
    const auto *ii_732 = buffer.data(ii + 732);
    const auto *ii_733 = buffer.data(ii + 733);
    const auto *ii_734 = buffer.data(ii + 734);
    const auto *ii_735 = buffer.data(ii + 735);
    const auto *ii_736 = buffer.data(ii + 736);
    const auto *ii_737 = buffer.data(ii + 737);
    const auto *ii_738 = buffer.data(ii + 738);
    const auto *ii_739 = buffer.data(ii + 739);
    const auto *ii_740 = buffer.data(ii + 740);
    const auto *ii_741 = buffer.data(ii + 741);
    const auto *ii_742 = buffer.data(ii + 742);
    const auto *ii_743 = buffer.data(ii + 743);
    const auto *ii_744 = buffer.data(ii + 744);
    const auto *ii_745 = buffer.data(ii + 745);
    const auto *ii_746 = buffer.data(ii + 746);
    const auto *ii_747 = buffer.data(ii + 747);
    const auto *ii_748 = buffer.data(ii + 748);
    const auto *ii_749 = buffer.data(ii + 749);
    const auto *ii_750 = buffer.data(ii + 750);
    const auto *ii_751 = buffer.data(ii + 751);
    const auto *ii_752 = buffer.data(ii + 752);

    const auto *li_1079 = buffer.data(li + 1079);
    const auto *li_1080 = buffer.data(li + 1080);
    const auto *li_1081 = buffer.data(li + 1081);
    const auto *li_1082 = buffer.data(li + 1082);
    const auto *li_1083 = buffer.data(li + 1083);
    const auto *li_1084 = buffer.data(li + 1084);
    const auto *li_1085 = buffer.data(li + 1085);
    const auto *li_1086 = buffer.data(li + 1086);
    const auto *li_1087 = buffer.data(li + 1087);
    const auto *li_1088 = buffer.data(li + 1088);
    const auto *li_1089 = buffer.data(li + 1089);
    const auto *li_1090 = buffer.data(li + 1090);
    const auto *li_1091 = buffer.data(li + 1091);
    const auto *li_1092 = buffer.data(li + 1092);
    const auto *li_1093 = buffer.data(li + 1093);
    const auto *li_1094 = buffer.data(li + 1094);
    const auto *li_1095 = buffer.data(li + 1095);
    const auto *li_1096 = buffer.data(li + 1096);
    const auto *li_1097 = buffer.data(li + 1097);
    const auto *li_1098 = buffer.data(li + 1098);
    const auto *li_1099 = buffer.data(li + 1099);
    const auto *li_1100 = buffer.data(li + 1100);
    const auto *li_1101 = buffer.data(li + 1101);
    const auto *li_1102 = buffer.data(li + 1102);
    const auto *li_1103 = buffer.data(li + 1103);
    const auto *li_1104 = buffer.data(li + 1104);
    const auto *li_1105 = buffer.data(li + 1105);
    const auto *li_1106 = buffer.data(li + 1106);
    const auto *li_1107 = buffer.data(li + 1107);
    const auto *li_1108 = buffer.data(li + 1108);
    const auto *li_1109 = buffer.data(li + 1109);
    const auto *li_1110 = buffer.data(li + 1110);
    const auto *li_1111 = buffer.data(li + 1111);
    const auto *li_1112 = buffer.data(li + 1112);
    const auto *li_1113 = buffer.data(li + 1113);
    const auto *li_1114 = buffer.data(li + 1114);
    const auto *li_1115 = buffer.data(li + 1115);
    const auto *li_1116 = buffer.data(li + 1116);
    const auto *li_1117 = buffer.data(li + 1117);
    const auto *li_1118 = buffer.data(li + 1118);
    const auto *li_1119 = buffer.data(li + 1119);
    const auto *li_1120 = buffer.data(li + 1120);
    const auto *li_1121 = buffer.data(li + 1121);
    const auto *li_1122 = buffer.data(li + 1122);
    const auto *li_1123 = buffer.data(li + 1123);
    const auto *li_1124 = buffer.data(li + 1124);
    const auto *li_1125 = buffer.data(li + 1125);
    const auto *li_1126 = buffer.data(li + 1126);
    const auto *li_1127 = buffer.data(li + 1127);
    const auto *li_1128 = buffer.data(li + 1128);
    const auto *li_1129 = buffer.data(li + 1129);
    const auto *li_1130 = buffer.data(li + 1130);
    const auto *li_1131 = buffer.data(li + 1131);
    const auto *li_1132 = buffer.data(li + 1132);
    const auto *li_1133 = buffer.data(li + 1133);
    const auto *li_1134 = buffer.data(li + 1134);
    const auto *li_1135 = buffer.data(li + 1135);
    const auto *li_1136 = buffer.data(li + 1136);
    const auto *li_1137 = buffer.data(li + 1137);
    const auto *li_1138 = buffer.data(li + 1138);
    const auto *li_1139 = buffer.data(li + 1139);
    const auto *li_1140 = buffer.data(li + 1140);
    const auto *li_1141 = buffer.data(li + 1141);
    const auto *li_1142 = buffer.data(li + 1142);
    const auto *li_1143 = buffer.data(li + 1143);
    const auto *li_1144 = buffer.data(li + 1144);
    const auto *li_1145 = buffer.data(li + 1145);
    const auto *li_1146 = buffer.data(li + 1146);
    const auto *li_1147 = buffer.data(li + 1147);
    const auto *li_1148 = buffer.data(li + 1148);
    const auto *li_1149 = buffer.data(li + 1149);
    const auto *li_1150 = buffer.data(li + 1150);
    const auto *li_1151 = buffer.data(li + 1151);
    const auto *li_1152 = buffer.data(li + 1152);
    const auto *li_1153 = buffer.data(li + 1153);
    const auto *li_1154 = buffer.data(li + 1154);
    const auto *li_1155 = buffer.data(li + 1155);
    const auto *li_1156 = buffer.data(li + 1156);
    const auto *li_1157 = buffer.data(li + 1157);
    const auto *li_1158 = buffer.data(li + 1158);
    const auto *li_1159 = buffer.data(li + 1159);
    const auto *li_1160 = buffer.data(li + 1160);
    const auto *li_1161 = buffer.data(li + 1161);
    const auto *li_1162 = buffer.data(li + 1162);
    const auto *li_1163 = buffer.data(li + 1163);
    const auto *li_1164 = buffer.data(li + 1164);
    const auto *li_1165 = buffer.data(li + 1165);
    const auto *li_1166 = buffer.data(li + 1166);
    const auto *li_1167 = buffer.data(li + 1167);
    const auto *li_1168 = buffer.data(li + 1168);
    const auto *li_1169 = buffer.data(li + 1169);
    const auto *li_1170 = buffer.data(li + 1170);
    const auto *li_1171 = buffer.data(li + 1171);
    const auto *li_1172 = buffer.data(li + 1172);
    const auto *li_1173 = buffer.data(li + 1173);
    const auto *li_1174 = buffer.data(li + 1174);
    const auto *li_1175 = buffer.data(li + 1175);
    const auto *li_1176 = buffer.data(li + 1176);
    const auto *li_1177 = buffer.data(li + 1177);
    const auto *li_1178 = buffer.data(li + 1178);
    const auto *li_1179 = buffer.data(li + 1179);
    const auto *li_1180 = buffer.data(li + 1180);
    const auto *li_1181 = buffer.data(li + 1181);
    const auto *li_1182 = buffer.data(li + 1182);
    const auto *li_1183 = buffer.data(li + 1183);
    const auto *li_1184 = buffer.data(li + 1184);
    const auto *li_1185 = buffer.data(li + 1185);
    const auto *li_1186 = buffer.data(li + 1186);
    const auto *li_1187 = buffer.data(li + 1187);
    const auto *li_1188 = buffer.data(li + 1188);
    const auto *li_1189 = buffer.data(li + 1189);
    const auto *li_1190 = buffer.data(li + 1190);
    const auto *li_1191 = buffer.data(li + 1191);
    const auto *li_1192 = buffer.data(li + 1192);
    const auto *li_1193 = buffer.data(li + 1193);
    const auto *li_1194 = buffer.data(li + 1194);
    const auto *li_1195 = buffer.data(li + 1195);
    const auto *li_1196 = buffer.data(li + 1196);
    const auto *li_1197 = buffer.data(li + 1197);
    const auto *li_1198 = buffer.data(li + 1198);
    const auto *li_1199 = buffer.data(li + 1199);
    const auto *li_1200 = buffer.data(li + 1200);
    const auto *li_1201 = buffer.data(li + 1201);
    const auto *li_1202 = buffer.data(li + 1202);
    const auto *li_1203 = buffer.data(li + 1203);
    const auto *li_1204 = buffer.data(li + 1204);
    const auto *li_1205 = buffer.data(li + 1205);
    const auto *li_1206 = buffer.data(li + 1206);
    const auto *li_1207 = buffer.data(li + 1207);
    const auto *li_1208 = buffer.data(li + 1208);
    const auto *li_1209 = buffer.data(li + 1209);
    const auto *li_1210 = buffer.data(li + 1210);
    const auto *li_1211 = buffer.data(li + 1211);
    const auto *li_1212 = buffer.data(li + 1212);
    const auto *li_1213 = buffer.data(li + 1213);
    const auto *li_1214 = buffer.data(li + 1214);
    const auto *li_1215 = buffer.data(li + 1215);
    const auto *li_1216 = buffer.data(li + 1216);
    const auto *li_1217 = buffer.data(li + 1217);
    const auto *li_1218 = buffer.data(li + 1218);
    const auto *li_1219 = buffer.data(li + 1219);
    const auto *li_1220 = buffer.data(li + 1220);
    const auto *li_1221 = buffer.data(li + 1221);
    const auto *li_1222 = buffer.data(li + 1222);
    const auto *li_1223 = buffer.data(li + 1223);
    const auto *li_1224 = buffer.data(li + 1224);
    const auto *li_1225 = buffer.data(li + 1225);
    const auto *li_1226 = buffer.data(li + 1226);
    const auto *li_1227 = buffer.data(li + 1227);
    const auto *li_1228 = buffer.data(li + 1228);

#pragma omp simd aligned(t_827, t_828, t_829, t_830, t_831, ii_603, ii_604, ii_605, ii_606, \
                         ii_607, li_1079, li_1080, li_1081, li_1082, \
                         li_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = -ii_603[k]
                   + f_0 * li_1079[k];

        t_828[k] = -ii_604[k]
                   + f_0 * li_1080[k];

        t_829[k] = -ii_605[k]
                   + f_0 * li_1081[k];

        t_830[k] = -ii_606[k]
                   + f_0 * li_1082[k];

        t_831[k] = -ii_607[k]
                   + f_0 * li_1083[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, t_836, ii_608, ii_609, ii_610, ii_611, \
                         ii_612, li_1084, li_1085, li_1086, li_1087, \
                         li_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = -ii_608[k]
                   + f_0 * li_1084[k];

        t_833[k] = -ii_609[k]
                   + f_0 * li_1085[k];

        t_834[k] = -ii_610[k]
                   + f_0 * li_1086[k];

        t_835[k] = -ii_611[k]
                   + f_0 * li_1087[k];

        t_836[k] = -ii_612[k]
                   + f_0 * li_1088[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, t_840, t_841, ii_613, ii_614, ii_615, ii_616, \
                         ii_617, li_1089, li_1090, li_1091, li_1092, \
                         li_1093 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = -ii_613[k]
                   + f_0 * li_1089[k];

        t_838[k] = -ii_614[k]
                   + f_0 * li_1090[k];

        t_839[k] = -ii_615[k]
                   + f_0 * li_1091[k];

        t_840[k] = -2.0 * ii_616[k]
                   + f_0 * li_1092[k];

        t_841[k] = -2.0 * ii_617[k]
                   + f_0 * li_1093[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, t_845, t_846, ii_618, ii_619, ii_620, ii_621, \
                         ii_622, li_1094, li_1095, li_1096, li_1097, \
                         li_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = -2.0 * ii_618[k]
                   + f_0 * li_1094[k];

        t_843[k] = -2.0 * ii_619[k]
                   + f_0 * li_1095[k];

        t_844[k] = -2.0 * ii_620[k]
                   + f_0 * li_1096[k];

        t_845[k] = -2.0 * ii_621[k]
                   + f_0 * li_1097[k];

        t_846[k] = -2.0 * ii_622[k]
                   + f_0 * li_1098[k];
    }

#pragma omp simd aligned(t_847, t_848, t_849, t_850, t_851, ii_623, ii_624, ii_625, ii_626, \
                         ii_627, li_1099, li_1100, li_1101, li_1102, \
                         li_1103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_847[k] = -2.0 * ii_623[k]
                   + f_0 * li_1099[k];

        t_848[k] = -2.0 * ii_624[k]
                   + f_0 * li_1100[k];

        t_849[k] = -2.0 * ii_625[k]
                   + f_0 * li_1101[k];

        t_850[k] = -2.0 * ii_626[k]
                   + f_0 * li_1102[k];

        t_851[k] = -2.0 * ii_627[k]
                   + f_0 * li_1103[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, t_856, ii_628, ii_629, ii_630, ii_631, \
                         ii_632, li_1104, li_1105, li_1106, li_1107, \
                         li_1108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = -2.0 * ii_628[k]
                   + f_0 * li_1104[k];

        t_853[k] = -2.0 * ii_629[k]
                   + f_0 * li_1105[k];

        t_854[k] = -2.0 * ii_630[k]
                   + f_0 * li_1106[k];

        t_855[k] = -2.0 * ii_631[k]
                   + f_0 * li_1107[k];

        t_856[k] = -2.0 * ii_632[k]
                   + f_0 * li_1108[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, t_860, t_861, ii_633, ii_634, ii_635, ii_636, \
                         ii_637, li_1109, li_1110, li_1111, li_1112, \
                         li_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = -2.0 * ii_633[k]
                   + f_0 * li_1109[k];

        t_858[k] = -2.0 * ii_634[k]
                   + f_0 * li_1110[k];

        t_859[k] = -2.0 * ii_635[k]
                   + f_0 * li_1111[k];

        t_860[k] = -2.0 * ii_636[k]
                   + f_0 * li_1112[k];

        t_861[k] = -2.0 * ii_637[k]
                   + f_0 * li_1113[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, t_866, ii_638, ii_639, ii_640, ii_641, \
                         ii_642, li_1114, li_1115, li_1116, li_1117, \
                         li_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = -2.0 * ii_638[k]
                   + f_0 * li_1114[k];

        t_863[k] = -2.0 * ii_639[k]
                   + f_0 * li_1115[k];

        t_864[k] = -2.0 * ii_640[k]
                   + f_0 * li_1116[k];

        t_865[k] = -2.0 * ii_641[k]
                   + f_0 * li_1117[k];

        t_866[k] = -2.0 * ii_642[k]
                   + f_0 * li_1118[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, t_870, t_871, ii_643, ii_644, ii_645, ii_646, \
                         ii_647, li_1119, li_1120, li_1121, li_1122, \
                         li_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = -2.0 * ii_643[k]
                   + f_0 * li_1119[k];

        t_868[k] = -3.0 * ii_644[k]
                   + f_0 * li_1120[k];

        t_869[k] = -3.0 * ii_645[k]
                   + f_0 * li_1121[k];

        t_870[k] = -3.0 * ii_646[k]
                   + f_0 * li_1122[k];

        t_871[k] = -3.0 * ii_647[k]
                   + f_0 * li_1123[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, t_875, t_876, ii_648, ii_649, ii_650, ii_651, \
                         ii_652, li_1124, li_1125, li_1126, li_1127, \
                         li_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = -3.0 * ii_648[k]
                   + f_0 * li_1124[k];

        t_873[k] = -3.0 * ii_649[k]
                   + f_0 * li_1125[k];

        t_874[k] = -3.0 * ii_650[k]
                   + f_0 * li_1126[k];

        t_875[k] = -3.0 * ii_651[k]
                   + f_0 * li_1127[k];

        t_876[k] = -3.0 * ii_652[k]
                   + f_0 * li_1128[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, t_880, t_881, ii_653, ii_654, ii_655, ii_656, \
                         ii_657, li_1129, li_1130, li_1131, li_1132, \
                         li_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = -3.0 * ii_653[k]
                   + f_0 * li_1129[k];

        t_878[k] = -3.0 * ii_654[k]
                   + f_0 * li_1130[k];

        t_879[k] = -3.0 * ii_655[k]
                   + f_0 * li_1131[k];

        t_880[k] = -3.0 * ii_656[k]
                   + f_0 * li_1132[k];

        t_881[k] = -3.0 * ii_657[k]
                   + f_0 * li_1133[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, t_886, ii_658, ii_659, ii_660, ii_661, \
                         ii_662, li_1134, li_1135, li_1136, li_1137, \
                         li_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = -3.0 * ii_658[k]
                   + f_0 * li_1134[k];

        t_883[k] = -3.0 * ii_659[k]
                   + f_0 * li_1135[k];

        t_884[k] = -3.0 * ii_660[k]
                   + f_0 * li_1136[k];

        t_885[k] = -3.0 * ii_661[k]
                   + f_0 * li_1137[k];

        t_886[k] = -3.0 * ii_662[k]
                   + f_0 * li_1138[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, t_890, t_891, ii_663, ii_664, ii_665, ii_666, \
                         ii_667, li_1139, li_1140, li_1141, li_1142, \
                         li_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = -3.0 * ii_663[k]
                   + f_0 * li_1139[k];

        t_888[k] = -3.0 * ii_664[k]
                   + f_0 * li_1140[k];

        t_889[k] = -3.0 * ii_665[k]
                   + f_0 * li_1141[k];

        t_890[k] = -3.0 * ii_666[k]
                   + f_0 * li_1142[k];

        t_891[k] = -3.0 * ii_667[k]
                   + f_0 * li_1143[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, t_895, t_896, ii_668, ii_669, ii_670, ii_671, \
                         ii_672, li_1144, li_1145, li_1146, li_1147, \
                         li_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = -3.0 * ii_668[k]
                   + f_0 * li_1144[k];

        t_893[k] = -3.0 * ii_669[k]
                   + f_0 * li_1145[k];

        t_894[k] = -3.0 * ii_670[k]
                   + f_0 * li_1146[k];

        t_895[k] = -3.0 * ii_671[k]
                   + f_0 * li_1147[k];

        t_896[k] = -4.0 * ii_672[k]
                   + f_0 * li_1148[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, t_900, t_901, ii_673, ii_674, ii_675, ii_676, \
                         ii_677, li_1149, li_1150, li_1151, li_1152, \
                         li_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = -4.0 * ii_673[k]
                   + f_0 * li_1149[k];

        t_898[k] = -4.0 * ii_674[k]
                   + f_0 * li_1150[k];

        t_899[k] = -4.0 * ii_675[k]
                   + f_0 * li_1151[k];

        t_900[k] = -4.0 * ii_676[k]
                   + f_0 * li_1152[k];

        t_901[k] = -4.0 * ii_677[k]
                   + f_0 * li_1153[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, t_905, t_906, ii_678, ii_679, ii_680, ii_681, \
                         ii_682, li_1154, li_1155, li_1156, li_1157, \
                         li_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = -4.0 * ii_678[k]
                   + f_0 * li_1154[k];

        t_903[k] = -4.0 * ii_679[k]
                   + f_0 * li_1155[k];

        t_904[k] = -4.0 * ii_680[k]
                   + f_0 * li_1156[k];

        t_905[k] = -4.0 * ii_681[k]
                   + f_0 * li_1157[k];

        t_906[k] = -4.0 * ii_682[k]
                   + f_0 * li_1158[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, t_910, t_911, ii_683, ii_684, ii_685, ii_686, \
                         ii_687, li_1159, li_1160, li_1161, li_1162, \
                         li_1163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = -4.0 * ii_683[k]
                   + f_0 * li_1159[k];

        t_908[k] = -4.0 * ii_684[k]
                   + f_0 * li_1160[k];

        t_909[k] = -4.0 * ii_685[k]
                   + f_0 * li_1161[k];

        t_910[k] = -4.0 * ii_686[k]
                   + f_0 * li_1162[k];

        t_911[k] = -4.0 * ii_687[k]
                   + f_0 * li_1163[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, t_915, t_916, ii_688, ii_689, ii_690, ii_691, \
                         ii_692, li_1164, li_1165, li_1166, li_1167, \
                         li_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = -4.0 * ii_688[k]
                   + f_0 * li_1164[k];

        t_913[k] = -4.0 * ii_689[k]
                   + f_0 * li_1165[k];

        t_914[k] = -4.0 * ii_690[k]
                   + f_0 * li_1166[k];

        t_915[k] = -4.0 * ii_691[k]
                   + f_0 * li_1167[k];

        t_916[k] = -4.0 * ii_692[k]
                   + f_0 * li_1168[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, t_920, t_921, ii_693, ii_694, ii_695, ii_696, \
                         ii_697, li_1169, li_1170, li_1171, li_1172, \
                         li_1173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = -4.0 * ii_693[k]
                   + f_0 * li_1169[k];

        t_918[k] = -4.0 * ii_694[k]
                   + f_0 * li_1170[k];

        t_919[k] = -4.0 * ii_695[k]
                   + f_0 * li_1171[k];

        t_920[k] = -4.0 * ii_696[k]
                   + f_0 * li_1172[k];

        t_921[k] = -4.0 * ii_697[k]
                   + f_0 * li_1173[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, ii_698, ii_699, ii_700, ii_701, \
                         ii_702, li_1174, li_1175, li_1176, li_1177, \
                         li_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = -4.0 * ii_698[k]
                   + f_0 * li_1174[k];

        t_923[k] = -4.0 * ii_699[k]
                   + f_0 * li_1175[k];

        t_924[k] = -5.0 * ii_700[k]
                   + f_0 * li_1176[k];

        t_925[k] = -5.0 * ii_701[k]
                   + f_0 * li_1177[k];

        t_926[k] = -5.0 * ii_702[k]
                   + f_0 * li_1178[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, t_931, ii_703, ii_704, ii_705, ii_706, \
                         ii_707, li_1179, li_1180, li_1181, li_1182, \
                         li_1183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = -5.0 * ii_703[k]
                   + f_0 * li_1179[k];

        t_928[k] = -5.0 * ii_704[k]
                   + f_0 * li_1180[k];

        t_929[k] = -5.0 * ii_705[k]
                   + f_0 * li_1181[k];

        t_930[k] = -5.0 * ii_706[k]
                   + f_0 * li_1182[k];

        t_931[k] = -5.0 * ii_707[k]
                   + f_0 * li_1183[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, t_935, t_936, ii_708, ii_709, ii_710, ii_711, \
                         ii_712, li_1184, li_1185, li_1186, li_1187, \
                         li_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = -5.0 * ii_708[k]
                   + f_0 * li_1184[k];

        t_933[k] = -5.0 * ii_709[k]
                   + f_0 * li_1185[k];

        t_934[k] = -5.0 * ii_710[k]
                   + f_0 * li_1186[k];

        t_935[k] = -5.0 * ii_711[k]
                   + f_0 * li_1187[k];

        t_936[k] = -5.0 * ii_712[k]
                   + f_0 * li_1188[k];
    }

#pragma omp simd aligned(t_937, t_938, t_939, t_940, t_941, ii_713, ii_714, ii_715, ii_716, \
                         ii_717, li_1189, li_1190, li_1191, li_1192, \
                         li_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_937[k] = -5.0 * ii_713[k]
                   + f_0 * li_1189[k];

        t_938[k] = -5.0 * ii_714[k]
                   + f_0 * li_1190[k];

        t_939[k] = -5.0 * ii_715[k]
                   + f_0 * li_1191[k];

        t_940[k] = -5.0 * ii_716[k]
                   + f_0 * li_1192[k];

        t_941[k] = -5.0 * ii_717[k]
                   + f_0 * li_1193[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, t_946, ii_718, ii_719, ii_720, ii_721, \
                         ii_722, li_1194, li_1195, li_1196, li_1197, \
                         li_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = -5.0 * ii_718[k]
                   + f_0 * li_1194[k];

        t_943[k] = -5.0 * ii_719[k]
                   + f_0 * li_1195[k];

        t_944[k] = -5.0 * ii_720[k]
                   + f_0 * li_1196[k];

        t_945[k] = -5.0 * ii_721[k]
                   + f_0 * li_1197[k];

        t_946[k] = -5.0 * ii_722[k]
                   + f_0 * li_1198[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, t_950, t_951, ii_723, ii_724, ii_725, ii_726, \
                         ii_727, li_1199, li_1200, li_1201, li_1202, \
                         li_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = -5.0 * ii_723[k]
                   + f_0 * li_1199[k];

        t_948[k] = -5.0 * ii_724[k]
                   + f_0 * li_1200[k];

        t_949[k] = -5.0 * ii_725[k]
                   + f_0 * li_1201[k];

        t_950[k] = -5.0 * ii_726[k]
                   + f_0 * li_1202[k];

        t_951[k] = -5.0 * ii_727[k]
                   + f_0 * li_1203[k];
    }

#pragma omp simd aligned(t_952, t_953, t_954, t_955, t_956, ii_728, ii_729, ii_730, ii_731, \
                         ii_732, li_1204, li_1205, li_1206, li_1207, \
                         li_1208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_952[k] = -6.0 * ii_728[k]
                   + f_0 * li_1204[k];

        t_953[k] = -6.0 * ii_729[k]
                   + f_0 * li_1205[k];

        t_954[k] = -6.0 * ii_730[k]
                   + f_0 * li_1206[k];

        t_955[k] = -6.0 * ii_731[k]
                   + f_0 * li_1207[k];

        t_956[k] = -6.0 * ii_732[k]
                   + f_0 * li_1208[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, t_961, ii_733, ii_734, ii_735, ii_736, \
                         ii_737, li_1209, li_1210, li_1211, li_1212, \
                         li_1213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = -6.0 * ii_733[k]
                   + f_0 * li_1209[k];

        t_958[k] = -6.0 * ii_734[k]
                   + f_0 * li_1210[k];

        t_959[k] = -6.0 * ii_735[k]
                   + f_0 * li_1211[k];

        t_960[k] = -6.0 * ii_736[k]
                   + f_0 * li_1212[k];

        t_961[k] = -6.0 * ii_737[k]
                   + f_0 * li_1213[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, t_966, ii_738, ii_739, ii_740, ii_741, \
                         ii_742, li_1214, li_1215, li_1216, li_1217, \
                         li_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = -6.0 * ii_738[k]
                   + f_0 * li_1214[k];

        t_963[k] = -6.0 * ii_739[k]
                   + f_0 * li_1215[k];

        t_964[k] = -6.0 * ii_740[k]
                   + f_0 * li_1216[k];

        t_965[k] = -6.0 * ii_741[k]
                   + f_0 * li_1217[k];

        t_966[k] = -6.0 * ii_742[k]
                   + f_0 * li_1218[k];
    }

#pragma omp simd aligned(t_967, t_968, t_969, t_970, t_971, ii_743, ii_744, ii_745, ii_746, \
                         ii_747, li_1219, li_1220, li_1221, li_1222, \
                         li_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_967[k] = -6.0 * ii_743[k]
                   + f_0 * li_1219[k];

        t_968[k] = -6.0 * ii_744[k]
                   + f_0 * li_1220[k];

        t_969[k] = -6.0 * ii_745[k]
                   + f_0 * li_1221[k];

        t_970[k] = -6.0 * ii_746[k]
                   + f_0 * li_1222[k];

        t_971[k] = -6.0 * ii_747[k]
                   + f_0 * li_1223[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, t_976, ii_748, ii_749, ii_750, ii_751, \
                         ii_752, li_1224, li_1225, li_1226, li_1227, \
                         li_1228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = -6.0 * ii_748[k]
                   + f_0 * li_1224[k];

        t_973[k] = -6.0 * ii_749[k]
                   + f_0 * li_1225[k];

        t_974[k] = -6.0 * ii_750[k]
                   + f_0 * li_1226[k];

        t_975[k] = -6.0 * ii_751[k]
                   + f_0 * li_1227[k];

        t_976[k] = -6.0 * ii_752[k]
                   + f_0 * li_1228[k];
    }
}

static auto
compute_prim_geom_10_ki_electron_repulsion_2_piece6(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t li,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *ii_753 = buffer.data(ii + 753);
    const auto *ii_754 = buffer.data(ii + 754);
    const auto *ii_755 = buffer.data(ii + 755);
    const auto *ii_756 = buffer.data(ii + 756);
    const auto *ii_757 = buffer.data(ii + 757);
    const auto *ii_758 = buffer.data(ii + 758);
    const auto *ii_759 = buffer.data(ii + 759);
    const auto *ii_760 = buffer.data(ii + 760);
    const auto *ii_761 = buffer.data(ii + 761);
    const auto *ii_762 = buffer.data(ii + 762);
    const auto *ii_763 = buffer.data(ii + 763);
    const auto *ii_764 = buffer.data(ii + 764);
    const auto *ii_765 = buffer.data(ii + 765);
    const auto *ii_766 = buffer.data(ii + 766);
    const auto *ii_767 = buffer.data(ii + 767);
    const auto *ii_768 = buffer.data(ii + 768);
    const auto *ii_769 = buffer.data(ii + 769);
    const auto *ii_770 = buffer.data(ii + 770);
    const auto *ii_771 = buffer.data(ii + 771);
    const auto *ii_772 = buffer.data(ii + 772);
    const auto *ii_773 = buffer.data(ii + 773);
    const auto *ii_774 = buffer.data(ii + 774);
    const auto *ii_775 = buffer.data(ii + 775);
    const auto *ii_776 = buffer.data(ii + 776);
    const auto *ii_777 = buffer.data(ii + 777);
    const auto *ii_778 = buffer.data(ii + 778);
    const auto *ii_779 = buffer.data(ii + 779);
    const auto *ii_780 = buffer.data(ii + 780);
    const auto *ii_781 = buffer.data(ii + 781);
    const auto *ii_782 = buffer.data(ii + 782);
    const auto *ii_783 = buffer.data(ii + 783);

    const auto *li_1229 = buffer.data(li + 1229);
    const auto *li_1230 = buffer.data(li + 1230);
    const auto *li_1231 = buffer.data(li + 1231);
    const auto *li_1232 = buffer.data(li + 1232);
    const auto *li_1233 = buffer.data(li + 1233);
    const auto *li_1234 = buffer.data(li + 1234);
    const auto *li_1235 = buffer.data(li + 1235);
    const auto *li_1236 = buffer.data(li + 1236);
    const auto *li_1237 = buffer.data(li + 1237);
    const auto *li_1238 = buffer.data(li + 1238);
    const auto *li_1239 = buffer.data(li + 1239);
    const auto *li_1240 = buffer.data(li + 1240);
    const auto *li_1241 = buffer.data(li + 1241);
    const auto *li_1242 = buffer.data(li + 1242);
    const auto *li_1243 = buffer.data(li + 1243);
    const auto *li_1244 = buffer.data(li + 1244);
    const auto *li_1245 = buffer.data(li + 1245);
    const auto *li_1246 = buffer.data(li + 1246);
    const auto *li_1247 = buffer.data(li + 1247);
    const auto *li_1248 = buffer.data(li + 1248);
    const auto *li_1249 = buffer.data(li + 1249);
    const auto *li_1250 = buffer.data(li + 1250);
    const auto *li_1251 = buffer.data(li + 1251);
    const auto *li_1252 = buffer.data(li + 1252);
    const auto *li_1253 = buffer.data(li + 1253);
    const auto *li_1254 = buffer.data(li + 1254);
    const auto *li_1255 = buffer.data(li + 1255);
    const auto *li_1256 = buffer.data(li + 1256);
    const auto *li_1257 = buffer.data(li + 1257);
    const auto *li_1258 = buffer.data(li + 1258);
    const auto *li_1259 = buffer.data(li + 1259);

#pragma omp simd aligned(t_977, t_978, t_979, t_980, t_981, ii_753, ii_754, ii_755, ii_756, \
                         ii_757, li_1229, li_1230, li_1231, li_1232, \
                         li_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = -6.0 * ii_753[k]
                   + f_0 * li_1229[k];

        t_978[k] = -6.0 * ii_754[k]
                   + f_0 * li_1230[k];

        t_979[k] = -6.0 * ii_755[k]
                   + f_0 * li_1231[k];

        t_980[k] = -7.0 * ii_756[k]
                   + f_0 * li_1232[k];

        t_981[k] = -7.0 * ii_757[k]
                   + f_0 * li_1233[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, t_985, t_986, ii_758, ii_759, ii_760, ii_761, \
                         ii_762, li_1234, li_1235, li_1236, li_1237, \
                         li_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = -7.0 * ii_758[k]
                   + f_0 * li_1234[k];

        t_983[k] = -7.0 * ii_759[k]
                   + f_0 * li_1235[k];

        t_984[k] = -7.0 * ii_760[k]
                   + f_0 * li_1236[k];

        t_985[k] = -7.0 * ii_761[k]
                   + f_0 * li_1237[k];

        t_986[k] = -7.0 * ii_762[k]
                   + f_0 * li_1238[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, t_990, t_991, ii_763, ii_764, ii_765, ii_766, \
                         ii_767, li_1239, li_1240, li_1241, li_1242, \
                         li_1243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = -7.0 * ii_763[k]
                   + f_0 * li_1239[k];

        t_988[k] = -7.0 * ii_764[k]
                   + f_0 * li_1240[k];

        t_989[k] = -7.0 * ii_765[k]
                   + f_0 * li_1241[k];

        t_990[k] = -7.0 * ii_766[k]
                   + f_0 * li_1242[k];

        t_991[k] = -7.0 * ii_767[k]
                   + f_0 * li_1243[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, t_996, ii_768, ii_769, ii_770, ii_771, \
                         ii_772, li_1244, li_1245, li_1246, li_1247, \
                         li_1248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = -7.0 * ii_768[k]
                   + f_0 * li_1244[k];

        t_993[k] = -7.0 * ii_769[k]
                   + f_0 * li_1245[k];

        t_994[k] = -7.0 * ii_770[k]
                   + f_0 * li_1246[k];

        t_995[k] = -7.0 * ii_771[k]
                   + f_0 * li_1247[k];

        t_996[k] = -7.0 * ii_772[k]
                   + f_0 * li_1248[k];
    }

#pragma omp simd aligned(t_997, t_998, t_999, t_1000, t_1001, ii_773, ii_774, ii_775, ii_776, \
                         ii_777, li_1249, li_1250, li_1251, li_1252, \
                         li_1253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_997[k] = -7.0 * ii_773[k]
                   + f_0 * li_1249[k];

        t_998[k] = -7.0 * ii_774[k]
                   + f_0 * li_1250[k];

        t_999[k] = -7.0 * ii_775[k]
                   + f_0 * li_1251[k];

        t_1000[k] = -7.0 * ii_776[k]
                    + f_0 * li_1252[k];

        t_1001[k] = -7.0 * ii_777[k]
                    + f_0 * li_1253[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, t_1005, t_1006, ii_778, ii_779, ii_780, \
                         ii_781, ii_782, li_1254, li_1255, li_1256, li_1257, \
                         li_1258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = -7.0 * ii_778[k]
                    + f_0 * li_1254[k];

        t_1003[k] = -7.0 * ii_779[k]
                    + f_0 * li_1255[k];

        t_1004[k] = -7.0 * ii_780[k]
                    + f_0 * li_1256[k];

        t_1005[k] = -7.0 * ii_781[k]
                    + f_0 * li_1257[k];

        t_1006[k] = -7.0 * ii_782[k]
                    + f_0 * li_1258[k];
    }

#pragma omp simd aligned(t_1007, ii_783, li_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = -7.0 * ii_783[k]
                    + f_0 * li_1259[k];
    }
}

auto
compute_prim_geom_10_ki_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ii, const size_t li,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ki_electron_repulsion_2_piece0(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_2_piece1(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_2_piece2(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_2_piece3(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_2_piece4(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_2_piece5(buffer, target, ii, li, ncols, alpha);

    compute_prim_geom_10_ki_electron_repulsion_2_piece6(buffer, target, ii, li, ncols, alpha);
}

}  // namespace simdt2ceri
