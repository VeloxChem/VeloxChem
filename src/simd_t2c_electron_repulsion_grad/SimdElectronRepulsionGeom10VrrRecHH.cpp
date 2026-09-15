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


#include "SimdElectronRepulsionGeom10VrrRecHH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_hh_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ih,
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

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_71 = buffer.data(gh + 71);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_74 = buffer.data(gh + 74);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_76 = buffer.data(gh + 76);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_115 = buffer.data(gh + 115);
    const auto *gh_116 = buffer.data(gh + 116);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_118 = buffer.data(gh + 118);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_124 = buffer.data(gh + 124);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_137 = buffer.data(gh + 137);
    const auto *gh_138 = buffer.data(gh + 138);
    const auto *gh_139 = buffer.data(gh + 139);
    const auto *gh_140 = buffer.data(gh + 140);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_142 = buffer.data(gh + 142);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_148 = buffer.data(gh + 148);
    const auto *gh_149 = buffer.data(gh + 149);

    const auto *ih_0 = buffer.data(ih + 0);
    const auto *ih_1 = buffer.data(ih + 1);
    const auto *ih_2 = buffer.data(ih + 2);
    const auto *ih_3 = buffer.data(ih + 3);
    const auto *ih_4 = buffer.data(ih + 4);
    const auto *ih_5 = buffer.data(ih + 5);
    const auto *ih_6 = buffer.data(ih + 6);
    const auto *ih_7 = buffer.data(ih + 7);
    const auto *ih_8 = buffer.data(ih + 8);
    const auto *ih_9 = buffer.data(ih + 9);
    const auto *ih_10 = buffer.data(ih + 10);
    const auto *ih_11 = buffer.data(ih + 11);
    const auto *ih_12 = buffer.data(ih + 12);
    const auto *ih_13 = buffer.data(ih + 13);
    const auto *ih_14 = buffer.data(ih + 14);
    const auto *ih_15 = buffer.data(ih + 15);
    const auto *ih_16 = buffer.data(ih + 16);
    const auto *ih_17 = buffer.data(ih + 17);
    const auto *ih_18 = buffer.data(ih + 18);
    const auto *ih_19 = buffer.data(ih + 19);
    const auto *ih_20 = buffer.data(ih + 20);
    const auto *ih_21 = buffer.data(ih + 21);
    const auto *ih_22 = buffer.data(ih + 22);
    const auto *ih_23 = buffer.data(ih + 23);
    const auto *ih_24 = buffer.data(ih + 24);
    const auto *ih_25 = buffer.data(ih + 25);
    const auto *ih_26 = buffer.data(ih + 26);
    const auto *ih_27 = buffer.data(ih + 27);
    const auto *ih_28 = buffer.data(ih + 28);
    const auto *ih_29 = buffer.data(ih + 29);
    const auto *ih_30 = buffer.data(ih + 30);
    const auto *ih_31 = buffer.data(ih + 31);
    const auto *ih_32 = buffer.data(ih + 32);
    const auto *ih_33 = buffer.data(ih + 33);
    const auto *ih_34 = buffer.data(ih + 34);
    const auto *ih_35 = buffer.data(ih + 35);
    const auto *ih_36 = buffer.data(ih + 36);
    const auto *ih_37 = buffer.data(ih + 37);
    const auto *ih_38 = buffer.data(ih + 38);
    const auto *ih_39 = buffer.data(ih + 39);
    const auto *ih_40 = buffer.data(ih + 40);
    const auto *ih_41 = buffer.data(ih + 41);
    const auto *ih_42 = buffer.data(ih + 42);
    const auto *ih_43 = buffer.data(ih + 43);
    const auto *ih_44 = buffer.data(ih + 44);
    const auto *ih_45 = buffer.data(ih + 45);
    const auto *ih_46 = buffer.data(ih + 46);
    const auto *ih_47 = buffer.data(ih + 47);
    const auto *ih_48 = buffer.data(ih + 48);
    const auto *ih_49 = buffer.data(ih + 49);
    const auto *ih_50 = buffer.data(ih + 50);
    const auto *ih_51 = buffer.data(ih + 51);
    const auto *ih_52 = buffer.data(ih + 52);
    const auto *ih_53 = buffer.data(ih + 53);
    const auto *ih_54 = buffer.data(ih + 54);
    const auto *ih_55 = buffer.data(ih + 55);
    const auto *ih_56 = buffer.data(ih + 56);
    const auto *ih_57 = buffer.data(ih + 57);
    const auto *ih_58 = buffer.data(ih + 58);
    const auto *ih_59 = buffer.data(ih + 59);
    const auto *ih_60 = buffer.data(ih + 60);
    const auto *ih_61 = buffer.data(ih + 61);
    const auto *ih_62 = buffer.data(ih + 62);
    const auto *ih_63 = buffer.data(ih + 63);
    const auto *ih_64 = buffer.data(ih + 64);
    const auto *ih_65 = buffer.data(ih + 65);
    const auto *ih_66 = buffer.data(ih + 66);
    const auto *ih_67 = buffer.data(ih + 67);
    const auto *ih_68 = buffer.data(ih + 68);
    const auto *ih_69 = buffer.data(ih + 69);
    const auto *ih_70 = buffer.data(ih + 70);
    const auto *ih_71 = buffer.data(ih + 71);
    const auto *ih_72 = buffer.data(ih + 72);
    const auto *ih_73 = buffer.data(ih + 73);
    const auto *ih_74 = buffer.data(ih + 74);
    const auto *ih_75 = buffer.data(ih + 75);
    const auto *ih_76 = buffer.data(ih + 76);
    const auto *ih_77 = buffer.data(ih + 77);
    const auto *ih_78 = buffer.data(ih + 78);
    const auto *ih_79 = buffer.data(ih + 79);
    const auto *ih_80 = buffer.data(ih + 80);
    const auto *ih_81 = buffer.data(ih + 81);
    const auto *ih_82 = buffer.data(ih + 82);
    const auto *ih_83 = buffer.data(ih + 83);
    const auto *ih_84 = buffer.data(ih + 84);
    const auto *ih_85 = buffer.data(ih + 85);
    const auto *ih_86 = buffer.data(ih + 86);
    const auto *ih_87 = buffer.data(ih + 87);
    const auto *ih_88 = buffer.data(ih + 88);
    const auto *ih_89 = buffer.data(ih + 89);
    const auto *ih_90 = buffer.data(ih + 90);
    const auto *ih_91 = buffer.data(ih + 91);
    const auto *ih_92 = buffer.data(ih + 92);
    const auto *ih_93 = buffer.data(ih + 93);
    const auto *ih_94 = buffer.data(ih + 94);
    const auto *ih_95 = buffer.data(ih + 95);
    const auto *ih_96 = buffer.data(ih + 96);
    const auto *ih_97 = buffer.data(ih + 97);
    const auto *ih_98 = buffer.data(ih + 98);
    const auto *ih_99 = buffer.data(ih + 99);
    const auto *ih_100 = buffer.data(ih + 100);
    const auto *ih_101 = buffer.data(ih + 101);
    const auto *ih_102 = buffer.data(ih + 102);
    const auto *ih_103 = buffer.data(ih + 103);
    const auto *ih_104 = buffer.data(ih + 104);
    const auto *ih_105 = buffer.data(ih + 105);
    const auto *ih_106 = buffer.data(ih + 106);
    const auto *ih_107 = buffer.data(ih + 107);
    const auto *ih_108 = buffer.data(ih + 108);
    const auto *ih_109 = buffer.data(ih + 109);
    const auto *ih_110 = buffer.data(ih + 110);
    const auto *ih_111 = buffer.data(ih + 111);
    const auto *ih_112 = buffer.data(ih + 112);
    const auto *ih_113 = buffer.data(ih + 113);
    const auto *ih_114 = buffer.data(ih + 114);
    const auto *ih_115 = buffer.data(ih + 115);
    const auto *ih_116 = buffer.data(ih + 116);
    const auto *ih_117 = buffer.data(ih + 117);
    const auto *ih_118 = buffer.data(ih + 118);
    const auto *ih_119 = buffer.data(ih + 119);
    const auto *ih_120 = buffer.data(ih + 120);
    const auto *ih_121 = buffer.data(ih + 121);
    const auto *ih_122 = buffer.data(ih + 122);
    const auto *ih_123 = buffer.data(ih + 123);
    const auto *ih_124 = buffer.data(ih + 124);
    const auto *ih_125 = buffer.data(ih + 125);
    const auto *ih_126 = buffer.data(ih + 126);
    const auto *ih_127 = buffer.data(ih + 127);
    const auto *ih_128 = buffer.data(ih + 128);
    const auto *ih_129 = buffer.data(ih + 129);
    const auto *ih_130 = buffer.data(ih + 130);
    const auto *ih_131 = buffer.data(ih + 131);
    const auto *ih_132 = buffer.data(ih + 132);
    const auto *ih_133 = buffer.data(ih + 133);
    const auto *ih_134 = buffer.data(ih + 134);
    const auto *ih_135 = buffer.data(ih + 135);
    const auto *ih_136 = buffer.data(ih + 136);
    const auto *ih_137 = buffer.data(ih + 137);
    const auto *ih_138 = buffer.data(ih + 138);
    const auto *ih_139 = buffer.data(ih + 139);
    const auto *ih_140 = buffer.data(ih + 140);
    const auto *ih_141 = buffer.data(ih + 141);
    const auto *ih_142 = buffer.data(ih + 142);
    const auto *ih_143 = buffer.data(ih + 143);
    const auto *ih_144 = buffer.data(ih + 144);
    const auto *ih_145 = buffer.data(ih + 145);
    const auto *ih_146 = buffer.data(ih + 146);
    const auto *ih_147 = buffer.data(ih + 147);
    const auto *ih_148 = buffer.data(ih + 148);
    const auto *ih_149 = buffer.data(ih + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, gh_0, gh_1, gh_2, gh_3, gh_4, ih_0, ih_1, \
                         ih_2, ih_3, ih_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -5.0 * gh_0[k]
                 + f_0 * ih_0[k];

        t_1[k] = -5.0 * gh_1[k]
                 + f_0 * ih_1[k];

        t_2[k] = -5.0 * gh_2[k]
                 + f_0 * ih_2[k];

        t_3[k] = -5.0 * gh_3[k]
                 + f_0 * ih_3[k];

        t_4[k] = -5.0 * gh_4[k]
                 + f_0 * ih_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, gh_5, gh_6, gh_7, gh_8, gh_9, ih_5, ih_6, \
                         ih_7, ih_8, ih_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -5.0 * gh_5[k]
                 + f_0 * ih_5[k];

        t_6[k] = -5.0 * gh_6[k]
                 + f_0 * ih_6[k];

        t_7[k] = -5.0 * gh_7[k]
                 + f_0 * ih_7[k];

        t_8[k] = -5.0 * gh_8[k]
                 + f_0 * ih_8[k];

        t_9[k] = -5.0 * gh_9[k]
                 + f_0 * ih_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, gh_10, gh_11, gh_12, gh_13, gh_14, \
                         ih_10, ih_11, ih_12, ih_13, ih_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -5.0 * gh_10[k]
                  + f_0 * ih_10[k];

        t_11[k] = -5.0 * gh_11[k]
                  + f_0 * ih_11[k];

        t_12[k] = -5.0 * gh_12[k]
                  + f_0 * ih_12[k];

        t_13[k] = -5.0 * gh_13[k]
                  + f_0 * ih_13[k];

        t_14[k] = -5.0 * gh_14[k]
                  + f_0 * ih_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, gh_15, gh_16, gh_17, gh_18, gh_19, \
                         ih_15, ih_16, ih_17, ih_18, ih_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -5.0 * gh_15[k]
                  + f_0 * ih_15[k];

        t_16[k] = -5.0 * gh_16[k]
                  + f_0 * ih_16[k];

        t_17[k] = -5.0 * gh_17[k]
                  + f_0 * ih_17[k];

        t_18[k] = -5.0 * gh_18[k]
                  + f_0 * ih_18[k];

        t_19[k] = -5.0 * gh_19[k]
                  + f_0 * ih_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, gh_20, gh_21, gh_22, gh_23, gh_24, \
                         ih_20, ih_21, ih_22, ih_23, ih_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -5.0 * gh_20[k]
                  + f_0 * ih_20[k];

        t_21[k] = -4.0 * gh_21[k]
                  + f_0 * ih_21[k];

        t_22[k] = -4.0 * gh_22[k]
                  + f_0 * ih_22[k];

        t_23[k] = -4.0 * gh_23[k]
                  + f_0 * ih_23[k];

        t_24[k] = -4.0 * gh_24[k]
                  + f_0 * ih_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, gh_25, gh_26, gh_27, gh_28, gh_29, \
                         ih_25, ih_26, ih_27, ih_28, ih_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -4.0 * gh_25[k]
                  + f_0 * ih_25[k];

        t_26[k] = -4.0 * gh_26[k]
                  + f_0 * ih_26[k];

        t_27[k] = -4.0 * gh_27[k]
                  + f_0 * ih_27[k];

        t_28[k] = -4.0 * gh_28[k]
                  + f_0 * ih_28[k];

        t_29[k] = -4.0 * gh_29[k]
                  + f_0 * ih_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, gh_30, gh_31, gh_32, gh_33, gh_34, \
                         ih_30, ih_31, ih_32, ih_33, ih_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -4.0 * gh_30[k]
                  + f_0 * ih_30[k];

        t_31[k] = -4.0 * gh_31[k]
                  + f_0 * ih_31[k];

        t_32[k] = -4.0 * gh_32[k]
                  + f_0 * ih_32[k];

        t_33[k] = -4.0 * gh_33[k]
                  + f_0 * ih_33[k];

        t_34[k] = -4.0 * gh_34[k]
                  + f_0 * ih_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gh_35, gh_36, gh_37, gh_38, gh_39, \
                         ih_35, ih_36, ih_37, ih_38, ih_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -4.0 * gh_35[k]
                  + f_0 * ih_35[k];

        t_36[k] = -4.0 * gh_36[k]
                  + f_0 * ih_36[k];

        t_37[k] = -4.0 * gh_37[k]
                  + f_0 * ih_37[k];

        t_38[k] = -4.0 * gh_38[k]
                  + f_0 * ih_38[k];

        t_39[k] = -4.0 * gh_39[k]
                  + f_0 * ih_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gh_40, gh_41, gh_42, gh_43, gh_44, \
                         ih_40, ih_41, ih_42, ih_43, ih_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -4.0 * gh_40[k]
                  + f_0 * ih_40[k];

        t_41[k] = -4.0 * gh_41[k]
                  + f_0 * ih_41[k];

        t_42[k] = -4.0 * gh_42[k]
                  + f_0 * ih_42[k];

        t_43[k] = -4.0 * gh_43[k]
                  + f_0 * ih_43[k];

        t_44[k] = -4.0 * gh_44[k]
                  + f_0 * ih_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gh_45, gh_46, gh_47, gh_48, gh_49, \
                         ih_45, ih_46, ih_47, ih_48, ih_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -4.0 * gh_45[k]
                  + f_0 * ih_45[k];

        t_46[k] = -4.0 * gh_46[k]
                  + f_0 * ih_46[k];

        t_47[k] = -4.0 * gh_47[k]
                  + f_0 * ih_47[k];

        t_48[k] = -4.0 * gh_48[k]
                  + f_0 * ih_48[k];

        t_49[k] = -4.0 * gh_49[k]
                  + f_0 * ih_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gh_50, gh_51, gh_52, gh_53, gh_54, \
                         ih_50, ih_51, ih_52, ih_53, ih_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * gh_50[k]
                  + f_0 * ih_50[k];

        t_51[k] = -4.0 * gh_51[k]
                  + f_0 * ih_51[k];

        t_52[k] = -4.0 * gh_52[k]
                  + f_0 * ih_52[k];

        t_53[k] = -4.0 * gh_53[k]
                  + f_0 * ih_53[k];

        t_54[k] = -4.0 * gh_54[k]
                  + f_0 * ih_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gh_55, gh_56, gh_57, gh_58, gh_59, \
                         ih_55, ih_56, ih_57, ih_58, ih_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -4.0 * gh_55[k]
                  + f_0 * ih_55[k];

        t_56[k] = -4.0 * gh_56[k]
                  + f_0 * ih_56[k];

        t_57[k] = -4.0 * gh_57[k]
                  + f_0 * ih_57[k];

        t_58[k] = -4.0 * gh_58[k]
                  + f_0 * ih_58[k];

        t_59[k] = -4.0 * gh_59[k]
                  + f_0 * ih_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, gh_60, gh_61, gh_62, gh_63, gh_64, \
                         ih_60, ih_61, ih_62, ih_63, ih_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -4.0 * gh_60[k]
                  + f_0 * ih_60[k];

        t_61[k] = -4.0 * gh_61[k]
                  + f_0 * ih_61[k];

        t_62[k] = -4.0 * gh_62[k]
                  + f_0 * ih_62[k];

        t_63[k] = -3.0 * gh_63[k]
                  + f_0 * ih_63[k];

        t_64[k] = -3.0 * gh_64[k]
                  + f_0 * ih_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gh_65, gh_66, gh_67, gh_68, gh_69, \
                         ih_65, ih_66, ih_67, ih_68, ih_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * gh_65[k]
                  + f_0 * ih_65[k];

        t_66[k] = -3.0 * gh_66[k]
                  + f_0 * ih_66[k];

        t_67[k] = -3.0 * gh_67[k]
                  + f_0 * ih_67[k];

        t_68[k] = -3.0 * gh_68[k]
                  + f_0 * ih_68[k];

        t_69[k] = -3.0 * gh_69[k]
                  + f_0 * ih_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gh_70, gh_71, gh_72, gh_73, gh_74, \
                         ih_70, ih_71, ih_72, ih_73, ih_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -3.0 * gh_70[k]
                  + f_0 * ih_70[k];

        t_71[k] = -3.0 * gh_71[k]
                  + f_0 * ih_71[k];

        t_72[k] = -3.0 * gh_72[k]
                  + f_0 * ih_72[k];

        t_73[k] = -3.0 * gh_73[k]
                  + f_0 * ih_73[k];

        t_74[k] = -3.0 * gh_74[k]
                  + f_0 * ih_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, gh_75, gh_76, gh_77, gh_78, gh_79, \
                         ih_75, ih_76, ih_77, ih_78, ih_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -3.0 * gh_75[k]
                  + f_0 * ih_75[k];

        t_76[k] = -3.0 * gh_76[k]
                  + f_0 * ih_76[k];

        t_77[k] = -3.0 * gh_77[k]
                  + f_0 * ih_77[k];

        t_78[k] = -3.0 * gh_78[k]
                  + f_0 * ih_78[k];

        t_79[k] = -3.0 * gh_79[k]
                  + f_0 * ih_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, gh_80, gh_81, gh_82, gh_83, gh_84, \
                         ih_80, ih_81, ih_82, ih_83, ih_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * gh_80[k]
                  + f_0 * ih_80[k];

        t_81[k] = -3.0 * gh_81[k]
                  + f_0 * ih_81[k];

        t_82[k] = -3.0 * gh_82[k]
                  + f_0 * ih_82[k];

        t_83[k] = -3.0 * gh_83[k]
                  + f_0 * ih_83[k];

        t_84[k] = -3.0 * gh_84[k]
                  + f_0 * ih_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gh_85, gh_86, gh_87, gh_88, gh_89, \
                         ih_85, ih_86, ih_87, ih_88, ih_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -3.0 * gh_85[k]
                  + f_0 * ih_85[k];

        t_86[k] = -3.0 * gh_86[k]
                  + f_0 * ih_86[k];

        t_87[k] = -3.0 * gh_87[k]
                  + f_0 * ih_87[k];

        t_88[k] = -3.0 * gh_88[k]
                  + f_0 * ih_88[k];

        t_89[k] = -3.0 * gh_89[k]
                  + f_0 * ih_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, gh_90, gh_91, gh_92, gh_93, gh_94, \
                         ih_90, ih_91, ih_92, ih_93, ih_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * gh_90[k]
                  + f_0 * ih_90[k];

        t_91[k] = -3.0 * gh_91[k]
                  + f_0 * ih_91[k];

        t_92[k] = -3.0 * gh_92[k]
                  + f_0 * ih_92[k];

        t_93[k] = -3.0 * gh_93[k]
                  + f_0 * ih_93[k];

        t_94[k] = -3.0 * gh_94[k]
                  + f_0 * ih_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, gh_95, gh_96, gh_97, gh_98, gh_99, \
                         ih_95, ih_96, ih_97, ih_98, ih_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * gh_95[k]
                  + f_0 * ih_95[k];

        t_96[k] = -3.0 * gh_96[k]
                  + f_0 * ih_96[k];

        t_97[k] = -3.0 * gh_97[k]
                  + f_0 * ih_97[k];

        t_98[k] = -3.0 * gh_98[k]
                  + f_0 * ih_98[k];

        t_99[k] = -3.0 * gh_99[k]
                  + f_0 * ih_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, gh_100, gh_101, gh_102, gh_103, \
                         gh_104, ih_100, ih_101, ih_102, ih_103, \
                         ih_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * gh_100[k]
                   + f_0 * ih_100[k];

        t_101[k] = -3.0 * gh_101[k]
                   + f_0 * ih_101[k];

        t_102[k] = -3.0 * gh_102[k]
                   + f_0 * ih_102[k];

        t_103[k] = -3.0 * gh_103[k]
                   + f_0 * ih_103[k];

        t_104[k] = -3.0 * gh_104[k]
                   + f_0 * ih_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, gh_105, gh_106, gh_107, gh_108, \
                         gh_109, ih_105, ih_106, ih_107, ih_108, \
                         ih_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * gh_105[k]
                   + f_0 * ih_105[k];

        t_106[k] = -3.0 * gh_106[k]
                   + f_0 * ih_106[k];

        t_107[k] = -3.0 * gh_107[k]
                   + f_0 * ih_107[k];

        t_108[k] = -3.0 * gh_108[k]
                   + f_0 * ih_108[k];

        t_109[k] = -3.0 * gh_109[k]
                   + f_0 * ih_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, gh_110, gh_111, gh_112, gh_113, \
                         gh_114, ih_110, ih_111, ih_112, ih_113, \
                         ih_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -3.0 * gh_110[k]
                   + f_0 * ih_110[k];

        t_111[k] = -3.0 * gh_111[k]
                   + f_0 * ih_111[k];

        t_112[k] = -3.0 * gh_112[k]
                   + f_0 * ih_112[k];

        t_113[k] = -3.0 * gh_113[k]
                   + f_0 * ih_113[k];

        t_114[k] = -3.0 * gh_114[k]
                   + f_0 * ih_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, gh_115, gh_116, gh_117, gh_118, \
                         gh_119, ih_115, ih_116, ih_117, ih_118, \
                         ih_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -3.0 * gh_115[k]
                   + f_0 * ih_115[k];

        t_116[k] = -3.0 * gh_116[k]
                   + f_0 * ih_116[k];

        t_117[k] = -3.0 * gh_117[k]
                   + f_0 * ih_117[k];

        t_118[k] = -3.0 * gh_118[k]
                   + f_0 * ih_118[k];

        t_119[k] = -3.0 * gh_119[k]
                   + f_0 * ih_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, gh_120, gh_121, gh_122, gh_123, \
                         gh_124, ih_120, ih_121, ih_122, ih_123, \
                         ih_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -3.0 * gh_120[k]
                   + f_0 * ih_120[k];

        t_121[k] = -3.0 * gh_121[k]
                   + f_0 * ih_121[k];

        t_122[k] = -3.0 * gh_122[k]
                   + f_0 * ih_122[k];

        t_123[k] = -3.0 * gh_123[k]
                   + f_0 * ih_123[k];

        t_124[k] = -3.0 * gh_124[k]
                   + f_0 * ih_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, gh_125, gh_126, gh_127, gh_128, \
                         gh_129, ih_125, ih_126, ih_127, ih_128, \
                         ih_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -3.0 * gh_125[k]
                   + f_0 * ih_125[k];

        t_126[k] = -2.0 * gh_126[k]
                   + f_0 * ih_126[k];

        t_127[k] = -2.0 * gh_127[k]
                   + f_0 * ih_127[k];

        t_128[k] = -2.0 * gh_128[k]
                   + f_0 * ih_128[k];

        t_129[k] = -2.0 * gh_129[k]
                   + f_0 * ih_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, gh_130, gh_131, gh_132, gh_133, \
                         gh_134, ih_130, ih_131, ih_132, ih_133, \
                         ih_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -2.0 * gh_130[k]
                   + f_0 * ih_130[k];

        t_131[k] = -2.0 * gh_131[k]
                   + f_0 * ih_131[k];

        t_132[k] = -2.0 * gh_132[k]
                   + f_0 * ih_132[k];

        t_133[k] = -2.0 * gh_133[k]
                   + f_0 * ih_133[k];

        t_134[k] = -2.0 * gh_134[k]
                   + f_0 * ih_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, gh_135, gh_136, gh_137, gh_138, \
                         gh_139, ih_135, ih_136, ih_137, ih_138, \
                         ih_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -2.0 * gh_135[k]
                   + f_0 * ih_135[k];

        t_136[k] = -2.0 * gh_136[k]
                   + f_0 * ih_136[k];

        t_137[k] = -2.0 * gh_137[k]
                   + f_0 * ih_137[k];

        t_138[k] = -2.0 * gh_138[k]
                   + f_0 * ih_138[k];

        t_139[k] = -2.0 * gh_139[k]
                   + f_0 * ih_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, gh_140, gh_141, gh_142, gh_143, \
                         gh_144, ih_140, ih_141, ih_142, ih_143, \
                         ih_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -2.0 * gh_140[k]
                   + f_0 * ih_140[k];

        t_141[k] = -2.0 * gh_141[k]
                   + f_0 * ih_141[k];

        t_142[k] = -2.0 * gh_142[k]
                   + f_0 * ih_142[k];

        t_143[k] = -2.0 * gh_143[k]
                   + f_0 * ih_143[k];

        t_144[k] = -2.0 * gh_144[k]
                   + f_0 * ih_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, gh_145, gh_146, gh_147, gh_148, \
                         gh_149, ih_145, ih_146, ih_147, ih_148, \
                         ih_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -2.0 * gh_145[k]
                   + f_0 * ih_145[k];

        t_146[k] = -2.0 * gh_146[k]
                   + f_0 * ih_146[k];

        t_147[k] = -2.0 * gh_147[k]
                   + f_0 * ih_147[k];

        t_148[k] = -2.0 * gh_148[k]
                   + f_0 * ih_148[k];

        t_149[k] = -2.0 * gh_149[k]
                   + f_0 * ih_149[k];
    }
}

static auto
compute_prim_geom_10_hh_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ih,
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

    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_190 = buffer.data(gh + 190);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_193 = buffer.data(gh + 193);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_195 = buffer.data(gh + 195);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_197 = buffer.data(gh + 197);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_199 = buffer.data(gh + 199);
    const auto *gh_200 = buffer.data(gh + 200);
    const auto *gh_201 = buffer.data(gh + 201);
    const auto *gh_202 = buffer.data(gh + 202);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_208 = buffer.data(gh + 208);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_212 = buffer.data(gh + 212);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_221 = buffer.data(gh + 221);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_223 = buffer.data(gh + 223);
    const auto *gh_224 = buffer.data(gh + 224);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_231 = buffer.data(gh + 231);
    const auto *gh_232 = buffer.data(gh + 232);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_234 = buffer.data(gh + 234);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_237 = buffer.data(gh + 237);
    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_241 = buffer.data(gh + 241);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_273 = buffer.data(gh + 273);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_275 = buffer.data(gh + 275);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_278 = buffer.data(gh + 278);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);
    const auto *gh_282 = buffer.data(gh + 282);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_287 = buffer.data(gh + 287);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_295 = buffer.data(gh + 295);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_298 = buffer.data(gh + 298);
    const auto *gh_299 = buffer.data(gh + 299);

    const auto *ih_150 = buffer.data(ih + 150);
    const auto *ih_151 = buffer.data(ih + 151);
    const auto *ih_152 = buffer.data(ih + 152);
    const auto *ih_153 = buffer.data(ih + 153);
    const auto *ih_154 = buffer.data(ih + 154);
    const auto *ih_155 = buffer.data(ih + 155);
    const auto *ih_156 = buffer.data(ih + 156);
    const auto *ih_157 = buffer.data(ih + 157);
    const auto *ih_158 = buffer.data(ih + 158);
    const auto *ih_159 = buffer.data(ih + 159);
    const auto *ih_160 = buffer.data(ih + 160);
    const auto *ih_161 = buffer.data(ih + 161);
    const auto *ih_162 = buffer.data(ih + 162);
    const auto *ih_163 = buffer.data(ih + 163);
    const auto *ih_164 = buffer.data(ih + 164);
    const auto *ih_165 = buffer.data(ih + 165);
    const auto *ih_166 = buffer.data(ih + 166);
    const auto *ih_167 = buffer.data(ih + 167);
    const auto *ih_168 = buffer.data(ih + 168);
    const auto *ih_169 = buffer.data(ih + 169);
    const auto *ih_170 = buffer.data(ih + 170);
    const auto *ih_171 = buffer.data(ih + 171);
    const auto *ih_172 = buffer.data(ih + 172);
    const auto *ih_173 = buffer.data(ih + 173);
    const auto *ih_174 = buffer.data(ih + 174);
    const auto *ih_175 = buffer.data(ih + 175);
    const auto *ih_176 = buffer.data(ih + 176);
    const auto *ih_177 = buffer.data(ih + 177);
    const auto *ih_178 = buffer.data(ih + 178);
    const auto *ih_179 = buffer.data(ih + 179);
    const auto *ih_180 = buffer.data(ih + 180);
    const auto *ih_181 = buffer.data(ih + 181);
    const auto *ih_182 = buffer.data(ih + 182);
    const auto *ih_183 = buffer.data(ih + 183);
    const auto *ih_184 = buffer.data(ih + 184);
    const auto *ih_185 = buffer.data(ih + 185);
    const auto *ih_186 = buffer.data(ih + 186);
    const auto *ih_187 = buffer.data(ih + 187);
    const auto *ih_188 = buffer.data(ih + 188);
    const auto *ih_189 = buffer.data(ih + 189);
    const auto *ih_190 = buffer.data(ih + 190);
    const auto *ih_191 = buffer.data(ih + 191);
    const auto *ih_192 = buffer.data(ih + 192);
    const auto *ih_193 = buffer.data(ih + 193);
    const auto *ih_194 = buffer.data(ih + 194);
    const auto *ih_195 = buffer.data(ih + 195);
    const auto *ih_196 = buffer.data(ih + 196);
    const auto *ih_197 = buffer.data(ih + 197);
    const auto *ih_198 = buffer.data(ih + 198);
    const auto *ih_199 = buffer.data(ih + 199);
    const auto *ih_200 = buffer.data(ih + 200);
    const auto *ih_201 = buffer.data(ih + 201);
    const auto *ih_202 = buffer.data(ih + 202);
    const auto *ih_203 = buffer.data(ih + 203);
    const auto *ih_204 = buffer.data(ih + 204);
    const auto *ih_205 = buffer.data(ih + 205);
    const auto *ih_206 = buffer.data(ih + 206);
    const auto *ih_207 = buffer.data(ih + 207);
    const auto *ih_208 = buffer.data(ih + 208);
    const auto *ih_209 = buffer.data(ih + 209);
    const auto *ih_210 = buffer.data(ih + 210);
    const auto *ih_211 = buffer.data(ih + 211);
    const auto *ih_212 = buffer.data(ih + 212);
    const auto *ih_213 = buffer.data(ih + 213);
    const auto *ih_214 = buffer.data(ih + 214);
    const auto *ih_215 = buffer.data(ih + 215);
    const auto *ih_216 = buffer.data(ih + 216);
    const auto *ih_217 = buffer.data(ih + 217);
    const auto *ih_218 = buffer.data(ih + 218);
    const auto *ih_219 = buffer.data(ih + 219);
    const auto *ih_220 = buffer.data(ih + 220);
    const auto *ih_221 = buffer.data(ih + 221);
    const auto *ih_222 = buffer.data(ih + 222);
    const auto *ih_223 = buffer.data(ih + 223);
    const auto *ih_224 = buffer.data(ih + 224);
    const auto *ih_225 = buffer.data(ih + 225);
    const auto *ih_226 = buffer.data(ih + 226);
    const auto *ih_227 = buffer.data(ih + 227);
    const auto *ih_228 = buffer.data(ih + 228);
    const auto *ih_229 = buffer.data(ih + 229);
    const auto *ih_230 = buffer.data(ih + 230);
    const auto *ih_231 = buffer.data(ih + 231);
    const auto *ih_232 = buffer.data(ih + 232);
    const auto *ih_233 = buffer.data(ih + 233);
    const auto *ih_234 = buffer.data(ih + 234);
    const auto *ih_235 = buffer.data(ih + 235);
    const auto *ih_236 = buffer.data(ih + 236);
    const auto *ih_237 = buffer.data(ih + 237);
    const auto *ih_238 = buffer.data(ih + 238);
    const auto *ih_239 = buffer.data(ih + 239);
    const auto *ih_240 = buffer.data(ih + 240);
    const auto *ih_241 = buffer.data(ih + 241);
    const auto *ih_242 = buffer.data(ih + 242);
    const auto *ih_243 = buffer.data(ih + 243);
    const auto *ih_244 = buffer.data(ih + 244);
    const auto *ih_245 = buffer.data(ih + 245);
    const auto *ih_246 = buffer.data(ih + 246);
    const auto *ih_247 = buffer.data(ih + 247);
    const auto *ih_248 = buffer.data(ih + 248);
    const auto *ih_249 = buffer.data(ih + 249);
    const auto *ih_250 = buffer.data(ih + 250);
    const auto *ih_251 = buffer.data(ih + 251);
    const auto *ih_252 = buffer.data(ih + 252);
    const auto *ih_253 = buffer.data(ih + 253);
    const auto *ih_254 = buffer.data(ih + 254);
    const auto *ih_255 = buffer.data(ih + 255);
    const auto *ih_256 = buffer.data(ih + 256);
    const auto *ih_257 = buffer.data(ih + 257);
    const auto *ih_258 = buffer.data(ih + 258);
    const auto *ih_259 = buffer.data(ih + 259);
    const auto *ih_260 = buffer.data(ih + 260);
    const auto *ih_261 = buffer.data(ih + 261);
    const auto *ih_262 = buffer.data(ih + 262);
    const auto *ih_263 = buffer.data(ih + 263);
    const auto *ih_264 = buffer.data(ih + 264);
    const auto *ih_265 = buffer.data(ih + 265);
    const auto *ih_266 = buffer.data(ih + 266);
    const auto *ih_267 = buffer.data(ih + 267);
    const auto *ih_268 = buffer.data(ih + 268);
    const auto *ih_269 = buffer.data(ih + 269);
    const auto *ih_270 = buffer.data(ih + 270);
    const auto *ih_271 = buffer.data(ih + 271);
    const auto *ih_272 = buffer.data(ih + 272);
    const auto *ih_273 = buffer.data(ih + 273);
    const auto *ih_274 = buffer.data(ih + 274);
    const auto *ih_275 = buffer.data(ih + 275);
    const auto *ih_276 = buffer.data(ih + 276);
    const auto *ih_277 = buffer.data(ih + 277);
    const auto *ih_278 = buffer.data(ih + 278);
    const auto *ih_279 = buffer.data(ih + 279);
    const auto *ih_280 = buffer.data(ih + 280);
    const auto *ih_281 = buffer.data(ih + 281);
    const auto *ih_282 = buffer.data(ih + 282);
    const auto *ih_283 = buffer.data(ih + 283);
    const auto *ih_284 = buffer.data(ih + 284);
    const auto *ih_285 = buffer.data(ih + 285);
    const auto *ih_286 = buffer.data(ih + 286);
    const auto *ih_287 = buffer.data(ih + 287);
    const auto *ih_288 = buffer.data(ih + 288);
    const auto *ih_289 = buffer.data(ih + 289);
    const auto *ih_290 = buffer.data(ih + 290);
    const auto *ih_291 = buffer.data(ih + 291);
    const auto *ih_292 = buffer.data(ih + 292);
    const auto *ih_293 = buffer.data(ih + 293);
    const auto *ih_294 = buffer.data(ih + 294);
    const auto *ih_295 = buffer.data(ih + 295);
    const auto *ih_296 = buffer.data(ih + 296);
    const auto *ih_297 = buffer.data(ih + 297);
    const auto *ih_298 = buffer.data(ih + 298);
    const auto *ih_299 = buffer.data(ih + 299);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, gh_150, gh_151, gh_152, gh_153, \
                         gh_154, ih_150, ih_151, ih_152, ih_153, \
                         ih_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -2.0 * gh_150[k]
                   + f_0 * ih_150[k];

        t_151[k] = -2.0 * gh_151[k]
                   + f_0 * ih_151[k];

        t_152[k] = -2.0 * gh_152[k]
                   + f_0 * ih_152[k];

        t_153[k] = -2.0 * gh_153[k]
                   + f_0 * ih_153[k];

        t_154[k] = -2.0 * gh_154[k]
                   + f_0 * ih_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, gh_155, gh_156, gh_157, gh_158, \
                         gh_159, ih_155, ih_156, ih_157, ih_158, \
                         ih_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -2.0 * gh_155[k]
                   + f_0 * ih_155[k];

        t_156[k] = -2.0 * gh_156[k]
                   + f_0 * ih_156[k];

        t_157[k] = -2.0 * gh_157[k]
                   + f_0 * ih_157[k];

        t_158[k] = -2.0 * gh_158[k]
                   + f_0 * ih_158[k];

        t_159[k] = -2.0 * gh_159[k]
                   + f_0 * ih_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, gh_160, gh_161, gh_162, gh_163, \
                         gh_164, ih_160, ih_161, ih_162, ih_163, \
                         ih_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -2.0 * gh_160[k]
                   + f_0 * ih_160[k];

        t_161[k] = -2.0 * gh_161[k]
                   + f_0 * ih_161[k];

        t_162[k] = -2.0 * gh_162[k]
                   + f_0 * ih_162[k];

        t_163[k] = -2.0 * gh_163[k]
                   + f_0 * ih_163[k];

        t_164[k] = -2.0 * gh_164[k]
                   + f_0 * ih_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, gh_165, gh_166, gh_167, gh_168, \
                         gh_169, ih_165, ih_166, ih_167, ih_168, \
                         ih_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -2.0 * gh_165[k]
                   + f_0 * ih_165[k];

        t_166[k] = -2.0 * gh_166[k]
                   + f_0 * ih_166[k];

        t_167[k] = -2.0 * gh_167[k]
                   + f_0 * ih_167[k];

        t_168[k] = -2.0 * gh_168[k]
                   + f_0 * ih_168[k];

        t_169[k] = -2.0 * gh_169[k]
                   + f_0 * ih_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, gh_170, gh_171, gh_172, gh_173, \
                         gh_174, ih_170, ih_171, ih_172, ih_173, \
                         ih_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * gh_170[k]
                   + f_0 * ih_170[k];

        t_171[k] = -2.0 * gh_171[k]
                   + f_0 * ih_171[k];

        t_172[k] = -2.0 * gh_172[k]
                   + f_0 * ih_172[k];

        t_173[k] = -2.0 * gh_173[k]
                   + f_0 * ih_173[k];

        t_174[k] = -2.0 * gh_174[k]
                   + f_0 * ih_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, gh_175, gh_176, gh_177, gh_178, \
                         gh_179, ih_175, ih_176, ih_177, ih_178, \
                         ih_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * gh_175[k]
                   + f_0 * ih_175[k];

        t_176[k] = -2.0 * gh_176[k]
                   + f_0 * ih_176[k];

        t_177[k] = -2.0 * gh_177[k]
                   + f_0 * ih_177[k];

        t_178[k] = -2.0 * gh_178[k]
                   + f_0 * ih_178[k];

        t_179[k] = -2.0 * gh_179[k]
                   + f_0 * ih_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gh_180, gh_181, gh_182, gh_183, \
                         gh_184, ih_180, ih_181, ih_182, ih_183, \
                         ih_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * gh_180[k]
                   + f_0 * ih_180[k];

        t_181[k] = -2.0 * gh_181[k]
                   + f_0 * ih_181[k];

        t_182[k] = -2.0 * gh_182[k]
                   + f_0 * ih_182[k];

        t_183[k] = -2.0 * gh_183[k]
                   + f_0 * ih_183[k];

        t_184[k] = -2.0 * gh_184[k]
                   + f_0 * ih_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gh_185, gh_186, gh_187, gh_188, \
                         gh_189, ih_185, ih_186, ih_187, ih_188, \
                         ih_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * gh_185[k]
                   + f_0 * ih_185[k];

        t_186[k] = -2.0 * gh_186[k]
                   + f_0 * ih_186[k];

        t_187[k] = -2.0 * gh_187[k]
                   + f_0 * ih_187[k];

        t_188[k] = -2.0 * gh_188[k]
                   + f_0 * ih_188[k];

        t_189[k] = -2.0 * gh_189[k]
                   + f_0 * ih_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gh_190, gh_191, gh_192, gh_193, \
                         gh_194, ih_190, ih_191, ih_192, ih_193, \
                         ih_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * gh_190[k]
                   + f_0 * ih_190[k];

        t_191[k] = -2.0 * gh_191[k]
                   + f_0 * ih_191[k];

        t_192[k] = -2.0 * gh_192[k]
                   + f_0 * ih_192[k];

        t_193[k] = -2.0 * gh_193[k]
                   + f_0 * ih_193[k];

        t_194[k] = -2.0 * gh_194[k]
                   + f_0 * ih_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gh_195, gh_196, gh_197, gh_198, \
                         gh_199, ih_195, ih_196, ih_197, ih_198, \
                         ih_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * gh_195[k]
                   + f_0 * ih_195[k];

        t_196[k] = -2.0 * gh_196[k]
                   + f_0 * ih_196[k];

        t_197[k] = -2.0 * gh_197[k]
                   + f_0 * ih_197[k];

        t_198[k] = -2.0 * gh_198[k]
                   + f_0 * ih_198[k];

        t_199[k] = -2.0 * gh_199[k]
                   + f_0 * ih_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gh_200, gh_201, gh_202, gh_203, \
                         gh_204, ih_200, ih_201, ih_202, ih_203, \
                         ih_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * gh_200[k]
                   + f_0 * ih_200[k];

        t_201[k] = -2.0 * gh_201[k]
                   + f_0 * ih_201[k];

        t_202[k] = -2.0 * gh_202[k]
                   + f_0 * ih_202[k];

        t_203[k] = -2.0 * gh_203[k]
                   + f_0 * ih_203[k];

        t_204[k] = -2.0 * gh_204[k]
                   + f_0 * ih_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gh_205, gh_206, gh_207, gh_208, \
                         gh_209, ih_205, ih_206, ih_207, ih_208, \
                         ih_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * gh_205[k]
                   + f_0 * ih_205[k];

        t_206[k] = -2.0 * gh_206[k]
                   + f_0 * ih_206[k];

        t_207[k] = -2.0 * gh_207[k]
                   + f_0 * ih_207[k];

        t_208[k] = -2.0 * gh_208[k]
                   + f_0 * ih_208[k];

        t_209[k] = -2.0 * gh_209[k]
                   + f_0 * ih_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, gh_210, gh_211, gh_212, gh_213, \
                         gh_214, ih_210, ih_211, ih_212, ih_213, \
                         ih_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -gh_210[k]
                   + f_0 * ih_210[k];

        t_211[k] = -gh_211[k]
                   + f_0 * ih_211[k];

        t_212[k] = -gh_212[k]
                   + f_0 * ih_212[k];

        t_213[k] = -gh_213[k]
                   + f_0 * ih_213[k];

        t_214[k] = -gh_214[k]
                   + f_0 * ih_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, gh_215, gh_216, gh_217, gh_218, \
                         gh_219, ih_215, ih_216, ih_217, ih_218, \
                         ih_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -gh_215[k]
                   + f_0 * ih_215[k];

        t_216[k] = -gh_216[k]
                   + f_0 * ih_216[k];

        t_217[k] = -gh_217[k]
                   + f_0 * ih_217[k];

        t_218[k] = -gh_218[k]
                   + f_0 * ih_218[k];

        t_219[k] = -gh_219[k]
                   + f_0 * ih_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, gh_220, gh_221, gh_222, gh_223, \
                         gh_224, ih_220, ih_221, ih_222, ih_223, \
                         ih_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -gh_220[k]
                   + f_0 * ih_220[k];

        t_221[k] = -gh_221[k]
                   + f_0 * ih_221[k];

        t_222[k] = -gh_222[k]
                   + f_0 * ih_222[k];

        t_223[k] = -gh_223[k]
                   + f_0 * ih_223[k];

        t_224[k] = -gh_224[k]
                   + f_0 * ih_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, gh_225, gh_226, gh_227, gh_228, \
                         gh_229, ih_225, ih_226, ih_227, ih_228, \
                         ih_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -gh_225[k]
                   + f_0 * ih_225[k];

        t_226[k] = -gh_226[k]
                   + f_0 * ih_226[k];

        t_227[k] = -gh_227[k]
                   + f_0 * ih_227[k];

        t_228[k] = -gh_228[k]
                   + f_0 * ih_228[k];

        t_229[k] = -gh_229[k]
                   + f_0 * ih_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, gh_230, gh_231, gh_232, gh_233, \
                         gh_234, ih_230, ih_231, ih_232, ih_233, \
                         ih_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -gh_230[k]
                   + f_0 * ih_230[k];

        t_231[k] = -gh_231[k]
                   + f_0 * ih_231[k];

        t_232[k] = -gh_232[k]
                   + f_0 * ih_232[k];

        t_233[k] = -gh_233[k]
                   + f_0 * ih_233[k];

        t_234[k] = -gh_234[k]
                   + f_0 * ih_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, gh_235, gh_236, gh_237, gh_238, \
                         gh_239, ih_235, ih_236, ih_237, ih_238, \
                         ih_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -gh_235[k]
                   + f_0 * ih_235[k];

        t_236[k] = -gh_236[k]
                   + f_0 * ih_236[k];

        t_237[k] = -gh_237[k]
                   + f_0 * ih_237[k];

        t_238[k] = -gh_238[k]
                   + f_0 * ih_238[k];

        t_239[k] = -gh_239[k]
                   + f_0 * ih_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, gh_240, gh_241, gh_242, gh_243, \
                         gh_244, ih_240, ih_241, ih_242, ih_243, \
                         ih_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -gh_240[k]
                   + f_0 * ih_240[k];

        t_241[k] = -gh_241[k]
                   + f_0 * ih_241[k];

        t_242[k] = -gh_242[k]
                   + f_0 * ih_242[k];

        t_243[k] = -gh_243[k]
                   + f_0 * ih_243[k];

        t_244[k] = -gh_244[k]
                   + f_0 * ih_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, gh_245, gh_246, gh_247, gh_248, \
                         gh_249, ih_245, ih_246, ih_247, ih_248, \
                         ih_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -gh_245[k]
                   + f_0 * ih_245[k];

        t_246[k] = -gh_246[k]
                   + f_0 * ih_246[k];

        t_247[k] = -gh_247[k]
                   + f_0 * ih_247[k];

        t_248[k] = -gh_248[k]
                   + f_0 * ih_248[k];

        t_249[k] = -gh_249[k]
                   + f_0 * ih_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, gh_250, gh_251, gh_252, gh_253, \
                         gh_254, ih_250, ih_251, ih_252, ih_253, \
                         ih_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -gh_250[k]
                   + f_0 * ih_250[k];

        t_251[k] = -gh_251[k]
                   + f_0 * ih_251[k];

        t_252[k] = -gh_252[k]
                   + f_0 * ih_252[k];

        t_253[k] = -gh_253[k]
                   + f_0 * ih_253[k];

        t_254[k] = -gh_254[k]
                   + f_0 * ih_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, gh_255, gh_256, gh_257, gh_258, \
                         gh_259, ih_255, ih_256, ih_257, ih_258, \
                         ih_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -gh_255[k]
                   + f_0 * ih_255[k];

        t_256[k] = -gh_256[k]
                   + f_0 * ih_256[k];

        t_257[k] = -gh_257[k]
                   + f_0 * ih_257[k];

        t_258[k] = -gh_258[k]
                   + f_0 * ih_258[k];

        t_259[k] = -gh_259[k]
                   + f_0 * ih_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, gh_260, gh_261, gh_262, gh_263, \
                         gh_264, ih_260, ih_261, ih_262, ih_263, \
                         ih_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -gh_260[k]
                   + f_0 * ih_260[k];

        t_261[k] = -gh_261[k]
                   + f_0 * ih_261[k];

        t_262[k] = -gh_262[k]
                   + f_0 * ih_262[k];

        t_263[k] = -gh_263[k]
                   + f_0 * ih_263[k];

        t_264[k] = -gh_264[k]
                   + f_0 * ih_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, gh_265, gh_266, gh_267, gh_268, \
                         gh_269, ih_265, ih_266, ih_267, ih_268, \
                         ih_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -gh_265[k]
                   + f_0 * ih_265[k];

        t_266[k] = -gh_266[k]
                   + f_0 * ih_266[k];

        t_267[k] = -gh_267[k]
                   + f_0 * ih_267[k];

        t_268[k] = -gh_268[k]
                   + f_0 * ih_268[k];

        t_269[k] = -gh_269[k]
                   + f_0 * ih_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, gh_270, gh_271, gh_272, gh_273, \
                         gh_274, ih_270, ih_271, ih_272, ih_273, \
                         ih_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -gh_270[k]
                   + f_0 * ih_270[k];

        t_271[k] = -gh_271[k]
                   + f_0 * ih_271[k];

        t_272[k] = -gh_272[k]
                   + f_0 * ih_272[k];

        t_273[k] = -gh_273[k]
                   + f_0 * ih_273[k];

        t_274[k] = -gh_274[k]
                   + f_0 * ih_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, gh_275, gh_276, gh_277, gh_278, \
                         gh_279, ih_275, ih_276, ih_277, ih_278, \
                         ih_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -gh_275[k]
                   + f_0 * ih_275[k];

        t_276[k] = -gh_276[k]
                   + f_0 * ih_276[k];

        t_277[k] = -gh_277[k]
                   + f_0 * ih_277[k];

        t_278[k] = -gh_278[k]
                   + f_0 * ih_278[k];

        t_279[k] = -gh_279[k]
                   + f_0 * ih_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, gh_280, gh_281, gh_282, gh_283, \
                         gh_284, ih_280, ih_281, ih_282, ih_283, \
                         ih_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -gh_280[k]
                   + f_0 * ih_280[k];

        t_281[k] = -gh_281[k]
                   + f_0 * ih_281[k];

        t_282[k] = -gh_282[k]
                   + f_0 * ih_282[k];

        t_283[k] = -gh_283[k]
                   + f_0 * ih_283[k];

        t_284[k] = -gh_284[k]
                   + f_0 * ih_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, gh_285, gh_286, gh_287, gh_288, \
                         gh_289, ih_285, ih_286, ih_287, ih_288, \
                         ih_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -gh_285[k]
                   + f_0 * ih_285[k];

        t_286[k] = -gh_286[k]
                   + f_0 * ih_286[k];

        t_287[k] = -gh_287[k]
                   + f_0 * ih_287[k];

        t_288[k] = -gh_288[k]
                   + f_0 * ih_288[k];

        t_289[k] = -gh_289[k]
                   + f_0 * ih_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, gh_290, gh_291, gh_292, gh_293, \
                         gh_294, ih_290, ih_291, ih_292, ih_293, \
                         ih_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -gh_290[k]
                   + f_0 * ih_290[k];

        t_291[k] = -gh_291[k]
                   + f_0 * ih_291[k];

        t_292[k] = -gh_292[k]
                   + f_0 * ih_292[k];

        t_293[k] = -gh_293[k]
                   + f_0 * ih_293[k];

        t_294[k] = -gh_294[k]
                   + f_0 * ih_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, gh_295, gh_296, gh_297, gh_298, \
                         gh_299, ih_295, ih_296, ih_297, ih_298, \
                         ih_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -gh_295[k]
                   + f_0 * ih_295[k];

        t_296[k] = -gh_296[k]
                   + f_0 * ih_296[k];

        t_297[k] = -gh_297[k]
                   + f_0 * ih_297[k];

        t_298[k] = -gh_298[k]
                   + f_0 * ih_298[k];

        t_299[k] = -gh_299[k]
                   + f_0 * ih_299[k];
    }
}

static auto
compute_prim_geom_10_hh_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ih,
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

    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_302 = buffer.data(gh + 302);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_304 = buffer.data(gh + 304);
    const auto *gh_305 = buffer.data(gh + 305);
    const auto *gh_306 = buffer.data(gh + 306);
    const auto *gh_307 = buffer.data(gh + 307);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

    const auto *ih_300 = buffer.data(ih + 300);
    const auto *ih_301 = buffer.data(ih + 301);
    const auto *ih_302 = buffer.data(ih + 302);
    const auto *ih_303 = buffer.data(ih + 303);
    const auto *ih_304 = buffer.data(ih + 304);
    const auto *ih_305 = buffer.data(ih + 305);
    const auto *ih_306 = buffer.data(ih + 306);
    const auto *ih_307 = buffer.data(ih + 307);
    const auto *ih_308 = buffer.data(ih + 308);
    const auto *ih_309 = buffer.data(ih + 309);
    const auto *ih_310 = buffer.data(ih + 310);
    const auto *ih_311 = buffer.data(ih + 311);
    const auto *ih_312 = buffer.data(ih + 312);
    const auto *ih_313 = buffer.data(ih + 313);
    const auto *ih_314 = buffer.data(ih + 314);
    const auto *ih_315 = buffer.data(ih + 315);
    const auto *ih_316 = buffer.data(ih + 316);
    const auto *ih_317 = buffer.data(ih + 317);
    const auto *ih_318 = buffer.data(ih + 318);
    const auto *ih_319 = buffer.data(ih + 319);
    const auto *ih_320 = buffer.data(ih + 320);
    const auto *ih_321 = buffer.data(ih + 321);
    const auto *ih_322 = buffer.data(ih + 322);
    const auto *ih_323 = buffer.data(ih + 323);
    const auto *ih_324 = buffer.data(ih + 324);
    const auto *ih_325 = buffer.data(ih + 325);
    const auto *ih_326 = buffer.data(ih + 326);
    const auto *ih_327 = buffer.data(ih + 327);
    const auto *ih_328 = buffer.data(ih + 328);
    const auto *ih_329 = buffer.data(ih + 329);
    const auto *ih_330 = buffer.data(ih + 330);
    const auto *ih_331 = buffer.data(ih + 331);
    const auto *ih_332 = buffer.data(ih + 332);
    const auto *ih_333 = buffer.data(ih + 333);
    const auto *ih_334 = buffer.data(ih + 334);
    const auto *ih_335 = buffer.data(ih + 335);
    const auto *ih_336 = buffer.data(ih + 336);
    const auto *ih_337 = buffer.data(ih + 337);
    const auto *ih_338 = buffer.data(ih + 338);
    const auto *ih_339 = buffer.data(ih + 339);
    const auto *ih_340 = buffer.data(ih + 340);
    const auto *ih_341 = buffer.data(ih + 341);
    const auto *ih_342 = buffer.data(ih + 342);
    const auto *ih_343 = buffer.data(ih + 343);
    const auto *ih_344 = buffer.data(ih + 344);
    const auto *ih_345 = buffer.data(ih + 345);
    const auto *ih_346 = buffer.data(ih + 346);
    const auto *ih_347 = buffer.data(ih + 347);
    const auto *ih_348 = buffer.data(ih + 348);
    const auto *ih_349 = buffer.data(ih + 349);
    const auto *ih_350 = buffer.data(ih + 350);
    const auto *ih_351 = buffer.data(ih + 351);
    const auto *ih_352 = buffer.data(ih + 352);
    const auto *ih_353 = buffer.data(ih + 353);
    const auto *ih_354 = buffer.data(ih + 354);
    const auto *ih_355 = buffer.data(ih + 355);
    const auto *ih_356 = buffer.data(ih + 356);
    const auto *ih_357 = buffer.data(ih + 357);
    const auto *ih_358 = buffer.data(ih + 358);
    const auto *ih_359 = buffer.data(ih + 359);
    const auto *ih_360 = buffer.data(ih + 360);
    const auto *ih_361 = buffer.data(ih + 361);
    const auto *ih_362 = buffer.data(ih + 362);
    const auto *ih_363 = buffer.data(ih + 363);
    const auto *ih_364 = buffer.data(ih + 364);
    const auto *ih_365 = buffer.data(ih + 365);
    const auto *ih_366 = buffer.data(ih + 366);
    const auto *ih_367 = buffer.data(ih + 367);
    const auto *ih_368 = buffer.data(ih + 368);
    const auto *ih_369 = buffer.data(ih + 369);
    const auto *ih_370 = buffer.data(ih + 370);
    const auto *ih_371 = buffer.data(ih + 371);
    const auto *ih_372 = buffer.data(ih + 372);
    const auto *ih_373 = buffer.data(ih + 373);
    const auto *ih_374 = buffer.data(ih + 374);
    const auto *ih_375 = buffer.data(ih + 375);
    const auto *ih_376 = buffer.data(ih + 376);
    const auto *ih_377 = buffer.data(ih + 377);
    const auto *ih_378 = buffer.data(ih + 378);
    const auto *ih_379 = buffer.data(ih + 379);
    const auto *ih_380 = buffer.data(ih + 380);
    const auto *ih_381 = buffer.data(ih + 381);
    const auto *ih_382 = buffer.data(ih + 382);
    const auto *ih_383 = buffer.data(ih + 383);
    const auto *ih_384 = buffer.data(ih + 384);
    const auto *ih_385 = buffer.data(ih + 385);
    const auto *ih_386 = buffer.data(ih + 386);
    const auto *ih_387 = buffer.data(ih + 387);
    const auto *ih_388 = buffer.data(ih + 388);
    const auto *ih_389 = buffer.data(ih + 389);
    const auto *ih_390 = buffer.data(ih + 390);
    const auto *ih_391 = buffer.data(ih + 391);
    const auto *ih_392 = buffer.data(ih + 392);
    const auto *ih_393 = buffer.data(ih + 393);
    const auto *ih_394 = buffer.data(ih + 394);
    const auto *ih_395 = buffer.data(ih + 395);
    const auto *ih_396 = buffer.data(ih + 396);
    const auto *ih_397 = buffer.data(ih + 397);
    const auto *ih_398 = buffer.data(ih + 398);
    const auto *ih_399 = buffer.data(ih + 399);
    const auto *ih_400 = buffer.data(ih + 400);
    const auto *ih_401 = buffer.data(ih + 401);
    const auto *ih_402 = buffer.data(ih + 402);
    const auto *ih_403 = buffer.data(ih + 403);
    const auto *ih_404 = buffer.data(ih + 404);
    const auto *ih_405 = buffer.data(ih + 405);
    const auto *ih_406 = buffer.data(ih + 406);
    const auto *ih_407 = buffer.data(ih + 407);
    const auto *ih_408 = buffer.data(ih + 408);
    const auto *ih_409 = buffer.data(ih + 409);
    const auto *ih_410 = buffer.data(ih + 410);
    const auto *ih_411 = buffer.data(ih + 411);
    const auto *ih_412 = buffer.data(ih + 412);
    const auto *ih_413 = buffer.data(ih + 413);
    const auto *ih_414 = buffer.data(ih + 414);
    const auto *ih_415 = buffer.data(ih + 415);
    const auto *ih_416 = buffer.data(ih + 416);
    const auto *ih_417 = buffer.data(ih + 417);
    const auto *ih_418 = buffer.data(ih + 418);
    const auto *ih_419 = buffer.data(ih + 419);
    const auto *ih_420 = buffer.data(ih + 420);
    const auto *ih_421 = buffer.data(ih + 421);
    const auto *ih_422 = buffer.data(ih + 422);
    const auto *ih_423 = buffer.data(ih + 423);
    const auto *ih_424 = buffer.data(ih + 424);
    const auto *ih_425 = buffer.data(ih + 425);
    const auto *ih_426 = buffer.data(ih + 426);
    const auto *ih_427 = buffer.data(ih + 427);
    const auto *ih_428 = buffer.data(ih + 428);
    const auto *ih_429 = buffer.data(ih + 429);
    const auto *ih_430 = buffer.data(ih + 430);
    const auto *ih_431 = buffer.data(ih + 431);
    const auto *ih_432 = buffer.data(ih + 432);
    const auto *ih_433 = buffer.data(ih + 433);
    const auto *ih_434 = buffer.data(ih + 434);
    const auto *ih_435 = buffer.data(ih + 435);
    const auto *ih_436 = buffer.data(ih + 436);
    const auto *ih_437 = buffer.data(ih + 437);
    const auto *ih_438 = buffer.data(ih + 438);
    const auto *ih_439 = buffer.data(ih + 439);
    const auto *ih_440 = buffer.data(ih + 440);

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, gh_300, gh_301, gh_302, gh_303, \
                         gh_304, ih_300, ih_301, ih_302, ih_303, \
                         ih_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -gh_300[k]
                   + f_0 * ih_300[k];

        t_301[k] = -gh_301[k]
                   + f_0 * ih_301[k];

        t_302[k] = -gh_302[k]
                   + f_0 * ih_302[k];

        t_303[k] = -gh_303[k]
                   + f_0 * ih_303[k];

        t_304[k] = -gh_304[k]
                   + f_0 * ih_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, gh_305, gh_306, gh_307, gh_308, \
                         gh_309, ih_305, ih_306, ih_307, ih_308, \
                         ih_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -gh_305[k]
                   + f_0 * ih_305[k];

        t_306[k] = -gh_306[k]
                   + f_0 * ih_306[k];

        t_307[k] = -gh_307[k]
                   + f_0 * ih_307[k];

        t_308[k] = -gh_308[k]
                   + f_0 * ih_308[k];

        t_309[k] = -gh_309[k]
                   + f_0 * ih_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, gh_310, gh_311, gh_312, gh_313, \
                         gh_314, ih_310, ih_311, ih_312, ih_313, \
                         ih_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -gh_310[k]
                   + f_0 * ih_310[k];

        t_311[k] = -gh_311[k]
                   + f_0 * ih_311[k];

        t_312[k] = -gh_312[k]
                   + f_0 * ih_312[k];

        t_313[k] = -gh_313[k]
                   + f_0 * ih_313[k];

        t_314[k] = -gh_314[k]
                   + f_0 * ih_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, t_321, t_322, ih_315, \
                         ih_316, ih_317, ih_318, ih_319, ih_320, ih_321, \
                         ih_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_0 * ih_315[k];

        t_316[k] = f_0 * ih_316[k];

        t_317[k] = f_0 * ih_317[k];

        t_318[k] = f_0 * ih_318[k];

        t_319[k] = f_0 * ih_319[k];

        t_320[k] = f_0 * ih_320[k];

        t_321[k] = f_0 * ih_321[k];

        t_322[k] = f_0 * ih_322[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, t_328, t_329, t_330, ih_323, \
                         ih_324, ih_325, ih_326, ih_327, ih_328, ih_329, \
                         ih_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_0 * ih_323[k];

        t_324[k] = f_0 * ih_324[k];

        t_325[k] = f_0 * ih_325[k];

        t_326[k] = f_0 * ih_326[k];

        t_327[k] = f_0 * ih_327[k];

        t_328[k] = f_0 * ih_328[k];

        t_329[k] = f_0 * ih_329[k];

        t_330[k] = f_0 * ih_330[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, t_336, t_337, t_338, ih_331, \
                         ih_332, ih_333, ih_334, ih_335, ih_336, ih_337, \
                         ih_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_0 * ih_331[k];

        t_332[k] = f_0 * ih_332[k];

        t_333[k] = f_0 * ih_333[k];

        t_334[k] = f_0 * ih_334[k];

        t_335[k] = f_0 * ih_335[k];

        t_336[k] = f_0 * ih_336[k];

        t_337[k] = f_0 * ih_337[k];

        t_338[k] = f_0 * ih_338[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, t_344, t_345, t_346, ih_339, \
                         ih_340, ih_341, ih_342, ih_343, ih_344, ih_345, \
                         ih_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_0 * ih_339[k];

        t_340[k] = f_0 * ih_340[k];

        t_341[k] = f_0 * ih_341[k];

        t_342[k] = f_0 * ih_342[k];

        t_343[k] = f_0 * ih_343[k];

        t_344[k] = f_0 * ih_344[k];

        t_345[k] = f_0 * ih_345[k];

        t_346[k] = f_0 * ih_346[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, t_352, t_353, t_354, ih_347, \
                         ih_348, ih_349, ih_350, ih_351, ih_352, ih_353, \
                         ih_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_0 * ih_347[k];

        t_348[k] = f_0 * ih_348[k];

        t_349[k] = f_0 * ih_349[k];

        t_350[k] = f_0 * ih_350[k];

        t_351[k] = f_0 * ih_351[k];

        t_352[k] = f_0 * ih_352[k];

        t_353[k] = f_0 * ih_353[k];

        t_354[k] = f_0 * ih_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, t_360, t_361, t_362, ih_355, \
                         ih_356, ih_357, ih_358, ih_359, ih_360, ih_361, \
                         ih_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_0 * ih_355[k];

        t_356[k] = f_0 * ih_356[k];

        t_357[k] = f_0 * ih_357[k];

        t_358[k] = f_0 * ih_358[k];

        t_359[k] = f_0 * ih_359[k];

        t_360[k] = f_0 * ih_360[k];

        t_361[k] = f_0 * ih_361[k];

        t_362[k] = f_0 * ih_362[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, t_368, t_369, t_370, ih_363, \
                         ih_364, ih_365, ih_366, ih_367, ih_368, ih_369, \
                         ih_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_0 * ih_363[k];

        t_364[k] = f_0 * ih_364[k];

        t_365[k] = f_0 * ih_365[k];

        t_366[k] = f_0 * ih_366[k];

        t_367[k] = f_0 * ih_367[k];

        t_368[k] = f_0 * ih_368[k];

        t_369[k] = f_0 * ih_369[k];

        t_370[k] = f_0 * ih_370[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, t_375, t_376, t_377, t_378, ih_371, \
                         ih_372, ih_373, ih_374, ih_375, ih_376, ih_377, \
                         ih_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_0 * ih_371[k];

        t_372[k] = f_0 * ih_372[k];

        t_373[k] = f_0 * ih_373[k];

        t_374[k] = f_0 * ih_374[k];

        t_375[k] = f_0 * ih_375[k];

        t_376[k] = f_0 * ih_376[k];

        t_377[k] = f_0 * ih_377[k];

        t_378[k] = f_0 * ih_378[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, t_384, t_385, t_386, ih_379, \
                         ih_380, ih_381, ih_382, ih_383, ih_384, ih_385, \
                         ih_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_0 * ih_379[k];

        t_380[k] = f_0 * ih_380[k];

        t_381[k] = f_0 * ih_381[k];

        t_382[k] = f_0 * ih_382[k];

        t_383[k] = f_0 * ih_383[k];

        t_384[k] = f_0 * ih_384[k];

        t_385[k] = f_0 * ih_385[k];

        t_386[k] = f_0 * ih_386[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, t_392, t_393, t_394, ih_387, \
                         ih_388, ih_389, ih_390, ih_391, ih_392, ih_393, \
                         ih_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_0 * ih_387[k];

        t_388[k] = f_0 * ih_388[k];

        t_389[k] = f_0 * ih_389[k];

        t_390[k] = f_0 * ih_390[k];

        t_391[k] = f_0 * ih_391[k];

        t_392[k] = f_0 * ih_392[k];

        t_393[k] = f_0 * ih_393[k];

        t_394[k] = f_0 * ih_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, t_400, t_401, t_402, ih_395, \
                         ih_396, ih_397, ih_398, ih_399, ih_400, ih_401, \
                         ih_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_0 * ih_395[k];

        t_396[k] = f_0 * ih_396[k];

        t_397[k] = f_0 * ih_397[k];

        t_398[k] = f_0 * ih_398[k];

        t_399[k] = f_0 * ih_399[k];

        t_400[k] = f_0 * ih_400[k];

        t_401[k] = f_0 * ih_401[k];

        t_402[k] = f_0 * ih_402[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, t_407, t_408, t_409, t_410, ih_403, \
                         ih_404, ih_405, ih_406, ih_407, ih_408, ih_409, \
                         ih_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_0 * ih_403[k];

        t_404[k] = f_0 * ih_404[k];

        t_405[k] = f_0 * ih_405[k];

        t_406[k] = f_0 * ih_406[k];

        t_407[k] = f_0 * ih_407[k];

        t_408[k] = f_0 * ih_408[k];

        t_409[k] = f_0 * ih_409[k];

        t_410[k] = f_0 * ih_410[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, t_415, t_416, t_417, t_418, ih_411, \
                         ih_412, ih_413, ih_414, ih_415, ih_416, ih_417, \
                         ih_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_0 * ih_411[k];

        t_412[k] = f_0 * ih_412[k];

        t_413[k] = f_0 * ih_413[k];

        t_414[k] = f_0 * ih_414[k];

        t_415[k] = f_0 * ih_415[k];

        t_416[k] = f_0 * ih_416[k];

        t_417[k] = f_0 * ih_417[k];

        t_418[k] = f_0 * ih_418[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, t_424, t_425, t_426, ih_419, \
                         ih_420, ih_421, ih_422, ih_423, ih_424, ih_425, \
                         ih_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_0 * ih_419[k];

        t_420[k] = f_0 * ih_420[k];

        t_421[k] = f_0 * ih_421[k];

        t_422[k] = f_0 * ih_422[k];

        t_423[k] = f_0 * ih_423[k];

        t_424[k] = f_0 * ih_424[k];

        t_425[k] = f_0 * ih_425[k];

        t_426[k] = f_0 * ih_426[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, t_431, t_432, t_433, t_434, ih_427, \
                         ih_428, ih_429, ih_430, ih_431, ih_432, ih_433, \
                         ih_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_0 * ih_427[k];

        t_428[k] = f_0 * ih_428[k];

        t_429[k] = f_0 * ih_429[k];

        t_430[k] = f_0 * ih_430[k];

        t_431[k] = f_0 * ih_431[k];

        t_432[k] = f_0 * ih_432[k];

        t_433[k] = f_0 * ih_433[k];

        t_434[k] = f_0 * ih_434[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, t_440, ih_435, ih_436, ih_437, \
                         ih_438, ih_439, ih_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_0 * ih_435[k];

        t_436[k] = f_0 * ih_436[k];

        t_437[k] = f_0 * ih_437[k];

        t_438[k] = f_0 * ih_438[k];

        t_439[k] = f_0 * ih_439[k];

        t_440[k] = f_0 * ih_440[k];
    }
}

auto
compute_prim_geom_10_hh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t gh, const size_t ih,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hh_electron_repulsion_0_piece0(buffer, target, gh, ih, ncols, alpha);

    compute_prim_geom_10_hh_electron_repulsion_0_piece1(buffer, target, gh, ih, ncols, alpha);

    compute_prim_geom_10_hh_electron_repulsion_0_piece2(buffer, target, gh, ih, ncols, alpha);
}

static auto
compute_prim_geom_10_hh_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ih,
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

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_71 = buffer.data(gh + 71);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_74 = buffer.data(gh + 74);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_76 = buffer.data(gh + 76);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);

    const auto *ih_21 = buffer.data(ih + 21);
    const auto *ih_22 = buffer.data(ih + 22);
    const auto *ih_23 = buffer.data(ih + 23);
    const auto *ih_24 = buffer.data(ih + 24);
    const auto *ih_25 = buffer.data(ih + 25);
    const auto *ih_26 = buffer.data(ih + 26);
    const auto *ih_27 = buffer.data(ih + 27);
    const auto *ih_28 = buffer.data(ih + 28);
    const auto *ih_29 = buffer.data(ih + 29);
    const auto *ih_30 = buffer.data(ih + 30);
    const auto *ih_31 = buffer.data(ih + 31);
    const auto *ih_32 = buffer.data(ih + 32);
    const auto *ih_33 = buffer.data(ih + 33);
    const auto *ih_34 = buffer.data(ih + 34);
    const auto *ih_35 = buffer.data(ih + 35);
    const auto *ih_36 = buffer.data(ih + 36);
    const auto *ih_37 = buffer.data(ih + 37);
    const auto *ih_38 = buffer.data(ih + 38);
    const auto *ih_39 = buffer.data(ih + 39);
    const auto *ih_40 = buffer.data(ih + 40);
    const auto *ih_41 = buffer.data(ih + 41);
    const auto *ih_63 = buffer.data(ih + 63);
    const auto *ih_64 = buffer.data(ih + 64);
    const auto *ih_65 = buffer.data(ih + 65);
    const auto *ih_66 = buffer.data(ih + 66);
    const auto *ih_67 = buffer.data(ih + 67);
    const auto *ih_68 = buffer.data(ih + 68);
    const auto *ih_69 = buffer.data(ih + 69);
    const auto *ih_70 = buffer.data(ih + 70);
    const auto *ih_71 = buffer.data(ih + 71);
    const auto *ih_72 = buffer.data(ih + 72);
    const auto *ih_73 = buffer.data(ih + 73);
    const auto *ih_74 = buffer.data(ih + 74);
    const auto *ih_75 = buffer.data(ih + 75);
    const auto *ih_76 = buffer.data(ih + 76);
    const auto *ih_77 = buffer.data(ih + 77);
    const auto *ih_78 = buffer.data(ih + 78);
    const auto *ih_79 = buffer.data(ih + 79);
    const auto *ih_80 = buffer.data(ih + 80);
    const auto *ih_81 = buffer.data(ih + 81);
    const auto *ih_82 = buffer.data(ih + 82);
    const auto *ih_83 = buffer.data(ih + 83);
    const auto *ih_84 = buffer.data(ih + 84);
    const auto *ih_85 = buffer.data(ih + 85);
    const auto *ih_86 = buffer.data(ih + 86);
    const auto *ih_87 = buffer.data(ih + 87);
    const auto *ih_88 = buffer.data(ih + 88);
    const auto *ih_89 = buffer.data(ih + 89);
    const auto *ih_90 = buffer.data(ih + 90);
    const auto *ih_91 = buffer.data(ih + 91);
    const auto *ih_92 = buffer.data(ih + 92);
    const auto *ih_93 = buffer.data(ih + 93);
    const auto *ih_94 = buffer.data(ih + 94);
    const auto *ih_95 = buffer.data(ih + 95);
    const auto *ih_96 = buffer.data(ih + 96);
    const auto *ih_97 = buffer.data(ih + 97);
    const auto *ih_98 = buffer.data(ih + 98);
    const auto *ih_99 = buffer.data(ih + 99);
    const auto *ih_100 = buffer.data(ih + 100);
    const auto *ih_101 = buffer.data(ih + 101);
    const auto *ih_102 = buffer.data(ih + 102);
    const auto *ih_103 = buffer.data(ih + 103);
    const auto *ih_104 = buffer.data(ih + 104);
    const auto *ih_126 = buffer.data(ih + 126);
    const auto *ih_127 = buffer.data(ih + 127);
    const auto *ih_128 = buffer.data(ih + 128);
    const auto *ih_129 = buffer.data(ih + 129);
    const auto *ih_130 = buffer.data(ih + 130);
    const auto *ih_131 = buffer.data(ih + 131);
    const auto *ih_132 = buffer.data(ih + 132);
    const auto *ih_133 = buffer.data(ih + 133);
    const auto *ih_134 = buffer.data(ih + 134);
    const auto *ih_135 = buffer.data(ih + 135);
    const auto *ih_136 = buffer.data(ih + 136);
    const auto *ih_137 = buffer.data(ih + 137);
    const auto *ih_138 = buffer.data(ih + 138);
    const auto *ih_139 = buffer.data(ih + 139);
    const auto *ih_140 = buffer.data(ih + 140);
    const auto *ih_141 = buffer.data(ih + 141);
    const auto *ih_142 = buffer.data(ih + 142);
    const auto *ih_143 = buffer.data(ih + 143);
    const auto *ih_144 = buffer.data(ih + 144);
    const auto *ih_145 = buffer.data(ih + 145);
    const auto *ih_146 = buffer.data(ih + 146);
    const auto *ih_147 = buffer.data(ih + 147);
    const auto *ih_148 = buffer.data(ih + 148);
    const auto *ih_149 = buffer.data(ih + 149);
    const auto *ih_150 = buffer.data(ih + 150);
    const auto *ih_151 = buffer.data(ih + 151);
    const auto *ih_152 = buffer.data(ih + 152);
    const auto *ih_153 = buffer.data(ih + 153);
    const auto *ih_154 = buffer.data(ih + 154);
    const auto *ih_155 = buffer.data(ih + 155);
    const auto *ih_156 = buffer.data(ih + 156);
    const auto *ih_157 = buffer.data(ih + 157);
    const auto *ih_158 = buffer.data(ih + 158);
    const auto *ih_159 = buffer.data(ih + 159);
    const auto *ih_160 = buffer.data(ih + 160);
    const auto *ih_161 = buffer.data(ih + 161);
    const auto *ih_162 = buffer.data(ih + 162);
    const auto *ih_163 = buffer.data(ih + 163);
    const auto *ih_164 = buffer.data(ih + 164);
    const auto *ih_165 = buffer.data(ih + 165);
    const auto *ih_166 = buffer.data(ih + 166);
    const auto *ih_167 = buffer.data(ih + 167);
    const auto *ih_168 = buffer.data(ih + 168);
    const auto *ih_169 = buffer.data(ih + 169);
    const auto *ih_170 = buffer.data(ih + 170);
    const auto *ih_171 = buffer.data(ih + 171);
    const auto *ih_172 = buffer.data(ih + 172);
    const auto *ih_173 = buffer.data(ih + 173);
    const auto *ih_174 = buffer.data(ih + 174);
    const auto *ih_175 = buffer.data(ih + 175);
    const auto *ih_176 = buffer.data(ih + 176);
    const auto *ih_177 = buffer.data(ih + 177);
    const auto *ih_178 = buffer.data(ih + 178);
    const auto *ih_179 = buffer.data(ih + 179);
    const auto *ih_180 = buffer.data(ih + 180);
    const auto *ih_181 = buffer.data(ih + 181);
    const auto *ih_182 = buffer.data(ih + 182);
    const auto *ih_183 = buffer.data(ih + 183);
    const auto *ih_184 = buffer.data(ih + 184);
    const auto *ih_185 = buffer.data(ih + 185);
    const auto *ih_186 = buffer.data(ih + 186);
    const auto *ih_187 = buffer.data(ih + 187);
    const auto *ih_188 = buffer.data(ih + 188);
    const auto *ih_210 = buffer.data(ih + 210);
    const auto *ih_211 = buffer.data(ih + 211);
    const auto *ih_212 = buffer.data(ih + 212);
    const auto *ih_213 = buffer.data(ih + 213);
    const auto *ih_214 = buffer.data(ih + 214);
    const auto *ih_215 = buffer.data(ih + 215);
    const auto *ih_216 = buffer.data(ih + 216);
    const auto *ih_217 = buffer.data(ih + 217);
    const auto *ih_218 = buffer.data(ih + 218);
    const auto *ih_219 = buffer.data(ih + 219);
    const auto *ih_220 = buffer.data(ih + 220);
    const auto *ih_221 = buffer.data(ih + 221);
    const auto *ih_222 = buffer.data(ih + 222);
    const auto *ih_223 = buffer.data(ih + 223);
    const auto *ih_224 = buffer.data(ih + 224);
    const auto *ih_225 = buffer.data(ih + 225);
    const auto *ih_226 = buffer.data(ih + 226);
    const auto *ih_227 = buffer.data(ih + 227);
    const auto *ih_228 = buffer.data(ih + 228);
    const auto *ih_229 = buffer.data(ih + 229);
    const auto *ih_230 = buffer.data(ih + 230);
    const auto *ih_231 = buffer.data(ih + 231);
    const auto *ih_232 = buffer.data(ih + 232);
    const auto *ih_233 = buffer.data(ih + 233);
    const auto *ih_234 = buffer.data(ih + 234);
    const auto *ih_235 = buffer.data(ih + 235);
    const auto *ih_236 = buffer.data(ih + 236);
    const auto *ih_237 = buffer.data(ih + 237);
    const auto *ih_238 = buffer.data(ih + 238);
    const auto *ih_239 = buffer.data(ih + 239);
    const auto *ih_240 = buffer.data(ih + 240);
    const auto *ih_241 = buffer.data(ih + 241);
    const auto *ih_242 = buffer.data(ih + 242);
    const auto *ih_243 = buffer.data(ih + 243);
    const auto *ih_244 = buffer.data(ih + 244);
    const auto *ih_245 = buffer.data(ih + 245);
    const auto *ih_246 = buffer.data(ih + 246);
    const auto *ih_247 = buffer.data(ih + 247);
    const auto *ih_248 = buffer.data(ih + 248);
    const auto *ih_249 = buffer.data(ih + 249);
    const auto *ih_250 = buffer.data(ih + 250);
    const auto *ih_251 = buffer.data(ih + 251);
    const auto *ih_252 = buffer.data(ih + 252);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ih_21, ih_22, ih_23, ih_24, \
                         ih_25, ih_26, ih_27, ih_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ih_21[k];

        t_1[k] = f_0 * ih_22[k];

        t_2[k] = f_0 * ih_23[k];

        t_3[k] = f_0 * ih_24[k];

        t_4[k] = f_0 * ih_25[k];

        t_5[k] = f_0 * ih_26[k];

        t_6[k] = f_0 * ih_27[k];

        t_7[k] = f_0 * ih_28[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ih_29, ih_30, ih_31, \
                         ih_32, ih_33, ih_34, ih_35, ih_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ih_29[k];

        t_9[k] = f_0 * ih_30[k];

        t_10[k] = f_0 * ih_31[k];

        t_11[k] = f_0 * ih_32[k];

        t_12[k] = f_0 * ih_33[k];

        t_13[k] = f_0 * ih_34[k];

        t_14[k] = f_0 * ih_35[k];

        t_15[k] = f_0 * ih_36[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, gh_0, gh_1, ih_37, ih_38, \
                         ih_39, ih_40, ih_41, ih_63, ih_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ih_37[k];

        t_17[k] = f_0 * ih_38[k];

        t_18[k] = f_0 * ih_39[k];

        t_19[k] = f_0 * ih_40[k];

        t_20[k] = f_0 * ih_41[k];

        t_21[k] = -gh_0[k]
                  + f_0 * ih_63[k];

        t_22[k] = -gh_1[k]
                  + f_0 * ih_64[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, gh_2, gh_3, gh_4, gh_5, gh_6, ih_65, \
                         ih_66, ih_67, ih_68, ih_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -gh_2[k]
                  + f_0 * ih_65[k];

        t_24[k] = -gh_3[k]
                  + f_0 * ih_66[k];

        t_25[k] = -gh_4[k]
                  + f_0 * ih_67[k];

        t_26[k] = -gh_5[k]
                  + f_0 * ih_68[k];

        t_27[k] = -gh_6[k]
                  + f_0 * ih_69[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, gh_7, gh_8, gh_9, gh_10, gh_11, ih_70, \
                         ih_71, ih_72, ih_73, ih_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -gh_7[k]
                  + f_0 * ih_70[k];

        t_29[k] = -gh_8[k]
                  + f_0 * ih_71[k];

        t_30[k] = -gh_9[k]
                  + f_0 * ih_72[k];

        t_31[k] = -gh_10[k]
                  + f_0 * ih_73[k];

        t_32[k] = -gh_11[k]
                  + f_0 * ih_74[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, gh_12, gh_13, gh_14, gh_15, gh_16, \
                         ih_75, ih_76, ih_77, ih_78, ih_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -gh_12[k]
                  + f_0 * ih_75[k];

        t_34[k] = -gh_13[k]
                  + f_0 * ih_76[k];

        t_35[k] = -gh_14[k]
                  + f_0 * ih_77[k];

        t_36[k] = -gh_15[k]
                  + f_0 * ih_78[k];

        t_37[k] = -gh_16[k]
                  + f_0 * ih_79[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, gh_17, gh_18, gh_19, gh_20, \
                         ih_80, ih_81, ih_82, ih_83, ih_84, ih_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -gh_17[k]
                  + f_0 * ih_80[k];

        t_39[k] = -gh_18[k]
                  + f_0 * ih_81[k];

        t_40[k] = -gh_19[k]
                  + f_0 * ih_82[k];

        t_41[k] = -gh_20[k]
                  + f_0 * ih_83[k];

        t_42[k] = f_0 * ih_84[k];

        t_43[k] = f_0 * ih_85[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, t_51, ih_86, ih_87, ih_88, \
                         ih_89, ih_90, ih_91, ih_92, ih_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * ih_86[k];

        t_45[k] = f_0 * ih_87[k];

        t_46[k] = f_0 * ih_88[k];

        t_47[k] = f_0 * ih_89[k];

        t_48[k] = f_0 * ih_90[k];

        t_49[k] = f_0 * ih_91[k];

        t_50[k] = f_0 * ih_92[k];

        t_51[k] = f_0 * ih_93[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, t_58, t_59, ih_94, ih_95, ih_96, \
                         ih_97, ih_98, ih_99, ih_100, ih_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * ih_94[k];

        t_53[k] = f_0 * ih_95[k];

        t_54[k] = f_0 * ih_96[k];

        t_55[k] = f_0 * ih_97[k];

        t_56[k] = f_0 * ih_98[k];

        t_57[k] = f_0 * ih_99[k];

        t_58[k] = f_0 * ih_100[k];

        t_59[k] = f_0 * ih_101[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, gh_21, gh_22, gh_23, ih_102, \
                         ih_103, ih_104, ih_126, ih_127, ih_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * ih_102[k];

        t_61[k] = f_0 * ih_103[k];

        t_62[k] = f_0 * ih_104[k];

        t_63[k] = -2.0 * gh_21[k]
                  + f_0 * ih_126[k];

        t_64[k] = -2.0 * gh_22[k]
                  + f_0 * ih_127[k];

        t_65[k] = -2.0 * gh_23[k]
                  + f_0 * ih_128[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, gh_24, gh_25, gh_26, gh_27, gh_28, \
                         ih_129, ih_130, ih_131, ih_132, ih_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -2.0 * gh_24[k]
                  + f_0 * ih_129[k];

        t_67[k] = -2.0 * gh_25[k]
                  + f_0 * ih_130[k];

        t_68[k] = -2.0 * gh_26[k]
                  + f_0 * ih_131[k];

        t_69[k] = -2.0 * gh_27[k]
                  + f_0 * ih_132[k];

        t_70[k] = -2.0 * gh_28[k]
                  + f_0 * ih_133[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, gh_29, gh_30, gh_31, gh_32, gh_33, \
                         ih_134, ih_135, ih_136, ih_137, ih_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -2.0 * gh_29[k]
                  + f_0 * ih_134[k];

        t_72[k] = -2.0 * gh_30[k]
                  + f_0 * ih_135[k];

        t_73[k] = -2.0 * gh_31[k]
                  + f_0 * ih_136[k];

        t_74[k] = -2.0 * gh_32[k]
                  + f_0 * ih_137[k];

        t_75[k] = -2.0 * gh_33[k]
                  + f_0 * ih_138[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, gh_34, gh_35, gh_36, gh_37, gh_38, \
                         ih_139, ih_140, ih_141, ih_142, ih_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * gh_34[k]
                  + f_0 * ih_139[k];

        t_77[k] = -2.0 * gh_35[k]
                  + f_0 * ih_140[k];

        t_78[k] = -2.0 * gh_36[k]
                  + f_0 * ih_141[k];

        t_79[k] = -2.0 * gh_37[k]
                  + f_0 * ih_142[k];

        t_80[k] = -2.0 * gh_38[k]
                  + f_0 * ih_143[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, gh_39, gh_40, gh_41, gh_42, gh_43, \
                         ih_144, ih_145, ih_146, ih_147, ih_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -2.0 * gh_39[k]
                  + f_0 * ih_144[k];

        t_82[k] = -2.0 * gh_40[k]
                  + f_0 * ih_145[k];

        t_83[k] = -2.0 * gh_41[k]
                  + f_0 * ih_146[k];

        t_84[k] = -gh_42[k]
                  + f_0 * ih_147[k];

        t_85[k] = -gh_43[k]
                  + f_0 * ih_148[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, gh_44, gh_45, gh_46, gh_47, gh_48, \
                         ih_149, ih_150, ih_151, ih_152, ih_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -gh_44[k]
                  + f_0 * ih_149[k];

        t_87[k] = -gh_45[k]
                  + f_0 * ih_150[k];

        t_88[k] = -gh_46[k]
                  + f_0 * ih_151[k];

        t_89[k] = -gh_47[k]
                  + f_0 * ih_152[k];

        t_90[k] = -gh_48[k]
                  + f_0 * ih_153[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, gh_49, gh_50, gh_51, gh_52, gh_53, \
                         ih_154, ih_155, ih_156, ih_157, ih_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = -gh_49[k]
                  + f_0 * ih_154[k];

        t_92[k] = -gh_50[k]
                  + f_0 * ih_155[k];

        t_93[k] = -gh_51[k]
                  + f_0 * ih_156[k];

        t_94[k] = -gh_52[k]
                  + f_0 * ih_157[k];

        t_95[k] = -gh_53[k]
                  + f_0 * ih_158[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, gh_54, gh_55, gh_56, gh_57, gh_58, \
                         ih_159, ih_160, ih_161, ih_162, ih_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -gh_54[k]
                  + f_0 * ih_159[k];

        t_97[k] = -gh_55[k]
                  + f_0 * ih_160[k];

        t_98[k] = -gh_56[k]
                  + f_0 * ih_161[k];

        t_99[k] = -gh_57[k]
                  + f_0 * ih_162[k];

        t_100[k] = -gh_58[k]
                   + f_0 * ih_163[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, gh_59, gh_60, gh_61, gh_62, \
                         ih_164, ih_165, ih_166, ih_167, ih_168, \
                         ih_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -gh_59[k]
                   + f_0 * ih_164[k];

        t_102[k] = -gh_60[k]
                   + f_0 * ih_165[k];

        t_103[k] = -gh_61[k]
                   + f_0 * ih_166[k];

        t_104[k] = -gh_62[k]
                   + f_0 * ih_167[k];

        t_105[k] = f_0 * ih_168[k];

        t_106[k] = f_0 * ih_169[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, t_114, ih_170, \
                         ih_171, ih_172, ih_173, ih_174, ih_175, ih_176, \
                         ih_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * ih_170[k];

        t_108[k] = f_0 * ih_171[k];

        t_109[k] = f_0 * ih_172[k];

        t_110[k] = f_0 * ih_173[k];

        t_111[k] = f_0 * ih_174[k];

        t_112[k] = f_0 * ih_175[k];

        t_113[k] = f_0 * ih_176[k];

        t_114[k] = f_0 * ih_177[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, t_121, t_122, ih_178, \
                         ih_179, ih_180, ih_181, ih_182, ih_183, ih_184, \
                         ih_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_0 * ih_178[k];

        t_116[k] = f_0 * ih_179[k];

        t_117[k] = f_0 * ih_180[k];

        t_118[k] = f_0 * ih_181[k];

        t_119[k] = f_0 * ih_182[k];

        t_120[k] = f_0 * ih_183[k];

        t_121[k] = f_0 * ih_184[k];

        t_122[k] = f_0 * ih_185[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, t_128, gh_63, gh_64, gh_65, \
                         ih_186, ih_187, ih_188, ih_210, ih_211, \
                         ih_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * ih_186[k];

        t_124[k] = f_0 * ih_187[k];

        t_125[k] = f_0 * ih_188[k];

        t_126[k] = -3.0 * gh_63[k]
                   + f_0 * ih_210[k];

        t_127[k] = -3.0 * gh_64[k]
                   + f_0 * ih_211[k];

        t_128[k] = -3.0 * gh_65[k]
                   + f_0 * ih_212[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, gh_66, gh_67, gh_68, gh_69, gh_70, \
                         ih_213, ih_214, ih_215, ih_216, ih_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -3.0 * gh_66[k]
                   + f_0 * ih_213[k];

        t_130[k] = -3.0 * gh_67[k]
                   + f_0 * ih_214[k];

        t_131[k] = -3.0 * gh_68[k]
                   + f_0 * ih_215[k];

        t_132[k] = -3.0 * gh_69[k]
                   + f_0 * ih_216[k];

        t_133[k] = -3.0 * gh_70[k]
                   + f_0 * ih_217[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, gh_71, gh_72, gh_73, gh_74, gh_75, \
                         ih_218, ih_219, ih_220, ih_221, ih_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -3.0 * gh_71[k]
                   + f_0 * ih_218[k];

        t_135[k] = -3.0 * gh_72[k]
                   + f_0 * ih_219[k];

        t_136[k] = -3.0 * gh_73[k]
                   + f_0 * ih_220[k];

        t_137[k] = -3.0 * gh_74[k]
                   + f_0 * ih_221[k];

        t_138[k] = -3.0 * gh_75[k]
                   + f_0 * ih_222[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, gh_76, gh_77, gh_78, gh_79, gh_80, \
                         ih_223, ih_224, ih_225, ih_226, ih_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -3.0 * gh_76[k]
                   + f_0 * ih_223[k];

        t_140[k] = -3.0 * gh_77[k]
                   + f_0 * ih_224[k];

        t_141[k] = -3.0 * gh_78[k]
                   + f_0 * ih_225[k];

        t_142[k] = -3.0 * gh_79[k]
                   + f_0 * ih_226[k];

        t_143[k] = -3.0 * gh_80[k]
                   + f_0 * ih_227[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, gh_81, gh_82, gh_83, gh_84, gh_85, \
                         ih_228, ih_229, ih_230, ih_231, ih_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -3.0 * gh_81[k]
                   + f_0 * ih_228[k];

        t_145[k] = -3.0 * gh_82[k]
                   + f_0 * ih_229[k];

        t_146[k] = -3.0 * gh_83[k]
                   + f_0 * ih_230[k];

        t_147[k] = -2.0 * gh_84[k]
                   + f_0 * ih_231[k];

        t_148[k] = -2.0 * gh_85[k]
                   + f_0 * ih_232[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, gh_86, gh_87, gh_88, gh_89, gh_90, \
                         ih_233, ih_234, ih_235, ih_236, ih_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * gh_86[k]
                   + f_0 * ih_233[k];

        t_150[k] = -2.0 * gh_87[k]
                   + f_0 * ih_234[k];

        t_151[k] = -2.0 * gh_88[k]
                   + f_0 * ih_235[k];

        t_152[k] = -2.0 * gh_89[k]
                   + f_0 * ih_236[k];

        t_153[k] = -2.0 * gh_90[k]
                   + f_0 * ih_237[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, gh_91, gh_92, gh_93, gh_94, gh_95, \
                         ih_238, ih_239, ih_240, ih_241, ih_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * gh_91[k]
                   + f_0 * ih_238[k];

        t_155[k] = -2.0 * gh_92[k]
                   + f_0 * ih_239[k];

        t_156[k] = -2.0 * gh_93[k]
                   + f_0 * ih_240[k];

        t_157[k] = -2.0 * gh_94[k]
                   + f_0 * ih_241[k];

        t_158[k] = -2.0 * gh_95[k]
                   + f_0 * ih_242[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, gh_96, gh_97, gh_98, gh_99, \
                         gh_100, ih_243, ih_244, ih_245, ih_246, \
                         ih_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * gh_96[k]
                   + f_0 * ih_243[k];

        t_160[k] = -2.0 * gh_97[k]
                   + f_0 * ih_244[k];

        t_161[k] = -2.0 * gh_98[k]
                   + f_0 * ih_245[k];

        t_162[k] = -2.0 * gh_99[k]
                   + f_0 * ih_246[k];

        t_163[k] = -2.0 * gh_100[k]
                   + f_0 * ih_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, gh_101, gh_102, gh_103, gh_104, \
                         gh_105, ih_248, ih_249, ih_250, ih_251, \
                         ih_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * gh_101[k]
                   + f_0 * ih_248[k];

        t_165[k] = -2.0 * gh_102[k]
                   + f_0 * ih_249[k];

        t_166[k] = -2.0 * gh_103[k]
                   + f_0 * ih_250[k];

        t_167[k] = -2.0 * gh_104[k]
                   + f_0 * ih_251[k];

        t_168[k] = -gh_105[k]
                   + f_0 * ih_252[k];
    }
}

static auto
compute_prim_geom_10_hh_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ih,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);

    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_115 = buffer.data(gh + 115);
    const auto *gh_116 = buffer.data(gh + 116);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_118 = buffer.data(gh + 118);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_124 = buffer.data(gh + 124);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_137 = buffer.data(gh + 137);
    const auto *gh_138 = buffer.data(gh + 138);
    const auto *gh_139 = buffer.data(gh + 139);
    const auto *gh_140 = buffer.data(gh + 140);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_142 = buffer.data(gh + 142);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_148 = buffer.data(gh + 148);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_190 = buffer.data(gh + 190);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_193 = buffer.data(gh + 193);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_195 = buffer.data(gh + 195);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_197 = buffer.data(gh + 197);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_199 = buffer.data(gh + 199);
    const auto *gh_200 = buffer.data(gh + 200);
    const auto *gh_201 = buffer.data(gh + 201);
    const auto *gh_202 = buffer.data(gh + 202);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_208 = buffer.data(gh + 208);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_212 = buffer.data(gh + 212);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_221 = buffer.data(gh + 221);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_223 = buffer.data(gh + 223);
    const auto *gh_224 = buffer.data(gh + 224);

    const auto *ih_253 = buffer.data(ih + 253);
    const auto *ih_254 = buffer.data(ih + 254);
    const auto *ih_255 = buffer.data(ih + 255);
    const auto *ih_256 = buffer.data(ih + 256);
    const auto *ih_257 = buffer.data(ih + 257);
    const auto *ih_258 = buffer.data(ih + 258);
    const auto *ih_259 = buffer.data(ih + 259);
    const auto *ih_260 = buffer.data(ih + 260);
    const auto *ih_261 = buffer.data(ih + 261);
    const auto *ih_262 = buffer.data(ih + 262);
    const auto *ih_263 = buffer.data(ih + 263);
    const auto *ih_264 = buffer.data(ih + 264);
    const auto *ih_265 = buffer.data(ih + 265);
    const auto *ih_266 = buffer.data(ih + 266);
    const auto *ih_267 = buffer.data(ih + 267);
    const auto *ih_268 = buffer.data(ih + 268);
    const auto *ih_269 = buffer.data(ih + 269);
    const auto *ih_270 = buffer.data(ih + 270);
    const auto *ih_271 = buffer.data(ih + 271);
    const auto *ih_272 = buffer.data(ih + 272);
    const auto *ih_273 = buffer.data(ih + 273);
    const auto *ih_274 = buffer.data(ih + 274);
    const auto *ih_275 = buffer.data(ih + 275);
    const auto *ih_276 = buffer.data(ih + 276);
    const auto *ih_277 = buffer.data(ih + 277);
    const auto *ih_278 = buffer.data(ih + 278);
    const auto *ih_279 = buffer.data(ih + 279);
    const auto *ih_280 = buffer.data(ih + 280);
    const auto *ih_281 = buffer.data(ih + 281);
    const auto *ih_282 = buffer.data(ih + 282);
    const auto *ih_283 = buffer.data(ih + 283);
    const auto *ih_284 = buffer.data(ih + 284);
    const auto *ih_285 = buffer.data(ih + 285);
    const auto *ih_286 = buffer.data(ih + 286);
    const auto *ih_287 = buffer.data(ih + 287);
    const auto *ih_288 = buffer.data(ih + 288);
    const auto *ih_289 = buffer.data(ih + 289);
    const auto *ih_290 = buffer.data(ih + 290);
    const auto *ih_291 = buffer.data(ih + 291);
    const auto *ih_292 = buffer.data(ih + 292);
    const auto *ih_293 = buffer.data(ih + 293);
    const auto *ih_315 = buffer.data(ih + 315);
    const auto *ih_316 = buffer.data(ih + 316);
    const auto *ih_317 = buffer.data(ih + 317);
    const auto *ih_318 = buffer.data(ih + 318);
    const auto *ih_319 = buffer.data(ih + 319);
    const auto *ih_320 = buffer.data(ih + 320);
    const auto *ih_321 = buffer.data(ih + 321);
    const auto *ih_322 = buffer.data(ih + 322);
    const auto *ih_323 = buffer.data(ih + 323);
    const auto *ih_324 = buffer.data(ih + 324);
    const auto *ih_325 = buffer.data(ih + 325);
    const auto *ih_326 = buffer.data(ih + 326);
    const auto *ih_327 = buffer.data(ih + 327);
    const auto *ih_328 = buffer.data(ih + 328);
    const auto *ih_329 = buffer.data(ih + 329);
    const auto *ih_330 = buffer.data(ih + 330);
    const auto *ih_331 = buffer.data(ih + 331);
    const auto *ih_332 = buffer.data(ih + 332);
    const auto *ih_333 = buffer.data(ih + 333);
    const auto *ih_334 = buffer.data(ih + 334);
    const auto *ih_335 = buffer.data(ih + 335);
    const auto *ih_336 = buffer.data(ih + 336);
    const auto *ih_337 = buffer.data(ih + 337);
    const auto *ih_338 = buffer.data(ih + 338);
    const auto *ih_339 = buffer.data(ih + 339);
    const auto *ih_340 = buffer.data(ih + 340);
    const auto *ih_341 = buffer.data(ih + 341);
    const auto *ih_342 = buffer.data(ih + 342);
    const auto *ih_343 = buffer.data(ih + 343);
    const auto *ih_344 = buffer.data(ih + 344);
    const auto *ih_345 = buffer.data(ih + 345);
    const auto *ih_346 = buffer.data(ih + 346);
    const auto *ih_347 = buffer.data(ih + 347);
    const auto *ih_348 = buffer.data(ih + 348);
    const auto *ih_349 = buffer.data(ih + 349);
    const auto *ih_350 = buffer.data(ih + 350);
    const auto *ih_351 = buffer.data(ih + 351);
    const auto *ih_352 = buffer.data(ih + 352);
    const auto *ih_353 = buffer.data(ih + 353);
    const auto *ih_354 = buffer.data(ih + 354);
    const auto *ih_355 = buffer.data(ih + 355);
    const auto *ih_356 = buffer.data(ih + 356);
    const auto *ih_357 = buffer.data(ih + 357);
    const auto *ih_358 = buffer.data(ih + 358);
    const auto *ih_359 = buffer.data(ih + 359);
    const auto *ih_360 = buffer.data(ih + 360);
    const auto *ih_361 = buffer.data(ih + 361);
    const auto *ih_362 = buffer.data(ih + 362);
    const auto *ih_363 = buffer.data(ih + 363);
    const auto *ih_364 = buffer.data(ih + 364);
    const auto *ih_365 = buffer.data(ih + 365);
    const auto *ih_366 = buffer.data(ih + 366);
    const auto *ih_367 = buffer.data(ih + 367);
    const auto *ih_368 = buffer.data(ih + 368);
    const auto *ih_369 = buffer.data(ih + 369);
    const auto *ih_370 = buffer.data(ih + 370);
    const auto *ih_371 = buffer.data(ih + 371);
    const auto *ih_372 = buffer.data(ih + 372);
    const auto *ih_373 = buffer.data(ih + 373);
    const auto *ih_374 = buffer.data(ih + 374);
    const auto *ih_375 = buffer.data(ih + 375);
    const auto *ih_376 = buffer.data(ih + 376);
    const auto *ih_377 = buffer.data(ih + 377);
    const auto *ih_378 = buffer.data(ih + 378);
    const auto *ih_379 = buffer.data(ih + 379);
    const auto *ih_380 = buffer.data(ih + 380);
    const auto *ih_381 = buffer.data(ih + 381);
    const auto *ih_382 = buffer.data(ih + 382);
    const auto *ih_383 = buffer.data(ih + 383);
    const auto *ih_384 = buffer.data(ih + 384);
    const auto *ih_385 = buffer.data(ih + 385);
    const auto *ih_386 = buffer.data(ih + 386);
    const auto *ih_387 = buffer.data(ih + 387);
    const auto *ih_388 = buffer.data(ih + 388);
    const auto *ih_389 = buffer.data(ih + 389);
    const auto *ih_390 = buffer.data(ih + 390);
    const auto *ih_391 = buffer.data(ih + 391);
    const auto *ih_392 = buffer.data(ih + 392);
    const auto *ih_393 = buffer.data(ih + 393);
    const auto *ih_394 = buffer.data(ih + 394);
    const auto *ih_395 = buffer.data(ih + 395);
    const auto *ih_396 = buffer.data(ih + 396);
    const auto *ih_397 = buffer.data(ih + 397);
    const auto *ih_398 = buffer.data(ih + 398);
    const auto *ih_399 = buffer.data(ih + 399);
    const auto *ih_400 = buffer.data(ih + 400);
    const auto *ih_401 = buffer.data(ih + 401);
    const auto *ih_402 = buffer.data(ih + 402);
    const auto *ih_403 = buffer.data(ih + 403);
    const auto *ih_404 = buffer.data(ih + 404);
    const auto *ih_405 = buffer.data(ih + 405);
    const auto *ih_406 = buffer.data(ih + 406);
    const auto *ih_407 = buffer.data(ih + 407);
    const auto *ih_408 = buffer.data(ih + 408);
    const auto *ih_409 = buffer.data(ih + 409);
    const auto *ih_410 = buffer.data(ih + 410);
    const auto *ih_411 = buffer.data(ih + 411);
    const auto *ih_412 = buffer.data(ih + 412);
    const auto *ih_413 = buffer.data(ih + 413);
    const auto *ih_414 = buffer.data(ih + 414);
    const auto *ih_415 = buffer.data(ih + 415);
    const auto *ih_416 = buffer.data(ih + 416);
    const auto *ih_417 = buffer.data(ih + 417);
    const auto *ih_418 = buffer.data(ih + 418);
    const auto *ih_419 = buffer.data(ih + 419);
    const auto *ih_441 = buffer.data(ih + 441);
    const auto *ih_442 = buffer.data(ih + 442);
    const auto *ih_443 = buffer.data(ih + 443);
    const auto *ih_444 = buffer.data(ih + 444);
    const auto *ih_445 = buffer.data(ih + 445);
    const auto *ih_446 = buffer.data(ih + 446);
    const auto *ih_447 = buffer.data(ih + 447);
    const auto *ih_448 = buffer.data(ih + 448);
    const auto *ih_449 = buffer.data(ih + 449);
    const auto *ih_450 = buffer.data(ih + 450);
    const auto *ih_451 = buffer.data(ih + 451);
    const auto *ih_452 = buffer.data(ih + 452);
    const auto *ih_453 = buffer.data(ih + 453);
    const auto *ih_454 = buffer.data(ih + 454);
    const auto *ih_455 = buffer.data(ih + 455);

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, gh_106, gh_107, gh_108, gh_109, \
                         gh_110, ih_253, ih_254, ih_255, ih_256, \
                         ih_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -gh_106[k]
                   + f_0 * ih_253[k];

        t_170[k] = -gh_107[k]
                   + f_0 * ih_254[k];

        t_171[k] = -gh_108[k]
                   + f_0 * ih_255[k];

        t_172[k] = -gh_109[k]
                   + f_0 * ih_256[k];

        t_173[k] = -gh_110[k]
                   + f_0 * ih_257[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, gh_111, gh_112, gh_113, gh_114, \
                         gh_115, ih_258, ih_259, ih_260, ih_261, \
                         ih_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -gh_111[k]
                   + f_0 * ih_258[k];

        t_175[k] = -gh_112[k]
                   + f_0 * ih_259[k];

        t_176[k] = -gh_113[k]
                   + f_0 * ih_260[k];

        t_177[k] = -gh_114[k]
                   + f_0 * ih_261[k];

        t_178[k] = -gh_115[k]
                   + f_0 * ih_262[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, gh_116, gh_117, gh_118, gh_119, \
                         gh_120, ih_263, ih_264, ih_265, ih_266, \
                         ih_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -gh_116[k]
                   + f_0 * ih_263[k];

        t_180[k] = -gh_117[k]
                   + f_0 * ih_264[k];

        t_181[k] = -gh_118[k]
                   + f_0 * ih_265[k];

        t_182[k] = -gh_119[k]
                   + f_0 * ih_266[k];

        t_183[k] = -gh_120[k]
                   + f_0 * ih_267[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, gh_121, gh_122, gh_123, gh_124, \
                         gh_125, ih_268, ih_269, ih_270, ih_271, \
                         ih_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -gh_121[k]
                   + f_0 * ih_268[k];

        t_185[k] = -gh_122[k]
                   + f_0 * ih_269[k];

        t_186[k] = -gh_123[k]
                   + f_0 * ih_270[k];

        t_187[k] = -gh_124[k]
                   + f_0 * ih_271[k];

        t_188[k] = -gh_125[k]
                   + f_0 * ih_272[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, t_196, ih_273, \
                         ih_274, ih_275, ih_276, ih_277, ih_278, ih_279, \
                         ih_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * ih_273[k];

        t_190[k] = f_0 * ih_274[k];

        t_191[k] = f_0 * ih_275[k];

        t_192[k] = f_0 * ih_276[k];

        t_193[k] = f_0 * ih_277[k];

        t_194[k] = f_0 * ih_278[k];

        t_195[k] = f_0 * ih_279[k];

        t_196[k] = f_0 * ih_280[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, t_203, t_204, ih_281, \
                         ih_282, ih_283, ih_284, ih_285, ih_286, ih_287, \
                         ih_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_0 * ih_281[k];

        t_198[k] = f_0 * ih_282[k];

        t_199[k] = f_0 * ih_283[k];

        t_200[k] = f_0 * ih_284[k];

        t_201[k] = f_0 * ih_285[k];

        t_202[k] = f_0 * ih_286[k];

        t_203[k] = f_0 * ih_287[k];

        t_204[k] = f_0 * ih_288[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, t_210, t_211, gh_126, gh_127, \
                         ih_289, ih_290, ih_291, ih_292, ih_293, ih_315, \
                         ih_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_0 * ih_289[k];

        t_206[k] = f_0 * ih_290[k];

        t_207[k] = f_0 * ih_291[k];

        t_208[k] = f_0 * ih_292[k];

        t_209[k] = f_0 * ih_293[k];

        t_210[k] = -4.0 * gh_126[k]
                   + f_0 * ih_315[k];

        t_211[k] = -4.0 * gh_127[k]
                   + f_0 * ih_316[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, gh_128, gh_129, gh_130, gh_131, \
                         gh_132, ih_317, ih_318, ih_319, ih_320, \
                         ih_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * gh_128[k]
                   + f_0 * ih_317[k];

        t_213[k] = -4.0 * gh_129[k]
                   + f_0 * ih_318[k];

        t_214[k] = -4.0 * gh_130[k]
                   + f_0 * ih_319[k];

        t_215[k] = -4.0 * gh_131[k]
                   + f_0 * ih_320[k];

        t_216[k] = -4.0 * gh_132[k]
                   + f_0 * ih_321[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, gh_133, gh_134, gh_135, gh_136, \
                         gh_137, ih_322, ih_323, ih_324, ih_325, \
                         ih_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * gh_133[k]
                   + f_0 * ih_322[k];

        t_218[k] = -4.0 * gh_134[k]
                   + f_0 * ih_323[k];

        t_219[k] = -4.0 * gh_135[k]
                   + f_0 * ih_324[k];

        t_220[k] = -4.0 * gh_136[k]
                   + f_0 * ih_325[k];

        t_221[k] = -4.0 * gh_137[k]
                   + f_0 * ih_326[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, gh_138, gh_139, gh_140, gh_141, \
                         gh_142, ih_327, ih_328, ih_329, ih_330, \
                         ih_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * gh_138[k]
                   + f_0 * ih_327[k];

        t_223[k] = -4.0 * gh_139[k]
                   + f_0 * ih_328[k];

        t_224[k] = -4.0 * gh_140[k]
                   + f_0 * ih_329[k];

        t_225[k] = -4.0 * gh_141[k]
                   + f_0 * ih_330[k];

        t_226[k] = -4.0 * gh_142[k]
                   + f_0 * ih_331[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, gh_143, gh_144, gh_145, gh_146, \
                         gh_147, ih_332, ih_333, ih_334, ih_335, \
                         ih_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = -4.0 * gh_143[k]
                   + f_0 * ih_332[k];

        t_228[k] = -4.0 * gh_144[k]
                   + f_0 * ih_333[k];

        t_229[k] = -4.0 * gh_145[k]
                   + f_0 * ih_334[k];

        t_230[k] = -4.0 * gh_146[k]
                   + f_0 * ih_335[k];

        t_231[k] = -3.0 * gh_147[k]
                   + f_0 * ih_336[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, gh_148, gh_149, gh_150, gh_151, \
                         gh_152, ih_337, ih_338, ih_339, ih_340, \
                         ih_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = -3.0 * gh_148[k]
                   + f_0 * ih_337[k];

        t_233[k] = -3.0 * gh_149[k]
                   + f_0 * ih_338[k];

        t_234[k] = -3.0 * gh_150[k]
                   + f_0 * ih_339[k];

        t_235[k] = -3.0 * gh_151[k]
                   + f_0 * ih_340[k];

        t_236[k] = -3.0 * gh_152[k]
                   + f_0 * ih_341[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, gh_153, gh_154, gh_155, gh_156, \
                         gh_157, ih_342, ih_343, ih_344, ih_345, \
                         ih_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = -3.0 * gh_153[k]
                   + f_0 * ih_342[k];

        t_238[k] = -3.0 * gh_154[k]
                   + f_0 * ih_343[k];

        t_239[k] = -3.0 * gh_155[k]
                   + f_0 * ih_344[k];

        t_240[k] = -3.0 * gh_156[k]
                   + f_0 * ih_345[k];

        t_241[k] = -3.0 * gh_157[k]
                   + f_0 * ih_346[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, gh_158, gh_159, gh_160, gh_161, \
                         gh_162, ih_347, ih_348, ih_349, ih_350, \
                         ih_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = -3.0 * gh_158[k]
                   + f_0 * ih_347[k];

        t_243[k] = -3.0 * gh_159[k]
                   + f_0 * ih_348[k];

        t_244[k] = -3.0 * gh_160[k]
                   + f_0 * ih_349[k];

        t_245[k] = -3.0 * gh_161[k]
                   + f_0 * ih_350[k];

        t_246[k] = -3.0 * gh_162[k]
                   + f_0 * ih_351[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, gh_163, gh_164, gh_165, gh_166, \
                         gh_167, ih_352, ih_353, ih_354, ih_355, \
                         ih_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = -3.0 * gh_163[k]
                   + f_0 * ih_352[k];

        t_248[k] = -3.0 * gh_164[k]
                   + f_0 * ih_353[k];

        t_249[k] = -3.0 * gh_165[k]
                   + f_0 * ih_354[k];

        t_250[k] = -3.0 * gh_166[k]
                   + f_0 * ih_355[k];

        t_251[k] = -3.0 * gh_167[k]
                   + f_0 * ih_356[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, gh_168, gh_169, gh_170, gh_171, \
                         gh_172, ih_357, ih_358, ih_359, ih_360, \
                         ih_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -2.0 * gh_168[k]
                   + f_0 * ih_357[k];

        t_253[k] = -2.0 * gh_169[k]
                   + f_0 * ih_358[k];

        t_254[k] = -2.0 * gh_170[k]
                   + f_0 * ih_359[k];

        t_255[k] = -2.0 * gh_171[k]
                   + f_0 * ih_360[k];

        t_256[k] = -2.0 * gh_172[k]
                   + f_0 * ih_361[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, gh_173, gh_174, gh_175, gh_176, \
                         gh_177, ih_362, ih_363, ih_364, ih_365, \
                         ih_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -2.0 * gh_173[k]
                   + f_0 * ih_362[k];

        t_258[k] = -2.0 * gh_174[k]
                   + f_0 * ih_363[k];

        t_259[k] = -2.0 * gh_175[k]
                   + f_0 * ih_364[k];

        t_260[k] = -2.0 * gh_176[k]
                   + f_0 * ih_365[k];

        t_261[k] = -2.0 * gh_177[k]
                   + f_0 * ih_366[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, gh_178, gh_179, gh_180, gh_181, \
                         gh_182, ih_367, ih_368, ih_369, ih_370, \
                         ih_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = -2.0 * gh_178[k]
                   + f_0 * ih_367[k];

        t_263[k] = -2.0 * gh_179[k]
                   + f_0 * ih_368[k];

        t_264[k] = -2.0 * gh_180[k]
                   + f_0 * ih_369[k];

        t_265[k] = -2.0 * gh_181[k]
                   + f_0 * ih_370[k];

        t_266[k] = -2.0 * gh_182[k]
                   + f_0 * ih_371[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, gh_183, gh_184, gh_185, gh_186, \
                         gh_187, ih_372, ih_373, ih_374, ih_375, \
                         ih_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = -2.0 * gh_183[k]
                   + f_0 * ih_372[k];

        t_268[k] = -2.0 * gh_184[k]
                   + f_0 * ih_373[k];

        t_269[k] = -2.0 * gh_185[k]
                   + f_0 * ih_374[k];

        t_270[k] = -2.0 * gh_186[k]
                   + f_0 * ih_375[k];

        t_271[k] = -2.0 * gh_187[k]
                   + f_0 * ih_376[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, gh_188, gh_189, gh_190, gh_191, \
                         gh_192, ih_377, ih_378, ih_379, ih_380, \
                         ih_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -2.0 * gh_188[k]
                   + f_0 * ih_377[k];

        t_273[k] = -gh_189[k]
                   + f_0 * ih_378[k];

        t_274[k] = -gh_190[k]
                   + f_0 * ih_379[k];

        t_275[k] = -gh_191[k]
                   + f_0 * ih_380[k];

        t_276[k] = -gh_192[k]
                   + f_0 * ih_381[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, gh_193, gh_194, gh_195, gh_196, \
                         gh_197, ih_382, ih_383, ih_384, ih_385, \
                         ih_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -gh_193[k]
                   + f_0 * ih_382[k];

        t_278[k] = -gh_194[k]
                   + f_0 * ih_383[k];

        t_279[k] = -gh_195[k]
                   + f_0 * ih_384[k];

        t_280[k] = -gh_196[k]
                   + f_0 * ih_385[k];

        t_281[k] = -gh_197[k]
                   + f_0 * ih_386[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, gh_198, gh_199, gh_200, gh_201, \
                         gh_202, ih_387, ih_388, ih_389, ih_390, \
                         ih_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = -gh_198[k]
                   + f_0 * ih_387[k];

        t_283[k] = -gh_199[k]
                   + f_0 * ih_388[k];

        t_284[k] = -gh_200[k]
                   + f_0 * ih_389[k];

        t_285[k] = -gh_201[k]
                   + f_0 * ih_390[k];

        t_286[k] = -gh_202[k]
                   + f_0 * ih_391[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, gh_203, gh_204, gh_205, gh_206, \
                         gh_207, ih_392, ih_393, ih_394, ih_395, \
                         ih_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = -gh_203[k]
                   + f_0 * ih_392[k];

        t_288[k] = -gh_204[k]
                   + f_0 * ih_393[k];

        t_289[k] = -gh_205[k]
                   + f_0 * ih_394[k];

        t_290[k] = -gh_206[k]
                   + f_0 * ih_395[k];

        t_291[k] = -gh_207[k]
                   + f_0 * ih_396[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, t_297, t_298, gh_208, gh_209, \
                         ih_397, ih_398, ih_399, ih_400, ih_401, ih_402, \
                         ih_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -gh_208[k]
                   + f_0 * ih_397[k];

        t_293[k] = -gh_209[k]
                   + f_0 * ih_398[k];

        t_294[k] = f_0 * ih_399[k];

        t_295[k] = f_0 * ih_400[k];

        t_296[k] = f_0 * ih_401[k];

        t_297[k] = f_0 * ih_402[k];

        t_298[k] = f_0 * ih_403[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, t_304, t_305, t_306, ih_404, \
                         ih_405, ih_406, ih_407, ih_408, ih_409, ih_410, \
                         ih_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_0 * ih_404[k];

        t_300[k] = f_0 * ih_405[k];

        t_301[k] = f_0 * ih_406[k];

        t_302[k] = f_0 * ih_407[k];

        t_303[k] = f_0 * ih_408[k];

        t_304[k] = f_0 * ih_409[k];

        t_305[k] = f_0 * ih_410[k];

        t_306[k] = f_0 * ih_411[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, t_312, t_313, t_314, ih_412, \
                         ih_413, ih_414, ih_415, ih_416, ih_417, ih_418, \
                         ih_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_0 * ih_412[k];

        t_308[k] = f_0 * ih_413[k];

        t_309[k] = f_0 * ih_414[k];

        t_310[k] = f_0 * ih_415[k];

        t_311[k] = f_0 * ih_416[k];

        t_312[k] = f_0 * ih_417[k];

        t_313[k] = f_0 * ih_418[k];

        t_314[k] = f_0 * ih_419[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, gh_210, gh_211, gh_212, gh_213, \
                         gh_214, ih_441, ih_442, ih_443, ih_444, \
                         ih_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -5.0 * gh_210[k]
                   + f_0 * ih_441[k];

        t_316[k] = -5.0 * gh_211[k]
                   + f_0 * ih_442[k];

        t_317[k] = -5.0 * gh_212[k]
                   + f_0 * ih_443[k];

        t_318[k] = -5.0 * gh_213[k]
                   + f_0 * ih_444[k];

        t_319[k] = -5.0 * gh_214[k]
                   + f_0 * ih_445[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, gh_215, gh_216, gh_217, gh_218, \
                         gh_219, ih_446, ih_447, ih_448, ih_449, \
                         ih_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -5.0 * gh_215[k]
                   + f_0 * ih_446[k];

        t_321[k] = -5.0 * gh_216[k]
                   + f_0 * ih_447[k];

        t_322[k] = -5.0 * gh_217[k]
                   + f_0 * ih_448[k];

        t_323[k] = -5.0 * gh_218[k]
                   + f_0 * ih_449[k];

        t_324[k] = -5.0 * gh_219[k]
                   + f_0 * ih_450[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, gh_220, gh_221, gh_222, gh_223, \
                         gh_224, ih_451, ih_452, ih_453, ih_454, \
                         ih_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -5.0 * gh_220[k]
                   + f_0 * ih_451[k];

        t_326[k] = -5.0 * gh_221[k]
                   + f_0 * ih_452[k];

        t_327[k] = -5.0 * gh_222[k]
                   + f_0 * ih_453[k];

        t_328[k] = -5.0 * gh_223[k]
                   + f_0 * ih_454[k];

        t_329[k] = -5.0 * gh_224[k]
                   + f_0 * ih_455[k];
    }
}

static auto
compute_prim_geom_10_hh_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ih,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_231 = buffer.data(gh + 231);
    const auto *gh_232 = buffer.data(gh + 232);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_234 = buffer.data(gh + 234);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_237 = buffer.data(gh + 237);
    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_241 = buffer.data(gh + 241);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_273 = buffer.data(gh + 273);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_275 = buffer.data(gh + 275);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_278 = buffer.data(gh + 278);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);
    const auto *gh_282 = buffer.data(gh + 282);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_287 = buffer.data(gh + 287);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_295 = buffer.data(gh + 295);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_298 = buffer.data(gh + 298);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_302 = buffer.data(gh + 302);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_304 = buffer.data(gh + 304);
    const auto *gh_305 = buffer.data(gh + 305);
    const auto *gh_306 = buffer.data(gh + 306);
    const auto *gh_307 = buffer.data(gh + 307);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

    const auto *ih_456 = buffer.data(ih + 456);
    const auto *ih_457 = buffer.data(ih + 457);
    const auto *ih_458 = buffer.data(ih + 458);
    const auto *ih_459 = buffer.data(ih + 459);
    const auto *ih_460 = buffer.data(ih + 460);
    const auto *ih_461 = buffer.data(ih + 461);
    const auto *ih_462 = buffer.data(ih + 462);
    const auto *ih_463 = buffer.data(ih + 463);
    const auto *ih_464 = buffer.data(ih + 464);
    const auto *ih_465 = buffer.data(ih + 465);
    const auto *ih_466 = buffer.data(ih + 466);
    const auto *ih_467 = buffer.data(ih + 467);
    const auto *ih_468 = buffer.data(ih + 468);
    const auto *ih_469 = buffer.data(ih + 469);
    const auto *ih_470 = buffer.data(ih + 470);
    const auto *ih_471 = buffer.data(ih + 471);
    const auto *ih_472 = buffer.data(ih + 472);
    const auto *ih_473 = buffer.data(ih + 473);
    const auto *ih_474 = buffer.data(ih + 474);
    const auto *ih_475 = buffer.data(ih + 475);
    const auto *ih_476 = buffer.data(ih + 476);
    const auto *ih_477 = buffer.data(ih + 477);
    const auto *ih_478 = buffer.data(ih + 478);
    const auto *ih_479 = buffer.data(ih + 479);
    const auto *ih_480 = buffer.data(ih + 480);
    const auto *ih_481 = buffer.data(ih + 481);
    const auto *ih_482 = buffer.data(ih + 482);
    const auto *ih_483 = buffer.data(ih + 483);
    const auto *ih_484 = buffer.data(ih + 484);
    const auto *ih_485 = buffer.data(ih + 485);
    const auto *ih_486 = buffer.data(ih + 486);
    const auto *ih_487 = buffer.data(ih + 487);
    const auto *ih_488 = buffer.data(ih + 488);
    const auto *ih_489 = buffer.data(ih + 489);
    const auto *ih_490 = buffer.data(ih + 490);
    const auto *ih_491 = buffer.data(ih + 491);
    const auto *ih_492 = buffer.data(ih + 492);
    const auto *ih_493 = buffer.data(ih + 493);
    const auto *ih_494 = buffer.data(ih + 494);
    const auto *ih_495 = buffer.data(ih + 495);
    const auto *ih_496 = buffer.data(ih + 496);
    const auto *ih_497 = buffer.data(ih + 497);
    const auto *ih_498 = buffer.data(ih + 498);
    const auto *ih_499 = buffer.data(ih + 499);
    const auto *ih_500 = buffer.data(ih + 500);
    const auto *ih_501 = buffer.data(ih + 501);
    const auto *ih_502 = buffer.data(ih + 502);
    const auto *ih_503 = buffer.data(ih + 503);
    const auto *ih_504 = buffer.data(ih + 504);
    const auto *ih_505 = buffer.data(ih + 505);
    const auto *ih_506 = buffer.data(ih + 506);
    const auto *ih_507 = buffer.data(ih + 507);
    const auto *ih_508 = buffer.data(ih + 508);
    const auto *ih_509 = buffer.data(ih + 509);
    const auto *ih_510 = buffer.data(ih + 510);
    const auto *ih_511 = buffer.data(ih + 511);
    const auto *ih_512 = buffer.data(ih + 512);
    const auto *ih_513 = buffer.data(ih + 513);
    const auto *ih_514 = buffer.data(ih + 514);
    const auto *ih_515 = buffer.data(ih + 515);
    const auto *ih_516 = buffer.data(ih + 516);
    const auto *ih_517 = buffer.data(ih + 517);
    const auto *ih_518 = buffer.data(ih + 518);
    const auto *ih_519 = buffer.data(ih + 519);
    const auto *ih_520 = buffer.data(ih + 520);
    const auto *ih_521 = buffer.data(ih + 521);
    const auto *ih_522 = buffer.data(ih + 522);
    const auto *ih_523 = buffer.data(ih + 523);
    const auto *ih_524 = buffer.data(ih + 524);
    const auto *ih_525 = buffer.data(ih + 525);
    const auto *ih_526 = buffer.data(ih + 526);
    const auto *ih_527 = buffer.data(ih + 527);
    const auto *ih_528 = buffer.data(ih + 528);
    const auto *ih_529 = buffer.data(ih + 529);
    const auto *ih_530 = buffer.data(ih + 530);
    const auto *ih_531 = buffer.data(ih + 531);
    const auto *ih_532 = buffer.data(ih + 532);
    const auto *ih_533 = buffer.data(ih + 533);
    const auto *ih_534 = buffer.data(ih + 534);
    const auto *ih_535 = buffer.data(ih + 535);
    const auto *ih_536 = buffer.data(ih + 536);
    const auto *ih_537 = buffer.data(ih + 537);
    const auto *ih_538 = buffer.data(ih + 538);
    const auto *ih_539 = buffer.data(ih + 539);
    const auto *ih_540 = buffer.data(ih + 540);
    const auto *ih_541 = buffer.data(ih + 541);
    const auto *ih_542 = buffer.data(ih + 542);
    const auto *ih_543 = buffer.data(ih + 543);
    const auto *ih_544 = buffer.data(ih + 544);
    const auto *ih_545 = buffer.data(ih + 545);
    const auto *ih_546 = buffer.data(ih + 546);
    const auto *ih_547 = buffer.data(ih + 547);
    const auto *ih_548 = buffer.data(ih + 548);
    const auto *ih_549 = buffer.data(ih + 549);
    const auto *ih_550 = buffer.data(ih + 550);
    const auto *ih_551 = buffer.data(ih + 551);
    const auto *ih_552 = buffer.data(ih + 552);
    const auto *ih_553 = buffer.data(ih + 553);
    const auto *ih_554 = buffer.data(ih + 554);
    const auto *ih_555 = buffer.data(ih + 555);
    const auto *ih_556 = buffer.data(ih + 556);
    const auto *ih_557 = buffer.data(ih + 557);
    const auto *ih_558 = buffer.data(ih + 558);
    const auto *ih_559 = buffer.data(ih + 559);
    const auto *ih_560 = buffer.data(ih + 560);
    const auto *ih_561 = buffer.data(ih + 561);
    const auto *ih_562 = buffer.data(ih + 562);
    const auto *ih_563 = buffer.data(ih + 563);
    const auto *ih_564 = buffer.data(ih + 564);
    const auto *ih_565 = buffer.data(ih + 565);
    const auto *ih_566 = buffer.data(ih + 566);

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, gh_225, gh_226, gh_227, gh_228, \
                         gh_229, ih_456, ih_457, ih_458, ih_459, \
                         ih_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -5.0 * gh_225[k]
                   + f_0 * ih_456[k];

        t_331[k] = -5.0 * gh_226[k]
                   + f_0 * ih_457[k];

        t_332[k] = -5.0 * gh_227[k]
                   + f_0 * ih_458[k];

        t_333[k] = -5.0 * gh_228[k]
                   + f_0 * ih_459[k];

        t_334[k] = -5.0 * gh_229[k]
                   + f_0 * ih_460[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, gh_230, gh_231, gh_232, gh_233, \
                         gh_234, ih_461, ih_462, ih_463, ih_464, \
                         ih_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -5.0 * gh_230[k]
                   + f_0 * ih_461[k];

        t_336[k] = -4.0 * gh_231[k]
                   + f_0 * ih_462[k];

        t_337[k] = -4.0 * gh_232[k]
                   + f_0 * ih_463[k];

        t_338[k] = -4.0 * gh_233[k]
                   + f_0 * ih_464[k];

        t_339[k] = -4.0 * gh_234[k]
                   + f_0 * ih_465[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, gh_235, gh_236, gh_237, gh_238, \
                         gh_239, ih_466, ih_467, ih_468, ih_469, \
                         ih_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -4.0 * gh_235[k]
                   + f_0 * ih_466[k];

        t_341[k] = -4.0 * gh_236[k]
                   + f_0 * ih_467[k];

        t_342[k] = -4.0 * gh_237[k]
                   + f_0 * ih_468[k];

        t_343[k] = -4.0 * gh_238[k]
                   + f_0 * ih_469[k];

        t_344[k] = -4.0 * gh_239[k]
                   + f_0 * ih_470[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, gh_240, gh_241, gh_242, gh_243, \
                         gh_244, ih_471, ih_472, ih_473, ih_474, \
                         ih_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -4.0 * gh_240[k]
                   + f_0 * ih_471[k];

        t_346[k] = -4.0 * gh_241[k]
                   + f_0 * ih_472[k];

        t_347[k] = -4.0 * gh_242[k]
                   + f_0 * ih_473[k];

        t_348[k] = -4.0 * gh_243[k]
                   + f_0 * ih_474[k];

        t_349[k] = -4.0 * gh_244[k]
                   + f_0 * ih_475[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, gh_245, gh_246, gh_247, gh_248, \
                         gh_249, ih_476, ih_477, ih_478, ih_479, \
                         ih_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -4.0 * gh_245[k]
                   + f_0 * ih_476[k];

        t_351[k] = -4.0 * gh_246[k]
                   + f_0 * ih_477[k];

        t_352[k] = -4.0 * gh_247[k]
                   + f_0 * ih_478[k];

        t_353[k] = -4.0 * gh_248[k]
                   + f_0 * ih_479[k];

        t_354[k] = -4.0 * gh_249[k]
                   + f_0 * ih_480[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, gh_250, gh_251, gh_252, gh_253, \
                         gh_254, ih_481, ih_482, ih_483, ih_484, \
                         ih_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -4.0 * gh_250[k]
                   + f_0 * ih_481[k];

        t_356[k] = -4.0 * gh_251[k]
                   + f_0 * ih_482[k];

        t_357[k] = -3.0 * gh_252[k]
                   + f_0 * ih_483[k];

        t_358[k] = -3.0 * gh_253[k]
                   + f_0 * ih_484[k];

        t_359[k] = -3.0 * gh_254[k]
                   + f_0 * ih_485[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, gh_255, gh_256, gh_257, gh_258, \
                         gh_259, ih_486, ih_487, ih_488, ih_489, \
                         ih_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -3.0 * gh_255[k]
                   + f_0 * ih_486[k];

        t_361[k] = -3.0 * gh_256[k]
                   + f_0 * ih_487[k];

        t_362[k] = -3.0 * gh_257[k]
                   + f_0 * ih_488[k];

        t_363[k] = -3.0 * gh_258[k]
                   + f_0 * ih_489[k];

        t_364[k] = -3.0 * gh_259[k]
                   + f_0 * ih_490[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, gh_260, gh_261, gh_262, gh_263, \
                         gh_264, ih_491, ih_492, ih_493, ih_494, \
                         ih_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -3.0 * gh_260[k]
                   + f_0 * ih_491[k];

        t_366[k] = -3.0 * gh_261[k]
                   + f_0 * ih_492[k];

        t_367[k] = -3.0 * gh_262[k]
                   + f_0 * ih_493[k];

        t_368[k] = -3.0 * gh_263[k]
                   + f_0 * ih_494[k];

        t_369[k] = -3.0 * gh_264[k]
                   + f_0 * ih_495[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, gh_265, gh_266, gh_267, gh_268, \
                         gh_269, ih_496, ih_497, ih_498, ih_499, \
                         ih_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -3.0 * gh_265[k]
                   + f_0 * ih_496[k];

        t_371[k] = -3.0 * gh_266[k]
                   + f_0 * ih_497[k];

        t_372[k] = -3.0 * gh_267[k]
                   + f_0 * ih_498[k];

        t_373[k] = -3.0 * gh_268[k]
                   + f_0 * ih_499[k];

        t_374[k] = -3.0 * gh_269[k]
                   + f_0 * ih_500[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, gh_270, gh_271, gh_272, gh_273, \
                         gh_274, ih_501, ih_502, ih_503, ih_504, \
                         ih_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -3.0 * gh_270[k]
                   + f_0 * ih_501[k];

        t_376[k] = -3.0 * gh_271[k]
                   + f_0 * ih_502[k];

        t_377[k] = -3.0 * gh_272[k]
                   + f_0 * ih_503[k];

        t_378[k] = -2.0 * gh_273[k]
                   + f_0 * ih_504[k];

        t_379[k] = -2.0 * gh_274[k]
                   + f_0 * ih_505[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, gh_275, gh_276, gh_277, gh_278, \
                         gh_279, ih_506, ih_507, ih_508, ih_509, \
                         ih_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -2.0 * gh_275[k]
                   + f_0 * ih_506[k];

        t_381[k] = -2.0 * gh_276[k]
                   + f_0 * ih_507[k];

        t_382[k] = -2.0 * gh_277[k]
                   + f_0 * ih_508[k];

        t_383[k] = -2.0 * gh_278[k]
                   + f_0 * ih_509[k];

        t_384[k] = -2.0 * gh_279[k]
                   + f_0 * ih_510[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, gh_280, gh_281, gh_282, gh_283, \
                         gh_284, ih_511, ih_512, ih_513, ih_514, \
                         ih_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -2.0 * gh_280[k]
                   + f_0 * ih_511[k];

        t_386[k] = -2.0 * gh_281[k]
                   + f_0 * ih_512[k];

        t_387[k] = -2.0 * gh_282[k]
                   + f_0 * ih_513[k];

        t_388[k] = -2.0 * gh_283[k]
                   + f_0 * ih_514[k];

        t_389[k] = -2.0 * gh_284[k]
                   + f_0 * ih_515[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, gh_285, gh_286, gh_287, gh_288, \
                         gh_289, ih_516, ih_517, ih_518, ih_519, \
                         ih_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -2.0 * gh_285[k]
                   + f_0 * ih_516[k];

        t_391[k] = -2.0 * gh_286[k]
                   + f_0 * ih_517[k];

        t_392[k] = -2.0 * gh_287[k]
                   + f_0 * ih_518[k];

        t_393[k] = -2.0 * gh_288[k]
                   + f_0 * ih_519[k];

        t_394[k] = -2.0 * gh_289[k]
                   + f_0 * ih_520[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, gh_290, gh_291, gh_292, gh_293, \
                         gh_294, ih_521, ih_522, ih_523, ih_524, \
                         ih_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -2.0 * gh_290[k]
                   + f_0 * ih_521[k];

        t_396[k] = -2.0 * gh_291[k]
                   + f_0 * ih_522[k];

        t_397[k] = -2.0 * gh_292[k]
                   + f_0 * ih_523[k];

        t_398[k] = -2.0 * gh_293[k]
                   + f_0 * ih_524[k];

        t_399[k] = -gh_294[k]
                   + f_0 * ih_525[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, gh_295, gh_296, gh_297, gh_298, \
                         gh_299, ih_526, ih_527, ih_528, ih_529, \
                         ih_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -gh_295[k]
                   + f_0 * ih_526[k];

        t_401[k] = -gh_296[k]
                   + f_0 * ih_527[k];

        t_402[k] = -gh_297[k]
                   + f_0 * ih_528[k];

        t_403[k] = -gh_298[k]
                   + f_0 * ih_529[k];

        t_404[k] = -gh_299[k]
                   + f_0 * ih_530[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, gh_300, gh_301, gh_302, gh_303, \
                         gh_304, ih_531, ih_532, ih_533, ih_534, \
                         ih_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -gh_300[k]
                   + f_0 * ih_531[k];

        t_406[k] = -gh_301[k]
                   + f_0 * ih_532[k];

        t_407[k] = -gh_302[k]
                   + f_0 * ih_533[k];

        t_408[k] = -gh_303[k]
                   + f_0 * ih_534[k];

        t_409[k] = -gh_304[k]
                   + f_0 * ih_535[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, gh_305, gh_306, gh_307, gh_308, \
                         gh_309, ih_536, ih_537, ih_538, ih_539, \
                         ih_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -gh_305[k]
                   + f_0 * ih_536[k];

        t_411[k] = -gh_306[k]
                   + f_0 * ih_537[k];

        t_412[k] = -gh_307[k]
                   + f_0 * ih_538[k];

        t_413[k] = -gh_308[k]
                   + f_0 * ih_539[k];

        t_414[k] = -gh_309[k]
                   + f_0 * ih_540[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, gh_310, gh_311, gh_312, gh_313, \
                         gh_314, ih_541, ih_542, ih_543, ih_544, \
                         ih_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -gh_310[k]
                   + f_0 * ih_541[k];

        t_416[k] = -gh_311[k]
                   + f_0 * ih_542[k];

        t_417[k] = -gh_312[k]
                   + f_0 * ih_543[k];

        t_418[k] = -gh_313[k]
                   + f_0 * ih_544[k];

        t_419[k] = -gh_314[k]
                   + f_0 * ih_545[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, t_425, t_426, t_427, ih_546, \
                         ih_547, ih_548, ih_549, ih_550, ih_551, ih_552, \
                         ih_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_0 * ih_546[k];

        t_421[k] = f_0 * ih_547[k];

        t_422[k] = f_0 * ih_548[k];

        t_423[k] = f_0 * ih_549[k];

        t_424[k] = f_0 * ih_550[k];

        t_425[k] = f_0 * ih_551[k];

        t_426[k] = f_0 * ih_552[k];

        t_427[k] = f_0 * ih_553[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, t_433, t_434, t_435, ih_554, \
                         ih_555, ih_556, ih_557, ih_558, ih_559, ih_560, \
                         ih_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_0 * ih_554[k];

        t_429[k] = f_0 * ih_555[k];

        t_430[k] = f_0 * ih_556[k];

        t_431[k] = f_0 * ih_557[k];

        t_432[k] = f_0 * ih_558[k];

        t_433[k] = f_0 * ih_559[k];

        t_434[k] = f_0 * ih_560[k];

        t_435[k] = f_0 * ih_561[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, ih_562, ih_563, ih_564, ih_565, \
                         ih_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_0 * ih_562[k];

        t_437[k] = f_0 * ih_563[k];

        t_438[k] = f_0 * ih_564[k];

        t_439[k] = f_0 * ih_565[k];

        t_440[k] = f_0 * ih_566[k];
    }
}

auto
compute_prim_geom_10_hh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t gh, const size_t ih,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hh_electron_repulsion_1_piece0(buffer, target, gh, ih, ncols, alpha);

    compute_prim_geom_10_hh_electron_repulsion_1_piece1(buffer, target, gh, ih, ncols, alpha);

    compute_prim_geom_10_hh_electron_repulsion_1_piece2(buffer, target, gh, ih, ncols, alpha);
}

static auto
compute_prim_geom_10_hh_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ih,
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

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_71 = buffer.data(gh + 71);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_74 = buffer.data(gh + 74);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_76 = buffer.data(gh + 76);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);

    const auto *ih_42 = buffer.data(ih + 42);
    const auto *ih_43 = buffer.data(ih + 43);
    const auto *ih_44 = buffer.data(ih + 44);
    const auto *ih_45 = buffer.data(ih + 45);
    const auto *ih_46 = buffer.data(ih + 46);
    const auto *ih_47 = buffer.data(ih + 47);
    const auto *ih_48 = buffer.data(ih + 48);
    const auto *ih_49 = buffer.data(ih + 49);
    const auto *ih_50 = buffer.data(ih + 50);
    const auto *ih_51 = buffer.data(ih + 51);
    const auto *ih_52 = buffer.data(ih + 52);
    const auto *ih_53 = buffer.data(ih + 53);
    const auto *ih_54 = buffer.data(ih + 54);
    const auto *ih_55 = buffer.data(ih + 55);
    const auto *ih_56 = buffer.data(ih + 56);
    const auto *ih_57 = buffer.data(ih + 57);
    const auto *ih_58 = buffer.data(ih + 58);
    const auto *ih_59 = buffer.data(ih + 59);
    const auto *ih_60 = buffer.data(ih + 60);
    const auto *ih_61 = buffer.data(ih + 61);
    const auto *ih_62 = buffer.data(ih + 62);
    const auto *ih_84 = buffer.data(ih + 84);
    const auto *ih_85 = buffer.data(ih + 85);
    const auto *ih_86 = buffer.data(ih + 86);
    const auto *ih_87 = buffer.data(ih + 87);
    const auto *ih_88 = buffer.data(ih + 88);
    const auto *ih_89 = buffer.data(ih + 89);
    const auto *ih_90 = buffer.data(ih + 90);
    const auto *ih_91 = buffer.data(ih + 91);
    const auto *ih_92 = buffer.data(ih + 92);
    const auto *ih_93 = buffer.data(ih + 93);
    const auto *ih_94 = buffer.data(ih + 94);
    const auto *ih_95 = buffer.data(ih + 95);
    const auto *ih_96 = buffer.data(ih + 96);
    const auto *ih_97 = buffer.data(ih + 97);
    const auto *ih_98 = buffer.data(ih + 98);
    const auto *ih_99 = buffer.data(ih + 99);
    const auto *ih_100 = buffer.data(ih + 100);
    const auto *ih_101 = buffer.data(ih + 101);
    const auto *ih_102 = buffer.data(ih + 102);
    const auto *ih_103 = buffer.data(ih + 103);
    const auto *ih_104 = buffer.data(ih + 104);
    const auto *ih_105 = buffer.data(ih + 105);
    const auto *ih_106 = buffer.data(ih + 106);
    const auto *ih_107 = buffer.data(ih + 107);
    const auto *ih_108 = buffer.data(ih + 108);
    const auto *ih_109 = buffer.data(ih + 109);
    const auto *ih_110 = buffer.data(ih + 110);
    const auto *ih_111 = buffer.data(ih + 111);
    const auto *ih_112 = buffer.data(ih + 112);
    const auto *ih_113 = buffer.data(ih + 113);
    const auto *ih_114 = buffer.data(ih + 114);
    const auto *ih_115 = buffer.data(ih + 115);
    const auto *ih_116 = buffer.data(ih + 116);
    const auto *ih_117 = buffer.data(ih + 117);
    const auto *ih_118 = buffer.data(ih + 118);
    const auto *ih_119 = buffer.data(ih + 119);
    const auto *ih_120 = buffer.data(ih + 120);
    const auto *ih_121 = buffer.data(ih + 121);
    const auto *ih_122 = buffer.data(ih + 122);
    const auto *ih_123 = buffer.data(ih + 123);
    const auto *ih_124 = buffer.data(ih + 124);
    const auto *ih_125 = buffer.data(ih + 125);
    const auto *ih_147 = buffer.data(ih + 147);
    const auto *ih_148 = buffer.data(ih + 148);
    const auto *ih_149 = buffer.data(ih + 149);
    const auto *ih_150 = buffer.data(ih + 150);
    const auto *ih_151 = buffer.data(ih + 151);
    const auto *ih_152 = buffer.data(ih + 152);
    const auto *ih_153 = buffer.data(ih + 153);
    const auto *ih_154 = buffer.data(ih + 154);
    const auto *ih_155 = buffer.data(ih + 155);
    const auto *ih_156 = buffer.data(ih + 156);
    const auto *ih_157 = buffer.data(ih + 157);
    const auto *ih_158 = buffer.data(ih + 158);
    const auto *ih_159 = buffer.data(ih + 159);
    const auto *ih_160 = buffer.data(ih + 160);
    const auto *ih_161 = buffer.data(ih + 161);
    const auto *ih_162 = buffer.data(ih + 162);
    const auto *ih_163 = buffer.data(ih + 163);
    const auto *ih_164 = buffer.data(ih + 164);
    const auto *ih_165 = buffer.data(ih + 165);
    const auto *ih_166 = buffer.data(ih + 166);
    const auto *ih_167 = buffer.data(ih + 167);
    const auto *ih_168 = buffer.data(ih + 168);
    const auto *ih_169 = buffer.data(ih + 169);
    const auto *ih_170 = buffer.data(ih + 170);
    const auto *ih_171 = buffer.data(ih + 171);
    const auto *ih_172 = buffer.data(ih + 172);
    const auto *ih_173 = buffer.data(ih + 173);
    const auto *ih_174 = buffer.data(ih + 174);
    const auto *ih_175 = buffer.data(ih + 175);
    const auto *ih_176 = buffer.data(ih + 176);
    const auto *ih_177 = buffer.data(ih + 177);
    const auto *ih_178 = buffer.data(ih + 178);
    const auto *ih_179 = buffer.data(ih + 179);
    const auto *ih_180 = buffer.data(ih + 180);
    const auto *ih_181 = buffer.data(ih + 181);
    const auto *ih_182 = buffer.data(ih + 182);
    const auto *ih_183 = buffer.data(ih + 183);
    const auto *ih_184 = buffer.data(ih + 184);
    const auto *ih_185 = buffer.data(ih + 185);
    const auto *ih_186 = buffer.data(ih + 186);
    const auto *ih_187 = buffer.data(ih + 187);
    const auto *ih_188 = buffer.data(ih + 188);
    const auto *ih_189 = buffer.data(ih + 189);
    const auto *ih_190 = buffer.data(ih + 190);
    const auto *ih_191 = buffer.data(ih + 191);
    const auto *ih_192 = buffer.data(ih + 192);
    const auto *ih_193 = buffer.data(ih + 193);
    const auto *ih_194 = buffer.data(ih + 194);
    const auto *ih_195 = buffer.data(ih + 195);
    const auto *ih_196 = buffer.data(ih + 196);
    const auto *ih_197 = buffer.data(ih + 197);
    const auto *ih_198 = buffer.data(ih + 198);
    const auto *ih_199 = buffer.data(ih + 199);
    const auto *ih_200 = buffer.data(ih + 200);
    const auto *ih_201 = buffer.data(ih + 201);
    const auto *ih_202 = buffer.data(ih + 202);
    const auto *ih_203 = buffer.data(ih + 203);
    const auto *ih_204 = buffer.data(ih + 204);
    const auto *ih_205 = buffer.data(ih + 205);
    const auto *ih_206 = buffer.data(ih + 206);
    const auto *ih_207 = buffer.data(ih + 207);
    const auto *ih_208 = buffer.data(ih + 208);
    const auto *ih_209 = buffer.data(ih + 209);
    const auto *ih_231 = buffer.data(ih + 231);
    const auto *ih_232 = buffer.data(ih + 232);
    const auto *ih_233 = buffer.data(ih + 233);
    const auto *ih_234 = buffer.data(ih + 234);
    const auto *ih_235 = buffer.data(ih + 235);
    const auto *ih_236 = buffer.data(ih + 236);
    const auto *ih_237 = buffer.data(ih + 237);
    const auto *ih_238 = buffer.data(ih + 238);
    const auto *ih_239 = buffer.data(ih + 239);
    const auto *ih_240 = buffer.data(ih + 240);
    const auto *ih_241 = buffer.data(ih + 241);
    const auto *ih_242 = buffer.data(ih + 242);
    const auto *ih_243 = buffer.data(ih + 243);
    const auto *ih_244 = buffer.data(ih + 244);
    const auto *ih_245 = buffer.data(ih + 245);
    const auto *ih_246 = buffer.data(ih + 246);
    const auto *ih_247 = buffer.data(ih + 247);
    const auto *ih_248 = buffer.data(ih + 248);
    const auto *ih_249 = buffer.data(ih + 249);
    const auto *ih_250 = buffer.data(ih + 250);
    const auto *ih_251 = buffer.data(ih + 251);
    const auto *ih_252 = buffer.data(ih + 252);
    const auto *ih_253 = buffer.data(ih + 253);
    const auto *ih_254 = buffer.data(ih + 254);
    const auto *ih_255 = buffer.data(ih + 255);
    const auto *ih_256 = buffer.data(ih + 256);
    const auto *ih_257 = buffer.data(ih + 257);
    const auto *ih_258 = buffer.data(ih + 258);
    const auto *ih_259 = buffer.data(ih + 259);
    const auto *ih_260 = buffer.data(ih + 260);
    const auto *ih_261 = buffer.data(ih + 261);
    const auto *ih_262 = buffer.data(ih + 262);
    const auto *ih_263 = buffer.data(ih + 263);
    const auto *ih_264 = buffer.data(ih + 264);
    const auto *ih_265 = buffer.data(ih + 265);
    const auto *ih_266 = buffer.data(ih + 266);
    const auto *ih_267 = buffer.data(ih + 267);
    const auto *ih_268 = buffer.data(ih + 268);
    const auto *ih_269 = buffer.data(ih + 269);
    const auto *ih_270 = buffer.data(ih + 270);
    const auto *ih_271 = buffer.data(ih + 271);
    const auto *ih_272 = buffer.data(ih + 272);
    const auto *ih_273 = buffer.data(ih + 273);
    const auto *ih_274 = buffer.data(ih + 274);
    const auto *ih_275 = buffer.data(ih + 275);
    const auto *ih_276 = buffer.data(ih + 276);
    const auto *ih_277 = buffer.data(ih + 277);
    const auto *ih_278 = buffer.data(ih + 278);
    const auto *ih_279 = buffer.data(ih + 279);
    const auto *ih_280 = buffer.data(ih + 280);
    const auto *ih_281 = buffer.data(ih + 281);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ih_42, ih_43, ih_44, ih_45, \
                         ih_46, ih_47, ih_48, ih_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ih_42[k];

        t_1[k] = f_0 * ih_43[k];

        t_2[k] = f_0 * ih_44[k];

        t_3[k] = f_0 * ih_45[k];

        t_4[k] = f_0 * ih_46[k];

        t_5[k] = f_0 * ih_47[k];

        t_6[k] = f_0 * ih_48[k];

        t_7[k] = f_0 * ih_49[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ih_50, ih_51, ih_52, \
                         ih_53, ih_54, ih_55, ih_56, ih_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ih_50[k];

        t_9[k] = f_0 * ih_51[k];

        t_10[k] = f_0 * ih_52[k];

        t_11[k] = f_0 * ih_53[k];

        t_12[k] = f_0 * ih_54[k];

        t_13[k] = f_0 * ih_55[k];

        t_14[k] = f_0 * ih_56[k];

        t_15[k] = f_0 * ih_57[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, ih_58, ih_59, ih_60, \
                         ih_61, ih_62, ih_84, ih_85, ih_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ih_58[k];

        t_17[k] = f_0 * ih_59[k];

        t_18[k] = f_0 * ih_60[k];

        t_19[k] = f_0 * ih_61[k];

        t_20[k] = f_0 * ih_62[k];

        t_21[k] = f_0 * ih_84[k];

        t_22[k] = f_0 * ih_85[k];

        t_23[k] = f_0 * ih_86[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, ih_87, ih_88, ih_89, \
                         ih_90, ih_91, ih_92, ih_93, ih_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * ih_87[k];

        t_25[k] = f_0 * ih_88[k];

        t_26[k] = f_0 * ih_89[k];

        t_27[k] = f_0 * ih_90[k];

        t_28[k] = f_0 * ih_91[k];

        t_29[k] = f_0 * ih_92[k];

        t_30[k] = f_0 * ih_93[k];

        t_31[k] = f_0 * ih_94[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, ih_95, ih_96, ih_97, \
                         ih_98, ih_99, ih_100, ih_101, ih_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ih_95[k];

        t_33[k] = f_0 * ih_96[k];

        t_34[k] = f_0 * ih_97[k];

        t_35[k] = f_0 * ih_98[k];

        t_36[k] = f_0 * ih_99[k];

        t_37[k] = f_0 * ih_100[k];

        t_38[k] = f_0 * ih_101[k];

        t_39[k] = f_0 * ih_102[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, gh_0, gh_1, gh_2, gh_3, ih_103, \
                         ih_104, ih_105, ih_106, ih_107, ih_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * ih_103[k];

        t_41[k] = f_0 * ih_104[k];

        t_42[k] = -gh_0[k]
                  + f_0 * ih_105[k];

        t_43[k] = -gh_1[k]
                  + f_0 * ih_106[k];

        t_44[k] = -gh_2[k]
                  + f_0 * ih_107[k];

        t_45[k] = -gh_3[k]
                  + f_0 * ih_108[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, gh_4, gh_5, gh_6, gh_7, gh_8, ih_109, \
                         ih_110, ih_111, ih_112, ih_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -gh_4[k]
                  + f_0 * ih_109[k];

        t_47[k] = -gh_5[k]
                  + f_0 * ih_110[k];

        t_48[k] = -gh_6[k]
                  + f_0 * ih_111[k];

        t_49[k] = -gh_7[k]
                  + f_0 * ih_112[k];

        t_50[k] = -gh_8[k]
                  + f_0 * ih_113[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, gh_9, gh_10, gh_11, gh_12, gh_13, \
                         ih_114, ih_115, ih_116, ih_117, ih_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -gh_9[k]
                  + f_0 * ih_114[k];

        t_52[k] = -gh_10[k]
                  + f_0 * ih_115[k];

        t_53[k] = -gh_11[k]
                  + f_0 * ih_116[k];

        t_54[k] = -gh_12[k]
                  + f_0 * ih_117[k];

        t_55[k] = -gh_13[k]
                  + f_0 * ih_118[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, gh_14, gh_15, gh_16, gh_17, gh_18, \
                         ih_119, ih_120, ih_121, ih_122, ih_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -gh_14[k]
                  + f_0 * ih_119[k];

        t_57[k] = -gh_15[k]
                  + f_0 * ih_120[k];

        t_58[k] = -gh_16[k]
                  + f_0 * ih_121[k];

        t_59[k] = -gh_17[k]
                  + f_0 * ih_122[k];

        t_60[k] = -gh_18[k]
                  + f_0 * ih_123[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, gh_19, gh_20, ih_124, \
                         ih_125, ih_147, ih_148, ih_149, ih_150, \
                         ih_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -gh_19[k]
                  + f_0 * ih_124[k];

        t_62[k] = -gh_20[k]
                  + f_0 * ih_125[k];

        t_63[k] = f_0 * ih_147[k];

        t_64[k] = f_0 * ih_148[k];

        t_65[k] = f_0 * ih_149[k];

        t_66[k] = f_0 * ih_150[k];

        t_67[k] = f_0 * ih_151[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, t_74, t_75, ih_152, ih_153, \
                         ih_154, ih_155, ih_156, ih_157, ih_158, \
                         ih_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * ih_152[k];

        t_69[k] = f_0 * ih_153[k];

        t_70[k] = f_0 * ih_154[k];

        t_71[k] = f_0 * ih_155[k];

        t_72[k] = f_0 * ih_156[k];

        t_73[k] = f_0 * ih_157[k];

        t_74[k] = f_0 * ih_158[k];

        t_75[k] = f_0 * ih_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, t_83, ih_160, ih_161, \
                         ih_162, ih_163, ih_164, ih_165, ih_166, \
                         ih_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_0 * ih_160[k];

        t_77[k] = f_0 * ih_161[k];

        t_78[k] = f_0 * ih_162[k];

        t_79[k] = f_0 * ih_163[k];

        t_80[k] = f_0 * ih_164[k];

        t_81[k] = f_0 * ih_165[k];

        t_82[k] = f_0 * ih_166[k];

        t_83[k] = f_0 * ih_167[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, gh_21, gh_22, gh_23, gh_24, gh_25, \
                         ih_168, ih_169, ih_170, ih_171, ih_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -gh_21[k]
                  + f_0 * ih_168[k];

        t_85[k] = -gh_22[k]
                  + f_0 * ih_169[k];

        t_86[k] = -gh_23[k]
                  + f_0 * ih_170[k];

        t_87[k] = -gh_24[k]
                  + f_0 * ih_171[k];

        t_88[k] = -gh_25[k]
                  + f_0 * ih_172[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, gh_26, gh_27, gh_28, gh_29, gh_30, \
                         ih_173, ih_174, ih_175, ih_176, ih_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -gh_26[k]
                  + f_0 * ih_173[k];

        t_90[k] = -gh_27[k]
                  + f_0 * ih_174[k];

        t_91[k] = -gh_28[k]
                  + f_0 * ih_175[k];

        t_92[k] = -gh_29[k]
                  + f_0 * ih_176[k];

        t_93[k] = -gh_30[k]
                  + f_0 * ih_177[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, gh_31, gh_32, gh_33, gh_34, gh_35, \
                         ih_178, ih_179, ih_180, ih_181, ih_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -gh_31[k]
                  + f_0 * ih_178[k];

        t_95[k] = -gh_32[k]
                  + f_0 * ih_179[k];

        t_96[k] = -gh_33[k]
                  + f_0 * ih_180[k];

        t_97[k] = -gh_34[k]
                  + f_0 * ih_181[k];

        t_98[k] = -gh_35[k]
                  + f_0 * ih_182[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, gh_36, gh_37, gh_38, gh_39, gh_40, \
                         ih_183, ih_184, ih_185, ih_186, ih_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -gh_36[k]
                  + f_0 * ih_183[k];

        t_100[k] = -gh_37[k]
                   + f_0 * ih_184[k];

        t_101[k] = -gh_38[k]
                   + f_0 * ih_185[k];

        t_102[k] = -gh_39[k]
                   + f_0 * ih_186[k];

        t_103[k] = -gh_40[k]
                   + f_0 * ih_187[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, gh_41, gh_42, gh_43, gh_44, gh_45, \
                         ih_188, ih_189, ih_190, ih_191, ih_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -gh_41[k]
                   + f_0 * ih_188[k];

        t_105[k] = -2.0 * gh_42[k]
                   + f_0 * ih_189[k];

        t_106[k] = -2.0 * gh_43[k]
                   + f_0 * ih_190[k];

        t_107[k] = -2.0 * gh_44[k]
                   + f_0 * ih_191[k];

        t_108[k] = -2.0 * gh_45[k]
                   + f_0 * ih_192[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, gh_46, gh_47, gh_48, gh_49, gh_50, \
                         ih_193, ih_194, ih_195, ih_196, ih_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -2.0 * gh_46[k]
                   + f_0 * ih_193[k];

        t_110[k] = -2.0 * gh_47[k]
                   + f_0 * ih_194[k];

        t_111[k] = -2.0 * gh_48[k]
                   + f_0 * ih_195[k];

        t_112[k] = -2.0 * gh_49[k]
                   + f_0 * ih_196[k];

        t_113[k] = -2.0 * gh_50[k]
                   + f_0 * ih_197[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, gh_51, gh_52, gh_53, gh_54, gh_55, \
                         ih_198, ih_199, ih_200, ih_201, ih_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -2.0 * gh_51[k]
                   + f_0 * ih_198[k];

        t_115[k] = -2.0 * gh_52[k]
                   + f_0 * ih_199[k];

        t_116[k] = -2.0 * gh_53[k]
                   + f_0 * ih_200[k];

        t_117[k] = -2.0 * gh_54[k]
                   + f_0 * ih_201[k];

        t_118[k] = -2.0 * gh_55[k]
                   + f_0 * ih_202[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, gh_56, gh_57, gh_58, gh_59, gh_60, \
                         ih_203, ih_204, ih_205, ih_206, ih_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -2.0 * gh_56[k]
                   + f_0 * ih_203[k];

        t_120[k] = -2.0 * gh_57[k]
                   + f_0 * ih_204[k];

        t_121[k] = -2.0 * gh_58[k]
                   + f_0 * ih_205[k];

        t_122[k] = -2.0 * gh_59[k]
                   + f_0 * ih_206[k];

        t_123[k] = -2.0 * gh_60[k]
                   + f_0 * ih_207[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, t_129, t_130, gh_61, gh_62, \
                         ih_208, ih_209, ih_231, ih_232, ih_233, ih_234, \
                         ih_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * gh_61[k]
                   + f_0 * ih_208[k];

        t_125[k] = -2.0 * gh_62[k]
                   + f_0 * ih_209[k];

        t_126[k] = f_0 * ih_231[k];

        t_127[k] = f_0 * ih_232[k];

        t_128[k] = f_0 * ih_233[k];

        t_129[k] = f_0 * ih_234[k];

        t_130[k] = f_0 * ih_235[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, t_136, t_137, t_138, ih_236, \
                         ih_237, ih_238, ih_239, ih_240, ih_241, ih_242, \
                         ih_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_0 * ih_236[k];

        t_132[k] = f_0 * ih_237[k];

        t_133[k] = f_0 * ih_238[k];

        t_134[k] = f_0 * ih_239[k];

        t_135[k] = f_0 * ih_240[k];

        t_136[k] = f_0 * ih_241[k];

        t_137[k] = f_0 * ih_242[k];

        t_138[k] = f_0 * ih_243[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, t_146, ih_244, \
                         ih_245, ih_246, ih_247, ih_248, ih_249, ih_250, \
                         ih_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_0 * ih_244[k];

        t_140[k] = f_0 * ih_245[k];

        t_141[k] = f_0 * ih_246[k];

        t_142[k] = f_0 * ih_247[k];

        t_143[k] = f_0 * ih_248[k];

        t_144[k] = f_0 * ih_249[k];

        t_145[k] = f_0 * ih_250[k];

        t_146[k] = f_0 * ih_251[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, gh_63, gh_64, gh_65, gh_66, gh_67, \
                         ih_252, ih_253, ih_254, ih_255, ih_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -gh_63[k]
                   + f_0 * ih_252[k];

        t_148[k] = -gh_64[k]
                   + f_0 * ih_253[k];

        t_149[k] = -gh_65[k]
                   + f_0 * ih_254[k];

        t_150[k] = -gh_66[k]
                   + f_0 * ih_255[k];

        t_151[k] = -gh_67[k]
                   + f_0 * ih_256[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, gh_68, gh_69, gh_70, gh_71, gh_72, \
                         ih_257, ih_258, ih_259, ih_260, ih_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -gh_68[k]
                   + f_0 * ih_257[k];

        t_153[k] = -gh_69[k]
                   + f_0 * ih_258[k];

        t_154[k] = -gh_70[k]
                   + f_0 * ih_259[k];

        t_155[k] = -gh_71[k]
                   + f_0 * ih_260[k];

        t_156[k] = -gh_72[k]
                   + f_0 * ih_261[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, gh_73, gh_74, gh_75, gh_76, gh_77, \
                         ih_262, ih_263, ih_264, ih_265, ih_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -gh_73[k]
                   + f_0 * ih_262[k];

        t_158[k] = -gh_74[k]
                   + f_0 * ih_263[k];

        t_159[k] = -gh_75[k]
                   + f_0 * ih_264[k];

        t_160[k] = -gh_76[k]
                   + f_0 * ih_265[k];

        t_161[k] = -gh_77[k]
                   + f_0 * ih_266[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, gh_78, gh_79, gh_80, gh_81, gh_82, \
                         ih_267, ih_268, ih_269, ih_270, ih_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -gh_78[k]
                   + f_0 * ih_267[k];

        t_163[k] = -gh_79[k]
                   + f_0 * ih_268[k];

        t_164[k] = -gh_80[k]
                   + f_0 * ih_269[k];

        t_165[k] = -gh_81[k]
                   + f_0 * ih_270[k];

        t_166[k] = -gh_82[k]
                   + f_0 * ih_271[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, gh_83, gh_84, gh_85, gh_86, gh_87, \
                         ih_272, ih_273, ih_274, ih_275, ih_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -gh_83[k]
                   + f_0 * ih_272[k];

        t_168[k] = -2.0 * gh_84[k]
                   + f_0 * ih_273[k];

        t_169[k] = -2.0 * gh_85[k]
                   + f_0 * ih_274[k];

        t_170[k] = -2.0 * gh_86[k]
                   + f_0 * ih_275[k];

        t_171[k] = -2.0 * gh_87[k]
                   + f_0 * ih_276[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, gh_88, gh_89, gh_90, gh_91, gh_92, \
                         ih_277, ih_278, ih_279, ih_280, ih_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -2.0 * gh_88[k]
                   + f_0 * ih_277[k];

        t_173[k] = -2.0 * gh_89[k]
                   + f_0 * ih_278[k];

        t_174[k] = -2.0 * gh_90[k]
                   + f_0 * ih_279[k];

        t_175[k] = -2.0 * gh_91[k]
                   + f_0 * ih_280[k];

        t_176[k] = -2.0 * gh_92[k]
                   + f_0 * ih_281[k];
    }
}

static auto
compute_prim_geom_10_hh_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ih,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);

    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_115 = buffer.data(gh + 115);
    const auto *gh_116 = buffer.data(gh + 116);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_118 = buffer.data(gh + 118);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_124 = buffer.data(gh + 124);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_137 = buffer.data(gh + 137);
    const auto *gh_138 = buffer.data(gh + 138);
    const auto *gh_139 = buffer.data(gh + 139);
    const auto *gh_140 = buffer.data(gh + 140);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_142 = buffer.data(gh + 142);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_148 = buffer.data(gh + 148);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_190 = buffer.data(gh + 190);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_193 = buffer.data(gh + 193);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_195 = buffer.data(gh + 195);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_197 = buffer.data(gh + 197);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_199 = buffer.data(gh + 199);
    const auto *gh_200 = buffer.data(gh + 200);
    const auto *gh_201 = buffer.data(gh + 201);
    const auto *gh_202 = buffer.data(gh + 202);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_208 = buffer.data(gh + 208);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);

    const auto *ih_282 = buffer.data(ih + 282);
    const auto *ih_283 = buffer.data(ih + 283);
    const auto *ih_284 = buffer.data(ih + 284);
    const auto *ih_285 = buffer.data(ih + 285);
    const auto *ih_286 = buffer.data(ih + 286);
    const auto *ih_287 = buffer.data(ih + 287);
    const auto *ih_288 = buffer.data(ih + 288);
    const auto *ih_289 = buffer.data(ih + 289);
    const auto *ih_290 = buffer.data(ih + 290);
    const auto *ih_291 = buffer.data(ih + 291);
    const auto *ih_292 = buffer.data(ih + 292);
    const auto *ih_293 = buffer.data(ih + 293);
    const auto *ih_294 = buffer.data(ih + 294);
    const auto *ih_295 = buffer.data(ih + 295);
    const auto *ih_296 = buffer.data(ih + 296);
    const auto *ih_297 = buffer.data(ih + 297);
    const auto *ih_298 = buffer.data(ih + 298);
    const auto *ih_299 = buffer.data(ih + 299);
    const auto *ih_300 = buffer.data(ih + 300);
    const auto *ih_301 = buffer.data(ih + 301);
    const auto *ih_302 = buffer.data(ih + 302);
    const auto *ih_303 = buffer.data(ih + 303);
    const auto *ih_304 = buffer.data(ih + 304);
    const auto *ih_305 = buffer.data(ih + 305);
    const auto *ih_306 = buffer.data(ih + 306);
    const auto *ih_307 = buffer.data(ih + 307);
    const auto *ih_308 = buffer.data(ih + 308);
    const auto *ih_309 = buffer.data(ih + 309);
    const auto *ih_310 = buffer.data(ih + 310);
    const auto *ih_311 = buffer.data(ih + 311);
    const auto *ih_312 = buffer.data(ih + 312);
    const auto *ih_313 = buffer.data(ih + 313);
    const auto *ih_314 = buffer.data(ih + 314);
    const auto *ih_336 = buffer.data(ih + 336);
    const auto *ih_337 = buffer.data(ih + 337);
    const auto *ih_338 = buffer.data(ih + 338);
    const auto *ih_339 = buffer.data(ih + 339);
    const auto *ih_340 = buffer.data(ih + 340);
    const auto *ih_341 = buffer.data(ih + 341);
    const auto *ih_342 = buffer.data(ih + 342);
    const auto *ih_343 = buffer.data(ih + 343);
    const auto *ih_344 = buffer.data(ih + 344);
    const auto *ih_345 = buffer.data(ih + 345);
    const auto *ih_346 = buffer.data(ih + 346);
    const auto *ih_347 = buffer.data(ih + 347);
    const auto *ih_348 = buffer.data(ih + 348);
    const auto *ih_349 = buffer.data(ih + 349);
    const auto *ih_350 = buffer.data(ih + 350);
    const auto *ih_351 = buffer.data(ih + 351);
    const auto *ih_352 = buffer.data(ih + 352);
    const auto *ih_353 = buffer.data(ih + 353);
    const auto *ih_354 = buffer.data(ih + 354);
    const auto *ih_355 = buffer.data(ih + 355);
    const auto *ih_356 = buffer.data(ih + 356);
    const auto *ih_357 = buffer.data(ih + 357);
    const auto *ih_358 = buffer.data(ih + 358);
    const auto *ih_359 = buffer.data(ih + 359);
    const auto *ih_360 = buffer.data(ih + 360);
    const auto *ih_361 = buffer.data(ih + 361);
    const auto *ih_362 = buffer.data(ih + 362);
    const auto *ih_363 = buffer.data(ih + 363);
    const auto *ih_364 = buffer.data(ih + 364);
    const auto *ih_365 = buffer.data(ih + 365);
    const auto *ih_366 = buffer.data(ih + 366);
    const auto *ih_367 = buffer.data(ih + 367);
    const auto *ih_368 = buffer.data(ih + 368);
    const auto *ih_369 = buffer.data(ih + 369);
    const auto *ih_370 = buffer.data(ih + 370);
    const auto *ih_371 = buffer.data(ih + 371);
    const auto *ih_372 = buffer.data(ih + 372);
    const auto *ih_373 = buffer.data(ih + 373);
    const auto *ih_374 = buffer.data(ih + 374);
    const auto *ih_375 = buffer.data(ih + 375);
    const auto *ih_376 = buffer.data(ih + 376);
    const auto *ih_377 = buffer.data(ih + 377);
    const auto *ih_378 = buffer.data(ih + 378);
    const auto *ih_379 = buffer.data(ih + 379);
    const auto *ih_380 = buffer.data(ih + 380);
    const auto *ih_381 = buffer.data(ih + 381);
    const auto *ih_382 = buffer.data(ih + 382);
    const auto *ih_383 = buffer.data(ih + 383);
    const auto *ih_384 = buffer.data(ih + 384);
    const auto *ih_385 = buffer.data(ih + 385);
    const auto *ih_386 = buffer.data(ih + 386);
    const auto *ih_387 = buffer.data(ih + 387);
    const auto *ih_388 = buffer.data(ih + 388);
    const auto *ih_389 = buffer.data(ih + 389);
    const auto *ih_390 = buffer.data(ih + 390);
    const auto *ih_391 = buffer.data(ih + 391);
    const auto *ih_392 = buffer.data(ih + 392);
    const auto *ih_393 = buffer.data(ih + 393);
    const auto *ih_394 = buffer.data(ih + 394);
    const auto *ih_395 = buffer.data(ih + 395);
    const auto *ih_396 = buffer.data(ih + 396);
    const auto *ih_397 = buffer.data(ih + 397);
    const auto *ih_398 = buffer.data(ih + 398);
    const auto *ih_399 = buffer.data(ih + 399);
    const auto *ih_400 = buffer.data(ih + 400);
    const auto *ih_401 = buffer.data(ih + 401);
    const auto *ih_402 = buffer.data(ih + 402);
    const auto *ih_403 = buffer.data(ih + 403);
    const auto *ih_404 = buffer.data(ih + 404);
    const auto *ih_405 = buffer.data(ih + 405);
    const auto *ih_406 = buffer.data(ih + 406);
    const auto *ih_407 = buffer.data(ih + 407);
    const auto *ih_408 = buffer.data(ih + 408);
    const auto *ih_409 = buffer.data(ih + 409);
    const auto *ih_410 = buffer.data(ih + 410);
    const auto *ih_411 = buffer.data(ih + 411);
    const auto *ih_412 = buffer.data(ih + 412);
    const auto *ih_413 = buffer.data(ih + 413);
    const auto *ih_414 = buffer.data(ih + 414);
    const auto *ih_415 = buffer.data(ih + 415);
    const auto *ih_416 = buffer.data(ih + 416);
    const auto *ih_417 = buffer.data(ih + 417);
    const auto *ih_418 = buffer.data(ih + 418);
    const auto *ih_419 = buffer.data(ih + 419);
    const auto *ih_420 = buffer.data(ih + 420);
    const auto *ih_421 = buffer.data(ih + 421);
    const auto *ih_422 = buffer.data(ih + 422);
    const auto *ih_423 = buffer.data(ih + 423);
    const auto *ih_424 = buffer.data(ih + 424);
    const auto *ih_425 = buffer.data(ih + 425);
    const auto *ih_426 = buffer.data(ih + 426);
    const auto *ih_427 = buffer.data(ih + 427);
    const auto *ih_428 = buffer.data(ih + 428);
    const auto *ih_429 = buffer.data(ih + 429);
    const auto *ih_430 = buffer.data(ih + 430);
    const auto *ih_431 = buffer.data(ih + 431);
    const auto *ih_432 = buffer.data(ih + 432);
    const auto *ih_433 = buffer.data(ih + 433);
    const auto *ih_434 = buffer.data(ih + 434);
    const auto *ih_435 = buffer.data(ih + 435);
    const auto *ih_436 = buffer.data(ih + 436);
    const auto *ih_437 = buffer.data(ih + 437);
    const auto *ih_438 = buffer.data(ih + 438);
    const auto *ih_439 = buffer.data(ih + 439);
    const auto *ih_440 = buffer.data(ih + 440);
    const auto *ih_462 = buffer.data(ih + 462);
    const auto *ih_463 = buffer.data(ih + 463);
    const auto *ih_464 = buffer.data(ih + 464);
    const auto *ih_465 = buffer.data(ih + 465);
    const auto *ih_466 = buffer.data(ih + 466);
    const auto *ih_467 = buffer.data(ih + 467);
    const auto *ih_468 = buffer.data(ih + 468);
    const auto *ih_469 = buffer.data(ih + 469);
    const auto *ih_470 = buffer.data(ih + 470);
    const auto *ih_471 = buffer.data(ih + 471);
    const auto *ih_472 = buffer.data(ih + 472);
    const auto *ih_473 = buffer.data(ih + 473);
    const auto *ih_474 = buffer.data(ih + 474);
    const auto *ih_475 = buffer.data(ih + 475);
    const auto *ih_476 = buffer.data(ih + 476);
    const auto *ih_477 = buffer.data(ih + 477);
    const auto *ih_478 = buffer.data(ih + 478);
    const auto *ih_479 = buffer.data(ih + 479);
    const auto *ih_480 = buffer.data(ih + 480);
    const auto *ih_481 = buffer.data(ih + 481);
    const auto *ih_482 = buffer.data(ih + 482);
    const auto *ih_483 = buffer.data(ih + 483);
    const auto *ih_484 = buffer.data(ih + 484);

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, gh_93, gh_94, gh_95, gh_96, gh_97, \
                         ih_282, ih_283, ih_284, ih_285, ih_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -2.0 * gh_93[k]
                   + f_0 * ih_282[k];

        t_178[k] = -2.0 * gh_94[k]
                   + f_0 * ih_283[k];

        t_179[k] = -2.0 * gh_95[k]
                   + f_0 * ih_284[k];

        t_180[k] = -2.0 * gh_96[k]
                   + f_0 * ih_285[k];

        t_181[k] = -2.0 * gh_97[k]
                   + f_0 * ih_286[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, gh_98, gh_99, gh_100, gh_101, \
                         gh_102, ih_287, ih_288, ih_289, ih_290, \
                         ih_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * gh_98[k]
                   + f_0 * ih_287[k];

        t_183[k] = -2.0 * gh_99[k]
                   + f_0 * ih_288[k];

        t_184[k] = -2.0 * gh_100[k]
                   + f_0 * ih_289[k];

        t_185[k] = -2.0 * gh_101[k]
                   + f_0 * ih_290[k];

        t_186[k] = -2.0 * gh_102[k]
                   + f_0 * ih_291[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, gh_103, gh_104, gh_105, gh_106, \
                         gh_107, ih_292, ih_293, ih_294, ih_295, \
                         ih_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * gh_103[k]
                   + f_0 * ih_292[k];

        t_188[k] = -2.0 * gh_104[k]
                   + f_0 * ih_293[k];

        t_189[k] = -3.0 * gh_105[k]
                   + f_0 * ih_294[k];

        t_190[k] = -3.0 * gh_106[k]
                   + f_0 * ih_295[k];

        t_191[k] = -3.0 * gh_107[k]
                   + f_0 * ih_296[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, gh_108, gh_109, gh_110, gh_111, \
                         gh_112, ih_297, ih_298, ih_299, ih_300, \
                         ih_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -3.0 * gh_108[k]
                   + f_0 * ih_297[k];

        t_193[k] = -3.0 * gh_109[k]
                   + f_0 * ih_298[k];

        t_194[k] = -3.0 * gh_110[k]
                   + f_0 * ih_299[k];

        t_195[k] = -3.0 * gh_111[k]
                   + f_0 * ih_300[k];

        t_196[k] = -3.0 * gh_112[k]
                   + f_0 * ih_301[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, gh_113, gh_114, gh_115, gh_116, \
                         gh_117, ih_302, ih_303, ih_304, ih_305, \
                         ih_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * gh_113[k]
                   + f_0 * ih_302[k];

        t_198[k] = -3.0 * gh_114[k]
                   + f_0 * ih_303[k];

        t_199[k] = -3.0 * gh_115[k]
                   + f_0 * ih_304[k];

        t_200[k] = -3.0 * gh_116[k]
                   + f_0 * ih_305[k];

        t_201[k] = -3.0 * gh_117[k]
                   + f_0 * ih_306[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, gh_118, gh_119, gh_120, gh_121, \
                         gh_122, ih_307, ih_308, ih_309, ih_310, \
                         ih_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * gh_118[k]
                   + f_0 * ih_307[k];

        t_203[k] = -3.0 * gh_119[k]
                   + f_0 * ih_308[k];

        t_204[k] = -3.0 * gh_120[k]
                   + f_0 * ih_309[k];

        t_205[k] = -3.0 * gh_121[k]
                   + f_0 * ih_310[k];

        t_206[k] = -3.0 * gh_122[k]
                   + f_0 * ih_311[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, gh_123, gh_124, gh_125, \
                         ih_312, ih_313, ih_314, ih_336, ih_337, \
                         ih_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * gh_123[k]
                   + f_0 * ih_312[k];

        t_208[k] = -3.0 * gh_124[k]
                   + f_0 * ih_313[k];

        t_209[k] = -3.0 * gh_125[k]
                   + f_0 * ih_314[k];

        t_210[k] = f_0 * ih_336[k];

        t_211[k] = f_0 * ih_337[k];

        t_212[k] = f_0 * ih_338[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, t_219, t_220, ih_339, \
                         ih_340, ih_341, ih_342, ih_343, ih_344, ih_345, \
                         ih_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_0 * ih_339[k];

        t_214[k] = f_0 * ih_340[k];

        t_215[k] = f_0 * ih_341[k];

        t_216[k] = f_0 * ih_342[k];

        t_217[k] = f_0 * ih_343[k];

        t_218[k] = f_0 * ih_344[k];

        t_219[k] = f_0 * ih_345[k];

        t_220[k] = f_0 * ih_346[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, t_226, t_227, t_228, ih_347, \
                         ih_348, ih_349, ih_350, ih_351, ih_352, ih_353, \
                         ih_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_0 * ih_347[k];

        t_222[k] = f_0 * ih_348[k];

        t_223[k] = f_0 * ih_349[k];

        t_224[k] = f_0 * ih_350[k];

        t_225[k] = f_0 * ih_351[k];

        t_226[k] = f_0 * ih_352[k];

        t_227[k] = f_0 * ih_353[k];

        t_228[k] = f_0 * ih_354[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, t_234, gh_126, gh_127, gh_128, \
                         gh_129, ih_355, ih_356, ih_357, ih_358, ih_359, \
                         ih_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * ih_355[k];

        t_230[k] = f_0 * ih_356[k];

        t_231[k] = -gh_126[k]
                   + f_0 * ih_357[k];

        t_232[k] = -gh_127[k]
                   + f_0 * ih_358[k];

        t_233[k] = -gh_128[k]
                   + f_0 * ih_359[k];

        t_234[k] = -gh_129[k]
                   + f_0 * ih_360[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, gh_130, gh_131, gh_132, gh_133, \
                         gh_134, ih_361, ih_362, ih_363, ih_364, \
                         ih_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -gh_130[k]
                   + f_0 * ih_361[k];

        t_236[k] = -gh_131[k]
                   + f_0 * ih_362[k];

        t_237[k] = -gh_132[k]
                   + f_0 * ih_363[k];

        t_238[k] = -gh_133[k]
                   + f_0 * ih_364[k];

        t_239[k] = -gh_134[k]
                   + f_0 * ih_365[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, gh_135, gh_136, gh_137, gh_138, \
                         gh_139, ih_366, ih_367, ih_368, ih_369, \
                         ih_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -gh_135[k]
                   + f_0 * ih_366[k];

        t_241[k] = -gh_136[k]
                   + f_0 * ih_367[k];

        t_242[k] = -gh_137[k]
                   + f_0 * ih_368[k];

        t_243[k] = -gh_138[k]
                   + f_0 * ih_369[k];

        t_244[k] = -gh_139[k]
                   + f_0 * ih_370[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, gh_140, gh_141, gh_142, gh_143, \
                         gh_144, ih_371, ih_372, ih_373, ih_374, \
                         ih_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -gh_140[k]
                   + f_0 * ih_371[k];

        t_246[k] = -gh_141[k]
                   + f_0 * ih_372[k];

        t_247[k] = -gh_142[k]
                   + f_0 * ih_373[k];

        t_248[k] = -gh_143[k]
                   + f_0 * ih_374[k];

        t_249[k] = -gh_144[k]
                   + f_0 * ih_375[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, gh_145, gh_146, gh_147, gh_148, \
                         gh_149, ih_376, ih_377, ih_378, ih_379, \
                         ih_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -gh_145[k]
                   + f_0 * ih_376[k];

        t_251[k] = -gh_146[k]
                   + f_0 * ih_377[k];

        t_252[k] = -2.0 * gh_147[k]
                   + f_0 * ih_378[k];

        t_253[k] = -2.0 * gh_148[k]
                   + f_0 * ih_379[k];

        t_254[k] = -2.0 * gh_149[k]
                   + f_0 * ih_380[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, gh_150, gh_151, gh_152, gh_153, \
                         gh_154, ih_381, ih_382, ih_383, ih_384, \
                         ih_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * gh_150[k]
                   + f_0 * ih_381[k];

        t_256[k] = -2.0 * gh_151[k]
                   + f_0 * ih_382[k];

        t_257[k] = -2.0 * gh_152[k]
                   + f_0 * ih_383[k];

        t_258[k] = -2.0 * gh_153[k]
                   + f_0 * ih_384[k];

        t_259[k] = -2.0 * gh_154[k]
                   + f_0 * ih_385[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, gh_155, gh_156, gh_157, gh_158, \
                         gh_159, ih_386, ih_387, ih_388, ih_389, \
                         ih_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * gh_155[k]
                   + f_0 * ih_386[k];

        t_261[k] = -2.0 * gh_156[k]
                   + f_0 * ih_387[k];

        t_262[k] = -2.0 * gh_157[k]
                   + f_0 * ih_388[k];

        t_263[k] = -2.0 * gh_158[k]
                   + f_0 * ih_389[k];

        t_264[k] = -2.0 * gh_159[k]
                   + f_0 * ih_390[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, gh_160, gh_161, gh_162, gh_163, \
                         gh_164, ih_391, ih_392, ih_393, ih_394, \
                         ih_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * gh_160[k]
                   + f_0 * ih_391[k];

        t_266[k] = -2.0 * gh_161[k]
                   + f_0 * ih_392[k];

        t_267[k] = -2.0 * gh_162[k]
                   + f_0 * ih_393[k];

        t_268[k] = -2.0 * gh_163[k]
                   + f_0 * ih_394[k];

        t_269[k] = -2.0 * gh_164[k]
                   + f_0 * ih_395[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, gh_165, gh_166, gh_167, gh_168, \
                         gh_169, ih_396, ih_397, ih_398, ih_399, \
                         ih_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * gh_165[k]
                   + f_0 * ih_396[k];

        t_271[k] = -2.0 * gh_166[k]
                   + f_0 * ih_397[k];

        t_272[k] = -2.0 * gh_167[k]
                   + f_0 * ih_398[k];

        t_273[k] = -3.0 * gh_168[k]
                   + f_0 * ih_399[k];

        t_274[k] = -3.0 * gh_169[k]
                   + f_0 * ih_400[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, gh_170, gh_171, gh_172, gh_173, \
                         gh_174, ih_401, ih_402, ih_403, ih_404, \
                         ih_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * gh_170[k]
                   + f_0 * ih_401[k];

        t_276[k] = -3.0 * gh_171[k]
                   + f_0 * ih_402[k];

        t_277[k] = -3.0 * gh_172[k]
                   + f_0 * ih_403[k];

        t_278[k] = -3.0 * gh_173[k]
                   + f_0 * ih_404[k];

        t_279[k] = -3.0 * gh_174[k]
                   + f_0 * ih_405[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, gh_175, gh_176, gh_177, gh_178, \
                         gh_179, ih_406, ih_407, ih_408, ih_409, \
                         ih_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * gh_175[k]
                   + f_0 * ih_406[k];

        t_281[k] = -3.0 * gh_176[k]
                   + f_0 * ih_407[k];

        t_282[k] = -3.0 * gh_177[k]
                   + f_0 * ih_408[k];

        t_283[k] = -3.0 * gh_178[k]
                   + f_0 * ih_409[k];

        t_284[k] = -3.0 * gh_179[k]
                   + f_0 * ih_410[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, gh_180, gh_181, gh_182, gh_183, \
                         gh_184, ih_411, ih_412, ih_413, ih_414, \
                         ih_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * gh_180[k]
                   + f_0 * ih_411[k];

        t_286[k] = -3.0 * gh_181[k]
                   + f_0 * ih_412[k];

        t_287[k] = -3.0 * gh_182[k]
                   + f_0 * ih_413[k];

        t_288[k] = -3.0 * gh_183[k]
                   + f_0 * ih_414[k];

        t_289[k] = -3.0 * gh_184[k]
                   + f_0 * ih_415[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, gh_185, gh_186, gh_187, gh_188, \
                         gh_189, ih_416, ih_417, ih_418, ih_419, \
                         ih_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * gh_185[k]
                   + f_0 * ih_416[k];

        t_291[k] = -3.0 * gh_186[k]
                   + f_0 * ih_417[k];

        t_292[k] = -3.0 * gh_187[k]
                   + f_0 * ih_418[k];

        t_293[k] = -3.0 * gh_188[k]
                   + f_0 * ih_419[k];

        t_294[k] = -4.0 * gh_189[k]
                   + f_0 * ih_420[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, gh_190, gh_191, gh_192, gh_193, \
                         gh_194, ih_421, ih_422, ih_423, ih_424, \
                         ih_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * gh_190[k]
                   + f_0 * ih_421[k];

        t_296[k] = -4.0 * gh_191[k]
                   + f_0 * ih_422[k];

        t_297[k] = -4.0 * gh_192[k]
                   + f_0 * ih_423[k];

        t_298[k] = -4.0 * gh_193[k]
                   + f_0 * ih_424[k];

        t_299[k] = -4.0 * gh_194[k]
                   + f_0 * ih_425[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, gh_195, gh_196, gh_197, gh_198, \
                         gh_199, ih_426, ih_427, ih_428, ih_429, \
                         ih_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * gh_195[k]
                   + f_0 * ih_426[k];

        t_301[k] = -4.0 * gh_196[k]
                   + f_0 * ih_427[k];

        t_302[k] = -4.0 * gh_197[k]
                   + f_0 * ih_428[k];

        t_303[k] = -4.0 * gh_198[k]
                   + f_0 * ih_429[k];

        t_304[k] = -4.0 * gh_199[k]
                   + f_0 * ih_430[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, gh_200, gh_201, gh_202, gh_203, \
                         gh_204, ih_431, ih_432, ih_433, ih_434, \
                         ih_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * gh_200[k]
                   + f_0 * ih_431[k];

        t_306[k] = -4.0 * gh_201[k]
                   + f_0 * ih_432[k];

        t_307[k] = -4.0 * gh_202[k]
                   + f_0 * ih_433[k];

        t_308[k] = -4.0 * gh_203[k]
                   + f_0 * ih_434[k];

        t_309[k] = -4.0 * gh_204[k]
                   + f_0 * ih_435[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, gh_205, gh_206, gh_207, gh_208, \
                         gh_209, ih_436, ih_437, ih_438, ih_439, \
                         ih_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -4.0 * gh_205[k]
                   + f_0 * ih_436[k];

        t_311[k] = -4.0 * gh_206[k]
                   + f_0 * ih_437[k];

        t_312[k] = -4.0 * gh_207[k]
                   + f_0 * ih_438[k];

        t_313[k] = -4.0 * gh_208[k]
                   + f_0 * ih_439[k];

        t_314[k] = -4.0 * gh_209[k]
                   + f_0 * ih_440[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, t_321, t_322, ih_462, \
                         ih_463, ih_464, ih_465, ih_466, ih_467, ih_468, \
                         ih_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_0 * ih_462[k];

        t_316[k] = f_0 * ih_463[k];

        t_317[k] = f_0 * ih_464[k];

        t_318[k] = f_0 * ih_465[k];

        t_319[k] = f_0 * ih_466[k];

        t_320[k] = f_0 * ih_467[k];

        t_321[k] = f_0 * ih_468[k];

        t_322[k] = f_0 * ih_469[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, t_328, t_329, t_330, ih_470, \
                         ih_471, ih_472, ih_473, ih_474, ih_475, ih_476, \
                         ih_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_0 * ih_470[k];

        t_324[k] = f_0 * ih_471[k];

        t_325[k] = f_0 * ih_472[k];

        t_326[k] = f_0 * ih_473[k];

        t_327[k] = f_0 * ih_474[k];

        t_328[k] = f_0 * ih_475[k];

        t_329[k] = f_0 * ih_476[k];

        t_330[k] = f_0 * ih_477[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, t_336, t_337, gh_210, gh_211, \
                         ih_478, ih_479, ih_480, ih_481, ih_482, ih_483, \
                         ih_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_0 * ih_478[k];

        t_332[k] = f_0 * ih_479[k];

        t_333[k] = f_0 * ih_480[k];

        t_334[k] = f_0 * ih_481[k];

        t_335[k] = f_0 * ih_482[k];

        t_336[k] = -gh_210[k]
                   + f_0 * ih_483[k];

        t_337[k] = -gh_211[k]
                   + f_0 * ih_484[k];
    }
}

static auto
compute_prim_geom_10_hh_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gh, const size_t ih,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gh_212 = buffer.data(gh + 212);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_221 = buffer.data(gh + 221);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_223 = buffer.data(gh + 223);
    const auto *gh_224 = buffer.data(gh + 224);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_231 = buffer.data(gh + 231);
    const auto *gh_232 = buffer.data(gh + 232);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_234 = buffer.data(gh + 234);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_237 = buffer.data(gh + 237);
    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_241 = buffer.data(gh + 241);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_273 = buffer.data(gh + 273);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_275 = buffer.data(gh + 275);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_278 = buffer.data(gh + 278);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);
    const auto *gh_282 = buffer.data(gh + 282);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_287 = buffer.data(gh + 287);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_295 = buffer.data(gh + 295);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_298 = buffer.data(gh + 298);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_302 = buffer.data(gh + 302);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_304 = buffer.data(gh + 304);
    const auto *gh_305 = buffer.data(gh + 305);
    const auto *gh_306 = buffer.data(gh + 306);
    const auto *gh_307 = buffer.data(gh + 307);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

    const auto *ih_485 = buffer.data(ih + 485);
    const auto *ih_486 = buffer.data(ih + 486);
    const auto *ih_487 = buffer.data(ih + 487);
    const auto *ih_488 = buffer.data(ih + 488);
    const auto *ih_489 = buffer.data(ih + 489);
    const auto *ih_490 = buffer.data(ih + 490);
    const auto *ih_491 = buffer.data(ih + 491);
    const auto *ih_492 = buffer.data(ih + 492);
    const auto *ih_493 = buffer.data(ih + 493);
    const auto *ih_494 = buffer.data(ih + 494);
    const auto *ih_495 = buffer.data(ih + 495);
    const auto *ih_496 = buffer.data(ih + 496);
    const auto *ih_497 = buffer.data(ih + 497);
    const auto *ih_498 = buffer.data(ih + 498);
    const auto *ih_499 = buffer.data(ih + 499);
    const auto *ih_500 = buffer.data(ih + 500);
    const auto *ih_501 = buffer.data(ih + 501);
    const auto *ih_502 = buffer.data(ih + 502);
    const auto *ih_503 = buffer.data(ih + 503);
    const auto *ih_504 = buffer.data(ih + 504);
    const auto *ih_505 = buffer.data(ih + 505);
    const auto *ih_506 = buffer.data(ih + 506);
    const auto *ih_507 = buffer.data(ih + 507);
    const auto *ih_508 = buffer.data(ih + 508);
    const auto *ih_509 = buffer.data(ih + 509);
    const auto *ih_510 = buffer.data(ih + 510);
    const auto *ih_511 = buffer.data(ih + 511);
    const auto *ih_512 = buffer.data(ih + 512);
    const auto *ih_513 = buffer.data(ih + 513);
    const auto *ih_514 = buffer.data(ih + 514);
    const auto *ih_515 = buffer.data(ih + 515);
    const auto *ih_516 = buffer.data(ih + 516);
    const auto *ih_517 = buffer.data(ih + 517);
    const auto *ih_518 = buffer.data(ih + 518);
    const auto *ih_519 = buffer.data(ih + 519);
    const auto *ih_520 = buffer.data(ih + 520);
    const auto *ih_521 = buffer.data(ih + 521);
    const auto *ih_522 = buffer.data(ih + 522);
    const auto *ih_523 = buffer.data(ih + 523);
    const auto *ih_524 = buffer.data(ih + 524);
    const auto *ih_525 = buffer.data(ih + 525);
    const auto *ih_526 = buffer.data(ih + 526);
    const auto *ih_527 = buffer.data(ih + 527);
    const auto *ih_528 = buffer.data(ih + 528);
    const auto *ih_529 = buffer.data(ih + 529);
    const auto *ih_530 = buffer.data(ih + 530);
    const auto *ih_531 = buffer.data(ih + 531);
    const auto *ih_532 = buffer.data(ih + 532);
    const auto *ih_533 = buffer.data(ih + 533);
    const auto *ih_534 = buffer.data(ih + 534);
    const auto *ih_535 = buffer.data(ih + 535);
    const auto *ih_536 = buffer.data(ih + 536);
    const auto *ih_537 = buffer.data(ih + 537);
    const auto *ih_538 = buffer.data(ih + 538);
    const auto *ih_539 = buffer.data(ih + 539);
    const auto *ih_540 = buffer.data(ih + 540);
    const auto *ih_541 = buffer.data(ih + 541);
    const auto *ih_542 = buffer.data(ih + 542);
    const auto *ih_543 = buffer.data(ih + 543);
    const auto *ih_544 = buffer.data(ih + 544);
    const auto *ih_545 = buffer.data(ih + 545);
    const auto *ih_546 = buffer.data(ih + 546);
    const auto *ih_547 = buffer.data(ih + 547);
    const auto *ih_548 = buffer.data(ih + 548);
    const auto *ih_549 = buffer.data(ih + 549);
    const auto *ih_550 = buffer.data(ih + 550);
    const auto *ih_551 = buffer.data(ih + 551);
    const auto *ih_552 = buffer.data(ih + 552);
    const auto *ih_553 = buffer.data(ih + 553);
    const auto *ih_554 = buffer.data(ih + 554);
    const auto *ih_555 = buffer.data(ih + 555);
    const auto *ih_556 = buffer.data(ih + 556);
    const auto *ih_557 = buffer.data(ih + 557);
    const auto *ih_558 = buffer.data(ih + 558);
    const auto *ih_559 = buffer.data(ih + 559);
    const auto *ih_560 = buffer.data(ih + 560);
    const auto *ih_561 = buffer.data(ih + 561);
    const auto *ih_562 = buffer.data(ih + 562);
    const auto *ih_563 = buffer.data(ih + 563);
    const auto *ih_564 = buffer.data(ih + 564);
    const auto *ih_565 = buffer.data(ih + 565);
    const auto *ih_566 = buffer.data(ih + 566);
    const auto *ih_567 = buffer.data(ih + 567);
    const auto *ih_568 = buffer.data(ih + 568);
    const auto *ih_569 = buffer.data(ih + 569);
    const auto *ih_570 = buffer.data(ih + 570);
    const auto *ih_571 = buffer.data(ih + 571);
    const auto *ih_572 = buffer.data(ih + 572);
    const auto *ih_573 = buffer.data(ih + 573);
    const auto *ih_574 = buffer.data(ih + 574);
    const auto *ih_575 = buffer.data(ih + 575);
    const auto *ih_576 = buffer.data(ih + 576);
    const auto *ih_577 = buffer.data(ih + 577);
    const auto *ih_578 = buffer.data(ih + 578);
    const auto *ih_579 = buffer.data(ih + 579);
    const auto *ih_580 = buffer.data(ih + 580);
    const auto *ih_581 = buffer.data(ih + 581);
    const auto *ih_582 = buffer.data(ih + 582);
    const auto *ih_583 = buffer.data(ih + 583);
    const auto *ih_584 = buffer.data(ih + 584);
    const auto *ih_585 = buffer.data(ih + 585);
    const auto *ih_586 = buffer.data(ih + 586);
    const auto *ih_587 = buffer.data(ih + 587);

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, gh_212, gh_213, gh_214, gh_215, \
                         gh_216, ih_485, ih_486, ih_487, ih_488, \
                         ih_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = -gh_212[k]
                   + f_0 * ih_485[k];

        t_339[k] = -gh_213[k]
                   + f_0 * ih_486[k];

        t_340[k] = -gh_214[k]
                   + f_0 * ih_487[k];

        t_341[k] = -gh_215[k]
                   + f_0 * ih_488[k];

        t_342[k] = -gh_216[k]
                   + f_0 * ih_489[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, gh_217, gh_218, gh_219, gh_220, \
                         gh_221, ih_490, ih_491, ih_492, ih_493, \
                         ih_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = -gh_217[k]
                   + f_0 * ih_490[k];

        t_344[k] = -gh_218[k]
                   + f_0 * ih_491[k];

        t_345[k] = -gh_219[k]
                   + f_0 * ih_492[k];

        t_346[k] = -gh_220[k]
                   + f_0 * ih_493[k];

        t_347[k] = -gh_221[k]
                   + f_0 * ih_494[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, gh_222, gh_223, gh_224, gh_225, \
                         gh_226, ih_495, ih_496, ih_497, ih_498, \
                         ih_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = -gh_222[k]
                   + f_0 * ih_495[k];

        t_349[k] = -gh_223[k]
                   + f_0 * ih_496[k];

        t_350[k] = -gh_224[k]
                   + f_0 * ih_497[k];

        t_351[k] = -gh_225[k]
                   + f_0 * ih_498[k];

        t_352[k] = -gh_226[k]
                   + f_0 * ih_499[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, gh_227, gh_228, gh_229, gh_230, \
                         gh_231, ih_500, ih_501, ih_502, ih_503, \
                         ih_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = -gh_227[k]
                   + f_0 * ih_500[k];

        t_354[k] = -gh_228[k]
                   + f_0 * ih_501[k];

        t_355[k] = -gh_229[k]
                   + f_0 * ih_502[k];

        t_356[k] = -gh_230[k]
                   + f_0 * ih_503[k];

        t_357[k] = -2.0 * gh_231[k]
                   + f_0 * ih_504[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, gh_232, gh_233, gh_234, gh_235, \
                         gh_236, ih_505, ih_506, ih_507, ih_508, \
                         ih_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = -2.0 * gh_232[k]
                   + f_0 * ih_505[k];

        t_359[k] = -2.0 * gh_233[k]
                   + f_0 * ih_506[k];

        t_360[k] = -2.0 * gh_234[k]
                   + f_0 * ih_507[k];

        t_361[k] = -2.0 * gh_235[k]
                   + f_0 * ih_508[k];

        t_362[k] = -2.0 * gh_236[k]
                   + f_0 * ih_509[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, gh_237, gh_238, gh_239, gh_240, \
                         gh_241, ih_510, ih_511, ih_512, ih_513, \
                         ih_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = -2.0 * gh_237[k]
                   + f_0 * ih_510[k];

        t_364[k] = -2.0 * gh_238[k]
                   + f_0 * ih_511[k];

        t_365[k] = -2.0 * gh_239[k]
                   + f_0 * ih_512[k];

        t_366[k] = -2.0 * gh_240[k]
                   + f_0 * ih_513[k];

        t_367[k] = -2.0 * gh_241[k]
                   + f_0 * ih_514[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, gh_242, gh_243, gh_244, gh_245, \
                         gh_246, ih_515, ih_516, ih_517, ih_518, \
                         ih_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = -2.0 * gh_242[k]
                   + f_0 * ih_515[k];

        t_369[k] = -2.0 * gh_243[k]
                   + f_0 * ih_516[k];

        t_370[k] = -2.0 * gh_244[k]
                   + f_0 * ih_517[k];

        t_371[k] = -2.0 * gh_245[k]
                   + f_0 * ih_518[k];

        t_372[k] = -2.0 * gh_246[k]
                   + f_0 * ih_519[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, gh_247, gh_248, gh_249, gh_250, \
                         gh_251, ih_520, ih_521, ih_522, ih_523, \
                         ih_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = -2.0 * gh_247[k]
                   + f_0 * ih_520[k];

        t_374[k] = -2.0 * gh_248[k]
                   + f_0 * ih_521[k];

        t_375[k] = -2.0 * gh_249[k]
                   + f_0 * ih_522[k];

        t_376[k] = -2.0 * gh_250[k]
                   + f_0 * ih_523[k];

        t_377[k] = -2.0 * gh_251[k]
                   + f_0 * ih_524[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, gh_252, gh_253, gh_254, gh_255, \
                         gh_256, ih_525, ih_526, ih_527, ih_528, \
                         ih_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = -3.0 * gh_252[k]
                   + f_0 * ih_525[k];

        t_379[k] = -3.0 * gh_253[k]
                   + f_0 * ih_526[k];

        t_380[k] = -3.0 * gh_254[k]
                   + f_0 * ih_527[k];

        t_381[k] = -3.0 * gh_255[k]
                   + f_0 * ih_528[k];

        t_382[k] = -3.0 * gh_256[k]
                   + f_0 * ih_529[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, gh_257, gh_258, gh_259, gh_260, \
                         gh_261, ih_530, ih_531, ih_532, ih_533, \
                         ih_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = -3.0 * gh_257[k]
                   + f_0 * ih_530[k];

        t_384[k] = -3.0 * gh_258[k]
                   + f_0 * ih_531[k];

        t_385[k] = -3.0 * gh_259[k]
                   + f_0 * ih_532[k];

        t_386[k] = -3.0 * gh_260[k]
                   + f_0 * ih_533[k];

        t_387[k] = -3.0 * gh_261[k]
                   + f_0 * ih_534[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, gh_262, gh_263, gh_264, gh_265, \
                         gh_266, ih_535, ih_536, ih_537, ih_538, \
                         ih_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = -3.0 * gh_262[k]
                   + f_0 * ih_535[k];

        t_389[k] = -3.0 * gh_263[k]
                   + f_0 * ih_536[k];

        t_390[k] = -3.0 * gh_264[k]
                   + f_0 * ih_537[k];

        t_391[k] = -3.0 * gh_265[k]
                   + f_0 * ih_538[k];

        t_392[k] = -3.0 * gh_266[k]
                   + f_0 * ih_539[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, gh_267, gh_268, gh_269, gh_270, \
                         gh_271, ih_540, ih_541, ih_542, ih_543, \
                         ih_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = -3.0 * gh_267[k]
                   + f_0 * ih_540[k];

        t_394[k] = -3.0 * gh_268[k]
                   + f_0 * ih_541[k];

        t_395[k] = -3.0 * gh_269[k]
                   + f_0 * ih_542[k];

        t_396[k] = -3.0 * gh_270[k]
                   + f_0 * ih_543[k];

        t_397[k] = -3.0 * gh_271[k]
                   + f_0 * ih_544[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, t_402, gh_272, gh_273, gh_274, gh_275, \
                         gh_276, ih_545, ih_546, ih_547, ih_548, \
                         ih_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = -3.0 * gh_272[k]
                   + f_0 * ih_545[k];

        t_399[k] = -4.0 * gh_273[k]
                   + f_0 * ih_546[k];

        t_400[k] = -4.0 * gh_274[k]
                   + f_0 * ih_547[k];

        t_401[k] = -4.0 * gh_275[k]
                   + f_0 * ih_548[k];

        t_402[k] = -4.0 * gh_276[k]
                   + f_0 * ih_549[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, t_407, gh_277, gh_278, gh_279, gh_280, \
                         gh_281, ih_550, ih_551, ih_552, ih_553, \
                         ih_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = -4.0 * gh_277[k]
                   + f_0 * ih_550[k];

        t_404[k] = -4.0 * gh_278[k]
                   + f_0 * ih_551[k];

        t_405[k] = -4.0 * gh_279[k]
                   + f_0 * ih_552[k];

        t_406[k] = -4.0 * gh_280[k]
                   + f_0 * ih_553[k];

        t_407[k] = -4.0 * gh_281[k]
                   + f_0 * ih_554[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, gh_282, gh_283, gh_284, gh_285, \
                         gh_286, ih_555, ih_556, ih_557, ih_558, \
                         ih_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = -4.0 * gh_282[k]
                   + f_0 * ih_555[k];

        t_409[k] = -4.0 * gh_283[k]
                   + f_0 * ih_556[k];

        t_410[k] = -4.0 * gh_284[k]
                   + f_0 * ih_557[k];

        t_411[k] = -4.0 * gh_285[k]
                   + f_0 * ih_558[k];

        t_412[k] = -4.0 * gh_286[k]
                   + f_0 * ih_559[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, gh_287, gh_288, gh_289, gh_290, \
                         gh_291, ih_560, ih_561, ih_562, ih_563, \
                         ih_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = -4.0 * gh_287[k]
                   + f_0 * ih_560[k];

        t_414[k] = -4.0 * gh_288[k]
                   + f_0 * ih_561[k];

        t_415[k] = -4.0 * gh_289[k]
                   + f_0 * ih_562[k];

        t_416[k] = -4.0 * gh_290[k]
                   + f_0 * ih_563[k];

        t_417[k] = -4.0 * gh_291[k]
                   + f_0 * ih_564[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, gh_292, gh_293, gh_294, gh_295, \
                         gh_296, ih_565, ih_566, ih_567, ih_568, \
                         ih_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = -4.0 * gh_292[k]
                   + f_0 * ih_565[k];

        t_419[k] = -4.0 * gh_293[k]
                   + f_0 * ih_566[k];

        t_420[k] = -5.0 * gh_294[k]
                   + f_0 * ih_567[k];

        t_421[k] = -5.0 * gh_295[k]
                   + f_0 * ih_568[k];

        t_422[k] = -5.0 * gh_296[k]
                   + f_0 * ih_569[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, gh_297, gh_298, gh_299, gh_300, \
                         gh_301, ih_570, ih_571, ih_572, ih_573, \
                         ih_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = -5.0 * gh_297[k]
                   + f_0 * ih_570[k];

        t_424[k] = -5.0 * gh_298[k]
                   + f_0 * ih_571[k];

        t_425[k] = -5.0 * gh_299[k]
                   + f_0 * ih_572[k];

        t_426[k] = -5.0 * gh_300[k]
                   + f_0 * ih_573[k];

        t_427[k] = -5.0 * gh_301[k]
                   + f_0 * ih_574[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, gh_302, gh_303, gh_304, gh_305, \
                         gh_306, ih_575, ih_576, ih_577, ih_578, \
                         ih_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = -5.0 * gh_302[k]
                   + f_0 * ih_575[k];

        t_429[k] = -5.0 * gh_303[k]
                   + f_0 * ih_576[k];

        t_430[k] = -5.0 * gh_304[k]
                   + f_0 * ih_577[k];

        t_431[k] = -5.0 * gh_305[k]
                   + f_0 * ih_578[k];

        t_432[k] = -5.0 * gh_306[k]
                   + f_0 * ih_579[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, gh_307, gh_308, gh_309, gh_310, \
                         gh_311, ih_580, ih_581, ih_582, ih_583, \
                         ih_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = -5.0 * gh_307[k]
                   + f_0 * ih_580[k];

        t_434[k] = -5.0 * gh_308[k]
                   + f_0 * ih_581[k];

        t_435[k] = -5.0 * gh_309[k]
                   + f_0 * ih_582[k];

        t_436[k] = -5.0 * gh_310[k]
                   + f_0 * ih_583[k];

        t_437[k] = -5.0 * gh_311[k]
                   + f_0 * ih_584[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, gh_312, gh_313, gh_314, ih_585, ih_586, \
                         ih_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = -5.0 * gh_312[k]
                   + f_0 * ih_585[k];

        t_439[k] = -5.0 * gh_313[k]
                   + f_0 * ih_586[k];

        t_440[k] = -5.0 * gh_314[k]
                   + f_0 * ih_587[k];
    }
}

auto
compute_prim_geom_10_hh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t gh, const size_t ih,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hh_electron_repulsion_2_piece0(buffer, target, gh, ih, ncols, alpha);

    compute_prim_geom_10_hh_electron_repulsion_2_piece1(buffer, target, gh, ih, ncols, alpha);

    compute_prim_geom_10_hh_electron_repulsion_2_piece2(buffer, target, gh, ih, ncols, alpha);
}

}  // namespace simdt2ceri
