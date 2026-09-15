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


#include "SimdElectronRepulsionGeom10VrrRecGH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_gh_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fh, const size_t hh,
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

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_79 = buffer.data(fh + 79);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_104 = buffer.data(fh + 104);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_108 = buffer.data(fh + 108);
    const auto *fh_109 = buffer.data(fh + 109);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_111 = buffer.data(fh + 111);
    const auto *fh_112 = buffer.data(fh + 112);
    const auto *fh_113 = buffer.data(fh + 113);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_124 = buffer.data(fh + 124);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_128 = buffer.data(fh + 128);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_130 = buffer.data(fh + 130);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_133 = buffer.data(fh + 133);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_137 = buffer.data(fh + 137);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_147 = buffer.data(fh + 147);
    const auto *fh_148 = buffer.data(fh + 148);
    const auto *fh_149 = buffer.data(fh + 149);

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_1 = buffer.data(hh + 1);
    const auto *hh_2 = buffer.data(hh + 2);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_4 = buffer.data(hh + 4);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_6 = buffer.data(hh + 6);
    const auto *hh_7 = buffer.data(hh + 7);
    const auto *hh_8 = buffer.data(hh + 8);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_10 = buffer.data(hh + 10);
    const auto *hh_11 = buffer.data(hh + 11);
    const auto *hh_12 = buffer.data(hh + 12);
    const auto *hh_13 = buffer.data(hh + 13);
    const auto *hh_14 = buffer.data(hh + 14);
    const auto *hh_15 = buffer.data(hh + 15);
    const auto *hh_16 = buffer.data(hh + 16);
    const auto *hh_17 = buffer.data(hh + 17);
    const auto *hh_18 = buffer.data(hh + 18);
    const auto *hh_19 = buffer.data(hh + 19);
    const auto *hh_20 = buffer.data(hh + 20);
    const auto *hh_21 = buffer.data(hh + 21);
    const auto *hh_22 = buffer.data(hh + 22);
    const auto *hh_23 = buffer.data(hh + 23);
    const auto *hh_24 = buffer.data(hh + 24);
    const auto *hh_25 = buffer.data(hh + 25);
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_29 = buffer.data(hh + 29);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_33 = buffer.data(hh + 33);
    const auto *hh_34 = buffer.data(hh + 34);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_43 = buffer.data(hh + 43);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_45 = buffer.data(hh + 45);
    const auto *hh_46 = buffer.data(hh + 46);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_48 = buffer.data(hh + 48);
    const auto *hh_49 = buffer.data(hh + 49);
    const auto *hh_50 = buffer.data(hh + 50);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_52 = buffer.data(hh + 52);
    const auto *hh_53 = buffer.data(hh + 53);
    const auto *hh_54 = buffer.data(hh + 54);
    const auto *hh_55 = buffer.data(hh + 55);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_57 = buffer.data(hh + 57);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_64 = buffer.data(hh + 64);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_67 = buffer.data(hh + 67);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_71 = buffer.data(hh + 71);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_74 = buffer.data(hh + 74);
    const auto *hh_75 = buffer.data(hh + 75);
    const auto *hh_76 = buffer.data(hh + 76);
    const auto *hh_77 = buffer.data(hh + 77);
    const auto *hh_78 = buffer.data(hh + 78);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_80 = buffer.data(hh + 80);
    const auto *hh_81 = buffer.data(hh + 81);
    const auto *hh_82 = buffer.data(hh + 82);
    const auto *hh_83 = buffer.data(hh + 83);
    const auto *hh_84 = buffer.data(hh + 84);
    const auto *hh_85 = buffer.data(hh + 85);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_88 = buffer.data(hh + 88);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_91 = buffer.data(hh + 91);
    const auto *hh_92 = buffer.data(hh + 92);
    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_94 = buffer.data(hh + 94);
    const auto *hh_95 = buffer.data(hh + 95);
    const auto *hh_96 = buffer.data(hh + 96);
    const auto *hh_97 = buffer.data(hh + 97);
    const auto *hh_98 = buffer.data(hh + 98);
    const auto *hh_99 = buffer.data(hh + 99);
    const auto *hh_100 = buffer.data(hh + 100);
    const auto *hh_101 = buffer.data(hh + 101);
    const auto *hh_102 = buffer.data(hh + 102);
    const auto *hh_103 = buffer.data(hh + 103);
    const auto *hh_104 = buffer.data(hh + 104);
    const auto *hh_105 = buffer.data(hh + 105);
    const auto *hh_106 = buffer.data(hh + 106);
    const auto *hh_107 = buffer.data(hh + 107);
    const auto *hh_108 = buffer.data(hh + 108);
    const auto *hh_109 = buffer.data(hh + 109);
    const auto *hh_110 = buffer.data(hh + 110);
    const auto *hh_111 = buffer.data(hh + 111);
    const auto *hh_112 = buffer.data(hh + 112);
    const auto *hh_113 = buffer.data(hh + 113);
    const auto *hh_114 = buffer.data(hh + 114);
    const auto *hh_115 = buffer.data(hh + 115);
    const auto *hh_116 = buffer.data(hh + 116);
    const auto *hh_117 = buffer.data(hh + 117);
    const auto *hh_118 = buffer.data(hh + 118);
    const auto *hh_119 = buffer.data(hh + 119);
    const auto *hh_120 = buffer.data(hh + 120);
    const auto *hh_121 = buffer.data(hh + 121);
    const auto *hh_122 = buffer.data(hh + 122);
    const auto *hh_123 = buffer.data(hh + 123);
    const auto *hh_124 = buffer.data(hh + 124);
    const auto *hh_125 = buffer.data(hh + 125);
    const auto *hh_126 = buffer.data(hh + 126);
    const auto *hh_127 = buffer.data(hh + 127);
    const auto *hh_128 = buffer.data(hh + 128);
    const auto *hh_129 = buffer.data(hh + 129);
    const auto *hh_130 = buffer.data(hh + 130);
    const auto *hh_131 = buffer.data(hh + 131);
    const auto *hh_132 = buffer.data(hh + 132);
    const auto *hh_133 = buffer.data(hh + 133);
    const auto *hh_134 = buffer.data(hh + 134);
    const auto *hh_135 = buffer.data(hh + 135);
    const auto *hh_136 = buffer.data(hh + 136);
    const auto *hh_137 = buffer.data(hh + 137);
    const auto *hh_138 = buffer.data(hh + 138);
    const auto *hh_139 = buffer.data(hh + 139);
    const auto *hh_140 = buffer.data(hh + 140);
    const auto *hh_141 = buffer.data(hh + 141);
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_146 = buffer.data(hh + 146);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_148 = buffer.data(hh + 148);
    const auto *hh_149 = buffer.data(hh + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, fh_0, fh_1, fh_2, fh_3, fh_4, hh_0, hh_1, \
                         hh_2, hh_3, hh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -4.0 * fh_0[k]
                 + f_0 * hh_0[k];

        t_1[k] = -4.0 * fh_1[k]
                 + f_0 * hh_1[k];

        t_2[k] = -4.0 * fh_2[k]
                 + f_0 * hh_2[k];

        t_3[k] = -4.0 * fh_3[k]
                 + f_0 * hh_3[k];

        t_4[k] = -4.0 * fh_4[k]
                 + f_0 * hh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, fh_5, fh_6, fh_7, fh_8, fh_9, hh_5, hh_6, \
                         hh_7, hh_8, hh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -4.0 * fh_5[k]
                 + f_0 * hh_5[k];

        t_6[k] = -4.0 * fh_6[k]
                 + f_0 * hh_6[k];

        t_7[k] = -4.0 * fh_7[k]
                 + f_0 * hh_7[k];

        t_8[k] = -4.0 * fh_8[k]
                 + f_0 * hh_8[k];

        t_9[k] = -4.0 * fh_9[k]
                 + f_0 * hh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, fh_10, fh_11, fh_12, fh_13, fh_14, \
                         hh_10, hh_11, hh_12, hh_13, hh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -4.0 * fh_10[k]
                  + f_0 * hh_10[k];

        t_11[k] = -4.0 * fh_11[k]
                  + f_0 * hh_11[k];

        t_12[k] = -4.0 * fh_12[k]
                  + f_0 * hh_12[k];

        t_13[k] = -4.0 * fh_13[k]
                  + f_0 * hh_13[k];

        t_14[k] = -4.0 * fh_14[k]
                  + f_0 * hh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, fh_15, fh_16, fh_17, fh_18, fh_19, \
                         hh_15, hh_16, hh_17, hh_18, hh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -4.0 * fh_15[k]
                  + f_0 * hh_15[k];

        t_16[k] = -4.0 * fh_16[k]
                  + f_0 * hh_16[k];

        t_17[k] = -4.0 * fh_17[k]
                  + f_0 * hh_17[k];

        t_18[k] = -4.0 * fh_18[k]
                  + f_0 * hh_18[k];

        t_19[k] = -4.0 * fh_19[k]
                  + f_0 * hh_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, fh_20, fh_21, fh_22, fh_23, fh_24, \
                         hh_20, hh_21, hh_22, hh_23, hh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -4.0 * fh_20[k]
                  + f_0 * hh_20[k];

        t_21[k] = -3.0 * fh_21[k]
                  + f_0 * hh_21[k];

        t_22[k] = -3.0 * fh_22[k]
                  + f_0 * hh_22[k];

        t_23[k] = -3.0 * fh_23[k]
                  + f_0 * hh_23[k];

        t_24[k] = -3.0 * fh_24[k]
                  + f_0 * hh_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, fh_25, fh_26, fh_27, fh_28, fh_29, \
                         hh_25, hh_26, hh_27, hh_28, hh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -3.0 * fh_25[k]
                  + f_0 * hh_25[k];

        t_26[k] = -3.0 * fh_26[k]
                  + f_0 * hh_26[k];

        t_27[k] = -3.0 * fh_27[k]
                  + f_0 * hh_27[k];

        t_28[k] = -3.0 * fh_28[k]
                  + f_0 * hh_28[k];

        t_29[k] = -3.0 * fh_29[k]
                  + f_0 * hh_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, fh_30, fh_31, fh_32, fh_33, fh_34, \
                         hh_30, hh_31, hh_32, hh_33, hh_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -3.0 * fh_30[k]
                  + f_0 * hh_30[k];

        t_31[k] = -3.0 * fh_31[k]
                  + f_0 * hh_31[k];

        t_32[k] = -3.0 * fh_32[k]
                  + f_0 * hh_32[k];

        t_33[k] = -3.0 * fh_33[k]
                  + f_0 * hh_33[k];

        t_34[k] = -3.0 * fh_34[k]
                  + f_0 * hh_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, fh_35, fh_36, fh_37, fh_38, fh_39, \
                         hh_35, hh_36, hh_37, hh_38, hh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -3.0 * fh_35[k]
                  + f_0 * hh_35[k];

        t_36[k] = -3.0 * fh_36[k]
                  + f_0 * hh_36[k];

        t_37[k] = -3.0 * fh_37[k]
                  + f_0 * hh_37[k];

        t_38[k] = -3.0 * fh_38[k]
                  + f_0 * hh_38[k];

        t_39[k] = -3.0 * fh_39[k]
                  + f_0 * hh_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, fh_40, fh_41, fh_42, fh_43, fh_44, \
                         hh_40, hh_41, hh_42, hh_43, hh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * fh_40[k]
                  + f_0 * hh_40[k];

        t_41[k] = -3.0 * fh_41[k]
                  + f_0 * hh_41[k];

        t_42[k] = -3.0 * fh_42[k]
                  + f_0 * hh_42[k];

        t_43[k] = -3.0 * fh_43[k]
                  + f_0 * hh_43[k];

        t_44[k] = -3.0 * fh_44[k]
                  + f_0 * hh_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, fh_45, fh_46, fh_47, fh_48, fh_49, \
                         hh_45, hh_46, hh_47, hh_48, hh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -3.0 * fh_45[k]
                  + f_0 * hh_45[k];

        t_46[k] = -3.0 * fh_46[k]
                  + f_0 * hh_46[k];

        t_47[k] = -3.0 * fh_47[k]
                  + f_0 * hh_47[k];

        t_48[k] = -3.0 * fh_48[k]
                  + f_0 * hh_48[k];

        t_49[k] = -3.0 * fh_49[k]
                  + f_0 * hh_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, fh_50, fh_51, fh_52, fh_53, fh_54, \
                         hh_50, hh_51, hh_52, hh_53, hh_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -3.0 * fh_50[k]
                  + f_0 * hh_50[k];

        t_51[k] = -3.0 * fh_51[k]
                  + f_0 * hh_51[k];

        t_52[k] = -3.0 * fh_52[k]
                  + f_0 * hh_52[k];

        t_53[k] = -3.0 * fh_53[k]
                  + f_0 * hh_53[k];

        t_54[k] = -3.0 * fh_54[k]
                  + f_0 * hh_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, fh_55, fh_56, fh_57, fh_58, fh_59, \
                         hh_55, hh_56, hh_57, hh_58, hh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -3.0 * fh_55[k]
                  + f_0 * hh_55[k];

        t_56[k] = -3.0 * fh_56[k]
                  + f_0 * hh_56[k];

        t_57[k] = -3.0 * fh_57[k]
                  + f_0 * hh_57[k];

        t_58[k] = -3.0 * fh_58[k]
                  + f_0 * hh_58[k];

        t_59[k] = -3.0 * fh_59[k]
                  + f_0 * hh_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, fh_60, fh_61, fh_62, fh_63, fh_64, \
                         hh_60, hh_61, hh_62, hh_63, hh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * fh_60[k]
                  + f_0 * hh_60[k];

        t_61[k] = -3.0 * fh_61[k]
                  + f_0 * hh_61[k];

        t_62[k] = -3.0 * fh_62[k]
                  + f_0 * hh_62[k];

        t_63[k] = -2.0 * fh_63[k]
                  + f_0 * hh_63[k];

        t_64[k] = -2.0 * fh_64[k]
                  + f_0 * hh_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, fh_65, fh_66, fh_67, fh_68, fh_69, \
                         hh_65, hh_66, hh_67, hh_68, hh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -2.0 * fh_65[k]
                  + f_0 * hh_65[k];

        t_66[k] = -2.0 * fh_66[k]
                  + f_0 * hh_66[k];

        t_67[k] = -2.0 * fh_67[k]
                  + f_0 * hh_67[k];

        t_68[k] = -2.0 * fh_68[k]
                  + f_0 * hh_68[k];

        t_69[k] = -2.0 * fh_69[k]
                  + f_0 * hh_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, fh_70, fh_71, fh_72, fh_73, fh_74, \
                         hh_70, hh_71, hh_72, hh_73, hh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * fh_70[k]
                  + f_0 * hh_70[k];

        t_71[k] = -2.0 * fh_71[k]
                  + f_0 * hh_71[k];

        t_72[k] = -2.0 * fh_72[k]
                  + f_0 * hh_72[k];

        t_73[k] = -2.0 * fh_73[k]
                  + f_0 * hh_73[k];

        t_74[k] = -2.0 * fh_74[k]
                  + f_0 * hh_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, fh_75, fh_76, fh_77, fh_78, fh_79, \
                         hh_75, hh_76, hh_77, hh_78, hh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * fh_75[k]
                  + f_0 * hh_75[k];

        t_76[k] = -2.0 * fh_76[k]
                  + f_0 * hh_76[k];

        t_77[k] = -2.0 * fh_77[k]
                  + f_0 * hh_77[k];

        t_78[k] = -2.0 * fh_78[k]
                  + f_0 * hh_78[k];

        t_79[k] = -2.0 * fh_79[k]
                  + f_0 * hh_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, fh_80, fh_81, fh_82, fh_83, fh_84, \
                         hh_80, hh_81, hh_82, hh_83, hh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -2.0 * fh_80[k]
                  + f_0 * hh_80[k];

        t_81[k] = -2.0 * fh_81[k]
                  + f_0 * hh_81[k];

        t_82[k] = -2.0 * fh_82[k]
                  + f_0 * hh_82[k];

        t_83[k] = -2.0 * fh_83[k]
                  + f_0 * hh_83[k];

        t_84[k] = -2.0 * fh_84[k]
                  + f_0 * hh_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, fh_85, fh_86, fh_87, fh_88, fh_89, \
                         hh_85, hh_86, hh_87, hh_88, hh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * fh_85[k]
                  + f_0 * hh_85[k];

        t_86[k] = -2.0 * fh_86[k]
                  + f_0 * hh_86[k];

        t_87[k] = -2.0 * fh_87[k]
                  + f_0 * hh_87[k];

        t_88[k] = -2.0 * fh_88[k]
                  + f_0 * hh_88[k];

        t_89[k] = -2.0 * fh_89[k]
                  + f_0 * hh_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, fh_90, fh_91, fh_92, fh_93, fh_94, \
                         hh_90, hh_91, hh_92, hh_93, hh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * fh_90[k]
                  + f_0 * hh_90[k];

        t_91[k] = -2.0 * fh_91[k]
                  + f_0 * hh_91[k];

        t_92[k] = -2.0 * fh_92[k]
                  + f_0 * hh_92[k];

        t_93[k] = -2.0 * fh_93[k]
                  + f_0 * hh_93[k];

        t_94[k] = -2.0 * fh_94[k]
                  + f_0 * hh_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, fh_95, fh_96, fh_97, fh_98, fh_99, \
                         hh_95, hh_96, hh_97, hh_98, hh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * fh_95[k]
                  + f_0 * hh_95[k];

        t_96[k] = -2.0 * fh_96[k]
                  + f_0 * hh_96[k];

        t_97[k] = -2.0 * fh_97[k]
                  + f_0 * hh_97[k];

        t_98[k] = -2.0 * fh_98[k]
                  + f_0 * hh_98[k];

        t_99[k] = -2.0 * fh_99[k]
                  + f_0 * hh_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, fh_100, fh_101, fh_102, fh_103, \
                         fh_104, hh_100, hh_101, hh_102, hh_103, \
                         hh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * fh_100[k]
                   + f_0 * hh_100[k];

        t_101[k] = -2.0 * fh_101[k]
                   + f_0 * hh_101[k];

        t_102[k] = -2.0 * fh_102[k]
                   + f_0 * hh_102[k];

        t_103[k] = -2.0 * fh_103[k]
                   + f_0 * hh_103[k];

        t_104[k] = -2.0 * fh_104[k]
                   + f_0 * hh_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, fh_105, fh_106, fh_107, fh_108, \
                         fh_109, hh_105, hh_106, hh_107, hh_108, \
                         hh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * fh_105[k]
                   + f_0 * hh_105[k];

        t_106[k] = -2.0 * fh_106[k]
                   + f_0 * hh_106[k];

        t_107[k] = -2.0 * fh_107[k]
                   + f_0 * hh_107[k];

        t_108[k] = -2.0 * fh_108[k]
                   + f_0 * hh_108[k];

        t_109[k] = -2.0 * fh_109[k]
                   + f_0 * hh_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, fh_110, fh_111, fh_112, fh_113, \
                         fh_114, hh_110, hh_111, hh_112, hh_113, \
                         hh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * fh_110[k]
                   + f_0 * hh_110[k];

        t_111[k] = -2.0 * fh_111[k]
                   + f_0 * hh_111[k];

        t_112[k] = -2.0 * fh_112[k]
                   + f_0 * hh_112[k];

        t_113[k] = -2.0 * fh_113[k]
                   + f_0 * hh_113[k];

        t_114[k] = -2.0 * fh_114[k]
                   + f_0 * hh_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, fh_115, fh_116, fh_117, fh_118, \
                         fh_119, hh_115, hh_116, hh_117, hh_118, \
                         hh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * fh_115[k]
                   + f_0 * hh_115[k];

        t_116[k] = -2.0 * fh_116[k]
                   + f_0 * hh_116[k];

        t_117[k] = -2.0 * fh_117[k]
                   + f_0 * hh_117[k];

        t_118[k] = -2.0 * fh_118[k]
                   + f_0 * hh_118[k];

        t_119[k] = -2.0 * fh_119[k]
                   + f_0 * hh_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, fh_120, fh_121, fh_122, fh_123, \
                         fh_124, hh_120, hh_121, hh_122, hh_123, \
                         hh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -2.0 * fh_120[k]
                   + f_0 * hh_120[k];

        t_121[k] = -2.0 * fh_121[k]
                   + f_0 * hh_121[k];

        t_122[k] = -2.0 * fh_122[k]
                   + f_0 * hh_122[k];

        t_123[k] = -2.0 * fh_123[k]
                   + f_0 * hh_123[k];

        t_124[k] = -2.0 * fh_124[k]
                   + f_0 * hh_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, fh_125, fh_126, fh_127, fh_128, \
                         fh_129, hh_125, hh_126, hh_127, hh_128, \
                         hh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -2.0 * fh_125[k]
                   + f_0 * hh_125[k];

        t_126[k] = -fh_126[k]
                   + f_0 * hh_126[k];

        t_127[k] = -fh_127[k]
                   + f_0 * hh_127[k];

        t_128[k] = -fh_128[k]
                   + f_0 * hh_128[k];

        t_129[k] = -fh_129[k]
                   + f_0 * hh_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, fh_130, fh_131, fh_132, fh_133, \
                         fh_134, hh_130, hh_131, hh_132, hh_133, \
                         hh_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -fh_130[k]
                   + f_0 * hh_130[k];

        t_131[k] = -fh_131[k]
                   + f_0 * hh_131[k];

        t_132[k] = -fh_132[k]
                   + f_0 * hh_132[k];

        t_133[k] = -fh_133[k]
                   + f_0 * hh_133[k];

        t_134[k] = -fh_134[k]
                   + f_0 * hh_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, fh_135, fh_136, fh_137, fh_138, \
                         fh_139, hh_135, hh_136, hh_137, hh_138, \
                         hh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -fh_135[k]
                   + f_0 * hh_135[k];

        t_136[k] = -fh_136[k]
                   + f_0 * hh_136[k];

        t_137[k] = -fh_137[k]
                   + f_0 * hh_137[k];

        t_138[k] = -fh_138[k]
                   + f_0 * hh_138[k];

        t_139[k] = -fh_139[k]
                   + f_0 * hh_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, fh_140, fh_141, fh_142, fh_143, \
                         fh_144, hh_140, hh_141, hh_142, hh_143, \
                         hh_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -fh_140[k]
                   + f_0 * hh_140[k];

        t_141[k] = -fh_141[k]
                   + f_0 * hh_141[k];

        t_142[k] = -fh_142[k]
                   + f_0 * hh_142[k];

        t_143[k] = -fh_143[k]
                   + f_0 * hh_143[k];

        t_144[k] = -fh_144[k]
                   + f_0 * hh_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, fh_145, fh_146, fh_147, fh_148, \
                         fh_149, hh_145, hh_146, hh_147, hh_148, \
                         hh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -fh_145[k]
                   + f_0 * hh_145[k];

        t_146[k] = -fh_146[k]
                   + f_0 * hh_146[k];

        t_147[k] = -fh_147[k]
                   + f_0 * hh_147[k];

        t_148[k] = -fh_148[k]
                   + f_0 * hh_148[k];

        t_149[k] = -fh_149[k]
                   + f_0 * hh_149[k];
    }
}

static auto
compute_prim_geom_10_gh_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fh, const size_t hh,
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
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);

    const auto *fh_150 = buffer.data(fh + 150);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_153 = buffer.data(fh + 153);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_157 = buffer.data(fh + 157);
    const auto *fh_158 = buffer.data(fh + 158);
    const auto *fh_159 = buffer.data(fh + 159);
    const auto *fh_160 = buffer.data(fh + 160);
    const auto *fh_161 = buffer.data(fh + 161);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_168 = buffer.data(fh + 168);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_170 = buffer.data(fh + 170);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_173 = buffer.data(fh + 173);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_177 = buffer.data(fh + 177);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_182 = buffer.data(fh + 182);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_190 = buffer.data(fh + 190);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_193 = buffer.data(fh + 193);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_197 = buffer.data(fh + 197);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_202 = buffer.data(fh + 202);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_151 = buffer.data(hh + 151);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_155 = buffer.data(hh + 155);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_158 = buffer.data(hh + 158);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_160 = buffer.data(hh + 160);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_169 = buffer.data(hh + 169);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_172 = buffer.data(hh + 172);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_178 = buffer.data(hh + 178);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_181 = buffer.data(hh + 181);
    const auto *hh_182 = buffer.data(hh + 182);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_189 = buffer.data(hh + 189);
    const auto *hh_190 = buffer.data(hh + 190);
    const auto *hh_191 = buffer.data(hh + 191);
    const auto *hh_192 = buffer.data(hh + 192);
    const auto *hh_193 = buffer.data(hh + 193);
    const auto *hh_194 = buffer.data(hh + 194);
    const auto *hh_195 = buffer.data(hh + 195);
    const auto *hh_196 = buffer.data(hh + 196);
    const auto *hh_197 = buffer.data(hh + 197);
    const auto *hh_198 = buffer.data(hh + 198);
    const auto *hh_199 = buffer.data(hh + 199);
    const auto *hh_200 = buffer.data(hh + 200);
    const auto *hh_201 = buffer.data(hh + 201);
    const auto *hh_202 = buffer.data(hh + 202);
    const auto *hh_203 = buffer.data(hh + 203);
    const auto *hh_204 = buffer.data(hh + 204);
    const auto *hh_205 = buffer.data(hh + 205);
    const auto *hh_206 = buffer.data(hh + 206);
    const auto *hh_207 = buffer.data(hh + 207);
    const auto *hh_208 = buffer.data(hh + 208);
    const auto *hh_209 = buffer.data(hh + 209);
    const auto *hh_210 = buffer.data(hh + 210);
    const auto *hh_211 = buffer.data(hh + 211);
    const auto *hh_212 = buffer.data(hh + 212);
    const auto *hh_213 = buffer.data(hh + 213);
    const auto *hh_214 = buffer.data(hh + 214);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_218 = buffer.data(hh + 218);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_221 = buffer.data(hh + 221);
    const auto *hh_222 = buffer.data(hh + 222);
    const auto *hh_223 = buffer.data(hh + 223);
    const auto *hh_224 = buffer.data(hh + 224);
    const auto *hh_225 = buffer.data(hh + 225);
    const auto *hh_226 = buffer.data(hh + 226);
    const auto *hh_227 = buffer.data(hh + 227);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_229 = buffer.data(hh + 229);
    const auto *hh_230 = buffer.data(hh + 230);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_232 = buffer.data(hh + 232);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_235 = buffer.data(hh + 235);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_238 = buffer.data(hh + 238);
    const auto *hh_239 = buffer.data(hh + 239);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_243 = buffer.data(hh + 243);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_246 = buffer.data(hh + 246);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_253 = buffer.data(hh + 253);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_256 = buffer.data(hh + 256);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_260 = buffer.data(hh + 260);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_262 = buffer.data(hh + 262);
    const auto *hh_263 = buffer.data(hh + 263);
    const auto *hh_264 = buffer.data(hh + 264);
    const auto *hh_265 = buffer.data(hh + 265);
    const auto *hh_266 = buffer.data(hh + 266);
    const auto *hh_267 = buffer.data(hh + 267);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_269 = buffer.data(hh + 269);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_271 = buffer.data(hh + 271);
    const auto *hh_272 = buffer.data(hh + 272);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_274 = buffer.data(hh + 274);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_277 = buffer.data(hh + 277);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_280 = buffer.data(hh + 280);
    const auto *hh_281 = buffer.data(hh + 281);
    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_283 = buffer.data(hh + 283);
    const auto *hh_284 = buffer.data(hh + 284);
    const auto *hh_285 = buffer.data(hh + 285);
    const auto *hh_286 = buffer.data(hh + 286);
    const auto *hh_287 = buffer.data(hh + 287);
    const auto *hh_288 = buffer.data(hh + 288);
    const auto *hh_289 = buffer.data(hh + 289);
    const auto *hh_290 = buffer.data(hh + 290);
    const auto *hh_291 = buffer.data(hh + 291);
    const auto *hh_292 = buffer.data(hh + 292);
    const auto *hh_293 = buffer.data(hh + 293);
    const auto *hh_294 = buffer.data(hh + 294);
    const auto *hh_295 = buffer.data(hh + 295);
    const auto *hh_296 = buffer.data(hh + 296);
    const auto *hh_297 = buffer.data(hh + 297);
    const auto *hh_298 = buffer.data(hh + 298);
    const auto *hh_299 = buffer.data(hh + 299);
    const auto *hh_300 = buffer.data(hh + 300);
    const auto *hh_301 = buffer.data(hh + 301);
    const auto *hh_302 = buffer.data(hh + 302);
    const auto *hh_303 = buffer.data(hh + 303);
    const auto *hh_304 = buffer.data(hh + 304);
    const auto *hh_305 = buffer.data(hh + 305);
    const auto *hh_306 = buffer.data(hh + 306);
    const auto *hh_307 = buffer.data(hh + 307);
    const auto *hh_308 = buffer.data(hh + 308);
    const auto *hh_309 = buffer.data(hh + 309);
    const auto *hh_310 = buffer.data(hh + 310);
    const auto *hh_311 = buffer.data(hh + 311);
    const auto *hh_312 = buffer.data(hh + 312);
    const auto *hh_313 = buffer.data(hh + 313);
    const auto *hh_314 = buffer.data(hh + 314);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, fh_150, fh_151, fh_152, fh_153, \
                         fh_154, hh_150, hh_151, hh_152, hh_153, \
                         hh_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -fh_150[k]
                   + f_0 * hh_150[k];

        t_151[k] = -fh_151[k]
                   + f_0 * hh_151[k];

        t_152[k] = -fh_152[k]
                   + f_0 * hh_152[k];

        t_153[k] = -fh_153[k]
                   + f_0 * hh_153[k];

        t_154[k] = -fh_154[k]
                   + f_0 * hh_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, fh_155, fh_156, fh_157, fh_158, \
                         fh_159, hh_155, hh_156, hh_157, hh_158, \
                         hh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -fh_155[k]
                   + f_0 * hh_155[k];

        t_156[k] = -fh_156[k]
                   + f_0 * hh_156[k];

        t_157[k] = -fh_157[k]
                   + f_0 * hh_157[k];

        t_158[k] = -fh_158[k]
                   + f_0 * hh_158[k];

        t_159[k] = -fh_159[k]
                   + f_0 * hh_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, fh_160, fh_161, fh_162, fh_163, \
                         fh_164, hh_160, hh_161, hh_162, hh_163, \
                         hh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -fh_160[k]
                   + f_0 * hh_160[k];

        t_161[k] = -fh_161[k]
                   + f_0 * hh_161[k];

        t_162[k] = -fh_162[k]
                   + f_0 * hh_162[k];

        t_163[k] = -fh_163[k]
                   + f_0 * hh_163[k];

        t_164[k] = -fh_164[k]
                   + f_0 * hh_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, fh_165, fh_166, fh_167, fh_168, \
                         fh_169, hh_165, hh_166, hh_167, hh_168, \
                         hh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -fh_165[k]
                   + f_0 * hh_165[k];

        t_166[k] = -fh_166[k]
                   + f_0 * hh_166[k];

        t_167[k] = -fh_167[k]
                   + f_0 * hh_167[k];

        t_168[k] = -fh_168[k]
                   + f_0 * hh_168[k];

        t_169[k] = -fh_169[k]
                   + f_0 * hh_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, fh_170, fh_171, fh_172, fh_173, \
                         fh_174, hh_170, hh_171, hh_172, hh_173, \
                         hh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -fh_170[k]
                   + f_0 * hh_170[k];

        t_171[k] = -fh_171[k]
                   + f_0 * hh_171[k];

        t_172[k] = -fh_172[k]
                   + f_0 * hh_172[k];

        t_173[k] = -fh_173[k]
                   + f_0 * hh_173[k];

        t_174[k] = -fh_174[k]
                   + f_0 * hh_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, fh_175, fh_176, fh_177, fh_178, \
                         fh_179, hh_175, hh_176, hh_177, hh_178, \
                         hh_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -fh_175[k]
                   + f_0 * hh_175[k];

        t_176[k] = -fh_176[k]
                   + f_0 * hh_176[k];

        t_177[k] = -fh_177[k]
                   + f_0 * hh_177[k];

        t_178[k] = -fh_178[k]
                   + f_0 * hh_178[k];

        t_179[k] = -fh_179[k]
                   + f_0 * hh_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, fh_180, fh_181, fh_182, fh_183, \
                         fh_184, hh_180, hh_181, hh_182, hh_183, \
                         hh_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -fh_180[k]
                   + f_0 * hh_180[k];

        t_181[k] = -fh_181[k]
                   + f_0 * hh_181[k];

        t_182[k] = -fh_182[k]
                   + f_0 * hh_182[k];

        t_183[k] = -fh_183[k]
                   + f_0 * hh_183[k];

        t_184[k] = -fh_184[k]
                   + f_0 * hh_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, fh_185, fh_186, fh_187, fh_188, \
                         fh_189, hh_185, hh_186, hh_187, hh_188, \
                         hh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -fh_185[k]
                   + f_0 * hh_185[k];

        t_186[k] = -fh_186[k]
                   + f_0 * hh_186[k];

        t_187[k] = -fh_187[k]
                   + f_0 * hh_187[k];

        t_188[k] = -fh_188[k]
                   + f_0 * hh_188[k];

        t_189[k] = -fh_189[k]
                   + f_0 * hh_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, fh_190, fh_191, fh_192, fh_193, \
                         fh_194, hh_190, hh_191, hh_192, hh_193, \
                         hh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -fh_190[k]
                   + f_0 * hh_190[k];

        t_191[k] = -fh_191[k]
                   + f_0 * hh_191[k];

        t_192[k] = -fh_192[k]
                   + f_0 * hh_192[k];

        t_193[k] = -fh_193[k]
                   + f_0 * hh_193[k];

        t_194[k] = -fh_194[k]
                   + f_0 * hh_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, fh_195, fh_196, fh_197, fh_198, \
                         fh_199, hh_195, hh_196, hh_197, hh_198, \
                         hh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -fh_195[k]
                   + f_0 * hh_195[k];

        t_196[k] = -fh_196[k]
                   + f_0 * hh_196[k];

        t_197[k] = -fh_197[k]
                   + f_0 * hh_197[k];

        t_198[k] = -fh_198[k]
                   + f_0 * hh_198[k];

        t_199[k] = -fh_199[k]
                   + f_0 * hh_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, fh_200, fh_201, fh_202, fh_203, \
                         fh_204, hh_200, hh_201, hh_202, hh_203, \
                         hh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -fh_200[k]
                   + f_0 * hh_200[k];

        t_201[k] = -fh_201[k]
                   + f_0 * hh_201[k];

        t_202[k] = -fh_202[k]
                   + f_0 * hh_202[k];

        t_203[k] = -fh_203[k]
                   + f_0 * hh_203[k];

        t_204[k] = -fh_204[k]
                   + f_0 * hh_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, fh_205, fh_206, fh_207, fh_208, \
                         fh_209, hh_205, hh_206, hh_207, hh_208, \
                         hh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -fh_205[k]
                   + f_0 * hh_205[k];

        t_206[k] = -fh_206[k]
                   + f_0 * hh_206[k];

        t_207[k] = -fh_207[k]
                   + f_0 * hh_207[k];

        t_208[k] = -fh_208[k]
                   + f_0 * hh_208[k];

        t_209[k] = -fh_209[k]
                   + f_0 * hh_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, hh_210, \
                         hh_211, hh_212, hh_213, hh_214, hh_215, hh_216, \
                         hh_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * hh_210[k];

        t_211[k] = f_0 * hh_211[k];

        t_212[k] = f_0 * hh_212[k];

        t_213[k] = f_0 * hh_213[k];

        t_214[k] = f_0 * hh_214[k];

        t_215[k] = f_0 * hh_215[k];

        t_216[k] = f_0 * hh_216[k];

        t_217[k] = f_0 * hh_217[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, t_224, t_225, hh_218, \
                         hh_219, hh_220, hh_221, hh_222, hh_223, hh_224, \
                         hh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * hh_218[k];

        t_219[k] = f_0 * hh_219[k];

        t_220[k] = f_0 * hh_220[k];

        t_221[k] = f_0 * hh_221[k];

        t_222[k] = f_0 * hh_222[k];

        t_223[k] = f_0 * hh_223[k];

        t_224[k] = f_0 * hh_224[k];

        t_225[k] = f_0 * hh_225[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, t_230, t_231, t_232, t_233, hh_226, \
                         hh_227, hh_228, hh_229, hh_230, hh_231, hh_232, \
                         hh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_0 * hh_226[k];

        t_227[k] = f_0 * hh_227[k];

        t_228[k] = f_0 * hh_228[k];

        t_229[k] = f_0 * hh_229[k];

        t_230[k] = f_0 * hh_230[k];

        t_231[k] = f_0 * hh_231[k];

        t_232[k] = f_0 * hh_232[k];

        t_233[k] = f_0 * hh_233[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, t_239, t_240, t_241, hh_234, \
                         hh_235, hh_236, hh_237, hh_238, hh_239, hh_240, \
                         hh_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_0 * hh_234[k];

        t_235[k] = f_0 * hh_235[k];

        t_236[k] = f_0 * hh_236[k];

        t_237[k] = f_0 * hh_237[k];

        t_238[k] = f_0 * hh_238[k];

        t_239[k] = f_0 * hh_239[k];

        t_240[k] = f_0 * hh_240[k];

        t_241[k] = f_0 * hh_241[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, t_247, t_248, t_249, hh_242, \
                         hh_243, hh_244, hh_245, hh_246, hh_247, hh_248, \
                         hh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_0 * hh_242[k];

        t_243[k] = f_0 * hh_243[k];

        t_244[k] = f_0 * hh_244[k];

        t_245[k] = f_0 * hh_245[k];

        t_246[k] = f_0 * hh_246[k];

        t_247[k] = f_0 * hh_247[k];

        t_248[k] = f_0 * hh_248[k];

        t_249[k] = f_0 * hh_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, t_255, t_256, t_257, hh_250, \
                         hh_251, hh_252, hh_253, hh_254, hh_255, hh_256, \
                         hh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_0 * hh_250[k];

        t_251[k] = f_0 * hh_251[k];

        t_252[k] = f_0 * hh_252[k];

        t_253[k] = f_0 * hh_253[k];

        t_254[k] = f_0 * hh_254[k];

        t_255[k] = f_0 * hh_255[k];

        t_256[k] = f_0 * hh_256[k];

        t_257[k] = f_0 * hh_257[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, t_263, t_264, t_265, hh_258, \
                         hh_259, hh_260, hh_261, hh_262, hh_263, hh_264, \
                         hh_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_0 * hh_258[k];

        t_259[k] = f_0 * hh_259[k];

        t_260[k] = f_0 * hh_260[k];

        t_261[k] = f_0 * hh_261[k];

        t_262[k] = f_0 * hh_262[k];

        t_263[k] = f_0 * hh_263[k];

        t_264[k] = f_0 * hh_264[k];

        t_265[k] = f_0 * hh_265[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, t_271, t_272, t_273, hh_266, \
                         hh_267, hh_268, hh_269, hh_270, hh_271, hh_272, \
                         hh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_0 * hh_266[k];

        t_267[k] = f_0 * hh_267[k];

        t_268[k] = f_0 * hh_268[k];

        t_269[k] = f_0 * hh_269[k];

        t_270[k] = f_0 * hh_270[k];

        t_271[k] = f_0 * hh_271[k];

        t_272[k] = f_0 * hh_272[k];

        t_273[k] = f_0 * hh_273[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, t_279, t_280, t_281, hh_274, \
                         hh_275, hh_276, hh_277, hh_278, hh_279, hh_280, \
                         hh_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_0 * hh_274[k];

        t_275[k] = f_0 * hh_275[k];

        t_276[k] = f_0 * hh_276[k];

        t_277[k] = f_0 * hh_277[k];

        t_278[k] = f_0 * hh_278[k];

        t_279[k] = f_0 * hh_279[k];

        t_280[k] = f_0 * hh_280[k];

        t_281[k] = f_0 * hh_281[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, t_287, t_288, t_289, hh_282, \
                         hh_283, hh_284, hh_285, hh_286, hh_287, hh_288, \
                         hh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_0 * hh_282[k];

        t_283[k] = f_0 * hh_283[k];

        t_284[k] = f_0 * hh_284[k];

        t_285[k] = f_0 * hh_285[k];

        t_286[k] = f_0 * hh_286[k];

        t_287[k] = f_0 * hh_287[k];

        t_288[k] = f_0 * hh_288[k];

        t_289[k] = f_0 * hh_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, t_295, t_296, t_297, hh_290, \
                         hh_291, hh_292, hh_293, hh_294, hh_295, hh_296, \
                         hh_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_0 * hh_290[k];

        t_291[k] = f_0 * hh_291[k];

        t_292[k] = f_0 * hh_292[k];

        t_293[k] = f_0 * hh_293[k];

        t_294[k] = f_0 * hh_294[k];

        t_295[k] = f_0 * hh_295[k];

        t_296[k] = f_0 * hh_296[k];

        t_297[k] = f_0 * hh_297[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, t_303, t_304, t_305, hh_298, \
                         hh_299, hh_300, hh_301, hh_302, hh_303, hh_304, \
                         hh_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_0 * hh_298[k];

        t_299[k] = f_0 * hh_299[k];

        t_300[k] = f_0 * hh_300[k];

        t_301[k] = f_0 * hh_301[k];

        t_302[k] = f_0 * hh_302[k];

        t_303[k] = f_0 * hh_303[k];

        t_304[k] = f_0 * hh_304[k];

        t_305[k] = f_0 * hh_305[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, t_312, t_313, hh_306, \
                         hh_307, hh_308, hh_309, hh_310, hh_311, hh_312, \
                         hh_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_0 * hh_306[k];

        t_307[k] = f_0 * hh_307[k];

        t_308[k] = f_0 * hh_308[k];

        t_309[k] = f_0 * hh_309[k];

        t_310[k] = f_0 * hh_310[k];

        t_311[k] = f_0 * hh_311[k];

        t_312[k] = f_0 * hh_312[k];

        t_313[k] = f_0 * hh_313[k];
    }

#pragma omp simd aligned(t_314, hh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_0 * hh_314[k];
    }
}

auto
compute_prim_geom_10_gh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t fh, const size_t hh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gh_electron_repulsion_0_piece0(buffer, target, fh, hh, ncols, alpha);

    compute_prim_geom_10_gh_electron_repulsion_0_piece1(buffer, target, fh, hh, ncols, alpha);
}

static auto
compute_prim_geom_10_gh_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fh, const size_t hh,
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

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_79 = buffer.data(fh + 79);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_104 = buffer.data(fh + 104);
    const auto *fh_105 = buffer.data(fh + 105);

    const auto *hh_21 = buffer.data(hh + 21);
    const auto *hh_22 = buffer.data(hh + 22);
    const auto *hh_23 = buffer.data(hh + 23);
    const auto *hh_24 = buffer.data(hh + 24);
    const auto *hh_25 = buffer.data(hh + 25);
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_29 = buffer.data(hh + 29);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_33 = buffer.data(hh + 33);
    const auto *hh_34 = buffer.data(hh + 34);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_64 = buffer.data(hh + 64);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_67 = buffer.data(hh + 67);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_71 = buffer.data(hh + 71);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_74 = buffer.data(hh + 74);
    const auto *hh_75 = buffer.data(hh + 75);
    const auto *hh_76 = buffer.data(hh + 76);
    const auto *hh_77 = buffer.data(hh + 77);
    const auto *hh_78 = buffer.data(hh + 78);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_80 = buffer.data(hh + 80);
    const auto *hh_81 = buffer.data(hh + 81);
    const auto *hh_82 = buffer.data(hh + 82);
    const auto *hh_83 = buffer.data(hh + 83);
    const auto *hh_84 = buffer.data(hh + 84);
    const auto *hh_85 = buffer.data(hh + 85);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_88 = buffer.data(hh + 88);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_91 = buffer.data(hh + 91);
    const auto *hh_92 = buffer.data(hh + 92);
    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_94 = buffer.data(hh + 94);
    const auto *hh_95 = buffer.data(hh + 95);
    const auto *hh_96 = buffer.data(hh + 96);
    const auto *hh_97 = buffer.data(hh + 97);
    const auto *hh_98 = buffer.data(hh + 98);
    const auto *hh_99 = buffer.data(hh + 99);
    const auto *hh_100 = buffer.data(hh + 100);
    const auto *hh_101 = buffer.data(hh + 101);
    const auto *hh_102 = buffer.data(hh + 102);
    const auto *hh_103 = buffer.data(hh + 103);
    const auto *hh_104 = buffer.data(hh + 104);
    const auto *hh_126 = buffer.data(hh + 126);
    const auto *hh_127 = buffer.data(hh + 127);
    const auto *hh_128 = buffer.data(hh + 128);
    const auto *hh_129 = buffer.data(hh + 129);
    const auto *hh_130 = buffer.data(hh + 130);
    const auto *hh_131 = buffer.data(hh + 131);
    const auto *hh_132 = buffer.data(hh + 132);
    const auto *hh_133 = buffer.data(hh + 133);
    const auto *hh_134 = buffer.data(hh + 134);
    const auto *hh_135 = buffer.data(hh + 135);
    const auto *hh_136 = buffer.data(hh + 136);
    const auto *hh_137 = buffer.data(hh + 137);
    const auto *hh_138 = buffer.data(hh + 138);
    const auto *hh_139 = buffer.data(hh + 139);
    const auto *hh_140 = buffer.data(hh + 140);
    const auto *hh_141 = buffer.data(hh + 141);
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_146 = buffer.data(hh + 146);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_148 = buffer.data(hh + 148);
    const auto *hh_149 = buffer.data(hh + 149);
    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_151 = buffer.data(hh + 151);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_155 = buffer.data(hh + 155);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_158 = buffer.data(hh + 158);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_160 = buffer.data(hh + 160);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_169 = buffer.data(hh + 169);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_172 = buffer.data(hh + 172);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_178 = buffer.data(hh + 178);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_181 = buffer.data(hh + 181);
    const auto *hh_182 = buffer.data(hh + 182);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_210 = buffer.data(hh + 210);
    const auto *hh_211 = buffer.data(hh + 211);
    const auto *hh_212 = buffer.data(hh + 212);
    const auto *hh_213 = buffer.data(hh + 213);
    const auto *hh_214 = buffer.data(hh + 214);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_218 = buffer.data(hh + 218);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_221 = buffer.data(hh + 221);
    const auto *hh_222 = buffer.data(hh + 222);
    const auto *hh_223 = buffer.data(hh + 223);
    const auto *hh_224 = buffer.data(hh + 224);
    const auto *hh_225 = buffer.data(hh + 225);
    const auto *hh_226 = buffer.data(hh + 226);
    const auto *hh_227 = buffer.data(hh + 227);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_229 = buffer.data(hh + 229);
    const auto *hh_230 = buffer.data(hh + 230);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_232 = buffer.data(hh + 232);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_235 = buffer.data(hh + 235);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_238 = buffer.data(hh + 238);
    const auto *hh_239 = buffer.data(hh + 239);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_243 = buffer.data(hh + 243);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_246 = buffer.data(hh + 246);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hh_21, hh_22, hh_23, hh_24, \
                         hh_25, hh_26, hh_27, hh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hh_21[k];

        t_1[k] = f_0 * hh_22[k];

        t_2[k] = f_0 * hh_23[k];

        t_3[k] = f_0 * hh_24[k];

        t_4[k] = f_0 * hh_25[k];

        t_5[k] = f_0 * hh_26[k];

        t_6[k] = f_0 * hh_27[k];

        t_7[k] = f_0 * hh_28[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hh_29, hh_30, hh_31, \
                         hh_32, hh_33, hh_34, hh_35, hh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hh_29[k];

        t_9[k] = f_0 * hh_30[k];

        t_10[k] = f_0 * hh_31[k];

        t_11[k] = f_0 * hh_32[k];

        t_12[k] = f_0 * hh_33[k];

        t_13[k] = f_0 * hh_34[k];

        t_14[k] = f_0 * hh_35[k];

        t_15[k] = f_0 * hh_36[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, fh_0, fh_1, hh_37, hh_38, \
                         hh_39, hh_40, hh_41, hh_63, hh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hh_37[k];

        t_17[k] = f_0 * hh_38[k];

        t_18[k] = f_0 * hh_39[k];

        t_19[k] = f_0 * hh_40[k];

        t_20[k] = f_0 * hh_41[k];

        t_21[k] = -fh_0[k]
                  + f_0 * hh_63[k];

        t_22[k] = -fh_1[k]
                  + f_0 * hh_64[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, fh_2, fh_3, fh_4, fh_5, fh_6, hh_65, \
                         hh_66, hh_67, hh_68, hh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -fh_2[k]
                  + f_0 * hh_65[k];

        t_24[k] = -fh_3[k]
                  + f_0 * hh_66[k];

        t_25[k] = -fh_4[k]
                  + f_0 * hh_67[k];

        t_26[k] = -fh_5[k]
                  + f_0 * hh_68[k];

        t_27[k] = -fh_6[k]
                  + f_0 * hh_69[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, fh_7, fh_8, fh_9, fh_10, fh_11, hh_70, \
                         hh_71, hh_72, hh_73, hh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -fh_7[k]
                  + f_0 * hh_70[k];

        t_29[k] = -fh_8[k]
                  + f_0 * hh_71[k];

        t_30[k] = -fh_9[k]
                  + f_0 * hh_72[k];

        t_31[k] = -fh_10[k]
                  + f_0 * hh_73[k];

        t_32[k] = -fh_11[k]
                  + f_0 * hh_74[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, fh_12, fh_13, fh_14, fh_15, fh_16, \
                         hh_75, hh_76, hh_77, hh_78, hh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -fh_12[k]
                  + f_0 * hh_75[k];

        t_34[k] = -fh_13[k]
                  + f_0 * hh_76[k];

        t_35[k] = -fh_14[k]
                  + f_0 * hh_77[k];

        t_36[k] = -fh_15[k]
                  + f_0 * hh_78[k];

        t_37[k] = -fh_16[k]
                  + f_0 * hh_79[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, fh_17, fh_18, fh_19, fh_20, \
                         hh_80, hh_81, hh_82, hh_83, hh_84, hh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -fh_17[k]
                  + f_0 * hh_80[k];

        t_39[k] = -fh_18[k]
                  + f_0 * hh_81[k];

        t_40[k] = -fh_19[k]
                  + f_0 * hh_82[k];

        t_41[k] = -fh_20[k]
                  + f_0 * hh_83[k];

        t_42[k] = f_0 * hh_84[k];

        t_43[k] = f_0 * hh_85[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, t_51, hh_86, hh_87, hh_88, \
                         hh_89, hh_90, hh_91, hh_92, hh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * hh_86[k];

        t_45[k] = f_0 * hh_87[k];

        t_46[k] = f_0 * hh_88[k];

        t_47[k] = f_0 * hh_89[k];

        t_48[k] = f_0 * hh_90[k];

        t_49[k] = f_0 * hh_91[k];

        t_50[k] = f_0 * hh_92[k];

        t_51[k] = f_0 * hh_93[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, t_58, t_59, hh_94, hh_95, hh_96, \
                         hh_97, hh_98, hh_99, hh_100, hh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * hh_94[k];

        t_53[k] = f_0 * hh_95[k];

        t_54[k] = f_0 * hh_96[k];

        t_55[k] = f_0 * hh_97[k];

        t_56[k] = f_0 * hh_98[k];

        t_57[k] = f_0 * hh_99[k];

        t_58[k] = f_0 * hh_100[k];

        t_59[k] = f_0 * hh_101[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, fh_21, fh_22, fh_23, hh_102, \
                         hh_103, hh_104, hh_126, hh_127, hh_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * hh_102[k];

        t_61[k] = f_0 * hh_103[k];

        t_62[k] = f_0 * hh_104[k];

        t_63[k] = -2.0 * fh_21[k]
                  + f_0 * hh_126[k];

        t_64[k] = -2.0 * fh_22[k]
                  + f_0 * hh_127[k];

        t_65[k] = -2.0 * fh_23[k]
                  + f_0 * hh_128[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, fh_24, fh_25, fh_26, fh_27, fh_28, \
                         hh_129, hh_130, hh_131, hh_132, hh_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -2.0 * fh_24[k]
                  + f_0 * hh_129[k];

        t_67[k] = -2.0 * fh_25[k]
                  + f_0 * hh_130[k];

        t_68[k] = -2.0 * fh_26[k]
                  + f_0 * hh_131[k];

        t_69[k] = -2.0 * fh_27[k]
                  + f_0 * hh_132[k];

        t_70[k] = -2.0 * fh_28[k]
                  + f_0 * hh_133[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, fh_29, fh_30, fh_31, fh_32, fh_33, \
                         hh_134, hh_135, hh_136, hh_137, hh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -2.0 * fh_29[k]
                  + f_0 * hh_134[k];

        t_72[k] = -2.0 * fh_30[k]
                  + f_0 * hh_135[k];

        t_73[k] = -2.0 * fh_31[k]
                  + f_0 * hh_136[k];

        t_74[k] = -2.0 * fh_32[k]
                  + f_0 * hh_137[k];

        t_75[k] = -2.0 * fh_33[k]
                  + f_0 * hh_138[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, fh_34, fh_35, fh_36, fh_37, fh_38, \
                         hh_139, hh_140, hh_141, hh_142, hh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * fh_34[k]
                  + f_0 * hh_139[k];

        t_77[k] = -2.0 * fh_35[k]
                  + f_0 * hh_140[k];

        t_78[k] = -2.0 * fh_36[k]
                  + f_0 * hh_141[k];

        t_79[k] = -2.0 * fh_37[k]
                  + f_0 * hh_142[k];

        t_80[k] = -2.0 * fh_38[k]
                  + f_0 * hh_143[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, fh_39, fh_40, fh_41, fh_42, fh_43, \
                         hh_144, hh_145, hh_146, hh_147, hh_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -2.0 * fh_39[k]
                  + f_0 * hh_144[k];

        t_82[k] = -2.0 * fh_40[k]
                  + f_0 * hh_145[k];

        t_83[k] = -2.0 * fh_41[k]
                  + f_0 * hh_146[k];

        t_84[k] = -fh_42[k]
                  + f_0 * hh_147[k];

        t_85[k] = -fh_43[k]
                  + f_0 * hh_148[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, fh_44, fh_45, fh_46, fh_47, fh_48, \
                         hh_149, hh_150, hh_151, hh_152, hh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -fh_44[k]
                  + f_0 * hh_149[k];

        t_87[k] = -fh_45[k]
                  + f_0 * hh_150[k];

        t_88[k] = -fh_46[k]
                  + f_0 * hh_151[k];

        t_89[k] = -fh_47[k]
                  + f_0 * hh_152[k];

        t_90[k] = -fh_48[k]
                  + f_0 * hh_153[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, fh_49, fh_50, fh_51, fh_52, fh_53, \
                         hh_154, hh_155, hh_156, hh_157, hh_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = -fh_49[k]
                  + f_0 * hh_154[k];

        t_92[k] = -fh_50[k]
                  + f_0 * hh_155[k];

        t_93[k] = -fh_51[k]
                  + f_0 * hh_156[k];

        t_94[k] = -fh_52[k]
                  + f_0 * hh_157[k];

        t_95[k] = -fh_53[k]
                  + f_0 * hh_158[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, fh_54, fh_55, fh_56, fh_57, fh_58, \
                         hh_159, hh_160, hh_161, hh_162, hh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -fh_54[k]
                  + f_0 * hh_159[k];

        t_97[k] = -fh_55[k]
                  + f_0 * hh_160[k];

        t_98[k] = -fh_56[k]
                  + f_0 * hh_161[k];

        t_99[k] = -fh_57[k]
                  + f_0 * hh_162[k];

        t_100[k] = -fh_58[k]
                   + f_0 * hh_163[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, fh_59, fh_60, fh_61, fh_62, \
                         hh_164, hh_165, hh_166, hh_167, hh_168, \
                         hh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -fh_59[k]
                   + f_0 * hh_164[k];

        t_102[k] = -fh_60[k]
                   + f_0 * hh_165[k];

        t_103[k] = -fh_61[k]
                   + f_0 * hh_166[k];

        t_104[k] = -fh_62[k]
                   + f_0 * hh_167[k];

        t_105[k] = f_0 * hh_168[k];

        t_106[k] = f_0 * hh_169[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, t_114, hh_170, \
                         hh_171, hh_172, hh_173, hh_174, hh_175, hh_176, \
                         hh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * hh_170[k];

        t_108[k] = f_0 * hh_171[k];

        t_109[k] = f_0 * hh_172[k];

        t_110[k] = f_0 * hh_173[k];

        t_111[k] = f_0 * hh_174[k];

        t_112[k] = f_0 * hh_175[k];

        t_113[k] = f_0 * hh_176[k];

        t_114[k] = f_0 * hh_177[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, t_121, t_122, hh_178, \
                         hh_179, hh_180, hh_181, hh_182, hh_183, hh_184, \
                         hh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_0 * hh_178[k];

        t_116[k] = f_0 * hh_179[k];

        t_117[k] = f_0 * hh_180[k];

        t_118[k] = f_0 * hh_181[k];

        t_119[k] = f_0 * hh_182[k];

        t_120[k] = f_0 * hh_183[k];

        t_121[k] = f_0 * hh_184[k];

        t_122[k] = f_0 * hh_185[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, t_128, fh_63, fh_64, fh_65, \
                         hh_186, hh_187, hh_188, hh_210, hh_211, \
                         hh_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * hh_186[k];

        t_124[k] = f_0 * hh_187[k];

        t_125[k] = f_0 * hh_188[k];

        t_126[k] = -3.0 * fh_63[k]
                   + f_0 * hh_210[k];

        t_127[k] = -3.0 * fh_64[k]
                   + f_0 * hh_211[k];

        t_128[k] = -3.0 * fh_65[k]
                   + f_0 * hh_212[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, fh_66, fh_67, fh_68, fh_69, fh_70, \
                         hh_213, hh_214, hh_215, hh_216, hh_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -3.0 * fh_66[k]
                   + f_0 * hh_213[k];

        t_130[k] = -3.0 * fh_67[k]
                   + f_0 * hh_214[k];

        t_131[k] = -3.0 * fh_68[k]
                   + f_0 * hh_215[k];

        t_132[k] = -3.0 * fh_69[k]
                   + f_0 * hh_216[k];

        t_133[k] = -3.0 * fh_70[k]
                   + f_0 * hh_217[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, fh_71, fh_72, fh_73, fh_74, fh_75, \
                         hh_218, hh_219, hh_220, hh_221, hh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -3.0 * fh_71[k]
                   + f_0 * hh_218[k];

        t_135[k] = -3.0 * fh_72[k]
                   + f_0 * hh_219[k];

        t_136[k] = -3.0 * fh_73[k]
                   + f_0 * hh_220[k];

        t_137[k] = -3.0 * fh_74[k]
                   + f_0 * hh_221[k];

        t_138[k] = -3.0 * fh_75[k]
                   + f_0 * hh_222[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, fh_76, fh_77, fh_78, fh_79, fh_80, \
                         hh_223, hh_224, hh_225, hh_226, hh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -3.0 * fh_76[k]
                   + f_0 * hh_223[k];

        t_140[k] = -3.0 * fh_77[k]
                   + f_0 * hh_224[k];

        t_141[k] = -3.0 * fh_78[k]
                   + f_0 * hh_225[k];

        t_142[k] = -3.0 * fh_79[k]
                   + f_0 * hh_226[k];

        t_143[k] = -3.0 * fh_80[k]
                   + f_0 * hh_227[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, fh_81, fh_82, fh_83, fh_84, fh_85, \
                         hh_228, hh_229, hh_230, hh_231, hh_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -3.0 * fh_81[k]
                   + f_0 * hh_228[k];

        t_145[k] = -3.0 * fh_82[k]
                   + f_0 * hh_229[k];

        t_146[k] = -3.0 * fh_83[k]
                   + f_0 * hh_230[k];

        t_147[k] = -2.0 * fh_84[k]
                   + f_0 * hh_231[k];

        t_148[k] = -2.0 * fh_85[k]
                   + f_0 * hh_232[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, fh_86, fh_87, fh_88, fh_89, fh_90, \
                         hh_233, hh_234, hh_235, hh_236, hh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * fh_86[k]
                   + f_0 * hh_233[k];

        t_150[k] = -2.0 * fh_87[k]
                   + f_0 * hh_234[k];

        t_151[k] = -2.0 * fh_88[k]
                   + f_0 * hh_235[k];

        t_152[k] = -2.0 * fh_89[k]
                   + f_0 * hh_236[k];

        t_153[k] = -2.0 * fh_90[k]
                   + f_0 * hh_237[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, fh_91, fh_92, fh_93, fh_94, fh_95, \
                         hh_238, hh_239, hh_240, hh_241, hh_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * fh_91[k]
                   + f_0 * hh_238[k];

        t_155[k] = -2.0 * fh_92[k]
                   + f_0 * hh_239[k];

        t_156[k] = -2.0 * fh_93[k]
                   + f_0 * hh_240[k];

        t_157[k] = -2.0 * fh_94[k]
                   + f_0 * hh_241[k];

        t_158[k] = -2.0 * fh_95[k]
                   + f_0 * hh_242[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, fh_96, fh_97, fh_98, fh_99, \
                         fh_100, hh_243, hh_244, hh_245, hh_246, \
                         hh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * fh_96[k]
                   + f_0 * hh_243[k];

        t_160[k] = -2.0 * fh_97[k]
                   + f_0 * hh_244[k];

        t_161[k] = -2.0 * fh_98[k]
                   + f_0 * hh_245[k];

        t_162[k] = -2.0 * fh_99[k]
                   + f_0 * hh_246[k];

        t_163[k] = -2.0 * fh_100[k]
                   + f_0 * hh_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, fh_101, fh_102, fh_103, fh_104, \
                         fh_105, hh_248, hh_249, hh_250, hh_251, \
                         hh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * fh_101[k]
                   + f_0 * hh_248[k];

        t_165[k] = -2.0 * fh_102[k]
                   + f_0 * hh_249[k];

        t_166[k] = -2.0 * fh_103[k]
                   + f_0 * hh_250[k];

        t_167[k] = -2.0 * fh_104[k]
                   + f_0 * hh_251[k];

        t_168[k] = -fh_105[k]
                   + f_0 * hh_252[k];
    }
}

static auto
compute_prim_geom_10_gh_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fh, const size_t hh,
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

    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_108 = buffer.data(fh + 108);
    const auto *fh_109 = buffer.data(fh + 109);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_111 = buffer.data(fh + 111);
    const auto *fh_112 = buffer.data(fh + 112);
    const auto *fh_113 = buffer.data(fh + 113);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_124 = buffer.data(fh + 124);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_128 = buffer.data(fh + 128);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_130 = buffer.data(fh + 130);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_133 = buffer.data(fh + 133);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_137 = buffer.data(fh + 137);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_147 = buffer.data(fh + 147);
    const auto *fh_148 = buffer.data(fh + 148);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_150 = buffer.data(fh + 150);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_153 = buffer.data(fh + 153);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_157 = buffer.data(fh + 157);
    const auto *fh_158 = buffer.data(fh + 158);
    const auto *fh_159 = buffer.data(fh + 159);
    const auto *fh_160 = buffer.data(fh + 160);
    const auto *fh_161 = buffer.data(fh + 161);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_168 = buffer.data(fh + 168);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_170 = buffer.data(fh + 170);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_173 = buffer.data(fh + 173);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_177 = buffer.data(fh + 177);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_182 = buffer.data(fh + 182);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_190 = buffer.data(fh + 190);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_193 = buffer.data(fh + 193);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_197 = buffer.data(fh + 197);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_202 = buffer.data(fh + 202);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *hh_253 = buffer.data(hh + 253);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_256 = buffer.data(hh + 256);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_260 = buffer.data(hh + 260);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_262 = buffer.data(hh + 262);
    const auto *hh_263 = buffer.data(hh + 263);
    const auto *hh_264 = buffer.data(hh + 264);
    const auto *hh_265 = buffer.data(hh + 265);
    const auto *hh_266 = buffer.data(hh + 266);
    const auto *hh_267 = buffer.data(hh + 267);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_269 = buffer.data(hh + 269);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_271 = buffer.data(hh + 271);
    const auto *hh_272 = buffer.data(hh + 272);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_274 = buffer.data(hh + 274);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_277 = buffer.data(hh + 277);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_280 = buffer.data(hh + 280);
    const auto *hh_281 = buffer.data(hh + 281);
    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_283 = buffer.data(hh + 283);
    const auto *hh_284 = buffer.data(hh + 284);
    const auto *hh_285 = buffer.data(hh + 285);
    const auto *hh_286 = buffer.data(hh + 286);
    const auto *hh_287 = buffer.data(hh + 287);
    const auto *hh_288 = buffer.data(hh + 288);
    const auto *hh_289 = buffer.data(hh + 289);
    const auto *hh_290 = buffer.data(hh + 290);
    const auto *hh_291 = buffer.data(hh + 291);
    const auto *hh_292 = buffer.data(hh + 292);
    const auto *hh_293 = buffer.data(hh + 293);
    const auto *hh_315 = buffer.data(hh + 315);
    const auto *hh_316 = buffer.data(hh + 316);
    const auto *hh_317 = buffer.data(hh + 317);
    const auto *hh_318 = buffer.data(hh + 318);
    const auto *hh_319 = buffer.data(hh + 319);
    const auto *hh_320 = buffer.data(hh + 320);
    const auto *hh_321 = buffer.data(hh + 321);
    const auto *hh_322 = buffer.data(hh + 322);
    const auto *hh_323 = buffer.data(hh + 323);
    const auto *hh_324 = buffer.data(hh + 324);
    const auto *hh_325 = buffer.data(hh + 325);
    const auto *hh_326 = buffer.data(hh + 326);
    const auto *hh_327 = buffer.data(hh + 327);
    const auto *hh_328 = buffer.data(hh + 328);
    const auto *hh_329 = buffer.data(hh + 329);
    const auto *hh_330 = buffer.data(hh + 330);
    const auto *hh_331 = buffer.data(hh + 331);
    const auto *hh_332 = buffer.data(hh + 332);
    const auto *hh_333 = buffer.data(hh + 333);
    const auto *hh_334 = buffer.data(hh + 334);
    const auto *hh_335 = buffer.data(hh + 335);
    const auto *hh_336 = buffer.data(hh + 336);
    const auto *hh_337 = buffer.data(hh + 337);
    const auto *hh_338 = buffer.data(hh + 338);
    const auto *hh_339 = buffer.data(hh + 339);
    const auto *hh_340 = buffer.data(hh + 340);
    const auto *hh_341 = buffer.data(hh + 341);
    const auto *hh_342 = buffer.data(hh + 342);
    const auto *hh_343 = buffer.data(hh + 343);
    const auto *hh_344 = buffer.data(hh + 344);
    const auto *hh_345 = buffer.data(hh + 345);
    const auto *hh_346 = buffer.data(hh + 346);
    const auto *hh_347 = buffer.data(hh + 347);
    const auto *hh_348 = buffer.data(hh + 348);
    const auto *hh_349 = buffer.data(hh + 349);
    const auto *hh_350 = buffer.data(hh + 350);
    const auto *hh_351 = buffer.data(hh + 351);
    const auto *hh_352 = buffer.data(hh + 352);
    const auto *hh_353 = buffer.data(hh + 353);
    const auto *hh_354 = buffer.data(hh + 354);
    const auto *hh_355 = buffer.data(hh + 355);
    const auto *hh_356 = buffer.data(hh + 356);
    const auto *hh_357 = buffer.data(hh + 357);
    const auto *hh_358 = buffer.data(hh + 358);
    const auto *hh_359 = buffer.data(hh + 359);
    const auto *hh_360 = buffer.data(hh + 360);
    const auto *hh_361 = buffer.data(hh + 361);
    const auto *hh_362 = buffer.data(hh + 362);
    const auto *hh_363 = buffer.data(hh + 363);
    const auto *hh_364 = buffer.data(hh + 364);
    const auto *hh_365 = buffer.data(hh + 365);
    const auto *hh_366 = buffer.data(hh + 366);
    const auto *hh_367 = buffer.data(hh + 367);
    const auto *hh_368 = buffer.data(hh + 368);
    const auto *hh_369 = buffer.data(hh + 369);
    const auto *hh_370 = buffer.data(hh + 370);
    const auto *hh_371 = buffer.data(hh + 371);
    const auto *hh_372 = buffer.data(hh + 372);
    const auto *hh_373 = buffer.data(hh + 373);
    const auto *hh_374 = buffer.data(hh + 374);
    const auto *hh_375 = buffer.data(hh + 375);
    const auto *hh_376 = buffer.data(hh + 376);
    const auto *hh_377 = buffer.data(hh + 377);
    const auto *hh_378 = buffer.data(hh + 378);
    const auto *hh_379 = buffer.data(hh + 379);
    const auto *hh_380 = buffer.data(hh + 380);
    const auto *hh_381 = buffer.data(hh + 381);
    const auto *hh_382 = buffer.data(hh + 382);
    const auto *hh_383 = buffer.data(hh + 383);
    const auto *hh_384 = buffer.data(hh + 384);
    const auto *hh_385 = buffer.data(hh + 385);
    const auto *hh_386 = buffer.data(hh + 386);
    const auto *hh_387 = buffer.data(hh + 387);
    const auto *hh_388 = buffer.data(hh + 388);
    const auto *hh_389 = buffer.data(hh + 389);
    const auto *hh_390 = buffer.data(hh + 390);
    const auto *hh_391 = buffer.data(hh + 391);
    const auto *hh_392 = buffer.data(hh + 392);
    const auto *hh_393 = buffer.data(hh + 393);
    const auto *hh_394 = buffer.data(hh + 394);
    const auto *hh_395 = buffer.data(hh + 395);
    const auto *hh_396 = buffer.data(hh + 396);
    const auto *hh_397 = buffer.data(hh + 397);
    const auto *hh_398 = buffer.data(hh + 398);
    const auto *hh_399 = buffer.data(hh + 399);
    const auto *hh_400 = buffer.data(hh + 400);
    const auto *hh_401 = buffer.data(hh + 401);
    const auto *hh_402 = buffer.data(hh + 402);
    const auto *hh_403 = buffer.data(hh + 403);
    const auto *hh_404 = buffer.data(hh + 404);
    const auto *hh_405 = buffer.data(hh + 405);
    const auto *hh_406 = buffer.data(hh + 406);
    const auto *hh_407 = buffer.data(hh + 407);
    const auto *hh_408 = buffer.data(hh + 408);
    const auto *hh_409 = buffer.data(hh + 409);
    const auto *hh_410 = buffer.data(hh + 410);
    const auto *hh_411 = buffer.data(hh + 411);
    const auto *hh_412 = buffer.data(hh + 412);
    const auto *hh_413 = buffer.data(hh + 413);
    const auto *hh_414 = buffer.data(hh + 414);
    const auto *hh_415 = buffer.data(hh + 415);
    const auto *hh_416 = buffer.data(hh + 416);
    const auto *hh_417 = buffer.data(hh + 417);
    const auto *hh_418 = buffer.data(hh + 418);
    const auto *hh_419 = buffer.data(hh + 419);

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, fh_106, fh_107, fh_108, fh_109, \
                         fh_110, hh_253, hh_254, hh_255, hh_256, \
                         hh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -fh_106[k]
                   + f_0 * hh_253[k];

        t_170[k] = -fh_107[k]
                   + f_0 * hh_254[k];

        t_171[k] = -fh_108[k]
                   + f_0 * hh_255[k];

        t_172[k] = -fh_109[k]
                   + f_0 * hh_256[k];

        t_173[k] = -fh_110[k]
                   + f_0 * hh_257[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, fh_111, fh_112, fh_113, fh_114, \
                         fh_115, hh_258, hh_259, hh_260, hh_261, \
                         hh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -fh_111[k]
                   + f_0 * hh_258[k];

        t_175[k] = -fh_112[k]
                   + f_0 * hh_259[k];

        t_176[k] = -fh_113[k]
                   + f_0 * hh_260[k];

        t_177[k] = -fh_114[k]
                   + f_0 * hh_261[k];

        t_178[k] = -fh_115[k]
                   + f_0 * hh_262[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, fh_116, fh_117, fh_118, fh_119, \
                         fh_120, hh_263, hh_264, hh_265, hh_266, \
                         hh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -fh_116[k]
                   + f_0 * hh_263[k];

        t_180[k] = -fh_117[k]
                   + f_0 * hh_264[k];

        t_181[k] = -fh_118[k]
                   + f_0 * hh_265[k];

        t_182[k] = -fh_119[k]
                   + f_0 * hh_266[k];

        t_183[k] = -fh_120[k]
                   + f_0 * hh_267[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, fh_121, fh_122, fh_123, fh_124, \
                         fh_125, hh_268, hh_269, hh_270, hh_271, \
                         hh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -fh_121[k]
                   + f_0 * hh_268[k];

        t_185[k] = -fh_122[k]
                   + f_0 * hh_269[k];

        t_186[k] = -fh_123[k]
                   + f_0 * hh_270[k];

        t_187[k] = -fh_124[k]
                   + f_0 * hh_271[k];

        t_188[k] = -fh_125[k]
                   + f_0 * hh_272[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, t_196, hh_273, \
                         hh_274, hh_275, hh_276, hh_277, hh_278, hh_279, \
                         hh_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * hh_273[k];

        t_190[k] = f_0 * hh_274[k];

        t_191[k] = f_0 * hh_275[k];

        t_192[k] = f_0 * hh_276[k];

        t_193[k] = f_0 * hh_277[k];

        t_194[k] = f_0 * hh_278[k];

        t_195[k] = f_0 * hh_279[k];

        t_196[k] = f_0 * hh_280[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, t_203, t_204, hh_281, \
                         hh_282, hh_283, hh_284, hh_285, hh_286, hh_287, \
                         hh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_0 * hh_281[k];

        t_198[k] = f_0 * hh_282[k];

        t_199[k] = f_0 * hh_283[k];

        t_200[k] = f_0 * hh_284[k];

        t_201[k] = f_0 * hh_285[k];

        t_202[k] = f_0 * hh_286[k];

        t_203[k] = f_0 * hh_287[k];

        t_204[k] = f_0 * hh_288[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, t_210, t_211, fh_126, fh_127, \
                         hh_289, hh_290, hh_291, hh_292, hh_293, hh_315, \
                         hh_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_0 * hh_289[k];

        t_206[k] = f_0 * hh_290[k];

        t_207[k] = f_0 * hh_291[k];

        t_208[k] = f_0 * hh_292[k];

        t_209[k] = f_0 * hh_293[k];

        t_210[k] = -4.0 * fh_126[k]
                   + f_0 * hh_315[k];

        t_211[k] = -4.0 * fh_127[k]
                   + f_0 * hh_316[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, fh_128, fh_129, fh_130, fh_131, \
                         fh_132, hh_317, hh_318, hh_319, hh_320, \
                         hh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -4.0 * fh_128[k]
                   + f_0 * hh_317[k];

        t_213[k] = -4.0 * fh_129[k]
                   + f_0 * hh_318[k];

        t_214[k] = -4.0 * fh_130[k]
                   + f_0 * hh_319[k];

        t_215[k] = -4.0 * fh_131[k]
                   + f_0 * hh_320[k];

        t_216[k] = -4.0 * fh_132[k]
                   + f_0 * hh_321[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, fh_133, fh_134, fh_135, fh_136, \
                         fh_137, hh_322, hh_323, hh_324, hh_325, \
                         hh_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = -4.0 * fh_133[k]
                   + f_0 * hh_322[k];

        t_218[k] = -4.0 * fh_134[k]
                   + f_0 * hh_323[k];

        t_219[k] = -4.0 * fh_135[k]
                   + f_0 * hh_324[k];

        t_220[k] = -4.0 * fh_136[k]
                   + f_0 * hh_325[k];

        t_221[k] = -4.0 * fh_137[k]
                   + f_0 * hh_326[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, fh_138, fh_139, fh_140, fh_141, \
                         fh_142, hh_327, hh_328, hh_329, hh_330, \
                         hh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -4.0 * fh_138[k]
                   + f_0 * hh_327[k];

        t_223[k] = -4.0 * fh_139[k]
                   + f_0 * hh_328[k];

        t_224[k] = -4.0 * fh_140[k]
                   + f_0 * hh_329[k];

        t_225[k] = -4.0 * fh_141[k]
                   + f_0 * hh_330[k];

        t_226[k] = -4.0 * fh_142[k]
                   + f_0 * hh_331[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, fh_143, fh_144, fh_145, fh_146, \
                         fh_147, hh_332, hh_333, hh_334, hh_335, \
                         hh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = -4.0 * fh_143[k]
                   + f_0 * hh_332[k];

        t_228[k] = -4.0 * fh_144[k]
                   + f_0 * hh_333[k];

        t_229[k] = -4.0 * fh_145[k]
                   + f_0 * hh_334[k];

        t_230[k] = -4.0 * fh_146[k]
                   + f_0 * hh_335[k];

        t_231[k] = -3.0 * fh_147[k]
                   + f_0 * hh_336[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, fh_148, fh_149, fh_150, fh_151, \
                         fh_152, hh_337, hh_338, hh_339, hh_340, \
                         hh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = -3.0 * fh_148[k]
                   + f_0 * hh_337[k];

        t_233[k] = -3.0 * fh_149[k]
                   + f_0 * hh_338[k];

        t_234[k] = -3.0 * fh_150[k]
                   + f_0 * hh_339[k];

        t_235[k] = -3.0 * fh_151[k]
                   + f_0 * hh_340[k];

        t_236[k] = -3.0 * fh_152[k]
                   + f_0 * hh_341[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, fh_153, fh_154, fh_155, fh_156, \
                         fh_157, hh_342, hh_343, hh_344, hh_345, \
                         hh_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = -3.0 * fh_153[k]
                   + f_0 * hh_342[k];

        t_238[k] = -3.0 * fh_154[k]
                   + f_0 * hh_343[k];

        t_239[k] = -3.0 * fh_155[k]
                   + f_0 * hh_344[k];

        t_240[k] = -3.0 * fh_156[k]
                   + f_0 * hh_345[k];

        t_241[k] = -3.0 * fh_157[k]
                   + f_0 * hh_346[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, fh_158, fh_159, fh_160, fh_161, \
                         fh_162, hh_347, hh_348, hh_349, hh_350, \
                         hh_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = -3.0 * fh_158[k]
                   + f_0 * hh_347[k];

        t_243[k] = -3.0 * fh_159[k]
                   + f_0 * hh_348[k];

        t_244[k] = -3.0 * fh_160[k]
                   + f_0 * hh_349[k];

        t_245[k] = -3.0 * fh_161[k]
                   + f_0 * hh_350[k];

        t_246[k] = -3.0 * fh_162[k]
                   + f_0 * hh_351[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, fh_163, fh_164, fh_165, fh_166, \
                         fh_167, hh_352, hh_353, hh_354, hh_355, \
                         hh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = -3.0 * fh_163[k]
                   + f_0 * hh_352[k];

        t_248[k] = -3.0 * fh_164[k]
                   + f_0 * hh_353[k];

        t_249[k] = -3.0 * fh_165[k]
                   + f_0 * hh_354[k];

        t_250[k] = -3.0 * fh_166[k]
                   + f_0 * hh_355[k];

        t_251[k] = -3.0 * fh_167[k]
                   + f_0 * hh_356[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, fh_168, fh_169, fh_170, fh_171, \
                         fh_172, hh_357, hh_358, hh_359, hh_360, \
                         hh_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -2.0 * fh_168[k]
                   + f_0 * hh_357[k];

        t_253[k] = -2.0 * fh_169[k]
                   + f_0 * hh_358[k];

        t_254[k] = -2.0 * fh_170[k]
                   + f_0 * hh_359[k];

        t_255[k] = -2.0 * fh_171[k]
                   + f_0 * hh_360[k];

        t_256[k] = -2.0 * fh_172[k]
                   + f_0 * hh_361[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, fh_173, fh_174, fh_175, fh_176, \
                         fh_177, hh_362, hh_363, hh_364, hh_365, \
                         hh_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = -2.0 * fh_173[k]
                   + f_0 * hh_362[k];

        t_258[k] = -2.0 * fh_174[k]
                   + f_0 * hh_363[k];

        t_259[k] = -2.0 * fh_175[k]
                   + f_0 * hh_364[k];

        t_260[k] = -2.0 * fh_176[k]
                   + f_0 * hh_365[k];

        t_261[k] = -2.0 * fh_177[k]
                   + f_0 * hh_366[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, fh_178, fh_179, fh_180, fh_181, \
                         fh_182, hh_367, hh_368, hh_369, hh_370, \
                         hh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = -2.0 * fh_178[k]
                   + f_0 * hh_367[k];

        t_263[k] = -2.0 * fh_179[k]
                   + f_0 * hh_368[k];

        t_264[k] = -2.0 * fh_180[k]
                   + f_0 * hh_369[k];

        t_265[k] = -2.0 * fh_181[k]
                   + f_0 * hh_370[k];

        t_266[k] = -2.0 * fh_182[k]
                   + f_0 * hh_371[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, fh_183, fh_184, fh_185, fh_186, \
                         fh_187, hh_372, hh_373, hh_374, hh_375, \
                         hh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = -2.0 * fh_183[k]
                   + f_0 * hh_372[k];

        t_268[k] = -2.0 * fh_184[k]
                   + f_0 * hh_373[k];

        t_269[k] = -2.0 * fh_185[k]
                   + f_0 * hh_374[k];

        t_270[k] = -2.0 * fh_186[k]
                   + f_0 * hh_375[k];

        t_271[k] = -2.0 * fh_187[k]
                   + f_0 * hh_376[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, fh_188, fh_189, fh_190, fh_191, \
                         fh_192, hh_377, hh_378, hh_379, hh_380, \
                         hh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -2.0 * fh_188[k]
                   + f_0 * hh_377[k];

        t_273[k] = -fh_189[k]
                   + f_0 * hh_378[k];

        t_274[k] = -fh_190[k]
                   + f_0 * hh_379[k];

        t_275[k] = -fh_191[k]
                   + f_0 * hh_380[k];

        t_276[k] = -fh_192[k]
                   + f_0 * hh_381[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, t_281, fh_193, fh_194, fh_195, fh_196, \
                         fh_197, hh_382, hh_383, hh_384, hh_385, \
                         hh_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -fh_193[k]
                   + f_0 * hh_382[k];

        t_278[k] = -fh_194[k]
                   + f_0 * hh_383[k];

        t_279[k] = -fh_195[k]
                   + f_0 * hh_384[k];

        t_280[k] = -fh_196[k]
                   + f_0 * hh_385[k];

        t_281[k] = -fh_197[k]
                   + f_0 * hh_386[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, fh_198, fh_199, fh_200, fh_201, \
                         fh_202, hh_387, hh_388, hh_389, hh_390, \
                         hh_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = -fh_198[k]
                   + f_0 * hh_387[k];

        t_283[k] = -fh_199[k]
                   + f_0 * hh_388[k];

        t_284[k] = -fh_200[k]
                   + f_0 * hh_389[k];

        t_285[k] = -fh_201[k]
                   + f_0 * hh_390[k];

        t_286[k] = -fh_202[k]
                   + f_0 * hh_391[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, t_291, fh_203, fh_204, fh_205, fh_206, \
                         fh_207, hh_392, hh_393, hh_394, hh_395, \
                         hh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = -fh_203[k]
                   + f_0 * hh_392[k];

        t_288[k] = -fh_204[k]
                   + f_0 * hh_393[k];

        t_289[k] = -fh_205[k]
                   + f_0 * hh_394[k];

        t_290[k] = -fh_206[k]
                   + f_0 * hh_395[k];

        t_291[k] = -fh_207[k]
                   + f_0 * hh_396[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, t_297, t_298, fh_208, fh_209, \
                         hh_397, hh_398, hh_399, hh_400, hh_401, hh_402, \
                         hh_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = -fh_208[k]
                   + f_0 * hh_397[k];

        t_293[k] = -fh_209[k]
                   + f_0 * hh_398[k];

        t_294[k] = f_0 * hh_399[k];

        t_295[k] = f_0 * hh_400[k];

        t_296[k] = f_0 * hh_401[k];

        t_297[k] = f_0 * hh_402[k];

        t_298[k] = f_0 * hh_403[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, t_304, t_305, t_306, hh_404, \
                         hh_405, hh_406, hh_407, hh_408, hh_409, hh_410, \
                         hh_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_0 * hh_404[k];

        t_300[k] = f_0 * hh_405[k];

        t_301[k] = f_0 * hh_406[k];

        t_302[k] = f_0 * hh_407[k];

        t_303[k] = f_0 * hh_408[k];

        t_304[k] = f_0 * hh_409[k];

        t_305[k] = f_0 * hh_410[k];

        t_306[k] = f_0 * hh_411[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, t_312, t_313, t_314, hh_412, \
                         hh_413, hh_414, hh_415, hh_416, hh_417, hh_418, \
                         hh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_0 * hh_412[k];

        t_308[k] = f_0 * hh_413[k];

        t_309[k] = f_0 * hh_414[k];

        t_310[k] = f_0 * hh_415[k];

        t_311[k] = f_0 * hh_416[k];

        t_312[k] = f_0 * hh_417[k];

        t_313[k] = f_0 * hh_418[k];

        t_314[k] = f_0 * hh_419[k];
    }
}

auto
compute_prim_geom_10_gh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t fh, const size_t hh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gh_electron_repulsion_1_piece0(buffer, target, fh, hh, ncols, alpha);

    compute_prim_geom_10_gh_electron_repulsion_1_piece1(buffer, target, fh, hh, ncols, alpha);
}

static auto
compute_prim_geom_10_gh_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fh, const size_t hh,
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

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_79 = buffer.data(fh + 79);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);

    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_43 = buffer.data(hh + 43);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_45 = buffer.data(hh + 45);
    const auto *hh_46 = buffer.data(hh + 46);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_48 = buffer.data(hh + 48);
    const auto *hh_49 = buffer.data(hh + 49);
    const auto *hh_50 = buffer.data(hh + 50);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_52 = buffer.data(hh + 52);
    const auto *hh_53 = buffer.data(hh + 53);
    const auto *hh_54 = buffer.data(hh + 54);
    const auto *hh_55 = buffer.data(hh + 55);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_57 = buffer.data(hh + 57);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_84 = buffer.data(hh + 84);
    const auto *hh_85 = buffer.data(hh + 85);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_88 = buffer.data(hh + 88);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_91 = buffer.data(hh + 91);
    const auto *hh_92 = buffer.data(hh + 92);
    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_94 = buffer.data(hh + 94);
    const auto *hh_95 = buffer.data(hh + 95);
    const auto *hh_96 = buffer.data(hh + 96);
    const auto *hh_97 = buffer.data(hh + 97);
    const auto *hh_98 = buffer.data(hh + 98);
    const auto *hh_99 = buffer.data(hh + 99);
    const auto *hh_100 = buffer.data(hh + 100);
    const auto *hh_101 = buffer.data(hh + 101);
    const auto *hh_102 = buffer.data(hh + 102);
    const auto *hh_103 = buffer.data(hh + 103);
    const auto *hh_104 = buffer.data(hh + 104);
    const auto *hh_105 = buffer.data(hh + 105);
    const auto *hh_106 = buffer.data(hh + 106);
    const auto *hh_107 = buffer.data(hh + 107);
    const auto *hh_108 = buffer.data(hh + 108);
    const auto *hh_109 = buffer.data(hh + 109);
    const auto *hh_110 = buffer.data(hh + 110);
    const auto *hh_111 = buffer.data(hh + 111);
    const auto *hh_112 = buffer.data(hh + 112);
    const auto *hh_113 = buffer.data(hh + 113);
    const auto *hh_114 = buffer.data(hh + 114);
    const auto *hh_115 = buffer.data(hh + 115);
    const auto *hh_116 = buffer.data(hh + 116);
    const auto *hh_117 = buffer.data(hh + 117);
    const auto *hh_118 = buffer.data(hh + 118);
    const auto *hh_119 = buffer.data(hh + 119);
    const auto *hh_120 = buffer.data(hh + 120);
    const auto *hh_121 = buffer.data(hh + 121);
    const auto *hh_122 = buffer.data(hh + 122);
    const auto *hh_123 = buffer.data(hh + 123);
    const auto *hh_124 = buffer.data(hh + 124);
    const auto *hh_125 = buffer.data(hh + 125);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_148 = buffer.data(hh + 148);
    const auto *hh_149 = buffer.data(hh + 149);
    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_151 = buffer.data(hh + 151);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_155 = buffer.data(hh + 155);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_158 = buffer.data(hh + 158);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_160 = buffer.data(hh + 160);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_169 = buffer.data(hh + 169);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_172 = buffer.data(hh + 172);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_178 = buffer.data(hh + 178);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_181 = buffer.data(hh + 181);
    const auto *hh_182 = buffer.data(hh + 182);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_189 = buffer.data(hh + 189);
    const auto *hh_190 = buffer.data(hh + 190);
    const auto *hh_191 = buffer.data(hh + 191);
    const auto *hh_192 = buffer.data(hh + 192);
    const auto *hh_193 = buffer.data(hh + 193);
    const auto *hh_194 = buffer.data(hh + 194);
    const auto *hh_195 = buffer.data(hh + 195);
    const auto *hh_196 = buffer.data(hh + 196);
    const auto *hh_197 = buffer.data(hh + 197);
    const auto *hh_198 = buffer.data(hh + 198);
    const auto *hh_199 = buffer.data(hh + 199);
    const auto *hh_200 = buffer.data(hh + 200);
    const auto *hh_201 = buffer.data(hh + 201);
    const auto *hh_202 = buffer.data(hh + 202);
    const auto *hh_203 = buffer.data(hh + 203);
    const auto *hh_204 = buffer.data(hh + 204);
    const auto *hh_205 = buffer.data(hh + 205);
    const auto *hh_206 = buffer.data(hh + 206);
    const auto *hh_207 = buffer.data(hh + 207);
    const auto *hh_208 = buffer.data(hh + 208);
    const auto *hh_209 = buffer.data(hh + 209);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_232 = buffer.data(hh + 232);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_235 = buffer.data(hh + 235);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_238 = buffer.data(hh + 238);
    const auto *hh_239 = buffer.data(hh + 239);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_243 = buffer.data(hh + 243);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_246 = buffer.data(hh + 246);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_253 = buffer.data(hh + 253);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_256 = buffer.data(hh + 256);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_260 = buffer.data(hh + 260);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_262 = buffer.data(hh + 262);
    const auto *hh_263 = buffer.data(hh + 263);
    const auto *hh_264 = buffer.data(hh + 264);
    const auto *hh_265 = buffer.data(hh + 265);
    const auto *hh_266 = buffer.data(hh + 266);
    const auto *hh_267 = buffer.data(hh + 267);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_269 = buffer.data(hh + 269);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_271 = buffer.data(hh + 271);
    const auto *hh_272 = buffer.data(hh + 272);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_274 = buffer.data(hh + 274);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_277 = buffer.data(hh + 277);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_280 = buffer.data(hh + 280);
    const auto *hh_281 = buffer.data(hh + 281);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hh_42, hh_43, hh_44, hh_45, \
                         hh_46, hh_47, hh_48, hh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hh_42[k];

        t_1[k] = f_0 * hh_43[k];

        t_2[k] = f_0 * hh_44[k];

        t_3[k] = f_0 * hh_45[k];

        t_4[k] = f_0 * hh_46[k];

        t_5[k] = f_0 * hh_47[k];

        t_6[k] = f_0 * hh_48[k];

        t_7[k] = f_0 * hh_49[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hh_50, hh_51, hh_52, \
                         hh_53, hh_54, hh_55, hh_56, hh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hh_50[k];

        t_9[k] = f_0 * hh_51[k];

        t_10[k] = f_0 * hh_52[k];

        t_11[k] = f_0 * hh_53[k];

        t_12[k] = f_0 * hh_54[k];

        t_13[k] = f_0 * hh_55[k];

        t_14[k] = f_0 * hh_56[k];

        t_15[k] = f_0 * hh_57[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, hh_58, hh_59, hh_60, \
                         hh_61, hh_62, hh_84, hh_85, hh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hh_58[k];

        t_17[k] = f_0 * hh_59[k];

        t_18[k] = f_0 * hh_60[k];

        t_19[k] = f_0 * hh_61[k];

        t_20[k] = f_0 * hh_62[k];

        t_21[k] = f_0 * hh_84[k];

        t_22[k] = f_0 * hh_85[k];

        t_23[k] = f_0 * hh_86[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, hh_87, hh_88, hh_89, \
                         hh_90, hh_91, hh_92, hh_93, hh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * hh_87[k];

        t_25[k] = f_0 * hh_88[k];

        t_26[k] = f_0 * hh_89[k];

        t_27[k] = f_0 * hh_90[k];

        t_28[k] = f_0 * hh_91[k];

        t_29[k] = f_0 * hh_92[k];

        t_30[k] = f_0 * hh_93[k];

        t_31[k] = f_0 * hh_94[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, hh_95, hh_96, hh_97, \
                         hh_98, hh_99, hh_100, hh_101, hh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * hh_95[k];

        t_33[k] = f_0 * hh_96[k];

        t_34[k] = f_0 * hh_97[k];

        t_35[k] = f_0 * hh_98[k];

        t_36[k] = f_0 * hh_99[k];

        t_37[k] = f_0 * hh_100[k];

        t_38[k] = f_0 * hh_101[k];

        t_39[k] = f_0 * hh_102[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, fh_0, fh_1, fh_2, fh_3, hh_103, \
                         hh_104, hh_105, hh_106, hh_107, hh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * hh_103[k];

        t_41[k] = f_0 * hh_104[k];

        t_42[k] = -fh_0[k]
                  + f_0 * hh_105[k];

        t_43[k] = -fh_1[k]
                  + f_0 * hh_106[k];

        t_44[k] = -fh_2[k]
                  + f_0 * hh_107[k];

        t_45[k] = -fh_3[k]
                  + f_0 * hh_108[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, fh_4, fh_5, fh_6, fh_7, fh_8, hh_109, \
                         hh_110, hh_111, hh_112, hh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -fh_4[k]
                  + f_0 * hh_109[k];

        t_47[k] = -fh_5[k]
                  + f_0 * hh_110[k];

        t_48[k] = -fh_6[k]
                  + f_0 * hh_111[k];

        t_49[k] = -fh_7[k]
                  + f_0 * hh_112[k];

        t_50[k] = -fh_8[k]
                  + f_0 * hh_113[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, fh_9, fh_10, fh_11, fh_12, fh_13, \
                         hh_114, hh_115, hh_116, hh_117, hh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -fh_9[k]
                  + f_0 * hh_114[k];

        t_52[k] = -fh_10[k]
                  + f_0 * hh_115[k];

        t_53[k] = -fh_11[k]
                  + f_0 * hh_116[k];

        t_54[k] = -fh_12[k]
                  + f_0 * hh_117[k];

        t_55[k] = -fh_13[k]
                  + f_0 * hh_118[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, fh_14, fh_15, fh_16, fh_17, fh_18, \
                         hh_119, hh_120, hh_121, hh_122, hh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -fh_14[k]
                  + f_0 * hh_119[k];

        t_57[k] = -fh_15[k]
                  + f_0 * hh_120[k];

        t_58[k] = -fh_16[k]
                  + f_0 * hh_121[k];

        t_59[k] = -fh_17[k]
                  + f_0 * hh_122[k];

        t_60[k] = -fh_18[k]
                  + f_0 * hh_123[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, fh_19, fh_20, hh_124, \
                         hh_125, hh_147, hh_148, hh_149, hh_150, \
                         hh_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -fh_19[k]
                  + f_0 * hh_124[k];

        t_62[k] = -fh_20[k]
                  + f_0 * hh_125[k];

        t_63[k] = f_0 * hh_147[k];

        t_64[k] = f_0 * hh_148[k];

        t_65[k] = f_0 * hh_149[k];

        t_66[k] = f_0 * hh_150[k];

        t_67[k] = f_0 * hh_151[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, t_74, t_75, hh_152, hh_153, \
                         hh_154, hh_155, hh_156, hh_157, hh_158, \
                         hh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * hh_152[k];

        t_69[k] = f_0 * hh_153[k];

        t_70[k] = f_0 * hh_154[k];

        t_71[k] = f_0 * hh_155[k];

        t_72[k] = f_0 * hh_156[k];

        t_73[k] = f_0 * hh_157[k];

        t_74[k] = f_0 * hh_158[k];

        t_75[k] = f_0 * hh_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, t_83, hh_160, hh_161, \
                         hh_162, hh_163, hh_164, hh_165, hh_166, \
                         hh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_0 * hh_160[k];

        t_77[k] = f_0 * hh_161[k];

        t_78[k] = f_0 * hh_162[k];

        t_79[k] = f_0 * hh_163[k];

        t_80[k] = f_0 * hh_164[k];

        t_81[k] = f_0 * hh_165[k];

        t_82[k] = f_0 * hh_166[k];

        t_83[k] = f_0 * hh_167[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, fh_21, fh_22, fh_23, fh_24, fh_25, \
                         hh_168, hh_169, hh_170, hh_171, hh_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -fh_21[k]
                  + f_0 * hh_168[k];

        t_85[k] = -fh_22[k]
                  + f_0 * hh_169[k];

        t_86[k] = -fh_23[k]
                  + f_0 * hh_170[k];

        t_87[k] = -fh_24[k]
                  + f_0 * hh_171[k];

        t_88[k] = -fh_25[k]
                  + f_0 * hh_172[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, fh_26, fh_27, fh_28, fh_29, fh_30, \
                         hh_173, hh_174, hh_175, hh_176, hh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -fh_26[k]
                  + f_0 * hh_173[k];

        t_90[k] = -fh_27[k]
                  + f_0 * hh_174[k];

        t_91[k] = -fh_28[k]
                  + f_0 * hh_175[k];

        t_92[k] = -fh_29[k]
                  + f_0 * hh_176[k];

        t_93[k] = -fh_30[k]
                  + f_0 * hh_177[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, fh_31, fh_32, fh_33, fh_34, fh_35, \
                         hh_178, hh_179, hh_180, hh_181, hh_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -fh_31[k]
                  + f_0 * hh_178[k];

        t_95[k] = -fh_32[k]
                  + f_0 * hh_179[k];

        t_96[k] = -fh_33[k]
                  + f_0 * hh_180[k];

        t_97[k] = -fh_34[k]
                  + f_0 * hh_181[k];

        t_98[k] = -fh_35[k]
                  + f_0 * hh_182[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, fh_36, fh_37, fh_38, fh_39, fh_40, \
                         hh_183, hh_184, hh_185, hh_186, hh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -fh_36[k]
                  + f_0 * hh_183[k];

        t_100[k] = -fh_37[k]
                   + f_0 * hh_184[k];

        t_101[k] = -fh_38[k]
                   + f_0 * hh_185[k];

        t_102[k] = -fh_39[k]
                   + f_0 * hh_186[k];

        t_103[k] = -fh_40[k]
                   + f_0 * hh_187[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, fh_41, fh_42, fh_43, fh_44, fh_45, \
                         hh_188, hh_189, hh_190, hh_191, hh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -fh_41[k]
                   + f_0 * hh_188[k];

        t_105[k] = -2.0 * fh_42[k]
                   + f_0 * hh_189[k];

        t_106[k] = -2.0 * fh_43[k]
                   + f_0 * hh_190[k];

        t_107[k] = -2.0 * fh_44[k]
                   + f_0 * hh_191[k];

        t_108[k] = -2.0 * fh_45[k]
                   + f_0 * hh_192[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, fh_46, fh_47, fh_48, fh_49, fh_50, \
                         hh_193, hh_194, hh_195, hh_196, hh_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -2.0 * fh_46[k]
                   + f_0 * hh_193[k];

        t_110[k] = -2.0 * fh_47[k]
                   + f_0 * hh_194[k];

        t_111[k] = -2.0 * fh_48[k]
                   + f_0 * hh_195[k];

        t_112[k] = -2.0 * fh_49[k]
                   + f_0 * hh_196[k];

        t_113[k] = -2.0 * fh_50[k]
                   + f_0 * hh_197[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, fh_51, fh_52, fh_53, fh_54, fh_55, \
                         hh_198, hh_199, hh_200, hh_201, hh_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -2.0 * fh_51[k]
                   + f_0 * hh_198[k];

        t_115[k] = -2.0 * fh_52[k]
                   + f_0 * hh_199[k];

        t_116[k] = -2.0 * fh_53[k]
                   + f_0 * hh_200[k];

        t_117[k] = -2.0 * fh_54[k]
                   + f_0 * hh_201[k];

        t_118[k] = -2.0 * fh_55[k]
                   + f_0 * hh_202[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, fh_56, fh_57, fh_58, fh_59, fh_60, \
                         hh_203, hh_204, hh_205, hh_206, hh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -2.0 * fh_56[k]
                   + f_0 * hh_203[k];

        t_120[k] = -2.0 * fh_57[k]
                   + f_0 * hh_204[k];

        t_121[k] = -2.0 * fh_58[k]
                   + f_0 * hh_205[k];

        t_122[k] = -2.0 * fh_59[k]
                   + f_0 * hh_206[k];

        t_123[k] = -2.0 * fh_60[k]
                   + f_0 * hh_207[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, t_129, t_130, fh_61, fh_62, \
                         hh_208, hh_209, hh_231, hh_232, hh_233, hh_234, \
                         hh_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * fh_61[k]
                   + f_0 * hh_208[k];

        t_125[k] = -2.0 * fh_62[k]
                   + f_0 * hh_209[k];

        t_126[k] = f_0 * hh_231[k];

        t_127[k] = f_0 * hh_232[k];

        t_128[k] = f_0 * hh_233[k];

        t_129[k] = f_0 * hh_234[k];

        t_130[k] = f_0 * hh_235[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, t_136, t_137, t_138, hh_236, \
                         hh_237, hh_238, hh_239, hh_240, hh_241, hh_242, \
                         hh_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_0 * hh_236[k];

        t_132[k] = f_0 * hh_237[k];

        t_133[k] = f_0 * hh_238[k];

        t_134[k] = f_0 * hh_239[k];

        t_135[k] = f_0 * hh_240[k];

        t_136[k] = f_0 * hh_241[k];

        t_137[k] = f_0 * hh_242[k];

        t_138[k] = f_0 * hh_243[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, t_146, hh_244, \
                         hh_245, hh_246, hh_247, hh_248, hh_249, hh_250, \
                         hh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_0 * hh_244[k];

        t_140[k] = f_0 * hh_245[k];

        t_141[k] = f_0 * hh_246[k];

        t_142[k] = f_0 * hh_247[k];

        t_143[k] = f_0 * hh_248[k];

        t_144[k] = f_0 * hh_249[k];

        t_145[k] = f_0 * hh_250[k];

        t_146[k] = f_0 * hh_251[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, fh_63, fh_64, fh_65, fh_66, fh_67, \
                         hh_252, hh_253, hh_254, hh_255, hh_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -fh_63[k]
                   + f_0 * hh_252[k];

        t_148[k] = -fh_64[k]
                   + f_0 * hh_253[k];

        t_149[k] = -fh_65[k]
                   + f_0 * hh_254[k];

        t_150[k] = -fh_66[k]
                   + f_0 * hh_255[k];

        t_151[k] = -fh_67[k]
                   + f_0 * hh_256[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, fh_68, fh_69, fh_70, fh_71, fh_72, \
                         hh_257, hh_258, hh_259, hh_260, hh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -fh_68[k]
                   + f_0 * hh_257[k];

        t_153[k] = -fh_69[k]
                   + f_0 * hh_258[k];

        t_154[k] = -fh_70[k]
                   + f_0 * hh_259[k];

        t_155[k] = -fh_71[k]
                   + f_0 * hh_260[k];

        t_156[k] = -fh_72[k]
                   + f_0 * hh_261[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, fh_73, fh_74, fh_75, fh_76, fh_77, \
                         hh_262, hh_263, hh_264, hh_265, hh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -fh_73[k]
                   + f_0 * hh_262[k];

        t_158[k] = -fh_74[k]
                   + f_0 * hh_263[k];

        t_159[k] = -fh_75[k]
                   + f_0 * hh_264[k];

        t_160[k] = -fh_76[k]
                   + f_0 * hh_265[k];

        t_161[k] = -fh_77[k]
                   + f_0 * hh_266[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, fh_78, fh_79, fh_80, fh_81, fh_82, \
                         hh_267, hh_268, hh_269, hh_270, hh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -fh_78[k]
                   + f_0 * hh_267[k];

        t_163[k] = -fh_79[k]
                   + f_0 * hh_268[k];

        t_164[k] = -fh_80[k]
                   + f_0 * hh_269[k];

        t_165[k] = -fh_81[k]
                   + f_0 * hh_270[k];

        t_166[k] = -fh_82[k]
                   + f_0 * hh_271[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, fh_83, fh_84, fh_85, fh_86, fh_87, \
                         hh_272, hh_273, hh_274, hh_275, hh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -fh_83[k]
                   + f_0 * hh_272[k];

        t_168[k] = -2.0 * fh_84[k]
                   + f_0 * hh_273[k];

        t_169[k] = -2.0 * fh_85[k]
                   + f_0 * hh_274[k];

        t_170[k] = -2.0 * fh_86[k]
                   + f_0 * hh_275[k];

        t_171[k] = -2.0 * fh_87[k]
                   + f_0 * hh_276[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, fh_88, fh_89, fh_90, fh_91, fh_92, \
                         hh_277, hh_278, hh_279, hh_280, hh_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -2.0 * fh_88[k]
                   + f_0 * hh_277[k];

        t_173[k] = -2.0 * fh_89[k]
                   + f_0 * hh_278[k];

        t_174[k] = -2.0 * fh_90[k]
                   + f_0 * hh_279[k];

        t_175[k] = -2.0 * fh_91[k]
                   + f_0 * hh_280[k];

        t_176[k] = -2.0 * fh_92[k]
                   + f_0 * hh_281[k];
    }
}

static auto
compute_prim_geom_10_gh_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t fh, const size_t hh,
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

    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_104 = buffer.data(fh + 104);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_108 = buffer.data(fh + 108);
    const auto *fh_109 = buffer.data(fh + 109);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_111 = buffer.data(fh + 111);
    const auto *fh_112 = buffer.data(fh + 112);
    const auto *fh_113 = buffer.data(fh + 113);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_124 = buffer.data(fh + 124);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_128 = buffer.data(fh + 128);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_130 = buffer.data(fh + 130);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_133 = buffer.data(fh + 133);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_137 = buffer.data(fh + 137);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_147 = buffer.data(fh + 147);
    const auto *fh_148 = buffer.data(fh + 148);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_150 = buffer.data(fh + 150);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_153 = buffer.data(fh + 153);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_157 = buffer.data(fh + 157);
    const auto *fh_158 = buffer.data(fh + 158);
    const auto *fh_159 = buffer.data(fh + 159);
    const auto *fh_160 = buffer.data(fh + 160);
    const auto *fh_161 = buffer.data(fh + 161);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_168 = buffer.data(fh + 168);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_170 = buffer.data(fh + 170);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_173 = buffer.data(fh + 173);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_177 = buffer.data(fh + 177);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_182 = buffer.data(fh + 182);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_190 = buffer.data(fh + 190);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_193 = buffer.data(fh + 193);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_197 = buffer.data(fh + 197);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_202 = buffer.data(fh + 202);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_283 = buffer.data(hh + 283);
    const auto *hh_284 = buffer.data(hh + 284);
    const auto *hh_285 = buffer.data(hh + 285);
    const auto *hh_286 = buffer.data(hh + 286);
    const auto *hh_287 = buffer.data(hh + 287);
    const auto *hh_288 = buffer.data(hh + 288);
    const auto *hh_289 = buffer.data(hh + 289);
    const auto *hh_290 = buffer.data(hh + 290);
    const auto *hh_291 = buffer.data(hh + 291);
    const auto *hh_292 = buffer.data(hh + 292);
    const auto *hh_293 = buffer.data(hh + 293);
    const auto *hh_294 = buffer.data(hh + 294);
    const auto *hh_295 = buffer.data(hh + 295);
    const auto *hh_296 = buffer.data(hh + 296);
    const auto *hh_297 = buffer.data(hh + 297);
    const auto *hh_298 = buffer.data(hh + 298);
    const auto *hh_299 = buffer.data(hh + 299);
    const auto *hh_300 = buffer.data(hh + 300);
    const auto *hh_301 = buffer.data(hh + 301);
    const auto *hh_302 = buffer.data(hh + 302);
    const auto *hh_303 = buffer.data(hh + 303);
    const auto *hh_304 = buffer.data(hh + 304);
    const auto *hh_305 = buffer.data(hh + 305);
    const auto *hh_306 = buffer.data(hh + 306);
    const auto *hh_307 = buffer.data(hh + 307);
    const auto *hh_308 = buffer.data(hh + 308);
    const auto *hh_309 = buffer.data(hh + 309);
    const auto *hh_310 = buffer.data(hh + 310);
    const auto *hh_311 = buffer.data(hh + 311);
    const auto *hh_312 = buffer.data(hh + 312);
    const auto *hh_313 = buffer.data(hh + 313);
    const auto *hh_314 = buffer.data(hh + 314);
    const auto *hh_336 = buffer.data(hh + 336);
    const auto *hh_337 = buffer.data(hh + 337);
    const auto *hh_338 = buffer.data(hh + 338);
    const auto *hh_339 = buffer.data(hh + 339);
    const auto *hh_340 = buffer.data(hh + 340);
    const auto *hh_341 = buffer.data(hh + 341);
    const auto *hh_342 = buffer.data(hh + 342);
    const auto *hh_343 = buffer.data(hh + 343);
    const auto *hh_344 = buffer.data(hh + 344);
    const auto *hh_345 = buffer.data(hh + 345);
    const auto *hh_346 = buffer.data(hh + 346);
    const auto *hh_347 = buffer.data(hh + 347);
    const auto *hh_348 = buffer.data(hh + 348);
    const auto *hh_349 = buffer.data(hh + 349);
    const auto *hh_350 = buffer.data(hh + 350);
    const auto *hh_351 = buffer.data(hh + 351);
    const auto *hh_352 = buffer.data(hh + 352);
    const auto *hh_353 = buffer.data(hh + 353);
    const auto *hh_354 = buffer.data(hh + 354);
    const auto *hh_355 = buffer.data(hh + 355);
    const auto *hh_356 = buffer.data(hh + 356);
    const auto *hh_357 = buffer.data(hh + 357);
    const auto *hh_358 = buffer.data(hh + 358);
    const auto *hh_359 = buffer.data(hh + 359);
    const auto *hh_360 = buffer.data(hh + 360);
    const auto *hh_361 = buffer.data(hh + 361);
    const auto *hh_362 = buffer.data(hh + 362);
    const auto *hh_363 = buffer.data(hh + 363);
    const auto *hh_364 = buffer.data(hh + 364);
    const auto *hh_365 = buffer.data(hh + 365);
    const auto *hh_366 = buffer.data(hh + 366);
    const auto *hh_367 = buffer.data(hh + 367);
    const auto *hh_368 = buffer.data(hh + 368);
    const auto *hh_369 = buffer.data(hh + 369);
    const auto *hh_370 = buffer.data(hh + 370);
    const auto *hh_371 = buffer.data(hh + 371);
    const auto *hh_372 = buffer.data(hh + 372);
    const auto *hh_373 = buffer.data(hh + 373);
    const auto *hh_374 = buffer.data(hh + 374);
    const auto *hh_375 = buffer.data(hh + 375);
    const auto *hh_376 = buffer.data(hh + 376);
    const auto *hh_377 = buffer.data(hh + 377);
    const auto *hh_378 = buffer.data(hh + 378);
    const auto *hh_379 = buffer.data(hh + 379);
    const auto *hh_380 = buffer.data(hh + 380);
    const auto *hh_381 = buffer.data(hh + 381);
    const auto *hh_382 = buffer.data(hh + 382);
    const auto *hh_383 = buffer.data(hh + 383);
    const auto *hh_384 = buffer.data(hh + 384);
    const auto *hh_385 = buffer.data(hh + 385);
    const auto *hh_386 = buffer.data(hh + 386);
    const auto *hh_387 = buffer.data(hh + 387);
    const auto *hh_388 = buffer.data(hh + 388);
    const auto *hh_389 = buffer.data(hh + 389);
    const auto *hh_390 = buffer.data(hh + 390);
    const auto *hh_391 = buffer.data(hh + 391);
    const auto *hh_392 = buffer.data(hh + 392);
    const auto *hh_393 = buffer.data(hh + 393);
    const auto *hh_394 = buffer.data(hh + 394);
    const auto *hh_395 = buffer.data(hh + 395);
    const auto *hh_396 = buffer.data(hh + 396);
    const auto *hh_397 = buffer.data(hh + 397);
    const auto *hh_398 = buffer.data(hh + 398);
    const auto *hh_399 = buffer.data(hh + 399);
    const auto *hh_400 = buffer.data(hh + 400);
    const auto *hh_401 = buffer.data(hh + 401);
    const auto *hh_402 = buffer.data(hh + 402);
    const auto *hh_403 = buffer.data(hh + 403);
    const auto *hh_404 = buffer.data(hh + 404);
    const auto *hh_405 = buffer.data(hh + 405);
    const auto *hh_406 = buffer.data(hh + 406);
    const auto *hh_407 = buffer.data(hh + 407);
    const auto *hh_408 = buffer.data(hh + 408);
    const auto *hh_409 = buffer.data(hh + 409);
    const auto *hh_410 = buffer.data(hh + 410);
    const auto *hh_411 = buffer.data(hh + 411);
    const auto *hh_412 = buffer.data(hh + 412);
    const auto *hh_413 = buffer.data(hh + 413);
    const auto *hh_414 = buffer.data(hh + 414);
    const auto *hh_415 = buffer.data(hh + 415);
    const auto *hh_416 = buffer.data(hh + 416);
    const auto *hh_417 = buffer.data(hh + 417);
    const auto *hh_418 = buffer.data(hh + 418);
    const auto *hh_419 = buffer.data(hh + 419);
    const auto *hh_420 = buffer.data(hh + 420);
    const auto *hh_421 = buffer.data(hh + 421);
    const auto *hh_422 = buffer.data(hh + 422);
    const auto *hh_423 = buffer.data(hh + 423);
    const auto *hh_424 = buffer.data(hh + 424);
    const auto *hh_425 = buffer.data(hh + 425);
    const auto *hh_426 = buffer.data(hh + 426);
    const auto *hh_427 = buffer.data(hh + 427);
    const auto *hh_428 = buffer.data(hh + 428);
    const auto *hh_429 = buffer.data(hh + 429);
    const auto *hh_430 = buffer.data(hh + 430);
    const auto *hh_431 = buffer.data(hh + 431);
    const auto *hh_432 = buffer.data(hh + 432);
    const auto *hh_433 = buffer.data(hh + 433);
    const auto *hh_434 = buffer.data(hh + 434);
    const auto *hh_435 = buffer.data(hh + 435);
    const auto *hh_436 = buffer.data(hh + 436);
    const auto *hh_437 = buffer.data(hh + 437);
    const auto *hh_438 = buffer.data(hh + 438);
    const auto *hh_439 = buffer.data(hh + 439);
    const auto *hh_440 = buffer.data(hh + 440);

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, fh_93, fh_94, fh_95, fh_96, fh_97, \
                         hh_282, hh_283, hh_284, hh_285, hh_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -2.0 * fh_93[k]
                   + f_0 * hh_282[k];

        t_178[k] = -2.0 * fh_94[k]
                   + f_0 * hh_283[k];

        t_179[k] = -2.0 * fh_95[k]
                   + f_0 * hh_284[k];

        t_180[k] = -2.0 * fh_96[k]
                   + f_0 * hh_285[k];

        t_181[k] = -2.0 * fh_97[k]
                   + f_0 * hh_286[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, fh_98, fh_99, fh_100, fh_101, \
                         fh_102, hh_287, hh_288, hh_289, hh_290, \
                         hh_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * fh_98[k]
                   + f_0 * hh_287[k];

        t_183[k] = -2.0 * fh_99[k]
                   + f_0 * hh_288[k];

        t_184[k] = -2.0 * fh_100[k]
                   + f_0 * hh_289[k];

        t_185[k] = -2.0 * fh_101[k]
                   + f_0 * hh_290[k];

        t_186[k] = -2.0 * fh_102[k]
                   + f_0 * hh_291[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, fh_103, fh_104, fh_105, fh_106, \
                         fh_107, hh_292, hh_293, hh_294, hh_295, \
                         hh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * fh_103[k]
                   + f_0 * hh_292[k];

        t_188[k] = -2.0 * fh_104[k]
                   + f_0 * hh_293[k];

        t_189[k] = -3.0 * fh_105[k]
                   + f_0 * hh_294[k];

        t_190[k] = -3.0 * fh_106[k]
                   + f_0 * hh_295[k];

        t_191[k] = -3.0 * fh_107[k]
                   + f_0 * hh_296[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, fh_108, fh_109, fh_110, fh_111, \
                         fh_112, hh_297, hh_298, hh_299, hh_300, \
                         hh_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -3.0 * fh_108[k]
                   + f_0 * hh_297[k];

        t_193[k] = -3.0 * fh_109[k]
                   + f_0 * hh_298[k];

        t_194[k] = -3.0 * fh_110[k]
                   + f_0 * hh_299[k];

        t_195[k] = -3.0 * fh_111[k]
                   + f_0 * hh_300[k];

        t_196[k] = -3.0 * fh_112[k]
                   + f_0 * hh_301[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, fh_113, fh_114, fh_115, fh_116, \
                         fh_117, hh_302, hh_303, hh_304, hh_305, \
                         hh_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -3.0 * fh_113[k]
                   + f_0 * hh_302[k];

        t_198[k] = -3.0 * fh_114[k]
                   + f_0 * hh_303[k];

        t_199[k] = -3.0 * fh_115[k]
                   + f_0 * hh_304[k];

        t_200[k] = -3.0 * fh_116[k]
                   + f_0 * hh_305[k];

        t_201[k] = -3.0 * fh_117[k]
                   + f_0 * hh_306[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, fh_118, fh_119, fh_120, fh_121, \
                         fh_122, hh_307, hh_308, hh_309, hh_310, \
                         hh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -3.0 * fh_118[k]
                   + f_0 * hh_307[k];

        t_203[k] = -3.0 * fh_119[k]
                   + f_0 * hh_308[k];

        t_204[k] = -3.0 * fh_120[k]
                   + f_0 * hh_309[k];

        t_205[k] = -3.0 * fh_121[k]
                   + f_0 * hh_310[k];

        t_206[k] = -3.0 * fh_122[k]
                   + f_0 * hh_311[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, fh_123, fh_124, fh_125, \
                         hh_312, hh_313, hh_314, hh_336, hh_337, \
                         hh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -3.0 * fh_123[k]
                   + f_0 * hh_312[k];

        t_208[k] = -3.0 * fh_124[k]
                   + f_0 * hh_313[k];

        t_209[k] = -3.0 * fh_125[k]
                   + f_0 * hh_314[k];

        t_210[k] = f_0 * hh_336[k];

        t_211[k] = f_0 * hh_337[k];

        t_212[k] = f_0 * hh_338[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, t_218, t_219, t_220, hh_339, \
                         hh_340, hh_341, hh_342, hh_343, hh_344, hh_345, \
                         hh_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_0 * hh_339[k];

        t_214[k] = f_0 * hh_340[k];

        t_215[k] = f_0 * hh_341[k];

        t_216[k] = f_0 * hh_342[k];

        t_217[k] = f_0 * hh_343[k];

        t_218[k] = f_0 * hh_344[k];

        t_219[k] = f_0 * hh_345[k];

        t_220[k] = f_0 * hh_346[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, t_226, t_227, t_228, hh_347, \
                         hh_348, hh_349, hh_350, hh_351, hh_352, hh_353, \
                         hh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_0 * hh_347[k];

        t_222[k] = f_0 * hh_348[k];

        t_223[k] = f_0 * hh_349[k];

        t_224[k] = f_0 * hh_350[k];

        t_225[k] = f_0 * hh_351[k];

        t_226[k] = f_0 * hh_352[k];

        t_227[k] = f_0 * hh_353[k];

        t_228[k] = f_0 * hh_354[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, t_234, fh_126, fh_127, fh_128, \
                         fh_129, hh_355, hh_356, hh_357, hh_358, hh_359, \
                         hh_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * hh_355[k];

        t_230[k] = f_0 * hh_356[k];

        t_231[k] = -fh_126[k]
                   + f_0 * hh_357[k];

        t_232[k] = -fh_127[k]
                   + f_0 * hh_358[k];

        t_233[k] = -fh_128[k]
                   + f_0 * hh_359[k];

        t_234[k] = -fh_129[k]
                   + f_0 * hh_360[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, fh_130, fh_131, fh_132, fh_133, \
                         fh_134, hh_361, hh_362, hh_363, hh_364, \
                         hh_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -fh_130[k]
                   + f_0 * hh_361[k];

        t_236[k] = -fh_131[k]
                   + f_0 * hh_362[k];

        t_237[k] = -fh_132[k]
                   + f_0 * hh_363[k];

        t_238[k] = -fh_133[k]
                   + f_0 * hh_364[k];

        t_239[k] = -fh_134[k]
                   + f_0 * hh_365[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, fh_135, fh_136, fh_137, fh_138, \
                         fh_139, hh_366, hh_367, hh_368, hh_369, \
                         hh_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -fh_135[k]
                   + f_0 * hh_366[k];

        t_241[k] = -fh_136[k]
                   + f_0 * hh_367[k];

        t_242[k] = -fh_137[k]
                   + f_0 * hh_368[k];

        t_243[k] = -fh_138[k]
                   + f_0 * hh_369[k];

        t_244[k] = -fh_139[k]
                   + f_0 * hh_370[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, fh_140, fh_141, fh_142, fh_143, \
                         fh_144, hh_371, hh_372, hh_373, hh_374, \
                         hh_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -fh_140[k]
                   + f_0 * hh_371[k];

        t_246[k] = -fh_141[k]
                   + f_0 * hh_372[k];

        t_247[k] = -fh_142[k]
                   + f_0 * hh_373[k];

        t_248[k] = -fh_143[k]
                   + f_0 * hh_374[k];

        t_249[k] = -fh_144[k]
                   + f_0 * hh_375[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, fh_145, fh_146, fh_147, fh_148, \
                         fh_149, hh_376, hh_377, hh_378, hh_379, \
                         hh_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -fh_145[k]
                   + f_0 * hh_376[k];

        t_251[k] = -fh_146[k]
                   + f_0 * hh_377[k];

        t_252[k] = -2.0 * fh_147[k]
                   + f_0 * hh_378[k];

        t_253[k] = -2.0 * fh_148[k]
                   + f_0 * hh_379[k];

        t_254[k] = -2.0 * fh_149[k]
                   + f_0 * hh_380[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, fh_150, fh_151, fh_152, fh_153, \
                         fh_154, hh_381, hh_382, hh_383, hh_384, \
                         hh_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * fh_150[k]
                   + f_0 * hh_381[k];

        t_256[k] = -2.0 * fh_151[k]
                   + f_0 * hh_382[k];

        t_257[k] = -2.0 * fh_152[k]
                   + f_0 * hh_383[k];

        t_258[k] = -2.0 * fh_153[k]
                   + f_0 * hh_384[k];

        t_259[k] = -2.0 * fh_154[k]
                   + f_0 * hh_385[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, fh_155, fh_156, fh_157, fh_158, \
                         fh_159, hh_386, hh_387, hh_388, hh_389, \
                         hh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * fh_155[k]
                   + f_0 * hh_386[k];

        t_261[k] = -2.0 * fh_156[k]
                   + f_0 * hh_387[k];

        t_262[k] = -2.0 * fh_157[k]
                   + f_0 * hh_388[k];

        t_263[k] = -2.0 * fh_158[k]
                   + f_0 * hh_389[k];

        t_264[k] = -2.0 * fh_159[k]
                   + f_0 * hh_390[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, fh_160, fh_161, fh_162, fh_163, \
                         fh_164, hh_391, hh_392, hh_393, hh_394, \
                         hh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * fh_160[k]
                   + f_0 * hh_391[k];

        t_266[k] = -2.0 * fh_161[k]
                   + f_0 * hh_392[k];

        t_267[k] = -2.0 * fh_162[k]
                   + f_0 * hh_393[k];

        t_268[k] = -2.0 * fh_163[k]
                   + f_0 * hh_394[k];

        t_269[k] = -2.0 * fh_164[k]
                   + f_0 * hh_395[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, fh_165, fh_166, fh_167, fh_168, \
                         fh_169, hh_396, hh_397, hh_398, hh_399, \
                         hh_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * fh_165[k]
                   + f_0 * hh_396[k];

        t_271[k] = -2.0 * fh_166[k]
                   + f_0 * hh_397[k];

        t_272[k] = -2.0 * fh_167[k]
                   + f_0 * hh_398[k];

        t_273[k] = -3.0 * fh_168[k]
                   + f_0 * hh_399[k];

        t_274[k] = -3.0 * fh_169[k]
                   + f_0 * hh_400[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, fh_170, fh_171, fh_172, fh_173, \
                         fh_174, hh_401, hh_402, hh_403, hh_404, \
                         hh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * fh_170[k]
                   + f_0 * hh_401[k];

        t_276[k] = -3.0 * fh_171[k]
                   + f_0 * hh_402[k];

        t_277[k] = -3.0 * fh_172[k]
                   + f_0 * hh_403[k];

        t_278[k] = -3.0 * fh_173[k]
                   + f_0 * hh_404[k];

        t_279[k] = -3.0 * fh_174[k]
                   + f_0 * hh_405[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, fh_175, fh_176, fh_177, fh_178, \
                         fh_179, hh_406, hh_407, hh_408, hh_409, \
                         hh_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * fh_175[k]
                   + f_0 * hh_406[k];

        t_281[k] = -3.0 * fh_176[k]
                   + f_0 * hh_407[k];

        t_282[k] = -3.0 * fh_177[k]
                   + f_0 * hh_408[k];

        t_283[k] = -3.0 * fh_178[k]
                   + f_0 * hh_409[k];

        t_284[k] = -3.0 * fh_179[k]
                   + f_0 * hh_410[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, fh_180, fh_181, fh_182, fh_183, \
                         fh_184, hh_411, hh_412, hh_413, hh_414, \
                         hh_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * fh_180[k]
                   + f_0 * hh_411[k];

        t_286[k] = -3.0 * fh_181[k]
                   + f_0 * hh_412[k];

        t_287[k] = -3.0 * fh_182[k]
                   + f_0 * hh_413[k];

        t_288[k] = -3.0 * fh_183[k]
                   + f_0 * hh_414[k];

        t_289[k] = -3.0 * fh_184[k]
                   + f_0 * hh_415[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, fh_185, fh_186, fh_187, fh_188, \
                         fh_189, hh_416, hh_417, hh_418, hh_419, \
                         hh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * fh_185[k]
                   + f_0 * hh_416[k];

        t_291[k] = -3.0 * fh_186[k]
                   + f_0 * hh_417[k];

        t_292[k] = -3.0 * fh_187[k]
                   + f_0 * hh_418[k];

        t_293[k] = -3.0 * fh_188[k]
                   + f_0 * hh_419[k];

        t_294[k] = -4.0 * fh_189[k]
                   + f_0 * hh_420[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, fh_190, fh_191, fh_192, fh_193, \
                         fh_194, hh_421, hh_422, hh_423, hh_424, \
                         hh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * fh_190[k]
                   + f_0 * hh_421[k];

        t_296[k] = -4.0 * fh_191[k]
                   + f_0 * hh_422[k];

        t_297[k] = -4.0 * fh_192[k]
                   + f_0 * hh_423[k];

        t_298[k] = -4.0 * fh_193[k]
                   + f_0 * hh_424[k];

        t_299[k] = -4.0 * fh_194[k]
                   + f_0 * hh_425[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, fh_195, fh_196, fh_197, fh_198, \
                         fh_199, hh_426, hh_427, hh_428, hh_429, \
                         hh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * fh_195[k]
                   + f_0 * hh_426[k];

        t_301[k] = -4.0 * fh_196[k]
                   + f_0 * hh_427[k];

        t_302[k] = -4.0 * fh_197[k]
                   + f_0 * hh_428[k];

        t_303[k] = -4.0 * fh_198[k]
                   + f_0 * hh_429[k];

        t_304[k] = -4.0 * fh_199[k]
                   + f_0 * hh_430[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, fh_200, fh_201, fh_202, fh_203, \
                         fh_204, hh_431, hh_432, hh_433, hh_434, \
                         hh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * fh_200[k]
                   + f_0 * hh_431[k];

        t_306[k] = -4.0 * fh_201[k]
                   + f_0 * hh_432[k];

        t_307[k] = -4.0 * fh_202[k]
                   + f_0 * hh_433[k];

        t_308[k] = -4.0 * fh_203[k]
                   + f_0 * hh_434[k];

        t_309[k] = -4.0 * fh_204[k]
                   + f_0 * hh_435[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, fh_205, fh_206, fh_207, fh_208, \
                         fh_209, hh_436, hh_437, hh_438, hh_439, \
                         hh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -4.0 * fh_205[k]
                   + f_0 * hh_436[k];

        t_311[k] = -4.0 * fh_206[k]
                   + f_0 * hh_437[k];

        t_312[k] = -4.0 * fh_207[k]
                   + f_0 * hh_438[k];

        t_313[k] = -4.0 * fh_208[k]
                   + f_0 * hh_439[k];

        t_314[k] = -4.0 * fh_209[k]
                   + f_0 * hh_440[k];
    }
}

auto
compute_prim_geom_10_gh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t fh, const size_t hh,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_gh_electron_repulsion_2_piece0(buffer, target, fh, hh, ncols, alpha);

    compute_prim_geom_10_gh_electron_repulsion_2_piece1(buffer, target, fh, hh, ncols, alpha);
}

}  // namespace simdt2ceri
