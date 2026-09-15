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


#include "SimdElectronRepulsionGeom10VrrRecFL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_fl_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dl, const size_t gl,
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

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);
    const auto *dl_3 = buffer.data(dl + 3);
    const auto *dl_4 = buffer.data(dl + 4);
    const auto *dl_5 = buffer.data(dl + 5);
    const auto *dl_6 = buffer.data(dl + 6);
    const auto *dl_7 = buffer.data(dl + 7);
    const auto *dl_8 = buffer.data(dl + 8);
    const auto *dl_9 = buffer.data(dl + 9);
    const auto *dl_10 = buffer.data(dl + 10);
    const auto *dl_11 = buffer.data(dl + 11);
    const auto *dl_12 = buffer.data(dl + 12);
    const auto *dl_13 = buffer.data(dl + 13);
    const auto *dl_14 = buffer.data(dl + 14);
    const auto *dl_15 = buffer.data(dl + 15);
    const auto *dl_16 = buffer.data(dl + 16);
    const auto *dl_17 = buffer.data(dl + 17);
    const auto *dl_18 = buffer.data(dl + 18);
    const auto *dl_19 = buffer.data(dl + 19);
    const auto *dl_20 = buffer.data(dl + 20);
    const auto *dl_21 = buffer.data(dl + 21);
    const auto *dl_22 = buffer.data(dl + 22);
    const auto *dl_23 = buffer.data(dl + 23);
    const auto *dl_24 = buffer.data(dl + 24);
    const auto *dl_25 = buffer.data(dl + 25);
    const auto *dl_26 = buffer.data(dl + 26);
    const auto *dl_27 = buffer.data(dl + 27);
    const auto *dl_28 = buffer.data(dl + 28);
    const auto *dl_29 = buffer.data(dl + 29);
    const auto *dl_30 = buffer.data(dl + 30);
    const auto *dl_31 = buffer.data(dl + 31);
    const auto *dl_32 = buffer.data(dl + 32);
    const auto *dl_33 = buffer.data(dl + 33);
    const auto *dl_34 = buffer.data(dl + 34);
    const auto *dl_35 = buffer.data(dl + 35);
    const auto *dl_36 = buffer.data(dl + 36);
    const auto *dl_37 = buffer.data(dl + 37);
    const auto *dl_38 = buffer.data(dl + 38);
    const auto *dl_39 = buffer.data(dl + 39);
    const auto *dl_40 = buffer.data(dl + 40);
    const auto *dl_41 = buffer.data(dl + 41);
    const auto *dl_42 = buffer.data(dl + 42);
    const auto *dl_43 = buffer.data(dl + 43);
    const auto *dl_44 = buffer.data(dl + 44);
    const auto *dl_45 = buffer.data(dl + 45);
    const auto *dl_46 = buffer.data(dl + 46);
    const auto *dl_47 = buffer.data(dl + 47);
    const auto *dl_48 = buffer.data(dl + 48);
    const auto *dl_49 = buffer.data(dl + 49);
    const auto *dl_50 = buffer.data(dl + 50);
    const auto *dl_51 = buffer.data(dl + 51);
    const auto *dl_52 = buffer.data(dl + 52);
    const auto *dl_53 = buffer.data(dl + 53);
    const auto *dl_54 = buffer.data(dl + 54);
    const auto *dl_55 = buffer.data(dl + 55);
    const auto *dl_56 = buffer.data(dl + 56);
    const auto *dl_57 = buffer.data(dl + 57);
    const auto *dl_58 = buffer.data(dl + 58);
    const auto *dl_59 = buffer.data(dl + 59);
    const auto *dl_60 = buffer.data(dl + 60);
    const auto *dl_61 = buffer.data(dl + 61);
    const auto *dl_62 = buffer.data(dl + 62);
    const auto *dl_63 = buffer.data(dl + 63);
    const auto *dl_64 = buffer.data(dl + 64);
    const auto *dl_65 = buffer.data(dl + 65);
    const auto *dl_66 = buffer.data(dl + 66);
    const auto *dl_67 = buffer.data(dl + 67);
    const auto *dl_68 = buffer.data(dl + 68);
    const auto *dl_69 = buffer.data(dl + 69);
    const auto *dl_70 = buffer.data(dl + 70);
    const auto *dl_71 = buffer.data(dl + 71);
    const auto *dl_72 = buffer.data(dl + 72);
    const auto *dl_73 = buffer.data(dl + 73);
    const auto *dl_74 = buffer.data(dl + 74);
    const auto *dl_75 = buffer.data(dl + 75);
    const auto *dl_76 = buffer.data(dl + 76);
    const auto *dl_77 = buffer.data(dl + 77);
    const auto *dl_78 = buffer.data(dl + 78);
    const auto *dl_79 = buffer.data(dl + 79);
    const auto *dl_80 = buffer.data(dl + 80);
    const auto *dl_81 = buffer.data(dl + 81);
    const auto *dl_82 = buffer.data(dl + 82);
    const auto *dl_83 = buffer.data(dl + 83);
    const auto *dl_84 = buffer.data(dl + 84);
    const auto *dl_85 = buffer.data(dl + 85);
    const auto *dl_86 = buffer.data(dl + 86);
    const auto *dl_87 = buffer.data(dl + 87);
    const auto *dl_88 = buffer.data(dl + 88);
    const auto *dl_89 = buffer.data(dl + 89);
    const auto *dl_90 = buffer.data(dl + 90);
    const auto *dl_91 = buffer.data(dl + 91);
    const auto *dl_92 = buffer.data(dl + 92);
    const auto *dl_93 = buffer.data(dl + 93);
    const auto *dl_94 = buffer.data(dl + 94);
    const auto *dl_95 = buffer.data(dl + 95);
    const auto *dl_96 = buffer.data(dl + 96);
    const auto *dl_97 = buffer.data(dl + 97);
    const auto *dl_98 = buffer.data(dl + 98);
    const auto *dl_99 = buffer.data(dl + 99);
    const auto *dl_100 = buffer.data(dl + 100);
    const auto *dl_101 = buffer.data(dl + 101);
    const auto *dl_102 = buffer.data(dl + 102);
    const auto *dl_103 = buffer.data(dl + 103);
    const auto *dl_104 = buffer.data(dl + 104);
    const auto *dl_105 = buffer.data(dl + 105);
    const auto *dl_106 = buffer.data(dl + 106);
    const auto *dl_107 = buffer.data(dl + 107);
    const auto *dl_108 = buffer.data(dl + 108);
    const auto *dl_109 = buffer.data(dl + 109);
    const auto *dl_110 = buffer.data(dl + 110);
    const auto *dl_111 = buffer.data(dl + 111);
    const auto *dl_112 = buffer.data(dl + 112);
    const auto *dl_113 = buffer.data(dl + 113);
    const auto *dl_114 = buffer.data(dl + 114);
    const auto *dl_115 = buffer.data(dl + 115);
    const auto *dl_116 = buffer.data(dl + 116);
    const auto *dl_117 = buffer.data(dl + 117);
    const auto *dl_118 = buffer.data(dl + 118);
    const auto *dl_119 = buffer.data(dl + 119);
    const auto *dl_120 = buffer.data(dl + 120);
    const auto *dl_121 = buffer.data(dl + 121);
    const auto *dl_122 = buffer.data(dl + 122);
    const auto *dl_123 = buffer.data(dl + 123);
    const auto *dl_124 = buffer.data(dl + 124);
    const auto *dl_125 = buffer.data(dl + 125);
    const auto *dl_126 = buffer.data(dl + 126);
    const auto *dl_127 = buffer.data(dl + 127);
    const auto *dl_128 = buffer.data(dl + 128);
    const auto *dl_129 = buffer.data(dl + 129);
    const auto *dl_130 = buffer.data(dl + 130);
    const auto *dl_131 = buffer.data(dl + 131);
    const auto *dl_132 = buffer.data(dl + 132);
    const auto *dl_133 = buffer.data(dl + 133);
    const auto *dl_134 = buffer.data(dl + 134);
    const auto *dl_135 = buffer.data(dl + 135);
    const auto *dl_136 = buffer.data(dl + 136);
    const auto *dl_137 = buffer.data(dl + 137);
    const auto *dl_138 = buffer.data(dl + 138);
    const auto *dl_139 = buffer.data(dl + 139);
    const auto *dl_140 = buffer.data(dl + 140);
    const auto *dl_141 = buffer.data(dl + 141);
    const auto *dl_142 = buffer.data(dl + 142);
    const auto *dl_143 = buffer.data(dl + 143);
    const auto *dl_144 = buffer.data(dl + 144);
    const auto *dl_145 = buffer.data(dl + 145);
    const auto *dl_146 = buffer.data(dl + 146);
    const auto *dl_147 = buffer.data(dl + 147);
    const auto *dl_148 = buffer.data(dl + 148);
    const auto *dl_149 = buffer.data(dl + 149);

    const auto *gl_0 = buffer.data(gl + 0);
    const auto *gl_1 = buffer.data(gl + 1);
    const auto *gl_2 = buffer.data(gl + 2);
    const auto *gl_3 = buffer.data(gl + 3);
    const auto *gl_4 = buffer.data(gl + 4);
    const auto *gl_5 = buffer.data(gl + 5);
    const auto *gl_6 = buffer.data(gl + 6);
    const auto *gl_7 = buffer.data(gl + 7);
    const auto *gl_8 = buffer.data(gl + 8);
    const auto *gl_9 = buffer.data(gl + 9);
    const auto *gl_10 = buffer.data(gl + 10);
    const auto *gl_11 = buffer.data(gl + 11);
    const auto *gl_12 = buffer.data(gl + 12);
    const auto *gl_13 = buffer.data(gl + 13);
    const auto *gl_14 = buffer.data(gl + 14);
    const auto *gl_15 = buffer.data(gl + 15);
    const auto *gl_16 = buffer.data(gl + 16);
    const auto *gl_17 = buffer.data(gl + 17);
    const auto *gl_18 = buffer.data(gl + 18);
    const auto *gl_19 = buffer.data(gl + 19);
    const auto *gl_20 = buffer.data(gl + 20);
    const auto *gl_21 = buffer.data(gl + 21);
    const auto *gl_22 = buffer.data(gl + 22);
    const auto *gl_23 = buffer.data(gl + 23);
    const auto *gl_24 = buffer.data(gl + 24);
    const auto *gl_25 = buffer.data(gl + 25);
    const auto *gl_26 = buffer.data(gl + 26);
    const auto *gl_27 = buffer.data(gl + 27);
    const auto *gl_28 = buffer.data(gl + 28);
    const auto *gl_29 = buffer.data(gl + 29);
    const auto *gl_30 = buffer.data(gl + 30);
    const auto *gl_31 = buffer.data(gl + 31);
    const auto *gl_32 = buffer.data(gl + 32);
    const auto *gl_33 = buffer.data(gl + 33);
    const auto *gl_34 = buffer.data(gl + 34);
    const auto *gl_35 = buffer.data(gl + 35);
    const auto *gl_36 = buffer.data(gl + 36);
    const auto *gl_37 = buffer.data(gl + 37);
    const auto *gl_38 = buffer.data(gl + 38);
    const auto *gl_39 = buffer.data(gl + 39);
    const auto *gl_40 = buffer.data(gl + 40);
    const auto *gl_41 = buffer.data(gl + 41);
    const auto *gl_42 = buffer.data(gl + 42);
    const auto *gl_43 = buffer.data(gl + 43);
    const auto *gl_44 = buffer.data(gl + 44);
    const auto *gl_45 = buffer.data(gl + 45);
    const auto *gl_46 = buffer.data(gl + 46);
    const auto *gl_47 = buffer.data(gl + 47);
    const auto *gl_48 = buffer.data(gl + 48);
    const auto *gl_49 = buffer.data(gl + 49);
    const auto *gl_50 = buffer.data(gl + 50);
    const auto *gl_51 = buffer.data(gl + 51);
    const auto *gl_52 = buffer.data(gl + 52);
    const auto *gl_53 = buffer.data(gl + 53);
    const auto *gl_54 = buffer.data(gl + 54);
    const auto *gl_55 = buffer.data(gl + 55);
    const auto *gl_56 = buffer.data(gl + 56);
    const auto *gl_57 = buffer.data(gl + 57);
    const auto *gl_58 = buffer.data(gl + 58);
    const auto *gl_59 = buffer.data(gl + 59);
    const auto *gl_60 = buffer.data(gl + 60);
    const auto *gl_61 = buffer.data(gl + 61);
    const auto *gl_62 = buffer.data(gl + 62);
    const auto *gl_63 = buffer.data(gl + 63);
    const auto *gl_64 = buffer.data(gl + 64);
    const auto *gl_65 = buffer.data(gl + 65);
    const auto *gl_66 = buffer.data(gl + 66);
    const auto *gl_67 = buffer.data(gl + 67);
    const auto *gl_68 = buffer.data(gl + 68);
    const auto *gl_69 = buffer.data(gl + 69);
    const auto *gl_70 = buffer.data(gl + 70);
    const auto *gl_71 = buffer.data(gl + 71);
    const auto *gl_72 = buffer.data(gl + 72);
    const auto *gl_73 = buffer.data(gl + 73);
    const auto *gl_74 = buffer.data(gl + 74);
    const auto *gl_75 = buffer.data(gl + 75);
    const auto *gl_76 = buffer.data(gl + 76);
    const auto *gl_77 = buffer.data(gl + 77);
    const auto *gl_78 = buffer.data(gl + 78);
    const auto *gl_79 = buffer.data(gl + 79);
    const auto *gl_80 = buffer.data(gl + 80);
    const auto *gl_81 = buffer.data(gl + 81);
    const auto *gl_82 = buffer.data(gl + 82);
    const auto *gl_83 = buffer.data(gl + 83);
    const auto *gl_84 = buffer.data(gl + 84);
    const auto *gl_85 = buffer.data(gl + 85);
    const auto *gl_86 = buffer.data(gl + 86);
    const auto *gl_87 = buffer.data(gl + 87);
    const auto *gl_88 = buffer.data(gl + 88);
    const auto *gl_89 = buffer.data(gl + 89);
    const auto *gl_90 = buffer.data(gl + 90);
    const auto *gl_91 = buffer.data(gl + 91);
    const auto *gl_92 = buffer.data(gl + 92);
    const auto *gl_93 = buffer.data(gl + 93);
    const auto *gl_94 = buffer.data(gl + 94);
    const auto *gl_95 = buffer.data(gl + 95);
    const auto *gl_96 = buffer.data(gl + 96);
    const auto *gl_97 = buffer.data(gl + 97);
    const auto *gl_98 = buffer.data(gl + 98);
    const auto *gl_99 = buffer.data(gl + 99);
    const auto *gl_100 = buffer.data(gl + 100);
    const auto *gl_101 = buffer.data(gl + 101);
    const auto *gl_102 = buffer.data(gl + 102);
    const auto *gl_103 = buffer.data(gl + 103);
    const auto *gl_104 = buffer.data(gl + 104);
    const auto *gl_105 = buffer.data(gl + 105);
    const auto *gl_106 = buffer.data(gl + 106);
    const auto *gl_107 = buffer.data(gl + 107);
    const auto *gl_108 = buffer.data(gl + 108);
    const auto *gl_109 = buffer.data(gl + 109);
    const auto *gl_110 = buffer.data(gl + 110);
    const auto *gl_111 = buffer.data(gl + 111);
    const auto *gl_112 = buffer.data(gl + 112);
    const auto *gl_113 = buffer.data(gl + 113);
    const auto *gl_114 = buffer.data(gl + 114);
    const auto *gl_115 = buffer.data(gl + 115);
    const auto *gl_116 = buffer.data(gl + 116);
    const auto *gl_117 = buffer.data(gl + 117);
    const auto *gl_118 = buffer.data(gl + 118);
    const auto *gl_119 = buffer.data(gl + 119);
    const auto *gl_120 = buffer.data(gl + 120);
    const auto *gl_121 = buffer.data(gl + 121);
    const auto *gl_122 = buffer.data(gl + 122);
    const auto *gl_123 = buffer.data(gl + 123);
    const auto *gl_124 = buffer.data(gl + 124);
    const auto *gl_125 = buffer.data(gl + 125);
    const auto *gl_126 = buffer.data(gl + 126);
    const auto *gl_127 = buffer.data(gl + 127);
    const auto *gl_128 = buffer.data(gl + 128);
    const auto *gl_129 = buffer.data(gl + 129);
    const auto *gl_130 = buffer.data(gl + 130);
    const auto *gl_131 = buffer.data(gl + 131);
    const auto *gl_132 = buffer.data(gl + 132);
    const auto *gl_133 = buffer.data(gl + 133);
    const auto *gl_134 = buffer.data(gl + 134);
    const auto *gl_135 = buffer.data(gl + 135);
    const auto *gl_136 = buffer.data(gl + 136);
    const auto *gl_137 = buffer.data(gl + 137);
    const auto *gl_138 = buffer.data(gl + 138);
    const auto *gl_139 = buffer.data(gl + 139);
    const auto *gl_140 = buffer.data(gl + 140);
    const auto *gl_141 = buffer.data(gl + 141);
    const auto *gl_142 = buffer.data(gl + 142);
    const auto *gl_143 = buffer.data(gl + 143);
    const auto *gl_144 = buffer.data(gl + 144);
    const auto *gl_145 = buffer.data(gl + 145);
    const auto *gl_146 = buffer.data(gl + 146);
    const auto *gl_147 = buffer.data(gl + 147);
    const auto *gl_148 = buffer.data(gl + 148);
    const auto *gl_149 = buffer.data(gl + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, dl_0, dl_1, dl_2, dl_3, dl_4, gl_0, gl_1, \
                         gl_2, gl_3, gl_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -3.0 * dl_0[k]
                 + f_0 * gl_0[k];

        t_1[k] = -3.0 * dl_1[k]
                 + f_0 * gl_1[k];

        t_2[k] = -3.0 * dl_2[k]
                 + f_0 * gl_2[k];

        t_3[k] = -3.0 * dl_3[k]
                 + f_0 * gl_3[k];

        t_4[k] = -3.0 * dl_4[k]
                 + f_0 * gl_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, dl_5, dl_6, dl_7, dl_8, dl_9, gl_5, gl_6, \
                         gl_7, gl_8, gl_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -3.0 * dl_5[k]
                 + f_0 * gl_5[k];

        t_6[k] = -3.0 * dl_6[k]
                 + f_0 * gl_6[k];

        t_7[k] = -3.0 * dl_7[k]
                 + f_0 * gl_7[k];

        t_8[k] = -3.0 * dl_8[k]
                 + f_0 * gl_8[k];

        t_9[k] = -3.0 * dl_9[k]
                 + f_0 * gl_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, dl_10, dl_11, dl_12, dl_13, dl_14, \
                         gl_10, gl_11, gl_12, gl_13, gl_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -3.0 * dl_10[k]
                  + f_0 * gl_10[k];

        t_11[k] = -3.0 * dl_11[k]
                  + f_0 * gl_11[k];

        t_12[k] = -3.0 * dl_12[k]
                  + f_0 * gl_12[k];

        t_13[k] = -3.0 * dl_13[k]
                  + f_0 * gl_13[k];

        t_14[k] = -3.0 * dl_14[k]
                  + f_0 * gl_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, dl_15, dl_16, dl_17, dl_18, dl_19, \
                         gl_15, gl_16, gl_17, gl_18, gl_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -3.0 * dl_15[k]
                  + f_0 * gl_15[k];

        t_16[k] = -3.0 * dl_16[k]
                  + f_0 * gl_16[k];

        t_17[k] = -3.0 * dl_17[k]
                  + f_0 * gl_17[k];

        t_18[k] = -3.0 * dl_18[k]
                  + f_0 * gl_18[k];

        t_19[k] = -3.0 * dl_19[k]
                  + f_0 * gl_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, dl_20, dl_21, dl_22, dl_23, dl_24, \
                         gl_20, gl_21, gl_22, gl_23, gl_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -3.0 * dl_20[k]
                  + f_0 * gl_20[k];

        t_21[k] = -3.0 * dl_21[k]
                  + f_0 * gl_21[k];

        t_22[k] = -3.0 * dl_22[k]
                  + f_0 * gl_22[k];

        t_23[k] = -3.0 * dl_23[k]
                  + f_0 * gl_23[k];

        t_24[k] = -3.0 * dl_24[k]
                  + f_0 * gl_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, dl_25, dl_26, dl_27, dl_28, dl_29, \
                         gl_25, gl_26, gl_27, gl_28, gl_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -3.0 * dl_25[k]
                  + f_0 * gl_25[k];

        t_26[k] = -3.0 * dl_26[k]
                  + f_0 * gl_26[k];

        t_27[k] = -3.0 * dl_27[k]
                  + f_0 * gl_27[k];

        t_28[k] = -3.0 * dl_28[k]
                  + f_0 * gl_28[k];

        t_29[k] = -3.0 * dl_29[k]
                  + f_0 * gl_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, dl_30, dl_31, dl_32, dl_33, dl_34, \
                         gl_30, gl_31, gl_32, gl_33, gl_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -3.0 * dl_30[k]
                  + f_0 * gl_30[k];

        t_31[k] = -3.0 * dl_31[k]
                  + f_0 * gl_31[k];

        t_32[k] = -3.0 * dl_32[k]
                  + f_0 * gl_32[k];

        t_33[k] = -3.0 * dl_33[k]
                  + f_0 * gl_33[k];

        t_34[k] = -3.0 * dl_34[k]
                  + f_0 * gl_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, dl_35, dl_36, dl_37, dl_38, dl_39, \
                         gl_35, gl_36, gl_37, gl_38, gl_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -3.0 * dl_35[k]
                  + f_0 * gl_35[k];

        t_36[k] = -3.0 * dl_36[k]
                  + f_0 * gl_36[k];

        t_37[k] = -3.0 * dl_37[k]
                  + f_0 * gl_37[k];

        t_38[k] = -3.0 * dl_38[k]
                  + f_0 * gl_38[k];

        t_39[k] = -3.0 * dl_39[k]
                  + f_0 * gl_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, dl_40, dl_41, dl_42, dl_43, dl_44, \
                         gl_40, gl_41, gl_42, gl_43, gl_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * dl_40[k]
                  + f_0 * gl_40[k];

        t_41[k] = -3.0 * dl_41[k]
                  + f_0 * gl_41[k];

        t_42[k] = -3.0 * dl_42[k]
                  + f_0 * gl_42[k];

        t_43[k] = -3.0 * dl_43[k]
                  + f_0 * gl_43[k];

        t_44[k] = -3.0 * dl_44[k]
                  + f_0 * gl_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, dl_45, dl_46, dl_47, dl_48, dl_49, \
                         gl_45, gl_46, gl_47, gl_48, gl_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * dl_45[k]
                  + f_0 * gl_45[k];

        t_46[k] = -2.0 * dl_46[k]
                  + f_0 * gl_46[k];

        t_47[k] = -2.0 * dl_47[k]
                  + f_0 * gl_47[k];

        t_48[k] = -2.0 * dl_48[k]
                  + f_0 * gl_48[k];

        t_49[k] = -2.0 * dl_49[k]
                  + f_0 * gl_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, dl_50, dl_51, dl_52, dl_53, dl_54, \
                         gl_50, gl_51, gl_52, gl_53, gl_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * dl_50[k]
                  + f_0 * gl_50[k];

        t_51[k] = -2.0 * dl_51[k]
                  + f_0 * gl_51[k];

        t_52[k] = -2.0 * dl_52[k]
                  + f_0 * gl_52[k];

        t_53[k] = -2.0 * dl_53[k]
                  + f_0 * gl_53[k];

        t_54[k] = -2.0 * dl_54[k]
                  + f_0 * gl_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, dl_55, dl_56, dl_57, dl_58, dl_59, \
                         gl_55, gl_56, gl_57, gl_58, gl_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * dl_55[k]
                  + f_0 * gl_55[k];

        t_56[k] = -2.0 * dl_56[k]
                  + f_0 * gl_56[k];

        t_57[k] = -2.0 * dl_57[k]
                  + f_0 * gl_57[k];

        t_58[k] = -2.0 * dl_58[k]
                  + f_0 * gl_58[k];

        t_59[k] = -2.0 * dl_59[k]
                  + f_0 * gl_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, dl_60, dl_61, dl_62, dl_63, dl_64, \
                         gl_60, gl_61, gl_62, gl_63, gl_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -2.0 * dl_60[k]
                  + f_0 * gl_60[k];

        t_61[k] = -2.0 * dl_61[k]
                  + f_0 * gl_61[k];

        t_62[k] = -2.0 * dl_62[k]
                  + f_0 * gl_62[k];

        t_63[k] = -2.0 * dl_63[k]
                  + f_0 * gl_63[k];

        t_64[k] = -2.0 * dl_64[k]
                  + f_0 * gl_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, dl_65, dl_66, dl_67, dl_68, dl_69, \
                         gl_65, gl_66, gl_67, gl_68, gl_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -2.0 * dl_65[k]
                  + f_0 * gl_65[k];

        t_66[k] = -2.0 * dl_66[k]
                  + f_0 * gl_66[k];

        t_67[k] = -2.0 * dl_67[k]
                  + f_0 * gl_67[k];

        t_68[k] = -2.0 * dl_68[k]
                  + f_0 * gl_68[k];

        t_69[k] = -2.0 * dl_69[k]
                  + f_0 * gl_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, dl_70, dl_71, dl_72, dl_73, dl_74, \
                         gl_70, gl_71, gl_72, gl_73, gl_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * dl_70[k]
                  + f_0 * gl_70[k];

        t_71[k] = -2.0 * dl_71[k]
                  + f_0 * gl_71[k];

        t_72[k] = -2.0 * dl_72[k]
                  + f_0 * gl_72[k];

        t_73[k] = -2.0 * dl_73[k]
                  + f_0 * gl_73[k];

        t_74[k] = -2.0 * dl_74[k]
                  + f_0 * gl_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, dl_75, dl_76, dl_77, dl_78, dl_79, \
                         gl_75, gl_76, gl_77, gl_78, gl_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * dl_75[k]
                  + f_0 * gl_75[k];

        t_76[k] = -2.0 * dl_76[k]
                  + f_0 * gl_76[k];

        t_77[k] = -2.0 * dl_77[k]
                  + f_0 * gl_77[k];

        t_78[k] = -2.0 * dl_78[k]
                  + f_0 * gl_78[k];

        t_79[k] = -2.0 * dl_79[k]
                  + f_0 * gl_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, dl_80, dl_81, dl_82, dl_83, dl_84, \
                         gl_80, gl_81, gl_82, gl_83, gl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -2.0 * dl_80[k]
                  + f_0 * gl_80[k];

        t_81[k] = -2.0 * dl_81[k]
                  + f_0 * gl_81[k];

        t_82[k] = -2.0 * dl_82[k]
                  + f_0 * gl_82[k];

        t_83[k] = -2.0 * dl_83[k]
                  + f_0 * gl_83[k];

        t_84[k] = -2.0 * dl_84[k]
                  + f_0 * gl_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, dl_85, dl_86, dl_87, dl_88, dl_89, \
                         gl_85, gl_86, gl_87, gl_88, gl_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * dl_85[k]
                  + f_0 * gl_85[k];

        t_86[k] = -2.0 * dl_86[k]
                  + f_0 * gl_86[k];

        t_87[k] = -2.0 * dl_87[k]
                  + f_0 * gl_87[k];

        t_88[k] = -2.0 * dl_88[k]
                  + f_0 * gl_88[k];

        t_89[k] = -2.0 * dl_89[k]
                  + f_0 * gl_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, dl_90, dl_91, dl_92, dl_93, dl_94, \
                         gl_90, gl_91, gl_92, gl_93, gl_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * dl_90[k]
                  + f_0 * gl_90[k];

        t_91[k] = -2.0 * dl_91[k]
                  + f_0 * gl_91[k];

        t_92[k] = -2.0 * dl_92[k]
                  + f_0 * gl_92[k];

        t_93[k] = -2.0 * dl_93[k]
                  + f_0 * gl_93[k];

        t_94[k] = -2.0 * dl_94[k]
                  + f_0 * gl_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, dl_95, dl_96, dl_97, dl_98, dl_99, \
                         gl_95, gl_96, gl_97, gl_98, gl_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * dl_95[k]
                  + f_0 * gl_95[k];

        t_96[k] = -2.0 * dl_96[k]
                  + f_0 * gl_96[k];

        t_97[k] = -2.0 * dl_97[k]
                  + f_0 * gl_97[k];

        t_98[k] = -2.0 * dl_98[k]
                  + f_0 * gl_98[k];

        t_99[k] = -2.0 * dl_99[k]
                  + f_0 * gl_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, dl_100, dl_101, dl_102, dl_103, \
                         dl_104, gl_100, gl_101, gl_102, gl_103, \
                         gl_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * dl_100[k]
                   + f_0 * gl_100[k];

        t_101[k] = -2.0 * dl_101[k]
                   + f_0 * gl_101[k];

        t_102[k] = -2.0 * dl_102[k]
                   + f_0 * gl_102[k];

        t_103[k] = -2.0 * dl_103[k]
                   + f_0 * gl_103[k];

        t_104[k] = -2.0 * dl_104[k]
                   + f_0 * gl_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, dl_105, dl_106, dl_107, dl_108, \
                         dl_109, gl_105, gl_106, gl_107, gl_108, \
                         gl_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * dl_105[k]
                   + f_0 * gl_105[k];

        t_106[k] = -2.0 * dl_106[k]
                   + f_0 * gl_106[k];

        t_107[k] = -2.0 * dl_107[k]
                   + f_0 * gl_107[k];

        t_108[k] = -2.0 * dl_108[k]
                   + f_0 * gl_108[k];

        t_109[k] = -2.0 * dl_109[k]
                   + f_0 * gl_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, dl_110, dl_111, dl_112, dl_113, \
                         dl_114, gl_110, gl_111, gl_112, gl_113, \
                         gl_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * dl_110[k]
                   + f_0 * gl_110[k];

        t_111[k] = -2.0 * dl_111[k]
                   + f_0 * gl_111[k];

        t_112[k] = -2.0 * dl_112[k]
                   + f_0 * gl_112[k];

        t_113[k] = -2.0 * dl_113[k]
                   + f_0 * gl_113[k];

        t_114[k] = -2.0 * dl_114[k]
                   + f_0 * gl_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, dl_115, dl_116, dl_117, dl_118, \
                         dl_119, gl_115, gl_116, gl_117, gl_118, \
                         gl_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * dl_115[k]
                   + f_0 * gl_115[k];

        t_116[k] = -2.0 * dl_116[k]
                   + f_0 * gl_116[k];

        t_117[k] = -2.0 * dl_117[k]
                   + f_0 * gl_117[k];

        t_118[k] = -2.0 * dl_118[k]
                   + f_0 * gl_118[k];

        t_119[k] = -2.0 * dl_119[k]
                   + f_0 * gl_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, dl_120, dl_121, dl_122, dl_123, \
                         dl_124, gl_120, gl_121, gl_122, gl_123, \
                         gl_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -2.0 * dl_120[k]
                   + f_0 * gl_120[k];

        t_121[k] = -2.0 * dl_121[k]
                   + f_0 * gl_121[k];

        t_122[k] = -2.0 * dl_122[k]
                   + f_0 * gl_122[k];

        t_123[k] = -2.0 * dl_123[k]
                   + f_0 * gl_123[k];

        t_124[k] = -2.0 * dl_124[k]
                   + f_0 * gl_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, dl_125, dl_126, dl_127, dl_128, \
                         dl_129, gl_125, gl_126, gl_127, gl_128, \
                         gl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -2.0 * dl_125[k]
                   + f_0 * gl_125[k];

        t_126[k] = -2.0 * dl_126[k]
                   + f_0 * gl_126[k];

        t_127[k] = -2.0 * dl_127[k]
                   + f_0 * gl_127[k];

        t_128[k] = -2.0 * dl_128[k]
                   + f_0 * gl_128[k];

        t_129[k] = -2.0 * dl_129[k]
                   + f_0 * gl_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, dl_130, dl_131, dl_132, dl_133, \
                         dl_134, gl_130, gl_131, gl_132, gl_133, \
                         gl_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -2.0 * dl_130[k]
                   + f_0 * gl_130[k];

        t_131[k] = -2.0 * dl_131[k]
                   + f_0 * gl_131[k];

        t_132[k] = -2.0 * dl_132[k]
                   + f_0 * gl_132[k];

        t_133[k] = -2.0 * dl_133[k]
                   + f_0 * gl_133[k];

        t_134[k] = -2.0 * dl_134[k]
                   + f_0 * gl_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, dl_135, dl_136, dl_137, dl_138, \
                         dl_139, gl_135, gl_136, gl_137, gl_138, \
                         gl_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -dl_135[k]
                   + f_0 * gl_135[k];

        t_136[k] = -dl_136[k]
                   + f_0 * gl_136[k];

        t_137[k] = -dl_137[k]
                   + f_0 * gl_137[k];

        t_138[k] = -dl_138[k]
                   + f_0 * gl_138[k];

        t_139[k] = -dl_139[k]
                   + f_0 * gl_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, dl_140, dl_141, dl_142, dl_143, \
                         dl_144, gl_140, gl_141, gl_142, gl_143, \
                         gl_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -dl_140[k]
                   + f_0 * gl_140[k];

        t_141[k] = -dl_141[k]
                   + f_0 * gl_141[k];

        t_142[k] = -dl_142[k]
                   + f_0 * gl_142[k];

        t_143[k] = -dl_143[k]
                   + f_0 * gl_143[k];

        t_144[k] = -dl_144[k]
                   + f_0 * gl_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, dl_145, dl_146, dl_147, dl_148, \
                         dl_149, gl_145, gl_146, gl_147, gl_148, \
                         gl_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -dl_145[k]
                   + f_0 * gl_145[k];

        t_146[k] = -dl_146[k]
                   + f_0 * gl_146[k];

        t_147[k] = -dl_147[k]
                   + f_0 * gl_147[k];

        t_148[k] = -dl_148[k]
                   + f_0 * gl_148[k];

        t_149[k] = -dl_149[k]
                   + f_0 * gl_149[k];
    }
}

static auto
compute_prim_geom_10_fl_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dl, const size_t gl,
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

    const auto *dl_150 = buffer.data(dl + 150);
    const auto *dl_151 = buffer.data(dl + 151);
    const auto *dl_152 = buffer.data(dl + 152);
    const auto *dl_153 = buffer.data(dl + 153);
    const auto *dl_154 = buffer.data(dl + 154);
    const auto *dl_155 = buffer.data(dl + 155);
    const auto *dl_156 = buffer.data(dl + 156);
    const auto *dl_157 = buffer.data(dl + 157);
    const auto *dl_158 = buffer.data(dl + 158);
    const auto *dl_159 = buffer.data(dl + 159);
    const auto *dl_160 = buffer.data(dl + 160);
    const auto *dl_161 = buffer.data(dl + 161);
    const auto *dl_162 = buffer.data(dl + 162);
    const auto *dl_163 = buffer.data(dl + 163);
    const auto *dl_164 = buffer.data(dl + 164);
    const auto *dl_165 = buffer.data(dl + 165);
    const auto *dl_166 = buffer.data(dl + 166);
    const auto *dl_167 = buffer.data(dl + 167);
    const auto *dl_168 = buffer.data(dl + 168);
    const auto *dl_169 = buffer.data(dl + 169);
    const auto *dl_170 = buffer.data(dl + 170);
    const auto *dl_171 = buffer.data(dl + 171);
    const auto *dl_172 = buffer.data(dl + 172);
    const auto *dl_173 = buffer.data(dl + 173);
    const auto *dl_174 = buffer.data(dl + 174);
    const auto *dl_175 = buffer.data(dl + 175);
    const auto *dl_176 = buffer.data(dl + 176);
    const auto *dl_177 = buffer.data(dl + 177);
    const auto *dl_178 = buffer.data(dl + 178);
    const auto *dl_179 = buffer.data(dl + 179);
    const auto *dl_180 = buffer.data(dl + 180);
    const auto *dl_181 = buffer.data(dl + 181);
    const auto *dl_182 = buffer.data(dl + 182);
    const auto *dl_183 = buffer.data(dl + 183);
    const auto *dl_184 = buffer.data(dl + 184);
    const auto *dl_185 = buffer.data(dl + 185);
    const auto *dl_186 = buffer.data(dl + 186);
    const auto *dl_187 = buffer.data(dl + 187);
    const auto *dl_188 = buffer.data(dl + 188);
    const auto *dl_189 = buffer.data(dl + 189);
    const auto *dl_190 = buffer.data(dl + 190);
    const auto *dl_191 = buffer.data(dl + 191);
    const auto *dl_192 = buffer.data(dl + 192);
    const auto *dl_193 = buffer.data(dl + 193);
    const auto *dl_194 = buffer.data(dl + 194);
    const auto *dl_195 = buffer.data(dl + 195);
    const auto *dl_196 = buffer.data(dl + 196);
    const auto *dl_197 = buffer.data(dl + 197);
    const auto *dl_198 = buffer.data(dl + 198);
    const auto *dl_199 = buffer.data(dl + 199);
    const auto *dl_200 = buffer.data(dl + 200);
    const auto *dl_201 = buffer.data(dl + 201);
    const auto *dl_202 = buffer.data(dl + 202);
    const auto *dl_203 = buffer.data(dl + 203);
    const auto *dl_204 = buffer.data(dl + 204);
    const auto *dl_205 = buffer.data(dl + 205);
    const auto *dl_206 = buffer.data(dl + 206);
    const auto *dl_207 = buffer.data(dl + 207);
    const auto *dl_208 = buffer.data(dl + 208);
    const auto *dl_209 = buffer.data(dl + 209);
    const auto *dl_210 = buffer.data(dl + 210);
    const auto *dl_211 = buffer.data(dl + 211);
    const auto *dl_212 = buffer.data(dl + 212);
    const auto *dl_213 = buffer.data(dl + 213);
    const auto *dl_214 = buffer.data(dl + 214);
    const auto *dl_215 = buffer.data(dl + 215);
    const auto *dl_216 = buffer.data(dl + 216);
    const auto *dl_217 = buffer.data(dl + 217);
    const auto *dl_218 = buffer.data(dl + 218);
    const auto *dl_219 = buffer.data(dl + 219);
    const auto *dl_220 = buffer.data(dl + 220);
    const auto *dl_221 = buffer.data(dl + 221);
    const auto *dl_222 = buffer.data(dl + 222);
    const auto *dl_223 = buffer.data(dl + 223);
    const auto *dl_224 = buffer.data(dl + 224);
    const auto *dl_225 = buffer.data(dl + 225);
    const auto *dl_226 = buffer.data(dl + 226);
    const auto *dl_227 = buffer.data(dl + 227);
    const auto *dl_228 = buffer.data(dl + 228);
    const auto *dl_229 = buffer.data(dl + 229);
    const auto *dl_230 = buffer.data(dl + 230);
    const auto *dl_231 = buffer.data(dl + 231);
    const auto *dl_232 = buffer.data(dl + 232);
    const auto *dl_233 = buffer.data(dl + 233);
    const auto *dl_234 = buffer.data(dl + 234);
    const auto *dl_235 = buffer.data(dl + 235);
    const auto *dl_236 = buffer.data(dl + 236);
    const auto *dl_237 = buffer.data(dl + 237);
    const auto *dl_238 = buffer.data(dl + 238);
    const auto *dl_239 = buffer.data(dl + 239);
    const auto *dl_240 = buffer.data(dl + 240);
    const auto *dl_241 = buffer.data(dl + 241);
    const auto *dl_242 = buffer.data(dl + 242);
    const auto *dl_243 = buffer.data(dl + 243);
    const auto *dl_244 = buffer.data(dl + 244);
    const auto *dl_245 = buffer.data(dl + 245);
    const auto *dl_246 = buffer.data(dl + 246);
    const auto *dl_247 = buffer.data(dl + 247);
    const auto *dl_248 = buffer.data(dl + 248);
    const auto *dl_249 = buffer.data(dl + 249);
    const auto *dl_250 = buffer.data(dl + 250);
    const auto *dl_251 = buffer.data(dl + 251);
    const auto *dl_252 = buffer.data(dl + 252);
    const auto *dl_253 = buffer.data(dl + 253);
    const auto *dl_254 = buffer.data(dl + 254);
    const auto *dl_255 = buffer.data(dl + 255);
    const auto *dl_256 = buffer.data(dl + 256);
    const auto *dl_257 = buffer.data(dl + 257);
    const auto *dl_258 = buffer.data(dl + 258);
    const auto *dl_259 = buffer.data(dl + 259);
    const auto *dl_260 = buffer.data(dl + 260);
    const auto *dl_261 = buffer.data(dl + 261);
    const auto *dl_262 = buffer.data(dl + 262);
    const auto *dl_263 = buffer.data(dl + 263);
    const auto *dl_264 = buffer.data(dl + 264);
    const auto *dl_265 = buffer.data(dl + 265);
    const auto *dl_266 = buffer.data(dl + 266);
    const auto *dl_267 = buffer.data(dl + 267);
    const auto *dl_268 = buffer.data(dl + 268);
    const auto *dl_269 = buffer.data(dl + 269);

    const auto *gl_150 = buffer.data(gl + 150);
    const auto *gl_151 = buffer.data(gl + 151);
    const auto *gl_152 = buffer.data(gl + 152);
    const auto *gl_153 = buffer.data(gl + 153);
    const auto *gl_154 = buffer.data(gl + 154);
    const auto *gl_155 = buffer.data(gl + 155);
    const auto *gl_156 = buffer.data(gl + 156);
    const auto *gl_157 = buffer.data(gl + 157);
    const auto *gl_158 = buffer.data(gl + 158);
    const auto *gl_159 = buffer.data(gl + 159);
    const auto *gl_160 = buffer.data(gl + 160);
    const auto *gl_161 = buffer.data(gl + 161);
    const auto *gl_162 = buffer.data(gl + 162);
    const auto *gl_163 = buffer.data(gl + 163);
    const auto *gl_164 = buffer.data(gl + 164);
    const auto *gl_165 = buffer.data(gl + 165);
    const auto *gl_166 = buffer.data(gl + 166);
    const auto *gl_167 = buffer.data(gl + 167);
    const auto *gl_168 = buffer.data(gl + 168);
    const auto *gl_169 = buffer.data(gl + 169);
    const auto *gl_170 = buffer.data(gl + 170);
    const auto *gl_171 = buffer.data(gl + 171);
    const auto *gl_172 = buffer.data(gl + 172);
    const auto *gl_173 = buffer.data(gl + 173);
    const auto *gl_174 = buffer.data(gl + 174);
    const auto *gl_175 = buffer.data(gl + 175);
    const auto *gl_176 = buffer.data(gl + 176);
    const auto *gl_177 = buffer.data(gl + 177);
    const auto *gl_178 = buffer.data(gl + 178);
    const auto *gl_179 = buffer.data(gl + 179);
    const auto *gl_180 = buffer.data(gl + 180);
    const auto *gl_181 = buffer.data(gl + 181);
    const auto *gl_182 = buffer.data(gl + 182);
    const auto *gl_183 = buffer.data(gl + 183);
    const auto *gl_184 = buffer.data(gl + 184);
    const auto *gl_185 = buffer.data(gl + 185);
    const auto *gl_186 = buffer.data(gl + 186);
    const auto *gl_187 = buffer.data(gl + 187);
    const auto *gl_188 = buffer.data(gl + 188);
    const auto *gl_189 = buffer.data(gl + 189);
    const auto *gl_190 = buffer.data(gl + 190);
    const auto *gl_191 = buffer.data(gl + 191);
    const auto *gl_192 = buffer.data(gl + 192);
    const auto *gl_193 = buffer.data(gl + 193);
    const auto *gl_194 = buffer.data(gl + 194);
    const auto *gl_195 = buffer.data(gl + 195);
    const auto *gl_196 = buffer.data(gl + 196);
    const auto *gl_197 = buffer.data(gl + 197);
    const auto *gl_198 = buffer.data(gl + 198);
    const auto *gl_199 = buffer.data(gl + 199);
    const auto *gl_200 = buffer.data(gl + 200);
    const auto *gl_201 = buffer.data(gl + 201);
    const auto *gl_202 = buffer.data(gl + 202);
    const auto *gl_203 = buffer.data(gl + 203);
    const auto *gl_204 = buffer.data(gl + 204);
    const auto *gl_205 = buffer.data(gl + 205);
    const auto *gl_206 = buffer.data(gl + 206);
    const auto *gl_207 = buffer.data(gl + 207);
    const auto *gl_208 = buffer.data(gl + 208);
    const auto *gl_209 = buffer.data(gl + 209);
    const auto *gl_210 = buffer.data(gl + 210);
    const auto *gl_211 = buffer.data(gl + 211);
    const auto *gl_212 = buffer.data(gl + 212);
    const auto *gl_213 = buffer.data(gl + 213);
    const auto *gl_214 = buffer.data(gl + 214);
    const auto *gl_215 = buffer.data(gl + 215);
    const auto *gl_216 = buffer.data(gl + 216);
    const auto *gl_217 = buffer.data(gl + 217);
    const auto *gl_218 = buffer.data(gl + 218);
    const auto *gl_219 = buffer.data(gl + 219);
    const auto *gl_220 = buffer.data(gl + 220);
    const auto *gl_221 = buffer.data(gl + 221);
    const auto *gl_222 = buffer.data(gl + 222);
    const auto *gl_223 = buffer.data(gl + 223);
    const auto *gl_224 = buffer.data(gl + 224);
    const auto *gl_225 = buffer.data(gl + 225);
    const auto *gl_226 = buffer.data(gl + 226);
    const auto *gl_227 = buffer.data(gl + 227);
    const auto *gl_228 = buffer.data(gl + 228);
    const auto *gl_229 = buffer.data(gl + 229);
    const auto *gl_230 = buffer.data(gl + 230);
    const auto *gl_231 = buffer.data(gl + 231);
    const auto *gl_232 = buffer.data(gl + 232);
    const auto *gl_233 = buffer.data(gl + 233);
    const auto *gl_234 = buffer.data(gl + 234);
    const auto *gl_235 = buffer.data(gl + 235);
    const auto *gl_236 = buffer.data(gl + 236);
    const auto *gl_237 = buffer.data(gl + 237);
    const auto *gl_238 = buffer.data(gl + 238);
    const auto *gl_239 = buffer.data(gl + 239);
    const auto *gl_240 = buffer.data(gl + 240);
    const auto *gl_241 = buffer.data(gl + 241);
    const auto *gl_242 = buffer.data(gl + 242);
    const auto *gl_243 = buffer.data(gl + 243);
    const auto *gl_244 = buffer.data(gl + 244);
    const auto *gl_245 = buffer.data(gl + 245);
    const auto *gl_246 = buffer.data(gl + 246);
    const auto *gl_247 = buffer.data(gl + 247);
    const auto *gl_248 = buffer.data(gl + 248);
    const auto *gl_249 = buffer.data(gl + 249);
    const auto *gl_250 = buffer.data(gl + 250);
    const auto *gl_251 = buffer.data(gl + 251);
    const auto *gl_252 = buffer.data(gl + 252);
    const auto *gl_253 = buffer.data(gl + 253);
    const auto *gl_254 = buffer.data(gl + 254);
    const auto *gl_255 = buffer.data(gl + 255);
    const auto *gl_256 = buffer.data(gl + 256);
    const auto *gl_257 = buffer.data(gl + 257);
    const auto *gl_258 = buffer.data(gl + 258);
    const auto *gl_259 = buffer.data(gl + 259);
    const auto *gl_260 = buffer.data(gl + 260);
    const auto *gl_261 = buffer.data(gl + 261);
    const auto *gl_262 = buffer.data(gl + 262);
    const auto *gl_263 = buffer.data(gl + 263);
    const auto *gl_264 = buffer.data(gl + 264);
    const auto *gl_265 = buffer.data(gl + 265);
    const auto *gl_266 = buffer.data(gl + 266);
    const auto *gl_267 = buffer.data(gl + 267);
    const auto *gl_268 = buffer.data(gl + 268);
    const auto *gl_269 = buffer.data(gl + 269);
    const auto *gl_270 = buffer.data(gl + 270);
    const auto *gl_271 = buffer.data(gl + 271);
    const auto *gl_272 = buffer.data(gl + 272);
    const auto *gl_273 = buffer.data(gl + 273);
    const auto *gl_274 = buffer.data(gl + 274);
    const auto *gl_275 = buffer.data(gl + 275);
    const auto *gl_276 = buffer.data(gl + 276);
    const auto *gl_277 = buffer.data(gl + 277);
    const auto *gl_278 = buffer.data(gl + 278);
    const auto *gl_279 = buffer.data(gl + 279);
    const auto *gl_280 = buffer.data(gl + 280);
    const auto *gl_281 = buffer.data(gl + 281);
    const auto *gl_282 = buffer.data(gl + 282);
    const auto *gl_283 = buffer.data(gl + 283);
    const auto *gl_284 = buffer.data(gl + 284);
    const auto *gl_285 = buffer.data(gl + 285);
    const auto *gl_286 = buffer.data(gl + 286);
    const auto *gl_287 = buffer.data(gl + 287);
    const auto *gl_288 = buffer.data(gl + 288);
    const auto *gl_289 = buffer.data(gl + 289);
    const auto *gl_290 = buffer.data(gl + 290);
    const auto *gl_291 = buffer.data(gl + 291);
    const auto *gl_292 = buffer.data(gl + 292);
    const auto *gl_293 = buffer.data(gl + 293);
    const auto *gl_294 = buffer.data(gl + 294);
    const auto *gl_295 = buffer.data(gl + 295);
    const auto *gl_296 = buffer.data(gl + 296);
    const auto *gl_297 = buffer.data(gl + 297);
    const auto *gl_298 = buffer.data(gl + 298);
    const auto *gl_299 = buffer.data(gl + 299);
    const auto *gl_300 = buffer.data(gl + 300);
    const auto *gl_301 = buffer.data(gl + 301);
    const auto *gl_302 = buffer.data(gl + 302);
    const auto *gl_303 = buffer.data(gl + 303);
    const auto *gl_304 = buffer.data(gl + 304);
    const auto *gl_305 = buffer.data(gl + 305);
    const auto *gl_306 = buffer.data(gl + 306);
    const auto *gl_307 = buffer.data(gl + 307);
    const auto *gl_308 = buffer.data(gl + 308);
    const auto *gl_309 = buffer.data(gl + 309);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, dl_150, dl_151, dl_152, dl_153, \
                         dl_154, gl_150, gl_151, gl_152, gl_153, \
                         gl_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -dl_150[k]
                   + f_0 * gl_150[k];

        t_151[k] = -dl_151[k]
                   + f_0 * gl_151[k];

        t_152[k] = -dl_152[k]
                   + f_0 * gl_152[k];

        t_153[k] = -dl_153[k]
                   + f_0 * gl_153[k];

        t_154[k] = -dl_154[k]
                   + f_0 * gl_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, dl_155, dl_156, dl_157, dl_158, \
                         dl_159, gl_155, gl_156, gl_157, gl_158, \
                         gl_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -dl_155[k]
                   + f_0 * gl_155[k];

        t_156[k] = -dl_156[k]
                   + f_0 * gl_156[k];

        t_157[k] = -dl_157[k]
                   + f_0 * gl_157[k];

        t_158[k] = -dl_158[k]
                   + f_0 * gl_158[k];

        t_159[k] = -dl_159[k]
                   + f_0 * gl_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, dl_160, dl_161, dl_162, dl_163, \
                         dl_164, gl_160, gl_161, gl_162, gl_163, \
                         gl_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -dl_160[k]
                   + f_0 * gl_160[k];

        t_161[k] = -dl_161[k]
                   + f_0 * gl_161[k];

        t_162[k] = -dl_162[k]
                   + f_0 * gl_162[k];

        t_163[k] = -dl_163[k]
                   + f_0 * gl_163[k];

        t_164[k] = -dl_164[k]
                   + f_0 * gl_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, dl_165, dl_166, dl_167, dl_168, \
                         dl_169, gl_165, gl_166, gl_167, gl_168, \
                         gl_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -dl_165[k]
                   + f_0 * gl_165[k];

        t_166[k] = -dl_166[k]
                   + f_0 * gl_166[k];

        t_167[k] = -dl_167[k]
                   + f_0 * gl_167[k];

        t_168[k] = -dl_168[k]
                   + f_0 * gl_168[k];

        t_169[k] = -dl_169[k]
                   + f_0 * gl_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, dl_170, dl_171, dl_172, dl_173, \
                         dl_174, gl_170, gl_171, gl_172, gl_173, \
                         gl_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -dl_170[k]
                   + f_0 * gl_170[k];

        t_171[k] = -dl_171[k]
                   + f_0 * gl_171[k];

        t_172[k] = -dl_172[k]
                   + f_0 * gl_172[k];

        t_173[k] = -dl_173[k]
                   + f_0 * gl_173[k];

        t_174[k] = -dl_174[k]
                   + f_0 * gl_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, dl_175, dl_176, dl_177, dl_178, \
                         dl_179, gl_175, gl_176, gl_177, gl_178, \
                         gl_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -dl_175[k]
                   + f_0 * gl_175[k];

        t_176[k] = -dl_176[k]
                   + f_0 * gl_176[k];

        t_177[k] = -dl_177[k]
                   + f_0 * gl_177[k];

        t_178[k] = -dl_178[k]
                   + f_0 * gl_178[k];

        t_179[k] = -dl_179[k]
                   + f_0 * gl_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, dl_180, dl_181, dl_182, dl_183, \
                         dl_184, gl_180, gl_181, gl_182, gl_183, \
                         gl_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -dl_180[k]
                   + f_0 * gl_180[k];

        t_181[k] = -dl_181[k]
                   + f_0 * gl_181[k];

        t_182[k] = -dl_182[k]
                   + f_0 * gl_182[k];

        t_183[k] = -dl_183[k]
                   + f_0 * gl_183[k];

        t_184[k] = -dl_184[k]
                   + f_0 * gl_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, dl_185, dl_186, dl_187, dl_188, \
                         dl_189, gl_185, gl_186, gl_187, gl_188, \
                         gl_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -dl_185[k]
                   + f_0 * gl_185[k];

        t_186[k] = -dl_186[k]
                   + f_0 * gl_186[k];

        t_187[k] = -dl_187[k]
                   + f_0 * gl_187[k];

        t_188[k] = -dl_188[k]
                   + f_0 * gl_188[k];

        t_189[k] = -dl_189[k]
                   + f_0 * gl_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, dl_190, dl_191, dl_192, dl_193, \
                         dl_194, gl_190, gl_191, gl_192, gl_193, \
                         gl_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -dl_190[k]
                   + f_0 * gl_190[k];

        t_191[k] = -dl_191[k]
                   + f_0 * gl_191[k];

        t_192[k] = -dl_192[k]
                   + f_0 * gl_192[k];

        t_193[k] = -dl_193[k]
                   + f_0 * gl_193[k];

        t_194[k] = -dl_194[k]
                   + f_0 * gl_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, dl_195, dl_196, dl_197, dl_198, \
                         dl_199, gl_195, gl_196, gl_197, gl_198, \
                         gl_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -dl_195[k]
                   + f_0 * gl_195[k];

        t_196[k] = -dl_196[k]
                   + f_0 * gl_196[k];

        t_197[k] = -dl_197[k]
                   + f_0 * gl_197[k];

        t_198[k] = -dl_198[k]
                   + f_0 * gl_198[k];

        t_199[k] = -dl_199[k]
                   + f_0 * gl_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, dl_200, dl_201, dl_202, dl_203, \
                         dl_204, gl_200, gl_201, gl_202, gl_203, \
                         gl_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -dl_200[k]
                   + f_0 * gl_200[k];

        t_201[k] = -dl_201[k]
                   + f_0 * gl_201[k];

        t_202[k] = -dl_202[k]
                   + f_0 * gl_202[k];

        t_203[k] = -dl_203[k]
                   + f_0 * gl_203[k];

        t_204[k] = -dl_204[k]
                   + f_0 * gl_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, dl_205, dl_206, dl_207, dl_208, \
                         dl_209, gl_205, gl_206, gl_207, gl_208, \
                         gl_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -dl_205[k]
                   + f_0 * gl_205[k];

        t_206[k] = -dl_206[k]
                   + f_0 * gl_206[k];

        t_207[k] = -dl_207[k]
                   + f_0 * gl_207[k];

        t_208[k] = -dl_208[k]
                   + f_0 * gl_208[k];

        t_209[k] = -dl_209[k]
                   + f_0 * gl_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, dl_210, dl_211, dl_212, dl_213, \
                         dl_214, gl_210, gl_211, gl_212, gl_213, \
                         gl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -dl_210[k]
                   + f_0 * gl_210[k];

        t_211[k] = -dl_211[k]
                   + f_0 * gl_211[k];

        t_212[k] = -dl_212[k]
                   + f_0 * gl_212[k];

        t_213[k] = -dl_213[k]
                   + f_0 * gl_213[k];

        t_214[k] = -dl_214[k]
                   + f_0 * gl_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, dl_215, dl_216, dl_217, dl_218, \
                         dl_219, gl_215, gl_216, gl_217, gl_218, \
                         gl_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -dl_215[k]
                   + f_0 * gl_215[k];

        t_216[k] = -dl_216[k]
                   + f_0 * gl_216[k];

        t_217[k] = -dl_217[k]
                   + f_0 * gl_217[k];

        t_218[k] = -dl_218[k]
                   + f_0 * gl_218[k];

        t_219[k] = -dl_219[k]
                   + f_0 * gl_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, dl_220, dl_221, dl_222, dl_223, \
                         dl_224, gl_220, gl_221, gl_222, gl_223, \
                         gl_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -dl_220[k]
                   + f_0 * gl_220[k];

        t_221[k] = -dl_221[k]
                   + f_0 * gl_221[k];

        t_222[k] = -dl_222[k]
                   + f_0 * gl_222[k];

        t_223[k] = -dl_223[k]
                   + f_0 * gl_223[k];

        t_224[k] = -dl_224[k]
                   + f_0 * gl_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, dl_225, dl_226, dl_227, dl_228, \
                         dl_229, gl_225, gl_226, gl_227, gl_228, \
                         gl_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -dl_225[k]
                   + f_0 * gl_225[k];

        t_226[k] = -dl_226[k]
                   + f_0 * gl_226[k];

        t_227[k] = -dl_227[k]
                   + f_0 * gl_227[k];

        t_228[k] = -dl_228[k]
                   + f_0 * gl_228[k];

        t_229[k] = -dl_229[k]
                   + f_0 * gl_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, dl_230, dl_231, dl_232, dl_233, \
                         dl_234, gl_230, gl_231, gl_232, gl_233, \
                         gl_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -dl_230[k]
                   + f_0 * gl_230[k];

        t_231[k] = -dl_231[k]
                   + f_0 * gl_231[k];

        t_232[k] = -dl_232[k]
                   + f_0 * gl_232[k];

        t_233[k] = -dl_233[k]
                   + f_0 * gl_233[k];

        t_234[k] = -dl_234[k]
                   + f_0 * gl_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, dl_235, dl_236, dl_237, dl_238, \
                         dl_239, gl_235, gl_236, gl_237, gl_238, \
                         gl_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -dl_235[k]
                   + f_0 * gl_235[k];

        t_236[k] = -dl_236[k]
                   + f_0 * gl_236[k];

        t_237[k] = -dl_237[k]
                   + f_0 * gl_237[k];

        t_238[k] = -dl_238[k]
                   + f_0 * gl_238[k];

        t_239[k] = -dl_239[k]
                   + f_0 * gl_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, dl_240, dl_241, dl_242, dl_243, \
                         dl_244, gl_240, gl_241, gl_242, gl_243, \
                         gl_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -dl_240[k]
                   + f_0 * gl_240[k];

        t_241[k] = -dl_241[k]
                   + f_0 * gl_241[k];

        t_242[k] = -dl_242[k]
                   + f_0 * gl_242[k];

        t_243[k] = -dl_243[k]
                   + f_0 * gl_243[k];

        t_244[k] = -dl_244[k]
                   + f_0 * gl_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, dl_245, dl_246, dl_247, dl_248, \
                         dl_249, gl_245, gl_246, gl_247, gl_248, \
                         gl_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -dl_245[k]
                   + f_0 * gl_245[k];

        t_246[k] = -dl_246[k]
                   + f_0 * gl_246[k];

        t_247[k] = -dl_247[k]
                   + f_0 * gl_247[k];

        t_248[k] = -dl_248[k]
                   + f_0 * gl_248[k];

        t_249[k] = -dl_249[k]
                   + f_0 * gl_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, dl_250, dl_251, dl_252, dl_253, \
                         dl_254, gl_250, gl_251, gl_252, gl_253, \
                         gl_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -dl_250[k]
                   + f_0 * gl_250[k];

        t_251[k] = -dl_251[k]
                   + f_0 * gl_251[k];

        t_252[k] = -dl_252[k]
                   + f_0 * gl_252[k];

        t_253[k] = -dl_253[k]
                   + f_0 * gl_253[k];

        t_254[k] = -dl_254[k]
                   + f_0 * gl_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, dl_255, dl_256, dl_257, dl_258, \
                         dl_259, gl_255, gl_256, gl_257, gl_258, \
                         gl_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -dl_255[k]
                   + f_0 * gl_255[k];

        t_256[k] = -dl_256[k]
                   + f_0 * gl_256[k];

        t_257[k] = -dl_257[k]
                   + f_0 * gl_257[k];

        t_258[k] = -dl_258[k]
                   + f_0 * gl_258[k];

        t_259[k] = -dl_259[k]
                   + f_0 * gl_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, dl_260, dl_261, dl_262, dl_263, \
                         dl_264, gl_260, gl_261, gl_262, gl_263, \
                         gl_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -dl_260[k]
                   + f_0 * gl_260[k];

        t_261[k] = -dl_261[k]
                   + f_0 * gl_261[k];

        t_262[k] = -dl_262[k]
                   + f_0 * gl_262[k];

        t_263[k] = -dl_263[k]
                   + f_0 * gl_263[k];

        t_264[k] = -dl_264[k]
                   + f_0 * gl_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, dl_265, dl_266, dl_267, dl_268, \
                         dl_269, gl_265, gl_266, gl_267, gl_268, \
                         gl_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -dl_265[k]
                   + f_0 * gl_265[k];

        t_266[k] = -dl_266[k]
                   + f_0 * gl_266[k];

        t_267[k] = -dl_267[k]
                   + f_0 * gl_267[k];

        t_268[k] = -dl_268[k]
                   + f_0 * gl_268[k];

        t_269[k] = -dl_269[k]
                   + f_0 * gl_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, t_276, t_277, gl_270, \
                         gl_271, gl_272, gl_273, gl_274, gl_275, gl_276, \
                         gl_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * gl_270[k];

        t_271[k] = f_0 * gl_271[k];

        t_272[k] = f_0 * gl_272[k];

        t_273[k] = f_0 * gl_273[k];

        t_274[k] = f_0 * gl_274[k];

        t_275[k] = f_0 * gl_275[k];

        t_276[k] = f_0 * gl_276[k];

        t_277[k] = f_0 * gl_277[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, t_284, t_285, gl_278, \
                         gl_279, gl_280, gl_281, gl_282, gl_283, gl_284, \
                         gl_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * gl_278[k];

        t_279[k] = f_0 * gl_279[k];

        t_280[k] = f_0 * gl_280[k];

        t_281[k] = f_0 * gl_281[k];

        t_282[k] = f_0 * gl_282[k];

        t_283[k] = f_0 * gl_283[k];

        t_284[k] = f_0 * gl_284[k];

        t_285[k] = f_0 * gl_285[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, t_291, t_292, t_293, gl_286, \
                         gl_287, gl_288, gl_289, gl_290, gl_291, gl_292, \
                         gl_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * gl_286[k];

        t_287[k] = f_0 * gl_287[k];

        t_288[k] = f_0 * gl_288[k];

        t_289[k] = f_0 * gl_289[k];

        t_290[k] = f_0 * gl_290[k];

        t_291[k] = f_0 * gl_291[k];

        t_292[k] = f_0 * gl_292[k];

        t_293[k] = f_0 * gl_293[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, t_299, t_300, t_301, gl_294, \
                         gl_295, gl_296, gl_297, gl_298, gl_299, gl_300, \
                         gl_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_0 * gl_294[k];

        t_295[k] = f_0 * gl_295[k];

        t_296[k] = f_0 * gl_296[k];

        t_297[k] = f_0 * gl_297[k];

        t_298[k] = f_0 * gl_298[k];

        t_299[k] = f_0 * gl_299[k];

        t_300[k] = f_0 * gl_300[k];

        t_301[k] = f_0 * gl_301[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, t_307, t_308, t_309, gl_302, \
                         gl_303, gl_304, gl_305, gl_306, gl_307, gl_308, \
                         gl_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * gl_302[k];

        t_303[k] = f_0 * gl_303[k];

        t_304[k] = f_0 * gl_304[k];

        t_305[k] = f_0 * gl_305[k];

        t_306[k] = f_0 * gl_306[k];

        t_307[k] = f_0 * gl_307[k];

        t_308[k] = f_0 * gl_308[k];

        t_309[k] = f_0 * gl_309[k];
    }
}

static auto
compute_prim_geom_10_fl_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gl, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gl_310 = buffer.data(gl + 310);
    const auto *gl_311 = buffer.data(gl + 311);
    const auto *gl_312 = buffer.data(gl + 312);
    const auto *gl_313 = buffer.data(gl + 313);
    const auto *gl_314 = buffer.data(gl + 314);
    const auto *gl_315 = buffer.data(gl + 315);
    const auto *gl_316 = buffer.data(gl + 316);
    const auto *gl_317 = buffer.data(gl + 317);
    const auto *gl_318 = buffer.data(gl + 318);
    const auto *gl_319 = buffer.data(gl + 319);
    const auto *gl_320 = buffer.data(gl + 320);
    const auto *gl_321 = buffer.data(gl + 321);
    const auto *gl_322 = buffer.data(gl + 322);
    const auto *gl_323 = buffer.data(gl + 323);
    const auto *gl_324 = buffer.data(gl + 324);
    const auto *gl_325 = buffer.data(gl + 325);
    const auto *gl_326 = buffer.data(gl + 326);
    const auto *gl_327 = buffer.data(gl + 327);
    const auto *gl_328 = buffer.data(gl + 328);
    const auto *gl_329 = buffer.data(gl + 329);
    const auto *gl_330 = buffer.data(gl + 330);
    const auto *gl_331 = buffer.data(gl + 331);
    const auto *gl_332 = buffer.data(gl + 332);
    const auto *gl_333 = buffer.data(gl + 333);
    const auto *gl_334 = buffer.data(gl + 334);
    const auto *gl_335 = buffer.data(gl + 335);
    const auto *gl_336 = buffer.data(gl + 336);
    const auto *gl_337 = buffer.data(gl + 337);
    const auto *gl_338 = buffer.data(gl + 338);
    const auto *gl_339 = buffer.data(gl + 339);
    const auto *gl_340 = buffer.data(gl + 340);
    const auto *gl_341 = buffer.data(gl + 341);
    const auto *gl_342 = buffer.data(gl + 342);
    const auto *gl_343 = buffer.data(gl + 343);
    const auto *gl_344 = buffer.data(gl + 344);
    const auto *gl_345 = buffer.data(gl + 345);
    const auto *gl_346 = buffer.data(gl + 346);
    const auto *gl_347 = buffer.data(gl + 347);
    const auto *gl_348 = buffer.data(gl + 348);
    const auto *gl_349 = buffer.data(gl + 349);
    const auto *gl_350 = buffer.data(gl + 350);
    const auto *gl_351 = buffer.data(gl + 351);
    const auto *gl_352 = buffer.data(gl + 352);
    const auto *gl_353 = buffer.data(gl + 353);
    const auto *gl_354 = buffer.data(gl + 354);
    const auto *gl_355 = buffer.data(gl + 355);
    const auto *gl_356 = buffer.data(gl + 356);
    const auto *gl_357 = buffer.data(gl + 357);
    const auto *gl_358 = buffer.data(gl + 358);
    const auto *gl_359 = buffer.data(gl + 359);
    const auto *gl_360 = buffer.data(gl + 360);
    const auto *gl_361 = buffer.data(gl + 361);
    const auto *gl_362 = buffer.data(gl + 362);
    const auto *gl_363 = buffer.data(gl + 363);
    const auto *gl_364 = buffer.data(gl + 364);
    const auto *gl_365 = buffer.data(gl + 365);
    const auto *gl_366 = buffer.data(gl + 366);
    const auto *gl_367 = buffer.data(gl + 367);
    const auto *gl_368 = buffer.data(gl + 368);
    const auto *gl_369 = buffer.data(gl + 369);
    const auto *gl_370 = buffer.data(gl + 370);
    const auto *gl_371 = buffer.data(gl + 371);
    const auto *gl_372 = buffer.data(gl + 372);
    const auto *gl_373 = buffer.data(gl + 373);
    const auto *gl_374 = buffer.data(gl + 374);
    const auto *gl_375 = buffer.data(gl + 375);
    const auto *gl_376 = buffer.data(gl + 376);
    const auto *gl_377 = buffer.data(gl + 377);
    const auto *gl_378 = buffer.data(gl + 378);
    const auto *gl_379 = buffer.data(gl + 379);
    const auto *gl_380 = buffer.data(gl + 380);
    const auto *gl_381 = buffer.data(gl + 381);
    const auto *gl_382 = buffer.data(gl + 382);
    const auto *gl_383 = buffer.data(gl + 383);
    const auto *gl_384 = buffer.data(gl + 384);
    const auto *gl_385 = buffer.data(gl + 385);
    const auto *gl_386 = buffer.data(gl + 386);
    const auto *gl_387 = buffer.data(gl + 387);
    const auto *gl_388 = buffer.data(gl + 388);
    const auto *gl_389 = buffer.data(gl + 389);
    const auto *gl_390 = buffer.data(gl + 390);
    const auto *gl_391 = buffer.data(gl + 391);
    const auto *gl_392 = buffer.data(gl + 392);
    const auto *gl_393 = buffer.data(gl + 393);
    const auto *gl_394 = buffer.data(gl + 394);
    const auto *gl_395 = buffer.data(gl + 395);
    const auto *gl_396 = buffer.data(gl + 396);
    const auto *gl_397 = buffer.data(gl + 397);
    const auto *gl_398 = buffer.data(gl + 398);
    const auto *gl_399 = buffer.data(gl + 399);
    const auto *gl_400 = buffer.data(gl + 400);
    const auto *gl_401 = buffer.data(gl + 401);
    const auto *gl_402 = buffer.data(gl + 402);
    const auto *gl_403 = buffer.data(gl + 403);
    const auto *gl_404 = buffer.data(gl + 404);
    const auto *gl_405 = buffer.data(gl + 405);
    const auto *gl_406 = buffer.data(gl + 406);
    const auto *gl_407 = buffer.data(gl + 407);
    const auto *gl_408 = buffer.data(gl + 408);
    const auto *gl_409 = buffer.data(gl + 409);
    const auto *gl_410 = buffer.data(gl + 410);
    const auto *gl_411 = buffer.data(gl + 411);
    const auto *gl_412 = buffer.data(gl + 412);
    const auto *gl_413 = buffer.data(gl + 413);
    const auto *gl_414 = buffer.data(gl + 414);
    const auto *gl_415 = buffer.data(gl + 415);
    const auto *gl_416 = buffer.data(gl + 416);
    const auto *gl_417 = buffer.data(gl + 417);
    const auto *gl_418 = buffer.data(gl + 418);
    const auto *gl_419 = buffer.data(gl + 419);
    const auto *gl_420 = buffer.data(gl + 420);
    const auto *gl_421 = buffer.data(gl + 421);
    const auto *gl_422 = buffer.data(gl + 422);
    const auto *gl_423 = buffer.data(gl + 423);
    const auto *gl_424 = buffer.data(gl + 424);
    const auto *gl_425 = buffer.data(gl + 425);
    const auto *gl_426 = buffer.data(gl + 426);
    const auto *gl_427 = buffer.data(gl + 427);
    const auto *gl_428 = buffer.data(gl + 428);
    const auto *gl_429 = buffer.data(gl + 429);
    const auto *gl_430 = buffer.data(gl + 430);
    const auto *gl_431 = buffer.data(gl + 431);
    const auto *gl_432 = buffer.data(gl + 432);
    const auto *gl_433 = buffer.data(gl + 433);
    const auto *gl_434 = buffer.data(gl + 434);
    const auto *gl_435 = buffer.data(gl + 435);
    const auto *gl_436 = buffer.data(gl + 436);
    const auto *gl_437 = buffer.data(gl + 437);
    const auto *gl_438 = buffer.data(gl + 438);
    const auto *gl_439 = buffer.data(gl + 439);
    const auto *gl_440 = buffer.data(gl + 440);
    const auto *gl_441 = buffer.data(gl + 441);
    const auto *gl_442 = buffer.data(gl + 442);
    const auto *gl_443 = buffer.data(gl + 443);
    const auto *gl_444 = buffer.data(gl + 444);
    const auto *gl_445 = buffer.data(gl + 445);
    const auto *gl_446 = buffer.data(gl + 446);
    const auto *gl_447 = buffer.data(gl + 447);
    const auto *gl_448 = buffer.data(gl + 448);
    const auto *gl_449 = buffer.data(gl + 449);

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, t_315, t_316, t_317, gl_310, \
                         gl_311, gl_312, gl_313, gl_314, gl_315, gl_316, \
                         gl_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_0 * gl_310[k];

        t_311[k] = f_0 * gl_311[k];

        t_312[k] = f_0 * gl_312[k];

        t_313[k] = f_0 * gl_313[k];

        t_314[k] = f_0 * gl_314[k];

        t_315[k] = f_0 * gl_315[k];

        t_316[k] = f_0 * gl_316[k];

        t_317[k] = f_0 * gl_317[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, t_323, t_324, t_325, gl_318, \
                         gl_319, gl_320, gl_321, gl_322, gl_323, gl_324, \
                         gl_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_0 * gl_318[k];

        t_319[k] = f_0 * gl_319[k];

        t_320[k] = f_0 * gl_320[k];

        t_321[k] = f_0 * gl_321[k];

        t_322[k] = f_0 * gl_322[k];

        t_323[k] = f_0 * gl_323[k];

        t_324[k] = f_0 * gl_324[k];

        t_325[k] = f_0 * gl_325[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, t_332, t_333, gl_326, \
                         gl_327, gl_328, gl_329, gl_330, gl_331, gl_332, \
                         gl_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_0 * gl_326[k];

        t_327[k] = f_0 * gl_327[k];

        t_328[k] = f_0 * gl_328[k];

        t_329[k] = f_0 * gl_329[k];

        t_330[k] = f_0 * gl_330[k];

        t_331[k] = f_0 * gl_331[k];

        t_332[k] = f_0 * gl_332[k];

        t_333[k] = f_0 * gl_333[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, t_339, t_340, t_341, gl_334, \
                         gl_335, gl_336, gl_337, gl_338, gl_339, gl_340, \
                         gl_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * gl_334[k];

        t_335[k] = f_0 * gl_335[k];

        t_336[k] = f_0 * gl_336[k];

        t_337[k] = f_0 * gl_337[k];

        t_338[k] = f_0 * gl_338[k];

        t_339[k] = f_0 * gl_339[k];

        t_340[k] = f_0 * gl_340[k];

        t_341[k] = f_0 * gl_341[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, t_347, t_348, t_349, gl_342, \
                         gl_343, gl_344, gl_345, gl_346, gl_347, gl_348, \
                         gl_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_0 * gl_342[k];

        t_343[k] = f_0 * gl_343[k];

        t_344[k] = f_0 * gl_344[k];

        t_345[k] = f_0 * gl_345[k];

        t_346[k] = f_0 * gl_346[k];

        t_347[k] = f_0 * gl_347[k];

        t_348[k] = f_0 * gl_348[k];

        t_349[k] = f_0 * gl_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, t_356, t_357, gl_350, \
                         gl_351, gl_352, gl_353, gl_354, gl_355, gl_356, \
                         gl_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_0 * gl_350[k];

        t_351[k] = f_0 * gl_351[k];

        t_352[k] = f_0 * gl_352[k];

        t_353[k] = f_0 * gl_353[k];

        t_354[k] = f_0 * gl_354[k];

        t_355[k] = f_0 * gl_355[k];

        t_356[k] = f_0 * gl_356[k];

        t_357[k] = f_0 * gl_357[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, t_364, t_365, gl_358, \
                         gl_359, gl_360, gl_361, gl_362, gl_363, gl_364, \
                         gl_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_0 * gl_358[k];

        t_359[k] = f_0 * gl_359[k];

        t_360[k] = f_0 * gl_360[k];

        t_361[k] = f_0 * gl_361[k];

        t_362[k] = f_0 * gl_362[k];

        t_363[k] = f_0 * gl_363[k];

        t_364[k] = f_0 * gl_364[k];

        t_365[k] = f_0 * gl_365[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, t_370, t_371, t_372, t_373, gl_366, \
                         gl_367, gl_368, gl_369, gl_370, gl_371, gl_372, \
                         gl_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_0 * gl_366[k];

        t_367[k] = f_0 * gl_367[k];

        t_368[k] = f_0 * gl_368[k];

        t_369[k] = f_0 * gl_369[k];

        t_370[k] = f_0 * gl_370[k];

        t_371[k] = f_0 * gl_371[k];

        t_372[k] = f_0 * gl_372[k];

        t_373[k] = f_0 * gl_373[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, t_378, t_379, t_380, t_381, gl_374, \
                         gl_375, gl_376, gl_377, gl_378, gl_379, gl_380, \
                         gl_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = f_0 * gl_374[k];

        t_375[k] = f_0 * gl_375[k];

        t_376[k] = f_0 * gl_376[k];

        t_377[k] = f_0 * gl_377[k];

        t_378[k] = f_0 * gl_378[k];

        t_379[k] = f_0 * gl_379[k];

        t_380[k] = f_0 * gl_380[k];

        t_381[k] = f_0 * gl_381[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, t_387, t_388, t_389, gl_382, \
                         gl_383, gl_384, gl_385, gl_386, gl_387, gl_388, \
                         gl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_0 * gl_382[k];

        t_383[k] = f_0 * gl_383[k];

        t_384[k] = f_0 * gl_384[k];

        t_385[k] = f_0 * gl_385[k];

        t_386[k] = f_0 * gl_386[k];

        t_387[k] = f_0 * gl_387[k];

        t_388[k] = f_0 * gl_388[k];

        t_389[k] = f_0 * gl_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, t_395, t_396, t_397, gl_390, \
                         gl_391, gl_392, gl_393, gl_394, gl_395, gl_396, \
                         gl_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_0 * gl_390[k];

        t_391[k] = f_0 * gl_391[k];

        t_392[k] = f_0 * gl_392[k];

        t_393[k] = f_0 * gl_393[k];

        t_394[k] = f_0 * gl_394[k];

        t_395[k] = f_0 * gl_395[k];

        t_396[k] = f_0 * gl_396[k];

        t_397[k] = f_0 * gl_397[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, t_402, t_403, t_404, t_405, gl_398, \
                         gl_399, gl_400, gl_401, gl_402, gl_403, gl_404, \
                         gl_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_0 * gl_398[k];

        t_399[k] = f_0 * gl_399[k];

        t_400[k] = f_0 * gl_400[k];

        t_401[k] = f_0 * gl_401[k];

        t_402[k] = f_0 * gl_402[k];

        t_403[k] = f_0 * gl_403[k];

        t_404[k] = f_0 * gl_404[k];

        t_405[k] = f_0 * gl_405[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, t_410, t_411, t_412, t_413, gl_406, \
                         gl_407, gl_408, gl_409, gl_410, gl_411, gl_412, \
                         gl_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = f_0 * gl_406[k];

        t_407[k] = f_0 * gl_407[k];

        t_408[k] = f_0 * gl_408[k];

        t_409[k] = f_0 * gl_409[k];

        t_410[k] = f_0 * gl_410[k];

        t_411[k] = f_0 * gl_411[k];

        t_412[k] = f_0 * gl_412[k];

        t_413[k] = f_0 * gl_413[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, t_419, t_420, t_421, gl_414, \
                         gl_415, gl_416, gl_417, gl_418, gl_419, gl_420, \
                         gl_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_0 * gl_414[k];

        t_415[k] = f_0 * gl_415[k];

        t_416[k] = f_0 * gl_416[k];

        t_417[k] = f_0 * gl_417[k];

        t_418[k] = f_0 * gl_418[k];

        t_419[k] = f_0 * gl_419[k];

        t_420[k] = f_0 * gl_420[k];

        t_421[k] = f_0 * gl_421[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, t_427, t_428, t_429, gl_422, \
                         gl_423, gl_424, gl_425, gl_426, gl_427, gl_428, \
                         gl_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_0 * gl_422[k];

        t_423[k] = f_0 * gl_423[k];

        t_424[k] = f_0 * gl_424[k];

        t_425[k] = f_0 * gl_425[k];

        t_426[k] = f_0 * gl_426[k];

        t_427[k] = f_0 * gl_427[k];

        t_428[k] = f_0 * gl_428[k];

        t_429[k] = f_0 * gl_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, t_435, t_436, t_437, gl_430, \
                         gl_431, gl_432, gl_433, gl_434, gl_435, gl_436, \
                         gl_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_0 * gl_430[k];

        t_431[k] = f_0 * gl_431[k];

        t_432[k] = f_0 * gl_432[k];

        t_433[k] = f_0 * gl_433[k];

        t_434[k] = f_0 * gl_434[k];

        t_435[k] = f_0 * gl_435[k];

        t_436[k] = f_0 * gl_436[k];

        t_437[k] = f_0 * gl_437[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, t_442, t_443, t_444, t_445, gl_438, \
                         gl_439, gl_440, gl_441, gl_442, gl_443, gl_444, \
                         gl_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_0 * gl_438[k];

        t_439[k] = f_0 * gl_439[k];

        t_440[k] = f_0 * gl_440[k];

        t_441[k] = f_0 * gl_441[k];

        t_442[k] = f_0 * gl_442[k];

        t_443[k] = f_0 * gl_443[k];

        t_444[k] = f_0 * gl_444[k];

        t_445[k] = f_0 * gl_445[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, gl_446, gl_447, gl_448, \
                         gl_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_0 * gl_446[k];

        t_447[k] = f_0 * gl_447[k];

        t_448[k] = f_0 * gl_448[k];

        t_449[k] = f_0 * gl_449[k];
    }
}

auto
compute_prim_geom_10_fl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t dl, const size_t gl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fl_electron_repulsion_0_piece0(buffer, target, dl, gl, ncols, alpha);

    compute_prim_geom_10_fl_electron_repulsion_0_piece1(buffer, target, dl, gl, ncols, alpha);

    compute_prim_geom_10_fl_electron_repulsion_0_piece2(buffer, target, gl, ncols, alpha);
}

static auto
compute_prim_geom_10_fl_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dl, const size_t gl,
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

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);
    const auto *dl_3 = buffer.data(dl + 3);
    const auto *dl_4 = buffer.data(dl + 4);
    const auto *dl_5 = buffer.data(dl + 5);
    const auto *dl_6 = buffer.data(dl + 6);
    const auto *dl_7 = buffer.data(dl + 7);
    const auto *dl_8 = buffer.data(dl + 8);
    const auto *dl_9 = buffer.data(dl + 9);
    const auto *dl_10 = buffer.data(dl + 10);
    const auto *dl_11 = buffer.data(dl + 11);
    const auto *dl_12 = buffer.data(dl + 12);
    const auto *dl_13 = buffer.data(dl + 13);
    const auto *dl_14 = buffer.data(dl + 14);
    const auto *dl_15 = buffer.data(dl + 15);
    const auto *dl_16 = buffer.data(dl + 16);
    const auto *dl_17 = buffer.data(dl + 17);
    const auto *dl_18 = buffer.data(dl + 18);
    const auto *dl_19 = buffer.data(dl + 19);
    const auto *dl_20 = buffer.data(dl + 20);
    const auto *dl_21 = buffer.data(dl + 21);
    const auto *dl_22 = buffer.data(dl + 22);
    const auto *dl_23 = buffer.data(dl + 23);
    const auto *dl_24 = buffer.data(dl + 24);
    const auto *dl_25 = buffer.data(dl + 25);
    const auto *dl_26 = buffer.data(dl + 26);
    const auto *dl_27 = buffer.data(dl + 27);
    const auto *dl_28 = buffer.data(dl + 28);
    const auto *dl_29 = buffer.data(dl + 29);
    const auto *dl_30 = buffer.data(dl + 30);
    const auto *dl_31 = buffer.data(dl + 31);
    const auto *dl_32 = buffer.data(dl + 32);
    const auto *dl_33 = buffer.data(dl + 33);
    const auto *dl_34 = buffer.data(dl + 34);
    const auto *dl_35 = buffer.data(dl + 35);
    const auto *dl_36 = buffer.data(dl + 36);
    const auto *dl_37 = buffer.data(dl + 37);
    const auto *dl_38 = buffer.data(dl + 38);
    const auto *dl_39 = buffer.data(dl + 39);
    const auto *dl_40 = buffer.data(dl + 40);
    const auto *dl_41 = buffer.data(dl + 41);
    const auto *dl_42 = buffer.data(dl + 42);
    const auto *dl_43 = buffer.data(dl + 43);
    const auto *dl_44 = buffer.data(dl + 44);
    const auto *dl_45 = buffer.data(dl + 45);
    const auto *dl_46 = buffer.data(dl + 46);
    const auto *dl_47 = buffer.data(dl + 47);
    const auto *dl_48 = buffer.data(dl + 48);
    const auto *dl_49 = buffer.data(dl + 49);
    const auto *dl_50 = buffer.data(dl + 50);
    const auto *dl_51 = buffer.data(dl + 51);
    const auto *dl_52 = buffer.data(dl + 52);
    const auto *dl_53 = buffer.data(dl + 53);
    const auto *dl_54 = buffer.data(dl + 54);
    const auto *dl_55 = buffer.data(dl + 55);
    const auto *dl_56 = buffer.data(dl + 56);
    const auto *dl_57 = buffer.data(dl + 57);
    const auto *dl_58 = buffer.data(dl + 58);
    const auto *dl_59 = buffer.data(dl + 59);
    const auto *dl_60 = buffer.data(dl + 60);
    const auto *dl_61 = buffer.data(dl + 61);
    const auto *dl_62 = buffer.data(dl + 62);
    const auto *dl_63 = buffer.data(dl + 63);
    const auto *dl_64 = buffer.data(dl + 64);
    const auto *dl_65 = buffer.data(dl + 65);
    const auto *dl_66 = buffer.data(dl + 66);
    const auto *dl_67 = buffer.data(dl + 67);
    const auto *dl_68 = buffer.data(dl + 68);
    const auto *dl_69 = buffer.data(dl + 69);
    const auto *dl_70 = buffer.data(dl + 70);
    const auto *dl_71 = buffer.data(dl + 71);
    const auto *dl_72 = buffer.data(dl + 72);
    const auto *dl_73 = buffer.data(dl + 73);
    const auto *dl_74 = buffer.data(dl + 74);
    const auto *dl_75 = buffer.data(dl + 75);
    const auto *dl_76 = buffer.data(dl + 76);
    const auto *dl_77 = buffer.data(dl + 77);
    const auto *dl_78 = buffer.data(dl + 78);
    const auto *dl_79 = buffer.data(dl + 79);
    const auto *dl_80 = buffer.data(dl + 80);
    const auto *dl_81 = buffer.data(dl + 81);
    const auto *dl_82 = buffer.data(dl + 82);
    const auto *dl_83 = buffer.data(dl + 83);
    const auto *dl_84 = buffer.data(dl + 84);
    const auto *dl_85 = buffer.data(dl + 85);
    const auto *dl_86 = buffer.data(dl + 86);
    const auto *dl_87 = buffer.data(dl + 87);
    const auto *dl_88 = buffer.data(dl + 88);

    const auto *gl_45 = buffer.data(gl + 45);
    const auto *gl_46 = buffer.data(gl + 46);
    const auto *gl_47 = buffer.data(gl + 47);
    const auto *gl_48 = buffer.data(gl + 48);
    const auto *gl_49 = buffer.data(gl + 49);
    const auto *gl_50 = buffer.data(gl + 50);
    const auto *gl_51 = buffer.data(gl + 51);
    const auto *gl_52 = buffer.data(gl + 52);
    const auto *gl_53 = buffer.data(gl + 53);
    const auto *gl_54 = buffer.data(gl + 54);
    const auto *gl_55 = buffer.data(gl + 55);
    const auto *gl_56 = buffer.data(gl + 56);
    const auto *gl_57 = buffer.data(gl + 57);
    const auto *gl_58 = buffer.data(gl + 58);
    const auto *gl_59 = buffer.data(gl + 59);
    const auto *gl_60 = buffer.data(gl + 60);
    const auto *gl_61 = buffer.data(gl + 61);
    const auto *gl_62 = buffer.data(gl + 62);
    const auto *gl_63 = buffer.data(gl + 63);
    const auto *gl_64 = buffer.data(gl + 64);
    const auto *gl_65 = buffer.data(gl + 65);
    const auto *gl_66 = buffer.data(gl + 66);
    const auto *gl_67 = buffer.data(gl + 67);
    const auto *gl_68 = buffer.data(gl + 68);
    const auto *gl_69 = buffer.data(gl + 69);
    const auto *gl_70 = buffer.data(gl + 70);
    const auto *gl_71 = buffer.data(gl + 71);
    const auto *gl_72 = buffer.data(gl + 72);
    const auto *gl_73 = buffer.data(gl + 73);
    const auto *gl_74 = buffer.data(gl + 74);
    const auto *gl_75 = buffer.data(gl + 75);
    const auto *gl_76 = buffer.data(gl + 76);
    const auto *gl_77 = buffer.data(gl + 77);
    const auto *gl_78 = buffer.data(gl + 78);
    const auto *gl_79 = buffer.data(gl + 79);
    const auto *gl_80 = buffer.data(gl + 80);
    const auto *gl_81 = buffer.data(gl + 81);
    const auto *gl_82 = buffer.data(gl + 82);
    const auto *gl_83 = buffer.data(gl + 83);
    const auto *gl_84 = buffer.data(gl + 84);
    const auto *gl_85 = buffer.data(gl + 85);
    const auto *gl_86 = buffer.data(gl + 86);
    const auto *gl_87 = buffer.data(gl + 87);
    const auto *gl_88 = buffer.data(gl + 88);
    const auto *gl_89 = buffer.data(gl + 89);
    const auto *gl_135 = buffer.data(gl + 135);
    const auto *gl_136 = buffer.data(gl + 136);
    const auto *gl_137 = buffer.data(gl + 137);
    const auto *gl_138 = buffer.data(gl + 138);
    const auto *gl_139 = buffer.data(gl + 139);
    const auto *gl_140 = buffer.data(gl + 140);
    const auto *gl_141 = buffer.data(gl + 141);
    const auto *gl_142 = buffer.data(gl + 142);
    const auto *gl_143 = buffer.data(gl + 143);
    const auto *gl_144 = buffer.data(gl + 144);
    const auto *gl_145 = buffer.data(gl + 145);
    const auto *gl_146 = buffer.data(gl + 146);
    const auto *gl_147 = buffer.data(gl + 147);
    const auto *gl_148 = buffer.data(gl + 148);
    const auto *gl_149 = buffer.data(gl + 149);
    const auto *gl_150 = buffer.data(gl + 150);
    const auto *gl_151 = buffer.data(gl + 151);
    const auto *gl_152 = buffer.data(gl + 152);
    const auto *gl_153 = buffer.data(gl + 153);
    const auto *gl_154 = buffer.data(gl + 154);
    const auto *gl_155 = buffer.data(gl + 155);
    const auto *gl_156 = buffer.data(gl + 156);
    const auto *gl_157 = buffer.data(gl + 157);
    const auto *gl_158 = buffer.data(gl + 158);
    const auto *gl_159 = buffer.data(gl + 159);
    const auto *gl_160 = buffer.data(gl + 160);
    const auto *gl_161 = buffer.data(gl + 161);
    const auto *gl_162 = buffer.data(gl + 162);
    const auto *gl_163 = buffer.data(gl + 163);
    const auto *gl_164 = buffer.data(gl + 164);
    const auto *gl_165 = buffer.data(gl + 165);
    const auto *gl_166 = buffer.data(gl + 166);
    const auto *gl_167 = buffer.data(gl + 167);
    const auto *gl_168 = buffer.data(gl + 168);
    const auto *gl_169 = buffer.data(gl + 169);
    const auto *gl_170 = buffer.data(gl + 170);
    const auto *gl_171 = buffer.data(gl + 171);
    const auto *gl_172 = buffer.data(gl + 172);
    const auto *gl_173 = buffer.data(gl + 173);
    const auto *gl_174 = buffer.data(gl + 174);
    const auto *gl_175 = buffer.data(gl + 175);
    const auto *gl_176 = buffer.data(gl + 176);
    const auto *gl_177 = buffer.data(gl + 177);
    const auto *gl_178 = buffer.data(gl + 178);
    const auto *gl_179 = buffer.data(gl + 179);
    const auto *gl_180 = buffer.data(gl + 180);
    const auto *gl_181 = buffer.data(gl + 181);
    const auto *gl_182 = buffer.data(gl + 182);
    const auto *gl_183 = buffer.data(gl + 183);
    const auto *gl_184 = buffer.data(gl + 184);
    const auto *gl_185 = buffer.data(gl + 185);
    const auto *gl_186 = buffer.data(gl + 186);
    const auto *gl_187 = buffer.data(gl + 187);
    const auto *gl_188 = buffer.data(gl + 188);
    const auto *gl_189 = buffer.data(gl + 189);
    const auto *gl_190 = buffer.data(gl + 190);
    const auto *gl_191 = buffer.data(gl + 191);
    const auto *gl_192 = buffer.data(gl + 192);
    const auto *gl_193 = buffer.data(gl + 193);
    const auto *gl_194 = buffer.data(gl + 194);
    const auto *gl_195 = buffer.data(gl + 195);
    const auto *gl_196 = buffer.data(gl + 196);
    const auto *gl_197 = buffer.data(gl + 197);
    const auto *gl_198 = buffer.data(gl + 198);
    const auto *gl_199 = buffer.data(gl + 199);
    const auto *gl_200 = buffer.data(gl + 200);
    const auto *gl_201 = buffer.data(gl + 201);
    const auto *gl_202 = buffer.data(gl + 202);
    const auto *gl_203 = buffer.data(gl + 203);
    const auto *gl_204 = buffer.data(gl + 204);
    const auto *gl_205 = buffer.data(gl + 205);
    const auto *gl_206 = buffer.data(gl + 206);
    const auto *gl_207 = buffer.data(gl + 207);
    const auto *gl_208 = buffer.data(gl + 208);
    const auto *gl_209 = buffer.data(gl + 209);
    const auto *gl_210 = buffer.data(gl + 210);
    const auto *gl_211 = buffer.data(gl + 211);
    const auto *gl_212 = buffer.data(gl + 212);
    const auto *gl_213 = buffer.data(gl + 213);
    const auto *gl_214 = buffer.data(gl + 214);
    const auto *gl_215 = buffer.data(gl + 215);
    const auto *gl_216 = buffer.data(gl + 216);
    const auto *gl_217 = buffer.data(gl + 217);
    const auto *gl_218 = buffer.data(gl + 218);
    const auto *gl_219 = buffer.data(gl + 219);
    const auto *gl_220 = buffer.data(gl + 220);
    const auto *gl_221 = buffer.data(gl + 221);
    const auto *gl_222 = buffer.data(gl + 222);
    const auto *gl_223 = buffer.data(gl + 223);
    const auto *gl_224 = buffer.data(gl + 224);
    const auto *gl_270 = buffer.data(gl + 270);
    const auto *gl_271 = buffer.data(gl + 271);
    const auto *gl_272 = buffer.data(gl + 272);
    const auto *gl_273 = buffer.data(gl + 273);
    const auto *gl_274 = buffer.data(gl + 274);
    const auto *gl_275 = buffer.data(gl + 275);
    const auto *gl_276 = buffer.data(gl + 276);
    const auto *gl_277 = buffer.data(gl + 277);
    const auto *gl_278 = buffer.data(gl + 278);
    const auto *gl_279 = buffer.data(gl + 279);
    const auto *gl_280 = buffer.data(gl + 280);
    const auto *gl_281 = buffer.data(gl + 281);
    const auto *gl_282 = buffer.data(gl + 282);
    const auto *gl_283 = buffer.data(gl + 283);
    const auto *gl_284 = buffer.data(gl + 284);
    const auto *gl_285 = buffer.data(gl + 285);
    const auto *gl_286 = buffer.data(gl + 286);
    const auto *gl_287 = buffer.data(gl + 287);
    const auto *gl_288 = buffer.data(gl + 288);
    const auto *gl_289 = buffer.data(gl + 289);
    const auto *gl_290 = buffer.data(gl + 290);
    const auto *gl_291 = buffer.data(gl + 291);
    const auto *gl_292 = buffer.data(gl + 292);
    const auto *gl_293 = buffer.data(gl + 293);
    const auto *gl_294 = buffer.data(gl + 294);
    const auto *gl_295 = buffer.data(gl + 295);
    const auto *gl_296 = buffer.data(gl + 296);
    const auto *gl_297 = buffer.data(gl + 297);
    const auto *gl_298 = buffer.data(gl + 298);
    const auto *gl_299 = buffer.data(gl + 299);
    const auto *gl_300 = buffer.data(gl + 300);
    const auto *gl_301 = buffer.data(gl + 301);
    const auto *gl_302 = buffer.data(gl + 302);
    const auto *gl_303 = buffer.data(gl + 303);
    const auto *gl_304 = buffer.data(gl + 304);
    const auto *gl_305 = buffer.data(gl + 305);
    const auto *gl_306 = buffer.data(gl + 306);
    const auto *gl_307 = buffer.data(gl + 307);
    const auto *gl_308 = buffer.data(gl + 308);
    const auto *gl_309 = buffer.data(gl + 309);
    const auto *gl_310 = buffer.data(gl + 310);
    const auto *gl_311 = buffer.data(gl + 311);
    const auto *gl_312 = buffer.data(gl + 312);
    const auto *gl_313 = buffer.data(gl + 313);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gl_45, gl_46, gl_47, gl_48, \
                         gl_49, gl_50, gl_51, gl_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gl_45[k];

        t_1[k] = f_0 * gl_46[k];

        t_2[k] = f_0 * gl_47[k];

        t_3[k] = f_0 * gl_48[k];

        t_4[k] = f_0 * gl_49[k];

        t_5[k] = f_0 * gl_50[k];

        t_6[k] = f_0 * gl_51[k];

        t_7[k] = f_0 * gl_52[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gl_53, gl_54, gl_55, \
                         gl_56, gl_57, gl_58, gl_59, gl_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gl_53[k];

        t_9[k] = f_0 * gl_54[k];

        t_10[k] = f_0 * gl_55[k];

        t_11[k] = f_0 * gl_56[k];

        t_12[k] = f_0 * gl_57[k];

        t_13[k] = f_0 * gl_58[k];

        t_14[k] = f_0 * gl_59[k];

        t_15[k] = f_0 * gl_60[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, gl_61, gl_62, gl_63, \
                         gl_64, gl_65, gl_66, gl_67, gl_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gl_61[k];

        t_17[k] = f_0 * gl_62[k];

        t_18[k] = f_0 * gl_63[k];

        t_19[k] = f_0 * gl_64[k];

        t_20[k] = f_0 * gl_65[k];

        t_21[k] = f_0 * gl_66[k];

        t_22[k] = f_0 * gl_67[k];

        t_23[k] = f_0 * gl_68[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, gl_69, gl_70, gl_71, \
                         gl_72, gl_73, gl_74, gl_75, gl_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * gl_69[k];

        t_25[k] = f_0 * gl_70[k];

        t_26[k] = f_0 * gl_71[k];

        t_27[k] = f_0 * gl_72[k];

        t_28[k] = f_0 * gl_73[k];

        t_29[k] = f_0 * gl_74[k];

        t_30[k] = f_0 * gl_75[k];

        t_31[k] = f_0 * gl_76[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, gl_77, gl_78, gl_79, \
                         gl_80, gl_81, gl_82, gl_83, gl_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * gl_77[k];

        t_33[k] = f_0 * gl_78[k];

        t_34[k] = f_0 * gl_79[k];

        t_35[k] = f_0 * gl_80[k];

        t_36[k] = f_0 * gl_81[k];

        t_37[k] = f_0 * gl_82[k];

        t_38[k] = f_0 * gl_83[k];

        t_39[k] = f_0 * gl_84[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, dl_0, dl_1, gl_85, gl_86, \
                         gl_87, gl_88, gl_89, gl_135, gl_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * gl_85[k];

        t_41[k] = f_0 * gl_86[k];

        t_42[k] = f_0 * gl_87[k];

        t_43[k] = f_0 * gl_88[k];

        t_44[k] = f_0 * gl_89[k];

        t_45[k] = -dl_0[k]
                  + f_0 * gl_135[k];

        t_46[k] = -dl_1[k]
                  + f_0 * gl_136[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, dl_2, dl_3, dl_4, dl_5, dl_6, gl_137, \
                         gl_138, gl_139, gl_140, gl_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -dl_2[k]
                  + f_0 * gl_137[k];

        t_48[k] = -dl_3[k]
                  + f_0 * gl_138[k];

        t_49[k] = -dl_4[k]
                  + f_0 * gl_139[k];

        t_50[k] = -dl_5[k]
                  + f_0 * gl_140[k];

        t_51[k] = -dl_6[k]
                  + f_0 * gl_141[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, dl_7, dl_8, dl_9, dl_10, dl_11, gl_142, \
                         gl_143, gl_144, gl_145, gl_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -dl_7[k]
                  + f_0 * gl_142[k];

        t_53[k] = -dl_8[k]
                  + f_0 * gl_143[k];

        t_54[k] = -dl_9[k]
                  + f_0 * gl_144[k];

        t_55[k] = -dl_10[k]
                  + f_0 * gl_145[k];

        t_56[k] = -dl_11[k]
                  + f_0 * gl_146[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, dl_12, dl_13, dl_14, dl_15, dl_16, \
                         gl_147, gl_148, gl_149, gl_150, gl_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -dl_12[k]
                  + f_0 * gl_147[k];

        t_58[k] = -dl_13[k]
                  + f_0 * gl_148[k];

        t_59[k] = -dl_14[k]
                  + f_0 * gl_149[k];

        t_60[k] = -dl_15[k]
                  + f_0 * gl_150[k];

        t_61[k] = -dl_16[k]
                  + f_0 * gl_151[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, dl_17, dl_18, dl_19, dl_20, dl_21, \
                         gl_152, gl_153, gl_154, gl_155, gl_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -dl_17[k]
                  + f_0 * gl_152[k];

        t_63[k] = -dl_18[k]
                  + f_0 * gl_153[k];

        t_64[k] = -dl_19[k]
                  + f_0 * gl_154[k];

        t_65[k] = -dl_20[k]
                  + f_0 * gl_155[k];

        t_66[k] = -dl_21[k]
                  + f_0 * gl_156[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, dl_22, dl_23, dl_24, dl_25, dl_26, \
                         gl_157, gl_158, gl_159, gl_160, gl_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -dl_22[k]
                  + f_0 * gl_157[k];

        t_68[k] = -dl_23[k]
                  + f_0 * gl_158[k];

        t_69[k] = -dl_24[k]
                  + f_0 * gl_159[k];

        t_70[k] = -dl_25[k]
                  + f_0 * gl_160[k];

        t_71[k] = -dl_26[k]
                  + f_0 * gl_161[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, dl_27, dl_28, dl_29, dl_30, dl_31, \
                         gl_162, gl_163, gl_164, gl_165, gl_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -dl_27[k]
                  + f_0 * gl_162[k];

        t_73[k] = -dl_28[k]
                  + f_0 * gl_163[k];

        t_74[k] = -dl_29[k]
                  + f_0 * gl_164[k];

        t_75[k] = -dl_30[k]
                  + f_0 * gl_165[k];

        t_76[k] = -dl_31[k]
                  + f_0 * gl_166[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, dl_32, dl_33, dl_34, dl_35, dl_36, \
                         gl_167, gl_168, gl_169, gl_170, gl_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -dl_32[k]
                  + f_0 * gl_167[k];

        t_78[k] = -dl_33[k]
                  + f_0 * gl_168[k];

        t_79[k] = -dl_34[k]
                  + f_0 * gl_169[k];

        t_80[k] = -dl_35[k]
                  + f_0 * gl_170[k];

        t_81[k] = -dl_36[k]
                  + f_0 * gl_171[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, dl_37, dl_38, dl_39, dl_40, dl_41, \
                         gl_172, gl_173, gl_174, gl_175, gl_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -dl_37[k]
                  + f_0 * gl_172[k];

        t_83[k] = -dl_38[k]
                  + f_0 * gl_173[k];

        t_84[k] = -dl_39[k]
                  + f_0 * gl_174[k];

        t_85[k] = -dl_40[k]
                  + f_0 * gl_175[k];

        t_86[k] = -dl_41[k]
                  + f_0 * gl_176[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, dl_42, dl_43, dl_44, gl_177, \
                         gl_178, gl_179, gl_180, gl_181, gl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -dl_42[k]
                  + f_0 * gl_177[k];

        t_88[k] = -dl_43[k]
                  + f_0 * gl_178[k];

        t_89[k] = -dl_44[k]
                  + f_0 * gl_179[k];

        t_90[k] = f_0 * gl_180[k];

        t_91[k] = f_0 * gl_181[k];

        t_92[k] = f_0 * gl_182[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, t_99, t_100, gl_183, gl_184, \
                         gl_185, gl_186, gl_187, gl_188, gl_189, \
                         gl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_0 * gl_183[k];

        t_94[k] = f_0 * gl_184[k];

        t_95[k] = f_0 * gl_185[k];

        t_96[k] = f_0 * gl_186[k];

        t_97[k] = f_0 * gl_187[k];

        t_98[k] = f_0 * gl_188[k];

        t_99[k] = f_0 * gl_189[k];

        t_100[k] = f_0 * gl_190[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, t_107, t_108, gl_191, \
                         gl_192, gl_193, gl_194, gl_195, gl_196, gl_197, \
                         gl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * gl_191[k];

        t_102[k] = f_0 * gl_192[k];

        t_103[k] = f_0 * gl_193[k];

        t_104[k] = f_0 * gl_194[k];

        t_105[k] = f_0 * gl_195[k];

        t_106[k] = f_0 * gl_196[k];

        t_107[k] = f_0 * gl_197[k];

        t_108[k] = f_0 * gl_198[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, t_114, t_115, t_116, gl_199, \
                         gl_200, gl_201, gl_202, gl_203, gl_204, gl_205, \
                         gl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_0 * gl_199[k];

        t_110[k] = f_0 * gl_200[k];

        t_111[k] = f_0 * gl_201[k];

        t_112[k] = f_0 * gl_202[k];

        t_113[k] = f_0 * gl_203[k];

        t_114[k] = f_0 * gl_204[k];

        t_115[k] = f_0 * gl_205[k];

        t_116[k] = f_0 * gl_206[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, t_122, t_123, t_124, gl_207, \
                         gl_208, gl_209, gl_210, gl_211, gl_212, gl_213, \
                         gl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * gl_207[k];

        t_118[k] = f_0 * gl_208[k];

        t_119[k] = f_0 * gl_209[k];

        t_120[k] = f_0 * gl_210[k];

        t_121[k] = f_0 * gl_211[k];

        t_122[k] = f_0 * gl_212[k];

        t_123[k] = f_0 * gl_213[k];

        t_124[k] = f_0 * gl_214[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, t_131, t_132, gl_215, \
                         gl_216, gl_217, gl_218, gl_219, gl_220, gl_221, \
                         gl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * gl_215[k];

        t_126[k] = f_0 * gl_216[k];

        t_127[k] = f_0 * gl_217[k];

        t_128[k] = f_0 * gl_218[k];

        t_129[k] = f_0 * gl_219[k];

        t_130[k] = f_0 * gl_220[k];

        t_131[k] = f_0 * gl_221[k];

        t_132[k] = f_0 * gl_222[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, t_138, dl_45, dl_46, dl_47, dl_48, \
                         gl_223, gl_224, gl_270, gl_271, gl_272, \
                         gl_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_0 * gl_223[k];

        t_134[k] = f_0 * gl_224[k];

        t_135[k] = -2.0 * dl_45[k]
                   + f_0 * gl_270[k];

        t_136[k] = -2.0 * dl_46[k]
                   + f_0 * gl_271[k];

        t_137[k] = -2.0 * dl_47[k]
                   + f_0 * gl_272[k];

        t_138[k] = -2.0 * dl_48[k]
                   + f_0 * gl_273[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, dl_49, dl_50, dl_51, dl_52, dl_53, \
                         gl_274, gl_275, gl_276, gl_277, gl_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -2.0 * dl_49[k]
                   + f_0 * gl_274[k];

        t_140[k] = -2.0 * dl_50[k]
                   + f_0 * gl_275[k];

        t_141[k] = -2.0 * dl_51[k]
                   + f_0 * gl_276[k];

        t_142[k] = -2.0 * dl_52[k]
                   + f_0 * gl_277[k];

        t_143[k] = -2.0 * dl_53[k]
                   + f_0 * gl_278[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, dl_54, dl_55, dl_56, dl_57, dl_58, \
                         gl_279, gl_280, gl_281, gl_282, gl_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -2.0 * dl_54[k]
                   + f_0 * gl_279[k];

        t_145[k] = -2.0 * dl_55[k]
                   + f_0 * gl_280[k];

        t_146[k] = -2.0 * dl_56[k]
                   + f_0 * gl_281[k];

        t_147[k] = -2.0 * dl_57[k]
                   + f_0 * gl_282[k];

        t_148[k] = -2.0 * dl_58[k]
                   + f_0 * gl_283[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, dl_59, dl_60, dl_61, dl_62, dl_63, \
                         gl_284, gl_285, gl_286, gl_287, gl_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -2.0 * dl_59[k]
                   + f_0 * gl_284[k];

        t_150[k] = -2.0 * dl_60[k]
                   + f_0 * gl_285[k];

        t_151[k] = -2.0 * dl_61[k]
                   + f_0 * gl_286[k];

        t_152[k] = -2.0 * dl_62[k]
                   + f_0 * gl_287[k];

        t_153[k] = -2.0 * dl_63[k]
                   + f_0 * gl_288[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, dl_64, dl_65, dl_66, dl_67, dl_68, \
                         gl_289, gl_290, gl_291, gl_292, gl_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -2.0 * dl_64[k]
                   + f_0 * gl_289[k];

        t_155[k] = -2.0 * dl_65[k]
                   + f_0 * gl_290[k];

        t_156[k] = -2.0 * dl_66[k]
                   + f_0 * gl_291[k];

        t_157[k] = -2.0 * dl_67[k]
                   + f_0 * gl_292[k];

        t_158[k] = -2.0 * dl_68[k]
                   + f_0 * gl_293[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, dl_69, dl_70, dl_71, dl_72, dl_73, \
                         gl_294, gl_295, gl_296, gl_297, gl_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -2.0 * dl_69[k]
                   + f_0 * gl_294[k];

        t_160[k] = -2.0 * dl_70[k]
                   + f_0 * gl_295[k];

        t_161[k] = -2.0 * dl_71[k]
                   + f_0 * gl_296[k];

        t_162[k] = -2.0 * dl_72[k]
                   + f_0 * gl_297[k];

        t_163[k] = -2.0 * dl_73[k]
                   + f_0 * gl_298[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, dl_74, dl_75, dl_76, dl_77, dl_78, \
                         gl_299, gl_300, gl_301, gl_302, gl_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -2.0 * dl_74[k]
                   + f_0 * gl_299[k];

        t_165[k] = -2.0 * dl_75[k]
                   + f_0 * gl_300[k];

        t_166[k] = -2.0 * dl_76[k]
                   + f_0 * gl_301[k];

        t_167[k] = -2.0 * dl_77[k]
                   + f_0 * gl_302[k];

        t_168[k] = -2.0 * dl_78[k]
                   + f_0 * gl_303[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, dl_79, dl_80, dl_81, dl_82, dl_83, \
                         gl_304, gl_305, gl_306, gl_307, gl_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -2.0 * dl_79[k]
                   + f_0 * gl_304[k];

        t_170[k] = -2.0 * dl_80[k]
                   + f_0 * gl_305[k];

        t_171[k] = -2.0 * dl_81[k]
                   + f_0 * gl_306[k];

        t_172[k] = -2.0 * dl_82[k]
                   + f_0 * gl_307[k];

        t_173[k] = -2.0 * dl_83[k]
                   + f_0 * gl_308[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, dl_84, dl_85, dl_86, dl_87, dl_88, \
                         gl_309, gl_310, gl_311, gl_312, gl_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -2.0 * dl_84[k]
                   + f_0 * gl_309[k];

        t_175[k] = -2.0 * dl_85[k]
                   + f_0 * gl_310[k];

        t_176[k] = -2.0 * dl_86[k]
                   + f_0 * gl_311[k];

        t_177[k] = -2.0 * dl_87[k]
                   + f_0 * gl_312[k];

        t_178[k] = -2.0 * dl_88[k]
                   + f_0 * gl_313[k];
    }
}

static auto
compute_prim_geom_10_fl_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dl, const size_t gl,
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

    const auto *dl_89 = buffer.data(dl + 89);
    const auto *dl_90 = buffer.data(dl + 90);
    const auto *dl_91 = buffer.data(dl + 91);
    const auto *dl_92 = buffer.data(dl + 92);
    const auto *dl_93 = buffer.data(dl + 93);
    const auto *dl_94 = buffer.data(dl + 94);
    const auto *dl_95 = buffer.data(dl + 95);
    const auto *dl_96 = buffer.data(dl + 96);
    const auto *dl_97 = buffer.data(dl + 97);
    const auto *dl_98 = buffer.data(dl + 98);
    const auto *dl_99 = buffer.data(dl + 99);
    const auto *dl_100 = buffer.data(dl + 100);
    const auto *dl_101 = buffer.data(dl + 101);
    const auto *dl_102 = buffer.data(dl + 102);
    const auto *dl_103 = buffer.data(dl + 103);
    const auto *dl_104 = buffer.data(dl + 104);
    const auto *dl_105 = buffer.data(dl + 105);
    const auto *dl_106 = buffer.data(dl + 106);
    const auto *dl_107 = buffer.data(dl + 107);
    const auto *dl_108 = buffer.data(dl + 108);
    const auto *dl_109 = buffer.data(dl + 109);
    const auto *dl_110 = buffer.data(dl + 110);
    const auto *dl_111 = buffer.data(dl + 111);
    const auto *dl_112 = buffer.data(dl + 112);
    const auto *dl_113 = buffer.data(dl + 113);
    const auto *dl_114 = buffer.data(dl + 114);
    const auto *dl_115 = buffer.data(dl + 115);
    const auto *dl_116 = buffer.data(dl + 116);
    const auto *dl_117 = buffer.data(dl + 117);
    const auto *dl_118 = buffer.data(dl + 118);
    const auto *dl_119 = buffer.data(dl + 119);
    const auto *dl_120 = buffer.data(dl + 120);
    const auto *dl_121 = buffer.data(dl + 121);
    const auto *dl_122 = buffer.data(dl + 122);
    const auto *dl_123 = buffer.data(dl + 123);
    const auto *dl_124 = buffer.data(dl + 124);
    const auto *dl_125 = buffer.data(dl + 125);
    const auto *dl_126 = buffer.data(dl + 126);
    const auto *dl_127 = buffer.data(dl + 127);
    const auto *dl_128 = buffer.data(dl + 128);
    const auto *dl_129 = buffer.data(dl + 129);
    const auto *dl_130 = buffer.data(dl + 130);
    const auto *dl_131 = buffer.data(dl + 131);
    const auto *dl_132 = buffer.data(dl + 132);
    const auto *dl_133 = buffer.data(dl + 133);
    const auto *dl_134 = buffer.data(dl + 134);
    const auto *dl_135 = buffer.data(dl + 135);
    const auto *dl_136 = buffer.data(dl + 136);
    const auto *dl_137 = buffer.data(dl + 137);
    const auto *dl_138 = buffer.data(dl + 138);
    const auto *dl_139 = buffer.data(dl + 139);
    const auto *dl_140 = buffer.data(dl + 140);
    const auto *dl_141 = buffer.data(dl + 141);
    const auto *dl_142 = buffer.data(dl + 142);
    const auto *dl_143 = buffer.data(dl + 143);
    const auto *dl_144 = buffer.data(dl + 144);
    const auto *dl_145 = buffer.data(dl + 145);
    const auto *dl_146 = buffer.data(dl + 146);
    const auto *dl_147 = buffer.data(dl + 147);
    const auto *dl_148 = buffer.data(dl + 148);
    const auto *dl_149 = buffer.data(dl + 149);
    const auto *dl_150 = buffer.data(dl + 150);
    const auto *dl_151 = buffer.data(dl + 151);
    const auto *dl_152 = buffer.data(dl + 152);
    const auto *dl_153 = buffer.data(dl + 153);
    const auto *dl_154 = buffer.data(dl + 154);
    const auto *dl_155 = buffer.data(dl + 155);
    const auto *dl_156 = buffer.data(dl + 156);
    const auto *dl_157 = buffer.data(dl + 157);
    const auto *dl_158 = buffer.data(dl + 158);
    const auto *dl_159 = buffer.data(dl + 159);
    const auto *dl_160 = buffer.data(dl + 160);
    const auto *dl_161 = buffer.data(dl + 161);
    const auto *dl_162 = buffer.data(dl + 162);
    const auto *dl_163 = buffer.data(dl + 163);
    const auto *dl_164 = buffer.data(dl + 164);
    const auto *dl_165 = buffer.data(dl + 165);
    const auto *dl_166 = buffer.data(dl + 166);
    const auto *dl_167 = buffer.data(dl + 167);
    const auto *dl_168 = buffer.data(dl + 168);
    const auto *dl_169 = buffer.data(dl + 169);
    const auto *dl_170 = buffer.data(dl + 170);
    const auto *dl_171 = buffer.data(dl + 171);
    const auto *dl_172 = buffer.data(dl + 172);
    const auto *dl_173 = buffer.data(dl + 173);
    const auto *dl_174 = buffer.data(dl + 174);
    const auto *dl_175 = buffer.data(dl + 175);
    const auto *dl_176 = buffer.data(dl + 176);
    const auto *dl_177 = buffer.data(dl + 177);
    const auto *dl_178 = buffer.data(dl + 178);
    const auto *dl_179 = buffer.data(dl + 179);
    const auto *dl_180 = buffer.data(dl + 180);
    const auto *dl_181 = buffer.data(dl + 181);
    const auto *dl_182 = buffer.data(dl + 182);
    const auto *dl_183 = buffer.data(dl + 183);
    const auto *dl_184 = buffer.data(dl + 184);
    const auto *dl_185 = buffer.data(dl + 185);
    const auto *dl_186 = buffer.data(dl + 186);
    const auto *dl_187 = buffer.data(dl + 187);
    const auto *dl_188 = buffer.data(dl + 188);
    const auto *dl_189 = buffer.data(dl + 189);
    const auto *dl_190 = buffer.data(dl + 190);
    const auto *dl_191 = buffer.data(dl + 191);
    const auto *dl_192 = buffer.data(dl + 192);
    const auto *dl_193 = buffer.data(dl + 193);
    const auto *dl_194 = buffer.data(dl + 194);
    const auto *dl_195 = buffer.data(dl + 195);
    const auto *dl_196 = buffer.data(dl + 196);
    const auto *dl_197 = buffer.data(dl + 197);
    const auto *dl_198 = buffer.data(dl + 198);
    const auto *dl_199 = buffer.data(dl + 199);
    const auto *dl_200 = buffer.data(dl + 200);
    const auto *dl_201 = buffer.data(dl + 201);
    const auto *dl_202 = buffer.data(dl + 202);
    const auto *dl_203 = buffer.data(dl + 203);
    const auto *dl_204 = buffer.data(dl + 204);
    const auto *dl_205 = buffer.data(dl + 205);
    const auto *dl_206 = buffer.data(dl + 206);
    const auto *dl_207 = buffer.data(dl + 207);
    const auto *dl_208 = buffer.data(dl + 208);
    const auto *dl_209 = buffer.data(dl + 209);

    const auto *gl_314 = buffer.data(gl + 314);
    const auto *gl_315 = buffer.data(gl + 315);
    const auto *gl_316 = buffer.data(gl + 316);
    const auto *gl_317 = buffer.data(gl + 317);
    const auto *gl_318 = buffer.data(gl + 318);
    const auto *gl_319 = buffer.data(gl + 319);
    const auto *gl_320 = buffer.data(gl + 320);
    const auto *gl_321 = buffer.data(gl + 321);
    const auto *gl_322 = buffer.data(gl + 322);
    const auto *gl_323 = buffer.data(gl + 323);
    const auto *gl_324 = buffer.data(gl + 324);
    const auto *gl_325 = buffer.data(gl + 325);
    const auto *gl_326 = buffer.data(gl + 326);
    const auto *gl_327 = buffer.data(gl + 327);
    const auto *gl_328 = buffer.data(gl + 328);
    const auto *gl_329 = buffer.data(gl + 329);
    const auto *gl_330 = buffer.data(gl + 330);
    const auto *gl_331 = buffer.data(gl + 331);
    const auto *gl_332 = buffer.data(gl + 332);
    const auto *gl_333 = buffer.data(gl + 333);
    const auto *gl_334 = buffer.data(gl + 334);
    const auto *gl_335 = buffer.data(gl + 335);
    const auto *gl_336 = buffer.data(gl + 336);
    const auto *gl_337 = buffer.data(gl + 337);
    const auto *gl_338 = buffer.data(gl + 338);
    const auto *gl_339 = buffer.data(gl + 339);
    const auto *gl_340 = buffer.data(gl + 340);
    const auto *gl_341 = buffer.data(gl + 341);
    const auto *gl_342 = buffer.data(gl + 342);
    const auto *gl_343 = buffer.data(gl + 343);
    const auto *gl_344 = buffer.data(gl + 344);
    const auto *gl_345 = buffer.data(gl + 345);
    const auto *gl_346 = buffer.data(gl + 346);
    const auto *gl_347 = buffer.data(gl + 347);
    const auto *gl_348 = buffer.data(gl + 348);
    const auto *gl_349 = buffer.data(gl + 349);
    const auto *gl_350 = buffer.data(gl + 350);
    const auto *gl_351 = buffer.data(gl + 351);
    const auto *gl_352 = buffer.data(gl + 352);
    const auto *gl_353 = buffer.data(gl + 353);
    const auto *gl_354 = buffer.data(gl + 354);
    const auto *gl_355 = buffer.data(gl + 355);
    const auto *gl_356 = buffer.data(gl + 356);
    const auto *gl_357 = buffer.data(gl + 357);
    const auto *gl_358 = buffer.data(gl + 358);
    const auto *gl_359 = buffer.data(gl + 359);
    const auto *gl_360 = buffer.data(gl + 360);
    const auto *gl_361 = buffer.data(gl + 361);
    const auto *gl_362 = buffer.data(gl + 362);
    const auto *gl_363 = buffer.data(gl + 363);
    const auto *gl_364 = buffer.data(gl + 364);
    const auto *gl_365 = buffer.data(gl + 365);
    const auto *gl_366 = buffer.data(gl + 366);
    const auto *gl_367 = buffer.data(gl + 367);
    const auto *gl_368 = buffer.data(gl + 368);
    const auto *gl_369 = buffer.data(gl + 369);
    const auto *gl_370 = buffer.data(gl + 370);
    const auto *gl_371 = buffer.data(gl + 371);
    const auto *gl_372 = buffer.data(gl + 372);
    const auto *gl_373 = buffer.data(gl + 373);
    const auto *gl_374 = buffer.data(gl + 374);
    const auto *gl_375 = buffer.data(gl + 375);
    const auto *gl_376 = buffer.data(gl + 376);
    const auto *gl_377 = buffer.data(gl + 377);
    const auto *gl_378 = buffer.data(gl + 378);
    const auto *gl_379 = buffer.data(gl + 379);
    const auto *gl_380 = buffer.data(gl + 380);
    const auto *gl_381 = buffer.data(gl + 381);
    const auto *gl_382 = buffer.data(gl + 382);
    const auto *gl_383 = buffer.data(gl + 383);
    const auto *gl_384 = buffer.data(gl + 384);
    const auto *gl_385 = buffer.data(gl + 385);
    const auto *gl_386 = buffer.data(gl + 386);
    const auto *gl_387 = buffer.data(gl + 387);
    const auto *gl_388 = buffer.data(gl + 388);
    const auto *gl_389 = buffer.data(gl + 389);
    const auto *gl_390 = buffer.data(gl + 390);
    const auto *gl_391 = buffer.data(gl + 391);
    const auto *gl_392 = buffer.data(gl + 392);
    const auto *gl_393 = buffer.data(gl + 393);
    const auto *gl_394 = buffer.data(gl + 394);
    const auto *gl_395 = buffer.data(gl + 395);
    const auto *gl_396 = buffer.data(gl + 396);
    const auto *gl_397 = buffer.data(gl + 397);
    const auto *gl_398 = buffer.data(gl + 398);
    const auto *gl_399 = buffer.data(gl + 399);
    const auto *gl_400 = buffer.data(gl + 400);
    const auto *gl_401 = buffer.data(gl + 401);
    const auto *gl_402 = buffer.data(gl + 402);
    const auto *gl_403 = buffer.data(gl + 403);
    const auto *gl_404 = buffer.data(gl + 404);
    const auto *gl_450 = buffer.data(gl + 450);
    const auto *gl_451 = buffer.data(gl + 451);
    const auto *gl_452 = buffer.data(gl + 452);
    const auto *gl_453 = buffer.data(gl + 453);
    const auto *gl_454 = buffer.data(gl + 454);
    const auto *gl_455 = buffer.data(gl + 455);
    const auto *gl_456 = buffer.data(gl + 456);
    const auto *gl_457 = buffer.data(gl + 457);
    const auto *gl_458 = buffer.data(gl + 458);
    const auto *gl_459 = buffer.data(gl + 459);
    const auto *gl_460 = buffer.data(gl + 460);
    const auto *gl_461 = buffer.data(gl + 461);
    const auto *gl_462 = buffer.data(gl + 462);
    const auto *gl_463 = buffer.data(gl + 463);
    const auto *gl_464 = buffer.data(gl + 464);
    const auto *gl_465 = buffer.data(gl + 465);
    const auto *gl_466 = buffer.data(gl + 466);
    const auto *gl_467 = buffer.data(gl + 467);
    const auto *gl_468 = buffer.data(gl + 468);
    const auto *gl_469 = buffer.data(gl + 469);
    const auto *gl_470 = buffer.data(gl + 470);
    const auto *gl_471 = buffer.data(gl + 471);
    const auto *gl_472 = buffer.data(gl + 472);
    const auto *gl_473 = buffer.data(gl + 473);
    const auto *gl_474 = buffer.data(gl + 474);
    const auto *gl_475 = buffer.data(gl + 475);
    const auto *gl_476 = buffer.data(gl + 476);
    const auto *gl_477 = buffer.data(gl + 477);
    const auto *gl_478 = buffer.data(gl + 478);
    const auto *gl_479 = buffer.data(gl + 479);
    const auto *gl_480 = buffer.data(gl + 480);
    const auto *gl_481 = buffer.data(gl + 481);
    const auto *gl_482 = buffer.data(gl + 482);
    const auto *gl_483 = buffer.data(gl + 483);
    const auto *gl_484 = buffer.data(gl + 484);
    const auto *gl_485 = buffer.data(gl + 485);
    const auto *gl_486 = buffer.data(gl + 486);
    const auto *gl_487 = buffer.data(gl + 487);
    const auto *gl_488 = buffer.data(gl + 488);
    const auto *gl_489 = buffer.data(gl + 489);
    const auto *gl_490 = buffer.data(gl + 490);
    const auto *gl_491 = buffer.data(gl + 491);
    const auto *gl_492 = buffer.data(gl + 492);
    const auto *gl_493 = buffer.data(gl + 493);
    const auto *gl_494 = buffer.data(gl + 494);
    const auto *gl_495 = buffer.data(gl + 495);
    const auto *gl_496 = buffer.data(gl + 496);
    const auto *gl_497 = buffer.data(gl + 497);
    const auto *gl_498 = buffer.data(gl + 498);
    const auto *gl_499 = buffer.data(gl + 499);
    const auto *gl_500 = buffer.data(gl + 500);
    const auto *gl_501 = buffer.data(gl + 501);
    const auto *gl_502 = buffer.data(gl + 502);
    const auto *gl_503 = buffer.data(gl + 503);
    const auto *gl_504 = buffer.data(gl + 504);
    const auto *gl_505 = buffer.data(gl + 505);
    const auto *gl_506 = buffer.data(gl + 506);
    const auto *gl_507 = buffer.data(gl + 507);
    const auto *gl_508 = buffer.data(gl + 508);
    const auto *gl_509 = buffer.data(gl + 509);
    const auto *gl_510 = buffer.data(gl + 510);
    const auto *gl_511 = buffer.data(gl + 511);
    const auto *gl_512 = buffer.data(gl + 512);
    const auto *gl_513 = buffer.data(gl + 513);
    const auto *gl_514 = buffer.data(gl + 514);
    const auto *gl_515 = buffer.data(gl + 515);
    const auto *gl_516 = buffer.data(gl + 516);
    const auto *gl_517 = buffer.data(gl + 517);
    const auto *gl_518 = buffer.data(gl + 518);
    const auto *gl_519 = buffer.data(gl + 519);
    const auto *gl_520 = buffer.data(gl + 520);
    const auto *gl_521 = buffer.data(gl + 521);
    const auto *gl_522 = buffer.data(gl + 522);
    const auto *gl_523 = buffer.data(gl + 523);
    const auto *gl_524 = buffer.data(gl + 524);

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, dl_89, dl_90, dl_91, dl_92, dl_93, \
                         gl_314, gl_315, gl_316, gl_317, gl_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -2.0 * dl_89[k]
                   + f_0 * gl_314[k];

        t_180[k] = -dl_90[k]
                   + f_0 * gl_315[k];

        t_181[k] = -dl_91[k]
                   + f_0 * gl_316[k];

        t_182[k] = -dl_92[k]
                   + f_0 * gl_317[k];

        t_183[k] = -dl_93[k]
                   + f_0 * gl_318[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, dl_94, dl_95, dl_96, dl_97, dl_98, \
                         gl_319, gl_320, gl_321, gl_322, gl_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -dl_94[k]
                   + f_0 * gl_319[k];

        t_185[k] = -dl_95[k]
                   + f_0 * gl_320[k];

        t_186[k] = -dl_96[k]
                   + f_0 * gl_321[k];

        t_187[k] = -dl_97[k]
                   + f_0 * gl_322[k];

        t_188[k] = -dl_98[k]
                   + f_0 * gl_323[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, dl_99, dl_100, dl_101, dl_102, \
                         dl_103, gl_324, gl_325, gl_326, gl_327, \
                         gl_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = -dl_99[k]
                   + f_0 * gl_324[k];

        t_190[k] = -dl_100[k]
                   + f_0 * gl_325[k];

        t_191[k] = -dl_101[k]
                   + f_0 * gl_326[k];

        t_192[k] = -dl_102[k]
                   + f_0 * gl_327[k];

        t_193[k] = -dl_103[k]
                   + f_0 * gl_328[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, dl_104, dl_105, dl_106, dl_107, \
                         dl_108, gl_329, gl_330, gl_331, gl_332, \
                         gl_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -dl_104[k]
                   + f_0 * gl_329[k];

        t_195[k] = -dl_105[k]
                   + f_0 * gl_330[k];

        t_196[k] = -dl_106[k]
                   + f_0 * gl_331[k];

        t_197[k] = -dl_107[k]
                   + f_0 * gl_332[k];

        t_198[k] = -dl_108[k]
                   + f_0 * gl_333[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, dl_109, dl_110, dl_111, dl_112, \
                         dl_113, gl_334, gl_335, gl_336, gl_337, \
                         gl_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = -dl_109[k]
                   + f_0 * gl_334[k];

        t_200[k] = -dl_110[k]
                   + f_0 * gl_335[k];

        t_201[k] = -dl_111[k]
                   + f_0 * gl_336[k];

        t_202[k] = -dl_112[k]
                   + f_0 * gl_337[k];

        t_203[k] = -dl_113[k]
                   + f_0 * gl_338[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, dl_114, dl_115, dl_116, dl_117, \
                         dl_118, gl_339, gl_340, gl_341, gl_342, \
                         gl_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -dl_114[k]
                   + f_0 * gl_339[k];

        t_205[k] = -dl_115[k]
                   + f_0 * gl_340[k];

        t_206[k] = -dl_116[k]
                   + f_0 * gl_341[k];

        t_207[k] = -dl_117[k]
                   + f_0 * gl_342[k];

        t_208[k] = -dl_118[k]
                   + f_0 * gl_343[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, dl_119, dl_120, dl_121, dl_122, \
                         dl_123, gl_344, gl_345, gl_346, gl_347, \
                         gl_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = -dl_119[k]
                   + f_0 * gl_344[k];

        t_210[k] = -dl_120[k]
                   + f_0 * gl_345[k];

        t_211[k] = -dl_121[k]
                   + f_0 * gl_346[k];

        t_212[k] = -dl_122[k]
                   + f_0 * gl_347[k];

        t_213[k] = -dl_123[k]
                   + f_0 * gl_348[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, dl_124, dl_125, dl_126, dl_127, \
                         dl_128, gl_349, gl_350, gl_351, gl_352, \
                         gl_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = -dl_124[k]
                   + f_0 * gl_349[k];

        t_215[k] = -dl_125[k]
                   + f_0 * gl_350[k];

        t_216[k] = -dl_126[k]
                   + f_0 * gl_351[k];

        t_217[k] = -dl_127[k]
                   + f_0 * gl_352[k];

        t_218[k] = -dl_128[k]
                   + f_0 * gl_353[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, dl_129, dl_130, dl_131, dl_132, \
                         dl_133, gl_354, gl_355, gl_356, gl_357, \
                         gl_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = -dl_129[k]
                   + f_0 * gl_354[k];

        t_220[k] = -dl_130[k]
                   + f_0 * gl_355[k];

        t_221[k] = -dl_131[k]
                   + f_0 * gl_356[k];

        t_222[k] = -dl_132[k]
                   + f_0 * gl_357[k];

        t_223[k] = -dl_133[k]
                   + f_0 * gl_358[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, t_229, t_230, dl_134, gl_359, \
                         gl_360, gl_361, gl_362, gl_363, gl_364, \
                         gl_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = -dl_134[k]
                   + f_0 * gl_359[k];

        t_225[k] = f_0 * gl_360[k];

        t_226[k] = f_0 * gl_361[k];

        t_227[k] = f_0 * gl_362[k];

        t_228[k] = f_0 * gl_363[k];

        t_229[k] = f_0 * gl_364[k];

        t_230[k] = f_0 * gl_365[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, t_236, t_237, t_238, gl_366, \
                         gl_367, gl_368, gl_369, gl_370, gl_371, gl_372, \
                         gl_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_0 * gl_366[k];

        t_232[k] = f_0 * gl_367[k];

        t_233[k] = f_0 * gl_368[k];

        t_234[k] = f_0 * gl_369[k];

        t_235[k] = f_0 * gl_370[k];

        t_236[k] = f_0 * gl_371[k];

        t_237[k] = f_0 * gl_372[k];

        t_238[k] = f_0 * gl_373[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, t_244, t_245, t_246, gl_374, \
                         gl_375, gl_376, gl_377, gl_378, gl_379, gl_380, \
                         gl_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_0 * gl_374[k];

        t_240[k] = f_0 * gl_375[k];

        t_241[k] = f_0 * gl_376[k];

        t_242[k] = f_0 * gl_377[k];

        t_243[k] = f_0 * gl_378[k];

        t_244[k] = f_0 * gl_379[k];

        t_245[k] = f_0 * gl_380[k];

        t_246[k] = f_0 * gl_381[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, t_252, t_253, t_254, gl_382, \
                         gl_383, gl_384, gl_385, gl_386, gl_387, gl_388, \
                         gl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_0 * gl_382[k];

        t_248[k] = f_0 * gl_383[k];

        t_249[k] = f_0 * gl_384[k];

        t_250[k] = f_0 * gl_385[k];

        t_251[k] = f_0 * gl_386[k];

        t_252[k] = f_0 * gl_387[k];

        t_253[k] = f_0 * gl_388[k];

        t_254[k] = f_0 * gl_389[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, t_260, t_261, t_262, gl_390, \
                         gl_391, gl_392, gl_393, gl_394, gl_395, gl_396, \
                         gl_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_0 * gl_390[k];

        t_256[k] = f_0 * gl_391[k];

        t_257[k] = f_0 * gl_392[k];

        t_258[k] = f_0 * gl_393[k];

        t_259[k] = f_0 * gl_394[k];

        t_260[k] = f_0 * gl_395[k];

        t_261[k] = f_0 * gl_396[k];

        t_262[k] = f_0 * gl_397[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, t_268, t_269, gl_398, gl_399, \
                         gl_400, gl_401, gl_402, gl_403, gl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_0 * gl_398[k];

        t_264[k] = f_0 * gl_399[k];

        t_265[k] = f_0 * gl_400[k];

        t_266[k] = f_0 * gl_401[k];

        t_267[k] = f_0 * gl_402[k];

        t_268[k] = f_0 * gl_403[k];

        t_269[k] = f_0 * gl_404[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, dl_135, dl_136, dl_137, dl_138, \
                         dl_139, gl_450, gl_451, gl_452, gl_453, \
                         gl_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -3.0 * dl_135[k]
                   + f_0 * gl_450[k];

        t_271[k] = -3.0 * dl_136[k]
                   + f_0 * gl_451[k];

        t_272[k] = -3.0 * dl_137[k]
                   + f_0 * gl_452[k];

        t_273[k] = -3.0 * dl_138[k]
                   + f_0 * gl_453[k];

        t_274[k] = -3.0 * dl_139[k]
                   + f_0 * gl_454[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, dl_140, dl_141, dl_142, dl_143, \
                         dl_144, gl_455, gl_456, gl_457, gl_458, \
                         gl_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -3.0 * dl_140[k]
                   + f_0 * gl_455[k];

        t_276[k] = -3.0 * dl_141[k]
                   + f_0 * gl_456[k];

        t_277[k] = -3.0 * dl_142[k]
                   + f_0 * gl_457[k];

        t_278[k] = -3.0 * dl_143[k]
                   + f_0 * gl_458[k];

        t_279[k] = -3.0 * dl_144[k]
                   + f_0 * gl_459[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, dl_145, dl_146, dl_147, dl_148, \
                         dl_149, gl_460, gl_461, gl_462, gl_463, \
                         gl_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -3.0 * dl_145[k]
                   + f_0 * gl_460[k];

        t_281[k] = -3.0 * dl_146[k]
                   + f_0 * gl_461[k];

        t_282[k] = -3.0 * dl_147[k]
                   + f_0 * gl_462[k];

        t_283[k] = -3.0 * dl_148[k]
                   + f_0 * gl_463[k];

        t_284[k] = -3.0 * dl_149[k]
                   + f_0 * gl_464[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, dl_150, dl_151, dl_152, dl_153, \
                         dl_154, gl_465, gl_466, gl_467, gl_468, \
                         gl_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -3.0 * dl_150[k]
                   + f_0 * gl_465[k];

        t_286[k] = -3.0 * dl_151[k]
                   + f_0 * gl_466[k];

        t_287[k] = -3.0 * dl_152[k]
                   + f_0 * gl_467[k];

        t_288[k] = -3.0 * dl_153[k]
                   + f_0 * gl_468[k];

        t_289[k] = -3.0 * dl_154[k]
                   + f_0 * gl_469[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, dl_155, dl_156, dl_157, dl_158, \
                         dl_159, gl_470, gl_471, gl_472, gl_473, \
                         gl_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -3.0 * dl_155[k]
                   + f_0 * gl_470[k];

        t_291[k] = -3.0 * dl_156[k]
                   + f_0 * gl_471[k];

        t_292[k] = -3.0 * dl_157[k]
                   + f_0 * gl_472[k];

        t_293[k] = -3.0 * dl_158[k]
                   + f_0 * gl_473[k];

        t_294[k] = -3.0 * dl_159[k]
                   + f_0 * gl_474[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, dl_160, dl_161, dl_162, dl_163, \
                         dl_164, gl_475, gl_476, gl_477, gl_478, \
                         gl_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -3.0 * dl_160[k]
                   + f_0 * gl_475[k];

        t_296[k] = -3.0 * dl_161[k]
                   + f_0 * gl_476[k];

        t_297[k] = -3.0 * dl_162[k]
                   + f_0 * gl_477[k];

        t_298[k] = -3.0 * dl_163[k]
                   + f_0 * gl_478[k];

        t_299[k] = -3.0 * dl_164[k]
                   + f_0 * gl_479[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, dl_165, dl_166, dl_167, dl_168, \
                         dl_169, gl_480, gl_481, gl_482, gl_483, \
                         gl_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -3.0 * dl_165[k]
                   + f_0 * gl_480[k];

        t_301[k] = -3.0 * dl_166[k]
                   + f_0 * gl_481[k];

        t_302[k] = -3.0 * dl_167[k]
                   + f_0 * gl_482[k];

        t_303[k] = -3.0 * dl_168[k]
                   + f_0 * gl_483[k];

        t_304[k] = -3.0 * dl_169[k]
                   + f_0 * gl_484[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, dl_170, dl_171, dl_172, dl_173, \
                         dl_174, gl_485, gl_486, gl_487, gl_488, \
                         gl_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -3.0 * dl_170[k]
                   + f_0 * gl_485[k];

        t_306[k] = -3.0 * dl_171[k]
                   + f_0 * gl_486[k];

        t_307[k] = -3.0 * dl_172[k]
                   + f_0 * gl_487[k];

        t_308[k] = -3.0 * dl_173[k]
                   + f_0 * gl_488[k];

        t_309[k] = -3.0 * dl_174[k]
                   + f_0 * gl_489[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, dl_175, dl_176, dl_177, dl_178, \
                         dl_179, gl_490, gl_491, gl_492, gl_493, \
                         gl_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * dl_175[k]
                   + f_0 * gl_490[k];

        t_311[k] = -3.0 * dl_176[k]
                   + f_0 * gl_491[k];

        t_312[k] = -3.0 * dl_177[k]
                   + f_0 * gl_492[k];

        t_313[k] = -3.0 * dl_178[k]
                   + f_0 * gl_493[k];

        t_314[k] = -3.0 * dl_179[k]
                   + f_0 * gl_494[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, dl_180, dl_181, dl_182, dl_183, \
                         dl_184, gl_495, gl_496, gl_497, gl_498, \
                         gl_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -2.0 * dl_180[k]
                   + f_0 * gl_495[k];

        t_316[k] = -2.0 * dl_181[k]
                   + f_0 * gl_496[k];

        t_317[k] = -2.0 * dl_182[k]
                   + f_0 * gl_497[k];

        t_318[k] = -2.0 * dl_183[k]
                   + f_0 * gl_498[k];

        t_319[k] = -2.0 * dl_184[k]
                   + f_0 * gl_499[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, dl_185, dl_186, dl_187, dl_188, \
                         dl_189, gl_500, gl_501, gl_502, gl_503, \
                         gl_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -2.0 * dl_185[k]
                   + f_0 * gl_500[k];

        t_321[k] = -2.0 * dl_186[k]
                   + f_0 * gl_501[k];

        t_322[k] = -2.0 * dl_187[k]
                   + f_0 * gl_502[k];

        t_323[k] = -2.0 * dl_188[k]
                   + f_0 * gl_503[k];

        t_324[k] = -2.0 * dl_189[k]
                   + f_0 * gl_504[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, dl_190, dl_191, dl_192, dl_193, \
                         dl_194, gl_505, gl_506, gl_507, gl_508, \
                         gl_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -2.0 * dl_190[k]
                   + f_0 * gl_505[k];

        t_326[k] = -2.0 * dl_191[k]
                   + f_0 * gl_506[k];

        t_327[k] = -2.0 * dl_192[k]
                   + f_0 * gl_507[k];

        t_328[k] = -2.0 * dl_193[k]
                   + f_0 * gl_508[k];

        t_329[k] = -2.0 * dl_194[k]
                   + f_0 * gl_509[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, dl_195, dl_196, dl_197, dl_198, \
                         dl_199, gl_510, gl_511, gl_512, gl_513, \
                         gl_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -2.0 * dl_195[k]
                   + f_0 * gl_510[k];

        t_331[k] = -2.0 * dl_196[k]
                   + f_0 * gl_511[k];

        t_332[k] = -2.0 * dl_197[k]
                   + f_0 * gl_512[k];

        t_333[k] = -2.0 * dl_198[k]
                   + f_0 * gl_513[k];

        t_334[k] = -2.0 * dl_199[k]
                   + f_0 * gl_514[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, dl_200, dl_201, dl_202, dl_203, \
                         dl_204, gl_515, gl_516, gl_517, gl_518, \
                         gl_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -2.0 * dl_200[k]
                   + f_0 * gl_515[k];

        t_336[k] = -2.0 * dl_201[k]
                   + f_0 * gl_516[k];

        t_337[k] = -2.0 * dl_202[k]
                   + f_0 * gl_517[k];

        t_338[k] = -2.0 * dl_203[k]
                   + f_0 * gl_518[k];

        t_339[k] = -2.0 * dl_204[k]
                   + f_0 * gl_519[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, dl_205, dl_206, dl_207, dl_208, \
                         dl_209, gl_520, gl_521, gl_522, gl_523, \
                         gl_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * dl_205[k]
                   + f_0 * gl_520[k];

        t_341[k] = -2.0 * dl_206[k]
                   + f_0 * gl_521[k];

        t_342[k] = -2.0 * dl_207[k]
                   + f_0 * gl_522[k];

        t_343[k] = -2.0 * dl_208[k]
                   + f_0 * gl_523[k];

        t_344[k] = -2.0 * dl_209[k]
                   + f_0 * gl_524[k];
    }
}

static auto
compute_prim_geom_10_fl_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dl, const size_t gl,
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

    const auto *dl_210 = buffer.data(dl + 210);
    const auto *dl_211 = buffer.data(dl + 211);
    const auto *dl_212 = buffer.data(dl + 212);
    const auto *dl_213 = buffer.data(dl + 213);
    const auto *dl_214 = buffer.data(dl + 214);
    const auto *dl_215 = buffer.data(dl + 215);
    const auto *dl_216 = buffer.data(dl + 216);
    const auto *dl_217 = buffer.data(dl + 217);
    const auto *dl_218 = buffer.data(dl + 218);
    const auto *dl_219 = buffer.data(dl + 219);
    const auto *dl_220 = buffer.data(dl + 220);
    const auto *dl_221 = buffer.data(dl + 221);
    const auto *dl_222 = buffer.data(dl + 222);
    const auto *dl_223 = buffer.data(dl + 223);
    const auto *dl_224 = buffer.data(dl + 224);
    const auto *dl_225 = buffer.data(dl + 225);
    const auto *dl_226 = buffer.data(dl + 226);
    const auto *dl_227 = buffer.data(dl + 227);
    const auto *dl_228 = buffer.data(dl + 228);
    const auto *dl_229 = buffer.data(dl + 229);
    const auto *dl_230 = buffer.data(dl + 230);
    const auto *dl_231 = buffer.data(dl + 231);
    const auto *dl_232 = buffer.data(dl + 232);
    const auto *dl_233 = buffer.data(dl + 233);
    const auto *dl_234 = buffer.data(dl + 234);
    const auto *dl_235 = buffer.data(dl + 235);
    const auto *dl_236 = buffer.data(dl + 236);
    const auto *dl_237 = buffer.data(dl + 237);
    const auto *dl_238 = buffer.data(dl + 238);
    const auto *dl_239 = buffer.data(dl + 239);
    const auto *dl_240 = buffer.data(dl + 240);
    const auto *dl_241 = buffer.data(dl + 241);
    const auto *dl_242 = buffer.data(dl + 242);
    const auto *dl_243 = buffer.data(dl + 243);
    const auto *dl_244 = buffer.data(dl + 244);
    const auto *dl_245 = buffer.data(dl + 245);
    const auto *dl_246 = buffer.data(dl + 246);
    const auto *dl_247 = buffer.data(dl + 247);
    const auto *dl_248 = buffer.data(dl + 248);
    const auto *dl_249 = buffer.data(dl + 249);
    const auto *dl_250 = buffer.data(dl + 250);
    const auto *dl_251 = buffer.data(dl + 251);
    const auto *dl_252 = buffer.data(dl + 252);
    const auto *dl_253 = buffer.data(dl + 253);
    const auto *dl_254 = buffer.data(dl + 254);
    const auto *dl_255 = buffer.data(dl + 255);
    const auto *dl_256 = buffer.data(dl + 256);
    const auto *dl_257 = buffer.data(dl + 257);
    const auto *dl_258 = buffer.data(dl + 258);
    const auto *dl_259 = buffer.data(dl + 259);
    const auto *dl_260 = buffer.data(dl + 260);
    const auto *dl_261 = buffer.data(dl + 261);
    const auto *dl_262 = buffer.data(dl + 262);
    const auto *dl_263 = buffer.data(dl + 263);
    const auto *dl_264 = buffer.data(dl + 264);
    const auto *dl_265 = buffer.data(dl + 265);
    const auto *dl_266 = buffer.data(dl + 266);
    const auto *dl_267 = buffer.data(dl + 267);
    const auto *dl_268 = buffer.data(dl + 268);
    const auto *dl_269 = buffer.data(dl + 269);

    const auto *gl_525 = buffer.data(gl + 525);
    const auto *gl_526 = buffer.data(gl + 526);
    const auto *gl_527 = buffer.data(gl + 527);
    const auto *gl_528 = buffer.data(gl + 528);
    const auto *gl_529 = buffer.data(gl + 529);
    const auto *gl_530 = buffer.data(gl + 530);
    const auto *gl_531 = buffer.data(gl + 531);
    const auto *gl_532 = buffer.data(gl + 532);
    const auto *gl_533 = buffer.data(gl + 533);
    const auto *gl_534 = buffer.data(gl + 534);
    const auto *gl_535 = buffer.data(gl + 535);
    const auto *gl_536 = buffer.data(gl + 536);
    const auto *gl_537 = buffer.data(gl + 537);
    const auto *gl_538 = buffer.data(gl + 538);
    const auto *gl_539 = buffer.data(gl + 539);
    const auto *gl_540 = buffer.data(gl + 540);
    const auto *gl_541 = buffer.data(gl + 541);
    const auto *gl_542 = buffer.data(gl + 542);
    const auto *gl_543 = buffer.data(gl + 543);
    const auto *gl_544 = buffer.data(gl + 544);
    const auto *gl_545 = buffer.data(gl + 545);
    const auto *gl_546 = buffer.data(gl + 546);
    const auto *gl_547 = buffer.data(gl + 547);
    const auto *gl_548 = buffer.data(gl + 548);
    const auto *gl_549 = buffer.data(gl + 549);
    const auto *gl_550 = buffer.data(gl + 550);
    const auto *gl_551 = buffer.data(gl + 551);
    const auto *gl_552 = buffer.data(gl + 552);
    const auto *gl_553 = buffer.data(gl + 553);
    const auto *gl_554 = buffer.data(gl + 554);
    const auto *gl_555 = buffer.data(gl + 555);
    const auto *gl_556 = buffer.data(gl + 556);
    const auto *gl_557 = buffer.data(gl + 557);
    const auto *gl_558 = buffer.data(gl + 558);
    const auto *gl_559 = buffer.data(gl + 559);
    const auto *gl_560 = buffer.data(gl + 560);
    const auto *gl_561 = buffer.data(gl + 561);
    const auto *gl_562 = buffer.data(gl + 562);
    const auto *gl_563 = buffer.data(gl + 563);
    const auto *gl_564 = buffer.data(gl + 564);
    const auto *gl_565 = buffer.data(gl + 565);
    const auto *gl_566 = buffer.data(gl + 566);
    const auto *gl_567 = buffer.data(gl + 567);
    const auto *gl_568 = buffer.data(gl + 568);
    const auto *gl_569 = buffer.data(gl + 569);
    const auto *gl_570 = buffer.data(gl + 570);
    const auto *gl_571 = buffer.data(gl + 571);
    const auto *gl_572 = buffer.data(gl + 572);
    const auto *gl_573 = buffer.data(gl + 573);
    const auto *gl_574 = buffer.data(gl + 574);
    const auto *gl_575 = buffer.data(gl + 575);
    const auto *gl_576 = buffer.data(gl + 576);
    const auto *gl_577 = buffer.data(gl + 577);
    const auto *gl_578 = buffer.data(gl + 578);
    const auto *gl_579 = buffer.data(gl + 579);
    const auto *gl_580 = buffer.data(gl + 580);
    const auto *gl_581 = buffer.data(gl + 581);
    const auto *gl_582 = buffer.data(gl + 582);
    const auto *gl_583 = buffer.data(gl + 583);
    const auto *gl_584 = buffer.data(gl + 584);
    const auto *gl_585 = buffer.data(gl + 585);
    const auto *gl_586 = buffer.data(gl + 586);
    const auto *gl_587 = buffer.data(gl + 587);
    const auto *gl_588 = buffer.data(gl + 588);
    const auto *gl_589 = buffer.data(gl + 589);
    const auto *gl_590 = buffer.data(gl + 590);
    const auto *gl_591 = buffer.data(gl + 591);
    const auto *gl_592 = buffer.data(gl + 592);
    const auto *gl_593 = buffer.data(gl + 593);
    const auto *gl_594 = buffer.data(gl + 594);
    const auto *gl_595 = buffer.data(gl + 595);
    const auto *gl_596 = buffer.data(gl + 596);
    const auto *gl_597 = buffer.data(gl + 597);
    const auto *gl_598 = buffer.data(gl + 598);
    const auto *gl_599 = buffer.data(gl + 599);
    const auto *gl_600 = buffer.data(gl + 600);
    const auto *gl_601 = buffer.data(gl + 601);
    const auto *gl_602 = buffer.data(gl + 602);
    const auto *gl_603 = buffer.data(gl + 603);
    const auto *gl_604 = buffer.data(gl + 604);
    const auto *gl_605 = buffer.data(gl + 605);
    const auto *gl_606 = buffer.data(gl + 606);
    const auto *gl_607 = buffer.data(gl + 607);
    const auto *gl_608 = buffer.data(gl + 608);
    const auto *gl_609 = buffer.data(gl + 609);
    const auto *gl_610 = buffer.data(gl + 610);
    const auto *gl_611 = buffer.data(gl + 611);
    const auto *gl_612 = buffer.data(gl + 612);
    const auto *gl_613 = buffer.data(gl + 613);
    const auto *gl_614 = buffer.data(gl + 614);
    const auto *gl_615 = buffer.data(gl + 615);
    const auto *gl_616 = buffer.data(gl + 616);
    const auto *gl_617 = buffer.data(gl + 617);
    const auto *gl_618 = buffer.data(gl + 618);
    const auto *gl_619 = buffer.data(gl + 619);
    const auto *gl_620 = buffer.data(gl + 620);
    const auto *gl_621 = buffer.data(gl + 621);
    const auto *gl_622 = buffer.data(gl + 622);
    const auto *gl_623 = buffer.data(gl + 623);
    const auto *gl_624 = buffer.data(gl + 624);
    const auto *gl_625 = buffer.data(gl + 625);
    const auto *gl_626 = buffer.data(gl + 626);
    const auto *gl_627 = buffer.data(gl + 627);
    const auto *gl_628 = buffer.data(gl + 628);
    const auto *gl_629 = buffer.data(gl + 629);

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, dl_210, dl_211, dl_212, dl_213, \
                         dl_214, gl_525, gl_526, gl_527, gl_528, \
                         gl_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * dl_210[k]
                   + f_0 * gl_525[k];

        t_346[k] = -2.0 * dl_211[k]
                   + f_0 * gl_526[k];

        t_347[k] = -2.0 * dl_212[k]
                   + f_0 * gl_527[k];

        t_348[k] = -2.0 * dl_213[k]
                   + f_0 * gl_528[k];

        t_349[k] = -2.0 * dl_214[k]
                   + f_0 * gl_529[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, dl_215, dl_216, dl_217, dl_218, \
                         dl_219, gl_530, gl_531, gl_532, gl_533, \
                         gl_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * dl_215[k]
                   + f_0 * gl_530[k];

        t_351[k] = -2.0 * dl_216[k]
                   + f_0 * gl_531[k];

        t_352[k] = -2.0 * dl_217[k]
                   + f_0 * gl_532[k];

        t_353[k] = -2.0 * dl_218[k]
                   + f_0 * gl_533[k];

        t_354[k] = -2.0 * dl_219[k]
                   + f_0 * gl_534[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, dl_220, dl_221, dl_222, dl_223, \
                         dl_224, gl_535, gl_536, gl_537, gl_538, \
                         gl_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * dl_220[k]
                   + f_0 * gl_535[k];

        t_356[k] = -2.0 * dl_221[k]
                   + f_0 * gl_536[k];

        t_357[k] = -2.0 * dl_222[k]
                   + f_0 * gl_537[k];

        t_358[k] = -2.0 * dl_223[k]
                   + f_0 * gl_538[k];

        t_359[k] = -2.0 * dl_224[k]
                   + f_0 * gl_539[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, dl_225, dl_226, dl_227, dl_228, \
                         dl_229, gl_540, gl_541, gl_542, gl_543, \
                         gl_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -dl_225[k]
                   + f_0 * gl_540[k];

        t_361[k] = -dl_226[k]
                   + f_0 * gl_541[k];

        t_362[k] = -dl_227[k]
                   + f_0 * gl_542[k];

        t_363[k] = -dl_228[k]
                   + f_0 * gl_543[k];

        t_364[k] = -dl_229[k]
                   + f_0 * gl_544[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, dl_230, dl_231, dl_232, dl_233, \
                         dl_234, gl_545, gl_546, gl_547, gl_548, \
                         gl_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -dl_230[k]
                   + f_0 * gl_545[k];

        t_366[k] = -dl_231[k]
                   + f_0 * gl_546[k];

        t_367[k] = -dl_232[k]
                   + f_0 * gl_547[k];

        t_368[k] = -dl_233[k]
                   + f_0 * gl_548[k];

        t_369[k] = -dl_234[k]
                   + f_0 * gl_549[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, dl_235, dl_236, dl_237, dl_238, \
                         dl_239, gl_550, gl_551, gl_552, gl_553, \
                         gl_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -dl_235[k]
                   + f_0 * gl_550[k];

        t_371[k] = -dl_236[k]
                   + f_0 * gl_551[k];

        t_372[k] = -dl_237[k]
                   + f_0 * gl_552[k];

        t_373[k] = -dl_238[k]
                   + f_0 * gl_553[k];

        t_374[k] = -dl_239[k]
                   + f_0 * gl_554[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, dl_240, dl_241, dl_242, dl_243, \
                         dl_244, gl_555, gl_556, gl_557, gl_558, \
                         gl_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -dl_240[k]
                   + f_0 * gl_555[k];

        t_376[k] = -dl_241[k]
                   + f_0 * gl_556[k];

        t_377[k] = -dl_242[k]
                   + f_0 * gl_557[k];

        t_378[k] = -dl_243[k]
                   + f_0 * gl_558[k];

        t_379[k] = -dl_244[k]
                   + f_0 * gl_559[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, dl_245, dl_246, dl_247, dl_248, \
                         dl_249, gl_560, gl_561, gl_562, gl_563, \
                         gl_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -dl_245[k]
                   + f_0 * gl_560[k];

        t_381[k] = -dl_246[k]
                   + f_0 * gl_561[k];

        t_382[k] = -dl_247[k]
                   + f_0 * gl_562[k];

        t_383[k] = -dl_248[k]
                   + f_0 * gl_563[k];

        t_384[k] = -dl_249[k]
                   + f_0 * gl_564[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, dl_250, dl_251, dl_252, dl_253, \
                         dl_254, gl_565, gl_566, gl_567, gl_568, \
                         gl_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -dl_250[k]
                   + f_0 * gl_565[k];

        t_386[k] = -dl_251[k]
                   + f_0 * gl_566[k];

        t_387[k] = -dl_252[k]
                   + f_0 * gl_567[k];

        t_388[k] = -dl_253[k]
                   + f_0 * gl_568[k];

        t_389[k] = -dl_254[k]
                   + f_0 * gl_569[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, dl_255, dl_256, dl_257, dl_258, \
                         dl_259, gl_570, gl_571, gl_572, gl_573, \
                         gl_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -dl_255[k]
                   + f_0 * gl_570[k];

        t_391[k] = -dl_256[k]
                   + f_0 * gl_571[k];

        t_392[k] = -dl_257[k]
                   + f_0 * gl_572[k];

        t_393[k] = -dl_258[k]
                   + f_0 * gl_573[k];

        t_394[k] = -dl_259[k]
                   + f_0 * gl_574[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, dl_260, dl_261, dl_262, dl_263, \
                         dl_264, gl_575, gl_576, gl_577, gl_578, \
                         gl_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -dl_260[k]
                   + f_0 * gl_575[k];

        t_396[k] = -dl_261[k]
                   + f_0 * gl_576[k];

        t_397[k] = -dl_262[k]
                   + f_0 * gl_577[k];

        t_398[k] = -dl_263[k]
                   + f_0 * gl_578[k];

        t_399[k] = -dl_264[k]
                   + f_0 * gl_579[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, dl_265, dl_266, dl_267, dl_268, \
                         dl_269, gl_580, gl_581, gl_582, gl_583, \
                         gl_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -dl_265[k]
                   + f_0 * gl_580[k];

        t_401[k] = -dl_266[k]
                   + f_0 * gl_581[k];

        t_402[k] = -dl_267[k]
                   + f_0 * gl_582[k];

        t_403[k] = -dl_268[k]
                   + f_0 * gl_583[k];

        t_404[k] = -dl_269[k]
                   + f_0 * gl_584[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, gl_585, \
                         gl_586, gl_587, gl_588, gl_589, gl_590, gl_591, \
                         gl_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * gl_585[k];

        t_406[k] = f_0 * gl_586[k];

        t_407[k] = f_0 * gl_587[k];

        t_408[k] = f_0 * gl_588[k];

        t_409[k] = f_0 * gl_589[k];

        t_410[k] = f_0 * gl_590[k];

        t_411[k] = f_0 * gl_591[k];

        t_412[k] = f_0 * gl_592[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, t_420, gl_593, \
                         gl_594, gl_595, gl_596, gl_597, gl_598, gl_599, \
                         gl_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * gl_593[k];

        t_414[k] = f_0 * gl_594[k];

        t_415[k] = f_0 * gl_595[k];

        t_416[k] = f_0 * gl_596[k];

        t_417[k] = f_0 * gl_597[k];

        t_418[k] = f_0 * gl_598[k];

        t_419[k] = f_0 * gl_599[k];

        t_420[k] = f_0 * gl_600[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, t_425, t_426, t_427, t_428, gl_601, \
                         gl_602, gl_603, gl_604, gl_605, gl_606, gl_607, \
                         gl_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_0 * gl_601[k];

        t_422[k] = f_0 * gl_602[k];

        t_423[k] = f_0 * gl_603[k];

        t_424[k] = f_0 * gl_604[k];

        t_425[k] = f_0 * gl_605[k];

        t_426[k] = f_0 * gl_606[k];

        t_427[k] = f_0 * gl_607[k];

        t_428[k] = f_0 * gl_608[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, t_434, t_435, t_436, gl_609, \
                         gl_610, gl_611, gl_612, gl_613, gl_614, gl_615, \
                         gl_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_0 * gl_609[k];

        t_430[k] = f_0 * gl_610[k];

        t_431[k] = f_0 * gl_611[k];

        t_432[k] = f_0 * gl_612[k];

        t_433[k] = f_0 * gl_613[k];

        t_434[k] = f_0 * gl_614[k];

        t_435[k] = f_0 * gl_615[k];

        t_436[k] = f_0 * gl_616[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, t_442, t_443, t_444, gl_617, \
                         gl_618, gl_619, gl_620, gl_621, gl_622, gl_623, \
                         gl_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_0 * gl_617[k];

        t_438[k] = f_0 * gl_618[k];

        t_439[k] = f_0 * gl_619[k];

        t_440[k] = f_0 * gl_620[k];

        t_441[k] = f_0 * gl_621[k];

        t_442[k] = f_0 * gl_622[k];

        t_443[k] = f_0 * gl_623[k];

        t_444[k] = f_0 * gl_624[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, gl_625, gl_626, gl_627, gl_628, \
                         gl_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_0 * gl_625[k];

        t_446[k] = f_0 * gl_626[k];

        t_447[k] = f_0 * gl_627[k];

        t_448[k] = f_0 * gl_628[k];

        t_449[k] = f_0 * gl_629[k];
    }
}

auto
compute_prim_geom_10_fl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t dl, const size_t gl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fl_electron_repulsion_1_piece0(buffer, target, dl, gl, ncols, alpha);

    compute_prim_geom_10_fl_electron_repulsion_1_piece1(buffer, target, dl, gl, ncols, alpha);

    compute_prim_geom_10_fl_electron_repulsion_1_piece2(buffer, target, dl, gl, ncols, alpha);
}

static auto
compute_prim_geom_10_fl_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dl, const size_t gl,
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

    const auto *dl_0 = buffer.data(dl + 0);
    const auto *dl_1 = buffer.data(dl + 1);
    const auto *dl_2 = buffer.data(dl + 2);
    const auto *dl_3 = buffer.data(dl + 3);
    const auto *dl_4 = buffer.data(dl + 4);
    const auto *dl_5 = buffer.data(dl + 5);
    const auto *dl_6 = buffer.data(dl + 6);
    const auto *dl_7 = buffer.data(dl + 7);
    const auto *dl_8 = buffer.data(dl + 8);
    const auto *dl_9 = buffer.data(dl + 9);
    const auto *dl_10 = buffer.data(dl + 10);
    const auto *dl_11 = buffer.data(dl + 11);
    const auto *dl_12 = buffer.data(dl + 12);
    const auto *dl_13 = buffer.data(dl + 13);
    const auto *dl_14 = buffer.data(dl + 14);
    const auto *dl_15 = buffer.data(dl + 15);
    const auto *dl_16 = buffer.data(dl + 16);
    const auto *dl_17 = buffer.data(dl + 17);
    const auto *dl_18 = buffer.data(dl + 18);
    const auto *dl_19 = buffer.data(dl + 19);
    const auto *dl_20 = buffer.data(dl + 20);
    const auto *dl_21 = buffer.data(dl + 21);
    const auto *dl_22 = buffer.data(dl + 22);
    const auto *dl_23 = buffer.data(dl + 23);
    const auto *dl_24 = buffer.data(dl + 24);
    const auto *dl_25 = buffer.data(dl + 25);
    const auto *dl_26 = buffer.data(dl + 26);
    const auto *dl_27 = buffer.data(dl + 27);
    const auto *dl_28 = buffer.data(dl + 28);
    const auto *dl_29 = buffer.data(dl + 29);
    const auto *dl_30 = buffer.data(dl + 30);
    const auto *dl_31 = buffer.data(dl + 31);
    const auto *dl_32 = buffer.data(dl + 32);
    const auto *dl_33 = buffer.data(dl + 33);
    const auto *dl_34 = buffer.data(dl + 34);
    const auto *dl_35 = buffer.data(dl + 35);
    const auto *dl_36 = buffer.data(dl + 36);
    const auto *dl_37 = buffer.data(dl + 37);
    const auto *dl_38 = buffer.data(dl + 38);
    const auto *dl_39 = buffer.data(dl + 39);
    const auto *dl_40 = buffer.data(dl + 40);
    const auto *dl_41 = buffer.data(dl + 41);
    const auto *dl_42 = buffer.data(dl + 42);
    const auto *dl_43 = buffer.data(dl + 43);
    const auto *dl_44 = buffer.data(dl + 44);
    const auto *dl_45 = buffer.data(dl + 45);
    const auto *dl_46 = buffer.data(dl + 46);
    const auto *dl_47 = buffer.data(dl + 47);
    const auto *dl_48 = buffer.data(dl + 48);
    const auto *dl_49 = buffer.data(dl + 49);
    const auto *dl_50 = buffer.data(dl + 50);
    const auto *dl_51 = buffer.data(dl + 51);
    const auto *dl_52 = buffer.data(dl + 52);
    const auto *dl_53 = buffer.data(dl + 53);
    const auto *dl_54 = buffer.data(dl + 54);
    const auto *dl_55 = buffer.data(dl + 55);
    const auto *dl_56 = buffer.data(dl + 56);
    const auto *dl_57 = buffer.data(dl + 57);
    const auto *dl_58 = buffer.data(dl + 58);
    const auto *dl_59 = buffer.data(dl + 59);

    const auto *gl_90 = buffer.data(gl + 90);
    const auto *gl_91 = buffer.data(gl + 91);
    const auto *gl_92 = buffer.data(gl + 92);
    const auto *gl_93 = buffer.data(gl + 93);
    const auto *gl_94 = buffer.data(gl + 94);
    const auto *gl_95 = buffer.data(gl + 95);
    const auto *gl_96 = buffer.data(gl + 96);
    const auto *gl_97 = buffer.data(gl + 97);
    const auto *gl_98 = buffer.data(gl + 98);
    const auto *gl_99 = buffer.data(gl + 99);
    const auto *gl_100 = buffer.data(gl + 100);
    const auto *gl_101 = buffer.data(gl + 101);
    const auto *gl_102 = buffer.data(gl + 102);
    const auto *gl_103 = buffer.data(gl + 103);
    const auto *gl_104 = buffer.data(gl + 104);
    const auto *gl_105 = buffer.data(gl + 105);
    const auto *gl_106 = buffer.data(gl + 106);
    const auto *gl_107 = buffer.data(gl + 107);
    const auto *gl_108 = buffer.data(gl + 108);
    const auto *gl_109 = buffer.data(gl + 109);
    const auto *gl_110 = buffer.data(gl + 110);
    const auto *gl_111 = buffer.data(gl + 111);
    const auto *gl_112 = buffer.data(gl + 112);
    const auto *gl_113 = buffer.data(gl + 113);
    const auto *gl_114 = buffer.data(gl + 114);
    const auto *gl_115 = buffer.data(gl + 115);
    const auto *gl_116 = buffer.data(gl + 116);
    const auto *gl_117 = buffer.data(gl + 117);
    const auto *gl_118 = buffer.data(gl + 118);
    const auto *gl_119 = buffer.data(gl + 119);
    const auto *gl_120 = buffer.data(gl + 120);
    const auto *gl_121 = buffer.data(gl + 121);
    const auto *gl_122 = buffer.data(gl + 122);
    const auto *gl_123 = buffer.data(gl + 123);
    const auto *gl_124 = buffer.data(gl + 124);
    const auto *gl_125 = buffer.data(gl + 125);
    const auto *gl_126 = buffer.data(gl + 126);
    const auto *gl_127 = buffer.data(gl + 127);
    const auto *gl_128 = buffer.data(gl + 128);
    const auto *gl_129 = buffer.data(gl + 129);
    const auto *gl_130 = buffer.data(gl + 130);
    const auto *gl_131 = buffer.data(gl + 131);
    const auto *gl_132 = buffer.data(gl + 132);
    const auto *gl_133 = buffer.data(gl + 133);
    const auto *gl_134 = buffer.data(gl + 134);
    const auto *gl_180 = buffer.data(gl + 180);
    const auto *gl_181 = buffer.data(gl + 181);
    const auto *gl_182 = buffer.data(gl + 182);
    const auto *gl_183 = buffer.data(gl + 183);
    const auto *gl_184 = buffer.data(gl + 184);
    const auto *gl_185 = buffer.data(gl + 185);
    const auto *gl_186 = buffer.data(gl + 186);
    const auto *gl_187 = buffer.data(gl + 187);
    const auto *gl_188 = buffer.data(gl + 188);
    const auto *gl_189 = buffer.data(gl + 189);
    const auto *gl_190 = buffer.data(gl + 190);
    const auto *gl_191 = buffer.data(gl + 191);
    const auto *gl_192 = buffer.data(gl + 192);
    const auto *gl_193 = buffer.data(gl + 193);
    const auto *gl_194 = buffer.data(gl + 194);
    const auto *gl_195 = buffer.data(gl + 195);
    const auto *gl_196 = buffer.data(gl + 196);
    const auto *gl_197 = buffer.data(gl + 197);
    const auto *gl_198 = buffer.data(gl + 198);
    const auto *gl_199 = buffer.data(gl + 199);
    const auto *gl_200 = buffer.data(gl + 200);
    const auto *gl_201 = buffer.data(gl + 201);
    const auto *gl_202 = buffer.data(gl + 202);
    const auto *gl_203 = buffer.data(gl + 203);
    const auto *gl_204 = buffer.data(gl + 204);
    const auto *gl_205 = buffer.data(gl + 205);
    const auto *gl_206 = buffer.data(gl + 206);
    const auto *gl_207 = buffer.data(gl + 207);
    const auto *gl_208 = buffer.data(gl + 208);
    const auto *gl_209 = buffer.data(gl + 209);
    const auto *gl_210 = buffer.data(gl + 210);
    const auto *gl_211 = buffer.data(gl + 211);
    const auto *gl_212 = buffer.data(gl + 212);
    const auto *gl_213 = buffer.data(gl + 213);
    const auto *gl_214 = buffer.data(gl + 214);
    const auto *gl_215 = buffer.data(gl + 215);
    const auto *gl_216 = buffer.data(gl + 216);
    const auto *gl_217 = buffer.data(gl + 217);
    const auto *gl_218 = buffer.data(gl + 218);
    const auto *gl_219 = buffer.data(gl + 219);
    const auto *gl_220 = buffer.data(gl + 220);
    const auto *gl_221 = buffer.data(gl + 221);
    const auto *gl_222 = buffer.data(gl + 222);
    const auto *gl_223 = buffer.data(gl + 223);
    const auto *gl_224 = buffer.data(gl + 224);
    const auto *gl_225 = buffer.data(gl + 225);
    const auto *gl_226 = buffer.data(gl + 226);
    const auto *gl_227 = buffer.data(gl + 227);
    const auto *gl_228 = buffer.data(gl + 228);
    const auto *gl_229 = buffer.data(gl + 229);
    const auto *gl_230 = buffer.data(gl + 230);
    const auto *gl_231 = buffer.data(gl + 231);
    const auto *gl_232 = buffer.data(gl + 232);
    const auto *gl_233 = buffer.data(gl + 233);
    const auto *gl_234 = buffer.data(gl + 234);
    const auto *gl_235 = buffer.data(gl + 235);
    const auto *gl_236 = buffer.data(gl + 236);
    const auto *gl_237 = buffer.data(gl + 237);
    const auto *gl_238 = buffer.data(gl + 238);
    const auto *gl_239 = buffer.data(gl + 239);
    const auto *gl_240 = buffer.data(gl + 240);
    const auto *gl_241 = buffer.data(gl + 241);
    const auto *gl_242 = buffer.data(gl + 242);
    const auto *gl_243 = buffer.data(gl + 243);
    const auto *gl_244 = buffer.data(gl + 244);
    const auto *gl_245 = buffer.data(gl + 245);
    const auto *gl_246 = buffer.data(gl + 246);
    const auto *gl_247 = buffer.data(gl + 247);
    const auto *gl_248 = buffer.data(gl + 248);
    const auto *gl_249 = buffer.data(gl + 249);
    const auto *gl_250 = buffer.data(gl + 250);
    const auto *gl_251 = buffer.data(gl + 251);
    const auto *gl_252 = buffer.data(gl + 252);
    const auto *gl_253 = buffer.data(gl + 253);
    const auto *gl_254 = buffer.data(gl + 254);
    const auto *gl_255 = buffer.data(gl + 255);
    const auto *gl_256 = buffer.data(gl + 256);
    const auto *gl_257 = buffer.data(gl + 257);
    const auto *gl_258 = buffer.data(gl + 258);
    const auto *gl_259 = buffer.data(gl + 259);
    const auto *gl_260 = buffer.data(gl + 260);
    const auto *gl_261 = buffer.data(gl + 261);
    const auto *gl_262 = buffer.data(gl + 262);
    const auto *gl_263 = buffer.data(gl + 263);
    const auto *gl_264 = buffer.data(gl + 264);
    const auto *gl_265 = buffer.data(gl + 265);
    const auto *gl_266 = buffer.data(gl + 266);
    const auto *gl_267 = buffer.data(gl + 267);
    const auto *gl_268 = buffer.data(gl + 268);
    const auto *gl_269 = buffer.data(gl + 269);
    const auto *gl_315 = buffer.data(gl + 315);
    const auto *gl_316 = buffer.data(gl + 316);
    const auto *gl_317 = buffer.data(gl + 317);
    const auto *gl_318 = buffer.data(gl + 318);
    const auto *gl_319 = buffer.data(gl + 319);
    const auto *gl_320 = buffer.data(gl + 320);
    const auto *gl_321 = buffer.data(gl + 321);
    const auto *gl_322 = buffer.data(gl + 322);
    const auto *gl_323 = buffer.data(gl + 323);
    const auto *gl_324 = buffer.data(gl + 324);
    const auto *gl_325 = buffer.data(gl + 325);
    const auto *gl_326 = buffer.data(gl + 326);
    const auto *gl_327 = buffer.data(gl + 327);
    const auto *gl_328 = buffer.data(gl + 328);
    const auto *gl_329 = buffer.data(gl + 329);
    const auto *gl_330 = buffer.data(gl + 330);
    const auto *gl_331 = buffer.data(gl + 331);
    const auto *gl_332 = buffer.data(gl + 332);
    const auto *gl_333 = buffer.data(gl + 333);
    const auto *gl_334 = buffer.data(gl + 334);
    const auto *gl_335 = buffer.data(gl + 335);
    const auto *gl_336 = buffer.data(gl + 336);
    const auto *gl_337 = buffer.data(gl + 337);
    const auto *gl_338 = buffer.data(gl + 338);
    const auto *gl_339 = buffer.data(gl + 339);
    const auto *gl_340 = buffer.data(gl + 340);
    const auto *gl_341 = buffer.data(gl + 341);
    const auto *gl_342 = buffer.data(gl + 342);
    const auto *gl_343 = buffer.data(gl + 343);
    const auto *gl_344 = buffer.data(gl + 344);
    const auto *gl_345 = buffer.data(gl + 345);
    const auto *gl_346 = buffer.data(gl + 346);
    const auto *gl_347 = buffer.data(gl + 347);
    const auto *gl_348 = buffer.data(gl + 348);
    const auto *gl_349 = buffer.data(gl + 349);
    const auto *gl_350 = buffer.data(gl + 350);
    const auto *gl_351 = buffer.data(gl + 351);
    const auto *gl_352 = buffer.data(gl + 352);
    const auto *gl_353 = buffer.data(gl + 353);
    const auto *gl_354 = buffer.data(gl + 354);
    const auto *gl_355 = buffer.data(gl + 355);
    const auto *gl_356 = buffer.data(gl + 356);
    const auto *gl_357 = buffer.data(gl + 357);
    const auto *gl_358 = buffer.data(gl + 358);
    const auto *gl_359 = buffer.data(gl + 359);
    const auto *gl_360 = buffer.data(gl + 360);
    const auto *gl_361 = buffer.data(gl + 361);
    const auto *gl_362 = buffer.data(gl + 362);
    const auto *gl_363 = buffer.data(gl + 363);
    const auto *gl_364 = buffer.data(gl + 364);
    const auto *gl_365 = buffer.data(gl + 365);
    const auto *gl_366 = buffer.data(gl + 366);
    const auto *gl_367 = buffer.data(gl + 367);
    const auto *gl_368 = buffer.data(gl + 368);
    const auto *gl_369 = buffer.data(gl + 369);
    const auto *gl_370 = buffer.data(gl + 370);
    const auto *gl_371 = buffer.data(gl + 371);
    const auto *gl_372 = buffer.data(gl + 372);
    const auto *gl_373 = buffer.data(gl + 373);
    const auto *gl_374 = buffer.data(gl + 374);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gl_90, gl_91, gl_92, gl_93, \
                         gl_94, gl_95, gl_96, gl_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gl_90[k];

        t_1[k] = f_0 * gl_91[k];

        t_2[k] = f_0 * gl_92[k];

        t_3[k] = f_0 * gl_93[k];

        t_4[k] = f_0 * gl_94[k];

        t_5[k] = f_0 * gl_95[k];

        t_6[k] = f_0 * gl_96[k];

        t_7[k] = f_0 * gl_97[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gl_98, gl_99, gl_100, \
                         gl_101, gl_102, gl_103, gl_104, gl_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gl_98[k];

        t_9[k] = f_0 * gl_99[k];

        t_10[k] = f_0 * gl_100[k];

        t_11[k] = f_0 * gl_101[k];

        t_12[k] = f_0 * gl_102[k];

        t_13[k] = f_0 * gl_103[k];

        t_14[k] = f_0 * gl_104[k];

        t_15[k] = f_0 * gl_105[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, gl_106, gl_107, \
                         gl_108, gl_109, gl_110, gl_111, gl_112, \
                         gl_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gl_106[k];

        t_17[k] = f_0 * gl_107[k];

        t_18[k] = f_0 * gl_108[k];

        t_19[k] = f_0 * gl_109[k];

        t_20[k] = f_0 * gl_110[k];

        t_21[k] = f_0 * gl_111[k];

        t_22[k] = f_0 * gl_112[k];

        t_23[k] = f_0 * gl_113[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, gl_114, gl_115, \
                         gl_116, gl_117, gl_118, gl_119, gl_120, \
                         gl_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * gl_114[k];

        t_25[k] = f_0 * gl_115[k];

        t_26[k] = f_0 * gl_116[k];

        t_27[k] = f_0 * gl_117[k];

        t_28[k] = f_0 * gl_118[k];

        t_29[k] = f_0 * gl_119[k];

        t_30[k] = f_0 * gl_120[k];

        t_31[k] = f_0 * gl_121[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, gl_122, gl_123, \
                         gl_124, gl_125, gl_126, gl_127, gl_128, \
                         gl_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * gl_122[k];

        t_33[k] = f_0 * gl_123[k];

        t_34[k] = f_0 * gl_124[k];

        t_35[k] = f_0 * gl_125[k];

        t_36[k] = f_0 * gl_126[k];

        t_37[k] = f_0 * gl_127[k];

        t_38[k] = f_0 * gl_128[k];

        t_39[k] = f_0 * gl_129[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, gl_130, gl_131, \
                         gl_132, gl_133, gl_134, gl_180, gl_181, \
                         gl_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * gl_130[k];

        t_41[k] = f_0 * gl_131[k];

        t_42[k] = f_0 * gl_132[k];

        t_43[k] = f_0 * gl_133[k];

        t_44[k] = f_0 * gl_134[k];

        t_45[k] = f_0 * gl_180[k];

        t_46[k] = f_0 * gl_181[k];

        t_47[k] = f_0 * gl_182[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, gl_183, gl_184, \
                         gl_185, gl_186, gl_187, gl_188, gl_189, \
                         gl_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * gl_183[k];

        t_49[k] = f_0 * gl_184[k];

        t_50[k] = f_0 * gl_185[k];

        t_51[k] = f_0 * gl_186[k];

        t_52[k] = f_0 * gl_187[k];

        t_53[k] = f_0 * gl_188[k];

        t_54[k] = f_0 * gl_189[k];

        t_55[k] = f_0 * gl_190[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, t_63, gl_191, gl_192, \
                         gl_193, gl_194, gl_195, gl_196, gl_197, \
                         gl_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * gl_191[k];

        t_57[k] = f_0 * gl_192[k];

        t_58[k] = f_0 * gl_193[k];

        t_59[k] = f_0 * gl_194[k];

        t_60[k] = f_0 * gl_195[k];

        t_61[k] = f_0 * gl_196[k];

        t_62[k] = f_0 * gl_197[k];

        t_63[k] = f_0 * gl_198[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, t_69, t_70, t_71, gl_199, gl_200, \
                         gl_201, gl_202, gl_203, gl_204, gl_205, \
                         gl_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * gl_199[k];

        t_65[k] = f_0 * gl_200[k];

        t_66[k] = f_0 * gl_201[k];

        t_67[k] = f_0 * gl_202[k];

        t_68[k] = f_0 * gl_203[k];

        t_69[k] = f_0 * gl_204[k];

        t_70[k] = f_0 * gl_205[k];

        t_71[k] = f_0 * gl_206[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, t_78, t_79, gl_207, gl_208, \
                         gl_209, gl_210, gl_211, gl_212, gl_213, \
                         gl_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * gl_207[k];

        t_73[k] = f_0 * gl_208[k];

        t_74[k] = f_0 * gl_209[k];

        t_75[k] = f_0 * gl_210[k];

        t_76[k] = f_0 * gl_211[k];

        t_77[k] = f_0 * gl_212[k];

        t_78[k] = f_0 * gl_213[k];

        t_79[k] = f_0 * gl_214[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, t_86, t_87, gl_215, gl_216, \
                         gl_217, gl_218, gl_219, gl_220, gl_221, \
                         gl_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * gl_215[k];

        t_81[k] = f_0 * gl_216[k];

        t_82[k] = f_0 * gl_217[k];

        t_83[k] = f_0 * gl_218[k];

        t_84[k] = f_0 * gl_219[k];

        t_85[k] = f_0 * gl_220[k];

        t_86[k] = f_0 * gl_221[k];

        t_87[k] = f_0 * gl_222[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, t_93, dl_0, dl_1, dl_2, dl_3, gl_223, \
                         gl_224, gl_225, gl_226, gl_227, gl_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * gl_223[k];

        t_89[k] = f_0 * gl_224[k];

        t_90[k] = -dl_0[k]
                  + f_0 * gl_225[k];

        t_91[k] = -dl_1[k]
                  + f_0 * gl_226[k];

        t_92[k] = -dl_2[k]
                  + f_0 * gl_227[k];

        t_93[k] = -dl_3[k]
                  + f_0 * gl_228[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, dl_4, dl_5, dl_6, dl_7, dl_8, gl_229, \
                         gl_230, gl_231, gl_232, gl_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -dl_4[k]
                  + f_0 * gl_229[k];

        t_95[k] = -dl_5[k]
                  + f_0 * gl_230[k];

        t_96[k] = -dl_6[k]
                  + f_0 * gl_231[k];

        t_97[k] = -dl_7[k]
                  + f_0 * gl_232[k];

        t_98[k] = -dl_8[k]
                  + f_0 * gl_233[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, dl_9, dl_10, dl_11, dl_12, dl_13, \
                         gl_234, gl_235, gl_236, gl_237, gl_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -dl_9[k]
                  + f_0 * gl_234[k];

        t_100[k] = -dl_10[k]
                   + f_0 * gl_235[k];

        t_101[k] = -dl_11[k]
                   + f_0 * gl_236[k];

        t_102[k] = -dl_12[k]
                   + f_0 * gl_237[k];

        t_103[k] = -dl_13[k]
                   + f_0 * gl_238[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, dl_14, dl_15, dl_16, dl_17, dl_18, \
                         gl_239, gl_240, gl_241, gl_242, gl_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -dl_14[k]
                   + f_0 * gl_239[k];

        t_105[k] = -dl_15[k]
                   + f_0 * gl_240[k];

        t_106[k] = -dl_16[k]
                   + f_0 * gl_241[k];

        t_107[k] = -dl_17[k]
                   + f_0 * gl_242[k];

        t_108[k] = -dl_18[k]
                   + f_0 * gl_243[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, dl_19, dl_20, dl_21, dl_22, dl_23, \
                         gl_244, gl_245, gl_246, gl_247, gl_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -dl_19[k]
                   + f_0 * gl_244[k];

        t_110[k] = -dl_20[k]
                   + f_0 * gl_245[k];

        t_111[k] = -dl_21[k]
                   + f_0 * gl_246[k];

        t_112[k] = -dl_22[k]
                   + f_0 * gl_247[k];

        t_113[k] = -dl_23[k]
                   + f_0 * gl_248[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, dl_24, dl_25, dl_26, dl_27, dl_28, \
                         gl_249, gl_250, gl_251, gl_252, gl_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -dl_24[k]
                   + f_0 * gl_249[k];

        t_115[k] = -dl_25[k]
                   + f_0 * gl_250[k];

        t_116[k] = -dl_26[k]
                   + f_0 * gl_251[k];

        t_117[k] = -dl_27[k]
                   + f_0 * gl_252[k];

        t_118[k] = -dl_28[k]
                   + f_0 * gl_253[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, dl_29, dl_30, dl_31, dl_32, dl_33, \
                         gl_254, gl_255, gl_256, gl_257, gl_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -dl_29[k]
                   + f_0 * gl_254[k];

        t_120[k] = -dl_30[k]
                   + f_0 * gl_255[k];

        t_121[k] = -dl_31[k]
                   + f_0 * gl_256[k];

        t_122[k] = -dl_32[k]
                   + f_0 * gl_257[k];

        t_123[k] = -dl_33[k]
                   + f_0 * gl_258[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, dl_34, dl_35, dl_36, dl_37, dl_38, \
                         gl_259, gl_260, gl_261, gl_262, gl_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -dl_34[k]
                   + f_0 * gl_259[k];

        t_125[k] = -dl_35[k]
                   + f_0 * gl_260[k];

        t_126[k] = -dl_36[k]
                   + f_0 * gl_261[k];

        t_127[k] = -dl_37[k]
                   + f_0 * gl_262[k];

        t_128[k] = -dl_38[k]
                   + f_0 * gl_263[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, dl_39, dl_40, dl_41, dl_42, dl_43, \
                         gl_264, gl_265, gl_266, gl_267, gl_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -dl_39[k]
                   + f_0 * gl_264[k];

        t_130[k] = -dl_40[k]
                   + f_0 * gl_265[k];

        t_131[k] = -dl_41[k]
                   + f_0 * gl_266[k];

        t_132[k] = -dl_42[k]
                   + f_0 * gl_267[k];

        t_133[k] = -dl_43[k]
                   + f_0 * gl_268[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, t_139, t_140, dl_44, gl_269, \
                         gl_315, gl_316, gl_317, gl_318, gl_319, \
                         gl_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -dl_44[k]
                   + f_0 * gl_269[k];

        t_135[k] = f_0 * gl_315[k];

        t_136[k] = f_0 * gl_316[k];

        t_137[k] = f_0 * gl_317[k];

        t_138[k] = f_0 * gl_318[k];

        t_139[k] = f_0 * gl_319[k];

        t_140[k] = f_0 * gl_320[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, t_146, t_147, t_148, gl_321, \
                         gl_322, gl_323, gl_324, gl_325, gl_326, gl_327, \
                         gl_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * gl_321[k];

        t_142[k] = f_0 * gl_322[k];

        t_143[k] = f_0 * gl_323[k];

        t_144[k] = f_0 * gl_324[k];

        t_145[k] = f_0 * gl_325[k];

        t_146[k] = f_0 * gl_326[k];

        t_147[k] = f_0 * gl_327[k];

        t_148[k] = f_0 * gl_328[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, t_154, t_155, t_156, gl_329, \
                         gl_330, gl_331, gl_332, gl_333, gl_334, gl_335, \
                         gl_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * gl_329[k];

        t_150[k] = f_0 * gl_330[k];

        t_151[k] = f_0 * gl_331[k];

        t_152[k] = f_0 * gl_332[k];

        t_153[k] = f_0 * gl_333[k];

        t_154[k] = f_0 * gl_334[k];

        t_155[k] = f_0 * gl_335[k];

        t_156[k] = f_0 * gl_336[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, t_162, t_163, t_164, gl_337, \
                         gl_338, gl_339, gl_340, gl_341, gl_342, gl_343, \
                         gl_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_0 * gl_337[k];

        t_158[k] = f_0 * gl_338[k];

        t_159[k] = f_0 * gl_339[k];

        t_160[k] = f_0 * gl_340[k];

        t_161[k] = f_0 * gl_341[k];

        t_162[k] = f_0 * gl_342[k];

        t_163[k] = f_0 * gl_343[k];

        t_164[k] = f_0 * gl_344[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, t_170, t_171, t_172, gl_345, \
                         gl_346, gl_347, gl_348, gl_349, gl_350, gl_351, \
                         gl_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_0 * gl_345[k];

        t_166[k] = f_0 * gl_346[k];

        t_167[k] = f_0 * gl_347[k];

        t_168[k] = f_0 * gl_348[k];

        t_169[k] = f_0 * gl_349[k];

        t_170[k] = f_0 * gl_350[k];

        t_171[k] = f_0 * gl_351[k];

        t_172[k] = f_0 * gl_352[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, gl_353, gl_354, \
                         gl_355, gl_356, gl_357, gl_358, gl_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * gl_353[k];

        t_174[k] = f_0 * gl_354[k];

        t_175[k] = f_0 * gl_355[k];

        t_176[k] = f_0 * gl_356[k];

        t_177[k] = f_0 * gl_357[k];

        t_178[k] = f_0 * gl_358[k];

        t_179[k] = f_0 * gl_359[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, dl_45, dl_46, dl_47, dl_48, dl_49, \
                         gl_360, gl_361, gl_362, gl_363, gl_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -dl_45[k]
                   + f_0 * gl_360[k];

        t_181[k] = -dl_46[k]
                   + f_0 * gl_361[k];

        t_182[k] = -dl_47[k]
                   + f_0 * gl_362[k];

        t_183[k] = -dl_48[k]
                   + f_0 * gl_363[k];

        t_184[k] = -dl_49[k]
                   + f_0 * gl_364[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, dl_50, dl_51, dl_52, dl_53, dl_54, \
                         gl_365, gl_366, gl_367, gl_368, gl_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -dl_50[k]
                   + f_0 * gl_365[k];

        t_186[k] = -dl_51[k]
                   + f_0 * gl_366[k];

        t_187[k] = -dl_52[k]
                   + f_0 * gl_367[k];

        t_188[k] = -dl_53[k]
                   + f_0 * gl_368[k];

        t_189[k] = -dl_54[k]
                   + f_0 * gl_369[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, dl_55, dl_56, dl_57, dl_58, dl_59, \
                         gl_370, gl_371, gl_372, gl_373, gl_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -dl_55[k]
                   + f_0 * gl_370[k];

        t_191[k] = -dl_56[k]
                   + f_0 * gl_371[k];

        t_192[k] = -dl_57[k]
                   + f_0 * gl_372[k];

        t_193[k] = -dl_58[k]
                   + f_0 * gl_373[k];

        t_194[k] = -dl_59[k]
                   + f_0 * gl_374[k];
    }
}

static auto
compute_prim_geom_10_fl_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dl, const size_t gl,
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

    const auto *dl_60 = buffer.data(dl + 60);
    const auto *dl_61 = buffer.data(dl + 61);
    const auto *dl_62 = buffer.data(dl + 62);
    const auto *dl_63 = buffer.data(dl + 63);
    const auto *dl_64 = buffer.data(dl + 64);
    const auto *dl_65 = buffer.data(dl + 65);
    const auto *dl_66 = buffer.data(dl + 66);
    const auto *dl_67 = buffer.data(dl + 67);
    const auto *dl_68 = buffer.data(dl + 68);
    const auto *dl_69 = buffer.data(dl + 69);
    const auto *dl_70 = buffer.data(dl + 70);
    const auto *dl_71 = buffer.data(dl + 71);
    const auto *dl_72 = buffer.data(dl + 72);
    const auto *dl_73 = buffer.data(dl + 73);
    const auto *dl_74 = buffer.data(dl + 74);
    const auto *dl_75 = buffer.data(dl + 75);
    const auto *dl_76 = buffer.data(dl + 76);
    const auto *dl_77 = buffer.data(dl + 77);
    const auto *dl_78 = buffer.data(dl + 78);
    const auto *dl_79 = buffer.data(dl + 79);
    const auto *dl_80 = buffer.data(dl + 80);
    const auto *dl_81 = buffer.data(dl + 81);
    const auto *dl_82 = buffer.data(dl + 82);
    const auto *dl_83 = buffer.data(dl + 83);
    const auto *dl_84 = buffer.data(dl + 84);
    const auto *dl_85 = buffer.data(dl + 85);
    const auto *dl_86 = buffer.data(dl + 86);
    const auto *dl_87 = buffer.data(dl + 87);
    const auto *dl_88 = buffer.data(dl + 88);
    const auto *dl_89 = buffer.data(dl + 89);
    const auto *dl_90 = buffer.data(dl + 90);
    const auto *dl_91 = buffer.data(dl + 91);
    const auto *dl_92 = buffer.data(dl + 92);
    const auto *dl_93 = buffer.data(dl + 93);
    const auto *dl_94 = buffer.data(dl + 94);
    const auto *dl_95 = buffer.data(dl + 95);
    const auto *dl_96 = buffer.data(dl + 96);
    const auto *dl_97 = buffer.data(dl + 97);
    const auto *dl_98 = buffer.data(dl + 98);
    const auto *dl_99 = buffer.data(dl + 99);
    const auto *dl_100 = buffer.data(dl + 100);
    const auto *dl_101 = buffer.data(dl + 101);
    const auto *dl_102 = buffer.data(dl + 102);
    const auto *dl_103 = buffer.data(dl + 103);
    const auto *dl_104 = buffer.data(dl + 104);
    const auto *dl_105 = buffer.data(dl + 105);
    const auto *dl_106 = buffer.data(dl + 106);
    const auto *dl_107 = buffer.data(dl + 107);
    const auto *dl_108 = buffer.data(dl + 108);
    const auto *dl_109 = buffer.data(dl + 109);
    const auto *dl_110 = buffer.data(dl + 110);
    const auto *dl_111 = buffer.data(dl + 111);
    const auto *dl_112 = buffer.data(dl + 112);
    const auto *dl_113 = buffer.data(dl + 113);
    const auto *dl_114 = buffer.data(dl + 114);
    const auto *dl_115 = buffer.data(dl + 115);
    const auto *dl_116 = buffer.data(dl + 116);
    const auto *dl_117 = buffer.data(dl + 117);
    const auto *dl_118 = buffer.data(dl + 118);
    const auto *dl_119 = buffer.data(dl + 119);
    const auto *dl_120 = buffer.data(dl + 120);
    const auto *dl_121 = buffer.data(dl + 121);
    const auto *dl_122 = buffer.data(dl + 122);
    const auto *dl_123 = buffer.data(dl + 123);
    const auto *dl_124 = buffer.data(dl + 124);
    const auto *dl_125 = buffer.data(dl + 125);
    const auto *dl_126 = buffer.data(dl + 126);
    const auto *dl_127 = buffer.data(dl + 127);
    const auto *dl_128 = buffer.data(dl + 128);
    const auto *dl_129 = buffer.data(dl + 129);
    const auto *dl_130 = buffer.data(dl + 130);
    const auto *dl_131 = buffer.data(dl + 131);
    const auto *dl_132 = buffer.data(dl + 132);
    const auto *dl_133 = buffer.data(dl + 133);
    const auto *dl_134 = buffer.data(dl + 134);
    const auto *dl_135 = buffer.data(dl + 135);
    const auto *dl_136 = buffer.data(dl + 136);
    const auto *dl_137 = buffer.data(dl + 137);
    const auto *dl_138 = buffer.data(dl + 138);
    const auto *dl_139 = buffer.data(dl + 139);
    const auto *dl_140 = buffer.data(dl + 140);
    const auto *dl_141 = buffer.data(dl + 141);
    const auto *dl_142 = buffer.data(dl + 142);
    const auto *dl_143 = buffer.data(dl + 143);
    const auto *dl_144 = buffer.data(dl + 144);
    const auto *dl_145 = buffer.data(dl + 145);
    const auto *dl_146 = buffer.data(dl + 146);
    const auto *dl_147 = buffer.data(dl + 147);
    const auto *dl_148 = buffer.data(dl + 148);
    const auto *dl_149 = buffer.data(dl + 149);
    const auto *dl_150 = buffer.data(dl + 150);
    const auto *dl_151 = buffer.data(dl + 151);
    const auto *dl_152 = buffer.data(dl + 152);
    const auto *dl_153 = buffer.data(dl + 153);
    const auto *dl_154 = buffer.data(dl + 154);
    const auto *dl_155 = buffer.data(dl + 155);
    const auto *dl_156 = buffer.data(dl + 156);
    const auto *dl_157 = buffer.data(dl + 157);
    const auto *dl_158 = buffer.data(dl + 158);
    const auto *dl_159 = buffer.data(dl + 159);
    const auto *dl_160 = buffer.data(dl + 160);
    const auto *dl_161 = buffer.data(dl + 161);
    const auto *dl_162 = buffer.data(dl + 162);
    const auto *dl_163 = buffer.data(dl + 163);
    const auto *dl_164 = buffer.data(dl + 164);
    const auto *dl_165 = buffer.data(dl + 165);
    const auto *dl_166 = buffer.data(dl + 166);
    const auto *dl_167 = buffer.data(dl + 167);
    const auto *dl_168 = buffer.data(dl + 168);
    const auto *dl_169 = buffer.data(dl + 169);
    const auto *dl_170 = buffer.data(dl + 170);
    const auto *dl_171 = buffer.data(dl + 171);
    const auto *dl_172 = buffer.data(dl + 172);
    const auto *dl_173 = buffer.data(dl + 173);
    const auto *dl_174 = buffer.data(dl + 174);
    const auto *dl_175 = buffer.data(dl + 175);
    const auto *dl_176 = buffer.data(dl + 176);

    const auto *gl_375 = buffer.data(gl + 375);
    const auto *gl_376 = buffer.data(gl + 376);
    const auto *gl_377 = buffer.data(gl + 377);
    const auto *gl_378 = buffer.data(gl + 378);
    const auto *gl_379 = buffer.data(gl + 379);
    const auto *gl_380 = buffer.data(gl + 380);
    const auto *gl_381 = buffer.data(gl + 381);
    const auto *gl_382 = buffer.data(gl + 382);
    const auto *gl_383 = buffer.data(gl + 383);
    const auto *gl_384 = buffer.data(gl + 384);
    const auto *gl_385 = buffer.data(gl + 385);
    const auto *gl_386 = buffer.data(gl + 386);
    const auto *gl_387 = buffer.data(gl + 387);
    const auto *gl_388 = buffer.data(gl + 388);
    const auto *gl_389 = buffer.data(gl + 389);
    const auto *gl_390 = buffer.data(gl + 390);
    const auto *gl_391 = buffer.data(gl + 391);
    const auto *gl_392 = buffer.data(gl + 392);
    const auto *gl_393 = buffer.data(gl + 393);
    const auto *gl_394 = buffer.data(gl + 394);
    const auto *gl_395 = buffer.data(gl + 395);
    const auto *gl_396 = buffer.data(gl + 396);
    const auto *gl_397 = buffer.data(gl + 397);
    const auto *gl_398 = buffer.data(gl + 398);
    const auto *gl_399 = buffer.data(gl + 399);
    const auto *gl_400 = buffer.data(gl + 400);
    const auto *gl_401 = buffer.data(gl + 401);
    const auto *gl_402 = buffer.data(gl + 402);
    const auto *gl_403 = buffer.data(gl + 403);
    const auto *gl_404 = buffer.data(gl + 404);
    const auto *gl_405 = buffer.data(gl + 405);
    const auto *gl_406 = buffer.data(gl + 406);
    const auto *gl_407 = buffer.data(gl + 407);
    const auto *gl_408 = buffer.data(gl + 408);
    const auto *gl_409 = buffer.data(gl + 409);
    const auto *gl_410 = buffer.data(gl + 410);
    const auto *gl_411 = buffer.data(gl + 411);
    const auto *gl_412 = buffer.data(gl + 412);
    const auto *gl_413 = buffer.data(gl + 413);
    const auto *gl_414 = buffer.data(gl + 414);
    const auto *gl_415 = buffer.data(gl + 415);
    const auto *gl_416 = buffer.data(gl + 416);
    const auto *gl_417 = buffer.data(gl + 417);
    const auto *gl_418 = buffer.data(gl + 418);
    const auto *gl_419 = buffer.data(gl + 419);
    const auto *gl_420 = buffer.data(gl + 420);
    const auto *gl_421 = buffer.data(gl + 421);
    const auto *gl_422 = buffer.data(gl + 422);
    const auto *gl_423 = buffer.data(gl + 423);
    const auto *gl_424 = buffer.data(gl + 424);
    const auto *gl_425 = buffer.data(gl + 425);
    const auto *gl_426 = buffer.data(gl + 426);
    const auto *gl_427 = buffer.data(gl + 427);
    const auto *gl_428 = buffer.data(gl + 428);
    const auto *gl_429 = buffer.data(gl + 429);
    const auto *gl_430 = buffer.data(gl + 430);
    const auto *gl_431 = buffer.data(gl + 431);
    const auto *gl_432 = buffer.data(gl + 432);
    const auto *gl_433 = buffer.data(gl + 433);
    const auto *gl_434 = buffer.data(gl + 434);
    const auto *gl_435 = buffer.data(gl + 435);
    const auto *gl_436 = buffer.data(gl + 436);
    const auto *gl_437 = buffer.data(gl + 437);
    const auto *gl_438 = buffer.data(gl + 438);
    const auto *gl_439 = buffer.data(gl + 439);
    const auto *gl_440 = buffer.data(gl + 440);
    const auto *gl_441 = buffer.data(gl + 441);
    const auto *gl_442 = buffer.data(gl + 442);
    const auto *gl_443 = buffer.data(gl + 443);
    const auto *gl_444 = buffer.data(gl + 444);
    const auto *gl_445 = buffer.data(gl + 445);
    const auto *gl_446 = buffer.data(gl + 446);
    const auto *gl_447 = buffer.data(gl + 447);
    const auto *gl_448 = buffer.data(gl + 448);
    const auto *gl_449 = buffer.data(gl + 449);
    const auto *gl_495 = buffer.data(gl + 495);
    const auto *gl_496 = buffer.data(gl + 496);
    const auto *gl_497 = buffer.data(gl + 497);
    const auto *gl_498 = buffer.data(gl + 498);
    const auto *gl_499 = buffer.data(gl + 499);
    const auto *gl_500 = buffer.data(gl + 500);
    const auto *gl_501 = buffer.data(gl + 501);
    const auto *gl_502 = buffer.data(gl + 502);
    const auto *gl_503 = buffer.data(gl + 503);
    const auto *gl_504 = buffer.data(gl + 504);
    const auto *gl_505 = buffer.data(gl + 505);
    const auto *gl_506 = buffer.data(gl + 506);
    const auto *gl_507 = buffer.data(gl + 507);
    const auto *gl_508 = buffer.data(gl + 508);
    const auto *gl_509 = buffer.data(gl + 509);
    const auto *gl_510 = buffer.data(gl + 510);
    const auto *gl_511 = buffer.data(gl + 511);
    const auto *gl_512 = buffer.data(gl + 512);
    const auto *gl_513 = buffer.data(gl + 513);
    const auto *gl_514 = buffer.data(gl + 514);
    const auto *gl_515 = buffer.data(gl + 515);
    const auto *gl_516 = buffer.data(gl + 516);
    const auto *gl_517 = buffer.data(gl + 517);
    const auto *gl_518 = buffer.data(gl + 518);
    const auto *gl_519 = buffer.data(gl + 519);
    const auto *gl_520 = buffer.data(gl + 520);
    const auto *gl_521 = buffer.data(gl + 521);
    const auto *gl_522 = buffer.data(gl + 522);
    const auto *gl_523 = buffer.data(gl + 523);
    const auto *gl_524 = buffer.data(gl + 524);
    const auto *gl_525 = buffer.data(gl + 525);
    const auto *gl_526 = buffer.data(gl + 526);
    const auto *gl_527 = buffer.data(gl + 527);
    const auto *gl_528 = buffer.data(gl + 528);
    const auto *gl_529 = buffer.data(gl + 529);
    const auto *gl_530 = buffer.data(gl + 530);
    const auto *gl_531 = buffer.data(gl + 531);
    const auto *gl_532 = buffer.data(gl + 532);
    const auto *gl_533 = buffer.data(gl + 533);
    const auto *gl_534 = buffer.data(gl + 534);
    const auto *gl_535 = buffer.data(gl + 535);
    const auto *gl_536 = buffer.data(gl + 536);
    const auto *gl_537 = buffer.data(gl + 537);
    const auto *gl_538 = buffer.data(gl + 538);
    const auto *gl_539 = buffer.data(gl + 539);
    const auto *gl_540 = buffer.data(gl + 540);
    const auto *gl_541 = buffer.data(gl + 541);
    const auto *gl_542 = buffer.data(gl + 542);
    const auto *gl_543 = buffer.data(gl + 543);
    const auto *gl_544 = buffer.data(gl + 544);
    const auto *gl_545 = buffer.data(gl + 545);
    const auto *gl_546 = buffer.data(gl + 546);
    const auto *gl_547 = buffer.data(gl + 547);
    const auto *gl_548 = buffer.data(gl + 548);
    const auto *gl_549 = buffer.data(gl + 549);
    const auto *gl_550 = buffer.data(gl + 550);
    const auto *gl_551 = buffer.data(gl + 551);
    const auto *gl_552 = buffer.data(gl + 552);
    const auto *gl_553 = buffer.data(gl + 553);
    const auto *gl_554 = buffer.data(gl + 554);
    const auto *gl_555 = buffer.data(gl + 555);
    const auto *gl_556 = buffer.data(gl + 556);
    const auto *gl_557 = buffer.data(gl + 557);
    const auto *gl_558 = buffer.data(gl + 558);
    const auto *gl_559 = buffer.data(gl + 559);
    const auto *gl_560 = buffer.data(gl + 560);
    const auto *gl_561 = buffer.data(gl + 561);
    const auto *gl_562 = buffer.data(gl + 562);
    const auto *gl_563 = buffer.data(gl + 563);
    const auto *gl_564 = buffer.data(gl + 564);
    const auto *gl_565 = buffer.data(gl + 565);
    const auto *gl_566 = buffer.data(gl + 566);
    const auto *gl_567 = buffer.data(gl + 567);
    const auto *gl_568 = buffer.data(gl + 568);
    const auto *gl_569 = buffer.data(gl + 569);
    const auto *gl_570 = buffer.data(gl + 570);
    const auto *gl_571 = buffer.data(gl + 571);
    const auto *gl_572 = buffer.data(gl + 572);
    const auto *gl_573 = buffer.data(gl + 573);
    const auto *gl_574 = buffer.data(gl + 574);
    const auto *gl_575 = buffer.data(gl + 575);
    const auto *gl_576 = buffer.data(gl + 576);
    const auto *gl_577 = buffer.data(gl + 577);
    const auto *gl_578 = buffer.data(gl + 578);
    const auto *gl_579 = buffer.data(gl + 579);
    const auto *gl_580 = buffer.data(gl + 580);
    const auto *gl_581 = buffer.data(gl + 581);

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, dl_60, dl_61, dl_62, dl_63, dl_64, \
                         gl_375, gl_376, gl_377, gl_378, gl_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -dl_60[k]
                   + f_0 * gl_375[k];

        t_196[k] = -dl_61[k]
                   + f_0 * gl_376[k];

        t_197[k] = -dl_62[k]
                   + f_0 * gl_377[k];

        t_198[k] = -dl_63[k]
                   + f_0 * gl_378[k];

        t_199[k] = -dl_64[k]
                   + f_0 * gl_379[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, dl_65, dl_66, dl_67, dl_68, dl_69, \
                         gl_380, gl_381, gl_382, gl_383, gl_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -dl_65[k]
                   + f_0 * gl_380[k];

        t_201[k] = -dl_66[k]
                   + f_0 * gl_381[k];

        t_202[k] = -dl_67[k]
                   + f_0 * gl_382[k];

        t_203[k] = -dl_68[k]
                   + f_0 * gl_383[k];

        t_204[k] = -dl_69[k]
                   + f_0 * gl_384[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, dl_70, dl_71, dl_72, dl_73, dl_74, \
                         gl_385, gl_386, gl_387, gl_388, gl_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -dl_70[k]
                   + f_0 * gl_385[k];

        t_206[k] = -dl_71[k]
                   + f_0 * gl_386[k];

        t_207[k] = -dl_72[k]
                   + f_0 * gl_387[k];

        t_208[k] = -dl_73[k]
                   + f_0 * gl_388[k];

        t_209[k] = -dl_74[k]
                   + f_0 * gl_389[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, dl_75, dl_76, dl_77, dl_78, dl_79, \
                         gl_390, gl_391, gl_392, gl_393, gl_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -dl_75[k]
                   + f_0 * gl_390[k];

        t_211[k] = -dl_76[k]
                   + f_0 * gl_391[k];

        t_212[k] = -dl_77[k]
                   + f_0 * gl_392[k];

        t_213[k] = -dl_78[k]
                   + f_0 * gl_393[k];

        t_214[k] = -dl_79[k]
                   + f_0 * gl_394[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, dl_80, dl_81, dl_82, dl_83, dl_84, \
                         gl_395, gl_396, gl_397, gl_398, gl_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -dl_80[k]
                   + f_0 * gl_395[k];

        t_216[k] = -dl_81[k]
                   + f_0 * gl_396[k];

        t_217[k] = -dl_82[k]
                   + f_0 * gl_397[k];

        t_218[k] = -dl_83[k]
                   + f_0 * gl_398[k];

        t_219[k] = -dl_84[k]
                   + f_0 * gl_399[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, dl_85, dl_86, dl_87, dl_88, dl_89, \
                         gl_400, gl_401, gl_402, gl_403, gl_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -dl_85[k]
                   + f_0 * gl_400[k];

        t_221[k] = -dl_86[k]
                   + f_0 * gl_401[k];

        t_222[k] = -dl_87[k]
                   + f_0 * gl_402[k];

        t_223[k] = -dl_88[k]
                   + f_0 * gl_403[k];

        t_224[k] = -dl_89[k]
                   + f_0 * gl_404[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, dl_90, dl_91, dl_92, dl_93, dl_94, \
                         gl_405, gl_406, gl_407, gl_408, gl_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * dl_90[k]
                   + f_0 * gl_405[k];

        t_226[k] = -2.0 * dl_91[k]
                   + f_0 * gl_406[k];

        t_227[k] = -2.0 * dl_92[k]
                   + f_0 * gl_407[k];

        t_228[k] = -2.0 * dl_93[k]
                   + f_0 * gl_408[k];

        t_229[k] = -2.0 * dl_94[k]
                   + f_0 * gl_409[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, dl_95, dl_96, dl_97, dl_98, dl_99, \
                         gl_410, gl_411, gl_412, gl_413, gl_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * dl_95[k]
                   + f_0 * gl_410[k];

        t_231[k] = -2.0 * dl_96[k]
                   + f_0 * gl_411[k];

        t_232[k] = -2.0 * dl_97[k]
                   + f_0 * gl_412[k];

        t_233[k] = -2.0 * dl_98[k]
                   + f_0 * gl_413[k];

        t_234[k] = -2.0 * dl_99[k]
                   + f_0 * gl_414[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, dl_100, dl_101, dl_102, dl_103, \
                         dl_104, gl_415, gl_416, gl_417, gl_418, \
                         gl_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * dl_100[k]
                   + f_0 * gl_415[k];

        t_236[k] = -2.0 * dl_101[k]
                   + f_0 * gl_416[k];

        t_237[k] = -2.0 * dl_102[k]
                   + f_0 * gl_417[k];

        t_238[k] = -2.0 * dl_103[k]
                   + f_0 * gl_418[k];

        t_239[k] = -2.0 * dl_104[k]
                   + f_0 * gl_419[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, dl_105, dl_106, dl_107, dl_108, \
                         dl_109, gl_420, gl_421, gl_422, gl_423, \
                         gl_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * dl_105[k]
                   + f_0 * gl_420[k];

        t_241[k] = -2.0 * dl_106[k]
                   + f_0 * gl_421[k];

        t_242[k] = -2.0 * dl_107[k]
                   + f_0 * gl_422[k];

        t_243[k] = -2.0 * dl_108[k]
                   + f_0 * gl_423[k];

        t_244[k] = -2.0 * dl_109[k]
                   + f_0 * gl_424[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, dl_110, dl_111, dl_112, dl_113, \
                         dl_114, gl_425, gl_426, gl_427, gl_428, \
                         gl_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * dl_110[k]
                   + f_0 * gl_425[k];

        t_246[k] = -2.0 * dl_111[k]
                   + f_0 * gl_426[k];

        t_247[k] = -2.0 * dl_112[k]
                   + f_0 * gl_427[k];

        t_248[k] = -2.0 * dl_113[k]
                   + f_0 * gl_428[k];

        t_249[k] = -2.0 * dl_114[k]
                   + f_0 * gl_429[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, dl_115, dl_116, dl_117, dl_118, \
                         dl_119, gl_430, gl_431, gl_432, gl_433, \
                         gl_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * dl_115[k]
                   + f_0 * gl_430[k];

        t_251[k] = -2.0 * dl_116[k]
                   + f_0 * gl_431[k];

        t_252[k] = -2.0 * dl_117[k]
                   + f_0 * gl_432[k];

        t_253[k] = -2.0 * dl_118[k]
                   + f_0 * gl_433[k];

        t_254[k] = -2.0 * dl_119[k]
                   + f_0 * gl_434[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, dl_120, dl_121, dl_122, dl_123, \
                         dl_124, gl_435, gl_436, gl_437, gl_438, \
                         gl_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * dl_120[k]
                   + f_0 * gl_435[k];

        t_256[k] = -2.0 * dl_121[k]
                   + f_0 * gl_436[k];

        t_257[k] = -2.0 * dl_122[k]
                   + f_0 * gl_437[k];

        t_258[k] = -2.0 * dl_123[k]
                   + f_0 * gl_438[k];

        t_259[k] = -2.0 * dl_124[k]
                   + f_0 * gl_439[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, dl_125, dl_126, dl_127, dl_128, \
                         dl_129, gl_440, gl_441, gl_442, gl_443, \
                         gl_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * dl_125[k]
                   + f_0 * gl_440[k];

        t_261[k] = -2.0 * dl_126[k]
                   + f_0 * gl_441[k];

        t_262[k] = -2.0 * dl_127[k]
                   + f_0 * gl_442[k];

        t_263[k] = -2.0 * dl_128[k]
                   + f_0 * gl_443[k];

        t_264[k] = -2.0 * dl_129[k]
                   + f_0 * gl_444[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, dl_130, dl_131, dl_132, dl_133, \
                         dl_134, gl_445, gl_446, gl_447, gl_448, \
                         gl_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * dl_130[k]
                   + f_0 * gl_445[k];

        t_266[k] = -2.0 * dl_131[k]
                   + f_0 * gl_446[k];

        t_267[k] = -2.0 * dl_132[k]
                   + f_0 * gl_447[k];

        t_268[k] = -2.0 * dl_133[k]
                   + f_0 * gl_448[k];

        t_269[k] = -2.0 * dl_134[k]
                   + f_0 * gl_449[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, t_276, t_277, gl_495, \
                         gl_496, gl_497, gl_498, gl_499, gl_500, gl_501, \
                         gl_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * gl_495[k];

        t_271[k] = f_0 * gl_496[k];

        t_272[k] = f_0 * gl_497[k];

        t_273[k] = f_0 * gl_498[k];

        t_274[k] = f_0 * gl_499[k];

        t_275[k] = f_0 * gl_500[k];

        t_276[k] = f_0 * gl_501[k];

        t_277[k] = f_0 * gl_502[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, t_284, t_285, gl_503, \
                         gl_504, gl_505, gl_506, gl_507, gl_508, gl_509, \
                         gl_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * gl_503[k];

        t_279[k] = f_0 * gl_504[k];

        t_280[k] = f_0 * gl_505[k];

        t_281[k] = f_0 * gl_506[k];

        t_282[k] = f_0 * gl_507[k];

        t_283[k] = f_0 * gl_508[k];

        t_284[k] = f_0 * gl_509[k];

        t_285[k] = f_0 * gl_510[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, t_291, t_292, t_293, gl_511, \
                         gl_512, gl_513, gl_514, gl_515, gl_516, gl_517, \
                         gl_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * gl_511[k];

        t_287[k] = f_0 * gl_512[k];

        t_288[k] = f_0 * gl_513[k];

        t_289[k] = f_0 * gl_514[k];

        t_290[k] = f_0 * gl_515[k];

        t_291[k] = f_0 * gl_516[k];

        t_292[k] = f_0 * gl_517[k];

        t_293[k] = f_0 * gl_518[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, t_299, t_300, t_301, gl_519, \
                         gl_520, gl_521, gl_522, gl_523, gl_524, gl_525, \
                         gl_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_0 * gl_519[k];

        t_295[k] = f_0 * gl_520[k];

        t_296[k] = f_0 * gl_521[k];

        t_297[k] = f_0 * gl_522[k];

        t_298[k] = f_0 * gl_523[k];

        t_299[k] = f_0 * gl_524[k];

        t_300[k] = f_0 * gl_525[k];

        t_301[k] = f_0 * gl_526[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, t_307, t_308, t_309, gl_527, \
                         gl_528, gl_529, gl_530, gl_531, gl_532, gl_533, \
                         gl_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * gl_527[k];

        t_303[k] = f_0 * gl_528[k];

        t_304[k] = f_0 * gl_529[k];

        t_305[k] = f_0 * gl_530[k];

        t_306[k] = f_0 * gl_531[k];

        t_307[k] = f_0 * gl_532[k];

        t_308[k] = f_0 * gl_533[k];

        t_309[k] = f_0 * gl_534[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, t_315, t_316, dl_135, dl_136, \
                         gl_535, gl_536, gl_537, gl_538, gl_539, gl_540, \
                         gl_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_0 * gl_535[k];

        t_311[k] = f_0 * gl_536[k];

        t_312[k] = f_0 * gl_537[k];

        t_313[k] = f_0 * gl_538[k];

        t_314[k] = f_0 * gl_539[k];

        t_315[k] = -dl_135[k]
                   + f_0 * gl_540[k];

        t_316[k] = -dl_136[k]
                   + f_0 * gl_541[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, dl_137, dl_138, dl_139, dl_140, \
                         dl_141, gl_542, gl_543, gl_544, gl_545, \
                         gl_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -dl_137[k]
                   + f_0 * gl_542[k];

        t_318[k] = -dl_138[k]
                   + f_0 * gl_543[k];

        t_319[k] = -dl_139[k]
                   + f_0 * gl_544[k];

        t_320[k] = -dl_140[k]
                   + f_0 * gl_545[k];

        t_321[k] = -dl_141[k]
                   + f_0 * gl_546[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, dl_142, dl_143, dl_144, dl_145, \
                         dl_146, gl_547, gl_548, gl_549, gl_550, \
                         gl_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -dl_142[k]
                   + f_0 * gl_547[k];

        t_323[k] = -dl_143[k]
                   + f_0 * gl_548[k];

        t_324[k] = -dl_144[k]
                   + f_0 * gl_549[k];

        t_325[k] = -dl_145[k]
                   + f_0 * gl_550[k];

        t_326[k] = -dl_146[k]
                   + f_0 * gl_551[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, dl_147, dl_148, dl_149, dl_150, \
                         dl_151, gl_552, gl_553, gl_554, gl_555, \
                         gl_556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -dl_147[k]
                   + f_0 * gl_552[k];

        t_328[k] = -dl_148[k]
                   + f_0 * gl_553[k];

        t_329[k] = -dl_149[k]
                   + f_0 * gl_554[k];

        t_330[k] = -dl_150[k]
                   + f_0 * gl_555[k];

        t_331[k] = -dl_151[k]
                   + f_0 * gl_556[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, dl_152, dl_153, dl_154, dl_155, \
                         dl_156, gl_557, gl_558, gl_559, gl_560, \
                         gl_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -dl_152[k]
                   + f_0 * gl_557[k];

        t_333[k] = -dl_153[k]
                   + f_0 * gl_558[k];

        t_334[k] = -dl_154[k]
                   + f_0 * gl_559[k];

        t_335[k] = -dl_155[k]
                   + f_0 * gl_560[k];

        t_336[k] = -dl_156[k]
                   + f_0 * gl_561[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, dl_157, dl_158, dl_159, dl_160, \
                         dl_161, gl_562, gl_563, gl_564, gl_565, \
                         gl_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -dl_157[k]
                   + f_0 * gl_562[k];

        t_338[k] = -dl_158[k]
                   + f_0 * gl_563[k];

        t_339[k] = -dl_159[k]
                   + f_0 * gl_564[k];

        t_340[k] = -dl_160[k]
                   + f_0 * gl_565[k];

        t_341[k] = -dl_161[k]
                   + f_0 * gl_566[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, dl_162, dl_163, dl_164, dl_165, \
                         dl_166, gl_567, gl_568, gl_569, gl_570, \
                         gl_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -dl_162[k]
                   + f_0 * gl_567[k];

        t_343[k] = -dl_163[k]
                   + f_0 * gl_568[k];

        t_344[k] = -dl_164[k]
                   + f_0 * gl_569[k];

        t_345[k] = -dl_165[k]
                   + f_0 * gl_570[k];

        t_346[k] = -dl_166[k]
                   + f_0 * gl_571[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, dl_167, dl_168, dl_169, dl_170, \
                         dl_171, gl_572, gl_573, gl_574, gl_575, \
                         gl_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -dl_167[k]
                   + f_0 * gl_572[k];

        t_348[k] = -dl_168[k]
                   + f_0 * gl_573[k];

        t_349[k] = -dl_169[k]
                   + f_0 * gl_574[k];

        t_350[k] = -dl_170[k]
                   + f_0 * gl_575[k];

        t_351[k] = -dl_171[k]
                   + f_0 * gl_576[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, dl_172, dl_173, dl_174, dl_175, \
                         dl_176, gl_577, gl_578, gl_579, gl_580, \
                         gl_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -dl_172[k]
                   + f_0 * gl_577[k];

        t_353[k] = -dl_173[k]
                   + f_0 * gl_578[k];

        t_354[k] = -dl_174[k]
                   + f_0 * gl_579[k];

        t_355[k] = -dl_175[k]
                   + f_0 * gl_580[k];

        t_356[k] = -dl_176[k]
                   + f_0 * gl_581[k];
    }
}

static auto
compute_prim_geom_10_fl_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t dl, const size_t gl,
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

    const auto *dl_177 = buffer.data(dl + 177);
    const auto *dl_178 = buffer.data(dl + 178);
    const auto *dl_179 = buffer.data(dl + 179);
    const auto *dl_180 = buffer.data(dl + 180);
    const auto *dl_181 = buffer.data(dl + 181);
    const auto *dl_182 = buffer.data(dl + 182);
    const auto *dl_183 = buffer.data(dl + 183);
    const auto *dl_184 = buffer.data(dl + 184);
    const auto *dl_185 = buffer.data(dl + 185);
    const auto *dl_186 = buffer.data(dl + 186);
    const auto *dl_187 = buffer.data(dl + 187);
    const auto *dl_188 = buffer.data(dl + 188);
    const auto *dl_189 = buffer.data(dl + 189);
    const auto *dl_190 = buffer.data(dl + 190);
    const auto *dl_191 = buffer.data(dl + 191);
    const auto *dl_192 = buffer.data(dl + 192);
    const auto *dl_193 = buffer.data(dl + 193);
    const auto *dl_194 = buffer.data(dl + 194);
    const auto *dl_195 = buffer.data(dl + 195);
    const auto *dl_196 = buffer.data(dl + 196);
    const auto *dl_197 = buffer.data(dl + 197);
    const auto *dl_198 = buffer.data(dl + 198);
    const auto *dl_199 = buffer.data(dl + 199);
    const auto *dl_200 = buffer.data(dl + 200);
    const auto *dl_201 = buffer.data(dl + 201);
    const auto *dl_202 = buffer.data(dl + 202);
    const auto *dl_203 = buffer.data(dl + 203);
    const auto *dl_204 = buffer.data(dl + 204);
    const auto *dl_205 = buffer.data(dl + 205);
    const auto *dl_206 = buffer.data(dl + 206);
    const auto *dl_207 = buffer.data(dl + 207);
    const auto *dl_208 = buffer.data(dl + 208);
    const auto *dl_209 = buffer.data(dl + 209);
    const auto *dl_210 = buffer.data(dl + 210);
    const auto *dl_211 = buffer.data(dl + 211);
    const auto *dl_212 = buffer.data(dl + 212);
    const auto *dl_213 = buffer.data(dl + 213);
    const auto *dl_214 = buffer.data(dl + 214);
    const auto *dl_215 = buffer.data(dl + 215);
    const auto *dl_216 = buffer.data(dl + 216);
    const auto *dl_217 = buffer.data(dl + 217);
    const auto *dl_218 = buffer.data(dl + 218);
    const auto *dl_219 = buffer.data(dl + 219);
    const auto *dl_220 = buffer.data(dl + 220);
    const auto *dl_221 = buffer.data(dl + 221);
    const auto *dl_222 = buffer.data(dl + 222);
    const auto *dl_223 = buffer.data(dl + 223);
    const auto *dl_224 = buffer.data(dl + 224);
    const auto *dl_225 = buffer.data(dl + 225);
    const auto *dl_226 = buffer.data(dl + 226);
    const auto *dl_227 = buffer.data(dl + 227);
    const auto *dl_228 = buffer.data(dl + 228);
    const auto *dl_229 = buffer.data(dl + 229);
    const auto *dl_230 = buffer.data(dl + 230);
    const auto *dl_231 = buffer.data(dl + 231);
    const auto *dl_232 = buffer.data(dl + 232);
    const auto *dl_233 = buffer.data(dl + 233);
    const auto *dl_234 = buffer.data(dl + 234);
    const auto *dl_235 = buffer.data(dl + 235);
    const auto *dl_236 = buffer.data(dl + 236);
    const auto *dl_237 = buffer.data(dl + 237);
    const auto *dl_238 = buffer.data(dl + 238);
    const auto *dl_239 = buffer.data(dl + 239);
    const auto *dl_240 = buffer.data(dl + 240);
    const auto *dl_241 = buffer.data(dl + 241);
    const auto *dl_242 = buffer.data(dl + 242);
    const auto *dl_243 = buffer.data(dl + 243);
    const auto *dl_244 = buffer.data(dl + 244);
    const auto *dl_245 = buffer.data(dl + 245);
    const auto *dl_246 = buffer.data(dl + 246);
    const auto *dl_247 = buffer.data(dl + 247);
    const auto *dl_248 = buffer.data(dl + 248);
    const auto *dl_249 = buffer.data(dl + 249);
    const auto *dl_250 = buffer.data(dl + 250);
    const auto *dl_251 = buffer.data(dl + 251);
    const auto *dl_252 = buffer.data(dl + 252);
    const auto *dl_253 = buffer.data(dl + 253);
    const auto *dl_254 = buffer.data(dl + 254);
    const auto *dl_255 = buffer.data(dl + 255);
    const auto *dl_256 = buffer.data(dl + 256);
    const auto *dl_257 = buffer.data(dl + 257);
    const auto *dl_258 = buffer.data(dl + 258);
    const auto *dl_259 = buffer.data(dl + 259);
    const auto *dl_260 = buffer.data(dl + 260);
    const auto *dl_261 = buffer.data(dl + 261);
    const auto *dl_262 = buffer.data(dl + 262);
    const auto *dl_263 = buffer.data(dl + 263);
    const auto *dl_264 = buffer.data(dl + 264);
    const auto *dl_265 = buffer.data(dl + 265);
    const auto *dl_266 = buffer.data(dl + 266);
    const auto *dl_267 = buffer.data(dl + 267);
    const auto *dl_268 = buffer.data(dl + 268);
    const auto *dl_269 = buffer.data(dl + 269);

    const auto *gl_582 = buffer.data(gl + 582);
    const auto *gl_583 = buffer.data(gl + 583);
    const auto *gl_584 = buffer.data(gl + 584);
    const auto *gl_585 = buffer.data(gl + 585);
    const auto *gl_586 = buffer.data(gl + 586);
    const auto *gl_587 = buffer.data(gl + 587);
    const auto *gl_588 = buffer.data(gl + 588);
    const auto *gl_589 = buffer.data(gl + 589);
    const auto *gl_590 = buffer.data(gl + 590);
    const auto *gl_591 = buffer.data(gl + 591);
    const auto *gl_592 = buffer.data(gl + 592);
    const auto *gl_593 = buffer.data(gl + 593);
    const auto *gl_594 = buffer.data(gl + 594);
    const auto *gl_595 = buffer.data(gl + 595);
    const auto *gl_596 = buffer.data(gl + 596);
    const auto *gl_597 = buffer.data(gl + 597);
    const auto *gl_598 = buffer.data(gl + 598);
    const auto *gl_599 = buffer.data(gl + 599);
    const auto *gl_600 = buffer.data(gl + 600);
    const auto *gl_601 = buffer.data(gl + 601);
    const auto *gl_602 = buffer.data(gl + 602);
    const auto *gl_603 = buffer.data(gl + 603);
    const auto *gl_604 = buffer.data(gl + 604);
    const auto *gl_605 = buffer.data(gl + 605);
    const auto *gl_606 = buffer.data(gl + 606);
    const auto *gl_607 = buffer.data(gl + 607);
    const auto *gl_608 = buffer.data(gl + 608);
    const auto *gl_609 = buffer.data(gl + 609);
    const auto *gl_610 = buffer.data(gl + 610);
    const auto *gl_611 = buffer.data(gl + 611);
    const auto *gl_612 = buffer.data(gl + 612);
    const auto *gl_613 = buffer.data(gl + 613);
    const auto *gl_614 = buffer.data(gl + 614);
    const auto *gl_615 = buffer.data(gl + 615);
    const auto *gl_616 = buffer.data(gl + 616);
    const auto *gl_617 = buffer.data(gl + 617);
    const auto *gl_618 = buffer.data(gl + 618);
    const auto *gl_619 = buffer.data(gl + 619);
    const auto *gl_620 = buffer.data(gl + 620);
    const auto *gl_621 = buffer.data(gl + 621);
    const auto *gl_622 = buffer.data(gl + 622);
    const auto *gl_623 = buffer.data(gl + 623);
    const auto *gl_624 = buffer.data(gl + 624);
    const auto *gl_625 = buffer.data(gl + 625);
    const auto *gl_626 = buffer.data(gl + 626);
    const auto *gl_627 = buffer.data(gl + 627);
    const auto *gl_628 = buffer.data(gl + 628);
    const auto *gl_629 = buffer.data(gl + 629);
    const auto *gl_630 = buffer.data(gl + 630);
    const auto *gl_631 = buffer.data(gl + 631);
    const auto *gl_632 = buffer.data(gl + 632);
    const auto *gl_633 = buffer.data(gl + 633);
    const auto *gl_634 = buffer.data(gl + 634);
    const auto *gl_635 = buffer.data(gl + 635);
    const auto *gl_636 = buffer.data(gl + 636);
    const auto *gl_637 = buffer.data(gl + 637);
    const auto *gl_638 = buffer.data(gl + 638);
    const auto *gl_639 = buffer.data(gl + 639);
    const auto *gl_640 = buffer.data(gl + 640);
    const auto *gl_641 = buffer.data(gl + 641);
    const auto *gl_642 = buffer.data(gl + 642);
    const auto *gl_643 = buffer.data(gl + 643);
    const auto *gl_644 = buffer.data(gl + 644);
    const auto *gl_645 = buffer.data(gl + 645);
    const auto *gl_646 = buffer.data(gl + 646);
    const auto *gl_647 = buffer.data(gl + 647);
    const auto *gl_648 = buffer.data(gl + 648);
    const auto *gl_649 = buffer.data(gl + 649);
    const auto *gl_650 = buffer.data(gl + 650);
    const auto *gl_651 = buffer.data(gl + 651);
    const auto *gl_652 = buffer.data(gl + 652);
    const auto *gl_653 = buffer.data(gl + 653);
    const auto *gl_654 = buffer.data(gl + 654);
    const auto *gl_655 = buffer.data(gl + 655);
    const auto *gl_656 = buffer.data(gl + 656);
    const auto *gl_657 = buffer.data(gl + 657);
    const auto *gl_658 = buffer.data(gl + 658);
    const auto *gl_659 = buffer.data(gl + 659);
    const auto *gl_660 = buffer.data(gl + 660);
    const auto *gl_661 = buffer.data(gl + 661);
    const auto *gl_662 = buffer.data(gl + 662);
    const auto *gl_663 = buffer.data(gl + 663);
    const auto *gl_664 = buffer.data(gl + 664);
    const auto *gl_665 = buffer.data(gl + 665);
    const auto *gl_666 = buffer.data(gl + 666);
    const auto *gl_667 = buffer.data(gl + 667);
    const auto *gl_668 = buffer.data(gl + 668);
    const auto *gl_669 = buffer.data(gl + 669);
    const auto *gl_670 = buffer.data(gl + 670);
    const auto *gl_671 = buffer.data(gl + 671);
    const auto *gl_672 = buffer.data(gl + 672);
    const auto *gl_673 = buffer.data(gl + 673);
    const auto *gl_674 = buffer.data(gl + 674);

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, dl_177, dl_178, dl_179, dl_180, \
                         dl_181, gl_582, gl_583, gl_584, gl_585, \
                         gl_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -dl_177[k]
                   + f_0 * gl_582[k];

        t_358[k] = -dl_178[k]
                   + f_0 * gl_583[k];

        t_359[k] = -dl_179[k]
                   + f_0 * gl_584[k];

        t_360[k] = -2.0 * dl_180[k]
                   + f_0 * gl_585[k];

        t_361[k] = -2.0 * dl_181[k]
                   + f_0 * gl_586[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, dl_182, dl_183, dl_184, dl_185, \
                         dl_186, gl_587, gl_588, gl_589, gl_590, \
                         gl_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -2.0 * dl_182[k]
                   + f_0 * gl_587[k];

        t_363[k] = -2.0 * dl_183[k]
                   + f_0 * gl_588[k];

        t_364[k] = -2.0 * dl_184[k]
                   + f_0 * gl_589[k];

        t_365[k] = -2.0 * dl_185[k]
                   + f_0 * gl_590[k];

        t_366[k] = -2.0 * dl_186[k]
                   + f_0 * gl_591[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, dl_187, dl_188, dl_189, dl_190, \
                         dl_191, gl_592, gl_593, gl_594, gl_595, \
                         gl_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -2.0 * dl_187[k]
                   + f_0 * gl_592[k];

        t_368[k] = -2.0 * dl_188[k]
                   + f_0 * gl_593[k];

        t_369[k] = -2.0 * dl_189[k]
                   + f_0 * gl_594[k];

        t_370[k] = -2.0 * dl_190[k]
                   + f_0 * gl_595[k];

        t_371[k] = -2.0 * dl_191[k]
                   + f_0 * gl_596[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, dl_192, dl_193, dl_194, dl_195, \
                         dl_196, gl_597, gl_598, gl_599, gl_600, \
                         gl_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -2.0 * dl_192[k]
                   + f_0 * gl_597[k];

        t_373[k] = -2.0 * dl_193[k]
                   + f_0 * gl_598[k];

        t_374[k] = -2.0 * dl_194[k]
                   + f_0 * gl_599[k];

        t_375[k] = -2.0 * dl_195[k]
                   + f_0 * gl_600[k];

        t_376[k] = -2.0 * dl_196[k]
                   + f_0 * gl_601[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, dl_197, dl_198, dl_199, dl_200, \
                         dl_201, gl_602, gl_603, gl_604, gl_605, \
                         gl_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -2.0 * dl_197[k]
                   + f_0 * gl_602[k];

        t_378[k] = -2.0 * dl_198[k]
                   + f_0 * gl_603[k];

        t_379[k] = -2.0 * dl_199[k]
                   + f_0 * gl_604[k];

        t_380[k] = -2.0 * dl_200[k]
                   + f_0 * gl_605[k];

        t_381[k] = -2.0 * dl_201[k]
                   + f_0 * gl_606[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, dl_202, dl_203, dl_204, dl_205, \
                         dl_206, gl_607, gl_608, gl_609, gl_610, \
                         gl_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -2.0 * dl_202[k]
                   + f_0 * gl_607[k];

        t_383[k] = -2.0 * dl_203[k]
                   + f_0 * gl_608[k];

        t_384[k] = -2.0 * dl_204[k]
                   + f_0 * gl_609[k];

        t_385[k] = -2.0 * dl_205[k]
                   + f_0 * gl_610[k];

        t_386[k] = -2.0 * dl_206[k]
                   + f_0 * gl_611[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, dl_207, dl_208, dl_209, dl_210, \
                         dl_211, gl_612, gl_613, gl_614, gl_615, \
                         gl_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -2.0 * dl_207[k]
                   + f_0 * gl_612[k];

        t_388[k] = -2.0 * dl_208[k]
                   + f_0 * gl_613[k];

        t_389[k] = -2.0 * dl_209[k]
                   + f_0 * gl_614[k];

        t_390[k] = -2.0 * dl_210[k]
                   + f_0 * gl_615[k];

        t_391[k] = -2.0 * dl_211[k]
                   + f_0 * gl_616[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, dl_212, dl_213, dl_214, dl_215, \
                         dl_216, gl_617, gl_618, gl_619, gl_620, \
                         gl_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -2.0 * dl_212[k]
                   + f_0 * gl_617[k];

        t_393[k] = -2.0 * dl_213[k]
                   + f_0 * gl_618[k];

        t_394[k] = -2.0 * dl_214[k]
                   + f_0 * gl_619[k];

        t_395[k] = -2.0 * dl_215[k]
                   + f_0 * gl_620[k];

        t_396[k] = -2.0 * dl_216[k]
                   + f_0 * gl_621[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, dl_217, dl_218, dl_219, dl_220, \
                         dl_221, gl_622, gl_623, gl_624, gl_625, \
                         gl_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -2.0 * dl_217[k]
                   + f_0 * gl_622[k];

        t_398[k] = -2.0 * dl_218[k]
                   + f_0 * gl_623[k];

        t_399[k] = -2.0 * dl_219[k]
                   + f_0 * gl_624[k];

        t_400[k] = -2.0 * dl_220[k]
                   + f_0 * gl_625[k];

        t_401[k] = -2.0 * dl_221[k]
                   + f_0 * gl_626[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, dl_222, dl_223, dl_224, dl_225, \
                         dl_226, gl_627, gl_628, gl_629, gl_630, \
                         gl_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -2.0 * dl_222[k]
                   + f_0 * gl_627[k];

        t_403[k] = -2.0 * dl_223[k]
                   + f_0 * gl_628[k];

        t_404[k] = -2.0 * dl_224[k]
                   + f_0 * gl_629[k];

        t_405[k] = -3.0 * dl_225[k]
                   + f_0 * gl_630[k];

        t_406[k] = -3.0 * dl_226[k]
                   + f_0 * gl_631[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, dl_227, dl_228, dl_229, dl_230, \
                         dl_231, gl_632, gl_633, gl_634, gl_635, \
                         gl_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -3.0 * dl_227[k]
                   + f_0 * gl_632[k];

        t_408[k] = -3.0 * dl_228[k]
                   + f_0 * gl_633[k];

        t_409[k] = -3.0 * dl_229[k]
                   + f_0 * gl_634[k];

        t_410[k] = -3.0 * dl_230[k]
                   + f_0 * gl_635[k];

        t_411[k] = -3.0 * dl_231[k]
                   + f_0 * gl_636[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, dl_232, dl_233, dl_234, dl_235, \
                         dl_236, gl_637, gl_638, gl_639, gl_640, \
                         gl_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -3.0 * dl_232[k]
                   + f_0 * gl_637[k];

        t_413[k] = -3.0 * dl_233[k]
                   + f_0 * gl_638[k];

        t_414[k] = -3.0 * dl_234[k]
                   + f_0 * gl_639[k];

        t_415[k] = -3.0 * dl_235[k]
                   + f_0 * gl_640[k];

        t_416[k] = -3.0 * dl_236[k]
                   + f_0 * gl_641[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, dl_237, dl_238, dl_239, dl_240, \
                         dl_241, gl_642, gl_643, gl_644, gl_645, \
                         gl_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -3.0 * dl_237[k]
                   + f_0 * gl_642[k];

        t_418[k] = -3.0 * dl_238[k]
                   + f_0 * gl_643[k];

        t_419[k] = -3.0 * dl_239[k]
                   + f_0 * gl_644[k];

        t_420[k] = -3.0 * dl_240[k]
                   + f_0 * gl_645[k];

        t_421[k] = -3.0 * dl_241[k]
                   + f_0 * gl_646[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, dl_242, dl_243, dl_244, dl_245, \
                         dl_246, gl_647, gl_648, gl_649, gl_650, \
                         gl_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = -3.0 * dl_242[k]
                   + f_0 * gl_647[k];

        t_423[k] = -3.0 * dl_243[k]
                   + f_0 * gl_648[k];

        t_424[k] = -3.0 * dl_244[k]
                   + f_0 * gl_649[k];

        t_425[k] = -3.0 * dl_245[k]
                   + f_0 * gl_650[k];

        t_426[k] = -3.0 * dl_246[k]
                   + f_0 * gl_651[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, t_431, dl_247, dl_248, dl_249, dl_250, \
                         dl_251, gl_652, gl_653, gl_654, gl_655, \
                         gl_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = -3.0 * dl_247[k]
                   + f_0 * gl_652[k];

        t_428[k] = -3.0 * dl_248[k]
                   + f_0 * gl_653[k];

        t_429[k] = -3.0 * dl_249[k]
                   + f_0 * gl_654[k];

        t_430[k] = -3.0 * dl_250[k]
                   + f_0 * gl_655[k];

        t_431[k] = -3.0 * dl_251[k]
                   + f_0 * gl_656[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, t_436, dl_252, dl_253, dl_254, dl_255, \
                         dl_256, gl_657, gl_658, gl_659, gl_660, \
                         gl_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = -3.0 * dl_252[k]
                   + f_0 * gl_657[k];

        t_433[k] = -3.0 * dl_253[k]
                   + f_0 * gl_658[k];

        t_434[k] = -3.0 * dl_254[k]
                   + f_0 * gl_659[k];

        t_435[k] = -3.0 * dl_255[k]
                   + f_0 * gl_660[k];

        t_436[k] = -3.0 * dl_256[k]
                   + f_0 * gl_661[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, t_441, dl_257, dl_258, dl_259, dl_260, \
                         dl_261, gl_662, gl_663, gl_664, gl_665, \
                         gl_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = -3.0 * dl_257[k]
                   + f_0 * gl_662[k];

        t_438[k] = -3.0 * dl_258[k]
                   + f_0 * gl_663[k];

        t_439[k] = -3.0 * dl_259[k]
                   + f_0 * gl_664[k];

        t_440[k] = -3.0 * dl_260[k]
                   + f_0 * gl_665[k];

        t_441[k] = -3.0 * dl_261[k]
                   + f_0 * gl_666[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, t_446, dl_262, dl_263, dl_264, dl_265, \
                         dl_266, gl_667, gl_668, gl_669, gl_670, \
                         gl_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = -3.0 * dl_262[k]
                   + f_0 * gl_667[k];

        t_443[k] = -3.0 * dl_263[k]
                   + f_0 * gl_668[k];

        t_444[k] = -3.0 * dl_264[k]
                   + f_0 * gl_669[k];

        t_445[k] = -3.0 * dl_265[k]
                   + f_0 * gl_670[k];

        t_446[k] = -3.0 * dl_266[k]
                   + f_0 * gl_671[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, dl_267, dl_268, dl_269, gl_672, gl_673, \
                         gl_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = -3.0 * dl_267[k]
                   + f_0 * gl_672[k];

        t_448[k] = -3.0 * dl_268[k]
                   + f_0 * gl_673[k];

        t_449[k] = -3.0 * dl_269[k]
                   + f_0 * gl_674[k];
    }
}

auto
compute_prim_geom_10_fl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t dl, const size_t gl,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fl_electron_repulsion_2_piece0(buffer, target, dl, gl, ncols, alpha);

    compute_prim_geom_10_fl_electron_repulsion_2_piece1(buffer, target, dl, gl, ncols, alpha);

    compute_prim_geom_10_fl_electron_repulsion_2_piece2(buffer, target, dl, gl, ncols, alpha);
}

}  // namespace simdt2ceri
