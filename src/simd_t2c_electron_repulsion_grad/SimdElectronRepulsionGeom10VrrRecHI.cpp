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


#include "SimdElectronRepulsionGeom10VrrRecHI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_hi_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_1 = buffer.data(gi + 1);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_41 = buffer.data(gi + 41);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_44 = buffer.data(gi + 44);
    const auto *gi_45 = buffer.data(gi + 45);
    const auto *gi_46 = buffer.data(gi + 46);
    const auto *gi_47 = buffer.data(gi + 47);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_50 = buffer.data(gi + 50);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_60 = buffer.data(gi + 60);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_63 = buffer.data(gi + 63);
    const auto *gi_64 = buffer.data(gi + 64);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_67 = buffer.data(gi + 67);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_69 = buffer.data(gi + 69);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_71 = buffer.data(gi + 71);
    const auto *gi_72 = buffer.data(gi + 72);
    const auto *gi_73 = buffer.data(gi + 73);
    const auto *gi_74 = buffer.data(gi + 74);
    const auto *gi_75 = buffer.data(gi + 75);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_77 = buffer.data(gi + 77);
    const auto *gi_78 = buffer.data(gi + 78);
    const auto *gi_79 = buffer.data(gi + 79);
    const auto *gi_80 = buffer.data(gi + 80);
    const auto *gi_81 = buffer.data(gi + 81);
    const auto *gi_82 = buffer.data(gi + 82);
    const auto *gi_83 = buffer.data(gi + 83);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_85 = buffer.data(gi + 85);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_88 = buffer.data(gi + 88);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_90 = buffer.data(gi + 90);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_92 = buffer.data(gi + 92);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_95 = buffer.data(gi + 95);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_97 = buffer.data(gi + 97);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_100 = buffer.data(gi + 100);
    const auto *gi_101 = buffer.data(gi + 101);
    const auto *gi_102 = buffer.data(gi + 102);
    const auto *gi_103 = buffer.data(gi + 103);
    const auto *gi_104 = buffer.data(gi + 104);
    const auto *gi_105 = buffer.data(gi + 105);
    const auto *gi_106 = buffer.data(gi + 106);
    const auto *gi_107 = buffer.data(gi + 107);
    const auto *gi_108 = buffer.data(gi + 108);
    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_110 = buffer.data(gi + 110);
    const auto *gi_111 = buffer.data(gi + 111);
    const auto *gi_112 = buffer.data(gi + 112);
    const auto *gi_113 = buffer.data(gi + 113);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_116 = buffer.data(gi + 116);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_119 = buffer.data(gi + 119);
    const auto *gi_120 = buffer.data(gi + 120);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_123 = buffer.data(gi + 123);
    const auto *gi_124 = buffer.data(gi + 124);
    const auto *gi_125 = buffer.data(gi + 125);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_127 = buffer.data(gi + 127);
    const auto *gi_128 = buffer.data(gi + 128);
    const auto *gi_129 = buffer.data(gi + 129);
    const auto *gi_130 = buffer.data(gi + 130);
    const auto *gi_131 = buffer.data(gi + 131);
    const auto *gi_132 = buffer.data(gi + 132);
    const auto *gi_133 = buffer.data(gi + 133);
    const auto *gi_134 = buffer.data(gi + 134);
    const auto *gi_135 = buffer.data(gi + 135);
    const auto *gi_136 = buffer.data(gi + 136);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_138 = buffer.data(gi + 138);
    const auto *gi_139 = buffer.data(gi + 139);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_141 = buffer.data(gi + 141);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_144 = buffer.data(gi + 144);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_147 = buffer.data(gi + 147);
    const auto *gi_148 = buffer.data(gi + 148);
    const auto *gi_149 = buffer.data(gi + 149);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, gi_0, gi_1, gi_2, gi_3, gi_4, ii_0, ii_1, \
                         ii_2, ii_3, ii_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -5.0 * gi_0[k]
                 + f_0 * ii_0[k];

        t_1[k] = -5.0 * gi_1[k]
                 + f_0 * ii_1[k];

        t_2[k] = -5.0 * gi_2[k]
                 + f_0 * ii_2[k];

        t_3[k] = -5.0 * gi_3[k]
                 + f_0 * ii_3[k];

        t_4[k] = -5.0 * gi_4[k]
                 + f_0 * ii_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, gi_5, gi_6, gi_7, gi_8, gi_9, ii_5, ii_6, \
                         ii_7, ii_8, ii_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -5.0 * gi_5[k]
                 + f_0 * ii_5[k];

        t_6[k] = -5.0 * gi_6[k]
                 + f_0 * ii_6[k];

        t_7[k] = -5.0 * gi_7[k]
                 + f_0 * ii_7[k];

        t_8[k] = -5.0 * gi_8[k]
                 + f_0 * ii_8[k];

        t_9[k] = -5.0 * gi_9[k]
                 + f_0 * ii_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, gi_10, gi_11, gi_12, gi_13, gi_14, \
                         ii_10, ii_11, ii_12, ii_13, ii_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -5.0 * gi_10[k]
                  + f_0 * ii_10[k];

        t_11[k] = -5.0 * gi_11[k]
                  + f_0 * ii_11[k];

        t_12[k] = -5.0 * gi_12[k]
                  + f_0 * ii_12[k];

        t_13[k] = -5.0 * gi_13[k]
                  + f_0 * ii_13[k];

        t_14[k] = -5.0 * gi_14[k]
                  + f_0 * ii_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, gi_15, gi_16, gi_17, gi_18, gi_19, \
                         ii_15, ii_16, ii_17, ii_18, ii_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -5.0 * gi_15[k]
                  + f_0 * ii_15[k];

        t_16[k] = -5.0 * gi_16[k]
                  + f_0 * ii_16[k];

        t_17[k] = -5.0 * gi_17[k]
                  + f_0 * ii_17[k];

        t_18[k] = -5.0 * gi_18[k]
                  + f_0 * ii_18[k];

        t_19[k] = -5.0 * gi_19[k]
                  + f_0 * ii_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, gi_20, gi_21, gi_22, gi_23, gi_24, \
                         ii_20, ii_21, ii_22, ii_23, ii_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -5.0 * gi_20[k]
                  + f_0 * ii_20[k];

        t_21[k] = -5.0 * gi_21[k]
                  + f_0 * ii_21[k];

        t_22[k] = -5.0 * gi_22[k]
                  + f_0 * ii_22[k];

        t_23[k] = -5.0 * gi_23[k]
                  + f_0 * ii_23[k];

        t_24[k] = -5.0 * gi_24[k]
                  + f_0 * ii_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, gi_25, gi_26, gi_27, gi_28, gi_29, \
                         ii_25, ii_26, ii_27, ii_28, ii_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -5.0 * gi_25[k]
                  + f_0 * ii_25[k];

        t_26[k] = -5.0 * gi_26[k]
                  + f_0 * ii_26[k];

        t_27[k] = -5.0 * gi_27[k]
                  + f_0 * ii_27[k];

        t_28[k] = -4.0 * gi_28[k]
                  + f_0 * ii_28[k];

        t_29[k] = -4.0 * gi_29[k]
                  + f_0 * ii_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, gi_30, gi_31, gi_32, gi_33, gi_34, \
                         ii_30, ii_31, ii_32, ii_33, ii_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -4.0 * gi_30[k]
                  + f_0 * ii_30[k];

        t_31[k] = -4.0 * gi_31[k]
                  + f_0 * ii_31[k];

        t_32[k] = -4.0 * gi_32[k]
                  + f_0 * ii_32[k];

        t_33[k] = -4.0 * gi_33[k]
                  + f_0 * ii_33[k];

        t_34[k] = -4.0 * gi_34[k]
                  + f_0 * ii_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gi_35, gi_36, gi_37, gi_38, gi_39, \
                         ii_35, ii_36, ii_37, ii_38, ii_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -4.0 * gi_35[k]
                  + f_0 * ii_35[k];

        t_36[k] = -4.0 * gi_36[k]
                  + f_0 * ii_36[k];

        t_37[k] = -4.0 * gi_37[k]
                  + f_0 * ii_37[k];

        t_38[k] = -4.0 * gi_38[k]
                  + f_0 * ii_38[k];

        t_39[k] = -4.0 * gi_39[k]
                  + f_0 * ii_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gi_40, gi_41, gi_42, gi_43, gi_44, \
                         ii_40, ii_41, ii_42, ii_43, ii_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -4.0 * gi_40[k]
                  + f_0 * ii_40[k];

        t_41[k] = -4.0 * gi_41[k]
                  + f_0 * ii_41[k];

        t_42[k] = -4.0 * gi_42[k]
                  + f_0 * ii_42[k];

        t_43[k] = -4.0 * gi_43[k]
                  + f_0 * ii_43[k];

        t_44[k] = -4.0 * gi_44[k]
                  + f_0 * ii_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gi_45, gi_46, gi_47, gi_48, gi_49, \
                         ii_45, ii_46, ii_47, ii_48, ii_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -4.0 * gi_45[k]
                  + f_0 * ii_45[k];

        t_46[k] = -4.0 * gi_46[k]
                  + f_0 * ii_46[k];

        t_47[k] = -4.0 * gi_47[k]
                  + f_0 * ii_47[k];

        t_48[k] = -4.0 * gi_48[k]
                  + f_0 * ii_48[k];

        t_49[k] = -4.0 * gi_49[k]
                  + f_0 * ii_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gi_50, gi_51, gi_52, gi_53, gi_54, \
                         ii_50, ii_51, ii_52, ii_53, ii_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * gi_50[k]
                  + f_0 * ii_50[k];

        t_51[k] = -4.0 * gi_51[k]
                  + f_0 * ii_51[k];

        t_52[k] = -4.0 * gi_52[k]
                  + f_0 * ii_52[k];

        t_53[k] = -4.0 * gi_53[k]
                  + f_0 * ii_53[k];

        t_54[k] = -4.0 * gi_54[k]
                  + f_0 * ii_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gi_55, gi_56, gi_57, gi_58, gi_59, \
                         ii_55, ii_56, ii_57, ii_58, ii_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -4.0 * gi_55[k]
                  + f_0 * ii_55[k];

        t_56[k] = -4.0 * gi_56[k]
                  + f_0 * ii_56[k];

        t_57[k] = -4.0 * gi_57[k]
                  + f_0 * ii_57[k];

        t_58[k] = -4.0 * gi_58[k]
                  + f_0 * ii_58[k];

        t_59[k] = -4.0 * gi_59[k]
                  + f_0 * ii_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, gi_60, gi_61, gi_62, gi_63, gi_64, \
                         ii_60, ii_61, ii_62, ii_63, ii_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -4.0 * gi_60[k]
                  + f_0 * ii_60[k];

        t_61[k] = -4.0 * gi_61[k]
                  + f_0 * ii_61[k];

        t_62[k] = -4.0 * gi_62[k]
                  + f_0 * ii_62[k];

        t_63[k] = -4.0 * gi_63[k]
                  + f_0 * ii_63[k];

        t_64[k] = -4.0 * gi_64[k]
                  + f_0 * ii_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gi_65, gi_66, gi_67, gi_68, gi_69, \
                         ii_65, ii_66, ii_67, ii_68, ii_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -4.0 * gi_65[k]
                  + f_0 * ii_65[k];

        t_66[k] = -4.0 * gi_66[k]
                  + f_0 * ii_66[k];

        t_67[k] = -4.0 * gi_67[k]
                  + f_0 * ii_67[k];

        t_68[k] = -4.0 * gi_68[k]
                  + f_0 * ii_68[k];

        t_69[k] = -4.0 * gi_69[k]
                  + f_0 * ii_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gi_70, gi_71, gi_72, gi_73, gi_74, \
                         ii_70, ii_71, ii_72, ii_73, ii_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -4.0 * gi_70[k]
                  + f_0 * ii_70[k];

        t_71[k] = -4.0 * gi_71[k]
                  + f_0 * ii_71[k];

        t_72[k] = -4.0 * gi_72[k]
                  + f_0 * ii_72[k];

        t_73[k] = -4.0 * gi_73[k]
                  + f_0 * ii_73[k];

        t_74[k] = -4.0 * gi_74[k]
                  + f_0 * ii_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, gi_75, gi_76, gi_77, gi_78, gi_79, \
                         ii_75, ii_76, ii_77, ii_78, ii_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -4.0 * gi_75[k]
                  + f_0 * ii_75[k];

        t_76[k] = -4.0 * gi_76[k]
                  + f_0 * ii_76[k];

        t_77[k] = -4.0 * gi_77[k]
                  + f_0 * ii_77[k];

        t_78[k] = -4.0 * gi_78[k]
                  + f_0 * ii_78[k];

        t_79[k] = -4.0 * gi_79[k]
                  + f_0 * ii_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, gi_80, gi_81, gi_82, gi_83, gi_84, \
                         ii_80, ii_81, ii_82, ii_83, ii_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -4.0 * gi_80[k]
                  + f_0 * ii_80[k];

        t_81[k] = -4.0 * gi_81[k]
                  + f_0 * ii_81[k];

        t_82[k] = -4.0 * gi_82[k]
                  + f_0 * ii_82[k];

        t_83[k] = -4.0 * gi_83[k]
                  + f_0 * ii_83[k];

        t_84[k] = -3.0 * gi_84[k]
                  + f_0 * ii_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gi_85, gi_86, gi_87, gi_88, gi_89, \
                         ii_85, ii_86, ii_87, ii_88, ii_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -3.0 * gi_85[k]
                  + f_0 * ii_85[k];

        t_86[k] = -3.0 * gi_86[k]
                  + f_0 * ii_86[k];

        t_87[k] = -3.0 * gi_87[k]
                  + f_0 * ii_87[k];

        t_88[k] = -3.0 * gi_88[k]
                  + f_0 * ii_88[k];

        t_89[k] = -3.0 * gi_89[k]
                  + f_0 * ii_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, gi_90, gi_91, gi_92, gi_93, gi_94, \
                         ii_90, ii_91, ii_92, ii_93, ii_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * gi_90[k]
                  + f_0 * ii_90[k];

        t_91[k] = -3.0 * gi_91[k]
                  + f_0 * ii_91[k];

        t_92[k] = -3.0 * gi_92[k]
                  + f_0 * ii_92[k];

        t_93[k] = -3.0 * gi_93[k]
                  + f_0 * ii_93[k];

        t_94[k] = -3.0 * gi_94[k]
                  + f_0 * ii_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, gi_95, gi_96, gi_97, gi_98, gi_99, \
                         ii_95, ii_96, ii_97, ii_98, ii_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * gi_95[k]
                  + f_0 * ii_95[k];

        t_96[k] = -3.0 * gi_96[k]
                  + f_0 * ii_96[k];

        t_97[k] = -3.0 * gi_97[k]
                  + f_0 * ii_97[k];

        t_98[k] = -3.0 * gi_98[k]
                  + f_0 * ii_98[k];

        t_99[k] = -3.0 * gi_99[k]
                  + f_0 * ii_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, gi_100, gi_101, gi_102, gi_103, \
                         gi_104, ii_100, ii_101, ii_102, ii_103, \
                         ii_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * gi_100[k]
                   + f_0 * ii_100[k];

        t_101[k] = -3.0 * gi_101[k]
                   + f_0 * ii_101[k];

        t_102[k] = -3.0 * gi_102[k]
                   + f_0 * ii_102[k];

        t_103[k] = -3.0 * gi_103[k]
                   + f_0 * ii_103[k];

        t_104[k] = -3.0 * gi_104[k]
                   + f_0 * ii_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, gi_105, gi_106, gi_107, gi_108, \
                         gi_109, ii_105, ii_106, ii_107, ii_108, \
                         ii_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * gi_105[k]
                   + f_0 * ii_105[k];

        t_106[k] = -3.0 * gi_106[k]
                   + f_0 * ii_106[k];

        t_107[k] = -3.0 * gi_107[k]
                   + f_0 * ii_107[k];

        t_108[k] = -3.0 * gi_108[k]
                   + f_0 * ii_108[k];

        t_109[k] = -3.0 * gi_109[k]
                   + f_0 * ii_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, gi_110, gi_111, gi_112, gi_113, \
                         gi_114, ii_110, ii_111, ii_112, ii_113, \
                         ii_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -3.0 * gi_110[k]
                   + f_0 * ii_110[k];

        t_111[k] = -3.0 * gi_111[k]
                   + f_0 * ii_111[k];

        t_112[k] = -3.0 * gi_112[k]
                   + f_0 * ii_112[k];

        t_113[k] = -3.0 * gi_113[k]
                   + f_0 * ii_113[k];

        t_114[k] = -3.0 * gi_114[k]
                   + f_0 * ii_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, gi_115, gi_116, gi_117, gi_118, \
                         gi_119, ii_115, ii_116, ii_117, ii_118, \
                         ii_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -3.0 * gi_115[k]
                   + f_0 * ii_115[k];

        t_116[k] = -3.0 * gi_116[k]
                   + f_0 * ii_116[k];

        t_117[k] = -3.0 * gi_117[k]
                   + f_0 * ii_117[k];

        t_118[k] = -3.0 * gi_118[k]
                   + f_0 * ii_118[k];

        t_119[k] = -3.0 * gi_119[k]
                   + f_0 * ii_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, gi_120, gi_121, gi_122, gi_123, \
                         gi_124, ii_120, ii_121, ii_122, ii_123, \
                         ii_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -3.0 * gi_120[k]
                   + f_0 * ii_120[k];

        t_121[k] = -3.0 * gi_121[k]
                   + f_0 * ii_121[k];

        t_122[k] = -3.0 * gi_122[k]
                   + f_0 * ii_122[k];

        t_123[k] = -3.0 * gi_123[k]
                   + f_0 * ii_123[k];

        t_124[k] = -3.0 * gi_124[k]
                   + f_0 * ii_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, gi_125, gi_126, gi_127, gi_128, \
                         gi_129, ii_125, ii_126, ii_127, ii_128, \
                         ii_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -3.0 * gi_125[k]
                   + f_0 * ii_125[k];

        t_126[k] = -3.0 * gi_126[k]
                   + f_0 * ii_126[k];

        t_127[k] = -3.0 * gi_127[k]
                   + f_0 * ii_127[k];

        t_128[k] = -3.0 * gi_128[k]
                   + f_0 * ii_128[k];

        t_129[k] = -3.0 * gi_129[k]
                   + f_0 * ii_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, gi_130, gi_131, gi_132, gi_133, \
                         gi_134, ii_130, ii_131, ii_132, ii_133, \
                         ii_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -3.0 * gi_130[k]
                   + f_0 * ii_130[k];

        t_131[k] = -3.0 * gi_131[k]
                   + f_0 * ii_131[k];

        t_132[k] = -3.0 * gi_132[k]
                   + f_0 * ii_132[k];

        t_133[k] = -3.0 * gi_133[k]
                   + f_0 * ii_133[k];

        t_134[k] = -3.0 * gi_134[k]
                   + f_0 * ii_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, gi_135, gi_136, gi_137, gi_138, \
                         gi_139, ii_135, ii_136, ii_137, ii_138, \
                         ii_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -3.0 * gi_135[k]
                   + f_0 * ii_135[k];

        t_136[k] = -3.0 * gi_136[k]
                   + f_0 * ii_136[k];

        t_137[k] = -3.0 * gi_137[k]
                   + f_0 * ii_137[k];

        t_138[k] = -3.0 * gi_138[k]
                   + f_0 * ii_138[k];

        t_139[k] = -3.0 * gi_139[k]
                   + f_0 * ii_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, gi_140, gi_141, gi_142, gi_143, \
                         gi_144, ii_140, ii_141, ii_142, ii_143, \
                         ii_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -3.0 * gi_140[k]
                   + f_0 * ii_140[k];

        t_141[k] = -3.0 * gi_141[k]
                   + f_0 * ii_141[k];

        t_142[k] = -3.0 * gi_142[k]
                   + f_0 * ii_142[k];

        t_143[k] = -3.0 * gi_143[k]
                   + f_0 * ii_143[k];

        t_144[k] = -3.0 * gi_144[k]
                   + f_0 * ii_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, gi_145, gi_146, gi_147, gi_148, \
                         gi_149, ii_145, ii_146, ii_147, ii_148, \
                         ii_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -3.0 * gi_145[k]
                   + f_0 * ii_145[k];

        t_146[k] = -3.0 * gi_146[k]
                   + f_0 * ii_146[k];

        t_147[k] = -3.0 * gi_147[k]
                   + f_0 * ii_147[k];

        t_148[k] = -3.0 * gi_148[k]
                   + f_0 * ii_148[k];

        t_149[k] = -3.0 * gi_149[k]
                   + f_0 * ii_149[k];
    }
}

static auto
compute_prim_geom_10_hi_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_151 = buffer.data(gi + 151);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);
    const auto *gi_158 = buffer.data(gi + 158);
    const auto *gi_159 = buffer.data(gi + 159);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_161 = buffer.data(gi + 161);
    const auto *gi_162 = buffer.data(gi + 162);
    const auto *gi_163 = buffer.data(gi + 163);
    const auto *gi_164 = buffer.data(gi + 164);
    const auto *gi_165 = buffer.data(gi + 165);
    const auto *gi_166 = buffer.data(gi + 166);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_170 = buffer.data(gi + 170);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_172 = buffer.data(gi + 172);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_175 = buffer.data(gi + 175);
    const auto *gi_176 = buffer.data(gi + 176);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_179 = buffer.data(gi + 179);
    const auto *gi_180 = buffer.data(gi + 180);
    const auto *gi_181 = buffer.data(gi + 181);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_184 = buffer.data(gi + 184);
    const auto *gi_185 = buffer.data(gi + 185);
    const auto *gi_186 = buffer.data(gi + 186);
    const auto *gi_187 = buffer.data(gi + 187);
    const auto *gi_188 = buffer.data(gi + 188);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_190 = buffer.data(gi + 190);
    const auto *gi_191 = buffer.data(gi + 191);
    const auto *gi_192 = buffer.data(gi + 192);
    const auto *gi_193 = buffer.data(gi + 193);
    const auto *gi_194 = buffer.data(gi + 194);
    const auto *gi_195 = buffer.data(gi + 195);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_197 = buffer.data(gi + 197);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_200 = buffer.data(gi + 200);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_203 = buffer.data(gi + 203);
    const auto *gi_204 = buffer.data(gi + 204);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_207 = buffer.data(gi + 207);
    const auto *gi_208 = buffer.data(gi + 208);
    const auto *gi_209 = buffer.data(gi + 209);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_211 = buffer.data(gi + 211);
    const auto *gi_212 = buffer.data(gi + 212);
    const auto *gi_213 = buffer.data(gi + 213);
    const auto *gi_214 = buffer.data(gi + 214);
    const auto *gi_215 = buffer.data(gi + 215);
    const auto *gi_216 = buffer.data(gi + 216);
    const auto *gi_217 = buffer.data(gi + 217);
    const auto *gi_218 = buffer.data(gi + 218);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_222 = buffer.data(gi + 222);
    const auto *gi_223 = buffer.data(gi + 223);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_225 = buffer.data(gi + 225);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_228 = buffer.data(gi + 228);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_231 = buffer.data(gi + 231);
    const auto *gi_232 = buffer.data(gi + 232);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_235 = buffer.data(gi + 235);
    const auto *gi_236 = buffer.data(gi + 236);
    const auto *gi_237 = buffer.data(gi + 237);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_239 = buffer.data(gi + 239);
    const auto *gi_240 = buffer.data(gi + 240);
    const auto *gi_241 = buffer.data(gi + 241);
    const auto *gi_242 = buffer.data(gi + 242);
    const auto *gi_243 = buffer.data(gi + 243);
    const auto *gi_244 = buffer.data(gi + 244);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_246 = buffer.data(gi + 246);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_250 = buffer.data(gi + 250);
    const auto *gi_251 = buffer.data(gi + 251);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_253 = buffer.data(gi + 253);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_256 = buffer.data(gi + 256);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_258 = buffer.data(gi + 258);
    const auto *gi_259 = buffer.data(gi + 259);
    const auto *gi_260 = buffer.data(gi + 260);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_263 = buffer.data(gi + 263);
    const auto *gi_264 = buffer.data(gi + 264);
    const auto *gi_265 = buffer.data(gi + 265);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_267 = buffer.data(gi + 267);
    const auto *gi_268 = buffer.data(gi + 268);
    const auto *gi_269 = buffer.data(gi + 269);
    const auto *gi_270 = buffer.data(gi + 270);
    const auto *gi_271 = buffer.data(gi + 271);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_273 = buffer.data(gi + 273);
    const auto *gi_274 = buffer.data(gi + 274);
    const auto *gi_275 = buffer.data(gi + 275);
    const auto *gi_276 = buffer.data(gi + 276);
    const auto *gi_277 = buffer.data(gi + 277);
    const auto *gi_278 = buffer.data(gi + 278);
    const auto *gi_279 = buffer.data(gi + 279);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_284 = buffer.data(gi + 284);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_287 = buffer.data(gi + 287);
    const auto *gi_288 = buffer.data(gi + 288);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_291 = buffer.data(gi + 291);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_293 = buffer.data(gi + 293);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_296 = buffer.data(gi + 296);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_299 = buffer.data(gi + 299);

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

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, gi_150, gi_151, gi_152, gi_153, \
                         gi_154, ii_150, ii_151, ii_152, ii_153, \
                         ii_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -3.0 * gi_150[k]
                   + f_0 * ii_150[k];

        t_151[k] = -3.0 * gi_151[k]
                   + f_0 * ii_151[k];

        t_152[k] = -3.0 * gi_152[k]
                   + f_0 * ii_152[k];

        t_153[k] = -3.0 * gi_153[k]
                   + f_0 * ii_153[k];

        t_154[k] = -3.0 * gi_154[k]
                   + f_0 * ii_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, gi_155, gi_156, gi_157, gi_158, \
                         gi_159, ii_155, ii_156, ii_157, ii_158, \
                         ii_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -3.0 * gi_155[k]
                   + f_0 * ii_155[k];

        t_156[k] = -3.0 * gi_156[k]
                   + f_0 * ii_156[k];

        t_157[k] = -3.0 * gi_157[k]
                   + f_0 * ii_157[k];

        t_158[k] = -3.0 * gi_158[k]
                   + f_0 * ii_158[k];

        t_159[k] = -3.0 * gi_159[k]
                   + f_0 * ii_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, gi_160, gi_161, gi_162, gi_163, \
                         gi_164, ii_160, ii_161, ii_162, ii_163, \
                         ii_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -3.0 * gi_160[k]
                   + f_0 * ii_160[k];

        t_161[k] = -3.0 * gi_161[k]
                   + f_0 * ii_161[k];

        t_162[k] = -3.0 * gi_162[k]
                   + f_0 * ii_162[k];

        t_163[k] = -3.0 * gi_163[k]
                   + f_0 * ii_163[k];

        t_164[k] = -3.0 * gi_164[k]
                   + f_0 * ii_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, gi_165, gi_166, gi_167, gi_168, \
                         gi_169, ii_165, ii_166, ii_167, ii_168, \
                         ii_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -3.0 * gi_165[k]
                   + f_0 * ii_165[k];

        t_166[k] = -3.0 * gi_166[k]
                   + f_0 * ii_166[k];

        t_167[k] = -3.0 * gi_167[k]
                   + f_0 * ii_167[k];

        t_168[k] = -2.0 * gi_168[k]
                   + f_0 * ii_168[k];

        t_169[k] = -2.0 * gi_169[k]
                   + f_0 * ii_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, gi_170, gi_171, gi_172, gi_173, \
                         gi_174, ii_170, ii_171, ii_172, ii_173, \
                         ii_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * gi_170[k]
                   + f_0 * ii_170[k];

        t_171[k] = -2.0 * gi_171[k]
                   + f_0 * ii_171[k];

        t_172[k] = -2.0 * gi_172[k]
                   + f_0 * ii_172[k];

        t_173[k] = -2.0 * gi_173[k]
                   + f_0 * ii_173[k];

        t_174[k] = -2.0 * gi_174[k]
                   + f_0 * ii_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, gi_175, gi_176, gi_177, gi_178, \
                         gi_179, ii_175, ii_176, ii_177, ii_178, \
                         ii_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * gi_175[k]
                   + f_0 * ii_175[k];

        t_176[k] = -2.0 * gi_176[k]
                   + f_0 * ii_176[k];

        t_177[k] = -2.0 * gi_177[k]
                   + f_0 * ii_177[k];

        t_178[k] = -2.0 * gi_178[k]
                   + f_0 * ii_178[k];

        t_179[k] = -2.0 * gi_179[k]
                   + f_0 * ii_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gi_180, gi_181, gi_182, gi_183, \
                         gi_184, ii_180, ii_181, ii_182, ii_183, \
                         ii_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -2.0 * gi_180[k]
                   + f_0 * ii_180[k];

        t_181[k] = -2.0 * gi_181[k]
                   + f_0 * ii_181[k];

        t_182[k] = -2.0 * gi_182[k]
                   + f_0 * ii_182[k];

        t_183[k] = -2.0 * gi_183[k]
                   + f_0 * ii_183[k];

        t_184[k] = -2.0 * gi_184[k]
                   + f_0 * ii_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gi_185, gi_186, gi_187, gi_188, \
                         gi_189, ii_185, ii_186, ii_187, ii_188, \
                         ii_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -2.0 * gi_185[k]
                   + f_0 * ii_185[k];

        t_186[k] = -2.0 * gi_186[k]
                   + f_0 * ii_186[k];

        t_187[k] = -2.0 * gi_187[k]
                   + f_0 * ii_187[k];

        t_188[k] = -2.0 * gi_188[k]
                   + f_0 * ii_188[k];

        t_189[k] = -2.0 * gi_189[k]
                   + f_0 * ii_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gi_190, gi_191, gi_192, gi_193, \
                         gi_194, ii_190, ii_191, ii_192, ii_193, \
                         ii_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -2.0 * gi_190[k]
                   + f_0 * ii_190[k];

        t_191[k] = -2.0 * gi_191[k]
                   + f_0 * ii_191[k];

        t_192[k] = -2.0 * gi_192[k]
                   + f_0 * ii_192[k];

        t_193[k] = -2.0 * gi_193[k]
                   + f_0 * ii_193[k];

        t_194[k] = -2.0 * gi_194[k]
                   + f_0 * ii_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gi_195, gi_196, gi_197, gi_198, \
                         gi_199, ii_195, ii_196, ii_197, ii_198, \
                         ii_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -2.0 * gi_195[k]
                   + f_0 * ii_195[k];

        t_196[k] = -2.0 * gi_196[k]
                   + f_0 * ii_196[k];

        t_197[k] = -2.0 * gi_197[k]
                   + f_0 * ii_197[k];

        t_198[k] = -2.0 * gi_198[k]
                   + f_0 * ii_198[k];

        t_199[k] = -2.0 * gi_199[k]
                   + f_0 * ii_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gi_200, gi_201, gi_202, gi_203, \
                         gi_204, ii_200, ii_201, ii_202, ii_203, \
                         ii_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * gi_200[k]
                   + f_0 * ii_200[k];

        t_201[k] = -2.0 * gi_201[k]
                   + f_0 * ii_201[k];

        t_202[k] = -2.0 * gi_202[k]
                   + f_0 * ii_202[k];

        t_203[k] = -2.0 * gi_203[k]
                   + f_0 * ii_203[k];

        t_204[k] = -2.0 * gi_204[k]
                   + f_0 * ii_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gi_205, gi_206, gi_207, gi_208, \
                         gi_209, ii_205, ii_206, ii_207, ii_208, \
                         ii_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * gi_205[k]
                   + f_0 * ii_205[k];

        t_206[k] = -2.0 * gi_206[k]
                   + f_0 * ii_206[k];

        t_207[k] = -2.0 * gi_207[k]
                   + f_0 * ii_207[k];

        t_208[k] = -2.0 * gi_208[k]
                   + f_0 * ii_208[k];

        t_209[k] = -2.0 * gi_209[k]
                   + f_0 * ii_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, gi_210, gi_211, gi_212, gi_213, \
                         gi_214, ii_210, ii_211, ii_212, ii_213, \
                         ii_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * gi_210[k]
                   + f_0 * ii_210[k];

        t_211[k] = -2.0 * gi_211[k]
                   + f_0 * ii_211[k];

        t_212[k] = -2.0 * gi_212[k]
                   + f_0 * ii_212[k];

        t_213[k] = -2.0 * gi_213[k]
                   + f_0 * ii_213[k];

        t_214[k] = -2.0 * gi_214[k]
                   + f_0 * ii_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, gi_215, gi_216, gi_217, gi_218, \
                         gi_219, ii_215, ii_216, ii_217, ii_218, \
                         ii_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * gi_215[k]
                   + f_0 * ii_215[k];

        t_216[k] = -2.0 * gi_216[k]
                   + f_0 * ii_216[k];

        t_217[k] = -2.0 * gi_217[k]
                   + f_0 * ii_217[k];

        t_218[k] = -2.0 * gi_218[k]
                   + f_0 * ii_218[k];

        t_219[k] = -2.0 * gi_219[k]
                   + f_0 * ii_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, gi_220, gi_221, gi_222, gi_223, \
                         gi_224, ii_220, ii_221, ii_222, ii_223, \
                         ii_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * gi_220[k]
                   + f_0 * ii_220[k];

        t_221[k] = -2.0 * gi_221[k]
                   + f_0 * ii_221[k];

        t_222[k] = -2.0 * gi_222[k]
                   + f_0 * ii_222[k];

        t_223[k] = -2.0 * gi_223[k]
                   + f_0 * ii_223[k];

        t_224[k] = -2.0 * gi_224[k]
                   + f_0 * ii_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, gi_225, gi_226, gi_227, gi_228, \
                         gi_229, ii_225, ii_226, ii_227, ii_228, \
                         ii_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -2.0 * gi_225[k]
                   + f_0 * ii_225[k];

        t_226[k] = -2.0 * gi_226[k]
                   + f_0 * ii_226[k];

        t_227[k] = -2.0 * gi_227[k]
                   + f_0 * ii_227[k];

        t_228[k] = -2.0 * gi_228[k]
                   + f_0 * ii_228[k];

        t_229[k] = -2.0 * gi_229[k]
                   + f_0 * ii_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, gi_230, gi_231, gi_232, gi_233, \
                         gi_234, ii_230, ii_231, ii_232, ii_233, \
                         ii_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -2.0 * gi_230[k]
                   + f_0 * ii_230[k];

        t_231[k] = -2.0 * gi_231[k]
                   + f_0 * ii_231[k];

        t_232[k] = -2.0 * gi_232[k]
                   + f_0 * ii_232[k];

        t_233[k] = -2.0 * gi_233[k]
                   + f_0 * ii_233[k];

        t_234[k] = -2.0 * gi_234[k]
                   + f_0 * ii_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, gi_235, gi_236, gi_237, gi_238, \
                         gi_239, ii_235, ii_236, ii_237, ii_238, \
                         ii_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -2.0 * gi_235[k]
                   + f_0 * ii_235[k];

        t_236[k] = -2.0 * gi_236[k]
                   + f_0 * ii_236[k];

        t_237[k] = -2.0 * gi_237[k]
                   + f_0 * ii_237[k];

        t_238[k] = -2.0 * gi_238[k]
                   + f_0 * ii_238[k];

        t_239[k] = -2.0 * gi_239[k]
                   + f_0 * ii_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, gi_240, gi_241, gi_242, gi_243, \
                         gi_244, ii_240, ii_241, ii_242, ii_243, \
                         ii_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -2.0 * gi_240[k]
                   + f_0 * ii_240[k];

        t_241[k] = -2.0 * gi_241[k]
                   + f_0 * ii_241[k];

        t_242[k] = -2.0 * gi_242[k]
                   + f_0 * ii_242[k];

        t_243[k] = -2.0 * gi_243[k]
                   + f_0 * ii_243[k];

        t_244[k] = -2.0 * gi_244[k]
                   + f_0 * ii_244[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, gi_245, gi_246, gi_247, gi_248, \
                         gi_249, ii_245, ii_246, ii_247, ii_248, \
                         ii_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -2.0 * gi_245[k]
                   + f_0 * ii_245[k];

        t_246[k] = -2.0 * gi_246[k]
                   + f_0 * ii_246[k];

        t_247[k] = -2.0 * gi_247[k]
                   + f_0 * ii_247[k];

        t_248[k] = -2.0 * gi_248[k]
                   + f_0 * ii_248[k];

        t_249[k] = -2.0 * gi_249[k]
                   + f_0 * ii_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, gi_250, gi_251, gi_252, gi_253, \
                         gi_254, ii_250, ii_251, ii_252, ii_253, \
                         ii_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * gi_250[k]
                   + f_0 * ii_250[k];

        t_251[k] = -2.0 * gi_251[k]
                   + f_0 * ii_251[k];

        t_252[k] = -2.0 * gi_252[k]
                   + f_0 * ii_252[k];

        t_253[k] = -2.0 * gi_253[k]
                   + f_0 * ii_253[k];

        t_254[k] = -2.0 * gi_254[k]
                   + f_0 * ii_254[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, gi_255, gi_256, gi_257, gi_258, \
                         gi_259, ii_255, ii_256, ii_257, ii_258, \
                         ii_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * gi_255[k]
                   + f_0 * ii_255[k];

        t_256[k] = -2.0 * gi_256[k]
                   + f_0 * ii_256[k];

        t_257[k] = -2.0 * gi_257[k]
                   + f_0 * ii_257[k];

        t_258[k] = -2.0 * gi_258[k]
                   + f_0 * ii_258[k];

        t_259[k] = -2.0 * gi_259[k]
                   + f_0 * ii_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, gi_260, gi_261, gi_262, gi_263, \
                         gi_264, ii_260, ii_261, ii_262, ii_263, \
                         ii_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -2.0 * gi_260[k]
                   + f_0 * ii_260[k];

        t_261[k] = -2.0 * gi_261[k]
                   + f_0 * ii_261[k];

        t_262[k] = -2.0 * gi_262[k]
                   + f_0 * ii_262[k];

        t_263[k] = -2.0 * gi_263[k]
                   + f_0 * ii_263[k];

        t_264[k] = -2.0 * gi_264[k]
                   + f_0 * ii_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, gi_265, gi_266, gi_267, gi_268, \
                         gi_269, ii_265, ii_266, ii_267, ii_268, \
                         ii_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -2.0 * gi_265[k]
                   + f_0 * ii_265[k];

        t_266[k] = -2.0 * gi_266[k]
                   + f_0 * ii_266[k];

        t_267[k] = -2.0 * gi_267[k]
                   + f_0 * ii_267[k];

        t_268[k] = -2.0 * gi_268[k]
                   + f_0 * ii_268[k];

        t_269[k] = -2.0 * gi_269[k]
                   + f_0 * ii_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, gi_270, gi_271, gi_272, gi_273, \
                         gi_274, ii_270, ii_271, ii_272, ii_273, \
                         ii_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -2.0 * gi_270[k]
                   + f_0 * ii_270[k];

        t_271[k] = -2.0 * gi_271[k]
                   + f_0 * ii_271[k];

        t_272[k] = -2.0 * gi_272[k]
                   + f_0 * ii_272[k];

        t_273[k] = -2.0 * gi_273[k]
                   + f_0 * ii_273[k];

        t_274[k] = -2.0 * gi_274[k]
                   + f_0 * ii_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, gi_275, gi_276, gi_277, gi_278, \
                         gi_279, ii_275, ii_276, ii_277, ii_278, \
                         ii_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -2.0 * gi_275[k]
                   + f_0 * ii_275[k];

        t_276[k] = -2.0 * gi_276[k]
                   + f_0 * ii_276[k];

        t_277[k] = -2.0 * gi_277[k]
                   + f_0 * ii_277[k];

        t_278[k] = -2.0 * gi_278[k]
                   + f_0 * ii_278[k];

        t_279[k] = -2.0 * gi_279[k]
                   + f_0 * ii_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, gi_280, gi_281, gi_282, gi_283, \
                         gi_284, ii_280, ii_281, ii_282, ii_283, \
                         ii_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -gi_280[k]
                   + f_0 * ii_280[k];

        t_281[k] = -gi_281[k]
                   + f_0 * ii_281[k];

        t_282[k] = -gi_282[k]
                   + f_0 * ii_282[k];

        t_283[k] = -gi_283[k]
                   + f_0 * ii_283[k];

        t_284[k] = -gi_284[k]
                   + f_0 * ii_284[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, gi_285, gi_286, gi_287, gi_288, \
                         gi_289, ii_285, ii_286, ii_287, ii_288, \
                         ii_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -gi_285[k]
                   + f_0 * ii_285[k];

        t_286[k] = -gi_286[k]
                   + f_0 * ii_286[k];

        t_287[k] = -gi_287[k]
                   + f_0 * ii_287[k];

        t_288[k] = -gi_288[k]
                   + f_0 * ii_288[k];

        t_289[k] = -gi_289[k]
                   + f_0 * ii_289[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, gi_290, gi_291, gi_292, gi_293, \
                         gi_294, ii_290, ii_291, ii_292, ii_293, \
                         ii_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -gi_290[k]
                   + f_0 * ii_290[k];

        t_291[k] = -gi_291[k]
                   + f_0 * ii_291[k];

        t_292[k] = -gi_292[k]
                   + f_0 * ii_292[k];

        t_293[k] = -gi_293[k]
                   + f_0 * ii_293[k];

        t_294[k] = -gi_294[k]
                   + f_0 * ii_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, gi_295, gi_296, gi_297, gi_298, \
                         gi_299, ii_295, ii_296, ii_297, ii_298, \
                         ii_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -gi_295[k]
                   + f_0 * ii_295[k];

        t_296[k] = -gi_296[k]
                   + f_0 * ii_296[k];

        t_297[k] = -gi_297[k]
                   + f_0 * ii_297[k];

        t_298[k] = -gi_298[k]
                   + f_0 * ii_298[k];

        t_299[k] = -gi_299[k]
                   + f_0 * ii_299[k];
    }
}

static auto
compute_prim_geom_10_hi_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_300 = buffer.data(gi + 300);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_307 = buffer.data(gi + 307);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_309 = buffer.data(gi + 309);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_312 = buffer.data(gi + 312);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_315 = buffer.data(gi + 315);
    const auto *gi_316 = buffer.data(gi + 316);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_319 = buffer.data(gi + 319);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_321 = buffer.data(gi + 321);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_323 = buffer.data(gi + 323);
    const auto *gi_324 = buffer.data(gi + 324);
    const auto *gi_325 = buffer.data(gi + 325);
    const auto *gi_326 = buffer.data(gi + 326);
    const auto *gi_327 = buffer.data(gi + 327);
    const auto *gi_328 = buffer.data(gi + 328);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_335 = buffer.data(gi + 335);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_337 = buffer.data(gi + 337);
    const auto *gi_338 = buffer.data(gi + 338);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_340 = buffer.data(gi + 340);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_343 = buffer.data(gi + 343);
    const auto *gi_344 = buffer.data(gi + 344);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_347 = buffer.data(gi + 347);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_349 = buffer.data(gi + 349);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_352 = buffer.data(gi + 352);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);
    const auto *gi_355 = buffer.data(gi + 355);
    const auto *gi_356 = buffer.data(gi + 356);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_363 = buffer.data(gi + 363);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_365 = buffer.data(gi + 365);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_368 = buffer.data(gi + 368);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_371 = buffer.data(gi + 371);
    const auto *gi_372 = buffer.data(gi + 372);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_375 = buffer.data(gi + 375);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_377 = buffer.data(gi + 377);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_379 = buffer.data(gi + 379);
    const auto *gi_380 = buffer.data(gi + 380);
    const auto *gi_381 = buffer.data(gi + 381);
    const auto *gi_382 = buffer.data(gi + 382);
    const auto *gi_383 = buffer.data(gi + 383);
    const auto *gi_384 = buffer.data(gi + 384);
    const auto *gi_385 = buffer.data(gi + 385);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_387 = buffer.data(gi + 387);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_389 = buffer.data(gi + 389);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_391 = buffer.data(gi + 391);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_393 = buffer.data(gi + 393);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_395 = buffer.data(gi + 395);
    const auto *gi_396 = buffer.data(gi + 396);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_399 = buffer.data(gi + 399);
    const auto *gi_400 = buffer.data(gi + 400);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_403 = buffer.data(gi + 403);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_405 = buffer.data(gi + 405);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_408 = buffer.data(gi + 408);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_411 = buffer.data(gi + 411);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

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

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, gi_300, gi_301, gi_302, gi_303, \
                         gi_304, ii_300, ii_301, ii_302, ii_303, \
                         ii_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -gi_300[k]
                   + f_0 * ii_300[k];

        t_301[k] = -gi_301[k]
                   + f_0 * ii_301[k];

        t_302[k] = -gi_302[k]
                   + f_0 * ii_302[k];

        t_303[k] = -gi_303[k]
                   + f_0 * ii_303[k];

        t_304[k] = -gi_304[k]
                   + f_0 * ii_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, gi_305, gi_306, gi_307, gi_308, \
                         gi_309, ii_305, ii_306, ii_307, ii_308, \
                         ii_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -gi_305[k]
                   + f_0 * ii_305[k];

        t_306[k] = -gi_306[k]
                   + f_0 * ii_306[k];

        t_307[k] = -gi_307[k]
                   + f_0 * ii_307[k];

        t_308[k] = -gi_308[k]
                   + f_0 * ii_308[k];

        t_309[k] = -gi_309[k]
                   + f_0 * ii_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, gi_310, gi_311, gi_312, gi_313, \
                         gi_314, ii_310, ii_311, ii_312, ii_313, \
                         ii_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -gi_310[k]
                   + f_0 * ii_310[k];

        t_311[k] = -gi_311[k]
                   + f_0 * ii_311[k];

        t_312[k] = -gi_312[k]
                   + f_0 * ii_312[k];

        t_313[k] = -gi_313[k]
                   + f_0 * ii_313[k];

        t_314[k] = -gi_314[k]
                   + f_0 * ii_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, gi_315, gi_316, gi_317, gi_318, \
                         gi_319, ii_315, ii_316, ii_317, ii_318, \
                         ii_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -gi_315[k]
                   + f_0 * ii_315[k];

        t_316[k] = -gi_316[k]
                   + f_0 * ii_316[k];

        t_317[k] = -gi_317[k]
                   + f_0 * ii_317[k];

        t_318[k] = -gi_318[k]
                   + f_0 * ii_318[k];

        t_319[k] = -gi_319[k]
                   + f_0 * ii_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, gi_320, gi_321, gi_322, gi_323, \
                         gi_324, ii_320, ii_321, ii_322, ii_323, \
                         ii_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -gi_320[k]
                   + f_0 * ii_320[k];

        t_321[k] = -gi_321[k]
                   + f_0 * ii_321[k];

        t_322[k] = -gi_322[k]
                   + f_0 * ii_322[k];

        t_323[k] = -gi_323[k]
                   + f_0 * ii_323[k];

        t_324[k] = -gi_324[k]
                   + f_0 * ii_324[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, gi_325, gi_326, gi_327, gi_328, \
                         gi_329, ii_325, ii_326, ii_327, ii_328, \
                         ii_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -gi_325[k]
                   + f_0 * ii_325[k];

        t_326[k] = -gi_326[k]
                   + f_0 * ii_326[k];

        t_327[k] = -gi_327[k]
                   + f_0 * ii_327[k];

        t_328[k] = -gi_328[k]
                   + f_0 * ii_328[k];

        t_329[k] = -gi_329[k]
                   + f_0 * ii_329[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, gi_330, gi_331, gi_332, gi_333, \
                         gi_334, ii_330, ii_331, ii_332, ii_333, \
                         ii_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -gi_330[k]
                   + f_0 * ii_330[k];

        t_331[k] = -gi_331[k]
                   + f_0 * ii_331[k];

        t_332[k] = -gi_332[k]
                   + f_0 * ii_332[k];

        t_333[k] = -gi_333[k]
                   + f_0 * ii_333[k];

        t_334[k] = -gi_334[k]
                   + f_0 * ii_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, gi_335, gi_336, gi_337, gi_338, \
                         gi_339, ii_335, ii_336, ii_337, ii_338, \
                         ii_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -gi_335[k]
                   + f_0 * ii_335[k];

        t_336[k] = -gi_336[k]
                   + f_0 * ii_336[k];

        t_337[k] = -gi_337[k]
                   + f_0 * ii_337[k];

        t_338[k] = -gi_338[k]
                   + f_0 * ii_338[k];

        t_339[k] = -gi_339[k]
                   + f_0 * ii_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, gi_340, gi_341, gi_342, gi_343, \
                         gi_344, ii_340, ii_341, ii_342, ii_343, \
                         ii_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -gi_340[k]
                   + f_0 * ii_340[k];

        t_341[k] = -gi_341[k]
                   + f_0 * ii_341[k];

        t_342[k] = -gi_342[k]
                   + f_0 * ii_342[k];

        t_343[k] = -gi_343[k]
                   + f_0 * ii_343[k];

        t_344[k] = -gi_344[k]
                   + f_0 * ii_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, gi_345, gi_346, gi_347, gi_348, \
                         gi_349, ii_345, ii_346, ii_347, ii_348, \
                         ii_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -gi_345[k]
                   + f_0 * ii_345[k];

        t_346[k] = -gi_346[k]
                   + f_0 * ii_346[k];

        t_347[k] = -gi_347[k]
                   + f_0 * ii_347[k];

        t_348[k] = -gi_348[k]
                   + f_0 * ii_348[k];

        t_349[k] = -gi_349[k]
                   + f_0 * ii_349[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, gi_350, gi_351, gi_352, gi_353, \
                         gi_354, ii_350, ii_351, ii_352, ii_353, \
                         ii_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -gi_350[k]
                   + f_0 * ii_350[k];

        t_351[k] = -gi_351[k]
                   + f_0 * ii_351[k];

        t_352[k] = -gi_352[k]
                   + f_0 * ii_352[k];

        t_353[k] = -gi_353[k]
                   + f_0 * ii_353[k];

        t_354[k] = -gi_354[k]
                   + f_0 * ii_354[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, gi_355, gi_356, gi_357, gi_358, \
                         gi_359, ii_355, ii_356, ii_357, ii_358, \
                         ii_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -gi_355[k]
                   + f_0 * ii_355[k];

        t_356[k] = -gi_356[k]
                   + f_0 * ii_356[k];

        t_357[k] = -gi_357[k]
                   + f_0 * ii_357[k];

        t_358[k] = -gi_358[k]
                   + f_0 * ii_358[k];

        t_359[k] = -gi_359[k]
                   + f_0 * ii_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, gi_360, gi_361, gi_362, gi_363, \
                         gi_364, ii_360, ii_361, ii_362, ii_363, \
                         ii_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -gi_360[k]
                   + f_0 * ii_360[k];

        t_361[k] = -gi_361[k]
                   + f_0 * ii_361[k];

        t_362[k] = -gi_362[k]
                   + f_0 * ii_362[k];

        t_363[k] = -gi_363[k]
                   + f_0 * ii_363[k];

        t_364[k] = -gi_364[k]
                   + f_0 * ii_364[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, gi_365, gi_366, gi_367, gi_368, \
                         gi_369, ii_365, ii_366, ii_367, ii_368, \
                         ii_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -gi_365[k]
                   + f_0 * ii_365[k];

        t_366[k] = -gi_366[k]
                   + f_0 * ii_366[k];

        t_367[k] = -gi_367[k]
                   + f_0 * ii_367[k];

        t_368[k] = -gi_368[k]
                   + f_0 * ii_368[k];

        t_369[k] = -gi_369[k]
                   + f_0 * ii_369[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, gi_370, gi_371, gi_372, gi_373, \
                         gi_374, ii_370, ii_371, ii_372, ii_373, \
                         ii_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -gi_370[k]
                   + f_0 * ii_370[k];

        t_371[k] = -gi_371[k]
                   + f_0 * ii_371[k];

        t_372[k] = -gi_372[k]
                   + f_0 * ii_372[k];

        t_373[k] = -gi_373[k]
                   + f_0 * ii_373[k];

        t_374[k] = -gi_374[k]
                   + f_0 * ii_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, gi_375, gi_376, gi_377, gi_378, \
                         gi_379, ii_375, ii_376, ii_377, ii_378, \
                         ii_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -gi_375[k]
                   + f_0 * ii_375[k];

        t_376[k] = -gi_376[k]
                   + f_0 * ii_376[k];

        t_377[k] = -gi_377[k]
                   + f_0 * ii_377[k];

        t_378[k] = -gi_378[k]
                   + f_0 * ii_378[k];

        t_379[k] = -gi_379[k]
                   + f_0 * ii_379[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, gi_380, gi_381, gi_382, gi_383, \
                         gi_384, ii_380, ii_381, ii_382, ii_383, \
                         ii_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -gi_380[k]
                   + f_0 * ii_380[k];

        t_381[k] = -gi_381[k]
                   + f_0 * ii_381[k];

        t_382[k] = -gi_382[k]
                   + f_0 * ii_382[k];

        t_383[k] = -gi_383[k]
                   + f_0 * ii_383[k];

        t_384[k] = -gi_384[k]
                   + f_0 * ii_384[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, gi_385, gi_386, gi_387, gi_388, \
                         gi_389, ii_385, ii_386, ii_387, ii_388, \
                         ii_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -gi_385[k]
                   + f_0 * ii_385[k];

        t_386[k] = -gi_386[k]
                   + f_0 * ii_386[k];

        t_387[k] = -gi_387[k]
                   + f_0 * ii_387[k];

        t_388[k] = -gi_388[k]
                   + f_0 * ii_388[k];

        t_389[k] = -gi_389[k]
                   + f_0 * ii_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, gi_390, gi_391, gi_392, gi_393, \
                         gi_394, ii_390, ii_391, ii_392, ii_393, \
                         ii_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -gi_390[k]
                   + f_0 * ii_390[k];

        t_391[k] = -gi_391[k]
                   + f_0 * ii_391[k];

        t_392[k] = -gi_392[k]
                   + f_0 * ii_392[k];

        t_393[k] = -gi_393[k]
                   + f_0 * ii_393[k];

        t_394[k] = -gi_394[k]
                   + f_0 * ii_394[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, gi_395, gi_396, gi_397, gi_398, \
                         gi_399, ii_395, ii_396, ii_397, ii_398, \
                         ii_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = -gi_395[k]
                   + f_0 * ii_395[k];

        t_396[k] = -gi_396[k]
                   + f_0 * ii_396[k];

        t_397[k] = -gi_397[k]
                   + f_0 * ii_397[k];

        t_398[k] = -gi_398[k]
                   + f_0 * ii_398[k];

        t_399[k] = -gi_399[k]
                   + f_0 * ii_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, gi_400, gi_401, gi_402, gi_403, \
                         gi_404, ii_400, ii_401, ii_402, ii_403, \
                         ii_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = -gi_400[k]
                   + f_0 * ii_400[k];

        t_401[k] = -gi_401[k]
                   + f_0 * ii_401[k];

        t_402[k] = -gi_402[k]
                   + f_0 * ii_402[k];

        t_403[k] = -gi_403[k]
                   + f_0 * ii_403[k];

        t_404[k] = -gi_404[k]
                   + f_0 * ii_404[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, gi_405, gi_406, gi_407, gi_408, \
                         gi_409, ii_405, ii_406, ii_407, ii_408, \
                         ii_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = -gi_405[k]
                   + f_0 * ii_405[k];

        t_406[k] = -gi_406[k]
                   + f_0 * ii_406[k];

        t_407[k] = -gi_407[k]
                   + f_0 * ii_407[k];

        t_408[k] = -gi_408[k]
                   + f_0 * ii_408[k];

        t_409[k] = -gi_409[k]
                   + f_0 * ii_409[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, gi_410, gi_411, gi_412, gi_413, \
                         gi_414, ii_410, ii_411, ii_412, ii_413, \
                         ii_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -gi_410[k]
                   + f_0 * ii_410[k];

        t_411[k] = -gi_411[k]
                   + f_0 * ii_411[k];

        t_412[k] = -gi_412[k]
                   + f_0 * ii_412[k];

        t_413[k] = -gi_413[k]
                   + f_0 * ii_413[k];

        t_414[k] = -gi_414[k]
                   + f_0 * ii_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, gi_415, gi_416, gi_417, gi_418, \
                         gi_419, ii_415, ii_416, ii_417, ii_418, \
                         ii_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = -gi_415[k]
                   + f_0 * ii_415[k];

        t_416[k] = -gi_416[k]
                   + f_0 * ii_416[k];

        t_417[k] = -gi_417[k]
                   + f_0 * ii_417[k];

        t_418[k] = -gi_418[k]
                   + f_0 * ii_418[k];

        t_419[k] = -gi_419[k]
                   + f_0 * ii_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, t_425, t_426, t_427, ii_420, \
                         ii_421, ii_422, ii_423, ii_424, ii_425, ii_426, \
                         ii_427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_0 * ii_420[k];

        t_421[k] = f_0 * ii_421[k];

        t_422[k] = f_0 * ii_422[k];

        t_423[k] = f_0 * ii_423[k];

        t_424[k] = f_0 * ii_424[k];

        t_425[k] = f_0 * ii_425[k];

        t_426[k] = f_0 * ii_426[k];

        t_427[k] = f_0 * ii_427[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, t_433, t_434, t_435, ii_428, \
                         ii_429, ii_430, ii_431, ii_432, ii_433, ii_434, \
                         ii_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_0 * ii_428[k];

        t_429[k] = f_0 * ii_429[k];

        t_430[k] = f_0 * ii_430[k];

        t_431[k] = f_0 * ii_431[k];

        t_432[k] = f_0 * ii_432[k];

        t_433[k] = f_0 * ii_433[k];

        t_434[k] = f_0 * ii_434[k];

        t_435[k] = f_0 * ii_435[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, t_441, t_442, t_443, ii_436, \
                         ii_437, ii_438, ii_439, ii_440, ii_441, ii_442, \
                         ii_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_0 * ii_436[k];

        t_437[k] = f_0 * ii_437[k];

        t_438[k] = f_0 * ii_438[k];

        t_439[k] = f_0 * ii_439[k];

        t_440[k] = f_0 * ii_440[k];

        t_441[k] = f_0 * ii_441[k];

        t_442[k] = f_0 * ii_442[k];

        t_443[k] = f_0 * ii_443[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, t_449, t_450, t_451, ii_444, \
                         ii_445, ii_446, ii_447, ii_448, ii_449, ii_450, \
                         ii_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_0 * ii_444[k];

        t_445[k] = f_0 * ii_445[k];

        t_446[k] = f_0 * ii_446[k];

        t_447[k] = f_0 * ii_447[k];

        t_448[k] = f_0 * ii_448[k];

        t_449[k] = f_0 * ii_449[k];

        t_450[k] = f_0 * ii_450[k];

        t_451[k] = f_0 * ii_451[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, t_457, t_458, t_459, ii_452, \
                         ii_453, ii_454, ii_455, ii_456, ii_457, ii_458, \
                         ii_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_0 * ii_452[k];

        t_453[k] = f_0 * ii_453[k];

        t_454[k] = f_0 * ii_454[k];

        t_455[k] = f_0 * ii_455[k];

        t_456[k] = f_0 * ii_456[k];

        t_457[k] = f_0 * ii_457[k];

        t_458[k] = f_0 * ii_458[k];

        t_459[k] = f_0 * ii_459[k];
    }
}

static auto
compute_prim_geom_10_hi_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t ii, const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, t_465, t_466, t_467, ii_460, \
                         ii_461, ii_462, ii_463, ii_464, ii_465, ii_466, \
                         ii_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_0 * ii_460[k];

        t_461[k] = f_0 * ii_461[k];

        t_462[k] = f_0 * ii_462[k];

        t_463[k] = f_0 * ii_463[k];

        t_464[k] = f_0 * ii_464[k];

        t_465[k] = f_0 * ii_465[k];

        t_466[k] = f_0 * ii_466[k];

        t_467[k] = f_0 * ii_467[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, t_471, t_472, t_473, t_474, t_475, ii_468, \
                         ii_469, ii_470, ii_471, ii_472, ii_473, ii_474, \
                         ii_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_0 * ii_468[k];

        t_469[k] = f_0 * ii_469[k];

        t_470[k] = f_0 * ii_470[k];

        t_471[k] = f_0 * ii_471[k];

        t_472[k] = f_0 * ii_472[k];

        t_473[k] = f_0 * ii_473[k];

        t_474[k] = f_0 * ii_474[k];

        t_475[k] = f_0 * ii_475[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, t_481, t_482, t_483, ii_476, \
                         ii_477, ii_478, ii_479, ii_480, ii_481, ii_482, \
                         ii_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_0 * ii_476[k];

        t_477[k] = f_0 * ii_477[k];

        t_478[k] = f_0 * ii_478[k];

        t_479[k] = f_0 * ii_479[k];

        t_480[k] = f_0 * ii_480[k];

        t_481[k] = f_0 * ii_481[k];

        t_482[k] = f_0 * ii_482[k];

        t_483[k] = f_0 * ii_483[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, t_489, t_490, t_491, ii_484, \
                         ii_485, ii_486, ii_487, ii_488, ii_489, ii_490, \
                         ii_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_0 * ii_484[k];

        t_485[k] = f_0 * ii_485[k];

        t_486[k] = f_0 * ii_486[k];

        t_487[k] = f_0 * ii_487[k];

        t_488[k] = f_0 * ii_488[k];

        t_489[k] = f_0 * ii_489[k];

        t_490[k] = f_0 * ii_490[k];

        t_491[k] = f_0 * ii_491[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, t_497, t_498, t_499, ii_492, \
                         ii_493, ii_494, ii_495, ii_496, ii_497, ii_498, \
                         ii_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_0 * ii_492[k];

        t_493[k] = f_0 * ii_493[k];

        t_494[k] = f_0 * ii_494[k];

        t_495[k] = f_0 * ii_495[k];

        t_496[k] = f_0 * ii_496[k];

        t_497[k] = f_0 * ii_497[k];

        t_498[k] = f_0 * ii_498[k];

        t_499[k] = f_0 * ii_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, t_505, t_506, t_507, ii_500, \
                         ii_501, ii_502, ii_503, ii_504, ii_505, ii_506, \
                         ii_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_0 * ii_500[k];

        t_501[k] = f_0 * ii_501[k];

        t_502[k] = f_0 * ii_502[k];

        t_503[k] = f_0 * ii_503[k];

        t_504[k] = f_0 * ii_504[k];

        t_505[k] = f_0 * ii_505[k];

        t_506[k] = f_0 * ii_506[k];

        t_507[k] = f_0 * ii_507[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, t_512, t_513, t_514, t_515, ii_508, \
                         ii_509, ii_510, ii_511, ii_512, ii_513, ii_514, \
                         ii_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_0 * ii_508[k];

        t_509[k] = f_0 * ii_509[k];

        t_510[k] = f_0 * ii_510[k];

        t_511[k] = f_0 * ii_511[k];

        t_512[k] = f_0 * ii_512[k];

        t_513[k] = f_0 * ii_513[k];

        t_514[k] = f_0 * ii_514[k];

        t_515[k] = f_0 * ii_515[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, t_520, t_521, t_522, t_523, ii_516, \
                         ii_517, ii_518, ii_519, ii_520, ii_521, ii_522, \
                         ii_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_0 * ii_516[k];

        t_517[k] = f_0 * ii_517[k];

        t_518[k] = f_0 * ii_518[k];

        t_519[k] = f_0 * ii_519[k];

        t_520[k] = f_0 * ii_520[k];

        t_521[k] = f_0 * ii_521[k];

        t_522[k] = f_0 * ii_522[k];

        t_523[k] = f_0 * ii_523[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, t_529, t_530, t_531, ii_524, \
                         ii_525, ii_526, ii_527, ii_528, ii_529, ii_530, \
                         ii_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_0 * ii_524[k];

        t_525[k] = f_0 * ii_525[k];

        t_526[k] = f_0 * ii_526[k];

        t_527[k] = f_0 * ii_527[k];

        t_528[k] = f_0 * ii_528[k];

        t_529[k] = f_0 * ii_529[k];

        t_530[k] = f_0 * ii_530[k];

        t_531[k] = f_0 * ii_531[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, t_537, t_538, t_539, ii_532, \
                         ii_533, ii_534, ii_535, ii_536, ii_537, ii_538, \
                         ii_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_0 * ii_532[k];

        t_533[k] = f_0 * ii_533[k];

        t_534[k] = f_0 * ii_534[k];

        t_535[k] = f_0 * ii_535[k];

        t_536[k] = f_0 * ii_536[k];

        t_537[k] = f_0 * ii_537[k];

        t_538[k] = f_0 * ii_538[k];

        t_539[k] = f_0 * ii_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, t_545, t_546, t_547, ii_540, \
                         ii_541, ii_542, ii_543, ii_544, ii_545, ii_546, \
                         ii_547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_0 * ii_540[k];

        t_541[k] = f_0 * ii_541[k];

        t_542[k] = f_0 * ii_542[k];

        t_543[k] = f_0 * ii_543[k];

        t_544[k] = f_0 * ii_544[k];

        t_545[k] = f_0 * ii_545[k];

        t_546[k] = f_0 * ii_546[k];

        t_547[k] = f_0 * ii_547[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, t_553, t_554, t_555, ii_548, \
                         ii_549, ii_550, ii_551, ii_552, ii_553, ii_554, \
                         ii_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_0 * ii_548[k];

        t_549[k] = f_0 * ii_549[k];

        t_550[k] = f_0 * ii_550[k];

        t_551[k] = f_0 * ii_551[k];

        t_552[k] = f_0 * ii_552[k];

        t_553[k] = f_0 * ii_553[k];

        t_554[k] = f_0 * ii_554[k];

        t_555[k] = f_0 * ii_555[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, t_561, t_562, t_563, ii_556, \
                         ii_557, ii_558, ii_559, ii_560, ii_561, ii_562, \
                         ii_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_0 * ii_556[k];

        t_557[k] = f_0 * ii_557[k];

        t_558[k] = f_0 * ii_558[k];

        t_559[k] = f_0 * ii_559[k];

        t_560[k] = f_0 * ii_560[k];

        t_561[k] = f_0 * ii_561[k];

        t_562[k] = f_0 * ii_562[k];

        t_563[k] = f_0 * ii_563[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, t_568, t_569, t_570, t_571, ii_564, \
                         ii_565, ii_566, ii_567, ii_568, ii_569, ii_570, \
                         ii_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_0 * ii_564[k];

        t_565[k] = f_0 * ii_565[k];

        t_566[k] = f_0 * ii_566[k];

        t_567[k] = f_0 * ii_567[k];

        t_568[k] = f_0 * ii_568[k];

        t_569[k] = f_0 * ii_569[k];

        t_570[k] = f_0 * ii_570[k];

        t_571[k] = f_0 * ii_571[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, t_577, t_578, t_579, ii_572, \
                         ii_573, ii_574, ii_575, ii_576, ii_577, ii_578, \
                         ii_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_0 * ii_572[k];

        t_573[k] = f_0 * ii_573[k];

        t_574[k] = f_0 * ii_574[k];

        t_575[k] = f_0 * ii_575[k];

        t_576[k] = f_0 * ii_576[k];

        t_577[k] = f_0 * ii_577[k];

        t_578[k] = f_0 * ii_578[k];

        t_579[k] = f_0 * ii_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, t_585, t_586, t_587, ii_580, \
                         ii_581, ii_582, ii_583, ii_584, ii_585, ii_586, \
                         ii_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_0 * ii_580[k];

        t_581[k] = f_0 * ii_581[k];

        t_582[k] = f_0 * ii_582[k];

        t_583[k] = f_0 * ii_583[k];

        t_584[k] = f_0 * ii_584[k];

        t_585[k] = f_0 * ii_585[k];

        t_586[k] = f_0 * ii_586[k];

        t_587[k] = f_0 * ii_587[k];
    }
}

auto
compute_prim_geom_10_hi_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t gi, const size_t ii,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hi_electron_repulsion_0_piece0(buffer, target, gi, ii, ncols, alpha);

    compute_prim_geom_10_hi_electron_repulsion_0_piece1(buffer, target, gi, ii, ncols, alpha);

    compute_prim_geom_10_hi_electron_repulsion_0_piece2(buffer, target, gi, ii, ncols, alpha);

    compute_prim_geom_10_hi_electron_repulsion_0_piece3(buffer, target, ii, ncols, alpha);
}

static auto
compute_prim_geom_10_hi_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_1 = buffer.data(gi + 1);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_41 = buffer.data(gi + 41);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_44 = buffer.data(gi + 44);
    const auto *gi_45 = buffer.data(gi + 45);
    const auto *gi_46 = buffer.data(gi + 46);
    const auto *gi_47 = buffer.data(gi + 47);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_50 = buffer.data(gi + 50);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_60 = buffer.data(gi + 60);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_63 = buffer.data(gi + 63);
    const auto *gi_64 = buffer.data(gi + 64);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_67 = buffer.data(gi + 67);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_69 = buffer.data(gi + 69);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_71 = buffer.data(gi + 71);
    const auto *gi_72 = buffer.data(gi + 72);
    const auto *gi_73 = buffer.data(gi + 73);
    const auto *gi_74 = buffer.data(gi + 74);
    const auto *gi_75 = buffer.data(gi + 75);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_77 = buffer.data(gi + 77);
    const auto *gi_78 = buffer.data(gi + 78);
    const auto *gi_79 = buffer.data(gi + 79);
    const auto *gi_80 = buffer.data(gi + 80);
    const auto *gi_81 = buffer.data(gi + 81);
    const auto *gi_82 = buffer.data(gi + 82);
    const auto *gi_83 = buffer.data(gi + 83);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_85 = buffer.data(gi + 85);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_88 = buffer.data(gi + 88);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_90 = buffer.data(gi + 90);

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
    const auto *ii_280 = buffer.data(ii + 280);
    const auto *ii_281 = buffer.data(ii + 281);
    const auto *ii_282 = buffer.data(ii + 282);
    const auto *ii_283 = buffer.data(ii + 283);
    const auto *ii_284 = buffer.data(ii + 284);
    const auto *ii_285 = buffer.data(ii + 285);
    const auto *ii_286 = buffer.data(ii + 286);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ii_28, ii_29, ii_30, ii_31, \
                         ii_32, ii_33, ii_34, ii_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ii_28[k];

        t_1[k] = f_0 * ii_29[k];

        t_2[k] = f_0 * ii_30[k];

        t_3[k] = f_0 * ii_31[k];

        t_4[k] = f_0 * ii_32[k];

        t_5[k] = f_0 * ii_33[k];

        t_6[k] = f_0 * ii_34[k];

        t_7[k] = f_0 * ii_35[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ii_36, ii_37, ii_38, \
                         ii_39, ii_40, ii_41, ii_42, ii_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ii_36[k];

        t_9[k] = f_0 * ii_37[k];

        t_10[k] = f_0 * ii_38[k];

        t_11[k] = f_0 * ii_39[k];

        t_12[k] = f_0 * ii_40[k];

        t_13[k] = f_0 * ii_41[k];

        t_14[k] = f_0 * ii_42[k];

        t_15[k] = f_0 * ii_43[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, ii_44, ii_45, ii_46, \
                         ii_47, ii_48, ii_49, ii_50, ii_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ii_44[k];

        t_17[k] = f_0 * ii_45[k];

        t_18[k] = f_0 * ii_46[k];

        t_19[k] = f_0 * ii_47[k];

        t_20[k] = f_0 * ii_48[k];

        t_21[k] = f_0 * ii_49[k];

        t_22[k] = f_0 * ii_50[k];

        t_23[k] = f_0 * ii_51[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, gi_0, gi_1, ii_52, ii_53, ii_54, \
                         ii_55, ii_84, ii_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * ii_52[k];

        t_25[k] = f_0 * ii_53[k];

        t_26[k] = f_0 * ii_54[k];

        t_27[k] = f_0 * ii_55[k];

        t_28[k] = -gi_0[k]
                  + f_0 * ii_84[k];

        t_29[k] = -gi_1[k]
                  + f_0 * ii_85[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, gi_2, gi_3, gi_4, gi_5, gi_6, ii_86, \
                         ii_87, ii_88, ii_89, ii_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -gi_2[k]
                  + f_0 * ii_86[k];

        t_31[k] = -gi_3[k]
                  + f_0 * ii_87[k];

        t_32[k] = -gi_4[k]
                  + f_0 * ii_88[k];

        t_33[k] = -gi_5[k]
                  + f_0 * ii_89[k];

        t_34[k] = -gi_6[k]
                  + f_0 * ii_90[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gi_7, gi_8, gi_9, gi_10, gi_11, ii_91, \
                         ii_92, ii_93, ii_94, ii_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -gi_7[k]
                  + f_0 * ii_91[k];

        t_36[k] = -gi_8[k]
                  + f_0 * ii_92[k];

        t_37[k] = -gi_9[k]
                  + f_0 * ii_93[k];

        t_38[k] = -gi_10[k]
                  + f_0 * ii_94[k];

        t_39[k] = -gi_11[k]
                  + f_0 * ii_95[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gi_12, gi_13, gi_14, gi_15, gi_16, \
                         ii_96, ii_97, ii_98, ii_99, ii_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -gi_12[k]
                  + f_0 * ii_96[k];

        t_41[k] = -gi_13[k]
                  + f_0 * ii_97[k];

        t_42[k] = -gi_14[k]
                  + f_0 * ii_98[k];

        t_43[k] = -gi_15[k]
                  + f_0 * ii_99[k];

        t_44[k] = -gi_16[k]
                  + f_0 * ii_100[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gi_17, gi_18, gi_19, gi_20, gi_21, \
                         ii_101, ii_102, ii_103, ii_104, ii_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -gi_17[k]
                  + f_0 * ii_101[k];

        t_46[k] = -gi_18[k]
                  + f_0 * ii_102[k];

        t_47[k] = -gi_19[k]
                  + f_0 * ii_103[k];

        t_48[k] = -gi_20[k]
                  + f_0 * ii_104[k];

        t_49[k] = -gi_21[k]
                  + f_0 * ii_105[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gi_22, gi_23, gi_24, gi_25, gi_26, \
                         ii_106, ii_107, ii_108, ii_109, ii_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -gi_22[k]
                  + f_0 * ii_106[k];

        t_51[k] = -gi_23[k]
                  + f_0 * ii_107[k];

        t_52[k] = -gi_24[k]
                  + f_0 * ii_108[k];

        t_53[k] = -gi_25[k]
                  + f_0 * ii_109[k];

        t_54[k] = -gi_26[k]
                  + f_0 * ii_110[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, gi_27, ii_111, ii_112, \
                         ii_113, ii_114, ii_115, ii_116, ii_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -gi_27[k]
                  + f_0 * ii_111[k];

        t_56[k] = f_0 * ii_112[k];

        t_57[k] = f_0 * ii_113[k];

        t_58[k] = f_0 * ii_114[k];

        t_59[k] = f_0 * ii_115[k];

        t_60[k] = f_0 * ii_116[k];

        t_61[k] = f_0 * ii_117[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, t_68, t_69, ii_118, ii_119, \
                         ii_120, ii_121, ii_122, ii_123, ii_124, \
                         ii_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * ii_118[k];

        t_63[k] = f_0 * ii_119[k];

        t_64[k] = f_0 * ii_120[k];

        t_65[k] = f_0 * ii_121[k];

        t_66[k] = f_0 * ii_122[k];

        t_67[k] = f_0 * ii_123[k];

        t_68[k] = f_0 * ii_124[k];

        t_69[k] = f_0 * ii_125[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, t_76, t_77, ii_126, ii_127, \
                         ii_128, ii_129, ii_130, ii_131, ii_132, \
                         ii_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * ii_126[k];

        t_71[k] = f_0 * ii_127[k];

        t_72[k] = f_0 * ii_128[k];

        t_73[k] = f_0 * ii_129[k];

        t_74[k] = f_0 * ii_130[k];

        t_75[k] = f_0 * ii_131[k];

        t_76[k] = f_0 * ii_132[k];

        t_77[k] = f_0 * ii_133[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, t_84, gi_28, ii_134, ii_135, \
                         ii_136, ii_137, ii_138, ii_139, ii_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * ii_134[k];

        t_79[k] = f_0 * ii_135[k];

        t_80[k] = f_0 * ii_136[k];

        t_81[k] = f_0 * ii_137[k];

        t_82[k] = f_0 * ii_138[k];

        t_83[k] = f_0 * ii_139[k];

        t_84[k] = -2.0 * gi_28[k]
                  + f_0 * ii_168[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gi_29, gi_30, gi_31, gi_32, gi_33, \
                         ii_169, ii_170, ii_171, ii_172, ii_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * gi_29[k]
                  + f_0 * ii_169[k];

        t_86[k] = -2.0 * gi_30[k]
                  + f_0 * ii_170[k];

        t_87[k] = -2.0 * gi_31[k]
                  + f_0 * ii_171[k];

        t_88[k] = -2.0 * gi_32[k]
                  + f_0 * ii_172[k];

        t_89[k] = -2.0 * gi_33[k]
                  + f_0 * ii_173[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, gi_34, gi_35, gi_36, gi_37, gi_38, \
                         ii_174, ii_175, ii_176, ii_177, ii_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * gi_34[k]
                  + f_0 * ii_174[k];

        t_91[k] = -2.0 * gi_35[k]
                  + f_0 * ii_175[k];

        t_92[k] = -2.0 * gi_36[k]
                  + f_0 * ii_176[k];

        t_93[k] = -2.0 * gi_37[k]
                  + f_0 * ii_177[k];

        t_94[k] = -2.0 * gi_38[k]
                  + f_0 * ii_178[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, gi_39, gi_40, gi_41, gi_42, gi_43, \
                         ii_179, ii_180, ii_181, ii_182, ii_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * gi_39[k]
                  + f_0 * ii_179[k];

        t_96[k] = -2.0 * gi_40[k]
                  + f_0 * ii_180[k];

        t_97[k] = -2.0 * gi_41[k]
                  + f_0 * ii_181[k];

        t_98[k] = -2.0 * gi_42[k]
                  + f_0 * ii_182[k];

        t_99[k] = -2.0 * gi_43[k]
                  + f_0 * ii_183[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, gi_44, gi_45, gi_46, gi_47, gi_48, \
                         ii_184, ii_185, ii_186, ii_187, ii_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * gi_44[k]
                   + f_0 * ii_184[k];

        t_101[k] = -2.0 * gi_45[k]
                   + f_0 * ii_185[k];

        t_102[k] = -2.0 * gi_46[k]
                   + f_0 * ii_186[k];

        t_103[k] = -2.0 * gi_47[k]
                   + f_0 * ii_187[k];

        t_104[k] = -2.0 * gi_48[k]
                   + f_0 * ii_188[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, gi_49, gi_50, gi_51, gi_52, gi_53, \
                         ii_189, ii_190, ii_191, ii_192, ii_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * gi_49[k]
                   + f_0 * ii_189[k];

        t_106[k] = -2.0 * gi_50[k]
                   + f_0 * ii_190[k];

        t_107[k] = -2.0 * gi_51[k]
                   + f_0 * ii_191[k];

        t_108[k] = -2.0 * gi_52[k]
                   + f_0 * ii_192[k];

        t_109[k] = -2.0 * gi_53[k]
                   + f_0 * ii_193[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, gi_54, gi_55, gi_56, gi_57, gi_58, \
                         ii_194, ii_195, ii_196, ii_197, ii_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * gi_54[k]
                   + f_0 * ii_194[k];

        t_111[k] = -2.0 * gi_55[k]
                   + f_0 * ii_195[k];

        t_112[k] = -gi_56[k]
                   + f_0 * ii_196[k];

        t_113[k] = -gi_57[k]
                   + f_0 * ii_197[k];

        t_114[k] = -gi_58[k]
                   + f_0 * ii_198[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, gi_59, gi_60, gi_61, gi_62, gi_63, \
                         ii_199, ii_200, ii_201, ii_202, ii_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -gi_59[k]
                   + f_0 * ii_199[k];

        t_116[k] = -gi_60[k]
                   + f_0 * ii_200[k];

        t_117[k] = -gi_61[k]
                   + f_0 * ii_201[k];

        t_118[k] = -gi_62[k]
                   + f_0 * ii_202[k];

        t_119[k] = -gi_63[k]
                   + f_0 * ii_203[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, gi_64, gi_65, gi_66, gi_67, gi_68, \
                         ii_204, ii_205, ii_206, ii_207, ii_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -gi_64[k]
                   + f_0 * ii_204[k];

        t_121[k] = -gi_65[k]
                   + f_0 * ii_205[k];

        t_122[k] = -gi_66[k]
                   + f_0 * ii_206[k];

        t_123[k] = -gi_67[k]
                   + f_0 * ii_207[k];

        t_124[k] = -gi_68[k]
                   + f_0 * ii_208[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, gi_69, gi_70, gi_71, gi_72, gi_73, \
                         ii_209, ii_210, ii_211, ii_212, ii_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -gi_69[k]
                   + f_0 * ii_209[k];

        t_126[k] = -gi_70[k]
                   + f_0 * ii_210[k];

        t_127[k] = -gi_71[k]
                   + f_0 * ii_211[k];

        t_128[k] = -gi_72[k]
                   + f_0 * ii_212[k];

        t_129[k] = -gi_73[k]
                   + f_0 * ii_213[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, gi_74, gi_75, gi_76, gi_77, gi_78, \
                         ii_214, ii_215, ii_216, ii_217, ii_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -gi_74[k]
                   + f_0 * ii_214[k];

        t_131[k] = -gi_75[k]
                   + f_0 * ii_215[k];

        t_132[k] = -gi_76[k]
                   + f_0 * ii_216[k];

        t_133[k] = -gi_77[k]
                   + f_0 * ii_217[k];

        t_134[k] = -gi_78[k]
                   + f_0 * ii_218[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, gi_79, gi_80, gi_81, gi_82, gi_83, \
                         ii_219, ii_220, ii_221, ii_222, ii_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -gi_79[k]
                   + f_0 * ii_219[k];

        t_136[k] = -gi_80[k]
                   + f_0 * ii_220[k];

        t_137[k] = -gi_81[k]
                   + f_0 * ii_221[k];

        t_138[k] = -gi_82[k]
                   + f_0 * ii_222[k];

        t_139[k] = -gi_83[k]
                   + f_0 * ii_223[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, ii_224, \
                         ii_225, ii_226, ii_227, ii_228, ii_229, ii_230, \
                         ii_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * ii_224[k];

        t_141[k] = f_0 * ii_225[k];

        t_142[k] = f_0 * ii_226[k];

        t_143[k] = f_0 * ii_227[k];

        t_144[k] = f_0 * ii_228[k];

        t_145[k] = f_0 * ii_229[k];

        t_146[k] = f_0 * ii_230[k];

        t_147[k] = f_0 * ii_231[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, t_154, t_155, ii_232, \
                         ii_233, ii_234, ii_235, ii_236, ii_237, ii_238, \
                         ii_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * ii_232[k];

        t_149[k] = f_0 * ii_233[k];

        t_150[k] = f_0 * ii_234[k];

        t_151[k] = f_0 * ii_235[k];

        t_152[k] = f_0 * ii_236[k];

        t_153[k] = f_0 * ii_237[k];

        t_154[k] = f_0 * ii_238[k];

        t_155[k] = f_0 * ii_239[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, t_162, t_163, ii_240, \
                         ii_241, ii_242, ii_243, ii_244, ii_245, ii_246, \
                         ii_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * ii_240[k];

        t_157[k] = f_0 * ii_241[k];

        t_158[k] = f_0 * ii_242[k];

        t_159[k] = f_0 * ii_243[k];

        t_160[k] = f_0 * ii_244[k];

        t_161[k] = f_0 * ii_245[k];

        t_162[k] = f_0 * ii_246[k];

        t_163[k] = f_0 * ii_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, gi_84, gi_85, ii_248, \
                         ii_249, ii_250, ii_251, ii_280, ii_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * ii_248[k];

        t_165[k] = f_0 * ii_249[k];

        t_166[k] = f_0 * ii_250[k];

        t_167[k] = f_0 * ii_251[k];

        t_168[k] = -3.0 * gi_84[k]
                   + f_0 * ii_280[k];

        t_169[k] = -3.0 * gi_85[k]
                   + f_0 * ii_281[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, gi_86, gi_87, gi_88, gi_89, gi_90, \
                         ii_282, ii_283, ii_284, ii_285, ii_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * gi_86[k]
                   + f_0 * ii_282[k];

        t_171[k] = -3.0 * gi_87[k]
                   + f_0 * ii_283[k];

        t_172[k] = -3.0 * gi_88[k]
                   + f_0 * ii_284[k];

        t_173[k] = -3.0 * gi_89[k]
                   + f_0 * ii_285[k];

        t_174[k] = -3.0 * gi_90[k]
                   + f_0 * ii_286[k];
    }
}

static auto
compute_prim_geom_10_hi_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_92 = buffer.data(gi + 92);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_95 = buffer.data(gi + 95);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_97 = buffer.data(gi + 97);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_100 = buffer.data(gi + 100);
    const auto *gi_101 = buffer.data(gi + 101);
    const auto *gi_102 = buffer.data(gi + 102);
    const auto *gi_103 = buffer.data(gi + 103);
    const auto *gi_104 = buffer.data(gi + 104);
    const auto *gi_105 = buffer.data(gi + 105);
    const auto *gi_106 = buffer.data(gi + 106);
    const auto *gi_107 = buffer.data(gi + 107);
    const auto *gi_108 = buffer.data(gi + 108);
    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_110 = buffer.data(gi + 110);
    const auto *gi_111 = buffer.data(gi + 111);
    const auto *gi_112 = buffer.data(gi + 112);
    const auto *gi_113 = buffer.data(gi + 113);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_116 = buffer.data(gi + 116);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_119 = buffer.data(gi + 119);
    const auto *gi_120 = buffer.data(gi + 120);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_123 = buffer.data(gi + 123);
    const auto *gi_124 = buffer.data(gi + 124);
    const auto *gi_125 = buffer.data(gi + 125);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_127 = buffer.data(gi + 127);
    const auto *gi_128 = buffer.data(gi + 128);
    const auto *gi_129 = buffer.data(gi + 129);
    const auto *gi_130 = buffer.data(gi + 130);
    const auto *gi_131 = buffer.data(gi + 131);
    const auto *gi_132 = buffer.data(gi + 132);
    const auto *gi_133 = buffer.data(gi + 133);
    const auto *gi_134 = buffer.data(gi + 134);
    const auto *gi_135 = buffer.data(gi + 135);
    const auto *gi_136 = buffer.data(gi + 136);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_138 = buffer.data(gi + 138);
    const auto *gi_139 = buffer.data(gi + 139);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_141 = buffer.data(gi + 141);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_144 = buffer.data(gi + 144);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_147 = buffer.data(gi + 147);
    const auto *gi_148 = buffer.data(gi + 148);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_151 = buffer.data(gi + 151);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);
    const auto *gi_158 = buffer.data(gi + 158);
    const auto *gi_159 = buffer.data(gi + 159);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_161 = buffer.data(gi + 161);
    const auto *gi_162 = buffer.data(gi + 162);
    const auto *gi_163 = buffer.data(gi + 163);
    const auto *gi_164 = buffer.data(gi + 164);
    const auto *gi_165 = buffer.data(gi + 165);
    const auto *gi_166 = buffer.data(gi + 166);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_170 = buffer.data(gi + 170);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_172 = buffer.data(gi + 172);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_175 = buffer.data(gi + 175);
    const auto *gi_176 = buffer.data(gi + 176);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_179 = buffer.data(gi + 179);
    const auto *gi_180 = buffer.data(gi + 180);
    const auto *gi_181 = buffer.data(gi + 181);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_184 = buffer.data(gi + 184);
    const auto *gi_185 = buffer.data(gi + 185);
    const auto *gi_186 = buffer.data(gi + 186);
    const auto *gi_187 = buffer.data(gi + 187);
    const auto *gi_188 = buffer.data(gi + 188);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_190 = buffer.data(gi + 190);
    const auto *gi_191 = buffer.data(gi + 191);
    const auto *gi_192 = buffer.data(gi + 192);
    const auto *gi_193 = buffer.data(gi + 193);
    const auto *gi_194 = buffer.data(gi + 194);
    const auto *gi_195 = buffer.data(gi + 195);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_197 = buffer.data(gi + 197);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_200 = buffer.data(gi + 200);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_203 = buffer.data(gi + 203);
    const auto *gi_204 = buffer.data(gi + 204);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_207 = buffer.data(gi + 207);
    const auto *gi_208 = buffer.data(gi + 208);
    const auto *gi_209 = buffer.data(gi + 209);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_211 = buffer.data(gi + 211);
    const auto *gi_212 = buffer.data(gi + 212);
    const auto *gi_213 = buffer.data(gi + 213);
    const auto *gi_214 = buffer.data(gi + 214);
    const auto *gi_215 = buffer.data(gi + 215);
    const auto *gi_216 = buffer.data(gi + 216);
    const auto *gi_217 = buffer.data(gi + 217);
    const auto *gi_218 = buffer.data(gi + 218);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_222 = buffer.data(gi + 222);

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

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, gi_91, gi_92, gi_93, gi_94, gi_95, \
                         ii_287, ii_288, ii_289, ii_290, ii_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * gi_91[k]
                   + f_0 * ii_287[k];

        t_176[k] = -3.0 * gi_92[k]
                   + f_0 * ii_288[k];

        t_177[k] = -3.0 * gi_93[k]
                   + f_0 * ii_289[k];

        t_178[k] = -3.0 * gi_94[k]
                   + f_0 * ii_290[k];

        t_179[k] = -3.0 * gi_95[k]
                   + f_0 * ii_291[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, gi_96, gi_97, gi_98, gi_99, \
                         gi_100, ii_292, ii_293, ii_294, ii_295, \
                         ii_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * gi_96[k]
                   + f_0 * ii_292[k];

        t_181[k] = -3.0 * gi_97[k]
                   + f_0 * ii_293[k];

        t_182[k] = -3.0 * gi_98[k]
                   + f_0 * ii_294[k];

        t_183[k] = -3.0 * gi_99[k]
                   + f_0 * ii_295[k];

        t_184[k] = -3.0 * gi_100[k]
                   + f_0 * ii_296[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, gi_101, gi_102, gi_103, gi_104, \
                         gi_105, ii_297, ii_298, ii_299, ii_300, \
                         ii_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * gi_101[k]
                   + f_0 * ii_297[k];

        t_186[k] = -3.0 * gi_102[k]
                   + f_0 * ii_298[k];

        t_187[k] = -3.0 * gi_103[k]
                   + f_0 * ii_299[k];

        t_188[k] = -3.0 * gi_104[k]
                   + f_0 * ii_300[k];

        t_189[k] = -3.0 * gi_105[k]
                   + f_0 * ii_301[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, gi_106, gi_107, gi_108, gi_109, \
                         gi_110, ii_302, ii_303, ii_304, ii_305, \
                         ii_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * gi_106[k]
                   + f_0 * ii_302[k];

        t_191[k] = -3.0 * gi_107[k]
                   + f_0 * ii_303[k];

        t_192[k] = -3.0 * gi_108[k]
                   + f_0 * ii_304[k];

        t_193[k] = -3.0 * gi_109[k]
                   + f_0 * ii_305[k];

        t_194[k] = -3.0 * gi_110[k]
                   + f_0 * ii_306[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, gi_111, gi_112, gi_113, gi_114, \
                         gi_115, ii_307, ii_308, ii_309, ii_310, \
                         ii_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * gi_111[k]
                   + f_0 * ii_307[k];

        t_196[k] = -2.0 * gi_112[k]
                   + f_0 * ii_308[k];

        t_197[k] = -2.0 * gi_113[k]
                   + f_0 * ii_309[k];

        t_198[k] = -2.0 * gi_114[k]
                   + f_0 * ii_310[k];

        t_199[k] = -2.0 * gi_115[k]
                   + f_0 * ii_311[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, gi_116, gi_117, gi_118, gi_119, \
                         gi_120, ii_312, ii_313, ii_314, ii_315, \
                         ii_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * gi_116[k]
                   + f_0 * ii_312[k];

        t_201[k] = -2.0 * gi_117[k]
                   + f_0 * ii_313[k];

        t_202[k] = -2.0 * gi_118[k]
                   + f_0 * ii_314[k];

        t_203[k] = -2.0 * gi_119[k]
                   + f_0 * ii_315[k];

        t_204[k] = -2.0 * gi_120[k]
                   + f_0 * ii_316[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, gi_121, gi_122, gi_123, gi_124, \
                         gi_125, ii_317, ii_318, ii_319, ii_320, \
                         ii_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * gi_121[k]
                   + f_0 * ii_317[k];

        t_206[k] = -2.0 * gi_122[k]
                   + f_0 * ii_318[k];

        t_207[k] = -2.0 * gi_123[k]
                   + f_0 * ii_319[k];

        t_208[k] = -2.0 * gi_124[k]
                   + f_0 * ii_320[k];

        t_209[k] = -2.0 * gi_125[k]
                   + f_0 * ii_321[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, gi_126, gi_127, gi_128, gi_129, \
                         gi_130, ii_322, ii_323, ii_324, ii_325, \
                         ii_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * gi_126[k]
                   + f_0 * ii_322[k];

        t_211[k] = -2.0 * gi_127[k]
                   + f_0 * ii_323[k];

        t_212[k] = -2.0 * gi_128[k]
                   + f_0 * ii_324[k];

        t_213[k] = -2.0 * gi_129[k]
                   + f_0 * ii_325[k];

        t_214[k] = -2.0 * gi_130[k]
                   + f_0 * ii_326[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, gi_131, gi_132, gi_133, gi_134, \
                         gi_135, ii_327, ii_328, ii_329, ii_330, \
                         ii_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * gi_131[k]
                   + f_0 * ii_327[k];

        t_216[k] = -2.0 * gi_132[k]
                   + f_0 * ii_328[k];

        t_217[k] = -2.0 * gi_133[k]
                   + f_0 * ii_329[k];

        t_218[k] = -2.0 * gi_134[k]
                   + f_0 * ii_330[k];

        t_219[k] = -2.0 * gi_135[k]
                   + f_0 * ii_331[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, gi_136, gi_137, gi_138, gi_139, \
                         gi_140, ii_332, ii_333, ii_334, ii_335, \
                         ii_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * gi_136[k]
                   + f_0 * ii_332[k];

        t_221[k] = -2.0 * gi_137[k]
                   + f_0 * ii_333[k];

        t_222[k] = -2.0 * gi_138[k]
                   + f_0 * ii_334[k];

        t_223[k] = -2.0 * gi_139[k]
                   + f_0 * ii_335[k];

        t_224[k] = -gi_140[k]
                   + f_0 * ii_336[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, gi_141, gi_142, gi_143, gi_144, \
                         gi_145, ii_337, ii_338, ii_339, ii_340, \
                         ii_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -gi_141[k]
                   + f_0 * ii_337[k];

        t_226[k] = -gi_142[k]
                   + f_0 * ii_338[k];

        t_227[k] = -gi_143[k]
                   + f_0 * ii_339[k];

        t_228[k] = -gi_144[k]
                   + f_0 * ii_340[k];

        t_229[k] = -gi_145[k]
                   + f_0 * ii_341[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, gi_146, gi_147, gi_148, gi_149, \
                         gi_150, ii_342, ii_343, ii_344, ii_345, \
                         ii_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -gi_146[k]
                   + f_0 * ii_342[k];

        t_231[k] = -gi_147[k]
                   + f_0 * ii_343[k];

        t_232[k] = -gi_148[k]
                   + f_0 * ii_344[k];

        t_233[k] = -gi_149[k]
                   + f_0 * ii_345[k];

        t_234[k] = -gi_150[k]
                   + f_0 * ii_346[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, gi_151, gi_152, gi_153, gi_154, \
                         gi_155, ii_347, ii_348, ii_349, ii_350, \
                         ii_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -gi_151[k]
                   + f_0 * ii_347[k];

        t_236[k] = -gi_152[k]
                   + f_0 * ii_348[k];

        t_237[k] = -gi_153[k]
                   + f_0 * ii_349[k];

        t_238[k] = -gi_154[k]
                   + f_0 * ii_350[k];

        t_239[k] = -gi_155[k]
                   + f_0 * ii_351[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, gi_156, gi_157, gi_158, gi_159, \
                         gi_160, ii_352, ii_353, ii_354, ii_355, \
                         ii_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -gi_156[k]
                   + f_0 * ii_352[k];

        t_241[k] = -gi_157[k]
                   + f_0 * ii_353[k];

        t_242[k] = -gi_158[k]
                   + f_0 * ii_354[k];

        t_243[k] = -gi_159[k]
                   + f_0 * ii_355[k];

        t_244[k] = -gi_160[k]
                   + f_0 * ii_356[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, gi_161, gi_162, gi_163, gi_164, \
                         gi_165, ii_357, ii_358, ii_359, ii_360, \
                         ii_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -gi_161[k]
                   + f_0 * ii_357[k];

        t_246[k] = -gi_162[k]
                   + f_0 * ii_358[k];

        t_247[k] = -gi_163[k]
                   + f_0 * ii_359[k];

        t_248[k] = -gi_164[k]
                   + f_0 * ii_360[k];

        t_249[k] = -gi_165[k]
                   + f_0 * ii_361[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, t_255, t_256, gi_166, gi_167, \
                         ii_362, ii_363, ii_364, ii_365, ii_366, ii_367, \
                         ii_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -gi_166[k]
                   + f_0 * ii_362[k];

        t_251[k] = -gi_167[k]
                   + f_0 * ii_363[k];

        t_252[k] = f_0 * ii_364[k];

        t_253[k] = f_0 * ii_365[k];

        t_254[k] = f_0 * ii_366[k];

        t_255[k] = f_0 * ii_367[k];

        t_256[k] = f_0 * ii_368[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, t_262, t_263, t_264, ii_369, \
                         ii_370, ii_371, ii_372, ii_373, ii_374, ii_375, \
                         ii_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_0 * ii_369[k];

        t_258[k] = f_0 * ii_370[k];

        t_259[k] = f_0 * ii_371[k];

        t_260[k] = f_0 * ii_372[k];

        t_261[k] = f_0 * ii_373[k];

        t_262[k] = f_0 * ii_374[k];

        t_263[k] = f_0 * ii_375[k];

        t_264[k] = f_0 * ii_376[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, t_271, t_272, ii_377, \
                         ii_378, ii_379, ii_380, ii_381, ii_382, ii_383, \
                         ii_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_0 * ii_377[k];

        t_266[k] = f_0 * ii_378[k];

        t_267[k] = f_0 * ii_379[k];

        t_268[k] = f_0 * ii_380[k];

        t_269[k] = f_0 * ii_381[k];

        t_270[k] = f_0 * ii_382[k];

        t_271[k] = f_0 * ii_383[k];

        t_272[k] = f_0 * ii_384[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, t_278, t_279, ii_385, ii_386, \
                         ii_387, ii_388, ii_389, ii_390, ii_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_0 * ii_385[k];

        t_274[k] = f_0 * ii_386[k];

        t_275[k] = f_0 * ii_387[k];

        t_276[k] = f_0 * ii_388[k];

        t_277[k] = f_0 * ii_389[k];

        t_278[k] = f_0 * ii_390[k];

        t_279[k] = f_0 * ii_391[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, gi_168, gi_169, gi_170, gi_171, \
                         gi_172, ii_420, ii_421, ii_422, ii_423, \
                         ii_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = -4.0 * gi_168[k]
                   + f_0 * ii_420[k];

        t_281[k] = -4.0 * gi_169[k]
                   + f_0 * ii_421[k];

        t_282[k] = -4.0 * gi_170[k]
                   + f_0 * ii_422[k];

        t_283[k] = -4.0 * gi_171[k]
                   + f_0 * ii_423[k];

        t_284[k] = -4.0 * gi_172[k]
                   + f_0 * ii_424[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, gi_173, gi_174, gi_175, gi_176, \
                         gi_177, ii_425, ii_426, ii_427, ii_428, \
                         ii_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = -4.0 * gi_173[k]
                   + f_0 * ii_425[k];

        t_286[k] = -4.0 * gi_174[k]
                   + f_0 * ii_426[k];

        t_287[k] = -4.0 * gi_175[k]
                   + f_0 * ii_427[k];

        t_288[k] = -4.0 * gi_176[k]
                   + f_0 * ii_428[k];

        t_289[k] = -4.0 * gi_177[k]
                   + f_0 * ii_429[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, gi_178, gi_179, gi_180, gi_181, \
                         gi_182, ii_430, ii_431, ii_432, ii_433, \
                         ii_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = -4.0 * gi_178[k]
                   + f_0 * ii_430[k];

        t_291[k] = -4.0 * gi_179[k]
                   + f_0 * ii_431[k];

        t_292[k] = -4.0 * gi_180[k]
                   + f_0 * ii_432[k];

        t_293[k] = -4.0 * gi_181[k]
                   + f_0 * ii_433[k];

        t_294[k] = -4.0 * gi_182[k]
                   + f_0 * ii_434[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, gi_183, gi_184, gi_185, gi_186, \
                         gi_187, ii_435, ii_436, ii_437, ii_438, \
                         ii_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = -4.0 * gi_183[k]
                   + f_0 * ii_435[k];

        t_296[k] = -4.0 * gi_184[k]
                   + f_0 * ii_436[k];

        t_297[k] = -4.0 * gi_185[k]
                   + f_0 * ii_437[k];

        t_298[k] = -4.0 * gi_186[k]
                   + f_0 * ii_438[k];

        t_299[k] = -4.0 * gi_187[k]
                   + f_0 * ii_439[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, gi_188, gi_189, gi_190, gi_191, \
                         gi_192, ii_440, ii_441, ii_442, ii_443, \
                         ii_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -4.0 * gi_188[k]
                   + f_0 * ii_440[k];

        t_301[k] = -4.0 * gi_189[k]
                   + f_0 * ii_441[k];

        t_302[k] = -4.0 * gi_190[k]
                   + f_0 * ii_442[k];

        t_303[k] = -4.0 * gi_191[k]
                   + f_0 * ii_443[k];

        t_304[k] = -4.0 * gi_192[k]
                   + f_0 * ii_444[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, gi_193, gi_194, gi_195, gi_196, \
                         gi_197, ii_445, ii_446, ii_447, ii_448, \
                         ii_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -4.0 * gi_193[k]
                   + f_0 * ii_445[k];

        t_306[k] = -4.0 * gi_194[k]
                   + f_0 * ii_446[k];

        t_307[k] = -4.0 * gi_195[k]
                   + f_0 * ii_447[k];

        t_308[k] = -3.0 * gi_196[k]
                   + f_0 * ii_448[k];

        t_309[k] = -3.0 * gi_197[k]
                   + f_0 * ii_449[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, gi_198, gi_199, gi_200, gi_201, \
                         gi_202, ii_450, ii_451, ii_452, ii_453, \
                         ii_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -3.0 * gi_198[k]
                   + f_0 * ii_450[k];

        t_311[k] = -3.0 * gi_199[k]
                   + f_0 * ii_451[k];

        t_312[k] = -3.0 * gi_200[k]
                   + f_0 * ii_452[k];

        t_313[k] = -3.0 * gi_201[k]
                   + f_0 * ii_453[k];

        t_314[k] = -3.0 * gi_202[k]
                   + f_0 * ii_454[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, gi_203, gi_204, gi_205, gi_206, \
                         gi_207, ii_455, ii_456, ii_457, ii_458, \
                         ii_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -3.0 * gi_203[k]
                   + f_0 * ii_455[k];

        t_316[k] = -3.0 * gi_204[k]
                   + f_0 * ii_456[k];

        t_317[k] = -3.0 * gi_205[k]
                   + f_0 * ii_457[k];

        t_318[k] = -3.0 * gi_206[k]
                   + f_0 * ii_458[k];

        t_319[k] = -3.0 * gi_207[k]
                   + f_0 * ii_459[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, gi_208, gi_209, gi_210, gi_211, \
                         gi_212, ii_460, ii_461, ii_462, ii_463, \
                         ii_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = -3.0 * gi_208[k]
                   + f_0 * ii_460[k];

        t_321[k] = -3.0 * gi_209[k]
                   + f_0 * ii_461[k];

        t_322[k] = -3.0 * gi_210[k]
                   + f_0 * ii_462[k];

        t_323[k] = -3.0 * gi_211[k]
                   + f_0 * ii_463[k];

        t_324[k] = -3.0 * gi_212[k]
                   + f_0 * ii_464[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, gi_213, gi_214, gi_215, gi_216, \
                         gi_217, ii_465, ii_466, ii_467, ii_468, \
                         ii_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = -3.0 * gi_213[k]
                   + f_0 * ii_465[k];

        t_326[k] = -3.0 * gi_214[k]
                   + f_0 * ii_466[k];

        t_327[k] = -3.0 * gi_215[k]
                   + f_0 * ii_467[k];

        t_328[k] = -3.0 * gi_216[k]
                   + f_0 * ii_468[k];

        t_329[k] = -3.0 * gi_217[k]
                   + f_0 * ii_469[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, gi_218, gi_219, gi_220, gi_221, \
                         gi_222, ii_470, ii_471, ii_472, ii_473, \
                         ii_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = -3.0 * gi_218[k]
                   + f_0 * ii_470[k];

        t_331[k] = -3.0 * gi_219[k]
                   + f_0 * ii_471[k];

        t_332[k] = -3.0 * gi_220[k]
                   + f_0 * ii_472[k];

        t_333[k] = -3.0 * gi_221[k]
                   + f_0 * ii_473[k];

        t_334[k] = -3.0 * gi_222[k]
                   + f_0 * ii_474[k];
    }
}

static auto
compute_prim_geom_10_hi_electron_repulsion_1_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_223 = buffer.data(gi + 223);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_225 = buffer.data(gi + 225);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_228 = buffer.data(gi + 228);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_231 = buffer.data(gi + 231);
    const auto *gi_232 = buffer.data(gi + 232);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_235 = buffer.data(gi + 235);
    const auto *gi_236 = buffer.data(gi + 236);
    const auto *gi_237 = buffer.data(gi + 237);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_239 = buffer.data(gi + 239);
    const auto *gi_240 = buffer.data(gi + 240);
    const auto *gi_241 = buffer.data(gi + 241);
    const auto *gi_242 = buffer.data(gi + 242);
    const auto *gi_243 = buffer.data(gi + 243);
    const auto *gi_244 = buffer.data(gi + 244);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_246 = buffer.data(gi + 246);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_250 = buffer.data(gi + 250);
    const auto *gi_251 = buffer.data(gi + 251);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_253 = buffer.data(gi + 253);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_256 = buffer.data(gi + 256);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_258 = buffer.data(gi + 258);
    const auto *gi_259 = buffer.data(gi + 259);
    const auto *gi_260 = buffer.data(gi + 260);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_263 = buffer.data(gi + 263);
    const auto *gi_264 = buffer.data(gi + 264);
    const auto *gi_265 = buffer.data(gi + 265);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_267 = buffer.data(gi + 267);
    const auto *gi_268 = buffer.data(gi + 268);
    const auto *gi_269 = buffer.data(gi + 269);
    const auto *gi_270 = buffer.data(gi + 270);
    const auto *gi_271 = buffer.data(gi + 271);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_273 = buffer.data(gi + 273);
    const auto *gi_274 = buffer.data(gi + 274);
    const auto *gi_275 = buffer.data(gi + 275);
    const auto *gi_276 = buffer.data(gi + 276);
    const auto *gi_277 = buffer.data(gi + 277);
    const auto *gi_278 = buffer.data(gi + 278);
    const auto *gi_279 = buffer.data(gi + 279);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_284 = buffer.data(gi + 284);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_287 = buffer.data(gi + 287);
    const auto *gi_288 = buffer.data(gi + 288);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_291 = buffer.data(gi + 291);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_293 = buffer.data(gi + 293);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_296 = buffer.data(gi + 296);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_299 = buffer.data(gi + 299);
    const auto *gi_300 = buffer.data(gi + 300);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_307 = buffer.data(gi + 307);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_309 = buffer.data(gi + 309);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_312 = buffer.data(gi + 312);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_315 = buffer.data(gi + 315);
    const auto *gi_316 = buffer.data(gi + 316);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_319 = buffer.data(gi + 319);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_321 = buffer.data(gi + 321);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_323 = buffer.data(gi + 323);
    const auto *gi_324 = buffer.data(gi + 324);
    const auto *gi_325 = buffer.data(gi + 325);
    const auto *gi_326 = buffer.data(gi + 326);
    const auto *gi_327 = buffer.data(gi + 327);
    const auto *gi_328 = buffer.data(gi + 328);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_335 = buffer.data(gi + 335);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_337 = buffer.data(gi + 337);
    const auto *gi_338 = buffer.data(gi + 338);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_340 = buffer.data(gi + 340);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_343 = buffer.data(gi + 343);
    const auto *gi_344 = buffer.data(gi + 344);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_347 = buffer.data(gi + 347);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_349 = buffer.data(gi + 349);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_352 = buffer.data(gi + 352);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);

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

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, gi_223, gi_224, gi_225, gi_226, \
                         gi_227, ii_475, ii_476, ii_477, ii_478, \
                         ii_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = -3.0 * gi_223[k]
                   + f_0 * ii_475[k];

        t_336[k] = -2.0 * gi_224[k]
                   + f_0 * ii_476[k];

        t_337[k] = -2.0 * gi_225[k]
                   + f_0 * ii_477[k];

        t_338[k] = -2.0 * gi_226[k]
                   + f_0 * ii_478[k];

        t_339[k] = -2.0 * gi_227[k]
                   + f_0 * ii_479[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, gi_228, gi_229, gi_230, gi_231, \
                         gi_232, ii_480, ii_481, ii_482, ii_483, \
                         ii_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = -2.0 * gi_228[k]
                   + f_0 * ii_480[k];

        t_341[k] = -2.0 * gi_229[k]
                   + f_0 * ii_481[k];

        t_342[k] = -2.0 * gi_230[k]
                   + f_0 * ii_482[k];

        t_343[k] = -2.0 * gi_231[k]
                   + f_0 * ii_483[k];

        t_344[k] = -2.0 * gi_232[k]
                   + f_0 * ii_484[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, gi_233, gi_234, gi_235, gi_236, \
                         gi_237, ii_485, ii_486, ii_487, ii_488, \
                         ii_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = -2.0 * gi_233[k]
                   + f_0 * ii_485[k];

        t_346[k] = -2.0 * gi_234[k]
                   + f_0 * ii_486[k];

        t_347[k] = -2.0 * gi_235[k]
                   + f_0 * ii_487[k];

        t_348[k] = -2.0 * gi_236[k]
                   + f_0 * ii_488[k];

        t_349[k] = -2.0 * gi_237[k]
                   + f_0 * ii_489[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, gi_238, gi_239, gi_240, gi_241, \
                         gi_242, ii_490, ii_491, ii_492, ii_493, \
                         ii_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = -2.0 * gi_238[k]
                   + f_0 * ii_490[k];

        t_351[k] = -2.0 * gi_239[k]
                   + f_0 * ii_491[k];

        t_352[k] = -2.0 * gi_240[k]
                   + f_0 * ii_492[k];

        t_353[k] = -2.0 * gi_241[k]
                   + f_0 * ii_493[k];

        t_354[k] = -2.0 * gi_242[k]
                   + f_0 * ii_494[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, gi_243, gi_244, gi_245, gi_246, \
                         gi_247, ii_495, ii_496, ii_497, ii_498, \
                         ii_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = -2.0 * gi_243[k]
                   + f_0 * ii_495[k];

        t_356[k] = -2.0 * gi_244[k]
                   + f_0 * ii_496[k];

        t_357[k] = -2.0 * gi_245[k]
                   + f_0 * ii_497[k];

        t_358[k] = -2.0 * gi_246[k]
                   + f_0 * ii_498[k];

        t_359[k] = -2.0 * gi_247[k]
                   + f_0 * ii_499[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, gi_248, gi_249, gi_250, gi_251, \
                         gi_252, ii_500, ii_501, ii_502, ii_503, \
                         ii_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -2.0 * gi_248[k]
                   + f_0 * ii_500[k];

        t_361[k] = -2.0 * gi_249[k]
                   + f_0 * ii_501[k];

        t_362[k] = -2.0 * gi_250[k]
                   + f_0 * ii_502[k];

        t_363[k] = -2.0 * gi_251[k]
                   + f_0 * ii_503[k];

        t_364[k] = -gi_252[k]
                   + f_0 * ii_504[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, gi_253, gi_254, gi_255, gi_256, \
                         gi_257, ii_505, ii_506, ii_507, ii_508, \
                         ii_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = -gi_253[k]
                   + f_0 * ii_505[k];

        t_366[k] = -gi_254[k]
                   + f_0 * ii_506[k];

        t_367[k] = -gi_255[k]
                   + f_0 * ii_507[k];

        t_368[k] = -gi_256[k]
                   + f_0 * ii_508[k];

        t_369[k] = -gi_257[k]
                   + f_0 * ii_509[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, gi_258, gi_259, gi_260, gi_261, \
                         gi_262, ii_510, ii_511, ii_512, ii_513, \
                         ii_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = -gi_258[k]
                   + f_0 * ii_510[k];

        t_371[k] = -gi_259[k]
                   + f_0 * ii_511[k];

        t_372[k] = -gi_260[k]
                   + f_0 * ii_512[k];

        t_373[k] = -gi_261[k]
                   + f_0 * ii_513[k];

        t_374[k] = -gi_262[k]
                   + f_0 * ii_514[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, gi_263, gi_264, gi_265, gi_266, \
                         gi_267, ii_515, ii_516, ii_517, ii_518, \
                         ii_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = -gi_263[k]
                   + f_0 * ii_515[k];

        t_376[k] = -gi_264[k]
                   + f_0 * ii_516[k];

        t_377[k] = -gi_265[k]
                   + f_0 * ii_517[k];

        t_378[k] = -gi_266[k]
                   + f_0 * ii_518[k];

        t_379[k] = -gi_267[k]
                   + f_0 * ii_519[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, gi_268, gi_269, gi_270, gi_271, \
                         gi_272, ii_520, ii_521, ii_522, ii_523, \
                         ii_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = -gi_268[k]
                   + f_0 * ii_520[k];

        t_381[k] = -gi_269[k]
                   + f_0 * ii_521[k];

        t_382[k] = -gi_270[k]
                   + f_0 * ii_522[k];

        t_383[k] = -gi_271[k]
                   + f_0 * ii_523[k];

        t_384[k] = -gi_272[k]
                   + f_0 * ii_524[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, gi_273, gi_274, gi_275, gi_276, \
                         gi_277, ii_525, ii_526, ii_527, ii_528, \
                         ii_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -gi_273[k]
                   + f_0 * ii_525[k];

        t_386[k] = -gi_274[k]
                   + f_0 * ii_526[k];

        t_387[k] = -gi_275[k]
                   + f_0 * ii_527[k];

        t_388[k] = -gi_276[k]
                   + f_0 * ii_528[k];

        t_389[k] = -gi_277[k]
                   + f_0 * ii_529[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, t_395, t_396, gi_278, gi_279, \
                         ii_530, ii_531, ii_532, ii_533, ii_534, ii_535, \
                         ii_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = -gi_278[k]
                   + f_0 * ii_530[k];

        t_391[k] = -gi_279[k]
                   + f_0 * ii_531[k];

        t_392[k] = f_0 * ii_532[k];

        t_393[k] = f_0 * ii_533[k];

        t_394[k] = f_0 * ii_534[k];

        t_395[k] = f_0 * ii_535[k];

        t_396[k] = f_0 * ii_536[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, t_402, t_403, t_404, ii_537, \
                         ii_538, ii_539, ii_540, ii_541, ii_542, ii_543, \
                         ii_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_0 * ii_537[k];

        t_398[k] = f_0 * ii_538[k];

        t_399[k] = f_0 * ii_539[k];

        t_400[k] = f_0 * ii_540[k];

        t_401[k] = f_0 * ii_541[k];

        t_402[k] = f_0 * ii_542[k];

        t_403[k] = f_0 * ii_543[k];

        t_404[k] = f_0 * ii_544[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, t_410, t_411, t_412, ii_545, \
                         ii_546, ii_547, ii_548, ii_549, ii_550, ii_551, \
                         ii_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_0 * ii_545[k];

        t_406[k] = f_0 * ii_546[k];

        t_407[k] = f_0 * ii_547[k];

        t_408[k] = f_0 * ii_548[k];

        t_409[k] = f_0 * ii_549[k];

        t_410[k] = f_0 * ii_550[k];

        t_411[k] = f_0 * ii_551[k];

        t_412[k] = f_0 * ii_552[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, t_419, ii_553, ii_554, \
                         ii_555, ii_556, ii_557, ii_558, ii_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_0 * ii_553[k];

        t_414[k] = f_0 * ii_554[k];

        t_415[k] = f_0 * ii_555[k];

        t_416[k] = f_0 * ii_556[k];

        t_417[k] = f_0 * ii_557[k];

        t_418[k] = f_0 * ii_558[k];

        t_419[k] = f_0 * ii_559[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, gi_280, gi_281, gi_282, gi_283, \
                         gi_284, ii_588, ii_589, ii_590, ii_591, \
                         ii_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = -5.0 * gi_280[k]
                   + f_0 * ii_588[k];

        t_421[k] = -5.0 * gi_281[k]
                   + f_0 * ii_589[k];

        t_422[k] = -5.0 * gi_282[k]
                   + f_0 * ii_590[k];

        t_423[k] = -5.0 * gi_283[k]
                   + f_0 * ii_591[k];

        t_424[k] = -5.0 * gi_284[k]
                   + f_0 * ii_592[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, gi_285, gi_286, gi_287, gi_288, \
                         gi_289, ii_593, ii_594, ii_595, ii_596, \
                         ii_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = -5.0 * gi_285[k]
                   + f_0 * ii_593[k];

        t_426[k] = -5.0 * gi_286[k]
                   + f_0 * ii_594[k];

        t_427[k] = -5.0 * gi_287[k]
                   + f_0 * ii_595[k];

        t_428[k] = -5.0 * gi_288[k]
                   + f_0 * ii_596[k];

        t_429[k] = -5.0 * gi_289[k]
                   + f_0 * ii_597[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, gi_290, gi_291, gi_292, gi_293, \
                         gi_294, ii_598, ii_599, ii_600, ii_601, \
                         ii_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -5.0 * gi_290[k]
                   + f_0 * ii_598[k];

        t_431[k] = -5.0 * gi_291[k]
                   + f_0 * ii_599[k];

        t_432[k] = -5.0 * gi_292[k]
                   + f_0 * ii_600[k];

        t_433[k] = -5.0 * gi_293[k]
                   + f_0 * ii_601[k];

        t_434[k] = -5.0 * gi_294[k]
                   + f_0 * ii_602[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, gi_295, gi_296, gi_297, gi_298, \
                         gi_299, ii_603, ii_604, ii_605, ii_606, \
                         ii_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -5.0 * gi_295[k]
                   + f_0 * ii_603[k];

        t_436[k] = -5.0 * gi_296[k]
                   + f_0 * ii_604[k];

        t_437[k] = -5.0 * gi_297[k]
                   + f_0 * ii_605[k];

        t_438[k] = -5.0 * gi_298[k]
                   + f_0 * ii_606[k];

        t_439[k] = -5.0 * gi_299[k]
                   + f_0 * ii_607[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, gi_300, gi_301, gi_302, gi_303, \
                         gi_304, ii_608, ii_609, ii_610, ii_611, \
                         ii_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = -5.0 * gi_300[k]
                   + f_0 * ii_608[k];

        t_441[k] = -5.0 * gi_301[k]
                   + f_0 * ii_609[k];

        t_442[k] = -5.0 * gi_302[k]
                   + f_0 * ii_610[k];

        t_443[k] = -5.0 * gi_303[k]
                   + f_0 * ii_611[k];

        t_444[k] = -5.0 * gi_304[k]
                   + f_0 * ii_612[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, gi_305, gi_306, gi_307, gi_308, \
                         gi_309, ii_613, ii_614, ii_615, ii_616, \
                         ii_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = -5.0 * gi_305[k]
                   + f_0 * ii_613[k];

        t_446[k] = -5.0 * gi_306[k]
                   + f_0 * ii_614[k];

        t_447[k] = -5.0 * gi_307[k]
                   + f_0 * ii_615[k];

        t_448[k] = -4.0 * gi_308[k]
                   + f_0 * ii_616[k];

        t_449[k] = -4.0 * gi_309[k]
                   + f_0 * ii_617[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, gi_310, gi_311, gi_312, gi_313, \
                         gi_314, ii_618, ii_619, ii_620, ii_621, \
                         ii_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -4.0 * gi_310[k]
                   + f_0 * ii_618[k];

        t_451[k] = -4.0 * gi_311[k]
                   + f_0 * ii_619[k];

        t_452[k] = -4.0 * gi_312[k]
                   + f_0 * ii_620[k];

        t_453[k] = -4.0 * gi_313[k]
                   + f_0 * ii_621[k];

        t_454[k] = -4.0 * gi_314[k]
                   + f_0 * ii_622[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, gi_315, gi_316, gi_317, gi_318, \
                         gi_319, ii_623, ii_624, ii_625, ii_626, \
                         ii_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = -4.0 * gi_315[k]
                   + f_0 * ii_623[k];

        t_456[k] = -4.0 * gi_316[k]
                   + f_0 * ii_624[k];

        t_457[k] = -4.0 * gi_317[k]
                   + f_0 * ii_625[k];

        t_458[k] = -4.0 * gi_318[k]
                   + f_0 * ii_626[k];

        t_459[k] = -4.0 * gi_319[k]
                   + f_0 * ii_627[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, gi_320, gi_321, gi_322, gi_323, \
                         gi_324, ii_628, ii_629, ii_630, ii_631, \
                         ii_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = -4.0 * gi_320[k]
                   + f_0 * ii_628[k];

        t_461[k] = -4.0 * gi_321[k]
                   + f_0 * ii_629[k];

        t_462[k] = -4.0 * gi_322[k]
                   + f_0 * ii_630[k];

        t_463[k] = -4.0 * gi_323[k]
                   + f_0 * ii_631[k];

        t_464[k] = -4.0 * gi_324[k]
                   + f_0 * ii_632[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, gi_325, gi_326, gi_327, gi_328, \
                         gi_329, ii_633, ii_634, ii_635, ii_636, \
                         ii_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = -4.0 * gi_325[k]
                   + f_0 * ii_633[k];

        t_466[k] = -4.0 * gi_326[k]
                   + f_0 * ii_634[k];

        t_467[k] = -4.0 * gi_327[k]
                   + f_0 * ii_635[k];

        t_468[k] = -4.0 * gi_328[k]
                   + f_0 * ii_636[k];

        t_469[k] = -4.0 * gi_329[k]
                   + f_0 * ii_637[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, gi_330, gi_331, gi_332, gi_333, \
                         gi_334, ii_638, ii_639, ii_640, ii_641, \
                         ii_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = -4.0 * gi_330[k]
                   + f_0 * ii_638[k];

        t_471[k] = -4.0 * gi_331[k]
                   + f_0 * ii_639[k];

        t_472[k] = -4.0 * gi_332[k]
                   + f_0 * ii_640[k];

        t_473[k] = -4.0 * gi_333[k]
                   + f_0 * ii_641[k];

        t_474[k] = -4.0 * gi_334[k]
                   + f_0 * ii_642[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, gi_335, gi_336, gi_337, gi_338, \
                         gi_339, ii_643, ii_644, ii_645, ii_646, \
                         ii_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = -4.0 * gi_335[k]
                   + f_0 * ii_643[k];

        t_476[k] = -3.0 * gi_336[k]
                   + f_0 * ii_644[k];

        t_477[k] = -3.0 * gi_337[k]
                   + f_0 * ii_645[k];

        t_478[k] = -3.0 * gi_338[k]
                   + f_0 * ii_646[k];

        t_479[k] = -3.0 * gi_339[k]
                   + f_0 * ii_647[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, gi_340, gi_341, gi_342, gi_343, \
                         gi_344, ii_648, ii_649, ii_650, ii_651, \
                         ii_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = -3.0 * gi_340[k]
                   + f_0 * ii_648[k];

        t_481[k] = -3.0 * gi_341[k]
                   + f_0 * ii_649[k];

        t_482[k] = -3.0 * gi_342[k]
                   + f_0 * ii_650[k];

        t_483[k] = -3.0 * gi_343[k]
                   + f_0 * ii_651[k];

        t_484[k] = -3.0 * gi_344[k]
                   + f_0 * ii_652[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, gi_345, gi_346, gi_347, gi_348, \
                         gi_349, ii_653, ii_654, ii_655, ii_656, \
                         ii_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -3.0 * gi_345[k]
                   + f_0 * ii_653[k];

        t_486[k] = -3.0 * gi_346[k]
                   + f_0 * ii_654[k];

        t_487[k] = -3.0 * gi_347[k]
                   + f_0 * ii_655[k];

        t_488[k] = -3.0 * gi_348[k]
                   + f_0 * ii_656[k];

        t_489[k] = -3.0 * gi_349[k]
                   + f_0 * ii_657[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, gi_350, gi_351, gi_352, gi_353, \
                         gi_354, ii_658, ii_659, ii_660, ii_661, \
                         ii_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = -3.0 * gi_350[k]
                   + f_0 * ii_658[k];

        t_491[k] = -3.0 * gi_351[k]
                   + f_0 * ii_659[k];

        t_492[k] = -3.0 * gi_352[k]
                   + f_0 * ii_660[k];

        t_493[k] = -3.0 * gi_353[k]
                   + f_0 * ii_661[k];

        t_494[k] = -3.0 * gi_354[k]
                   + f_0 * ii_662[k];
    }
}

static auto
compute_prim_geom_10_hi_electron_repulsion_1_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_355 = buffer.data(gi + 355);
    const auto *gi_356 = buffer.data(gi + 356);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_363 = buffer.data(gi + 363);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_365 = buffer.data(gi + 365);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_368 = buffer.data(gi + 368);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_371 = buffer.data(gi + 371);
    const auto *gi_372 = buffer.data(gi + 372);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_375 = buffer.data(gi + 375);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_377 = buffer.data(gi + 377);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_379 = buffer.data(gi + 379);
    const auto *gi_380 = buffer.data(gi + 380);
    const auto *gi_381 = buffer.data(gi + 381);
    const auto *gi_382 = buffer.data(gi + 382);
    const auto *gi_383 = buffer.data(gi + 383);
    const auto *gi_384 = buffer.data(gi + 384);
    const auto *gi_385 = buffer.data(gi + 385);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_387 = buffer.data(gi + 387);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_389 = buffer.data(gi + 389);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_391 = buffer.data(gi + 391);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_393 = buffer.data(gi + 393);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_395 = buffer.data(gi + 395);
    const auto *gi_396 = buffer.data(gi + 396);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_399 = buffer.data(gi + 399);
    const auto *gi_400 = buffer.data(gi + 400);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_403 = buffer.data(gi + 403);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_405 = buffer.data(gi + 405);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_408 = buffer.data(gi + 408);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_411 = buffer.data(gi + 411);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

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

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, gi_355, gi_356, gi_357, gi_358, \
                         gi_359, ii_663, ii_664, ii_665, ii_666, \
                         ii_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = -3.0 * gi_355[k]
                   + f_0 * ii_663[k];

        t_496[k] = -3.0 * gi_356[k]
                   + f_0 * ii_664[k];

        t_497[k] = -3.0 * gi_357[k]
                   + f_0 * ii_665[k];

        t_498[k] = -3.0 * gi_358[k]
                   + f_0 * ii_666[k];

        t_499[k] = -3.0 * gi_359[k]
                   + f_0 * ii_667[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, gi_360, gi_361, gi_362, gi_363, \
                         gi_364, ii_668, ii_669, ii_670, ii_671, \
                         ii_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = -3.0 * gi_360[k]
                   + f_0 * ii_668[k];

        t_501[k] = -3.0 * gi_361[k]
                   + f_0 * ii_669[k];

        t_502[k] = -3.0 * gi_362[k]
                   + f_0 * ii_670[k];

        t_503[k] = -3.0 * gi_363[k]
                   + f_0 * ii_671[k];

        t_504[k] = -2.0 * gi_364[k]
                   + f_0 * ii_672[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, gi_365, gi_366, gi_367, gi_368, \
                         gi_369, ii_673, ii_674, ii_675, ii_676, \
                         ii_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -2.0 * gi_365[k]
                   + f_0 * ii_673[k];

        t_506[k] = -2.0 * gi_366[k]
                   + f_0 * ii_674[k];

        t_507[k] = -2.0 * gi_367[k]
                   + f_0 * ii_675[k];

        t_508[k] = -2.0 * gi_368[k]
                   + f_0 * ii_676[k];

        t_509[k] = -2.0 * gi_369[k]
                   + f_0 * ii_677[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, gi_370, gi_371, gi_372, gi_373, \
                         gi_374, ii_678, ii_679, ii_680, ii_681, \
                         ii_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = -2.0 * gi_370[k]
                   + f_0 * ii_678[k];

        t_511[k] = -2.0 * gi_371[k]
                   + f_0 * ii_679[k];

        t_512[k] = -2.0 * gi_372[k]
                   + f_0 * ii_680[k];

        t_513[k] = -2.0 * gi_373[k]
                   + f_0 * ii_681[k];

        t_514[k] = -2.0 * gi_374[k]
                   + f_0 * ii_682[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, gi_375, gi_376, gi_377, gi_378, \
                         gi_379, ii_683, ii_684, ii_685, ii_686, \
                         ii_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = -2.0 * gi_375[k]
                   + f_0 * ii_683[k];

        t_516[k] = -2.0 * gi_376[k]
                   + f_0 * ii_684[k];

        t_517[k] = -2.0 * gi_377[k]
                   + f_0 * ii_685[k];

        t_518[k] = -2.0 * gi_378[k]
                   + f_0 * ii_686[k];

        t_519[k] = -2.0 * gi_379[k]
                   + f_0 * ii_687[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, gi_380, gi_381, gi_382, gi_383, \
                         gi_384, ii_688, ii_689, ii_690, ii_691, \
                         ii_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = -2.0 * gi_380[k]
                   + f_0 * ii_688[k];

        t_521[k] = -2.0 * gi_381[k]
                   + f_0 * ii_689[k];

        t_522[k] = -2.0 * gi_382[k]
                   + f_0 * ii_690[k];

        t_523[k] = -2.0 * gi_383[k]
                   + f_0 * ii_691[k];

        t_524[k] = -2.0 * gi_384[k]
                   + f_0 * ii_692[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, gi_385, gi_386, gi_387, gi_388, \
                         gi_389, ii_693, ii_694, ii_695, ii_696, \
                         ii_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -2.0 * gi_385[k]
                   + f_0 * ii_693[k];

        t_526[k] = -2.0 * gi_386[k]
                   + f_0 * ii_694[k];

        t_527[k] = -2.0 * gi_387[k]
                   + f_0 * ii_695[k];

        t_528[k] = -2.0 * gi_388[k]
                   + f_0 * ii_696[k];

        t_529[k] = -2.0 * gi_389[k]
                   + f_0 * ii_697[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, gi_390, gi_391, gi_392, gi_393, \
                         gi_394, ii_698, ii_699, ii_700, ii_701, \
                         ii_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = -2.0 * gi_390[k]
                   + f_0 * ii_698[k];

        t_531[k] = -2.0 * gi_391[k]
                   + f_0 * ii_699[k];

        t_532[k] = -gi_392[k]
                   + f_0 * ii_700[k];

        t_533[k] = -gi_393[k]
                   + f_0 * ii_701[k];

        t_534[k] = -gi_394[k]
                   + f_0 * ii_702[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, gi_395, gi_396, gi_397, gi_398, \
                         gi_399, ii_703, ii_704, ii_705, ii_706, \
                         ii_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = -gi_395[k]
                   + f_0 * ii_703[k];

        t_536[k] = -gi_396[k]
                   + f_0 * ii_704[k];

        t_537[k] = -gi_397[k]
                   + f_0 * ii_705[k];

        t_538[k] = -gi_398[k]
                   + f_0 * ii_706[k];

        t_539[k] = -gi_399[k]
                   + f_0 * ii_707[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, gi_400, gi_401, gi_402, gi_403, \
                         gi_404, ii_708, ii_709, ii_710, ii_711, \
                         ii_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = -gi_400[k]
                   + f_0 * ii_708[k];

        t_541[k] = -gi_401[k]
                   + f_0 * ii_709[k];

        t_542[k] = -gi_402[k]
                   + f_0 * ii_710[k];

        t_543[k] = -gi_403[k]
                   + f_0 * ii_711[k];

        t_544[k] = -gi_404[k]
                   + f_0 * ii_712[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, gi_405, gi_406, gi_407, gi_408, \
                         gi_409, ii_713, ii_714, ii_715, ii_716, \
                         ii_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = -gi_405[k]
                   + f_0 * ii_713[k];

        t_546[k] = -gi_406[k]
                   + f_0 * ii_714[k];

        t_547[k] = -gi_407[k]
                   + f_0 * ii_715[k];

        t_548[k] = -gi_408[k]
                   + f_0 * ii_716[k];

        t_549[k] = -gi_409[k]
                   + f_0 * ii_717[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, gi_410, gi_411, gi_412, gi_413, \
                         gi_414, ii_718, ii_719, ii_720, ii_721, \
                         ii_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = -gi_410[k]
                   + f_0 * ii_718[k];

        t_551[k] = -gi_411[k]
                   + f_0 * ii_719[k];

        t_552[k] = -gi_412[k]
                   + f_0 * ii_720[k];

        t_553[k] = -gi_413[k]
                   + f_0 * ii_721[k];

        t_554[k] = -gi_414[k]
                   + f_0 * ii_722[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, gi_415, gi_416, gi_417, gi_418, \
                         gi_419, ii_723, ii_724, ii_725, ii_726, \
                         ii_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = -gi_415[k]
                   + f_0 * ii_723[k];

        t_556[k] = -gi_416[k]
                   + f_0 * ii_724[k];

        t_557[k] = -gi_417[k]
                   + f_0 * ii_725[k];

        t_558[k] = -gi_418[k]
                   + f_0 * ii_726[k];

        t_559[k] = -gi_419[k]
                   + f_0 * ii_727[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, t_565, t_566, t_567, ii_728, \
                         ii_729, ii_730, ii_731, ii_732, ii_733, ii_734, \
                         ii_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_0 * ii_728[k];

        t_561[k] = f_0 * ii_729[k];

        t_562[k] = f_0 * ii_730[k];

        t_563[k] = f_0 * ii_731[k];

        t_564[k] = f_0 * ii_732[k];

        t_565[k] = f_0 * ii_733[k];

        t_566[k] = f_0 * ii_734[k];

        t_567[k] = f_0 * ii_735[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, t_572, t_573, t_574, t_575, ii_736, \
                         ii_737, ii_738, ii_739, ii_740, ii_741, ii_742, \
                         ii_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_0 * ii_736[k];

        t_569[k] = f_0 * ii_737[k];

        t_570[k] = f_0 * ii_738[k];

        t_571[k] = f_0 * ii_739[k];

        t_572[k] = f_0 * ii_740[k];

        t_573[k] = f_0 * ii_741[k];

        t_574[k] = f_0 * ii_742[k];

        t_575[k] = f_0 * ii_743[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, t_581, t_582, t_583, ii_744, \
                         ii_745, ii_746, ii_747, ii_748, ii_749, ii_750, \
                         ii_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_0 * ii_744[k];

        t_577[k] = f_0 * ii_745[k];

        t_578[k] = f_0 * ii_746[k];

        t_579[k] = f_0 * ii_747[k];

        t_580[k] = f_0 * ii_748[k];

        t_581[k] = f_0 * ii_749[k];

        t_582[k] = f_0 * ii_750[k];

        t_583[k] = f_0 * ii_751[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, ii_752, ii_753, ii_754, \
                         ii_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_0 * ii_752[k];

        t_585[k] = f_0 * ii_753[k];

        t_586[k] = f_0 * ii_754[k];

        t_587[k] = f_0 * ii_755[k];
    }
}

auto
compute_prim_geom_10_hi_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t gi, const size_t ii,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hi_electron_repulsion_1_piece0(buffer, target, gi, ii, ncols, alpha);

    compute_prim_geom_10_hi_electron_repulsion_1_piece1(buffer, target, gi, ii, ncols, alpha);

    compute_prim_geom_10_hi_electron_repulsion_1_piece2(buffer, target, gi, ii, ncols, alpha);

    compute_prim_geom_10_hi_electron_repulsion_1_piece3(buffer, target, gi, ii, ncols, alpha);
}

static auto
compute_prim_geom_10_hi_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_1 = buffer.data(gi + 1);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_41 = buffer.data(gi + 41);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_44 = buffer.data(gi + 44);
    const auto *gi_45 = buffer.data(gi + 45);
    const auto *gi_46 = buffer.data(gi + 46);
    const auto *gi_47 = buffer.data(gi + 47);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_50 = buffer.data(gi + 50);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_60 = buffer.data(gi + 60);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_63 = buffer.data(gi + 63);
    const auto *gi_64 = buffer.data(gi + 64);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_67 = buffer.data(gi + 67);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_69 = buffer.data(gi + 69);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_71 = buffer.data(gi + 71);
    const auto *gi_72 = buffer.data(gi + 72);
    const auto *gi_73 = buffer.data(gi + 73);
    const auto *gi_74 = buffer.data(gi + 74);
    const auto *gi_75 = buffer.data(gi + 75);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_77 = buffer.data(gi + 77);
    const auto *gi_78 = buffer.data(gi + 78);
    const auto *gi_79 = buffer.data(gi + 79);
    const auto *gi_80 = buffer.data(gi + 80);
    const auto *gi_81 = buffer.data(gi + 81);
    const auto *gi_82 = buffer.data(gi + 82);
    const auto *gi_83 = buffer.data(gi + 83);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, ii_56, ii_57, ii_58, ii_59, \
                         ii_60, ii_61, ii_62, ii_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ii_56[k];

        t_1[k] = f_0 * ii_57[k];

        t_2[k] = f_0 * ii_58[k];

        t_3[k] = f_0 * ii_59[k];

        t_4[k] = f_0 * ii_60[k];

        t_5[k] = f_0 * ii_61[k];

        t_6[k] = f_0 * ii_62[k];

        t_7[k] = f_0 * ii_63[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, ii_64, ii_65, ii_66, \
                         ii_67, ii_68, ii_69, ii_70, ii_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * ii_64[k];

        t_9[k] = f_0 * ii_65[k];

        t_10[k] = f_0 * ii_66[k];

        t_11[k] = f_0 * ii_67[k];

        t_12[k] = f_0 * ii_68[k];

        t_13[k] = f_0 * ii_69[k];

        t_14[k] = f_0 * ii_70[k];

        t_15[k] = f_0 * ii_71[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, ii_72, ii_73, ii_74, \
                         ii_75, ii_76, ii_77, ii_78, ii_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ii_72[k];

        t_17[k] = f_0 * ii_73[k];

        t_18[k] = f_0 * ii_74[k];

        t_19[k] = f_0 * ii_75[k];

        t_20[k] = f_0 * ii_76[k];

        t_21[k] = f_0 * ii_77[k];

        t_22[k] = f_0 * ii_78[k];

        t_23[k] = f_0 * ii_79[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, ii_80, ii_81, ii_82, \
                         ii_83, ii_112, ii_113, ii_114, ii_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * ii_80[k];

        t_25[k] = f_0 * ii_81[k];

        t_26[k] = f_0 * ii_82[k];

        t_27[k] = f_0 * ii_83[k];

        t_28[k] = f_0 * ii_112[k];

        t_29[k] = f_0 * ii_113[k];

        t_30[k] = f_0 * ii_114[k];

        t_31[k] = f_0 * ii_115[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, ii_116, ii_117, \
                         ii_118, ii_119, ii_120, ii_121, ii_122, \
                         ii_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * ii_116[k];

        t_33[k] = f_0 * ii_117[k];

        t_34[k] = f_0 * ii_118[k];

        t_35[k] = f_0 * ii_119[k];

        t_36[k] = f_0 * ii_120[k];

        t_37[k] = f_0 * ii_121[k];

        t_38[k] = f_0 * ii_122[k];

        t_39[k] = f_0 * ii_123[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, ii_124, ii_125, \
                         ii_126, ii_127, ii_128, ii_129, ii_130, \
                         ii_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * ii_124[k];

        t_41[k] = f_0 * ii_125[k];

        t_42[k] = f_0 * ii_126[k];

        t_43[k] = f_0 * ii_127[k];

        t_44[k] = f_0 * ii_128[k];

        t_45[k] = f_0 * ii_129[k];

        t_46[k] = f_0 * ii_130[k];

        t_47[k] = f_0 * ii_131[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, ii_132, ii_133, \
                         ii_134, ii_135, ii_136, ii_137, ii_138, \
                         ii_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * ii_132[k];

        t_49[k] = f_0 * ii_133[k];

        t_50[k] = f_0 * ii_134[k];

        t_51[k] = f_0 * ii_135[k];

        t_52[k] = f_0 * ii_136[k];

        t_53[k] = f_0 * ii_137[k];

        t_54[k] = f_0 * ii_138[k];

        t_55[k] = f_0 * ii_139[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, gi_0, gi_1, gi_2, gi_3, gi_4, ii_140, \
                         ii_141, ii_142, ii_143, ii_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -gi_0[k]
                  + f_0 * ii_140[k];

        t_57[k] = -gi_1[k]
                  + f_0 * ii_141[k];

        t_58[k] = -gi_2[k]
                  + f_0 * ii_142[k];

        t_59[k] = -gi_3[k]
                  + f_0 * ii_143[k];

        t_60[k] = -gi_4[k]
                  + f_0 * ii_144[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, gi_5, gi_6, gi_7, gi_8, gi_9, ii_145, \
                         ii_146, ii_147, ii_148, ii_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -gi_5[k]
                  + f_0 * ii_145[k];

        t_62[k] = -gi_6[k]
                  + f_0 * ii_146[k];

        t_63[k] = -gi_7[k]
                  + f_0 * ii_147[k];

        t_64[k] = -gi_8[k]
                  + f_0 * ii_148[k];

        t_65[k] = -gi_9[k]
                  + f_0 * ii_149[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, gi_10, gi_11, gi_12, gi_13, gi_14, \
                         ii_150, ii_151, ii_152, ii_153, ii_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -gi_10[k]
                  + f_0 * ii_150[k];

        t_67[k] = -gi_11[k]
                  + f_0 * ii_151[k];

        t_68[k] = -gi_12[k]
                  + f_0 * ii_152[k];

        t_69[k] = -gi_13[k]
                  + f_0 * ii_153[k];

        t_70[k] = -gi_14[k]
                  + f_0 * ii_154[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, gi_15, gi_16, gi_17, gi_18, gi_19, \
                         ii_155, ii_156, ii_157, ii_158, ii_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -gi_15[k]
                  + f_0 * ii_155[k];

        t_72[k] = -gi_16[k]
                  + f_0 * ii_156[k];

        t_73[k] = -gi_17[k]
                  + f_0 * ii_157[k];

        t_74[k] = -gi_18[k]
                  + f_0 * ii_158[k];

        t_75[k] = -gi_19[k]
                  + f_0 * ii_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, gi_20, gi_21, gi_22, gi_23, gi_24, \
                         ii_160, ii_161, ii_162, ii_163, ii_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -gi_20[k]
                  + f_0 * ii_160[k];

        t_77[k] = -gi_21[k]
                  + f_0 * ii_161[k];

        t_78[k] = -gi_22[k]
                  + f_0 * ii_162[k];

        t_79[k] = -gi_23[k]
                  + f_0 * ii_163[k];

        t_80[k] = -gi_24[k]
                  + f_0 * ii_164[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, gi_25, gi_26, gi_27, ii_165, \
                         ii_166, ii_167, ii_196, ii_197, ii_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -gi_25[k]
                  + f_0 * ii_165[k];

        t_82[k] = -gi_26[k]
                  + f_0 * ii_166[k];

        t_83[k] = -gi_27[k]
                  + f_0 * ii_167[k];

        t_84[k] = f_0 * ii_196[k];

        t_85[k] = f_0 * ii_197[k];

        t_86[k] = f_0 * ii_198[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, t_93, t_94, ii_199, ii_200, \
                         ii_201, ii_202, ii_203, ii_204, ii_205, \
                         ii_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_0 * ii_199[k];

        t_88[k] = f_0 * ii_200[k];

        t_89[k] = f_0 * ii_201[k];

        t_90[k] = f_0 * ii_202[k];

        t_91[k] = f_0 * ii_203[k];

        t_92[k] = f_0 * ii_204[k];

        t_93[k] = f_0 * ii_205[k];

        t_94[k] = f_0 * ii_206[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, t_101, t_102, ii_207, ii_208, \
                         ii_209, ii_210, ii_211, ii_212, ii_213, \
                         ii_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * ii_207[k];

        t_96[k] = f_0 * ii_208[k];

        t_97[k] = f_0 * ii_209[k];

        t_98[k] = f_0 * ii_210[k];

        t_99[k] = f_0 * ii_211[k];

        t_100[k] = f_0 * ii_212[k];

        t_101[k] = f_0 * ii_213[k];

        t_102[k] = f_0 * ii_214[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, t_109, t_110, ii_215, \
                         ii_216, ii_217, ii_218, ii_219, ii_220, ii_221, \
                         ii_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * ii_215[k];

        t_104[k] = f_0 * ii_216[k];

        t_105[k] = f_0 * ii_217[k];

        t_106[k] = f_0 * ii_218[k];

        t_107[k] = f_0 * ii_219[k];

        t_108[k] = f_0 * ii_220[k];

        t_109[k] = f_0 * ii_221[k];

        t_110[k] = f_0 * ii_222[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, gi_28, gi_29, gi_30, gi_31, \
                         ii_223, ii_224, ii_225, ii_226, ii_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * ii_223[k];

        t_112[k] = -gi_28[k]
                   + f_0 * ii_224[k];

        t_113[k] = -gi_29[k]
                   + f_0 * ii_225[k];

        t_114[k] = -gi_30[k]
                   + f_0 * ii_226[k];

        t_115[k] = -gi_31[k]
                   + f_0 * ii_227[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, gi_32, gi_33, gi_34, gi_35, gi_36, \
                         ii_228, ii_229, ii_230, ii_231, ii_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -gi_32[k]
                   + f_0 * ii_228[k];

        t_117[k] = -gi_33[k]
                   + f_0 * ii_229[k];

        t_118[k] = -gi_34[k]
                   + f_0 * ii_230[k];

        t_119[k] = -gi_35[k]
                   + f_0 * ii_231[k];

        t_120[k] = -gi_36[k]
                   + f_0 * ii_232[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, gi_37, gi_38, gi_39, gi_40, gi_41, \
                         ii_233, ii_234, ii_235, ii_236, ii_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -gi_37[k]
                   + f_0 * ii_233[k];

        t_122[k] = -gi_38[k]
                   + f_0 * ii_234[k];

        t_123[k] = -gi_39[k]
                   + f_0 * ii_235[k];

        t_124[k] = -gi_40[k]
                   + f_0 * ii_236[k];

        t_125[k] = -gi_41[k]
                   + f_0 * ii_237[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, gi_42, gi_43, gi_44, gi_45, gi_46, \
                         ii_238, ii_239, ii_240, ii_241, ii_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -gi_42[k]
                   + f_0 * ii_238[k];

        t_127[k] = -gi_43[k]
                   + f_0 * ii_239[k];

        t_128[k] = -gi_44[k]
                   + f_0 * ii_240[k];

        t_129[k] = -gi_45[k]
                   + f_0 * ii_241[k];

        t_130[k] = -gi_46[k]
                   + f_0 * ii_242[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, gi_47, gi_48, gi_49, gi_50, gi_51, \
                         ii_243, ii_244, ii_245, ii_246, ii_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -gi_47[k]
                   + f_0 * ii_243[k];

        t_132[k] = -gi_48[k]
                   + f_0 * ii_244[k];

        t_133[k] = -gi_49[k]
                   + f_0 * ii_245[k];

        t_134[k] = -gi_50[k]
                   + f_0 * ii_246[k];

        t_135[k] = -gi_51[k]
                   + f_0 * ii_247[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, gi_52, gi_53, gi_54, gi_55, gi_56, \
                         ii_248, ii_249, ii_250, ii_251, ii_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -gi_52[k]
                   + f_0 * ii_248[k];

        t_137[k] = -gi_53[k]
                   + f_0 * ii_249[k];

        t_138[k] = -gi_54[k]
                   + f_0 * ii_250[k];

        t_139[k] = -gi_55[k]
                   + f_0 * ii_251[k];

        t_140[k] = -2.0 * gi_56[k]
                   + f_0 * ii_252[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, gi_57, gi_58, gi_59, gi_60, gi_61, \
                         ii_253, ii_254, ii_255, ii_256, ii_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = -2.0 * gi_57[k]
                   + f_0 * ii_253[k];

        t_142[k] = -2.0 * gi_58[k]
                   + f_0 * ii_254[k];

        t_143[k] = -2.0 * gi_59[k]
                   + f_0 * ii_255[k];

        t_144[k] = -2.0 * gi_60[k]
                   + f_0 * ii_256[k];

        t_145[k] = -2.0 * gi_61[k]
                   + f_0 * ii_257[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, gi_62, gi_63, gi_64, gi_65, gi_66, \
                         ii_258, ii_259, ii_260, ii_261, ii_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -2.0 * gi_62[k]
                   + f_0 * ii_258[k];

        t_147[k] = -2.0 * gi_63[k]
                   + f_0 * ii_259[k];

        t_148[k] = -2.0 * gi_64[k]
                   + f_0 * ii_260[k];

        t_149[k] = -2.0 * gi_65[k]
                   + f_0 * ii_261[k];

        t_150[k] = -2.0 * gi_66[k]
                   + f_0 * ii_262[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, gi_67, gi_68, gi_69, gi_70, gi_71, \
                         ii_263, ii_264, ii_265, ii_266, ii_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = -2.0 * gi_67[k]
                   + f_0 * ii_263[k];

        t_152[k] = -2.0 * gi_68[k]
                   + f_0 * ii_264[k];

        t_153[k] = -2.0 * gi_69[k]
                   + f_0 * ii_265[k];

        t_154[k] = -2.0 * gi_70[k]
                   + f_0 * ii_266[k];

        t_155[k] = -2.0 * gi_71[k]
                   + f_0 * ii_267[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, gi_72, gi_73, gi_74, gi_75, gi_76, \
                         ii_268, ii_269, ii_270, ii_271, ii_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = -2.0 * gi_72[k]
                   + f_0 * ii_268[k];

        t_157[k] = -2.0 * gi_73[k]
                   + f_0 * ii_269[k];

        t_158[k] = -2.0 * gi_74[k]
                   + f_0 * ii_270[k];

        t_159[k] = -2.0 * gi_75[k]
                   + f_0 * ii_271[k];

        t_160[k] = -2.0 * gi_76[k]
                   + f_0 * ii_272[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, gi_77, gi_78, gi_79, gi_80, gi_81, \
                         ii_273, ii_274, ii_275, ii_276, ii_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = -2.0 * gi_77[k]
                   + f_0 * ii_273[k];

        t_162[k] = -2.0 * gi_78[k]
                   + f_0 * ii_274[k];

        t_163[k] = -2.0 * gi_79[k]
                   + f_0 * ii_275[k];

        t_164[k] = -2.0 * gi_80[k]
                   + f_0 * ii_276[k];

        t_165[k] = -2.0 * gi_81[k]
                   + f_0 * ii_277[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, t_172, gi_82, gi_83, \
                         ii_278, ii_279, ii_308, ii_309, ii_310, ii_311, \
                         ii_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = -2.0 * gi_82[k]
                   + f_0 * ii_278[k];

        t_167[k] = -2.0 * gi_83[k]
                   + f_0 * ii_279[k];

        t_168[k] = f_0 * ii_308[k];

        t_169[k] = f_0 * ii_309[k];

        t_170[k] = f_0 * ii_310[k];

        t_171[k] = f_0 * ii_311[k];

        t_172[k] = f_0 * ii_312[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, t_180, ii_313, \
                         ii_314, ii_315, ii_316, ii_317, ii_318, ii_319, \
                         ii_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * ii_313[k];

        t_174[k] = f_0 * ii_314[k];

        t_175[k] = f_0 * ii_315[k];

        t_176[k] = f_0 * ii_316[k];

        t_177[k] = f_0 * ii_317[k];

        t_178[k] = f_0 * ii_318[k];

        t_179[k] = f_0 * ii_319[k];

        t_180[k] = f_0 * ii_320[k];
    }
}

static auto
compute_prim_geom_10_hi_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_85 = buffer.data(gi + 85);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_88 = buffer.data(gi + 88);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_90 = buffer.data(gi + 90);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_92 = buffer.data(gi + 92);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_95 = buffer.data(gi + 95);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_97 = buffer.data(gi + 97);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_100 = buffer.data(gi + 100);
    const auto *gi_101 = buffer.data(gi + 101);
    const auto *gi_102 = buffer.data(gi + 102);
    const auto *gi_103 = buffer.data(gi + 103);
    const auto *gi_104 = buffer.data(gi + 104);
    const auto *gi_105 = buffer.data(gi + 105);
    const auto *gi_106 = buffer.data(gi + 106);
    const auto *gi_107 = buffer.data(gi + 107);
    const auto *gi_108 = buffer.data(gi + 108);
    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_110 = buffer.data(gi + 110);
    const auto *gi_111 = buffer.data(gi + 111);
    const auto *gi_112 = buffer.data(gi + 112);
    const auto *gi_113 = buffer.data(gi + 113);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_116 = buffer.data(gi + 116);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_119 = buffer.data(gi + 119);
    const auto *gi_120 = buffer.data(gi + 120);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_123 = buffer.data(gi + 123);
    const auto *gi_124 = buffer.data(gi + 124);
    const auto *gi_125 = buffer.data(gi + 125);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_127 = buffer.data(gi + 127);
    const auto *gi_128 = buffer.data(gi + 128);
    const auto *gi_129 = buffer.data(gi + 129);
    const auto *gi_130 = buffer.data(gi + 130);
    const auto *gi_131 = buffer.data(gi + 131);
    const auto *gi_132 = buffer.data(gi + 132);
    const auto *gi_133 = buffer.data(gi + 133);
    const auto *gi_134 = buffer.data(gi + 134);
    const auto *gi_135 = buffer.data(gi + 135);
    const auto *gi_136 = buffer.data(gi + 136);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_138 = buffer.data(gi + 138);
    const auto *gi_139 = buffer.data(gi + 139);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_141 = buffer.data(gi + 141);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_144 = buffer.data(gi + 144);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_147 = buffer.data(gi + 147);
    const auto *gi_148 = buffer.data(gi + 148);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_151 = buffer.data(gi + 151);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);
    const auto *gi_158 = buffer.data(gi + 158);
    const auto *gi_159 = buffer.data(gi + 159);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_161 = buffer.data(gi + 161);
    const auto *gi_162 = buffer.data(gi + 162);
    const auto *gi_163 = buffer.data(gi + 163);
    const auto *gi_164 = buffer.data(gi + 164);
    const auto *gi_165 = buffer.data(gi + 165);
    const auto *gi_166 = buffer.data(gi + 166);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_170 = buffer.data(gi + 170);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_172 = buffer.data(gi + 172);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_175 = buffer.data(gi + 175);
    const auto *gi_176 = buffer.data(gi + 176);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_179 = buffer.data(gi + 179);
    const auto *gi_180 = buffer.data(gi + 180);
    const auto *gi_181 = buffer.data(gi + 181);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_184 = buffer.data(gi + 184);
    const auto *gi_185 = buffer.data(gi + 185);
    const auto *gi_186 = buffer.data(gi + 186);
    const auto *gi_187 = buffer.data(gi + 187);
    const auto *gi_188 = buffer.data(gi + 188);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_190 = buffer.data(gi + 190);
    const auto *gi_191 = buffer.data(gi + 191);
    const auto *gi_192 = buffer.data(gi + 192);
    const auto *gi_193 = buffer.data(gi + 193);
    const auto *gi_194 = buffer.data(gi + 194);
    const auto *gi_195 = buffer.data(gi + 195);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_197 = buffer.data(gi + 197);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_200 = buffer.data(gi + 200);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_203 = buffer.data(gi + 203);
    const auto *gi_204 = buffer.data(gi + 204);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);

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

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, t_186, t_187, t_188, ii_321, \
                         ii_322, ii_323, ii_324, ii_325, ii_326, ii_327, \
                         ii_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_0 * ii_321[k];

        t_182[k] = f_0 * ii_322[k];

        t_183[k] = f_0 * ii_323[k];

        t_184[k] = f_0 * ii_324[k];

        t_185[k] = f_0 * ii_325[k];

        t_186[k] = f_0 * ii_326[k];

        t_187[k] = f_0 * ii_327[k];

        t_188[k] = f_0 * ii_328[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, ii_329, ii_330, \
                         ii_331, ii_332, ii_333, ii_334, ii_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * ii_329[k];

        t_190[k] = f_0 * ii_330[k];

        t_191[k] = f_0 * ii_331[k];

        t_192[k] = f_0 * ii_332[k];

        t_193[k] = f_0 * ii_333[k];

        t_194[k] = f_0 * ii_334[k];

        t_195[k] = f_0 * ii_335[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, gi_84, gi_85, gi_86, gi_87, gi_88, \
                         ii_336, ii_337, ii_338, ii_339, ii_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -gi_84[k]
                   + f_0 * ii_336[k];

        t_197[k] = -gi_85[k]
                   + f_0 * ii_337[k];

        t_198[k] = -gi_86[k]
                   + f_0 * ii_338[k];

        t_199[k] = -gi_87[k]
                   + f_0 * ii_339[k];

        t_200[k] = -gi_88[k]
                   + f_0 * ii_340[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, gi_89, gi_90, gi_91, gi_92, gi_93, \
                         ii_341, ii_342, ii_343, ii_344, ii_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -gi_89[k]
                   + f_0 * ii_341[k];

        t_202[k] = -gi_90[k]
                   + f_0 * ii_342[k];

        t_203[k] = -gi_91[k]
                   + f_0 * ii_343[k];

        t_204[k] = -gi_92[k]
                   + f_0 * ii_344[k];

        t_205[k] = -gi_93[k]
                   + f_0 * ii_345[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, gi_94, gi_95, gi_96, gi_97, gi_98, \
                         ii_346, ii_347, ii_348, ii_349, ii_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = -gi_94[k]
                   + f_0 * ii_346[k];

        t_207[k] = -gi_95[k]
                   + f_0 * ii_347[k];

        t_208[k] = -gi_96[k]
                   + f_0 * ii_348[k];

        t_209[k] = -gi_97[k]
                   + f_0 * ii_349[k];

        t_210[k] = -gi_98[k]
                   + f_0 * ii_350[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, gi_99, gi_100, gi_101, gi_102, \
                         gi_103, ii_351, ii_352, ii_353, ii_354, \
                         ii_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -gi_99[k]
                   + f_0 * ii_351[k];

        t_212[k] = -gi_100[k]
                   + f_0 * ii_352[k];

        t_213[k] = -gi_101[k]
                   + f_0 * ii_353[k];

        t_214[k] = -gi_102[k]
                   + f_0 * ii_354[k];

        t_215[k] = -gi_103[k]
                   + f_0 * ii_355[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, gi_104, gi_105, gi_106, gi_107, \
                         gi_108, ii_356, ii_357, ii_358, ii_359, \
                         ii_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = -gi_104[k]
                   + f_0 * ii_356[k];

        t_217[k] = -gi_105[k]
                   + f_0 * ii_357[k];

        t_218[k] = -gi_106[k]
                   + f_0 * ii_358[k];

        t_219[k] = -gi_107[k]
                   + f_0 * ii_359[k];

        t_220[k] = -gi_108[k]
                   + f_0 * ii_360[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, gi_109, gi_110, gi_111, gi_112, \
                         gi_113, ii_361, ii_362, ii_363, ii_364, \
                         ii_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = -gi_109[k]
                   + f_0 * ii_361[k];

        t_222[k] = -gi_110[k]
                   + f_0 * ii_362[k];

        t_223[k] = -gi_111[k]
                   + f_0 * ii_363[k];

        t_224[k] = -2.0 * gi_112[k]
                   + f_0 * ii_364[k];

        t_225[k] = -2.0 * gi_113[k]
                   + f_0 * ii_365[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, t_230, gi_114, gi_115, gi_116, gi_117, \
                         gi_118, ii_366, ii_367, ii_368, ii_369, \
                         ii_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = -2.0 * gi_114[k]
                   + f_0 * ii_366[k];

        t_227[k] = -2.0 * gi_115[k]
                   + f_0 * ii_367[k];

        t_228[k] = -2.0 * gi_116[k]
                   + f_0 * ii_368[k];

        t_229[k] = -2.0 * gi_117[k]
                   + f_0 * ii_369[k];

        t_230[k] = -2.0 * gi_118[k]
                   + f_0 * ii_370[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, gi_119, gi_120, gi_121, gi_122, \
                         gi_123, ii_371, ii_372, ii_373, ii_374, \
                         ii_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = -2.0 * gi_119[k]
                   + f_0 * ii_371[k];

        t_232[k] = -2.0 * gi_120[k]
                   + f_0 * ii_372[k];

        t_233[k] = -2.0 * gi_121[k]
                   + f_0 * ii_373[k];

        t_234[k] = -2.0 * gi_122[k]
                   + f_0 * ii_374[k];

        t_235[k] = -2.0 * gi_123[k]
                   + f_0 * ii_375[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, gi_124, gi_125, gi_126, gi_127, \
                         gi_128, ii_376, ii_377, ii_378, ii_379, \
                         ii_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = -2.0 * gi_124[k]
                   + f_0 * ii_376[k];

        t_237[k] = -2.0 * gi_125[k]
                   + f_0 * ii_377[k];

        t_238[k] = -2.0 * gi_126[k]
                   + f_0 * ii_378[k];

        t_239[k] = -2.0 * gi_127[k]
                   + f_0 * ii_379[k];

        t_240[k] = -2.0 * gi_128[k]
                   + f_0 * ii_380[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, gi_129, gi_130, gi_131, gi_132, \
                         gi_133, ii_381, ii_382, ii_383, ii_384, \
                         ii_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = -2.0 * gi_129[k]
                   + f_0 * ii_381[k];

        t_242[k] = -2.0 * gi_130[k]
                   + f_0 * ii_382[k];

        t_243[k] = -2.0 * gi_131[k]
                   + f_0 * ii_383[k];

        t_244[k] = -2.0 * gi_132[k]
                   + f_0 * ii_384[k];

        t_245[k] = -2.0 * gi_133[k]
                   + f_0 * ii_385[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, gi_134, gi_135, gi_136, gi_137, \
                         gi_138, ii_386, ii_387, ii_388, ii_389, \
                         ii_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = -2.0 * gi_134[k]
                   + f_0 * ii_386[k];

        t_247[k] = -2.0 * gi_135[k]
                   + f_0 * ii_387[k];

        t_248[k] = -2.0 * gi_136[k]
                   + f_0 * ii_388[k];

        t_249[k] = -2.0 * gi_137[k]
                   + f_0 * ii_389[k];

        t_250[k] = -2.0 * gi_138[k]
                   + f_0 * ii_390[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, gi_139, gi_140, gi_141, gi_142, \
                         gi_143, ii_391, ii_392, ii_393, ii_394, \
                         ii_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = -2.0 * gi_139[k]
                   + f_0 * ii_391[k];

        t_252[k] = -3.0 * gi_140[k]
                   + f_0 * ii_392[k];

        t_253[k] = -3.0 * gi_141[k]
                   + f_0 * ii_393[k];

        t_254[k] = -3.0 * gi_142[k]
                   + f_0 * ii_394[k];

        t_255[k] = -3.0 * gi_143[k]
                   + f_0 * ii_395[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, gi_144, gi_145, gi_146, gi_147, \
                         gi_148, ii_396, ii_397, ii_398, ii_399, \
                         ii_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = -3.0 * gi_144[k]
                   + f_0 * ii_396[k];

        t_257[k] = -3.0 * gi_145[k]
                   + f_0 * ii_397[k];

        t_258[k] = -3.0 * gi_146[k]
                   + f_0 * ii_398[k];

        t_259[k] = -3.0 * gi_147[k]
                   + f_0 * ii_399[k];

        t_260[k] = -3.0 * gi_148[k]
                   + f_0 * ii_400[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, gi_149, gi_150, gi_151, gi_152, \
                         gi_153, ii_401, ii_402, ii_403, ii_404, \
                         ii_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = -3.0 * gi_149[k]
                   + f_0 * ii_401[k];

        t_262[k] = -3.0 * gi_150[k]
                   + f_0 * ii_402[k];

        t_263[k] = -3.0 * gi_151[k]
                   + f_0 * ii_403[k];

        t_264[k] = -3.0 * gi_152[k]
                   + f_0 * ii_404[k];

        t_265[k] = -3.0 * gi_153[k]
                   + f_0 * ii_405[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, gi_154, gi_155, gi_156, gi_157, \
                         gi_158, ii_406, ii_407, ii_408, ii_409, \
                         ii_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = -3.0 * gi_154[k]
                   + f_0 * ii_406[k];

        t_267[k] = -3.0 * gi_155[k]
                   + f_0 * ii_407[k];

        t_268[k] = -3.0 * gi_156[k]
                   + f_0 * ii_408[k];

        t_269[k] = -3.0 * gi_157[k]
                   + f_0 * ii_409[k];

        t_270[k] = -3.0 * gi_158[k]
                   + f_0 * ii_410[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, gi_159, gi_160, gi_161, gi_162, \
                         gi_163, ii_411, ii_412, ii_413, ii_414, \
                         ii_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = -3.0 * gi_159[k]
                   + f_0 * ii_411[k];

        t_272[k] = -3.0 * gi_160[k]
                   + f_0 * ii_412[k];

        t_273[k] = -3.0 * gi_161[k]
                   + f_0 * ii_413[k];

        t_274[k] = -3.0 * gi_162[k]
                   + f_0 * ii_414[k];

        t_275[k] = -3.0 * gi_163[k]
                   + f_0 * ii_415[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, t_281, gi_164, gi_165, gi_166, \
                         gi_167, ii_416, ii_417, ii_418, ii_419, ii_448, \
                         ii_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = -3.0 * gi_164[k]
                   + f_0 * ii_416[k];

        t_277[k] = -3.0 * gi_165[k]
                   + f_0 * ii_417[k];

        t_278[k] = -3.0 * gi_166[k]
                   + f_0 * ii_418[k];

        t_279[k] = -3.0 * gi_167[k]
                   + f_0 * ii_419[k];

        t_280[k] = f_0 * ii_448[k];

        t_281[k] = f_0 * ii_449[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, t_286, t_287, t_288, t_289, ii_450, \
                         ii_451, ii_452, ii_453, ii_454, ii_455, ii_456, \
                         ii_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_0 * ii_450[k];

        t_283[k] = f_0 * ii_451[k];

        t_284[k] = f_0 * ii_452[k];

        t_285[k] = f_0 * ii_453[k];

        t_286[k] = f_0 * ii_454[k];

        t_287[k] = f_0 * ii_455[k];

        t_288[k] = f_0 * ii_456[k];

        t_289[k] = f_0 * ii_457[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, t_295, t_296, t_297, ii_458, \
                         ii_459, ii_460, ii_461, ii_462, ii_463, ii_464, \
                         ii_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_0 * ii_458[k];

        t_291[k] = f_0 * ii_459[k];

        t_292[k] = f_0 * ii_460[k];

        t_293[k] = f_0 * ii_461[k];

        t_294[k] = f_0 * ii_462[k];

        t_295[k] = f_0 * ii_463[k];

        t_296[k] = f_0 * ii_464[k];

        t_297[k] = f_0 * ii_465[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, t_303, t_304, t_305, ii_466, \
                         ii_467, ii_468, ii_469, ii_470, ii_471, ii_472, \
                         ii_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_0 * ii_466[k];

        t_299[k] = f_0 * ii_467[k];

        t_300[k] = f_0 * ii_468[k];

        t_301[k] = f_0 * ii_469[k];

        t_302[k] = f_0 * ii_470[k];

        t_303[k] = f_0 * ii_471[k];

        t_304[k] = f_0 * ii_472[k];

        t_305[k] = f_0 * ii_473[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, gi_168, gi_169, gi_170, \
                         gi_171, ii_474, ii_475, ii_476, ii_477, ii_478, \
                         ii_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_0 * ii_474[k];

        t_307[k] = f_0 * ii_475[k];

        t_308[k] = -gi_168[k]
                   + f_0 * ii_476[k];

        t_309[k] = -gi_169[k]
                   + f_0 * ii_477[k];

        t_310[k] = -gi_170[k]
                   + f_0 * ii_478[k];

        t_311[k] = -gi_171[k]
                   + f_0 * ii_479[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, gi_172, gi_173, gi_174, gi_175, \
                         gi_176, ii_480, ii_481, ii_482, ii_483, \
                         ii_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -gi_172[k]
                   + f_0 * ii_480[k];

        t_313[k] = -gi_173[k]
                   + f_0 * ii_481[k];

        t_314[k] = -gi_174[k]
                   + f_0 * ii_482[k];

        t_315[k] = -gi_175[k]
                   + f_0 * ii_483[k];

        t_316[k] = -gi_176[k]
                   + f_0 * ii_484[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, gi_177, gi_178, gi_179, gi_180, \
                         gi_181, ii_485, ii_486, ii_487, ii_488, \
                         ii_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = -gi_177[k]
                   + f_0 * ii_485[k];

        t_318[k] = -gi_178[k]
                   + f_0 * ii_486[k];

        t_319[k] = -gi_179[k]
                   + f_0 * ii_487[k];

        t_320[k] = -gi_180[k]
                   + f_0 * ii_488[k];

        t_321[k] = -gi_181[k]
                   + f_0 * ii_489[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, gi_182, gi_183, gi_184, gi_185, \
                         gi_186, ii_490, ii_491, ii_492, ii_493, \
                         ii_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = -gi_182[k]
                   + f_0 * ii_490[k];

        t_323[k] = -gi_183[k]
                   + f_0 * ii_491[k];

        t_324[k] = -gi_184[k]
                   + f_0 * ii_492[k];

        t_325[k] = -gi_185[k]
                   + f_0 * ii_493[k];

        t_326[k] = -gi_186[k]
                   + f_0 * ii_494[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, gi_187, gi_188, gi_189, gi_190, \
                         gi_191, ii_495, ii_496, ii_497, ii_498, \
                         ii_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = -gi_187[k]
                   + f_0 * ii_495[k];

        t_328[k] = -gi_188[k]
                   + f_0 * ii_496[k];

        t_329[k] = -gi_189[k]
                   + f_0 * ii_497[k];

        t_330[k] = -gi_190[k]
                   + f_0 * ii_498[k];

        t_331[k] = -gi_191[k]
                   + f_0 * ii_499[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, gi_192, gi_193, gi_194, gi_195, \
                         gi_196, ii_500, ii_501, ii_502, ii_503, \
                         ii_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = -gi_192[k]
                   + f_0 * ii_500[k];

        t_333[k] = -gi_193[k]
                   + f_0 * ii_501[k];

        t_334[k] = -gi_194[k]
                   + f_0 * ii_502[k];

        t_335[k] = -gi_195[k]
                   + f_0 * ii_503[k];

        t_336[k] = -2.0 * gi_196[k]
                   + f_0 * ii_504[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, gi_197, gi_198, gi_199, gi_200, \
                         gi_201, ii_505, ii_506, ii_507, ii_508, \
                         ii_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = -2.0 * gi_197[k]
                   + f_0 * ii_505[k];

        t_338[k] = -2.0 * gi_198[k]
                   + f_0 * ii_506[k];

        t_339[k] = -2.0 * gi_199[k]
                   + f_0 * ii_507[k];

        t_340[k] = -2.0 * gi_200[k]
                   + f_0 * ii_508[k];

        t_341[k] = -2.0 * gi_201[k]
                   + f_0 * ii_509[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, gi_202, gi_203, gi_204, gi_205, \
                         gi_206, ii_510, ii_511, ii_512, ii_513, \
                         ii_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = -2.0 * gi_202[k]
                   + f_0 * ii_510[k];

        t_343[k] = -2.0 * gi_203[k]
                   + f_0 * ii_511[k];

        t_344[k] = -2.0 * gi_204[k]
                   + f_0 * ii_512[k];

        t_345[k] = -2.0 * gi_205[k]
                   + f_0 * ii_513[k];

        t_346[k] = -2.0 * gi_206[k]
                   + f_0 * ii_514[k];
    }
}

static auto
compute_prim_geom_10_hi_electron_repulsion_2_piece2(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_207 = buffer.data(gi + 207);
    const auto *gi_208 = buffer.data(gi + 208);
    const auto *gi_209 = buffer.data(gi + 209);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_211 = buffer.data(gi + 211);
    const auto *gi_212 = buffer.data(gi + 212);
    const auto *gi_213 = buffer.data(gi + 213);
    const auto *gi_214 = buffer.data(gi + 214);
    const auto *gi_215 = buffer.data(gi + 215);
    const auto *gi_216 = buffer.data(gi + 216);
    const auto *gi_217 = buffer.data(gi + 217);
    const auto *gi_218 = buffer.data(gi + 218);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_222 = buffer.data(gi + 222);
    const auto *gi_223 = buffer.data(gi + 223);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_225 = buffer.data(gi + 225);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_228 = buffer.data(gi + 228);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_231 = buffer.data(gi + 231);
    const auto *gi_232 = buffer.data(gi + 232);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_235 = buffer.data(gi + 235);
    const auto *gi_236 = buffer.data(gi + 236);
    const auto *gi_237 = buffer.data(gi + 237);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_239 = buffer.data(gi + 239);
    const auto *gi_240 = buffer.data(gi + 240);
    const auto *gi_241 = buffer.data(gi + 241);
    const auto *gi_242 = buffer.data(gi + 242);
    const auto *gi_243 = buffer.data(gi + 243);
    const auto *gi_244 = buffer.data(gi + 244);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_246 = buffer.data(gi + 246);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_250 = buffer.data(gi + 250);
    const auto *gi_251 = buffer.data(gi + 251);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_253 = buffer.data(gi + 253);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_256 = buffer.data(gi + 256);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_258 = buffer.data(gi + 258);
    const auto *gi_259 = buffer.data(gi + 259);
    const auto *gi_260 = buffer.data(gi + 260);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_263 = buffer.data(gi + 263);
    const auto *gi_264 = buffer.data(gi + 264);
    const auto *gi_265 = buffer.data(gi + 265);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_267 = buffer.data(gi + 267);
    const auto *gi_268 = buffer.data(gi + 268);
    const auto *gi_269 = buffer.data(gi + 269);
    const auto *gi_270 = buffer.data(gi + 270);
    const auto *gi_271 = buffer.data(gi + 271);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_273 = buffer.data(gi + 273);
    const auto *gi_274 = buffer.data(gi + 274);
    const auto *gi_275 = buffer.data(gi + 275);
    const auto *gi_276 = buffer.data(gi + 276);
    const auto *gi_277 = buffer.data(gi + 277);
    const auto *gi_278 = buffer.data(gi + 278);
    const auto *gi_279 = buffer.data(gi + 279);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_284 = buffer.data(gi + 284);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_287 = buffer.data(gi + 287);
    const auto *gi_288 = buffer.data(gi + 288);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_291 = buffer.data(gi + 291);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_293 = buffer.data(gi + 293);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_296 = buffer.data(gi + 296);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_299 = buffer.data(gi + 299);
    const auto *gi_300 = buffer.data(gi + 300);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_307 = buffer.data(gi + 307);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_309 = buffer.data(gi + 309);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_312 = buffer.data(gi + 312);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_315 = buffer.data(gi + 315);
    const auto *gi_316 = buffer.data(gi + 316);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_319 = buffer.data(gi + 319);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_321 = buffer.data(gi + 321);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_323 = buffer.data(gi + 323);
    const auto *gi_324 = buffer.data(gi + 324);
    const auto *gi_325 = buffer.data(gi + 325);
    const auto *gi_326 = buffer.data(gi + 326);
    const auto *gi_327 = buffer.data(gi + 327);
    const auto *gi_328 = buffer.data(gi + 328);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_335 = buffer.data(gi + 335);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_337 = buffer.data(gi + 337);
    const auto *gi_338 = buffer.data(gi + 338);

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

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, gi_207, gi_208, gi_209, gi_210, \
                         gi_211, ii_515, ii_516, ii_517, ii_518, \
                         ii_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = -2.0 * gi_207[k]
                   + f_0 * ii_515[k];

        t_348[k] = -2.0 * gi_208[k]
                   + f_0 * ii_516[k];

        t_349[k] = -2.0 * gi_209[k]
                   + f_0 * ii_517[k];

        t_350[k] = -2.0 * gi_210[k]
                   + f_0 * ii_518[k];

        t_351[k] = -2.0 * gi_211[k]
                   + f_0 * ii_519[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, gi_212, gi_213, gi_214, gi_215, \
                         gi_216, ii_520, ii_521, ii_522, ii_523, \
                         ii_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = -2.0 * gi_212[k]
                   + f_0 * ii_520[k];

        t_353[k] = -2.0 * gi_213[k]
                   + f_0 * ii_521[k];

        t_354[k] = -2.0 * gi_214[k]
                   + f_0 * ii_522[k];

        t_355[k] = -2.0 * gi_215[k]
                   + f_0 * ii_523[k];

        t_356[k] = -2.0 * gi_216[k]
                   + f_0 * ii_524[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, gi_217, gi_218, gi_219, gi_220, \
                         gi_221, ii_525, ii_526, ii_527, ii_528, \
                         ii_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -2.0 * gi_217[k]
                   + f_0 * ii_525[k];

        t_358[k] = -2.0 * gi_218[k]
                   + f_0 * ii_526[k];

        t_359[k] = -2.0 * gi_219[k]
                   + f_0 * ii_527[k];

        t_360[k] = -2.0 * gi_220[k]
                   + f_0 * ii_528[k];

        t_361[k] = -2.0 * gi_221[k]
                   + f_0 * ii_529[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, gi_222, gi_223, gi_224, gi_225, \
                         gi_226, ii_530, ii_531, ii_532, ii_533, \
                         ii_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = -2.0 * gi_222[k]
                   + f_0 * ii_530[k];

        t_363[k] = -2.0 * gi_223[k]
                   + f_0 * ii_531[k];

        t_364[k] = -3.0 * gi_224[k]
                   + f_0 * ii_532[k];

        t_365[k] = -3.0 * gi_225[k]
                   + f_0 * ii_533[k];

        t_366[k] = -3.0 * gi_226[k]
                   + f_0 * ii_534[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, gi_227, gi_228, gi_229, gi_230, \
                         gi_231, ii_535, ii_536, ii_537, ii_538, \
                         ii_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = -3.0 * gi_227[k]
                   + f_0 * ii_535[k];

        t_368[k] = -3.0 * gi_228[k]
                   + f_0 * ii_536[k];

        t_369[k] = -3.0 * gi_229[k]
                   + f_0 * ii_537[k];

        t_370[k] = -3.0 * gi_230[k]
                   + f_0 * ii_538[k];

        t_371[k] = -3.0 * gi_231[k]
                   + f_0 * ii_539[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, gi_232, gi_233, gi_234, gi_235, \
                         gi_236, ii_540, ii_541, ii_542, ii_543, \
                         ii_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = -3.0 * gi_232[k]
                   + f_0 * ii_540[k];

        t_373[k] = -3.0 * gi_233[k]
                   + f_0 * ii_541[k];

        t_374[k] = -3.0 * gi_234[k]
                   + f_0 * ii_542[k];

        t_375[k] = -3.0 * gi_235[k]
                   + f_0 * ii_543[k];

        t_376[k] = -3.0 * gi_236[k]
                   + f_0 * ii_544[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, gi_237, gi_238, gi_239, gi_240, \
                         gi_241, ii_545, ii_546, ii_547, ii_548, \
                         ii_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -3.0 * gi_237[k]
                   + f_0 * ii_545[k];

        t_378[k] = -3.0 * gi_238[k]
                   + f_0 * ii_546[k];

        t_379[k] = -3.0 * gi_239[k]
                   + f_0 * ii_547[k];

        t_380[k] = -3.0 * gi_240[k]
                   + f_0 * ii_548[k];

        t_381[k] = -3.0 * gi_241[k]
                   + f_0 * ii_549[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, gi_242, gi_243, gi_244, gi_245, \
                         gi_246, ii_550, ii_551, ii_552, ii_553, \
                         ii_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -3.0 * gi_242[k]
                   + f_0 * ii_550[k];

        t_383[k] = -3.0 * gi_243[k]
                   + f_0 * ii_551[k];

        t_384[k] = -3.0 * gi_244[k]
                   + f_0 * ii_552[k];

        t_385[k] = -3.0 * gi_245[k]
                   + f_0 * ii_553[k];

        t_386[k] = -3.0 * gi_246[k]
                   + f_0 * ii_554[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, t_391, gi_247, gi_248, gi_249, gi_250, \
                         gi_251, ii_555, ii_556, ii_557, ii_558, \
                         ii_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = -3.0 * gi_247[k]
                   + f_0 * ii_555[k];

        t_388[k] = -3.0 * gi_248[k]
                   + f_0 * ii_556[k];

        t_389[k] = -3.0 * gi_249[k]
                   + f_0 * ii_557[k];

        t_390[k] = -3.0 * gi_250[k]
                   + f_0 * ii_558[k];

        t_391[k] = -3.0 * gi_251[k]
                   + f_0 * ii_559[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, t_396, gi_252, gi_253, gi_254, gi_255, \
                         gi_256, ii_560, ii_561, ii_562, ii_563, \
                         ii_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = -4.0 * gi_252[k]
                   + f_0 * ii_560[k];

        t_393[k] = -4.0 * gi_253[k]
                   + f_0 * ii_561[k];

        t_394[k] = -4.0 * gi_254[k]
                   + f_0 * ii_562[k];

        t_395[k] = -4.0 * gi_255[k]
                   + f_0 * ii_563[k];

        t_396[k] = -4.0 * gi_256[k]
                   + f_0 * ii_564[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, gi_257, gi_258, gi_259, gi_260, \
                         gi_261, ii_565, ii_566, ii_567, ii_568, \
                         ii_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = -4.0 * gi_257[k]
                   + f_0 * ii_565[k];

        t_398[k] = -4.0 * gi_258[k]
                   + f_0 * ii_566[k];

        t_399[k] = -4.0 * gi_259[k]
                   + f_0 * ii_567[k];

        t_400[k] = -4.0 * gi_260[k]
                   + f_0 * ii_568[k];

        t_401[k] = -4.0 * gi_261[k]
                   + f_0 * ii_569[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, gi_262, gi_263, gi_264, gi_265, \
                         gi_266, ii_570, ii_571, ii_572, ii_573, \
                         ii_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = -4.0 * gi_262[k]
                   + f_0 * ii_570[k];

        t_403[k] = -4.0 * gi_263[k]
                   + f_0 * ii_571[k];

        t_404[k] = -4.0 * gi_264[k]
                   + f_0 * ii_572[k];

        t_405[k] = -4.0 * gi_265[k]
                   + f_0 * ii_573[k];

        t_406[k] = -4.0 * gi_266[k]
                   + f_0 * ii_574[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, gi_267, gi_268, gi_269, gi_270, \
                         gi_271, ii_575, ii_576, ii_577, ii_578, \
                         ii_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -4.0 * gi_267[k]
                   + f_0 * ii_575[k];

        t_408[k] = -4.0 * gi_268[k]
                   + f_0 * ii_576[k];

        t_409[k] = -4.0 * gi_269[k]
                   + f_0 * ii_577[k];

        t_410[k] = -4.0 * gi_270[k]
                   + f_0 * ii_578[k];

        t_411[k] = -4.0 * gi_271[k]
                   + f_0 * ii_579[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, gi_272, gi_273, gi_274, gi_275, \
                         gi_276, ii_580, ii_581, ii_582, ii_583, \
                         ii_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = -4.0 * gi_272[k]
                   + f_0 * ii_580[k];

        t_413[k] = -4.0 * gi_273[k]
                   + f_0 * ii_581[k];

        t_414[k] = -4.0 * gi_274[k]
                   + f_0 * ii_582[k];

        t_415[k] = -4.0 * gi_275[k]
                   + f_0 * ii_583[k];

        t_416[k] = -4.0 * gi_276[k]
                   + f_0 * ii_584[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, t_421, t_422, gi_277, gi_278, gi_279, \
                         ii_585, ii_586, ii_587, ii_616, ii_617, \
                         ii_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -4.0 * gi_277[k]
                   + f_0 * ii_585[k];

        t_418[k] = -4.0 * gi_278[k]
                   + f_0 * ii_586[k];

        t_419[k] = -4.0 * gi_279[k]
                   + f_0 * ii_587[k];

        t_420[k] = f_0 * ii_616[k];

        t_421[k] = f_0 * ii_617[k];

        t_422[k] = f_0 * ii_618[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, t_428, t_429, t_430, ii_619, \
                         ii_620, ii_621, ii_622, ii_623, ii_624, ii_625, \
                         ii_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_0 * ii_619[k];

        t_424[k] = f_0 * ii_620[k];

        t_425[k] = f_0 * ii_621[k];

        t_426[k] = f_0 * ii_622[k];

        t_427[k] = f_0 * ii_623[k];

        t_428[k] = f_0 * ii_624[k];

        t_429[k] = f_0 * ii_625[k];

        t_430[k] = f_0 * ii_626[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, t_436, t_437, t_438, ii_627, \
                         ii_628, ii_629, ii_630, ii_631, ii_632, ii_633, \
                         ii_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_0 * ii_627[k];

        t_432[k] = f_0 * ii_628[k];

        t_433[k] = f_0 * ii_629[k];

        t_434[k] = f_0 * ii_630[k];

        t_435[k] = f_0 * ii_631[k];

        t_436[k] = f_0 * ii_632[k];

        t_437[k] = f_0 * ii_633[k];

        t_438[k] = f_0 * ii_634[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, t_444, t_445, t_446, ii_635, \
                         ii_636, ii_637, ii_638, ii_639, ii_640, ii_641, \
                         ii_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_0 * ii_635[k];

        t_440[k] = f_0 * ii_636[k];

        t_441[k] = f_0 * ii_637[k];

        t_442[k] = f_0 * ii_638[k];

        t_443[k] = f_0 * ii_639[k];

        t_444[k] = f_0 * ii_640[k];

        t_445[k] = f_0 * ii_641[k];

        t_446[k] = f_0 * ii_642[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, t_451, gi_280, gi_281, gi_282, gi_283, \
                         ii_643, ii_644, ii_645, ii_646, ii_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_0 * ii_643[k];

        t_448[k] = -gi_280[k]
                   + f_0 * ii_644[k];

        t_449[k] = -gi_281[k]
                   + f_0 * ii_645[k];

        t_450[k] = -gi_282[k]
                   + f_0 * ii_646[k];

        t_451[k] = -gi_283[k]
                   + f_0 * ii_647[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, gi_284, gi_285, gi_286, gi_287, \
                         gi_288, ii_648, ii_649, ii_650, ii_651, \
                         ii_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = -gi_284[k]
                   + f_0 * ii_648[k];

        t_453[k] = -gi_285[k]
                   + f_0 * ii_649[k];

        t_454[k] = -gi_286[k]
                   + f_0 * ii_650[k];

        t_455[k] = -gi_287[k]
                   + f_0 * ii_651[k];

        t_456[k] = -gi_288[k]
                   + f_0 * ii_652[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, gi_289, gi_290, gi_291, gi_292, \
                         gi_293, ii_653, ii_654, ii_655, ii_656, \
                         ii_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = -gi_289[k]
                   + f_0 * ii_653[k];

        t_458[k] = -gi_290[k]
                   + f_0 * ii_654[k];

        t_459[k] = -gi_291[k]
                   + f_0 * ii_655[k];

        t_460[k] = -gi_292[k]
                   + f_0 * ii_656[k];

        t_461[k] = -gi_293[k]
                   + f_0 * ii_657[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, gi_294, gi_295, gi_296, gi_297, \
                         gi_298, ii_658, ii_659, ii_660, ii_661, \
                         ii_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -gi_294[k]
                   + f_0 * ii_658[k];

        t_463[k] = -gi_295[k]
                   + f_0 * ii_659[k];

        t_464[k] = -gi_296[k]
                   + f_0 * ii_660[k];

        t_465[k] = -gi_297[k]
                   + f_0 * ii_661[k];

        t_466[k] = -gi_298[k]
                   + f_0 * ii_662[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, t_471, gi_299, gi_300, gi_301, gi_302, \
                         gi_303, ii_663, ii_664, ii_665, ii_666, \
                         ii_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = -gi_299[k]
                   + f_0 * ii_663[k];

        t_468[k] = -gi_300[k]
                   + f_0 * ii_664[k];

        t_469[k] = -gi_301[k]
                   + f_0 * ii_665[k];

        t_470[k] = -gi_302[k]
                   + f_0 * ii_666[k];

        t_471[k] = -gi_303[k]
                   + f_0 * ii_667[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, gi_304, gi_305, gi_306, gi_307, \
                         gi_308, ii_668, ii_669, ii_670, ii_671, \
                         ii_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = -gi_304[k]
                   + f_0 * ii_668[k];

        t_473[k] = -gi_305[k]
                   + f_0 * ii_669[k];

        t_474[k] = -gi_306[k]
                   + f_0 * ii_670[k];

        t_475[k] = -gi_307[k]
                   + f_0 * ii_671[k];

        t_476[k] = -2.0 * gi_308[k]
                   + f_0 * ii_672[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, t_481, gi_309, gi_310, gi_311, gi_312, \
                         gi_313, ii_673, ii_674, ii_675, ii_676, \
                         ii_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = -2.0 * gi_309[k]
                   + f_0 * ii_673[k];

        t_478[k] = -2.0 * gi_310[k]
                   + f_0 * ii_674[k];

        t_479[k] = -2.0 * gi_311[k]
                   + f_0 * ii_675[k];

        t_480[k] = -2.0 * gi_312[k]
                   + f_0 * ii_676[k];

        t_481[k] = -2.0 * gi_313[k]
                   + f_0 * ii_677[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, t_486, gi_314, gi_315, gi_316, gi_317, \
                         gi_318, ii_678, ii_679, ii_680, ii_681, \
                         ii_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -2.0 * gi_314[k]
                   + f_0 * ii_678[k];

        t_483[k] = -2.0 * gi_315[k]
                   + f_0 * ii_679[k];

        t_484[k] = -2.0 * gi_316[k]
                   + f_0 * ii_680[k];

        t_485[k] = -2.0 * gi_317[k]
                   + f_0 * ii_681[k];

        t_486[k] = -2.0 * gi_318[k]
                   + f_0 * ii_682[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, gi_319, gi_320, gi_321, gi_322, \
                         gi_323, ii_683, ii_684, ii_685, ii_686, \
                         ii_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = -2.0 * gi_319[k]
                   + f_0 * ii_683[k];

        t_488[k] = -2.0 * gi_320[k]
                   + f_0 * ii_684[k];

        t_489[k] = -2.0 * gi_321[k]
                   + f_0 * ii_685[k];

        t_490[k] = -2.0 * gi_322[k]
                   + f_0 * ii_686[k];

        t_491[k] = -2.0 * gi_323[k]
                   + f_0 * ii_687[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, gi_324, gi_325, gi_326, gi_327, \
                         gi_328, ii_688, ii_689, ii_690, ii_691, \
                         ii_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = -2.0 * gi_324[k]
                   + f_0 * ii_688[k];

        t_493[k] = -2.0 * gi_325[k]
                   + f_0 * ii_689[k];

        t_494[k] = -2.0 * gi_326[k]
                   + f_0 * ii_690[k];

        t_495[k] = -2.0 * gi_327[k]
                   + f_0 * ii_691[k];

        t_496[k] = -2.0 * gi_328[k]
                   + f_0 * ii_692[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, t_501, gi_329, gi_330, gi_331, gi_332, \
                         gi_333, ii_693, ii_694, ii_695, ii_696, \
                         ii_697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = -2.0 * gi_329[k]
                   + f_0 * ii_693[k];

        t_498[k] = -2.0 * gi_330[k]
                   + f_0 * ii_694[k];

        t_499[k] = -2.0 * gi_331[k]
                   + f_0 * ii_695[k];

        t_500[k] = -2.0 * gi_332[k]
                   + f_0 * ii_696[k];

        t_501[k] = -2.0 * gi_333[k]
                   + f_0 * ii_697[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, gi_334, gi_335, gi_336, gi_337, \
                         gi_338, ii_698, ii_699, ii_700, ii_701, \
                         ii_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = -2.0 * gi_334[k]
                   + f_0 * ii_698[k];

        t_503[k] = -2.0 * gi_335[k]
                   + f_0 * ii_699[k];

        t_504[k] = -3.0 * gi_336[k]
                   + f_0 * ii_700[k];

        t_505[k] = -3.0 * gi_337[k]
                   + f_0 * ii_701[k];

        t_506[k] = -3.0 * gi_338[k]
                   + f_0 * ii_702[k];
    }
}

static auto
compute_prim_geom_10_hi_electron_repulsion_2_piece3(CSimdMatrix &buffer, const size_t target,
                                                    const size_t gi, const size_t ii,
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

    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_340 = buffer.data(gi + 340);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_343 = buffer.data(gi + 343);
    const auto *gi_344 = buffer.data(gi + 344);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_347 = buffer.data(gi + 347);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_349 = buffer.data(gi + 349);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_352 = buffer.data(gi + 352);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);
    const auto *gi_355 = buffer.data(gi + 355);
    const auto *gi_356 = buffer.data(gi + 356);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_363 = buffer.data(gi + 363);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_365 = buffer.data(gi + 365);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_368 = buffer.data(gi + 368);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_371 = buffer.data(gi + 371);
    const auto *gi_372 = buffer.data(gi + 372);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_375 = buffer.data(gi + 375);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_377 = buffer.data(gi + 377);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_379 = buffer.data(gi + 379);
    const auto *gi_380 = buffer.data(gi + 380);
    const auto *gi_381 = buffer.data(gi + 381);
    const auto *gi_382 = buffer.data(gi + 382);
    const auto *gi_383 = buffer.data(gi + 383);
    const auto *gi_384 = buffer.data(gi + 384);
    const auto *gi_385 = buffer.data(gi + 385);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_387 = buffer.data(gi + 387);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_389 = buffer.data(gi + 389);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_391 = buffer.data(gi + 391);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_393 = buffer.data(gi + 393);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_395 = buffer.data(gi + 395);
    const auto *gi_396 = buffer.data(gi + 396);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_399 = buffer.data(gi + 399);
    const auto *gi_400 = buffer.data(gi + 400);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_403 = buffer.data(gi + 403);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_405 = buffer.data(gi + 405);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_408 = buffer.data(gi + 408);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_411 = buffer.data(gi + 411);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

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

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, gi_339, gi_340, gi_341, gi_342, \
                         gi_343, ii_703, ii_704, ii_705, ii_706, \
                         ii_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = -3.0 * gi_339[k]
                   + f_0 * ii_703[k];

        t_508[k] = -3.0 * gi_340[k]
                   + f_0 * ii_704[k];

        t_509[k] = -3.0 * gi_341[k]
                   + f_0 * ii_705[k];

        t_510[k] = -3.0 * gi_342[k]
                   + f_0 * ii_706[k];

        t_511[k] = -3.0 * gi_343[k]
                   + f_0 * ii_707[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, gi_344, gi_345, gi_346, gi_347, \
                         gi_348, ii_708, ii_709, ii_710, ii_711, \
                         ii_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = -3.0 * gi_344[k]
                   + f_0 * ii_708[k];

        t_513[k] = -3.0 * gi_345[k]
                   + f_0 * ii_709[k];

        t_514[k] = -3.0 * gi_346[k]
                   + f_0 * ii_710[k];

        t_515[k] = -3.0 * gi_347[k]
                   + f_0 * ii_711[k];

        t_516[k] = -3.0 * gi_348[k]
                   + f_0 * ii_712[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, gi_349, gi_350, gi_351, gi_352, \
                         gi_353, ii_713, ii_714, ii_715, ii_716, \
                         ii_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -3.0 * gi_349[k]
                   + f_0 * ii_713[k];

        t_518[k] = -3.0 * gi_350[k]
                   + f_0 * ii_714[k];

        t_519[k] = -3.0 * gi_351[k]
                   + f_0 * ii_715[k];

        t_520[k] = -3.0 * gi_352[k]
                   + f_0 * ii_716[k];

        t_521[k] = -3.0 * gi_353[k]
                   + f_0 * ii_717[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, gi_354, gi_355, gi_356, gi_357, \
                         gi_358, ii_718, ii_719, ii_720, ii_721, \
                         ii_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = -3.0 * gi_354[k]
                   + f_0 * ii_718[k];

        t_523[k] = -3.0 * gi_355[k]
                   + f_0 * ii_719[k];

        t_524[k] = -3.0 * gi_356[k]
                   + f_0 * ii_720[k];

        t_525[k] = -3.0 * gi_357[k]
                   + f_0 * ii_721[k];

        t_526[k] = -3.0 * gi_358[k]
                   + f_0 * ii_722[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, gi_359, gi_360, gi_361, gi_362, \
                         gi_363, ii_723, ii_724, ii_725, ii_726, \
                         ii_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = -3.0 * gi_359[k]
                   + f_0 * ii_723[k];

        t_528[k] = -3.0 * gi_360[k]
                   + f_0 * ii_724[k];

        t_529[k] = -3.0 * gi_361[k]
                   + f_0 * ii_725[k];

        t_530[k] = -3.0 * gi_362[k]
                   + f_0 * ii_726[k];

        t_531[k] = -3.0 * gi_363[k]
                   + f_0 * ii_727[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, t_536, gi_364, gi_365, gi_366, gi_367, \
                         gi_368, ii_728, ii_729, ii_730, ii_731, \
                         ii_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = -4.0 * gi_364[k]
                   + f_0 * ii_728[k];

        t_533[k] = -4.0 * gi_365[k]
                   + f_0 * ii_729[k];

        t_534[k] = -4.0 * gi_366[k]
                   + f_0 * ii_730[k];

        t_535[k] = -4.0 * gi_367[k]
                   + f_0 * ii_731[k];

        t_536[k] = -4.0 * gi_368[k]
                   + f_0 * ii_732[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, t_541, gi_369, gi_370, gi_371, gi_372, \
                         gi_373, ii_733, ii_734, ii_735, ii_736, \
                         ii_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = -4.0 * gi_369[k]
                   + f_0 * ii_733[k];

        t_538[k] = -4.0 * gi_370[k]
                   + f_0 * ii_734[k];

        t_539[k] = -4.0 * gi_371[k]
                   + f_0 * ii_735[k];

        t_540[k] = -4.0 * gi_372[k]
                   + f_0 * ii_736[k];

        t_541[k] = -4.0 * gi_373[k]
                   + f_0 * ii_737[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, t_546, gi_374, gi_375, gi_376, gi_377, \
                         gi_378, ii_738, ii_739, ii_740, ii_741, \
                         ii_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = -4.0 * gi_374[k]
                   + f_0 * ii_738[k];

        t_543[k] = -4.0 * gi_375[k]
                   + f_0 * ii_739[k];

        t_544[k] = -4.0 * gi_376[k]
                   + f_0 * ii_740[k];

        t_545[k] = -4.0 * gi_377[k]
                   + f_0 * ii_741[k];

        t_546[k] = -4.0 * gi_378[k]
                   + f_0 * ii_742[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, gi_379, gi_380, gi_381, gi_382, \
                         gi_383, ii_743, ii_744, ii_745, ii_746, \
                         ii_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = -4.0 * gi_379[k]
                   + f_0 * ii_743[k];

        t_548[k] = -4.0 * gi_380[k]
                   + f_0 * ii_744[k];

        t_549[k] = -4.0 * gi_381[k]
                   + f_0 * ii_745[k];

        t_550[k] = -4.0 * gi_382[k]
                   + f_0 * ii_746[k];

        t_551[k] = -4.0 * gi_383[k]
                   + f_0 * ii_747[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, t_556, gi_384, gi_385, gi_386, gi_387, \
                         gi_388, ii_748, ii_749, ii_750, ii_751, \
                         ii_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = -4.0 * gi_384[k]
                   + f_0 * ii_748[k];

        t_553[k] = -4.0 * gi_385[k]
                   + f_0 * ii_749[k];

        t_554[k] = -4.0 * gi_386[k]
                   + f_0 * ii_750[k];

        t_555[k] = -4.0 * gi_387[k]
                   + f_0 * ii_751[k];

        t_556[k] = -4.0 * gi_388[k]
                   + f_0 * ii_752[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, gi_389, gi_390, gi_391, gi_392, \
                         gi_393, ii_753, ii_754, ii_755, ii_756, \
                         ii_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = -4.0 * gi_389[k]
                   + f_0 * ii_753[k];

        t_558[k] = -4.0 * gi_390[k]
                   + f_0 * ii_754[k];

        t_559[k] = -4.0 * gi_391[k]
                   + f_0 * ii_755[k];

        t_560[k] = -5.0 * gi_392[k]
                   + f_0 * ii_756[k];

        t_561[k] = -5.0 * gi_393[k]
                   + f_0 * ii_757[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, gi_394, gi_395, gi_396, gi_397, \
                         gi_398, ii_758, ii_759, ii_760, ii_761, \
                         ii_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = -5.0 * gi_394[k]
                   + f_0 * ii_758[k];

        t_563[k] = -5.0 * gi_395[k]
                   + f_0 * ii_759[k];

        t_564[k] = -5.0 * gi_396[k]
                   + f_0 * ii_760[k];

        t_565[k] = -5.0 * gi_397[k]
                   + f_0 * ii_761[k];

        t_566[k] = -5.0 * gi_398[k]
                   + f_0 * ii_762[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, gi_399, gi_400, gi_401, gi_402, \
                         gi_403, ii_763, ii_764, ii_765, ii_766, \
                         ii_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -5.0 * gi_399[k]
                   + f_0 * ii_763[k];

        t_568[k] = -5.0 * gi_400[k]
                   + f_0 * ii_764[k];

        t_569[k] = -5.0 * gi_401[k]
                   + f_0 * ii_765[k];

        t_570[k] = -5.0 * gi_402[k]
                   + f_0 * ii_766[k];

        t_571[k] = -5.0 * gi_403[k]
                   + f_0 * ii_767[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, gi_404, gi_405, gi_406, gi_407, \
                         gi_408, ii_768, ii_769, ii_770, ii_771, \
                         ii_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = -5.0 * gi_404[k]
                   + f_0 * ii_768[k];

        t_573[k] = -5.0 * gi_405[k]
                   + f_0 * ii_769[k];

        t_574[k] = -5.0 * gi_406[k]
                   + f_0 * ii_770[k];

        t_575[k] = -5.0 * gi_407[k]
                   + f_0 * ii_771[k];

        t_576[k] = -5.0 * gi_408[k]
                   + f_0 * ii_772[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, gi_409, gi_410, gi_411, gi_412, \
                         gi_413, ii_773, ii_774, ii_775, ii_776, \
                         ii_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = -5.0 * gi_409[k]
                   + f_0 * ii_773[k];

        t_578[k] = -5.0 * gi_410[k]
                   + f_0 * ii_774[k];

        t_579[k] = -5.0 * gi_411[k]
                   + f_0 * ii_775[k];

        t_580[k] = -5.0 * gi_412[k]
                   + f_0 * ii_776[k];

        t_581[k] = -5.0 * gi_413[k]
                   + f_0 * ii_777[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, gi_414, gi_415, gi_416, gi_417, \
                         gi_418, ii_778, ii_779, ii_780, ii_781, \
                         ii_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = -5.0 * gi_414[k]
                   + f_0 * ii_778[k];

        t_583[k] = -5.0 * gi_415[k]
                   + f_0 * ii_779[k];

        t_584[k] = -5.0 * gi_416[k]
                   + f_0 * ii_780[k];

        t_585[k] = -5.0 * gi_417[k]
                   + f_0 * ii_781[k];

        t_586[k] = -5.0 * gi_418[k]
                   + f_0 * ii_782[k];
    }

#pragma omp simd aligned(t_587, gi_419, ii_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = -5.0 * gi_419[k]
                   + f_0 * ii_783[k];
    }
}

auto
compute_prim_geom_10_hi_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t gi, const size_t ii,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_hi_electron_repulsion_2_piece0(buffer, target, gi, ii, ncols, alpha);

    compute_prim_geom_10_hi_electron_repulsion_2_piece1(buffer, target, gi, ii, ncols, alpha);

    compute_prim_geom_10_hi_electron_repulsion_2_piece2(buffer, target, gi, ii, ncols, alpha);

    compute_prim_geom_10_hi_electron_repulsion_2_piece3(buffer, target, gi, ii, ncols, alpha);
}

}  // namespace simdt2ceri
