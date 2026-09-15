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


#include "SimdElectronRepulsionGeom10VrrRecFI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_fi_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t di, const size_t gi,
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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_77 = buffer.data(di + 77);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_82 = buffer.data(di + 82);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_86 = buffer.data(di + 86);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_88 = buffer.data(di + 88);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_91 = buffer.data(di + 91);
    const auto *di_92 = buffer.data(di + 92);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_97 = buffer.data(di + 97);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_100 = buffer.data(di + 100);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_103 = buffer.data(di + 103);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_112 = buffer.data(di + 112);
    const auto *di_113 = buffer.data(di + 113);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_116 = buffer.data(di + 116);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_119 = buffer.data(di + 119);
    const auto *di_120 = buffer.data(di + 120);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_123 = buffer.data(di + 123);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_125 = buffer.data(di + 125);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_127 = buffer.data(di + 127);
    const auto *di_128 = buffer.data(di + 128);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_131 = buffer.data(di + 131);
    const auto *di_132 = buffer.data(di + 132);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_141 = buffer.data(di + 141);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_144 = buffer.data(di + 144);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_147 = buffer.data(di + 147);
    const auto *di_148 = buffer.data(di + 148);
    const auto *di_149 = buffer.data(di + 149);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, di_0, di_1, di_2, di_3, di_4, gi_0, gi_1, \
                         gi_2, gi_3, gi_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -3.0 * di_0[k]
                 + f_0 * gi_0[k];

        t_1[k] = -3.0 * di_1[k]
                 + f_0 * gi_1[k];

        t_2[k] = -3.0 * di_2[k]
                 + f_0 * gi_2[k];

        t_3[k] = -3.0 * di_3[k]
                 + f_0 * gi_3[k];

        t_4[k] = -3.0 * di_4[k]
                 + f_0 * gi_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, di_5, di_6, di_7, di_8, di_9, gi_5, gi_6, \
                         gi_7, gi_8, gi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -3.0 * di_5[k]
                 + f_0 * gi_5[k];

        t_6[k] = -3.0 * di_6[k]
                 + f_0 * gi_6[k];

        t_7[k] = -3.0 * di_7[k]
                 + f_0 * gi_7[k];

        t_8[k] = -3.0 * di_8[k]
                 + f_0 * gi_8[k];

        t_9[k] = -3.0 * di_9[k]
                 + f_0 * gi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, di_10, di_11, di_12, di_13, di_14, \
                         gi_10, gi_11, gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -3.0 * di_10[k]
                  + f_0 * gi_10[k];

        t_11[k] = -3.0 * di_11[k]
                  + f_0 * gi_11[k];

        t_12[k] = -3.0 * di_12[k]
                  + f_0 * gi_12[k];

        t_13[k] = -3.0 * di_13[k]
                  + f_0 * gi_13[k];

        t_14[k] = -3.0 * di_14[k]
                  + f_0 * gi_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, di_15, di_16, di_17, di_18, di_19, \
                         gi_15, gi_16, gi_17, gi_18, gi_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -3.0 * di_15[k]
                  + f_0 * gi_15[k];

        t_16[k] = -3.0 * di_16[k]
                  + f_0 * gi_16[k];

        t_17[k] = -3.0 * di_17[k]
                  + f_0 * gi_17[k];

        t_18[k] = -3.0 * di_18[k]
                  + f_0 * gi_18[k];

        t_19[k] = -3.0 * di_19[k]
                  + f_0 * gi_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, di_20, di_21, di_22, di_23, di_24, \
                         gi_20, gi_21, gi_22, gi_23, gi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -3.0 * di_20[k]
                  + f_0 * gi_20[k];

        t_21[k] = -3.0 * di_21[k]
                  + f_0 * gi_21[k];

        t_22[k] = -3.0 * di_22[k]
                  + f_0 * gi_22[k];

        t_23[k] = -3.0 * di_23[k]
                  + f_0 * gi_23[k];

        t_24[k] = -3.0 * di_24[k]
                  + f_0 * gi_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, di_25, di_26, di_27, di_28, di_29, \
                         gi_25, gi_26, gi_27, gi_28, gi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -3.0 * di_25[k]
                  + f_0 * gi_25[k];

        t_26[k] = -3.0 * di_26[k]
                  + f_0 * gi_26[k];

        t_27[k] = -3.0 * di_27[k]
                  + f_0 * gi_27[k];

        t_28[k] = -2.0 * di_28[k]
                  + f_0 * gi_28[k];

        t_29[k] = -2.0 * di_29[k]
                  + f_0 * gi_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, di_30, di_31, di_32, di_33, di_34, \
                         gi_30, gi_31, gi_32, gi_33, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -2.0 * di_30[k]
                  + f_0 * gi_30[k];

        t_31[k] = -2.0 * di_31[k]
                  + f_0 * gi_31[k];

        t_32[k] = -2.0 * di_32[k]
                  + f_0 * gi_32[k];

        t_33[k] = -2.0 * di_33[k]
                  + f_0 * gi_33[k];

        t_34[k] = -2.0 * di_34[k]
                  + f_0 * gi_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, di_35, di_36, di_37, di_38, di_39, \
                         gi_35, gi_36, gi_37, gi_38, gi_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -2.0 * di_35[k]
                  + f_0 * gi_35[k];

        t_36[k] = -2.0 * di_36[k]
                  + f_0 * gi_36[k];

        t_37[k] = -2.0 * di_37[k]
                  + f_0 * gi_37[k];

        t_38[k] = -2.0 * di_38[k]
                  + f_0 * gi_38[k];

        t_39[k] = -2.0 * di_39[k]
                  + f_0 * gi_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, di_40, di_41, di_42, di_43, di_44, \
                         gi_40, gi_41, gi_42, gi_43, gi_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -2.0 * di_40[k]
                  + f_0 * gi_40[k];

        t_41[k] = -2.0 * di_41[k]
                  + f_0 * gi_41[k];

        t_42[k] = -2.0 * di_42[k]
                  + f_0 * gi_42[k];

        t_43[k] = -2.0 * di_43[k]
                  + f_0 * gi_43[k];

        t_44[k] = -2.0 * di_44[k]
                  + f_0 * gi_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, di_45, di_46, di_47, di_48, di_49, \
                         gi_45, gi_46, gi_47, gi_48, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * di_45[k]
                  + f_0 * gi_45[k];

        t_46[k] = -2.0 * di_46[k]
                  + f_0 * gi_46[k];

        t_47[k] = -2.0 * di_47[k]
                  + f_0 * gi_47[k];

        t_48[k] = -2.0 * di_48[k]
                  + f_0 * gi_48[k];

        t_49[k] = -2.0 * di_49[k]
                  + f_0 * gi_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, di_50, di_51, di_52, di_53, di_54, \
                         gi_50, gi_51, gi_52, gi_53, gi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * di_50[k]
                  + f_0 * gi_50[k];

        t_51[k] = -2.0 * di_51[k]
                  + f_0 * gi_51[k];

        t_52[k] = -2.0 * di_52[k]
                  + f_0 * gi_52[k];

        t_53[k] = -2.0 * di_53[k]
                  + f_0 * gi_53[k];

        t_54[k] = -2.0 * di_54[k]
                  + f_0 * gi_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, di_55, di_56, di_57, di_58, di_59, \
                         gi_55, gi_56, gi_57, gi_58, gi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * di_55[k]
                  + f_0 * gi_55[k];

        t_56[k] = -2.0 * di_56[k]
                  + f_0 * gi_56[k];

        t_57[k] = -2.0 * di_57[k]
                  + f_0 * gi_57[k];

        t_58[k] = -2.0 * di_58[k]
                  + f_0 * gi_58[k];

        t_59[k] = -2.0 * di_59[k]
                  + f_0 * gi_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, di_60, di_61, di_62, di_63, di_64, \
                         gi_60, gi_61, gi_62, gi_63, gi_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -2.0 * di_60[k]
                  + f_0 * gi_60[k];

        t_61[k] = -2.0 * di_61[k]
                  + f_0 * gi_61[k];

        t_62[k] = -2.0 * di_62[k]
                  + f_0 * gi_62[k];

        t_63[k] = -2.0 * di_63[k]
                  + f_0 * gi_63[k];

        t_64[k] = -2.0 * di_64[k]
                  + f_0 * gi_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, di_65, di_66, di_67, di_68, di_69, \
                         gi_65, gi_66, gi_67, gi_68, gi_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -2.0 * di_65[k]
                  + f_0 * gi_65[k];

        t_66[k] = -2.0 * di_66[k]
                  + f_0 * gi_66[k];

        t_67[k] = -2.0 * di_67[k]
                  + f_0 * gi_67[k];

        t_68[k] = -2.0 * di_68[k]
                  + f_0 * gi_68[k];

        t_69[k] = -2.0 * di_69[k]
                  + f_0 * gi_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, di_70, di_71, di_72, di_73, di_74, \
                         gi_70, gi_71, gi_72, gi_73, gi_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * di_70[k]
                  + f_0 * gi_70[k];

        t_71[k] = -2.0 * di_71[k]
                  + f_0 * gi_71[k];

        t_72[k] = -2.0 * di_72[k]
                  + f_0 * gi_72[k];

        t_73[k] = -2.0 * di_73[k]
                  + f_0 * gi_73[k];

        t_74[k] = -2.0 * di_74[k]
                  + f_0 * gi_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, di_75, di_76, di_77, di_78, di_79, \
                         gi_75, gi_76, gi_77, gi_78, gi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * di_75[k]
                  + f_0 * gi_75[k];

        t_76[k] = -2.0 * di_76[k]
                  + f_0 * gi_76[k];

        t_77[k] = -2.0 * di_77[k]
                  + f_0 * gi_77[k];

        t_78[k] = -2.0 * di_78[k]
                  + f_0 * gi_78[k];

        t_79[k] = -2.0 * di_79[k]
                  + f_0 * gi_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, di_80, di_81, di_82, di_83, di_84, \
                         gi_80, gi_81, gi_82, gi_83, gi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -2.0 * di_80[k]
                  + f_0 * gi_80[k];

        t_81[k] = -2.0 * di_81[k]
                  + f_0 * gi_81[k];

        t_82[k] = -2.0 * di_82[k]
                  + f_0 * gi_82[k];

        t_83[k] = -2.0 * di_83[k]
                  + f_0 * gi_83[k];

        t_84[k] = -di_84[k]
                  + f_0 * gi_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, di_85, di_86, di_87, di_88, di_89, \
                         gi_85, gi_86, gi_87, gi_88, gi_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -di_85[k]
                  + f_0 * gi_85[k];

        t_86[k] = -di_86[k]
                  + f_0 * gi_86[k];

        t_87[k] = -di_87[k]
                  + f_0 * gi_87[k];

        t_88[k] = -di_88[k]
                  + f_0 * gi_88[k];

        t_89[k] = -di_89[k]
                  + f_0 * gi_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, di_90, di_91, di_92, di_93, di_94, \
                         gi_90, gi_91, gi_92, gi_93, gi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -di_90[k]
                  + f_0 * gi_90[k];

        t_91[k] = -di_91[k]
                  + f_0 * gi_91[k];

        t_92[k] = -di_92[k]
                  + f_0 * gi_92[k];

        t_93[k] = -di_93[k]
                  + f_0 * gi_93[k];

        t_94[k] = -di_94[k]
                  + f_0 * gi_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, di_95, di_96, di_97, di_98, di_99, \
                         gi_95, gi_96, gi_97, gi_98, gi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -di_95[k]
                  + f_0 * gi_95[k];

        t_96[k] = -di_96[k]
                  + f_0 * gi_96[k];

        t_97[k] = -di_97[k]
                  + f_0 * gi_97[k];

        t_98[k] = -di_98[k]
                  + f_0 * gi_98[k];

        t_99[k] = -di_99[k]
                  + f_0 * gi_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, di_100, di_101, di_102, di_103, \
                         di_104, gi_100, gi_101, gi_102, gi_103, \
                         gi_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -di_100[k]
                   + f_0 * gi_100[k];

        t_101[k] = -di_101[k]
                   + f_0 * gi_101[k];

        t_102[k] = -di_102[k]
                   + f_0 * gi_102[k];

        t_103[k] = -di_103[k]
                   + f_0 * gi_103[k];

        t_104[k] = -di_104[k]
                   + f_0 * gi_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, di_105, di_106, di_107, di_108, \
                         di_109, gi_105, gi_106, gi_107, gi_108, \
                         gi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -di_105[k]
                   + f_0 * gi_105[k];

        t_106[k] = -di_106[k]
                   + f_0 * gi_106[k];

        t_107[k] = -di_107[k]
                   + f_0 * gi_107[k];

        t_108[k] = -di_108[k]
                   + f_0 * gi_108[k];

        t_109[k] = -di_109[k]
                   + f_0 * gi_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, di_110, di_111, di_112, di_113, \
                         di_114, gi_110, gi_111, gi_112, gi_113, \
                         gi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -di_110[k]
                   + f_0 * gi_110[k];

        t_111[k] = -di_111[k]
                   + f_0 * gi_111[k];

        t_112[k] = -di_112[k]
                   + f_0 * gi_112[k];

        t_113[k] = -di_113[k]
                   + f_0 * gi_113[k];

        t_114[k] = -di_114[k]
                   + f_0 * gi_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, di_115, di_116, di_117, di_118, \
                         di_119, gi_115, gi_116, gi_117, gi_118, \
                         gi_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -di_115[k]
                   + f_0 * gi_115[k];

        t_116[k] = -di_116[k]
                   + f_0 * gi_116[k];

        t_117[k] = -di_117[k]
                   + f_0 * gi_117[k];

        t_118[k] = -di_118[k]
                   + f_0 * gi_118[k];

        t_119[k] = -di_119[k]
                   + f_0 * gi_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, di_120, di_121, di_122, di_123, \
                         di_124, gi_120, gi_121, gi_122, gi_123, \
                         gi_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -di_120[k]
                   + f_0 * gi_120[k];

        t_121[k] = -di_121[k]
                   + f_0 * gi_121[k];

        t_122[k] = -di_122[k]
                   + f_0 * gi_122[k];

        t_123[k] = -di_123[k]
                   + f_0 * gi_123[k];

        t_124[k] = -di_124[k]
                   + f_0 * gi_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, di_125, di_126, di_127, di_128, \
                         di_129, gi_125, gi_126, gi_127, gi_128, \
                         gi_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -di_125[k]
                   + f_0 * gi_125[k];

        t_126[k] = -di_126[k]
                   + f_0 * gi_126[k];

        t_127[k] = -di_127[k]
                   + f_0 * gi_127[k];

        t_128[k] = -di_128[k]
                   + f_0 * gi_128[k];

        t_129[k] = -di_129[k]
                   + f_0 * gi_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, di_130, di_131, di_132, di_133, \
                         di_134, gi_130, gi_131, gi_132, gi_133, \
                         gi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -di_130[k]
                   + f_0 * gi_130[k];

        t_131[k] = -di_131[k]
                   + f_0 * gi_131[k];

        t_132[k] = -di_132[k]
                   + f_0 * gi_132[k];

        t_133[k] = -di_133[k]
                   + f_0 * gi_133[k];

        t_134[k] = -di_134[k]
                   + f_0 * gi_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, di_135, di_136, di_137, di_138, \
                         di_139, gi_135, gi_136, gi_137, gi_138, \
                         gi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -di_135[k]
                   + f_0 * gi_135[k];

        t_136[k] = -di_136[k]
                   + f_0 * gi_136[k];

        t_137[k] = -di_137[k]
                   + f_0 * gi_137[k];

        t_138[k] = -di_138[k]
                   + f_0 * gi_138[k];

        t_139[k] = -di_139[k]
                   + f_0 * gi_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, di_140, di_141, di_142, di_143, \
                         di_144, gi_140, gi_141, gi_142, gi_143, \
                         gi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -di_140[k]
                   + f_0 * gi_140[k];

        t_141[k] = -di_141[k]
                   + f_0 * gi_141[k];

        t_142[k] = -di_142[k]
                   + f_0 * gi_142[k];

        t_143[k] = -di_143[k]
                   + f_0 * gi_143[k];

        t_144[k] = -di_144[k]
                   + f_0 * gi_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, di_145, di_146, di_147, di_148, \
                         di_149, gi_145, gi_146, gi_147, gi_148, \
                         gi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -di_145[k]
                   + f_0 * gi_145[k];

        t_146[k] = -di_146[k]
                   + f_0 * gi_146[k];

        t_147[k] = -di_147[k]
                   + f_0 * gi_147[k];

        t_148[k] = -di_148[k]
                   + f_0 * gi_148[k];

        t_149[k] = -di_149[k]
                   + f_0 * gi_149[k];
    }
}

static auto
compute_prim_geom_10_fi_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t di, const size_t gi,
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

    const auto *di_150 = buffer.data(di + 150);
    const auto *di_151 = buffer.data(di + 151);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_153 = buffer.data(di + 153);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_156 = buffer.data(di + 156);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_159 = buffer.data(di + 159);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

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

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, di_150, di_151, di_152, di_153, \
                         di_154, gi_150, gi_151, gi_152, gi_153, \
                         gi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -di_150[k]
                   + f_0 * gi_150[k];

        t_151[k] = -di_151[k]
                   + f_0 * gi_151[k];

        t_152[k] = -di_152[k]
                   + f_0 * gi_152[k];

        t_153[k] = -di_153[k]
                   + f_0 * gi_153[k];

        t_154[k] = -di_154[k]
                   + f_0 * gi_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, di_155, di_156, di_157, di_158, \
                         di_159, gi_155, gi_156, gi_157, gi_158, \
                         gi_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -di_155[k]
                   + f_0 * gi_155[k];

        t_156[k] = -di_156[k]
                   + f_0 * gi_156[k];

        t_157[k] = -di_157[k]
                   + f_0 * gi_157[k];

        t_158[k] = -di_158[k]
                   + f_0 * gi_158[k];

        t_159[k] = -di_159[k]
                   + f_0 * gi_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, di_160, di_161, di_162, di_163, \
                         di_164, gi_160, gi_161, gi_162, gi_163, \
                         gi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -di_160[k]
                   + f_0 * gi_160[k];

        t_161[k] = -di_161[k]
                   + f_0 * gi_161[k];

        t_162[k] = -di_162[k]
                   + f_0 * gi_162[k];

        t_163[k] = -di_163[k]
                   + f_0 * gi_163[k];

        t_164[k] = -di_164[k]
                   + f_0 * gi_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, t_170, di_165, di_166, di_167, \
                         gi_165, gi_166, gi_167, gi_168, gi_169, \
                         gi_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -di_165[k]
                   + f_0 * gi_165[k];

        t_166[k] = -di_166[k]
                   + f_0 * gi_166[k];

        t_167[k] = -di_167[k]
                   + f_0 * gi_167[k];

        t_168[k] = f_0 * gi_168[k];

        t_169[k] = f_0 * gi_169[k];

        t_170[k] = f_0 * gi_170[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, t_175, t_176, t_177, t_178, gi_171, \
                         gi_172, gi_173, gi_174, gi_175, gi_176, gi_177, \
                         gi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_0 * gi_171[k];

        t_172[k] = f_0 * gi_172[k];

        t_173[k] = f_0 * gi_173[k];

        t_174[k] = f_0 * gi_174[k];

        t_175[k] = f_0 * gi_175[k];

        t_176[k] = f_0 * gi_176[k];

        t_177[k] = f_0 * gi_177[k];

        t_178[k] = f_0 * gi_178[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, t_184, t_185, t_186, gi_179, \
                         gi_180, gi_181, gi_182, gi_183, gi_184, gi_185, \
                         gi_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_0 * gi_179[k];

        t_180[k] = f_0 * gi_180[k];

        t_181[k] = f_0 * gi_181[k];

        t_182[k] = f_0 * gi_182[k];

        t_183[k] = f_0 * gi_183[k];

        t_184[k] = f_0 * gi_184[k];

        t_185[k] = f_0 * gi_185[k];

        t_186[k] = f_0 * gi_186[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, t_192, t_193, t_194, gi_187, \
                         gi_188, gi_189, gi_190, gi_191, gi_192, gi_193, \
                         gi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_0 * gi_187[k];

        t_188[k] = f_0 * gi_188[k];

        t_189[k] = f_0 * gi_189[k];

        t_190[k] = f_0 * gi_190[k];

        t_191[k] = f_0 * gi_191[k];

        t_192[k] = f_0 * gi_192[k];

        t_193[k] = f_0 * gi_193[k];

        t_194[k] = f_0 * gi_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, t_200, t_201, t_202, gi_195, \
                         gi_196, gi_197, gi_198, gi_199, gi_200, gi_201, \
                         gi_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_0 * gi_195[k];

        t_196[k] = f_0 * gi_196[k];

        t_197[k] = f_0 * gi_197[k];

        t_198[k] = f_0 * gi_198[k];

        t_199[k] = f_0 * gi_199[k];

        t_200[k] = f_0 * gi_200[k];

        t_201[k] = f_0 * gi_201[k];

        t_202[k] = f_0 * gi_202[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, t_209, t_210, gi_203, \
                         gi_204, gi_205, gi_206, gi_207, gi_208, gi_209, \
                         gi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * gi_203[k];

        t_204[k] = f_0 * gi_204[k];

        t_205[k] = f_0 * gi_205[k];

        t_206[k] = f_0 * gi_206[k];

        t_207[k] = f_0 * gi_207[k];

        t_208[k] = f_0 * gi_208[k];

        t_209[k] = f_0 * gi_209[k];

        t_210[k] = f_0 * gi_210[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, t_216, t_217, t_218, gi_211, \
                         gi_212, gi_213, gi_214, gi_215, gi_216, gi_217, \
                         gi_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_0 * gi_211[k];

        t_212[k] = f_0 * gi_212[k];

        t_213[k] = f_0 * gi_213[k];

        t_214[k] = f_0 * gi_214[k];

        t_215[k] = f_0 * gi_215[k];

        t_216[k] = f_0 * gi_216[k];

        t_217[k] = f_0 * gi_217[k];

        t_218[k] = f_0 * gi_218[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, t_224, t_225, t_226, gi_219, \
                         gi_220, gi_221, gi_222, gi_223, gi_224, gi_225, \
                         gi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_0 * gi_219[k];

        t_220[k] = f_0 * gi_220[k];

        t_221[k] = f_0 * gi_221[k];

        t_222[k] = f_0 * gi_222[k];

        t_223[k] = f_0 * gi_223[k];

        t_224[k] = f_0 * gi_224[k];

        t_225[k] = f_0 * gi_225[k];

        t_226[k] = f_0 * gi_226[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, t_232, t_233, t_234, gi_227, \
                         gi_228, gi_229, gi_230, gi_231, gi_232, gi_233, \
                         gi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_0 * gi_227[k];

        t_228[k] = f_0 * gi_228[k];

        t_229[k] = f_0 * gi_229[k];

        t_230[k] = f_0 * gi_230[k];

        t_231[k] = f_0 * gi_231[k];

        t_232[k] = f_0 * gi_232[k];

        t_233[k] = f_0 * gi_233[k];

        t_234[k] = f_0 * gi_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, t_240, t_241, t_242, gi_235, \
                         gi_236, gi_237, gi_238, gi_239, gi_240, gi_241, \
                         gi_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_0 * gi_235[k];

        t_236[k] = f_0 * gi_236[k];

        t_237[k] = f_0 * gi_237[k];

        t_238[k] = f_0 * gi_238[k];

        t_239[k] = f_0 * gi_239[k];

        t_240[k] = f_0 * gi_240[k];

        t_241[k] = f_0 * gi_241[k];

        t_242[k] = f_0 * gi_242[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, t_247, t_248, t_249, t_250, gi_243, \
                         gi_244, gi_245, gi_246, gi_247, gi_248, gi_249, \
                         gi_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_0 * gi_243[k];

        t_244[k] = f_0 * gi_244[k];

        t_245[k] = f_0 * gi_245[k];

        t_246[k] = f_0 * gi_246[k];

        t_247[k] = f_0 * gi_247[k];

        t_248[k] = f_0 * gi_248[k];

        t_249[k] = f_0 * gi_249[k];

        t_250[k] = f_0 * gi_250[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, t_256, t_257, t_258, gi_251, \
                         gi_252, gi_253, gi_254, gi_255, gi_256, gi_257, \
                         gi_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_0 * gi_251[k];

        t_252[k] = f_0 * gi_252[k];

        t_253[k] = f_0 * gi_253[k];

        t_254[k] = f_0 * gi_254[k];

        t_255[k] = f_0 * gi_255[k];

        t_256[k] = f_0 * gi_256[k];

        t_257[k] = f_0 * gi_257[k];

        t_258[k] = f_0 * gi_258[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, t_264, t_265, t_266, gi_259, \
                         gi_260, gi_261, gi_262, gi_263, gi_264, gi_265, \
                         gi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_0 * gi_259[k];

        t_260[k] = f_0 * gi_260[k];

        t_261[k] = f_0 * gi_261[k];

        t_262[k] = f_0 * gi_262[k];

        t_263[k] = f_0 * gi_263[k];

        t_264[k] = f_0 * gi_264[k];

        t_265[k] = f_0 * gi_265[k];

        t_266[k] = f_0 * gi_266[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, t_272, t_273, t_274, gi_267, \
                         gi_268, gi_269, gi_270, gi_271, gi_272, gi_273, \
                         gi_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_0 * gi_267[k];

        t_268[k] = f_0 * gi_268[k];

        t_269[k] = f_0 * gi_269[k];

        t_270[k] = f_0 * gi_270[k];

        t_271[k] = f_0 * gi_271[k];

        t_272[k] = f_0 * gi_272[k];

        t_273[k] = f_0 * gi_273[k];

        t_274[k] = f_0 * gi_274[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, gi_275, gi_276, gi_277, gi_278, \
                         gi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_0 * gi_275[k];

        t_276[k] = f_0 * gi_276[k];

        t_277[k] = f_0 * gi_277[k];

        t_278[k] = f_0 * gi_278[k];

        t_279[k] = f_0 * gi_279[k];
    }
}

auto
compute_prim_geom_10_fi_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t di, const size_t gi,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fi_electron_repulsion_0_piece0(buffer, target, di, gi, ncols, alpha);

    compute_prim_geom_10_fi_electron_repulsion_0_piece1(buffer, target, di, gi, ncols, alpha);
}

static auto
compute_prim_geom_10_fi_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t di, const size_t gi,
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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_77 = buffer.data(di + 77);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_82 = buffer.data(di + 82);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_86 = buffer.data(di + 86);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_88 = buffer.data(di + 88);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);

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
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_284 = buffer.data(gi + 284);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gi_28, gi_29, gi_30, gi_31, \
                         gi_32, gi_33, gi_34, gi_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gi_28[k];

        t_1[k] = f_0 * gi_29[k];

        t_2[k] = f_0 * gi_30[k];

        t_3[k] = f_0 * gi_31[k];

        t_4[k] = f_0 * gi_32[k];

        t_5[k] = f_0 * gi_33[k];

        t_6[k] = f_0 * gi_34[k];

        t_7[k] = f_0 * gi_35[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gi_36, gi_37, gi_38, \
                         gi_39, gi_40, gi_41, gi_42, gi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gi_36[k];

        t_9[k] = f_0 * gi_37[k];

        t_10[k] = f_0 * gi_38[k];

        t_11[k] = f_0 * gi_39[k];

        t_12[k] = f_0 * gi_40[k];

        t_13[k] = f_0 * gi_41[k];

        t_14[k] = f_0 * gi_42[k];

        t_15[k] = f_0 * gi_43[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, gi_44, gi_45, gi_46, \
                         gi_47, gi_48, gi_49, gi_50, gi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gi_44[k];

        t_17[k] = f_0 * gi_45[k];

        t_18[k] = f_0 * gi_46[k];

        t_19[k] = f_0 * gi_47[k];

        t_20[k] = f_0 * gi_48[k];

        t_21[k] = f_0 * gi_49[k];

        t_22[k] = f_0 * gi_50[k];

        t_23[k] = f_0 * gi_51[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, di_0, di_1, gi_52, gi_53, gi_54, \
                         gi_55, gi_84, gi_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * gi_52[k];

        t_25[k] = f_0 * gi_53[k];

        t_26[k] = f_0 * gi_54[k];

        t_27[k] = f_0 * gi_55[k];

        t_28[k] = -di_0[k]
                  + f_0 * gi_84[k];

        t_29[k] = -di_1[k]
                  + f_0 * gi_85[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, di_2, di_3, di_4, di_5, di_6, gi_86, \
                         gi_87, gi_88, gi_89, gi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -di_2[k]
                  + f_0 * gi_86[k];

        t_31[k] = -di_3[k]
                  + f_0 * gi_87[k];

        t_32[k] = -di_4[k]
                  + f_0 * gi_88[k];

        t_33[k] = -di_5[k]
                  + f_0 * gi_89[k];

        t_34[k] = -di_6[k]
                  + f_0 * gi_90[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, di_7, di_8, di_9, di_10, di_11, gi_91, \
                         gi_92, gi_93, gi_94, gi_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -di_7[k]
                  + f_0 * gi_91[k];

        t_36[k] = -di_8[k]
                  + f_0 * gi_92[k];

        t_37[k] = -di_9[k]
                  + f_0 * gi_93[k];

        t_38[k] = -di_10[k]
                  + f_0 * gi_94[k];

        t_39[k] = -di_11[k]
                  + f_0 * gi_95[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, di_12, di_13, di_14, di_15, di_16, \
                         gi_96, gi_97, gi_98, gi_99, gi_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -di_12[k]
                  + f_0 * gi_96[k];

        t_41[k] = -di_13[k]
                  + f_0 * gi_97[k];

        t_42[k] = -di_14[k]
                  + f_0 * gi_98[k];

        t_43[k] = -di_15[k]
                  + f_0 * gi_99[k];

        t_44[k] = -di_16[k]
                  + f_0 * gi_100[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, di_17, di_18, di_19, di_20, di_21, \
                         gi_101, gi_102, gi_103, gi_104, gi_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -di_17[k]
                  + f_0 * gi_101[k];

        t_46[k] = -di_18[k]
                  + f_0 * gi_102[k];

        t_47[k] = -di_19[k]
                  + f_0 * gi_103[k];

        t_48[k] = -di_20[k]
                  + f_0 * gi_104[k];

        t_49[k] = -di_21[k]
                  + f_0 * gi_105[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, di_22, di_23, di_24, di_25, di_26, \
                         gi_106, gi_107, gi_108, gi_109, gi_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -di_22[k]
                  + f_0 * gi_106[k];

        t_51[k] = -di_23[k]
                  + f_0 * gi_107[k];

        t_52[k] = -di_24[k]
                  + f_0 * gi_108[k];

        t_53[k] = -di_25[k]
                  + f_0 * gi_109[k];

        t_54[k] = -di_26[k]
                  + f_0 * gi_110[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, t_61, di_27, gi_111, gi_112, \
                         gi_113, gi_114, gi_115, gi_116, gi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -di_27[k]
                  + f_0 * gi_111[k];

        t_56[k] = f_0 * gi_112[k];

        t_57[k] = f_0 * gi_113[k];

        t_58[k] = f_0 * gi_114[k];

        t_59[k] = f_0 * gi_115[k];

        t_60[k] = f_0 * gi_116[k];

        t_61[k] = f_0 * gi_117[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, t_67, t_68, t_69, gi_118, gi_119, \
                         gi_120, gi_121, gi_122, gi_123, gi_124, \
                         gi_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * gi_118[k];

        t_63[k] = f_0 * gi_119[k];

        t_64[k] = f_0 * gi_120[k];

        t_65[k] = f_0 * gi_121[k];

        t_66[k] = f_0 * gi_122[k];

        t_67[k] = f_0 * gi_123[k];

        t_68[k] = f_0 * gi_124[k];

        t_69[k] = f_0 * gi_125[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, t_76, t_77, gi_126, gi_127, \
                         gi_128, gi_129, gi_130, gi_131, gi_132, \
                         gi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_0 * gi_126[k];

        t_71[k] = f_0 * gi_127[k];

        t_72[k] = f_0 * gi_128[k];

        t_73[k] = f_0 * gi_129[k];

        t_74[k] = f_0 * gi_130[k];

        t_75[k] = f_0 * gi_131[k];

        t_76[k] = f_0 * gi_132[k];

        t_77[k] = f_0 * gi_133[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, t_84, di_28, gi_134, gi_135, \
                         gi_136, gi_137, gi_138, gi_139, gi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * gi_134[k];

        t_79[k] = f_0 * gi_135[k];

        t_80[k] = f_0 * gi_136[k];

        t_81[k] = f_0 * gi_137[k];

        t_82[k] = f_0 * gi_138[k];

        t_83[k] = f_0 * gi_139[k];

        t_84[k] = -2.0 * di_28[k]
                  + f_0 * gi_168[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, di_29, di_30, di_31, di_32, di_33, \
                         gi_169, gi_170, gi_171, gi_172, gi_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -2.0 * di_29[k]
                  + f_0 * gi_169[k];

        t_86[k] = -2.0 * di_30[k]
                  + f_0 * gi_170[k];

        t_87[k] = -2.0 * di_31[k]
                  + f_0 * gi_171[k];

        t_88[k] = -2.0 * di_32[k]
                  + f_0 * gi_172[k];

        t_89[k] = -2.0 * di_33[k]
                  + f_0 * gi_173[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, di_34, di_35, di_36, di_37, di_38, \
                         gi_174, gi_175, gi_176, gi_177, gi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -2.0 * di_34[k]
                  + f_0 * gi_174[k];

        t_91[k] = -2.0 * di_35[k]
                  + f_0 * gi_175[k];

        t_92[k] = -2.0 * di_36[k]
                  + f_0 * gi_176[k];

        t_93[k] = -2.0 * di_37[k]
                  + f_0 * gi_177[k];

        t_94[k] = -2.0 * di_38[k]
                  + f_0 * gi_178[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, di_39, di_40, di_41, di_42, di_43, \
                         gi_179, gi_180, gi_181, gi_182, gi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -2.0 * di_39[k]
                  + f_0 * gi_179[k];

        t_96[k] = -2.0 * di_40[k]
                  + f_0 * gi_180[k];

        t_97[k] = -2.0 * di_41[k]
                  + f_0 * gi_181[k];

        t_98[k] = -2.0 * di_42[k]
                  + f_0 * gi_182[k];

        t_99[k] = -2.0 * di_43[k]
                  + f_0 * gi_183[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, di_44, di_45, di_46, di_47, di_48, \
                         gi_184, gi_185, gi_186, gi_187, gi_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * di_44[k]
                   + f_0 * gi_184[k];

        t_101[k] = -2.0 * di_45[k]
                   + f_0 * gi_185[k];

        t_102[k] = -2.0 * di_46[k]
                   + f_0 * gi_186[k];

        t_103[k] = -2.0 * di_47[k]
                   + f_0 * gi_187[k];

        t_104[k] = -2.0 * di_48[k]
                   + f_0 * gi_188[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, di_49, di_50, di_51, di_52, di_53, \
                         gi_189, gi_190, gi_191, gi_192, gi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * di_49[k]
                   + f_0 * gi_189[k];

        t_106[k] = -2.0 * di_50[k]
                   + f_0 * gi_190[k];

        t_107[k] = -2.0 * di_51[k]
                   + f_0 * gi_191[k];

        t_108[k] = -2.0 * di_52[k]
                   + f_0 * gi_192[k];

        t_109[k] = -2.0 * di_53[k]
                   + f_0 * gi_193[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, di_54, di_55, di_56, di_57, di_58, \
                         gi_194, gi_195, gi_196, gi_197, gi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * di_54[k]
                   + f_0 * gi_194[k];

        t_111[k] = -2.0 * di_55[k]
                   + f_0 * gi_195[k];

        t_112[k] = -di_56[k]
                   + f_0 * gi_196[k];

        t_113[k] = -di_57[k]
                   + f_0 * gi_197[k];

        t_114[k] = -di_58[k]
                   + f_0 * gi_198[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, di_59, di_60, di_61, di_62, di_63, \
                         gi_199, gi_200, gi_201, gi_202, gi_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -di_59[k]
                   + f_0 * gi_199[k];

        t_116[k] = -di_60[k]
                   + f_0 * gi_200[k];

        t_117[k] = -di_61[k]
                   + f_0 * gi_201[k];

        t_118[k] = -di_62[k]
                   + f_0 * gi_202[k];

        t_119[k] = -di_63[k]
                   + f_0 * gi_203[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, di_64, di_65, di_66, di_67, di_68, \
                         gi_204, gi_205, gi_206, gi_207, gi_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -di_64[k]
                   + f_0 * gi_204[k];

        t_121[k] = -di_65[k]
                   + f_0 * gi_205[k];

        t_122[k] = -di_66[k]
                   + f_0 * gi_206[k];

        t_123[k] = -di_67[k]
                   + f_0 * gi_207[k];

        t_124[k] = -di_68[k]
                   + f_0 * gi_208[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, di_69, di_70, di_71, di_72, di_73, \
                         gi_209, gi_210, gi_211, gi_212, gi_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -di_69[k]
                   + f_0 * gi_209[k];

        t_126[k] = -di_70[k]
                   + f_0 * gi_210[k];

        t_127[k] = -di_71[k]
                   + f_0 * gi_211[k];

        t_128[k] = -di_72[k]
                   + f_0 * gi_212[k];

        t_129[k] = -di_73[k]
                   + f_0 * gi_213[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, di_74, di_75, di_76, di_77, di_78, \
                         gi_214, gi_215, gi_216, gi_217, gi_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -di_74[k]
                   + f_0 * gi_214[k];

        t_131[k] = -di_75[k]
                   + f_0 * gi_215[k];

        t_132[k] = -di_76[k]
                   + f_0 * gi_216[k];

        t_133[k] = -di_77[k]
                   + f_0 * gi_217[k];

        t_134[k] = -di_78[k]
                   + f_0 * gi_218[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, di_79, di_80, di_81, di_82, di_83, \
                         gi_219, gi_220, gi_221, gi_222, gi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -di_79[k]
                   + f_0 * gi_219[k];

        t_136[k] = -di_80[k]
                   + f_0 * gi_220[k];

        t_137[k] = -di_81[k]
                   + f_0 * gi_221[k];

        t_138[k] = -di_82[k]
                   + f_0 * gi_222[k];

        t_139[k] = -di_83[k]
                   + f_0 * gi_223[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, gi_224, \
                         gi_225, gi_226, gi_227, gi_228, gi_229, gi_230, \
                         gi_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * gi_224[k];

        t_141[k] = f_0 * gi_225[k];

        t_142[k] = f_0 * gi_226[k];

        t_143[k] = f_0 * gi_227[k];

        t_144[k] = f_0 * gi_228[k];

        t_145[k] = f_0 * gi_229[k];

        t_146[k] = f_0 * gi_230[k];

        t_147[k] = f_0 * gi_231[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, t_154, t_155, gi_232, \
                         gi_233, gi_234, gi_235, gi_236, gi_237, gi_238, \
                         gi_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * gi_232[k];

        t_149[k] = f_0 * gi_233[k];

        t_150[k] = f_0 * gi_234[k];

        t_151[k] = f_0 * gi_235[k];

        t_152[k] = f_0 * gi_236[k];

        t_153[k] = f_0 * gi_237[k];

        t_154[k] = f_0 * gi_238[k];

        t_155[k] = f_0 * gi_239[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, t_162, t_163, gi_240, \
                         gi_241, gi_242, gi_243, gi_244, gi_245, gi_246, \
                         gi_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * gi_240[k];

        t_157[k] = f_0 * gi_241[k];

        t_158[k] = f_0 * gi_242[k];

        t_159[k] = f_0 * gi_243[k];

        t_160[k] = f_0 * gi_244[k];

        t_161[k] = f_0 * gi_245[k];

        t_162[k] = f_0 * gi_246[k];

        t_163[k] = f_0 * gi_247[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, di_84, di_85, gi_248, \
                         gi_249, gi_250, gi_251, gi_280, gi_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_0 * gi_248[k];

        t_165[k] = f_0 * gi_249[k];

        t_166[k] = f_0 * gi_250[k];

        t_167[k] = f_0 * gi_251[k];

        t_168[k] = -3.0 * di_84[k]
                   + f_0 * gi_280[k];

        t_169[k] = -3.0 * di_85[k]
                   + f_0 * gi_281[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, di_86, di_87, di_88, di_89, di_90, \
                         gi_282, gi_283, gi_284, gi_285, gi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -3.0 * di_86[k]
                   + f_0 * gi_282[k];

        t_171[k] = -3.0 * di_87[k]
                   + f_0 * gi_283[k];

        t_172[k] = -3.0 * di_88[k]
                   + f_0 * gi_284[k];

        t_173[k] = -3.0 * di_89[k]
                   + f_0 * gi_285[k];

        t_174[k] = -3.0 * di_90[k]
                   + f_0 * gi_286[k];
    }
}

static auto
compute_prim_geom_10_fi_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t di, const size_t gi,
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

    const auto *di_91 = buffer.data(di + 91);
    const auto *di_92 = buffer.data(di + 92);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_97 = buffer.data(di + 97);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_100 = buffer.data(di + 100);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_103 = buffer.data(di + 103);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_112 = buffer.data(di + 112);
    const auto *di_113 = buffer.data(di + 113);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_116 = buffer.data(di + 116);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_119 = buffer.data(di + 119);
    const auto *di_120 = buffer.data(di + 120);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_123 = buffer.data(di + 123);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_125 = buffer.data(di + 125);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_127 = buffer.data(di + 127);
    const auto *di_128 = buffer.data(di + 128);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_131 = buffer.data(di + 131);
    const auto *di_132 = buffer.data(di + 132);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_141 = buffer.data(di + 141);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_144 = buffer.data(di + 144);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_147 = buffer.data(di + 147);
    const auto *di_148 = buffer.data(di + 148);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_150 = buffer.data(di + 150);
    const auto *di_151 = buffer.data(di + 151);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_153 = buffer.data(di + 153);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_156 = buffer.data(di + 156);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_159 = buffer.data(di + 159);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

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

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, di_91, di_92, di_93, di_94, di_95, \
                         gi_287, gi_288, gi_289, gi_290, gi_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -3.0 * di_91[k]
                   + f_0 * gi_287[k];

        t_176[k] = -3.0 * di_92[k]
                   + f_0 * gi_288[k];

        t_177[k] = -3.0 * di_93[k]
                   + f_0 * gi_289[k];

        t_178[k] = -3.0 * di_94[k]
                   + f_0 * gi_290[k];

        t_179[k] = -3.0 * di_95[k]
                   + f_0 * gi_291[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, di_96, di_97, di_98, di_99, \
                         di_100, gi_292, gi_293, gi_294, gi_295, \
                         gi_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * di_96[k]
                   + f_0 * gi_292[k];

        t_181[k] = -3.0 * di_97[k]
                   + f_0 * gi_293[k];

        t_182[k] = -3.0 * di_98[k]
                   + f_0 * gi_294[k];

        t_183[k] = -3.0 * di_99[k]
                   + f_0 * gi_295[k];

        t_184[k] = -3.0 * di_100[k]
                   + f_0 * gi_296[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, di_101, di_102, di_103, di_104, \
                         di_105, gi_297, gi_298, gi_299, gi_300, \
                         gi_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * di_101[k]
                   + f_0 * gi_297[k];

        t_186[k] = -3.0 * di_102[k]
                   + f_0 * gi_298[k];

        t_187[k] = -3.0 * di_103[k]
                   + f_0 * gi_299[k];

        t_188[k] = -3.0 * di_104[k]
                   + f_0 * gi_300[k];

        t_189[k] = -3.0 * di_105[k]
                   + f_0 * gi_301[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, di_106, di_107, di_108, di_109, \
                         di_110, gi_302, gi_303, gi_304, gi_305, \
                         gi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -3.0 * di_106[k]
                   + f_0 * gi_302[k];

        t_191[k] = -3.0 * di_107[k]
                   + f_0 * gi_303[k];

        t_192[k] = -3.0 * di_108[k]
                   + f_0 * gi_304[k];

        t_193[k] = -3.0 * di_109[k]
                   + f_0 * gi_305[k];

        t_194[k] = -3.0 * di_110[k]
                   + f_0 * gi_306[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, di_111, di_112, di_113, di_114, \
                         di_115, gi_307, gi_308, gi_309, gi_310, \
                         gi_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -3.0 * di_111[k]
                   + f_0 * gi_307[k];

        t_196[k] = -2.0 * di_112[k]
                   + f_0 * gi_308[k];

        t_197[k] = -2.0 * di_113[k]
                   + f_0 * gi_309[k];

        t_198[k] = -2.0 * di_114[k]
                   + f_0 * gi_310[k];

        t_199[k] = -2.0 * di_115[k]
                   + f_0 * gi_311[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, di_116, di_117, di_118, di_119, \
                         di_120, gi_312, gi_313, gi_314, gi_315, \
                         gi_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -2.0 * di_116[k]
                   + f_0 * gi_312[k];

        t_201[k] = -2.0 * di_117[k]
                   + f_0 * gi_313[k];

        t_202[k] = -2.0 * di_118[k]
                   + f_0 * gi_314[k];

        t_203[k] = -2.0 * di_119[k]
                   + f_0 * gi_315[k];

        t_204[k] = -2.0 * di_120[k]
                   + f_0 * gi_316[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, di_121, di_122, di_123, di_124, \
                         di_125, gi_317, gi_318, gi_319, gi_320, \
                         gi_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -2.0 * di_121[k]
                   + f_0 * gi_317[k];

        t_206[k] = -2.0 * di_122[k]
                   + f_0 * gi_318[k];

        t_207[k] = -2.0 * di_123[k]
                   + f_0 * gi_319[k];

        t_208[k] = -2.0 * di_124[k]
                   + f_0 * gi_320[k];

        t_209[k] = -2.0 * di_125[k]
                   + f_0 * gi_321[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, di_126, di_127, di_128, di_129, \
                         di_130, gi_322, gi_323, gi_324, gi_325, \
                         gi_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -2.0 * di_126[k]
                   + f_0 * gi_322[k];

        t_211[k] = -2.0 * di_127[k]
                   + f_0 * gi_323[k];

        t_212[k] = -2.0 * di_128[k]
                   + f_0 * gi_324[k];

        t_213[k] = -2.0 * di_129[k]
                   + f_0 * gi_325[k];

        t_214[k] = -2.0 * di_130[k]
                   + f_0 * gi_326[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, di_131, di_132, di_133, di_134, \
                         di_135, gi_327, gi_328, gi_329, gi_330, \
                         gi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -2.0 * di_131[k]
                   + f_0 * gi_327[k];

        t_216[k] = -2.0 * di_132[k]
                   + f_0 * gi_328[k];

        t_217[k] = -2.0 * di_133[k]
                   + f_0 * gi_329[k];

        t_218[k] = -2.0 * di_134[k]
                   + f_0 * gi_330[k];

        t_219[k] = -2.0 * di_135[k]
                   + f_0 * gi_331[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, di_136, di_137, di_138, di_139, \
                         di_140, gi_332, gi_333, gi_334, gi_335, \
                         gi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -2.0 * di_136[k]
                   + f_0 * gi_332[k];

        t_221[k] = -2.0 * di_137[k]
                   + f_0 * gi_333[k];

        t_222[k] = -2.0 * di_138[k]
                   + f_0 * gi_334[k];

        t_223[k] = -2.0 * di_139[k]
                   + f_0 * gi_335[k];

        t_224[k] = -di_140[k]
                   + f_0 * gi_336[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, di_141, di_142, di_143, di_144, \
                         di_145, gi_337, gi_338, gi_339, gi_340, \
                         gi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -di_141[k]
                   + f_0 * gi_337[k];

        t_226[k] = -di_142[k]
                   + f_0 * gi_338[k];

        t_227[k] = -di_143[k]
                   + f_0 * gi_339[k];

        t_228[k] = -di_144[k]
                   + f_0 * gi_340[k];

        t_229[k] = -di_145[k]
                   + f_0 * gi_341[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, di_146, di_147, di_148, di_149, \
                         di_150, gi_342, gi_343, gi_344, gi_345, \
                         gi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -di_146[k]
                   + f_0 * gi_342[k];

        t_231[k] = -di_147[k]
                   + f_0 * gi_343[k];

        t_232[k] = -di_148[k]
                   + f_0 * gi_344[k];

        t_233[k] = -di_149[k]
                   + f_0 * gi_345[k];

        t_234[k] = -di_150[k]
                   + f_0 * gi_346[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, di_151, di_152, di_153, di_154, \
                         di_155, gi_347, gi_348, gi_349, gi_350, \
                         gi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -di_151[k]
                   + f_0 * gi_347[k];

        t_236[k] = -di_152[k]
                   + f_0 * gi_348[k];

        t_237[k] = -di_153[k]
                   + f_0 * gi_349[k];

        t_238[k] = -di_154[k]
                   + f_0 * gi_350[k];

        t_239[k] = -di_155[k]
                   + f_0 * gi_351[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, di_156, di_157, di_158, di_159, \
                         di_160, gi_352, gi_353, gi_354, gi_355, \
                         gi_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -di_156[k]
                   + f_0 * gi_352[k];

        t_241[k] = -di_157[k]
                   + f_0 * gi_353[k];

        t_242[k] = -di_158[k]
                   + f_0 * gi_354[k];

        t_243[k] = -di_159[k]
                   + f_0 * gi_355[k];

        t_244[k] = -di_160[k]
                   + f_0 * gi_356[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, di_161, di_162, di_163, di_164, \
                         di_165, gi_357, gi_358, gi_359, gi_360, \
                         gi_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -di_161[k]
                   + f_0 * gi_357[k];

        t_246[k] = -di_162[k]
                   + f_0 * gi_358[k];

        t_247[k] = -di_163[k]
                   + f_0 * gi_359[k];

        t_248[k] = -di_164[k]
                   + f_0 * gi_360[k];

        t_249[k] = -di_165[k]
                   + f_0 * gi_361[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, t_255, t_256, di_166, di_167, \
                         gi_362, gi_363, gi_364, gi_365, gi_366, gi_367, \
                         gi_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -di_166[k]
                   + f_0 * gi_362[k];

        t_251[k] = -di_167[k]
                   + f_0 * gi_363[k];

        t_252[k] = f_0 * gi_364[k];

        t_253[k] = f_0 * gi_365[k];

        t_254[k] = f_0 * gi_366[k];

        t_255[k] = f_0 * gi_367[k];

        t_256[k] = f_0 * gi_368[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, t_262, t_263, t_264, gi_369, \
                         gi_370, gi_371, gi_372, gi_373, gi_374, gi_375, \
                         gi_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_0 * gi_369[k];

        t_258[k] = f_0 * gi_370[k];

        t_259[k] = f_0 * gi_371[k];

        t_260[k] = f_0 * gi_372[k];

        t_261[k] = f_0 * gi_373[k];

        t_262[k] = f_0 * gi_374[k];

        t_263[k] = f_0 * gi_375[k];

        t_264[k] = f_0 * gi_376[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, t_271, t_272, gi_377, \
                         gi_378, gi_379, gi_380, gi_381, gi_382, gi_383, \
                         gi_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_0 * gi_377[k];

        t_266[k] = f_0 * gi_378[k];

        t_267[k] = f_0 * gi_379[k];

        t_268[k] = f_0 * gi_380[k];

        t_269[k] = f_0 * gi_381[k];

        t_270[k] = f_0 * gi_382[k];

        t_271[k] = f_0 * gi_383[k];

        t_272[k] = f_0 * gi_384[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, t_277, t_278, t_279, gi_385, gi_386, \
                         gi_387, gi_388, gi_389, gi_390, gi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_0 * gi_385[k];

        t_274[k] = f_0 * gi_386[k];

        t_275[k] = f_0 * gi_387[k];

        t_276[k] = f_0 * gi_388[k];

        t_277[k] = f_0 * gi_389[k];

        t_278[k] = f_0 * gi_390[k];

        t_279[k] = f_0 * gi_391[k];
    }
}

auto
compute_prim_geom_10_fi_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t di, const size_t gi,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fi_electron_repulsion_1_piece0(buffer, target, di, gi, ncols, alpha);

    compute_prim_geom_10_fi_electron_repulsion_1_piece1(buffer, target, di, gi, ncols, alpha);
}

static auto
compute_prim_geom_10_fi_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t di, const size_t gi,
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

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_1 = buffer.data(di + 1);
    const auto *di_2 = buffer.data(di + 2);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_7 = buffer.data(di + 7);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_16 = buffer.data(di + 16);
    const auto *di_17 = buffer.data(di + 17);
    const auto *di_18 = buffer.data(di + 18);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_23 = buffer.data(di + 23);
    const auto *di_24 = buffer.data(di + 24);
    const auto *di_25 = buffer.data(di + 25);
    const auto *di_26 = buffer.data(di + 26);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_28 = buffer.data(di + 28);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_30 = buffer.data(di + 30);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_32 = buffer.data(di + 32);
    const auto *di_33 = buffer.data(di + 33);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_35 = buffer.data(di + 35);
    const auto *di_36 = buffer.data(di + 36);
    const auto *di_37 = buffer.data(di + 37);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_39 = buffer.data(di + 39);
    const auto *di_40 = buffer.data(di + 40);
    const auto *di_41 = buffer.data(di + 41);
    const auto *di_42 = buffer.data(di + 42);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_44 = buffer.data(di + 44);
    const auto *di_45 = buffer.data(di + 45);
    const auto *di_46 = buffer.data(di + 46);
    const auto *di_47 = buffer.data(di + 47);
    const auto *di_48 = buffer.data(di + 48);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_50 = buffer.data(di + 50);
    const auto *di_51 = buffer.data(di + 51);
    const auto *di_52 = buffer.data(di + 52);
    const auto *di_53 = buffer.data(di + 53);
    const auto *di_54 = buffer.data(di + 54);
    const auto *di_55 = buffer.data(di + 55);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_57 = buffer.data(di + 57);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_59 = buffer.data(di + 59);
    const auto *di_60 = buffer.data(di + 60);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_62 = buffer.data(di + 62);
    const auto *di_63 = buffer.data(di + 63);
    const auto *di_64 = buffer.data(di + 64);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_66 = buffer.data(di + 66);
    const auto *di_67 = buffer.data(di + 67);
    const auto *di_68 = buffer.data(di + 68);
    const auto *di_69 = buffer.data(di + 69);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_71 = buffer.data(di + 71);
    const auto *di_72 = buffer.data(di + 72);
    const auto *di_73 = buffer.data(di + 73);
    const auto *di_74 = buffer.data(di + 74);
    const auto *di_75 = buffer.data(di + 75);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_77 = buffer.data(di + 77);
    const auto *di_78 = buffer.data(di + 78);
    const auto *di_79 = buffer.data(di + 79);
    const auto *di_80 = buffer.data(di + 80);
    const auto *di_81 = buffer.data(di + 81);
    const auto *di_82 = buffer.data(di + 82);
    const auto *di_83 = buffer.data(di + 83);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, gi_56, gi_57, gi_58, gi_59, \
                         gi_60, gi_61, gi_62, gi_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gi_56[k];

        t_1[k] = f_0 * gi_57[k];

        t_2[k] = f_0 * gi_58[k];

        t_3[k] = f_0 * gi_59[k];

        t_4[k] = f_0 * gi_60[k];

        t_5[k] = f_0 * gi_61[k];

        t_6[k] = f_0 * gi_62[k];

        t_7[k] = f_0 * gi_63[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, gi_64, gi_65, gi_66, \
                         gi_67, gi_68, gi_69, gi_70, gi_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gi_64[k];

        t_9[k] = f_0 * gi_65[k];

        t_10[k] = f_0 * gi_66[k];

        t_11[k] = f_0 * gi_67[k];

        t_12[k] = f_0 * gi_68[k];

        t_13[k] = f_0 * gi_69[k];

        t_14[k] = f_0 * gi_70[k];

        t_15[k] = f_0 * gi_71[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, gi_72, gi_73, gi_74, \
                         gi_75, gi_76, gi_77, gi_78, gi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * gi_72[k];

        t_17[k] = f_0 * gi_73[k];

        t_18[k] = f_0 * gi_74[k];

        t_19[k] = f_0 * gi_75[k];

        t_20[k] = f_0 * gi_76[k];

        t_21[k] = f_0 * gi_77[k];

        t_22[k] = f_0 * gi_78[k];

        t_23[k] = f_0 * gi_79[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, gi_80, gi_81, gi_82, \
                         gi_83, gi_112, gi_113, gi_114, gi_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * gi_80[k];

        t_25[k] = f_0 * gi_81[k];

        t_26[k] = f_0 * gi_82[k];

        t_27[k] = f_0 * gi_83[k];

        t_28[k] = f_0 * gi_112[k];

        t_29[k] = f_0 * gi_113[k];

        t_30[k] = f_0 * gi_114[k];

        t_31[k] = f_0 * gi_115[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, gi_116, gi_117, \
                         gi_118, gi_119, gi_120, gi_121, gi_122, \
                         gi_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * gi_116[k];

        t_33[k] = f_0 * gi_117[k];

        t_34[k] = f_0 * gi_118[k];

        t_35[k] = f_0 * gi_119[k];

        t_36[k] = f_0 * gi_120[k];

        t_37[k] = f_0 * gi_121[k];

        t_38[k] = f_0 * gi_122[k];

        t_39[k] = f_0 * gi_123[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, t_46, t_47, gi_124, gi_125, \
                         gi_126, gi_127, gi_128, gi_129, gi_130, \
                         gi_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * gi_124[k];

        t_41[k] = f_0 * gi_125[k];

        t_42[k] = f_0 * gi_126[k];

        t_43[k] = f_0 * gi_127[k];

        t_44[k] = f_0 * gi_128[k];

        t_45[k] = f_0 * gi_129[k];

        t_46[k] = f_0 * gi_130[k];

        t_47[k] = f_0 * gi_131[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, t_54, t_55, gi_132, gi_133, \
                         gi_134, gi_135, gi_136, gi_137, gi_138, \
                         gi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * gi_132[k];

        t_49[k] = f_0 * gi_133[k];

        t_50[k] = f_0 * gi_134[k];

        t_51[k] = f_0 * gi_135[k];

        t_52[k] = f_0 * gi_136[k];

        t_53[k] = f_0 * gi_137[k];

        t_54[k] = f_0 * gi_138[k];

        t_55[k] = f_0 * gi_139[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, di_0, di_1, di_2, di_3, di_4, gi_140, \
                         gi_141, gi_142, gi_143, gi_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -di_0[k]
                  + f_0 * gi_140[k];

        t_57[k] = -di_1[k]
                  + f_0 * gi_141[k];

        t_58[k] = -di_2[k]
                  + f_0 * gi_142[k];

        t_59[k] = -di_3[k]
                  + f_0 * gi_143[k];

        t_60[k] = -di_4[k]
                  + f_0 * gi_144[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, di_5, di_6, di_7, di_8, di_9, gi_145, \
                         gi_146, gi_147, gi_148, gi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -di_5[k]
                  + f_0 * gi_145[k];

        t_62[k] = -di_6[k]
                  + f_0 * gi_146[k];

        t_63[k] = -di_7[k]
                  + f_0 * gi_147[k];

        t_64[k] = -di_8[k]
                  + f_0 * gi_148[k];

        t_65[k] = -di_9[k]
                  + f_0 * gi_149[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, di_10, di_11, di_12, di_13, di_14, \
                         gi_150, gi_151, gi_152, gi_153, gi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -di_10[k]
                  + f_0 * gi_150[k];

        t_67[k] = -di_11[k]
                  + f_0 * gi_151[k];

        t_68[k] = -di_12[k]
                  + f_0 * gi_152[k];

        t_69[k] = -di_13[k]
                  + f_0 * gi_153[k];

        t_70[k] = -di_14[k]
                  + f_0 * gi_154[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, di_15, di_16, di_17, di_18, di_19, \
                         gi_155, gi_156, gi_157, gi_158, gi_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -di_15[k]
                  + f_0 * gi_155[k];

        t_72[k] = -di_16[k]
                  + f_0 * gi_156[k];

        t_73[k] = -di_17[k]
                  + f_0 * gi_157[k];

        t_74[k] = -di_18[k]
                  + f_0 * gi_158[k];

        t_75[k] = -di_19[k]
                  + f_0 * gi_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, di_20, di_21, di_22, di_23, di_24, \
                         gi_160, gi_161, gi_162, gi_163, gi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -di_20[k]
                  + f_0 * gi_160[k];

        t_77[k] = -di_21[k]
                  + f_0 * gi_161[k];

        t_78[k] = -di_22[k]
                  + f_0 * gi_162[k];

        t_79[k] = -di_23[k]
                  + f_0 * gi_163[k];

        t_80[k] = -di_24[k]
                  + f_0 * gi_164[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, di_25, di_26, di_27, gi_165, \
                         gi_166, gi_167, gi_196, gi_197, gi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -di_25[k]
                  + f_0 * gi_165[k];

        t_82[k] = -di_26[k]
                  + f_0 * gi_166[k];

        t_83[k] = -di_27[k]
                  + f_0 * gi_167[k];

        t_84[k] = f_0 * gi_196[k];

        t_85[k] = f_0 * gi_197[k];

        t_86[k] = f_0 * gi_198[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, t_92, t_93, t_94, gi_199, gi_200, \
                         gi_201, gi_202, gi_203, gi_204, gi_205, \
                         gi_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_0 * gi_199[k];

        t_88[k] = f_0 * gi_200[k];

        t_89[k] = f_0 * gi_201[k];

        t_90[k] = f_0 * gi_202[k];

        t_91[k] = f_0 * gi_203[k];

        t_92[k] = f_0 * gi_204[k];

        t_93[k] = f_0 * gi_205[k];

        t_94[k] = f_0 * gi_206[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, t_101, t_102, gi_207, gi_208, \
                         gi_209, gi_210, gi_211, gi_212, gi_213, \
                         gi_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * gi_207[k];

        t_96[k] = f_0 * gi_208[k];

        t_97[k] = f_0 * gi_209[k];

        t_98[k] = f_0 * gi_210[k];

        t_99[k] = f_0 * gi_211[k];

        t_100[k] = f_0 * gi_212[k];

        t_101[k] = f_0 * gi_213[k];

        t_102[k] = f_0 * gi_214[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, t_109, t_110, gi_215, \
                         gi_216, gi_217, gi_218, gi_219, gi_220, gi_221, \
                         gi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_0 * gi_215[k];

        t_104[k] = f_0 * gi_216[k];

        t_105[k] = f_0 * gi_217[k];

        t_106[k] = f_0 * gi_218[k];

        t_107[k] = f_0 * gi_219[k];

        t_108[k] = f_0 * gi_220[k];

        t_109[k] = f_0 * gi_221[k];

        t_110[k] = f_0 * gi_222[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, di_28, di_29, di_30, di_31, \
                         gi_223, gi_224, gi_225, gi_226, gi_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * gi_223[k];

        t_112[k] = -di_28[k]
                   + f_0 * gi_224[k];

        t_113[k] = -di_29[k]
                   + f_0 * gi_225[k];

        t_114[k] = -di_30[k]
                   + f_0 * gi_226[k];

        t_115[k] = -di_31[k]
                   + f_0 * gi_227[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, di_32, di_33, di_34, di_35, di_36, \
                         gi_228, gi_229, gi_230, gi_231, gi_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -di_32[k]
                   + f_0 * gi_228[k];

        t_117[k] = -di_33[k]
                   + f_0 * gi_229[k];

        t_118[k] = -di_34[k]
                   + f_0 * gi_230[k];

        t_119[k] = -di_35[k]
                   + f_0 * gi_231[k];

        t_120[k] = -di_36[k]
                   + f_0 * gi_232[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, di_37, di_38, di_39, di_40, di_41, \
                         gi_233, gi_234, gi_235, gi_236, gi_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -di_37[k]
                   + f_0 * gi_233[k];

        t_122[k] = -di_38[k]
                   + f_0 * gi_234[k];

        t_123[k] = -di_39[k]
                   + f_0 * gi_235[k];

        t_124[k] = -di_40[k]
                   + f_0 * gi_236[k];

        t_125[k] = -di_41[k]
                   + f_0 * gi_237[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, di_42, di_43, di_44, di_45, di_46, \
                         gi_238, gi_239, gi_240, gi_241, gi_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -di_42[k]
                   + f_0 * gi_238[k];

        t_127[k] = -di_43[k]
                   + f_0 * gi_239[k];

        t_128[k] = -di_44[k]
                   + f_0 * gi_240[k];

        t_129[k] = -di_45[k]
                   + f_0 * gi_241[k];

        t_130[k] = -di_46[k]
                   + f_0 * gi_242[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, di_47, di_48, di_49, di_50, di_51, \
                         gi_243, gi_244, gi_245, gi_246, gi_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -di_47[k]
                   + f_0 * gi_243[k];

        t_132[k] = -di_48[k]
                   + f_0 * gi_244[k];

        t_133[k] = -di_49[k]
                   + f_0 * gi_245[k];

        t_134[k] = -di_50[k]
                   + f_0 * gi_246[k];

        t_135[k] = -di_51[k]
                   + f_0 * gi_247[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, di_52, di_53, di_54, di_55, di_56, \
                         gi_248, gi_249, gi_250, gi_251, gi_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -di_52[k]
                   + f_0 * gi_248[k];

        t_137[k] = -di_53[k]
                   + f_0 * gi_249[k];

        t_138[k] = -di_54[k]
                   + f_0 * gi_250[k];

        t_139[k] = -di_55[k]
                   + f_0 * gi_251[k];

        t_140[k] = -2.0 * di_56[k]
                   + f_0 * gi_252[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, di_57, di_58, di_59, di_60, di_61, \
                         gi_253, gi_254, gi_255, gi_256, gi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = -2.0 * di_57[k]
                   + f_0 * gi_253[k];

        t_142[k] = -2.0 * di_58[k]
                   + f_0 * gi_254[k];

        t_143[k] = -2.0 * di_59[k]
                   + f_0 * gi_255[k];

        t_144[k] = -2.0 * di_60[k]
                   + f_0 * gi_256[k];

        t_145[k] = -2.0 * di_61[k]
                   + f_0 * gi_257[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, di_62, di_63, di_64, di_65, di_66, \
                         gi_258, gi_259, gi_260, gi_261, gi_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -2.0 * di_62[k]
                   + f_0 * gi_258[k];

        t_147[k] = -2.0 * di_63[k]
                   + f_0 * gi_259[k];

        t_148[k] = -2.0 * di_64[k]
                   + f_0 * gi_260[k];

        t_149[k] = -2.0 * di_65[k]
                   + f_0 * gi_261[k];

        t_150[k] = -2.0 * di_66[k]
                   + f_0 * gi_262[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, di_67, di_68, di_69, di_70, di_71, \
                         gi_263, gi_264, gi_265, gi_266, gi_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = -2.0 * di_67[k]
                   + f_0 * gi_263[k];

        t_152[k] = -2.0 * di_68[k]
                   + f_0 * gi_264[k];

        t_153[k] = -2.0 * di_69[k]
                   + f_0 * gi_265[k];

        t_154[k] = -2.0 * di_70[k]
                   + f_0 * gi_266[k];

        t_155[k] = -2.0 * di_71[k]
                   + f_0 * gi_267[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, di_72, di_73, di_74, di_75, di_76, \
                         gi_268, gi_269, gi_270, gi_271, gi_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = -2.0 * di_72[k]
                   + f_0 * gi_268[k];

        t_157[k] = -2.0 * di_73[k]
                   + f_0 * gi_269[k];

        t_158[k] = -2.0 * di_74[k]
                   + f_0 * gi_270[k];

        t_159[k] = -2.0 * di_75[k]
                   + f_0 * gi_271[k];

        t_160[k] = -2.0 * di_76[k]
                   + f_0 * gi_272[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, di_77, di_78, di_79, di_80, di_81, \
                         gi_273, gi_274, gi_275, gi_276, gi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = -2.0 * di_77[k]
                   + f_0 * gi_273[k];

        t_162[k] = -2.0 * di_78[k]
                   + f_0 * gi_274[k];

        t_163[k] = -2.0 * di_79[k]
                   + f_0 * gi_275[k];

        t_164[k] = -2.0 * di_80[k]
                   + f_0 * gi_276[k];

        t_165[k] = -2.0 * di_81[k]
                   + f_0 * gi_277[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, t_172, di_82, di_83, \
                         gi_278, gi_279, gi_308, gi_309, gi_310, gi_311, \
                         gi_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = -2.0 * di_82[k]
                   + f_0 * gi_278[k];

        t_167[k] = -2.0 * di_83[k]
                   + f_0 * gi_279[k];

        t_168[k] = f_0 * gi_308[k];

        t_169[k] = f_0 * gi_309[k];

        t_170[k] = f_0 * gi_310[k];

        t_171[k] = f_0 * gi_311[k];

        t_172[k] = f_0 * gi_312[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, t_179, t_180, gi_313, \
                         gi_314, gi_315, gi_316, gi_317, gi_318, gi_319, \
                         gi_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_0 * gi_313[k];

        t_174[k] = f_0 * gi_314[k];

        t_175[k] = f_0 * gi_315[k];

        t_176[k] = f_0 * gi_316[k];

        t_177[k] = f_0 * gi_317[k];

        t_178[k] = f_0 * gi_318[k];

        t_179[k] = f_0 * gi_319[k];

        t_180[k] = f_0 * gi_320[k];
    }
}

static auto
compute_prim_geom_10_fi_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t di, const size_t gi,
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

    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_86 = buffer.data(di + 86);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_88 = buffer.data(di + 88);
    const auto *di_89 = buffer.data(di + 89);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_91 = buffer.data(di + 91);
    const auto *di_92 = buffer.data(di + 92);
    const auto *di_93 = buffer.data(di + 93);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_95 = buffer.data(di + 95);
    const auto *di_96 = buffer.data(di + 96);
    const auto *di_97 = buffer.data(di + 97);
    const auto *di_98 = buffer.data(di + 98);
    const auto *di_99 = buffer.data(di + 99);
    const auto *di_100 = buffer.data(di + 100);
    const auto *di_101 = buffer.data(di + 101);
    const auto *di_102 = buffer.data(di + 102);
    const auto *di_103 = buffer.data(di + 103);
    const auto *di_104 = buffer.data(di + 104);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_106 = buffer.data(di + 106);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_112 = buffer.data(di + 112);
    const auto *di_113 = buffer.data(di + 113);
    const auto *di_114 = buffer.data(di + 114);
    const auto *di_115 = buffer.data(di + 115);
    const auto *di_116 = buffer.data(di + 116);
    const auto *di_117 = buffer.data(di + 117);
    const auto *di_118 = buffer.data(di + 118);
    const auto *di_119 = buffer.data(di + 119);
    const auto *di_120 = buffer.data(di + 120);
    const auto *di_121 = buffer.data(di + 121);
    const auto *di_122 = buffer.data(di + 122);
    const auto *di_123 = buffer.data(di + 123);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_125 = buffer.data(di + 125);
    const auto *di_126 = buffer.data(di + 126);
    const auto *di_127 = buffer.data(di + 127);
    const auto *di_128 = buffer.data(di + 128);
    const auto *di_129 = buffer.data(di + 129);
    const auto *di_130 = buffer.data(di + 130);
    const auto *di_131 = buffer.data(di + 131);
    const auto *di_132 = buffer.data(di + 132);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_141 = buffer.data(di + 141);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_143 = buffer.data(di + 143);
    const auto *di_144 = buffer.data(di + 144);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_146 = buffer.data(di + 146);
    const auto *di_147 = buffer.data(di + 147);
    const auto *di_148 = buffer.data(di + 148);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_150 = buffer.data(di + 150);
    const auto *di_151 = buffer.data(di + 151);
    const auto *di_152 = buffer.data(di + 152);
    const auto *di_153 = buffer.data(di + 153);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_155 = buffer.data(di + 155);
    const auto *di_156 = buffer.data(di + 156);
    const auto *di_157 = buffer.data(di + 157);
    const auto *di_158 = buffer.data(di + 158);
    const auto *di_159 = buffer.data(di + 159);
    const auto *di_160 = buffer.data(di + 160);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_166 = buffer.data(di + 166);
    const auto *di_167 = buffer.data(di + 167);

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

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, t_186, t_187, t_188, gi_321, \
                         gi_322, gi_323, gi_324, gi_325, gi_326, gi_327, \
                         gi_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_0 * gi_321[k];

        t_182[k] = f_0 * gi_322[k];

        t_183[k] = f_0 * gi_323[k];

        t_184[k] = f_0 * gi_324[k];

        t_185[k] = f_0 * gi_325[k];

        t_186[k] = f_0 * gi_326[k];

        t_187[k] = f_0 * gi_327[k];

        t_188[k] = f_0 * gi_328[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, t_194, t_195, gi_329, gi_330, \
                         gi_331, gi_332, gi_333, gi_334, gi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * gi_329[k];

        t_190[k] = f_0 * gi_330[k];

        t_191[k] = f_0 * gi_331[k];

        t_192[k] = f_0 * gi_332[k];

        t_193[k] = f_0 * gi_333[k];

        t_194[k] = f_0 * gi_334[k];

        t_195[k] = f_0 * gi_335[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, di_84, di_85, di_86, di_87, di_88, \
                         gi_336, gi_337, gi_338, gi_339, gi_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -di_84[k]
                   + f_0 * gi_336[k];

        t_197[k] = -di_85[k]
                   + f_0 * gi_337[k];

        t_198[k] = -di_86[k]
                   + f_0 * gi_338[k];

        t_199[k] = -di_87[k]
                   + f_0 * gi_339[k];

        t_200[k] = -di_88[k]
                   + f_0 * gi_340[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, di_89, di_90, di_91, di_92, di_93, \
                         gi_341, gi_342, gi_343, gi_344, gi_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -di_89[k]
                   + f_0 * gi_341[k];

        t_202[k] = -di_90[k]
                   + f_0 * gi_342[k];

        t_203[k] = -di_91[k]
                   + f_0 * gi_343[k];

        t_204[k] = -di_92[k]
                   + f_0 * gi_344[k];

        t_205[k] = -di_93[k]
                   + f_0 * gi_345[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, di_94, di_95, di_96, di_97, di_98, \
                         gi_346, gi_347, gi_348, gi_349, gi_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = -di_94[k]
                   + f_0 * gi_346[k];

        t_207[k] = -di_95[k]
                   + f_0 * gi_347[k];

        t_208[k] = -di_96[k]
                   + f_0 * gi_348[k];

        t_209[k] = -di_97[k]
                   + f_0 * gi_349[k];

        t_210[k] = -di_98[k]
                   + f_0 * gi_350[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, di_99, di_100, di_101, di_102, \
                         di_103, gi_351, gi_352, gi_353, gi_354, \
                         gi_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -di_99[k]
                   + f_0 * gi_351[k];

        t_212[k] = -di_100[k]
                   + f_0 * gi_352[k];

        t_213[k] = -di_101[k]
                   + f_0 * gi_353[k];

        t_214[k] = -di_102[k]
                   + f_0 * gi_354[k];

        t_215[k] = -di_103[k]
                   + f_0 * gi_355[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, di_104, di_105, di_106, di_107, \
                         di_108, gi_356, gi_357, gi_358, gi_359, \
                         gi_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = -di_104[k]
                   + f_0 * gi_356[k];

        t_217[k] = -di_105[k]
                   + f_0 * gi_357[k];

        t_218[k] = -di_106[k]
                   + f_0 * gi_358[k];

        t_219[k] = -di_107[k]
                   + f_0 * gi_359[k];

        t_220[k] = -di_108[k]
                   + f_0 * gi_360[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, di_109, di_110, di_111, di_112, \
                         di_113, gi_361, gi_362, gi_363, gi_364, \
                         gi_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = -di_109[k]
                   + f_0 * gi_361[k];

        t_222[k] = -di_110[k]
                   + f_0 * gi_362[k];

        t_223[k] = -di_111[k]
                   + f_0 * gi_363[k];

        t_224[k] = -2.0 * di_112[k]
                   + f_0 * gi_364[k];

        t_225[k] = -2.0 * di_113[k]
                   + f_0 * gi_365[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, t_230, di_114, di_115, di_116, di_117, \
                         di_118, gi_366, gi_367, gi_368, gi_369, \
                         gi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = -2.0 * di_114[k]
                   + f_0 * gi_366[k];

        t_227[k] = -2.0 * di_115[k]
                   + f_0 * gi_367[k];

        t_228[k] = -2.0 * di_116[k]
                   + f_0 * gi_368[k];

        t_229[k] = -2.0 * di_117[k]
                   + f_0 * gi_369[k];

        t_230[k] = -2.0 * di_118[k]
                   + f_0 * gi_370[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, di_119, di_120, di_121, di_122, \
                         di_123, gi_371, gi_372, gi_373, gi_374, \
                         gi_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = -2.0 * di_119[k]
                   + f_0 * gi_371[k];

        t_232[k] = -2.0 * di_120[k]
                   + f_0 * gi_372[k];

        t_233[k] = -2.0 * di_121[k]
                   + f_0 * gi_373[k];

        t_234[k] = -2.0 * di_122[k]
                   + f_0 * gi_374[k];

        t_235[k] = -2.0 * di_123[k]
                   + f_0 * gi_375[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, di_124, di_125, di_126, di_127, \
                         di_128, gi_376, gi_377, gi_378, gi_379, \
                         gi_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = -2.0 * di_124[k]
                   + f_0 * gi_376[k];

        t_237[k] = -2.0 * di_125[k]
                   + f_0 * gi_377[k];

        t_238[k] = -2.0 * di_126[k]
                   + f_0 * gi_378[k];

        t_239[k] = -2.0 * di_127[k]
                   + f_0 * gi_379[k];

        t_240[k] = -2.0 * di_128[k]
                   + f_0 * gi_380[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, di_129, di_130, di_131, di_132, \
                         di_133, gi_381, gi_382, gi_383, gi_384, \
                         gi_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = -2.0 * di_129[k]
                   + f_0 * gi_381[k];

        t_242[k] = -2.0 * di_130[k]
                   + f_0 * gi_382[k];

        t_243[k] = -2.0 * di_131[k]
                   + f_0 * gi_383[k];

        t_244[k] = -2.0 * di_132[k]
                   + f_0 * gi_384[k];

        t_245[k] = -2.0 * di_133[k]
                   + f_0 * gi_385[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, di_134, di_135, di_136, di_137, \
                         di_138, gi_386, gi_387, gi_388, gi_389, \
                         gi_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = -2.0 * di_134[k]
                   + f_0 * gi_386[k];

        t_247[k] = -2.0 * di_135[k]
                   + f_0 * gi_387[k];

        t_248[k] = -2.0 * di_136[k]
                   + f_0 * gi_388[k];

        t_249[k] = -2.0 * di_137[k]
                   + f_0 * gi_389[k];

        t_250[k] = -2.0 * di_138[k]
                   + f_0 * gi_390[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, di_139, di_140, di_141, di_142, \
                         di_143, gi_391, gi_392, gi_393, gi_394, \
                         gi_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = -2.0 * di_139[k]
                   + f_0 * gi_391[k];

        t_252[k] = -3.0 * di_140[k]
                   + f_0 * gi_392[k];

        t_253[k] = -3.0 * di_141[k]
                   + f_0 * gi_393[k];

        t_254[k] = -3.0 * di_142[k]
                   + f_0 * gi_394[k];

        t_255[k] = -3.0 * di_143[k]
                   + f_0 * gi_395[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, di_144, di_145, di_146, di_147, \
                         di_148, gi_396, gi_397, gi_398, gi_399, \
                         gi_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = -3.0 * di_144[k]
                   + f_0 * gi_396[k];

        t_257[k] = -3.0 * di_145[k]
                   + f_0 * gi_397[k];

        t_258[k] = -3.0 * di_146[k]
                   + f_0 * gi_398[k];

        t_259[k] = -3.0 * di_147[k]
                   + f_0 * gi_399[k];

        t_260[k] = -3.0 * di_148[k]
                   + f_0 * gi_400[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, di_149, di_150, di_151, di_152, \
                         di_153, gi_401, gi_402, gi_403, gi_404, \
                         gi_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = -3.0 * di_149[k]
                   + f_0 * gi_401[k];

        t_262[k] = -3.0 * di_150[k]
                   + f_0 * gi_402[k];

        t_263[k] = -3.0 * di_151[k]
                   + f_0 * gi_403[k];

        t_264[k] = -3.0 * di_152[k]
                   + f_0 * gi_404[k];

        t_265[k] = -3.0 * di_153[k]
                   + f_0 * gi_405[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, di_154, di_155, di_156, di_157, \
                         di_158, gi_406, gi_407, gi_408, gi_409, \
                         gi_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = -3.0 * di_154[k]
                   + f_0 * gi_406[k];

        t_267[k] = -3.0 * di_155[k]
                   + f_0 * gi_407[k];

        t_268[k] = -3.0 * di_156[k]
                   + f_0 * gi_408[k];

        t_269[k] = -3.0 * di_157[k]
                   + f_0 * gi_409[k];

        t_270[k] = -3.0 * di_158[k]
                   + f_0 * gi_410[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, di_159, di_160, di_161, di_162, \
                         di_163, gi_411, gi_412, gi_413, gi_414, \
                         gi_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = -3.0 * di_159[k]
                   + f_0 * gi_411[k];

        t_272[k] = -3.0 * di_160[k]
                   + f_0 * gi_412[k];

        t_273[k] = -3.0 * di_161[k]
                   + f_0 * gi_413[k];

        t_274[k] = -3.0 * di_162[k]
                   + f_0 * gi_414[k];

        t_275[k] = -3.0 * di_163[k]
                   + f_0 * gi_415[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, di_164, di_165, di_166, di_167, gi_416, \
                         gi_417, gi_418, gi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = -3.0 * di_164[k]
                   + f_0 * gi_416[k];

        t_277[k] = -3.0 * di_165[k]
                   + f_0 * gi_417[k];

        t_278[k] = -3.0 * di_166[k]
                   + f_0 * gi_418[k];

        t_279[k] = -3.0 * di_167[k]
                   + f_0 * gi_419[k];
    }
}

auto
compute_prim_geom_10_fi_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t di, const size_t gi,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_fi_electron_repulsion_2_piece0(buffer, target, di, gi, ncols, alpha);

    compute_prim_geom_10_fi_electron_repulsion_2_piece1(buffer, target, di, gi, ncols, alpha);
}

}  // namespace simdt2ceri
