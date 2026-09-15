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


#include "SimdElectronRepulsionGeom10VrrRecLD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_ld_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kd, const size_t md,
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

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_92 = buffer.data(kd + 92);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);

    const auto *md_0 = buffer.data(md + 0);
    const auto *md_1 = buffer.data(md + 1);
    const auto *md_2 = buffer.data(md + 2);
    const auto *md_3 = buffer.data(md + 3);
    const auto *md_4 = buffer.data(md + 4);
    const auto *md_5 = buffer.data(md + 5);
    const auto *md_6 = buffer.data(md + 6);
    const auto *md_7 = buffer.data(md + 7);
    const auto *md_8 = buffer.data(md + 8);
    const auto *md_9 = buffer.data(md + 9);
    const auto *md_10 = buffer.data(md + 10);
    const auto *md_11 = buffer.data(md + 11);
    const auto *md_12 = buffer.data(md + 12);
    const auto *md_13 = buffer.data(md + 13);
    const auto *md_14 = buffer.data(md + 14);
    const auto *md_15 = buffer.data(md + 15);
    const auto *md_16 = buffer.data(md + 16);
    const auto *md_17 = buffer.data(md + 17);
    const auto *md_18 = buffer.data(md + 18);
    const auto *md_19 = buffer.data(md + 19);
    const auto *md_20 = buffer.data(md + 20);
    const auto *md_21 = buffer.data(md + 21);
    const auto *md_22 = buffer.data(md + 22);
    const auto *md_23 = buffer.data(md + 23);
    const auto *md_24 = buffer.data(md + 24);
    const auto *md_25 = buffer.data(md + 25);
    const auto *md_26 = buffer.data(md + 26);
    const auto *md_27 = buffer.data(md + 27);
    const auto *md_28 = buffer.data(md + 28);
    const auto *md_29 = buffer.data(md + 29);
    const auto *md_30 = buffer.data(md + 30);
    const auto *md_31 = buffer.data(md + 31);
    const auto *md_32 = buffer.data(md + 32);
    const auto *md_33 = buffer.data(md + 33);
    const auto *md_34 = buffer.data(md + 34);
    const auto *md_35 = buffer.data(md + 35);
    const auto *md_36 = buffer.data(md + 36);
    const auto *md_37 = buffer.data(md + 37);
    const auto *md_38 = buffer.data(md + 38);
    const auto *md_39 = buffer.data(md + 39);
    const auto *md_40 = buffer.data(md + 40);
    const auto *md_41 = buffer.data(md + 41);
    const auto *md_42 = buffer.data(md + 42);
    const auto *md_43 = buffer.data(md + 43);
    const auto *md_44 = buffer.data(md + 44);
    const auto *md_45 = buffer.data(md + 45);
    const auto *md_46 = buffer.data(md + 46);
    const auto *md_47 = buffer.data(md + 47);
    const auto *md_48 = buffer.data(md + 48);
    const auto *md_49 = buffer.data(md + 49);
    const auto *md_50 = buffer.data(md + 50);
    const auto *md_51 = buffer.data(md + 51);
    const auto *md_52 = buffer.data(md + 52);
    const auto *md_53 = buffer.data(md + 53);
    const auto *md_54 = buffer.data(md + 54);
    const auto *md_55 = buffer.data(md + 55);
    const auto *md_56 = buffer.data(md + 56);
    const auto *md_57 = buffer.data(md + 57);
    const auto *md_58 = buffer.data(md + 58);
    const auto *md_59 = buffer.data(md + 59);
    const auto *md_60 = buffer.data(md + 60);
    const auto *md_61 = buffer.data(md + 61);
    const auto *md_62 = buffer.data(md + 62);
    const auto *md_63 = buffer.data(md + 63);
    const auto *md_64 = buffer.data(md + 64);
    const auto *md_65 = buffer.data(md + 65);
    const auto *md_66 = buffer.data(md + 66);
    const auto *md_67 = buffer.data(md + 67);
    const auto *md_68 = buffer.data(md + 68);
    const auto *md_69 = buffer.data(md + 69);
    const auto *md_70 = buffer.data(md + 70);
    const auto *md_71 = buffer.data(md + 71);
    const auto *md_72 = buffer.data(md + 72);
    const auto *md_73 = buffer.data(md + 73);
    const auto *md_74 = buffer.data(md + 74);
    const auto *md_75 = buffer.data(md + 75);
    const auto *md_76 = buffer.data(md + 76);
    const auto *md_77 = buffer.data(md + 77);
    const auto *md_78 = buffer.data(md + 78);
    const auto *md_79 = buffer.data(md + 79);
    const auto *md_80 = buffer.data(md + 80);
    const auto *md_81 = buffer.data(md + 81);
    const auto *md_82 = buffer.data(md + 82);
    const auto *md_83 = buffer.data(md + 83);
    const auto *md_84 = buffer.data(md + 84);
    const auto *md_85 = buffer.data(md + 85);
    const auto *md_86 = buffer.data(md + 86);
    const auto *md_87 = buffer.data(md + 87);
    const auto *md_88 = buffer.data(md + 88);
    const auto *md_89 = buffer.data(md + 89);
    const auto *md_90 = buffer.data(md + 90);
    const auto *md_91 = buffer.data(md + 91);
    const auto *md_92 = buffer.data(md + 92);
    const auto *md_93 = buffer.data(md + 93);
    const auto *md_94 = buffer.data(md + 94);
    const auto *md_95 = buffer.data(md + 95);
    const auto *md_96 = buffer.data(md + 96);
    const auto *md_97 = buffer.data(md + 97);
    const auto *md_98 = buffer.data(md + 98);
    const auto *md_99 = buffer.data(md + 99);
    const auto *md_100 = buffer.data(md + 100);
    const auto *md_101 = buffer.data(md + 101);
    const auto *md_102 = buffer.data(md + 102);
    const auto *md_103 = buffer.data(md + 103);
    const auto *md_104 = buffer.data(md + 104);
    const auto *md_105 = buffer.data(md + 105);
    const auto *md_106 = buffer.data(md + 106);
    const auto *md_107 = buffer.data(md + 107);
    const auto *md_108 = buffer.data(md + 108);
    const auto *md_109 = buffer.data(md + 109);
    const auto *md_110 = buffer.data(md + 110);
    const auto *md_111 = buffer.data(md + 111);
    const auto *md_112 = buffer.data(md + 112);
    const auto *md_113 = buffer.data(md + 113);
    const auto *md_114 = buffer.data(md + 114);
    const auto *md_115 = buffer.data(md + 115);
    const auto *md_116 = buffer.data(md + 116);
    const auto *md_117 = buffer.data(md + 117);
    const auto *md_118 = buffer.data(md + 118);
    const auto *md_119 = buffer.data(md + 119);
    const auto *md_120 = buffer.data(md + 120);
    const auto *md_121 = buffer.data(md + 121);
    const auto *md_122 = buffer.data(md + 122);
    const auto *md_123 = buffer.data(md + 123);
    const auto *md_124 = buffer.data(md + 124);
    const auto *md_125 = buffer.data(md + 125);
    const auto *md_126 = buffer.data(md + 126);
    const auto *md_127 = buffer.data(md + 127);
    const auto *md_128 = buffer.data(md + 128);
    const auto *md_129 = buffer.data(md + 129);
    const auto *md_130 = buffer.data(md + 130);
    const auto *md_131 = buffer.data(md + 131);
    const auto *md_132 = buffer.data(md + 132);
    const auto *md_133 = buffer.data(md + 133);
    const auto *md_134 = buffer.data(md + 134);
    const auto *md_135 = buffer.data(md + 135);
    const auto *md_136 = buffer.data(md + 136);
    const auto *md_137 = buffer.data(md + 137);
    const auto *md_138 = buffer.data(md + 138);
    const auto *md_139 = buffer.data(md + 139);
    const auto *md_140 = buffer.data(md + 140);
    const auto *md_141 = buffer.data(md + 141);
    const auto *md_142 = buffer.data(md + 142);
    const auto *md_143 = buffer.data(md + 143);
    const auto *md_144 = buffer.data(md + 144);
    const auto *md_145 = buffer.data(md + 145);
    const auto *md_146 = buffer.data(md + 146);
    const auto *md_147 = buffer.data(md + 147);
    const auto *md_148 = buffer.data(md + 148);
    const auto *md_149 = buffer.data(md + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, kd_0, kd_1, kd_2, kd_3, kd_4, md_0, md_1, \
                         md_2, md_3, md_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -8.0 * kd_0[k]
                 + f_0 * md_0[k];

        t_1[k] = -8.0 * kd_1[k]
                 + f_0 * md_1[k];

        t_2[k] = -8.0 * kd_2[k]
                 + f_0 * md_2[k];

        t_3[k] = -8.0 * kd_3[k]
                 + f_0 * md_3[k];

        t_4[k] = -8.0 * kd_4[k]
                 + f_0 * md_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, kd_5, kd_6, kd_7, kd_8, kd_9, md_5, md_6, \
                         md_7, md_8, md_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -8.0 * kd_5[k]
                 + f_0 * md_5[k];

        t_6[k] = -7.0 * kd_6[k]
                 + f_0 * md_6[k];

        t_7[k] = -7.0 * kd_7[k]
                 + f_0 * md_7[k];

        t_8[k] = -7.0 * kd_8[k]
                 + f_0 * md_8[k];

        t_9[k] = -7.0 * kd_9[k]
                 + f_0 * md_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, kd_10, kd_11, kd_12, kd_13, kd_14, \
                         md_10, md_11, md_12, md_13, md_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -7.0 * kd_10[k]
                  + f_0 * md_10[k];

        t_11[k] = -7.0 * kd_11[k]
                  + f_0 * md_11[k];

        t_12[k] = -7.0 * kd_12[k]
                  + f_0 * md_12[k];

        t_13[k] = -7.0 * kd_13[k]
                  + f_0 * md_13[k];

        t_14[k] = -7.0 * kd_14[k]
                  + f_0 * md_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, kd_15, kd_16, kd_17, kd_18, kd_19, \
                         md_15, md_16, md_17, md_18, md_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -7.0 * kd_15[k]
                  + f_0 * md_15[k];

        t_16[k] = -7.0 * kd_16[k]
                  + f_0 * md_16[k];

        t_17[k] = -7.0 * kd_17[k]
                  + f_0 * md_17[k];

        t_18[k] = -6.0 * kd_18[k]
                  + f_0 * md_18[k];

        t_19[k] = -6.0 * kd_19[k]
                  + f_0 * md_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, kd_20, kd_21, kd_22, kd_23, kd_24, \
                         md_20, md_21, md_22, md_23, md_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -6.0 * kd_20[k]
                  + f_0 * md_20[k];

        t_21[k] = -6.0 * kd_21[k]
                  + f_0 * md_21[k];

        t_22[k] = -6.0 * kd_22[k]
                  + f_0 * md_22[k];

        t_23[k] = -6.0 * kd_23[k]
                  + f_0 * md_23[k];

        t_24[k] = -6.0 * kd_24[k]
                  + f_0 * md_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, kd_25, kd_26, kd_27, kd_28, kd_29, \
                         md_25, md_26, md_27, md_28, md_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -6.0 * kd_25[k]
                  + f_0 * md_25[k];

        t_26[k] = -6.0 * kd_26[k]
                  + f_0 * md_26[k];

        t_27[k] = -6.0 * kd_27[k]
                  + f_0 * md_27[k];

        t_28[k] = -6.0 * kd_28[k]
                  + f_0 * md_28[k];

        t_29[k] = -6.0 * kd_29[k]
                  + f_0 * md_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, kd_30, kd_31, kd_32, kd_33, kd_34, \
                         md_30, md_31, md_32, md_33, md_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -6.0 * kd_30[k]
                  + f_0 * md_30[k];

        t_31[k] = -6.0 * kd_31[k]
                  + f_0 * md_31[k];

        t_32[k] = -6.0 * kd_32[k]
                  + f_0 * md_32[k];

        t_33[k] = -6.0 * kd_33[k]
                  + f_0 * md_33[k];

        t_34[k] = -6.0 * kd_34[k]
                  + f_0 * md_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, kd_35, kd_36, kd_37, kd_38, kd_39, \
                         md_35, md_36, md_37, md_38, md_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -6.0 * kd_35[k]
                  + f_0 * md_35[k];

        t_36[k] = -5.0 * kd_36[k]
                  + f_0 * md_36[k];

        t_37[k] = -5.0 * kd_37[k]
                  + f_0 * md_37[k];

        t_38[k] = -5.0 * kd_38[k]
                  + f_0 * md_38[k];

        t_39[k] = -5.0 * kd_39[k]
                  + f_0 * md_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, kd_40, kd_41, kd_42, kd_43, kd_44, \
                         md_40, md_41, md_42, md_43, md_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -5.0 * kd_40[k]
                  + f_0 * md_40[k];

        t_41[k] = -5.0 * kd_41[k]
                  + f_0 * md_41[k];

        t_42[k] = -5.0 * kd_42[k]
                  + f_0 * md_42[k];

        t_43[k] = -5.0 * kd_43[k]
                  + f_0 * md_43[k];

        t_44[k] = -5.0 * kd_44[k]
                  + f_0 * md_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, kd_45, kd_46, kd_47, kd_48, kd_49, \
                         md_45, md_46, md_47, md_48, md_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * kd_45[k]
                  + f_0 * md_45[k];

        t_46[k] = -5.0 * kd_46[k]
                  + f_0 * md_46[k];

        t_47[k] = -5.0 * kd_47[k]
                  + f_0 * md_47[k];

        t_48[k] = -5.0 * kd_48[k]
                  + f_0 * md_48[k];

        t_49[k] = -5.0 * kd_49[k]
                  + f_0 * md_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, kd_50, kd_51, kd_52, kd_53, kd_54, \
                         md_50, md_51, md_52, md_53, md_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -5.0 * kd_50[k]
                  + f_0 * md_50[k];

        t_51[k] = -5.0 * kd_51[k]
                  + f_0 * md_51[k];

        t_52[k] = -5.0 * kd_52[k]
                  + f_0 * md_52[k];

        t_53[k] = -5.0 * kd_53[k]
                  + f_0 * md_53[k];

        t_54[k] = -5.0 * kd_54[k]
                  + f_0 * md_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, kd_55, kd_56, kd_57, kd_58, kd_59, \
                         md_55, md_56, md_57, md_58, md_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -5.0 * kd_55[k]
                  + f_0 * md_55[k];

        t_56[k] = -5.0 * kd_56[k]
                  + f_0 * md_56[k];

        t_57[k] = -5.0 * kd_57[k]
                  + f_0 * md_57[k];

        t_58[k] = -5.0 * kd_58[k]
                  + f_0 * md_58[k];

        t_59[k] = -5.0 * kd_59[k]
                  + f_0 * md_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, kd_60, kd_61, kd_62, kd_63, kd_64, \
                         md_60, md_61, md_62, md_63, md_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -4.0 * kd_60[k]
                  + f_0 * md_60[k];

        t_61[k] = -4.0 * kd_61[k]
                  + f_0 * md_61[k];

        t_62[k] = -4.0 * kd_62[k]
                  + f_0 * md_62[k];

        t_63[k] = -4.0 * kd_63[k]
                  + f_0 * md_63[k];

        t_64[k] = -4.0 * kd_64[k]
                  + f_0 * md_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, kd_65, kd_66, kd_67, kd_68, kd_69, \
                         md_65, md_66, md_67, md_68, md_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -4.0 * kd_65[k]
                  + f_0 * md_65[k];

        t_66[k] = -4.0 * kd_66[k]
                  + f_0 * md_66[k];

        t_67[k] = -4.0 * kd_67[k]
                  + f_0 * md_67[k];

        t_68[k] = -4.0 * kd_68[k]
                  + f_0 * md_68[k];

        t_69[k] = -4.0 * kd_69[k]
                  + f_0 * md_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, kd_70, kd_71, kd_72, kd_73, kd_74, \
                         md_70, md_71, md_72, md_73, md_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -4.0 * kd_70[k]
                  + f_0 * md_70[k];

        t_71[k] = -4.0 * kd_71[k]
                  + f_0 * md_71[k];

        t_72[k] = -4.0 * kd_72[k]
                  + f_0 * md_72[k];

        t_73[k] = -4.0 * kd_73[k]
                  + f_0 * md_73[k];

        t_74[k] = -4.0 * kd_74[k]
                  + f_0 * md_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, kd_75, kd_76, kd_77, kd_78, kd_79, \
                         md_75, md_76, md_77, md_78, md_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -4.0 * kd_75[k]
                  + f_0 * md_75[k];

        t_76[k] = -4.0 * kd_76[k]
                  + f_0 * md_76[k];

        t_77[k] = -4.0 * kd_77[k]
                  + f_0 * md_77[k];

        t_78[k] = -4.0 * kd_78[k]
                  + f_0 * md_78[k];

        t_79[k] = -4.0 * kd_79[k]
                  + f_0 * md_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, kd_80, kd_81, kd_82, kd_83, kd_84, \
                         md_80, md_81, md_82, md_83, md_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -4.0 * kd_80[k]
                  + f_0 * md_80[k];

        t_81[k] = -4.0 * kd_81[k]
                  + f_0 * md_81[k];

        t_82[k] = -4.0 * kd_82[k]
                  + f_0 * md_82[k];

        t_83[k] = -4.0 * kd_83[k]
                  + f_0 * md_83[k];

        t_84[k] = -4.0 * kd_84[k]
                  + f_0 * md_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, kd_85, kd_86, kd_87, kd_88, kd_89, \
                         md_85, md_86, md_87, md_88, md_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * kd_85[k]
                  + f_0 * md_85[k];

        t_86[k] = -4.0 * kd_86[k]
                  + f_0 * md_86[k];

        t_87[k] = -4.0 * kd_87[k]
                  + f_0 * md_87[k];

        t_88[k] = -4.0 * kd_88[k]
                  + f_0 * md_88[k];

        t_89[k] = -4.0 * kd_89[k]
                  + f_0 * md_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, kd_90, kd_91, kd_92, kd_93, kd_94, \
                         md_90, md_91, md_92, md_93, md_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * kd_90[k]
                  + f_0 * md_90[k];

        t_91[k] = -3.0 * kd_91[k]
                  + f_0 * md_91[k];

        t_92[k] = -3.0 * kd_92[k]
                  + f_0 * md_92[k];

        t_93[k] = -3.0 * kd_93[k]
                  + f_0 * md_93[k];

        t_94[k] = -3.0 * kd_94[k]
                  + f_0 * md_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, kd_95, kd_96, kd_97, kd_98, kd_99, \
                         md_95, md_96, md_97, md_98, md_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * kd_95[k]
                  + f_0 * md_95[k];

        t_96[k] = -3.0 * kd_96[k]
                  + f_0 * md_96[k];

        t_97[k] = -3.0 * kd_97[k]
                  + f_0 * md_97[k];

        t_98[k] = -3.0 * kd_98[k]
                  + f_0 * md_98[k];

        t_99[k] = -3.0 * kd_99[k]
                  + f_0 * md_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, kd_100, kd_101, kd_102, kd_103, \
                         kd_104, md_100, md_101, md_102, md_103, \
                         md_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -3.0 * kd_100[k]
                   + f_0 * md_100[k];

        t_101[k] = -3.0 * kd_101[k]
                   + f_0 * md_101[k];

        t_102[k] = -3.0 * kd_102[k]
                   + f_0 * md_102[k];

        t_103[k] = -3.0 * kd_103[k]
                   + f_0 * md_103[k];

        t_104[k] = -3.0 * kd_104[k]
                   + f_0 * md_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, kd_105, kd_106, kd_107, kd_108, \
                         kd_109, md_105, md_106, md_107, md_108, \
                         md_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * kd_105[k]
                   + f_0 * md_105[k];

        t_106[k] = -3.0 * kd_106[k]
                   + f_0 * md_106[k];

        t_107[k] = -3.0 * kd_107[k]
                   + f_0 * md_107[k];

        t_108[k] = -3.0 * kd_108[k]
                   + f_0 * md_108[k];

        t_109[k] = -3.0 * kd_109[k]
                   + f_0 * md_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, kd_110, kd_111, kd_112, kd_113, \
                         kd_114, md_110, md_111, md_112, md_113, \
                         md_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -3.0 * kd_110[k]
                   + f_0 * md_110[k];

        t_111[k] = -3.0 * kd_111[k]
                   + f_0 * md_111[k];

        t_112[k] = -3.0 * kd_112[k]
                   + f_0 * md_112[k];

        t_113[k] = -3.0 * kd_113[k]
                   + f_0 * md_113[k];

        t_114[k] = -3.0 * kd_114[k]
                   + f_0 * md_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, kd_115, kd_116, kd_117, kd_118, \
                         kd_119, md_115, md_116, md_117, md_118, \
                         md_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -3.0 * kd_115[k]
                   + f_0 * md_115[k];

        t_116[k] = -3.0 * kd_116[k]
                   + f_0 * md_116[k];

        t_117[k] = -3.0 * kd_117[k]
                   + f_0 * md_117[k];

        t_118[k] = -3.0 * kd_118[k]
                   + f_0 * md_118[k];

        t_119[k] = -3.0 * kd_119[k]
                   + f_0 * md_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, kd_120, kd_121, kd_122, kd_123, \
                         kd_124, md_120, md_121, md_122, md_123, \
                         md_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -3.0 * kd_120[k]
                   + f_0 * md_120[k];

        t_121[k] = -3.0 * kd_121[k]
                   + f_0 * md_121[k];

        t_122[k] = -3.0 * kd_122[k]
                   + f_0 * md_122[k];

        t_123[k] = -3.0 * kd_123[k]
                   + f_0 * md_123[k];

        t_124[k] = -3.0 * kd_124[k]
                   + f_0 * md_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, kd_125, kd_126, kd_127, kd_128, \
                         kd_129, md_125, md_126, md_127, md_128, \
                         md_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -3.0 * kd_125[k]
                   + f_0 * md_125[k];

        t_126[k] = -2.0 * kd_126[k]
                   + f_0 * md_126[k];

        t_127[k] = -2.0 * kd_127[k]
                   + f_0 * md_127[k];

        t_128[k] = -2.0 * kd_128[k]
                   + f_0 * md_128[k];

        t_129[k] = -2.0 * kd_129[k]
                   + f_0 * md_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, kd_130, kd_131, kd_132, kd_133, \
                         kd_134, md_130, md_131, md_132, md_133, \
                         md_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -2.0 * kd_130[k]
                   + f_0 * md_130[k];

        t_131[k] = -2.0 * kd_131[k]
                   + f_0 * md_131[k];

        t_132[k] = -2.0 * kd_132[k]
                   + f_0 * md_132[k];

        t_133[k] = -2.0 * kd_133[k]
                   + f_0 * md_133[k];

        t_134[k] = -2.0 * kd_134[k]
                   + f_0 * md_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, kd_135, kd_136, kd_137, kd_138, \
                         kd_139, md_135, md_136, md_137, md_138, \
                         md_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -2.0 * kd_135[k]
                   + f_0 * md_135[k];

        t_136[k] = -2.0 * kd_136[k]
                   + f_0 * md_136[k];

        t_137[k] = -2.0 * kd_137[k]
                   + f_0 * md_137[k];

        t_138[k] = -2.0 * kd_138[k]
                   + f_0 * md_138[k];

        t_139[k] = -2.0 * kd_139[k]
                   + f_0 * md_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, kd_140, kd_141, kd_142, kd_143, \
                         kd_144, md_140, md_141, md_142, md_143, \
                         md_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -2.0 * kd_140[k]
                   + f_0 * md_140[k];

        t_141[k] = -2.0 * kd_141[k]
                   + f_0 * md_141[k];

        t_142[k] = -2.0 * kd_142[k]
                   + f_0 * md_142[k];

        t_143[k] = -2.0 * kd_143[k]
                   + f_0 * md_143[k];

        t_144[k] = -2.0 * kd_144[k]
                   + f_0 * md_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, kd_145, kd_146, kd_147, kd_148, \
                         kd_149, md_145, md_146, md_147, md_148, \
                         md_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -2.0 * kd_145[k]
                   + f_0 * md_145[k];

        t_146[k] = -2.0 * kd_146[k]
                   + f_0 * md_146[k];

        t_147[k] = -2.0 * kd_147[k]
                   + f_0 * md_147[k];

        t_148[k] = -2.0 * kd_148[k]
                   + f_0 * md_148[k];

        t_149[k] = -2.0 * kd_149[k]
                   + f_0 * md_149[k];
    }
}

static auto
compute_prim_geom_10_ld_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kd, const size_t md,
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

    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_163 = buffer.data(kd + 163);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_166 = buffer.data(kd + 166);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_168 = buffer.data(kd + 168);
    const auto *kd_169 = buffer.data(kd + 169);
    const auto *kd_170 = buffer.data(kd + 170);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_172 = buffer.data(kd + 172);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_174 = buffer.data(kd + 174);
    const auto *kd_175 = buffer.data(kd + 175);
    const auto *kd_176 = buffer.data(kd + 176);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_181 = buffer.data(kd + 181);
    const auto *kd_182 = buffer.data(kd + 182);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_187 = buffer.data(kd + 187);
    const auto *kd_188 = buffer.data(kd + 188);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_193 = buffer.data(kd + 193);
    const auto *kd_194 = buffer.data(kd + 194);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_199 = buffer.data(kd + 199);
    const auto *kd_200 = buffer.data(kd + 200);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_204 = buffer.data(kd + 204);
    const auto *kd_205 = buffer.data(kd + 205);
    const auto *kd_206 = buffer.data(kd + 206);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_210 = buffer.data(kd + 210);
    const auto *kd_211 = buffer.data(kd + 211);
    const auto *kd_212 = buffer.data(kd + 212);
    const auto *kd_213 = buffer.data(kd + 213);
    const auto *kd_214 = buffer.data(kd + 214);
    const auto *kd_215 = buffer.data(kd + 215);

    const auto *md_150 = buffer.data(md + 150);
    const auto *md_151 = buffer.data(md + 151);
    const auto *md_152 = buffer.data(md + 152);
    const auto *md_153 = buffer.data(md + 153);
    const auto *md_154 = buffer.data(md + 154);
    const auto *md_155 = buffer.data(md + 155);
    const auto *md_156 = buffer.data(md + 156);
    const auto *md_157 = buffer.data(md + 157);
    const auto *md_158 = buffer.data(md + 158);
    const auto *md_159 = buffer.data(md + 159);
    const auto *md_160 = buffer.data(md + 160);
    const auto *md_161 = buffer.data(md + 161);
    const auto *md_162 = buffer.data(md + 162);
    const auto *md_163 = buffer.data(md + 163);
    const auto *md_164 = buffer.data(md + 164);
    const auto *md_165 = buffer.data(md + 165);
    const auto *md_166 = buffer.data(md + 166);
    const auto *md_167 = buffer.data(md + 167);
    const auto *md_168 = buffer.data(md + 168);
    const auto *md_169 = buffer.data(md + 169);
    const auto *md_170 = buffer.data(md + 170);
    const auto *md_171 = buffer.data(md + 171);
    const auto *md_172 = buffer.data(md + 172);
    const auto *md_173 = buffer.data(md + 173);
    const auto *md_174 = buffer.data(md + 174);
    const auto *md_175 = buffer.data(md + 175);
    const auto *md_176 = buffer.data(md + 176);
    const auto *md_177 = buffer.data(md + 177);
    const auto *md_178 = buffer.data(md + 178);
    const auto *md_179 = buffer.data(md + 179);
    const auto *md_180 = buffer.data(md + 180);
    const auto *md_181 = buffer.data(md + 181);
    const auto *md_182 = buffer.data(md + 182);
    const auto *md_183 = buffer.data(md + 183);
    const auto *md_184 = buffer.data(md + 184);
    const auto *md_185 = buffer.data(md + 185);
    const auto *md_186 = buffer.data(md + 186);
    const auto *md_187 = buffer.data(md + 187);
    const auto *md_188 = buffer.data(md + 188);
    const auto *md_189 = buffer.data(md + 189);
    const auto *md_190 = buffer.data(md + 190);
    const auto *md_191 = buffer.data(md + 191);
    const auto *md_192 = buffer.data(md + 192);
    const auto *md_193 = buffer.data(md + 193);
    const auto *md_194 = buffer.data(md + 194);
    const auto *md_195 = buffer.data(md + 195);
    const auto *md_196 = buffer.data(md + 196);
    const auto *md_197 = buffer.data(md + 197);
    const auto *md_198 = buffer.data(md + 198);
    const auto *md_199 = buffer.data(md + 199);
    const auto *md_200 = buffer.data(md + 200);
    const auto *md_201 = buffer.data(md + 201);
    const auto *md_202 = buffer.data(md + 202);
    const auto *md_203 = buffer.data(md + 203);
    const auto *md_204 = buffer.data(md + 204);
    const auto *md_205 = buffer.data(md + 205);
    const auto *md_206 = buffer.data(md + 206);
    const auto *md_207 = buffer.data(md + 207);
    const auto *md_208 = buffer.data(md + 208);
    const auto *md_209 = buffer.data(md + 209);
    const auto *md_210 = buffer.data(md + 210);
    const auto *md_211 = buffer.data(md + 211);
    const auto *md_212 = buffer.data(md + 212);
    const auto *md_213 = buffer.data(md + 213);
    const auto *md_214 = buffer.data(md + 214);
    const auto *md_215 = buffer.data(md + 215);
    const auto *md_216 = buffer.data(md + 216);
    const auto *md_217 = buffer.data(md + 217);
    const auto *md_218 = buffer.data(md + 218);
    const auto *md_219 = buffer.data(md + 219);
    const auto *md_220 = buffer.data(md + 220);
    const auto *md_221 = buffer.data(md + 221);
    const auto *md_222 = buffer.data(md + 222);
    const auto *md_223 = buffer.data(md + 223);
    const auto *md_224 = buffer.data(md + 224);
    const auto *md_225 = buffer.data(md + 225);
    const auto *md_226 = buffer.data(md + 226);
    const auto *md_227 = buffer.data(md + 227);
    const auto *md_228 = buffer.data(md + 228);
    const auto *md_229 = buffer.data(md + 229);
    const auto *md_230 = buffer.data(md + 230);
    const auto *md_231 = buffer.data(md + 231);
    const auto *md_232 = buffer.data(md + 232);
    const auto *md_233 = buffer.data(md + 233);
    const auto *md_234 = buffer.data(md + 234);
    const auto *md_235 = buffer.data(md + 235);
    const auto *md_236 = buffer.data(md + 236);
    const auto *md_237 = buffer.data(md + 237);
    const auto *md_238 = buffer.data(md + 238);
    const auto *md_239 = buffer.data(md + 239);
    const auto *md_240 = buffer.data(md + 240);
    const auto *md_241 = buffer.data(md + 241);
    const auto *md_242 = buffer.data(md + 242);
    const auto *md_243 = buffer.data(md + 243);
    const auto *md_244 = buffer.data(md + 244);
    const auto *md_245 = buffer.data(md + 245);
    const auto *md_246 = buffer.data(md + 246);
    const auto *md_247 = buffer.data(md + 247);
    const auto *md_248 = buffer.data(md + 248);
    const auto *md_249 = buffer.data(md + 249);
    const auto *md_250 = buffer.data(md + 250);
    const auto *md_251 = buffer.data(md + 251);
    const auto *md_252 = buffer.data(md + 252);
    const auto *md_253 = buffer.data(md + 253);
    const auto *md_254 = buffer.data(md + 254);
    const auto *md_255 = buffer.data(md + 255);
    const auto *md_256 = buffer.data(md + 256);
    const auto *md_257 = buffer.data(md + 257);
    const auto *md_258 = buffer.data(md + 258);
    const auto *md_259 = buffer.data(md + 259);
    const auto *md_260 = buffer.data(md + 260);
    const auto *md_261 = buffer.data(md + 261);
    const auto *md_262 = buffer.data(md + 262);
    const auto *md_263 = buffer.data(md + 263);
    const auto *md_264 = buffer.data(md + 264);
    const auto *md_265 = buffer.data(md + 265);
    const auto *md_266 = buffer.data(md + 266);
    const auto *md_267 = buffer.data(md + 267);
    const auto *md_268 = buffer.data(md + 268);
    const auto *md_269 = buffer.data(md + 269);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, kd_150, kd_151, kd_152, kd_153, \
                         kd_154, md_150, md_151, md_152, md_153, \
                         md_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -2.0 * kd_150[k]
                   + f_0 * md_150[k];

        t_151[k] = -2.0 * kd_151[k]
                   + f_0 * md_151[k];

        t_152[k] = -2.0 * kd_152[k]
                   + f_0 * md_152[k];

        t_153[k] = -2.0 * kd_153[k]
                   + f_0 * md_153[k];

        t_154[k] = -2.0 * kd_154[k]
                   + f_0 * md_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, kd_155, kd_156, kd_157, kd_158, \
                         kd_159, md_155, md_156, md_157, md_158, \
                         md_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -2.0 * kd_155[k]
                   + f_0 * md_155[k];

        t_156[k] = -2.0 * kd_156[k]
                   + f_0 * md_156[k];

        t_157[k] = -2.0 * kd_157[k]
                   + f_0 * md_157[k];

        t_158[k] = -2.0 * kd_158[k]
                   + f_0 * md_158[k];

        t_159[k] = -2.0 * kd_159[k]
                   + f_0 * md_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, kd_160, kd_161, kd_162, kd_163, \
                         kd_164, md_160, md_161, md_162, md_163, \
                         md_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -2.0 * kd_160[k]
                   + f_0 * md_160[k];

        t_161[k] = -2.0 * kd_161[k]
                   + f_0 * md_161[k];

        t_162[k] = -2.0 * kd_162[k]
                   + f_0 * md_162[k];

        t_163[k] = -2.0 * kd_163[k]
                   + f_0 * md_163[k];

        t_164[k] = -2.0 * kd_164[k]
                   + f_0 * md_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, kd_165, kd_166, kd_167, kd_168, \
                         kd_169, md_165, md_166, md_167, md_168, \
                         md_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -2.0 * kd_165[k]
                   + f_0 * md_165[k];

        t_166[k] = -2.0 * kd_166[k]
                   + f_0 * md_166[k];

        t_167[k] = -2.0 * kd_167[k]
                   + f_0 * md_167[k];

        t_168[k] = -kd_168[k]
                   + f_0 * md_168[k];

        t_169[k] = -kd_169[k]
                   + f_0 * md_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, kd_170, kd_171, kd_172, kd_173, \
                         kd_174, md_170, md_171, md_172, md_173, \
                         md_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -kd_170[k]
                   + f_0 * md_170[k];

        t_171[k] = -kd_171[k]
                   + f_0 * md_171[k];

        t_172[k] = -kd_172[k]
                   + f_0 * md_172[k];

        t_173[k] = -kd_173[k]
                   + f_0 * md_173[k];

        t_174[k] = -kd_174[k]
                   + f_0 * md_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, kd_175, kd_176, kd_177, kd_178, \
                         kd_179, md_175, md_176, md_177, md_178, \
                         md_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -kd_175[k]
                   + f_0 * md_175[k];

        t_176[k] = -kd_176[k]
                   + f_0 * md_176[k];

        t_177[k] = -kd_177[k]
                   + f_0 * md_177[k];

        t_178[k] = -kd_178[k]
                   + f_0 * md_178[k];

        t_179[k] = -kd_179[k]
                   + f_0 * md_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, kd_180, kd_181, kd_182, kd_183, \
                         kd_184, md_180, md_181, md_182, md_183, \
                         md_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -kd_180[k]
                   + f_0 * md_180[k];

        t_181[k] = -kd_181[k]
                   + f_0 * md_181[k];

        t_182[k] = -kd_182[k]
                   + f_0 * md_182[k];

        t_183[k] = -kd_183[k]
                   + f_0 * md_183[k];

        t_184[k] = -kd_184[k]
                   + f_0 * md_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, kd_185, kd_186, kd_187, kd_188, \
                         kd_189, md_185, md_186, md_187, md_188, \
                         md_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -kd_185[k]
                   + f_0 * md_185[k];

        t_186[k] = -kd_186[k]
                   + f_0 * md_186[k];

        t_187[k] = -kd_187[k]
                   + f_0 * md_187[k];

        t_188[k] = -kd_188[k]
                   + f_0 * md_188[k];

        t_189[k] = -kd_189[k]
                   + f_0 * md_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, kd_190, kd_191, kd_192, kd_193, \
                         kd_194, md_190, md_191, md_192, md_193, \
                         md_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -kd_190[k]
                   + f_0 * md_190[k];

        t_191[k] = -kd_191[k]
                   + f_0 * md_191[k];

        t_192[k] = -kd_192[k]
                   + f_0 * md_192[k];

        t_193[k] = -kd_193[k]
                   + f_0 * md_193[k];

        t_194[k] = -kd_194[k]
                   + f_0 * md_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, kd_195, kd_196, kd_197, kd_198, \
                         kd_199, md_195, md_196, md_197, md_198, \
                         md_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -kd_195[k]
                   + f_0 * md_195[k];

        t_196[k] = -kd_196[k]
                   + f_0 * md_196[k];

        t_197[k] = -kd_197[k]
                   + f_0 * md_197[k];

        t_198[k] = -kd_198[k]
                   + f_0 * md_198[k];

        t_199[k] = -kd_199[k]
                   + f_0 * md_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, kd_200, kd_201, kd_202, kd_203, \
                         kd_204, md_200, md_201, md_202, md_203, \
                         md_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -kd_200[k]
                   + f_0 * md_200[k];

        t_201[k] = -kd_201[k]
                   + f_0 * md_201[k];

        t_202[k] = -kd_202[k]
                   + f_0 * md_202[k];

        t_203[k] = -kd_203[k]
                   + f_0 * md_203[k];

        t_204[k] = -kd_204[k]
                   + f_0 * md_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, kd_205, kd_206, kd_207, kd_208, \
                         kd_209, md_205, md_206, md_207, md_208, \
                         md_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -kd_205[k]
                   + f_0 * md_205[k];

        t_206[k] = -kd_206[k]
                   + f_0 * md_206[k];

        t_207[k] = -kd_207[k]
                   + f_0 * md_207[k];

        t_208[k] = -kd_208[k]
                   + f_0 * md_208[k];

        t_209[k] = -kd_209[k]
                   + f_0 * md_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, kd_210, kd_211, kd_212, kd_213, \
                         kd_214, md_210, md_211, md_212, md_213, \
                         md_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -kd_210[k]
                   + f_0 * md_210[k];

        t_211[k] = -kd_211[k]
                   + f_0 * md_211[k];

        t_212[k] = -kd_212[k]
                   + f_0 * md_212[k];

        t_213[k] = -kd_213[k]
                   + f_0 * md_213[k];

        t_214[k] = -kd_214[k]
                   + f_0 * md_214[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, t_220, t_221, kd_215, md_215, \
                         md_216, md_217, md_218, md_219, md_220, \
                         md_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -kd_215[k]
                   + f_0 * md_215[k];

        t_216[k] = f_0 * md_216[k];

        t_217[k] = f_0 * md_217[k];

        t_218[k] = f_0 * md_218[k];

        t_219[k] = f_0 * md_219[k];

        t_220[k] = f_0 * md_220[k];

        t_221[k] = f_0 * md_221[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, t_228, t_229, md_222, \
                         md_223, md_224, md_225, md_226, md_227, md_228, \
                         md_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * md_222[k];

        t_223[k] = f_0 * md_223[k];

        t_224[k] = f_0 * md_224[k];

        t_225[k] = f_0 * md_225[k];

        t_226[k] = f_0 * md_226[k];

        t_227[k] = f_0 * md_227[k];

        t_228[k] = f_0 * md_228[k];

        t_229[k] = f_0 * md_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, t_235, t_236, t_237, md_230, \
                         md_231, md_232, md_233, md_234, md_235, md_236, \
                         md_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_0 * md_230[k];

        t_231[k] = f_0 * md_231[k];

        t_232[k] = f_0 * md_232[k];

        t_233[k] = f_0 * md_233[k];

        t_234[k] = f_0 * md_234[k];

        t_235[k] = f_0 * md_235[k];

        t_236[k] = f_0 * md_236[k];

        t_237[k] = f_0 * md_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, t_243, t_244, t_245, md_238, \
                         md_239, md_240, md_241, md_242, md_243, md_244, \
                         md_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * md_238[k];

        t_239[k] = f_0 * md_239[k];

        t_240[k] = f_0 * md_240[k];

        t_241[k] = f_0 * md_241[k];

        t_242[k] = f_0 * md_242[k];

        t_243[k] = f_0 * md_243[k];

        t_244[k] = f_0 * md_244[k];

        t_245[k] = f_0 * md_245[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, t_253, md_246, \
                         md_247, md_248, md_249, md_250, md_251, md_252, \
                         md_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_0 * md_246[k];

        t_247[k] = f_0 * md_247[k];

        t_248[k] = f_0 * md_248[k];

        t_249[k] = f_0 * md_249[k];

        t_250[k] = f_0 * md_250[k];

        t_251[k] = f_0 * md_251[k];

        t_252[k] = f_0 * md_252[k];

        t_253[k] = f_0 * md_253[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, t_259, t_260, t_261, md_254, \
                         md_255, md_256, md_257, md_258, md_259, md_260, \
                         md_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_0 * md_254[k];

        t_255[k] = f_0 * md_255[k];

        t_256[k] = f_0 * md_256[k];

        t_257[k] = f_0 * md_257[k];

        t_258[k] = f_0 * md_258[k];

        t_259[k] = f_0 * md_259[k];

        t_260[k] = f_0 * md_260[k];

        t_261[k] = f_0 * md_261[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, t_267, t_268, t_269, md_262, \
                         md_263, md_264, md_265, md_266, md_267, md_268, \
                         md_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_0 * md_262[k];

        t_263[k] = f_0 * md_263[k];

        t_264[k] = f_0 * md_264[k];

        t_265[k] = f_0 * md_265[k];

        t_266[k] = f_0 * md_266[k];

        t_267[k] = f_0 * md_267[k];

        t_268[k] = f_0 * md_268[k];

        t_269[k] = f_0 * md_269[k];
    }
}

auto
compute_prim_geom_10_ld_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t kd, const size_t md,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ld_electron_repulsion_0_piece0(buffer, target, kd, md, ncols, alpha);

    compute_prim_geom_10_ld_electron_repulsion_0_piece1(buffer, target, kd, md, ncols, alpha);
}

static auto
compute_prim_geom_10_ld_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kd, const size_t md,
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

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_92 = buffer.data(kd + 92);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);

    const auto *md_6 = buffer.data(md + 6);
    const auto *md_7 = buffer.data(md + 7);
    const auto *md_8 = buffer.data(md + 8);
    const auto *md_9 = buffer.data(md + 9);
    const auto *md_10 = buffer.data(md + 10);
    const auto *md_11 = buffer.data(md + 11);
    const auto *md_18 = buffer.data(md + 18);
    const auto *md_19 = buffer.data(md + 19);
    const auto *md_20 = buffer.data(md + 20);
    const auto *md_21 = buffer.data(md + 21);
    const auto *md_22 = buffer.data(md + 22);
    const auto *md_23 = buffer.data(md + 23);
    const auto *md_24 = buffer.data(md + 24);
    const auto *md_25 = buffer.data(md + 25);
    const auto *md_26 = buffer.data(md + 26);
    const auto *md_27 = buffer.data(md + 27);
    const auto *md_28 = buffer.data(md + 28);
    const auto *md_29 = buffer.data(md + 29);
    const auto *md_36 = buffer.data(md + 36);
    const auto *md_37 = buffer.data(md + 37);
    const auto *md_38 = buffer.data(md + 38);
    const auto *md_39 = buffer.data(md + 39);
    const auto *md_40 = buffer.data(md + 40);
    const auto *md_41 = buffer.data(md + 41);
    const auto *md_42 = buffer.data(md + 42);
    const auto *md_43 = buffer.data(md + 43);
    const auto *md_44 = buffer.data(md + 44);
    const auto *md_45 = buffer.data(md + 45);
    const auto *md_46 = buffer.data(md + 46);
    const auto *md_47 = buffer.data(md + 47);
    const auto *md_48 = buffer.data(md + 48);
    const auto *md_49 = buffer.data(md + 49);
    const auto *md_50 = buffer.data(md + 50);
    const auto *md_51 = buffer.data(md + 51);
    const auto *md_52 = buffer.data(md + 52);
    const auto *md_53 = buffer.data(md + 53);
    const auto *md_60 = buffer.data(md + 60);
    const auto *md_61 = buffer.data(md + 61);
    const auto *md_62 = buffer.data(md + 62);
    const auto *md_63 = buffer.data(md + 63);
    const auto *md_64 = buffer.data(md + 64);
    const auto *md_65 = buffer.data(md + 65);
    const auto *md_66 = buffer.data(md + 66);
    const auto *md_67 = buffer.data(md + 67);
    const auto *md_68 = buffer.data(md + 68);
    const auto *md_69 = buffer.data(md + 69);
    const auto *md_70 = buffer.data(md + 70);
    const auto *md_71 = buffer.data(md + 71);
    const auto *md_72 = buffer.data(md + 72);
    const auto *md_73 = buffer.data(md + 73);
    const auto *md_74 = buffer.data(md + 74);
    const auto *md_75 = buffer.data(md + 75);
    const auto *md_76 = buffer.data(md + 76);
    const auto *md_77 = buffer.data(md + 77);
    const auto *md_78 = buffer.data(md + 78);
    const auto *md_79 = buffer.data(md + 79);
    const auto *md_80 = buffer.data(md + 80);
    const auto *md_81 = buffer.data(md + 81);
    const auto *md_82 = buffer.data(md + 82);
    const auto *md_83 = buffer.data(md + 83);
    const auto *md_90 = buffer.data(md + 90);
    const auto *md_91 = buffer.data(md + 91);
    const auto *md_92 = buffer.data(md + 92);
    const auto *md_93 = buffer.data(md + 93);
    const auto *md_94 = buffer.data(md + 94);
    const auto *md_95 = buffer.data(md + 95);
    const auto *md_96 = buffer.data(md + 96);
    const auto *md_97 = buffer.data(md + 97);
    const auto *md_98 = buffer.data(md + 98);
    const auto *md_99 = buffer.data(md + 99);
    const auto *md_100 = buffer.data(md + 100);
    const auto *md_101 = buffer.data(md + 101);
    const auto *md_102 = buffer.data(md + 102);
    const auto *md_103 = buffer.data(md + 103);
    const auto *md_104 = buffer.data(md + 104);
    const auto *md_105 = buffer.data(md + 105);
    const auto *md_106 = buffer.data(md + 106);
    const auto *md_107 = buffer.data(md + 107);
    const auto *md_108 = buffer.data(md + 108);
    const auto *md_109 = buffer.data(md + 109);
    const auto *md_110 = buffer.data(md + 110);
    const auto *md_111 = buffer.data(md + 111);
    const auto *md_112 = buffer.data(md + 112);
    const auto *md_113 = buffer.data(md + 113);
    const auto *md_114 = buffer.data(md + 114);
    const auto *md_115 = buffer.data(md + 115);
    const auto *md_116 = buffer.data(md + 116);
    const auto *md_117 = buffer.data(md + 117);
    const auto *md_118 = buffer.data(md + 118);
    const auto *md_119 = buffer.data(md + 119);
    const auto *md_126 = buffer.data(md + 126);
    const auto *md_127 = buffer.data(md + 127);
    const auto *md_128 = buffer.data(md + 128);
    const auto *md_129 = buffer.data(md + 129);
    const auto *md_130 = buffer.data(md + 130);
    const auto *md_131 = buffer.data(md + 131);
    const auto *md_132 = buffer.data(md + 132);
    const auto *md_133 = buffer.data(md + 133);
    const auto *md_134 = buffer.data(md + 134);
    const auto *md_135 = buffer.data(md + 135);
    const auto *md_136 = buffer.data(md + 136);
    const auto *md_137 = buffer.data(md + 137);
    const auto *md_138 = buffer.data(md + 138);
    const auto *md_139 = buffer.data(md + 139);
    const auto *md_140 = buffer.data(md + 140);
    const auto *md_141 = buffer.data(md + 141);
    const auto *md_142 = buffer.data(md + 142);
    const auto *md_143 = buffer.data(md + 143);
    const auto *md_144 = buffer.data(md + 144);
    const auto *md_145 = buffer.data(md + 145);
    const auto *md_146 = buffer.data(md + 146);
    const auto *md_147 = buffer.data(md + 147);
    const auto *md_148 = buffer.data(md + 148);
    const auto *md_149 = buffer.data(md + 149);
    const auto *md_150 = buffer.data(md + 150);
    const auto *md_151 = buffer.data(md + 151);
    const auto *md_152 = buffer.data(md + 152);
    const auto *md_153 = buffer.data(md + 153);
    const auto *md_154 = buffer.data(md + 154);
    const auto *md_155 = buffer.data(md + 155);
    const auto *md_156 = buffer.data(md + 156);
    const auto *md_157 = buffer.data(md + 157);
    const auto *md_158 = buffer.data(md + 158);
    const auto *md_159 = buffer.data(md + 159);
    const auto *md_160 = buffer.data(md + 160);
    const auto *md_161 = buffer.data(md + 161);
    const auto *md_168 = buffer.data(md + 168);
    const auto *md_169 = buffer.data(md + 169);
    const auto *md_170 = buffer.data(md + 170);
    const auto *md_171 = buffer.data(md + 171);
    const auto *md_172 = buffer.data(md + 172);
    const auto *md_173 = buffer.data(md + 173);
    const auto *md_174 = buffer.data(md + 174);
    const auto *md_175 = buffer.data(md + 175);
    const auto *md_176 = buffer.data(md + 176);
    const auto *md_177 = buffer.data(md + 177);
    const auto *md_178 = buffer.data(md + 178);
    const auto *md_179 = buffer.data(md + 179);
    const auto *md_180 = buffer.data(md + 180);
    const auto *md_181 = buffer.data(md + 181);
    const auto *md_182 = buffer.data(md + 182);
    const auto *md_183 = buffer.data(md + 183);
    const auto *md_184 = buffer.data(md + 184);
    const auto *md_185 = buffer.data(md + 185);
    const auto *md_186 = buffer.data(md + 186);
    const auto *md_187 = buffer.data(md + 187);
    const auto *md_188 = buffer.data(md + 188);
    const auto *md_189 = buffer.data(md + 189);
    const auto *md_190 = buffer.data(md + 190);
    const auto *md_191 = buffer.data(md + 191);
    const auto *md_192 = buffer.data(md + 192);
    const auto *md_193 = buffer.data(md + 193);
    const auto *md_194 = buffer.data(md + 194);
    const auto *md_195 = buffer.data(md + 195);
    const auto *md_196 = buffer.data(md + 196);
    const auto *md_197 = buffer.data(md + 197);
    const auto *md_198 = buffer.data(md + 198);
    const auto *md_199 = buffer.data(md + 199);
    const auto *md_200 = buffer.data(md + 200);
    const auto *md_201 = buffer.data(md + 201);
    const auto *md_202 = buffer.data(md + 202);
    const auto *md_203 = buffer.data(md + 203);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, kd_0, md_6, md_7, md_8, md_9, \
                         md_10, md_11, md_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * md_6[k];

        t_1[k] = f_0 * md_7[k];

        t_2[k] = f_0 * md_8[k];

        t_3[k] = f_0 * md_9[k];

        t_4[k] = f_0 * md_10[k];

        t_5[k] = f_0 * md_11[k];

        t_6[k] = -kd_0[k]
                 + f_0 * md_18[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, kd_1, kd_2, kd_3, kd_4, kd_5, md_19, \
                         md_20, md_21, md_22, md_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -kd_1[k]
                 + f_0 * md_19[k];

        t_8[k] = -kd_2[k]
                 + f_0 * md_20[k];

        t_9[k] = -kd_3[k]
                 + f_0 * md_21[k];

        t_10[k] = -kd_4[k]
                  + f_0 * md_22[k];

        t_11[k] = -kd_5[k]
                  + f_0 * md_23[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, t_18, kd_6, md_24, md_25, md_26, \
                         md_27, md_28, md_29, md_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * md_24[k];

        t_13[k] = f_0 * md_25[k];

        t_14[k] = f_0 * md_26[k];

        t_15[k] = f_0 * md_27[k];

        t_16[k] = f_0 * md_28[k];

        t_17[k] = f_0 * md_29[k];

        t_18[k] = -2.0 * kd_6[k]
                  + f_0 * md_36[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, kd_7, kd_8, kd_9, kd_10, kd_11, md_37, \
                         md_38, md_39, md_40, md_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -2.0 * kd_7[k]
                  + f_0 * md_37[k];

        t_20[k] = -2.0 * kd_8[k]
                  + f_0 * md_38[k];

        t_21[k] = -2.0 * kd_9[k]
                  + f_0 * md_39[k];

        t_22[k] = -2.0 * kd_10[k]
                  + f_0 * md_40[k];

        t_23[k] = -2.0 * kd_11[k]
                  + f_0 * md_41[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, kd_12, kd_13, kd_14, kd_15, kd_16, \
                         md_42, md_43, md_44, md_45, md_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -kd_12[k]
                  + f_0 * md_42[k];

        t_25[k] = -kd_13[k]
                  + f_0 * md_43[k];

        t_26[k] = -kd_14[k]
                  + f_0 * md_44[k];

        t_27[k] = -kd_15[k]
                  + f_0 * md_45[k];

        t_28[k] = -kd_16[k]
                  + f_0 * md_46[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, t_35, kd_17, md_47, md_48, md_49, \
                         md_50, md_51, md_52, md_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -kd_17[k]
                  + f_0 * md_47[k];

        t_30[k] = f_0 * md_48[k];

        t_31[k] = f_0 * md_49[k];

        t_32[k] = f_0 * md_50[k];

        t_33[k] = f_0 * md_51[k];

        t_34[k] = f_0 * md_52[k];

        t_35[k] = f_0 * md_53[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, kd_18, kd_19, kd_20, kd_21, kd_22, \
                         md_60, md_61, md_62, md_63, md_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -3.0 * kd_18[k]
                  + f_0 * md_60[k];

        t_37[k] = -3.0 * kd_19[k]
                  + f_0 * md_61[k];

        t_38[k] = -3.0 * kd_20[k]
                  + f_0 * md_62[k];

        t_39[k] = -3.0 * kd_21[k]
                  + f_0 * md_63[k];

        t_40[k] = -3.0 * kd_22[k]
                  + f_0 * md_64[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, kd_23, kd_24, kd_25, kd_26, kd_27, \
                         md_65, md_66, md_67, md_68, md_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -3.0 * kd_23[k]
                  + f_0 * md_65[k];

        t_42[k] = -2.0 * kd_24[k]
                  + f_0 * md_66[k];

        t_43[k] = -2.0 * kd_25[k]
                  + f_0 * md_67[k];

        t_44[k] = -2.0 * kd_26[k]
                  + f_0 * md_68[k];

        t_45[k] = -2.0 * kd_27[k]
                  + f_0 * md_69[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, kd_28, kd_29, kd_30, kd_31, kd_32, \
                         md_70, md_71, md_72, md_73, md_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -2.0 * kd_28[k]
                  + f_0 * md_70[k];

        t_47[k] = -2.0 * kd_29[k]
                  + f_0 * md_71[k];

        t_48[k] = -kd_30[k]
                  + f_0 * md_72[k];

        t_49[k] = -kd_31[k]
                  + f_0 * md_73[k];

        t_50[k] = -kd_32[k]
                  + f_0 * md_74[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, kd_33, kd_34, kd_35, md_75, \
                         md_76, md_77, md_78, md_79, md_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -kd_33[k]
                  + f_0 * md_75[k];

        t_52[k] = -kd_34[k]
                  + f_0 * md_76[k];

        t_53[k] = -kd_35[k]
                  + f_0 * md_77[k];

        t_54[k] = f_0 * md_78[k];

        t_55[k] = f_0 * md_79[k];

        t_56[k] = f_0 * md_80[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, kd_36, kd_37, kd_38, md_81, \
                         md_82, md_83, md_90, md_91, md_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_0 * md_81[k];

        t_58[k] = f_0 * md_82[k];

        t_59[k] = f_0 * md_83[k];

        t_60[k] = -4.0 * kd_36[k]
                  + f_0 * md_90[k];

        t_61[k] = -4.0 * kd_37[k]
                  + f_0 * md_91[k];

        t_62[k] = -4.0 * kd_38[k]
                  + f_0 * md_92[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, kd_39, kd_40, kd_41, kd_42, kd_43, \
                         md_93, md_94, md_95, md_96, md_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -4.0 * kd_39[k]
                  + f_0 * md_93[k];

        t_64[k] = -4.0 * kd_40[k]
                  + f_0 * md_94[k];

        t_65[k] = -4.0 * kd_41[k]
                  + f_0 * md_95[k];

        t_66[k] = -3.0 * kd_42[k]
                  + f_0 * md_96[k];

        t_67[k] = -3.0 * kd_43[k]
                  + f_0 * md_97[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, kd_44, kd_45, kd_46, kd_47, kd_48, \
                         md_98, md_99, md_100, md_101, md_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -3.0 * kd_44[k]
                  + f_0 * md_98[k];

        t_69[k] = -3.0 * kd_45[k]
                  + f_0 * md_99[k];

        t_70[k] = -3.0 * kd_46[k]
                  + f_0 * md_100[k];

        t_71[k] = -3.0 * kd_47[k]
                  + f_0 * md_101[k];

        t_72[k] = -2.0 * kd_48[k]
                  + f_0 * md_102[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, kd_49, kd_50, kd_51, kd_52, kd_53, \
                         md_103, md_104, md_105, md_106, md_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = -2.0 * kd_49[k]
                  + f_0 * md_103[k];

        t_74[k] = -2.0 * kd_50[k]
                  + f_0 * md_104[k];

        t_75[k] = -2.0 * kd_51[k]
                  + f_0 * md_105[k];

        t_76[k] = -2.0 * kd_52[k]
                  + f_0 * md_106[k];

        t_77[k] = -2.0 * kd_53[k]
                  + f_0 * md_107[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, kd_54, kd_55, kd_56, kd_57, kd_58, \
                         md_108, md_109, md_110, md_111, md_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -kd_54[k]
                  + f_0 * md_108[k];

        t_79[k] = -kd_55[k]
                  + f_0 * md_109[k];

        t_80[k] = -kd_56[k]
                  + f_0 * md_110[k];

        t_81[k] = -kd_57[k]
                  + f_0 * md_111[k];

        t_82[k] = -kd_58[k]
                  + f_0 * md_112[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, kd_59, md_113, md_114, \
                         md_115, md_116, md_117, md_118, md_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = -kd_59[k]
                  + f_0 * md_113[k];

        t_84[k] = f_0 * md_114[k];

        t_85[k] = f_0 * md_115[k];

        t_86[k] = f_0 * md_116[k];

        t_87[k] = f_0 * md_117[k];

        t_88[k] = f_0 * md_118[k];

        t_89[k] = f_0 * md_119[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, kd_60, kd_61, kd_62, kd_63, kd_64, \
                         md_126, md_127, md_128, md_129, md_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * kd_60[k]
                  + f_0 * md_126[k];

        t_91[k] = -5.0 * kd_61[k]
                  + f_0 * md_127[k];

        t_92[k] = -5.0 * kd_62[k]
                  + f_0 * md_128[k];

        t_93[k] = -5.0 * kd_63[k]
                  + f_0 * md_129[k];

        t_94[k] = -5.0 * kd_64[k]
                  + f_0 * md_130[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, kd_65, kd_66, kd_67, kd_68, kd_69, \
                         md_131, md_132, md_133, md_134, md_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * kd_65[k]
                  + f_0 * md_131[k];

        t_96[k] = -4.0 * kd_66[k]
                  + f_0 * md_132[k];

        t_97[k] = -4.0 * kd_67[k]
                  + f_0 * md_133[k];

        t_98[k] = -4.0 * kd_68[k]
                  + f_0 * md_134[k];

        t_99[k] = -4.0 * kd_69[k]
                  + f_0 * md_135[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, kd_70, kd_71, kd_72, kd_73, kd_74, \
                         md_136, md_137, md_138, md_139, md_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * kd_70[k]
                   + f_0 * md_136[k];

        t_101[k] = -4.0 * kd_71[k]
                   + f_0 * md_137[k];

        t_102[k] = -3.0 * kd_72[k]
                   + f_0 * md_138[k];

        t_103[k] = -3.0 * kd_73[k]
                   + f_0 * md_139[k];

        t_104[k] = -3.0 * kd_74[k]
                   + f_0 * md_140[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, kd_75, kd_76, kd_77, kd_78, kd_79, \
                         md_141, md_142, md_143, md_144, md_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * kd_75[k]
                   + f_0 * md_141[k];

        t_106[k] = -3.0 * kd_76[k]
                   + f_0 * md_142[k];

        t_107[k] = -3.0 * kd_77[k]
                   + f_0 * md_143[k];

        t_108[k] = -2.0 * kd_78[k]
                   + f_0 * md_144[k];

        t_109[k] = -2.0 * kd_79[k]
                   + f_0 * md_145[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, kd_80, kd_81, kd_82, kd_83, kd_84, \
                         md_146, md_147, md_148, md_149, md_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * kd_80[k]
                   + f_0 * md_146[k];

        t_111[k] = -2.0 * kd_81[k]
                   + f_0 * md_147[k];

        t_112[k] = -2.0 * kd_82[k]
                   + f_0 * md_148[k];

        t_113[k] = -2.0 * kd_83[k]
                   + f_0 * md_149[k];

        t_114[k] = -kd_84[k]
                   + f_0 * md_150[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, kd_85, kd_86, kd_87, kd_88, kd_89, \
                         md_151, md_152, md_153, md_154, md_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -kd_85[k]
                   + f_0 * md_151[k];

        t_116[k] = -kd_86[k]
                   + f_0 * md_152[k];

        t_117[k] = -kd_87[k]
                   + f_0 * md_153[k];

        t_118[k] = -kd_88[k]
                   + f_0 * md_154[k];

        t_119[k] = -kd_89[k]
                   + f_0 * md_155[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, t_126, kd_90, md_156, \
                         md_157, md_158, md_159, md_160, md_161, \
                         md_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_0 * md_156[k];

        t_121[k] = f_0 * md_157[k];

        t_122[k] = f_0 * md_158[k];

        t_123[k] = f_0 * md_159[k];

        t_124[k] = f_0 * md_160[k];

        t_125[k] = f_0 * md_161[k];

        t_126[k] = -6.0 * kd_90[k]
                   + f_0 * md_168[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, kd_91, kd_92, kd_93, kd_94, kd_95, \
                         md_169, md_170, md_171, md_172, md_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -6.0 * kd_91[k]
                   + f_0 * md_169[k];

        t_128[k] = -6.0 * kd_92[k]
                   + f_0 * md_170[k];

        t_129[k] = -6.0 * kd_93[k]
                   + f_0 * md_171[k];

        t_130[k] = -6.0 * kd_94[k]
                   + f_0 * md_172[k];

        t_131[k] = -6.0 * kd_95[k]
                   + f_0 * md_173[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, kd_96, kd_97, kd_98, kd_99, \
                         kd_100, md_174, md_175, md_176, md_177, \
                         md_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -5.0 * kd_96[k]
                   + f_0 * md_174[k];

        t_133[k] = -5.0 * kd_97[k]
                   + f_0 * md_175[k];

        t_134[k] = -5.0 * kd_98[k]
                   + f_0 * md_176[k];

        t_135[k] = -5.0 * kd_99[k]
                   + f_0 * md_177[k];

        t_136[k] = -5.0 * kd_100[k]
                   + f_0 * md_178[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, kd_101, kd_102, kd_103, kd_104, \
                         kd_105, md_179, md_180, md_181, md_182, \
                         md_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -5.0 * kd_101[k]
                   + f_0 * md_179[k];

        t_138[k] = -4.0 * kd_102[k]
                   + f_0 * md_180[k];

        t_139[k] = -4.0 * kd_103[k]
                   + f_0 * md_181[k];

        t_140[k] = -4.0 * kd_104[k]
                   + f_0 * md_182[k];

        t_141[k] = -4.0 * kd_105[k]
                   + f_0 * md_183[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, kd_106, kd_107, kd_108, kd_109, \
                         kd_110, md_184, md_185, md_186, md_187, \
                         md_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -4.0 * kd_106[k]
                   + f_0 * md_184[k];

        t_143[k] = -4.0 * kd_107[k]
                   + f_0 * md_185[k];

        t_144[k] = -3.0 * kd_108[k]
                   + f_0 * md_186[k];

        t_145[k] = -3.0 * kd_109[k]
                   + f_0 * md_187[k];

        t_146[k] = -3.0 * kd_110[k]
                   + f_0 * md_188[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, kd_111, kd_112, kd_113, kd_114, \
                         kd_115, md_189, md_190, md_191, md_192, \
                         md_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -3.0 * kd_111[k]
                   + f_0 * md_189[k];

        t_148[k] = -3.0 * kd_112[k]
                   + f_0 * md_190[k];

        t_149[k] = -3.0 * kd_113[k]
                   + f_0 * md_191[k];

        t_150[k] = -2.0 * kd_114[k]
                   + f_0 * md_192[k];

        t_151[k] = -2.0 * kd_115[k]
                   + f_0 * md_193[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, kd_116, kd_117, kd_118, kd_119, \
                         kd_120, md_194, md_195, md_196, md_197, \
                         md_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -2.0 * kd_116[k]
                   + f_0 * md_194[k];

        t_153[k] = -2.0 * kd_117[k]
                   + f_0 * md_195[k];

        t_154[k] = -2.0 * kd_118[k]
                   + f_0 * md_196[k];

        t_155[k] = -2.0 * kd_119[k]
                   + f_0 * md_197[k];

        t_156[k] = -kd_120[k]
                   + f_0 * md_198[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, kd_121, kd_122, kd_123, kd_124, \
                         kd_125, md_199, md_200, md_201, md_202, \
                         md_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -kd_121[k]
                   + f_0 * md_199[k];

        t_158[k] = -kd_122[k]
                   + f_0 * md_200[k];

        t_159[k] = -kd_123[k]
                   + f_0 * md_201[k];

        t_160[k] = -kd_124[k]
                   + f_0 * md_202[k];

        t_161[k] = -kd_125[k]
                   + f_0 * md_203[k];
    }
}

static auto
compute_prim_geom_10_ld_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kd, const size_t md,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_163 = buffer.data(kd + 163);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_166 = buffer.data(kd + 166);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_168 = buffer.data(kd + 168);
    const auto *kd_169 = buffer.data(kd + 169);
    const auto *kd_170 = buffer.data(kd + 170);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_172 = buffer.data(kd + 172);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_174 = buffer.data(kd + 174);
    const auto *kd_175 = buffer.data(kd + 175);
    const auto *kd_176 = buffer.data(kd + 176);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_181 = buffer.data(kd + 181);
    const auto *kd_182 = buffer.data(kd + 182);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_187 = buffer.data(kd + 187);
    const auto *kd_188 = buffer.data(kd + 188);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_193 = buffer.data(kd + 193);
    const auto *kd_194 = buffer.data(kd + 194);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_199 = buffer.data(kd + 199);
    const auto *kd_200 = buffer.data(kd + 200);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_204 = buffer.data(kd + 204);
    const auto *kd_205 = buffer.data(kd + 205);
    const auto *kd_206 = buffer.data(kd + 206);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_210 = buffer.data(kd + 210);
    const auto *kd_211 = buffer.data(kd + 211);
    const auto *kd_212 = buffer.data(kd + 212);
    const auto *kd_213 = buffer.data(kd + 213);
    const auto *kd_214 = buffer.data(kd + 214);
    const auto *kd_215 = buffer.data(kd + 215);

    const auto *md_204 = buffer.data(md + 204);
    const auto *md_205 = buffer.data(md + 205);
    const auto *md_206 = buffer.data(md + 206);
    const auto *md_207 = buffer.data(md + 207);
    const auto *md_208 = buffer.data(md + 208);
    const auto *md_209 = buffer.data(md + 209);
    const auto *md_216 = buffer.data(md + 216);
    const auto *md_217 = buffer.data(md + 217);
    const auto *md_218 = buffer.data(md + 218);
    const auto *md_219 = buffer.data(md + 219);
    const auto *md_220 = buffer.data(md + 220);
    const auto *md_221 = buffer.data(md + 221);
    const auto *md_222 = buffer.data(md + 222);
    const auto *md_223 = buffer.data(md + 223);
    const auto *md_224 = buffer.data(md + 224);
    const auto *md_225 = buffer.data(md + 225);
    const auto *md_226 = buffer.data(md + 226);
    const auto *md_227 = buffer.data(md + 227);
    const auto *md_228 = buffer.data(md + 228);
    const auto *md_229 = buffer.data(md + 229);
    const auto *md_230 = buffer.data(md + 230);
    const auto *md_231 = buffer.data(md + 231);
    const auto *md_232 = buffer.data(md + 232);
    const auto *md_233 = buffer.data(md + 233);
    const auto *md_234 = buffer.data(md + 234);
    const auto *md_235 = buffer.data(md + 235);
    const auto *md_236 = buffer.data(md + 236);
    const auto *md_237 = buffer.data(md + 237);
    const auto *md_238 = buffer.data(md + 238);
    const auto *md_239 = buffer.data(md + 239);
    const auto *md_240 = buffer.data(md + 240);
    const auto *md_241 = buffer.data(md + 241);
    const auto *md_242 = buffer.data(md + 242);
    const auto *md_243 = buffer.data(md + 243);
    const auto *md_244 = buffer.data(md + 244);
    const auto *md_245 = buffer.data(md + 245);
    const auto *md_246 = buffer.data(md + 246);
    const auto *md_247 = buffer.data(md + 247);
    const auto *md_248 = buffer.data(md + 248);
    const auto *md_249 = buffer.data(md + 249);
    const auto *md_250 = buffer.data(md + 250);
    const auto *md_251 = buffer.data(md + 251);
    const auto *md_252 = buffer.data(md + 252);
    const auto *md_253 = buffer.data(md + 253);
    const auto *md_254 = buffer.data(md + 254);
    const auto *md_255 = buffer.data(md + 255);
    const auto *md_256 = buffer.data(md + 256);
    const auto *md_257 = buffer.data(md + 257);
    const auto *md_258 = buffer.data(md + 258);
    const auto *md_259 = buffer.data(md + 259);
    const auto *md_260 = buffer.data(md + 260);
    const auto *md_261 = buffer.data(md + 261);
    const auto *md_262 = buffer.data(md + 262);
    const auto *md_263 = buffer.data(md + 263);
    const auto *md_270 = buffer.data(md + 270);
    const auto *md_271 = buffer.data(md + 271);
    const auto *md_272 = buffer.data(md + 272);
    const auto *md_273 = buffer.data(md + 273);
    const auto *md_274 = buffer.data(md + 274);
    const auto *md_275 = buffer.data(md + 275);
    const auto *md_276 = buffer.data(md + 276);
    const auto *md_277 = buffer.data(md + 277);
    const auto *md_278 = buffer.data(md + 278);
    const auto *md_279 = buffer.data(md + 279);
    const auto *md_280 = buffer.data(md + 280);
    const auto *md_281 = buffer.data(md + 281);
    const auto *md_282 = buffer.data(md + 282);
    const auto *md_283 = buffer.data(md + 283);
    const auto *md_284 = buffer.data(md + 284);
    const auto *md_285 = buffer.data(md + 285);
    const auto *md_286 = buffer.data(md + 286);
    const auto *md_287 = buffer.data(md + 287);
    const auto *md_288 = buffer.data(md + 288);
    const auto *md_289 = buffer.data(md + 289);
    const auto *md_290 = buffer.data(md + 290);
    const auto *md_291 = buffer.data(md + 291);
    const auto *md_292 = buffer.data(md + 292);
    const auto *md_293 = buffer.data(md + 293);
    const auto *md_294 = buffer.data(md + 294);
    const auto *md_295 = buffer.data(md + 295);
    const auto *md_296 = buffer.data(md + 296);
    const auto *md_297 = buffer.data(md + 297);
    const auto *md_298 = buffer.data(md + 298);
    const auto *md_299 = buffer.data(md + 299);
    const auto *md_300 = buffer.data(md + 300);
    const auto *md_301 = buffer.data(md + 301);
    const auto *md_302 = buffer.data(md + 302);
    const auto *md_303 = buffer.data(md + 303);
    const auto *md_304 = buffer.data(md + 304);
    const auto *md_305 = buffer.data(md + 305);
    const auto *md_306 = buffer.data(md + 306);
    const auto *md_307 = buffer.data(md + 307);
    const auto *md_308 = buffer.data(md + 308);
    const auto *md_309 = buffer.data(md + 309);
    const auto *md_310 = buffer.data(md + 310);
    const auto *md_311 = buffer.data(md + 311);
    const auto *md_312 = buffer.data(md + 312);
    const auto *md_313 = buffer.data(md + 313);
    const auto *md_314 = buffer.data(md + 314);
    const auto *md_315 = buffer.data(md + 315);
    const auto *md_316 = buffer.data(md + 316);
    const auto *md_317 = buffer.data(md + 317);
    const auto *md_318 = buffer.data(md + 318);
    const auto *md_319 = buffer.data(md + 319);
    const auto *md_320 = buffer.data(md + 320);
    const auto *md_321 = buffer.data(md + 321);
    const auto *md_322 = buffer.data(md + 322);
    const auto *md_323 = buffer.data(md + 323);

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, t_168, kd_126, md_204, \
                         md_205, md_206, md_207, md_208, md_209, \
                         md_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_0 * md_204[k];

        t_163[k] = f_0 * md_205[k];

        t_164[k] = f_0 * md_206[k];

        t_165[k] = f_0 * md_207[k];

        t_166[k] = f_0 * md_208[k];

        t_167[k] = f_0 * md_209[k];

        t_168[k] = -7.0 * kd_126[k]
                   + f_0 * md_216[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, kd_127, kd_128, kd_129, kd_130, \
                         kd_131, md_217, md_218, md_219, md_220, \
                         md_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -7.0 * kd_127[k]
                   + f_0 * md_217[k];

        t_170[k] = -7.0 * kd_128[k]
                   + f_0 * md_218[k];

        t_171[k] = -7.0 * kd_129[k]
                   + f_0 * md_219[k];

        t_172[k] = -7.0 * kd_130[k]
                   + f_0 * md_220[k];

        t_173[k] = -7.0 * kd_131[k]
                   + f_0 * md_221[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, kd_132, kd_133, kd_134, kd_135, \
                         kd_136, md_222, md_223, md_224, md_225, \
                         md_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = -6.0 * kd_132[k]
                   + f_0 * md_222[k];

        t_175[k] = -6.0 * kd_133[k]
                   + f_0 * md_223[k];

        t_176[k] = -6.0 * kd_134[k]
                   + f_0 * md_224[k];

        t_177[k] = -6.0 * kd_135[k]
                   + f_0 * md_225[k];

        t_178[k] = -6.0 * kd_136[k]
                   + f_0 * md_226[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, kd_137, kd_138, kd_139, kd_140, \
                         kd_141, md_227, md_228, md_229, md_230, \
                         md_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = -6.0 * kd_137[k]
                   + f_0 * md_227[k];

        t_180[k] = -5.0 * kd_138[k]
                   + f_0 * md_228[k];

        t_181[k] = -5.0 * kd_139[k]
                   + f_0 * md_229[k];

        t_182[k] = -5.0 * kd_140[k]
                   + f_0 * md_230[k];

        t_183[k] = -5.0 * kd_141[k]
                   + f_0 * md_231[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, kd_142, kd_143, kd_144, kd_145, \
                         kd_146, md_232, md_233, md_234, md_235, \
                         md_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = -5.0 * kd_142[k]
                   + f_0 * md_232[k];

        t_185[k] = -5.0 * kd_143[k]
                   + f_0 * md_233[k];

        t_186[k] = -4.0 * kd_144[k]
                   + f_0 * md_234[k];

        t_187[k] = -4.0 * kd_145[k]
                   + f_0 * md_235[k];

        t_188[k] = -4.0 * kd_146[k]
                   + f_0 * md_236[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, kd_147, kd_148, kd_149, kd_150, \
                         kd_151, md_237, md_238, md_239, md_240, \
                         md_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = -4.0 * kd_147[k]
                   + f_0 * md_237[k];

        t_190[k] = -4.0 * kd_148[k]
                   + f_0 * md_238[k];

        t_191[k] = -4.0 * kd_149[k]
                   + f_0 * md_239[k];

        t_192[k] = -3.0 * kd_150[k]
                   + f_0 * md_240[k];

        t_193[k] = -3.0 * kd_151[k]
                   + f_0 * md_241[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, kd_152, kd_153, kd_154, kd_155, \
                         kd_156, md_242, md_243, md_244, md_245, \
                         md_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -3.0 * kd_152[k]
                   + f_0 * md_242[k];

        t_195[k] = -3.0 * kd_153[k]
                   + f_0 * md_243[k];

        t_196[k] = -3.0 * kd_154[k]
                   + f_0 * md_244[k];

        t_197[k] = -3.0 * kd_155[k]
                   + f_0 * md_245[k];

        t_198[k] = -2.0 * kd_156[k]
                   + f_0 * md_246[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, kd_157, kd_158, kd_159, kd_160, \
                         kd_161, md_247, md_248, md_249, md_250, \
                         md_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = -2.0 * kd_157[k]
                   + f_0 * md_247[k];

        t_200[k] = -2.0 * kd_158[k]
                   + f_0 * md_248[k];

        t_201[k] = -2.0 * kd_159[k]
                   + f_0 * md_249[k];

        t_202[k] = -2.0 * kd_160[k]
                   + f_0 * md_250[k];

        t_203[k] = -2.0 * kd_161[k]
                   + f_0 * md_251[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, kd_162, kd_163, kd_164, kd_165, \
                         kd_166, md_252, md_253, md_254, md_255, \
                         md_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -kd_162[k]
                   + f_0 * md_252[k];

        t_205[k] = -kd_163[k]
                   + f_0 * md_253[k];

        t_206[k] = -kd_164[k]
                   + f_0 * md_254[k];

        t_207[k] = -kd_165[k]
                   + f_0 * md_255[k];

        t_208[k] = -kd_166[k]
                   + f_0 * md_256[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, t_214, t_215, kd_167, md_257, \
                         md_258, md_259, md_260, md_261, md_262, \
                         md_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = -kd_167[k]
                   + f_0 * md_257[k];

        t_210[k] = f_0 * md_258[k];

        t_211[k] = f_0 * md_259[k];

        t_212[k] = f_0 * md_260[k];

        t_213[k] = f_0 * md_261[k];

        t_214[k] = f_0 * md_262[k];

        t_215[k] = f_0 * md_263[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, kd_168, kd_169, kd_170, kd_171, \
                         kd_172, md_270, md_271, md_272, md_273, \
                         md_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = -8.0 * kd_168[k]
                   + f_0 * md_270[k];

        t_217[k] = -8.0 * kd_169[k]
                   + f_0 * md_271[k];

        t_218[k] = -8.0 * kd_170[k]
                   + f_0 * md_272[k];

        t_219[k] = -8.0 * kd_171[k]
                   + f_0 * md_273[k];

        t_220[k] = -8.0 * kd_172[k]
                   + f_0 * md_274[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, kd_173, kd_174, kd_175, kd_176, \
                         kd_177, md_275, md_276, md_277, md_278, \
                         md_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = -8.0 * kd_173[k]
                   + f_0 * md_275[k];

        t_222[k] = -7.0 * kd_174[k]
                   + f_0 * md_276[k];

        t_223[k] = -7.0 * kd_175[k]
                   + f_0 * md_277[k];

        t_224[k] = -7.0 * kd_176[k]
                   + f_0 * md_278[k];

        t_225[k] = -7.0 * kd_177[k]
                   + f_0 * md_279[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, t_230, kd_178, kd_179, kd_180, kd_181, \
                         kd_182, md_280, md_281, md_282, md_283, \
                         md_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = -7.0 * kd_178[k]
                   + f_0 * md_280[k];

        t_227[k] = -7.0 * kd_179[k]
                   + f_0 * md_281[k];

        t_228[k] = -6.0 * kd_180[k]
                   + f_0 * md_282[k];

        t_229[k] = -6.0 * kd_181[k]
                   + f_0 * md_283[k];

        t_230[k] = -6.0 * kd_182[k]
                   + f_0 * md_284[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, t_235, kd_183, kd_184, kd_185, kd_186, \
                         kd_187, md_285, md_286, md_287, md_288, \
                         md_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = -6.0 * kd_183[k]
                   + f_0 * md_285[k];

        t_232[k] = -6.0 * kd_184[k]
                   + f_0 * md_286[k];

        t_233[k] = -6.0 * kd_185[k]
                   + f_0 * md_287[k];

        t_234[k] = -5.0 * kd_186[k]
                   + f_0 * md_288[k];

        t_235[k] = -5.0 * kd_187[k]
                   + f_0 * md_289[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, kd_188, kd_189, kd_190, kd_191, \
                         kd_192, md_290, md_291, md_292, md_293, \
                         md_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = -5.0 * kd_188[k]
                   + f_0 * md_290[k];

        t_237[k] = -5.0 * kd_189[k]
                   + f_0 * md_291[k];

        t_238[k] = -5.0 * kd_190[k]
                   + f_0 * md_292[k];

        t_239[k] = -5.0 * kd_191[k]
                   + f_0 * md_293[k];

        t_240[k] = -4.0 * kd_192[k]
                   + f_0 * md_294[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, kd_193, kd_194, kd_195, kd_196, \
                         kd_197, md_295, md_296, md_297, md_298, \
                         md_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = -4.0 * kd_193[k]
                   + f_0 * md_295[k];

        t_242[k] = -4.0 * kd_194[k]
                   + f_0 * md_296[k];

        t_243[k] = -4.0 * kd_195[k]
                   + f_0 * md_297[k];

        t_244[k] = -4.0 * kd_196[k]
                   + f_0 * md_298[k];

        t_245[k] = -4.0 * kd_197[k]
                   + f_0 * md_299[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, kd_198, kd_199, kd_200, kd_201, \
                         kd_202, md_300, md_301, md_302, md_303, \
                         md_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = -3.0 * kd_198[k]
                   + f_0 * md_300[k];

        t_247[k] = -3.0 * kd_199[k]
                   + f_0 * md_301[k];

        t_248[k] = -3.0 * kd_200[k]
                   + f_0 * md_302[k];

        t_249[k] = -3.0 * kd_201[k]
                   + f_0 * md_303[k];

        t_250[k] = -3.0 * kd_202[k]
                   + f_0 * md_304[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, kd_203, kd_204, kd_205, kd_206, \
                         kd_207, md_305, md_306, md_307, md_308, \
                         md_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = -3.0 * kd_203[k]
                   + f_0 * md_305[k];

        t_252[k] = -2.0 * kd_204[k]
                   + f_0 * md_306[k];

        t_253[k] = -2.0 * kd_205[k]
                   + f_0 * md_307[k];

        t_254[k] = -2.0 * kd_206[k]
                   + f_0 * md_308[k];

        t_255[k] = -2.0 * kd_207[k]
                   + f_0 * md_309[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, kd_208, kd_209, kd_210, kd_211, \
                         kd_212, md_310, md_311, md_312, md_313, \
                         md_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = -2.0 * kd_208[k]
                   + f_0 * md_310[k];

        t_257[k] = -2.0 * kd_209[k]
                   + f_0 * md_311[k];

        t_258[k] = -kd_210[k]
                   + f_0 * md_312[k];

        t_259[k] = -kd_211[k]
                   + f_0 * md_313[k];

        t_260[k] = -kd_212[k]
                   + f_0 * md_314[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, t_266, kd_213, kd_214, kd_215, \
                         md_315, md_316, md_317, md_318, md_319, \
                         md_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = -kd_213[k]
                   + f_0 * md_315[k];

        t_262[k] = -kd_214[k]
                   + f_0 * md_316[k];

        t_263[k] = -kd_215[k]
                   + f_0 * md_317[k];

        t_264[k] = f_0 * md_318[k];

        t_265[k] = f_0 * md_319[k];

        t_266[k] = f_0 * md_320[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, md_321, md_322, md_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_0 * md_321[k];

        t_268[k] = f_0 * md_322[k];

        t_269[k] = f_0 * md_323[k];
    }
}

auto
compute_prim_geom_10_ld_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t kd, const size_t md,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ld_electron_repulsion_1_piece0(buffer, target, kd, md, ncols, alpha);

    compute_prim_geom_10_ld_electron_repulsion_1_piece1(buffer, target, kd, md, ncols, alpha);
}

static auto
compute_prim_geom_10_ld_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kd, const size_t md,
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

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_92 = buffer.data(kd + 92);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);

    const auto *md_12 = buffer.data(md + 12);
    const auto *md_13 = buffer.data(md + 13);
    const auto *md_14 = buffer.data(md + 14);
    const auto *md_15 = buffer.data(md + 15);
    const auto *md_16 = buffer.data(md + 16);
    const auto *md_17 = buffer.data(md + 17);
    const auto *md_24 = buffer.data(md + 24);
    const auto *md_25 = buffer.data(md + 25);
    const auto *md_26 = buffer.data(md + 26);
    const auto *md_27 = buffer.data(md + 27);
    const auto *md_28 = buffer.data(md + 28);
    const auto *md_29 = buffer.data(md + 29);
    const auto *md_30 = buffer.data(md + 30);
    const auto *md_31 = buffer.data(md + 31);
    const auto *md_32 = buffer.data(md + 32);
    const auto *md_33 = buffer.data(md + 33);
    const auto *md_34 = buffer.data(md + 34);
    const auto *md_35 = buffer.data(md + 35);
    const auto *md_42 = buffer.data(md + 42);
    const auto *md_43 = buffer.data(md + 43);
    const auto *md_44 = buffer.data(md + 44);
    const auto *md_45 = buffer.data(md + 45);
    const auto *md_46 = buffer.data(md + 46);
    const auto *md_47 = buffer.data(md + 47);
    const auto *md_48 = buffer.data(md + 48);
    const auto *md_49 = buffer.data(md + 49);
    const auto *md_50 = buffer.data(md + 50);
    const auto *md_51 = buffer.data(md + 51);
    const auto *md_52 = buffer.data(md + 52);
    const auto *md_53 = buffer.data(md + 53);
    const auto *md_54 = buffer.data(md + 54);
    const auto *md_55 = buffer.data(md + 55);
    const auto *md_56 = buffer.data(md + 56);
    const auto *md_57 = buffer.data(md + 57);
    const auto *md_58 = buffer.data(md + 58);
    const auto *md_59 = buffer.data(md + 59);
    const auto *md_66 = buffer.data(md + 66);
    const auto *md_67 = buffer.data(md + 67);
    const auto *md_68 = buffer.data(md + 68);
    const auto *md_69 = buffer.data(md + 69);
    const auto *md_70 = buffer.data(md + 70);
    const auto *md_71 = buffer.data(md + 71);
    const auto *md_72 = buffer.data(md + 72);
    const auto *md_73 = buffer.data(md + 73);
    const auto *md_74 = buffer.data(md + 74);
    const auto *md_75 = buffer.data(md + 75);
    const auto *md_76 = buffer.data(md + 76);
    const auto *md_77 = buffer.data(md + 77);
    const auto *md_78 = buffer.data(md + 78);
    const auto *md_79 = buffer.data(md + 79);
    const auto *md_80 = buffer.data(md + 80);
    const auto *md_81 = buffer.data(md + 81);
    const auto *md_82 = buffer.data(md + 82);
    const auto *md_83 = buffer.data(md + 83);
    const auto *md_84 = buffer.data(md + 84);
    const auto *md_85 = buffer.data(md + 85);
    const auto *md_86 = buffer.data(md + 86);
    const auto *md_87 = buffer.data(md + 87);
    const auto *md_88 = buffer.data(md + 88);
    const auto *md_89 = buffer.data(md + 89);
    const auto *md_96 = buffer.data(md + 96);
    const auto *md_97 = buffer.data(md + 97);
    const auto *md_98 = buffer.data(md + 98);
    const auto *md_99 = buffer.data(md + 99);
    const auto *md_100 = buffer.data(md + 100);
    const auto *md_101 = buffer.data(md + 101);
    const auto *md_102 = buffer.data(md + 102);
    const auto *md_103 = buffer.data(md + 103);
    const auto *md_104 = buffer.data(md + 104);
    const auto *md_105 = buffer.data(md + 105);
    const auto *md_106 = buffer.data(md + 106);
    const auto *md_107 = buffer.data(md + 107);
    const auto *md_108 = buffer.data(md + 108);
    const auto *md_109 = buffer.data(md + 109);
    const auto *md_110 = buffer.data(md + 110);
    const auto *md_111 = buffer.data(md + 111);
    const auto *md_112 = buffer.data(md + 112);
    const auto *md_113 = buffer.data(md + 113);
    const auto *md_114 = buffer.data(md + 114);
    const auto *md_115 = buffer.data(md + 115);
    const auto *md_116 = buffer.data(md + 116);
    const auto *md_117 = buffer.data(md + 117);
    const auto *md_118 = buffer.data(md + 118);
    const auto *md_119 = buffer.data(md + 119);
    const auto *md_120 = buffer.data(md + 120);
    const auto *md_121 = buffer.data(md + 121);
    const auto *md_122 = buffer.data(md + 122);
    const auto *md_123 = buffer.data(md + 123);
    const auto *md_124 = buffer.data(md + 124);
    const auto *md_125 = buffer.data(md + 125);
    const auto *md_132 = buffer.data(md + 132);
    const auto *md_133 = buffer.data(md + 133);
    const auto *md_134 = buffer.data(md + 134);
    const auto *md_135 = buffer.data(md + 135);
    const auto *md_136 = buffer.data(md + 136);
    const auto *md_137 = buffer.data(md + 137);
    const auto *md_138 = buffer.data(md + 138);
    const auto *md_139 = buffer.data(md + 139);
    const auto *md_140 = buffer.data(md + 140);
    const auto *md_141 = buffer.data(md + 141);
    const auto *md_142 = buffer.data(md + 142);
    const auto *md_143 = buffer.data(md + 143);
    const auto *md_144 = buffer.data(md + 144);
    const auto *md_145 = buffer.data(md + 145);
    const auto *md_146 = buffer.data(md + 146);
    const auto *md_147 = buffer.data(md + 147);
    const auto *md_148 = buffer.data(md + 148);
    const auto *md_149 = buffer.data(md + 149);
    const auto *md_150 = buffer.data(md + 150);
    const auto *md_151 = buffer.data(md + 151);
    const auto *md_152 = buffer.data(md + 152);
    const auto *md_153 = buffer.data(md + 153);
    const auto *md_154 = buffer.data(md + 154);
    const auto *md_155 = buffer.data(md + 155);
    const auto *md_156 = buffer.data(md + 156);
    const auto *md_157 = buffer.data(md + 157);
    const auto *md_158 = buffer.data(md + 158);
    const auto *md_159 = buffer.data(md + 159);
    const auto *md_160 = buffer.data(md + 160);
    const auto *md_161 = buffer.data(md + 161);
    const auto *md_162 = buffer.data(md + 162);
    const auto *md_163 = buffer.data(md + 163);
    const auto *md_164 = buffer.data(md + 164);
    const auto *md_165 = buffer.data(md + 165);
    const auto *md_166 = buffer.data(md + 166);
    const auto *md_167 = buffer.data(md + 167);
    const auto *md_174 = buffer.data(md + 174);
    const auto *md_175 = buffer.data(md + 175);
    const auto *md_176 = buffer.data(md + 176);
    const auto *md_177 = buffer.data(md + 177);
    const auto *md_178 = buffer.data(md + 178);
    const auto *md_179 = buffer.data(md + 179);
    const auto *md_180 = buffer.data(md + 180);
    const auto *md_181 = buffer.data(md + 181);
    const auto *md_182 = buffer.data(md + 182);
    const auto *md_183 = buffer.data(md + 183);
    const auto *md_184 = buffer.data(md + 184);
    const auto *md_185 = buffer.data(md + 185);
    const auto *md_186 = buffer.data(md + 186);
    const auto *md_187 = buffer.data(md + 187);
    const auto *md_188 = buffer.data(md + 188);
    const auto *md_189 = buffer.data(md + 189);
    const auto *md_190 = buffer.data(md + 190);
    const auto *md_191 = buffer.data(md + 191);
    const auto *md_192 = buffer.data(md + 192);
    const auto *md_193 = buffer.data(md + 193);
    const auto *md_194 = buffer.data(md + 194);
    const auto *md_195 = buffer.data(md + 195);
    const auto *md_196 = buffer.data(md + 196);
    const auto *md_197 = buffer.data(md + 197);
    const auto *md_198 = buffer.data(md + 198);
    const auto *md_199 = buffer.data(md + 199);
    const auto *md_200 = buffer.data(md + 200);
    const auto *md_201 = buffer.data(md + 201);
    const auto *md_202 = buffer.data(md + 202);
    const auto *md_203 = buffer.data(md + 203);
    const auto *md_204 = buffer.data(md + 204);
    const auto *md_205 = buffer.data(md + 205);
    const auto *md_206 = buffer.data(md + 206);
    const auto *md_207 = buffer.data(md + 207);
    const auto *md_208 = buffer.data(md + 208);
    const auto *md_209 = buffer.data(md + 209);
    const auto *md_210 = buffer.data(md + 210);
    const auto *md_211 = buffer.data(md + 211);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, md_12, md_13, md_14, md_15, \
                         md_16, md_17, md_24, md_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * md_12[k];

        t_1[k] = f_0 * md_13[k];

        t_2[k] = f_0 * md_14[k];

        t_3[k] = f_0 * md_15[k];

        t_4[k] = f_0 * md_16[k];

        t_5[k] = f_0 * md_17[k];

        t_6[k] = f_0 * md_24[k];

        t_7[k] = f_0 * md_25[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, kd_0, kd_1, md_26, md_27, md_28, \
                         md_29, md_30, md_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * md_26[k];

        t_9[k] = f_0 * md_27[k];

        t_10[k] = f_0 * md_28[k];

        t_11[k] = f_0 * md_29[k];

        t_12[k] = -kd_0[k]
                  + f_0 * md_30[k];

        t_13[k] = -kd_1[k]
                  + f_0 * md_31[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, kd_2, kd_3, kd_4, kd_5, md_32, \
                         md_33, md_34, md_35, md_42, md_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -kd_2[k]
                  + f_0 * md_32[k];

        t_15[k] = -kd_3[k]
                  + f_0 * md_33[k];

        t_16[k] = -kd_4[k]
                  + f_0 * md_34[k];

        t_17[k] = -kd_5[k]
                  + f_0 * md_35[k];

        t_18[k] = f_0 * md_42[k];

        t_19[k] = f_0 * md_43[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, kd_6, kd_7, md_44, md_45, md_46, \
                         md_47, md_48, md_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * md_44[k];

        t_21[k] = f_0 * md_45[k];

        t_22[k] = f_0 * md_46[k];

        t_23[k] = f_0 * md_47[k];

        t_24[k] = -kd_6[k]
                  + f_0 * md_48[k];

        t_25[k] = -kd_7[k]
                  + f_0 * md_49[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, kd_8, kd_9, kd_10, kd_11, kd_12, md_50, \
                         md_51, md_52, md_53, md_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -kd_8[k]
                  + f_0 * md_50[k];

        t_27[k] = -kd_9[k]
                  + f_0 * md_51[k];

        t_28[k] = -kd_10[k]
                  + f_0 * md_52[k];

        t_29[k] = -kd_11[k]
                  + f_0 * md_53[k];

        t_30[k] = -2.0 * kd_12[k]
                  + f_0 * md_54[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, kd_13, kd_14, kd_15, kd_16, kd_17, \
                         md_55, md_56, md_57, md_58, md_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -2.0 * kd_13[k]
                  + f_0 * md_55[k];

        t_32[k] = -2.0 * kd_14[k]
                  + f_0 * md_56[k];

        t_33[k] = -2.0 * kd_15[k]
                  + f_0 * md_57[k];

        t_34[k] = -2.0 * kd_16[k]
                  + f_0 * md_58[k];

        t_35[k] = -2.0 * kd_17[k]
                  + f_0 * md_59[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, t_42, kd_18, md_66, md_67, md_68, \
                         md_69, md_70, md_71, md_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * md_66[k];

        t_37[k] = f_0 * md_67[k];

        t_38[k] = f_0 * md_68[k];

        t_39[k] = f_0 * md_69[k];

        t_40[k] = f_0 * md_70[k];

        t_41[k] = f_0 * md_71[k];

        t_42[k] = -kd_18[k]
                  + f_0 * md_72[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, kd_19, kd_20, kd_21, kd_22, kd_23, \
                         md_73, md_74, md_75, md_76, md_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -kd_19[k]
                  + f_0 * md_73[k];

        t_44[k] = -kd_20[k]
                  + f_0 * md_74[k];

        t_45[k] = -kd_21[k]
                  + f_0 * md_75[k];

        t_46[k] = -kd_22[k]
                  + f_0 * md_76[k];

        t_47[k] = -kd_23[k]
                  + f_0 * md_77[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, kd_24, kd_25, kd_26, kd_27, kd_28, \
                         md_78, md_79, md_80, md_81, md_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -2.0 * kd_24[k]
                  + f_0 * md_78[k];

        t_49[k] = -2.0 * kd_25[k]
                  + f_0 * md_79[k];

        t_50[k] = -2.0 * kd_26[k]
                  + f_0 * md_80[k];

        t_51[k] = -2.0 * kd_27[k]
                  + f_0 * md_81[k];

        t_52[k] = -2.0 * kd_28[k]
                  + f_0 * md_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, kd_29, kd_30, kd_31, kd_32, kd_33, \
                         md_83, md_84, md_85, md_86, md_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -2.0 * kd_29[k]
                  + f_0 * md_83[k];

        t_54[k] = -3.0 * kd_30[k]
                  + f_0 * md_84[k];

        t_55[k] = -3.0 * kd_31[k]
                  + f_0 * md_85[k];

        t_56[k] = -3.0 * kd_32[k]
                  + f_0 * md_86[k];

        t_57[k] = -3.0 * kd_33[k]
                  + f_0 * md_87[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, t_63, t_64, kd_34, kd_35, md_88, md_89, \
                         md_96, md_97, md_98, md_99, md_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -3.0 * kd_34[k]
                  + f_0 * md_88[k];

        t_59[k] = -3.0 * kd_35[k]
                  + f_0 * md_89[k];

        t_60[k] = f_0 * md_96[k];

        t_61[k] = f_0 * md_97[k];

        t_62[k] = f_0 * md_98[k];

        t_63[k] = f_0 * md_99[k];

        t_64[k] = f_0 * md_100[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, kd_36, kd_37, kd_38, kd_39, md_101, \
                         md_102, md_103, md_104, md_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_0 * md_101[k];

        t_66[k] = -kd_36[k]
                  + f_0 * md_102[k];

        t_67[k] = -kd_37[k]
                  + f_0 * md_103[k];

        t_68[k] = -kd_38[k]
                  + f_0 * md_104[k];

        t_69[k] = -kd_39[k]
                  + f_0 * md_105[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, kd_40, kd_41, kd_42, kd_43, kd_44, \
                         md_106, md_107, md_108, md_109, md_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -kd_40[k]
                  + f_0 * md_106[k];

        t_71[k] = -kd_41[k]
                  + f_0 * md_107[k];

        t_72[k] = -2.0 * kd_42[k]
                  + f_0 * md_108[k];

        t_73[k] = -2.0 * kd_43[k]
                  + f_0 * md_109[k];

        t_74[k] = -2.0 * kd_44[k]
                  + f_0 * md_110[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, kd_45, kd_46, kd_47, kd_48, kd_49, \
                         md_111, md_112, md_113, md_114, md_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * kd_45[k]
                  + f_0 * md_111[k];

        t_76[k] = -2.0 * kd_46[k]
                  + f_0 * md_112[k];

        t_77[k] = -2.0 * kd_47[k]
                  + f_0 * md_113[k];

        t_78[k] = -3.0 * kd_48[k]
                  + f_0 * md_114[k];

        t_79[k] = -3.0 * kd_49[k]
                  + f_0 * md_115[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, kd_50, kd_51, kd_52, kd_53, kd_54, \
                         md_116, md_117, md_118, md_119, md_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * kd_50[k]
                  + f_0 * md_116[k];

        t_81[k] = -3.0 * kd_51[k]
                  + f_0 * md_117[k];

        t_82[k] = -3.0 * kd_52[k]
                  + f_0 * md_118[k];

        t_83[k] = -3.0 * kd_53[k]
                  + f_0 * md_119[k];

        t_84[k] = -4.0 * kd_54[k]
                  + f_0 * md_120[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, kd_55, kd_56, kd_57, kd_58, kd_59, \
                         md_121, md_122, md_123, md_124, md_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * kd_55[k]
                  + f_0 * md_121[k];

        t_86[k] = -4.0 * kd_56[k]
                  + f_0 * md_122[k];

        t_87[k] = -4.0 * kd_57[k]
                  + f_0 * md_123[k];

        t_88[k] = -4.0 * kd_58[k]
                  + f_0 * md_124[k];

        t_89[k] = -4.0 * kd_59[k]
                  + f_0 * md_125[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, kd_60, md_132, md_133, \
                         md_134, md_135, md_136, md_137, md_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * md_132[k];

        t_91[k] = f_0 * md_133[k];

        t_92[k] = f_0 * md_134[k];

        t_93[k] = f_0 * md_135[k];

        t_94[k] = f_0 * md_136[k];

        t_95[k] = f_0 * md_137[k];

        t_96[k] = -kd_60[k]
                  + f_0 * md_138[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, kd_61, kd_62, kd_63, kd_64, kd_65, \
                         md_139, md_140, md_141, md_142, md_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -kd_61[k]
                  + f_0 * md_139[k];

        t_98[k] = -kd_62[k]
                  + f_0 * md_140[k];

        t_99[k] = -kd_63[k]
                  + f_0 * md_141[k];

        t_100[k] = -kd_64[k]
                   + f_0 * md_142[k];

        t_101[k] = -kd_65[k]
                   + f_0 * md_143[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, kd_66, kd_67, kd_68, kd_69, kd_70, \
                         md_144, md_145, md_146, md_147, md_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -2.0 * kd_66[k]
                   + f_0 * md_144[k];

        t_103[k] = -2.0 * kd_67[k]
                   + f_0 * md_145[k];

        t_104[k] = -2.0 * kd_68[k]
                   + f_0 * md_146[k];

        t_105[k] = -2.0 * kd_69[k]
                   + f_0 * md_147[k];

        t_106[k] = -2.0 * kd_70[k]
                   + f_0 * md_148[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, kd_71, kd_72, kd_73, kd_74, kd_75, \
                         md_149, md_150, md_151, md_152, md_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -2.0 * kd_71[k]
                   + f_0 * md_149[k];

        t_108[k] = -3.0 * kd_72[k]
                   + f_0 * md_150[k];

        t_109[k] = -3.0 * kd_73[k]
                   + f_0 * md_151[k];

        t_110[k] = -3.0 * kd_74[k]
                   + f_0 * md_152[k];

        t_111[k] = -3.0 * kd_75[k]
                   + f_0 * md_153[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, kd_76, kd_77, kd_78, kd_79, kd_80, \
                         md_154, md_155, md_156, md_157, md_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -3.0 * kd_76[k]
                   + f_0 * md_154[k];

        t_113[k] = -3.0 * kd_77[k]
                   + f_0 * md_155[k];

        t_114[k] = -4.0 * kd_78[k]
                   + f_0 * md_156[k];

        t_115[k] = -4.0 * kd_79[k]
                   + f_0 * md_157[k];

        t_116[k] = -4.0 * kd_80[k]
                   + f_0 * md_158[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, kd_81, kd_82, kd_83, kd_84, kd_85, \
                         md_159, md_160, md_161, md_162, md_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -4.0 * kd_81[k]
                   + f_0 * md_159[k];

        t_118[k] = -4.0 * kd_82[k]
                   + f_0 * md_160[k];

        t_119[k] = -4.0 * kd_83[k]
                   + f_0 * md_161[k];

        t_120[k] = -5.0 * kd_84[k]
                   + f_0 * md_162[k];

        t_121[k] = -5.0 * kd_85[k]
                   + f_0 * md_163[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, t_127, kd_86, kd_87, kd_88, kd_89, \
                         md_164, md_165, md_166, md_167, md_174, \
                         md_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -5.0 * kd_86[k]
                   + f_0 * md_164[k];

        t_123[k] = -5.0 * kd_87[k]
                   + f_0 * md_165[k];

        t_124[k] = -5.0 * kd_88[k]
                   + f_0 * md_166[k];

        t_125[k] = -5.0 * kd_89[k]
                   + f_0 * md_167[k];

        t_126[k] = f_0 * md_174[k];

        t_127[k] = f_0 * md_175[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, t_133, kd_90, kd_91, md_176, \
                         md_177, md_178, md_179, md_180, md_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_0 * md_176[k];

        t_129[k] = f_0 * md_177[k];

        t_130[k] = f_0 * md_178[k];

        t_131[k] = f_0 * md_179[k];

        t_132[k] = -kd_90[k]
                   + f_0 * md_180[k];

        t_133[k] = -kd_91[k]
                   + f_0 * md_181[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, kd_92, kd_93, kd_94, kd_95, kd_96, \
                         md_182, md_183, md_184, md_185, md_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -kd_92[k]
                   + f_0 * md_182[k];

        t_135[k] = -kd_93[k]
                   + f_0 * md_183[k];

        t_136[k] = -kd_94[k]
                   + f_0 * md_184[k];

        t_137[k] = -kd_95[k]
                   + f_0 * md_185[k];

        t_138[k] = -2.0 * kd_96[k]
                   + f_0 * md_186[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, kd_97, kd_98, kd_99, kd_100, \
                         kd_101, md_187, md_188, md_189, md_190, \
                         md_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -2.0 * kd_97[k]
                   + f_0 * md_187[k];

        t_140[k] = -2.0 * kd_98[k]
                   + f_0 * md_188[k];

        t_141[k] = -2.0 * kd_99[k]
                   + f_0 * md_189[k];

        t_142[k] = -2.0 * kd_100[k]
                   + f_0 * md_190[k];

        t_143[k] = -2.0 * kd_101[k]
                   + f_0 * md_191[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, kd_102, kd_103, kd_104, kd_105, \
                         kd_106, md_192, md_193, md_194, md_195, \
                         md_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -3.0 * kd_102[k]
                   + f_0 * md_192[k];

        t_145[k] = -3.0 * kd_103[k]
                   + f_0 * md_193[k];

        t_146[k] = -3.0 * kd_104[k]
                   + f_0 * md_194[k];

        t_147[k] = -3.0 * kd_105[k]
                   + f_0 * md_195[k];

        t_148[k] = -3.0 * kd_106[k]
                   + f_0 * md_196[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, kd_107, kd_108, kd_109, kd_110, \
                         kd_111, md_197, md_198, md_199, md_200, \
                         md_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -3.0 * kd_107[k]
                   + f_0 * md_197[k];

        t_150[k] = -4.0 * kd_108[k]
                   + f_0 * md_198[k];

        t_151[k] = -4.0 * kd_109[k]
                   + f_0 * md_199[k];

        t_152[k] = -4.0 * kd_110[k]
                   + f_0 * md_200[k];

        t_153[k] = -4.0 * kd_111[k]
                   + f_0 * md_201[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, kd_112, kd_113, kd_114, kd_115, \
                         kd_116, md_202, md_203, md_204, md_205, \
                         md_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = -4.0 * kd_112[k]
                   + f_0 * md_202[k];

        t_155[k] = -4.0 * kd_113[k]
                   + f_0 * md_203[k];

        t_156[k] = -5.0 * kd_114[k]
                   + f_0 * md_204[k];

        t_157[k] = -5.0 * kd_115[k]
                   + f_0 * md_205[k];

        t_158[k] = -5.0 * kd_116[k]
                   + f_0 * md_206[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, kd_117, kd_118, kd_119, kd_120, \
                         kd_121, md_207, md_208, md_209, md_210, \
                         md_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = -5.0 * kd_117[k]
                   + f_0 * md_207[k];

        t_160[k] = -5.0 * kd_118[k]
                   + f_0 * md_208[k];

        t_161[k] = -5.0 * kd_119[k]
                   + f_0 * md_209[k];

        t_162[k] = -6.0 * kd_120[k]
                   + f_0 * md_210[k];

        t_163[k] = -6.0 * kd_121[k]
                   + f_0 * md_211[k];
    }
}

static auto
compute_prim_geom_10_ld_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t kd, const size_t md,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);
    const auto *kd_141 = buffer.data(kd + 141);
    const auto *kd_142 = buffer.data(kd + 142);
    const auto *kd_143 = buffer.data(kd + 143);
    const auto *kd_144 = buffer.data(kd + 144);
    const auto *kd_145 = buffer.data(kd + 145);
    const auto *kd_146 = buffer.data(kd + 146);
    const auto *kd_147 = buffer.data(kd + 147);
    const auto *kd_148 = buffer.data(kd + 148);
    const auto *kd_149 = buffer.data(kd + 149);
    const auto *kd_150 = buffer.data(kd + 150);
    const auto *kd_151 = buffer.data(kd + 151);
    const auto *kd_152 = buffer.data(kd + 152);
    const auto *kd_153 = buffer.data(kd + 153);
    const auto *kd_154 = buffer.data(kd + 154);
    const auto *kd_155 = buffer.data(kd + 155);
    const auto *kd_156 = buffer.data(kd + 156);
    const auto *kd_157 = buffer.data(kd + 157);
    const auto *kd_158 = buffer.data(kd + 158);
    const auto *kd_159 = buffer.data(kd + 159);
    const auto *kd_160 = buffer.data(kd + 160);
    const auto *kd_161 = buffer.data(kd + 161);
    const auto *kd_162 = buffer.data(kd + 162);
    const auto *kd_163 = buffer.data(kd + 163);
    const auto *kd_164 = buffer.data(kd + 164);
    const auto *kd_165 = buffer.data(kd + 165);
    const auto *kd_166 = buffer.data(kd + 166);
    const auto *kd_167 = buffer.data(kd + 167);
    const auto *kd_168 = buffer.data(kd + 168);
    const auto *kd_169 = buffer.data(kd + 169);
    const auto *kd_170 = buffer.data(kd + 170);
    const auto *kd_171 = buffer.data(kd + 171);
    const auto *kd_172 = buffer.data(kd + 172);
    const auto *kd_173 = buffer.data(kd + 173);
    const auto *kd_174 = buffer.data(kd + 174);
    const auto *kd_175 = buffer.data(kd + 175);
    const auto *kd_176 = buffer.data(kd + 176);
    const auto *kd_177 = buffer.data(kd + 177);
    const auto *kd_178 = buffer.data(kd + 178);
    const auto *kd_179 = buffer.data(kd + 179);
    const auto *kd_180 = buffer.data(kd + 180);
    const auto *kd_181 = buffer.data(kd + 181);
    const auto *kd_182 = buffer.data(kd + 182);
    const auto *kd_183 = buffer.data(kd + 183);
    const auto *kd_184 = buffer.data(kd + 184);
    const auto *kd_185 = buffer.data(kd + 185);
    const auto *kd_186 = buffer.data(kd + 186);
    const auto *kd_187 = buffer.data(kd + 187);
    const auto *kd_188 = buffer.data(kd + 188);
    const auto *kd_189 = buffer.data(kd + 189);
    const auto *kd_190 = buffer.data(kd + 190);
    const auto *kd_191 = buffer.data(kd + 191);
    const auto *kd_192 = buffer.data(kd + 192);
    const auto *kd_193 = buffer.data(kd + 193);
    const auto *kd_194 = buffer.data(kd + 194);
    const auto *kd_195 = buffer.data(kd + 195);
    const auto *kd_196 = buffer.data(kd + 196);
    const auto *kd_197 = buffer.data(kd + 197);
    const auto *kd_198 = buffer.data(kd + 198);
    const auto *kd_199 = buffer.data(kd + 199);
    const auto *kd_200 = buffer.data(kd + 200);
    const auto *kd_201 = buffer.data(kd + 201);
    const auto *kd_202 = buffer.data(kd + 202);
    const auto *kd_203 = buffer.data(kd + 203);
    const auto *kd_204 = buffer.data(kd + 204);
    const auto *kd_205 = buffer.data(kd + 205);
    const auto *kd_206 = buffer.data(kd + 206);
    const auto *kd_207 = buffer.data(kd + 207);
    const auto *kd_208 = buffer.data(kd + 208);
    const auto *kd_209 = buffer.data(kd + 209);
    const auto *kd_210 = buffer.data(kd + 210);
    const auto *kd_211 = buffer.data(kd + 211);
    const auto *kd_212 = buffer.data(kd + 212);
    const auto *kd_213 = buffer.data(kd + 213);
    const auto *kd_214 = buffer.data(kd + 214);
    const auto *kd_215 = buffer.data(kd + 215);

    const auto *md_212 = buffer.data(md + 212);
    const auto *md_213 = buffer.data(md + 213);
    const auto *md_214 = buffer.data(md + 214);
    const auto *md_215 = buffer.data(md + 215);
    const auto *md_222 = buffer.data(md + 222);
    const auto *md_223 = buffer.data(md + 223);
    const auto *md_224 = buffer.data(md + 224);
    const auto *md_225 = buffer.data(md + 225);
    const auto *md_226 = buffer.data(md + 226);
    const auto *md_227 = buffer.data(md + 227);
    const auto *md_228 = buffer.data(md + 228);
    const auto *md_229 = buffer.data(md + 229);
    const auto *md_230 = buffer.data(md + 230);
    const auto *md_231 = buffer.data(md + 231);
    const auto *md_232 = buffer.data(md + 232);
    const auto *md_233 = buffer.data(md + 233);
    const auto *md_234 = buffer.data(md + 234);
    const auto *md_235 = buffer.data(md + 235);
    const auto *md_236 = buffer.data(md + 236);
    const auto *md_237 = buffer.data(md + 237);
    const auto *md_238 = buffer.data(md + 238);
    const auto *md_239 = buffer.data(md + 239);
    const auto *md_240 = buffer.data(md + 240);
    const auto *md_241 = buffer.data(md + 241);
    const auto *md_242 = buffer.data(md + 242);
    const auto *md_243 = buffer.data(md + 243);
    const auto *md_244 = buffer.data(md + 244);
    const auto *md_245 = buffer.data(md + 245);
    const auto *md_246 = buffer.data(md + 246);
    const auto *md_247 = buffer.data(md + 247);
    const auto *md_248 = buffer.data(md + 248);
    const auto *md_249 = buffer.data(md + 249);
    const auto *md_250 = buffer.data(md + 250);
    const auto *md_251 = buffer.data(md + 251);
    const auto *md_252 = buffer.data(md + 252);
    const auto *md_253 = buffer.data(md + 253);
    const auto *md_254 = buffer.data(md + 254);
    const auto *md_255 = buffer.data(md + 255);
    const auto *md_256 = buffer.data(md + 256);
    const auto *md_257 = buffer.data(md + 257);
    const auto *md_258 = buffer.data(md + 258);
    const auto *md_259 = buffer.data(md + 259);
    const auto *md_260 = buffer.data(md + 260);
    const auto *md_261 = buffer.data(md + 261);
    const auto *md_262 = buffer.data(md + 262);
    const auto *md_263 = buffer.data(md + 263);
    const auto *md_264 = buffer.data(md + 264);
    const auto *md_265 = buffer.data(md + 265);
    const auto *md_266 = buffer.data(md + 266);
    const auto *md_267 = buffer.data(md + 267);
    const auto *md_268 = buffer.data(md + 268);
    const auto *md_269 = buffer.data(md + 269);
    const auto *md_276 = buffer.data(md + 276);
    const auto *md_277 = buffer.data(md + 277);
    const auto *md_278 = buffer.data(md + 278);
    const auto *md_279 = buffer.data(md + 279);
    const auto *md_280 = buffer.data(md + 280);
    const auto *md_281 = buffer.data(md + 281);
    const auto *md_282 = buffer.data(md + 282);
    const auto *md_283 = buffer.data(md + 283);
    const auto *md_284 = buffer.data(md + 284);
    const auto *md_285 = buffer.data(md + 285);
    const auto *md_286 = buffer.data(md + 286);
    const auto *md_287 = buffer.data(md + 287);
    const auto *md_288 = buffer.data(md + 288);
    const auto *md_289 = buffer.data(md + 289);
    const auto *md_290 = buffer.data(md + 290);
    const auto *md_291 = buffer.data(md + 291);
    const auto *md_292 = buffer.data(md + 292);
    const auto *md_293 = buffer.data(md + 293);
    const auto *md_294 = buffer.data(md + 294);
    const auto *md_295 = buffer.data(md + 295);
    const auto *md_296 = buffer.data(md + 296);
    const auto *md_297 = buffer.data(md + 297);
    const auto *md_298 = buffer.data(md + 298);
    const auto *md_299 = buffer.data(md + 299);
    const auto *md_300 = buffer.data(md + 300);
    const auto *md_301 = buffer.data(md + 301);
    const auto *md_302 = buffer.data(md + 302);
    const auto *md_303 = buffer.data(md + 303);
    const auto *md_304 = buffer.data(md + 304);
    const auto *md_305 = buffer.data(md + 305);
    const auto *md_306 = buffer.data(md + 306);
    const auto *md_307 = buffer.data(md + 307);
    const auto *md_308 = buffer.data(md + 308);
    const auto *md_309 = buffer.data(md + 309);
    const auto *md_310 = buffer.data(md + 310);
    const auto *md_311 = buffer.data(md + 311);
    const auto *md_312 = buffer.data(md + 312);
    const auto *md_313 = buffer.data(md + 313);
    const auto *md_314 = buffer.data(md + 314);
    const auto *md_315 = buffer.data(md + 315);
    const auto *md_316 = buffer.data(md + 316);
    const auto *md_317 = buffer.data(md + 317);
    const auto *md_318 = buffer.data(md + 318);
    const auto *md_319 = buffer.data(md + 319);
    const auto *md_320 = buffer.data(md + 320);
    const auto *md_321 = buffer.data(md + 321);
    const auto *md_322 = buffer.data(md + 322);
    const auto *md_323 = buffer.data(md + 323);
    const auto *md_324 = buffer.data(md + 324);
    const auto *md_325 = buffer.data(md + 325);
    const auto *md_326 = buffer.data(md + 326);
    const auto *md_327 = buffer.data(md + 327);
    const auto *md_328 = buffer.data(md + 328);
    const auto *md_329 = buffer.data(md + 329);

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, kd_122, kd_123, kd_124, \
                         kd_125, md_212, md_213, md_214, md_215, md_222, \
                         md_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -6.0 * kd_122[k]
                   + f_0 * md_212[k];

        t_165[k] = -6.0 * kd_123[k]
                   + f_0 * md_213[k];

        t_166[k] = -6.0 * kd_124[k]
                   + f_0 * md_214[k];

        t_167[k] = -6.0 * kd_125[k]
                   + f_0 * md_215[k];

        t_168[k] = f_0 * md_222[k];

        t_169[k] = f_0 * md_223[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, t_175, kd_126, kd_127, md_224, \
                         md_225, md_226, md_227, md_228, md_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_0 * md_224[k];

        t_171[k] = f_0 * md_225[k];

        t_172[k] = f_0 * md_226[k];

        t_173[k] = f_0 * md_227[k];

        t_174[k] = -kd_126[k]
                   + f_0 * md_228[k];

        t_175[k] = -kd_127[k]
                   + f_0 * md_229[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, t_180, kd_128, kd_129, kd_130, kd_131, \
                         kd_132, md_230, md_231, md_232, md_233, \
                         md_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = -kd_128[k]
                   + f_0 * md_230[k];

        t_177[k] = -kd_129[k]
                   + f_0 * md_231[k];

        t_178[k] = -kd_130[k]
                   + f_0 * md_232[k];

        t_179[k] = -kd_131[k]
                   + f_0 * md_233[k];

        t_180[k] = -2.0 * kd_132[k]
                   + f_0 * md_234[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, kd_133, kd_134, kd_135, kd_136, \
                         kd_137, md_235, md_236, md_237, md_238, \
                         md_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = -2.0 * kd_133[k]
                   + f_0 * md_235[k];

        t_182[k] = -2.0 * kd_134[k]
                   + f_0 * md_236[k];

        t_183[k] = -2.0 * kd_135[k]
                   + f_0 * md_237[k];

        t_184[k] = -2.0 * kd_136[k]
                   + f_0 * md_238[k];

        t_185[k] = -2.0 * kd_137[k]
                   + f_0 * md_239[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, kd_138, kd_139, kd_140, kd_141, \
                         kd_142, md_240, md_241, md_242, md_243, \
                         md_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = -3.0 * kd_138[k]
                   + f_0 * md_240[k];

        t_187[k] = -3.0 * kd_139[k]
                   + f_0 * md_241[k];

        t_188[k] = -3.0 * kd_140[k]
                   + f_0 * md_242[k];

        t_189[k] = -3.0 * kd_141[k]
                   + f_0 * md_243[k];

        t_190[k] = -3.0 * kd_142[k]
                   + f_0 * md_244[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, kd_143, kd_144, kd_145, kd_146, \
                         kd_147, md_245, md_246, md_247, md_248, \
                         md_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = -3.0 * kd_143[k]
                   + f_0 * md_245[k];

        t_192[k] = -4.0 * kd_144[k]
                   + f_0 * md_246[k];

        t_193[k] = -4.0 * kd_145[k]
                   + f_0 * md_247[k];

        t_194[k] = -4.0 * kd_146[k]
                   + f_0 * md_248[k];

        t_195[k] = -4.0 * kd_147[k]
                   + f_0 * md_249[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, kd_148, kd_149, kd_150, kd_151, \
                         kd_152, md_250, md_251, md_252, md_253, \
                         md_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -4.0 * kd_148[k]
                   + f_0 * md_250[k];

        t_197[k] = -4.0 * kd_149[k]
                   + f_0 * md_251[k];

        t_198[k] = -5.0 * kd_150[k]
                   + f_0 * md_252[k];

        t_199[k] = -5.0 * kd_151[k]
                   + f_0 * md_253[k];

        t_200[k] = -5.0 * kd_152[k]
                   + f_0 * md_254[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, kd_153, kd_154, kd_155, kd_156, \
                         kd_157, md_255, md_256, md_257, md_258, \
                         md_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -5.0 * kd_153[k]
                   + f_0 * md_255[k];

        t_202[k] = -5.0 * kd_154[k]
                   + f_0 * md_256[k];

        t_203[k] = -5.0 * kd_155[k]
                   + f_0 * md_257[k];

        t_204[k] = -6.0 * kd_156[k]
                   + f_0 * md_258[k];

        t_205[k] = -6.0 * kd_157[k]
                   + f_0 * md_259[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, kd_158, kd_159, kd_160, kd_161, \
                         kd_162, md_260, md_261, md_262, md_263, \
                         md_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = -6.0 * kd_158[k]
                   + f_0 * md_260[k];

        t_207[k] = -6.0 * kd_159[k]
                   + f_0 * md_261[k];

        t_208[k] = -6.0 * kd_160[k]
                   + f_0 * md_262[k];

        t_209[k] = -6.0 * kd_161[k]
                   + f_0 * md_263[k];

        t_210[k] = -7.0 * kd_162[k]
                   + f_0 * md_264[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, kd_163, kd_164, kd_165, kd_166, \
                         kd_167, md_265, md_266, md_267, md_268, \
                         md_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -7.0 * kd_163[k]
                   + f_0 * md_265[k];

        t_212[k] = -7.0 * kd_164[k]
                   + f_0 * md_266[k];

        t_213[k] = -7.0 * kd_165[k]
                   + f_0 * md_267[k];

        t_214[k] = -7.0 * kd_166[k]
                   + f_0 * md_268[k];

        t_215[k] = -7.0 * kd_167[k]
                   + f_0 * md_269[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, t_221, t_222, kd_168, md_276, \
                         md_277, md_278, md_279, md_280, md_281, \
                         md_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_0 * md_276[k];

        t_217[k] = f_0 * md_277[k];

        t_218[k] = f_0 * md_278[k];

        t_219[k] = f_0 * md_279[k];

        t_220[k] = f_0 * md_280[k];

        t_221[k] = f_0 * md_281[k];

        t_222[k] = -kd_168[k]
                   + f_0 * md_282[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, kd_169, kd_170, kd_171, kd_172, \
                         kd_173, md_283, md_284, md_285, md_286, \
                         md_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = -kd_169[k]
                   + f_0 * md_283[k];

        t_224[k] = -kd_170[k]
                   + f_0 * md_284[k];

        t_225[k] = -kd_171[k]
                   + f_0 * md_285[k];

        t_226[k] = -kd_172[k]
                   + f_0 * md_286[k];

        t_227[k] = -kd_173[k]
                   + f_0 * md_287[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, kd_174, kd_175, kd_176, kd_177, \
                         kd_178, md_288, md_289, md_290, md_291, \
                         md_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = -2.0 * kd_174[k]
                   + f_0 * md_288[k];

        t_229[k] = -2.0 * kd_175[k]
                   + f_0 * md_289[k];

        t_230[k] = -2.0 * kd_176[k]
                   + f_0 * md_290[k];

        t_231[k] = -2.0 * kd_177[k]
                   + f_0 * md_291[k];

        t_232[k] = -2.0 * kd_178[k]
                   + f_0 * md_292[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, t_237, kd_179, kd_180, kd_181, kd_182, \
                         kd_183, md_293, md_294, md_295, md_296, \
                         md_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = -2.0 * kd_179[k]
                   + f_0 * md_293[k];

        t_234[k] = -3.0 * kd_180[k]
                   + f_0 * md_294[k];

        t_235[k] = -3.0 * kd_181[k]
                   + f_0 * md_295[k];

        t_236[k] = -3.0 * kd_182[k]
                   + f_0 * md_296[k];

        t_237[k] = -3.0 * kd_183[k]
                   + f_0 * md_297[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, kd_184, kd_185, kd_186, kd_187, \
                         kd_188, md_298, md_299, md_300, md_301, \
                         md_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = -3.0 * kd_184[k]
                   + f_0 * md_298[k];

        t_239[k] = -3.0 * kd_185[k]
                   + f_0 * md_299[k];

        t_240[k] = -4.0 * kd_186[k]
                   + f_0 * md_300[k];

        t_241[k] = -4.0 * kd_187[k]
                   + f_0 * md_301[k];

        t_242[k] = -4.0 * kd_188[k]
                   + f_0 * md_302[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, t_247, kd_189, kd_190, kd_191, kd_192, \
                         kd_193, md_303, md_304, md_305, md_306, \
                         md_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = -4.0 * kd_189[k]
                   + f_0 * md_303[k];

        t_244[k] = -4.0 * kd_190[k]
                   + f_0 * md_304[k];

        t_245[k] = -4.0 * kd_191[k]
                   + f_0 * md_305[k];

        t_246[k] = -5.0 * kd_192[k]
                   + f_0 * md_306[k];

        t_247[k] = -5.0 * kd_193[k]
                   + f_0 * md_307[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, t_252, kd_194, kd_195, kd_196, kd_197, \
                         kd_198, md_308, md_309, md_310, md_311, \
                         md_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = -5.0 * kd_194[k]
                   + f_0 * md_308[k];

        t_249[k] = -5.0 * kd_195[k]
                   + f_0 * md_309[k];

        t_250[k] = -5.0 * kd_196[k]
                   + f_0 * md_310[k];

        t_251[k] = -5.0 * kd_197[k]
                   + f_0 * md_311[k];

        t_252[k] = -6.0 * kd_198[k]
                   + f_0 * md_312[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, kd_199, kd_200, kd_201, kd_202, \
                         kd_203, md_313, md_314, md_315, md_316, \
                         md_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = -6.0 * kd_199[k]
                   + f_0 * md_313[k];

        t_254[k] = -6.0 * kd_200[k]
                   + f_0 * md_314[k];

        t_255[k] = -6.0 * kd_201[k]
                   + f_0 * md_315[k];

        t_256[k] = -6.0 * kd_202[k]
                   + f_0 * md_316[k];

        t_257[k] = -6.0 * kd_203[k]
                   + f_0 * md_317[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, kd_204, kd_205, kd_206, kd_207, \
                         kd_208, md_318, md_319, md_320, md_321, \
                         md_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -7.0 * kd_204[k]
                   + f_0 * md_318[k];

        t_259[k] = -7.0 * kd_205[k]
                   + f_0 * md_319[k];

        t_260[k] = -7.0 * kd_206[k]
                   + f_0 * md_320[k];

        t_261[k] = -7.0 * kd_207[k]
                   + f_0 * md_321[k];

        t_262[k] = -7.0 * kd_208[k]
                   + f_0 * md_322[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, kd_209, kd_210, kd_211, kd_212, \
                         kd_213, md_323, md_324, md_325, md_326, \
                         md_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = -7.0 * kd_209[k]
                   + f_0 * md_323[k];

        t_264[k] = -8.0 * kd_210[k]
                   + f_0 * md_324[k];

        t_265[k] = -8.0 * kd_211[k]
                   + f_0 * md_325[k];

        t_266[k] = -8.0 * kd_212[k]
                   + f_0 * md_326[k];

        t_267[k] = -8.0 * kd_213[k]
                   + f_0 * md_327[k];
    }

#pragma omp simd aligned(t_268, t_269, kd_214, kd_215, md_328, md_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = -8.0 * kd_214[k]
                   + f_0 * md_328[k];

        t_269[k] = -8.0 * kd_215[k]
                   + f_0 * md_329[k];
    }
}

auto
compute_prim_geom_10_ld_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t kd, const size_t md,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_ld_electron_repulsion_2_piece0(buffer, target, kd, md, ncols, alpha);

    compute_prim_geom_10_ld_electron_repulsion_2_piece1(buffer, target, kd, md, ncols, alpha);
}

}  // namespace simdt2ceri
