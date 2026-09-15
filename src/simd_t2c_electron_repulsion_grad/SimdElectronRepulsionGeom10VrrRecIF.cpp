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


#include "SimdElectronRepulsionGeom10VrrRecIF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_geom_10_if_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hf, const size_t kf,
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

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_102 = buffer.data(hf + 102);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_104 = buffer.data(hf + 104);
    const auto *hf_105 = buffer.data(hf + 105);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_107 = buffer.data(hf + 107);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_111 = buffer.data(hf + 111);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_113 = buffer.data(hf + 113);
    const auto *hf_114 = buffer.data(hf + 114);
    const auto *hf_115 = buffer.data(hf + 115);
    const auto *hf_116 = buffer.data(hf + 116);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_121 = buffer.data(hf + 121);
    const auto *hf_122 = buffer.data(hf + 122);
    const auto *hf_123 = buffer.data(hf + 123);
    const auto *hf_124 = buffer.data(hf + 124);
    const auto *hf_125 = buffer.data(hf + 125);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_127 = buffer.data(hf + 127);
    const auto *hf_128 = buffer.data(hf + 128);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_130 = buffer.data(hf + 130);
    const auto *hf_131 = buffer.data(hf + 131);
    const auto *hf_132 = buffer.data(hf + 132);
    const auto *hf_133 = buffer.data(hf + 133);
    const auto *hf_134 = buffer.data(hf + 134);
    const auto *hf_135 = buffer.data(hf + 135);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_137 = buffer.data(hf + 137);
    const auto *hf_138 = buffer.data(hf + 138);
    const auto *hf_139 = buffer.data(hf + 139);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_141 = buffer.data(hf + 141);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_143 = buffer.data(hf + 143);
    const auto *hf_144 = buffer.data(hf + 144);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_146 = buffer.data(hf + 146);
    const auto *hf_147 = buffer.data(hf + 147);
    const auto *hf_148 = buffer.data(hf + 148);
    const auto *hf_149 = buffer.data(hf + 149);

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_1 = buffer.data(kf + 1);
    const auto *kf_2 = buffer.data(kf + 2);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_4 = buffer.data(kf + 4);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_12 = buffer.data(kf + 12);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_15 = buffer.data(kf + 15);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_53 = buffer.data(kf + 53);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_65 = buffer.data(kf + 65);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_93 = buffer.data(kf + 93);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_98 = buffer.data(kf + 98);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_104 = buffer.data(kf + 104);
    const auto *kf_105 = buffer.data(kf + 105);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_113 = buffer.data(kf + 113);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_115 = buffer.data(kf + 115);
    const auto *kf_116 = buffer.data(kf + 116);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_121 = buffer.data(kf + 121);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_123 = buffer.data(kf + 123);
    const auto *kf_124 = buffer.data(kf + 124);
    const auto *kf_125 = buffer.data(kf + 125);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_131 = buffer.data(kf + 131);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_135 = buffer.data(kf + 135);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_143 = buffer.data(kf + 143);
    const auto *kf_144 = buffer.data(kf + 144);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_147 = buffer.data(kf + 147);
    const auto *kf_148 = buffer.data(kf + 148);
    const auto *kf_149 = buffer.data(kf + 149);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, hf_0, hf_1, hf_2, hf_3, hf_4, kf_0, kf_1, \
                         kf_2, kf_3, kf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -6.0 * hf_0[k]
                 + f_0 * kf_0[k];

        t_1[k] = -6.0 * hf_1[k]
                 + f_0 * kf_1[k];

        t_2[k] = -6.0 * hf_2[k]
                 + f_0 * kf_2[k];

        t_3[k] = -6.0 * hf_3[k]
                 + f_0 * kf_3[k];

        t_4[k] = -6.0 * hf_4[k]
                 + f_0 * kf_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, hf_5, hf_6, hf_7, hf_8, hf_9, kf_5, kf_6, \
                         kf_7, kf_8, kf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -6.0 * hf_5[k]
                 + f_0 * kf_5[k];

        t_6[k] = -6.0 * hf_6[k]
                 + f_0 * kf_6[k];

        t_7[k] = -6.0 * hf_7[k]
                 + f_0 * kf_7[k];

        t_8[k] = -6.0 * hf_8[k]
                 + f_0 * kf_8[k];

        t_9[k] = -6.0 * hf_9[k]
                 + f_0 * kf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, hf_10, hf_11, hf_12, hf_13, hf_14, \
                         kf_10, kf_11, kf_12, kf_13, kf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -5.0 * hf_10[k]
                  + f_0 * kf_10[k];

        t_11[k] = -5.0 * hf_11[k]
                  + f_0 * kf_11[k];

        t_12[k] = -5.0 * hf_12[k]
                  + f_0 * kf_12[k];

        t_13[k] = -5.0 * hf_13[k]
                  + f_0 * kf_13[k];

        t_14[k] = -5.0 * hf_14[k]
                  + f_0 * kf_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, hf_15, hf_16, hf_17, hf_18, hf_19, \
                         kf_15, kf_16, kf_17, kf_18, kf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -5.0 * hf_15[k]
                  + f_0 * kf_15[k];

        t_16[k] = -5.0 * hf_16[k]
                  + f_0 * kf_16[k];

        t_17[k] = -5.0 * hf_17[k]
                  + f_0 * kf_17[k];

        t_18[k] = -5.0 * hf_18[k]
                  + f_0 * kf_18[k];

        t_19[k] = -5.0 * hf_19[k]
                  + f_0 * kf_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, hf_20, hf_21, hf_22, hf_23, hf_24, \
                         kf_20, kf_21, kf_22, kf_23, kf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -5.0 * hf_20[k]
                  + f_0 * kf_20[k];

        t_21[k] = -5.0 * hf_21[k]
                  + f_0 * kf_21[k];

        t_22[k] = -5.0 * hf_22[k]
                  + f_0 * kf_22[k];

        t_23[k] = -5.0 * hf_23[k]
                  + f_0 * kf_23[k];

        t_24[k] = -5.0 * hf_24[k]
                  + f_0 * kf_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, hf_25, hf_26, hf_27, hf_28, hf_29, \
                         kf_25, kf_26, kf_27, kf_28, kf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -5.0 * hf_25[k]
                  + f_0 * kf_25[k];

        t_26[k] = -5.0 * hf_26[k]
                  + f_0 * kf_26[k];

        t_27[k] = -5.0 * hf_27[k]
                  + f_0 * kf_27[k];

        t_28[k] = -5.0 * hf_28[k]
                  + f_0 * kf_28[k];

        t_29[k] = -5.0 * hf_29[k]
                  + f_0 * kf_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, hf_30, hf_31, hf_32, hf_33, hf_34, \
                         kf_30, kf_31, kf_32, kf_33, kf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -4.0 * hf_30[k]
                  + f_0 * kf_30[k];

        t_31[k] = -4.0 * hf_31[k]
                  + f_0 * kf_31[k];

        t_32[k] = -4.0 * hf_32[k]
                  + f_0 * kf_32[k];

        t_33[k] = -4.0 * hf_33[k]
                  + f_0 * kf_33[k];

        t_34[k] = -4.0 * hf_34[k]
                  + f_0 * kf_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, hf_35, hf_36, hf_37, hf_38, hf_39, \
                         kf_35, kf_36, kf_37, kf_38, kf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -4.0 * hf_35[k]
                  + f_0 * kf_35[k];

        t_36[k] = -4.0 * hf_36[k]
                  + f_0 * kf_36[k];

        t_37[k] = -4.0 * hf_37[k]
                  + f_0 * kf_37[k];

        t_38[k] = -4.0 * hf_38[k]
                  + f_0 * kf_38[k];

        t_39[k] = -4.0 * hf_39[k]
                  + f_0 * kf_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hf_40, hf_41, hf_42, hf_43, hf_44, \
                         kf_40, kf_41, kf_42, kf_43, kf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -4.0 * hf_40[k]
                  + f_0 * kf_40[k];

        t_41[k] = -4.0 * hf_41[k]
                  + f_0 * kf_41[k];

        t_42[k] = -4.0 * hf_42[k]
                  + f_0 * kf_42[k];

        t_43[k] = -4.0 * hf_43[k]
                  + f_0 * kf_43[k];

        t_44[k] = -4.0 * hf_44[k]
                  + f_0 * kf_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hf_45, hf_46, hf_47, hf_48, hf_49, \
                         kf_45, kf_46, kf_47, kf_48, kf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -4.0 * hf_45[k]
                  + f_0 * kf_45[k];

        t_46[k] = -4.0 * hf_46[k]
                  + f_0 * kf_46[k];

        t_47[k] = -4.0 * hf_47[k]
                  + f_0 * kf_47[k];

        t_48[k] = -4.0 * hf_48[k]
                  + f_0 * kf_48[k];

        t_49[k] = -4.0 * hf_49[k]
                  + f_0 * kf_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hf_50, hf_51, hf_52, hf_53, hf_54, \
                         kf_50, kf_51, kf_52, kf_53, kf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * hf_50[k]
                  + f_0 * kf_50[k];

        t_51[k] = -4.0 * hf_51[k]
                  + f_0 * kf_51[k];

        t_52[k] = -4.0 * hf_52[k]
                  + f_0 * kf_52[k];

        t_53[k] = -4.0 * hf_53[k]
                  + f_0 * kf_53[k];

        t_54[k] = -4.0 * hf_54[k]
                  + f_0 * kf_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hf_55, hf_56, hf_57, hf_58, hf_59, \
                         kf_55, kf_56, kf_57, kf_58, kf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -4.0 * hf_55[k]
                  + f_0 * kf_55[k];

        t_56[k] = -4.0 * hf_56[k]
                  + f_0 * kf_56[k];

        t_57[k] = -4.0 * hf_57[k]
                  + f_0 * kf_57[k];

        t_58[k] = -4.0 * hf_58[k]
                  + f_0 * kf_58[k];

        t_59[k] = -4.0 * hf_59[k]
                  + f_0 * kf_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, hf_60, hf_61, hf_62, hf_63, hf_64, \
                         kf_60, kf_61, kf_62, kf_63, kf_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * hf_60[k]
                  + f_0 * kf_60[k];

        t_61[k] = -3.0 * hf_61[k]
                  + f_0 * kf_61[k];

        t_62[k] = -3.0 * hf_62[k]
                  + f_0 * kf_62[k];

        t_63[k] = -3.0 * hf_63[k]
                  + f_0 * kf_63[k];

        t_64[k] = -3.0 * hf_64[k]
                  + f_0 * kf_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hf_65, hf_66, hf_67, hf_68, hf_69, \
                         kf_65, kf_66, kf_67, kf_68, kf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * hf_65[k]
                  + f_0 * kf_65[k];

        t_66[k] = -3.0 * hf_66[k]
                  + f_0 * kf_66[k];

        t_67[k] = -3.0 * hf_67[k]
                  + f_0 * kf_67[k];

        t_68[k] = -3.0 * hf_68[k]
                  + f_0 * kf_68[k];

        t_69[k] = -3.0 * hf_69[k]
                  + f_0 * kf_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hf_70, hf_71, hf_72, hf_73, hf_74, \
                         kf_70, kf_71, kf_72, kf_73, kf_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -3.0 * hf_70[k]
                  + f_0 * kf_70[k];

        t_71[k] = -3.0 * hf_71[k]
                  + f_0 * kf_71[k];

        t_72[k] = -3.0 * hf_72[k]
                  + f_0 * kf_72[k];

        t_73[k] = -3.0 * hf_73[k]
                  + f_0 * kf_73[k];

        t_74[k] = -3.0 * hf_74[k]
                  + f_0 * kf_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, hf_75, hf_76, hf_77, hf_78, hf_79, \
                         kf_75, kf_76, kf_77, kf_78, kf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -3.0 * hf_75[k]
                  + f_0 * kf_75[k];

        t_76[k] = -3.0 * hf_76[k]
                  + f_0 * kf_76[k];

        t_77[k] = -3.0 * hf_77[k]
                  + f_0 * kf_77[k];

        t_78[k] = -3.0 * hf_78[k]
                  + f_0 * kf_78[k];

        t_79[k] = -3.0 * hf_79[k]
                  + f_0 * kf_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, hf_80, hf_81, hf_82, hf_83, hf_84, \
                         kf_80, kf_81, kf_82, kf_83, kf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * hf_80[k]
                  + f_0 * kf_80[k];

        t_81[k] = -3.0 * hf_81[k]
                  + f_0 * kf_81[k];

        t_82[k] = -3.0 * hf_82[k]
                  + f_0 * kf_82[k];

        t_83[k] = -3.0 * hf_83[k]
                  + f_0 * kf_83[k];

        t_84[k] = -3.0 * hf_84[k]
                  + f_0 * kf_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hf_85, hf_86, hf_87, hf_88, hf_89, \
                         kf_85, kf_86, kf_87, kf_88, kf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -3.0 * hf_85[k]
                  + f_0 * kf_85[k];

        t_86[k] = -3.0 * hf_86[k]
                  + f_0 * kf_86[k];

        t_87[k] = -3.0 * hf_87[k]
                  + f_0 * kf_87[k];

        t_88[k] = -3.0 * hf_88[k]
                  + f_0 * kf_88[k];

        t_89[k] = -3.0 * hf_89[k]
                  + f_0 * kf_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, hf_90, hf_91, hf_92, hf_93, hf_94, \
                         kf_90, kf_91, kf_92, kf_93, kf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -3.0 * hf_90[k]
                  + f_0 * kf_90[k];

        t_91[k] = -3.0 * hf_91[k]
                  + f_0 * kf_91[k];

        t_92[k] = -3.0 * hf_92[k]
                  + f_0 * kf_92[k];

        t_93[k] = -3.0 * hf_93[k]
                  + f_0 * kf_93[k];

        t_94[k] = -3.0 * hf_94[k]
                  + f_0 * kf_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, hf_95, hf_96, hf_97, hf_98, hf_99, \
                         kf_95, kf_96, kf_97, kf_98, kf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -3.0 * hf_95[k]
                  + f_0 * kf_95[k];

        t_96[k] = -3.0 * hf_96[k]
                  + f_0 * kf_96[k];

        t_97[k] = -3.0 * hf_97[k]
                  + f_0 * kf_97[k];

        t_98[k] = -3.0 * hf_98[k]
                  + f_0 * kf_98[k];

        t_99[k] = -3.0 * hf_99[k]
                  + f_0 * kf_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, hf_100, hf_101, hf_102, hf_103, \
                         hf_104, kf_100, kf_101, kf_102, kf_103, \
                         kf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -2.0 * hf_100[k]
                   + f_0 * kf_100[k];

        t_101[k] = -2.0 * hf_101[k]
                   + f_0 * kf_101[k];

        t_102[k] = -2.0 * hf_102[k]
                   + f_0 * kf_102[k];

        t_103[k] = -2.0 * hf_103[k]
                   + f_0 * kf_103[k];

        t_104[k] = -2.0 * hf_104[k]
                   + f_0 * kf_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, hf_105, hf_106, hf_107, hf_108, \
                         hf_109, kf_105, kf_106, kf_107, kf_108, \
                         kf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -2.0 * hf_105[k]
                   + f_0 * kf_105[k];

        t_106[k] = -2.0 * hf_106[k]
                   + f_0 * kf_106[k];

        t_107[k] = -2.0 * hf_107[k]
                   + f_0 * kf_107[k];

        t_108[k] = -2.0 * hf_108[k]
                   + f_0 * kf_108[k];

        t_109[k] = -2.0 * hf_109[k]
                   + f_0 * kf_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, hf_110, hf_111, hf_112, hf_113, \
                         hf_114, kf_110, kf_111, kf_112, kf_113, \
                         kf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * hf_110[k]
                   + f_0 * kf_110[k];

        t_111[k] = -2.0 * hf_111[k]
                   + f_0 * kf_111[k];

        t_112[k] = -2.0 * hf_112[k]
                   + f_0 * kf_112[k];

        t_113[k] = -2.0 * hf_113[k]
                   + f_0 * kf_113[k];

        t_114[k] = -2.0 * hf_114[k]
                   + f_0 * kf_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, hf_115, hf_116, hf_117, hf_118, \
                         hf_119, kf_115, kf_116, kf_117, kf_118, \
                         kf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -2.0 * hf_115[k]
                   + f_0 * kf_115[k];

        t_116[k] = -2.0 * hf_116[k]
                   + f_0 * kf_116[k];

        t_117[k] = -2.0 * hf_117[k]
                   + f_0 * kf_117[k];

        t_118[k] = -2.0 * hf_118[k]
                   + f_0 * kf_118[k];

        t_119[k] = -2.0 * hf_119[k]
                   + f_0 * kf_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, hf_120, hf_121, hf_122, hf_123, \
                         hf_124, kf_120, kf_121, kf_122, kf_123, \
                         kf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -2.0 * hf_120[k]
                   + f_0 * kf_120[k];

        t_121[k] = -2.0 * hf_121[k]
                   + f_0 * kf_121[k];

        t_122[k] = -2.0 * hf_122[k]
                   + f_0 * kf_122[k];

        t_123[k] = -2.0 * hf_123[k]
                   + f_0 * kf_123[k];

        t_124[k] = -2.0 * hf_124[k]
                   + f_0 * kf_124[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, hf_125, hf_126, hf_127, hf_128, \
                         hf_129, kf_125, kf_126, kf_127, kf_128, \
                         kf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -2.0 * hf_125[k]
                   + f_0 * kf_125[k];

        t_126[k] = -2.0 * hf_126[k]
                   + f_0 * kf_126[k];

        t_127[k] = -2.0 * hf_127[k]
                   + f_0 * kf_127[k];

        t_128[k] = -2.0 * hf_128[k]
                   + f_0 * kf_128[k];

        t_129[k] = -2.0 * hf_129[k]
                   + f_0 * kf_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, hf_130, hf_131, hf_132, hf_133, \
                         hf_134, kf_130, kf_131, kf_132, kf_133, \
                         kf_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -2.0 * hf_130[k]
                   + f_0 * kf_130[k];

        t_131[k] = -2.0 * hf_131[k]
                   + f_0 * kf_131[k];

        t_132[k] = -2.0 * hf_132[k]
                   + f_0 * kf_132[k];

        t_133[k] = -2.0 * hf_133[k]
                   + f_0 * kf_133[k];

        t_134[k] = -2.0 * hf_134[k]
                   + f_0 * kf_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, hf_135, hf_136, hf_137, hf_138, \
                         hf_139, kf_135, kf_136, kf_137, kf_138, \
                         kf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = -2.0 * hf_135[k]
                   + f_0 * kf_135[k];

        t_136[k] = -2.0 * hf_136[k]
                   + f_0 * kf_136[k];

        t_137[k] = -2.0 * hf_137[k]
                   + f_0 * kf_137[k];

        t_138[k] = -2.0 * hf_138[k]
                   + f_0 * kf_138[k];

        t_139[k] = -2.0 * hf_139[k]
                   + f_0 * kf_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, hf_140, hf_141, hf_142, hf_143, \
                         hf_144, kf_140, kf_141, kf_142, kf_143, \
                         kf_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -2.0 * hf_140[k]
                   + f_0 * kf_140[k];

        t_141[k] = -2.0 * hf_141[k]
                   + f_0 * kf_141[k];

        t_142[k] = -2.0 * hf_142[k]
                   + f_0 * kf_142[k];

        t_143[k] = -2.0 * hf_143[k]
                   + f_0 * kf_143[k];

        t_144[k] = -2.0 * hf_144[k]
                   + f_0 * kf_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, hf_145, hf_146, hf_147, hf_148, \
                         hf_149, kf_145, kf_146, kf_147, kf_148, \
                         kf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = -2.0 * hf_145[k]
                   + f_0 * kf_145[k];

        t_146[k] = -2.0 * hf_146[k]
                   + f_0 * kf_146[k];

        t_147[k] = -2.0 * hf_147[k]
                   + f_0 * kf_147[k];

        t_148[k] = -2.0 * hf_148[k]
                   + f_0 * kf_148[k];

        t_149[k] = -2.0 * hf_149[k]
                   + f_0 * kf_149[k];
    }
}

static auto
compute_prim_geom_10_if_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hf, const size_t kf,
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

    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_152 = buffer.data(hf + 152);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_154 = buffer.data(hf + 154);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_160 = buffer.data(hf + 160);
    const auto *hf_161 = buffer.data(hf + 161);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_163 = buffer.data(hf + 163);
    const auto *hf_164 = buffer.data(hf + 164);
    const auto *hf_165 = buffer.data(hf + 165);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_171 = buffer.data(hf + 171);
    const auto *hf_172 = buffer.data(hf + 172);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_174 = buffer.data(hf + 174);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_181 = buffer.data(hf + 181);
    const auto *hf_182 = buffer.data(hf + 182);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_184 = buffer.data(hf + 184);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_190 = buffer.data(hf + 190);
    const auto *hf_191 = buffer.data(hf + 191);
    const auto *hf_192 = buffer.data(hf + 192);
    const auto *hf_193 = buffer.data(hf + 193);
    const auto *hf_194 = buffer.data(hf + 194);
    const auto *hf_195 = buffer.data(hf + 195);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_201 = buffer.data(hf + 201);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_204 = buffer.data(hf + 204);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_151 = buffer.data(kf + 151);
    const auto *kf_152 = buffer.data(kf + 152);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_154 = buffer.data(kf + 154);
    const auto *kf_155 = buffer.data(kf + 155);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_157 = buffer.data(kf + 157);
    const auto *kf_158 = buffer.data(kf + 158);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_161 = buffer.data(kf + 161);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_163 = buffer.data(kf + 163);
    const auto *kf_164 = buffer.data(kf + 164);
    const auto *kf_165 = buffer.data(kf + 165);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_171 = buffer.data(kf + 171);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_173 = buffer.data(kf + 173);
    const auto *kf_174 = buffer.data(kf + 174);
    const auto *kf_175 = buffer.data(kf + 175);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_181 = buffer.data(kf + 181);
    const auto *kf_182 = buffer.data(kf + 182);
    const auto *kf_183 = buffer.data(kf + 183);
    const auto *kf_184 = buffer.data(kf + 184);
    const auto *kf_185 = buffer.data(kf + 185);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_191 = buffer.data(kf + 191);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_193 = buffer.data(kf + 193);
    const auto *kf_194 = buffer.data(kf + 194);
    const auto *kf_195 = buffer.data(kf + 195);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_201 = buffer.data(kf + 201);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_203 = buffer.data(kf + 203);
    const auto *kf_204 = buffer.data(kf + 204);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_207 = buffer.data(kf + 207);
    const auto *kf_208 = buffer.data(kf + 208);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_211 = buffer.data(kf + 211);
    const auto *kf_212 = buffer.data(kf + 212);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_214 = buffer.data(kf + 214);
    const auto *kf_215 = buffer.data(kf + 215);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_217 = buffer.data(kf + 217);
    const auto *kf_218 = buffer.data(kf + 218);
    const auto *kf_219 = buffer.data(kf + 219);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_221 = buffer.data(kf + 221);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_223 = buffer.data(kf + 223);
    const auto *kf_224 = buffer.data(kf + 224);
    const auto *kf_225 = buffer.data(kf + 225);
    const auto *kf_226 = buffer.data(kf + 226);
    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_231 = buffer.data(kf + 231);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_233 = buffer.data(kf + 233);
    const auto *kf_234 = buffer.data(kf + 234);
    const auto *kf_235 = buffer.data(kf + 235);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);
    const auto *kf_238 = buffer.data(kf + 238);
    const auto *kf_239 = buffer.data(kf + 239);
    const auto *kf_240 = buffer.data(kf + 240);
    const auto *kf_241 = buffer.data(kf + 241);
    const auto *kf_242 = buffer.data(kf + 242);
    const auto *kf_243 = buffer.data(kf + 243);
    const auto *kf_244 = buffer.data(kf + 244);
    const auto *kf_245 = buffer.data(kf + 245);
    const auto *kf_246 = buffer.data(kf + 246);
    const auto *kf_247 = buffer.data(kf + 247);
    const auto *kf_248 = buffer.data(kf + 248);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_250 = buffer.data(kf + 250);
    const auto *kf_251 = buffer.data(kf + 251);
    const auto *kf_252 = buffer.data(kf + 252);
    const auto *kf_253 = buffer.data(kf + 253);
    const auto *kf_254 = buffer.data(kf + 254);
    const auto *kf_255 = buffer.data(kf + 255);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_257 = buffer.data(kf + 257);
    const auto *kf_258 = buffer.data(kf + 258);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_260 = buffer.data(kf + 260);
    const auto *kf_261 = buffer.data(kf + 261);
    const auto *kf_262 = buffer.data(kf + 262);
    const auto *kf_263 = buffer.data(kf + 263);
    const auto *kf_264 = buffer.data(kf + 264);
    const auto *kf_265 = buffer.data(kf + 265);
    const auto *kf_266 = buffer.data(kf + 266);
    const auto *kf_267 = buffer.data(kf + 267);
    const auto *kf_268 = buffer.data(kf + 268);
    const auto *kf_269 = buffer.data(kf + 269);
    const auto *kf_270 = buffer.data(kf + 270);
    const auto *kf_271 = buffer.data(kf + 271);
    const auto *kf_272 = buffer.data(kf + 272);
    const auto *kf_273 = buffer.data(kf + 273);
    const auto *kf_274 = buffer.data(kf + 274);
    const auto *kf_275 = buffer.data(kf + 275);
    const auto *kf_276 = buffer.data(kf + 276);
    const auto *kf_277 = buffer.data(kf + 277);
    const auto *kf_278 = buffer.data(kf + 278);
    const auto *kf_279 = buffer.data(kf + 279);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, hf_150, hf_151, hf_152, hf_153, \
                         hf_154, kf_150, kf_151, kf_152, kf_153, \
                         kf_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -hf_150[k]
                   + f_0 * kf_150[k];

        t_151[k] = -hf_151[k]
                   + f_0 * kf_151[k];

        t_152[k] = -hf_152[k]
                   + f_0 * kf_152[k];

        t_153[k] = -hf_153[k]
                   + f_0 * kf_153[k];

        t_154[k] = -hf_154[k]
                   + f_0 * kf_154[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, hf_155, hf_156, hf_157, hf_158, \
                         hf_159, kf_155, kf_156, kf_157, kf_158, \
                         kf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -hf_155[k]
                   + f_0 * kf_155[k];

        t_156[k] = -hf_156[k]
                   + f_0 * kf_156[k];

        t_157[k] = -hf_157[k]
                   + f_0 * kf_157[k];

        t_158[k] = -hf_158[k]
                   + f_0 * kf_158[k];

        t_159[k] = -hf_159[k]
                   + f_0 * kf_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, hf_160, hf_161, hf_162, hf_163, \
                         hf_164, kf_160, kf_161, kf_162, kf_163, \
                         kf_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -hf_160[k]
                   + f_0 * kf_160[k];

        t_161[k] = -hf_161[k]
                   + f_0 * kf_161[k];

        t_162[k] = -hf_162[k]
                   + f_0 * kf_162[k];

        t_163[k] = -hf_163[k]
                   + f_0 * kf_163[k];

        t_164[k] = -hf_164[k]
                   + f_0 * kf_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, hf_165, hf_166, hf_167, hf_168, \
                         hf_169, kf_165, kf_166, kf_167, kf_168, \
                         kf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -hf_165[k]
                   + f_0 * kf_165[k];

        t_166[k] = -hf_166[k]
                   + f_0 * kf_166[k];

        t_167[k] = -hf_167[k]
                   + f_0 * kf_167[k];

        t_168[k] = -hf_168[k]
                   + f_0 * kf_168[k];

        t_169[k] = -hf_169[k]
                   + f_0 * kf_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hf_170, hf_171, hf_172, hf_173, \
                         hf_174, kf_170, kf_171, kf_172, kf_173, \
                         kf_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -hf_170[k]
                   + f_0 * kf_170[k];

        t_171[k] = -hf_171[k]
                   + f_0 * kf_171[k];

        t_172[k] = -hf_172[k]
                   + f_0 * kf_172[k];

        t_173[k] = -hf_173[k]
                   + f_0 * kf_173[k];

        t_174[k] = -hf_174[k]
                   + f_0 * kf_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hf_175, hf_176, hf_177, hf_178, \
                         hf_179, kf_175, kf_176, kf_177, kf_178, \
                         kf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -hf_175[k]
                   + f_0 * kf_175[k];

        t_176[k] = -hf_176[k]
                   + f_0 * kf_176[k];

        t_177[k] = -hf_177[k]
                   + f_0 * kf_177[k];

        t_178[k] = -hf_178[k]
                   + f_0 * kf_178[k];

        t_179[k] = -hf_179[k]
                   + f_0 * kf_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hf_180, hf_181, hf_182, hf_183, \
                         hf_184, kf_180, kf_181, kf_182, kf_183, \
                         kf_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -hf_180[k]
                   + f_0 * kf_180[k];

        t_181[k] = -hf_181[k]
                   + f_0 * kf_181[k];

        t_182[k] = -hf_182[k]
                   + f_0 * kf_182[k];

        t_183[k] = -hf_183[k]
                   + f_0 * kf_183[k];

        t_184[k] = -hf_184[k]
                   + f_0 * kf_184[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hf_185, hf_186, hf_187, hf_188, \
                         hf_189, kf_185, kf_186, kf_187, kf_188, \
                         kf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -hf_185[k]
                   + f_0 * kf_185[k];

        t_186[k] = -hf_186[k]
                   + f_0 * kf_186[k];

        t_187[k] = -hf_187[k]
                   + f_0 * kf_187[k];

        t_188[k] = -hf_188[k]
                   + f_0 * kf_188[k];

        t_189[k] = -hf_189[k]
                   + f_0 * kf_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hf_190, hf_191, hf_192, hf_193, \
                         hf_194, kf_190, kf_191, kf_192, kf_193, \
                         kf_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -hf_190[k]
                   + f_0 * kf_190[k];

        t_191[k] = -hf_191[k]
                   + f_0 * kf_191[k];

        t_192[k] = -hf_192[k]
                   + f_0 * kf_192[k];

        t_193[k] = -hf_193[k]
                   + f_0 * kf_193[k];

        t_194[k] = -hf_194[k]
                   + f_0 * kf_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hf_195, hf_196, hf_197, hf_198, \
                         hf_199, kf_195, kf_196, kf_197, kf_198, \
                         kf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -hf_195[k]
                   + f_0 * kf_195[k];

        t_196[k] = -hf_196[k]
                   + f_0 * kf_196[k];

        t_197[k] = -hf_197[k]
                   + f_0 * kf_197[k];

        t_198[k] = -hf_198[k]
                   + f_0 * kf_198[k];

        t_199[k] = -hf_199[k]
                   + f_0 * kf_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hf_200, hf_201, hf_202, hf_203, \
                         hf_204, kf_200, kf_201, kf_202, kf_203, \
                         kf_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -hf_200[k]
                   + f_0 * kf_200[k];

        t_201[k] = -hf_201[k]
                   + f_0 * kf_201[k];

        t_202[k] = -hf_202[k]
                   + f_0 * kf_202[k];

        t_203[k] = -hf_203[k]
                   + f_0 * kf_203[k];

        t_204[k] = -hf_204[k]
                   + f_0 * kf_204[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hf_205, hf_206, hf_207, hf_208, \
                         hf_209, kf_205, kf_206, kf_207, kf_208, \
                         kf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -hf_205[k]
                   + f_0 * kf_205[k];

        t_206[k] = -hf_206[k]
                   + f_0 * kf_206[k];

        t_207[k] = -hf_207[k]
                   + f_0 * kf_207[k];

        t_208[k] = -hf_208[k]
                   + f_0 * kf_208[k];

        t_209[k] = -hf_209[k]
                   + f_0 * kf_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, kf_210, \
                         kf_211, kf_212, kf_213, kf_214, kf_215, kf_216, \
                         kf_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * kf_210[k];

        t_211[k] = f_0 * kf_211[k];

        t_212[k] = f_0 * kf_212[k];

        t_213[k] = f_0 * kf_213[k];

        t_214[k] = f_0 * kf_214[k];

        t_215[k] = f_0 * kf_215[k];

        t_216[k] = f_0 * kf_216[k];

        t_217[k] = f_0 * kf_217[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, t_224, t_225, kf_218, \
                         kf_219, kf_220, kf_221, kf_222, kf_223, kf_224, \
                         kf_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * kf_218[k];

        t_219[k] = f_0 * kf_219[k];

        t_220[k] = f_0 * kf_220[k];

        t_221[k] = f_0 * kf_221[k];

        t_222[k] = f_0 * kf_222[k];

        t_223[k] = f_0 * kf_223[k];

        t_224[k] = f_0 * kf_224[k];

        t_225[k] = f_0 * kf_225[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, t_230, t_231, t_232, t_233, kf_226, \
                         kf_227, kf_228, kf_229, kf_230, kf_231, kf_232, \
                         kf_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_0 * kf_226[k];

        t_227[k] = f_0 * kf_227[k];

        t_228[k] = f_0 * kf_228[k];

        t_229[k] = f_0 * kf_229[k];

        t_230[k] = f_0 * kf_230[k];

        t_231[k] = f_0 * kf_231[k];

        t_232[k] = f_0 * kf_232[k];

        t_233[k] = f_0 * kf_233[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, t_239, t_240, t_241, kf_234, \
                         kf_235, kf_236, kf_237, kf_238, kf_239, kf_240, \
                         kf_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_0 * kf_234[k];

        t_235[k] = f_0 * kf_235[k];

        t_236[k] = f_0 * kf_236[k];

        t_237[k] = f_0 * kf_237[k];

        t_238[k] = f_0 * kf_238[k];

        t_239[k] = f_0 * kf_239[k];

        t_240[k] = f_0 * kf_240[k];

        t_241[k] = f_0 * kf_241[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, t_247, t_248, t_249, kf_242, \
                         kf_243, kf_244, kf_245, kf_246, kf_247, kf_248, \
                         kf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_0 * kf_242[k];

        t_243[k] = f_0 * kf_243[k];

        t_244[k] = f_0 * kf_244[k];

        t_245[k] = f_0 * kf_245[k];

        t_246[k] = f_0 * kf_246[k];

        t_247[k] = f_0 * kf_247[k];

        t_248[k] = f_0 * kf_248[k];

        t_249[k] = f_0 * kf_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, t_255, t_256, t_257, kf_250, \
                         kf_251, kf_252, kf_253, kf_254, kf_255, kf_256, \
                         kf_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_0 * kf_250[k];

        t_251[k] = f_0 * kf_251[k];

        t_252[k] = f_0 * kf_252[k];

        t_253[k] = f_0 * kf_253[k];

        t_254[k] = f_0 * kf_254[k];

        t_255[k] = f_0 * kf_255[k];

        t_256[k] = f_0 * kf_256[k];

        t_257[k] = f_0 * kf_257[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, t_263, t_264, t_265, kf_258, \
                         kf_259, kf_260, kf_261, kf_262, kf_263, kf_264, \
                         kf_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_0 * kf_258[k];

        t_259[k] = f_0 * kf_259[k];

        t_260[k] = f_0 * kf_260[k];

        t_261[k] = f_0 * kf_261[k];

        t_262[k] = f_0 * kf_262[k];

        t_263[k] = f_0 * kf_263[k];

        t_264[k] = f_0 * kf_264[k];

        t_265[k] = f_0 * kf_265[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, t_271, t_272, t_273, kf_266, \
                         kf_267, kf_268, kf_269, kf_270, kf_271, kf_272, \
                         kf_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_0 * kf_266[k];

        t_267[k] = f_0 * kf_267[k];

        t_268[k] = f_0 * kf_268[k];

        t_269[k] = f_0 * kf_269[k];

        t_270[k] = f_0 * kf_270[k];

        t_271[k] = f_0 * kf_271[k];

        t_272[k] = f_0 * kf_272[k];

        t_273[k] = f_0 * kf_273[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, t_279, kf_274, kf_275, kf_276, \
                         kf_277, kf_278, kf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_0 * kf_274[k];

        t_275[k] = f_0 * kf_275[k];

        t_276[k] = f_0 * kf_276[k];

        t_277[k] = f_0 * kf_277[k];

        t_278[k] = f_0 * kf_278[k];

        t_279[k] = f_0 * kf_279[k];
    }
}

auto
compute_prim_geom_10_if_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t hf, const size_t kf,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_if_electron_repulsion_0_piece0(buffer, target, hf, kf, ncols, alpha);

    compute_prim_geom_10_if_electron_repulsion_0_piece1(buffer, target, hf, kf, ncols, alpha);
}

static auto
compute_prim_geom_10_if_electron_repulsion_1_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hf, const size_t kf,
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

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_102 = buffer.data(hf + 102);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_104 = buffer.data(hf + 104);
    const auto *hf_105 = buffer.data(hf + 105);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_107 = buffer.data(hf + 107);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_111 = buffer.data(hf + 111);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_113 = buffer.data(hf + 113);
    const auto *hf_114 = buffer.data(hf + 114);
    const auto *hf_115 = buffer.data(hf + 115);
    const auto *hf_116 = buffer.data(hf + 116);

    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_12 = buffer.data(kf + 12);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_15 = buffer.data(kf + 15);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_65 = buffer.data(kf + 65);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_104 = buffer.data(kf + 104);
    const auto *kf_105 = buffer.data(kf + 105);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_113 = buffer.data(kf + 113);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_115 = buffer.data(kf + 115);
    const auto *kf_116 = buffer.data(kf + 116);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_121 = buffer.data(kf + 121);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_123 = buffer.data(kf + 123);
    const auto *kf_124 = buffer.data(kf + 124);
    const auto *kf_125 = buffer.data(kf + 125);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_131 = buffer.data(kf + 131);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_135 = buffer.data(kf + 135);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_151 = buffer.data(kf + 151);
    const auto *kf_152 = buffer.data(kf + 152);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_154 = buffer.data(kf + 154);
    const auto *kf_155 = buffer.data(kf + 155);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_157 = buffer.data(kf + 157);
    const auto *kf_158 = buffer.data(kf + 158);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_161 = buffer.data(kf + 161);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_163 = buffer.data(kf + 163);
    const auto *kf_164 = buffer.data(kf + 164);
    const auto *kf_165 = buffer.data(kf + 165);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_171 = buffer.data(kf + 171);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_173 = buffer.data(kf + 173);
    const auto *kf_174 = buffer.data(kf + 174);
    const auto *kf_175 = buffer.data(kf + 175);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_181 = buffer.data(kf + 181);
    const auto *kf_182 = buffer.data(kf + 182);
    const auto *kf_183 = buffer.data(kf + 183);
    const auto *kf_184 = buffer.data(kf + 184);
    const auto *kf_185 = buffer.data(kf + 185);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_191 = buffer.data(kf + 191);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_193 = buffer.data(kf + 193);
    const auto *kf_194 = buffer.data(kf + 194);
    const auto *kf_195 = buffer.data(kf + 195);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_211 = buffer.data(kf + 211);
    const auto *kf_212 = buffer.data(kf + 212);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_214 = buffer.data(kf + 214);
    const auto *kf_215 = buffer.data(kf + 215);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_217 = buffer.data(kf + 217);
    const auto *kf_218 = buffer.data(kf + 218);
    const auto *kf_219 = buffer.data(kf + 219);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_221 = buffer.data(kf + 221);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_223 = buffer.data(kf + 223);
    const auto *kf_224 = buffer.data(kf + 224);
    const auto *kf_225 = buffer.data(kf + 225);
    const auto *kf_226 = buffer.data(kf + 226);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kf_10, kf_11, kf_12, kf_13, \
                         kf_14, kf_15, kf_16, kf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kf_10[k];

        t_1[k] = f_0 * kf_11[k];

        t_2[k] = f_0 * kf_12[k];

        t_3[k] = f_0 * kf_13[k];

        t_4[k] = f_0 * kf_14[k];

        t_5[k] = f_0 * kf_15[k];

        t_6[k] = f_0 * kf_16[k];

        t_7[k] = f_0 * kf_17[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, hf_0, hf_1, hf_2, hf_3, kf_18, \
                         kf_19, kf_30, kf_31, kf_32, kf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kf_18[k];

        t_9[k] = f_0 * kf_19[k];

        t_10[k] = -hf_0[k]
                  + f_0 * kf_30[k];

        t_11[k] = -hf_1[k]
                  + f_0 * kf_31[k];

        t_12[k] = -hf_2[k]
                  + f_0 * kf_32[k];

        t_13[k] = -hf_3[k]
                  + f_0 * kf_33[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, hf_4, hf_5, hf_6, hf_7, hf_8, kf_34, \
                         kf_35, kf_36, kf_37, kf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -hf_4[k]
                  + f_0 * kf_34[k];

        t_15[k] = -hf_5[k]
                  + f_0 * kf_35[k];

        t_16[k] = -hf_6[k]
                  + f_0 * kf_36[k];

        t_17[k] = -hf_7[k]
                  + f_0 * kf_37[k];

        t_18[k] = -hf_8[k]
                  + f_0 * kf_38[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, hf_9, kf_39, kf_40, kf_41, \
                         kf_42, kf_43, kf_44, kf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -hf_9[k]
                  + f_0 * kf_39[k];

        t_20[k] = f_0 * kf_40[k];

        t_21[k] = f_0 * kf_41[k];

        t_22[k] = f_0 * kf_42[k];

        t_23[k] = f_0 * kf_43[k];

        t_24[k] = f_0 * kf_44[k];

        t_25[k] = f_0 * kf_45[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, hf_10, hf_11, kf_46, kf_47, \
                         kf_48, kf_49, kf_60, kf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * kf_46[k];

        t_27[k] = f_0 * kf_47[k];

        t_28[k] = f_0 * kf_48[k];

        t_29[k] = f_0 * kf_49[k];

        t_30[k] = -2.0 * hf_10[k]
                  + f_0 * kf_60[k];

        t_31[k] = -2.0 * hf_11[k]
                  + f_0 * kf_61[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, hf_12, hf_13, hf_14, hf_15, hf_16, \
                         kf_62, kf_63, kf_64, kf_65, kf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -2.0 * hf_12[k]
                  + f_0 * kf_62[k];

        t_33[k] = -2.0 * hf_13[k]
                  + f_0 * kf_63[k];

        t_34[k] = -2.0 * hf_14[k]
                  + f_0 * kf_64[k];

        t_35[k] = -2.0 * hf_15[k]
                  + f_0 * kf_65[k];

        t_36[k] = -2.0 * hf_16[k]
                  + f_0 * kf_66[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, hf_17, hf_18, hf_19, hf_20, hf_21, \
                         kf_67, kf_68, kf_69, kf_70, kf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -2.0 * hf_17[k]
                  + f_0 * kf_67[k];

        t_38[k] = -2.0 * hf_18[k]
                  + f_0 * kf_68[k];

        t_39[k] = -2.0 * hf_19[k]
                  + f_0 * kf_69[k];

        t_40[k] = -hf_20[k]
                  + f_0 * kf_70[k];

        t_41[k] = -hf_21[k]
                  + f_0 * kf_71[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, hf_22, hf_23, hf_24, hf_25, hf_26, \
                         kf_72, kf_73, kf_74, kf_75, kf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -hf_22[k]
                  + f_0 * kf_72[k];

        t_43[k] = -hf_23[k]
                  + f_0 * kf_73[k];

        t_44[k] = -hf_24[k]
                  + f_0 * kf_74[k];

        t_45[k] = -hf_25[k]
                  + f_0 * kf_75[k];

        t_46[k] = -hf_26[k]
                  + f_0 * kf_76[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, hf_27, hf_28, hf_29, kf_77, \
                         kf_78, kf_79, kf_80, kf_81, kf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -hf_27[k]
                  + f_0 * kf_77[k];

        t_48[k] = -hf_28[k]
                  + f_0 * kf_78[k];

        t_49[k] = -hf_29[k]
                  + f_0 * kf_79[k];

        t_50[k] = f_0 * kf_80[k];

        t_51[k] = f_0 * kf_81[k];

        t_52[k] = f_0 * kf_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, kf_83, kf_84, kf_85, kf_86, \
                         kf_87, kf_88, kf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * kf_83[k];

        t_54[k] = f_0 * kf_84[k];

        t_55[k] = f_0 * kf_85[k];

        t_56[k] = f_0 * kf_86[k];

        t_57[k] = f_0 * kf_87[k];

        t_58[k] = f_0 * kf_88[k];

        t_59[k] = f_0 * kf_89[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, hf_30, hf_31, hf_32, hf_33, hf_34, \
                         kf_100, kf_101, kf_102, kf_103, kf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * hf_30[k]
                  + f_0 * kf_100[k];

        t_61[k] = -3.0 * hf_31[k]
                  + f_0 * kf_101[k];

        t_62[k] = -3.0 * hf_32[k]
                  + f_0 * kf_102[k];

        t_63[k] = -3.0 * hf_33[k]
                  + f_0 * kf_103[k];

        t_64[k] = -3.0 * hf_34[k]
                  + f_0 * kf_104[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, hf_35, hf_36, hf_37, hf_38, hf_39, \
                         kf_105, kf_106, kf_107, kf_108, kf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * hf_35[k]
                  + f_0 * kf_105[k];

        t_66[k] = -3.0 * hf_36[k]
                  + f_0 * kf_106[k];

        t_67[k] = -3.0 * hf_37[k]
                  + f_0 * kf_107[k];

        t_68[k] = -3.0 * hf_38[k]
                  + f_0 * kf_108[k];

        t_69[k] = -3.0 * hf_39[k]
                  + f_0 * kf_109[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, hf_40, hf_41, hf_42, hf_43, hf_44, \
                         kf_110, kf_111, kf_112, kf_113, kf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * hf_40[k]
                  + f_0 * kf_110[k];

        t_71[k] = -2.0 * hf_41[k]
                  + f_0 * kf_111[k];

        t_72[k] = -2.0 * hf_42[k]
                  + f_0 * kf_112[k];

        t_73[k] = -2.0 * hf_43[k]
                  + f_0 * kf_113[k];

        t_74[k] = -2.0 * hf_44[k]
                  + f_0 * kf_114[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, hf_45, hf_46, hf_47, hf_48, hf_49, \
                         kf_115, kf_116, kf_117, kf_118, kf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * hf_45[k]
                  + f_0 * kf_115[k];

        t_76[k] = -2.0 * hf_46[k]
                  + f_0 * kf_116[k];

        t_77[k] = -2.0 * hf_47[k]
                  + f_0 * kf_117[k];

        t_78[k] = -2.0 * hf_48[k]
                  + f_0 * kf_118[k];

        t_79[k] = -2.0 * hf_49[k]
                  + f_0 * kf_119[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, hf_50, hf_51, hf_52, hf_53, hf_54, \
                         kf_120, kf_121, kf_122, kf_123, kf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -hf_50[k]
                  + f_0 * kf_120[k];

        t_81[k] = -hf_51[k]
                  + f_0 * kf_121[k];

        t_82[k] = -hf_52[k]
                  + f_0 * kf_122[k];

        t_83[k] = -hf_53[k]
                  + f_0 * kf_123[k];

        t_84[k] = -hf_54[k]
                  + f_0 * kf_124[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, hf_55, hf_56, hf_57, hf_58, hf_59, \
                         kf_125, kf_126, kf_127, kf_128, kf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -hf_55[k]
                  + f_0 * kf_125[k];

        t_86[k] = -hf_56[k]
                  + f_0 * kf_126[k];

        t_87[k] = -hf_57[k]
                  + f_0 * kf_127[k];

        t_88[k] = -hf_58[k]
                  + f_0 * kf_128[k];

        t_89[k] = -hf_59[k]
                  + f_0 * kf_129[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, kf_130, kf_131, \
                         kf_132, kf_133, kf_134, kf_135, kf_136, \
                         kf_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * kf_130[k];

        t_91[k] = f_0 * kf_131[k];

        t_92[k] = f_0 * kf_132[k];

        t_93[k] = f_0 * kf_133[k];

        t_94[k] = f_0 * kf_134[k];

        t_95[k] = f_0 * kf_135[k];

        t_96[k] = f_0 * kf_136[k];

        t_97[k] = f_0 * kf_137[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, hf_60, hf_61, hf_62, hf_63, \
                         kf_138, kf_139, kf_150, kf_151, kf_152, \
                         kf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * kf_138[k];

        t_99[k] = f_0 * kf_139[k];

        t_100[k] = -4.0 * hf_60[k]
                   + f_0 * kf_150[k];

        t_101[k] = -4.0 * hf_61[k]
                   + f_0 * kf_151[k];

        t_102[k] = -4.0 * hf_62[k]
                   + f_0 * kf_152[k];

        t_103[k] = -4.0 * hf_63[k]
                   + f_0 * kf_153[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, hf_64, hf_65, hf_66, hf_67, hf_68, \
                         kf_154, kf_155, kf_156, kf_157, kf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -4.0 * hf_64[k]
                   + f_0 * kf_154[k];

        t_105[k] = -4.0 * hf_65[k]
                   + f_0 * kf_155[k];

        t_106[k] = -4.0 * hf_66[k]
                   + f_0 * kf_156[k];

        t_107[k] = -4.0 * hf_67[k]
                   + f_0 * kf_157[k];

        t_108[k] = -4.0 * hf_68[k]
                   + f_0 * kf_158[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, hf_69, hf_70, hf_71, hf_72, hf_73, \
                         kf_159, kf_160, kf_161, kf_162, kf_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -4.0 * hf_69[k]
                   + f_0 * kf_159[k];

        t_110[k] = -3.0 * hf_70[k]
                   + f_0 * kf_160[k];

        t_111[k] = -3.0 * hf_71[k]
                   + f_0 * kf_161[k];

        t_112[k] = -3.0 * hf_72[k]
                   + f_0 * kf_162[k];

        t_113[k] = -3.0 * hf_73[k]
                   + f_0 * kf_163[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, hf_74, hf_75, hf_76, hf_77, hf_78, \
                         kf_164, kf_165, kf_166, kf_167, kf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -3.0 * hf_74[k]
                   + f_0 * kf_164[k];

        t_115[k] = -3.0 * hf_75[k]
                   + f_0 * kf_165[k];

        t_116[k] = -3.0 * hf_76[k]
                   + f_0 * kf_166[k];

        t_117[k] = -3.0 * hf_77[k]
                   + f_0 * kf_167[k];

        t_118[k] = -3.0 * hf_78[k]
                   + f_0 * kf_168[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, hf_79, hf_80, hf_81, hf_82, hf_83, \
                         kf_169, kf_170, kf_171, kf_172, kf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -3.0 * hf_79[k]
                   + f_0 * kf_169[k];

        t_120[k] = -2.0 * hf_80[k]
                   + f_0 * kf_170[k];

        t_121[k] = -2.0 * hf_81[k]
                   + f_0 * kf_171[k];

        t_122[k] = -2.0 * hf_82[k]
                   + f_0 * kf_172[k];

        t_123[k] = -2.0 * hf_83[k]
                   + f_0 * kf_173[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, hf_84, hf_85, hf_86, hf_87, hf_88, \
                         kf_174, kf_175, kf_176, kf_177, kf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * hf_84[k]
                   + f_0 * kf_174[k];

        t_125[k] = -2.0 * hf_85[k]
                   + f_0 * kf_175[k];

        t_126[k] = -2.0 * hf_86[k]
                   + f_0 * kf_176[k];

        t_127[k] = -2.0 * hf_87[k]
                   + f_0 * kf_177[k];

        t_128[k] = -2.0 * hf_88[k]
                   + f_0 * kf_178[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, hf_89, hf_90, hf_91, hf_92, hf_93, \
                         kf_179, kf_180, kf_181, kf_182, kf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -2.0 * hf_89[k]
                   + f_0 * kf_179[k];

        t_130[k] = -hf_90[k]
                   + f_0 * kf_180[k];

        t_131[k] = -hf_91[k]
                   + f_0 * kf_181[k];

        t_132[k] = -hf_92[k]
                   + f_0 * kf_182[k];

        t_133[k] = -hf_93[k]
                   + f_0 * kf_183[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, hf_94, hf_95, hf_96, hf_97, hf_98, \
                         kf_184, kf_185, kf_186, kf_187, kf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -hf_94[k]
                   + f_0 * kf_184[k];

        t_135[k] = -hf_95[k]
                   + f_0 * kf_185[k];

        t_136[k] = -hf_96[k]
                   + f_0 * kf_186[k];

        t_137[k] = -hf_97[k]
                   + f_0 * kf_187[k];

        t_138[k] = -hf_98[k]
                   + f_0 * kf_188[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, hf_99, kf_189, \
                         kf_190, kf_191, kf_192, kf_193, kf_194, \
                         kf_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -hf_99[k]
                   + f_0 * kf_189[k];

        t_140[k] = f_0 * kf_190[k];

        t_141[k] = f_0 * kf_191[k];

        t_142[k] = f_0 * kf_192[k];

        t_143[k] = f_0 * kf_193[k];

        t_144[k] = f_0 * kf_194[k];

        t_145[k] = f_0 * kf_195[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, t_151, hf_100, hf_101, kf_196, \
                         kf_197, kf_198, kf_199, kf_210, kf_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_0 * kf_196[k];

        t_147[k] = f_0 * kf_197[k];

        t_148[k] = f_0 * kf_198[k];

        t_149[k] = f_0 * kf_199[k];

        t_150[k] = -5.0 * hf_100[k]
                   + f_0 * kf_210[k];

        t_151[k] = -5.0 * hf_101[k]
                   + f_0 * kf_211[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, hf_102, hf_103, hf_104, hf_105, \
                         hf_106, kf_212, kf_213, kf_214, kf_215, \
                         kf_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = -5.0 * hf_102[k]
                   + f_0 * kf_212[k];

        t_153[k] = -5.0 * hf_103[k]
                   + f_0 * kf_213[k];

        t_154[k] = -5.0 * hf_104[k]
                   + f_0 * kf_214[k];

        t_155[k] = -5.0 * hf_105[k]
                   + f_0 * kf_215[k];

        t_156[k] = -5.0 * hf_106[k]
                   + f_0 * kf_216[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, hf_107, hf_108, hf_109, hf_110, \
                         hf_111, kf_217, kf_218, kf_219, kf_220, \
                         kf_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = -5.0 * hf_107[k]
                   + f_0 * kf_217[k];

        t_158[k] = -5.0 * hf_108[k]
                   + f_0 * kf_218[k];

        t_159[k] = -5.0 * hf_109[k]
                   + f_0 * kf_219[k];

        t_160[k] = -4.0 * hf_110[k]
                   + f_0 * kf_220[k];

        t_161[k] = -4.0 * hf_111[k]
                   + f_0 * kf_221[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, hf_112, hf_113, hf_114, hf_115, \
                         hf_116, kf_222, kf_223, kf_224, kf_225, \
                         kf_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -4.0 * hf_112[k]
                   + f_0 * kf_222[k];

        t_163[k] = -4.0 * hf_113[k]
                   + f_0 * kf_223[k];

        t_164[k] = -4.0 * hf_114[k]
                   + f_0 * kf_224[k];

        t_165[k] = -4.0 * hf_115[k]
                   + f_0 * kf_225[k];

        t_166[k] = -4.0 * hf_116[k]
                   + f_0 * kf_226[k];
    }
}

static auto
compute_prim_geom_10_if_electron_repulsion_1_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hf, const size_t kf,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_121 = buffer.data(hf + 121);
    const auto *hf_122 = buffer.data(hf + 122);
    const auto *hf_123 = buffer.data(hf + 123);
    const auto *hf_124 = buffer.data(hf + 124);
    const auto *hf_125 = buffer.data(hf + 125);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_127 = buffer.data(hf + 127);
    const auto *hf_128 = buffer.data(hf + 128);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_130 = buffer.data(hf + 130);
    const auto *hf_131 = buffer.data(hf + 131);
    const auto *hf_132 = buffer.data(hf + 132);
    const auto *hf_133 = buffer.data(hf + 133);
    const auto *hf_134 = buffer.data(hf + 134);
    const auto *hf_135 = buffer.data(hf + 135);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_137 = buffer.data(hf + 137);
    const auto *hf_138 = buffer.data(hf + 138);
    const auto *hf_139 = buffer.data(hf + 139);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_141 = buffer.data(hf + 141);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_143 = buffer.data(hf + 143);
    const auto *hf_144 = buffer.data(hf + 144);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_146 = buffer.data(hf + 146);
    const auto *hf_147 = buffer.data(hf + 147);
    const auto *hf_148 = buffer.data(hf + 148);
    const auto *hf_149 = buffer.data(hf + 149);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_152 = buffer.data(hf + 152);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_154 = buffer.data(hf + 154);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_160 = buffer.data(hf + 160);
    const auto *hf_161 = buffer.data(hf + 161);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_163 = buffer.data(hf + 163);
    const auto *hf_164 = buffer.data(hf + 164);
    const auto *hf_165 = buffer.data(hf + 165);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_171 = buffer.data(hf + 171);
    const auto *hf_172 = buffer.data(hf + 172);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_174 = buffer.data(hf + 174);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_181 = buffer.data(hf + 181);
    const auto *hf_182 = buffer.data(hf + 182);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_184 = buffer.data(hf + 184);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_190 = buffer.data(hf + 190);
    const auto *hf_191 = buffer.data(hf + 191);
    const auto *hf_192 = buffer.data(hf + 192);
    const auto *hf_193 = buffer.data(hf + 193);
    const auto *hf_194 = buffer.data(hf + 194);
    const auto *hf_195 = buffer.data(hf + 195);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_201 = buffer.data(hf + 201);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_204 = buffer.data(hf + 204);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_231 = buffer.data(kf + 231);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_233 = buffer.data(kf + 233);
    const auto *kf_234 = buffer.data(kf + 234);
    const auto *kf_235 = buffer.data(kf + 235);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);
    const auto *kf_238 = buffer.data(kf + 238);
    const auto *kf_239 = buffer.data(kf + 239);
    const auto *kf_240 = buffer.data(kf + 240);
    const auto *kf_241 = buffer.data(kf + 241);
    const auto *kf_242 = buffer.data(kf + 242);
    const auto *kf_243 = buffer.data(kf + 243);
    const auto *kf_244 = buffer.data(kf + 244);
    const auto *kf_245 = buffer.data(kf + 245);
    const auto *kf_246 = buffer.data(kf + 246);
    const auto *kf_247 = buffer.data(kf + 247);
    const auto *kf_248 = buffer.data(kf + 248);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_250 = buffer.data(kf + 250);
    const auto *kf_251 = buffer.data(kf + 251);
    const auto *kf_252 = buffer.data(kf + 252);
    const auto *kf_253 = buffer.data(kf + 253);
    const auto *kf_254 = buffer.data(kf + 254);
    const auto *kf_255 = buffer.data(kf + 255);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_257 = buffer.data(kf + 257);
    const auto *kf_258 = buffer.data(kf + 258);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_260 = buffer.data(kf + 260);
    const auto *kf_261 = buffer.data(kf + 261);
    const auto *kf_262 = buffer.data(kf + 262);
    const auto *kf_263 = buffer.data(kf + 263);
    const auto *kf_264 = buffer.data(kf + 264);
    const auto *kf_265 = buffer.data(kf + 265);
    const auto *kf_266 = buffer.data(kf + 266);
    const auto *kf_267 = buffer.data(kf + 267);
    const auto *kf_268 = buffer.data(kf + 268);
    const auto *kf_269 = buffer.data(kf + 269);
    const auto *kf_280 = buffer.data(kf + 280);
    const auto *kf_281 = buffer.data(kf + 281);
    const auto *kf_282 = buffer.data(kf + 282);
    const auto *kf_283 = buffer.data(kf + 283);
    const auto *kf_284 = buffer.data(kf + 284);
    const auto *kf_285 = buffer.data(kf + 285);
    const auto *kf_286 = buffer.data(kf + 286);
    const auto *kf_287 = buffer.data(kf + 287);
    const auto *kf_288 = buffer.data(kf + 288);
    const auto *kf_289 = buffer.data(kf + 289);
    const auto *kf_290 = buffer.data(kf + 290);
    const auto *kf_291 = buffer.data(kf + 291);
    const auto *kf_292 = buffer.data(kf + 292);
    const auto *kf_293 = buffer.data(kf + 293);
    const auto *kf_294 = buffer.data(kf + 294);
    const auto *kf_295 = buffer.data(kf + 295);
    const auto *kf_296 = buffer.data(kf + 296);
    const auto *kf_297 = buffer.data(kf + 297);
    const auto *kf_298 = buffer.data(kf + 298);
    const auto *kf_299 = buffer.data(kf + 299);
    const auto *kf_300 = buffer.data(kf + 300);
    const auto *kf_301 = buffer.data(kf + 301);
    const auto *kf_302 = buffer.data(kf + 302);
    const auto *kf_303 = buffer.data(kf + 303);
    const auto *kf_304 = buffer.data(kf + 304);
    const auto *kf_305 = buffer.data(kf + 305);
    const auto *kf_306 = buffer.data(kf + 306);
    const auto *kf_307 = buffer.data(kf + 307);
    const auto *kf_308 = buffer.data(kf + 308);
    const auto *kf_309 = buffer.data(kf + 309);
    const auto *kf_310 = buffer.data(kf + 310);
    const auto *kf_311 = buffer.data(kf + 311);
    const auto *kf_312 = buffer.data(kf + 312);
    const auto *kf_313 = buffer.data(kf + 313);
    const auto *kf_314 = buffer.data(kf + 314);
    const auto *kf_315 = buffer.data(kf + 315);
    const auto *kf_316 = buffer.data(kf + 316);
    const auto *kf_317 = buffer.data(kf + 317);
    const auto *kf_318 = buffer.data(kf + 318);
    const auto *kf_319 = buffer.data(kf + 319);
    const auto *kf_320 = buffer.data(kf + 320);
    const auto *kf_321 = buffer.data(kf + 321);
    const auto *kf_322 = buffer.data(kf + 322);
    const auto *kf_323 = buffer.data(kf + 323);
    const auto *kf_324 = buffer.data(kf + 324);
    const auto *kf_325 = buffer.data(kf + 325);
    const auto *kf_326 = buffer.data(kf + 326);
    const auto *kf_327 = buffer.data(kf + 327);
    const auto *kf_328 = buffer.data(kf + 328);
    const auto *kf_329 = buffer.data(kf + 329);
    const auto *kf_330 = buffer.data(kf + 330);
    const auto *kf_331 = buffer.data(kf + 331);
    const auto *kf_332 = buffer.data(kf + 332);
    const auto *kf_333 = buffer.data(kf + 333);
    const auto *kf_334 = buffer.data(kf + 334);
    const auto *kf_335 = buffer.data(kf + 335);
    const auto *kf_336 = buffer.data(kf + 336);
    const auto *kf_337 = buffer.data(kf + 337);
    const auto *kf_338 = buffer.data(kf + 338);
    const auto *kf_339 = buffer.data(kf + 339);
    const auto *kf_340 = buffer.data(kf + 340);
    const auto *kf_341 = buffer.data(kf + 341);
    const auto *kf_342 = buffer.data(kf + 342);
    const auto *kf_343 = buffer.data(kf + 343);
    const auto *kf_344 = buffer.data(kf + 344);
    const auto *kf_345 = buffer.data(kf + 345);
    const auto *kf_346 = buffer.data(kf + 346);
    const auto *kf_347 = buffer.data(kf + 347);
    const auto *kf_348 = buffer.data(kf + 348);
    const auto *kf_349 = buffer.data(kf + 349);

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, hf_117, hf_118, hf_119, hf_120, \
                         hf_121, kf_227, kf_228, kf_229, kf_230, \
                         kf_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -4.0 * hf_117[k]
                   + f_0 * kf_227[k];

        t_168[k] = -4.0 * hf_118[k]
                   + f_0 * kf_228[k];

        t_169[k] = -4.0 * hf_119[k]
                   + f_0 * kf_229[k];

        t_170[k] = -3.0 * hf_120[k]
                   + f_0 * kf_230[k];

        t_171[k] = -3.0 * hf_121[k]
                   + f_0 * kf_231[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, hf_122, hf_123, hf_124, hf_125, \
                         hf_126, kf_232, kf_233, kf_234, kf_235, \
                         kf_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -3.0 * hf_122[k]
                   + f_0 * kf_232[k];

        t_173[k] = -3.0 * hf_123[k]
                   + f_0 * kf_233[k];

        t_174[k] = -3.0 * hf_124[k]
                   + f_0 * kf_234[k];

        t_175[k] = -3.0 * hf_125[k]
                   + f_0 * kf_235[k];

        t_176[k] = -3.0 * hf_126[k]
                   + f_0 * kf_236[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, hf_127, hf_128, hf_129, hf_130, \
                         hf_131, kf_237, kf_238, kf_239, kf_240, \
                         kf_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -3.0 * hf_127[k]
                   + f_0 * kf_237[k];

        t_178[k] = -3.0 * hf_128[k]
                   + f_0 * kf_238[k];

        t_179[k] = -3.0 * hf_129[k]
                   + f_0 * kf_239[k];

        t_180[k] = -2.0 * hf_130[k]
                   + f_0 * kf_240[k];

        t_181[k] = -2.0 * hf_131[k]
                   + f_0 * kf_241[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, hf_132, hf_133, hf_134, hf_135, \
                         hf_136, kf_242, kf_243, kf_244, kf_245, \
                         kf_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = -2.0 * hf_132[k]
                   + f_0 * kf_242[k];

        t_183[k] = -2.0 * hf_133[k]
                   + f_0 * kf_243[k];

        t_184[k] = -2.0 * hf_134[k]
                   + f_0 * kf_244[k];

        t_185[k] = -2.0 * hf_135[k]
                   + f_0 * kf_245[k];

        t_186[k] = -2.0 * hf_136[k]
                   + f_0 * kf_246[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, hf_137, hf_138, hf_139, hf_140, \
                         hf_141, kf_247, kf_248, kf_249, kf_250, \
                         kf_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = -2.0 * hf_137[k]
                   + f_0 * kf_247[k];

        t_188[k] = -2.0 * hf_138[k]
                   + f_0 * kf_248[k];

        t_189[k] = -2.0 * hf_139[k]
                   + f_0 * kf_249[k];

        t_190[k] = -hf_140[k]
                   + f_0 * kf_250[k];

        t_191[k] = -hf_141[k]
                   + f_0 * kf_251[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, hf_142, hf_143, hf_144, hf_145, \
                         hf_146, kf_252, kf_253, kf_254, kf_255, \
                         kf_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -hf_142[k]
                   + f_0 * kf_252[k];

        t_193[k] = -hf_143[k]
                   + f_0 * kf_253[k];

        t_194[k] = -hf_144[k]
                   + f_0 * kf_254[k];

        t_195[k] = -hf_145[k]
                   + f_0 * kf_255[k];

        t_196[k] = -hf_146[k]
                   + f_0 * kf_256[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, t_202, hf_147, hf_148, hf_149, \
                         kf_257, kf_258, kf_259, kf_260, kf_261, \
                         kf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = -hf_147[k]
                   + f_0 * kf_257[k];

        t_198[k] = -hf_148[k]
                   + f_0 * kf_258[k];

        t_199[k] = -hf_149[k]
                   + f_0 * kf_259[k];

        t_200[k] = f_0 * kf_260[k];

        t_201[k] = f_0 * kf_261[k];

        t_202[k] = f_0 * kf_262[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, t_208, t_209, kf_263, kf_264, \
                         kf_265, kf_266, kf_267, kf_268, kf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * kf_263[k];

        t_204[k] = f_0 * kf_264[k];

        t_205[k] = f_0 * kf_265[k];

        t_206[k] = f_0 * kf_266[k];

        t_207[k] = f_0 * kf_267[k];

        t_208[k] = f_0 * kf_268[k];

        t_209[k] = f_0 * kf_269[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, hf_150, hf_151, hf_152, hf_153, \
                         hf_154, kf_280, kf_281, kf_282, kf_283, \
                         kf_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -6.0 * hf_150[k]
                   + f_0 * kf_280[k];

        t_211[k] = -6.0 * hf_151[k]
                   + f_0 * kf_281[k];

        t_212[k] = -6.0 * hf_152[k]
                   + f_0 * kf_282[k];

        t_213[k] = -6.0 * hf_153[k]
                   + f_0 * kf_283[k];

        t_214[k] = -6.0 * hf_154[k]
                   + f_0 * kf_284[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, hf_155, hf_156, hf_157, hf_158, \
                         hf_159, kf_285, kf_286, kf_287, kf_288, \
                         kf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -6.0 * hf_155[k]
                   + f_0 * kf_285[k];

        t_216[k] = -6.0 * hf_156[k]
                   + f_0 * kf_286[k];

        t_217[k] = -6.0 * hf_157[k]
                   + f_0 * kf_287[k];

        t_218[k] = -6.0 * hf_158[k]
                   + f_0 * kf_288[k];

        t_219[k] = -6.0 * hf_159[k]
                   + f_0 * kf_289[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, hf_160, hf_161, hf_162, hf_163, \
                         hf_164, kf_290, kf_291, kf_292, kf_293, \
                         kf_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = -5.0 * hf_160[k]
                   + f_0 * kf_290[k];

        t_221[k] = -5.0 * hf_161[k]
                   + f_0 * kf_291[k];

        t_222[k] = -5.0 * hf_162[k]
                   + f_0 * kf_292[k];

        t_223[k] = -5.0 * hf_163[k]
                   + f_0 * kf_293[k];

        t_224[k] = -5.0 * hf_164[k]
                   + f_0 * kf_294[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, hf_165, hf_166, hf_167, hf_168, \
                         hf_169, kf_295, kf_296, kf_297, kf_298, \
                         kf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -5.0 * hf_165[k]
                   + f_0 * kf_295[k];

        t_226[k] = -5.0 * hf_166[k]
                   + f_0 * kf_296[k];

        t_227[k] = -5.0 * hf_167[k]
                   + f_0 * kf_297[k];

        t_228[k] = -5.0 * hf_168[k]
                   + f_0 * kf_298[k];

        t_229[k] = -5.0 * hf_169[k]
                   + f_0 * kf_299[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, hf_170, hf_171, hf_172, hf_173, \
                         hf_174, kf_300, kf_301, kf_302, kf_303, \
                         kf_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -4.0 * hf_170[k]
                   + f_0 * kf_300[k];

        t_231[k] = -4.0 * hf_171[k]
                   + f_0 * kf_301[k];

        t_232[k] = -4.0 * hf_172[k]
                   + f_0 * kf_302[k];

        t_233[k] = -4.0 * hf_173[k]
                   + f_0 * kf_303[k];

        t_234[k] = -4.0 * hf_174[k]
                   + f_0 * kf_304[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, hf_175, hf_176, hf_177, hf_178, \
                         hf_179, kf_305, kf_306, kf_307, kf_308, \
                         kf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = -4.0 * hf_175[k]
                   + f_0 * kf_305[k];

        t_236[k] = -4.0 * hf_176[k]
                   + f_0 * kf_306[k];

        t_237[k] = -4.0 * hf_177[k]
                   + f_0 * kf_307[k];

        t_238[k] = -4.0 * hf_178[k]
                   + f_0 * kf_308[k];

        t_239[k] = -4.0 * hf_179[k]
                   + f_0 * kf_309[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, hf_180, hf_181, hf_182, hf_183, \
                         hf_184, kf_310, kf_311, kf_312, kf_313, \
                         kf_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -3.0 * hf_180[k]
                   + f_0 * kf_310[k];

        t_241[k] = -3.0 * hf_181[k]
                   + f_0 * kf_311[k];

        t_242[k] = -3.0 * hf_182[k]
                   + f_0 * kf_312[k];

        t_243[k] = -3.0 * hf_183[k]
                   + f_0 * kf_313[k];

        t_244[k] = -3.0 * hf_184[k]
                   + f_0 * kf_314[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, hf_185, hf_186, hf_187, hf_188, \
                         hf_189, kf_315, kf_316, kf_317, kf_318, \
                         kf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = -3.0 * hf_185[k]
                   + f_0 * kf_315[k];

        t_246[k] = -3.0 * hf_186[k]
                   + f_0 * kf_316[k];

        t_247[k] = -3.0 * hf_187[k]
                   + f_0 * kf_317[k];

        t_248[k] = -3.0 * hf_188[k]
                   + f_0 * kf_318[k];

        t_249[k] = -3.0 * hf_189[k]
                   + f_0 * kf_319[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, hf_190, hf_191, hf_192, hf_193, \
                         hf_194, kf_320, kf_321, kf_322, kf_323, \
                         kf_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -2.0 * hf_190[k]
                   + f_0 * kf_320[k];

        t_251[k] = -2.0 * hf_191[k]
                   + f_0 * kf_321[k];

        t_252[k] = -2.0 * hf_192[k]
                   + f_0 * kf_322[k];

        t_253[k] = -2.0 * hf_193[k]
                   + f_0 * kf_323[k];

        t_254[k] = -2.0 * hf_194[k]
                   + f_0 * kf_324[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, hf_195, hf_196, hf_197, hf_198, \
                         hf_199, kf_325, kf_326, kf_327, kf_328, \
                         kf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -2.0 * hf_195[k]
                   + f_0 * kf_325[k];

        t_256[k] = -2.0 * hf_196[k]
                   + f_0 * kf_326[k];

        t_257[k] = -2.0 * hf_197[k]
                   + f_0 * kf_327[k];

        t_258[k] = -2.0 * hf_198[k]
                   + f_0 * kf_328[k];

        t_259[k] = -2.0 * hf_199[k]
                   + f_0 * kf_329[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, hf_200, hf_201, hf_202, hf_203, \
                         hf_204, kf_330, kf_331, kf_332, kf_333, \
                         kf_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = -hf_200[k]
                   + f_0 * kf_330[k];

        t_261[k] = -hf_201[k]
                   + f_0 * kf_331[k];

        t_262[k] = -hf_202[k]
                   + f_0 * kf_332[k];

        t_263[k] = -hf_203[k]
                   + f_0 * kf_333[k];

        t_264[k] = -hf_204[k]
                   + f_0 * kf_334[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, hf_205, hf_206, hf_207, hf_208, \
                         hf_209, kf_335, kf_336, kf_337, kf_338, \
                         kf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = -hf_205[k]
                   + f_0 * kf_335[k];

        t_266[k] = -hf_206[k]
                   + f_0 * kf_336[k];

        t_267[k] = -hf_207[k]
                   + f_0 * kf_337[k];

        t_268[k] = -hf_208[k]
                   + f_0 * kf_338[k];

        t_269[k] = -hf_209[k]
                   + f_0 * kf_339[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, t_276, t_277, kf_340, \
                         kf_341, kf_342, kf_343, kf_344, kf_345, kf_346, \
                         kf_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * kf_340[k];

        t_271[k] = f_0 * kf_341[k];

        t_272[k] = f_0 * kf_342[k];

        t_273[k] = f_0 * kf_343[k];

        t_274[k] = f_0 * kf_344[k];

        t_275[k] = f_0 * kf_345[k];

        t_276[k] = f_0 * kf_346[k];

        t_277[k] = f_0 * kf_347[k];
    }

#pragma omp simd aligned(t_278, t_279, kf_348, kf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * kf_348[k];

        t_279[k] = f_0 * kf_349[k];
    }
}

auto
compute_prim_geom_10_if_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t hf, const size_t kf,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_if_electron_repulsion_1_piece0(buffer, target, hf, kf, ncols, alpha);

    compute_prim_geom_10_if_electron_repulsion_1_piece1(buffer, target, hf, kf, ncols, alpha);
}

static auto
compute_prim_geom_10_if_electron_repulsion_2_piece0(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hf, const size_t kf,
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

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_102 = buffer.data(hf + 102);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_104 = buffer.data(hf + 104);
    const auto *hf_105 = buffer.data(hf + 105);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_107 = buffer.data(hf + 107);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);

    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_53 = buffer.data(kf + 53);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_93 = buffer.data(kf + 93);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_98 = buffer.data(kf + 98);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_113 = buffer.data(kf + 113);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_115 = buffer.data(kf + 115);
    const auto *kf_116 = buffer.data(kf + 116);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_121 = buffer.data(kf + 121);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_123 = buffer.data(kf + 123);
    const auto *kf_124 = buffer.data(kf + 124);
    const auto *kf_125 = buffer.data(kf + 125);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_131 = buffer.data(kf + 131);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_135 = buffer.data(kf + 135);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_143 = buffer.data(kf + 143);
    const auto *kf_144 = buffer.data(kf + 144);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_147 = buffer.data(kf + 147);
    const auto *kf_148 = buffer.data(kf + 148);
    const auto *kf_149 = buffer.data(kf + 149);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_161 = buffer.data(kf + 161);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_163 = buffer.data(kf + 163);
    const auto *kf_164 = buffer.data(kf + 164);
    const auto *kf_165 = buffer.data(kf + 165);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_171 = buffer.data(kf + 171);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_173 = buffer.data(kf + 173);
    const auto *kf_174 = buffer.data(kf + 174);
    const auto *kf_175 = buffer.data(kf + 175);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_181 = buffer.data(kf + 181);
    const auto *kf_182 = buffer.data(kf + 182);
    const auto *kf_183 = buffer.data(kf + 183);
    const auto *kf_184 = buffer.data(kf + 184);
    const auto *kf_185 = buffer.data(kf + 185);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_191 = buffer.data(kf + 191);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_193 = buffer.data(kf + 193);
    const auto *kf_194 = buffer.data(kf + 194);
    const auto *kf_195 = buffer.data(kf + 195);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_201 = buffer.data(kf + 201);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_203 = buffer.data(kf + 203);
    const auto *kf_204 = buffer.data(kf + 204);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_207 = buffer.data(kf + 207);
    const auto *kf_208 = buffer.data(kf + 208);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_221 = buffer.data(kf + 221);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_223 = buffer.data(kf + 223);
    const auto *kf_224 = buffer.data(kf + 224);
    const auto *kf_225 = buffer.data(kf + 225);
    const auto *kf_226 = buffer.data(kf + 226);
    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_231 = buffer.data(kf + 231);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_233 = buffer.data(kf + 233);
    const auto *kf_234 = buffer.data(kf + 234);
    const auto *kf_235 = buffer.data(kf + 235);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);
    const auto *kf_238 = buffer.data(kf + 238);
    const auto *kf_239 = buffer.data(kf + 239);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, kf_20, kf_21, kf_22, kf_23, \
                         kf_24, kf_25, kf_26, kf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kf_20[k];

        t_1[k] = f_0 * kf_21[k];

        t_2[k] = f_0 * kf_22[k];

        t_3[k] = f_0 * kf_23[k];

        t_4[k] = f_0 * kf_24[k];

        t_5[k] = f_0 * kf_25[k];

        t_6[k] = f_0 * kf_26[k];

        t_7[k] = f_0 * kf_27[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, kf_28, kf_29, kf_40, \
                         kf_41, kf_42, kf_43, kf_44, kf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * kf_28[k];

        t_9[k] = f_0 * kf_29[k];

        t_10[k] = f_0 * kf_40[k];

        t_11[k] = f_0 * kf_41[k];

        t_12[k] = f_0 * kf_42[k];

        t_13[k] = f_0 * kf_43[k];

        t_14[k] = f_0 * kf_44[k];

        t_15[k] = f_0 * kf_45[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, hf_0, hf_1, kf_46, kf_47, kf_48, \
                         kf_49, kf_50, kf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * kf_46[k];

        t_17[k] = f_0 * kf_47[k];

        t_18[k] = f_0 * kf_48[k];

        t_19[k] = f_0 * kf_49[k];

        t_20[k] = -hf_0[k]
                  + f_0 * kf_50[k];

        t_21[k] = -hf_1[k]
                  + f_0 * kf_51[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, hf_2, hf_3, hf_4, hf_5, hf_6, kf_52, \
                         kf_53, kf_54, kf_55, kf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -hf_2[k]
                  + f_0 * kf_52[k];

        t_23[k] = -hf_3[k]
                  + f_0 * kf_53[k];

        t_24[k] = -hf_4[k]
                  + f_0 * kf_54[k];

        t_25[k] = -hf_5[k]
                  + f_0 * kf_55[k];

        t_26[k] = -hf_6[k]
                  + f_0 * kf_56[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, hf_7, hf_8, hf_9, kf_57, kf_58, \
                         kf_59, kf_70, kf_71, kf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -hf_7[k]
                  + f_0 * kf_57[k];

        t_28[k] = -hf_8[k]
                  + f_0 * kf_58[k];

        t_29[k] = -hf_9[k]
                  + f_0 * kf_59[k];

        t_30[k] = f_0 * kf_70[k];

        t_31[k] = f_0 * kf_71[k];

        t_32[k] = f_0 * kf_72[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, t_39, kf_73, kf_74, kf_75, kf_76, \
                         kf_77, kf_78, kf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * kf_73[k];

        t_34[k] = f_0 * kf_74[k];

        t_35[k] = f_0 * kf_75[k];

        t_36[k] = f_0 * kf_76[k];

        t_37[k] = f_0 * kf_77[k];

        t_38[k] = f_0 * kf_78[k];

        t_39[k] = f_0 * kf_79[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, hf_10, hf_11, hf_12, hf_13, hf_14, \
                         kf_80, kf_81, kf_82, kf_83, kf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -hf_10[k]
                  + f_0 * kf_80[k];

        t_41[k] = -hf_11[k]
                  + f_0 * kf_81[k];

        t_42[k] = -hf_12[k]
                  + f_0 * kf_82[k];

        t_43[k] = -hf_13[k]
                  + f_0 * kf_83[k];

        t_44[k] = -hf_14[k]
                  + f_0 * kf_84[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, hf_15, hf_16, hf_17, hf_18, hf_19, \
                         kf_85, kf_86, kf_87, kf_88, kf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -hf_15[k]
                  + f_0 * kf_85[k];

        t_46[k] = -hf_16[k]
                  + f_0 * kf_86[k];

        t_47[k] = -hf_17[k]
                  + f_0 * kf_87[k];

        t_48[k] = -hf_18[k]
                  + f_0 * kf_88[k];

        t_49[k] = -hf_19[k]
                  + f_0 * kf_89[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, hf_20, hf_21, hf_22, hf_23, hf_24, \
                         kf_90, kf_91, kf_92, kf_93, kf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * hf_20[k]
                  + f_0 * kf_90[k];

        t_51[k] = -2.0 * hf_21[k]
                  + f_0 * kf_91[k];

        t_52[k] = -2.0 * hf_22[k]
                  + f_0 * kf_92[k];

        t_53[k] = -2.0 * hf_23[k]
                  + f_0 * kf_93[k];

        t_54[k] = -2.0 * hf_24[k]
                  + f_0 * kf_94[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, hf_25, hf_26, hf_27, hf_28, hf_29, \
                         kf_95, kf_96, kf_97, kf_98, kf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * hf_25[k]
                  + f_0 * kf_95[k];

        t_56[k] = -2.0 * hf_26[k]
                  + f_0 * kf_96[k];

        t_57[k] = -2.0 * hf_27[k]
                  + f_0 * kf_97[k];

        t_58[k] = -2.0 * hf_28[k]
                  + f_0 * kf_98[k];

        t_59[k] = -2.0 * hf_29[k]
                  + f_0 * kf_99[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, t_66, t_67, kf_110, kf_111, \
                         kf_112, kf_113, kf_114, kf_115, kf_116, \
                         kf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * kf_110[k];

        t_61[k] = f_0 * kf_111[k];

        t_62[k] = f_0 * kf_112[k];

        t_63[k] = f_0 * kf_113[k];

        t_64[k] = f_0 * kf_114[k];

        t_65[k] = f_0 * kf_115[k];

        t_66[k] = f_0 * kf_116[k];

        t_67[k] = f_0 * kf_117[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, hf_30, hf_31, hf_32, hf_33, \
                         kf_118, kf_119, kf_120, kf_121, kf_122, \
                         kf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * kf_118[k];

        t_69[k] = f_0 * kf_119[k];

        t_70[k] = -hf_30[k]
                  + f_0 * kf_120[k];

        t_71[k] = -hf_31[k]
                  + f_0 * kf_121[k];

        t_72[k] = -hf_32[k]
                  + f_0 * kf_122[k];

        t_73[k] = -hf_33[k]
                  + f_0 * kf_123[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, hf_34, hf_35, hf_36, hf_37, hf_38, \
                         kf_124, kf_125, kf_126, kf_127, kf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = -hf_34[k]
                  + f_0 * kf_124[k];

        t_75[k] = -hf_35[k]
                  + f_0 * kf_125[k];

        t_76[k] = -hf_36[k]
                  + f_0 * kf_126[k];

        t_77[k] = -hf_37[k]
                  + f_0 * kf_127[k];

        t_78[k] = -hf_38[k]
                  + f_0 * kf_128[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, hf_39, hf_40, hf_41, hf_42, hf_43, \
                         kf_129, kf_130, kf_131, kf_132, kf_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -hf_39[k]
                  + f_0 * kf_129[k];

        t_80[k] = -2.0 * hf_40[k]
                  + f_0 * kf_130[k];

        t_81[k] = -2.0 * hf_41[k]
                  + f_0 * kf_131[k];

        t_82[k] = -2.0 * hf_42[k]
                  + f_0 * kf_132[k];

        t_83[k] = -2.0 * hf_43[k]
                  + f_0 * kf_133[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, hf_44, hf_45, hf_46, hf_47, hf_48, \
                         kf_134, kf_135, kf_136, kf_137, kf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -2.0 * hf_44[k]
                  + f_0 * kf_134[k];

        t_85[k] = -2.0 * hf_45[k]
                  + f_0 * kf_135[k];

        t_86[k] = -2.0 * hf_46[k]
                  + f_0 * kf_136[k];

        t_87[k] = -2.0 * hf_47[k]
                  + f_0 * kf_137[k];

        t_88[k] = -2.0 * hf_48[k]
                  + f_0 * kf_138[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, hf_49, hf_50, hf_51, hf_52, hf_53, \
                         kf_139, kf_140, kf_141, kf_142, kf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -2.0 * hf_49[k]
                  + f_0 * kf_139[k];

        t_90[k] = -3.0 * hf_50[k]
                  + f_0 * kf_140[k];

        t_91[k] = -3.0 * hf_51[k]
                  + f_0 * kf_141[k];

        t_92[k] = -3.0 * hf_52[k]
                  + f_0 * kf_142[k];

        t_93[k] = -3.0 * hf_53[k]
                  + f_0 * kf_143[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, hf_54, hf_55, hf_56, hf_57, hf_58, \
                         kf_144, kf_145, kf_146, kf_147, kf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -3.0 * hf_54[k]
                  + f_0 * kf_144[k];

        t_95[k] = -3.0 * hf_55[k]
                  + f_0 * kf_145[k];

        t_96[k] = -3.0 * hf_56[k]
                  + f_0 * kf_146[k];

        t_97[k] = -3.0 * hf_57[k]
                  + f_0 * kf_147[k];

        t_98[k] = -3.0 * hf_58[k]
                  + f_0 * kf_148[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, t_105, hf_59, kf_149, \
                         kf_160, kf_161, kf_162, kf_163, kf_164, \
                         kf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -3.0 * hf_59[k]
                  + f_0 * kf_149[k];

        t_100[k] = f_0 * kf_160[k];

        t_101[k] = f_0 * kf_161[k];

        t_102[k] = f_0 * kf_162[k];

        t_103[k] = f_0 * kf_163[k];

        t_104[k] = f_0 * kf_164[k];

        t_105[k] = f_0 * kf_165[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, hf_60, hf_61, kf_166, \
                         kf_167, kf_168, kf_169, kf_170, kf_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * kf_166[k];

        t_107[k] = f_0 * kf_167[k];

        t_108[k] = f_0 * kf_168[k];

        t_109[k] = f_0 * kf_169[k];

        t_110[k] = -hf_60[k]
                   + f_0 * kf_170[k];

        t_111[k] = -hf_61[k]
                   + f_0 * kf_171[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, hf_62, hf_63, hf_64, hf_65, hf_66, \
                         kf_172, kf_173, kf_174, kf_175, kf_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -hf_62[k]
                   + f_0 * kf_172[k];

        t_113[k] = -hf_63[k]
                   + f_0 * kf_173[k];

        t_114[k] = -hf_64[k]
                   + f_0 * kf_174[k];

        t_115[k] = -hf_65[k]
                   + f_0 * kf_175[k];

        t_116[k] = -hf_66[k]
                   + f_0 * kf_176[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, hf_67, hf_68, hf_69, hf_70, hf_71, \
                         kf_177, kf_178, kf_179, kf_180, kf_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -hf_67[k]
                   + f_0 * kf_177[k];

        t_118[k] = -hf_68[k]
                   + f_0 * kf_178[k];

        t_119[k] = -hf_69[k]
                   + f_0 * kf_179[k];

        t_120[k] = -2.0 * hf_70[k]
                   + f_0 * kf_180[k];

        t_121[k] = -2.0 * hf_71[k]
                   + f_0 * kf_181[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, hf_72, hf_73, hf_74, hf_75, hf_76, \
                         kf_182, kf_183, kf_184, kf_185, kf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * hf_72[k]
                   + f_0 * kf_182[k];

        t_123[k] = -2.0 * hf_73[k]
                   + f_0 * kf_183[k];

        t_124[k] = -2.0 * hf_74[k]
                   + f_0 * kf_184[k];

        t_125[k] = -2.0 * hf_75[k]
                   + f_0 * kf_185[k];

        t_126[k] = -2.0 * hf_76[k]
                   + f_0 * kf_186[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, hf_77, hf_78, hf_79, hf_80, hf_81, \
                         kf_187, kf_188, kf_189, kf_190, kf_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * hf_77[k]
                   + f_0 * kf_187[k];

        t_128[k] = -2.0 * hf_78[k]
                   + f_0 * kf_188[k];

        t_129[k] = -2.0 * hf_79[k]
                   + f_0 * kf_189[k];

        t_130[k] = -3.0 * hf_80[k]
                   + f_0 * kf_190[k];

        t_131[k] = -3.0 * hf_81[k]
                   + f_0 * kf_191[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, hf_82, hf_83, hf_84, hf_85, hf_86, \
                         kf_192, kf_193, kf_194, kf_195, kf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -3.0 * hf_82[k]
                   + f_0 * kf_192[k];

        t_133[k] = -3.0 * hf_83[k]
                   + f_0 * kf_193[k];

        t_134[k] = -3.0 * hf_84[k]
                   + f_0 * kf_194[k];

        t_135[k] = -3.0 * hf_85[k]
                   + f_0 * kf_195[k];

        t_136[k] = -3.0 * hf_86[k]
                   + f_0 * kf_196[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, hf_87, hf_88, hf_89, hf_90, hf_91, \
                         kf_197, kf_198, kf_199, kf_200, kf_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * hf_87[k]
                   + f_0 * kf_197[k];

        t_138[k] = -3.0 * hf_88[k]
                   + f_0 * kf_198[k];

        t_139[k] = -3.0 * hf_89[k]
                   + f_0 * kf_199[k];

        t_140[k] = -4.0 * hf_90[k]
                   + f_0 * kf_200[k];

        t_141[k] = -4.0 * hf_91[k]
                   + f_0 * kf_201[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, hf_92, hf_93, hf_94, hf_95, hf_96, \
                         kf_202, kf_203, kf_204, kf_205, kf_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -4.0 * hf_92[k]
                   + f_0 * kf_202[k];

        t_143[k] = -4.0 * hf_93[k]
                   + f_0 * kf_203[k];

        t_144[k] = -4.0 * hf_94[k]
                   + f_0 * kf_204[k];

        t_145[k] = -4.0 * hf_95[k]
                   + f_0 * kf_205[k];

        t_146[k] = -4.0 * hf_96[k]
                   + f_0 * kf_206[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, t_152, hf_97, hf_98, hf_99, \
                         kf_207, kf_208, kf_209, kf_220, kf_221, \
                         kf_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -4.0 * hf_97[k]
                   + f_0 * kf_207[k];

        t_148[k] = -4.0 * hf_98[k]
                   + f_0 * kf_208[k];

        t_149[k] = -4.0 * hf_99[k]
                   + f_0 * kf_209[k];

        t_150[k] = f_0 * kf_220[k];

        t_151[k] = f_0 * kf_221[k];

        t_152[k] = f_0 * kf_222[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, t_158, t_159, kf_223, kf_224, \
                         kf_225, kf_226, kf_227, kf_228, kf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * kf_223[k];

        t_154[k] = f_0 * kf_224[k];

        t_155[k] = f_0 * kf_225[k];

        t_156[k] = f_0 * kf_226[k];

        t_157[k] = f_0 * kf_227[k];

        t_158[k] = f_0 * kf_228[k];

        t_159[k] = f_0 * kf_229[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, hf_100, hf_101, hf_102, hf_103, \
                         hf_104, kf_230, kf_231, kf_232, kf_233, \
                         kf_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -hf_100[k]
                   + f_0 * kf_230[k];

        t_161[k] = -hf_101[k]
                   + f_0 * kf_231[k];

        t_162[k] = -hf_102[k]
                   + f_0 * kf_232[k];

        t_163[k] = -hf_103[k]
                   + f_0 * kf_233[k];

        t_164[k] = -hf_104[k]
                   + f_0 * kf_234[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, hf_105, hf_106, hf_107, hf_108, \
                         hf_109, kf_235, kf_236, kf_237, kf_238, \
                         kf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -hf_105[k]
                   + f_0 * kf_235[k];

        t_166[k] = -hf_106[k]
                   + f_0 * kf_236[k];

        t_167[k] = -hf_107[k]
                   + f_0 * kf_237[k];

        t_168[k] = -hf_108[k]
                   + f_0 * kf_238[k];

        t_169[k] = -hf_109[k]
                   + f_0 * kf_239[k];
    }
}

static auto
compute_prim_geom_10_if_electron_repulsion_2_piece1(CSimdMatrix &buffer, const size_t target,
                                                    const size_t hf, const size_t kf,
                                                    const size_t ncols,
                                                    const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_111 = buffer.data(hf + 111);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_113 = buffer.data(hf + 113);
    const auto *hf_114 = buffer.data(hf + 114);
    const auto *hf_115 = buffer.data(hf + 115);
    const auto *hf_116 = buffer.data(hf + 116);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_121 = buffer.data(hf + 121);
    const auto *hf_122 = buffer.data(hf + 122);
    const auto *hf_123 = buffer.data(hf + 123);
    const auto *hf_124 = buffer.data(hf + 124);
    const auto *hf_125 = buffer.data(hf + 125);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_127 = buffer.data(hf + 127);
    const auto *hf_128 = buffer.data(hf + 128);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_130 = buffer.data(hf + 130);
    const auto *hf_131 = buffer.data(hf + 131);
    const auto *hf_132 = buffer.data(hf + 132);
    const auto *hf_133 = buffer.data(hf + 133);
    const auto *hf_134 = buffer.data(hf + 134);
    const auto *hf_135 = buffer.data(hf + 135);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_137 = buffer.data(hf + 137);
    const auto *hf_138 = buffer.data(hf + 138);
    const auto *hf_139 = buffer.data(hf + 139);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_141 = buffer.data(hf + 141);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_143 = buffer.data(hf + 143);
    const auto *hf_144 = buffer.data(hf + 144);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_146 = buffer.data(hf + 146);
    const auto *hf_147 = buffer.data(hf + 147);
    const auto *hf_148 = buffer.data(hf + 148);
    const auto *hf_149 = buffer.data(hf + 149);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_152 = buffer.data(hf + 152);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_154 = buffer.data(hf + 154);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_160 = buffer.data(hf + 160);
    const auto *hf_161 = buffer.data(hf + 161);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_163 = buffer.data(hf + 163);
    const auto *hf_164 = buffer.data(hf + 164);
    const auto *hf_165 = buffer.data(hf + 165);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_171 = buffer.data(hf + 171);
    const auto *hf_172 = buffer.data(hf + 172);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_174 = buffer.data(hf + 174);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_181 = buffer.data(hf + 181);
    const auto *hf_182 = buffer.data(hf + 182);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_184 = buffer.data(hf + 184);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_190 = buffer.data(hf + 190);
    const auto *hf_191 = buffer.data(hf + 191);
    const auto *hf_192 = buffer.data(hf + 192);
    const auto *hf_193 = buffer.data(hf + 193);
    const auto *hf_194 = buffer.data(hf + 194);
    const auto *hf_195 = buffer.data(hf + 195);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_201 = buffer.data(hf + 201);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_204 = buffer.data(hf + 204);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *kf_240 = buffer.data(kf + 240);
    const auto *kf_241 = buffer.data(kf + 241);
    const auto *kf_242 = buffer.data(kf + 242);
    const auto *kf_243 = buffer.data(kf + 243);
    const auto *kf_244 = buffer.data(kf + 244);
    const auto *kf_245 = buffer.data(kf + 245);
    const auto *kf_246 = buffer.data(kf + 246);
    const auto *kf_247 = buffer.data(kf + 247);
    const auto *kf_248 = buffer.data(kf + 248);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_250 = buffer.data(kf + 250);
    const auto *kf_251 = buffer.data(kf + 251);
    const auto *kf_252 = buffer.data(kf + 252);
    const auto *kf_253 = buffer.data(kf + 253);
    const auto *kf_254 = buffer.data(kf + 254);
    const auto *kf_255 = buffer.data(kf + 255);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_257 = buffer.data(kf + 257);
    const auto *kf_258 = buffer.data(kf + 258);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_260 = buffer.data(kf + 260);
    const auto *kf_261 = buffer.data(kf + 261);
    const auto *kf_262 = buffer.data(kf + 262);
    const auto *kf_263 = buffer.data(kf + 263);
    const auto *kf_264 = buffer.data(kf + 264);
    const auto *kf_265 = buffer.data(kf + 265);
    const auto *kf_266 = buffer.data(kf + 266);
    const auto *kf_267 = buffer.data(kf + 267);
    const auto *kf_268 = buffer.data(kf + 268);
    const auto *kf_269 = buffer.data(kf + 269);
    const auto *kf_270 = buffer.data(kf + 270);
    const auto *kf_271 = buffer.data(kf + 271);
    const auto *kf_272 = buffer.data(kf + 272);
    const auto *kf_273 = buffer.data(kf + 273);
    const auto *kf_274 = buffer.data(kf + 274);
    const auto *kf_275 = buffer.data(kf + 275);
    const auto *kf_276 = buffer.data(kf + 276);
    const auto *kf_277 = buffer.data(kf + 277);
    const auto *kf_278 = buffer.data(kf + 278);
    const auto *kf_279 = buffer.data(kf + 279);
    const auto *kf_290 = buffer.data(kf + 290);
    const auto *kf_291 = buffer.data(kf + 291);
    const auto *kf_292 = buffer.data(kf + 292);
    const auto *kf_293 = buffer.data(kf + 293);
    const auto *kf_294 = buffer.data(kf + 294);
    const auto *kf_295 = buffer.data(kf + 295);
    const auto *kf_296 = buffer.data(kf + 296);
    const auto *kf_297 = buffer.data(kf + 297);
    const auto *kf_298 = buffer.data(kf + 298);
    const auto *kf_299 = buffer.data(kf + 299);
    const auto *kf_300 = buffer.data(kf + 300);
    const auto *kf_301 = buffer.data(kf + 301);
    const auto *kf_302 = buffer.data(kf + 302);
    const auto *kf_303 = buffer.data(kf + 303);
    const auto *kf_304 = buffer.data(kf + 304);
    const auto *kf_305 = buffer.data(kf + 305);
    const auto *kf_306 = buffer.data(kf + 306);
    const auto *kf_307 = buffer.data(kf + 307);
    const auto *kf_308 = buffer.data(kf + 308);
    const auto *kf_309 = buffer.data(kf + 309);
    const auto *kf_310 = buffer.data(kf + 310);
    const auto *kf_311 = buffer.data(kf + 311);
    const auto *kf_312 = buffer.data(kf + 312);
    const auto *kf_313 = buffer.data(kf + 313);
    const auto *kf_314 = buffer.data(kf + 314);
    const auto *kf_315 = buffer.data(kf + 315);
    const auto *kf_316 = buffer.data(kf + 316);
    const auto *kf_317 = buffer.data(kf + 317);
    const auto *kf_318 = buffer.data(kf + 318);
    const auto *kf_319 = buffer.data(kf + 319);
    const auto *kf_320 = buffer.data(kf + 320);
    const auto *kf_321 = buffer.data(kf + 321);
    const auto *kf_322 = buffer.data(kf + 322);
    const auto *kf_323 = buffer.data(kf + 323);
    const auto *kf_324 = buffer.data(kf + 324);
    const auto *kf_325 = buffer.data(kf + 325);
    const auto *kf_326 = buffer.data(kf + 326);
    const auto *kf_327 = buffer.data(kf + 327);
    const auto *kf_328 = buffer.data(kf + 328);
    const auto *kf_329 = buffer.data(kf + 329);
    const auto *kf_330 = buffer.data(kf + 330);
    const auto *kf_331 = buffer.data(kf + 331);
    const auto *kf_332 = buffer.data(kf + 332);
    const auto *kf_333 = buffer.data(kf + 333);
    const auto *kf_334 = buffer.data(kf + 334);
    const auto *kf_335 = buffer.data(kf + 335);
    const auto *kf_336 = buffer.data(kf + 336);
    const auto *kf_337 = buffer.data(kf + 337);
    const auto *kf_338 = buffer.data(kf + 338);
    const auto *kf_339 = buffer.data(kf + 339);
    const auto *kf_340 = buffer.data(kf + 340);
    const auto *kf_341 = buffer.data(kf + 341);
    const auto *kf_342 = buffer.data(kf + 342);
    const auto *kf_343 = buffer.data(kf + 343);
    const auto *kf_344 = buffer.data(kf + 344);
    const auto *kf_345 = buffer.data(kf + 345);
    const auto *kf_346 = buffer.data(kf + 346);
    const auto *kf_347 = buffer.data(kf + 347);
    const auto *kf_348 = buffer.data(kf + 348);
    const auto *kf_349 = buffer.data(kf + 349);
    const auto *kf_350 = buffer.data(kf + 350);
    const auto *kf_351 = buffer.data(kf + 351);
    const auto *kf_352 = buffer.data(kf + 352);
    const auto *kf_353 = buffer.data(kf + 353);
    const auto *kf_354 = buffer.data(kf + 354);
    const auto *kf_355 = buffer.data(kf + 355);
    const auto *kf_356 = buffer.data(kf + 356);
    const auto *kf_357 = buffer.data(kf + 357);
    const auto *kf_358 = buffer.data(kf + 358);
    const auto *kf_359 = buffer.data(kf + 359);

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, hf_110, hf_111, hf_112, hf_113, \
                         hf_114, kf_240, kf_241, kf_242, kf_243, \
                         kf_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -2.0 * hf_110[k]
                   + f_0 * kf_240[k];

        t_171[k] = -2.0 * hf_111[k]
                   + f_0 * kf_241[k];

        t_172[k] = -2.0 * hf_112[k]
                   + f_0 * kf_242[k];

        t_173[k] = -2.0 * hf_113[k]
                   + f_0 * kf_243[k];

        t_174[k] = -2.0 * hf_114[k]
                   + f_0 * kf_244[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, hf_115, hf_116, hf_117, hf_118, \
                         hf_119, kf_245, kf_246, kf_247, kf_248, \
                         kf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = -2.0 * hf_115[k]
                   + f_0 * kf_245[k];

        t_176[k] = -2.0 * hf_116[k]
                   + f_0 * kf_246[k];

        t_177[k] = -2.0 * hf_117[k]
                   + f_0 * kf_247[k];

        t_178[k] = -2.0 * hf_118[k]
                   + f_0 * kf_248[k];

        t_179[k] = -2.0 * hf_119[k]
                   + f_0 * kf_249[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, hf_120, hf_121, hf_122, hf_123, \
                         hf_124, kf_250, kf_251, kf_252, kf_253, \
                         kf_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -3.0 * hf_120[k]
                   + f_0 * kf_250[k];

        t_181[k] = -3.0 * hf_121[k]
                   + f_0 * kf_251[k];

        t_182[k] = -3.0 * hf_122[k]
                   + f_0 * kf_252[k];

        t_183[k] = -3.0 * hf_123[k]
                   + f_0 * kf_253[k];

        t_184[k] = -3.0 * hf_124[k]
                   + f_0 * kf_254[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, hf_125, hf_126, hf_127, hf_128, \
                         hf_129, kf_255, kf_256, kf_257, kf_258, \
                         kf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -3.0 * hf_125[k]
                   + f_0 * kf_255[k];

        t_186[k] = -3.0 * hf_126[k]
                   + f_0 * kf_256[k];

        t_187[k] = -3.0 * hf_127[k]
                   + f_0 * kf_257[k];

        t_188[k] = -3.0 * hf_128[k]
                   + f_0 * kf_258[k];

        t_189[k] = -3.0 * hf_129[k]
                   + f_0 * kf_259[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, hf_130, hf_131, hf_132, hf_133, \
                         hf_134, kf_260, kf_261, kf_262, kf_263, \
                         kf_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -4.0 * hf_130[k]
                   + f_0 * kf_260[k];

        t_191[k] = -4.0 * hf_131[k]
                   + f_0 * kf_261[k];

        t_192[k] = -4.0 * hf_132[k]
                   + f_0 * kf_262[k];

        t_193[k] = -4.0 * hf_133[k]
                   + f_0 * kf_263[k];

        t_194[k] = -4.0 * hf_134[k]
                   + f_0 * kf_264[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, hf_135, hf_136, hf_137, hf_138, \
                         hf_139, kf_265, kf_266, kf_267, kf_268, \
                         kf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -4.0 * hf_135[k]
                   + f_0 * kf_265[k];

        t_196[k] = -4.0 * hf_136[k]
                   + f_0 * kf_266[k];

        t_197[k] = -4.0 * hf_137[k]
                   + f_0 * kf_267[k];

        t_198[k] = -4.0 * hf_138[k]
                   + f_0 * kf_268[k];

        t_199[k] = -4.0 * hf_139[k]
                   + f_0 * kf_269[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, hf_140, hf_141, hf_142, hf_143, \
                         hf_144, kf_270, kf_271, kf_272, kf_273, \
                         kf_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -5.0 * hf_140[k]
                   + f_0 * kf_270[k];

        t_201[k] = -5.0 * hf_141[k]
                   + f_0 * kf_271[k];

        t_202[k] = -5.0 * hf_142[k]
                   + f_0 * kf_272[k];

        t_203[k] = -5.0 * hf_143[k]
                   + f_0 * kf_273[k];

        t_204[k] = -5.0 * hf_144[k]
                   + f_0 * kf_274[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, hf_145, hf_146, hf_147, hf_148, \
                         hf_149, kf_275, kf_276, kf_277, kf_278, \
                         kf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -5.0 * hf_145[k]
                   + f_0 * kf_275[k];

        t_206[k] = -5.0 * hf_146[k]
                   + f_0 * kf_276[k];

        t_207[k] = -5.0 * hf_147[k]
                   + f_0 * kf_277[k];

        t_208[k] = -5.0 * hf_148[k]
                   + f_0 * kf_278[k];

        t_209[k] = -5.0 * hf_149[k]
                   + f_0 * kf_279[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, t_217, kf_290, \
                         kf_291, kf_292, kf_293, kf_294, kf_295, kf_296, \
                         kf_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * kf_290[k];

        t_211[k] = f_0 * kf_291[k];

        t_212[k] = f_0 * kf_292[k];

        t_213[k] = f_0 * kf_293[k];

        t_214[k] = f_0 * kf_294[k];

        t_215[k] = f_0 * kf_295[k];

        t_216[k] = f_0 * kf_296[k];

        t_217[k] = f_0 * kf_297[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, hf_150, hf_151, hf_152, \
                         hf_153, kf_298, kf_299, kf_300, kf_301, kf_302, \
                         kf_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * kf_298[k];

        t_219[k] = f_0 * kf_299[k];

        t_220[k] = -hf_150[k]
                   + f_0 * kf_300[k];

        t_221[k] = -hf_151[k]
                   + f_0 * kf_301[k];

        t_222[k] = -hf_152[k]
                   + f_0 * kf_302[k];

        t_223[k] = -hf_153[k]
                   + f_0 * kf_303[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, hf_154, hf_155, hf_156, hf_157, \
                         hf_158, kf_304, kf_305, kf_306, kf_307, \
                         kf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = -hf_154[k]
                   + f_0 * kf_304[k];

        t_225[k] = -hf_155[k]
                   + f_0 * kf_305[k];

        t_226[k] = -hf_156[k]
                   + f_0 * kf_306[k];

        t_227[k] = -hf_157[k]
                   + f_0 * kf_307[k];

        t_228[k] = -hf_158[k]
                   + f_0 * kf_308[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, hf_159, hf_160, hf_161, hf_162, \
                         hf_163, kf_309, kf_310, kf_311, kf_312, \
                         kf_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = -hf_159[k]
                   + f_0 * kf_309[k];

        t_230[k] = -2.0 * hf_160[k]
                   + f_0 * kf_310[k];

        t_231[k] = -2.0 * hf_161[k]
                   + f_0 * kf_311[k];

        t_232[k] = -2.0 * hf_162[k]
                   + f_0 * kf_312[k];

        t_233[k] = -2.0 * hf_163[k]
                   + f_0 * kf_313[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, hf_164, hf_165, hf_166, hf_167, \
                         hf_168, kf_314, kf_315, kf_316, kf_317, \
                         kf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = -2.0 * hf_164[k]
                   + f_0 * kf_314[k];

        t_235[k] = -2.0 * hf_165[k]
                   + f_0 * kf_315[k];

        t_236[k] = -2.0 * hf_166[k]
                   + f_0 * kf_316[k];

        t_237[k] = -2.0 * hf_167[k]
                   + f_0 * kf_317[k];

        t_238[k] = -2.0 * hf_168[k]
                   + f_0 * kf_318[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, hf_169, hf_170, hf_171, hf_172, \
                         hf_173, kf_319, kf_320, kf_321, kf_322, \
                         kf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = -2.0 * hf_169[k]
                   + f_0 * kf_319[k];

        t_240[k] = -3.0 * hf_170[k]
                   + f_0 * kf_320[k];

        t_241[k] = -3.0 * hf_171[k]
                   + f_0 * kf_321[k];

        t_242[k] = -3.0 * hf_172[k]
                   + f_0 * kf_322[k];

        t_243[k] = -3.0 * hf_173[k]
                   + f_0 * kf_323[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, t_248, hf_174, hf_175, hf_176, hf_177, \
                         hf_178, kf_324, kf_325, kf_326, kf_327, \
                         kf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = -3.0 * hf_174[k]
                   + f_0 * kf_324[k];

        t_245[k] = -3.0 * hf_175[k]
                   + f_0 * kf_325[k];

        t_246[k] = -3.0 * hf_176[k]
                   + f_0 * kf_326[k];

        t_247[k] = -3.0 * hf_177[k]
                   + f_0 * kf_327[k];

        t_248[k] = -3.0 * hf_178[k]
                   + f_0 * kf_328[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, hf_179, hf_180, hf_181, hf_182, \
                         hf_183, kf_329, kf_330, kf_331, kf_332, \
                         kf_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = -3.0 * hf_179[k]
                   + f_0 * kf_329[k];

        t_250[k] = -4.0 * hf_180[k]
                   + f_0 * kf_330[k];

        t_251[k] = -4.0 * hf_181[k]
                   + f_0 * kf_331[k];

        t_252[k] = -4.0 * hf_182[k]
                   + f_0 * kf_332[k];

        t_253[k] = -4.0 * hf_183[k]
                   + f_0 * kf_333[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, hf_184, hf_185, hf_186, hf_187, \
                         hf_188, kf_334, kf_335, kf_336, kf_337, \
                         kf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = -4.0 * hf_184[k]
                   + f_0 * kf_334[k];

        t_255[k] = -4.0 * hf_185[k]
                   + f_0 * kf_335[k];

        t_256[k] = -4.0 * hf_186[k]
                   + f_0 * kf_336[k];

        t_257[k] = -4.0 * hf_187[k]
                   + f_0 * kf_337[k];

        t_258[k] = -4.0 * hf_188[k]
                   + f_0 * kf_338[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, hf_189, hf_190, hf_191, hf_192, \
                         hf_193, kf_339, kf_340, kf_341, kf_342, \
                         kf_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = -4.0 * hf_189[k]
                   + f_0 * kf_339[k];

        t_260[k] = -5.0 * hf_190[k]
                   + f_0 * kf_340[k];

        t_261[k] = -5.0 * hf_191[k]
                   + f_0 * kf_341[k];

        t_262[k] = -5.0 * hf_192[k]
                   + f_0 * kf_342[k];

        t_263[k] = -5.0 * hf_193[k]
                   + f_0 * kf_343[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, hf_194, hf_195, hf_196, hf_197, \
                         hf_198, kf_344, kf_345, kf_346, kf_347, \
                         kf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = -5.0 * hf_194[k]
                   + f_0 * kf_344[k];

        t_265[k] = -5.0 * hf_195[k]
                   + f_0 * kf_345[k];

        t_266[k] = -5.0 * hf_196[k]
                   + f_0 * kf_346[k];

        t_267[k] = -5.0 * hf_197[k]
                   + f_0 * kf_347[k];

        t_268[k] = -5.0 * hf_198[k]
                   + f_0 * kf_348[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, t_273, hf_199, hf_200, hf_201, hf_202, \
                         hf_203, kf_349, kf_350, kf_351, kf_352, \
                         kf_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = -5.0 * hf_199[k]
                   + f_0 * kf_349[k];

        t_270[k] = -6.0 * hf_200[k]
                   + f_0 * kf_350[k];

        t_271[k] = -6.0 * hf_201[k]
                   + f_0 * kf_351[k];

        t_272[k] = -6.0 * hf_202[k]
                   + f_0 * kf_352[k];

        t_273[k] = -6.0 * hf_203[k]
                   + f_0 * kf_353[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, hf_204, hf_205, hf_206, hf_207, \
                         hf_208, kf_354, kf_355, kf_356, kf_357, \
                         kf_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = -6.0 * hf_204[k]
                   + f_0 * kf_354[k];

        t_275[k] = -6.0 * hf_205[k]
                   + f_0 * kf_355[k];

        t_276[k] = -6.0 * hf_206[k]
                   + f_0 * kf_356[k];

        t_277[k] = -6.0 * hf_207[k]
                   + f_0 * kf_357[k];

        t_278[k] = -6.0 * hf_208[k]
                   + f_0 * kf_358[k];
    }

#pragma omp simd aligned(t_279, hf_209, kf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = -6.0 * hf_209[k]
                   + f_0 * kf_359[k];
    }
}

auto
compute_prim_geom_10_if_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t hf, const size_t kf,
                                             const size_t ncols, const double alpha) -> void
{
    compute_prim_geom_10_if_electron_repulsion_2_piece0(buffer, target, hf, kf, ncols, alpha);

    compute_prim_geom_10_if_electron_repulsion_2_piece1(buffer, target, hf, kf, ncols, alpha);
}

}  // namespace simdt2ceri
