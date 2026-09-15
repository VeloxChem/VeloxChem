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


#include "SimdElectronRepulsionGeom10VrrRecGF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_gf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ff, const size_t hf,
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

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);
    const auto *ff_54 = buffer.data(ff + 54);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_58 = buffer.data(ff + 58);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_62 = buffer.data(ff + 62);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_64 = buffer.data(ff + 64);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_70 = buffer.data(ff + 70);
    const auto *ff_71 = buffer.data(ff + 71);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_73 = buffer.data(ff + 73);
    const auto *ff_74 = buffer.data(ff + 74);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_80 = buffer.data(ff + 80);
    const auto *ff_81 = buffer.data(ff + 81);
    const auto *ff_82 = buffer.data(ff + 82);
    const auto *ff_83 = buffer.data(ff + 83);
    const auto *ff_84 = buffer.data(ff + 84);
    const auto *ff_85 = buffer.data(ff + 85);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_91 = buffer.data(ff + 91);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_94 = buffer.data(ff + 94);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_98 = buffer.data(ff + 98);
    const auto *ff_99 = buffer.data(ff + 99);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ff_0, ff_1, ff_2, ff_3, ff_4, hf_0, hf_1, \
                         hf_2, hf_3, hf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -4.0 * ff_0[k]
                 + f_0 * hf_0[k];

        t_1[k] = -4.0 * ff_1[k]
                 + f_0 * hf_1[k];

        t_2[k] = -4.0 * ff_2[k]
                 + f_0 * hf_2[k];

        t_3[k] = -4.0 * ff_3[k]
                 + f_0 * hf_3[k];

        t_4[k] = -4.0 * ff_4[k]
                 + f_0 * hf_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ff_5, ff_6, ff_7, ff_8, ff_9, hf_5, hf_6, \
                         hf_7, hf_8, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -4.0 * ff_5[k]
                 + f_0 * hf_5[k];

        t_6[k] = -4.0 * ff_6[k]
                 + f_0 * hf_6[k];

        t_7[k] = -4.0 * ff_7[k]
                 + f_0 * hf_7[k];

        t_8[k] = -4.0 * ff_8[k]
                 + f_0 * hf_8[k];

        t_9[k] = -4.0 * ff_9[k]
                 + f_0 * hf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ff_10, ff_11, ff_12, ff_13, ff_14, \
                         hf_10, hf_11, hf_12, hf_13, hf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -3.0 * ff_10[k]
                  + f_0 * hf_10[k];

        t_11[k] = -3.0 * ff_11[k]
                  + f_0 * hf_11[k];

        t_12[k] = -3.0 * ff_12[k]
                  + f_0 * hf_12[k];

        t_13[k] = -3.0 * ff_13[k]
                  + f_0 * hf_13[k];

        t_14[k] = -3.0 * ff_14[k]
                  + f_0 * hf_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ff_15, ff_16, ff_17, ff_18, ff_19, \
                         hf_15, hf_16, hf_17, hf_18, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -3.0 * ff_15[k]
                  + f_0 * hf_15[k];

        t_16[k] = -3.0 * ff_16[k]
                  + f_0 * hf_16[k];

        t_17[k] = -3.0 * ff_17[k]
                  + f_0 * hf_17[k];

        t_18[k] = -3.0 * ff_18[k]
                  + f_0 * hf_18[k];

        t_19[k] = -3.0 * ff_19[k]
                  + f_0 * hf_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ff_20, ff_21, ff_22, ff_23, ff_24, \
                         hf_20, hf_21, hf_22, hf_23, hf_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -3.0 * ff_20[k]
                  + f_0 * hf_20[k];

        t_21[k] = -3.0 * ff_21[k]
                  + f_0 * hf_21[k];

        t_22[k] = -3.0 * ff_22[k]
                  + f_0 * hf_22[k];

        t_23[k] = -3.0 * ff_23[k]
                  + f_0 * hf_23[k];

        t_24[k] = -3.0 * ff_24[k]
                  + f_0 * hf_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ff_25, ff_26, ff_27, ff_28, ff_29, \
                         hf_25, hf_26, hf_27, hf_28, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -3.0 * ff_25[k]
                  + f_0 * hf_25[k];

        t_26[k] = -3.0 * ff_26[k]
                  + f_0 * hf_26[k];

        t_27[k] = -3.0 * ff_27[k]
                  + f_0 * hf_27[k];

        t_28[k] = -3.0 * ff_28[k]
                  + f_0 * hf_28[k];

        t_29[k] = -3.0 * ff_29[k]
                  + f_0 * hf_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ff_30, ff_31, ff_32, ff_33, ff_34, \
                         hf_30, hf_31, hf_32, hf_33, hf_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -2.0 * ff_30[k]
                  + f_0 * hf_30[k];

        t_31[k] = -2.0 * ff_31[k]
                  + f_0 * hf_31[k];

        t_32[k] = -2.0 * ff_32[k]
                  + f_0 * hf_32[k];

        t_33[k] = -2.0 * ff_33[k]
                  + f_0 * hf_33[k];

        t_34[k] = -2.0 * ff_34[k]
                  + f_0 * hf_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ff_35, ff_36, ff_37, ff_38, ff_39, \
                         hf_35, hf_36, hf_37, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -2.0 * ff_35[k]
                  + f_0 * hf_35[k];

        t_36[k] = -2.0 * ff_36[k]
                  + f_0 * hf_36[k];

        t_37[k] = -2.0 * ff_37[k]
                  + f_0 * hf_37[k];

        t_38[k] = -2.0 * ff_38[k]
                  + f_0 * hf_38[k];

        t_39[k] = -2.0 * ff_39[k]
                  + f_0 * hf_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ff_40, ff_41, ff_42, ff_43, ff_44, \
                         hf_40, hf_41, hf_42, hf_43, hf_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -2.0 * ff_40[k]
                  + f_0 * hf_40[k];

        t_41[k] = -2.0 * ff_41[k]
                  + f_0 * hf_41[k];

        t_42[k] = -2.0 * ff_42[k]
                  + f_0 * hf_42[k];

        t_43[k] = -2.0 * ff_43[k]
                  + f_0 * hf_43[k];

        t_44[k] = -2.0 * ff_44[k]
                  + f_0 * hf_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ff_45, ff_46, ff_47, ff_48, ff_49, \
                         hf_45, hf_46, hf_47, hf_48, hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * ff_45[k]
                  + f_0 * hf_45[k];

        t_46[k] = -2.0 * ff_46[k]
                  + f_0 * hf_46[k];

        t_47[k] = -2.0 * ff_47[k]
                  + f_0 * hf_47[k];

        t_48[k] = -2.0 * ff_48[k]
                  + f_0 * hf_48[k];

        t_49[k] = -2.0 * ff_49[k]
                  + f_0 * hf_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ff_50, ff_51, ff_52, ff_53, ff_54, \
                         hf_50, hf_51, hf_52, hf_53, hf_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * ff_50[k]
                  + f_0 * hf_50[k];

        t_51[k] = -2.0 * ff_51[k]
                  + f_0 * hf_51[k];

        t_52[k] = -2.0 * ff_52[k]
                  + f_0 * hf_52[k];

        t_53[k] = -2.0 * ff_53[k]
                  + f_0 * hf_53[k];

        t_54[k] = -2.0 * ff_54[k]
                  + f_0 * hf_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ff_55, ff_56, ff_57, ff_58, ff_59, \
                         hf_55, hf_56, hf_57, hf_58, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * ff_55[k]
                  + f_0 * hf_55[k];

        t_56[k] = -2.0 * ff_56[k]
                  + f_0 * hf_56[k];

        t_57[k] = -2.0 * ff_57[k]
                  + f_0 * hf_57[k];

        t_58[k] = -2.0 * ff_58[k]
                  + f_0 * hf_58[k];

        t_59[k] = -2.0 * ff_59[k]
                  + f_0 * hf_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ff_60, ff_61, ff_62, ff_63, ff_64, \
                         hf_60, hf_61, hf_62, hf_63, hf_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -ff_60[k]
                  + f_0 * hf_60[k];

        t_61[k] = -ff_61[k]
                  + f_0 * hf_61[k];

        t_62[k] = -ff_62[k]
                  + f_0 * hf_62[k];

        t_63[k] = -ff_63[k]
                  + f_0 * hf_63[k];

        t_64[k] = -ff_64[k]
                  + f_0 * hf_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ff_65, ff_66, ff_67, ff_68, ff_69, \
                         hf_65, hf_66, hf_67, hf_68, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -ff_65[k]
                  + f_0 * hf_65[k];

        t_66[k] = -ff_66[k]
                  + f_0 * hf_66[k];

        t_67[k] = -ff_67[k]
                  + f_0 * hf_67[k];

        t_68[k] = -ff_68[k]
                  + f_0 * hf_68[k];

        t_69[k] = -ff_69[k]
                  + f_0 * hf_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ff_70, ff_71, ff_72, ff_73, ff_74, \
                         hf_70, hf_71, hf_72, hf_73, hf_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -ff_70[k]
                  + f_0 * hf_70[k];

        t_71[k] = -ff_71[k]
                  + f_0 * hf_71[k];

        t_72[k] = -ff_72[k]
                  + f_0 * hf_72[k];

        t_73[k] = -ff_73[k]
                  + f_0 * hf_73[k];

        t_74[k] = -ff_74[k]
                  + f_0 * hf_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ff_75, ff_76, ff_77, ff_78, ff_79, \
                         hf_75, hf_76, hf_77, hf_78, hf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -ff_75[k]
                  + f_0 * hf_75[k];

        t_76[k] = -ff_76[k]
                  + f_0 * hf_76[k];

        t_77[k] = -ff_77[k]
                  + f_0 * hf_77[k];

        t_78[k] = -ff_78[k]
                  + f_0 * hf_78[k];

        t_79[k] = -ff_79[k]
                  + f_0 * hf_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ff_80, ff_81, ff_82, ff_83, ff_84, \
                         hf_80, hf_81, hf_82, hf_83, hf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -ff_80[k]
                  + f_0 * hf_80[k];

        t_81[k] = -ff_81[k]
                  + f_0 * hf_81[k];

        t_82[k] = -ff_82[k]
                  + f_0 * hf_82[k];

        t_83[k] = -ff_83[k]
                  + f_0 * hf_83[k];

        t_84[k] = -ff_84[k]
                  + f_0 * hf_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ff_85, ff_86, ff_87, ff_88, ff_89, \
                         hf_85, hf_86, hf_87, hf_88, hf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -ff_85[k]
                  + f_0 * hf_85[k];

        t_86[k] = -ff_86[k]
                  + f_0 * hf_86[k];

        t_87[k] = -ff_87[k]
                  + f_0 * hf_87[k];

        t_88[k] = -ff_88[k]
                  + f_0 * hf_88[k];

        t_89[k] = -ff_89[k]
                  + f_0 * hf_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ff_90, ff_91, ff_92, ff_93, ff_94, \
                         hf_90, hf_91, hf_92, hf_93, hf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -ff_90[k]
                  + f_0 * hf_90[k];

        t_91[k] = -ff_91[k]
                  + f_0 * hf_91[k];

        t_92[k] = -ff_92[k]
                  + f_0 * hf_92[k];

        t_93[k] = -ff_93[k]
                  + f_0 * hf_93[k];

        t_94[k] = -ff_94[k]
                  + f_0 * hf_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ff_95, ff_96, ff_97, ff_98, ff_99, \
                         hf_95, hf_96, hf_97, hf_98, hf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -ff_95[k]
                  + f_0 * hf_95[k];

        t_96[k] = -ff_96[k]
                  + f_0 * hf_96[k];

        t_97[k] = -ff_97[k]
                  + f_0 * hf_97[k];

        t_98[k] = -ff_98[k]
                  + f_0 * hf_98[k];

        t_99[k] = -ff_99[k]
                  + f_0 * hf_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, t_105, t_106, t_107, hf_100, \
                         hf_101, hf_102, hf_103, hf_104, hf_105, hf_106, \
                         hf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_0 * hf_100[k];

        t_101[k] = f_0 * hf_101[k];

        t_102[k] = f_0 * hf_102[k];

        t_103[k] = f_0 * hf_103[k];

        t_104[k] = f_0 * hf_104[k];

        t_105[k] = f_0 * hf_105[k];

        t_106[k] = f_0 * hf_106[k];

        t_107[k] = f_0 * hf_107[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, t_114, t_115, hf_108, \
                         hf_109, hf_110, hf_111, hf_112, hf_113, hf_114, \
                         hf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * hf_108[k];

        t_109[k] = f_0 * hf_109[k];

        t_110[k] = f_0 * hf_110[k];

        t_111[k] = f_0 * hf_111[k];

        t_112[k] = f_0 * hf_112[k];

        t_113[k] = f_0 * hf_113[k];

        t_114[k] = f_0 * hf_114[k];

        t_115[k] = f_0 * hf_115[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, t_120, t_121, t_122, t_123, hf_116, \
                         hf_117, hf_118, hf_119, hf_120, hf_121, hf_122, \
                         hf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_0 * hf_116[k];

        t_117[k] = f_0 * hf_117[k];

        t_118[k] = f_0 * hf_118[k];

        t_119[k] = f_0 * hf_119[k];

        t_120[k] = f_0 * hf_120[k];

        t_121[k] = f_0 * hf_121[k];

        t_122[k] = f_0 * hf_122[k];

        t_123[k] = f_0 * hf_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, t_129, t_130, t_131, hf_124, \
                         hf_125, hf_126, hf_127, hf_128, hf_129, hf_130, \
                         hf_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_0 * hf_124[k];

        t_125[k] = f_0 * hf_125[k];

        t_126[k] = f_0 * hf_126[k];

        t_127[k] = f_0 * hf_127[k];

        t_128[k] = f_0 * hf_128[k];

        t_129[k] = f_0 * hf_129[k];

        t_130[k] = f_0 * hf_130[k];

        t_131[k] = f_0 * hf_131[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, t_138, t_139, hf_132, \
                         hf_133, hf_134, hf_135, hf_136, hf_137, hf_138, \
                         hf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_0 * hf_132[k];

        t_133[k] = f_0 * hf_133[k];

        t_134[k] = f_0 * hf_134[k];

        t_135[k] = f_0 * hf_135[k];

        t_136[k] = f_0 * hf_136[k];

        t_137[k] = f_0 * hf_137[k];

        t_138[k] = f_0 * hf_138[k];

        t_139[k] = f_0 * hf_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, t_146, t_147, hf_140, \
                         hf_141, hf_142, hf_143, hf_144, hf_145, hf_146, \
                         hf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * hf_140[k];

        t_141[k] = f_0 * hf_141[k];

        t_142[k] = f_0 * hf_142[k];

        t_143[k] = f_0 * hf_143[k];

        t_144[k] = f_0 * hf_144[k];

        t_145[k] = f_0 * hf_145[k];

        t_146[k] = f_0 * hf_146[k];

        t_147[k] = f_0 * hf_147[k];
    }

#pragma omp simd aligned(t_148, t_149, hf_148, hf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * hf_148[k];

        t_149[k] = f_0 * hf_149[k];
    }
}

auto
compute_prim_geom_10_gf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ff, const size_t hf,
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

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);
    const auto *ff_54 = buffer.data(ff + 54);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_58 = buffer.data(ff + 58);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_62 = buffer.data(ff + 62);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_64 = buffer.data(ff + 64);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_70 = buffer.data(ff + 70);
    const auto *ff_71 = buffer.data(ff + 71);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_73 = buffer.data(ff + 73);
    const auto *ff_74 = buffer.data(ff + 74);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_80 = buffer.data(ff + 80);
    const auto *ff_81 = buffer.data(ff + 81);
    const auto *ff_82 = buffer.data(ff + 82);
    const auto *ff_83 = buffer.data(ff + 83);
    const auto *ff_84 = buffer.data(ff + 84);
    const auto *ff_85 = buffer.data(ff + 85);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_91 = buffer.data(ff + 91);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_94 = buffer.data(ff + 94);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_98 = buffer.data(ff + 98);
    const auto *ff_99 = buffer.data(ff + 99);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hf_10, hf_11, hf_12, hf_13, \
                         hf_14, hf_15, hf_16, hf_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_10[k];

        t_1[k] = f_0 * hf_11[k];

        t_2[k] = f_0 * hf_12[k];

        t_3[k] = f_0 * hf_13[k];

        t_4[k] = f_0 * hf_14[k];

        t_5[k] = f_0 * hf_15[k];

        t_6[k] = f_0 * hf_16[k];

        t_7[k] = f_0 * hf_17[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, ff_0, ff_1, ff_2, ff_3, hf_18, \
                         hf_19, hf_30, hf_31, hf_32, hf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hf_18[k];

        t_9[k] = f_0 * hf_19[k];

        t_10[k] = -ff_0[k]
                  + f_0 * hf_30[k];

        t_11[k] = -ff_1[k]
                  + f_0 * hf_31[k];

        t_12[k] = -ff_2[k]
                  + f_0 * hf_32[k];

        t_13[k] = -ff_3[k]
                  + f_0 * hf_33[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, ff_4, ff_5, ff_6, ff_7, ff_8, hf_34, \
                         hf_35, hf_36, hf_37, hf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -ff_4[k]
                  + f_0 * hf_34[k];

        t_15[k] = -ff_5[k]
                  + f_0 * hf_35[k];

        t_16[k] = -ff_6[k]
                  + f_0 * hf_36[k];

        t_17[k] = -ff_7[k]
                  + f_0 * hf_37[k];

        t_18[k] = -ff_8[k]
                  + f_0 * hf_38[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, t_24, t_25, ff_9, hf_39, hf_40, hf_41, \
                         hf_42, hf_43, hf_44, hf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -ff_9[k]
                  + f_0 * hf_39[k];

        t_20[k] = f_0 * hf_40[k];

        t_21[k] = f_0 * hf_41[k];

        t_22[k] = f_0 * hf_42[k];

        t_23[k] = f_0 * hf_43[k];

        t_24[k] = f_0 * hf_44[k];

        t_25[k] = f_0 * hf_45[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, t_31, ff_10, ff_11, hf_46, hf_47, \
                         hf_48, hf_49, hf_60, hf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hf_46[k];

        t_27[k] = f_0 * hf_47[k];

        t_28[k] = f_0 * hf_48[k];

        t_29[k] = f_0 * hf_49[k];

        t_30[k] = -2.0 * ff_10[k]
                  + f_0 * hf_60[k];

        t_31[k] = -2.0 * ff_11[k]
                  + f_0 * hf_61[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, ff_12, ff_13, ff_14, ff_15, ff_16, \
                         hf_62, hf_63, hf_64, hf_65, hf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -2.0 * ff_12[k]
                  + f_0 * hf_62[k];

        t_33[k] = -2.0 * ff_13[k]
                  + f_0 * hf_63[k];

        t_34[k] = -2.0 * ff_14[k]
                  + f_0 * hf_64[k];

        t_35[k] = -2.0 * ff_15[k]
                  + f_0 * hf_65[k];

        t_36[k] = -2.0 * ff_16[k]
                  + f_0 * hf_66[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, ff_17, ff_18, ff_19, ff_20, ff_21, \
                         hf_67, hf_68, hf_69, hf_70, hf_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -2.0 * ff_17[k]
                  + f_0 * hf_67[k];

        t_38[k] = -2.0 * ff_18[k]
                  + f_0 * hf_68[k];

        t_39[k] = -2.0 * ff_19[k]
                  + f_0 * hf_69[k];

        t_40[k] = -ff_20[k]
                  + f_0 * hf_70[k];

        t_41[k] = -ff_21[k]
                  + f_0 * hf_71[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, ff_22, ff_23, ff_24, ff_25, ff_26, \
                         hf_72, hf_73, hf_74, hf_75, hf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -ff_22[k]
                  + f_0 * hf_72[k];

        t_43[k] = -ff_23[k]
                  + f_0 * hf_73[k];

        t_44[k] = -ff_24[k]
                  + f_0 * hf_74[k];

        t_45[k] = -ff_25[k]
                  + f_0 * hf_75[k];

        t_46[k] = -ff_26[k]
                  + f_0 * hf_76[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, t_52, ff_27, ff_28, ff_29, hf_77, \
                         hf_78, hf_79, hf_80, hf_81, hf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -ff_27[k]
                  + f_0 * hf_77[k];

        t_48[k] = -ff_28[k]
                  + f_0 * hf_78[k];

        t_49[k] = -ff_29[k]
                  + f_0 * hf_79[k];

        t_50[k] = f_0 * hf_80[k];

        t_51[k] = f_0 * hf_81[k];

        t_52[k] = f_0 * hf_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, hf_83, hf_84, hf_85, hf_86, \
                         hf_87, hf_88, hf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * hf_83[k];

        t_54[k] = f_0 * hf_84[k];

        t_55[k] = f_0 * hf_85[k];

        t_56[k] = f_0 * hf_86[k];

        t_57[k] = f_0 * hf_87[k];

        t_58[k] = f_0 * hf_88[k];

        t_59[k] = f_0 * hf_89[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ff_30, ff_31, ff_32, ff_33, ff_34, \
                         hf_100, hf_101, hf_102, hf_103, hf_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -3.0 * ff_30[k]
                  + f_0 * hf_100[k];

        t_61[k] = -3.0 * ff_31[k]
                  + f_0 * hf_101[k];

        t_62[k] = -3.0 * ff_32[k]
                  + f_0 * hf_102[k];

        t_63[k] = -3.0 * ff_33[k]
                  + f_0 * hf_103[k];

        t_64[k] = -3.0 * ff_34[k]
                  + f_0 * hf_104[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ff_35, ff_36, ff_37, ff_38, ff_39, \
                         hf_105, hf_106, hf_107, hf_108, hf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -3.0 * ff_35[k]
                  + f_0 * hf_105[k];

        t_66[k] = -3.0 * ff_36[k]
                  + f_0 * hf_106[k];

        t_67[k] = -3.0 * ff_37[k]
                  + f_0 * hf_107[k];

        t_68[k] = -3.0 * ff_38[k]
                  + f_0 * hf_108[k];

        t_69[k] = -3.0 * ff_39[k]
                  + f_0 * hf_109[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ff_40, ff_41, ff_42, ff_43, ff_44, \
                         hf_110, hf_111, hf_112, hf_113, hf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -2.0 * ff_40[k]
                  + f_0 * hf_110[k];

        t_71[k] = -2.0 * ff_41[k]
                  + f_0 * hf_111[k];

        t_72[k] = -2.0 * ff_42[k]
                  + f_0 * hf_112[k];

        t_73[k] = -2.0 * ff_43[k]
                  + f_0 * hf_113[k];

        t_74[k] = -2.0 * ff_44[k]
                  + f_0 * hf_114[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ff_45, ff_46, ff_47, ff_48, ff_49, \
                         hf_115, hf_116, hf_117, hf_118, hf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * ff_45[k]
                  + f_0 * hf_115[k];

        t_76[k] = -2.0 * ff_46[k]
                  + f_0 * hf_116[k];

        t_77[k] = -2.0 * ff_47[k]
                  + f_0 * hf_117[k];

        t_78[k] = -2.0 * ff_48[k]
                  + f_0 * hf_118[k];

        t_79[k] = -2.0 * ff_49[k]
                  + f_0 * hf_119[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ff_50, ff_51, ff_52, ff_53, ff_54, \
                         hf_120, hf_121, hf_122, hf_123, hf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -ff_50[k]
                  + f_0 * hf_120[k];

        t_81[k] = -ff_51[k]
                  + f_0 * hf_121[k];

        t_82[k] = -ff_52[k]
                  + f_0 * hf_122[k];

        t_83[k] = -ff_53[k]
                  + f_0 * hf_123[k];

        t_84[k] = -ff_54[k]
                  + f_0 * hf_124[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ff_55, ff_56, ff_57, ff_58, ff_59, \
                         hf_125, hf_126, hf_127, hf_128, hf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -ff_55[k]
                  + f_0 * hf_125[k];

        t_86[k] = -ff_56[k]
                  + f_0 * hf_126[k];

        t_87[k] = -ff_57[k]
                  + f_0 * hf_127[k];

        t_88[k] = -ff_58[k]
                  + f_0 * hf_128[k];

        t_89[k] = -ff_59[k]
                  + f_0 * hf_129[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, hf_130, hf_131, \
                         hf_132, hf_133, hf_134, hf_135, hf_136, \
                         hf_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * hf_130[k];

        t_91[k] = f_0 * hf_131[k];

        t_92[k] = f_0 * hf_132[k];

        t_93[k] = f_0 * hf_133[k];

        t_94[k] = f_0 * hf_134[k];

        t_95[k] = f_0 * hf_135[k];

        t_96[k] = f_0 * hf_136[k];

        t_97[k] = f_0 * hf_137[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, ff_60, ff_61, ff_62, ff_63, \
                         hf_138, hf_139, hf_150, hf_151, hf_152, \
                         hf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * hf_138[k];

        t_99[k] = f_0 * hf_139[k];

        t_100[k] = -4.0 * ff_60[k]
                   + f_0 * hf_150[k];

        t_101[k] = -4.0 * ff_61[k]
                   + f_0 * hf_151[k];

        t_102[k] = -4.0 * ff_62[k]
                   + f_0 * hf_152[k];

        t_103[k] = -4.0 * ff_63[k]
                   + f_0 * hf_153[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, ff_64, ff_65, ff_66, ff_67, ff_68, \
                         hf_154, hf_155, hf_156, hf_157, hf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -4.0 * ff_64[k]
                   + f_0 * hf_154[k];

        t_105[k] = -4.0 * ff_65[k]
                   + f_0 * hf_155[k];

        t_106[k] = -4.0 * ff_66[k]
                   + f_0 * hf_156[k];

        t_107[k] = -4.0 * ff_67[k]
                   + f_0 * hf_157[k];

        t_108[k] = -4.0 * ff_68[k]
                   + f_0 * hf_158[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, ff_69, ff_70, ff_71, ff_72, ff_73, \
                         hf_159, hf_160, hf_161, hf_162, hf_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -4.0 * ff_69[k]
                   + f_0 * hf_159[k];

        t_110[k] = -3.0 * ff_70[k]
                   + f_0 * hf_160[k];

        t_111[k] = -3.0 * ff_71[k]
                   + f_0 * hf_161[k];

        t_112[k] = -3.0 * ff_72[k]
                   + f_0 * hf_162[k];

        t_113[k] = -3.0 * ff_73[k]
                   + f_0 * hf_163[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, ff_74, ff_75, ff_76, ff_77, ff_78, \
                         hf_164, hf_165, hf_166, hf_167, hf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -3.0 * ff_74[k]
                   + f_0 * hf_164[k];

        t_115[k] = -3.0 * ff_75[k]
                   + f_0 * hf_165[k];

        t_116[k] = -3.0 * ff_76[k]
                   + f_0 * hf_166[k];

        t_117[k] = -3.0 * ff_77[k]
                   + f_0 * hf_167[k];

        t_118[k] = -3.0 * ff_78[k]
                   + f_0 * hf_168[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, ff_79, ff_80, ff_81, ff_82, ff_83, \
                         hf_169, hf_170, hf_171, hf_172, hf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -3.0 * ff_79[k]
                   + f_0 * hf_169[k];

        t_120[k] = -2.0 * ff_80[k]
                   + f_0 * hf_170[k];

        t_121[k] = -2.0 * ff_81[k]
                   + f_0 * hf_171[k];

        t_122[k] = -2.0 * ff_82[k]
                   + f_0 * hf_172[k];

        t_123[k] = -2.0 * ff_83[k]
                   + f_0 * hf_173[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, ff_84, ff_85, ff_86, ff_87, ff_88, \
                         hf_174, hf_175, hf_176, hf_177, hf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * ff_84[k]
                   + f_0 * hf_174[k];

        t_125[k] = -2.0 * ff_85[k]
                   + f_0 * hf_175[k];

        t_126[k] = -2.0 * ff_86[k]
                   + f_0 * hf_176[k];

        t_127[k] = -2.0 * ff_87[k]
                   + f_0 * hf_177[k];

        t_128[k] = -2.0 * ff_88[k]
                   + f_0 * hf_178[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, ff_89, ff_90, ff_91, ff_92, ff_93, \
                         hf_179, hf_180, hf_181, hf_182, hf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -2.0 * ff_89[k]
                   + f_0 * hf_179[k];

        t_130[k] = -ff_90[k]
                   + f_0 * hf_180[k];

        t_131[k] = -ff_91[k]
                   + f_0 * hf_181[k];

        t_132[k] = -ff_92[k]
                   + f_0 * hf_182[k];

        t_133[k] = -ff_93[k]
                   + f_0 * hf_183[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, ff_94, ff_95, ff_96, ff_97, ff_98, \
                         hf_184, hf_185, hf_186, hf_187, hf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -ff_94[k]
                   + f_0 * hf_184[k];

        t_135[k] = -ff_95[k]
                   + f_0 * hf_185[k];

        t_136[k] = -ff_96[k]
                   + f_0 * hf_186[k];

        t_137[k] = -ff_97[k]
                   + f_0 * hf_187[k];

        t_138[k] = -ff_98[k]
                   + f_0 * hf_188[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, t_144, t_145, ff_99, hf_189, \
                         hf_190, hf_191, hf_192, hf_193, hf_194, \
                         hf_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -ff_99[k]
                   + f_0 * hf_189[k];

        t_140[k] = f_0 * hf_190[k];

        t_141[k] = f_0 * hf_191[k];

        t_142[k] = f_0 * hf_192[k];

        t_143[k] = f_0 * hf_193[k];

        t_144[k] = f_0 * hf_194[k];

        t_145[k] = f_0 * hf_195[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, hf_196, hf_197, hf_198, \
                         hf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_0 * hf_196[k];

        t_147[k] = f_0 * hf_197[k];

        t_148[k] = f_0 * hf_198[k];

        t_149[k] = f_0 * hf_199[k];
    }
}

auto
compute_prim_geom_10_gf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ff, const size_t hf,
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

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_43 = buffer.data(ff + 43);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_46 = buffer.data(ff + 46);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);
    const auto *ff_54 = buffer.data(ff + 54);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_58 = buffer.data(ff + 58);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_62 = buffer.data(ff + 62);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_64 = buffer.data(ff + 64);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_70 = buffer.data(ff + 70);
    const auto *ff_71 = buffer.data(ff + 71);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_73 = buffer.data(ff + 73);
    const auto *ff_74 = buffer.data(ff + 74);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_80 = buffer.data(ff + 80);
    const auto *ff_81 = buffer.data(ff + 81);
    const auto *ff_82 = buffer.data(ff + 82);
    const auto *ff_83 = buffer.data(ff + 83);
    const auto *ff_84 = buffer.data(ff + 84);
    const auto *ff_85 = buffer.data(ff + 85);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_91 = buffer.data(ff + 91);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_94 = buffer.data(ff + 94);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_98 = buffer.data(ff + 98);
    const auto *ff_99 = buffer.data(ff + 99);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, hf_20, hf_21, hf_22, hf_23, \
                         hf_24, hf_25, hf_26, hf_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_20[k];

        t_1[k] = f_0 * hf_21[k];

        t_2[k] = f_0 * hf_22[k];

        t_3[k] = f_0 * hf_23[k];

        t_4[k] = f_0 * hf_24[k];

        t_5[k] = f_0 * hf_25[k];

        t_6[k] = f_0 * hf_26[k];

        t_7[k] = f_0 * hf_27[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, hf_28, hf_29, hf_40, \
                         hf_41, hf_42, hf_43, hf_44, hf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * hf_28[k];

        t_9[k] = f_0 * hf_29[k];

        t_10[k] = f_0 * hf_40[k];

        t_11[k] = f_0 * hf_41[k];

        t_12[k] = f_0 * hf_42[k];

        t_13[k] = f_0 * hf_43[k];

        t_14[k] = f_0 * hf_44[k];

        t_15[k] = f_0 * hf_45[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, ff_0, ff_1, hf_46, hf_47, hf_48, \
                         hf_49, hf_50, hf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * hf_46[k];

        t_17[k] = f_0 * hf_47[k];

        t_18[k] = f_0 * hf_48[k];

        t_19[k] = f_0 * hf_49[k];

        t_20[k] = -ff_0[k]
                  + f_0 * hf_50[k];

        t_21[k] = -ff_1[k]
                  + f_0 * hf_51[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ff_2, ff_3, ff_4, ff_5, ff_6, hf_52, \
                         hf_53, hf_54, hf_55, hf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -ff_2[k]
                  + f_0 * hf_52[k];

        t_23[k] = -ff_3[k]
                  + f_0 * hf_53[k];

        t_24[k] = -ff_4[k]
                  + f_0 * hf_54[k];

        t_25[k] = -ff_5[k]
                  + f_0 * hf_55[k];

        t_26[k] = -ff_6[k]
                  + f_0 * hf_56[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, t_32, ff_7, ff_8, ff_9, hf_57, hf_58, \
                         hf_59, hf_70, hf_71, hf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -ff_7[k]
                  + f_0 * hf_57[k];

        t_28[k] = -ff_8[k]
                  + f_0 * hf_58[k];

        t_29[k] = -ff_9[k]
                  + f_0 * hf_59[k];

        t_30[k] = f_0 * hf_70[k];

        t_31[k] = f_0 * hf_71[k];

        t_32[k] = f_0 * hf_72[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, t_39, hf_73, hf_74, hf_75, hf_76, \
                         hf_77, hf_78, hf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * hf_73[k];

        t_34[k] = f_0 * hf_74[k];

        t_35[k] = f_0 * hf_75[k];

        t_36[k] = f_0 * hf_76[k];

        t_37[k] = f_0 * hf_77[k];

        t_38[k] = f_0 * hf_78[k];

        t_39[k] = f_0 * hf_79[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ff_10, ff_11, ff_12, ff_13, ff_14, \
                         hf_80, hf_81, hf_82, hf_83, hf_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -ff_10[k]
                  + f_0 * hf_80[k];

        t_41[k] = -ff_11[k]
                  + f_0 * hf_81[k];

        t_42[k] = -ff_12[k]
                  + f_0 * hf_82[k];

        t_43[k] = -ff_13[k]
                  + f_0 * hf_83[k];

        t_44[k] = -ff_14[k]
                  + f_0 * hf_84[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ff_15, ff_16, ff_17, ff_18, ff_19, \
                         hf_85, hf_86, hf_87, hf_88, hf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -ff_15[k]
                  + f_0 * hf_85[k];

        t_46[k] = -ff_16[k]
                  + f_0 * hf_86[k];

        t_47[k] = -ff_17[k]
                  + f_0 * hf_87[k];

        t_48[k] = -ff_18[k]
                  + f_0 * hf_88[k];

        t_49[k] = -ff_19[k]
                  + f_0 * hf_89[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ff_20, ff_21, ff_22, ff_23, ff_24, \
                         hf_90, hf_91, hf_92, hf_93, hf_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * ff_20[k]
                  + f_0 * hf_90[k];

        t_51[k] = -2.0 * ff_21[k]
                  + f_0 * hf_91[k];

        t_52[k] = -2.0 * ff_22[k]
                  + f_0 * hf_92[k];

        t_53[k] = -2.0 * ff_23[k]
                  + f_0 * hf_93[k];

        t_54[k] = -2.0 * ff_24[k]
                  + f_0 * hf_94[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ff_25, ff_26, ff_27, ff_28, ff_29, \
                         hf_95, hf_96, hf_97, hf_98, hf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * ff_25[k]
                  + f_0 * hf_95[k];

        t_56[k] = -2.0 * ff_26[k]
                  + f_0 * hf_96[k];

        t_57[k] = -2.0 * ff_27[k]
                  + f_0 * hf_97[k];

        t_58[k] = -2.0 * ff_28[k]
                  + f_0 * hf_98[k];

        t_59[k] = -2.0 * ff_29[k]
                  + f_0 * hf_99[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, t_66, t_67, hf_110, hf_111, \
                         hf_112, hf_113, hf_114, hf_115, hf_116, \
                         hf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * hf_110[k];

        t_61[k] = f_0 * hf_111[k];

        t_62[k] = f_0 * hf_112[k];

        t_63[k] = f_0 * hf_113[k];

        t_64[k] = f_0 * hf_114[k];

        t_65[k] = f_0 * hf_115[k];

        t_66[k] = f_0 * hf_116[k];

        t_67[k] = f_0 * hf_117[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, ff_30, ff_31, ff_32, ff_33, \
                         hf_118, hf_119, hf_120, hf_121, hf_122, \
                         hf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * hf_118[k];

        t_69[k] = f_0 * hf_119[k];

        t_70[k] = -ff_30[k]
                  + f_0 * hf_120[k];

        t_71[k] = -ff_31[k]
                  + f_0 * hf_121[k];

        t_72[k] = -ff_32[k]
                  + f_0 * hf_122[k];

        t_73[k] = -ff_33[k]
                  + f_0 * hf_123[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, ff_34, ff_35, ff_36, ff_37, ff_38, \
                         hf_124, hf_125, hf_126, hf_127, hf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = -ff_34[k]
                  + f_0 * hf_124[k];

        t_75[k] = -ff_35[k]
                  + f_0 * hf_125[k];

        t_76[k] = -ff_36[k]
                  + f_0 * hf_126[k];

        t_77[k] = -ff_37[k]
                  + f_0 * hf_127[k];

        t_78[k] = -ff_38[k]
                  + f_0 * hf_128[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, ff_39, ff_40, ff_41, ff_42, ff_43, \
                         hf_129, hf_130, hf_131, hf_132, hf_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -ff_39[k]
                  + f_0 * hf_129[k];

        t_80[k] = -2.0 * ff_40[k]
                  + f_0 * hf_130[k];

        t_81[k] = -2.0 * ff_41[k]
                  + f_0 * hf_131[k];

        t_82[k] = -2.0 * ff_42[k]
                  + f_0 * hf_132[k];

        t_83[k] = -2.0 * ff_43[k]
                  + f_0 * hf_133[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ff_44, ff_45, ff_46, ff_47, ff_48, \
                         hf_134, hf_135, hf_136, hf_137, hf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -2.0 * ff_44[k]
                  + f_0 * hf_134[k];

        t_85[k] = -2.0 * ff_45[k]
                  + f_0 * hf_135[k];

        t_86[k] = -2.0 * ff_46[k]
                  + f_0 * hf_136[k];

        t_87[k] = -2.0 * ff_47[k]
                  + f_0 * hf_137[k];

        t_88[k] = -2.0 * ff_48[k]
                  + f_0 * hf_138[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ff_49, ff_50, ff_51, ff_52, ff_53, \
                         hf_139, hf_140, hf_141, hf_142, hf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -2.0 * ff_49[k]
                  + f_0 * hf_139[k];

        t_90[k] = -3.0 * ff_50[k]
                  + f_0 * hf_140[k];

        t_91[k] = -3.0 * ff_51[k]
                  + f_0 * hf_141[k];

        t_92[k] = -3.0 * ff_52[k]
                  + f_0 * hf_142[k];

        t_93[k] = -3.0 * ff_53[k]
                  + f_0 * hf_143[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ff_54, ff_55, ff_56, ff_57, ff_58, \
                         hf_144, hf_145, hf_146, hf_147, hf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -3.0 * ff_54[k]
                  + f_0 * hf_144[k];

        t_95[k] = -3.0 * ff_55[k]
                  + f_0 * hf_145[k];

        t_96[k] = -3.0 * ff_56[k]
                  + f_0 * hf_146[k];

        t_97[k] = -3.0 * ff_57[k]
                  + f_0 * hf_147[k];

        t_98[k] = -3.0 * ff_58[k]
                  + f_0 * hf_148[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, t_104, t_105, ff_59, hf_149, \
                         hf_160, hf_161, hf_162, hf_163, hf_164, \
                         hf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -3.0 * ff_59[k]
                  + f_0 * hf_149[k];

        t_100[k] = f_0 * hf_160[k];

        t_101[k] = f_0 * hf_161[k];

        t_102[k] = f_0 * hf_162[k];

        t_103[k] = f_0 * hf_163[k];

        t_104[k] = f_0 * hf_164[k];

        t_105[k] = f_0 * hf_165[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, ff_60, ff_61, hf_166, \
                         hf_167, hf_168, hf_169, hf_170, hf_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * hf_166[k];

        t_107[k] = f_0 * hf_167[k];

        t_108[k] = f_0 * hf_168[k];

        t_109[k] = f_0 * hf_169[k];

        t_110[k] = -ff_60[k]
                   + f_0 * hf_170[k];

        t_111[k] = -ff_61[k]
                   + f_0 * hf_171[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, ff_62, ff_63, ff_64, ff_65, ff_66, \
                         hf_172, hf_173, hf_174, hf_175, hf_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -ff_62[k]
                   + f_0 * hf_172[k];

        t_113[k] = -ff_63[k]
                   + f_0 * hf_173[k];

        t_114[k] = -ff_64[k]
                   + f_0 * hf_174[k];

        t_115[k] = -ff_65[k]
                   + f_0 * hf_175[k];

        t_116[k] = -ff_66[k]
                   + f_0 * hf_176[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, ff_67, ff_68, ff_69, ff_70, ff_71, \
                         hf_177, hf_178, hf_179, hf_180, hf_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -ff_67[k]
                   + f_0 * hf_177[k];

        t_118[k] = -ff_68[k]
                   + f_0 * hf_178[k];

        t_119[k] = -ff_69[k]
                   + f_0 * hf_179[k];

        t_120[k] = -2.0 * ff_70[k]
                   + f_0 * hf_180[k];

        t_121[k] = -2.0 * ff_71[k]
                   + f_0 * hf_181[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, ff_72, ff_73, ff_74, ff_75, ff_76, \
                         hf_182, hf_183, hf_184, hf_185, hf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -2.0 * ff_72[k]
                   + f_0 * hf_182[k];

        t_123[k] = -2.0 * ff_73[k]
                   + f_0 * hf_183[k];

        t_124[k] = -2.0 * ff_74[k]
                   + f_0 * hf_184[k];

        t_125[k] = -2.0 * ff_75[k]
                   + f_0 * hf_185[k];

        t_126[k] = -2.0 * ff_76[k]
                   + f_0 * hf_186[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, ff_77, ff_78, ff_79, ff_80, ff_81, \
                         hf_187, hf_188, hf_189, hf_190, hf_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -2.0 * ff_77[k]
                   + f_0 * hf_187[k];

        t_128[k] = -2.0 * ff_78[k]
                   + f_0 * hf_188[k];

        t_129[k] = -2.0 * ff_79[k]
                   + f_0 * hf_189[k];

        t_130[k] = -3.0 * ff_80[k]
                   + f_0 * hf_190[k];

        t_131[k] = -3.0 * ff_81[k]
                   + f_0 * hf_191[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, ff_82, ff_83, ff_84, ff_85, ff_86, \
                         hf_192, hf_193, hf_194, hf_195, hf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = -3.0 * ff_82[k]
                   + f_0 * hf_192[k];

        t_133[k] = -3.0 * ff_83[k]
                   + f_0 * hf_193[k];

        t_134[k] = -3.0 * ff_84[k]
                   + f_0 * hf_194[k];

        t_135[k] = -3.0 * ff_85[k]
                   + f_0 * hf_195[k];

        t_136[k] = -3.0 * ff_86[k]
                   + f_0 * hf_196[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, ff_87, ff_88, ff_89, ff_90, ff_91, \
                         hf_197, hf_198, hf_199, hf_200, hf_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -3.0 * ff_87[k]
                   + f_0 * hf_197[k];

        t_138[k] = -3.0 * ff_88[k]
                   + f_0 * hf_198[k];

        t_139[k] = -3.0 * ff_89[k]
                   + f_0 * hf_199[k];

        t_140[k] = -4.0 * ff_90[k]
                   + f_0 * hf_200[k];

        t_141[k] = -4.0 * ff_91[k]
                   + f_0 * hf_201[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, t_146, ff_92, ff_93, ff_94, ff_95, ff_96, \
                         hf_202, hf_203, hf_204, hf_205, hf_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = -4.0 * ff_92[k]
                   + f_0 * hf_202[k];

        t_143[k] = -4.0 * ff_93[k]
                   + f_0 * hf_203[k];

        t_144[k] = -4.0 * ff_94[k]
                   + f_0 * hf_204[k];

        t_145[k] = -4.0 * ff_95[k]
                   + f_0 * hf_205[k];

        t_146[k] = -4.0 * ff_96[k]
                   + f_0 * hf_206[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, ff_97, ff_98, ff_99, hf_207, hf_208, \
                         hf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -4.0 * ff_97[k]
                   + f_0 * hf_207[k];

        t_148[k] = -4.0 * ff_98[k]
                   + f_0 * hf_208[k];

        t_149[k] = -4.0 * ff_99[k]
                   + f_0 * hf_209[k];
    }
}

}  // namespace simdt2ceri
