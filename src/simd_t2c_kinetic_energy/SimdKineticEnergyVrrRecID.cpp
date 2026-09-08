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


#include "SimdKineticEnergyVrrRecID.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_id_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gd_s, const size_t gd,
                                 const size_t hp, const size_t hd, const size_t is_s,
                                 const size_t id_s, const size_t is, const size_t ip,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 4.0 * beta / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = beta / p;
    const auto f_8 = 3.0 * beta / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 2.0 * beta / p;
    const auto f_11 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_1 = buffer.data(gd_s + 1);
    const auto *gd_s_2 = buffer.data(gd_s + 2);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_23 = buffer.data(gd_s + 23);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);

    const auto *is_s_0 = buffer.data(is_s + 0);
    const auto *is_s_3 = buffer.data(is_s + 3);
    const auto *is_s_4 = buffer.data(is_s + 4);
    const auto *is_s_5 = buffer.data(is_s + 5);
    const auto *is_s_8 = buffer.data(is_s + 8);
    const auto *is_s_9 = buffer.data(is_s + 9);
    const auto *is_s_13 = buffer.data(is_s + 13);
    const auto *is_s_16 = buffer.data(is_s + 16);
    const auto *is_s_18 = buffer.data(is_s + 18);
    const auto *is_s_19 = buffer.data(is_s + 19);
    const auto *is_s_20 = buffer.data(is_s + 20);
    const auto *is_s_22 = buffer.data(is_s + 22);

    const auto *id_s_0 = buffer.data(id_s + 0);
    const auto *id_s_1 = buffer.data(id_s + 1);
    const auto *id_s_2 = buffer.data(id_s + 2);
    const auto *id_s_3 = buffer.data(id_s + 3);
    const auto *id_s_4 = buffer.data(id_s + 4);
    const auto *id_s_5 = buffer.data(id_s + 5);
    const auto *id_s_6 = buffer.data(id_s + 6);
    const auto *id_s_7 = buffer.data(id_s + 7);
    const auto *id_s_8 = buffer.data(id_s + 8);
    const auto *id_s_9 = buffer.data(id_s + 9);
    const auto *id_s_10 = buffer.data(id_s + 10);
    const auto *id_s_11 = buffer.data(id_s + 11);
    const auto *id_s_12 = buffer.data(id_s + 12);
    const auto *id_s_13 = buffer.data(id_s + 13);
    const auto *id_s_14 = buffer.data(id_s + 14);
    const auto *id_s_15 = buffer.data(id_s + 15);
    const auto *id_s_16 = buffer.data(id_s + 16);
    const auto *id_s_17 = buffer.data(id_s + 17);
    const auto *id_s_18 = buffer.data(id_s + 18);
    const auto *id_s_19 = buffer.data(id_s + 19);
    const auto *id_s_20 = buffer.data(id_s + 20);
    const auto *id_s_21 = buffer.data(id_s + 21);
    const auto *id_s_22 = buffer.data(id_s + 22);
    const auto *id_s_23 = buffer.data(id_s + 23);
    const auto *id_s_24 = buffer.data(id_s + 24);
    const auto *id_s_25 = buffer.data(id_s + 25);
    const auto *id_s_26 = buffer.data(id_s + 26);
    const auto *id_s_27 = buffer.data(id_s + 27);
    const auto *id_s_28 = buffer.data(id_s + 28);
    const auto *id_s_29 = buffer.data(id_s + 29);
    const auto *id_s_30 = buffer.data(id_s + 30);
    const auto *id_s_31 = buffer.data(id_s + 31);
    const auto *id_s_32 = buffer.data(id_s + 32);
    const auto *id_s_33 = buffer.data(id_s + 33);
    const auto *id_s_34 = buffer.data(id_s + 34);
    const auto *id_s_35 = buffer.data(id_s + 35);
    const auto *id_s_36 = buffer.data(id_s + 36);
    const auto *id_s_37 = buffer.data(id_s + 37);
    const auto *id_s_38 = buffer.data(id_s + 38);
    const auto *id_s_39 = buffer.data(id_s + 39);
    const auto *id_s_40 = buffer.data(id_s + 40);
    const auto *id_s_41 = buffer.data(id_s + 41);
    const auto *id_s_42 = buffer.data(id_s + 42);
    const auto *id_s_43 = buffer.data(id_s + 43);
    const auto *id_s_44 = buffer.data(id_s + 44);
    const auto *id_s_45 = buffer.data(id_s + 45);
    const auto *id_s_46 = buffer.data(id_s + 46);
    const auto *id_s_47 = buffer.data(id_s + 47);
    const auto *id_s_48 = buffer.data(id_s + 48);
    const auto *id_s_49 = buffer.data(id_s + 49);
    const auto *id_s_50 = buffer.data(id_s + 50);
    const auto *id_s_51 = buffer.data(id_s + 51);
    const auto *id_s_52 = buffer.data(id_s + 52);
    const auto *id_s_53 = buffer.data(id_s + 53);
    const auto *id_s_54 = buffer.data(id_s + 54);
    const auto *id_s_55 = buffer.data(id_s + 55);
    const auto *id_s_56 = buffer.data(id_s + 56);
    const auto *id_s_57 = buffer.data(id_s + 57);
    const auto *id_s_58 = buffer.data(id_s + 58);
    const auto *id_s_59 = buffer.data(id_s + 59);
    const auto *id_s_60 = buffer.data(id_s + 60);
    const auto *id_s_61 = buffer.data(id_s + 61);
    const auto *id_s_62 = buffer.data(id_s + 62);
    const auto *id_s_63 = buffer.data(id_s + 63);
    const auto *id_s_64 = buffer.data(id_s + 64);
    const auto *id_s_65 = buffer.data(id_s + 65);
    const auto *id_s_66 = buffer.data(id_s + 66);
    const auto *id_s_67 = buffer.data(id_s + 67);
    const auto *id_s_68 = buffer.data(id_s + 68);
    const auto *id_s_69 = buffer.data(id_s + 69);
    const auto *id_s_70 = buffer.data(id_s + 70);
    const auto *id_s_71 = buffer.data(id_s + 71);
    const auto *id_s_72 = buffer.data(id_s + 72);
    const auto *id_s_73 = buffer.data(id_s + 73);
    const auto *id_s_74 = buffer.data(id_s + 74);
    const auto *id_s_75 = buffer.data(id_s + 75);
    const auto *id_s_76 = buffer.data(id_s + 76);
    const auto *id_s_77 = buffer.data(id_s + 77);
    const auto *id_s_78 = buffer.data(id_s + 78);
    const auto *id_s_79 = buffer.data(id_s + 79);
    const auto *id_s_80 = buffer.data(id_s + 80);
    const auto *id_s_81 = buffer.data(id_s + 81);
    const auto *id_s_82 = buffer.data(id_s + 82);
    const auto *id_s_83 = buffer.data(id_s + 83);
    const auto *id_s_84 = buffer.data(id_s + 84);
    const auto *id_s_85 = buffer.data(id_s + 85);
    const auto *id_s_86 = buffer.data(id_s + 86);
    const auto *id_s_87 = buffer.data(id_s + 87);
    const auto *id_s_88 = buffer.data(id_s + 88);
    const auto *id_s_89 = buffer.data(id_s + 89);
    const auto *id_s_90 = buffer.data(id_s + 90);
    const auto *id_s_91 = buffer.data(id_s + 91);
    const auto *id_s_92 = buffer.data(id_s + 92);
    const auto *id_s_93 = buffer.data(id_s + 93);
    const auto *id_s_94 = buffer.data(id_s + 94);
    const auto *id_s_95 = buffer.data(id_s + 95);
    const auto *id_s_96 = buffer.data(id_s + 96);
    const auto *id_s_97 = buffer.data(id_s + 97);
    const auto *id_s_98 = buffer.data(id_s + 98);
    const auto *id_s_99 = buffer.data(id_s + 99);
    const auto *id_s_100 = buffer.data(id_s + 100);
    const auto *id_s_101 = buffer.data(id_s + 101);
    const auto *id_s_102 = buffer.data(id_s + 102);
    const auto *id_s_103 = buffer.data(id_s + 103);
    const auto *id_s_104 = buffer.data(id_s + 104);
    const auto *id_s_105 = buffer.data(id_s + 105);
    const auto *id_s_106 = buffer.data(id_s + 106);
    const auto *id_s_107 = buffer.data(id_s + 107);
    const auto *id_s_108 = buffer.data(id_s + 108);
    const auto *id_s_109 = buffer.data(id_s + 109);
    const auto *id_s_110 = buffer.data(id_s + 110);
    const auto *id_s_111 = buffer.data(id_s + 111);
    const auto *id_s_112 = buffer.data(id_s + 112);
    const auto *id_s_113 = buffer.data(id_s + 113);
    const auto *id_s_114 = buffer.data(id_s + 114);
    const auto *id_s_115 = buffer.data(id_s + 115);
    const auto *id_s_116 = buffer.data(id_s + 116);
    const auto *id_s_117 = buffer.data(id_s + 117);
    const auto *id_s_118 = buffer.data(id_s + 118);
    const auto *id_s_119 = buffer.data(id_s + 119);
    const auto *id_s_120 = buffer.data(id_s + 120);
    const auto *id_s_121 = buffer.data(id_s + 121);
    const auto *id_s_122 = buffer.data(id_s + 122);
    const auto *id_s_123 = buffer.data(id_s + 123);
    const auto *id_s_124 = buffer.data(id_s + 124);
    const auto *id_s_125 = buffer.data(id_s + 125);
    const auto *id_s_126 = buffer.data(id_s + 126);
    const auto *id_s_127 = buffer.data(id_s + 127);
    const auto *id_s_128 = buffer.data(id_s + 128);
    const auto *id_s_129 = buffer.data(id_s + 129);
    const auto *id_s_130 = buffer.data(id_s + 130);
    const auto *id_s_131 = buffer.data(id_s + 131);
    const auto *id_s_132 = buffer.data(id_s + 132);
    const auto *id_s_133 = buffer.data(id_s + 133);
    const auto *id_s_134 = buffer.data(id_s + 134);
    const auto *id_s_135 = buffer.data(id_s + 135);
    const auto *id_s_136 = buffer.data(id_s + 136);
    const auto *id_s_137 = buffer.data(id_s + 137);
    const auto *id_s_138 = buffer.data(id_s + 138);
    const auto *id_s_139 = buffer.data(id_s + 139);
    const auto *id_s_140 = buffer.data(id_s + 140);
    const auto *id_s_141 = buffer.data(id_s + 141);
    const auto *id_s_142 = buffer.data(id_s + 142);
    const auto *id_s_143 = buffer.data(id_s + 143);
    const auto *id_s_144 = buffer.data(id_s + 144);
    const auto *id_s_145 = buffer.data(id_s + 145);
    const auto *id_s_146 = buffer.data(id_s + 146);
    const auto *id_s_147 = buffer.data(id_s + 147);
    const auto *id_s_148 = buffer.data(id_s + 148);
    const auto *id_s_149 = buffer.data(id_s + 149);
    const auto *id_s_150 = buffer.data(id_s + 150);
    const auto *id_s_151 = buffer.data(id_s + 151);
    const auto *id_s_152 = buffer.data(id_s + 152);
    const auto *id_s_153 = buffer.data(id_s + 153);
    const auto *id_s_154 = buffer.data(id_s + 154);
    const auto *id_s_155 = buffer.data(id_s + 155);
    const auto *id_s_156 = buffer.data(id_s + 156);
    const auto *id_s_157 = buffer.data(id_s + 157);
    const auto *id_s_158 = buffer.data(id_s + 158);
    const auto *id_s_159 = buffer.data(id_s + 159);
    const auto *id_s_160 = buffer.data(id_s + 160);
    const auto *id_s_161 = buffer.data(id_s + 161);
    const auto *id_s_162 = buffer.data(id_s + 162);
    const auto *id_s_163 = buffer.data(id_s + 163);
    const auto *id_s_164 = buffer.data(id_s + 164);
    const auto *id_s_165 = buffer.data(id_s + 165);
    const auto *id_s_166 = buffer.data(id_s + 166);
    const auto *id_s_167 = buffer.data(id_s + 167);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_22 = buffer.data(is + 22);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_39 = buffer.data(ip + 39);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);
    const auto *ip_42 = buffer.data(ip + 42);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_45 = buffer.data(ip + 45);
    const auto *ip_46 = buffer.data(ip + 46);
    const auto *ip_47 = buffer.data(ip + 47);
    const auto *ip_48 = buffer.data(ip + 48);
    const auto *ip_49 = buffer.data(ip + 49);
    const auto *ip_50 = buffer.data(ip + 50);
    const auto *ip_51 = buffer.data(ip + 51);
    const auto *ip_52 = buffer.data(ip + 52);
    const auto *ip_53 = buffer.data(ip + 53);
    const auto *ip_54 = buffer.data(ip + 54);
    const auto *ip_55 = buffer.data(ip + 55);
    const auto *ip_56 = buffer.data(ip + 56);
    const auto *ip_57 = buffer.data(ip + 57);
    const auto *ip_58 = buffer.data(ip + 58);
    const auto *ip_59 = buffer.data(ip + 59);
    const auto *ip_60 = buffer.data(ip + 60);
    const auto *ip_61 = buffer.data(ip + 61);
    const auto *ip_62 = buffer.data(ip + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, hp_0, is_s_0, id_s_0, id_s_1, \
                         id_s_2, id_s_3, is_0, ip_0, ip_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 - f_1 * is_s_0[k]
                 + f_2 * id_s_0[k]
                 + f_3 * is_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_2 * id_s_1[k]
                 + pb_y[k] * ip_0[k];

        t_2[k] = f_2 * id_s_2[k]
                 + pb_z[k] * ip_0[k];

        t_3[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_3[k]
                 + f_3 * is_0[k]
                 + pb_y[k] * ip_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_y, pb_y, pb_z, hd_0, is_s_0, id_s_4, id_s_5, \
                         id_s_6, is_0, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * id_s_4[k]
                 + pb_y[k] * ip_2[k];

        t_5[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_5[k]
                 + f_3 * is_0[k]
                 + pb_z[k] * ip_2[k];

        t_6[k] = pa_y[k] * hd_0[k]
                 + f_2 * id_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, pb_z, gd_s_2, gd_2, hp_3, hd_5, id_s_7, \
                         id_s_8, id_s_9, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_4 * hp_3[k]
                 + f_2 * id_s_7[k]
                 + pb_x[k] * ip_4[k];

        t_8[k] = f_2 * id_s_8[k]
                 + pb_z[k] * ip_3[k];

        t_9[k] = -f_5 * gd_s_2[k]
                 + f_6 * gd_2[k]
                 + pa_x[k] * hd_5[k]
                 + f_2 * id_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_y, pb_z, hd_0, hd_2, id_s_10, \
                         id_s_11, id_s_12, id_s_13, ip_4, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * id_s_10[k]
                  + pb_z[k] * ip_4[k];

        t_11[k] = pa_y[k] * hd_2[k]
                  + f_2 * id_s_11[k];

        t_12[k] = pa_z[k] * hd_0[k]
                  + f_2 * id_s_12[k];

        t_13[k] = f_2 * id_s_13[k]
                  + pb_y[k] * ip_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_x, pb_y, hp_4, hd_1, id_s_14, id_s_15, \
                         id_s_16, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * hp_4[k]
                  + f_2 * id_s_14[k]
                  + pb_x[k] * ip_6[k];

        t_15[k] = pa_z[k] * hd_1[k]
                  + f_2 * id_s_15[k];

        t_16[k] = f_2 * id_s_16[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pa_y, gd_s_0, gd_s_4, gd_0, gd_4, hd_3, hd_8, \
                         id_s_17, id_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_5 * gd_s_4[k]
                  + f_6 * gd_4[k]
                  + pa_x[k] * hd_8[k]
                  + f_2 * id_s_17[k];

        t_18[k] = -f_7 * gd_s_0[k]
                  + f_3 * gd_0[k]
                  + pa_y[k] * hd_3[k]
                  + f_2 * id_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pb_x, pb_z, gd_s_6, gd_6, hp_5, hd_11, \
                         id_s_19, id_s_20, id_s_21, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_6 * hp_5[k]
                  + f_2 * id_s_19[k]
                  + pb_x[k] * ip_8[k];

        t_20[k] = f_2 * id_s_20[k]
                  + pb_z[k] * ip_7[k];

        t_21[k] = -f_8 * gd_s_6[k]
                  + f_9 * gd_6[k]
                  + pa_x[k] * hd_11[k]
                  + f_2 * id_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pb_z, hd_6, is_s_3, id_s_22, id_s_23, \
                         id_s_24, is_3, ip_8, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * id_s_22[k]
                  + pb_z[k] * ip_8[k];

        t_23[k] = -f_1 * is_s_3[k]
                  + f_2 * id_s_23[k]
                  + f_3 * is_3[k]
                  + pb_z[k] * ip_9[k];

        t_24[k] = pa_y[k] * hd_6[k]
                  + f_2 * id_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pa_z, pb_y, hp_4, hd_4, hd_5, hd_7, \
                         id_s_25, id_s_26, id_s_27, id_s_28, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * hd_4[k]
                  + f_2 * id_s_25[k];

        t_26[k] = pa_y[k] * hd_7[k]
                  + f_2 * id_s_26[k];

        t_27[k] = pa_z[k] * hd_5[k]
                  + f_2 * id_s_27[k];

        t_28[k] = f_3 * hp_4[k]
                  + f_2 * id_s_28[k]
                  + pb_y[k] * ip_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_y, gd_s_0, gd_0, hd_6, hd_8, \
                         id_s_29, id_s_30, id_s_31, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * hd_8[k]
                  + f_2 * id_s_29[k];

        t_30[k] = -f_7 * gd_s_0[k]
                  + f_3 * gd_0[k]
                  + pa_z[k] * hd_6[k]
                  + f_2 * id_s_30[k];

        t_31[k] = f_2 * id_s_31[k]
                  + pb_y[k] * ip_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, hp_7, is_s_4, id_s_32, id_s_33, \
                         id_s_34, is_4, ip_12, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_6 * hp_7[k]
                  + f_2 * id_s_32[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = -f_1 * is_s_4[k]
                  + f_2 * id_s_33[k]
                  + f_3 * is_4[k]
                  + pb_y[k] * ip_12[k];

        t_34[k] = f_2 * id_s_34[k]
                  + pb_y[k] * ip_13[k];
    }

#pragma omp simd aligned(t_35, t_36, pa_x, pa_y, gd_s_1, gd_s_8, gd_1, gd_8, hd_9, hd_14, \
                         id_s_35, id_s_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_8 * gd_s_8[k]
                  + f_9 * gd_8[k]
                  + pa_x[k] * hd_14[k]
                  + f_2 * id_s_35[k];

        t_36[k] = -f_10 * gd_s_1[k]
                  + f_11 * gd_1[k]
                  + pa_y[k] * hd_9[k]
                  + f_2 * id_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_x, pb_x, pb_z, gd_s_9, gd_9, hp_8, hd_17, \
                         id_s_37, id_s_38, id_s_39, ip_14, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * hp_8[k]
                  + f_2 * id_s_37[k]
                  + pb_x[k] * ip_15[k];

        t_38[k] = f_2 * id_s_38[k]
                  + pb_z[k] * ip_14[k];

        t_39[k] = -f_10 * gd_s_9[k]
                  + f_11 * gd_9[k]
                  + pa_x[k] * hd_17[k]
                  + f_2 * id_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_z, hd_9, hd_10, is_s_5, id_s_40, \
                         id_s_41, id_s_42, id_s_43, is_5, ip_15, \
                         ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * id_s_40[k]
                  + pb_z[k] * ip_15[k];

        t_41[k] = -f_1 * is_s_5[k]
                  + f_2 * id_s_41[k]
                  + f_3 * is_5[k]
                  + pb_z[k] * ip_16[k];

        t_42[k] = pa_z[k] * hd_9[k]
                  + f_2 * id_s_42[k];

        t_43[k] = pa_z[k] * hd_10[k]
                  + f_2 * id_s_43[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_z, pb_x, pb_y, hp_6, hp_9, hd_11, id_s_44, \
                         id_s_45, id_s_46, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * hp_9[k]
                  + f_2 * id_s_44[k]
                  + pb_x[k] * ip_17[k];

        t_45[k] = pa_z[k] * hd_11[k]
                  + f_2 * id_s_45[k];

        t_46[k] = f_11 * hp_6[k]
                  + f_2 * id_s_46[k]
                  + pb_y[k] * ip_17[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pa_y, pb_x, gd_s_10, gd_10, hp_10, hd_12, \
                         hd_18, id_s_47, id_s_48, id_s_49, ip_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -f_10 * gd_s_10[k]
                  + f_11 * gd_10[k]
                  + pa_x[k] * hd_18[k]
                  + f_2 * id_s_47[k];

        t_48[k] = pa_y[k] * hd_12[k]
                  + f_2 * id_s_48[k];

        t_49[k] = f_9 * hp_10[k]
                  + f_2 * id_s_49[k]
                  + pb_x[k] * ip_18[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pa_y, pb_y, gd_s_11, gd_11, hp_7, hd_13, \
                         hd_20, id_s_50, id_s_51, id_s_52, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * hd_13[k]
                  + f_2 * id_s_50[k];

        t_51[k] = -f_10 * gd_s_11[k]
                  + f_11 * gd_11[k]
                  + pa_x[k] * hd_20[k]
                  + f_2 * id_s_51[k];

        t_52[k] = f_3 * hp_7[k]
                  + f_2 * id_s_52[k]
                  + pb_y[k] * ip_19[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pa_z, pb_y, gd_s_3, gd_3, hd_12, hd_14, \
                         id_s_53, id_s_54, id_s_55, ip_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_y[k] * hd_14[k]
                  + f_2 * id_s_53[k];

        t_54[k] = -f_10 * gd_s_3[k]
                  + f_11 * gd_3[k]
                  + pa_z[k] * hd_12[k]
                  + f_2 * id_s_54[k];

        t_55[k] = f_2 * id_s_55[k]
                  + pb_y[k] * ip_20[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pb_x, pb_y, hp_12, is_s_8, id_s_56, id_s_57, \
                         id_s_58, is_8, ip_21, ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_9 * hp_12[k]
                  + f_2 * id_s_56[k]
                  + pb_x[k] * ip_22[k];

        t_57[k] = -f_1 * is_s_8[k]
                  + f_2 * id_s_57[k]
                  + f_3 * is_8[k]
                  + pb_y[k] * ip_21[k];

        t_58[k] = f_2 * id_s_58[k]
                  + pb_y[k] * ip_22[k];
    }

#pragma omp simd aligned(t_59, t_60, pa_x, pa_y, gd_s_5, gd_s_12, gd_5, gd_12, hd_15, hd_23, \
                         id_s_59, id_s_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_10 * gd_s_12[k]
                  + f_11 * gd_12[k]
                  + pa_x[k] * hd_23[k]
                  + f_2 * id_s_59[k];

        t_60[k] = -f_8 * gd_s_5[k]
                  + f_9 * gd_5[k]
                  + pa_y[k] * hd_15[k]
                  + f_2 * id_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_x, pb_x, pb_z, gd_s_14, gd_14, hp_13, hd_26, \
                         id_s_61, id_s_62, id_s_63, ip_23, ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_11 * hp_13[k]
                  + f_2 * id_s_61[k]
                  + pb_x[k] * ip_24[k];

        t_62[k] = f_2 * id_s_62[k]
                  + pb_z[k] * ip_23[k];

        t_63[k] = -f_7 * gd_s_14[k]
                  + f_3 * gd_14[k]
                  + pa_x[k] * hd_26[k]
                  + f_2 * id_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_z, hd_15, hd_16, is_s_9, id_s_64, \
                         id_s_65, id_s_66, id_s_67, is_9, ip_24, \
                         ip_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * id_s_64[k]
                  + pb_z[k] * ip_24[k];

        t_65[k] = -f_1 * is_s_9[k]
                  + f_2 * id_s_65[k]
                  + f_3 * is_9[k]
                  + pb_z[k] * ip_25[k];

        t_66[k] = pa_z[k] * hd_15[k]
                  + f_2 * id_s_66[k];

        t_67[k] = pa_z[k] * hd_16[k]
                  + f_2 * id_s_67[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_z, pb_x, pb_y, hp_9, hp_14, hd_17, id_s_68, \
                         id_s_69, id_s_70, ip_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_11 * hp_14[k]
                  + f_2 * id_s_68[k]
                  + pb_x[k] * ip_26[k];

        t_69[k] = pa_z[k] * hd_17[k]
                  + f_2 * id_s_69[k];

        t_70[k] = f_9 * hp_9[k]
                  + f_2 * id_s_70[k]
                  + pb_y[k] * ip_26[k];
    }

#pragma omp simd aligned(t_71, t_72, pa_x, pa_y, gd_s_7, gd_s_16, gd_7, gd_16, hd_19, hd_27, \
                         id_s_71, id_s_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -f_7 * gd_s_16[k]
                  + f_3 * gd_16[k]
                  + pa_x[k] * hd_27[k]
                  + f_2 * id_s_71[k];

        t_72[k] = -f_7 * gd_s_7[k]
                  + f_3 * gd_7[k]
                  + pa_y[k] * hd_19[k]
                  + f_2 * id_s_72[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pb_x, gd_s_17, gd_17, hp_15, hp_16, hd_28, \
                         id_s_73, id_s_74, id_s_75, ip_27, ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_11 * hp_15[k]
                  + f_2 * id_s_73[k]
                  + pb_x[k] * ip_27[k];

        t_74[k] = f_11 * hp_16[k]
                  + f_2 * id_s_74[k]
                  + pb_x[k] * ip_28[k];

        t_75[k] = -f_7 * gd_s_17[k]
                  + f_3 * gd_17[k]
                  + pa_x[k] * hd_28[k]
                  + f_2 * id_s_75[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_x, pa_y, pb_y, gd_s_18, gd_18, hp_11, hd_21, \
                         hd_29, id_s_76, id_s_77, id_s_78, ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_11 * hp_11[k]
                  + f_2 * id_s_76[k]
                  + pb_y[k] * ip_28[k];

        t_77[k] = -f_7 * gd_s_18[k]
                  + f_3 * gd_18[k]
                  + pa_x[k] * hd_29[k]
                  + f_2 * id_s_77[k];

        t_78[k] = pa_y[k] * hd_21[k]
                  + f_2 * id_s_78[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pa_x, pa_y, pb_x, gd_s_19, gd_19, hp_17, hd_22, \
                         hd_30, id_s_79, id_s_80, id_s_81, ip_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * hp_17[k]
                  + f_2 * id_s_79[k]
                  + pb_x[k] * ip_29[k];

        t_80[k] = pa_y[k] * hd_22[k]
                  + f_2 * id_s_80[k];

        t_81[k] = -f_7 * gd_s_19[k]
                  + f_3 * gd_19[k]
                  + pa_x[k] * hd_30[k]
                  + f_2 * id_s_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pa_z, pb_y, gd_s_7, gd_7, hp_12, hd_21, \
                         hd_23, id_s_82, id_s_83, id_s_84, ip_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * hp_12[k]
                  + f_2 * id_s_82[k]
                  + pb_y[k] * ip_30[k];

        t_83[k] = pa_y[k] * hd_23[k]
                  + f_2 * id_s_83[k];

        t_84[k] = -f_8 * gd_s_7[k]
                  + f_9 * gd_7[k]
                  + pa_z[k] * hd_21[k]
                  + f_2 * id_s_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, pb_y, hp_18, is_s_13, id_s_85, id_s_86, \
                         id_s_87, id_s_88, is_13, ip_31, ip_32, ip_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_2 * id_s_85[k]
                  + pb_y[k] * ip_31[k];

        t_86[k] = f_11 * hp_18[k]
                  + f_2 * id_s_86[k]
                  + pb_x[k] * ip_33[k];

        t_87[k] = -f_1 * is_s_13[k]
                  + f_2 * id_s_87[k]
                  + f_3 * is_13[k]
                  + pb_y[k] * ip_32[k];

        t_88[k] = f_2 * id_s_88[k]
                  + pb_y[k] * ip_33[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pb_x, gd_s_23, gd_23, hp_19, hp_20, hd_33, \
                         hd_34, id_s_89, id_s_90, id_s_91, ip_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -f_7 * gd_s_23[k]
                  + f_3 * gd_23[k]
                  + pa_x[k] * hd_33[k]
                  + f_2 * id_s_89[k];

        t_90[k] = f_11 * hp_19[k]
                  + pa_x[k] * hd_34[k]
                  + f_2 * id_s_90[k];

        t_91[k] = f_3 * hp_20[k]
                  + f_2 * id_s_91[k]
                  + pb_x[k] * ip_35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pb_z, hd_35, hd_36, id_s_92, id_s_93, \
                         id_s_94, id_s_95, ip_34, ip_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_2 * id_s_92[k]
                  + pb_z[k] * ip_34[k];

        t_93[k] = pa_x[k] * hd_35[k]
                  + f_2 * id_s_93[k];

        t_94[k] = f_2 * id_s_94[k]
                  + pb_z[k] * ip_35[k];

        t_95[k] = pa_x[k] * hd_36[k]
                  + f_2 * id_s_95[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pa_z, pb_x, hp_22, hd_24, hd_25, hd_37, \
                         id_s_96, id_s_97, id_s_98, id_s_99, ip_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * hd_24[k]
                  + f_2 * id_s_96[k];

        t_97[k] = pa_z[k] * hd_25[k]
                  + f_2 * id_s_97[k];

        t_98[k] = f_3 * hp_22[k]
                  + f_2 * id_s_98[k]
                  + pb_x[k] * ip_36[k];

        t_99[k] = pa_x[k] * hd_37[k]
                  + f_2 * id_s_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_x, hp_23, hp_24, hd_38, hd_39, \
                         hd_40, id_s_100, id_s_101, id_s_102, id_s_103, \
                         ip_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * hd_38[k]
                   + f_2 * id_s_100[k];

        t_101[k] = pa_x[k] * hd_39[k]
                   + f_2 * id_s_101[k];

        t_102[k] = f_11 * hp_23[k]
                   + pa_x[k] * hd_40[k]
                   + f_2 * id_s_102[k];

        t_103[k] = f_3 * hp_24[k]
                   + f_2 * id_s_103[k]
                   + pb_x[k] * ip_37[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pb_x, hp_25, hd_41, hd_42, hd_43, \
                         id_s_104, id_s_105, id_s_106, id_s_107, \
                         ip_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * hp_25[k]
                   + f_2 * id_s_104[k]
                   + pb_x[k] * ip_38[k];

        t_105[k] = pa_x[k] * hd_41[k]
                   + f_2 * id_s_105[k];

        t_106[k] = pa_x[k] * hd_42[k]
                   + f_2 * id_s_106[k];

        t_107[k] = pa_x[k] * hd_43[k]
                   + f_2 * id_s_107[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_x, pb_x, hp_26, hp_27, hp_28, hd_44, \
                         id_s_108, id_s_109, id_s_110, ip_39, ip_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_11 * hp_26[k]
                   + pa_x[k] * hd_44[k]
                   + f_2 * id_s_108[k];

        t_109[k] = f_3 * hp_27[k]
                   + f_2 * id_s_109[k]
                   + pb_x[k] * ip_39[k];

        t_110[k] = f_3 * hp_28[k]
                   + f_2 * id_s_110[k]
                   + pb_x[k] * ip_40[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pa_y, hd_31, hd_45, hd_46, hd_47, \
                         id_s_111, id_s_112, id_s_113, id_s_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pa_x[k] * hd_45[k]
                   + f_2 * id_s_111[k];

        t_112[k] = pa_x[k] * hd_46[k]
                   + f_2 * id_s_112[k];

        t_113[k] = pa_x[k] * hd_47[k]
                   + f_2 * id_s_113[k];

        t_114[k] = pa_y[k] * hd_31[k]
                   + f_2 * id_s_114[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pa_x, pa_y, pb_x, hp_29, hd_32, hd_48, \
                         hd_49, id_s_115, id_s_116, id_s_117, id_s_118, \
                         ip_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_3 * hp_29[k]
                   + f_2 * id_s_115[k]
                   + pb_x[k] * ip_41[k];

        t_116[k] = pa_y[k] * hd_32[k]
                   + f_2 * id_s_116[k];

        t_117[k] = pa_x[k] * hd_48[k]
                   + f_2 * id_s_117[k];

        t_118[k] = pa_x[k] * hd_49[k]
                   + f_2 * id_s_118[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pa_x, pb_y, hp_31, hd_50, hd_51, id_s_119, \
                         id_s_120, id_s_121, ip_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pa_x[k] * hd_50[k]
                   + f_2 * id_s_119[k];

        t_120[k] = f_11 * hp_31[k]
                   + pa_x[k] * hd_51[k]
                   + f_2 * id_s_120[k];

        t_121[k] = f_2 * id_s_121[k]
                   + pb_y[k] * ip_42[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_x, pb_x, pb_y, hp_33, hd_52, hd_53, \
                         id_s_122, id_s_123, id_s_124, id_s_125, \
                         ip_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_3 * hp_33[k]
                   + f_2 * id_s_122[k]
                   + pb_x[k] * ip_43[k];

        t_123[k] = pa_x[k] * hd_52[k]
                   + f_2 * id_s_123[k];

        t_124[k] = f_2 * id_s_124[k]
                   + pb_y[k] * ip_43[k];

        t_125[k] = pa_x[k] * hd_53[k]
                   + f_2 * id_s_125[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pb_x, pb_y, hp_20, is_s_16, id_s_126, \
                         id_s_127, id_s_128, id_s_129, is_16, ip_44, ip_45, \
                         ip_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -f_1 * is_s_16[k]
                   + f_2 * id_s_126[k]
                   + f_3 * is_16[k]
                   + pb_x[k] * ip_44[k];

        t_127[k] = f_2 * id_s_127[k]
                   + pb_x[k] * ip_45[k];

        t_128[k] = f_2 * id_s_128[k]
                   + pb_x[k] * ip_46[k];

        t_129[k] = f_0 * hp_20[k]
                   - f_1 * is_s_16[k]
                   + f_2 * id_s_129[k]
                   + f_3 * is_16[k]
                   + pb_y[k] * ip_45[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pa_z, pb_z, hd_34, is_s_16, id_s_130, id_s_131, \
                         id_s_132, is_16, ip_45, ip_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_2 * id_s_130[k]
                   + pb_z[k] * ip_45[k];

        t_131[k] = -f_1 * is_s_16[k]
                   + f_2 * id_s_131[k]
                   + f_3 * is_16[k]
                   + pb_z[k] * ip_46[k];

        t_132[k] = pa_z[k] * hd_34[k]
                   + f_2 * id_s_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_z, pb_x, pb_y, hp_22, hd_35, id_s_133, \
                         id_s_134, id_s_135, id_s_136, ip_47, ip_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_2 * id_s_133[k]
                   + pb_x[k] * ip_47[k];

        t_134[k] = f_2 * id_s_134[k]
                   + pb_x[k] * ip_48[k];

        t_135[k] = pa_z[k] * hd_35[k]
                   + f_2 * id_s_135[k];

        t_136[k] = f_4 * hp_22[k]
                   + f_2 * id_s_136[k]
                   + pb_y[k] * ip_48[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pb_x, gd_s_16, gd_16, hd_39, is_s_18, \
                         id_s_137, id_s_138, id_s_139, is_18, ip_49, \
                         ip_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -f_5 * gd_s_16[k]
                   + f_6 * gd_16[k]
                   + pa_y[k] * hd_39[k]
                   + f_2 * id_s_137[k];

        t_138[k] = -f_1 * is_s_18[k]
                   + f_2 * id_s_138[k]
                   + f_3 * is_18[k]
                   + pb_x[k] * ip_49[k];

        t_139[k] = f_2 * id_s_139[k]
                   + pb_x[k] * ip_50[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_z, pb_x, pb_y, gd_s_14, gd_14, hp_25, hd_37, \
                         id_s_140, id_s_141, id_s_142, ip_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_2 * id_s_140[k]
                   + pb_x[k] * ip_51[k];

        t_141[k] = -f_7 * gd_s_14[k]
                   + f_3 * gd_14[k]
                   + pa_z[k] * hd_37[k]
                   + f_2 * id_s_141[k];

        t_142[k] = f_6 * hp_25[k]
                   + f_2 * id_s_142[k]
                   + pb_y[k] * ip_51[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_y, pb_x, gd_s_18, gd_18, hd_43, is_s_19, \
                         id_s_143, id_s_144, id_s_145, is_19, ip_52, \
                         ip_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = -f_8 * gd_s_18[k]
                   + f_9 * gd_18[k]
                   + pa_y[k] * hd_43[k]
                   + f_2 * id_s_143[k];

        t_144[k] = -f_1 * is_s_19[k]
                   + f_2 * id_s_144[k]
                   + f_3 * is_19[k]
                   + pb_x[k] * ip_52[k];

        t_145[k] = f_2 * id_s_145[k]
                   + pb_x[k] * ip_53[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_z, pb_x, pb_y, gd_s_15, gd_15, hp_28, hd_41, \
                         id_s_146, id_s_147, id_s_148, ip_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_2 * id_s_146[k]
                   + pb_x[k] * ip_54[k];

        t_147[k] = -f_10 * gd_s_15[k]
                   + f_11 * gd_15[k]
                   + pa_z[k] * hd_41[k]
                   + f_2 * id_s_147[k];

        t_148[k] = f_9 * hp_28[k]
                   + f_2 * id_s_148[k]
                   + pb_y[k] * ip_54[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pa_y, pb_x, gd_s_20, gd_20, hd_47, is_s_20, \
                         id_s_149, id_s_150, id_s_151, is_20, ip_55, \
                         ip_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -f_10 * gd_s_20[k]
                   + f_11 * gd_20[k]
                   + pa_y[k] * hd_47[k]
                   + f_2 * id_s_149[k];

        t_150[k] = -f_1 * is_s_20[k]
                   + f_2 * id_s_150[k]
                   + f_3 * is_20[k]
                   + pb_x[k] * ip_55[k];

        t_151[k] = f_2 * id_s_151[k]
                   + pb_x[k] * ip_56[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pa_z, pb_x, pb_y, gd_s_17, gd_17, hp_30, hd_45, \
                         id_s_152, id_s_153, id_s_154, ip_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_2 * id_s_152[k]
                   + pb_x[k] * ip_57[k];

        t_153[k] = -f_8 * gd_s_17[k]
                   + f_9 * gd_17[k]
                   + pa_z[k] * hd_45[k]
                   + f_2 * id_s_153[k];

        t_154[k] = f_11 * hp_30[k]
                   + f_2 * id_s_154[k]
                   + pb_y[k] * ip_57[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_y, pb_x, gd_s_23, gd_23, hd_50, hd_51, \
                         id_s_155, id_s_156, id_s_157, id_s_158, ip_58, \
                         ip_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = -f_7 * gd_s_23[k]
                   + f_3 * gd_23[k]
                   + pa_y[k] * hd_50[k]
                   + f_2 * id_s_155[k];

        t_156[k] = pa_y[k] * hd_51[k]
                   + f_2 * id_s_156[k];

        t_157[k] = f_2 * id_s_157[k]
                   + pb_x[k] * ip_58[k];

        t_158[k] = f_2 * id_s_158[k]
                   + pb_x[k] * ip_59[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_y, pb_y, hp_32, hp_33, hd_52, hd_53, \
                         id_s_159, id_s_160, id_s_161, ip_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_11 * hp_32[k]
                   + pa_y[k] * hd_52[k]
                   + f_2 * id_s_159[k];

        t_160[k] = f_3 * hp_33[k]
                   + f_2 * id_s_160[k]
                   + pb_y[k] * ip_59[k];

        t_161[k] = pa_y[k] * hd_53[k]
                   + f_2 * id_s_161[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_x, pb_y, is_s_22, id_s_162, id_s_163, \
                         id_s_164, id_s_165, is_22, ip_60, ip_61, \
                         ip_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -f_1 * is_s_22[k]
                   + f_2 * id_s_162[k]
                   + f_3 * is_22[k]
                   + pb_x[k] * ip_60[k];

        t_163[k] = f_2 * id_s_163[k]
                   + pb_x[k] * ip_61[k];

        t_164[k] = f_2 * id_s_164[k]
                   + pb_x[k] * ip_62[k];

        t_165[k] = -f_1 * is_s_22[k]
                   + f_2 * id_s_165[k]
                   + f_3 * is_22[k]
                   + pb_y[k] * ip_61[k];
    }

#pragma omp simd aligned(t_166, t_167, pb_y, pb_z, hp_33, is_s_22, id_s_166, id_s_167, is_22, \
                         ip_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_2 * id_s_166[k]
                   + pb_y[k] * ip_62[k];

        t_167[k] = f_0 * hp_33[k]
                   - f_1 * is_s_22[k]
                   + f_2 * id_s_167[k]
                   + f_3 * is_22[k]
                   + pb_z[k] * ip_62[k];
    }
}

auto
compute_prim_id_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gd_s, const size_t gd,
                                 const size_t hp, const size_t hd, const size_t is_s,
                                 const size_t id_s, const size_t is, const size_t ip,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 4.0 * beta / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = beta / p;
    const auto f_8 = 3.0 * beta / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 2.0 * beta / p;
    const auto f_11 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_22 = buffer.data(gd_s + 22);
    const auto *gd_s_24 = buffer.data(gd_s + 24);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_27 = buffer.data(gd_s + 27);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_32 = buffer.data(gd_s + 32);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_32 = buffer.data(gd + 32);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_69 = buffer.data(hd + 69);

    const auto *is_s_0 = buffer.data(is_s + 0);
    const auto *is_s_3 = buffer.data(is_s + 3);
    const auto *is_s_4 = buffer.data(is_s + 4);
    const auto *is_s_5 = buffer.data(is_s + 5);
    const auto *is_s_7 = buffer.data(is_s + 7);
    const auto *is_s_8 = buffer.data(is_s + 8);
    const auto *is_s_11 = buffer.data(is_s + 11);
    const auto *is_s_12 = buffer.data(is_s + 12);
    const auto *is_s_14 = buffer.data(is_s + 14);
    const auto *is_s_15 = buffer.data(is_s + 15);
    const auto *is_s_16 = buffer.data(is_s + 16);
    const auto *is_s_18 = buffer.data(is_s + 18);

    const auto *id_s_0 = buffer.data(id_s + 0);
    const auto *id_s_1 = buffer.data(id_s + 1);
    const auto *id_s_2 = buffer.data(id_s + 2);
    const auto *id_s_3 = buffer.data(id_s + 3);
    const auto *id_s_4 = buffer.data(id_s + 4);
    const auto *id_s_5 = buffer.data(id_s + 5);
    const auto *id_s_6 = buffer.data(id_s + 6);
    const auto *id_s_7 = buffer.data(id_s + 7);
    const auto *id_s_8 = buffer.data(id_s + 8);
    const auto *id_s_9 = buffer.data(id_s + 9);
    const auto *id_s_10 = buffer.data(id_s + 10);
    const auto *id_s_11 = buffer.data(id_s + 11);
    const auto *id_s_12 = buffer.data(id_s + 12);
    const auto *id_s_13 = buffer.data(id_s + 13);
    const auto *id_s_14 = buffer.data(id_s + 14);
    const auto *id_s_15 = buffer.data(id_s + 15);
    const auto *id_s_16 = buffer.data(id_s + 16);
    const auto *id_s_17 = buffer.data(id_s + 17);
    const auto *id_s_18 = buffer.data(id_s + 18);
    const auto *id_s_19 = buffer.data(id_s + 19);
    const auto *id_s_20 = buffer.data(id_s + 20);
    const auto *id_s_21 = buffer.data(id_s + 21);
    const auto *id_s_22 = buffer.data(id_s + 22);
    const auto *id_s_23 = buffer.data(id_s + 23);
    const auto *id_s_24 = buffer.data(id_s + 24);
    const auto *id_s_25 = buffer.data(id_s + 25);
    const auto *id_s_26 = buffer.data(id_s + 26);
    const auto *id_s_27 = buffer.data(id_s + 27);
    const auto *id_s_28 = buffer.data(id_s + 28);
    const auto *id_s_29 = buffer.data(id_s + 29);
    const auto *id_s_30 = buffer.data(id_s + 30);
    const auto *id_s_31 = buffer.data(id_s + 31);
    const auto *id_s_32 = buffer.data(id_s + 32);
    const auto *id_s_33 = buffer.data(id_s + 33);
    const auto *id_s_34 = buffer.data(id_s + 34);
    const auto *id_s_35 = buffer.data(id_s + 35);
    const auto *id_s_36 = buffer.data(id_s + 36);
    const auto *id_s_37 = buffer.data(id_s + 37);
    const auto *id_s_38 = buffer.data(id_s + 38);
    const auto *id_s_39 = buffer.data(id_s + 39);
    const auto *id_s_40 = buffer.data(id_s + 40);
    const auto *id_s_41 = buffer.data(id_s + 41);
    const auto *id_s_42 = buffer.data(id_s + 42);
    const auto *id_s_43 = buffer.data(id_s + 43);
    const auto *id_s_44 = buffer.data(id_s + 44);
    const auto *id_s_45 = buffer.data(id_s + 45);
    const auto *id_s_46 = buffer.data(id_s + 46);
    const auto *id_s_47 = buffer.data(id_s + 47);
    const auto *id_s_48 = buffer.data(id_s + 48);
    const auto *id_s_49 = buffer.data(id_s + 49);
    const auto *id_s_50 = buffer.data(id_s + 50);
    const auto *id_s_51 = buffer.data(id_s + 51);
    const auto *id_s_52 = buffer.data(id_s + 52);
    const auto *id_s_53 = buffer.data(id_s + 53);
    const auto *id_s_54 = buffer.data(id_s + 54);
    const auto *id_s_55 = buffer.data(id_s + 55);
    const auto *id_s_56 = buffer.data(id_s + 56);
    const auto *id_s_57 = buffer.data(id_s + 57);
    const auto *id_s_58 = buffer.data(id_s + 58);
    const auto *id_s_59 = buffer.data(id_s + 59);
    const auto *id_s_60 = buffer.data(id_s + 60);
    const auto *id_s_61 = buffer.data(id_s + 61);
    const auto *id_s_62 = buffer.data(id_s + 62);
    const auto *id_s_63 = buffer.data(id_s + 63);
    const auto *id_s_64 = buffer.data(id_s + 64);
    const auto *id_s_65 = buffer.data(id_s + 65);
    const auto *id_s_66 = buffer.data(id_s + 66);
    const auto *id_s_67 = buffer.data(id_s + 67);
    const auto *id_s_68 = buffer.data(id_s + 68);
    const auto *id_s_69 = buffer.data(id_s + 69);
    const auto *id_s_70 = buffer.data(id_s + 70);
    const auto *id_s_71 = buffer.data(id_s + 71);
    const auto *id_s_72 = buffer.data(id_s + 72);
    const auto *id_s_73 = buffer.data(id_s + 73);
    const auto *id_s_74 = buffer.data(id_s + 74);
    const auto *id_s_75 = buffer.data(id_s + 75);
    const auto *id_s_76 = buffer.data(id_s + 76);
    const auto *id_s_77 = buffer.data(id_s + 77);
    const auto *id_s_78 = buffer.data(id_s + 78);
    const auto *id_s_79 = buffer.data(id_s + 79);
    const auto *id_s_80 = buffer.data(id_s + 80);
    const auto *id_s_81 = buffer.data(id_s + 81);
    const auto *id_s_82 = buffer.data(id_s + 82);
    const auto *id_s_83 = buffer.data(id_s + 83);
    const auto *id_s_84 = buffer.data(id_s + 84);
    const auto *id_s_85 = buffer.data(id_s + 85);
    const auto *id_s_86 = buffer.data(id_s + 86);
    const auto *id_s_87 = buffer.data(id_s + 87);
    const auto *id_s_88 = buffer.data(id_s + 88);
    const auto *id_s_89 = buffer.data(id_s + 89);
    const auto *id_s_90 = buffer.data(id_s + 90);
    const auto *id_s_91 = buffer.data(id_s + 91);
    const auto *id_s_92 = buffer.data(id_s + 92);
    const auto *id_s_93 = buffer.data(id_s + 93);
    const auto *id_s_94 = buffer.data(id_s + 94);
    const auto *id_s_95 = buffer.data(id_s + 95);
    const auto *id_s_96 = buffer.data(id_s + 96);
    const auto *id_s_97 = buffer.data(id_s + 97);
    const auto *id_s_98 = buffer.data(id_s + 98);
    const auto *id_s_99 = buffer.data(id_s + 99);
    const auto *id_s_100 = buffer.data(id_s + 100);
    const auto *id_s_101 = buffer.data(id_s + 101);
    const auto *id_s_102 = buffer.data(id_s + 102);
    const auto *id_s_103 = buffer.data(id_s + 103);
    const auto *id_s_104 = buffer.data(id_s + 104);
    const auto *id_s_105 = buffer.data(id_s + 105);
    const auto *id_s_106 = buffer.data(id_s + 106);
    const auto *id_s_107 = buffer.data(id_s + 107);
    const auto *id_s_108 = buffer.data(id_s + 108);
    const auto *id_s_109 = buffer.data(id_s + 109);
    const auto *id_s_110 = buffer.data(id_s + 110);
    const auto *id_s_112 = buffer.data(id_s + 112);
    const auto *id_s_113 = buffer.data(id_s + 113);
    const auto *id_s_114 = buffer.data(id_s + 114);
    const auto *id_s_115 = buffer.data(id_s + 115);
    const auto *id_s_116 = buffer.data(id_s + 116);
    const auto *id_s_117 = buffer.data(id_s + 117);
    const auto *id_s_118 = buffer.data(id_s + 118);
    const auto *id_s_119 = buffer.data(id_s + 119);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_39 = buffer.data(ip + 39);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is_s_0, id_s_0, id_s_1, \
                         id_s_2, is_0, ip_0, ip_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 - f_1 * is_s_0[k]
                 + f_2 * id_s_0[k]
                 + f_3 * is_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_2 * id_s_1[k]
                 + pb_z[k] * ip_0[k];

        t_2[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_2[k]
                 + f_3 * is_0[k]
                 + pb_y[k] * ip_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_x, pb_z, hp_3, hd_0, is_s_0, id_s_3, id_s_4, \
                         id_s_5, is_0, ip_2, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_3[k]
                 + f_3 * is_0[k]
                 + pb_z[k] * ip_2[k];

        t_4[k] = pa_y[k] * hd_0[k]
                 + f_2 * id_s_4[k];

        t_5[k] = f_4 * hp_3[k]
                 + f_2 * id_s_5[k]
                 + pb_x[k] * ip_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pa_z, gd_s_4, gd_4, hd_0, hd_2, hd_4, \
                         id_s_6, id_s_7, id_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_5 * gd_s_4[k]
                 + f_6 * gd_4[k]
                 + pa_x[k] * hd_4[k]
                 + f_2 * id_s_6[k];

        t_7[k] = pa_y[k] * hd_2[k]
                 + f_2 * id_s_7[k];

        t_8[k] = pa_z[k] * hd_0[k]
                 + f_2 * id_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, pb_y, gd_s_6, gd_6, hp_4, hd_7, id_s_9, \
                         id_s_10, id_s_11, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * hp_4[k]
                 + f_2 * id_s_9[k]
                 + pb_x[k] * ip_4[k];

        t_10[k] = f_2 * id_s_10[k]
                  + pb_y[k] * ip_4[k];

        t_11[k] = -f_5 * gd_s_6[k]
                  + f_6 * gd_6[k]
                  + pa_x[k] * hd_7[k]
                  + f_2 * id_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, pa_y, pb_x, gd_s_0, gd_0, hp_5, hd_3, id_s_12, id_s_13, \
                         ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_7 * gd_s_0[k]
                  + f_3 * gd_0[k]
                  + pa_y[k] * hd_3[k]
                  + f_2 * id_s_12[k];

        t_13[k] = f_6 * hp_5[k]
                  + f_2 * id_s_13[k]
                  + pb_x[k] * ip_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_z, gd_s_8, gd_8, hd_4, hd_9, is_s_3, \
                         id_s_14, id_s_15, id_s_16, is_3, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -f_8 * gd_s_8[k]
                  + f_9 * gd_8[k]
                  + pa_x[k] * hd_9[k]
                  + f_2 * id_s_14[k];

        t_15[k] = -f_1 * is_s_3[k]
                  + f_2 * id_s_15[k]
                  + f_3 * is_3[k]
                  + pb_z[k] * ip_6[k];

        t_16[k] = pa_z[k] * hd_4[k]
                  + f_2 * id_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_y, pa_z, pb_y, gd_s_0, gd_0, hp_4, hd_6, hd_7, \
                         id_s_17, id_s_18, id_s_19, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * hp_4[k]
                  + f_2 * id_s_17[k]
                  + pb_y[k] * ip_7[k];

        t_18[k] = pa_y[k] * hd_7[k]
                  + f_2 * id_s_18[k];

        t_19[k] = -f_7 * gd_s_0[k]
                  + f_3 * gd_0[k]
                  + pa_z[k] * hd_6[k]
                  + f_2 * id_s_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_y, hp_9, is_s_4, id_s_20, id_s_21, \
                         id_s_22, is_4, ip_8, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * hp_9[k]
                  + f_2 * id_s_20[k]
                  + pb_x[k] * ip_9[k];

        t_21[k] = -f_1 * is_s_4[k]
                  + f_2 * id_s_21[k]
                  + f_3 * is_4[k]
                  + pb_y[k] * ip_8[k];

        t_22[k] = f_2 * id_s_22[k]
                  + pb_y[k] * ip_9[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_x, pa_y, gd_s_3, gd_s_10, gd_3, gd_10, hd_8, hd_16, \
                         id_s_23, id_s_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_8 * gd_s_10[k]
                  + f_9 * gd_10[k]
                  + pa_x[k] * hd_16[k]
                  + f_2 * id_s_23[k];

        t_24[k] = -f_10 * gd_s_3[k]
                  + f_11 * gd_3[k]
                  + pa_y[k] * hd_8[k]
                  + f_2 * id_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pa_x, pb_x, gd_s_12, gd_12, hp_10, hd_18, id_s_25, \
                         id_s_26, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_9 * hp_10[k]
                  + f_2 * id_s_25[k]
                  + pb_x[k] * ip_10[k];

        t_26[k] = -f_10 * gd_s_12[k]
                  + f_11 * gd_12[k]
                  + pa_x[k] * hd_18[k]
                  + f_2 * id_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_z, pb_z, hd_8, hd_9, is_s_5, id_s_27, id_s_28, \
                         id_s_29, is_5, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * is_s_5[k]
                  + f_2 * id_s_27[k]
                  + f_3 * is_5[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_8[k]
                  + f_2 * id_s_28[k];

        t_29[k] = pa_z[k] * hd_9[k]
                  + f_2 * id_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_y, pb_y, gd_s_13, gd_13, hp_7, hd_14, \
                         hd_23, id_s_30, id_s_31, id_s_32, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_11 * hp_7[k]
                  + f_2 * id_s_30[k]
                  + pb_y[k] * ip_12[k];

        t_31[k] = -f_10 * gd_s_13[k]
                  + f_11 * gd_13[k]
                  + pa_x[k] * hd_23[k]
                  + f_2 * id_s_31[k];

        t_32[k] = pa_y[k] * hd_14[k]
                  + f_2 * id_s_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_x, pa_y, pb_y, gd_s_14, gd_14, hp_9, hd_16, \
                         hd_25, id_s_33, id_s_34, id_s_35, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -f_10 * gd_s_14[k]
                  + f_11 * gd_14[k]
                  + pa_x[k] * hd_25[k]
                  + f_2 * id_s_33[k];

        t_34[k] = f_3 * hp_9[k]
                  + f_2 * id_s_34[k]
                  + pb_y[k] * ip_13[k];

        t_35[k] = pa_y[k] * hd_16[k]
                  + f_2 * id_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, pa_z, pb_x, gd_s_5, gd_5, hp_15, hd_14, id_s_36, id_s_37, \
                         ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_10 * gd_s_5[k]
                  + f_11 * gd_5[k]
                  + pa_z[k] * hd_14[k]
                  + f_2 * id_s_36[k];

        t_37[k] = f_9 * hp_15[k]
                  + f_2 * id_s_37[k]
                  + pb_x[k] * ip_15[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pb_y, gd_s_16, gd_16, hd_30, is_s_7, id_s_38, \
                         id_s_39, id_s_40, is_7, ip_14, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -f_1 * is_s_7[k]
                  + f_2 * id_s_38[k]
                  + f_3 * is_7[k]
                  + pb_y[k] * ip_14[k];

        t_39[k] = f_2 * id_s_39[k]
                  + pb_y[k] * ip_15[k];

        t_40[k] = -f_10 * gd_s_16[k]
                  + f_11 * gd_16[k]
                  + pa_x[k] * hd_30[k]
                  + f_2 * id_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, pa_y, pb_x, gd_s_7, gd_7, hp_16, hd_17, id_s_41, id_s_42, \
                         ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -f_8 * gd_s_7[k]
                  + f_9 * gd_7[k]
                  + pa_y[k] * hd_17[k]
                  + f_2 * id_s_41[k];

        t_42[k] = f_11 * hp_16[k]
                  + f_2 * id_s_42[k]
                  + pb_x[k] * ip_16[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pa_z, pb_z, gd_s_18, gd_18, hd_17, hd_32, \
                         is_s_8, id_s_43, id_s_44, id_s_45, is_8, \
                         ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -f_7 * gd_s_18[k]
                  + f_3 * gd_18[k]
                  + pa_x[k] * hd_32[k]
                  + f_2 * id_s_43[k];

        t_44[k] = -f_1 * is_s_8[k]
                  + f_2 * id_s_44[k]
                  + f_3 * is_8[k]
                  + pb_z[k] * ip_17[k];

        t_45[k] = pa_z[k] * hd_17[k]
                  + f_2 * id_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_x, pa_z, pb_y, gd_s_22, gd_22, hp_12, hd_18, \
                         hd_36, id_s_46, id_s_47, id_s_48, ip_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * hd_18[k]
                  + f_2 * id_s_46[k];

        t_47[k] = f_9 * hp_12[k]
                  + f_2 * id_s_47[k]
                  + pb_y[k] * ip_18[k];

        t_48[k] = -f_7 * gd_s_22[k]
                  + f_3 * gd_22[k]
                  + pa_x[k] * hd_36[k]
                  + f_2 * id_s_48[k];
    }

#pragma omp simd aligned(t_49, t_50, pa_x, pa_y, gd_s_9, gd_s_24, gd_9, gd_24, hd_24, hd_38, \
                         id_s_49, id_s_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = -f_7 * gd_s_9[k]
                  + f_3 * gd_9[k]
                  + pa_y[k] * hd_24[k]
                  + f_2 * id_s_49[k];

        t_50[k] = -f_7 * gd_s_24[k]
                  + f_3 * gd_24[k]
                  + pa_x[k] * hd_38[k]
                  + f_2 * id_s_50[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pa_x, pa_y, pb_y, gd_s_26, gd_26, hp_13, hd_28, \
                         hd_40, id_s_51, id_s_52, id_s_53, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_11 * hp_13[k]
                  + f_2 * id_s_51[k]
                  + pb_y[k] * ip_19[k];

        t_52[k] = -f_7 * gd_s_26[k]
                  + f_3 * gd_26[k]
                  + pa_x[k] * hd_40[k]
                  + f_2 * id_s_52[k];

        t_53[k] = pa_y[k] * hd_28[k]
                  + f_2 * id_s_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_x, pa_y, pb_y, gd_s_27, gd_27, hp_15, hd_30, \
                         hd_42, id_s_54, id_s_55, id_s_56, ip_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = -f_7 * gd_s_27[k]
                  + f_3 * gd_27[k]
                  + pa_x[k] * hd_42[k]
                  + f_2 * id_s_54[k];

        t_55[k] = f_3 * hp_15[k]
                  + f_2 * id_s_55[k]
                  + pb_y[k] * ip_20[k];

        t_56[k] = pa_y[k] * hd_30[k]
                  + f_2 * id_s_56[k];
    }

#pragma omp simd aligned(t_57, t_58, pa_z, pb_x, gd_s_9, gd_9, hp_17, hd_28, id_s_57, id_s_58, \
                         ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -f_8 * gd_s_9[k]
                  + f_9 * gd_9[k]
                  + pa_z[k] * hd_28[k]
                  + f_2 * id_s_57[k];

        t_58[k] = f_11 * hp_17[k]
                  + f_2 * id_s_58[k]
                  + pb_x[k] * ip_22[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_x, pb_y, gd_s_32, gd_32, hd_46, is_s_11, \
                         id_s_59, id_s_60, id_s_61, is_11, ip_21, \
                         ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_1 * is_s_11[k]
                  + f_2 * id_s_59[k]
                  + f_3 * is_11[k]
                  + pb_y[k] * ip_21[k];

        t_60[k] = f_2 * id_s_60[k]
                  + pb_y[k] * ip_22[k];

        t_61[k] = -f_7 * gd_s_32[k]
                  + f_3 * gd_32[k]
                  + pa_x[k] * hd_46[k]
                  + f_2 * id_s_61[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pb_x, hp_18, hp_19, hd_47, hd_49, \
                         hd_50, id_s_62, id_s_63, id_s_64, id_s_65, \
                         ip_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_11 * hp_18[k]
                  + pa_x[k] * hd_47[k]
                  + f_2 * id_s_62[k];

        t_63[k] = f_3 * hp_19[k]
                  + f_2 * id_s_63[k]
                  + pb_x[k] * ip_23[k];

        t_64[k] = pa_x[k] * hd_49[k]
                  + f_2 * id_s_64[k];

        t_65[k] = pa_x[k] * hd_50[k]
                  + f_2 * id_s_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pa_z, hp_22, hd_31, hd_52, hd_53, \
                         hd_54, id_s_66, id_s_67, id_s_68, id_s_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * hd_31[k]
                  + f_2 * id_s_66[k];

        t_67[k] = pa_x[k] * hd_52[k]
                  + f_2 * id_s_67[k];

        t_68[k] = pa_x[k] * hd_53[k]
                  + f_2 * id_s_68[k];

        t_69[k] = f_11 * hp_22[k]
                  + pa_x[k] * hd_54[k]
                  + f_2 * id_s_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, hp_24, hd_55, hd_56, hd_57, hd_58, \
                         id_s_70, id_s_71, id_s_72, id_s_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_x[k] * hd_55[k]
                  + f_2 * id_s_70[k];

        t_71[k] = pa_x[k] * hd_56[k]
                  + f_2 * id_s_71[k];

        t_72[k] = pa_x[k] * hd_57[k]
                  + f_2 * id_s_72[k];

        t_73[k] = f_11 * hp_24[k]
                  + pa_x[k] * hd_58[k]
                  + f_2 * id_s_73[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pa_y, hd_44, hd_59, hd_60, hd_61, \
                         id_s_74, id_s_75, id_s_76, id_s_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_x[k] * hd_59[k]
                  + f_2 * id_s_74[k];

        t_75[k] = pa_x[k] * hd_60[k]
                  + f_2 * id_s_75[k];

        t_76[k] = pa_x[k] * hd_61[k]
                  + f_2 * id_s_76[k];

        t_77[k] = pa_y[k] * hd_44[k]
                  + f_2 * id_s_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_x, hp_27, hp_29, hd_62, hd_63, \
                         hd_65, id_s_78, id_s_79, id_s_80, id_s_81, \
                         ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_x[k] * hd_62[k]
                  + f_2 * id_s_78[k];

        t_79[k] = pa_x[k] * hd_63[k]
                  + f_2 * id_s_79[k];

        t_80[k] = f_11 * hp_27[k]
                  + pa_x[k] * hd_65[k]
                  + f_2 * id_s_80[k];

        t_81[k] = f_3 * hp_29[k]
                  + f_2 * id_s_81[k]
                  + pb_x[k] * ip_24[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pb_x, hd_67, hd_69, is_s_12, id_s_82, \
                         id_s_83, id_s_84, id_s_85, is_12, ip_25, \
                         ip_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_x[k] * hd_67[k]
                  + f_2 * id_s_82[k];

        t_83[k] = pa_x[k] * hd_69[k]
                  + f_2 * id_s_83[k];

        t_84[k] = -f_1 * is_s_12[k]
                  + f_2 * id_s_84[k]
                  + f_3 * is_12[k]
                  + pb_x[k] * ip_25[k];

        t_85[k] = f_2 * id_s_85[k]
                  + pb_x[k] * ip_26[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_y, pb_z, hp_19, is_s_12, id_s_86, id_s_87, \
                         id_s_88, is_12, ip_26, ip_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_0 * hp_19[k]
                  - f_1 * is_s_12[k]
                  + f_2 * id_s_86[k]
                  + f_3 * is_12[k]
                  + pb_y[k] * ip_26[k];

        t_87[k] = f_2 * id_s_87[k]
                  + pb_z[k] * ip_26[k];

        t_88[k] = -f_1 * is_s_12[k]
                  + f_2 * id_s_88[k]
                  + f_3 * is_12[k]
                  + pb_z[k] * ip_27[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_z, pb_x, pb_y, hp_21, hd_49, id_s_89, id_s_90, \
                         id_s_91, ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_2 * id_s_89[k]
                  + pb_x[k] * ip_28[k];

        t_90[k] = pa_z[k] * hd_49[k]
                  + f_2 * id_s_90[k];

        t_91[k] = f_4 * hp_21[k]
                  + f_2 * id_s_91[k]
                  + pb_y[k] * ip_28[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pa_y, pb_x, gd_s_22, gd_22, hd_53, is_s_14, \
                         id_s_92, id_s_93, id_s_94, is_14, ip_29, \
                         ip_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -f_5 * gd_s_22[k]
                  + f_6 * gd_22[k]
                  + pa_y[k] * hd_53[k]
                  + f_2 * id_s_92[k];

        t_93[k] = -f_1 * is_s_14[k]
                  + f_2 * id_s_93[k]
                  + f_3 * is_14[k]
                  + pb_x[k] * ip_29[k];

        t_94[k] = f_2 * id_s_94[k]
                  + pb_x[k] * ip_30[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pa_z, pb_x, pb_y, gd_s_18, gd_18, hp_23, hd_51, \
                         id_s_95, id_s_96, id_s_97, ip_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_2 * id_s_95[k]
                  + pb_x[k] * ip_31[k];

        t_96[k] = -f_7 * gd_s_18[k]
                  + f_3 * gd_18[k]
                  + pa_z[k] * hd_51[k]
                  + f_2 * id_s_96[k];

        t_97[k] = f_6 * hp_23[k]
                  + f_2 * id_s_97[k]
                  + pb_y[k] * ip_31[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_y, pb_x, gd_s_26, gd_26, hd_57, is_s_15, \
                         id_s_98, id_s_99, id_s_100, is_15, ip_32, \
                         ip_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = -f_8 * gd_s_26[k]
                  + f_9 * gd_26[k]
                  + pa_y[k] * hd_57[k]
                  + f_2 * id_s_98[k];

        t_99[k] = -f_1 * is_s_15[k]
                  + f_2 * id_s_99[k]
                  + f_3 * is_15[k]
                  + pb_x[k] * ip_32[k];

        t_100[k] = f_2 * id_s_100[k]
                   + pb_x[k] * ip_33[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_z, pb_x, pb_y, gd_s_20, gd_20, hp_25, hd_55, \
                         id_s_101, id_s_102, id_s_103, ip_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_2 * id_s_101[k]
                   + pb_x[k] * ip_34[k];

        t_102[k] = -f_10 * gd_s_20[k]
                   + f_11 * gd_20[k]
                   + pa_z[k] * hd_55[k]
                   + f_2 * id_s_102[k];

        t_103[k] = f_9 * hp_25[k]
                   + f_2 * id_s_103[k]
                   + pb_y[k] * ip_34[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pb_x, gd_s_29, gd_29, hd_61, is_s_16, \
                         id_s_104, id_s_105, id_s_106, is_16, ip_35, \
                         ip_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -f_10 * gd_s_29[k]
                   + f_11 * gd_29[k]
                   + pa_y[k] * hd_61[k]
                   + f_2 * id_s_104[k];

        t_105[k] = -f_1 * is_s_16[k]
                   + f_2 * id_s_105[k]
                   + f_3 * is_16[k]
                   + pb_x[k] * ip_35[k];

        t_106[k] = f_2 * id_s_106[k]
                   + pb_x[k] * ip_36[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pa_z, pb_x, pb_y, gd_s_24, gd_24, hp_26, hd_59, \
                         id_s_107, id_s_108, id_s_109, ip_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_2 * id_s_107[k]
                   + pb_x[k] * ip_37[k];

        t_108[k] = -f_8 * gd_s_24[k]
                   + f_9 * gd_24[k]
                   + pa_z[k] * hd_59[k]
                   + f_2 * id_s_108[k];

        t_109[k] = f_11 * hp_26[k]
                   + f_2 * id_s_109[k]
                   + pb_y[k] * ip_37[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_y, pb_y, gd_s_32, gd_32, hp_28, hp_29, hd_64, \
                         hd_67, id_s_110, id_s_112, id_s_113, ip_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -f_7 * gd_s_32[k]
                   + f_3 * gd_32[k]
                   + pa_y[k] * hd_64[k]
                   + f_2 * id_s_110[k];

        t_111[k] = f_11 * hp_28[k]
                   + pa_y[k] * hd_67[k]
                   + f_2 * id_s_112[k];

        t_112[k] = f_3 * hp_29[k]
                   + f_2 * id_s_113[k]
                   + pb_y[k] * ip_38[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_y, pb_x, hd_69, is_s_18, id_s_114, id_s_115, \
                         id_s_116, is_17, ip_39, ip_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_y[k] * hd_69[k]
                   + f_2 * id_s_114[k];

        t_114[k] = -f_1 * is_s_18[k]
                   + f_2 * id_s_115[k]
                   + f_3 * is_17[k]
                   + pb_x[k] * ip_39[k];

        t_115[k] = f_2 * id_s_116[k]
                   + pb_x[k] * ip_41[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pb_y, pb_z, hp_29, is_s_18, id_s_117, id_s_118, \
                         id_s_119, is_17, ip_40, ip_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -f_1 * is_s_18[k]
                   + f_2 * id_s_117[k]
                   + f_3 * is_17[k]
                   + pb_y[k] * ip_40[k];

        t_117[k] = f_2 * id_s_118[k]
                   + pb_y[k] * ip_41[k];

        t_118[k] = f_0 * hp_29[k]
                   - f_1 * is_s_18[k]
                   + f_2 * id_s_119[k]
                   + f_3 * is_17[k]
                   + pb_z[k] * ip_41[k];
    }
}

auto
compute_prim_id_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gd_s, const size_t gd,
                                 const size_t hp, const size_t hd, const size_t is_s,
                                 const size_t id_s, const size_t is, const size_t ip,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 4.0 * beta / p;
    const auto f_5 = 2.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 3.0 * beta / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 2.0 * beta / p;
    const auto f_10 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_14 = buffer.data(gd_s + 14);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_19 = buffer.data(gd_s + 19);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_25 = buffer.data(gd_s + 25);
    const auto *gd_s_26 = buffer.data(gd_s + 26);
    const auto *gd_s_28 = buffer.data(gd_s + 28);
    const auto *gd_s_29 = buffer.data(gd_s + 29);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_31 = buffer.data(gd_s + 31);
    const auto *gd_s_35 = buffer.data(gd_s + 35);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_35 = buffer.data(gd + 35);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_59 = buffer.data(hd + 59);

    const auto *is_s_0 = buffer.data(is_s + 0);
    const auto *is_s_2 = buffer.data(is_s + 2);
    const auto *is_s_3 = buffer.data(is_s + 3);
    const auto *is_s_4 = buffer.data(is_s + 4);
    const auto *is_s_5 = buffer.data(is_s + 5);
    const auto *is_s_6 = buffer.data(is_s + 6);
    const auto *is_s_7 = buffer.data(is_s + 7);
    const auto *is_s_8 = buffer.data(is_s + 8);
    const auto *is_s_10 = buffer.data(is_s + 10);
    const auto *is_s_11 = buffer.data(is_s + 11);
    const auto *is_s_12 = buffer.data(is_s + 12);
    const auto *is_s_14 = buffer.data(is_s + 14);

    const auto *id_s_0 = buffer.data(id_s + 0);
    const auto *id_s_1 = buffer.data(id_s + 1);
    const auto *id_s_2 = buffer.data(id_s + 2);
    const auto *id_s_3 = buffer.data(id_s + 3);
    const auto *id_s_4 = buffer.data(id_s + 4);
    const auto *id_s_5 = buffer.data(id_s + 5);
    const auto *id_s_6 = buffer.data(id_s + 6);
    const auto *id_s_7 = buffer.data(id_s + 7);
    const auto *id_s_8 = buffer.data(id_s + 8);
    const auto *id_s_9 = buffer.data(id_s + 9);
    const auto *id_s_10 = buffer.data(id_s + 10);
    const auto *id_s_11 = buffer.data(id_s + 11);
    const auto *id_s_12 = buffer.data(id_s + 12);
    const auto *id_s_13 = buffer.data(id_s + 13);
    const auto *id_s_14 = buffer.data(id_s + 14);
    const auto *id_s_15 = buffer.data(id_s + 15);
    const auto *id_s_16 = buffer.data(id_s + 16);
    const auto *id_s_17 = buffer.data(id_s + 17);
    const auto *id_s_18 = buffer.data(id_s + 18);
    const auto *id_s_19 = buffer.data(id_s + 19);
    const auto *id_s_20 = buffer.data(id_s + 20);
    const auto *id_s_21 = buffer.data(id_s + 21);
    const auto *id_s_22 = buffer.data(id_s + 22);
    const auto *id_s_23 = buffer.data(id_s + 23);
    const auto *id_s_24 = buffer.data(id_s + 24);
    const auto *id_s_25 = buffer.data(id_s + 25);
    const auto *id_s_26 = buffer.data(id_s + 26);
    const auto *id_s_27 = buffer.data(id_s + 27);
    const auto *id_s_28 = buffer.data(id_s + 28);
    const auto *id_s_29 = buffer.data(id_s + 29);
    const auto *id_s_30 = buffer.data(id_s + 30);
    const auto *id_s_31 = buffer.data(id_s + 31);
    const auto *id_s_32 = buffer.data(id_s + 32);
    const auto *id_s_33 = buffer.data(id_s + 33);
    const auto *id_s_34 = buffer.data(id_s + 34);
    const auto *id_s_35 = buffer.data(id_s + 35);
    const auto *id_s_36 = buffer.data(id_s + 36);
    const auto *id_s_37 = buffer.data(id_s + 37);
    const auto *id_s_38 = buffer.data(id_s + 38);
    const auto *id_s_39 = buffer.data(id_s + 39);
    const auto *id_s_40 = buffer.data(id_s + 40);
    const auto *id_s_41 = buffer.data(id_s + 41);
    const auto *id_s_42 = buffer.data(id_s + 42);
    const auto *id_s_43 = buffer.data(id_s + 43);
    const auto *id_s_44 = buffer.data(id_s + 44);
    const auto *id_s_45 = buffer.data(id_s + 45);
    const auto *id_s_46 = buffer.data(id_s + 46);
    const auto *id_s_47 = buffer.data(id_s + 47);
    const auto *id_s_48 = buffer.data(id_s + 48);
    const auto *id_s_49 = buffer.data(id_s + 49);
    const auto *id_s_50 = buffer.data(id_s + 50);
    const auto *id_s_51 = buffer.data(id_s + 51);
    const auto *id_s_52 = buffer.data(id_s + 52);
    const auto *id_s_53 = buffer.data(id_s + 53);
    const auto *id_s_54 = buffer.data(id_s + 54);
    const auto *id_s_55 = buffer.data(id_s + 55);
    const auto *id_s_56 = buffer.data(id_s + 56);
    const auto *id_s_57 = buffer.data(id_s + 57);
    const auto *id_s_59 = buffer.data(id_s + 59);
    const auto *id_s_60 = buffer.data(id_s + 60);
    const auto *id_s_61 = buffer.data(id_s + 61);
    const auto *id_s_62 = buffer.data(id_s + 62);
    const auto *id_s_63 = buffer.data(id_s + 63);
    const auto *id_s_64 = buffer.data(id_s + 64);
    const auto *id_s_65 = buffer.data(id_s + 65);
    const auto *id_s_66 = buffer.data(id_s + 66);
    const auto *id_s_67 = buffer.data(id_s + 67);
    const auto *id_s_68 = buffer.data(id_s + 68);
    const auto *id_s_69 = buffer.data(id_s + 69);
    const auto *id_s_70 = buffer.data(id_s + 70);
    const auto *id_s_71 = buffer.data(id_s + 71);
    const auto *id_s_72 = buffer.data(id_s + 72);
    const auto *id_s_73 = buffer.data(id_s + 73);
    const auto *id_s_74 = buffer.data(id_s + 74);
    const auto *id_s_75 = buffer.data(id_s + 75);
    const auto *id_s_76 = buffer.data(id_s + 76);
    const auto *id_s_77 = buffer.data(id_s + 77);
    const auto *id_s_79 = buffer.data(id_s + 79);
    const auto *id_s_81 = buffer.data(id_s + 81);
    const auto *id_s_82 = buffer.data(id_s + 82);
    const auto *id_s_83 = buffer.data(id_s + 83);
    const auto *id_s_84 = buffer.data(id_s + 84);
    const auto *id_s_85 = buffer.data(id_s + 85);
    const auto *id_s_86 = buffer.data(id_s + 86);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is_s_0, id_s_0, id_s_1, \
                         id_s_2, is_0, ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 - f_1 * is_s_0[k]
                 + f_2 * id_s_0[k]
                 + f_3 * is_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_1[k]
                 + f_3 * is_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_2[k]
                 + f_3 * is_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, gd_s_4, gd_4, hd_0, hd_2, hd_4, \
                         id_s_3, id_s_4, id_s_5, id_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * hd_0[k]
                 + f_2 * id_s_3[k];

        t_4[k] = -f_4 * gd_s_4[k]
                 + f_5 * gd_4[k]
                 + pa_x[k] * hd_4[k]
                 + f_2 * id_s_4[k];

        t_5[k] = pa_y[k] * hd_2[k]
                 + f_2 * id_s_5[k];

        t_6[k] = pa_z[k] * hd_0[k]
                 + f_2 * id_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pb_y, gd_s_0, gd_s_7, gd_0, gd_7, hd_3, \
                         hd_8, id_s_7, id_s_8, id_s_9, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * id_s_7[k]
                 + pb_y[k] * ip_3[k];

        t_8[k] = -f_4 * gd_s_7[k]
                 + f_5 * gd_7[k]
                 + pa_x[k] * hd_8[k]
                 + f_2 * id_s_8[k];

        t_9[k] = -f_6 * gd_s_0[k]
                 + f_3 * gd_0[k]
                 + pa_y[k] * hd_3[k]
                 + f_2 * id_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_z, gd_s_9, gd_9, hd_4, hd_10, \
                         is_s_2, id_s_10, id_s_11, id_s_12, is_2, \
                         ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_7 * gd_s_9[k]
                  + f_8 * gd_9[k]
                  + pa_x[k] * hd_10[k]
                  + f_2 * id_s_10[k];

        t_11[k] = -f_1 * is_s_2[k]
                  + f_2 * id_s_11[k]
                  + f_3 * is_2[k]
                  + pb_z[k] * ip_4[k];

        t_12[k] = pa_z[k] * hd_4[k]
                  + f_2 * id_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_y, gd_s_0, gd_0, hd_6, hd_8, is_s_3, \
                         id_s_13, id_s_14, id_s_15, is_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * hd_8[k]
                  + f_2 * id_s_13[k];

        t_14[k] = -f_6 * gd_s_0[k]
                  + f_3 * gd_0[k]
                  + pa_z[k] * hd_6[k]
                  + f_2 * id_s_14[k];

        t_15[k] = -f_1 * is_s_3[k]
                  + f_2 * id_s_15[k]
                  + f_3 * is_3[k]
                  + pb_y[k] * ip_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_y, gd_s_3, gd_s_14, gd_3, gd_14, \
                         hd_9, hd_17, id_s_16, id_s_17, id_s_18, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * id_s_16[k]
                  + pb_y[k] * ip_6[k];

        t_17[k] = -f_7 * gd_s_14[k]
                  + f_8 * gd_14[k]
                  + pa_x[k] * hd_17[k]
                  + f_2 * id_s_17[k];

        t_18[k] = -f_9 * gd_s_3[k]
                  + f_10 * gd_3[k]
                  + pa_y[k] * hd_9[k]
                  + f_2 * id_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_z, pb_z, gd_s_16, gd_16, hd_9, hd_19, \
                         is_s_4, id_s_19, id_s_20, id_s_21, is_4, \
                         ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_9 * gd_s_16[k]
                  + f_10 * gd_16[k]
                  + pa_x[k] * hd_19[k]
                  + f_2 * id_s_19[k];

        t_20[k] = -f_1 * is_s_4[k]
                  + f_2 * id_s_20[k]
                  + f_3 * is_4[k]
                  + pb_z[k] * ip_7[k];

        t_21[k] = pa_z[k] * hd_9[k]
                  + f_2 * id_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_z, gd_s_18, gd_s_19, gd_18, gd_19, hd_10, \
                         hd_23, hd_25, id_s_22, id_s_23, id_s_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_z[k] * hd_10[k]
                  + f_2 * id_s_22[k];

        t_23[k] = -f_9 * gd_s_18[k]
                  + f_10 * gd_18[k]
                  + pa_x[k] * hd_23[k]
                  + f_2 * id_s_23[k];

        t_24[k] = -f_9 * gd_s_19[k]
                  + f_10 * gd_19[k]
                  + pa_x[k] * hd_25[k]
                  + f_2 * id_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_y, gd_s_6, gd_6, hd_14, hd_17, \
                         is_s_5, id_s_25, id_s_26, id_s_27, is_5, \
                         ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_y[k] * hd_17[k]
                  + f_2 * id_s_25[k];

        t_26[k] = -f_9 * gd_s_6[k]
                  + f_10 * gd_6[k]
                  + pa_z[k] * hd_14[k]
                  + f_2 * id_s_26[k];

        t_27[k] = -f_1 * is_s_5[k]
                  + f_2 * id_s_27[k]
                  + f_3 * is_5[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pa_y, pb_y, gd_s_8, gd_s_21, gd_8, gd_21, \
                         hd_18, hd_30, id_s_28, id_s_29, id_s_30, \
                         ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_2 * id_s_28[k]
                  + pb_y[k] * ip_9[k];

        t_29[k] = -f_9 * gd_s_21[k]
                  + f_10 * gd_21[k]
                  + pa_x[k] * hd_30[k]
                  + f_2 * id_s_29[k];

        t_30[k] = -f_7 * gd_s_8[k]
                  + f_8 * gd_8[k]
                  + pa_y[k] * hd_18[k]
                  + f_2 * id_s_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_x, pa_z, pb_z, gd_s_23, gd_23, hd_18, hd_32, \
                         is_s_6, id_s_31, id_s_32, id_s_33, is_6, \
                         ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_6 * gd_s_23[k]
                  + f_3 * gd_23[k]
                  + pa_x[k] * hd_32[k]
                  + f_2 * id_s_31[k];

        t_32[k] = -f_1 * is_s_6[k]
                  + f_2 * id_s_32[k]
                  + f_3 * is_6[k]
                  + pb_z[k] * ip_10[k];

        t_33[k] = pa_z[k] * hd_18[k]
                  + f_2 * id_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_y, pa_z, gd_s_13, gd_s_26, gd_13, gd_26, \
                         hd_19, hd_24, hd_34, id_s_34, id_s_35, \
                         id_s_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * hd_19[k]
                  + f_2 * id_s_34[k];

        t_35[k] = -f_6 * gd_s_26[k]
                  + f_3 * gd_26[k]
                  + pa_x[k] * hd_34[k]
                  + f_2 * id_s_35[k];

        t_36[k] = -f_6 * gd_s_13[k]
                  + f_3 * gd_13[k]
                  + pa_y[k] * hd_24[k]
                  + f_2 * id_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_x, gd_s_28, gd_s_29, gd_s_30, gd_28, gd_29, \
                         gd_30, hd_36, hd_37, hd_38, id_s_37, id_s_38, \
                         id_s_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_6 * gd_s_28[k]
                  + f_3 * gd_28[k]
                  + pa_x[k] * hd_36[k]
                  + f_2 * id_s_37[k];

        t_38[k] = -f_6 * gd_s_29[k]
                  + f_3 * gd_29[k]
                  + pa_x[k] * hd_37[k]
                  + f_2 * id_s_38[k];

        t_39[k] = -f_6 * gd_s_30[k]
                  + f_3 * gd_30[k]
                  + pa_x[k] * hd_38[k]
                  + f_2 * id_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_y, gd_s_13, gd_13, hd_27, hd_30, \
                         is_s_7, id_s_40, id_s_41, id_s_42, is_7, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * hd_30[k]
                  + f_2 * id_s_40[k];

        t_41[k] = -f_7 * gd_s_13[k]
                  + f_8 * gd_13[k]
                  + pa_z[k] * hd_27[k]
                  + f_2 * id_s_41[k];

        t_42[k] = -f_1 * is_s_7[k]
                  + f_2 * id_s_42[k]
                  + f_3 * is_7[k]
                  + pb_y[k] * ip_11[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pb_y, gd_s_35, gd_35, hp_10, hd_40, hd_41, \
                         id_s_43, id_s_44, id_s_45, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * id_s_43[k]
                  + pb_y[k] * ip_12[k];

        t_44[k] = -f_6 * gd_s_35[k]
                  + f_3 * gd_35[k]
                  + pa_x[k] * hd_40[k]
                  + f_2 * id_s_44[k];

        t_45[k] = f_10 * hp_10[k]
                  + pa_x[k] * hd_41[k]
                  + f_2 * id_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pa_z, hp_13, hp_15, hd_31, hd_43, \
                         hd_47, hd_50, id_s_46, id_s_47, id_s_48, \
                         id_s_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_x[k] * hd_43[k]
                  + f_2 * id_s_46[k];

        t_47[k] = pa_z[k] * hd_31[k]
                  + f_2 * id_s_47[k];

        t_48[k] = f_10 * hp_13[k]
                  + pa_x[k] * hd_47[k]
                  + f_2 * id_s_48[k];

        t_49[k] = f_10 * hp_15[k]
                  + pa_x[k] * hd_50[k]
                  + f_2 * id_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pb_x, hp_18, hd_55, hd_59, is_s_8, id_s_50, \
                         id_s_51, id_s_52, is_8, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_10 * hp_18[k]
                  + pa_x[k] * hd_55[k]
                  + f_2 * id_s_50[k];

        t_51[k] = pa_x[k] * hd_59[k]
                  + f_2 * id_s_51[k];

        t_52[k] = -f_1 * is_s_8[k]
                  + f_2 * id_s_52[k]
                  + f_3 * is_8[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, pb_z, hp_11, is_s_8, id_s_53, id_s_54, \
                         id_s_55, is_8, ip_14, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_2 * id_s_53[k]
                  + pb_x[k] * ip_14[k];

        t_54[k] = f_0 * hp_11[k]
                  - f_1 * is_s_8[k]
                  + f_2 * id_s_54[k]
                  + f_3 * is_8[k]
                  + pb_y[k] * ip_14[k];

        t_55[k] = -f_1 * is_s_8[k]
                  + f_2 * id_s_55[k]
                  + f_3 * is_8[k]
                  + pb_z[k] * ip_15[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_y, pa_z, pb_x, gd_s_26, gd_26, hd_43, hd_46, \
                         id_s_56, id_s_57, id_s_59, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_2 * id_s_56[k]
                  + pb_x[k] * ip_16[k];

        t_57[k] = pa_z[k] * hd_43[k]
                  + f_2 * id_s_57[k];

        t_58[k] = -f_4 * gd_s_26[k]
                  + f_5 * gd_26[k]
                  + pa_y[k] * hd_46[k]
                  + f_2 * id_s_59[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_x, is_s_10, id_s_60, id_s_61, id_s_62, is_10, \
                         ip_17, ip_18, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_1 * is_s_10[k]
                  + f_2 * id_s_60[k]
                  + f_3 * is_10[k]
                  + pb_x[k] * ip_17[k];

        t_60[k] = f_2 * id_s_61[k]
                  + pb_x[k] * ip_18[k];

        t_61[k] = f_2 * id_s_62[k]
                  + pb_x[k] * ip_19[k];
    }

#pragma omp simd aligned(t_62, t_63, pa_z, pb_y, gd_s_23, gd_23, hp_14, hd_45, id_s_63, \
                         id_s_64, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -f_6 * gd_s_23[k]
                  + f_3 * gd_23[k]
                  + pa_z[k] * hd_45[k]
                  + f_2 * id_s_63[k];

        t_63[k] = f_5 * hp_14[k]
                  + f_2 * id_s_64[k]
                  + pb_y[k] * ip_19[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_y, pb_x, gd_s_29, gd_29, hd_49, is_s_11, \
                         id_s_65, id_s_66, id_s_67, is_11, ip_20, \
                         ip_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = -f_7 * gd_s_29[k]
                  + f_8 * gd_29[k]
                  + pa_y[k] * hd_49[k]
                  + f_2 * id_s_65[k];

        t_65[k] = -f_1 * is_s_11[k]
                  + f_2 * id_s_66[k]
                  + f_3 * is_11[k]
                  + pb_x[k] * ip_20[k];

        t_66[k] = f_2 * id_s_67[k]
                  + pb_x[k] * ip_21[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_z, pb_x, pb_y, gd_s_25, gd_25, hp_16, hd_48, \
                         id_s_68, id_s_69, id_s_70, ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_2 * id_s_68[k]
                  + pb_x[k] * ip_22[k];

        t_68[k] = -f_9 * gd_s_25[k]
                  + f_10 * gd_25[k]
                  + pa_z[k] * hd_48[k]
                  + f_2 * id_s_69[k];

        t_69[k] = f_8 * hp_16[k]
                  + f_2 * id_s_70[k]
                  + pb_y[k] * ip_22[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pb_x, gd_s_31, gd_31, hd_52, is_s_12, \
                         id_s_71, id_s_72, id_s_73, is_12, ip_23, \
                         ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -f_9 * gd_s_31[k]
                  + f_10 * gd_31[k]
                  + pa_y[k] * hd_52[k]
                  + f_2 * id_s_71[k];

        t_71[k] = -f_1 * is_s_12[k]
                  + f_2 * id_s_72[k]
                  + f_3 * is_12[k]
                  + pb_x[k] * ip_23[k];

        t_72[k] = f_2 * id_s_73[k]
                  + pb_x[k] * ip_24[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_z, pb_x, pb_y, gd_s_28, gd_28, hp_17, hd_51, \
                         id_s_74, id_s_75, id_s_76, ip_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_2 * id_s_74[k]
                  + pb_x[k] * ip_25[k];

        t_74[k] = -f_7 * gd_s_28[k]
                  + f_8 * gd_28[k]
                  + pa_z[k] * hd_51[k]
                  + f_2 * id_s_75[k];

        t_75[k] = f_10 * hp_17[k]
                  + f_2 * id_s_76[k]
                  + pb_y[k] * ip_25[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_y, gd_s_35, gd_35, hp_19, hd_54, hd_57, hd_59, \
                         id_s_77, id_s_79, id_s_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -f_6 * gd_s_35[k]
                  + f_3 * gd_35[k]
                  + pa_y[k] * hd_54[k]
                  + f_2 * id_s_77[k];

        t_77[k] = f_10 * hp_19[k]
                  + pa_y[k] * hd_57[k]
                  + f_2 * id_s_79[k];

        t_78[k] = pa_y[k] * hd_59[k]
                  + f_2 * id_s_81[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_x, pb_y, is_s_14, id_s_82, id_s_83, \
                         id_s_84, id_s_85, is_13, ip_26, ip_27, ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -f_1 * is_s_14[k]
                  + f_2 * id_s_82[k]
                  + f_3 * is_13[k]
                  + pb_x[k] * ip_26[k];

        t_80[k] = f_2 * id_s_83[k]
                  + pb_x[k] * ip_28[k];

        t_81[k] = -f_1 * is_s_14[k]
                  + f_2 * id_s_84[k]
                  + f_3 * is_13[k]
                  + pb_y[k] * ip_27[k];

        t_82[k] = f_2 * id_s_85[k]
                  + pb_y[k] * ip_28[k];
    }

#pragma omp simd aligned(t_83, pb_z, hp_20, is_s_14, id_s_86, is_13, \
                         ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * hp_20[k]
                  - f_1 * is_s_14[k]
                  + f_2 * id_s_86[k]
                  + f_3 * is_13[k]
                  + pb_z[k] * ip_28[k];
    }
}

auto
compute_prim_id_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gd_s, const size_t gd,
                                 const size_t hp, const size_t hd, const size_t is_s,
                                 const size_t id_s, const size_t is, const size_t ip,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 4.0 * beta / p;
    const auto f_5 = 2.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 3.0 * beta / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 2.0 * beta / p;
    const auto f_10 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_27 = buffer.data(gd_s + 27);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_27 = buffer.data(gd + 27);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_20 = buffer.data(hp + 20);

    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_37 = buffer.data(hd + 37);

    const auto *is_s_0 = buffer.data(is_s + 0);
    const auto *is_s_2 = buffer.data(is_s + 2);
    const auto *is_s_3 = buffer.data(is_s + 3);
    const auto *is_s_4 = buffer.data(is_s + 4);
    const auto *is_s_5 = buffer.data(is_s + 5);
    const auto *is_s_6 = buffer.data(is_s + 6);
    const auto *is_s_7 = buffer.data(is_s + 7);
    const auto *is_s_8 = buffer.data(is_s + 8);
    const auto *is_s_10 = buffer.data(is_s + 10);
    const auto *is_s_11 = buffer.data(is_s + 11);
    const auto *is_s_12 = buffer.data(is_s + 12);
    const auto *is_s_14 = buffer.data(is_s + 14);

    const auto *id_s_0 = buffer.data(id_s + 0);
    const auto *id_s_1 = buffer.data(id_s + 1);
    const auto *id_s_2 = buffer.data(id_s + 2);
    const auto *id_s_3 = buffer.data(id_s + 3);
    const auto *id_s_5 = buffer.data(id_s + 5);
    const auto *id_s_6 = buffer.data(id_s + 6);
    const auto *id_s_7 = buffer.data(id_s + 7);
    const auto *id_s_8 = buffer.data(id_s + 8);
    const auto *id_s_9 = buffer.data(id_s + 9);
    const auto *id_s_10 = buffer.data(id_s + 10);
    const auto *id_s_11 = buffer.data(id_s + 11);
    const auto *id_s_12 = buffer.data(id_s + 12);
    const auto *id_s_13 = buffer.data(id_s + 13);
    const auto *id_s_14 = buffer.data(id_s + 14);
    const auto *id_s_15 = buffer.data(id_s + 15);
    const auto *id_s_16 = buffer.data(id_s + 16);
    const auto *id_s_17 = buffer.data(id_s + 17);
    const auto *id_s_18 = buffer.data(id_s + 18);
    const auto *id_s_19 = buffer.data(id_s + 19);
    const auto *id_s_20 = buffer.data(id_s + 20);
    const auto *id_s_21 = buffer.data(id_s + 21);
    const auto *id_s_22 = buffer.data(id_s + 22);
    const auto *id_s_23 = buffer.data(id_s + 23);
    const auto *id_s_24 = buffer.data(id_s + 24);
    const auto *id_s_25 = buffer.data(id_s + 25);
    const auto *id_s_26 = buffer.data(id_s + 26);
    const auto *id_s_27 = buffer.data(id_s + 27);
    const auto *id_s_30 = buffer.data(id_s + 30);
    const auto *id_s_31 = buffer.data(id_s + 31);
    const auto *id_s_32 = buffer.data(id_s + 32);
    const auto *id_s_33 = buffer.data(id_s + 33);
    const auto *id_s_34 = buffer.data(id_s + 34);
    const auto *id_s_37 = buffer.data(id_s + 37);
    const auto *id_s_38 = buffer.data(id_s + 38);
    const auto *id_s_39 = buffer.data(id_s + 39);
    const auto *id_s_40 = buffer.data(id_s + 40);
    const auto *id_s_41 = buffer.data(id_s + 41);
    const auto *id_s_42 = buffer.data(id_s + 42);
    const auto *id_s_43 = buffer.data(id_s + 43);
    const auto *id_s_44 = buffer.data(id_s + 44);
    const auto *id_s_45 = buffer.data(id_s + 45);
    const auto *id_s_46 = buffer.data(id_s + 46);
    const auto *id_s_47 = buffer.data(id_s + 47);
    const auto *id_s_48 = buffer.data(id_s + 48);
    const auto *id_s_49 = buffer.data(id_s + 49);
    const auto *id_s_50 = buffer.data(id_s + 50);
    const auto *id_s_51 = buffer.data(id_s + 51);
    const auto *id_s_52 = buffer.data(id_s + 52);
    const auto *id_s_53 = buffer.data(id_s + 53);
    const auto *id_s_54 = buffer.data(id_s + 54);
    const auto *id_s_55 = buffer.data(id_s + 55);
    const auto *id_s_60 = buffer.data(id_s + 60);
    const auto *id_s_61 = buffer.data(id_s + 61);
    const auto *id_s_62 = buffer.data(id_s + 62);
    const auto *id_s_63 = buffer.data(id_s + 63);
    const auto *id_s_64 = buffer.data(id_s + 64);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is_s_0, id_s_0, id_s_1, \
                         id_s_2, is_0, ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 - f_1 * is_s_0[k]
                 + f_2 * id_s_0[k]
                 + f_3 * is_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_1[k]
                 + f_3 * is_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_2[k]
                 + f_3 * is_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, gd_s_4, gd_s_6, gd_4, gd_6, hd_4, hd_7, \
                         id_s_3, id_s_5, id_s_6, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * gd_s_4[k]
                 + f_5 * gd_4[k]
                 + pa_x[k] * hd_4[k]
                 + f_2 * id_s_3[k];

        t_4[k] = f_2 * id_s_5[k]
                 + pb_y[k] * ip_3[k];

        t_5[k] = -f_4 * gd_s_6[k]
                 + f_5 * gd_6[k]
                 + pa_x[k] * hd_7[k]
                 + f_2 * id_s_6[k];
    }

#pragma omp simd aligned(t_6, t_7, pa_x, pa_y, gd_s_0, gd_s_8, gd_0, gd_8, hd_3, hd_9, id_s_7, \
                         id_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_6 * gd_s_0[k]
                 + f_3 * gd_0[k]
                 + pa_y[k] * hd_3[k]
                 + f_2 * id_s_7[k];

        t_7[k] = -f_7 * gd_s_8[k]
                 + f_8 * gd_8[k]
                 + pa_x[k] * hd_9[k]
                 + f_2 * id_s_8[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_z, pb_z, gd_s_0, gd_0, hd_5, is_s_2, id_s_9, id_s_10, \
                         is_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_1 * is_s_2[k]
                 + f_2 * id_s_9[k]
                 + f_3 * is_2[k]
                 + pb_z[k] * ip_4[k];

        t_9[k] = -f_6 * gd_s_0[k]
                 + f_3 * gd_0[k]
                 + pa_z[k] * hd_5[k]
                 + f_2 * id_s_10[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_y, gd_s_11, gd_11, hd_14, is_s_3, id_s_11, \
                         id_s_12, id_s_13, is_3, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * is_s_3[k]
                  + f_2 * id_s_11[k]
                  + f_3 * is_3[k]
                  + pb_y[k] * ip_5[k];

        t_11[k] = f_2 * id_s_12[k]
                  + pb_y[k] * ip_6[k];

        t_12[k] = -f_7 * gd_s_11[k]
                  + f_8 * gd_11[k]
                  + pa_x[k] * hd_14[k]
                  + f_2 * id_s_13[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, gd_s_3, gd_s_12, gd_3, gd_12, hd_8, hd_16, \
                         id_s_14, id_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_9 * gd_s_3[k]
                  + f_10 * gd_3[k]
                  + pa_y[k] * hd_8[k]
                  + f_2 * id_s_14[k];

        t_14[k] = -f_9 * gd_s_12[k]
                  + f_10 * gd_12[k]
                  + pa_x[k] * hd_16[k]
                  + f_2 * id_s_15[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_z, pb_z, gd_s_5, gd_5, hd_11, is_s_4, id_s_16, \
                         id_s_17, is_4, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * is_s_4[k]
                  + f_2 * id_s_16[k]
                  + f_3 * is_4[k]
                  + pb_z[k] * ip_7[k];

        t_16[k] = -f_9 * gd_s_5[k]
                  + f_10 * gd_5[k]
                  + pa_z[k] * hd_11[k]
                  + f_2 * id_s_17[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pb_y, gd_s_13, gd_13, hd_21, is_s_5, id_s_18, \
                         id_s_19, id_s_20, is_5, ip_8, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_1 * is_s_5[k]
                  + f_2 * id_s_18[k]
                  + f_3 * is_5[k]
                  + pb_y[k] * ip_8[k];

        t_18[k] = f_2 * id_s_19[k]
                  + pb_y[k] * ip_9[k];

        t_19[k] = -f_9 * gd_s_13[k]
                  + f_10 * gd_13[k]
                  + pa_x[k] * hd_21[k]
                  + f_2 * id_s_20[k];
    }

#pragma omp simd aligned(t_20, t_21, pa_x, pa_y, gd_s_7, gd_s_15, gd_7, gd_15, hd_15, hd_22, \
                         id_s_21, id_s_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_7 * gd_s_7[k]
                  + f_8 * gd_7[k]
                  + pa_y[k] * hd_15[k]
                  + f_2 * id_s_21[k];

        t_21[k] = -f_6 * gd_s_15[k]
                  + f_3 * gd_15[k]
                  + pa_x[k] * hd_22[k]
                  + f_2 * id_s_22[k];
    }

#pragma omp simd aligned(t_22, t_23, pa_z, pb_z, gd_s_10, gd_10, hd_18, is_s_6, id_s_23, \
                         id_s_24, is_6, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_1 * is_s_6[k]
                  + f_2 * id_s_23[k]
                  + f_3 * is_6[k]
                  + pb_z[k] * ip_10[k];

        t_23[k] = -f_7 * gd_s_10[k]
                  + f_8 * gd_10[k]
                  + pa_z[k] * hd_18[k]
                  + f_2 * id_s_24[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_y, gd_s_27, gd_27, hd_23, is_s_7, id_s_25, \
                         id_s_26, id_s_27, is_7, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * is_s_7[k]
                  + f_2 * id_s_25[k]
                  + f_3 * is_7[k]
                  + pb_y[k] * ip_11[k];

        t_25[k] = f_2 * id_s_26[k]
                  + pb_y[k] * ip_12[k];

        t_26[k] = -f_6 * gd_s_27[k]
                  + f_3 * gd_27[k]
                  + pa_x[k] * hd_23[k]
                  + f_2 * id_s_27[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, pb_y, hp_11, is_s_8, id_s_30, id_s_31, \
                         id_s_32, is_8, ip_13, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * is_s_8[k]
                  + f_2 * id_s_30[k]
                  + f_3 * is_8[k]
                  + pb_x[k] * ip_13[k];

        t_28[k] = f_2 * id_s_31[k]
                  + pb_x[k] * ip_14[k];

        t_29[k] = f_0 * hp_11[k]
                  - f_1 * is_s_8[k]
                  + f_2 * id_s_32[k]
                  + f_3 * is_8[k]
                  + pb_y[k] * ip_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, pb_z, gd_s_18, gd_18, hd_29, is_s_8, \
                         id_s_33, id_s_34, id_s_37, is_8, ip_15, \
                         ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_1 * is_s_8[k]
                  + f_2 * id_s_33[k]
                  + f_3 * is_8[k]
                  + pb_z[k] * ip_15[k];

        t_31[k] = f_2 * id_s_34[k]
                  + pb_x[k] * ip_16[k];

        t_32[k] = -f_4 * gd_s_18[k]
                  + f_5 * gd_18[k]
                  + pa_y[k] * hd_29[k]
                  + f_2 * id_s_37[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, is_s_10, id_s_38, id_s_39, id_s_40, is_10, \
                         ip_17, ip_18, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -f_1 * is_s_10[k]
                  + f_2 * id_s_38[k]
                  + f_3 * is_10[k]
                  + pb_x[k] * ip_17[k];

        t_34[k] = f_2 * id_s_39[k]
                  + pb_x[k] * ip_18[k];

        t_35[k] = f_2 * id_s_40[k]
                  + pb_x[k] * ip_19[k];
    }

#pragma omp simd aligned(t_36, t_37, pa_z, pb_y, gd_s_15, gd_15, hp_14, hd_28, id_s_41, \
                         id_s_42, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_6 * gd_s_15[k]
                  + f_3 * gd_15[k]
                  + pa_z[k] * hd_28[k]
                  + f_2 * id_s_41[k];

        t_37[k] = f_5 * hp_14[k]
                  + f_2 * id_s_42[k]
                  + pb_y[k] * ip_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pb_x, gd_s_21, gd_21, hd_32, is_s_11, \
                         id_s_43, id_s_44, id_s_45, is_11, ip_20, \
                         ip_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -f_7 * gd_s_21[k]
                  + f_8 * gd_21[k]
                  + pa_y[k] * hd_32[k]
                  + f_2 * id_s_43[k];

        t_39[k] = -f_1 * is_s_11[k]
                  + f_2 * id_s_44[k]
                  + f_3 * is_11[k]
                  + pb_x[k] * ip_20[k];

        t_40[k] = f_2 * id_s_45[k]
                  + pb_x[k] * ip_21[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_z, pb_x, pb_y, gd_s_17, gd_17, hp_16, hd_31, \
                         id_s_46, id_s_47, id_s_48, ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_2 * id_s_46[k]
                  + pb_x[k] * ip_22[k];

        t_42[k] = -f_9 * gd_s_17[k]
                  + f_10 * gd_17[k]
                  + pa_z[k] * hd_31[k]
                  + f_2 * id_s_47[k];

        t_43[k] = f_8 * hp_16[k]
                  + f_2 * id_s_48[k]
                  + pb_y[k] * ip_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pb_x, gd_s_23, gd_23, hd_35, is_s_12, \
                         id_s_49, id_s_50, id_s_51, is_12, ip_23, \
                         ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = -f_9 * gd_s_23[k]
                  + f_10 * gd_23[k]
                  + pa_y[k] * hd_35[k]
                  + f_2 * id_s_49[k];

        t_45[k] = -f_1 * is_s_12[k]
                  + f_2 * id_s_50[k]
                  + f_3 * is_12[k]
                  + pb_x[k] * ip_23[k];

        t_46[k] = f_2 * id_s_51[k]
                  + pb_x[k] * ip_24[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_z, pb_x, pb_y, gd_s_20, gd_20, hp_17, hd_34, \
                         id_s_52, id_s_53, id_s_54, ip_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_2 * id_s_52[k]
                  + pb_x[k] * ip_25[k];

        t_48[k] = -f_7 * gd_s_20[k]
                  + f_8 * gd_20[k]
                  + pa_z[k] * hd_34[k]
                  + f_2 * id_s_53[k];

        t_49[k] = f_10 * hp_17[k]
                  + f_2 * id_s_54[k]
                  + pb_y[k] * ip_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_y, pb_x, gd_s_27, gd_27, hd_37, is_s_14, \
                         id_s_55, id_s_60, id_s_61, is_13, ip_26, \
                         ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_6 * gd_s_27[k]
                  + f_3 * gd_27[k]
                  + pa_y[k] * hd_37[k]
                  + f_2 * id_s_55[k];

        t_51[k] = -f_1 * is_s_14[k]
                  + f_2 * id_s_60[k]
                  + f_3 * is_13[k]
                  + pb_x[k] * ip_26[k];

        t_52[k] = f_2 * id_s_61[k]
                  + pb_x[k] * ip_28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pb_z, hp_20, is_s_14, id_s_62, id_s_63, \
                         id_s_64, is_13, ip_27, ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_1 * is_s_14[k]
                  + f_2 * id_s_62[k]
                  + f_3 * is_13[k]
                  + pb_y[k] * ip_27[k];

        t_54[k] = f_2 * id_s_63[k]
                  + pb_y[k] * ip_28[k];

        t_55[k] = f_0 * hp_20[k]
                  - f_1 * is_s_14[k]
                  + f_2 * id_s_64[k]
                  + f_3 * is_13[k]
                  + pb_z[k] * ip_28[k];
    }
}

auto
compute_prim_id_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gd_s, const size_t gd,
                                 const size_t hp, const size_t hd, const size_t is_s,
                                 const size_t id_s, const size_t is, const size_t ip,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 4.0 * beta / p;
    const auto f_5 = 2.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 3.0 * beta / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 2.0 * beta / p;
    const auto f_10 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_4 = buffer.data(gd_s + 4);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_6 = buffer.data(gd_s + 6);
    const auto *gd_s_7 = buffer.data(gd_s + 7);
    const auto *gd_s_8 = buffer.data(gd_s + 8);
    const auto *gd_s_10 = buffer.data(gd_s + 10);
    const auto *gd_s_11 = buffer.data(gd_s + 11);
    const auto *gd_s_12 = buffer.data(gd_s + 12);
    const auto *gd_s_13 = buffer.data(gd_s + 13);
    const auto *gd_s_15 = buffer.data(gd_s + 15);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_20 = buffer.data(gd_s + 20);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_27 = buffer.data(gd_s + 27);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_26 = buffer.data(gd + 26);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_20 = buffer.data(hp + 20);

    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);

    const auto *is_s_0 = buffer.data(is_s + 0);
    const auto *is_s_2 = buffer.data(is_s + 2);
    const auto *is_s_3 = buffer.data(is_s + 3);
    const auto *is_s_4 = buffer.data(is_s + 4);
    const auto *is_s_5 = buffer.data(is_s + 5);
    const auto *is_s_6 = buffer.data(is_s + 6);
    const auto *is_s_7 = buffer.data(is_s + 7);
    const auto *is_s_8 = buffer.data(is_s + 8);
    const auto *is_s_10 = buffer.data(is_s + 10);
    const auto *is_s_11 = buffer.data(is_s + 11);
    const auto *is_s_12 = buffer.data(is_s + 12);
    const auto *is_s_14 = buffer.data(is_s + 14);

    const auto *id_s_0 = buffer.data(id_s + 0);
    const auto *id_s_1 = buffer.data(id_s + 1);
    const auto *id_s_2 = buffer.data(id_s + 2);
    const auto *id_s_3 = buffer.data(id_s + 3);
    const auto *id_s_5 = buffer.data(id_s + 5);
    const auto *id_s_6 = buffer.data(id_s + 6);
    const auto *id_s_7 = buffer.data(id_s + 7);
    const auto *id_s_8 = buffer.data(id_s + 8);
    const auto *id_s_9 = buffer.data(id_s + 9);
    const auto *id_s_10 = buffer.data(id_s + 10);
    const auto *id_s_11 = buffer.data(id_s + 11);
    const auto *id_s_12 = buffer.data(id_s + 12);
    const auto *id_s_13 = buffer.data(id_s + 13);
    const auto *id_s_14 = buffer.data(id_s + 14);
    const auto *id_s_15 = buffer.data(id_s + 15);
    const auto *id_s_16 = buffer.data(id_s + 16);
    const auto *id_s_17 = buffer.data(id_s + 17);
    const auto *id_s_18 = buffer.data(id_s + 18);
    const auto *id_s_19 = buffer.data(id_s + 19);
    const auto *id_s_20 = buffer.data(id_s + 20);
    const auto *id_s_21 = buffer.data(id_s + 21);
    const auto *id_s_22 = buffer.data(id_s + 22);
    const auto *id_s_23 = buffer.data(id_s + 23);
    const auto *id_s_24 = buffer.data(id_s + 24);
    const auto *id_s_25 = buffer.data(id_s + 25);
    const auto *id_s_26 = buffer.data(id_s + 26);
    const auto *id_s_27 = buffer.data(id_s + 27);
    const auto *id_s_30 = buffer.data(id_s + 30);
    const auto *id_s_31 = buffer.data(id_s + 31);
    const auto *id_s_32 = buffer.data(id_s + 32);
    const auto *id_s_33 = buffer.data(id_s + 33);
    const auto *id_s_34 = buffer.data(id_s + 34);
    const auto *id_s_37 = buffer.data(id_s + 37);
    const auto *id_s_38 = buffer.data(id_s + 38);
    const auto *id_s_39 = buffer.data(id_s + 39);
    const auto *id_s_40 = buffer.data(id_s + 40);
    const auto *id_s_41 = buffer.data(id_s + 41);
    const auto *id_s_42 = buffer.data(id_s + 42);
    const auto *id_s_43 = buffer.data(id_s + 43);
    const auto *id_s_44 = buffer.data(id_s + 44);
    const auto *id_s_45 = buffer.data(id_s + 45);
    const auto *id_s_46 = buffer.data(id_s + 46);
    const auto *id_s_47 = buffer.data(id_s + 47);
    const auto *id_s_48 = buffer.data(id_s + 48);
    const auto *id_s_49 = buffer.data(id_s + 49);
    const auto *id_s_50 = buffer.data(id_s + 50);
    const auto *id_s_51 = buffer.data(id_s + 51);
    const auto *id_s_52 = buffer.data(id_s + 52);
    const auto *id_s_53 = buffer.data(id_s + 53);
    const auto *id_s_54 = buffer.data(id_s + 54);
    const auto *id_s_55 = buffer.data(id_s + 55);
    const auto *id_s_60 = buffer.data(id_s + 60);
    const auto *id_s_61 = buffer.data(id_s + 61);
    const auto *id_s_62 = buffer.data(id_s + 62);
    const auto *id_s_63 = buffer.data(id_s + 63);
    const auto *id_s_64 = buffer.data(id_s + 64);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is_s_0, id_s_0, id_s_1, \
                         id_s_2, is_0, ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 - f_1 * is_s_0[k]
                 + f_2 * id_s_0[k]
                 + f_3 * is_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_1[k]
                 + f_3 * is_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = -f_1 * is_s_0[k]
                 + f_2 * id_s_2[k]
                 + f_3 * is_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, gd_s_4, gd_s_6, gd_4, gd_6, hd_4, hd_7, \
                         id_s_3, id_s_5, id_s_6, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * gd_s_4[k]
                 + f_5 * gd_4[k]
                 + pa_x[k] * hd_4[k]
                 + f_2 * id_s_3[k];

        t_4[k] = f_2 * id_s_5[k]
                 + pb_y[k] * ip_3[k];

        t_5[k] = -f_4 * gd_s_6[k]
                 + f_5 * gd_6[k]
                 + pa_x[k] * hd_7[k]
                 + f_2 * id_s_6[k];
    }

#pragma omp simd aligned(t_6, t_7, pa_x, pa_y, gd_s_0, gd_s_8, gd_0, gd_8, hd_3, hd_9, id_s_7, \
                         id_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_6 * gd_s_0[k]
                 + f_3 * gd_0[k]
                 + pa_y[k] * hd_3[k]
                 + f_2 * id_s_7[k];

        t_7[k] = -f_7 * gd_s_8[k]
                 + f_8 * gd_8[k]
                 + pa_x[k] * hd_9[k]
                 + f_2 * id_s_8[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_z, pb_z, gd_s_0, gd_0, hd_5, is_s_2, id_s_9, id_s_10, \
                         is_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_1 * is_s_2[k]
                 + f_2 * id_s_9[k]
                 + f_3 * is_2[k]
                 + pb_z[k] * ip_4[k];

        t_9[k] = -f_6 * gd_s_0[k]
                 + f_3 * gd_0[k]
                 + pa_z[k] * hd_5[k]
                 + f_2 * id_s_10[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_y, gd_s_11, gd_11, hd_14, is_s_3, id_s_11, \
                         id_s_12, id_s_13, is_3, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * is_s_3[k]
                  + f_2 * id_s_11[k]
                  + f_3 * is_3[k]
                  + pb_y[k] * ip_5[k];

        t_11[k] = f_2 * id_s_12[k]
                  + pb_y[k] * ip_6[k];

        t_12[k] = -f_7 * gd_s_11[k]
                  + f_8 * gd_11[k]
                  + pa_x[k] * hd_14[k]
                  + f_2 * id_s_13[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, gd_s_3, gd_s_12, gd_3, gd_12, hd_8, hd_16, \
                         id_s_14, id_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_9 * gd_s_3[k]
                  + f_10 * gd_3[k]
                  + pa_y[k] * hd_8[k]
                  + f_2 * id_s_14[k];

        t_14[k] = -f_9 * gd_s_12[k]
                  + f_10 * gd_12[k]
                  + pa_x[k] * hd_16[k]
                  + f_2 * id_s_15[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_z, pb_z, gd_s_5, gd_5, hd_11, is_s_4, id_s_16, \
                         id_s_17, is_4, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * is_s_4[k]
                  + f_2 * id_s_16[k]
                  + f_3 * is_4[k]
                  + pb_z[k] * ip_7[k];

        t_16[k] = -f_9 * gd_s_5[k]
                  + f_10 * gd_5[k]
                  + pa_z[k] * hd_11[k]
                  + f_2 * id_s_17[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pb_y, gd_s_13, gd_13, hd_21, is_s_5, id_s_18, \
                         id_s_19, id_s_20, is_5, ip_8, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_1 * is_s_5[k]
                  + f_2 * id_s_18[k]
                  + f_3 * is_5[k]
                  + pb_y[k] * ip_8[k];

        t_18[k] = f_2 * id_s_19[k]
                  + pb_y[k] * ip_9[k];

        t_19[k] = -f_9 * gd_s_13[k]
                  + f_10 * gd_13[k]
                  + pa_x[k] * hd_21[k]
                  + f_2 * id_s_20[k];
    }

#pragma omp simd aligned(t_20, t_21, pa_x, pa_y, gd_s_7, gd_s_15, gd_7, gd_15, hd_15, hd_22, \
                         id_s_21, id_s_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_7 * gd_s_7[k]
                  + f_8 * gd_7[k]
                  + pa_y[k] * hd_15[k]
                  + f_2 * id_s_21[k];

        t_21[k] = -f_6 * gd_s_15[k]
                  + f_3 * gd_15[k]
                  + pa_x[k] * hd_22[k]
                  + f_2 * id_s_22[k];
    }

#pragma omp simd aligned(t_22, t_23, pa_z, pb_z, gd_s_10, gd_10, hd_18, is_s_6, id_s_23, \
                         id_s_24, is_6, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_1 * is_s_6[k]
                  + f_2 * id_s_23[k]
                  + f_3 * is_6[k]
                  + pb_z[k] * ip_10[k];

        t_23[k] = -f_7 * gd_s_10[k]
                  + f_8 * gd_10[k]
                  + pa_z[k] * hd_18[k]
                  + f_2 * id_s_24[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_y, gd_s_27, gd_26, hd_23, is_s_7, id_s_25, \
                         id_s_26, id_s_27, is_7, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * is_s_7[k]
                  + f_2 * id_s_25[k]
                  + f_3 * is_7[k]
                  + pb_y[k] * ip_11[k];

        t_25[k] = f_2 * id_s_26[k]
                  + pb_y[k] * ip_12[k];

        t_26[k] = -f_6 * gd_s_27[k]
                  + f_3 * gd_26[k]
                  + pa_x[k] * hd_23[k]
                  + f_2 * id_s_27[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, pb_y, hp_11, is_s_8, id_s_30, id_s_31, \
                         id_s_32, is_8, ip_13, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * is_s_8[k]
                  + f_2 * id_s_30[k]
                  + f_3 * is_8[k]
                  + pb_x[k] * ip_13[k];

        t_28[k] = f_2 * id_s_31[k]
                  + pb_x[k] * ip_14[k];

        t_29[k] = f_0 * hp_11[k]
                  - f_1 * is_s_8[k]
                  + f_2 * id_s_32[k]
                  + f_3 * is_8[k]
                  + pb_y[k] * ip_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, pb_z, gd_s_18, gd_18, hd_29, is_s_8, \
                         id_s_33, id_s_34, id_s_37, is_8, ip_15, \
                         ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_1 * is_s_8[k]
                  + f_2 * id_s_33[k]
                  + f_3 * is_8[k]
                  + pb_z[k] * ip_15[k];

        t_31[k] = f_2 * id_s_34[k]
                  + pb_x[k] * ip_16[k];

        t_32[k] = -f_4 * gd_s_18[k]
                  + f_5 * gd_18[k]
                  + pa_y[k] * hd_29[k]
                  + f_2 * id_s_37[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, is_s_10, id_s_38, id_s_39, id_s_40, is_10, \
                         ip_17, ip_18, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -f_1 * is_s_10[k]
                  + f_2 * id_s_38[k]
                  + f_3 * is_10[k]
                  + pb_x[k] * ip_17[k];

        t_34[k] = f_2 * id_s_39[k]
                  + pb_x[k] * ip_18[k];

        t_35[k] = f_2 * id_s_40[k]
                  + pb_x[k] * ip_19[k];
    }

#pragma omp simd aligned(t_36, t_37, pa_z, pb_y, gd_s_15, gd_15, hp_14, hd_28, id_s_41, \
                         id_s_42, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_6 * gd_s_15[k]
                  + f_3 * gd_15[k]
                  + pa_z[k] * hd_28[k]
                  + f_2 * id_s_41[k];

        t_37[k] = f_5 * hp_14[k]
                  + f_2 * id_s_42[k]
                  + pb_y[k] * ip_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pb_x, gd_s_21, gd_21, hd_32, is_s_11, \
                         id_s_43, id_s_44, id_s_45, is_11, ip_20, \
                         ip_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -f_7 * gd_s_21[k]
                  + f_8 * gd_21[k]
                  + pa_y[k] * hd_32[k]
                  + f_2 * id_s_43[k];

        t_39[k] = -f_1 * is_s_11[k]
                  + f_2 * id_s_44[k]
                  + f_3 * is_11[k]
                  + pb_x[k] * ip_20[k];

        t_40[k] = f_2 * id_s_45[k]
                  + pb_x[k] * ip_21[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_z, pb_x, pb_y, gd_s_17, gd_17, hp_16, hd_31, \
                         id_s_46, id_s_47, id_s_48, ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_2 * id_s_46[k]
                  + pb_x[k] * ip_22[k];

        t_42[k] = -f_9 * gd_s_17[k]
                  + f_10 * gd_17[k]
                  + pa_z[k] * hd_31[k]
                  + f_2 * id_s_47[k];

        t_43[k] = f_8 * hp_16[k]
                  + f_2 * id_s_48[k]
                  + pb_y[k] * ip_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pb_x, gd_s_23, gd_22, hd_35, is_s_12, \
                         id_s_49, id_s_50, id_s_51, is_12, ip_23, \
                         ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = -f_9 * gd_s_23[k]
                  + f_10 * gd_22[k]
                  + pa_y[k] * hd_35[k]
                  + f_2 * id_s_49[k];

        t_45[k] = -f_1 * is_s_12[k]
                  + f_2 * id_s_50[k]
                  + f_3 * is_12[k]
                  + pb_x[k] * ip_23[k];

        t_46[k] = f_2 * id_s_51[k]
                  + pb_x[k] * ip_24[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_z, pb_x, pb_y, gd_s_20, gd_20, hp_17, hd_34, \
                         id_s_52, id_s_53, id_s_54, ip_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_2 * id_s_52[k]
                  + pb_x[k] * ip_25[k];

        t_48[k] = -f_7 * gd_s_20[k]
                  + f_8 * gd_20[k]
                  + pa_z[k] * hd_34[k]
                  + f_2 * id_s_53[k];

        t_49[k] = f_10 * hp_17[k]
                  + f_2 * id_s_54[k]
                  + pb_y[k] * ip_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_y, pb_x, gd_s_27, gd_26, hd_36, is_s_14, \
                         id_s_55, id_s_60, id_s_61, is_13, ip_26, \
                         ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_6 * gd_s_27[k]
                  + f_3 * gd_26[k]
                  + pa_y[k] * hd_36[k]
                  + f_2 * id_s_55[k];

        t_51[k] = -f_1 * is_s_14[k]
                  + f_2 * id_s_60[k]
                  + f_3 * is_13[k]
                  + pb_x[k] * ip_26[k];

        t_52[k] = f_2 * id_s_61[k]
                  + pb_x[k] * ip_28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pb_z, hp_20, is_s_14, id_s_62, id_s_63, \
                         id_s_64, is_13, ip_27, ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_1 * is_s_14[k]
                  + f_2 * id_s_62[k]
                  + f_3 * is_13[k]
                  + pb_y[k] * ip_27[k];

        t_54[k] = f_2 * id_s_63[k]
                  + pb_y[k] * ip_28[k];

        t_55[k] = f_0 * hp_20[k]
                  - f_1 * is_s_14[k]
                  + f_2 * id_s_64[k]
                  + f_3 * is_13[k]
                  + pb_z[k] * ip_28[k];
    }
}

}  // namespace simdkin
