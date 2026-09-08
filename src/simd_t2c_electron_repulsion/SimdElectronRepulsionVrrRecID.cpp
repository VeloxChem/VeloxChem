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


#include "SimdElectronRepulsionVrrRecID.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_id_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 2.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.0 / p;
    const auto f_8 = 1.5 / alpha;
    const auto f_9 = 1.5 * beta / (alpha * p);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 1.5 / p;

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

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
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
    const auto *hp_21 = buffer.data(hp + 21);
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
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, hp_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = pb_y[k] * ip_0[k];

        t_2[k] = pb_z[k] * ip_0[k];

        t_3[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_4[k] = pb_y[k] * ip_2[k];

        t_5[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_y, pb_x, pb_z, hp_1, hp_3, hd_0, \
                         hd_1, hd_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * hd_0[k];

        t_7[k] = f_3 * hp_3[k]
                 + pb_x[k] * ip_4[k];

        t_8[k] = pb_z[k] * ip_3[k];

        t_9[k] = f_4 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_10[k] = pb_z[k] * ip_4[k];

        t_11[k] = pa_y[k] * hd_2[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, hp_2, hp_4, \
                         hd_0, hd_1, hd_2, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * hd_0[k];

        t_13[k] = pb_y[k] * ip_5[k];

        t_14[k] = f_3 * hp_4[k]
                  + pb_x[k] * ip_6[k];

        t_15[k] = pa_z[k] * hd_1[k];

        t_16[k] = pb_y[k] * ip_6[k];

        t_17[k] = f_4 * hp_2[k]
                  + pa_z[k] * hd_2[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_x, pb_z, gd0_0, gd1_0, hp_5, hd_3, ip_7, \
                         ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_0[k]
                  - f_6 * gd1_0[k]
                  + pa_y[k] * hd_3[k];

        t_19[k] = f_7 * hp_5[k]
                  + pb_x[k] * ip_8[k];

        t_20[k] = pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_y, pb_z, gd0_4, gd1_4, hd_6, hd_11, \
                         is0_1, is1_1, ip_8, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * gd0_4[k]
                  - f_9 * gd1_4[k]
                  + pa_x[k] * hd_11[k];

        t_22[k] = pb_z[k] * ip_8[k];

        t_23[k] = f_1 * is0_1[k]
                  - f_2 * is1_1[k]
                  + pb_z[k] * ip_9[k];

        t_24[k] = pa_y[k] * hd_6[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, hp_4, hd_4, hd_5, \
                         hd_7, hd_8, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * hd_4[k];

        t_26[k] = pa_y[k] * hd_7[k];

        t_27[k] = pa_z[k] * hd_5[k];

        t_28[k] = f_10 * hp_4[k]
                  + pb_y[k] * ip_10[k];

        t_29[k] = pa_y[k] * hd_8[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pb_x, pb_y, gd0_0, gd1_0, hp_9, hd_6, \
                         is0_2, is1_2, ip_11, ip_12, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * gd0_0[k]
                  - f_6 * gd1_0[k]
                  + pa_z[k] * hd_6[k];

        t_31[k] = pb_y[k] * ip_11[k];

        t_32[k] = f_7 * hp_9[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_12[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_y, pb_y, gd0_1, gd0_6, gd1_1, gd1_6, hd_9, \
                         hd_16, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * ip_13[k];

        t_35[k] = f_8 * gd0_6[k]
                  - f_9 * gd1_6[k]
                  + pa_x[k] * hd_16[k];

        t_36[k] = f_11 * gd0_1[k]
                  - f_12 * gd1_1[k]
                  + pa_y[k] * hd_9[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pb_x, pb_z, gd0_7, gd1_7, hp_10, hd_19, \
                         ip_14, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_13 * hp_10[k]
                  + pb_x[k] * ip_15[k];

        t_38[k] = pb_z[k] * ip_14[k];

        t_39[k] = f_11 * gd0_7[k]
                  - f_12 * gd1_7[k]
                  + pa_x[k] * hd_19[k];

        t_40[k] = pb_z[k] * ip_15[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_z, pb_x, pb_z, hp_12, hd_9, hd_10, \
                         hd_11, is0_3, is1_3, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_16[k];

        t_42[k] = pa_z[k] * hd_9[k];

        t_43[k] = pa_z[k] * hd_10[k];

        t_44[k] = f_13 * hp_12[k]
                  + pb_x[k] * ip_17[k];

        t_45[k] = pa_z[k] * hd_11[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_x, pb_y, hp_6, hp_7, hp_13, \
                         hd_12, hd_13, ip_17, ip_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * hp_7[k]
                  + pb_y[k] * ip_17[k];

        t_47[k] = f_4 * hp_6[k]
                  + pa_z[k] * hd_12[k];

        t_48[k] = pa_y[k] * hd_13[k];

        t_49[k] = f_13 * hp_13[k]
                  + pb_x[k] * ip_18[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_y, hp_8, hp_9, hd_14, hd_15, hd_16, \
                         ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * hd_14[k];

        t_51[k] = f_4 * hp_8[k]
                  + pa_y[k] * hd_15[k];

        t_52[k] = f_10 * hp_9[k]
                  + pb_y[k] * ip_19[k];

        t_53[k] = pa_y[k] * hd_16[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, gd0_2, gd1_2, hp_16, hd_13, \
                         is0_4, is1_4, ip_20, ip_21, ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * gd0_2[k]
                  - f_12 * gd1_2[k]
                  + pa_z[k] * hd_13[k];

        t_55[k] = pb_y[k] * ip_20[k];

        t_56[k] = f_13 * hp_16[k]
                  + pb_x[k] * ip_22[k];

        t_57[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_21[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pa_y, pb_y, gd0_3, gd0_8, gd1_3, gd1_8, \
                         hd_17, hd_25, ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_y[k] * ip_22[k];

        t_59[k] = f_11 * gd0_8[k]
                  - f_12 * gd1_8[k]
                  + pa_x[k] * hd_25[k];

        t_60[k] = f_8 * gd0_3[k]
                  - f_9 * gd1_3[k]
                  + pa_y[k] * hd_17[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_x, pb_z, gd0_9, gd1_9, hp_17, hd_28, \
                         ip_23, ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_4 * hp_17[k]
                  + pb_x[k] * ip_24[k];

        t_62[k] = pb_z[k] * ip_23[k];

        t_63[k] = f_5 * gd0_9[k]
                  - f_6 * gd1_9[k]
                  + pa_x[k] * hd_28[k];

        t_64[k] = pb_z[k] * ip_24[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_z, pb_x, pb_z, hp_18, hd_17, hd_18, \
                         hd_19, is0_5, is1_5, ip_25, ip_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_25[k];

        t_66[k] = pa_z[k] * hd_17[k];

        t_67[k] = pa_z[k] * hd_18[k];

        t_68[k] = f_4 * hp_18[k]
                  + pb_x[k] * ip_26[k];

        t_69[k] = pa_z[k] * hd_19[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pa_z, pb_y, gd0_5, gd1_5, hp_11, hp_12, \
                         hd_20, hd_21, ip_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * hp_12[k]
                  + pb_y[k] * ip_26[k];

        t_71[k] = f_4 * hp_11[k]
                  + pa_z[k] * hd_20[k];

        t_72[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_21[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_x, pb_y, gd0_11, gd1_11, hp_14, \
                         hp_19, hp_20, hd_29, ip_27, ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_4 * hp_19[k]
                  + pb_x[k] * ip_27[k];

        t_74[k] = f_4 * hp_20[k]
                  + pb_x[k] * ip_28[k];

        t_75[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_x[k] * hd_29[k];

        t_76[k] = f_4 * hp_14[k]
                  + pb_y[k] * ip_28[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_x, pa_y, pb_x, gd0_12, gd1_12, hp_21, \
                         hd_22, hd_23, hd_30, ip_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_x[k] * hd_30[k];

        t_78[k] = pa_y[k] * hd_22[k];

        t_79[k] = f_4 * hp_21[k]
                  + pb_x[k] * ip_29[k];

        t_80[k] = pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pa_z, pb_y, gd0_5, gd1_5, hp_15, hp_16, \
                         hd_22, hd_24, hd_25, ip_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * hp_15[k]
                  + pa_y[k] * hd_24[k];

        t_82[k] = f_10 * hp_16[k]
                  + pb_y[k] * ip_30[k];

        t_83[k] = pa_y[k] * hd_25[k];

        t_84[k] = f_8 * gd0_5[k]
                  - f_9 * gd1_5[k]
                  + pa_z[k] * hd_22[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, pb_y, hp_22, is0_6, is1_6, ip_31, \
                         ip_32, ip_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_y[k] * ip_31[k];

        t_86[k] = f_4 * hp_22[k]
                  + pb_x[k] * ip_33[k];

        t_87[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_32[k];

        t_88[k] = pb_y[k] * ip_33[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_x, pb_x, pb_z, gd0_14, gd1_14, hp_23, \
                         hp_24, hd_33, hd_34, ip_34, ip_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_5 * gd0_14[k]
                  - f_6 * gd1_14[k]
                  + pa_x[k] * hd_33[k];

        t_90[k] = f_4 * hp_23[k]
                  + pa_x[k] * hd_34[k];

        t_91[k] = f_10 * hp_24[k]
                  + pb_x[k] * ip_35[k];

        t_92[k] = pb_z[k] * ip_34[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_x, pa_z, pb_z, hd_26, hd_27, hd_35, \
                         hd_36, ip_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_x[k] * hd_35[k];

        t_94[k] = pb_z[k] * ip_35[k];

        t_95[k] = pa_x[k] * hd_36[k];

        t_96[k] = pa_z[k] * hd_26[k];

        t_97[k] = pa_z[k] * hd_27[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_x, pb_x, hp_26, hp_27, hd_37, \
                         hd_38, hd_39, hd_40, ip_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_10 * hp_26[k]
                  + pb_x[k] * ip_36[k];

        t_99[k] = pa_x[k] * hd_37[k];

        t_100[k] = pa_x[k] * hd_38[k];

        t_101[k] = pa_x[k] * hd_39[k];

        t_102[k] = f_4 * hp_27[k]
                   + pa_x[k] * hd_40[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_x, pb_x, hp_28, hp_29, hd_41, \
                         hd_42, hd_43, ip_37, ip_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_10 * hp_28[k]
                   + pb_x[k] * ip_37[k];

        t_104[k] = f_10 * hp_29[k]
                   + pb_x[k] * ip_38[k];

        t_105[k] = pa_x[k] * hd_41[k];

        t_106[k] = pa_x[k] * hd_42[k];

        t_107[k] = pa_x[k] * hd_43[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_x, pb_x, hp_30, hp_31, hp_32, \
                         hd_44, hd_45, hd_46, ip_39, ip_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_4 * hp_30[k]
                   + pa_x[k] * hd_44[k];

        t_109[k] = f_10 * hp_31[k]
                   + pb_x[k] * ip_39[k];

        t_110[k] = f_10 * hp_32[k]
                   + pb_x[k] * ip_40[k];

        t_111[k] = pa_x[k] * hd_45[k];

        t_112[k] = pa_x[k] * hd_46[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, t_118, pa_x, pa_y, pb_x, hp_33, \
                         hd_31, hd_32, hd_47, hd_48, hd_49, ip_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_x[k] * hd_47[k];

        t_114[k] = pa_y[k] * hd_31[k];

        t_115[k] = f_10 * hp_33[k]
                   + pb_x[k] * ip_41[k];

        t_116[k] = pa_y[k] * hd_32[k];

        t_117[k] = pa_x[k] * hd_48[k];

        t_118[k] = pa_x[k] * hd_49[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, t_124, pa_x, pb_x, pb_y, hp_35, \
                         hp_37, hd_50, hd_51, hd_52, ip_42, ip_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pa_x[k] * hd_50[k];

        t_120[k] = f_4 * hp_35[k]
                   + pa_x[k] * hd_51[k];

        t_121[k] = pb_y[k] * ip_42[k];

        t_122[k] = f_10 * hp_37[k]
                   + pb_x[k] * ip_43[k];

        t_123[k] = pa_x[k] * hd_52[k];

        t_124[k] = pb_y[k] * ip_43[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pa_x, pb_x, pb_y, hp_24, hd_53, \
                         is0_7, is1_7, ip_44, ip_45, ip_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_x[k] * hd_53[k];

        t_126[k] = f_1 * is0_7[k]
                   - f_2 * is1_7[k]
                   + pb_x[k] * ip_44[k];

        t_127[k] = pb_x[k] * ip_45[k];

        t_128[k] = pb_x[k] * ip_46[k];

        t_129[k] = f_0 * hp_24[k]
                   + f_1 * is0_7[k]
                   - f_2 * is1_7[k]
                   + pb_y[k] * ip_45[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pa_z, pb_x, pb_z, hd_34, is0_7, \
                         is1_7, ip_45, ip_46, ip_47, ip_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_z[k] * ip_45[k];

        t_131[k] = f_1 * is0_7[k]
                   - f_2 * is1_7[k]
                   + pb_z[k] * ip_46[k];

        t_132[k] = pa_z[k] * hd_34[k];

        t_133[k] = pb_x[k] * ip_47[k];

        t_134[k] = pb_x[k] * ip_48[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_x, pb_y, hp_25, hp_26, hd_35, \
                         hd_36, is0_8, is1_8, ip_48, ip_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = pa_z[k] * hd_35[k];

        t_136[k] = f_3 * hp_26[k]
                   + pb_y[k] * ip_48[k];

        t_137[k] = f_4 * hp_25[k]
                   + pa_z[k] * hd_36[k];

        t_138[k] = f_1 * is0_8[k]
                   - f_2 * is1_8[k]
                   + pb_x[k] * ip_49[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_z, pb_x, pb_y, gd0_9, gd1_9, hp_29, \
                         hd_37, ip_50, ip_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_x[k] * ip_50[k];

        t_140[k] = pb_x[k] * ip_51[k];

        t_141[k] = f_5 * gd0_9[k]
                   - f_6 * gd1_9[k]
                   + pa_z[k] * hd_37[k];

        t_142[k] = f_7 * hp_29[k]
                   + pb_y[k] * ip_51[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_y, pb_x, gd0_12, gd1_12, hd_43, is0_9, \
                         is1_9, ip_52, ip_53, ip_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_8 * gd0_12[k]
                   - f_9 * gd1_12[k]
                   + pa_y[k] * hd_43[k];

        t_144[k] = f_1 * is0_9[k]
                   - f_2 * is1_9[k]
                   + pb_x[k] * ip_52[k];

        t_145[k] = pb_x[k] * ip_53[k];

        t_146[k] = pb_x[k] * ip_54[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_y, pa_z, pb_y, gd0_10, gd0_13, gd1_10, \
                         gd1_13, hp_32, hd_41, hd_47, ip_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_11 * gd0_10[k]
                   - f_12 * gd1_10[k]
                   + pa_z[k] * hd_41[k];

        t_148[k] = f_13 * hp_32[k]
                   + pb_y[k] * ip_54[k];

        t_149[k] = f_11 * gd0_13[k]
                   - f_12 * gd1_13[k]
                   + pa_y[k] * hd_47[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_z, pb_x, gd0_11, gd1_11, hd_45, \
                         is0_10, is1_10, ip_55, ip_56, ip_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * is0_10[k]
                   - f_2 * is1_10[k]
                   + pb_x[k] * ip_55[k];

        t_151[k] = pb_x[k] * ip_56[k];

        t_152[k] = pb_x[k] * ip_57[k];

        t_153[k] = f_8 * gd0_11[k]
                   - f_9 * gd1_11[k]
                   + pa_z[k] * hd_45[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, pa_y, pb_x, pb_y, gd0_14, gd1_14, \
                         hp_34, hd_50, hd_51, ip_57, ip_58, ip_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_4 * hp_34[k]
                   + pb_y[k] * ip_57[k];

        t_155[k] = f_5 * gd0_14[k]
                   - f_6 * gd1_14[k]
                   + pa_y[k] * hd_50[k];

        t_156[k] = pa_y[k] * hd_51[k];

        t_157[k] = pb_x[k] * ip_58[k];

        t_158[k] = pb_x[k] * ip_59[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_y, pb_x, pb_y, hp_36, hp_37, hd_52, \
                         hd_53, is0_11, is1_11, ip_59, ip_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_4 * hp_36[k]
                   + pa_y[k] * hd_52[k];

        t_160[k] = f_10 * hp_37[k]
                   + pb_y[k] * ip_59[k];

        t_161[k] = pa_y[k] * hd_53[k];

        t_162[k] = f_1 * is0_11[k]
                   - f_2 * is1_11[k]
                   + pb_x[k] * ip_60[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pb_x, pb_y, pb_z, hp_37, is0_11, \
                         is1_11, ip_61, ip_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pb_x[k] * ip_61[k];

        t_164[k] = pb_x[k] * ip_62[k];

        t_165[k] = f_1 * is0_11[k]
                   - f_2 * is1_11[k]
                   + pb_y[k] * ip_61[k];

        t_166[k] = pb_y[k] * ip_62[k];

        t_167[k] = f_0 * hp_37[k]
                   + f_1 * is0_11[k]
                   - f_2 * is1_11[k]
                   + pb_z[k] * ip_62[k];
    }
}

auto
compute_prim_id_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 2.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.0 / p;
    const auto f_8 = 1.5 / alpha;
    const auto f_9 = 1.5 * beta / (alpha * p);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_24 = buffer.data(gd1 + 24);
    const auto *gd1_26 = buffer.data(gd1 + 26);
    const auto *gd1_29 = buffer.data(gd1 + 29);
    const auto *gd1_32 = buffer.data(gd1 + 32);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, \
                         is1_0, ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = pb_z[k] * ip_0[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_3[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_4[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_x, hp_1, hp_3, hp_4, hd_0, \
                         hd_1, hd_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * hp_3[k]
                 + pb_x[k] * ip_3[k];

        t_6[k] = f_4 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_7[k] = pa_y[k] * hd_2[k];

        t_8[k] = pa_z[k] * hd_0[k];

        t_9[k] = f_3 * hp_4[k]
                 + pb_x[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_x, gd0_0, gd1_0, hp_2, hp_5, \
                         hd_1, hd_2, hd_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * hd_1[k];

        t_11[k] = f_4 * hp_2[k]
                  + pa_z[k] * hd_2[k];

        t_12[k] = f_5 * gd0_0[k]
                  - f_6 * gd1_0[k]
                  + pa_y[k] * hd_3[k];

        t_13[k] = f_7 * hp_5[k]
                  + pb_x[k] * ip_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_z, gd0_4, gd1_8, hd_4, hd_8, is0_1, \
                         is1_1, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * gd0_4[k]
                  - f_9 * gd1_8[k]
                  + pa_x[k] * hd_8[k];

        t_15[k] = f_1 * is0_1[k]
                  - f_2 * is1_1[k]
                  + pb_z[k] * ip_6[k];

        t_16[k] = pa_z[k] * hd_4[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, pb_x, pb_y, gd0_0, gd1_0, hp_4, \
                         hp_9, hd_5, hd_6, ip_7, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_10 * hp_4[k]
                  + pb_y[k] * ip_7[k];

        t_18[k] = pa_y[k] * hd_6[k];

        t_19[k] = f_5 * gd0_0[k]
                  - f_6 * gd1_0[k]
                  + pa_z[k] * hd_5[k];

        t_20[k] = f_7 * hp_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_y, gd0_1, gd0_6, gd1_3, gd1_12, \
                         hd_7, hd_12, is0_2, is1_2, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_8[k];

        t_22[k] = f_8 * gd0_6[k]
                  - f_9 * gd1_12[k]
                  + pa_x[k] * hd_12[k];

        t_23[k] = f_11 * gd0_1[k]
                  - f_12 * gd1_3[k]
                  + pa_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_x, pb_z, gd0_7, gd1_14, hp_10, hd_14, \
                         is0_3, is1_3, ip_10, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_13 * hp_10[k]
                  + pb_x[k] * ip_10[k];

        t_25[k] = f_11 * gd0_7[k]
                  - f_12 * gd1_14[k]
                  + pa_x[k] * hd_14[k];

        t_26[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, hp_6, hp_7, hd_7, \
                         hd_8, hd_9, hd_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_z[k] * hd_7[k];

        t_28[k] = pa_z[k] * hd_8[k];

        t_29[k] = f_4 * hp_7[k]
                  + pb_y[k] * ip_12[k];

        t_30[k] = f_4 * hp_6[k]
                  + pa_z[k] * hd_9[k];

        t_31[k] = pa_y[k] * hd_10[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pa_z, pb_y, gd0_2, gd1_5, hp_8, hp_9, \
                         hd_10, hd_11, hd_12, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_4 * hp_8[k]
                  + pa_y[k] * hd_11[k];

        t_33[k] = f_10 * hp_9[k]
                  + pb_y[k] * ip_13[k];

        t_34[k] = pa_y[k] * hd_12[k];

        t_35[k] = f_11 * gd0_2[k]
                  - f_12 * gd1_5[k]
                  + pa_z[k] * hd_10[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_x, pb_y, gd0_8, gd1_16, hp_15, hd_19, \
                         is0_4, is1_4, ip_14, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_13 * hp_15[k]
                  + pb_x[k] * ip_15[k];

        t_37[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_14[k];

        t_38[k] = f_11 * gd0_8[k]
                  - f_12 * gd1_16[k]
                  + pa_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_y, pb_x, gd0_3, gd0_9, gd1_7, gd1_18, \
                         hp_16, hd_13, hd_21, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_8 * gd0_3[k]
                  - f_9 * gd1_7[k]
                  + pa_y[k] * hd_13[k];

        t_40[k] = f_4 * hp_16[k]
                  + pb_x[k] * ip_16[k];

        t_41[k] = f_5 * gd0_9[k]
                  - f_6 * gd1_18[k]
                  + pa_x[k] * hd_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_y, pb_z, hp_12, hd_13, hd_14, is0_5, \
                         is1_5, ip_17, ip_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_17[k];

        t_43[k] = pa_z[k] * hd_13[k];

        t_44[k] = pa_z[k] * hd_14[k];

        t_45[k] = f_13 * hp_12[k]
                  + pb_y[k] * ip_18[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_x, pa_y, pa_z, gd0_5, gd0_11, gd1_10, gd1_24, \
                         hp_11, hd_15, hd_16, hd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * hp_11[k]
                  + pa_z[k] * hd_15[k];

        t_47[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_10[k]
                  + pa_y[k] * hd_16[k];

        t_48[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_24[k]
                  + pa_x[k] * hd_22[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pa_y, pb_y, gd0_12, gd1_26, hp_13, \
                         hp_14, hd_17, hd_18, hd_23, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_4 * hp_13[k]
                  + pb_y[k] * ip_19[k];

        t_50[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_26[k]
                  + pa_x[k] * hd_23[k];

        t_51[k] = pa_y[k] * hd_17[k];

        t_52[k] = f_4 * hp_14[k]
                  + pa_y[k] * hd_18[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pa_z, pb_x, pb_y, gd0_5, gd1_10, hp_15, \
                         hp_17, hd_17, hd_19, ip_20, ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * hp_15[k]
                  + pb_y[k] * ip_20[k];

        t_54[k] = pa_y[k] * hd_19[k];

        t_55[k] = f_8 * gd0_5[k]
                  - f_9 * gd1_10[k]
                  + pa_z[k] * hd_17[k];

        t_56[k] = f_4 * hp_17[k]
                  + pb_x[k] * ip_22[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_x, pb_y, gd0_14, gd1_32, hp_18, hd_25, hd_26, \
                         is0_6, is1_6, ip_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_21[k];

        t_58[k] = f_5 * gd0_14[k]
                  - f_6 * gd1_32[k]
                  + pa_x[k] * hd_25[k];

        t_59[k] = f_4 * hp_18[k]
                  + pa_x[k] * hd_26[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_x, pa_z, pb_x, hp_19, hd_20, \
                         hd_27, hd_28, hd_30, hd_31, ip_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_10 * hp_19[k]
                  + pb_x[k] * ip_23[k];

        t_61[k] = pa_x[k] * hd_27[k];

        t_62[k] = pa_x[k] * hd_28[k];

        t_63[k] = pa_z[k] * hd_20[k];

        t_64[k] = pa_x[k] * hd_30[k];

        t_65[k] = pa_x[k] * hd_31[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, pa_x, hp_22, hp_24, hd_32, hd_33, \
                         hd_34, hd_35, hd_36, hd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * hp_22[k]
                  + pa_x[k] * hd_32[k];

        t_67[k] = pa_x[k] * hd_33[k];

        t_68[k] = pa_x[k] * hd_34[k];

        t_69[k] = pa_x[k] * hd_35[k];

        t_70[k] = f_4 * hp_24[k]
                  + pa_x[k] * hd_36[k];

        t_71[k] = pa_x[k] * hd_37[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, t_77, pa_x, pa_y, hp_27, hd_24, hd_38, \
                         hd_39, hd_40, hd_41, hd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_x[k] * hd_38[k];

        t_73[k] = pa_x[k] * hd_39[k];

        t_74[k] = pa_y[k] * hd_24[k];

        t_75[k] = pa_x[k] * hd_40[k];

        t_76[k] = pa_x[k] * hd_41[k];

        t_77[k] = f_4 * hp_27[k]
                  + pa_x[k] * hd_43[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_x, hp_29, hd_44, hd_45, is0_7, \
                         is1_7, ip_24, ip_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * hp_29[k]
                  + pb_x[k] * ip_24[k];

        t_79[k] = pa_x[k] * hd_44[k];

        t_80[k] = pa_x[k] * hd_45[k];

        t_81[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_x[k] * ip_25[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_z, pb_y, pb_z, hp_19, hd_26, hd_27, \
                         is0_7, is1_7, ip_26, ip_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * hp_19[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_26[k];

        t_83[k] = pb_z[k] * ip_26[k];

        t_84[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_z[k] * ip_27[k];

        t_85[k] = pa_z[k] * hd_26[k];

        t_86[k] = pa_z[k] * hd_27[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_z, pb_x, pb_y, hp_20, hp_21, hd_28, is0_8, \
                         is1_8, ip_28, ip_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * hp_21[k]
                  + pb_y[k] * ip_28[k];

        t_88[k] = f_4 * hp_20[k]
                  + pa_z[k] * hd_28[k];

        t_89[k] = f_1 * is0_8[k]
                  - f_2 * is1_8[k]
                  + pb_x[k] * ip_29[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pa_y, pa_z, pb_y, gd0_9, gd0_12, gd1_18, gd1_26, \
                         hp_23, hd_29, hd_35, ip_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_5 * gd0_9[k]
                  - f_6 * gd1_18[k]
                  + pa_z[k] * hd_29[k];

        t_91[k] = f_7 * hp_23[k]
                  + pb_y[k] * ip_30[k];

        t_92[k] = f_8 * gd0_12[k]
                  - f_9 * gd1_26[k]
                  + pa_y[k] * hd_35[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pa_z, pb_x, pb_y, gd0_10, gd1_20, hp_25, hd_33, \
                         is0_9, is1_9, ip_31, ip_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_31[k];

        t_94[k] = f_11 * gd0_10[k]
                  - f_12 * gd1_20[k]
                  + pa_z[k] * hd_33[k];

        t_95[k] = f_13 * hp_25[k]
                  + pb_y[k] * ip_32[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_y, pa_z, pb_x, gd0_11, gd0_13, gd1_24, gd1_29, \
                         hd_37, hd_39, is0_10, is1_10, ip_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_11 * gd0_13[k]
                  - f_12 * gd1_29[k]
                  + pa_y[k] * hd_39[k];

        t_97[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_33[k];

        t_98[k] = f_8 * gd0_11[k]
                  - f_9 * gd1_24[k]
                  + pa_z[k] * hd_37[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_y, pb_y, gd0_14, gd1_32, hp_26, hp_28, \
                         hd_42, hd_43, hd_44, ip_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_4 * hp_26[k]
                  + pb_y[k] * ip_34[k];

        t_100[k] = f_5 * gd0_14[k]
                   - f_6 * gd1_32[k]
                   + pa_y[k] * hd_42[k];

        t_101[k] = pa_y[k] * hd_43[k];

        t_102[k] = f_4 * hp_28[k]
                   + pa_y[k] * hd_44[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_y, pb_x, pb_y, hp_29, hd_45, \
                         is0_11, is1_11, ip_35, ip_36, ip_37, ip_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_10 * hp_29[k]
                   + pb_y[k] * ip_35[k];

        t_104[k] = pa_y[k] * hd_45[k];

        t_105[k] = f_1 * is0_11[k]
                   - f_2 * is1_11[k]
                   + pb_x[k] * ip_36[k];

        t_106[k] = f_1 * is0_11[k]
                   - f_2 * is1_11[k]
                   + pb_y[k] * ip_37[k];

        t_107[k] = pb_y[k] * ip_38[k];
    }

#pragma omp simd aligned(t_108, pb_z, hp_29, is0_11, is1_11, ip_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * hp_29[k]
                   + f_1 * is0_11[k]
                   - f_2 * is1_11[k]
                   + pb_z[k] * ip_38[k];
    }
}

auto
compute_prim_id_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is0_0, is1_0, ip_0, ip_1, \
                         ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, gd0_0, gd0_4, gd1_0, gd1_4, hd_3, \
                         hd_6, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_6[k];

        t_5[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, gd0_0, gd0_6, gd1_0, gd1_6, hd_4, \
                         hd_10, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_7[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];

        t_8[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, gd0_1, gd0_7, gd1_1, gd1_7, hd_5, \
                         hd_12, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_5[k];

        t_10[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_7[k]
                  + pa_x[k] * hd_12[k];

        t_11[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_2, gd1_8, hd_8, \
                         hd_16, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_8[k];

        t_13[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_14[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_3, gd1_9, \
                         hd_11, hd_17, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_11[k];

        t_16[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_17[k];

        t_17[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, gd0_5, gd0_14, gd1_5, gd1_14, \
                         hd_14, hd_18, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_14[k];

        t_19[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_20[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_18[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pb_y, pb_z, hp_8, is0_7, is0_8, is1_7, \
                         is1_8, ip_9, ip_10, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_x[k] * ip_9[k];

        t_22[k] = f_0 * hp_8[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_10[k];

        t_23[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_z[k] * ip_11[k];

        t_24[k] = f_1 * is0_8[k]
                  - f_2 * is1_8[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_9, gd1_12, \
                         hd_22, hd_25, is0_9, is1_9, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_22[k];

        t_26[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_25[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_10, gd1_13, \
                         hd_24, hd_28, is0_10, is1_10, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_24[k];

        t_29[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_28[k];

        t_30[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_11, gd1_14, \
                         hd_27, hd_29, is0_11, is1_11, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_27[k];

        t_32[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_29[k];

        t_33[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_15[k];
    }

#pragma omp simd aligned(t_34, t_35, pb_y, pb_z, hp_14, is0_11, is1_11, ip_16, \
                         ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_y[k] * ip_16[k];

        t_35[k] = f_0 * hp_14[k]
                  + f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_23 = buffer.data(gd1 + 23);
    const auto *gd1_26 = buffer.data(gd1 + 26);
    const auto *gd1_29 = buffer.data(gd1 + 29);
    const auto *gd1_30 = buffer.data(gd1 + 30);
    const auto *gd1_32 = buffer.data(gd1 + 32);
    const auto *gd1_35 = buffer.data(gd1 + 35);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_47 = buffer.data(hd + 47);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is0_0, is1_0, ip_0, ip_1, \
                         ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, gd0_0, gd0_4, gd1_0, gd1_10, hd_3, \
                         hd_10, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_10[k]
                 + pa_x[k] * hd_10[k];

        t_5[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, gd0_0, gd0_6, gd1_0, gd1_16, hd_6, \
                         hd_16, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_6[k];

        t_7[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];

        t_8[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_16[k]
                 + pa_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, gd0_1, gd0_7, gd1_3, gd1_18, hd_9, \
                         hd_18, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_3[k]
                 + pa_y[k] * hd_9[k];

        t_10[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_18[k]
                  + pa_x[k] * hd_18[k];

        t_11[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_6, gd1_21, \
                         hd_14, hd_27, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_6[k]
                  + pa_z[k] * hd_14[k];

        t_13[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_14[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_21[k]
                  + pa_x[k] * hd_27[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_9, gd1_23, \
                         hd_17, hd_29, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_9[k]
                  + pa_y[k] * hd_17[k];

        t_16[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_23[k]
                  + pa_x[k] * hd_29[k];

        t_17[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, gd0_5, gd0_14, gd1_14, gd1_35, \
                         hd_25, hd_33, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_14[k]
                  + pa_z[k] * hd_25[k];

        t_19[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_20[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_35[k]
                  + pa_x[k] * hd_33[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pb_y, pb_z, hp_8, is0_7, is0_8, is1_7, \
                         is1_8, ip_9, ip_10, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_x[k] * ip_9[k];

        t_22[k] = f_0 * hp_8[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_10[k];

        t_23[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_z[k] * ip_11[k];

        t_24[k] = f_1 * is0_8[k]
                  - f_2 * is1_8[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_23, gd1_30, \
                         hd_38, hd_42, is0_9, is1_9, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_23[k]
                  + pa_z[k] * hd_38[k];

        t_26[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_30[k]
                  + pa_y[k] * hd_42[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_26, gd1_32, \
                         hd_41, hd_45, is0_10, is1_10, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_26[k]
                  + pa_z[k] * hd_41[k];

        t_29[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_32[k]
                  + pa_y[k] * hd_45[k];

        t_30[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_29, gd1_35, \
                         hd_44, hd_47, is0_11, is1_11, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_29[k]
                  + pa_z[k] * hd_44[k];

        t_32[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_35[k]
                  + pa_y[k] * hd_47[k];

        t_33[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_15[k];
    }

#pragma omp simd aligned(t_34, t_35, pb_y, pb_z, hp_14, is0_11, is1_11, ip_16, \
                         ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_y[k] * ip_16[k];

        t_35[k] = f_0 * hp_14[k]
                  + f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_23 = buffer.data(gd0 + 23);
    const auto *gd0_26 = buffer.data(gd0 + 26);
    const auto *gd0_29 = buffer.data(gd0 + 29);
    const auto *gd0_30 = buffer.data(gd0 + 30);
    const auto *gd0_32 = buffer.data(gd0 + 32);
    const auto *gd0_35 = buffer.data(gd0 + 35);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_23 = buffer.data(gd1 + 23);
    const auto *gd1_24 = buffer.data(gd1 + 24);
    const auto *gd1_26 = buffer.data(gd1 + 26);
    const auto *gd1_29 = buffer.data(gd1 + 29);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
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
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, gd0_0, gd1_0, hp_1, hp_2, \
                         hd_0, hd_1, hd_2, hd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_5[k] = pa_y[k] * hd_2[k];

        t_6[k] = pa_z[k] * hd_0[k];

        t_7[k] = pa_z[k] * hd_1[k];

        t_8[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_9[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_z, gd0_10, gd1_8, hd_4, \
                         hd_6, hd_8, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * gd0_10[k]
                  - f_7 * gd1_8[k]
                  + pa_x[k] * hd_8[k];

        t_11[k] = f_1 * is0_1[k]
                  - f_2 * is1_1[k]
                  + pb_z[k] * ip_3[k];

        t_12[k] = pa_z[k] * hd_4[k];

        t_13[k] = pa_y[k] * hd_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_y, gd0_0, gd0_16, gd1_0, gd1_12, \
                         hd_5, hd_12, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * gd0_0[k]
                  - f_5 * gd1_0[k]
                  + pa_z[k] * hd_5[k];

        t_15[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_4[k];

        t_16[k] = f_6 * gd0_16[k]
                  - f_7 * gd1_12[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, gd0_3, gd0_18, gd1_3, gd1_14, \
                         hd_7, hd_14, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * gd0_3[k]
                  - f_9 * gd1_3[k]
                  + pa_y[k] * hd_7[k];

        t_18[k] = f_8 * gd0_18[k]
                  - f_9 * gd1_14[k]
                  + pa_x[k] * hd_14[k];

        t_19[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, hp_3, hp_4, hd_7, hd_8, \
                         hd_9, hd_11, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * hd_7[k];

        t_21[k] = pa_z[k] * hd_8[k];

        t_22[k] = f_3 * hp_3[k]
                  + pa_z[k] * hd_9[k];

        t_23[k] = f_3 * hp_4[k]
                  + pa_y[k] * hd_11[k];

        t_24[k] = pa_y[k] * hd_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_z, pb_y, gd0_6, gd0_21, gd1_5, gd1_16, \
                         hd_10, hd_19, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * gd0_6[k]
                  - f_9 * gd1_5[k]
                  + pa_z[k] * hd_10[k];

        t_26[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_27[k] = f_8 * gd0_21[k]
                  - f_9 * gd1_16[k]
                  + pa_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pa_y, pb_z, gd0_9, gd0_23, gd1_7, gd1_18, \
                         hd_13, hd_21, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_6 * gd0_9[k]
                  - f_7 * gd1_7[k]
                  + pa_y[k] * hd_13[k];

        t_29[k] = f_4 * gd0_23[k]
                  - f_5 * gd1_18[k]
                  + pa_x[k] * hd_21[k];

        t_30[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, gd0_14, gd1_10, hp_5, hd_13, \
                         hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_z[k] * hd_13[k];

        t_32[k] = pa_z[k] * hd_14[k];

        t_33[k] = f_3 * hp_5[k]
                  + pa_z[k] * hd_15[k];

        t_34[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_10[k]
                  + pa_y[k] * hd_16[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pa_y, gd0_29, gd0_30, gd1_23, gd1_24, \
                         hp_6, hd_18, hd_19, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_4 * gd0_29[k]
                  - f_5 * gd1_23[k]
                  + pa_x[k] * hd_22[k];

        t_36[k] = f_4 * gd0_30[k]
                  - f_5 * gd1_24[k]
                  + pa_x[k] * hd_23[k];

        t_37[k] = f_3 * hp_6[k]
                  + pa_y[k] * hd_18[k];

        t_38[k] = pa_y[k] * hd_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_z, pb_y, gd0_14, gd0_35, gd1_10, gd1_29, \
                         hd_17, hd_25, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_6 * gd0_14[k]
                  - f_7 * gd1_10[k]
                  + pa_z[k] * hd_17[k];

        t_40[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_41[k] = f_4 * gd0_35[k]
                  - f_5 * gd1_29[k]
                  + pa_x[k] * hd_25[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_x, pa_z, hp_7, hp_10, hp_11, hd_20, \
                         hd_26, hd_27, hd_31, hd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * hp_7[k]
                  + pa_x[k] * hd_26[k];

        t_43[k] = pa_x[k] * hd_27[k];

        t_44[k] = pa_z[k] * hd_20[k];

        t_45[k] = f_3 * hp_10[k]
                  + pa_x[k] * hd_31[k];

        t_46[k] = f_3 * hp_11[k]
                  + pa_x[k] * hd_34[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pb_x, pb_y, hp_8, hp_12, hd_39, hd_41, \
                         is0_7, is1_7, ip_9, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * hp_12[k]
                  + pa_x[k] * hd_39[k];

        t_48[k] = pa_x[k] * hd_41[k];

        t_49[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_x[k] * ip_9[k];

        t_50[k] = f_0 * hp_8[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_10[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_z, hp_9, hd_26, hd_27, hd_28, is0_7, \
                         is1_7, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_z[k] * ip_11[k];

        t_52[k] = pa_z[k] * hd_26[k];

        t_53[k] = pa_z[k] * hd_27[k];

        t_54[k] = f_3 * hp_9[k]
                  + pa_z[k] * hd_28[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pa_z, pb_x, gd0_23, gd0_30, gd1_18, gd1_24, \
                         hd_29, hd_33, is0_8, is1_8, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * is0_8[k]
                  - f_2 * is1_8[k]
                  + pb_x[k] * ip_12[k];

        t_56[k] = f_4 * gd0_23[k]
                  - f_5 * gd1_18[k]
                  + pa_z[k] * hd_29[k];

        t_57[k] = f_6 * gd0_30[k]
                  - f_7 * gd1_24[k]
                  + pa_y[k] * hd_33[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_y, pa_z, pb_x, gd0_26, gd0_32, gd1_20, gd1_26, \
                         hd_32, hd_36, is0_9, is1_9, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_13[k];

        t_59[k] = f_8 * gd0_26[k]
                  - f_9 * gd1_20[k]
                  + pa_z[k] * hd_32[k];

        t_60[k] = f_8 * gd0_32[k]
                  - f_9 * gd1_26[k]
                  + pa_y[k] * hd_36[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_y, pa_z, pb_x, gd0_29, gd0_35, gd1_23, gd1_29, \
                         hd_35, hd_38, is0_10, is1_10, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_14[k];

        t_62[k] = f_6 * gd0_29[k]
                  - f_7 * gd1_23[k]
                  + pa_z[k] * hd_35[k];

        t_63[k] = f_4 * gd0_35[k]
                  - f_5 * gd1_29[k]
                  + pa_y[k] * hd_38[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, hp_13, hd_40, hd_41, \
                         is0_11, is1_11, ip_15, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * hp_13[k]
                  + pa_y[k] * hd_40[k];

        t_65[k] = pa_y[k] * hd_41[k];

        t_66[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_15[k];

        t_67[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_y[k] * ip_16[k];
    }

#pragma omp simd aligned(t_68, pb_z, hp_14, is0_11, is1_11, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * hp_14[k]
                  + f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is0_0, is1_0, ip_0, ip_1, \
                         ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, gd0_0, gd0_4, gd1_0, gd1_6, hd_3, \
                         hd_6, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_5[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, gd0_0, gd0_6, gd1_0, gd1_10, hd_4, \
                         hd_10, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_7[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];

        t_8[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_10[k]
                 + pa_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, gd0_1, gd0_7, gd1_3, gd1_11, hd_5, \
                         hd_12, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_3[k]
                 + pa_y[k] * hd_5[k];

        t_10[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_11[k]
                  + pa_x[k] * hd_12[k];

        t_11[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_4, gd1_12, \
                         hd_8, hd_16, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_8[k];

        t_13[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_14[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_12[k]
                  + pa_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_5, gd1_14, \
                         hd_11, hd_17, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_16[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_17[k];

        t_17[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, gd0_5, gd0_14, gd1_8, gd1_23, \
                         hd_14, hd_18, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_8[k]
                  + pa_z[k] * hd_14[k];

        t_19[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_20[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_23[k]
                  + pa_x[k] * hd_18[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pb_y, pb_z, hp_8, is0_7, is0_8, is1_7, \
                         is1_8, ip_9, ip_10, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_x[k] * ip_9[k];

        t_22[k] = f_0 * hp_8[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_10[k];

        t_23[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_z[k] * ip_11[k];

        t_24[k] = f_1 * is0_8[k]
                  - f_2 * is1_8[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_14, gd1_19, \
                         hd_22, hd_25, is0_9, is1_9, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_14[k]
                  + pa_z[k] * hd_22[k];

        t_26[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_19[k]
                  + pa_y[k] * hd_25[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_16, gd1_20, \
                         hd_24, hd_28, is0_10, is1_10, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_16[k]
                  + pa_z[k] * hd_24[k];

        t_29[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_20[k]
                  + pa_y[k] * hd_28[k];

        t_30[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_18, gd1_23, \
                         hd_27, hd_29, is0_11, is1_11, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_18[k]
                  + pa_z[k] * hd_27[k];

        t_32[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_23[k]
                  + pa_y[k] * hd_29[k];

        t_33[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_15[k];
    }

#pragma omp simd aligned(t_34, t_35, pb_y, pb_z, hp_14, is0_11, is1_11, ip_16, \
                         ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_y[k] * ip_16[k];

        t_35[k] = f_0 * hp_14[k]
                  + f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_20 = buffer.data(gd0 + 20);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_21 = buffer.data(gd1 + 21);
    const auto *gd1_23 = buffer.data(gd1 + 23);
    const auto *gd1_26 = buffer.data(gd1 + 26);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is0_0, is1_0, ip_0, ip_1, \
                         ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, gd0_0, gd0_6, gd1_0, gd1_7, hd_3, \
                         hd_7, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_4[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_7[k]
                 + pa_x[k] * hd_7[k];

        t_5[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, gd0_0, gd0_10, gd1_0, gd1_11, hd_4, \
                         hd_11, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_7[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];

        t_8[k] = f_5 * gd0_10[k]
                 - f_6 * gd1_11[k]
                 + pa_x[k] * hd_11[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, gd0_3, gd0_11, gd1_3, gd1_12, \
                         hd_6, hd_13, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * gd0_3[k]
                 - f_8 * gd1_3[k]
                 + pa_y[k] * hd_6[k];

        t_10[k] = f_7 * gd0_11[k]
                  - f_8 * gd1_12[k]
                  + pa_x[k] * hd_13[k];

        t_11[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, gd0_4, gd0_12, gd1_4, gd1_13, \
                         hd_9, hd_17, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_4[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_9[k];

        t_13[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_14[k] = f_7 * gd0_12[k]
                  - f_8 * gd1_13[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, gd0_5, gd0_14, gd1_6, gd1_15, \
                         hd_12, hd_18, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_6[k]
                  + pa_y[k] * hd_12[k];

        t_16[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_15[k]
                  + pa_x[k] * hd_18[k];

        t_17[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, gd0_8, gd0_23, gd1_9, gd1_26, \
                         hd_15, hd_19, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_8[k]
                  - f_6 * gd1_9[k]
                  + pa_z[k] * hd_15[k];

        t_19[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_20[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_26[k]
                  + pa_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pb_y, pb_z, hp_8, is0_7, is0_8, is1_7, \
                         is1_8, ip_9, ip_10, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_x[k] * ip_9[k];

        t_22[k] = f_0 * hp_8[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_10[k];

        t_23[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_z[k] * ip_11[k];

        t_24[k] = f_1 * is0_8[k]
                  - f_2 * is1_8[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, gd0_14, gd0_19, gd1_15, gd1_21, \
                         hd_23, hd_27, is0_9, is1_9, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_15[k]
                  + pa_z[k] * hd_23[k];

        t_26[k] = f_5 * gd0_19[k]
                  - f_6 * gd1_21[k]
                  + pa_y[k] * hd_27[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, pb_x, gd0_16, gd0_20, gd1_17, gd1_23, \
                         hd_26, hd_30, is0_10, is1_10, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_7 * gd0_16[k]
                  - f_8 * gd1_17[k]
                  + pa_z[k] * hd_26[k];

        t_29[k] = f_7 * gd0_20[k]
                  - f_8 * gd1_23[k]
                  + pa_y[k] * hd_30[k];

        t_30[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, gd0_18, gd0_23, gd1_20, gd1_26, \
                         hd_29, hd_32, is0_11, is1_11, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_5 * gd0_18[k]
                  - f_6 * gd1_20[k]
                  + pa_z[k] * hd_29[k];

        t_32[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_26[k]
                  + pa_y[k] * hd_32[k];

        t_33[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_15[k];
    }

#pragma omp simd aligned(t_34, t_35, pb_y, pb_z, hp_14, is0_11, is1_11, ip_16, \
                         ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_y[k] * ip_16[k];

        t_35[k] = f_0 * hp_14[k]
                  + f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_20 = buffer.data(gd0 + 20);
    const auto *gd0_21 = buffer.data(gd0 + 21);
    const auto *gd0_23 = buffer.data(gd0 + 23);
    const auto *gd0_26 = buffer.data(gd0 + 26);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is0_0, is1_0, ip_0, ip_1, \
                         ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, gd0_0, gd0_7, gd1_0, gd1_6, hd_3, \
                         hd_6, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_4[k] = f_5 * gd0_7[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_5[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, gd0_0, gd0_11, gd1_0, gd1_10, hd_4, \
                         hd_10, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_7[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];

        t_8[k] = f_5 * gd0_11[k]
                 - f_6 * gd1_10[k]
                 + pa_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, gd0_3, gd0_12, gd1_3, gd1_11, \
                         hd_5, hd_12, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * gd0_3[k]
                 - f_8 * gd1_3[k]
                 + pa_y[k] * hd_5[k];

        t_10[k] = f_7 * gd0_12[k]
                  - f_8 * gd1_11[k]
                  + pa_x[k] * hd_12[k];

        t_11[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, gd0_4, gd0_13, gd1_4, gd1_12, \
                         hd_8, hd_16, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_4[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_8[k];

        t_13[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_14[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_12[k]
                  + pa_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, gd0_6, gd0_15, gd1_5, gd1_14, \
                         hd_11, hd_17, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * gd0_6[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_16[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_17[k];

        t_17[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, gd0_9, gd0_26, gd1_8, gd1_23, \
                         hd_14, hd_18, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_9[k]
                  - f_6 * gd1_8[k]
                  + pa_z[k] * hd_14[k];

        t_19[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_20[k] = f_3 * gd0_26[k]
                  - f_4 * gd1_23[k]
                  + pa_x[k] * hd_18[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pb_y, pb_z, hp_8, is0_7, is0_8, is1_7, \
                         is1_8, ip_9, ip_10, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_x[k] * ip_9[k];

        t_22[k] = f_0 * hp_8[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_10[k];

        t_23[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_z[k] * ip_11[k];

        t_24[k] = f_1 * is0_8[k]
                  - f_2 * is1_8[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, gd0_15, gd0_21, gd1_14, gd1_19, \
                         hd_22, hd_25, is0_9, is1_9, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_14[k]
                  + pa_z[k] * hd_22[k];

        t_26[k] = f_5 * gd0_21[k]
                  - f_6 * gd1_19[k]
                  + pa_y[k] * hd_25[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, pb_x, gd0_17, gd0_23, gd1_16, gd1_20, \
                         hd_24, hd_28, is0_10, is1_10, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_7 * gd0_17[k]
                  - f_8 * gd1_16[k]
                  + pa_z[k] * hd_24[k];

        t_29[k] = f_7 * gd0_23[k]
                  - f_8 * gd1_20[k]
                  + pa_y[k] * hd_28[k];

        t_30[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, gd0_20, gd0_26, gd1_18, gd1_23, \
                         hd_27, hd_29, is0_11, is1_11, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_5 * gd0_20[k]
                  - f_6 * gd1_18[k]
                  + pa_z[k] * hd_27[k];

        t_32[k] = f_3 * gd0_26[k]
                  - f_4 * gd1_23[k]
                  + pa_y[k] * hd_29[k];

        t_33[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_15[k];
    }

#pragma omp simd aligned(t_34, t_35, pb_y, pb_z, hp_14, is0_11, is1_11, ip_16, \
                         ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_y[k] * ip_16[k];

        t_35[k] = f_0 * hp_14[k]
                  + f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_20 = buffer.data(gd0 + 20);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is0_0, is1_0, ip_0, ip_1, \
                         ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, gd0_0, gd0_6, gd1_0, gd1_6, hd_3, \
                         hd_6, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_4[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_5[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, gd0_0, gd0_10, gd1_0, gd1_10, hd_4, \
                         hd_10, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_7[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];

        t_8[k] = f_5 * gd0_10[k]
                 - f_6 * gd1_10[k]
                 + pa_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, gd0_3, gd0_11, gd1_3, gd1_11, \
                         hd_5, hd_12, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * gd0_3[k]
                 - f_8 * gd1_3[k]
                 + pa_y[k] * hd_5[k];

        t_10[k] = f_7 * gd0_11[k]
                  - f_8 * gd1_11[k]
                  + pa_x[k] * hd_12[k];

        t_11[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, gd0_4, gd0_12, gd1_4, gd1_12, \
                         hd_8, hd_16, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_4[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_8[k];

        t_13[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_14[k] = f_7 * gd0_12[k]
                  - f_8 * gd1_12[k]
                  + pa_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, gd0_5, gd0_14, gd1_5, gd1_14, \
                         hd_11, hd_17, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_16[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_17[k];

        t_17[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, gd0_8, gd0_23, gd1_8, gd1_23, \
                         hd_14, hd_18, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_8[k]
                  - f_6 * gd1_8[k]
                  + pa_z[k] * hd_14[k];

        t_19[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_20[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pa_x[k] * hd_18[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pb_y, pb_z, hp_8, is0_7, is0_8, is1_7, \
                         is1_8, ip_9, ip_10, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_x[k] * ip_9[k];

        t_22[k] = f_0 * hp_8[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_10[k];

        t_23[k] = f_1 * is0_7[k]
                  - f_2 * is1_7[k]
                  + pb_z[k] * ip_11[k];

        t_24[k] = f_1 * is0_8[k]
                  - f_2 * is1_8[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_x, gd0_14, gd0_19, gd1_14, gd1_19, \
                         hd_22, hd_25, is0_9, is1_9, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_z[k] * hd_22[k];

        t_26[k] = f_5 * gd0_19[k]
                  - f_6 * gd1_19[k]
                  + pa_y[k] * hd_25[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, pb_x, gd0_16, gd0_20, gd1_16, gd1_20, \
                         hd_24, hd_28, is0_10, is1_10, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_7 * gd0_16[k]
                  - f_8 * gd1_16[k]
                  + pa_z[k] * hd_24[k];

        t_29[k] = f_7 * gd0_20[k]
                  - f_8 * gd1_20[k]
                  + pa_y[k] * hd_28[k];

        t_30[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, gd0_18, gd0_23, gd1_18, gd1_23, \
                         hd_27, hd_29, is0_11, is1_11, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_5 * gd0_18[k]
                  - f_6 * gd1_18[k]
                  + pa_z[k] * hd_27[k];

        t_32[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pa_y[k] * hd_29[k];

        t_33[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_15[k];
    }

#pragma omp simd aligned(t_34, t_35, pb_y, pb_z, hp_14, is0_11, is1_11, ip_16, \
                         ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_y[k] * ip_16[k];

        t_35[k] = f_0 * hp_14[k]
                  + f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gd0, const size_t gd1,
                                     const size_t hp, const size_t hd, const size_t is0,
                                     const size_t is1, const size_t ip, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 2.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.0 / p;
    const auto f_8 = 1.5 / alpha;
    const auto f_9 = 1.5 * beta / (alpha * p);
    const auto f_10 = 1.0 / alpha;
    const auto f_11 = beta / (alpha * p);
    const auto f_12 = 1.5 / p;
    const auto f_13 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
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
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_13 = buffer.data(is1 + 13);
    const auto *is1_16 = buffer.data(is1 + 16);
    const auto *is1_18 = buffer.data(is1 + 18);
    const auto *is1_19 = buffer.data(is1 + 19);
    const auto *is1_20 = buffer.data(is1 + 20);
    const auto *is1_22 = buffer.data(is1 + 22);

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
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, pb_x, hp_1, hp_3, hp_4, hd_0, hd_1, \
                         ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * hp_3[k]
                 + pb_x[k] * ip_3[k];

        t_5[k] = f_4 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_6[k] = pa_z[k] * hd_0[k];

        t_7[k] = f_3 * hp_4[k]
                 + pb_x[k] * ip_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_y, pa_z, pb_x, gd0_0, gd1_0, hp_2, hp_5, hd_2, \
                         hd_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_9[k] = f_5 * gd0_0[k]
                 - f_6 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_10[k] = f_7 * hp_5[k]
                  + pb_x[k] * ip_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_z, pb_z, gd0_0, gd0_4, gd1_0, gd1_4, hd_4, \
                         hd_6, is0_1, is1_3, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_8 * gd0_4[k]
                  - f_9 * gd1_4[k]
                  + pa_x[k] * hd_6[k];

        t_12[k] = f_1 * is0_1[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_6[k];

        t_13[k] = f_5 * gd0_0[k]
                  - f_6 * gd1_0[k]
                  + pa_z[k] * hd_4[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pb_x, pb_y, gd0_6, gd1_6, hp_6, hd_8, is0_2, \
                         is1_4, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * hp_6[k]
                  + pb_x[k] * ip_8[k];

        t_15[k] = f_1 * is0_2[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_7[k];

        t_16[k] = f_8 * gd0_6[k]
                  - f_9 * gd1_6[k]
                  + pa_x[k] * hd_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_x, gd0_1, gd0_7, gd1_1, gd1_7, hp_7, \
                         hd_5, hd_10, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_10 * gd0_1[k]
                  - f_11 * gd1_1[k]
                  + pa_y[k] * hd_5[k];

        t_18[k] = f_12 * hp_7[k]
                  + pb_x[k] * ip_9[k];

        t_19[k] = f_10 * gd0_7[k]
                  - f_11 * gd1_7[k]
                  + pa_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pa_z, pb_x, pb_z, gd0_2, gd1_2, hp_8, \
                         hd_7, is0_3, is1_5, ip_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * is0_3[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_10[k];

        t_21[k] = pa_y[k] * hd_7[k];

        t_22[k] = f_10 * gd0_2[k]
                  - f_11 * gd1_2[k]
                  + pa_z[k] * hd_7[k];

        t_23[k] = f_12 * hp_8[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_y, pb_y, gd0_3, gd0_8, gd1_3, gd1_8, hd_9, \
                         hd_13, is0_4, is1_8, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_1 * is0_4[k]
                  - f_2 * is1_8[k]
                  + pb_y[k] * ip_11[k];

        t_25[k] = f_10 * gd0_8[k]
                  - f_11 * gd1_8[k]
                  + pa_x[k] * hd_13[k];

        t_26[k] = f_8 * gd0_3[k]
                  - f_9 * gd1_3[k]
                  + pa_y[k] * hd_9[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pb_x, pb_z, gd0_9, gd1_9, hp_9, hd_14, is0_5, \
                         is1_9, ip_13, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * hp_9[k]
                  + pb_x[k] * ip_13[k];

        t_28[k] = f_5 * gd0_9[k]
                  - f_6 * gd1_9[k]
                  + pa_x[k] * hd_14[k];

        t_29[k] = f_1 * is0_5[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, gd0_5, gd0_11, gd0_12, gd1_5, \
                         gd1_11, gd1_12, hd_11, hd_12, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_31[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_x[k] * hd_15[k];

        t_32[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_x[k] * hd_16[k];

        t_33[k] = pa_y[k] * hd_12[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_x, pb_y, gd0_5, gd1_5, hp_10, hd_12, \
                         is0_6, is1_13, ip_15, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_8 * gd0_5[k]
                  - f_9 * gd1_5[k]
                  + pa_z[k] * hd_12[k];

        t_35[k] = f_4 * hp_10[k]
                  + pb_x[k] * ip_16[k];

        t_36[k] = f_1 * is0_6[k]
                  - f_2 * is1_13[k]
                  + pb_y[k] * ip_15[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_x, pb_x, gd0_14, gd1_14, hp_11, \
                         hp_12, hd_17, hd_18, hd_19, hd_22, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * gd0_14[k]
                  - f_6 * gd1_14[k]
                  + pa_x[k] * hd_17[k];

        t_38[k] = f_4 * hp_11[k]
                  + pa_x[k] * hd_18[k];

        t_39[k] = f_13 * hp_12[k]
                  + pb_x[k] * ip_17[k];

        t_40[k] = pa_x[k] * hd_19[k];

        t_41[k] = pa_x[k] * hd_22[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_x, pb_x, hp_21, hp_23, hd_23, \
                         hd_24, hd_25, hd_27, hd_29, ip_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_x[k] * hd_23[k];

        t_43[k] = pa_x[k] * hd_24[k];

        t_44[k] = pa_x[k] * hd_25[k];

        t_45[k] = f_4 * hp_21[k]
                  + pa_x[k] * hd_27[k];

        t_46[k] = f_13 * hp_23[k]
                  + pb_x[k] * ip_18[k];

        t_47[k] = pa_x[k] * hd_29[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, pb_z, hp_12, hd_19, is0_9, \
                         is1_16, ip_19, ip_20, ip_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * is0_9[k]
                  - f_2 * is1_16[k]
                  + pb_x[k] * ip_19[k];

        t_49[k] = f_0 * hp_12[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_16[k]
                  + pb_y[k] * ip_20[k];

        t_50[k] = f_1 * is0_9[k]
                  - f_2 * is1_16[k]
                  + pb_z[k] * ip_21[k];

        t_51[k] = pa_z[k] * hd_19[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_z, pb_x, pb_y, hp_13, hp_14, hd_20, is0_10, \
                         is1_18, ip_22, ip_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * hp_14[k]
                  + pb_y[k] * ip_22[k];

        t_53[k] = f_4 * hp_13[k]
                  + pa_z[k] * hd_20[k];

        t_54[k] = f_1 * is0_10[k]
                  - f_2 * is1_18[k]
                  + pb_x[k] * ip_23[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pa_z, pb_y, gd0_9, gd0_12, gd1_9, gd1_12, \
                         hp_16, hd_21, hd_23, ip_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_5 * gd0_9[k]
                  - f_6 * gd1_9[k]
                  + pa_z[k] * hd_21[k];

        t_56[k] = f_7 * hp_16[k]
                  + pb_y[k] * ip_25[k];

        t_57[k] = f_8 * gd0_12[k]
                  - f_9 * gd1_12[k]
                  + pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_z, pb_x, pb_y, gd0_10, gd1_10, hp_18, hd_22, \
                         is0_11, is1_19, ip_26, ip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * is0_11[k]
                  - f_2 * is1_19[k]
                  + pb_x[k] * ip_26[k];

        t_59[k] = f_10 * gd0_10[k]
                  - f_11 * gd1_10[k]
                  + pa_z[k] * hd_22[k];

        t_60[k] = f_12 * hp_18[k]
                  + pb_y[k] * ip_28[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_y, pa_z, pb_x, gd0_11, gd0_13, gd1_11, gd1_13, \
                         hd_24, hd_25, is0_12, is1_20, ip_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_10 * gd0_13[k]
                  - f_11 * gd1_13[k]
                  + pa_y[k] * hd_25[k];

        t_62[k] = f_1 * is0_12[k]
                  - f_2 * is1_20[k]
                  + pb_x[k] * ip_29[k];

        t_63[k] = f_8 * gd0_11[k]
                  - f_9 * gd1_11[k]
                  + pa_z[k] * hd_24[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_y, gd0_14, gd1_14, hp_20, hp_22, \
                         hp_23, hd_26, hd_28, ip_31, ip_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_4 * hp_20[k]
                  + pb_y[k] * ip_31[k];

        t_65[k] = f_5 * gd0_14[k]
                  - f_6 * gd1_14[k]
                  + pa_y[k] * hd_26[k];

        t_66[k] = f_4 * hp_22[k]
                  + pa_y[k] * hd_28[k];

        t_67[k] = f_13 * hp_23[k]
                  + pb_y[k] * ip_33[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pb_x, pb_y, pb_z, hp_23, hd_29, is0_14, \
                         is1_22, ip_34, ip_35, ip_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pa_y[k] * hd_29[k];

        t_69[k] = f_1 * is0_14[k]
                  - f_2 * is1_22[k]
                  + pb_x[k] * ip_34[k];

        t_70[k] = f_1 * is0_14[k]
                  - f_2 * is1_22[k]
                  + pb_y[k] * ip_35[k];

        t_71[k] = f_0 * hp_23[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_22[k]
                  + pb_z[k] * ip_36[k];
    }
}

auto
compute_prim_id_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 0.5 / p;
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.5 / p;
    const auto f_13 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
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
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, gd0_0, gd1_0, hp_1, hp_2, \
                         hd_0, hd_1, hd_2, hd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_5[k] = pa_y[k] * hd_2[k];

        t_6[k] = pa_z[k] * hd_0[k];

        t_7[k] = pa_z[k] * hd_1[k];

        t_8[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_9[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_z, gd0_4, gd1_4, hd_4, hd_8, is0_1, \
                         is1_1, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * gd0_4[k]
                  - f_7 * gd1_4[k]
                  + pa_x[k] * hd_8[k];

        t_11[k] = f_1 * is0_1[k]
                  - f_2 * is1_1[k]
                  + pb_z[k] * ip_6[k];

        t_12[k] = pa_z[k] * hd_4[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pa_z, pb_y, gd0_0, gd1_0, hp_3, hd_5, \
                         hd_6, is0_2, is1_2, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * hp_3[k]
                  + pb_y[k] * ip_7[k];

        t_14[k] = pa_y[k] * hd_6[k];

        t_15[k] = f_4 * gd0_0[k]
                  - f_5 * gd1_0[k]
                  + pa_z[k] * hd_5[k];

        t_16[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, gd0_1, gd0_6, gd0_7, gd1_1, gd1_6, \
                         gd1_7, hd_7, hd_12, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_6 * gd0_6[k]
                  - f_7 * gd1_6[k]
                  + pa_x[k] * hd_12[k];

        t_18[k] = f_9 * gd0_1[k]
                  - f_10 * gd1_1[k]
                  + pa_y[k] * hd_7[k];

        t_19[k] = f_9 * gd0_7[k]
                  - f_10 * gd1_7[k]
                  + pa_x[k] * hd_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_z, pb_y, pb_z, hp_5, hd_7, hd_8, is0_3, \
                         is1_3, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_11[k];

        t_21[k] = pa_z[k] * hd_7[k];

        t_22[k] = pa_z[k] * hd_8[k];

        t_23[k] = f_3 * hp_5[k]
                  + pb_y[k] * ip_12[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pa_z, pb_y, hp_4, hp_6, hp_7, \
                         hd_9, hd_10, hd_11, hd_12, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_3 * hp_4[k]
                  + pa_z[k] * hd_9[k];

        t_25[k] = pa_y[k] * hd_10[k];

        t_26[k] = f_3 * hp_6[k]
                  + pa_y[k] * hd_11[k];

        t_27[k] = f_8 * hp_7[k]
                  + pb_y[k] * ip_13[k];

        t_28[k] = pa_y[k] * hd_12[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_2, gd1_8, \
                         hd_10, hd_19, is0_4, is1_4, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_9 * gd0_2[k]
                  - f_10 * gd1_2[k]
                  + pa_z[k] * hd_10[k];

        t_30[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_14[k];

        t_31[k] = f_9 * gd0_8[k]
                  - f_10 * gd1_8[k]
                  + pa_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_3, gd1_9, \
                         hd_13, hd_21, is0_5, is1_5, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_6 * gd0_3[k]
                  - f_7 * gd1_3[k]
                  + pa_y[k] * hd_13[k];

        t_33[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_9[k]
                  + pa_x[k] * hd_21[k];

        t_34[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_z, pb_y, hp_8, hp_9, hd_13, hd_14, hd_15, \
                         ip_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_z[k] * hd_13[k];

        t_36[k] = pa_z[k] * hd_14[k];

        t_37[k] = f_11 * hp_9[k]
                  + pb_y[k] * ip_18[k];

        t_38[k] = f_3 * hp_8[k]
                  + pa_z[k] * hd_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_y, pb_y, gd0_5, gd0_11, gd1_5, gd1_11, \
                         hp_10, hd_16, hd_22, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_4 * gd0_5[k]
                  - f_5 * gd1_5[k]
                  + pa_y[k] * hd_16[k];

        t_40[k] = f_4 * gd0_11[k]
                  - f_5 * gd1_11[k]
                  + pa_x[k] * hd_22[k];

        t_41[k] = f_3 * hp_10[k]
                  + pb_y[k] * ip_19[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, gd0_12, gd1_12, hp_11, \
                         hp_12, hd_17, hd_18, hd_23, ip_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_4 * gd0_12[k]
                  - f_5 * gd1_12[k]
                  + pa_x[k] * hd_23[k];

        t_43[k] = pa_y[k] * hd_17[k];

        t_44[k] = f_3 * hp_11[k]
                  + pa_y[k] * hd_18[k];

        t_45[k] = f_8 * hp_12[k]
                  + pb_y[k] * ip_20[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_y, pa_z, pb_y, gd0_5, gd1_5, hd_17, hd_19, \
                         is0_6, is1_6, ip_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_y[k] * hd_19[k];

        t_47[k] = f_6 * gd0_5[k]
                  - f_7 * gd1_5[k]
                  + pa_z[k] * hd_17[k];

        t_48[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_21[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_x, pa_z, gd0_14, gd1_14, hp_13, \
                         hd_20, hd_25, hd_26, hd_27, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_14[k]
                  + pa_x[k] * hd_25[k];

        t_50[k] = f_3 * hp_13[k]
                  + pa_x[k] * hd_26[k];

        t_51[k] = pa_x[k] * hd_27[k];

        t_52[k] = pa_x[k] * hd_28[k];

        t_53[k] = pa_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pa_x, hp_17, hd_30, hd_31, hd_32, \
                         hd_33, hd_34, hd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_x[k] * hd_30[k];

        t_55[k] = pa_x[k] * hd_31[k];

        t_56[k] = f_3 * hp_17[k]
                  + pa_x[k] * hd_32[k];

        t_57[k] = pa_x[k] * hd_33[k];

        t_58[k] = pa_x[k] * hd_34[k];

        t_59[k] = pa_x[k] * hd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_x, pa_y, hp_19, hd_24, hd_36, \
                         hd_37, hd_38, hd_39, hd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * hp_19[k]
                  + pa_x[k] * hd_36[k];

        t_61[k] = pa_x[k] * hd_37[k];

        t_62[k] = pa_x[k] * hd_38[k];

        t_63[k] = pa_x[k] * hd_39[k];

        t_64[k] = pa_y[k] * hd_24[k];

        t_65[k] = pa_x[k] * hd_40[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_x, pb_x, hp_22, hd_41, hd_43, hd_44, \
                         hd_45, is0_7, is1_9, ip_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_x[k] * hd_41[k];

        t_67[k] = f_3 * hp_22[k]
                  + pa_x[k] * hd_43[k];

        t_68[k] = pa_x[k] * hd_44[k];

        t_69[k] = pa_x[k] * hd_45[k];

        t_70[k] = f_1 * is0_7[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_z, pb_y, pb_z, hp_14, hd_26, hd_27, is0_7, \
                         is1_9, ip_26, ip_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_0 * hp_14[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_26[k];

        t_72[k] = f_1 * is0_7[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_27[k];

        t_73[k] = pa_z[k] * hd_26[k];

        t_74[k] = pa_z[k] * hd_27[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_z, pb_x, pb_y, hp_15, hp_16, hd_28, is0_8, \
                         is1_10, ip_28, ip_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_12 * hp_16[k]
                  + pb_y[k] * ip_28[k];

        t_76[k] = f_3 * hp_15[k]
                  + pa_z[k] * hd_28[k];

        t_77[k] = f_1 * is0_8[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_29[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_y, pa_z, pb_y, gd0_9, gd0_12, gd1_9, gd1_12, \
                         hp_18, hd_29, hd_35, ip_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_9[k]
                  + pa_z[k] * hd_29[k];

        t_79[k] = f_13 * hp_18[k]
                  + pb_y[k] * ip_30[k];

        t_80[k] = f_6 * gd0_12[k]
                  - f_7 * gd1_12[k]
                  + pa_y[k] * hd_35[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pb_x, pb_y, gd0_10, gd1_10, hp_20, hd_33, \
                         is0_9, is1_11, ip_31, ip_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_1 * is0_9[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_31[k];

        t_82[k] = f_9 * gd0_10[k]
                  - f_10 * gd1_10[k]
                  + pa_z[k] * hd_33[k];

        t_83[k] = f_11 * hp_20[k]
                  + pb_y[k] * ip_32[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pa_z, pb_x, gd0_11, gd0_13, gd1_11, gd1_13, \
                         hd_37, hd_39, is0_10, is1_12, ip_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_9 * gd0_13[k]
                  - f_10 * gd1_13[k]
                  + pa_y[k] * hd_39[k];

        t_85[k] = f_1 * is0_10[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_33[k];

        t_86[k] = f_6 * gd0_11[k]
                  - f_7 * gd1_11[k]
                  + pa_z[k] * hd_37[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pb_y, gd0_14, gd1_14, hp_21, hp_23, \
                         hd_42, hd_43, hd_44, ip_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * hp_21[k]
                  + pb_y[k] * ip_34[k];

        t_88[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_14[k]
                  + pa_y[k] * hd_42[k];

        t_89[k] = pa_y[k] * hd_43[k];

        t_90[k] = f_3 * hp_23[k]
                  + pa_y[k] * hd_44[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pb_x, pb_y, hp_24, hd_45, is0_11, \
                         is1_14, ip_35, ip_36, ip_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_8 * hp_24[k]
                  + pb_y[k] * ip_35[k];

        t_92[k] = pa_y[k] * hd_45[k];

        t_93[k] = f_1 * is0_11[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_36[k];

        t_94[k] = f_1 * is0_11[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_37[k];
    }

#pragma omp simd aligned(t_95, pb_z, hp_24, is0_11, is1_14, ip_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_0 * hp_24[k]
                  + f_1 * is0_11[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_38[k];
    }
}

auto
compute_prim_id_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.5 / p;
    const auto f_11 = 2.0 / p;
    const auto f_12 = 1.5 / p;
    const auto f_13 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_20 = buffer.data(gd1 + 20);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_14 = buffer.data(is1 + 14);
    const auto *is1_16 = buffer.data(is1 + 16);
    const auto *is1_17 = buffer.data(is1 + 17);
    const auto *is1_18 = buffer.data(is1 + 18);
    const auto *is1_20 = buffer.data(is1 + 20);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, gd0_0, gd1_0, hp_1, hp_2, hd_0, hd_1, \
                         hd_2, hd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_5[k] = pa_z[k] * hd_0[k];

        t_6[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_7[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_z, pb_z, gd0_0, gd0_4, gd1_0, gd1_6, hd_4, \
                         hd_6, is0_1, is1_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * gd0_4[k]
                 - f_7 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_9[k] = f_1 * is0_1[k]
                 - f_2 * is1_3[k]
                 + pb_z[k] * ip_4[k];

        t_10[k] = f_4 * gd0_0[k]
                  - f_5 * gd1_0[k]
                  + pa_z[k] * hd_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, pb_y, gd0_1, gd0_6, gd1_3, gd1_8, hd_5, \
                         hd_8, is0_2, is1_4, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * is0_2[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_5[k];

        t_12[k] = f_6 * gd0_6[k]
                  - f_7 * gd1_8[k]
                  + pa_x[k] * hd_8[k];

        t_13[k] = f_8 * gd0_1[k]
                  - f_9 * gd1_3[k]
                  + pa_y[k] * hd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_z, gd0_7, gd1_9, hd_7, hd_10, is0_3, \
                         is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * gd0_7[k]
                  - f_9 * gd1_9[k]
                  + pa_x[k] * hd_10[k];

        t_15[k] = f_1 * is0_3[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_16[k] = pa_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_4, gd1_10, \
                         hd_7, hd_13, is0_4, is1_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * gd0_2[k]
                  - f_9 * gd1_4[k]
                  + pa_z[k] * hd_7[k];

        t_18[k] = f_1 * is0_4[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_8[k];

        t_19[k] = f_8 * gd0_8[k]
                  - f_9 * gd1_10[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_5, gd1_12, \
                         hd_9, hd_14, is0_5, is1_8, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * gd0_3[k]
                  - f_7 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_21[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_12[k]
                  + pa_x[k] * hd_14[k];

        t_22[k] = f_1 * is0_5[k]
                  - f_2 * is1_8[k]
                  + pb_z[k] * ip_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, gd0_5, gd0_11, gd0_12, gd1_7, \
                         gd1_15, gd1_16, hd_11, hd_12, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * gd0_5[k]
                  - f_5 * gd1_7[k]
                  + pa_y[k] * hd_11[k];

        t_24[k] = f_4 * gd0_11[k]
                  - f_5 * gd1_15[k]
                  + pa_x[k] * hd_15[k];

        t_25[k] = f_4 * gd0_12[k]
                  - f_5 * gd1_16[k]
                  + pa_x[k] * hd_16[k];

        t_26[k] = pa_y[k] * hd_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_z, pb_y, gd0_5, gd0_14, gd1_7, gd1_20, \
                         hd_12, hd_17, is0_6, is1_11, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * gd0_5[k]
                  - f_7 * gd1_7[k]
                  + pa_z[k] * hd_12[k];

        t_28[k] = f_1 * is0_6[k]
                  - f_2 * is1_11[k]
                  + pb_y[k] * ip_11[k];

        t_29[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_20[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_x, hp_10, hd_18, hd_19, hd_22, \
                         hd_23, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hp_10[k]
                  + pa_x[k] * hd_18[k];

        t_31[k] = pa_x[k] * hd_19[k];

        t_32[k] = pa_x[k] * hd_22[k];

        t_33[k] = pa_x[k] * hd_23[k];

        t_34[k] = pa_x[k] * hd_24[k];

        t_35[k] = pa_x[k] * hd_25[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_x, pb_y, hp_11, hp_19, hd_27, hd_29, \
                         is0_9, is1_14, ip_13, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * hp_19[k]
                  + pa_x[k] * hd_27[k];

        t_37[k] = pa_x[k] * hd_29[k];

        t_38[k] = f_1 * is0_9[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_13[k];

        t_39[k] = f_0 * hp_11[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_14[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_y, pb_z, hp_12, hp_13, hd_19, hd_20, \
                         is0_9, is1_14, ip_15, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * is0_9[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_15[k];

        t_41[k] = pa_z[k] * hd_19[k];

        t_42[k] = f_10 * hp_13[k]
                  + pb_y[k] * ip_16[k];

        t_43[k] = f_3 * hp_12[k]
                  + pa_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_z, pb_x, pb_y, gd0_9, gd1_12, hp_15, hd_21, \
                         is0_10, is1_16, ip_17, ip_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * is0_10[k]
                  - f_2 * is1_16[k]
                  + pb_x[k] * ip_17[k];

        t_45[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_12[k]
                  + pa_z[k] * hd_21[k];

        t_46[k] = f_11 * hp_15[k]
                  + pb_y[k] * ip_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pa_z, pb_x, gd0_10, gd0_12, gd1_14, gd1_16, \
                         hd_22, hd_23, is0_11, is1_17, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * gd0_12[k]
                  - f_7 * gd1_16[k]
                  + pa_y[k] * hd_23[k];

        t_48[k] = f_1 * is0_11[k]
                  - f_2 * is1_17[k]
                  + pb_x[k] * ip_19[k];

        t_49[k] = f_8 * gd0_10[k]
                  - f_9 * gd1_14[k]
                  + pa_z[k] * hd_22[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_y, pb_x, pb_y, gd0_13, gd1_17, hp_17, hd_25, \
                         is0_12, is1_18, ip_20, ip_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_12 * hp_17[k]
                  + pb_y[k] * ip_20[k];

        t_51[k] = f_8 * gd0_13[k]
                  - f_9 * gd1_17[k]
                  + pa_y[k] * hd_25[k];

        t_52[k] = f_1 * is0_12[k]
                  - f_2 * is1_18[k]
                  + pb_x[k] * ip_21[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pa_z, pb_y, gd0_11, gd0_14, gd1_15, gd1_20, \
                         hp_18, hd_24, hd_26, ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_6 * gd0_11[k]
                  - f_7 * gd1_15[k]
                  + pa_z[k] * hd_24[k];

        t_54[k] = f_3 * hp_18[k]
                  + pb_y[k] * ip_22[k];

        t_55[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_20[k]
                  + pa_y[k] * hd_26[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pb_x, pb_y, hp_20, hp_21, hd_28, hd_29, \
                         is0_14, is1_20, ip_23, ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_3 * hp_20[k]
                  + pa_y[k] * hd_28[k];

        t_57[k] = f_13 * hp_21[k]
                  + pb_y[k] * ip_23[k];

        t_58[k] = pa_y[k] * hd_29[k];

        t_59[k] = f_1 * is0_14[k]
                  - f_2 * is1_20[k]
                  + pb_x[k] * ip_24[k];
    }

#pragma omp simd aligned(t_60, t_61, pb_y, pb_z, hp_21, is0_14, is1_20, ip_25, \
                         ip_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * is0_14[k]
                  - f_2 * is1_20[k]
                  + pb_y[k] * ip_25[k];

        t_61[k] = f_0 * hp_21[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_20[k]
                  + pb_z[k] * ip_26[k];
    }
}

auto
compute_prim_id_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_53 = buffer.data(hd + 53);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, gd0_0, gd0_4, gd1_0, gd1_4, hd_0, \
                         hd_3, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_6[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_6, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_8[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_6[k];

        t_9[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_1, gd0_6, gd0_7, gd1_1, gd1_6, \
                         gd1_7, hd_9, hd_16, hd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd0_6[k]
                  - f_6 * gd1_6[k]
                  + pa_x[k] * hd_16[k];

        t_11[k] = f_7 * gd0_1[k]
                  - f_8 * gd1_1[k]
                  + pa_y[k] * hd_9[k];

        t_12[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_7[k]
                  + pa_x[k] * hd_18[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, gd0_2, gd1_2, hd_14, is0_3, \
                         is0_4, is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_14[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_14[k];

        t_15[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_17, hd_28, hd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_28[k];

        t_17[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_17[k];

        t_18[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_30[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, gd0_5, gd1_5, hd_26, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_20[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_26[k];

        t_21[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, gd0_14, gd1_14, hd_36, hd_38, \
                         hd_53, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_36[k];

        t_23[k] = pa_x[k] * hd_38[k];

        t_24[k] = pa_x[k] * hd_53[k];

        t_25[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hp_8, hd_38, is0_9, is1_9, ip_10, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_38[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_9, gd1_12, \
                         hd_41, hd_45, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_30[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_41[k];

        t_31[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_45[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_10, gd1_13, \
                         hd_44, hd_48, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_44[k];

        t_34[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_48[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_11, gd1_14, \
                         hd_47, hd_50, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_36[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_47[k];

        t_37[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_50[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, hp_14, hd_53, is0_14, \
                         is1_14, ip_15, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hd_53[k];

        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_41[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_23 = buffer.data(gd1 + 23);
    const auto *gd1_24 = buffer.data(gd1 + 24);
    const auto *gd1_26 = buffer.data(gd1 + 26);
    const auto *gd1_29 = buffer.data(gd1 + 29);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
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
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, gd0_0, gd1_0, hp_1, hp_2, \
                         hd_0, hd_1, hd_2, hd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_5[k] = pa_y[k] * hd_2[k];

        t_6[k] = pa_z[k] * hd_0[k];

        t_7[k] = pa_z[k] * hd_1[k];

        t_8[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_9[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_z, gd0_4, gd1_8, hd_4, \
                         hd_6, hd_8, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * gd0_4[k]
                  - f_7 * gd1_8[k]
                  + pa_x[k] * hd_8[k];

        t_11[k] = f_1 * is0_1[k]
                  - f_2 * is1_1[k]
                  + pb_z[k] * ip_3[k];

        t_12[k] = pa_z[k] * hd_4[k];

        t_13[k] = pa_y[k] * hd_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_y, gd0_0, gd0_6, gd1_0, gd1_12, \
                         hd_5, hd_12, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * gd0_0[k]
                  - f_5 * gd1_0[k]
                  + pa_z[k] * hd_5[k];

        t_15[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_4[k];

        t_16[k] = f_6 * gd0_6[k]
                  - f_7 * gd1_12[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, gd0_1, gd0_7, gd1_3, gd1_14, \
                         hd_7, hd_14, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * gd0_1[k]
                  - f_9 * gd1_3[k]
                  + pa_y[k] * hd_7[k];

        t_18[k] = f_8 * gd0_7[k]
                  - f_9 * gd1_14[k]
                  + pa_x[k] * hd_14[k];

        t_19[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, hp_3, hp_4, hd_7, hd_8, \
                         hd_9, hd_11, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * hd_7[k];

        t_21[k] = pa_z[k] * hd_8[k];

        t_22[k] = f_3 * hp_3[k]
                  + pa_z[k] * hd_9[k];

        t_23[k] = f_3 * hp_4[k]
                  + pa_y[k] * hd_11[k];

        t_24[k] = pa_y[k] * hd_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_5, gd1_16, \
                         hd_10, hd_19, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * gd0_2[k]
                  - f_9 * gd1_5[k]
                  + pa_z[k] * hd_10[k];

        t_26[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_27[k] = f_8 * gd0_8[k]
                  - f_9 * gd1_16[k]
                  + pa_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_7, gd1_18, \
                         hd_13, hd_21, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_6 * gd0_3[k]
                  - f_7 * gd1_7[k]
                  + pa_y[k] * hd_13[k];

        t_29[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_18[k]
                  + pa_x[k] * hd_21[k];

        t_30[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, gd0_5, gd1_10, hp_5, hd_13, \
                         hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_z[k] * hd_13[k];

        t_32[k] = pa_z[k] * hd_14[k];

        t_33[k] = f_3 * hp_5[k]
                  + pa_z[k] * hd_15[k];

        t_34[k] = f_4 * gd0_5[k]
                  - f_5 * gd1_10[k]
                  + pa_y[k] * hd_16[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pa_y, gd0_11, gd0_12, gd1_23, gd1_24, \
                         hp_6, hd_18, hd_19, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_4 * gd0_11[k]
                  - f_5 * gd1_23[k]
                  + pa_x[k] * hd_22[k];

        t_36[k] = f_4 * gd0_12[k]
                  - f_5 * gd1_24[k]
                  + pa_x[k] * hd_23[k];

        t_37[k] = f_3 * hp_6[k]
                  + pa_y[k] * hd_18[k];

        t_38[k] = pa_y[k] * hd_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pa_z, pb_y, gd0_5, gd0_14, gd1_10, gd1_29, \
                         hd_17, hd_25, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_6 * gd0_5[k]
                  - f_7 * gd1_10[k]
                  + pa_z[k] * hd_17[k];

        t_40[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_41[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_29[k]
                  + pa_x[k] * hd_25[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_x, pa_z, hp_7, hp_10, hp_11, hd_20, \
                         hd_26, hd_27, hd_31, hd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * hp_7[k]
                  + pa_x[k] * hd_26[k];

        t_43[k] = pa_x[k] * hd_27[k];

        t_44[k] = pa_z[k] * hd_20[k];

        t_45[k] = f_3 * hp_10[k]
                  + pa_x[k] * hd_31[k];

        t_46[k] = f_3 * hp_11[k]
                  + pa_x[k] * hd_34[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pb_x, pb_y, hp_8, hp_12, hd_39, hd_41, \
                         is0_9, is1_9, ip_9, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * hp_12[k]
                  + pa_x[k] * hd_39[k];

        t_48[k] = pa_x[k] * hd_41[k];

        t_49[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];

        t_50[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_z, hp_9, hd_26, hd_27, hd_28, is0_9, \
                         is1_9, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_52[k] = pa_z[k] * hd_26[k];

        t_53[k] = pa_z[k] * hd_27[k];

        t_54[k] = f_3 * hp_9[k]
                  + pa_z[k] * hd_28[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_18, gd1_24, \
                         hd_29, hd_33, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_56[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_18[k]
                  + pa_z[k] * hd_29[k];

        t_57[k] = f_6 * gd0_12[k]
                  - f_7 * gd1_24[k]
                  + pa_y[k] * hd_33[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_20, gd1_26, \
                         hd_32, hd_36, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_59[k] = f_8 * gd0_10[k]
                  - f_9 * gd1_20[k]
                  + pa_z[k] * hd_32[k];

        t_60[k] = f_8 * gd0_13[k]
                  - f_9 * gd1_26[k]
                  + pa_y[k] * hd_36[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_23, gd1_29, \
                         hd_35, hd_38, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_62[k] = f_6 * gd0_11[k]
                  - f_7 * gd1_23[k]
                  + pa_z[k] * hd_35[k];

        t_63[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_29[k]
                  + pa_y[k] * hd_38[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, hp_13, hd_40, hd_41, \
                         is0_14, is1_14, ip_15, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * hp_13[k]
                  + pa_y[k] * hd_40[k];

        t_65[k] = pa_y[k] * hd_41[k];

        t_66[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_67[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];
    }

#pragma omp simd aligned(t_68, pb_z, hp_14, is0_14, is1_14, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_20 = buffer.data(gd0 + 20);
    const auto *gd0_23 = buffer.data(gd0 + 23);
    const auto *gd0_24 = buffer.data(gd0 + 24);
    const auto *gd0_26 = buffer.data(gd0 + 26);
    const auto *gd0_29 = buffer.data(gd0 + 29);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_20 = buffer.data(gd1 + 20);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, gd0_0, gd1_0, hp_1, hp_2, hd_0, hd_1, \
                         hd_2, hd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_5[k] = pa_z[k] * hd_0[k];

        t_6[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_7[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_z, pb_z, gd0_0, gd0_8, gd1_0, gd1_6, hd_4, \
                         hd_6, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * gd0_8[k]
                 - f_7 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_9[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_10[k] = f_4 * gd0_0[k]
                  - f_5 * gd1_0[k]
                  + pa_z[k] * hd_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, pb_y, gd0_3, gd0_12, gd1_3, gd1_8, \
                         hd_5, hd_8, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_4[k];

        t_12[k] = f_6 * gd0_12[k]
                  - f_7 * gd1_8[k]
                  + pa_x[k] * hd_8[k];

        t_13[k] = f_8 * gd0_3[k]
                  - f_9 * gd1_3[k]
                  + pa_y[k] * hd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_z, gd0_14, gd1_9, hd_7, hd_10, \
                         is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * gd0_14[k]
                  - f_9 * gd1_9[k]
                  + pa_x[k] * hd_10[k];

        t_15[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_16[k] = pa_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_z, pb_y, gd0_5, gd0_16, gd1_4, gd1_10, \
                         hd_7, hd_13, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * gd0_5[k]
                  - f_9 * gd1_4[k]
                  + pa_z[k] * hd_7[k];

        t_18[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_19[k] = f_8 * gd0_16[k]
                  - f_9 * gd1_10[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pa_y, pb_z, gd0_7, gd0_18, gd1_5, gd1_12, \
                         hd_9, hd_14, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * gd0_7[k]
                  - f_7 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_21[k] = f_4 * gd0_18[k]
                  - f_5 * gd1_12[k]
                  + pa_x[k] * hd_14[k];

        t_22[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, gd0_10, gd0_23, gd0_24, gd1_7, \
                         gd1_15, gd1_16, hd_11, hd_12, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * gd0_10[k]
                  - f_5 * gd1_7[k]
                  + pa_y[k] * hd_11[k];

        t_24[k] = f_4 * gd0_23[k]
                  - f_5 * gd1_15[k]
                  + pa_x[k] * hd_15[k];

        t_25[k] = f_4 * gd0_24[k]
                  - f_5 * gd1_16[k]
                  + pa_x[k] * hd_16[k];

        t_26[k] = pa_y[k] * hd_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_z, pb_y, gd0_10, gd0_29, gd1_7, gd1_20, \
                         hd_12, hd_17, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * gd0_10[k]
                  - f_7 * gd1_7[k]
                  + pa_z[k] * hd_12[k];

        t_28[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_29[k] = f_4 * gd0_29[k]
                  - f_5 * gd1_20[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_x, hp_7, hd_18, hd_19, hd_22, \
                         hd_23, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hp_7[k]
                  + pa_x[k] * hd_18[k];

        t_31[k] = pa_x[k] * hd_19[k];

        t_32[k] = pa_x[k] * hd_22[k];

        t_33[k] = pa_x[k] * hd_23[k];

        t_34[k] = pa_x[k] * hd_24[k];

        t_35[k] = pa_x[k] * hd_25[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_x, pb_y, hp_8, hp_12, hd_27, hd_29, \
                         is0_9, is1_9, ip_9, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * hp_12[k]
                  + pa_x[k] * hd_27[k];

        t_37[k] = pa_x[k] * hd_29[k];

        t_38[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];

        t_39[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_x, pb_z, hp_9, hd_19, hd_20, is0_9, \
                         is0_10, is1_9, is1_10, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_41[k] = pa_z[k] * hd_19[k];

        t_42[k] = f_3 * hp_9[k]
                  + pa_z[k] * hd_20[k];

        t_43[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_x, gd0_18, gd0_24, gd1_12, gd1_16, \
                         hd_21, hd_23, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_4 * gd0_18[k]
                  - f_5 * gd1_12[k]
                  + pa_z[k] * hd_21[k];

        t_45[k] = f_6 * gd0_24[k]
                  - f_7 * gd1_16[k]
                  + pa_y[k] * hd_23[k];

        t_46[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pa_z, pb_x, gd0_20, gd0_26, gd1_14, gd1_17, \
                         hd_22, hd_25, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_8 * gd0_20[k]
                  - f_9 * gd1_14[k]
                  + pa_z[k] * hd_22[k];

        t_48[k] = f_8 * gd0_26[k]
                  - f_9 * gd1_17[k]
                  + pa_y[k] * hd_25[k];

        t_49[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pa_z, gd0_23, gd0_29, gd1_15, gd1_20, \
                         hp_13, hd_24, hd_26, hd_28, hd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * gd0_23[k]
                  - f_7 * gd1_15[k]
                  + pa_z[k] * hd_24[k];

        t_51[k] = f_4 * gd0_29[k]
                  - f_5 * gd1_20[k]
                  + pa_y[k] * hd_26[k];

        t_52[k] = f_3 * hp_13[k]
                  + pa_y[k] * hd_28[k];

        t_53[k] = pa_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, pb_y, pb_z, hp_14, is0_14, is1_14, ip_15, \
                         ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_55[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_56[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_7 = buffer.data(is0 + 7);
    const auto *is0_8 = buffer.data(is0 + 8);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hp_0, is0_0, is1_0, ip_0, ip_1, \
                         ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, gd0_0, gd0_4, gd1_0, gd1_4, hd_3, \
                         hd_6, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_6[k];

        t_5[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_y, gd0_0, gd0_6, gd1_0, gd1_6, hd_4, \
                         hd_10, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_7[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];

        t_8[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_z, gd0_1, gd0_7, gd1_1, gd1_7, hd_5, \
                         hd_12, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_5[k];

        t_10[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_7[k]
                  + pa_x[k] * hd_12[k];

        t_11[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_2, gd1_8, hd_8, \
                         hd_16, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_8[k];

        t_13[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_14[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_16[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_3, gd1_9, \
                         hd_11, hd_17, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_11[k];

        t_16[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_17[k];

        t_17[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_z, pb_y, gd0_5, gd0_14, gd1_5, gd1_14, \
                         hd_14, hd_18, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_14[k];

        t_19[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_20[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_18[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pb_x, pb_y, hp_8, hd_20, hd_32, is0_7, \
                         is1_9, ip_9, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_x[k] * hd_20[k];

        t_22[k] = pa_x[k] * hd_32[k];

        t_23[k] = f_1 * is0_7[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];

        t_24[k] = f_0 * hp_8[k]
                  + f_1 * is0_7[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_z, pb_x, pb_z, gd0_9, gd1_9, hd_22, is0_7, \
                         is0_8, is1_9, is1_10, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_1 * is0_7[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_26[k] = f_1 * is0_8[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_27[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_22[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, pb_x, gd0_10, gd0_12, gd1_10, gd1_12, \
                         hd_24, hd_25, is0_9, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_25[k];

        t_29[k] = f_1 * is0_9[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_30[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_24[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, gd0_11, gd0_13, gd1_11, gd1_13, \
                         hd_27, hd_28, is0_10, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_28[k];

        t_32[k] = f_1 * is0_10[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_33[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_27[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pb_x, pb_y, gd0_14, gd1_14, hd_29, \
                         hd_32, is0_11, is1_14, ip_15, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_29[k];

        t_35[k] = pa_y[k] * hd_32[k];

        t_36[k] = f_1 * is0_11[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_37[k] = f_1 * is0_11[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];
    }

#pragma omp simd aligned(t_38, pb_z, hp_14, is0_11, is1_14, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * hp_14[k]
                  + f_1 * is0_11[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_35 = buffer.data(hd + 35);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, gd0_0, gd0_4, gd1_0, gd1_6, hd_0, \
                         hd_3, hd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_6[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_7[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_8[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_9[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_1, gd0_6, gd0_7, gd1_3, gd1_10, \
                         gd1_11, hd_6, hd_11, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd0_6[k]
                  - f_6 * gd1_10[k]
                  + pa_x[k] * hd_11[k];

        t_11[k] = f_7 * gd0_1[k]
                  - f_8 * gd1_3[k]
                  + pa_y[k] * hd_6[k];

        t_12[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_11[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, gd0_2, gd1_4, hd_9, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_14[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_9[k];

        t_15[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_5, gd1_12, \
                         gd1_14, hd_12, hd_17, hd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_12[k]
                  + pa_x[k] * hd_17[k];

        t_17[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_12[k];

        t_18[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, gd0_5, gd1_8, hd_15, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_20[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_8[k]
                  + pa_z[k] * hd_15[k];

        t_21[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, gd0_14, gd1_23, hd_19, hd_21, \
                         hd_35, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_23[k]
                  + pa_x[k] * hd_19[k];

        t_23[k] = pa_x[k] * hd_21[k];

        t_24[k] = pa_x[k] * hd_35[k];

        t_25[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hp_8, hd_21, is0_9, is1_9, ip_10, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_21[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_14, gd1_19, \
                         hd_23, hd_27, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_30[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_14[k]
                  + pa_z[k] * hd_23[k];

        t_31[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_19[k]
                  + pa_y[k] * hd_27[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_16, gd1_20, \
                         hd_26, hd_30, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_16[k]
                  + pa_z[k] * hd_26[k];

        t_34[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_20[k]
                  + pa_y[k] * hd_30[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_18, gd1_23, \
                         hd_29, hd_32, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_36[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_18[k]
                  + pa_z[k] * hd_29[k];

        t_37[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_23[k]
                  + pa_y[k] * hd_32[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, hp_14, hd_35, is0_14, \
                         is1_14, ip_15, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hd_35[k];

        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_41[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_20 = buffer.data(gd0 + 20);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, gd0_0, gd0_6, gd1_0, gd1_6, \
                         hp_2, hd_0, hd_2, hd_3, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_6[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_7[k] = f_6 * gd0_6[k]
                 - f_7 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_9[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_10[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, gd0_3, gd0_10, gd0_11, gd1_3, gd1_10, \
                         gd1_11, hd_5, hd_10, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_6 * gd0_10[k]
                  - f_7 * gd1_10[k]
                  + pa_x[k] * hd_10[k];

        t_12[k] = f_8 * gd0_3[k]
                  - f_9 * gd1_3[k]
                  + pa_y[k] * hd_5[k];

        t_13[k] = f_8 * gd0_11[k]
                  - f_9 * gd1_11[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_y, pb_z, gd0_4, gd1_4, hd_8, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_15[k] = f_8 * gd0_4[k]
                  - f_9 * gd1_4[k]
                  + pa_z[k] * hd_8[k];

        t_16[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, gd0_5, gd0_12, gd0_14, gd1_5, gd1_12, \
                         gd1_14, hd_11, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * gd0_12[k]
                  - f_9 * gd1_12[k]
                  + pa_x[k] * hd_16[k];

        t_18[k] = f_6 * gd0_5[k]
                  - f_7 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_19[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_14[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, pb_z, gd0_8, gd1_8, hd_14, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_21[k] = f_6 * gd0_8[k]
                  - f_7 * gd1_8[k]
                  + pa_z[k] * hd_14[k];

        t_22[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, gd0_23, gd1_23, hd_18, hd_20, \
                         hd_32, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * gd0_23[k]
                  - f_5 * gd1_23[k]
                  + pa_x[k] * hd_18[k];

        t_24[k] = pa_x[k] * hd_20[k];

        t_25[k] = pa_x[k] * hd_32[k];

        t_26[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_z, pb_y, pb_z, hp_8, hp_9, hd_20, hd_21, \
                         is0_9, is1_9, ip_10, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_28[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_29[k] = pa_z[k] * hd_20[k];

        t_30[k] = f_3 * hp_9[k]
                  + pa_z[k] * hd_21[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, gd0_14, gd0_19, gd1_14, gd1_19, \
                         hd_22, hd_25, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_32[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_14[k]
                  + pa_z[k] * hd_22[k];

        t_33[k] = f_6 * gd0_19[k]
                  - f_7 * gd1_19[k]
                  + pa_y[k] * hd_25[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_x, gd0_16, gd0_20, gd1_16, gd1_20, \
                         hd_24, hd_28, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_35[k] = f_8 * gd0_16[k]
                  - f_9 * gd1_16[k]
                  + pa_z[k] * hd_24[k];

        t_36[k] = f_8 * gd0_20[k]
                  - f_9 * gd1_20[k]
                  + pa_y[k] * hd_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, gd0_18, gd0_23, gd1_18, gd1_23, \
                         hd_27, hd_29, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_38[k] = f_6 * gd0_18[k]
                  - f_7 * gd1_18[k]
                  + pa_z[k] * hd_27[k];

        t_39[k] = f_4 * gd0_23[k]
                  - f_5 * gd1_23[k]
                  + pa_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_x, pb_y, hp_13, hd_31, hd_32, \
                         is0_14, is1_14, ip_15, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_3 * hp_13[k]
                  + pa_y[k] * hd_31[k];

        t_41[k] = pa_y[k] * hd_32[k];

        t_42[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_43[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];
    }

#pragma omp simd aligned(t_44, pb_z, hp_14, is0_14, is1_14, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_20 = buffer.data(gd0 + 20);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_20 = buffer.data(gd1 + 20);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
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
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_26 = buffer.data(hd + 26);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, gd0_0, gd0_6, gd1_0, gd1_6, hd_0, \
                         hd_3, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_8[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_9[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_3, gd0_10, gd0_11, gd1_3, gd1_8, \
                         gd1_9, hd_5, hd_8, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd0_10[k]
                  - f_6 * gd1_8[k]
                  + pa_x[k] * hd_8[k];

        t_11[k] = f_7 * gd0_3[k]
                  - f_8 * gd1_3[k]
                  + pa_y[k] * hd_5[k];

        t_12[k] = f_7 * gd0_11[k]
                  - f_8 * gd1_9[k]
                  + pa_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, gd0_4, gd1_4, hd_7, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_14[k] = f_7 * gd0_4[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_7[k];

        t_15[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, gd0_5, gd0_12, gd0_14, gd1_5, gd1_10, \
                         gd1_12, hd_9, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gd0_12[k]
                  - f_8 * gd1_10[k]
                  + pa_x[k] * hd_12[k];

        t_17[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_18[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, gd0_8, gd1_7, hd_11, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_20[k] = f_5 * gd0_8[k]
                  - f_6 * gd1_7[k]
                  + pa_z[k] * hd_11[k];

        t_21[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, gd0_23, gd1_20, hd_14, hd_16, \
                         hd_26, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_20[k]
                  + pa_x[k] * hd_14[k];

        t_23[k] = pa_x[k] * hd_16[k];

        t_24[k] = pa_x[k] * hd_26[k];

        t_25[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hp_8, hd_16, is0_9, is1_9, ip_10, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_16[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, gd0_14, gd0_19, gd1_12, gd1_16, \
                         hd_18, hd_20, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_30[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_12[k]
                  + pa_z[k] * hd_18[k];

        t_31[k] = f_5 * gd0_19[k]
                  - f_6 * gd1_16[k]
                  + pa_y[k] * hd_20[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_x, gd0_16, gd0_20, gd1_14, gd1_17, \
                         hd_19, hd_22, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_7 * gd0_16[k]
                  - f_8 * gd1_14[k]
                  + pa_z[k] * hd_19[k];

        t_34[k] = f_7 * gd0_20[k]
                  - f_8 * gd1_17[k]
                  + pa_y[k] * hd_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, gd0_18, gd0_23, gd1_15, gd1_20, \
                         hd_21, hd_23, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_36[k] = f_5 * gd0_18[k]
                  - f_6 * gd1_15[k]
                  + pa_z[k] * hd_21[k];

        t_37[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_20[k]
                  + pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, hp_14, hd_26, is0_14, \
                         is1_14, ip_15, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hd_26[k];

        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_41[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, gd0_0, gd0_4, gd1_0, gd1_4, hd_0, \
                         hd_3, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_6[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_8[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_9[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_1, gd0_6, gd0_7, gd1_1, gd1_6, \
                         gd1_7, hd_5, hd_10, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd0_6[k]
                  - f_6 * gd1_6[k]
                  + pa_x[k] * hd_10[k];

        t_11[k] = f_7 * gd0_1[k]
                  - f_8 * gd1_1[k]
                  + pa_y[k] * hd_5[k];

        t_12[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_7[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, gd0_2, gd1_2, hd_8, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_14[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_8[k];

        t_15[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_11, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_16[k];

        t_17[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_11[k];

        t_18[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, gd0_5, gd1_5, hd_14, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_20[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_14[k];

        t_21[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, gd0_14, gd1_14, hd_18, hd_20, \
                         hd_32, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_18[k];

        t_23[k] = pa_x[k] * hd_20[k];

        t_24[k] = pa_x[k] * hd_32[k];

        t_25[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hp_8, hd_20, is0_9, is1_9, ip_10, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_9, gd1_12, \
                         hd_22, hd_25, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_30[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_22[k];

        t_31[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_10, gd1_13, \
                         hd_24, hd_28, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_24[k];

        t_34[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_28[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_11, gd1_14, \
                         hd_27, hd_29, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_36[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_27[k];

        t_37[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, hp_14, hd_32, is0_14, \
                         is1_14, ip_15, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hd_32[k];

        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_41[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, gd0_0, gd0_4, gd1_0, gd1_6, hd_0, \
                         hd_3, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_6[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_8[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_9[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_1, gd0_6, gd0_7, gd1_3, gd1_10, \
                         gd1_11, hd_5, hd_10, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd0_6[k]
                  - f_6 * gd1_10[k]
                  + pa_x[k] * hd_10[k];

        t_11[k] = f_7 * gd0_1[k]
                  - f_8 * gd1_3[k]
                  + pa_y[k] * hd_5[k];

        t_12[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_11[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, gd0_2, gd1_4, hd_8, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_14[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_8[k];

        t_15[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_5, gd1_12, \
                         gd1_14, hd_11, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_12[k]
                  + pa_x[k] * hd_16[k];

        t_17[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_18[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, gd0_5, gd1_8, hd_14, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_20[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_8[k]
                  + pa_z[k] * hd_14[k];

        t_21[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, gd0_14, gd1_23, hd_18, hd_20, \
                         hd_32, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_23[k]
                  + pa_x[k] * hd_18[k];

        t_23[k] = pa_x[k] * hd_20[k];

        t_24[k] = pa_x[k] * hd_32[k];

        t_25[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hp_8, hd_20, is0_9, is1_9, ip_10, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_14, gd1_19, \
                         hd_22, hd_25, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_30[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_14[k]
                  + pa_z[k] * hd_22[k];

        t_31[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_19[k]
                  + pa_y[k] * hd_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_16, gd1_20, \
                         hd_24, hd_28, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_16[k]
                  + pa_z[k] * hd_24[k];

        t_34[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_20[k]
                  + pa_y[k] * hd_28[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_18, gd1_23, \
                         hd_27, hd_29, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_36[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_18[k]
                  + pa_z[k] * hd_27[k];

        t_37[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_23[k]
                  + pa_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, hp_14, hd_32, is0_14, \
                         is1_14, ip_15, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hd_32[k];

        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_41[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_18 = buffer.data(gd0 + 18);
    const auto *gd0_19 = buffer.data(gd0 + 19);
    const auto *gd0_20 = buffer.data(gd0 + 20);
    const auto *gd0_23 = buffer.data(gd0 + 23);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_18 = buffer.data(gd1 + 18);
    const auto *gd1_19 = buffer.data(gd1 + 19);
    const auto *gd1_20 = buffer.data(gd1 + 20);
    const auto *gd1_23 = buffer.data(gd1 + 23);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, gd0_0, gd0_6, gd1_0, gd1_6, hd_0, \
                         hd_3, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_8[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_9[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_3, gd0_10, gd0_11, gd1_3, gd1_10, \
                         gd1_11, hd_5, hd_10, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd0_10[k]
                  - f_6 * gd1_10[k]
                  + pa_x[k] * hd_10[k];

        t_11[k] = f_7 * gd0_3[k]
                  - f_8 * gd1_3[k]
                  + pa_y[k] * hd_5[k];

        t_12[k] = f_7 * gd0_11[k]
                  - f_8 * gd1_11[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, gd0_4, gd1_4, hd_8, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_14[k] = f_7 * gd0_4[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_8[k];

        t_15[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, gd0_5, gd0_12, gd0_14, gd1_5, gd1_12, \
                         gd1_14, hd_11, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gd0_12[k]
                  - f_8 * gd1_12[k]
                  + pa_x[k] * hd_16[k];

        t_17[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_18[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, gd0_8, gd1_8, hd_14, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_20[k] = f_5 * gd0_8[k]
                  - f_6 * gd1_8[k]
                  + pa_z[k] * hd_14[k];

        t_21[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, gd0_23, gd1_23, hd_18, hd_20, \
                         hd_32, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pa_x[k] * hd_18[k];

        t_23[k] = pa_x[k] * hd_20[k];

        t_24[k] = pa_x[k] * hd_32[k];

        t_25[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hp_8, hd_20, is0_9, is1_9, ip_10, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, gd0_14, gd0_19, gd1_14, gd1_19, \
                         hd_22, hd_25, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_30[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_z[k] * hd_22[k];

        t_31[k] = f_5 * gd0_19[k]
                  - f_6 * gd1_19[k]
                  + pa_y[k] * hd_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_x, gd0_16, gd0_20, gd1_16, gd1_20, \
                         hd_24, hd_28, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_7 * gd0_16[k]
                  - f_8 * gd1_16[k]
                  + pa_z[k] * hd_24[k];

        t_34[k] = f_7 * gd0_20[k]
                  - f_8 * gd1_20[k]
                  + pa_y[k] * hd_28[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, gd0_18, gd0_23, gd1_18, gd1_23, \
                         hd_27, hd_29, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_36[k] = f_5 * gd0_18[k]
                  - f_6 * gd1_18[k]
                  + pa_z[k] * hd_27[k];

        t_37[k] = f_3 * gd0_23[k]
                  - f_4 * gd1_23[k]
                  + pa_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, hp_14, hd_32, is0_14, \
                         is1_14, ip_15, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hd_32[k];

        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_41[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_14 = buffer.data(is1 + 14);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gd0_0, gd1_0, hp_0, hd_0, hd_1, \
                         is0_0, is1_0, ip_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = pa_y[k] * hd_0[k];

        t_2[k] = pa_z[k] * hd_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, gd0_0, gd0_4, gd0_6, gd1_0, gd1_4, gd1_6, \
                         hd_2, hd_4, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, gd0_1, gd0_2, gd0_7, gd1_1, \
                         gd1_2, gd1_7, hd_3, hd_5, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_8[k] = f_7 * gd0_7[k]
                 - f_8 * gd1_7[k]
                 + pa_x[k] * hd_8[k];

        t_9[k] = pa_y[k] * hd_5[k];

        t_10[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_7, hd_11, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_11[k];

        t_12[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];

        t_13[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_y, gd0_5, gd0_11, gd0_12, gd1_5, \
                         gd1_11, gd1_12, hd_9, hd_10, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_15[k] = f_3 * gd0_11[k]
                  - f_4 * gd1_11[k]
                  + pa_x[k] * hd_13[k];

        t_16[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_14[k];

        t_17[k] = pa_y[k] * hd_10[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pa_z, gd0_5, gd0_14, gd1_5, \
                         gd1_14, hd_10, hd_15, hd_16, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_10[k];

        t_19[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_15[k];

        t_20[k] = pa_x[k] * hd_16[k];

        t_21[k] = pa_x[k] * hd_18[k];

        t_22[k] = pa_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pa_z, pb_y, hp_1, hd_16, hd_20, \
                         hd_21, hd_23, is0_9, is1_9, ip_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * hd_20[k];

        t_24[k] = pa_x[k] * hd_21[k];

        t_25[k] = pa_x[k] * hd_23[k];

        t_26[k] = f_0 * hp_1[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_1[k];

        t_27[k] = pa_z[k] * hd_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, gd0_9, gd0_10, gd0_12, gd1_9, gd1_10, \
                         gd1_12, hd_17, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_17[k];

        t_29[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_19[k];

        t_30[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, gd0_11, gd0_13, gd0_14, gd1_11, \
                         gd1_13, gd1_14, hd_20, hd_21, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_21[k];

        t_32[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_20[k];

        t_33[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_22[k];

        t_34[k] = pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_35, pb_z, hp_2, is0_14, is1_14, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * hp_2[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_2[k];
    }
}

auto
compute_prim_id_electron_repulsion_23(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_18 = buffer.data(is1 + 18);
    const auto *is1_24 = buffer.data(is1 + 24);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_27 = buffer.data(ip + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gd0_0, gd1_0, hp_0, hd_0, hd_1, \
                         is0_0, is1_0, ip_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = pa_y[k] * hd_0[k];

        t_2[k] = pa_z[k] * hd_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, gd0_0, gd0_4, gd0_6, gd1_0, gd1_4, gd1_6, \
                         hd_2, hd_4, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, gd0_1, gd0_2, gd0_7, gd1_1, \
                         gd1_2, gd1_7, hd_3, hd_5, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_8[k] = f_7 * gd0_7[k]
                 - f_8 * gd1_7[k]
                 + pa_x[k] * hd_8[k];

        t_9[k] = pa_y[k] * hd_5[k];

        t_10[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_7, hd_11, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_11[k];

        t_12[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];

        t_13[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_y, gd0_5, gd0_11, gd0_12, gd1_5, \
                         gd1_11, gd1_12, hd_9, hd_10, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_15[k] = f_3 * gd0_11[k]
                  - f_4 * gd1_11[k]
                  + pa_x[k] * hd_13[k];

        t_16[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_14[k];

        t_17[k] = pa_y[k] * hd_10[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pa_z, gd0_5, gd0_14, gd1_5, \
                         gd1_14, hd_10, hd_15, hd_16, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_10[k];

        t_19[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_15[k];

        t_20[k] = pa_x[k] * hd_16[k];

        t_21[k] = pa_x[k] * hd_18[k];

        t_22[k] = pa_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pa_z, pb_y, hp_1, hd_16, hd_20, \
                         hd_21, hd_23, is0_9, is1_18, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * hd_20[k];

        t_24[k] = pa_x[k] * hd_21[k];

        t_25[k] = pa_x[k] * hd_23[k];

        t_26[k] = f_0 * hp_1[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_18[k]
                  + pb_y[k] * ip_14[k];

        t_27[k] = pa_z[k] * hd_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, gd0_9, gd0_10, gd0_12, gd1_9, gd1_10, \
                         gd1_12, hd_17, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_17[k];

        t_29[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_19[k];

        t_30[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, gd0_11, gd0_13, gd0_14, gd1_11, \
                         gd1_13, gd1_14, hd_20, hd_21, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_21[k];

        t_32[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_20[k];

        t_33[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_22[k];

        t_34[k] = pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_35, pb_z, hp_2, is0_14, is1_24, ip_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * hp_2[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_24[k]
                  + pb_z[k] * ip_27[k];
    }
}

auto
compute_prim_id_electron_repulsion_24(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_18 = buffer.data(is0 + 18);
    const auto *is0_24 = buffer.data(is0 + 24);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_14 = buffer.data(is1 + 14);
    const auto *is1_20 = buffer.data(is1 + 20);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, gd0_0, gd0_4, gd1_0, gd1_4, \
                         hd_0, hd_1, hd_2, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];

        t_6[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_7[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_y, gd0_1, gd0_6, gd0_7, gd1_1, gd1_6, \
                         gd1_7, hd_3, hd_5, hd_6, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_9[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_10[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_7[k]
                  + pa_x[k] * hd_8[k];

        t_11[k] = pa_y[k] * hd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_y, pa_z, gd0_2, gd0_3, gd0_8, gd1_2, \
                         gd1_3, gd1_8, hd_5, hd_7, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_5[k];

        t_13[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_11[k];

        t_14[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, gd0_5, gd0_9, gd0_11, gd1_5, gd1_9, \
                         gd1_11, hd_9, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_12[k];

        t_16[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_17[k] = f_3 * gd0_11[k]
                  - f_4 * gd1_11[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pa_z, gd0_5, gd0_12, gd0_14, \
                         gd1_5, gd1_12, gd1_14, hd_10, hd_14, hd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_14[k];

        t_19[k] = pa_y[k] * hd_10[k];

        t_20[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_10[k];

        t_21[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pa_x, hd_16, hd_18, hd_19, hd_20, \
                         hd_21, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_x[k] * hd_16[k];

        t_23[k] = pa_x[k] * hd_18[k];

        t_24[k] = pa_x[k] * hd_19[k];

        t_25[k] = pa_x[k] * hd_20[k];

        t_26[k] = pa_x[k] * hd_21[k];

        t_27[k] = pa_x[k] * hd_23[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_z, pb_x, pb_y, pb_z, hp_1, hd_16, is0_18, \
                         is1_14, ip_3, ip_4, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * is0_18[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_3[k];

        t_29[k] = f_0 * hp_1[k]
                  + f_1 * is0_18[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_4[k];

        t_30[k] = f_1 * is0_18[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_5[k];

        t_31[k] = pa_z[k] * hd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, gd0_9, gd0_10, gd0_12, gd1_9, gd1_10, \
                         gd1_12, hd_17, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_17[k];

        t_33[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_19[k];

        t_34[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_y, pa_z, gd0_11, gd0_13, gd0_14, gd1_11, \
                         gd1_13, gd1_14, hd_20, hd_21, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_21[k];

        t_36[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_20[k];

        t_37[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_22[k];

        t_38[k] = pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, pb_y, pb_z, hp_2, is0_24, is1_20, ip_6, ip_7, \
                         ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * is0_24[k]
                  - f_2 * is1_20[k]
                  + pb_x[k] * ip_6[k];

        t_40[k] = f_1 * is0_24[k]
                  - f_2 * is1_20[k]
                  + pb_y[k] * ip_7[k];

        t_41[k] = f_0 * hp_2[k]
                  + f_1 * is0_24[k]
                  - f_2 * is1_20[k]
                  + pb_z[k] * ip_8[k];
    }
}

auto
compute_prim_id_electron_repulsion_25(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_14 = buffer.data(is1 + 14);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_32 = buffer.data(ip + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gd0_0, gd1_0, hp_0, hd_0, hd_1, \
                         is0_0, is1_0, ip_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = pa_y[k] * hd_0[k];

        t_2[k] = pa_z[k] * hd_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, gd0_0, gd0_4, gd0_6, gd1_0, gd1_4, gd1_6, \
                         hd_2, hd_4, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, gd0_1, gd0_2, gd0_7, gd1_1, \
                         gd1_2, gd1_7, hd_3, hd_5, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_8[k] = f_7 * gd0_7[k]
                 - f_8 * gd1_7[k]
                 + pa_x[k] * hd_8[k];

        t_9[k] = pa_y[k] * hd_5[k];

        t_10[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_7, hd_11, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_11[k];

        t_12[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];

        t_13[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_y, gd0_5, gd0_11, gd0_12, gd1_5, \
                         gd1_11, gd1_12, hd_9, hd_10, hd_13, hd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_15[k] = f_3 * gd0_11[k]
                  - f_4 * gd1_11[k]
                  + pa_x[k] * hd_13[k];

        t_16[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_14[k];

        t_17[k] = pa_y[k] * hd_10[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_x, pa_z, gd0_5, gd0_14, gd1_5, \
                         gd1_14, hd_10, hd_15, hd_16, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_10[k];

        t_19[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_15[k];

        t_20[k] = pa_x[k] * hd_16[k];

        t_21[k] = pa_x[k] * hd_18[k];

        t_22[k] = pa_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pa_x, pa_z, pb_y, hp_1, hd_16, hd_20, \
                         hd_21, hd_23, is0_9, is1_9, ip_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * hd_20[k];

        t_24[k] = pa_x[k] * hd_21[k];

        t_25[k] = pa_x[k] * hd_23[k];

        t_26[k] = f_0 * hp_1[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_20[k];

        t_27[k] = pa_z[k] * hd_16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_y, pa_z, gd0_9, gd0_10, gd0_12, gd1_9, gd1_10, \
                         gd1_12, hd_17, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_17[k];

        t_29[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_19[k];

        t_30[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, gd0_11, gd0_13, gd0_14, gd1_11, \
                         gd1_13, gd1_14, hd_20, hd_21, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_21[k];

        t_32[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_20[k];

        t_33[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_22[k];

        t_34[k] = pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_35, pb_z, hp_2, is0_14, is1_14, ip_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * hp_2[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_32[k];
    }
}

auto
compute_prim_id_electron_repulsion_26(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.5 / p;
    const auto f_11 = 2.0 / p;
    const auto f_12 = 1.5 / p;
    const auto f_13 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_7 = buffer.data(is1 + 7);
    const auto *is1_8 = buffer.data(is1 + 8);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_16 = buffer.data(is1 + 16);
    const auto *is1_18 = buffer.data(is1 + 18);
    const auto *is1_19 = buffer.data(is1 + 19);
    const auto *is1_20 = buffer.data(is1 + 20);
    const auto *is1_22 = buffer.data(is1 + 22);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, gd0_0, gd1_0, hp_1, hp_2, hd_0, hd_1, \
                         hd_2, hd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_5[k] = pa_z[k] * hd_0[k];

        t_6[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_7[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_z, pb_z, gd0_0, gd0_4, gd1_0, gd1_4, hd_4, \
                         hd_6, is0_1, is1_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * gd0_4[k]
                 - f_7 * gd1_4[k]
                 + pa_x[k] * hd_6[k];

        t_9[k] = f_1 * is0_1[k]
                 - f_2 * is1_3[k]
                 + pb_z[k] * ip_4[k];

        t_10[k] = f_4 * gd0_0[k]
                  - f_5 * gd1_0[k]
                  + pa_z[k] * hd_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, pb_y, gd0_1, gd0_6, gd1_1, gd1_6, hd_5, \
                         hd_8, is0_2, is1_4, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * is0_2[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_5[k];

        t_12[k] = f_6 * gd0_6[k]
                  - f_7 * gd1_6[k]
                  + pa_x[k] * hd_8[k];

        t_13[k] = f_8 * gd0_1[k]
                  - f_9 * gd1_1[k]
                  + pa_y[k] * hd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_z, gd0_7, gd1_7, hd_7, hd_10, is0_3, \
                         is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * gd0_7[k]
                  - f_9 * gd1_7[k]
                  + pa_x[k] * hd_10[k];

        t_15[k] = f_1 * is0_3[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_16[k] = pa_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_2, gd1_8, hd_7, \
                         hd_13, is0_4, is1_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * gd0_2[k]
                  - f_9 * gd1_2[k]
                  + pa_z[k] * hd_7[k];

        t_18[k] = f_1 * is0_4[k]
                  - f_2 * is1_7[k]
                  + pb_y[k] * ip_8[k];

        t_19[k] = f_8 * gd0_8[k]
                  - f_9 * gd1_8[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_3, gd1_9, hd_9, \
                         hd_14, is0_5, is1_8, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * gd0_3[k]
                  - f_7 * gd1_3[k]
                  + pa_y[k] * hd_9[k];

        t_21[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_9[k]
                  + pa_x[k] * hd_14[k];

        t_22[k] = f_1 * is0_5[k]
                  - f_2 * is1_8[k]
                  + pb_z[k] * ip_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, gd0_5, gd0_11, gd0_12, gd1_5, \
                         gd1_11, gd1_12, hd_11, hd_12, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * gd0_5[k]
                  - f_5 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_24[k] = f_4 * gd0_11[k]
                  - f_5 * gd1_11[k]
                  + pa_x[k] * hd_15[k];

        t_25[k] = f_4 * gd0_12[k]
                  - f_5 * gd1_12[k]
                  + pa_x[k] * hd_16[k];

        t_26[k] = pa_y[k] * hd_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_z, pb_y, gd0_5, gd0_14, gd1_5, gd1_14, \
                         hd_12, hd_17, is0_6, is1_11, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * gd0_5[k]
                  - f_7 * gd1_5[k]
                  + pa_z[k] * hd_12[k];

        t_28[k] = f_1 * is0_6[k]
                  - f_2 * is1_11[k]
                  + pb_y[k] * ip_11[k];

        t_29[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_14[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_x, hp_6, hd_18, hd_19, hd_22, \
                         hd_23, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hp_6[k]
                  + pa_x[k] * hd_18[k];

        t_31[k] = pa_x[k] * hd_19[k];

        t_32[k] = pa_x[k] * hd_22[k];

        t_33[k] = pa_x[k] * hd_23[k];

        t_34[k] = pa_x[k] * hd_24[k];

        t_35[k] = pa_x[k] * hd_25[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_x, pb_y, hp_7, hp_13, hd_27, hd_29, \
                         is0_9, is1_16, ip_13, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * hp_13[k]
                  + pa_x[k] * hd_27[k];

        t_37[k] = pa_x[k] * hd_29[k];

        t_38[k] = f_1 * is0_9[k]
                  - f_2 * is1_16[k]
                  + pb_x[k] * ip_13[k];

        t_39[k] = f_0 * hp_7[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_16[k]
                  + pb_y[k] * ip_14[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_y, pb_z, hp_8, hp_9, hd_19, hd_20, \
                         is0_9, is1_16, ip_15, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * is0_9[k]
                  - f_2 * is1_16[k]
                  + pb_z[k] * ip_15[k];

        t_41[k] = pa_z[k] * hd_19[k];

        t_42[k] = f_10 * hp_9[k]
                  + pb_y[k] * ip_16[k];

        t_43[k] = f_3 * hp_8[k]
                  + pa_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_z, pb_x, pb_y, gd0_9, gd1_9, hp_10, hd_21, \
                         is0_10, is1_18, ip_17, ip_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_1 * is0_10[k]
                  - f_2 * is1_18[k]
                  + pb_x[k] * ip_17[k];

        t_45[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_9[k]
                  + pa_z[k] * hd_21[k];

        t_46[k] = f_11 * hp_10[k]
                  + pb_y[k] * ip_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pa_z, pb_x, gd0_10, gd0_12, gd1_10, gd1_12, \
                         hd_22, hd_23, is0_11, is1_19, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * gd0_12[k]
                  - f_7 * gd1_12[k]
                  + pa_y[k] * hd_23[k];

        t_48[k] = f_1 * is0_11[k]
                  - f_2 * is1_19[k]
                  + pb_x[k] * ip_19[k];

        t_49[k] = f_8 * gd0_10[k]
                  - f_9 * gd1_10[k]
                  + pa_z[k] * hd_22[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_y, pb_x, pb_y, gd0_13, gd1_13, hp_11, hd_25, \
                         is0_12, is1_20, ip_20, ip_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_12 * hp_11[k]
                  + pb_y[k] * ip_20[k];

        t_51[k] = f_8 * gd0_13[k]
                  - f_9 * gd1_13[k]
                  + pa_y[k] * hd_25[k];

        t_52[k] = f_1 * is0_12[k]
                  - f_2 * is1_20[k]
                  + pb_x[k] * ip_21[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pa_z, pb_y, gd0_11, gd0_14, gd1_11, gd1_14, \
                         hp_12, hd_24, hd_26, ip_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_6 * gd0_11[k]
                  - f_7 * gd1_11[k]
                  + pa_z[k] * hd_24[k];

        t_54[k] = f_3 * hp_12[k]
                  + pb_y[k] * ip_22[k];

        t_55[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_14[k]
                  + pa_y[k] * hd_26[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pb_x, pb_y, hp_14, hp_15, hd_28, hd_29, \
                         is0_14, is1_22, ip_23, ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_3 * hp_14[k]
                  + pa_y[k] * hd_28[k];

        t_57[k] = f_13 * hp_15[k]
                  + pb_y[k] * ip_23[k];

        t_58[k] = pa_y[k] * hd_29[k];

        t_59[k] = f_1 * is0_14[k]
                  - f_2 * is1_22[k]
                  + pb_x[k] * ip_24[k];
    }

#pragma omp simd aligned(t_60, t_61, pb_y, pb_z, hp_15, is0_14, is1_22, ip_25, \
                         ip_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * is0_14[k]
                  - f_2 * is1_22[k]
                  + pb_y[k] * ip_25[k];

        t_61[k] = f_0 * hp_15[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_22[k]
                  + pb_z[k] * ip_26[k];
    }
}

auto
compute_prim_id_electron_repulsion_27(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_8 = buffer.data(hp + 8);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_16 = buffer.data(is0 + 16);
    const auto *is0_22 = buffer.data(is0 + 22);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_13 = buffer.data(is1 + 13);
    const auto *is1_19 = buffer.data(is1 + 19);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, gd0_0, gd0_4, gd1_0, gd1_4, \
                         hd_0, hd_1, hd_2, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];

        t_6[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_7[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_y, gd0_1, gd0_6, gd0_7, gd1_1, gd1_6, \
                         gd1_7, hd_3, hd_5, hd_6, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_9[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_10[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_7[k]
                  + pa_x[k] * hd_8[k];

        t_11[k] = pa_y[k] * hd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_y, pa_z, gd0_2, gd0_3, gd0_8, gd1_2, \
                         gd1_3, gd1_8, hd_5, hd_7, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_5[k];

        t_13[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_11[k];

        t_14[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, gd0_5, gd0_9, gd0_11, gd1_5, gd1_9, \
                         gd1_11, hd_9, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_12[k];

        t_16[k] = f_3 * gd0_5[k]
                  - f_4 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_17[k] = f_3 * gd0_11[k]
                  - f_4 * gd1_11[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pa_z, gd0_5, gd0_12, gd0_14, \
                         gd1_5, gd1_12, gd1_14, hd_10, hd_14, hd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_14[k];

        t_19[k] = pa_y[k] * hd_10[k];

        t_20[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_10[k];

        t_21[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pa_x, hd_16, hd_18, hd_19, hd_20, \
                         hd_21, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_x[k] * hd_16[k];

        t_23[k] = pa_x[k] * hd_18[k];

        t_24[k] = pa_x[k] * hd_19[k];

        t_25[k] = pa_x[k] * hd_20[k];

        t_26[k] = pa_x[k] * hd_21[k];

        t_27[k] = pa_x[k] * hd_23[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_z, pb_x, pb_y, pb_z, hp_4, hd_16, is0_16, \
                         is1_13, ip_3, ip_4, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * is0_16[k]
                  - f_2 * is1_13[k]
                  + pb_x[k] * ip_3[k];

        t_29[k] = f_0 * hp_4[k]
                  + f_1 * is0_16[k]
                  - f_2 * is1_13[k]
                  + pb_y[k] * ip_4[k];

        t_30[k] = f_1 * is0_16[k]
                  - f_2 * is1_13[k]
                  + pb_z[k] * ip_5[k];

        t_31[k] = pa_z[k] * hd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, gd0_9, gd0_10, gd0_12, gd1_9, gd1_10, \
                         gd1_12, hd_17, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_17[k];

        t_33[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_19[k];

        t_34[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_y, pa_z, gd0_11, gd0_13, gd0_14, gd1_11, \
                         gd1_13, gd1_14, hd_20, hd_21, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_21[k];

        t_36[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_20[k];

        t_37[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_22[k];

        t_38[k] = pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, pb_y, pb_z, hp_8, is0_22, is1_19, ip_6, ip_7, \
                         ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * is0_22[k]
                  - f_2 * is1_19[k]
                  + pb_x[k] * ip_6[k];

        t_40[k] = f_1 * is0_22[k]
                  - f_2 * is1_19[k]
                  + pb_y[k] * ip_7[k];

        t_41[k] = f_0 * hp_8[k]
                  + f_1 * is0_22[k]
                  - f_2 * is1_19[k]
                  + pb_z[k] * ip_8[k];
    }
}

auto
compute_prim_id_electron_repulsion_28(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_14 = buffer.data(is1 + 14);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gd0_0, gd1_0, hp_0, hd_0, hd_1, \
                         is0_0, is1_0, ip_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = pa_y[k] * hd_0[k];

        t_2[k] = pa_z[k] * hd_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, gd0_0, gd0_4, gd0_6, gd1_0, gd1_4, gd1_6, \
                         hd_2, hd_4, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pa_z, gd0_1, gd0_2, gd0_7, gd1_1, gd1_2, \
                         gd1_7, hd_3, hd_5, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_8[k] = f_7 * gd0_7[k]
                 - f_8 * gd1_7[k]
                 + pa_x[k] * hd_8[k];

        t_9[k] = f_7 * gd0_2[k]
                 - f_8 * gd1_2[k]
                 + pa_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_7, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_10[k];

        t_11[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];

        t_12[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_z, gd0_5, gd0_14, gd1_5, gd1_14, \
                         hd_9, hd_12, hd_13, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_9[k];

        t_14[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_12[k];

        t_15[k] = pa_x[k] * hd_13[k];

        t_16[k] = pa_x[k] * hd_20[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, gd0_9, gd1_9, hp_1, hd_13, hd_14, \
                         is0_9, is1_9, ip_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * hp_1[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_1[k];

        t_18[k] = pa_z[k] * hd_13[k];

        t_19[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pa_z, gd0_10, gd0_12, gd0_13, gd1_10, gd1_12, \
                         gd1_13, hd_15, hd_16, hd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_16[k];

        t_21[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_15[k];

        t_22[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_18[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, gd0_11, gd0_14, gd1_11, gd1_14, hd_17, \
                         hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_17[k];

        t_24[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_19[k];

        t_25[k] = pa_y[k] * hd_20[k];
    }

#pragma omp simd aligned(t_26, pb_z, hp_2, is0_14, is1_14, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_2[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_2[k];
    }
}

auto
compute_prim_id_electron_repulsion_29(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_14 = buffer.data(is1 + 14);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gd0_0, gd1_0, hp_0, hd_0, hd_1, \
                         is0_0, is1_0, ip_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = pa_y[k] * hd_0[k];

        t_2[k] = pa_z[k] * hd_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, gd0_0, gd0_4, gd0_6, gd1_0, gd1_4, gd1_6, \
                         hd_2, hd_4, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pa_z, gd0_1, gd0_2, gd0_7, gd1_1, gd1_2, \
                         gd1_7, hd_3, hd_5, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_8[k] = f_7 * gd0_7[k]
                 - f_8 * gd1_7[k]
                 + pa_x[k] * hd_8[k];

        t_9[k] = f_7 * gd0_2[k]
                 - f_8 * gd1_2[k]
                 + pa_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_7, hd_11, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_11[k];

        t_11[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];

        t_12[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_z, gd0_5, gd0_14, gd1_5, gd1_14, \
                         hd_10, hd_15, hd_16, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_10[k];

        t_14[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_15[k];

        t_15[k] = pa_x[k] * hd_16[k];

        t_16[k] = pa_x[k] * hd_23[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, gd0_9, gd1_9, hp_1, hd_16, hd_17, \
                         is0_9, is1_9, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * hp_1[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_18[k] = pa_z[k] * hd_16[k];

        t_19[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pa_z, gd0_10, gd0_12, gd0_13, gd1_10, gd1_12, \
                         gd1_13, hd_18, hd_19, hd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_19[k];

        t_21[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_18[k];

        t_22[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_21[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, gd0_11, gd0_14, gd1_11, gd1_14, hd_20, \
                         hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_20[k];

        t_24[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_22[k];

        t_25[k] = pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_26, pb_z, hp_2, is0_14, is1_14, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_2[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_30(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
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
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, gd0_0, gd1_0, hp_1, hp_2, \
                         hd_0, hd_1, hd_2, hd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_5[k] = pa_y[k] * hd_2[k];

        t_6[k] = pa_z[k] * hd_0[k];

        t_7[k] = pa_z[k] * hd_1[k];

        t_8[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_9[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, pa_z, pb_z, gd0_4, gd1_4, hd_4, \
                         hd_6, hd_8, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * gd0_4[k]
                  - f_7 * gd1_4[k]
                  + pa_x[k] * hd_8[k];

        t_11[k] = f_1 * is0_1[k]
                  - f_2 * is1_1[k]
                  + pb_z[k] * ip_3[k];

        t_12[k] = pa_z[k] * hd_4[k];

        t_13[k] = pa_y[k] * hd_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_y, gd0_0, gd0_6, gd1_0, gd1_6, hd_5, \
                         hd_12, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * gd0_0[k]
                  - f_5 * gd1_0[k]
                  + pa_z[k] * hd_5[k];

        t_15[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_4[k];

        t_16[k] = f_6 * gd0_6[k]
                  - f_7 * gd1_6[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, pb_z, gd0_1, gd0_7, gd1_1, gd1_7, hd_7, \
                         hd_14, is0_3, is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * gd0_1[k]
                  - f_9 * gd1_1[k]
                  + pa_y[k] * hd_7[k];

        t_18[k] = f_8 * gd0_7[k]
                  - f_9 * gd1_7[k]
                  + pa_x[k] * hd_14[k];

        t_19[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pa_y, pa_z, hp_3, hp_4, hd_7, \
                         hd_8, hd_9, hd_10, hd_11, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * hd_7[k];

        t_21[k] = pa_z[k] * hd_8[k];

        t_22[k] = f_3 * hp_3[k]
                  + pa_z[k] * hd_9[k];

        t_23[k] = pa_y[k] * hd_10[k];

        t_24[k] = f_3 * hp_4[k]
                  + pa_y[k] * hd_11[k];

        t_25[k] = pa_y[k] * hd_12[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_2, gd1_8, \
                         hd_10, hd_19, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * gd0_2[k]
                  - f_9 * gd1_2[k]
                  + pa_z[k] * hd_10[k];

        t_27[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_28[k] = f_8 * gd0_8[k]
                  - f_9 * gd1_8[k]
                  + pa_x[k] * hd_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_3, gd1_9, \
                         hd_13, hd_21, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_6 * gd0_3[k]
                  - f_7 * gd1_3[k]
                  + pa_y[k] * hd_13[k];

        t_30[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_9[k]
                  + pa_x[k] * hd_21[k];

        t_31[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_y, pa_z, gd0_5, gd1_5, hp_5, hd_13, hd_14, \
                         hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pa_z[k] * hd_13[k];

        t_33[k] = pa_z[k] * hd_14[k];

        t_34[k] = f_3 * hp_5[k]
                  + pa_z[k] * hd_15[k];

        t_35[k] = f_4 * gd0_5[k]
                  - f_5 * gd1_5[k]
                  + pa_y[k] * hd_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pa_y, gd0_11, gd0_12, gd1_11, gd1_12, \
                         hp_6, hd_17, hd_18, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_4 * gd0_11[k]
                  - f_5 * gd1_11[k]
                  + pa_x[k] * hd_22[k];

        t_37[k] = f_4 * gd0_12[k]
                  - f_5 * gd1_12[k]
                  + pa_x[k] * hd_23[k];

        t_38[k] = pa_y[k] * hd_17[k];

        t_39[k] = f_3 * hp_6[k]
                  + pa_y[k] * hd_18[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pa_z, pb_y, gd0_5, gd1_5, hd_17, hd_19, \
                         is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * hd_19[k];

        t_41[k] = f_6 * gd0_5[k]
                  - f_7 * gd1_5[k]
                  + pa_z[k] * hd_17[k];

        t_42[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_x, pa_z, gd0_14, gd1_14, hp_7, \
                         hp_10, hd_20, hd_25, hd_26, hd_27, hd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_14[k]
                  + pa_x[k] * hd_25[k];

        t_44[k] = f_3 * hp_7[k]
                  + pa_x[k] * hd_26[k];

        t_45[k] = pa_x[k] * hd_27[k];

        t_46[k] = pa_z[k] * hd_20[k];

        t_47[k] = f_3 * hp_10[k]
                  + pa_x[k] * hd_31[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, t_53, pa_x, hp_11, hp_12, hd_32, hd_33, \
                         hd_34, hd_35, hd_36, hd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_x[k] * hd_32[k];

        t_49[k] = pa_x[k] * hd_33[k];

        t_50[k] = f_3 * hp_11[k]
                  + pa_x[k] * hd_34[k];

        t_51[k] = pa_x[k] * hd_35[k];

        t_52[k] = pa_x[k] * hd_36[k];

        t_53[k] = f_3 * hp_12[k]
                  + pa_x[k] * hd_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_x, pb_y, pb_z, hp_8, hd_41, is0_9, \
                         is1_9, ip_9, ip_10, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_x[k] * hd_41[k];

        t_55[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];

        t_56[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_57[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_x, hp_9, hd_26, hd_27, hd_28, \
                         is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pa_z[k] * hd_26[k];

        t_59[k] = pa_z[k] * hd_27[k];

        t_60[k] = f_3 * hp_9[k]
                  + pa_z[k] * hd_28[k];

        t_61[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_9, gd1_12, \
                         hd_29, hd_33, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_9[k]
                  + pa_z[k] * hd_29[k];

        t_63[k] = f_6 * gd0_12[k]
                  - f_7 * gd1_12[k]
                  + pa_y[k] * hd_33[k];

        t_64[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_10, gd1_13, \
                         hd_32, hd_36, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_8 * gd0_10[k]
                  - f_9 * gd1_10[k]
                  + pa_z[k] * hd_32[k];

        t_66[k] = f_8 * gd0_13[k]
                  - f_9 * gd1_13[k]
                  + pa_y[k] * hd_36[k];

        t_67[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pa_z, gd0_11, gd0_14, gd1_11, gd1_14, \
                         hp_13, hd_35, hd_38, hd_40, hd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_6 * gd0_11[k]
                  - f_7 * gd1_11[k]
                  + pa_z[k] * hd_35[k];

        t_69[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_14[k]
                  + pa_y[k] * hd_38[k];

        t_70[k] = f_3 * hp_13[k]
                  + pa_y[k] * hd_40[k];

        t_71[k] = pa_y[k] * hd_41[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, pb_y, pb_z, hp_14, is0_14, is1_14, ip_15, \
                         ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_73[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_74[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_31(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_20 = buffer.data(gd1 + 20);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pa_z, gd0_0, gd1_0, hp_1, hp_2, hd_0, hd_1, \
                         hd_2, hd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * hp_1[k]
                 + pa_y[k] * hd_1[k];

        t_5[k] = pa_z[k] * hd_0[k];

        t_6[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_7[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_z, pb_z, gd0_0, gd0_4, gd1_0, gd1_6, hd_4, \
                         hd_6, is0_1, is1_1, ip_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * gd0_4[k]
                 - f_7 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_9[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_10[k] = f_4 * gd0_0[k]
                  - f_5 * gd1_0[k]
                  + pa_z[k] * hd_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, pb_y, gd0_1, gd0_6, gd1_3, gd1_8, hd_5, \
                         hd_8, is0_2, is1_2, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_4[k];

        t_12[k] = f_6 * gd0_6[k]
                  - f_7 * gd1_8[k]
                  + pa_x[k] * hd_8[k];

        t_13[k] = f_8 * gd0_1[k]
                  - f_9 * gd1_3[k]
                  + pa_y[k] * hd_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_z, gd0_7, gd1_9, hd_7, hd_10, is0_3, \
                         is1_3, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * gd0_7[k]
                  - f_9 * gd1_9[k]
                  + pa_x[k] * hd_10[k];

        t_15[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_16[k] = pa_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_z, pb_y, gd0_2, gd0_8, gd1_4, gd1_10, \
                         hd_7, hd_13, is0_4, is1_4, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * gd0_2[k]
                  - f_9 * gd1_4[k]
                  + pa_z[k] * hd_7[k];

        t_18[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];

        t_19[k] = f_8 * gd0_8[k]
                  - f_9 * gd1_10[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pa_y, pb_z, gd0_3, gd0_9, gd1_5, gd1_12, \
                         hd_9, hd_14, is0_5, is1_5, ip_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * gd0_3[k]
                  - f_7 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_21[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_12[k]
                  + pa_x[k] * hd_14[k];

        t_22[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, gd0_5, gd0_11, gd0_12, gd1_7, \
                         gd1_15, gd1_16, hd_11, hd_12, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * gd0_5[k]
                  - f_5 * gd1_7[k]
                  + pa_y[k] * hd_11[k];

        t_24[k] = f_4 * gd0_11[k]
                  - f_5 * gd1_15[k]
                  + pa_x[k] * hd_15[k];

        t_25[k] = f_4 * gd0_12[k]
                  - f_5 * gd1_16[k]
                  + pa_x[k] * hd_16[k];

        t_26[k] = pa_y[k] * hd_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_z, pb_y, gd0_5, gd0_14, gd1_7, gd1_20, \
                         hd_12, hd_17, is0_6, is1_6, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_6 * gd0_5[k]
                  - f_7 * gd1_7[k]
                  + pa_z[k] * hd_12[k];

        t_28[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];

        t_29[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_20[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pa_x, hp_7, hd_18, hd_19, hd_22, \
                         hd_23, hd_24, hd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * hp_7[k]
                  + pa_x[k] * hd_18[k];

        t_31[k] = pa_x[k] * hd_19[k];

        t_32[k] = pa_x[k] * hd_22[k];

        t_33[k] = pa_x[k] * hd_23[k];

        t_34[k] = pa_x[k] * hd_24[k];

        t_35[k] = pa_x[k] * hd_25[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_x, pb_y, hp_8, hp_12, hd_27, hd_29, \
                         is0_9, is1_9, ip_9, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_3 * hp_12[k]
                  + pa_x[k] * hd_27[k];

        t_37[k] = pa_x[k] * hd_29[k];

        t_38[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];

        t_39[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_x, pb_z, hp_9, hd_19, hd_20, is0_9, \
                         is0_10, is1_9, is1_10, ip_11, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_41[k] = pa_z[k] * hd_19[k];

        t_42[k] = f_3 * hp_9[k]
                  + pa_z[k] * hd_20[k];

        t_43[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_12, gd1_16, \
                         hd_21, hd_23, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_12[k]
                  + pa_z[k] * hd_21[k];

        t_45[k] = f_6 * gd0_12[k]
                  - f_7 * gd1_16[k]
                  + pa_y[k] * hd_23[k];

        t_46[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_14, gd1_17, \
                         hd_22, hd_25, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_8 * gd0_10[k]
                  - f_9 * gd1_14[k]
                  + pa_z[k] * hd_22[k];

        t_48[k] = f_8 * gd0_13[k]
                  - f_9 * gd1_17[k]
                  + pa_y[k] * hd_25[k];

        t_49[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pa_z, gd0_11, gd0_14, gd1_15, gd1_20, \
                         hp_13, hd_24, hd_26, hd_28, hd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * gd0_11[k]
                  - f_7 * gd1_15[k]
                  + pa_z[k] * hd_24[k];

        t_51[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_20[k]
                  + pa_y[k] * hd_26[k];

        t_52[k] = f_3 * hp_13[k]
                  + pa_y[k] * hd_28[k];

        t_53[k] = pa_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, pb_y, pb_z, hp_14, is0_14, is1_14, ip_15, \
                         ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_55[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_56[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_32(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_20 = buffer.data(gd0 + 20);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_8 = buffer.data(hp + 8);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_14 = buffer.data(is1 + 14);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, gd0_0, gd0_6, gd1_0, gd1_4, \
                         hd_0, hd_1, hd_2, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_7[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_1, gd1_6, \
                         gd1_7, hd_3, hd_5, hd_6, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * gd0_8[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_9[k] = f_7 * gd0_3[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_10[k] = f_7 * gd0_9[k]
                  - f_8 * gd1_7[k]
                  + pa_x[k] * hd_8[k];

        t_11[k] = pa_y[k] * hd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_y, pa_z, gd0_4, gd0_5, gd0_10, gd1_2, \
                         gd1_3, gd1_8, hd_5, hd_7, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * gd0_4[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_5[k];

        t_13[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_11[k];

        t_14[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_y, gd0_7, gd0_12, gd0_15, gd1_5, gd1_9, \
                         gd1_11, hd_9, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_12[k];

        t_16[k] = f_3 * gd0_7[k]
                  - f_4 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_17[k] = f_3 * gd0_15[k]
                  - f_4 * gd1_11[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pa_z, gd0_7, gd0_16, gd0_20, \
                         gd1_5, gd1_12, gd1_14, hd_10, hd_14, hd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * gd0_16[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_14[k];

        t_19[k] = pa_y[k] * hd_10[k];

        t_20[k] = f_5 * gd0_7[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_10[k];

        t_21[k] = f_3 * gd0_20[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pa_x, hd_16, hd_18, hd_19, hd_20, \
                         hd_21, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_x[k] * hd_16[k];

        t_23[k] = pa_x[k] * hd_18[k];

        t_24[k] = pa_x[k] * hd_19[k];

        t_25[k] = pa_x[k] * hd_20[k];

        t_26[k] = pa_x[k] * hd_21[k];

        t_27[k] = pa_x[k] * hd_23[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_z, pb_x, pb_y, pb_z, hp_4, hd_16, is0_9, \
                         is1_9, ip_3, ip_4, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_3[k];

        t_29[k] = f_0 * hp_4[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_4[k];

        t_30[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_5[k];

        t_31[k] = pa_z[k] * hd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, gd0_12, gd0_14, gd0_16, gd1_9, gd1_10, \
                         gd1_12, hd_17, hd_18, hd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_17[k];

        t_33[k] = f_5 * gd0_16[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_19[k];

        t_34[k] = f_7 * gd0_14[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_y, pa_z, gd0_15, gd0_17, gd0_20, gd1_11, \
                         gd1_13, gd1_14, hd_20, hd_21, hd_22, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_7 * gd0_17[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_21[k];

        t_36[k] = f_5 * gd0_15[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_20[k];

        t_37[k] = f_3 * gd0_20[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_22[k];

        t_38[k] = pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, pb_y, pb_z, hp_8, is0_14, is1_14, ip_6, ip_7, \
                         ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_6[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_7[k];

        t_41[k] = f_0 * hp_8[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_8[k];
    }
}

auto
compute_prim_id_electron_repulsion_33(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_14 = buffer.data(is1 + 14);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_17 = buffer.data(ip + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pa_z, pb_x, gd0_0, gd1_0, hp_0, hd_0, hd_1, \
                         is0_0, is1_0, ip_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = pa_y[k] * hd_0[k];

        t_2[k] = pa_z[k] * hd_0[k];

        t_3[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_z, gd0_0, gd0_4, gd0_6, gd1_0, gd1_4, gd1_6, \
                         hd_2, hd_4, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pa_z, gd0_1, gd0_2, gd0_7, gd1_1, gd1_2, \
                         gd1_7, hd_3, hd_5, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_7 * gd0_1[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_8[k] = f_7 * gd0_7[k]
                 - f_8 * gd1_7[k]
                 + pa_x[k] * hd_8[k];

        t_9[k] = f_7 * gd0_2[k]
                 - f_8 * gd1_2[k]
                 + pa_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_7, hd_10, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_10[k];

        t_11[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];

        t_12[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_11[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_z, gd0_5, gd0_14, gd1_5, gd1_14, \
                         hd_9, hd_12, hd_13, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_9[k];

        t_14[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_12[k];

        t_15[k] = pa_x[k] * hd_13[k];

        t_16[k] = pa_x[k] * hd_20[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_z, pb_y, gd0_9, gd1_9, hp_1, hd_13, hd_14, \
                         is0_9, is1_9, ip_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * hp_1[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_18[k] = pa_z[k] * hd_13[k];

        t_19[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pa_z, gd0_10, gd0_12, gd0_13, gd1_10, gd1_12, \
                         gd1_13, hd_15, hd_16, hd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_16[k];

        t_21[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_15[k];

        t_22[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_18[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, gd0_11, gd0_14, gd1_11, gd1_14, hd_17, \
                         hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_17[k];

        t_24[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_19[k];

        t_25[k] = pa_y[k] * hd_20[k];
    }

#pragma omp simd aligned(t_26, pb_z, hp_2, is0_14, is1_14, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_2[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_34(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 1.5 / alpha;
    const auto f_7 = 1.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_20 = buffer.data(gd1 + 20);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, gd0_0, gd0_4, gd1_0, gd1_6, \
                         hp_2, hd_0, hd_2, hd_3, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * hp_2[k]
                 + pa_z[k] * hd_2[k];

        t_6[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_7[k] = f_6 * gd0_4[k]
                 - f_7 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_9[k] = f_4 * gd0_0[k]
                 - f_5 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_10[k] = f_1 * is0_2[k]
                  - f_2 * is1_2[k]
                  + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, gd0_1, gd0_6, gd0_7, gd1_3, gd1_8, \
                         gd1_9, hd_5, hd_10, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_6 * gd0_6[k]
                  - f_7 * gd1_8[k]
                  + pa_x[k] * hd_10[k];

        t_12[k] = f_8 * gd0_1[k]
                  - f_9 * gd1_3[k]
                  + pa_y[k] * hd_5[k];

        t_13[k] = f_8 * gd0_7[k]
                  - f_9 * gd1_9[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_y, pb_z, gd0_2, gd1_4, hd_8, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_15[k] = f_8 * gd0_2[k]
                  - f_9 * gd1_4[k]
                  + pa_z[k] * hd_8[k];

        t_16[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_5, gd1_10, \
                         gd1_12, hd_11, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * gd0_8[k]
                  - f_9 * gd1_10[k]
                  + pa_x[k] * hd_16[k];

        t_18[k] = f_6 * gd0_3[k]
                  - f_7 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_19[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_12[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_z, pb_y, pb_z, gd0_5, gd1_7, hd_14, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_21[k] = f_6 * gd0_5[k]
                  - f_7 * gd1_7[k]
                  + pa_z[k] * hd_14[k];

        t_22[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, gd0_14, gd1_20, hd_18, hd_20, \
                         hd_32, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_20[k]
                  + pa_x[k] * hd_18[k];

        t_24[k] = pa_x[k] * hd_20[k];

        t_25[k] = pa_x[k] * hd_32[k];

        t_26[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_z, pb_y, pb_z, hp_8, hp_9, hd_20, hd_21, \
                         is0_9, is1_9, ip_10, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_28[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_29[k] = pa_z[k] * hd_20[k];

        t_30[k] = f_3 * hp_9[k]
                  + pa_z[k] * hd_21[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_12, gd1_16, \
                         hd_22, hd_25, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_32[k] = f_4 * gd0_9[k]
                  - f_5 * gd1_12[k]
                  + pa_z[k] * hd_22[k];

        t_33[k] = f_6 * gd0_12[k]
                  - f_7 * gd1_16[k]
                  + pa_y[k] * hd_25[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_14, gd1_17, \
                         hd_24, hd_28, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_35[k] = f_8 * gd0_10[k]
                  - f_9 * gd1_14[k]
                  + pa_z[k] * hd_24[k];

        t_36[k] = f_8 * gd0_13[k]
                  - f_9 * gd1_17[k]
                  + pa_y[k] * hd_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_15, gd1_20, \
                         hd_27, hd_29, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_38[k] = f_6 * gd0_11[k]
                  - f_7 * gd1_15[k]
                  + pa_z[k] * hd_27[k];

        t_39[k] = f_4 * gd0_14[k]
                  - f_5 * gd1_20[k]
                  + pa_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pb_x, pb_y, hp_13, hd_31, hd_32, \
                         is0_14, is1_14, ip_15, ip_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_3 * hp_13[k]
                  + pa_y[k] * hd_31[k];

        t_41[k] = pa_y[k] * hd_32[k];

        t_42[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_43[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];
    }

#pragma omp simd aligned(t_44, pb_z, hp_14, is0_14, is1_14, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_35(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_20 = buffer.data(gd0 + 20);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_20 = buffer.data(gd1 + 20);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
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
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_26 = buffer.data(hd + 26);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, gd0_0, gd0_6, gd1_0, gd1_6, hd_0, \
                         hd_3, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_8[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_9[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_5, hd_8, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd0_8[k]
                  - f_6 * gd1_8[k]
                  + pa_x[k] * hd_8[k];

        t_11[k] = f_7 * gd0_3[k]
                  - f_8 * gd1_3[k]
                  + pa_y[k] * hd_5[k];

        t_12[k] = f_7 * gd0_9[k]
                  - f_8 * gd1_9[k]
                  + pa_x[k] * hd_10[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, gd0_4, gd1_4, hd_7, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_14[k] = f_7 * gd0_4[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_7[k];

        t_15[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, gd0_5, gd0_10, gd0_12, gd1_5, gd1_10, \
                         gd1_12, hd_9, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_x[k] * hd_12[k];

        t_17[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_9[k];

        t_18[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, gd0_7, gd1_7, hd_11, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_20[k] = f_5 * gd0_7[k]
                  - f_6 * gd1_7[k]
                  + pa_z[k] * hd_11[k];

        t_21[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, gd0_20, gd1_20, hd_14, hd_16, \
                         hd_26, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gd0_20[k]
                  - f_4 * gd1_20[k]
                  + pa_x[k] * hd_14[k];

        t_23[k] = pa_x[k] * hd_16[k];

        t_24[k] = pa_x[k] * hd_26[k];

        t_25[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hp_8, hd_16, is0_9, is1_9, ip_10, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_16[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, gd0_12, gd0_16, gd1_12, gd1_16, \
                         hd_18, hd_20, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_30[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_12[k]
                  + pa_z[k] * hd_18[k];

        t_31[k] = f_5 * gd0_16[k]
                  - f_6 * gd1_16[k]
                  + pa_y[k] * hd_20[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_x, gd0_14, gd0_17, gd1_14, gd1_17, \
                         hd_19, hd_22, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_7 * gd0_14[k]
                  - f_8 * gd1_14[k]
                  + pa_z[k] * hd_19[k];

        t_34[k] = f_7 * gd0_17[k]
                  - f_8 * gd1_17[k]
                  + pa_y[k] * hd_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, gd0_15, gd0_20, gd1_15, gd1_20, \
                         hd_21, hd_23, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_36[k] = f_5 * gd0_15[k]
                  - f_6 * gd1_15[k]
                  + pa_z[k] * hd_21[k];

        t_37[k] = f_3 * gd0_20[k]
                  - f_4 * gd1_20[k]
                  + pa_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, hp_14, hd_26, is0_14, \
                         is1_14, ip_15, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hd_26[k];

        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_41[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_36(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_20 = buffer.data(gd0 + 20);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_1 = buffer.data(gd1 + 1);
    const auto *gd1_2 = buffer.data(gd1 + 2);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_11 = buffer.data(gd1 + 11);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_13 = buffer.data(gd1 + 13);
    const auto *gd1_14 = buffer.data(gd1 + 14);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_8 = buffer.data(hp + 8);

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

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_14 = buffer.data(is1 + 14);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, gd0_0, gd0_6, gd1_0, gd1_4, \
                         hd_0, hd_1, hd_2, hd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_1[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_4[k]
                 + pa_x[k] * hd_4[k];

        t_7[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_2[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_1, gd1_6, gd1_7, \
                         hd_3, hd_6, hd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * gd0_8[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];

        t_9[k] = f_7 * gd0_3[k]
                 - f_8 * gd1_1[k]
                 + pa_y[k] * hd_3[k];

        t_10[k] = f_7 * gd0_9[k]
                  - f_8 * gd1_7[k]
                  + pa_x[k] * hd_8[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_x, pa_y, pa_z, gd0_4, gd0_5, gd0_10, gd1_2, \
                         gd1_3, gd1_8, hd_5, hd_7, hd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_7 * gd0_4[k]
                  - f_8 * gd1_2[k]
                  + pa_z[k] * hd_5[k];

        t_12[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_8[k]
                  + pa_x[k] * hd_10[k];

        t_13[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_3[k]
                  + pa_y[k] * hd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, gd0_7, gd0_12, gd0_20, gd1_5, \
                         gd1_9, gd1_14, hd_9, hd_11, hd_12, hd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_9[k]
                  + pa_x[k] * hd_11[k];

        t_15[k] = f_5 * gd0_7[k]
                  - f_6 * gd1_5[k]
                  + pa_z[k] * hd_9[k];

        t_16[k] = f_3 * gd0_20[k]
                  - f_4 * gd1_14[k]
                  + pa_x[k] * hd_12[k];

        t_17[k] = pa_x[k] * hd_13[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pb_x, pb_y, pb_z, hp_4, hd_20, is0_9, \
                         is1_9, ip_3, ip_4, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pa_x[k] * hd_20[k];

        t_19[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_3[k];

        t_20[k] = f_0 * hp_4[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_4[k];

        t_21[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_5[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pa_z, gd0_12, gd0_14, gd0_16, gd1_9, \
                         gd1_10, gd1_12, hd_13, hd_14, hd_15, hd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_z[k] * hd_13[k];

        t_23[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_9[k]
                  + pa_z[k] * hd_14[k];

        t_24[k] = f_5 * gd0_16[k]
                  - f_6 * gd1_12[k]
                  + pa_y[k] * hd_16[k];

        t_25[k] = f_7 * gd0_14[k]
                  - f_8 * gd1_10[k]
                  + pa_z[k] * hd_15[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, gd0_15, gd0_17, gd0_20, gd1_11, \
                         gd1_13, gd1_14, hd_17, hd_18, hd_19, hd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * gd0_17[k]
                  - f_8 * gd1_13[k]
                  + pa_y[k] * hd_18[k];

        t_27[k] = f_5 * gd0_15[k]
                  - f_6 * gd1_11[k]
                  + pa_z[k] * hd_17[k];

        t_28[k] = f_3 * gd0_20[k]
                  - f_4 * gd1_14[k]
                  + pa_y[k] * hd_19[k];

        t_29[k] = pa_y[k] * hd_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, pb_y, pb_z, hp_8, is0_14, is1_14, ip_6, ip_7, \
                         ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_6[k];

        t_31[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_7[k];

        t_32[k] = f_0 * hp_8[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_8[k];
    }
}

auto
compute_prim_id_electron_repulsion_37(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_1 = buffer.data(gd0 + 1);
    const auto *gd0_2 = buffer.data(gd0 + 2);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_11 = buffer.data(gd0 + 11);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_13 = buffer.data(gd0 + 13);
    const auto *gd0_14 = buffer.data(gd0 + 14);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_20 = buffer.data(gd1 + 20);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, gd0_0, gd0_4, gd1_0, gd1_6, hd_0, \
                         hd_3, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_6[k] = f_5 * gd0_4[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_8[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_9[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_1, gd0_6, gd0_7, gd1_3, gd1_8, \
                         gd1_9, hd_5, hd_10, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd0_6[k]
                  - f_6 * gd1_8[k]
                  + pa_x[k] * hd_10[k];

        t_11[k] = f_7 * gd0_1[k]
                  - f_8 * gd1_3[k]
                  + pa_y[k] * hd_5[k];

        t_12[k] = f_7 * gd0_7[k]
                  - f_8 * gd1_9[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, gd0_2, gd1_4, hd_8, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_14[k] = f_7 * gd0_2[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_8[k];

        t_15[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_5, gd1_10, \
                         gd1_12, hd_11, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gd0_8[k]
                  - f_8 * gd1_10[k]
                  + pa_x[k] * hd_16[k];

        t_17[k] = f_5 * gd0_3[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_18[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, gd0_5, gd1_7, hd_14, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_20[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_7[k]
                  + pa_z[k] * hd_14[k];

        t_21[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, gd0_14, gd1_20, hd_18, hd_20, \
                         hd_32, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_20[k]
                  + pa_x[k] * hd_18[k];

        t_23[k] = pa_x[k] * hd_20[k];

        t_24[k] = pa_x[k] * hd_32[k];

        t_25[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hp_8, hd_20, is0_9, is1_9, ip_10, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, gd0_9, gd0_12, gd1_12, gd1_16, \
                         hd_22, hd_25, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_30[k] = f_3 * gd0_9[k]
                  - f_4 * gd1_12[k]
                  + pa_z[k] * hd_22[k];

        t_31[k] = f_5 * gd0_12[k]
                  - f_6 * gd1_16[k]
                  + pa_y[k] * hd_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_x, gd0_10, gd0_13, gd1_14, gd1_17, \
                         hd_24, hd_28, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_14[k]
                  + pa_z[k] * hd_24[k];

        t_34[k] = f_7 * gd0_13[k]
                  - f_8 * gd1_17[k]
                  + pa_y[k] * hd_28[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, gd0_11, gd0_14, gd1_15, gd1_20, \
                         hd_27, hd_29, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_36[k] = f_5 * gd0_11[k]
                  - f_6 * gd1_15[k]
                  + pa_z[k] * hd_27[k];

        t_37[k] = f_3 * gd0_14[k]
                  - f_4 * gd1_20[k]
                  + pa_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, hp_14, hd_32, is0_14, \
                         is1_14, ip_15, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hd_32[k];

        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_41[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

auto
compute_prim_id_electron_repulsion_38(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t gd0, const size_t gd1,
                                      const size_t hp, const size_t hd, const size_t is0,
                                      const size_t is1, const size_t ip, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 1.5 / alpha;
    const auto f_6 = 1.5 * beta / (alpha * p);
    const auto f_7 = 1.0 / alpha;
    const auto f_8 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gd0_0 = buffer.data(gd0 + 0);
    const auto *gd0_3 = buffer.data(gd0 + 3);
    const auto *gd0_4 = buffer.data(gd0 + 4);
    const auto *gd0_5 = buffer.data(gd0 + 5);
    const auto *gd0_6 = buffer.data(gd0 + 6);
    const auto *gd0_7 = buffer.data(gd0 + 7);
    const auto *gd0_8 = buffer.data(gd0 + 8);
    const auto *gd0_9 = buffer.data(gd0 + 9);
    const auto *gd0_10 = buffer.data(gd0 + 10);
    const auto *gd0_12 = buffer.data(gd0 + 12);
    const auto *gd0_14 = buffer.data(gd0 + 14);
    const auto *gd0_15 = buffer.data(gd0 + 15);
    const auto *gd0_16 = buffer.data(gd0 + 16);
    const auto *gd0_17 = buffer.data(gd0 + 17);
    const auto *gd0_20 = buffer.data(gd0 + 20);

    const auto *gd1_0 = buffer.data(gd1 + 0);
    const auto *gd1_3 = buffer.data(gd1 + 3);
    const auto *gd1_4 = buffer.data(gd1 + 4);
    const auto *gd1_5 = buffer.data(gd1 + 5);
    const auto *gd1_6 = buffer.data(gd1 + 6);
    const auto *gd1_7 = buffer.data(gd1 + 7);
    const auto *gd1_8 = buffer.data(gd1 + 8);
    const auto *gd1_9 = buffer.data(gd1 + 9);
    const auto *gd1_10 = buffer.data(gd1 + 10);
    const auto *gd1_12 = buffer.data(gd1 + 12);
    const auto *gd1_14 = buffer.data(gd1 + 14);
    const auto *gd1_15 = buffer.data(gd1 + 15);
    const auto *gd1_16 = buffer.data(gd1 + 16);
    const auto *gd1_17 = buffer.data(gd1 + 17);
    const auto *gd1_20 = buffer.data(gd1 + 20);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_14 = buffer.data(hp + 14);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_32 = buffer.data(hd + 32);

    const auto *is0_0 = buffer.data(is0 + 0);
    const auto *is0_1 = buffer.data(is0 + 1);
    const auto *is0_2 = buffer.data(is0 + 2);
    const auto *is0_3 = buffer.data(is0 + 3);
    const auto *is0_4 = buffer.data(is0 + 4);
    const auto *is0_5 = buffer.data(is0 + 5);
    const auto *is0_6 = buffer.data(is0 + 6);
    const auto *is0_9 = buffer.data(is0 + 9);
    const auto *is0_10 = buffer.data(is0 + 10);
    const auto *is0_11 = buffer.data(is0 + 11);
    const auto *is0_12 = buffer.data(is0 + 12);
    const auto *is0_14 = buffer.data(is0 + 14);

    const auto *is1_0 = buffer.data(is1 + 0);
    const auto *is1_1 = buffer.data(is1 + 1);
    const auto *is1_2 = buffer.data(is1 + 2);
    const auto *is1_3 = buffer.data(is1 + 3);
    const auto *is1_4 = buffer.data(is1 + 4);
    const auto *is1_5 = buffer.data(is1 + 5);
    const auto *is1_6 = buffer.data(is1 + 6);
    const auto *is1_9 = buffer.data(is1 + 9);
    const auto *is1_10 = buffer.data(is1 + 10);
    const auto *is1_11 = buffer.data(is1 + 11);
    const auto *is1_12 = buffer.data(is1 + 12);
    const auto *is1_14 = buffer.data(is1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_x, pb_y, pb_z, hp_0, hd_0, is0_0, is1_0, \
                         ip_0, ip_1, ip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hp_0[k]
                 + f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_x[k] * ip_0[k];

        t_1[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_y[k] * ip_1[k];

        t_2[k] = f_1 * is0_0[k]
                 - f_2 * is1_0[k]
                 + pb_z[k] * ip_2[k];

        t_3[k] = pa_y[k] * hd_0[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pa_y, pa_z, gd0_0, gd0_6, gd1_0, gd1_6, hd_0, \
                         hd_3, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = pa_z[k] * hd_0[k];

        t_5[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_y[k] * hd_3[k];

        t_6[k] = f_5 * gd0_6[k]
                 - f_6 * gd1_6[k]
                 + pa_x[k] * hd_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_z, pb_y, pb_z, gd0_0, gd1_0, hd_4, is0_1, is0_2, \
                         is1_1, is1_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * is0_1[k]
                 - f_2 * is1_1[k]
                 + pb_z[k] * ip_3[k];

        t_8[k] = f_3 * gd0_0[k]
                 - f_4 * gd1_0[k]
                 + pa_z[k] * hd_4[k];

        t_9[k] = f_1 * is0_2[k]
                 - f_2 * is1_2[k]
                 + pb_y[k] * ip_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, gd0_3, gd0_8, gd0_9, gd1_3, gd1_8, \
                         gd1_9, hd_5, hd_10, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * gd0_8[k]
                  - f_6 * gd1_8[k]
                  + pa_x[k] * hd_10[k];

        t_11[k] = f_7 * gd0_3[k]
                  - f_8 * gd1_3[k]
                  + pa_y[k] * hd_5[k];

        t_12[k] = f_7 * gd0_9[k]
                  - f_8 * gd1_9[k]
                  + pa_x[k] * hd_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_z, pb_y, pb_z, gd0_4, gd1_4, hd_8, is0_3, is0_4, \
                         is1_3, is1_4, ip_5, ip_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * is0_3[k]
                  - f_2 * is1_3[k]
                  + pb_z[k] * ip_5[k];

        t_14[k] = f_7 * gd0_4[k]
                  - f_8 * gd1_4[k]
                  + pa_z[k] * hd_8[k];

        t_15[k] = f_1 * is0_4[k]
                  - f_2 * is1_4[k]
                  + pb_y[k] * ip_6[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, gd0_5, gd0_10, gd0_12, gd1_5, gd1_10, \
                         gd1_12, hd_11, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * gd0_10[k]
                  - f_8 * gd1_10[k]
                  + pa_x[k] * hd_16[k];

        t_17[k] = f_5 * gd0_5[k]
                  - f_6 * gd1_5[k]
                  + pa_y[k] * hd_11[k];

        t_18[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_12[k]
                  + pa_x[k] * hd_17[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_y, pb_z, gd0_7, gd1_7, hd_14, is0_5, \
                         is0_6, is1_5, is1_6, ip_7, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * is0_5[k]
                  - f_2 * is1_5[k]
                  + pb_z[k] * ip_7[k];

        t_20[k] = f_5 * gd0_7[k]
                  - f_6 * gd1_7[k]
                  + pa_z[k] * hd_14[k];

        t_21[k] = f_1 * is0_6[k]
                  - f_2 * is1_6[k]
                  + pb_y[k] * ip_8[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, gd0_20, gd1_20, hd_18, hd_20, \
                         hd_32, is0_9, is1_9, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * gd0_20[k]
                  - f_4 * gd1_20[k]
                  + pa_x[k] * hd_18[k];

        t_23[k] = pa_x[k] * hd_20[k];

        t_24[k] = pa_x[k] * hd_32[k];

        t_25[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_x[k] * ip_9[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, pb_z, hp_8, hd_20, is0_9, is1_9, ip_10, \
                         ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * hp_8[k]
                  + f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_y[k] * ip_10[k];

        t_27[k] = f_1 * is0_9[k]
                  - f_2 * is1_9[k]
                  + pb_z[k] * ip_11[k];

        t_28[k] = pa_z[k] * hd_20[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_x, gd0_12, gd0_16, gd1_12, gd1_16, \
                         hd_22, hd_25, is0_10, is1_10, ip_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_1 * is0_10[k]
                  - f_2 * is1_10[k]
                  + pb_x[k] * ip_12[k];

        t_30[k] = f_3 * gd0_12[k]
                  - f_4 * gd1_12[k]
                  + pa_z[k] * hd_22[k];

        t_31[k] = f_5 * gd0_16[k]
                  - f_6 * gd1_16[k]
                  + pa_y[k] * hd_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_x, gd0_14, gd0_17, gd1_14, gd1_17, \
                         hd_24, hd_28, is0_11, is1_11, ip_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * is0_11[k]
                  - f_2 * is1_11[k]
                  + pb_x[k] * ip_13[k];

        t_33[k] = f_7 * gd0_14[k]
                  - f_8 * gd1_14[k]
                  + pa_z[k] * hd_24[k];

        t_34[k] = f_7 * gd0_17[k]
                  - f_8 * gd1_17[k]
                  + pa_y[k] * hd_28[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pa_z, pb_x, gd0_15, gd0_20, gd1_15, gd1_20, \
                         hd_27, hd_29, is0_12, is1_12, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * is0_12[k]
                  - f_2 * is1_12[k]
                  + pb_x[k] * ip_14[k];

        t_36[k] = f_5 * gd0_15[k]
                  - f_6 * gd1_15[k]
                  + pa_z[k] * hd_27[k];

        t_37[k] = f_3 * gd0_20[k]
                  - f_4 * gd1_20[k]
                  + pa_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_x, pb_y, pb_z, hp_14, hd_32, is0_14, \
                         is1_14, ip_15, ip_16, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hd_32[k];

        t_39[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_x[k] * ip_15[k];

        t_40[k] = f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_y[k] * ip_16[k];

        t_41[k] = f_0 * hp_14[k]
                  + f_1 * is0_14[k]
                  - f_2 * is1_14[k]
                  + pb_z[k] * ip_17[k];
    }
}

}  // namespace simdt2ceri
