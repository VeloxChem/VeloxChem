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


#include "SimdKineticEnergyVrrRecFG.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_fg_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * alpha / p;
    const auto f_8 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_s_2 = buffer.data(pg_s + 2);
    const auto *pg_s_5 = buffer.data(pg_s + 5);

    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_5 = buffer.data(pg + 5);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);
    const auto *df_28 = buffer.data(df + 28);
    const auto *df_29 = buffer.data(df + 29);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_31 = buffer.data(df + 31);
    const auto *df_32 = buffer.data(df + 32);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_25 = buffer.data(fd_s + 25);
    const auto *fd_s_26 = buffer.data(fd_s + 26);
    const auto *fd_s_27 = buffer.data(fd_s + 27);
    const auto *fd_s_28 = buffer.data(fd_s + 28);
    const auto *fd_s_29 = buffer.data(fd_s + 29);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_23 = buffer.data(fg_s + 23);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_29 = buffer.data(fg_s + 29);
    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_31 = buffer.data(fg_s + 31);
    const auto *fg_s_32 = buffer.data(fg_s + 32);
    const auto *fg_s_33 = buffer.data(fg_s + 33);
    const auto *fg_s_34 = buffer.data(fg_s + 34);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_36 = buffer.data(fg_s + 36);
    const auto *fg_s_37 = buffer.data(fg_s + 37);
    const auto *fg_s_38 = buffer.data(fg_s + 38);
    const auto *fg_s_39 = buffer.data(fg_s + 39);
    const auto *fg_s_40 = buffer.data(fg_s + 40);
    const auto *fg_s_41 = buffer.data(fg_s + 41);
    const auto *fg_s_42 = buffer.data(fg_s + 42);
    const auto *fg_s_43 = buffer.data(fg_s + 43);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_45 = buffer.data(fg_s + 45);
    const auto *fg_s_46 = buffer.data(fg_s + 46);
    const auto *fg_s_47 = buffer.data(fg_s + 47);
    const auto *fg_s_48 = buffer.data(fg_s + 48);
    const auto *fg_s_49 = buffer.data(fg_s + 49);
    const auto *fg_s_50 = buffer.data(fg_s + 50);
    const auto *fg_s_51 = buffer.data(fg_s + 51);
    const auto *fg_s_52 = buffer.data(fg_s + 52);
    const auto *fg_s_53 = buffer.data(fg_s + 53);
    const auto *fg_s_54 = buffer.data(fg_s + 54);
    const auto *fg_s_55 = buffer.data(fg_s + 55);
    const auto *fg_s_56 = buffer.data(fg_s + 56);
    const auto *fg_s_57 = buffer.data(fg_s + 57);
    const auto *fg_s_58 = buffer.data(fg_s + 58);
    const auto *fg_s_59 = buffer.data(fg_s + 59);
    const auto *fg_s_60 = buffer.data(fg_s + 60);
    const auto *fg_s_61 = buffer.data(fg_s + 61);
    const auto *fg_s_62 = buffer.data(fg_s + 62);
    const auto *fg_s_63 = buffer.data(fg_s + 63);
    const auto *fg_s_64 = buffer.data(fg_s + 64);
    const auto *fg_s_65 = buffer.data(fg_s + 65);
    const auto *fg_s_66 = buffer.data(fg_s + 66);
    const auto *fg_s_67 = buffer.data(fg_s + 67);
    const auto *fg_s_68 = buffer.data(fg_s + 68);
    const auto *fg_s_69 = buffer.data(fg_s + 69);
    const auto *fg_s_70 = buffer.data(fg_s + 70);
    const auto *fg_s_71 = buffer.data(fg_s + 71);
    const auto *fg_s_72 = buffer.data(fg_s + 72);
    const auto *fg_s_73 = buffer.data(fg_s + 73);
    const auto *fg_s_74 = buffer.data(fg_s + 74);
    const auto *fg_s_75 = buffer.data(fg_s + 75);
    const auto *fg_s_76 = buffer.data(fg_s + 76);
    const auto *fg_s_77 = buffer.data(fg_s + 77);
    const auto *fg_s_78 = buffer.data(fg_s + 78);
    const auto *fg_s_79 = buffer.data(fg_s + 79);
    const auto *fg_s_80 = buffer.data(fg_s + 80);
    const auto *fg_s_81 = buffer.data(fg_s + 81);
    const auto *fg_s_82 = buffer.data(fg_s + 82);
    const auto *fg_s_83 = buffer.data(fg_s + 83);
    const auto *fg_s_84 = buffer.data(fg_s + 84);
    const auto *fg_s_85 = buffer.data(fg_s + 85);
    const auto *fg_s_86 = buffer.data(fg_s + 86);
    const auto *fg_s_87 = buffer.data(fg_s + 87);
    const auto *fg_s_88 = buffer.data(fg_s + 88);
    const auto *fg_s_89 = buffer.data(fg_s + 89);
    const auto *fg_s_90 = buffer.data(fg_s + 90);
    const auto *fg_s_91 = buffer.data(fg_s + 91);
    const auto *fg_s_92 = buffer.data(fg_s + 92);
    const auto *fg_s_93 = buffer.data(fg_s + 93);
    const auto *fg_s_94 = buffer.data(fg_s + 94);
    const auto *fg_s_95 = buffer.data(fg_s + 95);
    const auto *fg_s_96 = buffer.data(fg_s + 96);
    const auto *fg_s_97 = buffer.data(fg_s + 97);
    const auto *fg_s_98 = buffer.data(fg_s + 98);
    const auto *fg_s_99 = buffer.data(fg_s + 99);
    const auto *fg_s_100 = buffer.data(fg_s + 100);
    const auto *fg_s_101 = buffer.data(fg_s + 101);
    const auto *fg_s_102 = buffer.data(fg_s + 102);
    const auto *fg_s_103 = buffer.data(fg_s + 103);
    const auto *fg_s_104 = buffer.data(fg_s + 104);
    const auto *fg_s_105 = buffer.data(fg_s + 105);
    const auto *fg_s_106 = buffer.data(fg_s + 106);
    const auto *fg_s_107 = buffer.data(fg_s + 107);
    const auto *fg_s_108 = buffer.data(fg_s + 108);
    const auto *fg_s_109 = buffer.data(fg_s + 109);
    const auto *fg_s_110 = buffer.data(fg_s + 110);
    const auto *fg_s_111 = buffer.data(fg_s + 111);
    const auto *fg_s_112 = buffer.data(fg_s + 112);
    const auto *fg_s_113 = buffer.data(fg_s + 113);
    const auto *fg_s_114 = buffer.data(fg_s + 114);
    const auto *fg_s_115 = buffer.data(fg_s + 115);
    const auto *fg_s_116 = buffer.data(fg_s + 116);
    const auto *fg_s_117 = buffer.data(fg_s + 117);
    const auto *fg_s_118 = buffer.data(fg_s + 118);
    const auto *fg_s_119 = buffer.data(fg_s + 119);
    const auto *fg_s_120 = buffer.data(fg_s + 120);
    const auto *fg_s_121 = buffer.data(fg_s + 121);
    const auto *fg_s_122 = buffer.data(fg_s + 122);
    const auto *fg_s_123 = buffer.data(fg_s + 123);
    const auto *fg_s_124 = buffer.data(fg_s + 124);
    const auto *fg_s_125 = buffer.data(fg_s + 125);
    const auto *fg_s_126 = buffer.data(fg_s + 126);
    const auto *fg_s_127 = buffer.data(fg_s + 127);
    const auto *fg_s_128 = buffer.data(fg_s + 128);
    const auto *fg_s_129 = buffer.data(fg_s + 129);
    const auto *fg_s_130 = buffer.data(fg_s + 130);
    const auto *fg_s_131 = buffer.data(fg_s + 131);
    const auto *fg_s_132 = buffer.data(fg_s + 132);
    const auto *fg_s_133 = buffer.data(fg_s + 133);
    const auto *fg_s_134 = buffer.data(fg_s + 134);
    const auto *fg_s_135 = buffer.data(fg_s + 135);
    const auto *fg_s_136 = buffer.data(fg_s + 136);
    const auto *fg_s_137 = buffer.data(fg_s + 137);
    const auto *fg_s_138 = buffer.data(fg_s + 138);
    const auto *fg_s_139 = buffer.data(fg_s + 139);
    const auto *fg_s_140 = buffer.data(fg_s + 140);
    const auto *fg_s_141 = buffer.data(fg_s + 141);
    const auto *fg_s_142 = buffer.data(fg_s + 142);
    const auto *fg_s_143 = buffer.data(fg_s + 143);
    const auto *fg_s_144 = buffer.data(fg_s + 144);
    const auto *fg_s_145 = buffer.data(fg_s + 145);
    const auto *fg_s_146 = buffer.data(fg_s + 146);
    const auto *fg_s_147 = buffer.data(fg_s + 147);
    const auto *fg_s_148 = buffer.data(fg_s + 148);
    const auto *fg_s_149 = buffer.data(fg_s + 149);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fg_s_3, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_2 * fg_s_1[k]
                 + pb_y[k] * ff_0[k];

        t_2[k] = f_2 * fg_s_2[k]
                 + pb_z[k] * ff_0[k];

        t_3[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_3[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, df_3, fd_s_0, fg_s_4, fg_s_5, \
                         fg_s_6, fd_0, ff_2, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * fg_s_4[k]
                 + pb_y[k] * ff_2[k];

        t_5[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_5[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_6[k] = f_0 * df_3[k]
                 + f_2 * fg_s_6[k]
                 + pb_x[k] * ff_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, df_4, fg_s_7, fg_s_8, fg_s_9, ff_3, \
                         ff_4, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * fg_s_7[k]
                 + pb_z[k] * ff_3[k];

        t_8[k] = f_2 * fg_s_8[k]
                 + pb_y[k] * ff_4[k];

        t_9[k] = f_0 * df_4[k]
                 + f_2 * fg_s_9[k]
                 + pb_x[k] * ff_7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, fd_s_1, fd_s_2, fg_s_10, fg_s_11, \
                         fg_s_12, fd_1, fd_2, ff_5, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * fd_s_1[k]
                  + f_2 * fg_s_10[k]
                  + f_0 * fd_1[k]
                  + pb_y[k] * ff_5[k];

        t_11[k] = f_2 * fg_s_11[k]
                  + pb_z[k] * ff_5[k];

        t_12[k] = -f_3 * fd_s_2[k]
                  + f_2 * fg_s_12[k]
                  + f_4 * fd_2[k]
                  + pb_y[k] * ff_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, dg_0, fd_s_2, fg_s_13, fg_s_14, \
                         fg_s_15, fd_2, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * fg_s_13[k]
                  + pb_y[k] * ff_7[k];

        t_14[k] = -f_1 * fd_s_2[k]
                  + f_2 * fg_s_14[k]
                  + f_0 * fd_2[k]
                  + pb_z[k] * ff_7[k];

        t_15[k] = pa_y[k] * dg_0[k]
                  + f_2 * fg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, df_0, df_1, dg_1, fg_s_16, \
                         fg_s_17, fg_s_18, fg_s_19, ff_8, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * df_0[k]
                  + f_2 * fg_s_16[k]
                  + pb_y[k] * ff_8[k];

        t_17[k] = f_2 * fg_s_17[k]
                  + pb_z[k] * ff_8[k];

        t_18[k] = f_5 * df_1[k]
                  + pa_y[k] * dg_1[k]
                  + f_2 * fg_s_18[k];

        t_19[k] = f_2 * fg_s_19[k]
                  + pb_z[k] * ff_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_x, pb_z, df_6, dg_2, fg_s_20, fg_s_21, \
                         fg_s_22, ff_10, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * dg_2[k]
                  + f_2 * fg_s_20[k];

        t_21[k] = f_5 * df_6[k]
                  + f_2 * fg_s_21[k]
                  + pb_x[k] * ff_11[k];

        t_22[k] = f_2 * fg_s_22[k]
                  + pb_z[k] * ff_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, pb_x, pg_s_2, pg_2, df_7, dg_4, dg_10, \
                         fg_s_23, fg_s_24, fg_s_25, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * df_7[k]
                  + f_2 * fg_s_23[k]
                  + pb_x[k] * ff_13[k];

        t_24[k] = pa_y[k] * dg_4[k]
                  + f_2 * fg_s_24[k];

        t_25[k] = -f_6 * pg_s_2[k]
                  + f_4 * pg_2[k]
                  + pa_x[k] * dg_10[k]
                  + f_2 * fg_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_y, pb_z, df_4, fd_s_4, fg_s_26, fg_s_27, \
                         fg_s_28, fd_4, ff_11, ff_12, ff_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * fg_s_26[k]
                  + pb_z[k] * ff_11[k];

        t_27[k] = -f_3 * fd_s_4[k]
                  + f_2 * fg_s_27[k]
                  + f_4 * fd_4[k]
                  + pb_z[k] * ff_12[k];

        t_28[k] = f_4 * df_4[k]
                  + f_2 * fg_s_28[k]
                  + pb_y[k] * ff_14[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pb_y, pb_z, df_0, dg_0, dg_6, \
                         fg_s_29, fg_s_30, fg_s_31, fg_s_32, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * dg_6[k]
                  + f_2 * fg_s_29[k];

        t_30[k] = pa_z[k] * dg_0[k]
                  + f_2 * fg_s_30[k];

        t_31[k] = f_2 * fg_s_31[k]
                  + pb_y[k] * ff_15[k];

        t_32[k] = f_4 * df_0[k]
                  + f_2 * fg_s_32[k]
                  + pb_z[k] * ff_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_y, df_2, dg_1, dg_2, dg_3, fg_s_33, \
                         fg_s_34, fg_s_35, fg_s_36, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * dg_1[k]
                  + f_2 * fg_s_33[k];

        t_34[k] = f_2 * fg_s_34[k]
                  + pb_y[k] * ff_16[k];

        t_35[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_2[k]
                  + f_2 * fg_s_35[k];

        t_36[k] = pa_z[k] * dg_3[k]
                  + f_2 * fg_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, df_10, df_11, fg_s_37, fg_s_38, \
                         fg_s_39, ff_17, ff_18, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * df_10[k]
                  + f_2 * fg_s_37[k]
                  + pb_x[k] * ff_18[k];

        t_38[k] = f_2 * fg_s_38[k]
                  + pb_y[k] * ff_17[k];

        t_39[k] = f_5 * df_11[k]
                  + f_2 * fg_s_39[k]
                  + pb_x[k] * ff_20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_z, pb_y, dg_5, fd_s_7, fd_s_8, fg_s_40, fg_s_41, \
                         fg_s_42, fd_7, fd_8, ff_18, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_z[k] * dg_5[k]
                  + f_2 * fg_s_40[k];

        t_41[k] = -f_7 * fd_s_7[k]
                  + f_2 * fg_s_41[k]
                  + f_5 * fd_7[k]
                  + pb_y[k] * ff_18[k];

        t_42[k] = -f_3 * fd_s_8[k]
                  + f_2 * fg_s_42[k]
                  + f_4 * fd_8[k]
                  + pb_y[k] * ff_19[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pb_y, pg_s_5, pg_5, df_12, dg_15, dg_16, \
                         fg_s_43, fg_s_44, fg_s_45, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * fg_s_43[k]
                  + pb_y[k] * ff_20[k];

        t_44[k] = -f_6 * pg_s_5[k]
                  + f_4 * pg_5[k]
                  + pa_x[k] * dg_15[k]
                  + f_2 * fg_s_44[k];

        t_45[k] = f_8 * df_12[k]
                  + pa_x[k] * dg_16[k]
                  + f_2 * fg_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_y, pb_z, df_5, df_14, dg_18, \
                         fg_s_46, fg_s_47, fg_s_48, fg_s_49, ff_21, \
                         ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_5 * df_5[k]
                  + f_2 * fg_s_46[k]
                  + pb_y[k] * ff_21[k];

        t_47[k] = f_2 * fg_s_47[k]
                  + pb_z[k] * ff_21[k];

        t_48[k] = f_5 * df_14[k]
                  + pa_x[k] * dg_18[k]
                  + f_2 * fg_s_48[k];

        t_49[k] = f_2 * fg_s_49[k]
                  + pb_z[k] * ff_22[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pb_x, pb_z, df_15, df_16, dg_20, fg_s_50, \
                         fg_s_51, fg_s_52, ff_23, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_5 * df_15[k]
                  + pa_x[k] * dg_20[k]
                  + f_2 * fg_s_50[k];

        t_51[k] = f_4 * df_16[k]
                  + f_2 * fg_s_51[k]
                  + pb_x[k] * ff_24[k];

        t_52[k] = f_2 * fg_s_52[k]
                  + pb_z[k] * ff_23[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_x, pb_x, df_18, df_19, dg_21, fg_s_53, fg_s_54, \
                         fg_s_55, ff_25, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_4 * df_18[k]
                  + f_2 * fg_s_53[k]
                  + pb_x[k] * ff_25[k];

        t_54[k] = f_4 * df_19[k]
                  + f_2 * fg_s_54[k]
                  + pb_x[k] * ff_26[k];

        t_55[k] = pa_x[k] * dg_21[k]
                  + f_2 * fg_s_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_z, dg_22, dg_23, dg_24, fg_s_56, \
                         fg_s_57, fg_s_58, fg_s_59, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_2 * fg_s_56[k]
                  + pb_z[k] * ff_24[k];

        t_57[k] = pa_x[k] * dg_22[k]
                  + f_2 * fg_s_57[k];

        t_58[k] = pa_x[k] * dg_23[k]
                  + f_2 * fg_s_58[k];

        t_59[k] = pa_x[k] * dg_24[k]
                  + f_2 * fg_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pa_z, dg_7, dg_8, dg_11, dg_12, \
                         fg_s_60, fg_s_61, fg_s_62, fg_s_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * dg_11[k]
                  + f_2 * fg_s_60[k];

        t_61[k] = pa_z[k] * dg_7[k]
                  + f_2 * fg_s_61[k];

        t_62[k] = pa_y[k] * dg_12[k]
                  + f_2 * fg_s_62[k];

        t_63[k] = pa_z[k] * dg_8[k]
                  + f_2 * fg_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_y, pa_z, pb_y, df_9, dg_9, dg_13, fg_s_64, \
                         fg_s_65, fg_s_66, ff_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_4 * df_9[k]
                  + f_2 * fg_s_64[k]
                  + pb_y[k] * ff_27[k];

        t_65[k] = pa_y[k] * dg_13[k]
                  + f_2 * fg_s_65[k];

        t_66[k] = pa_z[k] * dg_9[k]
                  + f_2 * fg_s_66[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_y, pb_x, df_21, df_22, dg_14, fg_s_67, fg_s_68, \
                         fg_s_69, ff_28, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_4 * df_21[k]
                  + f_2 * fg_s_67[k]
                  + pb_x[k] * ff_28[k];

        t_68[k] = f_4 * df_22[k]
                  + f_2 * fg_s_68[k]
                  + pb_x[k] * ff_29[k];

        t_69[k] = pa_y[k] * dg_14[k]
                  + f_2 * fg_s_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_x, dg_25, dg_26, dg_27, dg_28, \
                         dg_29, fg_s_70, fg_s_71, fg_s_72, fg_s_73, \
                         fg_s_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_x[k] * dg_25[k]
                  + f_2 * fg_s_70[k];

        t_71[k] = pa_x[k] * dg_26[k]
                  + f_2 * fg_s_71[k];

        t_72[k] = pa_x[k] * dg_27[k]
                  + f_2 * fg_s_72[k];

        t_73[k] = pa_x[k] * dg_28[k]
                  + f_2 * fg_s_73[k];

        t_74[k] = pa_x[k] * dg_29[k]
                  + f_2 * fg_s_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_x, pb_y, pb_z, df_8, df_24, dg_30, fg_s_75, \
                         fg_s_76, fg_s_77, ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * df_24[k]
                  + pa_x[k] * dg_30[k]
                  + f_2 * fg_s_75[k];

        t_76[k] = f_2 * fg_s_76[k]
                  + pb_y[k] * ff_30[k];

        t_77[k] = f_5 * df_8[k]
                  + f_2 * fg_s_77[k]
                  + pb_z[k] * ff_30[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_x, pb_y, df_27, df_28, dg_33, dg_35, fg_s_78, \
                         fg_s_79, fg_s_80, ff_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_5 * df_27[k]
                  + pa_x[k] * dg_33[k]
                  + f_2 * fg_s_78[k];

        t_79[k] = f_2 * fg_s_79[k]
                  + pb_y[k] * ff_31[k];

        t_80[k] = f_5 * df_28[k]
                  + pa_x[k] * dg_35[k]
                  + f_2 * fg_s_80[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_x, pb_y, df_29, df_30, fg_s_81, fg_s_82, \
                         fg_s_83, ff_32, ff_33, ff_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * df_29[k]
                  + f_2 * fg_s_81[k]
                  + pb_x[k] * ff_33[k];

        t_82[k] = f_4 * df_30[k]
                  + f_2 * fg_s_82[k]
                  + pb_x[k] * ff_34[k];

        t_83[k] = f_2 * fg_s_83[k]
                  + pb_y[k] * ff_32[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pb_x, df_32, dg_36, dg_37, dg_38, \
                         fg_s_84, fg_s_85, fg_s_86, fg_s_87, ff_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * df_32[k]
                  + f_2 * fg_s_84[k]
                  + pb_x[k] * ff_35[k];

        t_85[k] = pa_x[k] * dg_36[k]
                  + f_2 * fg_s_85[k];

        t_86[k] = pa_x[k] * dg_37[k]
                  + f_2 * fg_s_86[k];

        t_87[k] = pa_x[k] * dg_38[k]
                  + f_2 * fg_s_87[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pa_x, pb_x, pb_y, dg_39, fd_s_13, fg_s_88, fg_s_89, \
                         fg_s_90, fd_13, ff_35, ff_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * fg_s_88[k]
                  + pb_y[k] * ff_35[k];

        t_89[k] = pa_x[k] * dg_39[k]
                  + f_2 * fg_s_89[k];

        t_90[k] = -f_1 * fd_s_13[k]
                  + f_2 * fg_s_90[k]
                  + f_0 * fd_13[k]
                  + pb_x[k] * ff_36[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_x, pb_z, fd_s_14, fd_s_15, fg_s_91, fg_s_92, \
                         fg_s_93, fd_14, fd_15, ff_36, ff_37, ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = -f_7 * fd_s_14[k]
                  + f_2 * fg_s_91[k]
                  + f_5 * fd_14[k]
                  + pb_x[k] * ff_37[k];

        t_92[k] = f_2 * fg_s_92[k]
                  + pb_z[k] * ff_36[k];

        t_93[k] = -f_3 * fd_s_15[k]
                  + f_2 * fg_s_93[k]
                  + f_4 * fd_15[k]
                  + pb_x[k] * ff_38[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_x, pb_z, fd_s_16, fg_s_94, fg_s_95, \
                         fg_s_96, fg_s_97, fd_16, ff_37, ff_39, ff_40, \
                         ff_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_2 * fg_s_94[k]
                  + pb_z[k] * ff_37[k];

        t_95[k] = -f_3 * fd_s_16[k]
                  + f_2 * fg_s_95[k]
                  + f_4 * fd_16[k]
                  + pb_x[k] * ff_39[k];

        t_96[k] = f_2 * fg_s_96[k]
                  + pb_x[k] * ff_40[k];

        t_97[k] = f_2 * fg_s_97[k]
                  + pb_x[k] * ff_41[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pb_x, pb_y, df_16, fd_s_15, fg_s_98, fg_s_99, \
                         fg_s_100, fd_15, ff_40, ff_42, ff_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_2 * fg_s_98[k]
                  + pb_x[k] * ff_42[k];

        t_99[k] = f_2 * fg_s_99[k]
                  + pb_x[k] * ff_43[k];

        t_100[k] = f_0 * df_16[k]
                   - f_1 * fd_s_15[k]
                   + f_2 * fg_s_100[k]
                   + f_0 * fd_15[k]
                   + pb_y[k] * ff_40[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pb_y, pb_z, df_19, fd_s_15, fg_s_101, fg_s_102, \
                         fg_s_103, fd_15, ff_40, ff_41, ff_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_2 * fg_s_101[k]
                   + pb_z[k] * ff_40[k];

        t_102[k] = -f_3 * fd_s_15[k]
                   + f_2 * fg_s_102[k]
                   + f_4 * fd_15[k]
                   + pb_z[k] * ff_41[k];

        t_103[k] = f_0 * df_19[k]
                   + f_2 * fg_s_103[k]
                   + pb_y[k] * ff_43[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pb_z, dg_16, dg_17, fd_s_16, fg_s_104, \
                         fg_s_105, fg_s_106, fd_16, ff_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -f_1 * fd_s_16[k]
                   + f_2 * fg_s_104[k]
                   + f_0 * fd_16[k]
                   + pb_z[k] * ff_43[k];

        t_105[k] = pa_z[k] * dg_16[k]
                   + f_2 * fg_s_105[k];

        t_106[k] = pa_z[k] * dg_17[k]
                   + f_2 * fg_s_106[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pa_z, pb_x, df_13, dg_18, dg_19, fd_s_17, \
                         fg_s_107, fg_s_108, fg_s_109, fd_17, ff_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -f_7 * fd_s_17[k]
                   + f_2 * fg_s_107[k]
                   + f_5 * fd_17[k]
                   + pb_x[k] * ff_44[k];

        t_108[k] = pa_z[k] * dg_18[k]
                   + f_2 * fg_s_108[k];

        t_109[k] = f_4 * df_13[k]
                   + pa_z[k] * dg_19[k]
                   + f_2 * fg_s_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_x, fd_s_20, fg_s_110, fg_s_111, \
                         fg_s_112, fg_s_113, fd_19, ff_45, ff_46, ff_47, \
                         ff_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -f_3 * fd_s_20[k]
                   + f_2 * fg_s_110[k]
                   + f_4 * fd_19[k]
                   + pb_x[k] * ff_45[k];

        t_111[k] = f_2 * fg_s_111[k]
                   + pb_x[k] * ff_46[k];

        t_112[k] = f_2 * fg_s_112[k]
                   + pb_x[k] * ff_47[k];

        t_113[k] = f_2 * fg_s_113[k]
                   + pb_x[k] * ff_48[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_z, pb_x, pb_z, df_16, dg_21, fg_s_114, \
                         fg_s_115, fg_s_116, ff_46, ff_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_2 * fg_s_114[k]
                   + pb_x[k] * ff_49[k];

        t_115[k] = pa_z[k] * dg_21[k]
                   + f_2 * fg_s_115[k];

        t_116[k] = f_4 * df_16[k]
                   + f_2 * fg_s_116[k]
                   + pb_z[k] * ff_46[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_y, pa_z, pb_y, pg_s_5, pg_5, df_17, df_23, \
                         dg_22, dg_29, fg_s_117, fg_s_118, fg_s_119, \
                         ff_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * df_17[k]
                   + pa_z[k] * dg_22[k]
                   + f_2 * fg_s_117[k];

        t_118[k] = f_5 * df_23[k]
                   + f_2 * fg_s_118[k]
                   + pb_y[k] * ff_49[k];

        t_119[k] = -f_6 * pg_s_5[k]
                   + f_4 * pg_5[k]
                   + pa_y[k] * dg_29[k]
                   + f_2 * fg_s_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, df_24, df_25, dg_30, dg_31, dg_32, \
                         dg_33, fg_s_120, fg_s_121, fg_s_122, \
                         fg_s_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_y[k] * dg_30[k]
                   + f_2 * fg_s_120[k];

        t_121[k] = f_4 * df_24[k]
                   + pa_y[k] * dg_31[k]
                   + f_2 * fg_s_121[k];

        t_122[k] = pa_y[k] * dg_32[k]
                   + f_2 * fg_s_122[k];

        t_123[k] = f_5 * df_25[k]
                   + pa_y[k] * dg_33[k]
                   + f_2 * fg_s_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pb_x, df_26, dg_34, dg_35, \
                         fg_s_124, fg_s_125, fg_s_126, fg_s_127, ff_50, \
                         ff_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_4 * df_26[k]
                   + pa_y[k] * dg_34[k]
                   + f_2 * fg_s_124[k];

        t_125[k] = pa_y[k] * dg_35[k]
                   + f_2 * fg_s_125[k];

        t_126[k] = f_2 * fg_s_126[k]
                   + pb_x[k] * ff_50[k];

        t_127[k] = f_2 * fg_s_127[k]
                   + pb_x[k] * ff_51[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pa_y, pb_x, df_29, dg_36, fg_s_128, fg_s_129, \
                         fg_s_130, ff_52, ff_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_2 * fg_s_128[k]
                   + pb_x[k] * ff_52[k];

        t_129[k] = f_2 * fg_s_129[k]
                   + pb_x[k] * ff_53[k];

        t_130[k] = f_8 * df_29[k]
                   + pa_y[k] * dg_36[k]
                   + f_2 * fg_s_130[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, pa_y, pb_y, pb_z, df_20, df_31, df_32, dg_38, \
                         fg_s_131, fg_s_132, fg_s_133, ff_50, ff_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_5 * df_20[k]
                   + f_2 * fg_s_131[k]
                   + pb_z[k] * ff_50[k];

        t_132[k] = f_5 * df_31[k]
                   + pa_y[k] * dg_38[k]
                   + f_2 * fg_s_132[k];

        t_133[k] = f_4 * df_32[k]
                   + f_2 * fg_s_133[k]
                   + pb_y[k] * ff_53[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pa_y, pb_x, pb_y, dg_39, fd_s_25, fg_s_134, \
                         fg_s_135, fg_s_136, fd_22, ff_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_y[k] * dg_39[k]
                   + f_2 * fg_s_134[k];

        t_135[k] = -f_1 * fd_s_25[k]
                   + f_2 * fg_s_135[k]
                   + f_0 * fd_22[k]
                   + pb_x[k] * ff_54[k];

        t_136[k] = f_2 * fg_s_136[k]
                   + pb_y[k] * ff_54[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pb_x, pb_y, fd_s_26, fd_s_27, fg_s_137, \
                         fg_s_138, fg_s_139, fd_23, fd_24, ff_55, \
                         ff_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -f_7 * fd_s_26[k]
                   + f_2 * fg_s_137[k]
                   + f_5 * fd_23[k]
                   + pb_x[k] * ff_55[k];

        t_138[k] = -f_3 * fd_s_27[k]
                   + f_2 * fg_s_138[k]
                   + f_4 * fd_24[k]
                   + pb_x[k] * ff_56[k];

        t_139[k] = f_2 * fg_s_139[k]
                   + pb_y[k] * ff_55[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pb_x, fd_s_29, fg_s_140, fg_s_141, \
                         fg_s_142, fg_s_143, fd_26, ff_57, ff_58, ff_59, \
                         ff_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -f_3 * fd_s_29[k]
                   + f_2 * fg_s_140[k]
                   + f_4 * fd_26[k]
                   + pb_x[k] * ff_57[k];

        t_141[k] = f_2 * fg_s_141[k]
                   + pb_x[k] * ff_58[k];

        t_142[k] = f_2 * fg_s_142[k]
                   + pb_x[k] * ff_59[k];

        t_143[k] = f_2 * fg_s_143[k]
                   + pb_x[k] * ff_60[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pb_x, pb_y, fd_s_27, fd_s_28, fg_s_144, \
                         fg_s_145, fg_s_146, fd_24, fd_25, ff_58, ff_59, \
                         ff_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_2 * fg_s_144[k]
                   + pb_x[k] * ff_61[k];

        t_145[k] = -f_1 * fd_s_27[k]
                   + f_2 * fg_s_145[k]
                   + f_0 * fd_24[k]
                   + pb_y[k] * ff_58[k];

        t_146[k] = -f_7 * fd_s_28[k]
                   + f_2 * fg_s_146[k]
                   + f_5 * fd_25[k]
                   + pb_y[k] * ff_59[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pb_y, pb_z, df_32, fd_s_29, fg_s_147, fg_s_148, \
                         fg_s_149, fd_26, ff_60, ff_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -f_3 * fd_s_29[k]
                   + f_2 * fg_s_147[k]
                   + f_4 * fd_26[k]
                   + pb_y[k] * ff_60[k];

        t_148[k] = f_2 * fg_s_148[k]
                   + pb_y[k] * ff_61[k];

        t_149[k] = f_0 * df_32[k]
                   - f_1 * fd_s_29[k]
                   + f_2 * fg_s_149[k]
                   + f_0 * fd_26[k]
                   + pb_z[k] * ff_61[k];
    }
}

auto
compute_prim_fg_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * alpha / p;
    const auto f_8 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_s_6 = buffer.data(pg_s + 6);
    const auto *pg_s_15 = buffer.data(pg_s + 15);

    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_15 = buffer.data(pg + 15);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_51 = buffer.data(dg + 51);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_29 = buffer.data(fg_s + 29);
    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_31 = buffer.data(fg_s + 31);
    const auto *fg_s_33 = buffer.data(fg_s + 33);
    const auto *fg_s_34 = buffer.data(fg_s + 34);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_36 = buffer.data(fg_s + 36);
    const auto *fg_s_37 = buffer.data(fg_s + 37);
    const auto *fg_s_38 = buffer.data(fg_s + 38);
    const auto *fg_s_39 = buffer.data(fg_s + 39);
    const auto *fg_s_40 = buffer.data(fg_s + 40);
    const auto *fg_s_41 = buffer.data(fg_s + 41);
    const auto *fg_s_42 = buffer.data(fg_s + 42);
    const auto *fg_s_43 = buffer.data(fg_s + 43);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_45 = buffer.data(fg_s + 45);
    const auto *fg_s_46 = buffer.data(fg_s + 46);
    const auto *fg_s_48 = buffer.data(fg_s + 48);
    const auto *fg_s_51 = buffer.data(fg_s + 51);
    const auto *fg_s_52 = buffer.data(fg_s + 52);
    const auto *fg_s_53 = buffer.data(fg_s + 53);
    const auto *fg_s_54 = buffer.data(fg_s + 54);
    const auto *fg_s_55 = buffer.data(fg_s + 55);
    const auto *fg_s_56 = buffer.data(fg_s + 56);
    const auto *fg_s_57 = buffer.data(fg_s + 57);
    const auto *fg_s_58 = buffer.data(fg_s + 58);
    const auto *fg_s_59 = buffer.data(fg_s + 59);
    const auto *fg_s_60 = buffer.data(fg_s + 60);
    const auto *fg_s_61 = buffer.data(fg_s + 61);
    const auto *fg_s_62 = buffer.data(fg_s + 62);
    const auto *fg_s_63 = buffer.data(fg_s + 63);
    const auto *fg_s_64 = buffer.data(fg_s + 64);
    const auto *fg_s_65 = buffer.data(fg_s + 65);
    const auto *fg_s_66 = buffer.data(fg_s + 66);
    const auto *fg_s_67 = buffer.data(fg_s + 67);
    const auto *fg_s_68 = buffer.data(fg_s + 68);
    const auto *fg_s_69 = buffer.data(fg_s + 69);
    const auto *fg_s_71 = buffer.data(fg_s + 71);
    const auto *fg_s_74 = buffer.data(fg_s + 74);
    const auto *fg_s_75 = buffer.data(fg_s + 75);
    const auto *fg_s_76 = buffer.data(fg_s + 76);
    const auto *fg_s_77 = buffer.data(fg_s + 77);
    const auto *fg_s_78 = buffer.data(fg_s + 78);
    const auto *fg_s_79 = buffer.data(fg_s + 79);
    const auto *fg_s_86 = buffer.data(fg_s + 86);
    const auto *fg_s_87 = buffer.data(fg_s + 87);
    const auto *fg_s_88 = buffer.data(fg_s + 88);
    const auto *fg_s_89 = buffer.data(fg_s + 89);
    const auto *fg_s_90 = buffer.data(fg_s + 90);
    const auto *fg_s_91 = buffer.data(fg_s + 91);
    const auto *fg_s_92 = buffer.data(fg_s + 92);
    const auto *fg_s_93 = buffer.data(fg_s + 93);
    const auto *fg_s_94 = buffer.data(fg_s + 94);
    const auto *fg_s_95 = buffer.data(fg_s + 95);
    const auto *fg_s_96 = buffer.data(fg_s + 96);
    const auto *fg_s_97 = buffer.data(fg_s + 97);
    const auto *fg_s_98 = buffer.data(fg_s + 98);
    const auto *fg_s_99 = buffer.data(fg_s + 99);
    const auto *fg_s_100 = buffer.data(fg_s + 100);
    const auto *fg_s_101 = buffer.data(fg_s + 101);
    const auto *fg_s_102 = buffer.data(fg_s + 102);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fg_s_3, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_2 * fg_s_1[k]
                 + pb_y[k] * ff_0[k];

        t_2[k] = f_2 * fg_s_2[k]
                 + pb_z[k] * ff_0[k];

        t_3[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_3[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_z, df_3, df_5, fd_s_0, fg_s_4, fg_s_5, \
                         fg_s_6, fd_0, ff_2, ff_3, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_4[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_5[k] = f_0 * df_3[k]
                 + f_2 * fg_s_5[k]
                 + pb_x[k] * ff_3[k];

        t_6[k] = f_0 * df_5[k]
                 + f_2 * fg_s_6[k]
                 + pb_x[k] * ff_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, fd_s_1, fd_s_2, fg_s_7, fg_s_8, fg_s_9, fd_1, \
                         fd_2, ff_3, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_1 * fd_s_1[k]
                 + f_2 * fg_s_7[k]
                 + f_0 * fd_1[k]
                 + pb_y[k] * ff_3[k];

        t_8[k] = -f_3 * fd_s_2[k]
                 + f_2 * fg_s_8[k]
                 + f_4 * fd_2[k]
                 + pb_y[k] * ff_4[k];

        t_9[k] = f_2 * fg_s_9[k]
                 + pb_y[k] * ff_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_y, pb_z, df_0, dg_0, fd_s_2, fg_s_10, \
                         fg_s_11, fg_s_12, fd_2, ff_5, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * fd_s_2[k]
                  + f_2 * fg_s_10[k]
                  + f_0 * fd_2[k]
                  + pb_z[k] * ff_5[k];

        t_11[k] = pa_y[k] * dg_0[k]
                  + f_2 * fg_s_11[k];

        t_12[k] = f_4 * df_0[k]
                  + f_2 * fg_s_12[k]
                  + pb_y[k] * ff_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, df_1, df_7, dg_3, dg_4, fg_s_13, \
                         fg_s_14, fg_s_15, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * df_1[k]
                  + pa_y[k] * dg_3[k]
                  + f_2 * fg_s_13[k];

        t_14[k] = pa_y[k] * dg_4[k]
                  + f_2 * fg_s_14[k];

        t_15[k] = f_5 * df_7[k]
                  + f_2 * fg_s_15[k]
                  + pb_x[k] * ff_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_z, pg_s_6, pg_6, dg_11, fd_s_3, fg_s_16, \
                         fg_s_17, fg_s_18, fd_3, ff_7, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_6 * pg_s_6[k]
                  + f_4 * pg_6[k]
                  + pa_x[k] * dg_11[k]
                  + f_2 * fg_s_16[k];

        t_17[k] = f_2 * fg_s_17[k]
                  + pb_z[k] * ff_7[k];

        t_18[k] = -f_3 * fd_s_3[k]
                  + f_2 * fg_s_18[k]
                  + f_4 * fd_3[k]
                  + pb_z[k] * ff_8[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_y, df_5, dg_0, dg_7, fg_s_19, \
                         fg_s_20, fg_s_21, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_4 * df_5[k]
                  + f_2 * fg_s_19[k]
                  + pb_y[k] * ff_9[k];

        t_20[k] = pa_y[k] * dg_7[k]
                  + f_2 * fg_s_20[k];

        t_21[k] = pa_z[k] * dg_0[k]
                  + f_2 * fg_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_z, pb_x, pb_z, df_0, df_2, df_9, dg_4, fg_s_22, \
                         fg_s_24, fg_s_25, ff_10, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_4 * df_0[k]
                  + f_2 * fg_s_22[k]
                  + pb_z[k] * ff_10[k];

        t_23[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_4[k]
                  + f_2 * fg_s_24[k];

        t_24[k] = f_5 * df_9[k]
                  + f_2 * fg_s_25[k]
                  + pb_x[k] * ff_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_y, fd_s_5, fd_s_6, fg_s_26, fg_s_27, fg_s_28, \
                         fd_4, fd_5, ff_11, ff_12, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_7 * fd_s_5[k]
                  + f_2 * fg_s_26[k]
                  + f_5 * fd_4[k]
                  + pb_y[k] * ff_11[k];

        t_26[k] = -f_3 * fd_s_6[k]
                  + f_2 * fg_s_27[k]
                  + f_4 * fd_5[k]
                  + pb_y[k] * ff_12[k];

        t_27[k] = f_2 * fg_s_28[k]
                  + pb_y[k] * ff_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pb_y, pg_s_15, pg_15, df_6, df_10, dg_19, \
                         dg_20, fg_s_29, fg_s_30, fg_s_31, ff_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_6 * pg_s_15[k]
                  + f_4 * pg_15[k]
                  + pa_x[k] * dg_19[k]
                  + f_2 * fg_s_29[k];

        t_29[k] = f_8 * df_10[k]
                  + pa_x[k] * dg_20[k]
                  + f_2 * fg_s_30[k];

        t_30[k] = f_5 * df_6[k]
                  + f_2 * fg_s_31[k]
                  + pb_y[k] * ff_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_x, pb_x, df_12, df_13, df_14, dg_22, dg_24, \
                         fg_s_33, fg_s_34, fg_s_35, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_5 * df_12[k]
                  + pa_x[k] * dg_22[k]
                  + f_2 * fg_s_33[k];

        t_32[k] = f_5 * df_13[k]
                  + pa_x[k] * dg_24[k]
                  + f_2 * fg_s_34[k];

        t_33[k] = f_4 * df_14[k]
                  + f_2 * fg_s_35[k]
                  + pb_x[k] * ff_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, dg_28, dg_30, dg_31, dg_32, fg_s_36, \
                         fg_s_37, fg_s_38, fg_s_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_x[k] * dg_28[k]
                  + f_2 * fg_s_36[k];

        t_35[k] = pa_x[k] * dg_30[k]
                  + f_2 * fg_s_37[k];

        t_36[k] = pa_x[k] * dg_31[k]
                  + f_2 * fg_s_38[k];

        t_37[k] = pa_x[k] * dg_32[k]
                  + f_2 * fg_s_39[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pa_y, pa_z, dg_9, dg_15, dg_16, dg_34, \
                         fg_s_40, fg_s_41, fg_s_42, fg_s_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * dg_15[k]
                  + f_2 * fg_s_40[k];

        t_39[k] = pa_z[k] * dg_9[k]
                  + f_2 * fg_s_41[k];

        t_40[k] = pa_y[k] * dg_16[k]
                  + f_2 * fg_s_42[k];

        t_41[k] = pa_x[k] * dg_34[k]
                  + f_2 * fg_s_43[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_z, df_8, df_20, dg_35, dg_36, dg_38, \
                         fg_s_44, fg_s_45, fg_s_46, fg_s_48, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_x[k] * dg_35[k]
                  + f_2 * fg_s_44[k];

        t_43[k] = pa_x[k] * dg_36[k]
                  + f_2 * fg_s_45[k];

        t_44[k] = f_8 * df_20[k]
                  + pa_x[k] * dg_38[k]
                  + f_2 * fg_s_46[k];

        t_45[k] = f_5 * df_8[k]
                  + f_2 * fg_s_48[k]
                  + pb_z[k] * ff_16[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_x, df_23, df_27, dg_43, dg_47, \
                         dg_48, fg_s_51, fg_s_52, fg_s_53, fg_s_54, \
                         ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_5 * df_23[k]
                  + pa_x[k] * dg_43[k]
                  + f_2 * fg_s_51[k];

        t_47[k] = f_4 * df_27[k]
                  + f_2 * fg_s_52[k]
                  + pb_x[k] * ff_17[k];

        t_48[k] = pa_x[k] * dg_47[k]
                  + f_2 * fg_s_53[k];

        t_49[k] = pa_x[k] * dg_48[k]
                  + f_2 * fg_s_54[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pb_x, dg_49, dg_51, fd_s_9, fg_s_55, fg_s_56, \
                         fg_s_57, fd_6, ff_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_x[k] * dg_49[k]
                  + f_2 * fg_s_55[k];

        t_51[k] = pa_x[k] * dg_51[k]
                  + f_2 * fg_s_56[k];

        t_52[k] = -f_1 * fd_s_9[k]
                  + f_2 * fg_s_57[k]
                  + f_0 * fd_6[k]
                  + pb_x[k] * ff_18[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, fd_s_10, fd_s_11, fd_s_12, fg_s_58, fg_s_59, \
                         fg_s_60, fd_7, fd_8, fd_9, ff_19, ff_20, \
                         ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_7 * fd_s_10[k]
                  + f_2 * fg_s_58[k]
                  + f_5 * fd_7[k]
                  + pb_x[k] * ff_19[k];

        t_54[k] = -f_3 * fd_s_11[k]
                  + f_2 * fg_s_59[k]
                  + f_4 * fd_8[k]
                  + pb_x[k] * ff_20[k];

        t_55[k] = -f_3 * fd_s_12[k]
                  + f_2 * fg_s_60[k]
                  + f_4 * fd_9[k]
                  + pb_x[k] * ff_21[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pb_y, df_14, fd_s_11, fg_s_61, fg_s_62, \
                         fg_s_63, fg_s_64, fd_8, ff_22, ff_24, ff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_2 * fg_s_61[k]
                  + pb_x[k] * ff_22[k];

        t_57[k] = f_2 * fg_s_62[k]
                  + pb_x[k] * ff_24[k];

        t_58[k] = f_2 * fg_s_63[k]
                  + pb_x[k] * ff_25[k];

        t_59[k] = f_0 * df_14[k]
                  - f_1 * fd_s_11[k]
                  + f_2 * fg_s_64[k]
                  + f_0 * fd_8[k]
                  + pb_y[k] * ff_22[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_y, pb_z, df_17, fd_s_11, fg_s_65, fg_s_66, \
                         fg_s_67, fd_8, ff_22, ff_23, ff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * fg_s_65[k]
                  + pb_z[k] * ff_22[k];

        t_61[k] = -f_3 * fd_s_11[k]
                  + f_2 * fg_s_66[k]
                  + f_4 * fd_8[k]
                  + pb_z[k] * ff_23[k];

        t_62[k] = f_0 * df_17[k]
                  + f_2 * fg_s_67[k]
                  + pb_y[k] * ff_25[k];
    }

#pragma omp simd aligned(t_63, t_64, pb_x, pb_z, fd_s_12, fd_s_13, fg_s_68, fg_s_69, fd_9, \
                         fd_10, ff_25, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -f_1 * fd_s_12[k]
                  + f_2 * fg_s_68[k]
                  + f_0 * fd_9[k]
                  + pb_z[k] * ff_25[k];

        t_64[k] = -f_7 * fd_s_13[k]
                  + f_2 * fg_s_69[k]
                  + f_5 * fd_10[k]
                  + pb_x[k] * ff_26[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_z, pb_x, dg_28, fd_s_16, fg_s_71, fg_s_74, \
                         fg_s_75, fd_11, ff_27, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -f_3 * fd_s_16[k]
                  + f_2 * fg_s_71[k]
                  + f_4 * fd_11[k]
                  + pb_x[k] * ff_27[k];

        t_66[k] = f_2 * fg_s_74[k]
                  + pb_x[k] * ff_29[k];

        t_67[k] = pa_z[k] * dg_28[k]
                  + f_2 * fg_s_75[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_z, pb_y, pb_z, df_14, df_15, df_19, dg_30, \
                         fg_s_76, fg_s_77, fg_s_78, ff_28, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * df_14[k]
                  + f_2 * fg_s_76[k]
                  + pb_z[k] * ff_28[k];

        t_69[k] = f_5 * df_15[k]
                  + pa_z[k] * dg_30[k]
                  + f_2 * fg_s_77[k];

        t_70[k] = f_5 * df_19[k]
                  + f_2 * fg_s_78[k]
                  + pb_y[k] * ff_29[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_y, pb_z, pg_s_15, pg_15, df_18, df_24, dg_37, \
                         dg_47, fg_s_79, fg_s_86, fg_s_87, ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -f_6 * pg_s_15[k]
                  + f_4 * pg_15[k]
                  + pa_y[k] * dg_37[k]
                  + f_2 * fg_s_79[k];

        t_72[k] = f_8 * df_24[k]
                  + pa_y[k] * dg_47[k]
                  + f_2 * fg_s_86[k];

        t_73[k] = f_5 * df_18[k]
                  + f_2 * fg_s_87[k]
                  + pb_z[k] * ff_30[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_y, pb_y, df_26, df_27, dg_49, dg_51, fg_s_88, \
                         fg_s_89, fg_s_90, ff_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_5 * df_26[k]
                  + pa_y[k] * dg_49[k]
                  + f_2 * fg_s_88[k];

        t_75[k] = f_4 * df_27[k]
                  + f_2 * fg_s_89[k]
                  + pb_y[k] * ff_31[k];

        t_76[k] = pa_y[k] * dg_51[k]
                  + f_2 * fg_s_90[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_x, fd_s_21, fd_s_22, fd_s_23, fg_s_91, fg_s_92, \
                         fg_s_93, fd_12, fd_13, fd_14, ff_32, ff_33, \
                         ff_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -f_1 * fd_s_21[k]
                  + f_2 * fg_s_91[k]
                  + f_0 * fd_12[k]
                  + pb_x[k] * ff_32[k];

        t_78[k] = -f_7 * fd_s_22[k]
                  + f_2 * fg_s_92[k]
                  + f_5 * fd_13[k]
                  + pb_x[k] * ff_33[k];

        t_79[k] = -f_3 * fd_s_23[k]
                  + f_2 * fg_s_93[k]
                  + f_4 * fd_14[k]
                  + pb_x[k] * ff_34[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pb_x, fd_s_25, fg_s_94, fg_s_95, fg_s_96, \
                         fg_s_97, fd_16, ff_35, ff_36, ff_37, ff_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -f_3 * fd_s_25[k]
                  + f_2 * fg_s_94[k]
                  + f_4 * fd_16[k]
                  + pb_x[k] * ff_35[k];

        t_81[k] = f_2 * fg_s_95[k]
                  + pb_x[k] * ff_36[k];

        t_82[k] = f_2 * fg_s_96[k]
                  + pb_x[k] * ff_37[k];

        t_83[k] = f_2 * fg_s_97[k]
                  + pb_x[k] * ff_39[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_y, fd_s_23, fd_s_24, fd_s_25, fg_s_98, fg_s_99, \
                         fg_s_100, fd_14, fd_15, fd_16, ff_36, ff_37, \
                         ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -f_1 * fd_s_23[k]
                  + f_2 * fg_s_98[k]
                  + f_0 * fd_14[k]
                  + pb_y[k] * ff_36[k];

        t_85[k] = -f_7 * fd_s_24[k]
                  + f_2 * fg_s_99[k]
                  + f_5 * fd_15[k]
                  + pb_y[k] * ff_37[k];

        t_86[k] = -f_3 * fd_s_25[k]
                  + f_2 * fg_s_100[k]
                  + f_4 * fd_16[k]
                  + pb_y[k] * ff_38[k];
    }

#pragma omp simd aligned(t_87, t_88, pb_y, pb_z, df_27, fd_s_25, fg_s_101, fg_s_102, fd_16, \
                         ff_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_2 * fg_s_101[k]
                  + pb_y[k] * ff_39[k];

        t_88[k] = f_0 * df_27[k]
                  - f_1 * fd_s_25[k]
                  + f_2 * fg_s_102[k]
                  + f_0 * fd_16[k]
                  + pb_z[k] * ff_39[k];
    }
}

auto
compute_prim_fg_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_s_5 = buffer.data(pg_s + 5);
    const auto *pg_s_14 = buffer.data(pg_s + 14);

    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_9 = buffer.data(pg + 9);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_24 = buffer.data(df + 24);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_39 = buffer.data(dg + 39);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_23 = buffer.data(fg_s + 23);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_34 = buffer.data(fg_s + 34);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_36 = buffer.data(fg_s + 36);
    const auto *fg_s_37 = buffer.data(fg_s + 37);
    const auto *fg_s_38 = buffer.data(fg_s + 38);
    const auto *fg_s_39 = buffer.data(fg_s + 39);
    const auto *fg_s_40 = buffer.data(fg_s + 40);
    const auto *fg_s_41 = buffer.data(fg_s + 41);
    const auto *fg_s_42 = buffer.data(fg_s + 42);
    const auto *fg_s_43 = buffer.data(fg_s + 43);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_45 = buffer.data(fg_s + 45);
    const auto *fg_s_46 = buffer.data(fg_s + 46);
    const auto *fg_s_47 = buffer.data(fg_s + 47);
    const auto *fg_s_49 = buffer.data(fg_s + 49);
    const auto *fg_s_52 = buffer.data(fg_s + 52);
    const auto *fg_s_53 = buffer.data(fg_s + 53);
    const auto *fg_s_57 = buffer.data(fg_s + 57);
    const auto *fg_s_64 = buffer.data(fg_s + 64);
    const auto *fg_s_68 = buffer.data(fg_s + 68);
    const auto *fg_s_69 = buffer.data(fg_s + 69);
    const auto *fg_s_70 = buffer.data(fg_s + 70);
    const auto *fg_s_71 = buffer.data(fg_s + 71);
    const auto *fg_s_72 = buffer.data(fg_s + 72);
    const auto *fg_s_73 = buffer.data(fg_s + 73);
    const auto *fg_s_74 = buffer.data(fg_s + 74);
    const auto *fg_s_75 = buffer.data(fg_s + 75);
    const auto *fg_s_76 = buffer.data(fg_s + 76);
    const auto *fg_s_77 = buffer.data(fg_s + 77);
    const auto *fg_s_78 = buffer.data(fg_s + 78);
    const auto *fg_s_79 = buffer.data(fg_s + 79);
    const auto *fg_s_80 = buffer.data(fg_s + 80);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fg_s_3, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_2 * fg_s_1[k]
                 + pb_y[k] * ff_0[k];

        t_2[k] = f_2 * fg_s_2[k]
                 + pb_z[k] * ff_0[k];

        t_3[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_3[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, fd_s_0, fd_s_1, fg_s_4, fg_s_5, fd_0, fd_1, \
                         ff_2, ff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_4[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_5[k] = -f_1 * fd_s_1[k]
                 + f_2 * fg_s_5[k]
                 + f_0 * fd_1[k]
                 + pb_y[k] * ff_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd_s_2, fg_s_6, fg_s_7, \
                         fg_s_8, fg_s_9, fd_2, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_3 * fd_s_2[k]
                 + f_2 * fg_s_6[k]
                 + f_4 * fd_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = f_2 * fg_s_7[k]
                 + pb_y[k] * ff_5[k];

        t_8[k] = -f_1 * fd_s_2[k]
                 + f_2 * fg_s_8[k]
                 + f_0 * fd_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_z, pg_s_5, pg_3, dg_10, fd_s_3, fg_s_11, \
                         fg_s_12, fg_s_13, fd_3, ff_6, ff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_5 * pg_s_5[k]
                  + f_4 * pg_3[k]
                  + pa_x[k] * dg_10[k]
                  + f_2 * fg_s_11[k];

        t_11[k] = f_2 * fg_s_12[k]
                  + pb_z[k] * ff_6[k];

        t_12[k] = -f_3 * fd_s_3[k]
                  + f_2 * fg_s_13[k]
                  + f_4 * fd_3[k]
                  + pb_z[k] * ff_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_y, dg_0, dg_8, fd_s_5, fg_s_14, \
                         fg_s_15, fg_s_19, fd_4, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * dg_8[k]
                  + f_2 * fg_s_14[k];

        t_14[k] = pa_z[k] * dg_0[k]
                  + f_2 * fg_s_15[k];

        t_15[k] = -f_6 * fd_s_5[k]
                  + f_2 * fg_s_19[k]
                  + f_7 * fd_4[k]
                  + pb_y[k] * ff_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_y, pg_s_14, pg_9, dg_12, fd_s_6, fg_s_20, \
                         fg_s_21, fg_s_22, fd_5, ff_9, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_3 * fd_s_6[k]
                  + f_2 * fg_s_20[k]
                  + f_4 * fd_5[k]
                  + pb_y[k] * ff_9[k];

        t_17[k] = f_2 * fg_s_21[k]
                  + pb_y[k] * ff_10[k];

        t_18[k] = -f_5 * pg_s_14[k]
                  + f_4 * pg_9[k]
                  + pa_x[k] * dg_12[k]
                  + f_2 * fg_s_22[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, df_8, df_17, dg_13, dg_20, dg_28, \
                         dg_39, fg_s_23, fg_s_27, fg_s_28, fg_s_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_8 * df_8[k]
                  + pa_x[k] * dg_13[k]
                  + f_2 * fg_s_23[k];

        t_20[k] = pa_x[k] * dg_20[k]
                  + f_2 * fg_s_27[k];

        t_21[k] = f_8 * df_17[k]
                  + pa_x[k] * dg_28[k]
                  + f_2 * fg_s_28[k];

        t_22[k] = pa_x[k] * dg_39[k]
                  + f_2 * fg_s_34[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_x, fd_s_9, fd_s_10, fd_s_11, fg_s_35, fg_s_36, \
                         fg_s_37, fd_6, fd_7, fd_8, ff_11, ff_12, \
                         ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_1 * fd_s_9[k]
                  + f_2 * fg_s_35[k]
                  + f_0 * fd_6[k]
                  + pb_x[k] * ff_11[k];

        t_24[k] = -f_6 * fd_s_10[k]
                  + f_2 * fg_s_36[k]
                  + f_7 * fd_7[k]
                  + pb_x[k] * ff_12[k];

        t_25[k] = -f_3 * fd_s_11[k]
                  + f_2 * fg_s_37[k]
                  + f_4 * fd_8[k]
                  + pb_x[k] * ff_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, fd_s_12, fg_s_38, fg_s_39, fg_s_40, \
                         fg_s_41, fd_9, ff_14, ff_15, ff_17, ff_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -f_3 * fd_s_12[k]
                  + f_2 * fg_s_38[k]
                  + f_4 * fd_9[k]
                  + pb_x[k] * ff_14[k];

        t_27[k] = f_2 * fg_s_39[k]
                  + pb_x[k] * ff_15[k];

        t_28[k] = f_2 * fg_s_40[k]
                  + pb_x[k] * ff_17[k];

        t_29[k] = f_2 * fg_s_41[k]
                  + pb_x[k] * ff_18[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, pb_z, df_12, fd_s_11, fg_s_42, fg_s_43, \
                         fg_s_44, fd_8, ff_15, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * df_12[k]
                  - f_1 * fd_s_11[k]
                  + f_2 * fg_s_42[k]
                  + f_0 * fd_8[k]
                  + pb_y[k] * ff_15[k];

        t_31[k] = f_2 * fg_s_43[k]
                  + pb_z[k] * ff_15[k];

        t_32[k] = -f_3 * fd_s_11[k]
                  + f_2 * fg_s_44[k]
                  + f_4 * fd_8[k]
                  + pb_z[k] * ff_16[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_y, pb_z, df_15, fd_s_12, fd_s_13, fg_s_45, \
                         fg_s_46, fg_s_47, fd_9, fd_10, ff_18, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * df_15[k]
                  + f_2 * fg_s_45[k]
                  + pb_y[k] * ff_18[k];

        t_34[k] = -f_1 * fd_s_12[k]
                  + f_2 * fg_s_46[k]
                  + f_0 * fd_9[k]
                  + pb_z[k] * ff_18[k];

        t_35[k] = -f_6 * fd_s_13[k]
                  + f_2 * fg_s_47[k]
                  + f_7 * fd_10[k]
                  + pb_x[k] * ff_19[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_z, pb_x, dg_20, fd_s_16, fg_s_49, fg_s_52, \
                         fg_s_53, fd_11, ff_20, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_3 * fd_s_16[k]
                  + f_2 * fg_s_49[k]
                  + f_4 * fd_11[k]
                  + pb_x[k] * ff_20[k];

        t_37[k] = f_2 * fg_s_52[k]
                  + pb_x[k] * ff_21[k];

        t_38[k] = pa_z[k] * dg_20[k]
                  + f_2 * fg_s_53[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pg_s_14, pg_9, df_21, dg_27, dg_35, dg_39, \
                         fg_s_57, fg_s_64, fg_s_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_5 * pg_s_14[k]
                  + f_4 * pg_9[k]
                  + pa_y[k] * dg_27[k]
                  + f_2 * fg_s_57[k];

        t_40[k] = f_8 * df_21[k]
                  + pa_y[k] * dg_35[k]
                  + f_2 * fg_s_64[k];

        t_41[k] = pa_y[k] * dg_39[k]
                  + f_2 * fg_s_68[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, fd_s_21, fd_s_22, fd_s_23, fg_s_69, fg_s_70, \
                         fg_s_71, fd_12, fd_13, fd_14, ff_22, ff_23, \
                         ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -f_1 * fd_s_21[k]
                  + f_2 * fg_s_69[k]
                  + f_0 * fd_12[k]
                  + pb_x[k] * ff_22[k];

        t_43[k] = -f_6 * fd_s_22[k]
                  + f_2 * fg_s_70[k]
                  + f_7 * fd_13[k]
                  + pb_x[k] * ff_23[k];

        t_44[k] = -f_3 * fd_s_23[k]
                  + f_2 * fg_s_71[k]
                  + f_4 * fd_14[k]
                  + pb_x[k] * ff_24[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_x, fd_s_25, fg_s_72, fg_s_73, fg_s_74, \
                         fg_s_75, fd_16, ff_25, ff_26, ff_27, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -f_3 * fd_s_25[k]
                  + f_2 * fg_s_72[k]
                  + f_4 * fd_16[k]
                  + pb_x[k] * ff_25[k];

        t_46[k] = f_2 * fg_s_73[k]
                  + pb_x[k] * ff_26[k];

        t_47[k] = f_2 * fg_s_74[k]
                  + pb_x[k] * ff_27[k];

        t_48[k] = f_2 * fg_s_75[k]
                  + pb_x[k] * ff_29[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_y, fd_s_23, fd_s_24, fd_s_25, fg_s_76, fg_s_77, \
                         fg_s_78, fd_14, fd_15, fd_16, ff_26, ff_27, \
                         ff_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = -f_1 * fd_s_23[k]
                  + f_2 * fg_s_76[k]
                  + f_0 * fd_14[k]
                  + pb_y[k] * ff_26[k];

        t_50[k] = -f_6 * fd_s_24[k]
                  + f_2 * fg_s_77[k]
                  + f_7 * fd_15[k]
                  + pb_y[k] * ff_27[k];

        t_51[k] = -f_3 * fd_s_25[k]
                  + f_2 * fg_s_78[k]
                  + f_4 * fd_16[k]
                  + pb_y[k] * ff_28[k];
    }

#pragma omp simd aligned(t_52, t_53, pb_y, pb_z, df_24, fd_s_25, fg_s_79, fg_s_80, fd_16, \
                         ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * fg_s_79[k]
                  + pb_y[k] * ff_29[k];

        t_53[k] = f_0 * df_24[k]
                  - f_1 * fd_s_25[k]
                  + f_2 * fg_s_80[k]
                  + f_0 * fd_16[k]
                  + pb_z[k] * ff_29[k];
    }
}

auto
compute_prim_fg_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = 2.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_s_0 = buffer.data(pg_s + 0);
    const auto *pg_s_2 = buffer.data(pg_s + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);
    const auto *fd_s_25 = buffer.data(fd_s + 25);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_23 = buffer.data(fg_s + 23);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_29 = buffer.data(fg_s + 29);
    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_31 = buffer.data(fg_s + 31);
    const auto *fg_s_32 = buffer.data(fg_s + 32);
    const auto *fg_s_33 = buffer.data(fg_s + 33);
    const auto *fg_s_34 = buffer.data(fg_s + 34);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_36 = buffer.data(fg_s + 36);
    const auto *fg_s_37 = buffer.data(fg_s + 37);
    const auto *fg_s_38 = buffer.data(fg_s + 38);
    const auto *fg_s_39 = buffer.data(fg_s + 39);
    const auto *fg_s_40 = buffer.data(fg_s + 40);
    const auto *fg_s_41 = buffer.data(fg_s + 41);
    const auto *fg_s_42 = buffer.data(fg_s + 42);
    const auto *fg_s_43 = buffer.data(fg_s + 43);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_45 = buffer.data(fg_s + 45);
    const auto *fg_s_46 = buffer.data(fg_s + 46);
    const auto *fg_s_47 = buffer.data(fg_s + 47);
    const auto *fg_s_48 = buffer.data(fg_s + 48);
    const auto *fg_s_49 = buffer.data(fg_s + 49);
    const auto *fg_s_50 = buffer.data(fg_s + 50);
    const auto *fg_s_51 = buffer.data(fg_s + 51);
    const auto *fg_s_52 = buffer.data(fg_s + 52);
    const auto *fg_s_53 = buffer.data(fg_s + 53);
    const auto *fg_s_54 = buffer.data(fg_s + 54);
    const auto *fg_s_55 = buffer.data(fg_s + 55);
    const auto *fg_s_56 = buffer.data(fg_s + 56);
    const auto *fg_s_57 = buffer.data(fg_s + 57);
    const auto *fg_s_58 = buffer.data(fg_s + 58);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_44 = buffer.data(ff + 44);
    const auto *ff_45 = buffer.data(ff + 45);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_49 = buffer.data(ff + 49);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_51 = buffer.data(ff + 51);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_53 = buffer.data(ff + 53);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fd_0, ff_0, ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_1[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];

        t_2[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_2[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, df_3, df_4, fd_s_1, fg_s_3, fg_s_4, \
                         fg_s_5, fd_1, ff_3, ff_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * df_3[k]
                 + f_2 * fg_s_3[k]
                 + pb_x[k] * ff_3[k];

        t_4[k] = f_0 * df_4[k]
                 + f_2 * fg_s_4[k]
                 + pb_x[k] * ff_4[k];

        t_5[k] = -f_1 * fd_s_1[k]
                 + f_2 * fg_s_5[k]
                 + f_0 * fd_1[k]
                 + pb_y[k] * ff_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_y, pb_y, pb_z, df_0, dg_0, fd_s_2, fg_s_6, fg_s_7, \
                         fg_s_8, fd_2, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * fd_s_2[k]
                 + f_2 * fg_s_6[k]
                 + f_0 * fd_2[k]
                 + pb_z[k] * ff_4[k];

        t_7[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_7[k];

        t_8[k] = f_4 * df_0[k]
                 + f_2 * fg_s_8[k]
                 + pb_y[k] * ff_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pb_x, pg_s_0, pg_0, df_1, df_6, dg_1, \
                         dg_3, fg_s_9, fg_s_10, fg_s_11, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * df_1[k]
                 + pa_y[k] * dg_1[k]
                 + f_2 * fg_s_9[k];

        t_10[k] = f_5 * df_6[k]
                  + f_2 * fg_s_10[k]
                  + pb_x[k] * ff_6[k];

        t_11[k] = -f_6 * pg_s_0[k]
                  + f_4 * pg_0[k]
                  + pa_x[k] * dg_3[k]
                  + f_2 * fg_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_z, pb_z, df_0, df_2, dg_0, dg_2, fg_s_12, \
                         fg_s_13, fg_s_14, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * dg_0[k]
                  + f_2 * fg_s_12[k];

        t_13[k] = f_4 * df_0[k]
                  + f_2 * fg_s_13[k]
                  + pb_z[k] * ff_9[k];

        t_14[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_2[k]
                  + f_2 * fg_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, pg_s_2, pg_2, df_9, df_10, dg_4, dg_5, \
                         fg_s_15, fg_s_16, fg_s_17, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * df_9[k]
                  + f_2 * fg_s_15[k]
                  + pb_x[k] * ff_12[k];

        t_16[k] = -f_6 * pg_s_2[k]
                  + f_4 * pg_2[k]
                  + pa_x[k] * dg_4[k]
                  + f_2 * fg_s_16[k];

        t_17[k] = f_7 * df_10[k]
                  + pa_x[k] * dg_5[k]
                  + f_2 * fg_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pb_x, pb_y, df_5, df_11, df_12, dg_6, \
                         fg_s_18, fg_s_19, fg_s_20, ff_13, ff_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * df_5[k]
                  + f_2 * fg_s_18[k]
                  + pb_y[k] * ff_13[k];

        t_19[k] = f_5 * df_11[k]
                  + pa_x[k] * dg_6[k]
                  + f_2 * fg_s_19[k];

        t_20[k] = f_4 * df_12[k]
                  + f_2 * fg_s_20[k]
                  + pb_x[k] * ff_15[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pb_z, df_7, df_20, dg_7, dg_9, dg_11, \
                         fg_s_21, fg_s_22, fg_s_23, fg_s_24, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_x[k] * dg_7[k]
                  + f_2 * fg_s_21[k];

        t_22[k] = pa_x[k] * dg_9[k]
                  + f_2 * fg_s_22[k];

        t_23[k] = f_7 * df_20[k]
                  + pa_x[k] * dg_11[k]
                  + f_2 * fg_s_23[k];

        t_24[k] = f_5 * df_7[k]
                  + f_2 * fg_s_24[k]
                  + pb_z[k] * ff_21[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pb_x, df_23, df_27, dg_13, dg_16, fg_s_25, \
                         fg_s_26, fg_s_27, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_5 * df_23[k]
                  + pa_x[k] * dg_13[k]
                  + f_2 * fg_s_25[k];

        t_26[k] = f_4 * df_27[k]
                  + f_2 * fg_s_26[k]
                  + pb_x[k] * ff_26[k];

        t_27[k] = pa_x[k] * dg_16[k]
                  + f_2 * fg_s_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, fd_s_13, fd_s_14, fd_s_15, fg_s_28, fg_s_29, \
                         fg_s_30, fd_13, fd_14, fd_15, ff_27, ff_28, \
                         ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_1 * fd_s_13[k]
                  + f_2 * fg_s_28[k]
                  + f_0 * fd_13[k]
                  + pb_x[k] * ff_27[k];

        t_29[k] = -f_8 * fd_s_14[k]
                  + f_2 * fg_s_29[k]
                  + f_5 * fd_14[k]
                  + pb_x[k] * ff_28[k];

        t_30[k] = -f_3 * fd_s_15[k]
                  + f_2 * fg_s_30[k]
                  + f_4 * fd_15[k]
                  + pb_x[k] * ff_29[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_x, pb_z, fd_s_16, fg_s_31, fg_s_32, fd_16, ff_28, \
                         ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_2 * fg_s_31[k]
                  + pb_z[k] * ff_28[k];

        t_32[k] = -f_3 * fd_s_16[k]
                  + f_2 * fg_s_32[k]
                  + f_4 * fd_16[k]
                  + pb_x[k] * ff_30[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_y, pb_z, df_12, df_15, fd_s_15, fg_s_33, \
                         fg_s_34, fg_s_35, fd_15, ff_31, ff_32, ff_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * df_12[k]
                  - f_1 * fd_s_15[k]
                  + f_2 * fg_s_33[k]
                  + f_0 * fd_15[k]
                  + pb_y[k] * ff_31[k];

        t_34[k] = -f_3 * fd_s_15[k]
                  + f_2 * fg_s_34[k]
                  + f_4 * fd_15[k]
                  + pb_z[k] * ff_32[k];

        t_35[k] = f_0 * df_15[k]
                  + f_2 * fg_s_35[k]
                  + pb_y[k] * ff_34[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_z, pb_x, pb_z, dg_7, fd_s_16, fd_s_18, fg_s_36, \
                         fg_s_37, fg_s_38, fd_16, fd_18, ff_34, ff_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_1 * fd_s_16[k]
                  + f_2 * fg_s_36[k]
                  + f_0 * fd_16[k]
                  + pb_z[k] * ff_34[k];

        t_37[k] = -f_3 * fd_s_18[k]
                  + f_2 * fg_s_37[k]
                  + f_4 * fd_18[k]
                  + pb_x[k] * ff_35[k];

        t_38[k] = pa_z[k] * dg_7[k]
                  + f_2 * fg_s_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_z, pb_y, pb_z, df_12, df_13, df_19, dg_8, \
                         fg_s_39, fg_s_40, fg_s_41, ff_36, ff_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_4 * df_12[k]
                  + f_2 * fg_s_39[k]
                  + pb_z[k] * ff_36[k];

        t_40[k] = f_5 * df_13[k]
                  + pa_z[k] * dg_8[k]
                  + f_2 * fg_s_40[k];

        t_41[k] = f_5 * df_19[k]
                  + f_2 * fg_s_41[k]
                  + pb_y[k] * ff_39[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_y, pg_s_2, pg_2, df_21, df_24, dg_10, dg_12, \
                         dg_14, fg_s_42, fg_s_43, fg_s_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = -f_6 * pg_s_2[k]
                  + f_4 * pg_2[k]
                  + pa_y[k] * dg_10[k]
                  + f_2 * fg_s_42[k];

        t_43[k] = f_5 * df_21[k]
                  + pa_y[k] * dg_12[k]
                  + f_2 * fg_s_43[k];

        t_44[k] = f_7 * df_24[k]
                  + pa_y[k] * dg_14[k]
                  + f_2 * fg_s_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, df_16, df_26, df_27, dg_15, \
                         fg_s_45, fg_s_46, fg_s_47, ff_41, ff_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_5 * df_16[k]
                  + f_2 * fg_s_45[k]
                  + pb_z[k] * ff_41[k];

        t_46[k] = f_5 * df_26[k]
                  + pa_y[k] * dg_15[k]
                  + f_2 * fg_s_46[k];

        t_47[k] = f_4 * df_27[k]
                  + f_2 * fg_s_47[k]
                  + pb_y[k] * ff_44[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pa_y, pb_x, pb_y, dg_16, fd_s_21, fg_s_48, fg_s_49, \
                         fg_s_50, fd_21, ff_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * dg_16[k]
                  + f_2 * fg_s_48[k];

        t_49[k] = -f_1 * fd_s_21[k]
                  + f_2 * fg_s_49[k]
                  + f_0 * fd_21[k]
                  + pb_x[k] * ff_45[k];

        t_50[k] = f_2 * fg_s_50[k]
                  + pb_y[k] * ff_45[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_x, pb_y, fd_s_22, fd_s_23, fg_s_51, fg_s_52, \
                         fg_s_53, fd_22, fd_23, ff_47, ff_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -f_8 * fd_s_22[k]
                  + f_2 * fg_s_51[k]
                  + f_5 * fd_22[k]
                  + pb_x[k] * ff_47[k];

        t_52[k] = -f_3 * fd_s_23[k]
                  + f_2 * fg_s_52[k]
                  + f_4 * fd_23[k]
                  + pb_x[k] * ff_48[k];

        t_53[k] = f_2 * fg_s_53[k]
                  + pb_y[k] * ff_47[k];
    }

#pragma omp simd aligned(t_54, t_55, pb_x, pb_y, fd_s_23, fd_s_25, fg_s_54, fg_s_55, fd_23, \
                         fd_25, ff_49, ff_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = -f_3 * fd_s_25[k]
                  + f_2 * fg_s_54[k]
                  + f_4 * fd_25[k]
                  + pb_x[k] * ff_49[k];

        t_55[k] = -f_1 * fd_s_23[k]
                  + f_2 * fg_s_55[k]
                  + f_0 * fd_23[k]
                  + pb_y[k] * ff_50[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pb_y, pb_z, df_27, fd_s_24, fd_s_25, fg_s_56, \
                         fg_s_57, fg_s_58, fd_24, fd_25, ff_51, ff_52, \
                         ff_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -f_8 * fd_s_24[k]
                  + f_2 * fg_s_56[k]
                  + f_5 * fd_24[k]
                  + pb_y[k] * ff_51[k];

        t_57[k] = -f_3 * fd_s_25[k]
                  + f_2 * fg_s_57[k]
                  + f_4 * fd_25[k]
                  + pb_y[k] * ff_52[k];

        t_58[k] = f_0 * df_27[k]
                  - f_1 * fd_s_25[k]
                  + f_2 * fg_s_58[k]
                  + f_0 * fd_25[k]
                  + pb_z[k] * ff_53[k];
    }
}

auto
compute_prim_fg_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * alpha / p;
    const auto f_8 = 2.0 / p;

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

    const auto *pg_s_2 = buffer.data(pg_s + 2);
    const auto *pg_s_6 = buffer.data(pg_s + 6);

    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_6 = buffer.data(pg + 6);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_38 = buffer.data(dg + 38);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_23 = buffer.data(fg_s + 23);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_29 = buffer.data(fg_s + 29);
    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_31 = buffer.data(fg_s + 31);
    const auto *fg_s_32 = buffer.data(fg_s + 32);
    const auto *fg_s_33 = buffer.data(fg_s + 33);
    const auto *fg_s_34 = buffer.data(fg_s + 34);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_36 = buffer.data(fg_s + 36);
    const auto *fg_s_37 = buffer.data(fg_s + 37);
    const auto *fg_s_38 = buffer.data(fg_s + 38);
    const auto *fg_s_39 = buffer.data(fg_s + 39);
    const auto *fg_s_40 = buffer.data(fg_s + 40);
    const auto *fg_s_41 = buffer.data(fg_s + 41);
    const auto *fg_s_42 = buffer.data(fg_s + 42);
    const auto *fg_s_43 = buffer.data(fg_s + 43);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_45 = buffer.data(fg_s + 45);
    const auto *fg_s_46 = buffer.data(fg_s + 46);
    const auto *fg_s_47 = buffer.data(fg_s + 47);
    const auto *fg_s_48 = buffer.data(fg_s + 48);
    const auto *fg_s_49 = buffer.data(fg_s + 49);
    const auto *fg_s_50 = buffer.data(fg_s + 50);
    const auto *fg_s_51 = buffer.data(fg_s + 51);
    const auto *fg_s_52 = buffer.data(fg_s + 52);
    const auto *fg_s_53 = buffer.data(fg_s + 53);
    const auto *fg_s_54 = buffer.data(fg_s + 54);
    const auto *fg_s_55 = buffer.data(fg_s + 55);
    const auto *fg_s_56 = buffer.data(fg_s + 56);
    const auto *fg_s_57 = buffer.data(fg_s + 57);
    const auto *fg_s_58 = buffer.data(fg_s + 58);
    const auto *fg_s_59 = buffer.data(fg_s + 59);
    const auto *fg_s_60 = buffer.data(fg_s + 60);
    const auto *fg_s_61 = buffer.data(fg_s + 61);
    const auto *fg_s_62 = buffer.data(fg_s + 62);
    const auto *fg_s_63 = buffer.data(fg_s + 63);
    const auto *fg_s_64 = buffer.data(fg_s + 64);
    const auto *fg_s_65 = buffer.data(fg_s + 65);
    const auto *fg_s_66 = buffer.data(fg_s + 66);
    const auto *fg_s_67 = buffer.data(fg_s + 67);
    const auto *fg_s_68 = buffer.data(fg_s + 68);
    const auto *fg_s_69 = buffer.data(fg_s + 69);
    const auto *fg_s_70 = buffer.data(fg_s + 70);
    const auto *fg_s_71 = buffer.data(fg_s + 71);
    const auto *fg_s_72 = buffer.data(fg_s + 72);
    const auto *fg_s_73 = buffer.data(fg_s + 73);
    const auto *fg_s_74 = buffer.data(fg_s + 74);
    const auto *fg_s_75 = buffer.data(fg_s + 75);
    const auto *fg_s_76 = buffer.data(fg_s + 76);
    const auto *fg_s_77 = buffer.data(fg_s + 77);
    const auto *fg_s_78 = buffer.data(fg_s + 78);
    const auto *fg_s_79 = buffer.data(fg_s + 79);
    const auto *fg_s_80 = buffer.data(fg_s + 80);
    const auto *fg_s_81 = buffer.data(fg_s + 81);
    const auto *fg_s_82 = buffer.data(fg_s + 82);
    const auto *fg_s_83 = buffer.data(fg_s + 83);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
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
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_40 = buffer.data(ff + 40);
    const auto *ff_41 = buffer.data(ff + 41);
    const auto *ff_42 = buffer.data(ff + 42);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fg_s_3, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_2 * fg_s_1[k]
                 + pb_y[k] * ff_0[k];

        t_2[k] = f_2 * fg_s_2[k]
                 + pb_z[k] * ff_0[k];

        t_3[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_3[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, fd_s_0, fd_s_1, fg_s_4, fg_s_5, fd_0, fd_1, \
                         ff_2, ff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_4[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_5[k] = -f_1 * fd_s_1[k]
                 + f_2 * fg_s_5[k]
                 + f_0 * fd_1[k]
                 + pb_y[k] * ff_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_y, pb_y, pb_z, dg_0, fd_s_2, fg_s_6, fg_s_7, \
                         fg_s_8, fd_2, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_3 * fd_s_2[k]
                 + f_2 * fg_s_6[k]
                 + f_4 * fd_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = -f_1 * fd_s_2[k]
                 + f_2 * fg_s_7[k]
                 + f_0 * fd_2[k]
                 + pb_z[k] * ff_5[k];

        t_8[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pg_s_2, pg_2, df_1, dg_3, dg_4, dg_9, \
                         fg_s_9, fg_s_10, fg_s_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * df_1[k]
                 + pa_y[k] * dg_3[k]
                 + f_2 * fg_s_9[k];

        t_10[k] = pa_y[k] * dg_4[k]
                  + f_2 * fg_s_10[k];

        t_11[k] = -f_6 * pg_s_2[k]
                  + f_4 * pg_2[k]
                  + pa_x[k] * dg_9[k]
                  + f_2 * fg_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pb_y, pb_z, df_4, dg_6, fd_s_4, fg_s_12, \
                         fg_s_13, fg_s_14, fd_4, ff_8, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_3 * fd_s_4[k]
                  + f_2 * fg_s_12[k]
                  + f_4 * fd_4[k]
                  + pb_z[k] * ff_8[k];

        t_13[k] = f_4 * df_4[k]
                  + f_2 * fg_s_13[k]
                  + pb_y[k] * ff_9[k];

        t_14[k] = pa_y[k] * dg_6[k]
                  + f_2 * fg_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_z, pb_z, df_0, df_2, dg_0, dg_4, fg_s_15, \
                         fg_s_16, fg_s_17, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_z[k] * dg_0[k]
                  + f_2 * fg_s_15[k];

        t_16[k] = f_4 * df_0[k]
                  + f_2 * fg_s_16[k]
                  + pb_z[k] * ff_10[k];

        t_17[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_4[k]
                  + f_2 * fg_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, pb_y, fd_s_6, fd_s_7, fg_s_18, fg_s_19, fd_6, fd_7, \
                         ff_11, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_7 * fd_s_6[k]
                  + f_2 * fg_s_18[k]
                  + f_5 * fd_6[k]
                  + pb_y[k] * ff_11[k];

        t_19[k] = -f_3 * fd_s_7[k]
                  + f_2 * fg_s_19[k]
                  + f_4 * fd_7[k]
                  + pb_y[k] * ff_12[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_x, pg_s_6, pg_6, df_8, df_9, dg_13, dg_14, \
                         dg_15, fg_s_20, fg_s_21, fg_s_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_6 * pg_s_6[k]
                  + f_4 * pg_6[k]
                  + pa_x[k] * dg_13[k]
                  + f_2 * fg_s_20[k];

        t_21[k] = f_8 * df_8[k]
                  + pa_x[k] * dg_14[k]
                  + f_2 * fg_s_21[k];

        t_22[k] = f_5 * df_9[k]
                  + pa_x[k] * dg_15[k]
                  + f_2 * fg_s_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, df_10, df_11, dg_16, dg_18, \
                         dg_20, fg_s_23, fg_s_24, fg_s_25, fg_s_26, \
                         ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * df_10[k]
                  + pa_x[k] * dg_16[k]
                  + f_2 * fg_s_23[k];

        t_24[k] = f_4 * df_11[k]
                  + f_2 * fg_s_24[k]
                  + pb_x[k] * ff_16[k];

        t_25[k] = pa_x[k] * dg_18[k]
                  + f_2 * fg_s_25[k];

        t_26[k] = pa_x[k] * dg_20[k]
                  + f_2 * fg_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pa_y, pa_z, dg_8, dg_11, dg_21, dg_22, \
                         fg_s_27, fg_s_28, fg_s_29, fg_s_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_x[k] * dg_21[k]
                  + f_2 * fg_s_27[k];

        t_28[k] = pa_x[k] * dg_22[k]
                  + f_2 * fg_s_28[k];

        t_29[k] = pa_y[k] * dg_11[k]
                  + f_2 * fg_s_29[k];

        t_30[k] = pa_z[k] * dg_8[k]
                  + f_2 * fg_s_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pa_y, dg_12, dg_24, dg_25, dg_26, \
                         fg_s_31, fg_s_32, fg_s_33, fg_s_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_y[k] * dg_12[k]
                  + f_2 * fg_s_31[k];

        t_32[k] = pa_x[k] * dg_24[k]
                  + f_2 * fg_s_32[k];

        t_33[k] = pa_x[k] * dg_25[k]
                  + f_2 * fg_s_33[k];

        t_34[k] = pa_x[k] * dg_26[k]
                  + f_2 * fg_s_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pb_z, df_6, df_16, df_19, dg_28, dg_31, \
                         fg_s_35, fg_s_36, fg_s_37, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * df_16[k]
                  + pa_x[k] * dg_28[k]
                  + f_2 * fg_s_35[k];

        t_36[k] = f_5 * df_6[k]
                  + f_2 * fg_s_36[k]
                  + pb_z[k] * ff_17[k];

        t_37[k] = f_5 * df_19[k]
                  + pa_x[k] * dg_31[k]
                  + f_2 * fg_s_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pb_x, df_23, dg_34, dg_35, dg_36, \
                         fg_s_38, fg_s_39, fg_s_40, fg_s_41, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_4 * df_23[k]
                  + f_2 * fg_s_38[k]
                  + pb_x[k] * ff_19[k];

        t_39[k] = pa_x[k] * dg_34[k]
                  + f_2 * fg_s_39[k];

        t_40[k] = pa_x[k] * dg_35[k]
                  + f_2 * fg_s_40[k];

        t_41[k] = pa_x[k] * dg_36[k]
                  + f_2 * fg_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_x, pb_x, dg_38, fd_s_12, fd_s_13, fg_s_42, \
                         fg_s_43, fg_s_44, fd_10, fd_11, ff_20, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_x[k] * dg_38[k]
                  + f_2 * fg_s_42[k];

        t_43[k] = -f_1 * fd_s_12[k]
                  + f_2 * fg_s_43[k]
                  + f_0 * fd_10[k]
                  + pb_x[k] * ff_20[k];

        t_44[k] = -f_7 * fd_s_13[k]
                  + f_2 * fg_s_44[k]
                  + f_5 * fd_11[k]
                  + pb_x[k] * ff_21[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, pb_z, fd_s_14, fd_s_15, fg_s_45, fg_s_46, \
                         fg_s_47, fd_12, fd_13, ff_21, ff_22, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -f_3 * fd_s_14[k]
                  + f_2 * fg_s_45[k]
                  + f_4 * fd_12[k]
                  + pb_x[k] * ff_22[k];

        t_46[k] = f_2 * fg_s_46[k]
                  + pb_z[k] * ff_21[k];

        t_47[k] = -f_3 * fd_s_15[k]
                  + f_2 * fg_s_47[k]
                  + f_4 * fd_13[k]
                  + pb_x[k] * ff_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_y, df_11, fd_s_14, fg_s_48, fg_s_49, \
                         fg_s_50, fg_s_51, fd_12, ff_24, ff_26, ff_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_2 * fg_s_48[k]
                  + pb_x[k] * ff_24[k];

        t_49[k] = f_2 * fg_s_49[k]
                  + pb_x[k] * ff_26[k];

        t_50[k] = f_2 * fg_s_50[k]
                  + pb_x[k] * ff_27[k];

        t_51[k] = f_0 * df_11[k]
                  - f_1 * fd_s_14[k]
                  + f_2 * fg_s_51[k]
                  + f_0 * fd_12[k]
                  + pb_y[k] * ff_24[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_y, pb_z, df_13, fd_s_14, fg_s_52, fg_s_53, \
                         fg_s_54, fd_12, ff_24, ff_25, ff_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * fg_s_52[k]
                  + pb_z[k] * ff_24[k];

        t_53[k] = -f_3 * fd_s_14[k]
                  + f_2 * fg_s_53[k]
                  + f_4 * fd_12[k]
                  + pb_z[k] * ff_25[k];

        t_54[k] = f_0 * df_13[k]
                  + f_2 * fg_s_54[k]
                  + pb_y[k] * ff_27[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pb_z, fd_s_15, fd_s_17, fg_s_55, fg_s_56, \
                         fg_s_57, fd_13, fd_15, ff_27, ff_28, ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -f_1 * fd_s_15[k]
                  + f_2 * fg_s_55[k]
                  + f_0 * fd_13[k]
                  + pb_z[k] * ff_27[k];

        t_56[k] = -f_3 * fd_s_17[k]
                  + f_2 * fg_s_56[k]
                  + f_4 * fd_15[k]
                  + pb_x[k] * ff_28[k];

        t_57[k] = f_2 * fg_s_57[k]
                  + pb_x[k] * ff_30[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_z, pb_z, df_11, df_12, dg_18, dg_20, fg_s_58, \
                         fg_s_59, fg_s_60, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pa_z[k] * dg_18[k]
                  + f_2 * fg_s_58[k];

        t_59[k] = f_4 * df_11[k]
                  + f_2 * fg_s_59[k]
                  + pb_z[k] * ff_29[k];

        t_60[k] = f_5 * df_12[k]
                  + pa_z[k] * dg_20[k]
                  + f_2 * fg_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_y, pb_y, pg_s_6, pg_6, df_15, df_17, dg_27, \
                         dg_30, fg_s_61, fg_s_62, fg_s_63, ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * df_15[k]
                  + f_2 * fg_s_61[k]
                  + pb_y[k] * ff_30[k];

        t_62[k] = -f_6 * pg_s_6[k]
                  + f_4 * pg_6[k]
                  + pa_y[k] * dg_27[k]
                  + f_2 * fg_s_62[k];

        t_63[k] = f_5 * df_17[k]
                  + pa_y[k] * dg_30[k]
                  + f_2 * fg_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_y, pb_x, pb_z, df_14, df_20, dg_34, fg_s_64, \
                         fg_s_65, fg_s_66, ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * fg_s_64[k]
                  + pb_x[k] * ff_32[k];

        t_65[k] = f_8 * df_20[k]
                  + pa_y[k] * dg_34[k]
                  + f_2 * fg_s_65[k];

        t_66[k] = f_5 * df_14[k]
                  + f_2 * fg_s_66[k]
                  + pb_z[k] * ff_32[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_y, pb_y, df_22, df_23, dg_36, dg_38, fg_s_67, \
                         fg_s_68, fg_s_69, ff_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_5 * df_22[k]
                  + pa_y[k] * dg_36[k]
                  + f_2 * fg_s_67[k];

        t_68[k] = f_4 * df_23[k]
                  + f_2 * fg_s_68[k]
                  + pb_y[k] * ff_34[k];

        t_69[k] = pa_y[k] * dg_38[k]
                  + f_2 * fg_s_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_x, pb_y, fd_s_20, fd_s_21, fg_s_70, fg_s_71, \
                         fg_s_72, fd_18, fd_19, ff_35, ff_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -f_1 * fd_s_20[k]
                  + f_2 * fg_s_70[k]
                  + f_0 * fd_18[k]
                  + pb_x[k] * ff_35[k];

        t_71[k] = f_2 * fg_s_71[k]
                  + pb_y[k] * ff_35[k];

        t_72[k] = -f_7 * fd_s_21[k]
                  + f_2 * fg_s_72[k]
                  + f_5 * fd_19[k]
                  + pb_x[k] * ff_36[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pb_x, pb_y, fd_s_22, fd_s_24, fg_s_73, fg_s_74, \
                         fg_s_75, fd_20, fd_22, ff_36, ff_37, ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = -f_3 * fd_s_22[k]
                  + f_2 * fg_s_73[k]
                  + f_4 * fd_20[k]
                  + pb_x[k] * ff_37[k];

        t_74[k] = f_2 * fg_s_74[k]
                  + pb_y[k] * ff_36[k];

        t_75[k] = -f_3 * fd_s_24[k]
                  + f_2 * fg_s_75[k]
                  + f_4 * fd_22[k]
                  + pb_x[k] * ff_38[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_x, pb_y, fd_s_22, fg_s_76, fg_s_77, \
                         fg_s_78, fg_s_79, fd_20, ff_39, ff_40, ff_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * fg_s_76[k]
                  + pb_x[k] * ff_39[k];

        t_77[k] = f_2 * fg_s_77[k]
                  + pb_x[k] * ff_40[k];

        t_78[k] = f_2 * fg_s_78[k]
                  + pb_x[k] * ff_42[k];

        t_79[k] = -f_1 * fd_s_22[k]
                  + f_2 * fg_s_79[k]
                  + f_0 * fd_20[k]
                  + pb_y[k] * ff_39[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pb_y, fd_s_23, fd_s_24, fg_s_80, fg_s_81, fg_s_82, \
                         fd_21, fd_22, ff_40, ff_41, ff_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -f_7 * fd_s_23[k]
                  + f_2 * fg_s_80[k]
                  + f_5 * fd_21[k]
                  + pb_y[k] * ff_40[k];

        t_81[k] = -f_3 * fd_s_24[k]
                  + f_2 * fg_s_81[k]
                  + f_4 * fd_22[k]
                  + pb_y[k] * ff_41[k];

        t_82[k] = f_2 * fg_s_82[k]
                  + pb_y[k] * ff_42[k];
    }

#pragma omp simd aligned(t_83, pb_z, df_23, fd_s_24, fg_s_83, fd_22, \
                         ff_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * df_23[k]
                  - f_1 * fd_s_24[k]
                  + f_2 * fg_s_83[k]
                  + f_0 * fd_22[k]
                  + pb_z[k] * ff_42[k];
    }
}

auto
compute_prim_fg_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 * alpha / p;
    const auto f_8 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_s_3 = buffer.data(pg_s + 3);
    const auto *pg_s_9 = buffer.data(pg_s + 9);

    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_9 = buffer.data(pg + 9);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_31 = buffer.data(dg + 31);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_21 = buffer.data(fd_s + 21);
    const auto *fd_s_22 = buffer.data(fd_s + 22);
    const auto *fd_s_23 = buffer.data(fd_s + 23);
    const auto *fd_s_24 = buffer.data(fd_s + 24);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_29 = buffer.data(fg_s + 29);
    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_31 = buffer.data(fg_s + 31);
    const auto *fg_s_32 = buffer.data(fg_s + 32);
    const auto *fg_s_33 = buffer.data(fg_s + 33);
    const auto *fg_s_34 = buffer.data(fg_s + 34);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_36 = buffer.data(fg_s + 36);
    const auto *fg_s_37 = buffer.data(fg_s + 37);
    const auto *fg_s_38 = buffer.data(fg_s + 38);
    const auto *fg_s_39 = buffer.data(fg_s + 39);
    const auto *fg_s_40 = buffer.data(fg_s + 40);
    const auto *fg_s_41 = buffer.data(fg_s + 41);
    const auto *fg_s_42 = buffer.data(fg_s + 42);
    const auto *fg_s_43 = buffer.data(fg_s + 43);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_45 = buffer.data(fg_s + 45);
    const auto *fg_s_46 = buffer.data(fg_s + 46);
    const auto *fg_s_47 = buffer.data(fg_s + 47);
    const auto *fg_s_48 = buffer.data(fg_s + 48);
    const auto *fg_s_49 = buffer.data(fg_s + 49);
    const auto *fg_s_50 = buffer.data(fg_s + 50);
    const auto *fg_s_51 = buffer.data(fg_s + 51);
    const auto *fg_s_52 = buffer.data(fg_s + 52);
    const auto *fg_s_53 = buffer.data(fg_s + 53);
    const auto *fg_s_54 = buffer.data(fg_s + 54);
    const auto *fg_s_55 = buffer.data(fg_s + 55);
    const auto *fg_s_56 = buffer.data(fg_s + 56);
    const auto *fg_s_57 = buffer.data(fg_s + 57);
    const auto *fg_s_58 = buffer.data(fg_s + 58);
    const auto *fg_s_59 = buffer.data(fg_s + 59);
    const auto *fg_s_60 = buffer.data(fg_s + 60);
    const auto *fg_s_61 = buffer.data(fg_s + 61);
    const auto *fg_s_62 = buffer.data(fg_s + 62);
    const auto *fg_s_63 = buffer.data(fg_s + 63);
    const auto *fg_s_64 = buffer.data(fg_s + 64);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_7 = buffer.data(ff + 7);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_14 = buffer.data(ff + 14);
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
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fg_s_3, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_2 * fg_s_1[k]
                 + pb_y[k] * ff_0[k];

        t_2[k] = f_2 * fg_s_2[k]
                 + pb_z[k] * ff_0[k];

        t_3[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_3[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, fd_s_0, fd_s_1, fg_s_4, fg_s_5, fd_0, fd_1, \
                         ff_2, ff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_4[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_5[k] = -f_1 * fd_s_1[k]
                 + f_2 * fg_s_5[k]
                 + f_0 * fd_1[k]
                 + pb_y[k] * ff_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, dg_0, fd_s_2, fg_s_6, fg_s_7, \
                         fg_s_8, fg_s_9, fd_2, ff_4, ff_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_3 * fd_s_2[k]
                 + f_2 * fg_s_6[k]
                 + f_4 * fd_2[k]
                 + pb_y[k] * ff_4[k];

        t_7[k] = f_2 * fg_s_7[k]
                 + pb_y[k] * ff_5[k];

        t_8[k] = -f_1 * fd_s_2[k]
                 + f_2 * fg_s_8[k]
                 + f_0 * fd_2[k]
                 + pb_z[k] * ff_5[k];

        t_9[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_z, pg_s_3, pg_3, dg_7, fd_s_4, fg_s_10, \
                         fg_s_11, fg_s_12, fd_4, ff_7, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_5 * pg_s_3[k]
                  + f_4 * pg_3[k]
                  + pa_x[k] * dg_7[k]
                  + f_2 * fg_s_10[k];

        t_11[k] = f_2 * fg_s_11[k]
                  + pb_z[k] * ff_7[k];

        t_12[k] = -f_3 * fd_s_4[k]
                  + f_2 * fg_s_12[k]
                  + f_4 * fd_4[k]
                  + pb_z[k] * ff_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, df_2, dg_0, dg_4, dg_6, fg_s_13, \
                         fg_s_14, fg_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * dg_6[k]
                  + f_2 * fg_s_13[k];

        t_14[k] = pa_z[k] * dg_0[k]
                  + f_2 * fg_s_14[k];

        t_15[k] = f_6 * df_2[k]
                  + pa_z[k] * dg_4[k]
                  + f_2 * fg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_y, fd_s_6, fd_s_7, fg_s_16, fg_s_17, fg_s_18, \
                         fd_6, fd_7, ff_10, ff_11, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_7 * fd_s_6[k]
                  + f_2 * fg_s_16[k]
                  + f_6 * fd_6[k]
                  + pb_y[k] * ff_10[k];

        t_17[k] = -f_3 * fd_s_7[k]
                  + f_2 * fg_s_17[k]
                  + f_4 * fd_7[k]
                  + pb_y[k] * ff_11[k];

        t_18[k] = f_2 * fg_s_18[k]
                  + pb_y[k] * ff_12[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pg_s_9, pg_9, df_7, df_8, dg_9, dg_10, dg_11, \
                         fg_s_19, fg_s_20, fg_s_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_5 * pg_s_9[k]
                  + f_4 * pg_9[k]
                  + pa_x[k] * dg_9[k]
                  + f_2 * fg_s_19[k];

        t_20[k] = f_8 * df_7[k]
                  + pa_x[k] * dg_10[k]
                  + f_2 * fg_s_20[k];

        t_21[k] = f_6 * df_8[k]
                  + pa_x[k] * dg_11[k]
                  + f_2 * fg_s_22[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pb_x, df_10, df_14, dg_15, dg_21, fg_s_24, \
                         fg_s_25, fg_s_26, ff_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_4 * df_10[k]
                  + f_2 * fg_s_24[k]
                  + pb_x[k] * ff_14[k];

        t_23[k] = pa_x[k] * dg_15[k]
                  + f_2 * fg_s_25[k];

        t_24[k] = f_8 * df_14[k]
                  + pa_x[k] * dg_21[k]
                  + f_2 * fg_s_26[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pb_x, df_17, df_21, dg_24, dg_31, fg_s_27, \
                         fg_s_28, fg_s_29, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * df_17[k]
                  + pa_x[k] * dg_24[k]
                  + f_2 * fg_s_27[k];

        t_26[k] = f_4 * df_21[k]
                  + f_2 * fg_s_28[k]
                  + pb_x[k] * ff_16[k];

        t_27[k] = pa_x[k] * dg_31[k]
                  + f_2 * fg_s_29[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, fd_s_12, fd_s_13, fd_s_14, fg_s_30, fg_s_31, \
                         fg_s_32, fd_10, fd_11, fd_12, ff_17, ff_18, \
                         ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_1 * fd_s_12[k]
                  + f_2 * fg_s_30[k]
                  + f_0 * fd_10[k]
                  + pb_x[k] * ff_17[k];

        t_29[k] = -f_7 * fd_s_13[k]
                  + f_2 * fg_s_31[k]
                  + f_6 * fd_11[k]
                  + pb_x[k] * ff_18[k];

        t_30[k] = -f_3 * fd_s_14[k]
                  + f_2 * fg_s_32[k]
                  + f_4 * fd_12[k]
                  + pb_x[k] * ff_19[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_x, fd_s_15, fg_s_33, fg_s_34, fg_s_35, \
                         fg_s_36, fd_13, ff_20, ff_21, ff_23, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_3 * fd_s_15[k]
                  + f_2 * fg_s_33[k]
                  + f_4 * fd_13[k]
                  + pb_x[k] * ff_20[k];

        t_32[k] = f_2 * fg_s_34[k]
                  + pb_x[k] * ff_21[k];

        t_33[k] = f_2 * fg_s_35[k]
                  + pb_x[k] * ff_23[k];

        t_34[k] = f_2 * fg_s_36[k]
                  + pb_x[k] * ff_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pb_z, df_10, fd_s_14, fg_s_37, fg_s_38, \
                         fg_s_39, fd_12, ff_21, ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * df_10[k]
                  - f_1 * fd_s_14[k]
                  + f_2 * fg_s_37[k]
                  + f_0 * fd_12[k]
                  + pb_y[k] * ff_21[k];

        t_36[k] = f_2 * fg_s_38[k]
                  + pb_z[k] * ff_21[k];

        t_37[k] = -f_3 * fd_s_14[k]
                  + f_2 * fg_s_39[k]
                  + f_4 * fd_12[k]
                  + pb_z[k] * ff_22[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, pb_y, pb_z, df_12, fd_s_15, fd_s_17, fg_s_40, \
                         fg_s_41, fg_s_42, fd_13, fd_15, ff_24, ff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * df_12[k]
                  + f_2 * fg_s_40[k]
                  + pb_y[k] * ff_24[k];

        t_39[k] = -f_1 * fd_s_15[k]
                  + f_2 * fg_s_41[k]
                  + f_0 * fd_13[k]
                  + pb_z[k] * ff_24[k];

        t_40[k] = -f_3 * fd_s_17[k]
                  + f_2 * fg_s_42[k]
                  + f_4 * fd_15[k]
                  + pb_x[k] * ff_25[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_z, pb_x, df_11, dg_15, dg_17, fg_s_43, fg_s_44, \
                         fg_s_45, ff_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_2 * fg_s_43[k]
                  + pb_x[k] * ff_27[k];

        t_42[k] = pa_z[k] * dg_15[k]
                  + f_2 * fg_s_44[k];

        t_43[k] = f_6 * df_11[k]
                  + pa_z[k] * dg_17[k]
                  + f_2 * fg_s_45[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_y, pb_x, pg_s_9, pg_9, df_15, dg_20, dg_23, \
                         fg_s_46, fg_s_47, fg_s_48, ff_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = -f_5 * pg_s_9[k]
                  + f_4 * pg_9[k]
                  + pa_y[k] * dg_20[k]
                  + f_2 * fg_s_46[k];

        t_45[k] = f_6 * df_15[k]
                  + pa_y[k] * dg_23[k]
                  + f_2 * fg_s_47[k];

        t_46[k] = f_2 * fg_s_48[k]
                  + pb_x[k] * ff_28[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pb_y, df_18, df_20, df_21, dg_27, dg_29, \
                         fg_s_49, fg_s_50, fg_s_51, ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_8 * df_18[k]
                  + pa_y[k] * dg_27[k]
                  + f_2 * fg_s_49[k];

        t_48[k] = f_6 * df_20[k]
                  + pa_y[k] * dg_29[k]
                  + f_2 * fg_s_50[k];

        t_49[k] = f_4 * df_21[k]
                  + f_2 * fg_s_51[k]
                  + pb_y[k] * ff_30[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_y, pb_x, dg_31, fd_s_20, fd_s_21, fg_s_52, \
                         fg_s_53, fg_s_54, fd_17, fd_18, ff_31, ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * dg_31[k]
                  + f_2 * fg_s_52[k];

        t_51[k] = -f_1 * fd_s_20[k]
                  + f_2 * fg_s_53[k]
                  + f_0 * fd_17[k]
                  + pb_x[k] * ff_31[k];

        t_52[k] = -f_7 * fd_s_21[k]
                  + f_2 * fg_s_54[k]
                  + f_6 * fd_18[k]
                  + pb_x[k] * ff_32[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, fd_s_22, fd_s_24, fg_s_55, fg_s_56, fg_s_57, \
                         fd_19, fd_21, ff_33, ff_34, ff_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_3 * fd_s_22[k]
                  + f_2 * fg_s_55[k]
                  + f_4 * fd_19[k]
                  + pb_x[k] * ff_33[k];

        t_54[k] = -f_3 * fd_s_24[k]
                  + f_2 * fg_s_56[k]
                  + f_4 * fd_21[k]
                  + pb_x[k] * ff_34[k];

        t_55[k] = f_2 * fg_s_57[k]
                  + pb_x[k] * ff_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pb_x, pb_y, fd_s_22, fg_s_58, fg_s_59, fg_s_60, \
                         fd_19, ff_35, ff_36, ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_2 * fg_s_58[k]
                  + pb_x[k] * ff_36[k];

        t_57[k] = f_2 * fg_s_59[k]
                  + pb_x[k] * ff_38[k];

        t_58[k] = -f_1 * fd_s_22[k]
                  + f_2 * fg_s_60[k]
                  + f_0 * fd_19[k]
                  + pb_y[k] * ff_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_y, fd_s_23, fd_s_24, fg_s_61, fg_s_62, fg_s_63, \
                         fd_20, fd_21, ff_36, ff_37, ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_7 * fd_s_23[k]
                  + f_2 * fg_s_61[k]
                  + f_6 * fd_20[k]
                  + pb_y[k] * ff_36[k];

        t_60[k] = -f_3 * fd_s_24[k]
                  + f_2 * fg_s_62[k]
                  + f_4 * fd_21[k]
                  + pb_y[k] * ff_37[k];

        t_61[k] = f_2 * fg_s_63[k]
                  + pb_y[k] * ff_38[k];
    }

#pragma omp simd aligned(t_62, pb_z, df_21, fd_s_24, fg_s_64, fd_21, \
                         ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * df_21[k]
                  - f_1 * fd_s_24[k]
                  + f_2 * fg_s_64[k]
                  + f_0 * fd_21[k]
                  + pb_z[k] * ff_38[k];
    }
}

auto
compute_prim_fg_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;

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

    const auto *pg_s_0 = buffer.data(pg_s + 0);
    const auto *pg_s_2 = buffer.data(pg_s + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_20 = buffer.data(fd_s + 20);
    const auto *fd_s_22 = buffer.data(fd_s + 22);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_23 = buffer.data(fg_s + 23);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_22 = buffer.data(fd + 22);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_37 = buffer.data(ff + 37);
    const auto *ff_38 = buffer.data(ff + 38);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fd_0, ff_0, ff_1, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_1[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];

        t_2[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_2[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, pg_s_0, pg_0, dg_0, dg_2, fg_s_3, \
                         fg_s_4, fg_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_3[k];

        t_4[k] = -f_5 * pg_s_0[k]
                 + f_4 * pg_0[k]
                 + pa_x[k] * dg_2[k]
                 + f_2 * fg_s_4[k];

        t_5[k] = pa_z[k] * dg_0[k]
                 + f_2 * fg_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pg_s_2, pg_2, df_1, dg_1, dg_3, dg_4, \
                         fg_s_6, fg_s_7, fg_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * df_1[k]
                 + pa_z[k] * dg_1[k]
                 + f_2 * fg_s_6[k];

        t_7[k] = -f_5 * pg_s_2[k]
                 + f_4 * pg_2[k]
                 + pa_x[k] * dg_3[k]
                 + f_2 * fg_s_7[k];

        t_8[k] = pa_x[k] * dg_4[k]
                 + f_2 * fg_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, dg_6, dg_10, fd_s_12, fg_s_9, fg_s_10, \
                         fg_s_11, fd_12, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_x[k] * dg_6[k]
                 + f_2 * fg_s_9[k];

        t_10[k] = pa_x[k] * dg_10[k]
                  + f_2 * fg_s_10[k];

        t_11[k] = -f_1 * fd_s_12[k]
                  + f_2 * fg_s_11[k]
                  + f_0 * fd_12[k]
                  + pb_x[k] * ff_17[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, pb_y, pb_z, df_5, fd_s_13, fg_s_12, fg_s_13, \
                         fg_s_14, fd_13, ff_18, ff_19, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_3 * fd_s_13[k]
                  + f_2 * fg_s_12[k]
                  + f_4 * fd_13[k]
                  + pb_x[k] * ff_18[k];

        t_13[k] = f_0 * df_5[k]
                  - f_1 * fd_s_13[k]
                  + f_2 * fg_s_13[k]
                  + f_0 * fd_13[k]
                  + pb_y[k] * ff_19[k];

        t_14[k] = -f_3 * fd_s_13[k]
                  + f_2 * fg_s_14[k]
                  + f_4 * fd_13[k]
                  + pb_z[k] * ff_20[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_y, pa_z, pg_s_2, pg_2, df_6, dg_4, dg_5, dg_7, \
                         fg_s_15, fg_s_16, fg_s_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_z[k] * dg_4[k]
                  + f_2 * fg_s_15[k];

        t_16[k] = f_6 * df_6[k]
                  + pa_z[k] * dg_5[k]
                  + f_2 * fg_s_16[k];

        t_17[k] = -f_5 * pg_s_2[k]
                  + f_4 * pg_2[k]
                  + pa_y[k] * dg_7[k]
                  + f_2 * fg_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, df_9, df_10, dg_8, dg_9, dg_10, fg_s_18, \
                         fg_s_19, fg_s_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * df_9[k]
                  + pa_y[k] * dg_8[k]
                  + f_2 * fg_s_18[k];

        t_19[k] = f_6 * df_10[k]
                  + pa_y[k] * dg_9[k]
                  + f_2 * fg_s_19[k];

        t_20[k] = pa_y[k] * dg_10[k]
                  + f_2 * fg_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, fd_s_19, fd_s_20, fd_s_22, fg_s_21, fg_s_22, \
                         fg_s_23, fd_19, fd_20, fd_22, ff_31, ff_33, \
                         ff_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * fd_s_19[k]
                  + f_2 * fg_s_21[k]
                  + f_0 * fd_19[k]
                  + pb_x[k] * ff_31[k];

        t_22[k] = -f_3 * fd_s_20[k]
                  + f_2 * fg_s_22[k]
                  + f_4 * fd_20[k]
                  + pb_x[k] * ff_33[k];

        t_23[k] = -f_3 * fd_s_22[k]
                  + f_2 * fg_s_23[k]
                  + f_4 * fd_22[k]
                  + pb_x[k] * ff_34[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pb_z, df_11, fd_s_20, fd_s_22, fg_s_24, \
                         fg_s_25, fg_s_26, fd_20, fd_22, ff_35, ff_37, \
                         ff_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * fd_s_20[k]
                  + f_2 * fg_s_24[k]
                  + f_0 * fd_20[k]
                  + pb_y[k] * ff_35[k];

        t_25[k] = -f_3 * fd_s_22[k]
                  + f_2 * fg_s_25[k]
                  + f_4 * fd_22[k]
                  + pb_y[k] * ff_37[k];

        t_26[k] = f_0 * df_11[k]
                  - f_1 * fd_s_22[k]
                  + f_2 * fg_s_26[k]
                  + f_0 * fd_22[k]
                  + pb_z[k] * ff_38[k];
    }
}

auto
compute_prim_fg_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.0 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = 2.0 * alpha / p;

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

    const auto *pg_s_1 = buffer.data(pg_s + 1);
    const auto *pg_s_4 = buffer.data(pg_s + 4);

    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_19 = buffer.data(dg + 19);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_23 = buffer.data(fg_s + 23);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_29 = buffer.data(fg_s + 29);
    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_31 = buffer.data(fg_s + 31);
    const auto *fg_s_32 = buffer.data(fg_s + 32);
    const auto *fg_s_33 = buffer.data(fg_s + 33);
    const auto *fg_s_34 = buffer.data(fg_s + 34);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_36 = buffer.data(fg_s + 36);
    const auto *fg_s_37 = buffer.data(fg_s + 37);
    const auto *fg_s_38 = buffer.data(fg_s + 38);
    const auto *fg_s_39 = buffer.data(fg_s + 39);
    const auto *fg_s_40 = buffer.data(fg_s + 40);
    const auto *fg_s_41 = buffer.data(fg_s + 41);
    const auto *fg_s_42 = buffer.data(fg_s + 42);
    const auto *fg_s_43 = buffer.data(fg_s + 43);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_45 = buffer.data(fg_s + 45);
    const auto *fg_s_46 = buffer.data(fg_s + 46);
    const auto *fg_s_47 = buffer.data(fg_s + 47);
    const auto *fg_s_48 = buffer.data(fg_s + 48);
    const auto *fg_s_49 = buffer.data(fg_s + 49);
    const auto *fg_s_50 = buffer.data(fg_s + 50);
    const auto *fg_s_51 = buffer.data(fg_s + 51);
    const auto *fg_s_52 = buffer.data(fg_s + 52);
    const auto *fg_s_53 = buffer.data(fg_s + 53);
    const auto *fg_s_54 = buffer.data(fg_s + 54);
    const auto *fg_s_55 = buffer.data(fg_s + 55);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
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
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_34 = buffer.data(ff + 34);
    const auto *ff_35 = buffer.data(ff + 35);
    const auto *ff_36 = buffer.data(ff + 36);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fg_s_3, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_2 * fg_s_1[k]
                 + pb_y[k] * ff_0[k];

        t_2[k] = f_2 * fg_s_2[k]
                 + pb_z[k] * ff_0[k];

        t_3[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_3[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, fd_s_0, fd_s_1, fg_s_4, fg_s_5, fd_0, fd_1, \
                         ff_2, ff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_4[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_5[k] = -f_1 * fd_s_1[k]
                 + f_2 * fg_s_5[k]
                 + f_0 * fd_1[k]
                 + pb_y[k] * ff_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_y, pb_z, df_1, dg_0, dg_2, fd_s_2, fg_s_6, fg_s_7, \
                         fg_s_8, fd_2, ff_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * fd_s_2[k]
                 + f_2 * fg_s_6[k]
                 + f_0 * fd_2[k]
                 + pb_z[k] * ff_4[k];

        t_7[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_7[k];

        t_8[k] = f_5 * df_1[k]
                 + pa_y[k] * dg_2[k]
                 + f_2 * fg_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, pb_z, pg_s_1, pg_1, df_0, dg_0, dg_4, \
                         fg_s_9, fg_s_10, fg_s_11, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_6 * pg_s_1[k]
                 + f_4 * pg_1[k]
                 + pa_x[k] * dg_4[k]
                 + f_2 * fg_s_9[k];

        t_10[k] = pa_z[k] * dg_0[k]
                  + f_2 * fg_s_10[k];

        t_11[k] = f_4 * df_0[k]
                  + f_2 * fg_s_11[k]
                  + pb_z[k] * ff_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pa_z, pg_s_4, pg_4, df_2, df_7, dg_3, dg_5, \
                         dg_6, fg_s_12, fg_s_13, fg_s_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_3[k]
                  + f_2 * fg_s_12[k];

        t_13[k] = -f_6 * pg_s_4[k]
                  + f_4 * pg_4[k]
                  + pa_x[k] * dg_5[k]
                  + f_2 * fg_s_13[k];

        t_14[k] = f_7 * df_7[k]
                  + pa_x[k] * dg_6[k]
                  + f_2 * fg_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, df_8, df_9, dg_7, dg_8, dg_11, \
                         fg_s_15, fg_s_16, fg_s_17, fg_s_18, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_5 * df_8[k]
                  + pa_x[k] * dg_7[k]
                  + f_2 * fg_s_15[k];

        t_16[k] = f_4 * df_9[k]
                  + f_2 * fg_s_16[k]
                  + pb_x[k] * ff_12[k];

        t_17[k] = pa_x[k] * dg_8[k]
                  + f_2 * fg_s_17[k];

        t_18[k] = pa_x[k] * dg_11[k]
                  + f_2 * fg_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pb_z, df_5, df_14, df_17, dg_13, dg_15, \
                         fg_s_19, fg_s_20, fg_s_21, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * df_14[k]
                  + pa_x[k] * dg_13[k]
                  + f_2 * fg_s_19[k];

        t_20[k] = f_5 * df_5[k]
                  + f_2 * fg_s_20[k]
                  + pb_z[k] * ff_13[k];

        t_21[k] = f_5 * df_17[k]
                  + pa_x[k] * dg_15[k]
                  + f_2 * fg_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pb_x, df_20, dg_19, fd_s_9, fg_s_22, fg_s_23, \
                         fg_s_24, fd_9, ff_15, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_4 * df_20[k]
                  + f_2 * fg_s_22[k]
                  + pb_x[k] * ff_15[k];

        t_23[k] = pa_x[k] * dg_19[k]
                  + f_2 * fg_s_23[k];

        t_24[k] = -f_1 * fd_s_9[k]
                  + f_2 * fg_s_24[k]
                  + f_0 * fd_9[k]
                  + pb_x[k] * ff_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_x, fd_s_10, fd_s_11, fg_s_25, fg_s_26, fg_s_27, \
                         fd_10, fd_11, ff_17, ff_18, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_3 * fd_s_10[k]
                  + f_2 * fg_s_25[k]
                  + f_4 * fd_10[k]
                  + pb_x[k] * ff_17[k];

        t_26[k] = -f_3 * fd_s_11[k]
                  + f_2 * fg_s_26[k]
                  + f_4 * fd_11[k]
                  + pb_x[k] * ff_18[k];

        t_27[k] = f_2 * fg_s_27[k]
                  + pb_x[k] * ff_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_y, pb_z, df_9, fd_s_10, fg_s_28, fg_s_29, \
                         fg_s_30, fd_10, ff_19, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * df_9[k]
                  - f_1 * fd_s_10[k]
                  + f_2 * fg_s_28[k]
                  + f_0 * fd_10[k]
                  + pb_y[k] * ff_19[k];

        t_29[k] = f_2 * fg_s_29[k]
                  + pb_z[k] * ff_19[k];

        t_30[k] = -f_3 * fd_s_10[k]
                  + f_2 * fg_s_30[k]
                  + f_4 * fd_10[k]
                  + pb_z[k] * ff_20[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_x, pb_y, pb_z, df_11, fd_s_11, fd_s_13, fg_s_31, \
                         fg_s_32, fg_s_33, fd_11, fd_13, ff_21, ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * df_11[k]
                  + f_2 * fg_s_31[k]
                  + pb_y[k] * ff_21[k];

        t_32[k] = -f_1 * fd_s_11[k]
                  + f_2 * fg_s_32[k]
                  + f_0 * fd_11[k]
                  + pb_z[k] * ff_21[k];

        t_33[k] = -f_3 * fd_s_13[k]
                  + f_2 * fg_s_33[k]
                  + f_4 * fd_13[k]
                  + pb_x[k] * ff_22[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_z, df_9, df_10, dg_8, dg_10, fg_s_34, \
                         fg_s_35, fg_s_36, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * dg_8[k]
                  + f_2 * fg_s_34[k];

        t_35[k] = f_4 * df_9[k]
                  + f_2 * fg_s_35[k]
                  + pb_z[k] * ff_23[k];

        t_36[k] = f_5 * df_10[k]
                  + pa_z[k] * dg_10[k]
                  + f_2 * fg_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pb_y, pg_s_4, pg_4, df_13, df_15, dg_12, \
                         dg_14, fg_s_37, fg_s_38, fg_s_39, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * df_13[k]
                  + f_2 * fg_s_37[k]
                  + pb_y[k] * ff_24[k];

        t_38[k] = -f_6 * pg_s_4[k]
                  + f_4 * pg_4[k]
                  + pa_y[k] * dg_12[k]
                  + f_2 * fg_s_38[k];

        t_39[k] = f_5 * df_15[k]
                  + pa_y[k] * dg_14[k]
                  + f_2 * fg_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pb_z, df_12, df_18, df_19, dg_16, dg_17, \
                         fg_s_40, fg_s_41, fg_s_42, ff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * df_18[k]
                  + pa_y[k] * dg_16[k]
                  + f_2 * fg_s_40[k];

        t_41[k] = f_5 * df_12[k]
                  + f_2 * fg_s_41[k]
                  + pb_z[k] * ff_26[k];

        t_42[k] = f_5 * df_19[k]
                  + pa_y[k] * dg_17[k]
                  + f_2 * fg_s_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_y, pb_x, pb_y, df_20, dg_19, fd_s_16, fg_s_43, \
                         fg_s_44, fg_s_45, fd_16, ff_28, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_4 * df_20[k]
                  + f_2 * fg_s_43[k]
                  + pb_y[k] * ff_28[k];

        t_44[k] = pa_y[k] * dg_19[k]
                  + f_2 * fg_s_44[k];

        t_45[k] = -f_1 * fd_s_16[k]
                  + f_2 * fg_s_45[k]
                  + f_0 * fd_16[k]
                  + pb_x[k] * ff_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pb_y, fd_s_17, fd_s_19, fg_s_46, fg_s_47, \
                         fg_s_48, fd_17, fd_19, ff_29, ff_31, ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * fg_s_46[k]
                  + pb_y[k] * ff_29[k];

        t_47[k] = -f_3 * fd_s_17[k]
                  + f_2 * fg_s_47[k]
                  + f_4 * fd_17[k]
                  + pb_x[k] * ff_31[k];

        t_48[k] = -f_3 * fd_s_19[k]
                  + f_2 * fg_s_48[k]
                  + f_4 * fd_19[k]
                  + pb_x[k] * ff_32[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, pb_y, fd_s_17, fg_s_49, fg_s_50, fg_s_51, \
                         fd_17, ff_33, ff_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_2 * fg_s_49[k]
                  + pb_x[k] * ff_33[k];

        t_50[k] = f_2 * fg_s_50[k]
                  + pb_x[k] * ff_36[k];

        t_51[k] = -f_1 * fd_s_17[k]
                  + f_2 * fg_s_51[k]
                  + f_0 * fd_17[k]
                  + pb_y[k] * ff_33[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_y, fd_s_18, fd_s_19, fg_s_52, fg_s_53, fg_s_54, \
                         fd_18, fd_19, ff_34, ff_35, ff_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -f_8 * fd_s_18[k]
                  + f_2 * fg_s_52[k]
                  + f_5 * fd_18[k]
                  + pb_y[k] * ff_34[k];

        t_53[k] = -f_3 * fd_s_19[k]
                  + f_2 * fg_s_53[k]
                  + f_4 * fd_19[k]
                  + pb_y[k] * ff_35[k];

        t_54[k] = f_2 * fg_s_54[k]
                  + pb_y[k] * ff_36[k];
    }

#pragma omp simd aligned(t_55, pb_z, df_20, fd_s_19, fg_s_55, fd_19, \
                         ff_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * df_20[k]
                  - f_1 * fd_s_19[k]
                  + f_2 * fg_s_55[k]
                  + f_0 * fd_19[k]
                  + pb_z[k] * ff_36[k];
    }
}

auto
compute_prim_fg_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 2.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_s_1 = buffer.data(pg_s + 1);
    const auto *pg_s_4 = buffer.data(pg_s + 4);

    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_27 = buffer.data(dg + 27);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_23 = buffer.data(fg_s + 23);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_29 = buffer.data(fg_s + 29);
    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_31 = buffer.data(fg_s + 31);
    const auto *fg_s_32 = buffer.data(fg_s + 32);
    const auto *fg_s_33 = buffer.data(fg_s + 33);
    const auto *fg_s_34 = buffer.data(fg_s + 34);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_36 = buffer.data(fg_s + 36);
    const auto *fg_s_37 = buffer.data(fg_s + 37);
    const auto *fg_s_38 = buffer.data(fg_s + 38);
    const auto *fg_s_39 = buffer.data(fg_s + 39);
    const auto *fg_s_40 = buffer.data(fg_s + 40);
    const auto *fg_s_41 = buffer.data(fg_s + 41);
    const auto *fg_s_42 = buffer.data(fg_s + 42);
    const auto *fg_s_43 = buffer.data(fg_s + 43);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_45 = buffer.data(fg_s + 45);
    const auto *fg_s_46 = buffer.data(fg_s + 46);
    const auto *fg_s_47 = buffer.data(fg_s + 47);
    const auto *fg_s_48 = buffer.data(fg_s + 48);
    const auto *fg_s_49 = buffer.data(fg_s + 49);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_4 = buffer.data(ff + 4);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_32 = buffer.data(ff + 32);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fg_s_3, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_2 * fg_s_1[k]
                 + pb_y[k] * ff_0[k];

        t_2[k] = f_2 * fg_s_2[k]
                 + pb_z[k] * ff_0[k];

        t_3[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_3[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, fd_s_0, fd_s_1, fg_s_4, fg_s_5, fd_0, fd_1, \
                         ff_2, ff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_4[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_5[k] = -f_1 * fd_s_1[k]
                 + f_2 * fg_s_5[k]
                 + f_0 * fd_1[k]
                 + pb_y[k] * ff_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_z, pg_s_1, pg_1, dg_0, dg_6, fd_s_2, \
                         fg_s_6, fg_s_7, fg_s_8, fd_2, ff_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * fd_s_2[k]
                 + f_2 * fg_s_6[k]
                 + f_0 * fd_2[k]
                 + pb_z[k] * ff_4[k];

        t_7[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_7[k];

        t_8[k] = -f_5 * pg_s_1[k]
                 + f_4 * pg_1[k]
                 + pa_x[k] * dg_6[k]
                 + f_2 * fg_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_y, pa_z, pg_s_4, pg_4, dg_0, dg_5, dg_8, \
                         fg_s_9, fg_s_10, fg_s_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_y[k] * dg_5[k]
                 + f_2 * fg_s_9[k];

        t_10[k] = pa_z[k] * dg_0[k]
                  + f_2 * fg_s_10[k];

        t_11[k] = -f_5 * pg_s_4[k]
                  + f_4 * pg_4[k]
                  + pa_x[k] * dg_8[k]
                  + f_2 * fg_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_x, df_6, df_7, df_8, dg_9, dg_10, fg_s_12, \
                         fg_s_13, fg_s_14, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_6 * df_6[k]
                  + pa_x[k] * dg_9[k]
                  + f_2 * fg_s_12[k];

        t_13[k] = f_7 * df_7[k]
                  + pa_x[k] * dg_10[k]
                  + f_2 * fg_s_13[k];

        t_14[k] = f_4 * df_8[k]
                  + f_2 * fg_s_14[k]
                  + pb_x[k] * ff_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, df_12, df_15, dg_12, dg_18, dg_21, fg_s_15, \
                         fg_s_16, fg_s_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_x[k] * dg_12[k]
                  + f_2 * fg_s_15[k];

        t_16[k] = f_6 * df_12[k]
                  + pa_x[k] * dg_18[k]
                  + f_2 * fg_s_16[k];

        t_17[k] = f_7 * df_15[k]
                  + pa_x[k] * dg_21[k]
                  + f_2 * fg_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pb_x, df_18, dg_27, fd_s_9, fg_s_18, fg_s_19, \
                         fg_s_20, fd_9, ff_12, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_4 * df_18[k]
                  + f_2 * fg_s_18[k]
                  + pb_x[k] * ff_12[k];

        t_19[k] = pa_x[k] * dg_27[k]
                  + f_2 * fg_s_19[k];

        t_20[k] = -f_1 * fd_s_9[k]
                  + f_2 * fg_s_20[k]
                  + f_0 * fd_9[k]
                  + pb_x[k] * ff_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, fd_s_10, fd_s_11, fg_s_21, fg_s_22, fg_s_23, \
                         fd_10, fd_11, ff_14, ff_15, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_3 * fd_s_10[k]
                  + f_2 * fg_s_21[k]
                  + f_4 * fd_10[k]
                  + pb_x[k] * ff_14[k];

        t_22[k] = -f_3 * fd_s_11[k]
                  + f_2 * fg_s_22[k]
                  + f_4 * fd_11[k]
                  + pb_x[k] * ff_15[k];

        t_23[k] = f_2 * fg_s_23[k]
                  + pb_x[k] * ff_16[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pb_y, pb_z, df_8, fd_s_10, fg_s_24, fg_s_25, \
                         fg_s_26, fd_10, ff_16, ff_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_2 * fg_s_24[k]
                  + pb_x[k] * ff_18[k];

        t_25[k] = f_0 * df_8[k]
                  - f_1 * fd_s_10[k]
                  + f_2 * fg_s_25[k]
                  + f_0 * fd_10[k]
                  + pb_y[k] * ff_16[k];

        t_26[k] = f_2 * fg_s_26[k]
                  + pb_z[k] * ff_16[k];
    }

#pragma omp simd aligned(t_27, t_28, pb_z, fd_s_10, fd_s_11, fg_s_27, fg_s_28, fd_10, fd_11, \
                         ff_17, ff_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_3 * fd_s_10[k]
                  + f_2 * fg_s_27[k]
                  + f_4 * fd_10[k]
                  + pb_z[k] * ff_17[k];

        t_28[k] = -f_1 * fd_s_11[k]
                  + f_2 * fg_s_28[k]
                  + f_0 * fd_11[k]
                  + pb_z[k] * ff_18[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_z, pb_x, dg_12, fd_s_13, fg_s_29, fg_s_30, \
                         fg_s_31, fd_13, ff_19, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_3 * fd_s_13[k]
                  + f_2 * fg_s_29[k]
                  + f_4 * fd_13[k]
                  + pb_x[k] * ff_19[k];

        t_30[k] = f_2 * fg_s_30[k]
                  + pb_x[k] * ff_21[k];

        t_31[k] = pa_z[k] * dg_12[k]
                  + f_2 * fg_s_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pb_x, pg_s_4, pg_4, df_13, dg_17, dg_20, \
                         fg_s_32, fg_s_33, fg_s_34, ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -f_5 * pg_s_4[k]
                  + f_4 * pg_4[k]
                  + pa_y[k] * dg_17[k]
                  + f_2 * fg_s_32[k];

        t_33[k] = f_7 * df_13[k]
                  + pa_y[k] * dg_20[k]
                  + f_2 * fg_s_33[k];

        t_34[k] = f_2 * fg_s_34[k]
                  + pb_x[k] * ff_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pb_y, df_16, df_17, df_18, dg_24, dg_25, \
                         fg_s_35, fg_s_36, fg_s_37, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * df_16[k]
                  + pa_y[k] * dg_24[k]
                  + f_2 * fg_s_35[k];

        t_36[k] = f_7 * df_17[k]
                  + pa_y[k] * dg_25[k]
                  + f_2 * fg_s_36[k];

        t_37[k] = f_4 * df_18[k]
                  + f_2 * fg_s_37[k]
                  + pb_y[k] * ff_24[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pb_x, pb_y, dg_27, fd_s_16, fg_s_38, fg_s_39, \
                         fg_s_40, fd_15, ff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * dg_27[k]
                  + f_2 * fg_s_38[k];

        t_39[k] = -f_1 * fd_s_16[k]
                  + f_2 * fg_s_39[k]
                  + f_0 * fd_15[k]
                  + pb_x[k] * ff_25[k];

        t_40[k] = f_2 * fg_s_40[k]
                  + pb_y[k] * ff_25[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pb_x, fd_s_17, fd_s_19, fg_s_41, fg_s_42, fg_s_43, \
                         fd_16, fd_18, ff_27, ff_28, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -f_3 * fd_s_17[k]
                  + f_2 * fg_s_41[k]
                  + f_4 * fd_16[k]
                  + pb_x[k] * ff_27[k];

        t_42[k] = -f_3 * fd_s_19[k]
                  + f_2 * fg_s_42[k]
                  + f_4 * fd_18[k]
                  + pb_x[k] * ff_28[k];

        t_43[k] = f_2 * fg_s_43[k]
                  + pb_x[k] * ff_29[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, pb_y, fd_s_17, fd_s_18, fg_s_44, fg_s_45, \
                         fg_s_46, fd_16, fd_17, ff_29, ff_30, ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_2 * fg_s_44[k]
                  + pb_x[k] * ff_32[k];

        t_45[k] = -f_1 * fd_s_17[k]
                  + f_2 * fg_s_45[k]
                  + f_0 * fd_16[k]
                  + pb_y[k] * ff_29[k];

        t_46[k] = -f_8 * fd_s_18[k]
                  + f_2 * fg_s_46[k]
                  + f_7 * fd_17[k]
                  + pb_y[k] * ff_30[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pb_y, pb_z, df_18, fd_s_19, fg_s_47, fg_s_48, \
                         fg_s_49, fd_18, ff_31, ff_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -f_3 * fd_s_19[k]
                  + f_2 * fg_s_47[k]
                  + f_4 * fd_18[k]
                  + pb_y[k] * ff_31[k];

        t_48[k] = f_2 * fg_s_48[k]
                  + pb_y[k] * ff_32[k];

        t_49[k] = f_0 * df_18[k]
                  - f_1 * fd_s_19[k]
                  + f_2 * fg_s_49[k]
                  + f_0 * fd_18[k]
                  + pb_z[k] * ff_32[k];
    }
}

auto
compute_prim_fg_kinetic_energy_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg_s, const size_t pg,
                                 const size_t df, const size_t dg, const size_t fd_s,
                                 const size_t fg_s, const size_t fd, const size_t ff,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_s_0 = buffer.data(pg_s + 0);
    const auto *pg_s_2 = buffer.data(pg_s + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_2 = buffer.data(pg + 2);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_19 = buffer.data(fd_s + 19);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);

#pragma omp simd aligned(t_0, t_1, t_2, pa_y, pb_x, pb_z, df_0, dg_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_1[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_1[k];

        t_2[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_z, pg_s_0, pg_0, df_1, dg_0, dg_1, dg_2, \
                         fg_s_3, fg_s_4, fg_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_5 * pg_s_0[k]
                 + f_4 * pg_0[k]
                 + pa_x[k] * dg_2[k]
                 + f_2 * fg_s_3[k];

        t_4[k] = pa_z[k] * dg_0[k]
                 + f_2 * fg_s_4[k];

        t_5[k] = f_6 * df_1[k]
                 + pa_z[k] * dg_1[k]
                 + f_2 * fg_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pg_s_2, pg_2, dg_3, dg_4, dg_6, dg_10, \
                         fg_s_6, fg_s_7, fg_s_8, fg_s_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_5 * pg_s_2[k]
                 + f_4 * pg_2[k]
                 + pa_x[k] * dg_3[k]
                 + f_2 * fg_s_6[k];

        t_7[k] = pa_x[k] * dg_4[k]
                 + f_2 * fg_s_7[k];

        t_8[k] = pa_x[k] * dg_6[k]
                 + f_2 * fg_s_8[k];

        t_9[k] = pa_x[k] * dg_10[k]
                 + f_2 * fg_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_z, pb_y, pb_z, df_4, dg_4, fd_s_11, fg_s_10, \
                         fg_s_11, fg_s_12, fd_11, ff_11, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * df_4[k]
                  - f_1 * fd_s_11[k]
                  + f_2 * fg_s_10[k]
                  + f_0 * fd_11[k]
                  + pb_y[k] * ff_11[k];

        t_11[k] = -f_3 * fd_s_11[k]
                  + f_2 * fg_s_11[k]
                  + f_4 * fd_11[k]
                  + pb_z[k] * ff_12[k];

        t_12[k] = pa_z[k] * dg_4[k]
                  + f_2 * fg_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pg_s_2, pg_2, df_5, df_7, dg_5, dg_7, \
                         dg_8, fg_s_13, fg_s_14, fg_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * df_5[k]
                  + pa_z[k] * dg_5[k]
                  + f_2 * fg_s_13[k];

        t_14[k] = -f_5 * pg_s_2[k]
                  + f_4 * pg_2[k]
                  + pa_y[k] * dg_7[k]
                  + f_2 * fg_s_14[k];

        t_15[k] = f_7 * df_7[k]
                  + pa_y[k] * dg_8[k]
                  + f_2 * fg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_y, pb_y, df_8, dg_9, dg_10, fd_s_18, fg_s_16, \
                         fg_s_17, fg_s_18, fd_18, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_6 * df_8[k]
                  + pa_y[k] * dg_9[k]
                  + f_2 * fg_s_16[k];

        t_17[k] = pa_y[k] * dg_10[k]
                  + f_2 * fg_s_17[k];

        t_18[k] = -f_1 * fd_s_18[k]
                  + f_2 * fg_s_18[k]
                  + f_0 * fd_18[k]
                  + pb_y[k] * ff_19[k];
    }

#pragma omp simd aligned(t_19, t_20, pb_y, pb_z, df_9, fd_s_19, fg_s_19, fg_s_20, fd_19, \
                         ff_20, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_3 * fd_s_19[k]
                  + f_2 * fg_s_19[k]
                  + f_4 * fd_19[k]
                  + pb_y[k] * ff_20[k];

        t_20[k] = f_0 * df_9[k]
                  - f_1 * fd_s_19[k]
                  + f_2 * fg_s_20[k]
                  + f_0 * fd_19[k]
                  + pb_z[k] * ff_21[k];
    }
}

auto
compute_prim_fg_kinetic_energy_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pg_s, const size_t pg,
                                  const size_t df, const size_t dg, const size_t fd_s,
                                  const size_t fg_s, const size_t fd, const size_t ff,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_s_1 = buffer.data(pg_s + 1);
    const auto *pg_s_4 = buffer.data(pg_s + 4);

    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_15 = buffer.data(dg + 15);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_23 = buffer.data(fg_s + 23);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_29 = buffer.data(fg_s + 29);
    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_31 = buffer.data(fg_s + 31);
    const auto *fg_s_32 = buffer.data(fg_s + 32);
    const auto *fg_s_33 = buffer.data(fg_s + 33);
    const auto *fg_s_34 = buffer.data(fg_s + 34);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_2 * fg_s_1[k]
                 + pb_z[k] * ff_0[k];

        t_2[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_2[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_z, pg_s_1, pg_1, dg_0, dg_3, fd_s_0, \
                         fg_s_3, fg_s_4, fg_s_5, fd_0, ff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_3[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_4[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_4[k];

        t_5[k] = -f_5 * pg_s_1[k]
                 + f_4 * pg_1[k]
                 + pa_x[k] * dg_3[k]
                 + f_2 * fg_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pg_s_4, pg_4, df_1, dg_0, dg_2, dg_4, \
                         fg_s_6, fg_s_7, fg_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_z[k] * dg_0[k]
                 + f_2 * fg_s_6[k];

        t_7[k] = f_6 * df_1[k]
                 + pa_z[k] * dg_2[k]
                 + f_2 * fg_s_7[k];

        t_8[k] = -f_5 * pg_s_4[k]
                 + f_4 * pg_4[k]
                 + pa_x[k] * dg_4[k]
                 + f_2 * fg_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_x, pb_x, df_5, df_6, dg_5, dg_6, dg_9, \
                         fg_s_9, fg_s_10, fg_s_11, fg_s_12, ff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * df_5[k]
                 + pa_x[k] * dg_5[k]
                 + f_2 * fg_s_9[k];

        t_10[k] = f_4 * df_6[k]
                  + f_2 * fg_s_10[k]
                  + pb_x[k] * ff_10[k];

        t_11[k] = pa_x[k] * dg_6[k]
                  + f_2 * fg_s_11[k];

        t_12[k] = pa_x[k] * dg_9[k]
                  + f_2 * fg_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_x, df_10, df_13, dg_11, dg_15, fg_s_13, \
                         fg_s_14, fg_s_15, ff_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * df_10[k]
                  + pa_x[k] * dg_11[k]
                  + f_2 * fg_s_13[k];

        t_14[k] = f_4 * df_13[k]
                  + f_2 * fg_s_14[k]
                  + pb_x[k] * ff_13[k];

        t_15[k] = pa_x[k] * dg_15[k]
                  + f_2 * fg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_y, df_6, fd_s_8, fd_s_9, fg_s_16, fg_s_17, \
                         fg_s_18, fd_8, fd_9, ff_14, ff_15, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_1 * fd_s_8[k]
                  + f_2 * fg_s_16[k]
                  + f_0 * fd_8[k]
                  + pb_x[k] * ff_14[k];

        t_17[k] = -f_3 * fd_s_9[k]
                  + f_2 * fg_s_17[k]
                  + f_4 * fd_9[k]
                  + pb_x[k] * ff_15[k];

        t_18[k] = f_0 * df_6[k]
                  - f_1 * fd_s_9[k]
                  + f_2 * fg_s_18[k]
                  + f_0 * fd_9[k]
                  + pb_y[k] * ff_16[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_z, dg_6, fd_s_9, fg_s_19, fg_s_20, \
                         fg_s_21, fd_9, ff_16, ff_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * fg_s_19[k]
                  + pb_z[k] * ff_16[k];

        t_20[k] = -f_3 * fd_s_9[k]
                  + f_2 * fg_s_20[k]
                  + f_4 * fd_9[k]
                  + pb_z[k] * ff_17[k];

        t_21[k] = pa_z[k] * dg_6[k]
                  + f_2 * fg_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, pg_s_4, pg_4, df_7, df_11, dg_8, dg_10, \
                         dg_12, fg_s_22, fg_s_23, fg_s_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_6 * df_7[k]
                  + pa_z[k] * dg_8[k]
                  + f_2 * fg_s_22[k];

        t_23[k] = -f_5 * pg_s_4[k]
                  + f_4 * pg_4[k]
                  + pa_y[k] * dg_10[k]
                  + f_2 * fg_s_23[k];

        t_24[k] = f_7 * df_11[k]
                  + pa_y[k] * dg_12[k]
                  + f_2 * fg_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pb_y, df_12, df_13, dg_13, dg_15, fg_s_25, \
                         fg_s_26, fg_s_27, ff_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_6 * df_12[k]
                  + pa_y[k] * dg_13[k]
                  + f_2 * fg_s_25[k];

        t_26[k] = f_4 * df_13[k]
                  + f_2 * fg_s_26[k]
                  + pb_y[k] * ff_23[k];

        t_27[k] = pa_y[k] * dg_15[k]
                  + f_2 * fg_s_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_x, fd_s_14, fd_s_15, fd_s_16, fg_s_28, fg_s_29, \
                         fg_s_30, fd_14, fd_15, fd_16, ff_24, ff_26, \
                         ff_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_1 * fd_s_14[k]
                  + f_2 * fg_s_28[k]
                  + f_0 * fd_14[k]
                  + pb_x[k] * ff_24[k];

        t_29[k] = -f_3 * fd_s_15[k]
                  + f_2 * fg_s_29[k]
                  + f_4 * fd_15[k]
                  + pb_x[k] * ff_26[k];

        t_30[k] = -f_3 * fd_s_16[k]
                  + f_2 * fg_s_30[k]
                  + f_4 * fd_16[k]
                  + pb_x[k] * ff_27[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_y, fd_s_15, fd_s_16, fg_s_31, fg_s_32, fg_s_33, \
                         fd_15, fd_16, ff_28, ff_29, ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_1 * fd_s_15[k]
                  + f_2 * fg_s_31[k]
                  + f_0 * fd_15[k]
                  + pb_y[k] * ff_28[k];

        t_32[k] = -f_3 * fd_s_16[k]
                  + f_2 * fg_s_32[k]
                  + f_4 * fd_16[k]
                  + pb_y[k] * ff_29[k];

        t_33[k] = f_2 * fg_s_33[k]
                  + pb_y[k] * ff_30[k];
    }

#pragma omp simd aligned(t_34, pb_z, df_13, fd_s_16, fg_s_34, fd_16, \
                         ff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * df_13[k]
                  - f_1 * fd_s_16[k]
                  + f_2 * fg_s_34[k]
                  + f_0 * fd_16[k]
                  + pb_z[k] * ff_30[k];
    }
}

auto
compute_prim_fg_kinetic_energy_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                  const size_t pb, const size_t pg_s, const size_t pg,
                                  const size_t df, const size_t dg, const size_t fd_s,
                                  const size_t fg_s, const size_t fd, const size_t ff,
                                  const size_t ncols, const double alpha, const double beta,
                                  const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = beta / p;
    const auto f_6 = 2.0 / p;
    const auto f_7 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_s_1 = buffer.data(pg_s + 1);
    const auto *pg_s_4 = buffer.data(pg_s + 4);

    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_18 = buffer.data(dg + 18);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_1 = buffer.data(fg_s + 1);
    const auto *fg_s_2 = buffer.data(fg_s + 2);
    const auto *fg_s_3 = buffer.data(fg_s + 3);
    const auto *fg_s_4 = buffer.data(fg_s + 4);
    const auto *fg_s_5 = buffer.data(fg_s + 5);
    const auto *fg_s_6 = buffer.data(fg_s + 6);
    const auto *fg_s_7 = buffer.data(fg_s + 7);
    const auto *fg_s_8 = buffer.data(fg_s + 8);
    const auto *fg_s_9 = buffer.data(fg_s + 9);
    const auto *fg_s_10 = buffer.data(fg_s + 10);
    const auto *fg_s_11 = buffer.data(fg_s + 11);
    const auto *fg_s_12 = buffer.data(fg_s + 12);
    const auto *fg_s_13 = buffer.data(fg_s + 13);
    const auto *fg_s_14 = buffer.data(fg_s + 14);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_16 = buffer.data(fg_s + 16);
    const auto *fg_s_17 = buffer.data(fg_s + 17);
    const auto *fg_s_18 = buffer.data(fg_s + 18);
    const auto *fg_s_19 = buffer.data(fg_s + 19);
    const auto *fg_s_20 = buffer.data(fg_s + 20);
    const auto *fg_s_21 = buffer.data(fg_s + 21);
    const auto *fg_s_22 = buffer.data(fg_s + 22);
    const auto *fg_s_23 = buffer.data(fg_s + 23);
    const auto *fg_s_24 = buffer.data(fg_s + 24);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_26 = buffer.data(fg_s + 26);
    const auto *fg_s_27 = buffer.data(fg_s + 27);
    const auto *fg_s_28 = buffer.data(fg_s + 28);
    const auto *fg_s_29 = buffer.data(fg_s + 29);
    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_31 = buffer.data(fg_s + 31);
    const auto *fg_s_32 = buffer.data(fg_s + 32);
    const auto *fg_s_33 = buffer.data(fg_s + 33);
    const auto *fg_s_34 = buffer.data(fg_s + 34);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_36 = buffer.data(fg_s + 36);
    const auto *fg_s_37 = buffer.data(fg_s + 37);
    const auto *fg_s_38 = buffer.data(fg_s + 38);
    const auto *fg_s_39 = buffer.data(fg_s + 39);
    const auto *fg_s_40 = buffer.data(fg_s + 40);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_12 = buffer.data(ff + 12);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_14 = buffer.data(ff + 14);
    const auto *ff_15 = buffer.data(ff + 15);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_21 = buffer.data(ff + 21);
    const auto *ff_23 = buffer.data(ff + 23);
    const auto *ff_24 = buffer.data(ff + 24);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_26 = buffer.data(ff + 26);
    const auto *ff_27 = buffer.data(ff + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, df_0, fd_s_0, fg_s_0, fg_s_1, \
                         fg_s_2, fg_s_3, fd_0, ff_0, ff_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * df_0[k]
                 - f_1 * fd_s_0[k]
                 + f_2 * fg_s_0[k]
                 + f_0 * fd_0[k]
                 + pb_x[k] * ff_0[k];

        t_1[k] = f_2 * fg_s_1[k]
                 + pb_y[k] * ff_0[k];

        t_2[k] = f_2 * fg_s_2[k]
                 + pb_z[k] * ff_0[k];

        t_3[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_3[k]
                 + f_4 * fd_0[k]
                 + pb_y[k] * ff_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_y, pb_z, dg_0, fd_s_0, fd_s_1, fg_s_4, fg_s_5, \
                         fg_s_6, fd_0, fd_1, ff_2, ff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_4[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_5[k] = -f_1 * fd_s_1[k]
                 + f_2 * fg_s_5[k]
                 + f_0 * fd_1[k]
                 + pb_z[k] * ff_3[k];

        t_6[k] = pa_y[k] * dg_0[k]
                 + f_2 * fg_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_z, pg_s_1, pg_s_4, pg_1, pg_4, dg_0, dg_3, \
                         dg_4, fg_s_7, fg_s_8, fg_s_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_5 * pg_s_1[k]
                 + f_4 * pg_1[k]
                 + pa_x[k] * dg_3[k]
                 + f_2 * fg_s_7[k];

        t_8[k] = pa_z[k] * dg_0[k]
                 + f_2 * fg_s_8[k];

        t_9[k] = -f_5 * pg_s_4[k]
                 + f_4 * pg_4[k]
                 + pa_x[k] * dg_4[k]
                 + f_2 * fg_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_x, df_4, df_5, df_6, dg_5, dg_6, fg_s_10, \
                         fg_s_11, fg_s_12, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_6 * df_4[k]
                  + pa_x[k] * dg_5[k]
                  + f_2 * fg_s_10[k];

        t_11[k] = f_7 * df_5[k]
                  + pa_x[k] * dg_6[k]
                  + f_2 * fg_s_11[k];

        t_12[k] = f_4 * df_6[k]
                  + f_2 * fg_s_12[k]
                  + pb_x[k] * ff_9[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, df_9, df_10, dg_8, dg_12, dg_13, fg_s_13, \
                         fg_s_14, fg_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_x[k] * dg_8[k]
                  + f_2 * fg_s_13[k];

        t_14[k] = f_6 * df_9[k]
                  + pa_x[k] * dg_12[k]
                  + f_2 * fg_s_14[k];

        t_15[k] = f_7 * df_10[k]
                  + pa_x[k] * dg_13[k]
                  + f_2 * fg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_x, df_13, dg_18, fd_s_8, fg_s_16, fg_s_17, \
                         fg_s_18, fd_8, ff_11, ff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * df_13[k]
                  + f_2 * fg_s_16[k]
                  + pb_x[k] * ff_11[k];

        t_17[k] = pa_x[k] * dg_18[k]
                  + f_2 * fg_s_17[k];

        t_18[k] = -f_1 * fd_s_8[k]
                  + f_2 * fg_s_18[k]
                  + f_0 * fd_8[k]
                  + pb_x[k] * ff_12[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pb_y, pb_z, df_6, fd_s_9, fg_s_19, \
                         fg_s_20, fg_s_21, fg_s_22, fd_9, ff_13, \
                         ff_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_3 * fd_s_9[k]
                  + f_2 * fg_s_19[k]
                  + f_4 * fd_9[k]
                  + pb_x[k] * ff_13[k];

        t_20[k] = f_2 * fg_s_20[k]
                  + pb_x[k] * ff_14[k];

        t_21[k] = f_0 * df_6[k]
                  - f_1 * fd_s_9[k]
                  + f_2 * fg_s_21[k]
                  + f_0 * fd_9[k]
                  + pb_y[k] * ff_14[k];

        t_22[k] = f_2 * fg_s_22[k]
                  + pb_z[k] * ff_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_z, pb_z, dg_8, fd_s_9, fd_s_10, fg_s_23, \
                         fg_s_24, fg_s_25, fd_9, fd_10, ff_15, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_3 * fd_s_9[k]
                  + f_2 * fg_s_23[k]
                  + f_4 * fd_9[k]
                  + pb_z[k] * ff_15[k];

        t_24[k] = -f_1 * fd_s_10[k]
                  + f_2 * fg_s_24[k]
                  + f_0 * fd_10[k]
                  + pb_z[k] * ff_16[k];

        t_25[k] = pa_z[k] * dg_8[k]
                  + f_2 * fg_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pg_s_4, pg_4, df_11, df_12, dg_11, dg_15, \
                         dg_16, fg_s_26, fg_s_27, fg_s_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -f_5 * pg_s_4[k]
                  + f_4 * pg_4[k]
                  + pa_y[k] * dg_11[k]
                  + f_2 * fg_s_26[k];

        t_27[k] = f_6 * df_11[k]
                  + pa_y[k] * dg_15[k]
                  + f_2 * fg_s_27[k];

        t_28[k] = f_7 * df_12[k]
                  + pa_y[k] * dg_16[k]
                  + f_2 * fg_s_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pb_x, pb_y, df_13, dg_18, fd_s_14, fg_s_29, \
                         fg_s_30, fg_s_31, fd_14, ff_20, ff_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_4 * df_13[k]
                  + f_2 * fg_s_29[k]
                  + pb_y[k] * ff_20[k];

        t_30[k] = pa_y[k] * dg_18[k]
                  + f_2 * fg_s_30[k];

        t_31[k] = -f_1 * fd_s_14[k]
                  + f_2 * fg_s_31[k]
                  + f_0 * fd_14[k]
                  + pb_x[k] * ff_21[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, fd_s_15, fd_s_16, fg_s_32, fg_s_33, \
                         fg_s_34, fd_15, fd_16, ff_21, ff_23, ff_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * fg_s_32[k]
                  + pb_y[k] * ff_21[k];

        t_33[k] = -f_3 * fd_s_15[k]
                  + f_2 * fg_s_33[k]
                  + f_4 * fd_15[k]
                  + pb_x[k] * ff_23[k];

        t_34[k] = -f_3 * fd_s_16[k]
                  + f_2 * fg_s_34[k]
                  + f_4 * fd_16[k]
                  + pb_x[k] * ff_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, pb_y, fd_s_15, fg_s_35, fg_s_36, fg_s_37, \
                         fd_15, ff_25, ff_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_2 * fg_s_35[k]
                  + pb_x[k] * ff_25[k];

        t_36[k] = f_2 * fg_s_36[k]
                  + pb_x[k] * ff_27[k];

        t_37[k] = -f_1 * fd_s_15[k]
                  + f_2 * fg_s_37[k]
                  + f_0 * fd_15[k]
                  + pb_y[k] * ff_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_y, pb_z, df_13, fd_s_16, fg_s_38, fg_s_39, \
                         fg_s_40, fd_16, ff_26, ff_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -f_3 * fd_s_16[k]
                  + f_2 * fg_s_38[k]
                  + f_4 * fd_16[k]
                  + pb_y[k] * ff_26[k];

        t_39[k] = f_2 * fg_s_39[k]
                  + pb_y[k] * ff_27[k];

        t_40[k] = f_0 * df_13[k]
                  - f_1 * fd_s_16[k]
                  + f_2 * fg_s_40[k]
                  + f_0 * fd_16[k]
                  + pb_z[k] * ff_27[k];
    }
}

}  // namespace simdkin
