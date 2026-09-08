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

    const auto *pg_s_25 = buffer.data(pg_s + 25);
    const auto *pg_s_44 = buffer.data(pg_s + 44);

    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_44 = buffer.data(pg + 44);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_27 = buffer.data(df + 27);
    const auto *df_29 = buffer.data(df + 29);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_31 = buffer.data(df + 31);
    const auto *df_33 = buffer.data(df + 33);
    const auto *df_35 = buffer.data(df + 35);
    const auto *df_36 = buffer.data(df + 36);
    const auto *df_37 = buffer.data(df + 37);
    const auto *df_38 = buffer.data(df + 38);
    const auto *df_39 = buffer.data(df + 39);
    const auto *df_46 = buffer.data(df + 46);
    const auto *df_47 = buffer.data(df + 47);
    const auto *df_48 = buffer.data(df + 48);
    const auto *df_49 = buffer.data(df + 49);
    const auto *df_50 = buffer.data(df + 50);
    const auto *df_51 = buffer.data(df + 51);
    const auto *df_52 = buffer.data(df + 52);
    const auto *df_53 = buffer.data(df + 53);
    const auto *df_55 = buffer.data(df + 55);
    const auto *df_56 = buffer.data(df + 56);
    const auto *df_57 = buffer.data(df + 57);
    const auto *df_58 = buffer.data(df + 58);
    const auto *df_59 = buffer.data(df + 59);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_79 = buffer.data(dg + 79);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_36 = buffer.data(fd_s + 36);
    const auto *fd_s_37 = buffer.data(fd_s + 37);
    const auto *fd_s_39 = buffer.data(fd_s + 39);
    const auto *fd_s_41 = buffer.data(fd_s + 41);
    const auto *fd_s_44 = buffer.data(fd_s + 44);
    const auto *fd_s_47 = buffer.data(fd_s + 47);
    const auto *fd_s_54 = buffer.data(fd_s + 54);
    const auto *fd_s_56 = buffer.data(fd_s + 56);
    const auto *fd_s_57 = buffer.data(fd_s + 57);
    const auto *fd_s_58 = buffer.data(fd_s + 58);
    const auto *fd_s_59 = buffer.data(fd_s + 59);

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
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_37 = buffer.data(fd + 37);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_41 = buffer.data(fd + 41);
    const auto *fd_44 = buffer.data(fd + 44);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_54 = buffer.data(fd + 54);
    const auto *fd_56 = buffer.data(fd + 56);
    const auto *fd_57 = buffer.data(fd + 57);
    const auto *fd_58 = buffer.data(fd + 58);
    const auto *fd_59 = buffer.data(fd + 59);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_98 = buffer.data(ff + 98);
    const auto *ff_99 = buffer.data(ff + 99);

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

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, df_6, fd_s_0, fg_s_4, fg_s_5, \
                         fg_s_6, fd_0, ff_2, ff_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * fg_s_4[k]
                 + pb_y[k] * ff_2[k];

        t_5[k] = -f_3 * fd_s_0[k]
                 + f_2 * fg_s_5[k]
                 + f_4 * fd_0[k]
                 + pb_z[k] * ff_2[k];

        t_6[k] = f_0 * df_6[k]
                 + f_2 * fg_s_6[k]
                 + pb_x[k] * ff_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, df_9, fg_s_7, fg_s_8, fg_s_9, ff_3, \
                         ff_5, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * fg_s_7[k]
                 + pb_z[k] * ff_3[k];

        t_8[k] = f_2 * fg_s_8[k]
                 + pb_y[k] * ff_5[k];

        t_9[k] = f_0 * df_9[k]
                 + f_2 * fg_s_9[k]
                 + pb_x[k] * ff_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, fd_s_3, fd_s_5, fg_s_10, fg_s_11, \
                         fg_s_12, fd_3, fd_5, ff_6, ff_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * fd_s_3[k]
                  + f_2 * fg_s_10[k]
                  + f_0 * fd_3[k]
                  + pb_y[k] * ff_6[k];

        t_11[k] = f_2 * fg_s_11[k]
                  + pb_z[k] * ff_6[k];

        t_12[k] = -f_3 * fd_s_5[k]
                  + f_2 * fg_s_12[k]
                  + f_4 * fd_5[k]
                  + pb_y[k] * ff_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, dg_0, fd_s_5, fg_s_13, fg_s_14, \
                         fg_s_15, fd_5, ff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * fg_s_13[k]
                  + pb_y[k] * ff_9[k];

        t_14[k] = -f_1 * fd_s_5[k]
                  + f_2 * fg_s_14[k]
                  + f_0 * fd_5[k]
                  + pb_z[k] * ff_9[k];

        t_15[k] = pa_y[k] * dg_0[k]
                  + f_2 * fg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, df_0, df_1, dg_3, fg_s_16, \
                         fg_s_17, fg_s_18, fg_s_19, ff_10, ff_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * df_0[k]
                  + f_2 * fg_s_16[k]
                  + pb_y[k] * ff_10[k];

        t_17[k] = f_2 * fg_s_17[k]
                  + pb_z[k] * ff_10[k];

        t_18[k] = f_5 * df_1[k]
                  + pa_y[k] * dg_3[k]
                  + f_2 * fg_s_18[k];

        t_19[k] = f_2 * fg_s_19[k]
                  + pb_z[k] * ff_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_x, pb_z, df_16, dg_5, fg_s_20, fg_s_21, \
                         fg_s_22, ff_13, ff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * dg_5[k]
                  + f_2 * fg_s_20[k];

        t_21[k] = f_5 * df_16[k]
                  + f_2 * fg_s_21[k]
                  + pb_x[k] * ff_16[k];

        t_22[k] = f_2 * fg_s_22[k]
                  + pb_z[k] * ff_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, pb_x, pg_s_25, pg_25, df_18, dg_9, \
                         dg_25, fg_s_23, fg_s_24, fg_s_25, ff_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_5 * df_18[k]
                  + f_2 * fg_s_23[k]
                  + pb_x[k] * ff_18[k];

        t_24[k] = pa_y[k] * dg_9[k]
                  + f_2 * fg_s_24[k];

        t_25[k] = -f_6 * pg_s_25[k]
                  + f_4 * pg_25[k]
                  + pa_x[k] * dg_25[k]
                  + f_2 * fg_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_y, pb_z, df_9, fd_s_9, fg_s_26, fg_s_27, \
                         fg_s_28, fd_9, ff_16, ff_17, ff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * fg_s_26[k]
                  + pb_z[k] * ff_16[k];

        t_27[k] = -f_3 * fd_s_9[k]
                  + f_2 * fg_s_27[k]
                  + f_4 * fd_9[k]
                  + pb_z[k] * ff_17[k];

        t_28[k] = f_4 * df_9[k]
                  + f_2 * fg_s_28[k]
                  + pb_y[k] * ff_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pb_y, pb_z, df_0, dg_0, dg_14, \
                         fg_s_29, fg_s_30, fg_s_31, fg_s_32, ff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * dg_14[k]
                  + f_2 * fg_s_29[k];

        t_30[k] = pa_z[k] * dg_0[k]
                  + f_2 * fg_s_30[k];

        t_31[k] = f_2 * fg_s_31[k]
                  + pb_y[k] * ff_20[k];

        t_32[k] = f_4 * df_0[k]
                  + f_2 * fg_s_32[k]
                  + pb_z[k] * ff_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_y, df_2, dg_3, dg_5, dg_6, fg_s_33, \
                         fg_s_34, fg_s_35, fg_s_36, ff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * dg_3[k]
                  + f_2 * fg_s_33[k];

        t_34[k] = f_2 * fg_s_34[k]
                  + pb_y[k] * ff_22[k];

        t_35[k] = f_5 * df_2[k]
                  + pa_z[k] * dg_5[k]
                  + f_2 * fg_s_35[k];

        t_36[k] = pa_z[k] * dg_6[k]
                  + f_2 * fg_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, df_27, df_29, fg_s_37, fg_s_38, \
                         fg_s_39, ff_25, ff_27, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_5 * df_27[k]
                  + f_2 * fg_s_37[k]
                  + pb_x[k] * ff_27[k];

        t_38[k] = f_2 * fg_s_38[k]
                  + pb_y[k] * ff_25[k];

        t_39[k] = f_5 * df_29[k]
                  + f_2 * fg_s_39[k]
                  + pb_x[k] * ff_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_z, pb_y, dg_10, fd_s_16, fd_s_17, fg_s_40, \
                         fg_s_41, fg_s_42, fd_16, fd_17, ff_27, ff_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_z[k] * dg_10[k]
                  + f_2 * fg_s_40[k];

        t_41[k] = -f_7 * fd_s_16[k]
                  + f_2 * fg_s_41[k]
                  + f_5 * fd_16[k]
                  + pb_y[k] * ff_27[k];

        t_42[k] = -f_3 * fd_s_17[k]
                  + f_2 * fg_s_42[k]
                  + f_4 * fd_17[k]
                  + pb_y[k] * ff_28[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pb_y, pg_s_44, pg_44, df_30, dg_44, dg_45, \
                         fg_s_43, fg_s_44, fg_s_45, ff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * fg_s_43[k]
                  + pb_y[k] * ff_29[k];

        t_44[k] = -f_6 * pg_s_44[k]
                  + f_4 * pg_44[k]
                  + pa_x[k] * dg_44[k]
                  + f_2 * fg_s_44[k];

        t_45[k] = f_8 * df_30[k]
                  + pa_x[k] * dg_45[k]
                  + f_2 * fg_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_y, pb_z, df_10, df_33, dg_48, \
                         fg_s_46, fg_s_47, fg_s_48, fg_s_49, ff_30, \
                         ff_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_5 * df_10[k]
                  + f_2 * fg_s_46[k]
                  + pb_y[k] * ff_30[k];

        t_47[k] = f_2 * fg_s_47[k]
                  + pb_z[k] * ff_30[k];

        t_48[k] = f_5 * df_33[k]
                  + pa_x[k] * dg_48[k]
                  + f_2 * fg_s_48[k];

        t_49[k] = f_2 * fg_s_49[k]
                  + pb_z[k] * ff_31[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pb_x, pb_z, df_35, df_36, dg_50, fg_s_50, \
                         fg_s_51, fg_s_52, ff_33, ff_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_5 * df_35[k]
                  + pa_x[k] * dg_50[k]
                  + f_2 * fg_s_50[k];

        t_51[k] = f_4 * df_36[k]
                  + f_2 * fg_s_51[k]
                  + pb_x[k] * ff_36[k];

        t_52[k] = f_2 * fg_s_52[k]
                  + pb_z[k] * ff_33[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_x, pb_x, df_38, df_39, dg_55, fg_s_53, fg_s_54, \
                         fg_s_55, ff_38, ff_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_4 * df_38[k]
                  + f_2 * fg_s_53[k]
                  + pb_x[k] * ff_38[k];

        t_54[k] = f_4 * df_39[k]
                  + f_2 * fg_s_54[k]
                  + pb_x[k] * ff_39[k];

        t_55[k] = pa_x[k] * dg_55[k]
                  + f_2 * fg_s_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_z, dg_57, dg_58, dg_59, fg_s_56, \
                         fg_s_57, fg_s_58, fg_s_59, ff_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_2 * fg_s_56[k]
                  + pb_z[k] * ff_36[k];

        t_57[k] = pa_x[k] * dg_57[k]
                  + f_2 * fg_s_57[k];

        t_58[k] = pa_x[k] * dg_58[k]
                  + f_2 * fg_s_58[k];

        t_59[k] = pa_x[k] * dg_59[k]
                  + f_2 * fg_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pa_z, dg_16, dg_18, dg_30, dg_32, \
                         fg_s_60, fg_s_61, fg_s_62, fg_s_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * dg_30[k]
                  + f_2 * fg_s_60[k];

        t_61[k] = pa_z[k] * dg_16[k]
                  + f_2 * fg_s_61[k];

        t_62[k] = pa_y[k] * dg_32[k]
                  + f_2 * fg_s_62[k];

        t_63[k] = pa_z[k] * dg_18[k]
                  + f_2 * fg_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_y, pa_z, pb_y, df_22, dg_21, dg_35, fg_s_64, \
                         fg_s_65, fg_s_66, ff_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_4 * df_22[k]
                  + f_2 * fg_s_64[k]
                  + pb_y[k] * ff_42[k];

        t_65[k] = pa_y[k] * dg_35[k]
                  + f_2 * fg_s_65[k];

        t_66[k] = pa_z[k] * dg_21[k]
                  + f_2 * fg_s_66[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_y, pb_x, df_47, df_48, dg_39, fg_s_67, fg_s_68, \
                         fg_s_69, ff_47, ff_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_4 * df_47[k]
                  + f_2 * fg_s_67[k]
                  + pb_x[k] * ff_47[k];

        t_68[k] = f_4 * df_48[k]
                  + f_2 * fg_s_68[k]
                  + pb_x[k] * ff_48[k];

        t_69[k] = pa_y[k] * dg_39[k]
                  + f_2 * fg_s_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_x, dg_70, dg_71, dg_72, dg_73, \
                         dg_74, fg_s_70, fg_s_71, fg_s_72, fg_s_73, \
                         fg_s_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_x[k] * dg_70[k]
                  + f_2 * fg_s_70[k];

        t_71[k] = pa_x[k] * dg_71[k]
                  + f_2 * fg_s_71[k];

        t_72[k] = pa_x[k] * dg_72[k]
                  + f_2 * fg_s_72[k];

        t_73[k] = pa_x[k] * dg_73[k]
                  + f_2 * fg_s_73[k];

        t_74[k] = pa_x[k] * dg_74[k]
                  + f_2 * fg_s_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_x, pb_y, pb_z, df_20, df_50, dg_75, fg_s_75, \
                         fg_s_76, fg_s_77, ff_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * df_50[k]
                  + pa_x[k] * dg_75[k]
                  + f_2 * fg_s_75[k];

        t_76[k] = f_2 * fg_s_76[k]
                  + pb_y[k] * ff_50[k];

        t_77[k] = f_5 * df_20[k]
                  + f_2 * fg_s_77[k]
                  + pb_z[k] * ff_50[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_x, pb_y, df_53, df_55, dg_78, dg_80, fg_s_78, \
                         fg_s_79, fg_s_80, ff_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_5 * df_53[k]
                  + pa_x[k] * dg_78[k]
                  + f_2 * fg_s_78[k];

        t_79[k] = f_2 * fg_s_79[k]
                  + pb_y[k] * ff_52[k];

        t_80[k] = f_5 * df_55[k]
                  + pa_x[k] * dg_80[k]
                  + f_2 * fg_s_80[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_x, pb_y, df_56, df_57, fg_s_81, fg_s_82, \
                         fg_s_83, ff_55, ff_56, ff_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * df_56[k]
                  + f_2 * fg_s_81[k]
                  + pb_x[k] * ff_56[k];

        t_82[k] = f_4 * df_57[k]
                  + f_2 * fg_s_82[k]
                  + pb_x[k] * ff_57[k];

        t_83[k] = f_2 * fg_s_83[k]
                  + pb_y[k] * ff_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pb_x, df_59, dg_85, dg_86, dg_87, \
                         fg_s_84, fg_s_85, fg_s_86, fg_s_87, ff_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_4 * df_59[k]
                  + f_2 * fg_s_84[k]
                  + pb_x[k] * ff_59[k];

        t_85[k] = pa_x[k] * dg_85[k]
                  + f_2 * fg_s_85[k];

        t_86[k] = pa_x[k] * dg_86[k]
                  + f_2 * fg_s_86[k];

        t_87[k] = pa_x[k] * dg_87[k]
                  + f_2 * fg_s_87[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pa_x, pb_x, pb_y, dg_89, fd_s_36, fg_s_88, fg_s_89, \
                         fg_s_90, fd_36, ff_59, ff_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * fg_s_88[k]
                  + pb_y[k] * ff_59[k];

        t_89[k] = pa_x[k] * dg_89[k]
                  + f_2 * fg_s_89[k];

        t_90[k] = -f_1 * fd_s_36[k]
                  + f_2 * fg_s_90[k]
                  + f_0 * fd_36[k]
                  + pb_x[k] * ff_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_x, pb_z, fd_s_37, fd_s_39, fg_s_91, fg_s_92, \
                         fg_s_93, fd_37, fd_39, ff_60, ff_61, ff_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = -f_7 * fd_s_37[k]
                  + f_2 * fg_s_91[k]
                  + f_5 * fd_37[k]
                  + pb_x[k] * ff_61[k];

        t_92[k] = f_2 * fg_s_92[k]
                  + pb_z[k] * ff_60[k];

        t_93[k] = -f_3 * fd_s_39[k]
                  + f_2 * fg_s_93[k]
                  + f_4 * fd_39[k]
                  + pb_x[k] * ff_63[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_x, pb_z, fd_s_41, fg_s_94, fg_s_95, \
                         fg_s_96, fg_s_97, fd_41, ff_61, ff_65, ff_66, \
                         ff_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_2 * fg_s_94[k]
                  + pb_z[k] * ff_61[k];

        t_95[k] = -f_3 * fd_s_41[k]
                  + f_2 * fg_s_95[k]
                  + f_4 * fd_41[k]
                  + pb_x[k] * ff_65[k];

        t_96[k] = f_2 * fg_s_96[k]
                  + pb_x[k] * ff_66[k];

        t_97[k] = f_2 * fg_s_97[k]
                  + pb_x[k] * ff_67[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pb_x, pb_y, df_36, fd_s_39, fg_s_98, fg_s_99, \
                         fg_s_100, fd_39, ff_66, ff_68, ff_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_2 * fg_s_98[k]
                  + pb_x[k] * ff_68[k];

        t_99[k] = f_2 * fg_s_99[k]
                  + pb_x[k] * ff_69[k];

        t_100[k] = f_0 * df_36[k]
                   - f_1 * fd_s_39[k]
                   + f_2 * fg_s_100[k]
                   + f_0 * fd_39[k]
                   + pb_y[k] * ff_66[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pb_y, pb_z, df_39, fd_s_39, fg_s_101, fg_s_102, \
                         fg_s_103, fd_39, ff_66, ff_67, ff_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_2 * fg_s_101[k]
                   + pb_z[k] * ff_66[k];

        t_102[k] = -f_3 * fd_s_39[k]
                   + f_2 * fg_s_102[k]
                   + f_4 * fd_39[k]
                   + pb_z[k] * ff_67[k];

        t_103[k] = f_0 * df_39[k]
                   + f_2 * fg_s_103[k]
                   + pb_y[k] * ff_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pb_z, dg_45, dg_46, fd_s_41, fg_s_104, \
                         fg_s_105, fg_s_106, fd_41, ff_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -f_1 * fd_s_41[k]
                   + f_2 * fg_s_104[k]
                   + f_0 * fd_41[k]
                   + pb_z[k] * ff_69[k];

        t_105[k] = pa_z[k] * dg_45[k]
                   + f_2 * fg_s_105[k];

        t_106[k] = pa_z[k] * dg_46[k]
                   + f_2 * fg_s_106[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pa_z, pb_x, df_31, dg_48, dg_49, fd_s_44, \
                         fg_s_107, fg_s_108, fg_s_109, fd_44, ff_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -f_7 * fd_s_44[k]
                   + f_2 * fg_s_107[k]
                   + f_5 * fd_44[k]
                   + pb_x[k] * ff_72[k];

        t_108[k] = pa_z[k] * dg_48[k]
                   + f_2 * fg_s_108[k];

        t_109[k] = f_4 * df_31[k]
                   + pa_z[k] * dg_49[k]
                   + f_2 * fg_s_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_x, fd_s_47, fg_s_110, fg_s_111, \
                         fg_s_112, fg_s_113, fd_47, ff_75, ff_76, ff_77, \
                         ff_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -f_3 * fd_s_47[k]
                   + f_2 * fg_s_110[k]
                   + f_4 * fd_47[k]
                   + pb_x[k] * ff_75[k];

        t_111[k] = f_2 * fg_s_111[k]
                   + pb_x[k] * ff_76[k];

        t_112[k] = f_2 * fg_s_112[k]
                   + pb_x[k] * ff_77[k];

        t_113[k] = f_2 * fg_s_113[k]
                   + pb_x[k] * ff_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_z, pb_x, pb_z, df_36, dg_55, fg_s_114, \
                         fg_s_115, fg_s_116, ff_76, ff_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_2 * fg_s_114[k]
                   + pb_x[k] * ff_79[k];

        t_115[k] = pa_z[k] * dg_55[k]
                   + f_2 * fg_s_115[k];

        t_116[k] = f_4 * df_36[k]
                   + f_2 * fg_s_116[k]
                   + pb_z[k] * ff_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_y, pa_z, pb_y, pg_s_44, pg_44, df_37, df_49, \
                         dg_57, dg_74, fg_s_117, fg_s_118, fg_s_119, \
                         ff_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * df_37[k]
                   + pa_z[k] * dg_57[k]
                   + f_2 * fg_s_117[k];

        t_118[k] = f_5 * df_49[k]
                   + f_2 * fg_s_118[k]
                   + pb_y[k] * ff_79[k];

        t_119[k] = -f_6 * pg_s_44[k]
                   + f_4 * pg_44[k]
                   + pa_y[k] * dg_74[k]
                   + f_2 * fg_s_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, df_50, df_51, dg_75, dg_76, dg_77, \
                         dg_78, fg_s_120, fg_s_121, fg_s_122, \
                         fg_s_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_y[k] * dg_75[k]
                   + f_2 * fg_s_120[k];

        t_121[k] = f_4 * df_50[k]
                   + pa_y[k] * dg_76[k]
                   + f_2 * fg_s_121[k];

        t_122[k] = pa_y[k] * dg_77[k]
                   + f_2 * fg_s_122[k];

        t_123[k] = f_5 * df_51[k]
                   + pa_y[k] * dg_78[k]
                   + f_2 * fg_s_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pb_x, df_52, dg_79, dg_80, \
                         fg_s_124, fg_s_125, fg_s_126, fg_s_127, ff_86, \
                         ff_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_4 * df_52[k]
                   + pa_y[k] * dg_79[k]
                   + f_2 * fg_s_124[k];

        t_125[k] = pa_y[k] * dg_80[k]
                   + f_2 * fg_s_125[k];

        t_126[k] = f_2 * fg_s_126[k]
                   + pb_x[k] * ff_86[k];

        t_127[k] = f_2 * fg_s_127[k]
                   + pb_x[k] * ff_87[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pa_y, pb_x, df_56, dg_85, fg_s_128, fg_s_129, \
                         fg_s_130, ff_88, ff_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_2 * fg_s_128[k]
                   + pb_x[k] * ff_88[k];

        t_129[k] = f_2 * fg_s_129[k]
                   + pb_x[k] * ff_89[k];

        t_130[k] = f_8 * df_56[k]
                   + pa_y[k] * dg_85[k]
                   + f_2 * fg_s_130[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, pa_y, pb_y, pb_z, df_46, df_58, df_59, dg_87, \
                         fg_s_131, fg_s_132, fg_s_133, ff_86, ff_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_5 * df_46[k]
                   + f_2 * fg_s_131[k]
                   + pb_z[k] * ff_86[k];

        t_132[k] = f_5 * df_58[k]
                   + pa_y[k] * dg_87[k]
                   + f_2 * fg_s_132[k];

        t_133[k] = f_4 * df_59[k]
                   + f_2 * fg_s_133[k]
                   + pb_y[k] * ff_89[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pa_y, pb_x, pb_y, dg_89, fd_s_54, fg_s_134, \
                         fg_s_135, fg_s_136, fd_54, ff_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_y[k] * dg_89[k]
                   + f_2 * fg_s_134[k];

        t_135[k] = -f_1 * fd_s_54[k]
                   + f_2 * fg_s_135[k]
                   + f_0 * fd_54[k]
                   + pb_x[k] * ff_90[k];

        t_136[k] = f_2 * fg_s_136[k]
                   + pb_y[k] * ff_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pb_x, pb_y, fd_s_56, fd_s_57, fg_s_137, \
                         fg_s_138, fg_s_139, fd_56, fd_57, ff_92, \
                         ff_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = -f_7 * fd_s_56[k]
                   + f_2 * fg_s_137[k]
                   + f_5 * fd_56[k]
                   + pb_x[k] * ff_92[k];

        t_138[k] = -f_3 * fd_s_57[k]
                   + f_2 * fg_s_138[k]
                   + f_4 * fd_57[k]
                   + pb_x[k] * ff_93[k];

        t_139[k] = f_2 * fg_s_139[k]
                   + pb_y[k] * ff_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pb_x, fd_s_59, fg_s_140, fg_s_141, \
                         fg_s_142, fg_s_143, fd_59, ff_95, ff_96, ff_97, \
                         ff_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -f_3 * fd_s_59[k]
                   + f_2 * fg_s_140[k]
                   + f_4 * fd_59[k]
                   + pb_x[k] * ff_95[k];

        t_141[k] = f_2 * fg_s_141[k]
                   + pb_x[k] * ff_96[k];

        t_142[k] = f_2 * fg_s_142[k]
                   + pb_x[k] * ff_97[k];

        t_143[k] = f_2 * fg_s_143[k]
                   + pb_x[k] * ff_98[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pb_x, pb_y, fd_s_57, fd_s_58, fg_s_144, \
                         fg_s_145, fg_s_146, fd_57, fd_58, ff_96, ff_97, \
                         ff_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_2 * fg_s_144[k]
                   + pb_x[k] * ff_99[k];

        t_145[k] = -f_1 * fd_s_57[k]
                   + f_2 * fg_s_145[k]
                   + f_0 * fd_57[k]
                   + pb_y[k] * ff_96[k];

        t_146[k] = -f_7 * fd_s_58[k]
                   + f_2 * fg_s_146[k]
                   + f_5 * fd_58[k]
                   + pb_y[k] * ff_97[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pb_y, pb_z, df_59, fd_s_59, fg_s_147, fg_s_148, \
                         fg_s_149, fd_59, ff_98, ff_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -f_3 * fd_s_59[k]
                   + f_2 * fg_s_147[k]
                   + f_4 * fd_59[k]
                   + pb_y[k] * ff_98[k];

        t_148[k] = f_2 * fg_s_148[k]
                   + pb_y[k] * ff_99[k];

        t_149[k] = f_0 * df_59[k]
                   - f_1 * fd_s_59[k]
                   + f_2 * fg_s_149[k]
                   + f_0 * fd_59[k]
                   + pb_z[k] * ff_99[k];
    }
}

}  // namespace simdkin
