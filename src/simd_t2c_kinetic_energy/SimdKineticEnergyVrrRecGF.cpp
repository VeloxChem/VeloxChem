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


#include "SimdKineticEnergyVrrRecGF.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_gf_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t df_s, const size_t df,
                                 const size_t fd, const size_t ff, const size_t gp_s,
                                 const size_t gf_s, const size_t gp, const size_t gd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 2.0 * beta / p;
    const auto f_7 = alpha / p;
    const auto f_8 = beta / p;

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

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_16 = buffer.data(df_s + 16);
    const auto *df_s_29 = buffer.data(df_s + 29);
    const auto *df_s_36 = buffer.data(df_s + 36);
    const auto *df_s_49 = buffer.data(df_s + 49);
    const auto *df_s_59 = buffer.data(df_s + 59);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_29 = buffer.data(df + 29);
    const auto *df_36 = buffer.data(df + 36);
    const auto *df_49 = buffer.data(df + 49);
    const auto *df_59 = buffer.data(df + 59);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_41 = buffer.data(fd + 41);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_46 = buffer.data(fd + 46);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_52 = buffer.data(fd + 52);
    const auto *fd_53 = buffer.data(fd + 53);
    const auto *fd_54 = buffer.data(fd + 54);
    const auto *fd_57 = buffer.data(fd + 57);
    const auto *fd_59 = buffer.data(fd + 59);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_3 = buffer.data(ff + 3);
    const auto *ff_5 = buffer.data(ff + 5);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_11 = buffer.data(ff + 11);
    const auto *ff_13 = buffer.data(ff + 13);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_25 = buffer.data(ff + 25);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_31 = buffer.data(ff + 31);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_91 = buffer.data(ff + 91);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_99 = buffer.data(ff + 99);

    const auto *gp_s_0 = buffer.data(gp_s + 0);
    const auto *gp_s_1 = buffer.data(gp_s + 1);
    const auto *gp_s_2 = buffer.data(gp_s + 2);
    const auto *gp_s_8 = buffer.data(gp_s + 8);
    const auto *gp_s_11 = buffer.data(gp_s + 11);
    const auto *gp_s_16 = buffer.data(gp_s + 16);
    const auto *gp_s_17 = buffer.data(gp_s + 17);
    const auto *gp_s_30 = buffer.data(gp_s + 30);
    const auto *gp_s_31 = buffer.data(gp_s + 31);
    const auto *gp_s_32 = buffer.data(gp_s + 32);
    const auto *gp_s_35 = buffer.data(gp_s + 35);
    const auto *gp_s_36 = buffer.data(gp_s + 36);
    const auto *gp_s_37 = buffer.data(gp_s + 37);
    const auto *gp_s_38 = buffer.data(gp_s + 38);
    const auto *gp_s_42 = buffer.data(gp_s + 42);
    const auto *gp_s_43 = buffer.data(gp_s + 43);
    const auto *gp_s_44 = buffer.data(gp_s + 44);

    const auto *gf_s_0 = buffer.data(gf_s + 0);
    const auto *gf_s_1 = buffer.data(gf_s + 1);
    const auto *gf_s_2 = buffer.data(gf_s + 2);
    const auto *gf_s_3 = buffer.data(gf_s + 3);
    const auto *gf_s_4 = buffer.data(gf_s + 4);
    const auto *gf_s_5 = buffer.data(gf_s + 5);
    const auto *gf_s_6 = buffer.data(gf_s + 6);
    const auto *gf_s_7 = buffer.data(gf_s + 7);
    const auto *gf_s_8 = buffer.data(gf_s + 8);
    const auto *gf_s_9 = buffer.data(gf_s + 9);
    const auto *gf_s_10 = buffer.data(gf_s + 10);
    const auto *gf_s_11 = buffer.data(gf_s + 11);
    const auto *gf_s_12 = buffer.data(gf_s + 12);
    const auto *gf_s_13 = buffer.data(gf_s + 13);
    const auto *gf_s_14 = buffer.data(gf_s + 14);
    const auto *gf_s_15 = buffer.data(gf_s + 15);
    const auto *gf_s_16 = buffer.data(gf_s + 16);
    const auto *gf_s_17 = buffer.data(gf_s + 17);
    const auto *gf_s_18 = buffer.data(gf_s + 18);
    const auto *gf_s_19 = buffer.data(gf_s + 19);
    const auto *gf_s_20 = buffer.data(gf_s + 20);
    const auto *gf_s_21 = buffer.data(gf_s + 21);
    const auto *gf_s_22 = buffer.data(gf_s + 22);
    const auto *gf_s_23 = buffer.data(gf_s + 23);
    const auto *gf_s_24 = buffer.data(gf_s + 24);
    const auto *gf_s_25 = buffer.data(gf_s + 25);
    const auto *gf_s_26 = buffer.data(gf_s + 26);
    const auto *gf_s_27 = buffer.data(gf_s + 27);
    const auto *gf_s_28 = buffer.data(gf_s + 28);
    const auto *gf_s_29 = buffer.data(gf_s + 29);
    const auto *gf_s_30 = buffer.data(gf_s + 30);
    const auto *gf_s_31 = buffer.data(gf_s + 31);
    const auto *gf_s_32 = buffer.data(gf_s + 32);
    const auto *gf_s_33 = buffer.data(gf_s + 33);
    const auto *gf_s_34 = buffer.data(gf_s + 34);
    const auto *gf_s_35 = buffer.data(gf_s + 35);
    const auto *gf_s_36 = buffer.data(gf_s + 36);
    const auto *gf_s_37 = buffer.data(gf_s + 37);
    const auto *gf_s_38 = buffer.data(gf_s + 38);
    const auto *gf_s_39 = buffer.data(gf_s + 39);
    const auto *gf_s_40 = buffer.data(gf_s + 40);
    const auto *gf_s_41 = buffer.data(gf_s + 41);
    const auto *gf_s_42 = buffer.data(gf_s + 42);
    const auto *gf_s_43 = buffer.data(gf_s + 43);
    const auto *gf_s_44 = buffer.data(gf_s + 44);
    const auto *gf_s_45 = buffer.data(gf_s + 45);
    const auto *gf_s_46 = buffer.data(gf_s + 46);
    const auto *gf_s_47 = buffer.data(gf_s + 47);
    const auto *gf_s_48 = buffer.data(gf_s + 48);
    const auto *gf_s_49 = buffer.data(gf_s + 49);
    const auto *gf_s_50 = buffer.data(gf_s + 50);
    const auto *gf_s_51 = buffer.data(gf_s + 51);
    const auto *gf_s_52 = buffer.data(gf_s + 52);
    const auto *gf_s_53 = buffer.data(gf_s + 53);
    const auto *gf_s_54 = buffer.data(gf_s + 54);
    const auto *gf_s_55 = buffer.data(gf_s + 55);
    const auto *gf_s_56 = buffer.data(gf_s + 56);
    const auto *gf_s_57 = buffer.data(gf_s + 57);
    const auto *gf_s_58 = buffer.data(gf_s + 58);
    const auto *gf_s_59 = buffer.data(gf_s + 59);
    const auto *gf_s_60 = buffer.data(gf_s + 60);
    const auto *gf_s_61 = buffer.data(gf_s + 61);
    const auto *gf_s_62 = buffer.data(gf_s + 62);
    const auto *gf_s_63 = buffer.data(gf_s + 63);
    const auto *gf_s_64 = buffer.data(gf_s + 64);
    const auto *gf_s_65 = buffer.data(gf_s + 65);
    const auto *gf_s_66 = buffer.data(gf_s + 66);
    const auto *gf_s_67 = buffer.data(gf_s + 67);
    const auto *gf_s_68 = buffer.data(gf_s + 68);
    const auto *gf_s_69 = buffer.data(gf_s + 69);
    const auto *gf_s_70 = buffer.data(gf_s + 70);
    const auto *gf_s_71 = buffer.data(gf_s + 71);
    const auto *gf_s_72 = buffer.data(gf_s + 72);
    const auto *gf_s_73 = buffer.data(gf_s + 73);
    const auto *gf_s_74 = buffer.data(gf_s + 74);
    const auto *gf_s_75 = buffer.data(gf_s + 75);
    const auto *gf_s_76 = buffer.data(gf_s + 76);
    const auto *gf_s_77 = buffer.data(gf_s + 77);
    const auto *gf_s_78 = buffer.data(gf_s + 78);
    const auto *gf_s_79 = buffer.data(gf_s + 79);
    const auto *gf_s_80 = buffer.data(gf_s + 80);
    const auto *gf_s_81 = buffer.data(gf_s + 81);
    const auto *gf_s_82 = buffer.data(gf_s + 82);
    const auto *gf_s_83 = buffer.data(gf_s + 83);
    const auto *gf_s_84 = buffer.data(gf_s + 84);
    const auto *gf_s_85 = buffer.data(gf_s + 85);
    const auto *gf_s_86 = buffer.data(gf_s + 86);
    const auto *gf_s_87 = buffer.data(gf_s + 87);
    const auto *gf_s_88 = buffer.data(gf_s + 88);
    const auto *gf_s_89 = buffer.data(gf_s + 89);
    const auto *gf_s_90 = buffer.data(gf_s + 90);
    const auto *gf_s_91 = buffer.data(gf_s + 91);
    const auto *gf_s_92 = buffer.data(gf_s + 92);
    const auto *gf_s_93 = buffer.data(gf_s + 93);
    const auto *gf_s_94 = buffer.data(gf_s + 94);
    const auto *gf_s_95 = buffer.data(gf_s + 95);
    const auto *gf_s_96 = buffer.data(gf_s + 96);
    const auto *gf_s_97 = buffer.data(gf_s + 97);
    const auto *gf_s_98 = buffer.data(gf_s + 98);
    const auto *gf_s_99 = buffer.data(gf_s + 99);
    const auto *gf_s_100 = buffer.data(gf_s + 100);
    const auto *gf_s_101 = buffer.data(gf_s + 101);
    const auto *gf_s_102 = buffer.data(gf_s + 102);
    const auto *gf_s_103 = buffer.data(gf_s + 103);
    const auto *gf_s_104 = buffer.data(gf_s + 104);
    const auto *gf_s_105 = buffer.data(gf_s + 105);
    const auto *gf_s_106 = buffer.data(gf_s + 106);
    const auto *gf_s_107 = buffer.data(gf_s + 107);
    const auto *gf_s_108 = buffer.data(gf_s + 108);
    const auto *gf_s_109 = buffer.data(gf_s + 109);
    const auto *gf_s_110 = buffer.data(gf_s + 110);
    const auto *gf_s_111 = buffer.data(gf_s + 111);
    const auto *gf_s_112 = buffer.data(gf_s + 112);
    const auto *gf_s_113 = buffer.data(gf_s + 113);
    const auto *gf_s_114 = buffer.data(gf_s + 114);
    const auto *gf_s_115 = buffer.data(gf_s + 115);
    const auto *gf_s_116 = buffer.data(gf_s + 116);
    const auto *gf_s_117 = buffer.data(gf_s + 117);
    const auto *gf_s_118 = buffer.data(gf_s + 118);
    const auto *gf_s_119 = buffer.data(gf_s + 119);
    const auto *gf_s_120 = buffer.data(gf_s + 120);
    const auto *gf_s_121 = buffer.data(gf_s + 121);
    const auto *gf_s_122 = buffer.data(gf_s + 122);
    const auto *gf_s_123 = buffer.data(gf_s + 123);
    const auto *gf_s_124 = buffer.data(gf_s + 124);
    const auto *gf_s_125 = buffer.data(gf_s + 125);
    const auto *gf_s_126 = buffer.data(gf_s + 126);
    const auto *gf_s_127 = buffer.data(gf_s + 127);
    const auto *gf_s_128 = buffer.data(gf_s + 128);
    const auto *gf_s_129 = buffer.data(gf_s + 129);
    const auto *gf_s_130 = buffer.data(gf_s + 130);
    const auto *gf_s_131 = buffer.data(gf_s + 131);
    const auto *gf_s_132 = buffer.data(gf_s + 132);
    const auto *gf_s_133 = buffer.data(gf_s + 133);
    const auto *gf_s_134 = buffer.data(gf_s + 134);
    const auto *gf_s_135 = buffer.data(gf_s + 135);
    const auto *gf_s_136 = buffer.data(gf_s + 136);
    const auto *gf_s_137 = buffer.data(gf_s + 137);
    const auto *gf_s_138 = buffer.data(gf_s + 138);
    const auto *gf_s_139 = buffer.data(gf_s + 139);
    const auto *gf_s_140 = buffer.data(gf_s + 140);
    const auto *gf_s_141 = buffer.data(gf_s + 141);
    const auto *gf_s_142 = buffer.data(gf_s + 142);
    const auto *gf_s_143 = buffer.data(gf_s + 143);
    const auto *gf_s_144 = buffer.data(gf_s + 144);
    const auto *gf_s_145 = buffer.data(gf_s + 145);
    const auto *gf_s_146 = buffer.data(gf_s + 146);
    const auto *gf_s_147 = buffer.data(gf_s + 147);
    const auto *gf_s_148 = buffer.data(gf_s + 148);
    const auto *gf_s_149 = buffer.data(gf_s + 149);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_61 = buffer.data(gd + 61);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_64 = buffer.data(gd + 64);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_68 = buffer.data(gd + 68);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_73 = buffer.data(gd + 73);
    const auto *gd_74 = buffer.data(gd + 74);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_86 = buffer.data(gd + 86);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_88 = buffer.data(gd + 88);
    const auto *gd_89 = buffer.data(gd + 89);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, fd_0, gp_s_0, gf_s_0, gf_s_1, \
                         gf_s_2, gp_0, gd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fd_0[k]
                 - f_1 * gp_s_0[k]
                 + f_2 * gf_s_0[k]
                 + f_3 * gp_0[k]
                 + pb_x[k] * gd_0[k];

        t_1[k] = f_2 * gf_s_1[k]
                 + pb_y[k] * gd_0[k];

        t_2[k] = f_2 * gf_s_2[k]
                 + pb_z[k] * gd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, fd_3, fd_5, gf_s_3, gf_s_4, gf_s_5, gd_2, \
                         gd_3, gd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * fd_3[k]
                 + f_2 * gf_s_3[k]
                 + pb_x[k] * gd_3[k];

        t_4[k] = f_2 * gf_s_4[k]
                 + pb_y[k] * gd_2[k];

        t_5[k] = f_0 * fd_5[k]
                 + f_2 * gf_s_5[k]
                 + pb_x[k] * gd_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pb_z, gp_s_1, gp_s_2, gf_s_6, gf_s_7, \
                         gf_s_8, gf_s_9, gp_1, gp_2, gd_3, gd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * gp_s_1[k]
                 + f_2 * gf_s_6[k]
                 + f_3 * gp_1[k]
                 + pb_y[k] * gd_3[k];

        t_7[k] = f_2 * gf_s_7[k]
                 + pb_z[k] * gd_3[k];

        t_8[k] = f_2 * gf_s_8[k]
                 + pb_y[k] * gd_5[k];

        t_9[k] = -f_1 * gp_s_2[k]
                 + f_2 * gf_s_9[k]
                 + f_3 * gp_2[k]
                 + pb_z[k] * gd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_y, pb_z, fd_0, ff_0, gf_s_10, gf_s_11, \
                         gf_s_12, gd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * ff_0[k]
                  + f_2 * gf_s_10[k];

        t_11[k] = f_4 * fd_0[k]
                  + f_2 * gf_s_11[k]
                  + pb_y[k] * gd_6[k];

        t_12[k] = f_2 * gf_s_12[k]
                  + pb_z[k] * gd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, pb_z, fd_9, ff_5, gf_s_13, gf_s_14, \
                         gf_s_15, gd_7, gd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * fd_9[k]
                  + f_2 * gf_s_13[k]
                  + pb_x[k] * gd_9[k];

        t_14[k] = f_2 * gf_s_14[k]
                  + pb_z[k] * gd_7[k];

        t_15[k] = pa_y[k] * ff_5[k]
                  + f_2 * gf_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_y, pb_z, df_s_16, df_16, fd_5, ff_16, \
                         gf_s_16, gf_s_17, gf_s_18, gd_9, gd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_6 * df_s_16[k]
                  + f_3 * df_16[k]
                  + pa_x[k] * ff_16[k]
                  + f_2 * gf_s_16[k];

        t_17[k] = f_2 * gf_s_17[k]
                  + pb_z[k] * gd_9[k];

        t_18[k] = f_4 * fd_5[k]
                  + f_2 * gf_s_18[k]
                  + pb_y[k] * gd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_y, pb_z, fd_0, ff_0, ff_9, \
                         gf_s_19, gf_s_20, gf_s_21, gf_s_22, gd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * ff_9[k]
                  + f_2 * gf_s_19[k];

        t_20[k] = pa_z[k] * ff_0[k]
                  + f_2 * gf_s_20[k];

        t_21[k] = f_2 * gf_s_21[k]
                  + pb_y[k] * gd_12[k];

        t_22[k] = f_4 * fd_0[k]
                  + f_2 * gf_s_22[k]
                  + pb_z[k] * gd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_x, pb_y, fd_17, ff_3, ff_6, gf_s_23, \
                         gf_s_24, gf_s_25, gf_s_26, gd_14, gd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_z[k] * ff_3[k]
                  + f_2 * gf_s_23[k];

        t_24[k] = f_2 * gf_s_24[k]
                  + pb_y[k] * gd_14[k];

        t_25[k] = f_5 * fd_17[k]
                  + f_2 * gf_s_25[k]
                  + pb_x[k] * gd_17[k];

        t_26[k] = pa_z[k] * ff_6[k]
                  + f_2 * gf_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pb_y, df_s_29, df_29, ff_29, gp_s_8, gf_s_27, \
                         gf_s_28, gf_s_29, gp_8, gd_16, gd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_7 * gp_s_8[k]
                  + f_2 * gf_s_27[k]
                  + f_4 * gp_8[k]
                  + pb_y[k] * gd_16[k];

        t_28[k] = f_2 * gf_s_28[k]
                  + pb_y[k] * gd_17[k];

        t_29[k] = -f_6 * df_s_29[k]
                  + f_3 * df_29[k]
                  + pa_x[k] * ff_29[k]
                  + f_2 * gf_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_y, pb_z, df_s_0, df_0, fd_6, ff_10, \
                         gf_s_30, gf_s_31, gf_s_32, gd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_8 * df_s_0[k]
                  + f_4 * df_0[k]
                  + pa_y[k] * ff_10[k]
                  + f_2 * gf_s_30[k];

        t_31[k] = f_3 * fd_6[k]
                  + f_2 * gf_s_31[k]
                  + pb_y[k] * gd_18[k];

        t_32[k] = f_2 * gf_s_32[k]
                  + pb_z[k] * gd_18[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_z, fd_21, fd_23, gf_s_33, gf_s_34, \
                         gf_s_35, gd_19, gd_21, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * fd_21[k]
                  + f_2 * gf_s_33[k]
                  + pb_x[k] * gd_21[k];

        t_34[k] = f_2 * gf_s_34[k]
                  + pb_z[k] * gd_19[k];

        t_35[k] = f_3 * fd_23[k]
                  + f_2 * gf_s_35[k]
                  + pb_x[k] * gd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_y, pb_z, df_s_36, df_36, fd_11, ff_36, \
                         gf_s_36, gf_s_37, gf_s_38, gd_21, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_8 * df_s_36[k]
                  + f_4 * df_36[k]
                  + pa_x[k] * ff_36[k]
                  + f_2 * gf_s_36[k];

        t_37[k] = f_2 * gf_s_37[k]
                  + pb_z[k] * gd_21[k];

        t_38[k] = f_3 * fd_11[k]
                  + f_2 * gf_s_38[k]
                  + pb_y[k] * gd_23[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pa_z, pb_z, ff_11, ff_20, gp_s_11, gf_s_39, \
                         gf_s_40, gf_s_41, gp_11, gd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_1 * gp_s_11[k]
                  + f_2 * gf_s_39[k]
                  + f_3 * gp_11[k]
                  + pb_z[k] * gd_23[k];

        t_40[k] = pa_y[k] * ff_20[k]
                  + f_2 * gf_s_40[k];

        t_41[k] = pa_z[k] * ff_11[k]
                  + f_2 * gf_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pa_z, pb_x, fd_28, ff_13, ff_22, ff_25, \
                         gf_s_42, gf_s_43, gf_s_44, gf_s_45, gd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * ff_22[k]
                  + f_2 * gf_s_42[k];

        t_43[k] = pa_z[k] * ff_13[k]
                  + f_2 * gf_s_43[k];

        t_44[k] = f_3 * fd_28[k]
                  + f_2 * gf_s_44[k]
                  + pb_x[k] * gd_28[k];

        t_45[k] = pa_y[k] * ff_25[k]
                  + f_2 * gf_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_z, pb_y, pb_z, fd_9, fd_17, ff_16, gf_s_46, \
                         gf_s_47, gf_s_48, gd_27, gd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * ff_16[k]
                  + f_2 * gf_s_46[k];

        t_47[k] = f_4 * fd_9[k]
                  + f_2 * gf_s_47[k]
                  + pb_z[k] * gd_27[k];

        t_48[k] = f_4 * fd_17[k]
                  + f_2 * gf_s_48[k]
                  + pb_y[k] * gd_29[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pa_z, pb_y, df_s_0, df_0, ff_20, ff_29, \
                         gf_s_49, gf_s_50, gf_s_51, gd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_y[k] * ff_29[k]
                  + f_2 * gf_s_49[k];

        t_50[k] = -f_8 * df_s_0[k]
                  + f_4 * df_0[k]
                  + pa_z[k] * ff_20[k]
                  + f_2 * gf_s_50[k];

        t_51[k] = f_2 * gf_s_51[k]
                  + pb_y[k] * gd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_y, pb_z, fd_12, fd_33, gf_s_52, gf_s_53, \
                         gf_s_54, gd_30, gd_32, gd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * fd_12[k]
                  + f_2 * gf_s_52[k]
                  + pb_z[k] * gd_30[k];

        t_53[k] = f_3 * fd_33[k]
                  + f_2 * gf_s_53[k]
                  + pb_x[k] * gd_33[k];

        t_54[k] = f_2 * gf_s_54[k]
                  + pb_y[k] * gd_32[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pb_y, fd_35, gp_s_16, gp_s_17, gf_s_55, \
                         gf_s_56, gf_s_57, gp_16, gp_17, gd_33, gd_34, \
                         gd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_3 * fd_35[k]
                  + f_2 * gf_s_55[k]
                  + pb_x[k] * gd_35[k];

        t_56[k] = -f_1 * gp_s_16[k]
                  + f_2 * gf_s_56[k]
                  + f_3 * gp_16[k]
                  + pb_y[k] * gd_33[k];

        t_57[k] = -f_7 * gp_s_17[k]
                  + f_2 * gf_s_57[k]
                  + f_4 * gp_17[k]
                  + pb_y[k] * gd_34[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pb_y, df_s_59, df_59, fd_36, ff_59, ff_60, \
                         gf_s_58, gf_s_59, gf_s_60, gd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * gf_s_58[k]
                  + pb_y[k] * gd_35[k];

        t_59[k] = -f_8 * df_s_59[k]
                  + f_4 * df_59[k]
                  + pa_x[k] * ff_59[k]
                  + f_2 * gf_s_59[k];

        t_60[k] = f_5 * fd_36[k]
                  + pa_x[k] * ff_60[k]
                  + f_2 * gf_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_x, pb_y, pb_z, fd_18, fd_39, gf_s_61, \
                         gf_s_62, gf_s_63, gf_s_64, gd_36, gd_37, \
                         gd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * fd_18[k]
                  + f_2 * gf_s_61[k]
                  + pb_y[k] * gd_36[k];

        t_62[k] = f_2 * gf_s_62[k]
                  + pb_z[k] * gd_36[k];

        t_63[k] = f_4 * fd_39[k]
                  + f_2 * gf_s_63[k]
                  + pb_x[k] * gd_39[k];

        t_64[k] = f_2 * gf_s_64[k]
                  + pb_z[k] * gd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_x, pb_z, fd_41, ff_66, ff_68, \
                         gf_s_65, gf_s_66, gf_s_67, gf_s_68, gd_39, \
                         gd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_4 * fd_41[k]
                  + f_2 * gf_s_65[k]
                  + pb_x[k] * gd_41[k];

        t_66[k] = pa_x[k] * ff_66[k]
                  + f_2 * gf_s_66[k];

        t_67[k] = f_2 * gf_s_67[k]
                  + pb_z[k] * gd_39[k];

        t_68[k] = pa_x[k] * ff_68[k]
                  + f_2 * gf_s_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pa_z, pb_z, fd_18, ff_30, ff_31, ff_69, \
                         gf_s_69, gf_s_70, gf_s_71, gf_s_72, gd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_x[k] * ff_69[k]
                  + f_2 * gf_s_69[k];

        t_70[k] = pa_z[k] * ff_30[k]
                  + f_2 * gf_s_70[k];

        t_71[k] = pa_z[k] * ff_31[k]
                  + f_2 * gf_s_71[k];

        t_72[k] = f_4 * fd_18[k]
                  + f_2 * gf_s_72[k]
                  + pb_z[k] * gd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_z, pb_x, fd_46, fd_47, ff_33, gf_s_73, gf_s_74, \
                         gf_s_75, gd_46, gd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_z[k] * ff_33[k]
                  + f_2 * gf_s_73[k];

        t_74[k] = f_4 * fd_46[k]
                  + f_2 * gf_s_74[k]
                  + pb_x[k] * gd_46[k];

        t_75[k] = f_4 * fd_47[k]
                  + f_2 * gf_s_75[k]
                  + pb_x[k] * gd_47[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_x, ff_76, ff_77, ff_78, ff_79, gf_s_76, \
                         gf_s_77, gf_s_78, gf_s_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pa_x[k] * ff_76[k]
                  + f_2 * gf_s_76[k];

        t_77[k] = pa_x[k] * ff_77[k]
                  + f_2 * gf_s_77[k];

        t_78[k] = pa_x[k] * ff_78[k]
                  + f_2 * gf_s_78[k];

        t_79[k] = pa_x[k] * ff_79[k]
                  + f_2 * gf_s_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pa_y, pb_y, fd_30, ff_50, ff_52, gf_s_80, gf_s_81, \
                         gf_s_82, gd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_y[k] * ff_50[k]
                  + f_2 * gf_s_80[k];

        t_81[k] = f_4 * fd_30[k]
                  + f_2 * gf_s_81[k]
                  + pb_y[k] * gd_48[k];

        t_82[k] = pa_y[k] * ff_52[k]
                  + f_2 * gf_s_82[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pa_y, pb_x, fd_51, fd_52, ff_55, gf_s_83, gf_s_84, \
                         gf_s_85, gd_51, gd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_4 * fd_51[k]
                  + f_2 * gf_s_83[k]
                  + pb_x[k] * gd_51[k];

        t_84[k] = f_4 * fd_52[k]
                  + f_2 * gf_s_84[k]
                  + pb_x[k] * gd_52[k];

        t_85[k] = pa_y[k] * ff_55[k]
                  + f_2 * gf_s_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, ff_86, ff_87, ff_88, ff_89, gf_s_86, \
                         gf_s_87, gf_s_88, gf_s_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_x[k] * ff_86[k]
                  + f_2 * gf_s_86[k];

        t_87[k] = pa_x[k] * ff_87[k]
                  + f_2 * gf_s_87[k];

        t_88[k] = pa_x[k] * ff_88[k]
                  + f_2 * gf_s_88[k];

        t_89[k] = pa_x[k] * ff_89[k]
                  + f_2 * gf_s_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pa_x, pb_y, pb_z, fd_30, fd_54, ff_90, gf_s_90, \
                         gf_s_91, gf_s_92, gd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_5 * fd_54[k]
                  + pa_x[k] * ff_90[k]
                  + f_2 * gf_s_90[k];

        t_91[k] = f_2 * gf_s_91[k]
                  + pb_y[k] * gd_54[k];

        t_92[k] = f_5 * fd_30[k]
                  + f_2 * gf_s_92[k]
                  + pb_z[k] * gd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_x, pb_y, fd_57, fd_59, gf_s_93, gf_s_94, \
                         gf_s_95, gd_56, gd_57, gd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_4 * fd_57[k]
                  + f_2 * gf_s_93[k]
                  + pb_x[k] * gd_57[k];

        t_94[k] = f_2 * gf_s_94[k]
                  + pb_y[k] * gd_56[k];

        t_95[k] = f_4 * fd_59[k]
                  + f_2 * gf_s_95[k]
                  + pb_x[k] * gd_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pb_y, ff_96, ff_97, ff_99, gf_s_96, \
                         gf_s_97, gf_s_98, gf_s_99, gd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_x[k] * ff_96[k]
                  + f_2 * gf_s_96[k];

        t_97[k] = pa_x[k] * ff_97[k]
                  + f_2 * gf_s_97[k];

        t_98[k] = f_2 * gf_s_98[k]
                  + pb_y[k] * gd_59[k];

        t_99[k] = pa_x[k] * ff_99[k]
                  + f_2 * gf_s_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pb_x, pb_z, gp_s_30, gp_s_31, gf_s_100, \
                         gf_s_101, gf_s_102, gp_30, gp_31, gd_60, \
                         gd_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -f_1 * gp_s_30[k]
                   + f_2 * gf_s_100[k]
                   + f_3 * gp_30[k]
                   + pb_x[k] * gd_60[k];

        t_101[k] = -f_7 * gp_s_31[k]
                   + f_2 * gf_s_101[k]
                   + f_4 * gp_31[k]
                   + pb_x[k] * gd_61[k];

        t_102[k] = f_2 * gf_s_102[k]
                   + pb_z[k] * gd_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pb_x, pb_y, fd_39, gp_s_31, gf_s_103, \
                         gf_s_104, gf_s_105, gf_s_106, gp_31, gd_63, gd_64, \
                         gd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_2 * gf_s_103[k]
                   + pb_x[k] * gd_63[k];

        t_104[k] = f_2 * gf_s_104[k]
                   + pb_x[k] * gd_64[k];

        t_105[k] = f_2 * gf_s_105[k]
                   + pb_x[k] * gd_65[k];

        t_106[k] = f_0 * fd_39[k]
                   - f_1 * gp_s_31[k]
                   + f_2 * gf_s_106[k]
                   + f_3 * gp_31[k]
                   + pb_y[k] * gd_63[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_y, pb_z, fd_41, gp_s_32, gf_s_107, gf_s_108, \
                         gf_s_109, gp_32, gd_63, gd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_2 * gf_s_107[k]
                   + pb_z[k] * gd_63[k];

        t_108[k] = f_0 * fd_41[k]
                   + f_2 * gf_s_108[k]
                   + pb_y[k] * gd_65[k];

        t_109[k] = -f_1 * gp_s_32[k]
                   + f_2 * gf_s_109[k]
                   + f_3 * gp_32[k]
                   + pb_z[k] * gd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pb_x, ff_60, ff_61, gp_s_35, \
                         gf_s_110, gf_s_111, gf_s_112, gf_s_113, gp_35, gd_68, \
                         gd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * ff_60[k]
                   + f_2 * gf_s_110[k];

        t_111[k] = pa_z[k] * ff_61[k]
                   + f_2 * gf_s_111[k];

        t_112[k] = -f_7 * gp_s_35[k]
                   + f_2 * gf_s_112[k]
                   + f_4 * gp_35[k]
                   + pb_x[k] * gd_68[k];

        t_113[k] = f_2 * gf_s_113[k]
                   + pb_x[k] * gd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pb_x, pb_z, fd_39, ff_66, gf_s_114, \
                         gf_s_115, gf_s_116, gf_s_117, gd_69, gd_70, \
                         gd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_2 * gf_s_114[k]
                   + pb_x[k] * gd_70[k];

        t_115[k] = f_2 * gf_s_115[k]
                   + pb_x[k] * gd_71[k];

        t_116[k] = pa_z[k] * ff_66[k]
                   + f_2 * gf_s_116[k];

        t_117[k] = f_4 * fd_39[k]
                   + f_2 * gf_s_117[k]
                   + pb_z[k] * gd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, pa_y, pb_y, df_s_49, df_49, fd_47, ff_79, gf_s_118, \
                         gf_s_119, gd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_5 * fd_47[k]
                   + f_2 * gf_s_118[k]
                   + pb_y[k] * gd_71[k];

        t_119[k] = -f_6 * df_s_49[k]
                   + f_3 * df_49[k]
                   + pa_y[k] * ff_79[k]
                   + f_2 * gf_s_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pb_x, gp_s_36, gp_s_37, gp_s_38, gf_s_120, \
                         gf_s_121, gf_s_122, gp_36, gp_37, gp_38, gd_72, gd_73, \
                         gd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -f_1 * gp_s_36[k]
                   + f_2 * gf_s_120[k]
                   + f_3 * gp_36[k]
                   + pb_x[k] * gd_72[k];

        t_121[k] = -f_7 * gp_s_37[k]
                   + f_2 * gf_s_121[k]
                   + f_4 * gp_37[k]
                   + pb_x[k] * gd_73[k];

        t_122[k] = -f_7 * gp_s_38[k]
                   + f_2 * gf_s_122[k]
                   + f_4 * gp_38[k]
                   + pb_x[k] * gd_74[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_z, pb_x, df_s_36, df_36, ff_76, \
                         gf_s_123, gf_s_124, gf_s_125, gf_s_126, gd_75, gd_76, \
                         gd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_2 * gf_s_123[k]
                   + pb_x[k] * gd_75[k];

        t_124[k] = f_2 * gf_s_124[k]
                   + pb_x[k] * gd_76[k];

        t_125[k] = f_2 * gf_s_125[k]
                   + pb_x[k] * gd_77[k];

        t_126[k] = -f_8 * df_s_36[k]
                   + f_4 * df_36[k]
                   + pa_z[k] * ff_76[k]
                   + f_2 * gf_s_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pa_y, pb_y, pb_z, df_s_59, df_59, fd_45, fd_53, \
                         ff_89, gf_s_127, gf_s_128, gf_s_129, gd_75, \
                         gd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_3 * fd_45[k]
                   + f_2 * gf_s_127[k]
                   + pb_z[k] * gd_75[k];

        t_128[k] = f_3 * fd_53[k]
                   + f_2 * gf_s_128[k]
                   + pb_y[k] * gd_77[k];

        t_129[k] = -f_8 * df_s_59[k]
                   + f_4 * df_59[k]
                   + pa_y[k] * ff_89[k]
                   + f_2 * gf_s_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_y, pb_x, fd_54, ff_90, ff_91, ff_92, \
                         gf_s_130, gf_s_131, gf_s_132, gf_s_133, \
                         gd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pa_y[k] * ff_90[k]
                   + f_2 * gf_s_130[k];

        t_131[k] = f_4 * fd_54[k]
                   + pa_y[k] * ff_91[k]
                   + f_2 * gf_s_131[k];

        t_132[k] = pa_y[k] * ff_92[k]
                   + f_2 * gf_s_132[k];

        t_133[k] = f_2 * gf_s_133[k]
                   + pb_x[k] * gd_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pa_y, pb_x, fd_57, ff_96, gf_s_134, gf_s_135, \
                         gf_s_136, gd_82, gd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_2 * gf_s_134[k]
                   + pb_x[k] * gd_82[k];

        t_135[k] = f_2 * gf_s_135[k]
                   + pb_x[k] * gd_83[k];

        t_136[k] = f_5 * fd_57[k]
                   + pa_y[k] * ff_96[k]
                   + f_2 * gf_s_136[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pb_y, pb_z, fd_51, fd_59, ff_99, gf_s_137, \
                         gf_s_138, gf_s_139, gd_81, gd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_5 * fd_51[k]
                   + f_2 * gf_s_137[k]
                   + pb_z[k] * gd_81[k];

        t_138[k] = f_4 * fd_59[k]
                   + f_2 * gf_s_138[k]
                   + pb_y[k] * gd_83[k];

        t_139[k] = pa_y[k] * ff_99[k]
                   + f_2 * gf_s_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pb_x, pb_y, gp_s_42, gp_s_44, gf_s_140, \
                         gf_s_141, gf_s_142, gp_42, gp_44, gd_84, \
                         gd_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = -f_1 * gp_s_42[k]
                   + f_2 * gf_s_140[k]
                   + f_3 * gp_42[k]
                   + pb_x[k] * gd_84[k];

        t_141[k] = f_2 * gf_s_141[k]
                   + pb_y[k] * gd_84[k];

        t_142[k] = -f_7 * gp_s_44[k]
                   + f_2 * gf_s_142[k]
                   + f_4 * gp_44[k]
                   + pb_x[k] * gd_86[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pb_x, pb_y, gp_s_43, gf_s_143, gf_s_144, \
                         gf_s_145, gf_s_146, gp_43, gd_87, gd_88, \
                         gd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_2 * gf_s_143[k]
                   + pb_x[k] * gd_87[k];

        t_144[k] = f_2 * gf_s_144[k]
                   + pb_x[k] * gd_88[k];

        t_145[k] = f_2 * gf_s_145[k]
                   + pb_x[k] * gd_89[k];

        t_146[k] = -f_1 * gp_s_43[k]
                   + f_2 * gf_s_146[k]
                   + f_3 * gp_43[k]
                   + pb_y[k] * gd_87[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pb_y, pb_z, fd_59, gp_s_44, gf_s_147, gf_s_148, \
                         gf_s_149, gp_44, gd_88, gd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -f_7 * gp_s_44[k]
                   + f_2 * gf_s_147[k]
                   + f_4 * gp_44[k]
                   + pb_y[k] * gd_88[k];

        t_148[k] = f_2 * gf_s_148[k]
                   + pb_y[k] * gd_89[k];

        t_149[k] = f_0 * fd_59[k]
                   - f_1 * gp_s_44[k]
                   + f_2 * gf_s_149[k]
                   + f_3 * gp_44[k]
                   + pb_z[k] * gd_89[k];
    }
}

}  // namespace simdkin
