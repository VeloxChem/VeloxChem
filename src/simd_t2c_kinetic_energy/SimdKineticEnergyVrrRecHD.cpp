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


#include "SimdKineticEnergyVrrRecHD.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_hd_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 * beta / p;
    const auto f_6 = 1.5 / p;
    const auto f_7 = beta / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_15 = buffer.data(fd_s + 15);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_15 = buffer.data(fd + 15);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);

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
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_8 = buffer.data(hs_s + 8);
    const auto *hs_s_11 = buffer.data(hs_s + 11);
    const auto *hs_s_13 = buffer.data(hs_s + 13);
    const auto *hs_s_14 = buffer.data(hs_s + 14);
    const auto *hs_s_16 = buffer.data(hs_s + 16);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_4 = buffer.data(hd_s + 4);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_21 = buffer.data(hd_s + 21);
    const auto *hd_s_22 = buffer.data(hd_s + 22);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_24 = buffer.data(hd_s + 24);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_28 = buffer.data(hd_s + 28);
    const auto *hd_s_29 = buffer.data(hd_s + 29);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_38 = buffer.data(hd_s + 38);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_42 = buffer.data(hd_s + 42);
    const auto *hd_s_43 = buffer.data(hd_s + 43);
    const auto *hd_s_44 = buffer.data(hd_s + 44);
    const auto *hd_s_45 = buffer.data(hd_s + 45);
    const auto *hd_s_46 = buffer.data(hd_s + 46);
    const auto *hd_s_47 = buffer.data(hd_s + 47);
    const auto *hd_s_48 = buffer.data(hd_s + 48);
    const auto *hd_s_49 = buffer.data(hd_s + 49);
    const auto *hd_s_50 = buffer.data(hd_s + 50);
    const auto *hd_s_51 = buffer.data(hd_s + 51);
    const auto *hd_s_52 = buffer.data(hd_s + 52);
    const auto *hd_s_53 = buffer.data(hd_s + 53);
    const auto *hd_s_54 = buffer.data(hd_s + 54);
    const auto *hd_s_55 = buffer.data(hd_s + 55);
    const auto *hd_s_56 = buffer.data(hd_s + 56);
    const auto *hd_s_57 = buffer.data(hd_s + 57);
    const auto *hd_s_58 = buffer.data(hd_s + 58);
    const auto *hd_s_59 = buffer.data(hd_s + 59);
    const auto *hd_s_60 = buffer.data(hd_s + 60);
    const auto *hd_s_61 = buffer.data(hd_s + 61);
    const auto *hd_s_62 = buffer.data(hd_s + 62);
    const auto *hd_s_63 = buffer.data(hd_s + 63);
    const auto *hd_s_64 = buffer.data(hd_s + 64);
    const auto *hd_s_65 = buffer.data(hd_s + 65);
    const auto *hd_s_66 = buffer.data(hd_s + 66);
    const auto *hd_s_67 = buffer.data(hd_s + 67);
    const auto *hd_s_68 = buffer.data(hd_s + 68);
    const auto *hd_s_69 = buffer.data(hd_s + 69);
    const auto *hd_s_70 = buffer.data(hd_s + 70);
    const auto *hd_s_71 = buffer.data(hd_s + 71);
    const auto *hd_s_72 = buffer.data(hd_s + 72);
    const auto *hd_s_73 = buffer.data(hd_s + 73);
    const auto *hd_s_74 = buffer.data(hd_s + 74);
    const auto *hd_s_75 = buffer.data(hd_s + 75);
    const auto *hd_s_76 = buffer.data(hd_s + 76);
    const auto *hd_s_77 = buffer.data(hd_s + 77);
    const auto *hd_s_78 = buffer.data(hd_s + 78);
    const auto *hd_s_79 = buffer.data(hd_s + 79);
    const auto *hd_s_80 = buffer.data(hd_s + 80);
    const auto *hd_s_81 = buffer.data(hd_s + 81);
    const auto *hd_s_82 = buffer.data(hd_s + 82);
    const auto *hd_s_83 = buffer.data(hd_s + 83);
    const auto *hd_s_84 = buffer.data(hd_s + 84);
    const auto *hd_s_85 = buffer.data(hd_s + 85);
    const auto *hd_s_86 = buffer.data(hd_s + 86);
    const auto *hd_s_87 = buffer.data(hd_s + 87);
    const auto *hd_s_88 = buffer.data(hd_s + 88);
    const auto *hd_s_89 = buffer.data(hd_s + 89);
    const auto *hd_s_90 = buffer.data(hd_s + 90);
    const auto *hd_s_91 = buffer.data(hd_s + 91);
    const auto *hd_s_92 = buffer.data(hd_s + 92);
    const auto *hd_s_93 = buffer.data(hd_s + 93);
    const auto *hd_s_94 = buffer.data(hd_s + 94);
    const auto *hd_s_95 = buffer.data(hd_s + 95);
    const auto *hd_s_96 = buffer.data(hd_s + 96);
    const auto *hd_s_97 = buffer.data(hd_s + 97);
    const auto *hd_s_98 = buffer.data(hd_s + 98);
    const auto *hd_s_99 = buffer.data(hd_s + 99);
    const auto *hd_s_100 = buffer.data(hd_s + 100);
    const auto *hd_s_101 = buffer.data(hd_s + 101);
    const auto *hd_s_102 = buffer.data(hd_s + 102);
    const auto *hd_s_103 = buffer.data(hd_s + 103);
    const auto *hd_s_104 = buffer.data(hd_s + 104);
    const auto *hd_s_105 = buffer.data(hd_s + 105);
    const auto *hd_s_106 = buffer.data(hd_s + 106);
    const auto *hd_s_107 = buffer.data(hd_s + 107);
    const auto *hd_s_108 = buffer.data(hd_s + 108);
    const auto *hd_s_109 = buffer.data(hd_s + 109);
    const auto *hd_s_110 = buffer.data(hd_s + 110);
    const auto *hd_s_111 = buffer.data(hd_s + 111);
    const auto *hd_s_112 = buffer.data(hd_s + 112);
    const auto *hd_s_113 = buffer.data(hd_s + 113);
    const auto *hd_s_114 = buffer.data(hd_s + 114);
    const auto *hd_s_115 = buffer.data(hd_s + 115);
    const auto *hd_s_116 = buffer.data(hd_s + 116);
    const auto *hd_s_117 = buffer.data(hd_s + 117);
    const auto *hd_s_118 = buffer.data(hd_s + 118);
    const auto *hd_s_119 = buffer.data(hd_s + 119);
    const auto *hd_s_120 = buffer.data(hd_s + 120);
    const auto *hd_s_121 = buffer.data(hd_s + 121);
    const auto *hd_s_122 = buffer.data(hd_s + 122);
    const auto *hd_s_123 = buffer.data(hd_s + 123);
    const auto *hd_s_124 = buffer.data(hd_s + 124);
    const auto *hd_s_125 = buffer.data(hd_s + 125);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_16 = buffer.data(hs + 16);

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
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hd_s_3, hs_0, hp_0, hp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_2 * hd_s_1[k]
                 + pb_y[k] * hp_0[k];

        t_2[k] = f_2 * hd_s_2[k]
                 + pb_z[k] * hp_0[k];

        t_3[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_3[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_y, pb_y, pb_z, gd_0, hs_s_0, hd_s_4, hd_s_5, \
                         hd_s_6, hs_0, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hd_s_4[k]
                 + pb_y[k] * hp_2[k];

        t_5[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_5[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];

        t_6[k] = pa_y[k] * gd_0[k]
                 + f_2 * hd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_x, pb_z, fd_s_2, fd_2, gp_3, gd_5, hd_s_7, \
                         hd_s_8, hd_s_9, hp_3, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_4 * gp_3[k]
                 + f_2 * hd_s_7[k]
                 + pb_x[k] * hp_4[k];

        t_8[k] = f_2 * hd_s_8[k]
                 + pb_z[k] * hp_3[k];

        t_9[k] = -f_5 * fd_s_2[k]
                 + f_6 * fd_2[k]
                 + pa_x[k] * gd_5[k]
                 + f_2 * hd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_y, pb_z, gd_0, gd_2, hd_s_10, \
                         hd_s_11, hd_s_12, hd_s_13, hp_4, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * hd_s_10[k]
                  + pb_z[k] * hp_4[k];

        t_11[k] = pa_y[k] * gd_2[k]
                  + f_2 * hd_s_11[k];

        t_12[k] = pa_z[k] * gd_0[k]
                  + f_2 * hd_s_12[k];

        t_13[k] = f_2 * hd_s_13[k]
                  + pb_y[k] * hp_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_z, pb_x, pb_y, gp_4, gd_1, hd_s_14, hd_s_15, \
                         hd_s_16, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_4 * gp_4[k]
                  + f_2 * hd_s_14[k]
                  + pb_x[k] * hp_6[k];

        t_15[k] = pa_z[k] * gd_1[k]
                  + f_2 * hd_s_15[k];

        t_16[k] = f_2 * hd_s_16[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pa_y, fd_s_0, fd_s_4, fd_0, fd_4, gd_3, gd_8, \
                         hd_s_17, hd_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_5 * fd_s_4[k]
                  + f_6 * fd_4[k]
                  + pa_x[k] * gd_8[k]
                  + f_2 * hd_s_17[k];

        t_18[k] = -f_7 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_y[k] * gd_3[k]
                  + f_2 * hd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pb_x, pb_z, fd_s_5, fd_5, gp_5, gd_11, \
                         hd_s_19, hd_s_20, hd_s_21, hp_7, hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_6 * gp_5[k]
                  + f_2 * hd_s_19[k]
                  + pb_x[k] * hp_8[k];

        t_20[k] = f_2 * hd_s_20[k]
                  + pb_z[k] * hp_7[k];

        t_21[k] = -f_8 * fd_s_5[k]
                  + f_9 * fd_5[k]
                  + pa_x[k] * gd_11[k]
                  + f_2 * hd_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pb_z, gd_6, hs_s_3, hd_s_22, hd_s_23, \
                         hd_s_24, hs_3, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * hd_s_22[k]
                  + pb_z[k] * hp_8[k];

        t_23[k] = -f_1 * hs_s_3[k]
                  + f_2 * hd_s_23[k]
                  + f_3 * hs_3[k]
                  + pb_z[k] * hp_9[k];

        t_24[k] = pa_y[k] * gd_6[k]
                  + f_2 * hd_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, pa_z, pb_y, gp_4, gd_4, gd_5, gd_7, \
                         hd_s_25, hd_s_26, hd_s_27, hd_s_28, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * gd_4[k]
                  + f_2 * hd_s_25[k];

        t_26[k] = pa_y[k] * gd_7[k]
                  + f_2 * hd_s_26[k];

        t_27[k] = pa_z[k] * gd_5[k]
                  + f_2 * hd_s_27[k];

        t_28[k] = f_3 * gp_4[k]
                  + f_2 * hd_s_28[k]
                  + pb_y[k] * hp_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_y, fd_s_0, fd_0, gd_6, gd_8, \
                         hd_s_29, hd_s_30, hd_s_31, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gd_8[k]
                  + f_2 * hd_s_29[k];

        t_30[k] = -f_7 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_z[k] * gd_6[k]
                  + f_2 * hd_s_30[k];

        t_31[k] = f_2 * hd_s_31[k]
                  + pb_y[k] * hp_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, gp_7, hs_s_4, hd_s_32, hd_s_33, \
                         hd_s_34, hs_4, hp_12, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_6 * gp_7[k]
                  + f_2 * hd_s_32[k]
                  + pb_x[k] * hp_13[k];

        t_33[k] = -f_1 * hs_s_4[k]
                  + f_2 * hd_s_33[k]
                  + f_3 * hs_4[k]
                  + pb_y[k] * hp_12[k];

        t_34[k] = f_2 * hd_s_34[k]
                  + pb_y[k] * hp_13[k];
    }

#pragma omp simd aligned(t_35, t_36, pa_x, pa_y, fd_s_1, fd_s_6, fd_1, fd_6, gd_9, gd_14, \
                         hd_s_35, hd_s_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_8 * fd_s_6[k]
                  + f_9 * fd_6[k]
                  + pa_x[k] * gd_14[k]
                  + f_2 * hd_s_35[k];

        t_36[k] = -f_8 * fd_s_1[k]
                  + f_9 * fd_1[k]
                  + pa_y[k] * gd_9[k]
                  + f_2 * hd_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_x, pb_x, pb_z, fd_s_8, fd_8, gp_8, gd_17, \
                         hd_s_37, hd_s_38, hd_s_39, hp_14, hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * gp_8[k]
                  + f_2 * hd_s_37[k]
                  + pb_x[k] * hp_15[k];

        t_38[k] = f_2 * hd_s_38[k]
                  + pb_z[k] * hp_14[k];

        t_39[k] = -f_7 * fd_s_8[k]
                  + f_3 * fd_8[k]
                  + pa_x[k] * gd_17[k]
                  + f_2 * hd_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_z, pb_z, gd_9, gd_10, hs_s_5, hd_s_40, \
                         hd_s_41, hd_s_42, hd_s_43, hs_5, hp_15, \
                         hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * hd_s_40[k]
                  + pb_z[k] * hp_15[k];

        t_41[k] = -f_1 * hs_s_5[k]
                  + f_2 * hd_s_41[k]
                  + f_3 * hs_5[k]
                  + pb_z[k] * hp_16[k];

        t_42[k] = pa_z[k] * gd_9[k]
                  + f_2 * hd_s_42[k];

        t_43[k] = pa_z[k] * gd_10[k]
                  + f_2 * hd_s_43[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_z, pb_x, pb_y, gp_6, gp_9, gd_11, hd_s_44, \
                         hd_s_45, hd_s_46, hp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * gp_9[k]
                  + f_2 * hd_s_44[k]
                  + pb_x[k] * hp_17[k];

        t_45[k] = pa_z[k] * gd_11[k]
                  + f_2 * hd_s_45[k];

        t_46[k] = f_9 * gp_6[k]
                  + f_2 * hd_s_46[k]
                  + pb_y[k] * hp_17[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pa_y, pb_x, fd_s_10, fd_10, gp_10, gd_12, \
                         gd_18, hd_s_47, hd_s_48, hd_s_49, hp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -f_7 * fd_s_10[k]
                  + f_3 * fd_10[k]
                  + pa_x[k] * gd_18[k]
                  + f_2 * hd_s_47[k];

        t_48[k] = pa_y[k] * gd_12[k]
                  + f_2 * hd_s_48[k];

        t_49[k] = f_9 * gp_10[k]
                  + f_2 * hd_s_49[k]
                  + pb_x[k] * hp_18[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pa_y, pb_y, fd_s_11, fd_11, gp_7, gd_13, \
                         gd_19, hd_s_50, hd_s_51, hd_s_52, hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * gd_13[k]
                  + f_2 * hd_s_50[k];

        t_51[k] = -f_7 * fd_s_11[k]
                  + f_3 * fd_11[k]
                  + pa_x[k] * gd_19[k]
                  + f_2 * hd_s_51[k];

        t_52[k] = f_3 * gp_7[k]
                  + f_2 * hd_s_52[k]
                  + pb_y[k] * hp_19[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pa_z, pb_y, fd_s_3, fd_3, gd_12, gd_14, \
                         hd_s_53, hd_s_54, hd_s_55, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_y[k] * gd_14[k]
                  + f_2 * hd_s_53[k];

        t_54[k] = -f_8 * fd_s_3[k]
                  + f_9 * fd_3[k]
                  + pa_z[k] * gd_12[k]
                  + f_2 * hd_s_54[k];

        t_55[k] = f_2 * hd_s_55[k]
                  + pb_y[k] * hp_20[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pb_x, pb_y, gp_11, hs_s_8, hd_s_56, hd_s_57, \
                         hd_s_58, hs_8, hp_21, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_9 * gp_11[k]
                  + f_2 * hd_s_56[k]
                  + pb_x[k] * hp_22[k];

        t_57[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_57[k]
                  + f_3 * hs_8[k]
                  + pb_y[k] * hp_21[k];

        t_58[k] = f_2 * hd_s_58[k]
                  + pb_y[k] * hp_22[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_x, pb_x, fd_s_15, fd_15, gp_12, gp_13, gd_22, \
                         gd_23, hd_s_59, hd_s_60, hd_s_61, hp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_7 * fd_s_15[k]
                  + f_3 * fd_15[k]
                  + pa_x[k] * gd_22[k]
                  + f_2 * hd_s_59[k];

        t_60[k] = f_9 * gp_12[k]
                  + pa_x[k] * gd_23[k]
                  + f_2 * hd_s_60[k];

        t_61[k] = f_3 * gp_13[k]
                  + f_2 * hd_s_61[k]
                  + pb_x[k] * hp_24[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pb_z, gd_24, gd_25, hd_s_62, hd_s_63, \
                         hd_s_64, hd_s_65, hp_23, hp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_2 * hd_s_62[k]
                  + pb_z[k] * hp_23[k];

        t_63[k] = pa_x[k] * gd_24[k]
                  + f_2 * hd_s_63[k];

        t_64[k] = f_2 * hd_s_64[k]
                  + pb_z[k] * hp_24[k];

        t_65[k] = pa_x[k] * gd_25[k]
                  + f_2 * hd_s_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pa_z, pb_x, gp_15, gd_15, gd_16, gd_26, \
                         hd_s_66, hd_s_67, hd_s_68, hd_s_69, hp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * gd_15[k]
                  + f_2 * hd_s_66[k];

        t_67[k] = pa_z[k] * gd_16[k]
                  + f_2 * hd_s_67[k];

        t_68[k] = f_3 * gp_15[k]
                  + f_2 * hd_s_68[k]
                  + pb_x[k] * hp_25[k];

        t_69[k] = pa_x[k] * gd_26[k]
                  + f_2 * hd_s_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pb_x, gp_16, gp_17, gd_27, gd_28, \
                         gd_29, hd_s_70, hd_s_71, hd_s_72, hd_s_73, \
                         hp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_x[k] * gd_27[k]
                  + f_2 * hd_s_70[k];

        t_71[k] = pa_x[k] * gd_28[k]
                  + f_2 * hd_s_71[k];

        t_72[k] = f_9 * gp_16[k]
                  + pa_x[k] * gd_29[k]
                  + f_2 * hd_s_72[k];

        t_73[k] = f_3 * gp_17[k]
                  + f_2 * hd_s_73[k]
                  + pb_x[k] * hp_26[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pb_x, gp_18, gd_30, gd_31, gd_32, \
                         hd_s_74, hd_s_75, hd_s_76, hd_s_77, hp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * gp_18[k]
                  + f_2 * hd_s_74[k]
                  + pb_x[k] * hp_27[k];

        t_75[k] = pa_x[k] * gd_30[k]
                  + f_2 * hd_s_75[k];

        t_76[k] = pa_x[k] * gd_31[k]
                  + f_2 * hd_s_76[k];

        t_77[k] = pa_x[k] * gd_32[k]
                  + f_2 * hd_s_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pa_y, pb_x, gp_19, gd_20, gd_21, gd_33, \
                         hd_s_78, hd_s_79, hd_s_80, hd_s_81, hp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * gd_20[k]
                  + f_2 * hd_s_78[k];

        t_79[k] = f_3 * gp_19[k]
                  + f_2 * hd_s_79[k]
                  + pb_x[k] * hp_28[k];

        t_80[k] = pa_y[k] * gd_21[k]
                  + f_2 * hd_s_80[k];

        t_81[k] = pa_x[k] * gd_33[k]
                  + f_2 * hd_s_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pb_y, gp_21, gd_34, gd_35, gd_36, \
                         hd_s_82, hd_s_83, hd_s_84, hd_s_85, hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_x[k] * gd_34[k]
                  + f_2 * hd_s_82[k];

        t_83[k] = pa_x[k] * gd_35[k]
                  + f_2 * hd_s_83[k];

        t_84[k] = f_9 * gp_21[k]
                  + pa_x[k] * gd_36[k]
                  + f_2 * hd_s_84[k];

        t_85[k] = f_2 * hd_s_85[k]
                  + pb_y[k] * hp_29[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pb_x, pb_y, gp_23, gd_37, gd_38, \
                         hd_s_86, hd_s_87, hd_s_88, hd_s_89, hp_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_3 * gp_23[k]
                  + f_2 * hd_s_86[k]
                  + pb_x[k] * hp_30[k];

        t_87[k] = pa_x[k] * gd_37[k]
                  + f_2 * hd_s_87[k];

        t_88[k] = f_2 * hd_s_88[k]
                  + pb_y[k] * hp_30[k];

        t_89[k] = pa_x[k] * gd_38[k]
                  + f_2 * hd_s_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_x, pb_y, gp_13, hs_s_11, hd_s_90, hd_s_91, \
                         hd_s_92, hd_s_93, hs_11, hp_31, hp_32, hp_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_90[k]
                  + f_3 * hs_11[k]
                  + pb_x[k] * hp_31[k];

        t_91[k] = f_2 * hd_s_91[k]
                  + pb_x[k] * hp_32[k];

        t_92[k] = f_2 * hd_s_92[k]
                  + pb_x[k] * hp_33[k];

        t_93[k] = f_0 * gp_13[k]
                  - f_1 * hs_s_11[k]
                  + f_2 * hd_s_93[k]
                  + f_3 * hs_11[k]
                  + pb_y[k] * hp_32[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_z, pb_z, gd_23, hs_s_11, hd_s_94, hd_s_95, \
                         hd_s_96, hs_11, hp_32, hp_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_2 * hd_s_94[k]
                  + pb_z[k] * hp_32[k];

        t_95[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_95[k]
                  + f_3 * hs_11[k]
                  + pb_z[k] * hp_33[k];

        t_96[k] = pa_z[k] * gd_23[k]
                  + f_2 * hd_s_96[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_z, pb_x, pb_y, gp_15, gd_24, hd_s_97, \
                         hd_s_98, hd_s_99, hd_s_100, hp_34, hp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_2 * hd_s_97[k]
                  + pb_x[k] * hp_34[k];

        t_98[k] = f_2 * hd_s_98[k]
                  + pb_x[k] * hp_35[k];

        t_99[k] = pa_z[k] * gd_24[k]
                  + f_2 * hd_s_99[k];

        t_100[k] = f_4 * gp_15[k]
                   + f_2 * hd_s_100[k]
                   + pb_y[k] * hp_35[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_y, pb_x, fd_s_10, fd_10, gd_28, hs_s_13, \
                         hd_s_101, hd_s_102, hd_s_103, hs_13, hp_36, \
                         hp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -f_5 * fd_s_10[k]
                   + f_6 * fd_10[k]
                   + pa_y[k] * gd_28[k]
                   + f_2 * hd_s_101[k];

        t_102[k] = -f_1 * hs_s_13[k]
                   + f_2 * hd_s_102[k]
                   + f_3 * hs_13[k]
                   + pb_x[k] * hp_36[k];

        t_103[k] = f_2 * hd_s_103[k]
                   + pb_x[k] * hp_37[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pb_x, pb_y, fd_s_8, fd_8, gp_18, gd_26, \
                         hd_s_104, hd_s_105, hd_s_106, hp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_2 * hd_s_104[k]
                   + pb_x[k] * hp_38[k];

        t_105[k] = -f_7 * fd_s_8[k]
                   + f_3 * fd_8[k]
                   + pa_z[k] * gd_26[k]
                   + f_2 * hd_s_105[k];

        t_106[k] = f_6 * gp_18[k]
                   + f_2 * hd_s_106[k]
                   + pb_y[k] * hp_38[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pa_y, pb_x, fd_s_12, fd_12, gd_32, hs_s_14, \
                         hd_s_107, hd_s_108, hd_s_109, hs_14, hp_39, \
                         hp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -f_8 * fd_s_12[k]
                   + f_9 * fd_12[k]
                   + pa_y[k] * gd_32[k]
                   + f_2 * hd_s_107[k];

        t_108[k] = -f_1 * hs_s_14[k]
                   + f_2 * hd_s_108[k]
                   + f_3 * hs_14[k]
                   + pb_x[k] * hp_39[k];

        t_109[k] = f_2 * hd_s_109[k]
                   + pb_x[k] * hp_40[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pb_x, pb_y, fd_s_9, fd_9, gp_20, gd_30, \
                         hd_s_110, hd_s_111, hd_s_112, hp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_2 * hd_s_110[k]
                   + pb_x[k] * hp_41[k];

        t_111[k] = -f_8 * fd_s_9[k]
                   + f_9 * fd_9[k]
                   + pa_z[k] * gd_30[k]
                   + f_2 * hd_s_111[k];

        t_112[k] = f_9 * gp_20[k]
                   + f_2 * hd_s_112[k]
                   + pb_y[k] * hp_41[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_y, pb_x, fd_s_15, fd_15, gd_35, gd_36, \
                         hd_s_113, hd_s_114, hd_s_115, hd_s_116, hp_42, \
                         hp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = -f_7 * fd_s_15[k]
                   + f_3 * fd_15[k]
                   + pa_y[k] * gd_35[k]
                   + f_2 * hd_s_113[k];

        t_114[k] = pa_y[k] * gd_36[k]
                   + f_2 * hd_s_114[k];

        t_115[k] = f_2 * hd_s_115[k]
                   + pb_x[k] * hp_42[k];

        t_116[k] = f_2 * hd_s_116[k]
                   + pb_x[k] * hp_43[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_y, pb_y, gp_22, gp_23, gd_37, gd_38, \
                         hd_s_117, hd_s_118, hd_s_119, hp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_9 * gp_22[k]
                   + pa_y[k] * gd_37[k]
                   + f_2 * hd_s_117[k];

        t_118[k] = f_3 * gp_23[k]
                   + f_2 * hd_s_118[k]
                   + pb_y[k] * hp_43[k];

        t_119[k] = pa_y[k] * gd_38[k]
                   + f_2 * hd_s_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_x, pb_y, hs_s_16, hd_s_120, hd_s_121, \
                         hd_s_122, hd_s_123, hs_16, hp_44, hp_45, \
                         hp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -f_1 * hs_s_16[k]
                   + f_2 * hd_s_120[k]
                   + f_3 * hs_16[k]
                   + pb_x[k] * hp_44[k];

        t_121[k] = f_2 * hd_s_121[k]
                   + pb_x[k] * hp_45[k];

        t_122[k] = f_2 * hd_s_122[k]
                   + pb_x[k] * hp_46[k];

        t_123[k] = -f_1 * hs_s_16[k]
                   + f_2 * hd_s_123[k]
                   + f_3 * hs_16[k]
                   + pb_y[k] * hp_45[k];
    }

#pragma omp simd aligned(t_124, t_125, pb_y, pb_z, gp_23, hs_s_16, hd_s_124, hd_s_125, hs_16, \
                         hp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * hd_s_124[k]
                   + pb_y[k] * hp_46[k];

        t_125[k] = f_0 * gp_23[k]
                   - f_1 * hs_s_16[k]
                   + f_2 * hd_s_125[k]
                   + f_3 * hs_16[k]
                   + pb_z[k] * hp_46[k];
    }
}

auto
compute_prim_hd_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 * beta / p;
    const auto f_6 = 1.5 / p;
    const auto f_7 = beta / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_19 = buffer.data(fd_s + 19);
    const auto *fd_s_22 = buffer.data(fd_s + 22);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_22 = buffer.data(fd + 22);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_47 = buffer.data(gd + 47);

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_7 = buffer.data(hs_s + 7);
    const auto *hs_s_8 = buffer.data(hs_s + 8);
    const auto *hs_s_10 = buffer.data(hs_s + 10);
    const auto *hs_s_11 = buffer.data(hs_s + 11);
    const auto *hs_s_13 = buffer.data(hs_s + 13);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_4 = buffer.data(hd_s + 4);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_21 = buffer.data(hd_s + 21);
    const auto *hd_s_22 = buffer.data(hd_s + 22);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_24 = buffer.data(hd_s + 24);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_28 = buffer.data(hd_s + 28);
    const auto *hd_s_29 = buffer.data(hd_s + 29);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_38 = buffer.data(hd_s + 38);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_42 = buffer.data(hd_s + 42);
    const auto *hd_s_43 = buffer.data(hd_s + 43);
    const auto *hd_s_44 = buffer.data(hd_s + 44);
    const auto *hd_s_45 = buffer.data(hd_s + 45);
    const auto *hd_s_46 = buffer.data(hd_s + 46);
    const auto *hd_s_47 = buffer.data(hd_s + 47);
    const auto *hd_s_48 = buffer.data(hd_s + 48);
    const auto *hd_s_49 = buffer.data(hd_s + 49);
    const auto *hd_s_50 = buffer.data(hd_s + 50);
    const auto *hd_s_51 = buffer.data(hd_s + 51);
    const auto *hd_s_52 = buffer.data(hd_s + 52);
    const auto *hd_s_53 = buffer.data(hd_s + 53);
    const auto *hd_s_54 = buffer.data(hd_s + 54);
    const auto *hd_s_55 = buffer.data(hd_s + 55);
    const auto *hd_s_56 = buffer.data(hd_s + 56);
    const auto *hd_s_57 = buffer.data(hd_s + 57);
    const auto *hd_s_58 = buffer.data(hd_s + 58);
    const auto *hd_s_59 = buffer.data(hd_s + 59);
    const auto *hd_s_60 = buffer.data(hd_s + 60);
    const auto *hd_s_61 = buffer.data(hd_s + 61);
    const auto *hd_s_62 = buffer.data(hd_s + 62);
    const auto *hd_s_63 = buffer.data(hd_s + 63);
    const auto *hd_s_64 = buffer.data(hd_s + 64);
    const auto *hd_s_65 = buffer.data(hd_s + 65);
    const auto *hd_s_66 = buffer.data(hd_s + 66);
    const auto *hd_s_67 = buffer.data(hd_s + 67);
    const auto *hd_s_68 = buffer.data(hd_s + 68);
    const auto *hd_s_69 = buffer.data(hd_s + 69);
    const auto *hd_s_70 = buffer.data(hd_s + 70);
    const auto *hd_s_71 = buffer.data(hd_s + 71);
    const auto *hd_s_72 = buffer.data(hd_s + 72);
    const auto *hd_s_73 = buffer.data(hd_s + 73);
    const auto *hd_s_74 = buffer.data(hd_s + 74);
    const auto *hd_s_75 = buffer.data(hd_s + 75);
    const auto *hd_s_76 = buffer.data(hd_s + 76);
    const auto *hd_s_77 = buffer.data(hd_s + 77);
    const auto *hd_s_78 = buffer.data(hd_s + 78);
    const auto *hd_s_79 = buffer.data(hd_s + 79);
    const auto *hd_s_81 = buffer.data(hd_s + 81);
    const auto *hd_s_82 = buffer.data(hd_s + 82);
    const auto *hd_s_83 = buffer.data(hd_s + 83);
    const auto *hd_s_84 = buffer.data(hd_s + 84);
    const auto *hd_s_85 = buffer.data(hd_s + 85);
    const auto *hd_s_86 = buffer.data(hd_s + 86);
    const auto *hd_s_87 = buffer.data(hd_s + 87);
    const auto *hd_s_88 = buffer.data(hd_s + 88);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hs_0, hp_0, hp_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = f_2 * hd_s_1[k]
                 + pb_z[k] * hp_0[k];

        t_2[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_2[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pb_x, pb_z, gp_3, gd_0, hs_s_0, hd_s_3, hd_s_4, \
                         hd_s_5, hs_0, hp_2, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_3[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];

        t_4[k] = pa_y[k] * gd_0[k]
                 + f_2 * hd_s_4[k];

        t_5[k] = f_4 * gp_3[k]
                 + f_2 * hd_s_5[k]
                 + pb_x[k] * hp_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pa_z, fd_s_4, fd_4, gd_0, gd_2, gd_4, \
                         hd_s_6, hd_s_7, hd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_5 * fd_s_4[k]
                 + f_6 * fd_4[k]
                 + pa_x[k] * gd_4[k]
                 + f_2 * hd_s_6[k];

        t_7[k] = pa_y[k] * gd_2[k]
                 + f_2 * hd_s_7[k];

        t_8[k] = pa_z[k] * gd_0[k]
                 + f_2 * hd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_x, pb_y, fd_s_6, fd_6, gp_4, gd_7, hd_s_9, \
                         hd_s_10, hd_s_11, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * gp_4[k]
                 + f_2 * hd_s_9[k]
                 + pb_x[k] * hp_4[k];

        t_10[k] = f_2 * hd_s_10[k]
                  + pb_y[k] * hp_4[k];

        t_11[k] = -f_5 * fd_s_6[k]
                  + f_6 * fd_6[k]
                  + pa_x[k] * gd_7[k]
                  + f_2 * hd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, pa_y, pb_x, fd_s_0, fd_0, gp_5, gd_3, hd_s_12, hd_s_13, \
                         hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_7 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_y[k] * gd_3[k]
                  + f_2 * hd_s_12[k];

        t_13[k] = f_6 * gp_5[k]
                  + f_2 * hd_s_13[k]
                  + pb_x[k] * hp_5[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_z, pb_z, fd_s_8, fd_8, gd_4, gd_9, hs_s_3, \
                         hd_s_14, hd_s_15, hd_s_16, hs_3, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -f_8 * fd_s_8[k]
                  + f_9 * fd_8[k]
                  + pa_x[k] * gd_9[k]
                  + f_2 * hd_s_14[k];

        t_15[k] = -f_1 * hs_s_3[k]
                  + f_2 * hd_s_15[k]
                  + f_3 * hs_3[k]
                  + pb_z[k] * hp_6[k];

        t_16[k] = pa_z[k] * gd_4[k]
                  + f_2 * hd_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_y, pa_z, pb_y, fd_s_0, fd_0, gp_4, gd_6, gd_7, \
                         hd_s_17, hd_s_18, hd_s_19, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * gp_4[k]
                  + f_2 * hd_s_17[k]
                  + pb_y[k] * hp_7[k];

        t_18[k] = pa_y[k] * gd_7[k]
                  + f_2 * hd_s_18[k];

        t_19[k] = -f_7 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_z[k] * gd_6[k]
                  + f_2 * hd_s_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_y, gp_9, hs_s_4, hd_s_20, hd_s_21, \
                         hd_s_22, hs_4, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_6 * gp_9[k]
                  + f_2 * hd_s_20[k]
                  + pb_x[k] * hp_9[k];

        t_21[k] = -f_1 * hs_s_4[k]
                  + f_2 * hd_s_21[k]
                  + f_3 * hs_4[k]
                  + pb_y[k] * hp_8[k];

        t_22[k] = f_2 * hd_s_22[k]
                  + pb_y[k] * hp_9[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_x, pa_y, fd_s_3, fd_s_10, fd_3, fd_10, gd_8, gd_16, \
                         hd_s_23, hd_s_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_8 * fd_s_10[k]
                  + f_9 * fd_10[k]
                  + pa_x[k] * gd_16[k]
                  + f_2 * hd_s_23[k];

        t_24[k] = -f_8 * fd_s_3[k]
                  + f_9 * fd_3[k]
                  + pa_y[k] * gd_8[k]
                  + f_2 * hd_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pa_x, pb_x, fd_s_12, fd_12, gp_10, gd_18, hd_s_25, \
                         hd_s_26, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_9 * gp_10[k]
                  + f_2 * hd_s_25[k]
                  + pb_x[k] * hp_10[k];

        t_26[k] = -f_7 * fd_s_12[k]
                  + f_3 * fd_12[k]
                  + pa_x[k] * gd_18[k]
                  + f_2 * hd_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_z, pb_z, gd_8, gd_9, hs_s_5, hd_s_27, hd_s_28, \
                         hd_s_29, hs_5, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * hs_s_5[k]
                  + f_2 * hd_s_27[k]
                  + f_3 * hs_5[k]
                  + pb_z[k] * hp_11[k];

        t_28[k] = pa_z[k] * gd_8[k]
                  + f_2 * hd_s_28[k];

        t_29[k] = pa_z[k] * gd_9[k]
                  + f_2 * hd_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_x, pa_y, pb_y, fd_s_16, fd_16, gp_7, gd_14, \
                         gd_22, hd_s_30, hd_s_31, hd_s_32, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * gp_7[k]
                  + f_2 * hd_s_30[k]
                  + pb_y[k] * hp_12[k];

        t_31[k] = -f_7 * fd_s_16[k]
                  + f_3 * fd_16[k]
                  + pa_x[k] * gd_22[k]
                  + f_2 * hd_s_31[k];

        t_32[k] = pa_y[k] * gd_14[k]
                  + f_2 * hd_s_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_x, pa_y, pb_y, fd_s_17, fd_17, gp_9, gd_16, \
                         gd_24, hd_s_33, hd_s_34, hd_s_35, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -f_7 * fd_s_17[k]
                  + f_3 * fd_17[k]
                  + pa_x[k] * gd_24[k]
                  + f_2 * hd_s_33[k];

        t_34[k] = f_3 * gp_9[k]
                  + f_2 * hd_s_34[k]
                  + pb_y[k] * hp_13[k];

        t_35[k] = pa_y[k] * gd_16[k]
                  + f_2 * hd_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, pa_z, pb_x, fd_s_5, fd_5, gp_11, gd_14, hd_s_36, hd_s_37, \
                         hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_8 * fd_s_5[k]
                  + f_9 * fd_5[k]
                  + pa_z[k] * gd_14[k]
                  + f_2 * hd_s_36[k];

        t_37[k] = f_9 * gp_11[k]
                  + f_2 * hd_s_37[k]
                  + pb_x[k] * hp_15[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pb_y, fd_s_22, fd_22, gd_28, hs_s_7, hd_s_38, \
                         hd_s_39, hd_s_40, hs_7, hp_14, hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -f_1 * hs_s_7[k]
                  + f_2 * hd_s_38[k]
                  + f_3 * hs_7[k]
                  + pb_y[k] * hp_14[k];

        t_39[k] = f_2 * hd_s_39[k]
                  + pb_y[k] * hp_15[k];

        t_40[k] = -f_7 * fd_s_22[k]
                  + f_3 * fd_22[k]
                  + pa_x[k] * gd_28[k]
                  + f_2 * hd_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_x, gp_12, gp_13, gd_29, gd_31, \
                         gd_32, hd_s_41, hd_s_42, hd_s_43, hd_s_44, \
                         hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_9 * gp_12[k]
                  + pa_x[k] * gd_29[k]
                  + f_2 * hd_s_41[k];

        t_42[k] = f_3 * gp_13[k]
                  + f_2 * hd_s_42[k]
                  + pb_x[k] * hp_16[k];

        t_43[k] = pa_x[k] * gd_31[k]
                  + f_2 * hd_s_43[k];

        t_44[k] = pa_x[k] * gd_32[k]
                  + f_2 * hd_s_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pa_z, gp_16, gd_17, gd_34, gd_35, \
                         gd_36, hd_s_45, hd_s_46, hd_s_47, hd_s_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * gd_17[k]
                  + f_2 * hd_s_45[k];

        t_46[k] = pa_x[k] * gd_34[k]
                  + f_2 * hd_s_46[k];

        t_47[k] = pa_x[k] * gd_35[k]
                  + f_2 * hd_s_47[k];

        t_48[k] = f_9 * gp_16[k]
                  + pa_x[k] * gd_36[k]
                  + f_2 * hd_s_48[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pa_y, gd_26, gd_37, gd_38, gd_39, \
                         hd_s_49, hd_s_50, hd_s_51, hd_s_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_x[k] * gd_37[k]
                  + f_2 * hd_s_49[k];

        t_50[k] = pa_x[k] * gd_38[k]
                  + f_2 * hd_s_50[k];

        t_51[k] = pa_x[k] * gd_39[k]
                  + f_2 * hd_s_51[k];

        t_52[k] = pa_y[k] * gd_26[k]
                  + f_2 * hd_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_x, pb_x, gp_19, gp_21, gd_40, gd_41, \
                         gd_43, hd_s_53, hd_s_54, hd_s_55, hd_s_56, \
                         hp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_x[k] * gd_40[k]
                  + f_2 * hd_s_53[k];

        t_54[k] = pa_x[k] * gd_41[k]
                  + f_2 * hd_s_54[k];

        t_55[k] = f_9 * gp_19[k]
                  + pa_x[k] * gd_43[k]
                  + f_2 * hd_s_55[k];

        t_56[k] = f_3 * gp_21[k]
                  + f_2 * hd_s_56[k]
                  + pb_x[k] * hp_17[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_x, pb_x, gd_45, gd_47, hs_s_8, hd_s_57, \
                         hd_s_58, hd_s_59, hd_s_60, hs_8, hp_18, \
                         hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_x[k] * gd_45[k]
                  + f_2 * hd_s_57[k];

        t_58[k] = pa_x[k] * gd_47[k]
                  + f_2 * hd_s_58[k];

        t_59[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_59[k]
                  + f_3 * hs_8[k]
                  + pb_x[k] * hp_18[k];

        t_60[k] = f_2 * hd_s_60[k]
                  + pb_x[k] * hp_19[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_y, pb_z, gp_13, hs_s_8, hd_s_61, hd_s_62, \
                         hd_s_63, hs_8, hp_19, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * gp_13[k]
                  - f_1 * hs_s_8[k]
                  + f_2 * hd_s_61[k]
                  + f_3 * hs_8[k]
                  + pb_y[k] * hp_19[k];

        t_62[k] = f_2 * hd_s_62[k]
                  + pb_z[k] * hp_19[k];

        t_63[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_63[k]
                  + f_3 * hs_8[k]
                  + pb_z[k] * hp_20[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_z, pb_x, pb_y, gp_15, gd_31, hd_s_64, hd_s_65, \
                         hd_s_66, hp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * hd_s_64[k]
                  + pb_x[k] * hp_21[k];

        t_65[k] = pa_z[k] * gd_31[k]
                  + f_2 * hd_s_65[k];

        t_66[k] = f_4 * gp_15[k]
                  + f_2 * hd_s_66[k]
                  + pb_y[k] * hp_21[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_y, pb_x, fd_s_16, fd_16, gd_35, hs_s_10, \
                         hd_s_67, hd_s_68, hd_s_69, hs_10, hp_22, \
                         hp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -f_5 * fd_s_16[k]
                  + f_6 * fd_16[k]
                  + pa_y[k] * gd_35[k]
                  + f_2 * hd_s_67[k];

        t_68[k] = -f_1 * hs_s_10[k]
                  + f_2 * hd_s_68[k]
                  + f_3 * hs_10[k]
                  + pb_x[k] * hp_22[k];

        t_69[k] = f_2 * hd_s_69[k]
                  + pb_x[k] * hp_23[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pb_x, pb_y, fd_s_12, fd_12, gp_17, gd_33, \
                         hd_s_70, hd_s_71, hd_s_72, hp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_2 * hd_s_70[k]
                  + pb_x[k] * hp_24[k];

        t_71[k] = -f_7 * fd_s_12[k]
                  + f_3 * fd_12[k]
                  + pa_z[k] * gd_33[k]
                  + f_2 * hd_s_71[k];

        t_72[k] = f_6 * gp_17[k]
                  + f_2 * hd_s_72[k]
                  + pb_y[k] * hp_24[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_y, pb_x, fd_s_19, fd_19, gd_39, hs_s_11, \
                         hd_s_73, hd_s_74, hd_s_75, hs_11, hp_25, \
                         hp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = -f_8 * fd_s_19[k]
                  + f_9 * fd_19[k]
                  + pa_y[k] * gd_39[k]
                  + f_2 * hd_s_73[k];

        t_74[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_74[k]
                  + f_3 * hs_11[k]
                  + pb_x[k] * hp_25[k];

        t_75[k] = f_2 * hd_s_75[k]
                  + pb_x[k] * hp_26[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_z, pb_x, pb_y, fd_s_14, fd_14, gp_18, gd_37, \
                         hd_s_76, hd_s_77, hd_s_78, hp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * hd_s_76[k]
                  + pb_x[k] * hp_27[k];

        t_77[k] = -f_8 * fd_s_14[k]
                  + f_9 * fd_14[k]
                  + pa_z[k] * gd_37[k]
                  + f_2 * hd_s_77[k];

        t_78[k] = f_9 * gp_18[k]
                  + f_2 * hd_s_78[k]
                  + pb_y[k] * hp_27[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pa_y, pb_y, fd_s_22, fd_22, gp_20, gp_21, gd_42, \
                         gd_45, hd_s_79, hd_s_81, hd_s_82, hp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -f_7 * fd_s_22[k]
                  + f_3 * fd_22[k]
                  + pa_y[k] * gd_42[k]
                  + f_2 * hd_s_79[k];

        t_80[k] = f_9 * gp_20[k]
                  + pa_y[k] * gd_45[k]
                  + f_2 * hd_s_81[k];

        t_81[k] = f_3 * gp_21[k]
                  + f_2 * hd_s_82[k]
                  + pb_y[k] * hp_28[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pb_x, gd_47, hs_s_13, hd_s_83, hd_s_84, \
                         hd_s_85, hs_12, hp_29, hp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_y[k] * gd_47[k]
                  + f_2 * hd_s_83[k];

        t_83[k] = -f_1 * hs_s_13[k]
                  + f_2 * hd_s_84[k]
                  + f_3 * hs_12[k]
                  + pb_x[k] * hp_29[k];

        t_84[k] = f_2 * hd_s_85[k]
                  + pb_x[k] * hp_31[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pb_y, pb_z, gp_21, hs_s_13, hd_s_86, hd_s_87, \
                         hd_s_88, hs_12, hp_30, hp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -f_1 * hs_s_13[k]
                  + f_2 * hd_s_86[k]
                  + f_3 * hs_12[k]
                  + pb_y[k] * hp_30[k];

        t_86[k] = f_2 * hd_s_87[k]
                  + pb_y[k] * hp_31[k];

        t_87[k] = f_0 * gp_21[k]
                  - f_1 * hs_s_13[k]
                  + f_2 * hd_s_88[k]
                  + f_3 * hs_12[k]
                  + pb_z[k] * hp_31[k];
    }
}

auto
compute_prim_hd_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 * beta / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 1.0 / p;

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

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_17 = buffer.data(fd_s + 17);
    const auto *fd_s_18 = buffer.data(fd_s + 18);
    const auto *fd_s_22 = buffer.data(fd_s + 22);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_22 = buffer.data(fd + 22);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_40 = buffer.data(gd + 40);

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_2 = buffer.data(hs_s + 2);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_6 = buffer.data(hs_s + 6);
    const auto *hs_s_8 = buffer.data(hs_s + 8);
    const auto *hs_s_9 = buffer.data(hs_s + 9);
    const auto *hs_s_11 = buffer.data(hs_s + 11);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_4 = buffer.data(hd_s + 4);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_21 = buffer.data(hd_s + 21);
    const auto *hd_s_22 = buffer.data(hd_s + 22);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_24 = buffer.data(hd_s + 24);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_28 = buffer.data(hd_s + 28);
    const auto *hd_s_29 = buffer.data(hd_s + 29);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_38 = buffer.data(hd_s + 38);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_43 = buffer.data(hd_s + 43);
    const auto *hd_s_44 = buffer.data(hd_s + 44);
    const auto *hd_s_45 = buffer.data(hd_s + 45);
    const auto *hd_s_46 = buffer.data(hd_s + 46);
    const auto *hd_s_47 = buffer.data(hd_s + 47);
    const auto *hd_s_48 = buffer.data(hd_s + 48);
    const auto *hd_s_49 = buffer.data(hd_s + 49);
    const auto *hd_s_50 = buffer.data(hd_s + 50);
    const auto *hd_s_51 = buffer.data(hd_s + 51);
    const auto *hd_s_52 = buffer.data(hd_s + 52);
    const auto *hd_s_53 = buffer.data(hd_s + 53);
    const auto *hd_s_54 = buffer.data(hd_s + 54);
    const auto *hd_s_55 = buffer.data(hd_s + 55);
    const auto *hd_s_57 = buffer.data(hd_s + 57);
    const auto *hd_s_59 = buffer.data(hd_s + 59);
    const auto *hd_s_60 = buffer.data(hd_s + 60);
    const auto *hd_s_61 = buffer.data(hd_s + 61);
    const auto *hd_s_62 = buffer.data(hd_s + 62);
    const auto *hd_s_63 = buffer.data(hd_s + 63);
    const auto *hd_s_64 = buffer.data(hd_s + 64);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hs_0, hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_1[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_2[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, fd_s_4, fd_4, gd_0, gd_2, gd_4, \
                         hd_s_3, hd_s_4, hd_s_5, hd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * gd_0[k]
                 + f_2 * hd_s_3[k];

        t_4[k] = -f_4 * fd_s_4[k]
                 + f_5 * fd_4[k]
                 + pa_x[k] * gd_4[k]
                 + f_2 * hd_s_4[k];

        t_5[k] = pa_y[k] * gd_2[k]
                 + f_2 * hd_s_5[k];

        t_6[k] = pa_z[k] * gd_0[k]
                 + f_2 * hd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pb_y, fd_s_0, fd_s_7, fd_0, fd_7, gd_3, \
                         gd_8, hd_s_7, hd_s_8, hd_s_9, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * hd_s_7[k]
                 + pb_y[k] * hp_3[k];

        t_8[k] = -f_4 * fd_s_7[k]
                 + f_5 * fd_7[k]
                 + pa_x[k] * gd_8[k]
                 + f_2 * hd_s_8[k];

        t_9[k] = -f_6 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_y[k] * gd_3[k]
                 + f_2 * hd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_z, fd_s_9, fd_9, gd_4, gd_10, \
                         hs_s_2, hd_s_10, hd_s_11, hd_s_12, hs_2, \
                         hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_7 * fd_s_9[k]
                  + f_8 * fd_9[k]
                  + pa_x[k] * gd_10[k]
                  + f_2 * hd_s_10[k];

        t_11[k] = -f_1 * hs_s_2[k]
                  + f_2 * hd_s_11[k]
                  + f_3 * hs_2[k]
                  + pb_z[k] * hp_4[k];

        t_12[k] = pa_z[k] * gd_4[k]
                  + f_2 * hd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_y, fd_s_0, fd_0, gd_6, gd_8, hs_s_3, \
                         hd_s_13, hd_s_14, hd_s_15, hs_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * gd_8[k]
                  + f_2 * hd_s_13[k];

        t_14[k] = -f_6 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_z[k] * gd_6[k]
                  + f_2 * hd_s_14[k];

        t_15[k] = -f_1 * hs_s_3[k]
                  + f_2 * hd_s_15[k]
                  + f_3 * hs_3[k]
                  + pb_y[k] * hp_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_y, fd_s_3, fd_s_11, fd_3, fd_11, \
                         gd_9, gd_17, hd_s_16, hd_s_17, hd_s_18, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * hd_s_16[k]
                  + pb_y[k] * hp_6[k];

        t_17[k] = -f_7 * fd_s_11[k]
                  + f_8 * fd_11[k]
                  + pa_x[k] * gd_17[k]
                  + f_2 * hd_s_17[k];

        t_18[k] = -f_7 * fd_s_3[k]
                  + f_8 * fd_3[k]
                  + pa_y[k] * gd_9[k]
                  + f_2 * hd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_z, pb_z, fd_s_13, fd_13, gd_9, gd_19, \
                         hs_s_4, hd_s_19, hd_s_20, hd_s_21, hs_4, \
                         hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_6 * fd_s_13[k]
                  + f_3 * fd_13[k]
                  + pa_x[k] * gd_19[k]
                  + f_2 * hd_s_19[k];

        t_20[k] = -f_1 * hs_s_4[k]
                  + f_2 * hd_s_20[k]
                  + f_3 * hs_4[k]
                  + pb_z[k] * hp_7[k];

        t_21[k] = pa_z[k] * gd_9[k]
                  + f_2 * hd_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_z, fd_s_16, fd_s_17, fd_16, fd_17, gd_10, \
                         gd_21, gd_22, hd_s_22, hd_s_23, hd_s_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_z[k] * gd_10[k]
                  + f_2 * hd_s_22[k];

        t_23[k] = -f_6 * fd_s_16[k]
                  + f_3 * fd_16[k]
                  + pa_x[k] * gd_21[k]
                  + f_2 * hd_s_23[k];

        t_24[k] = -f_6 * fd_s_17[k]
                  + f_3 * fd_17[k]
                  + pa_x[k] * gd_22[k]
                  + f_2 * hd_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_y, pa_z, pb_y, fd_s_6, fd_6, gd_14, gd_17, \
                         hs_s_5, hd_s_25, hd_s_26, hd_s_27, hs_5, \
                         hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_y[k] * gd_17[k]
                  + f_2 * hd_s_25[k];

        t_26[k] = -f_7 * fd_s_6[k]
                  + f_8 * fd_6[k]
                  + pa_z[k] * gd_14[k]
                  + f_2 * hd_s_26[k];

        t_27[k] = -f_1 * hs_s_5[k]
                  + f_2 * hd_s_27[k]
                  + f_3 * hs_5[k]
                  + pb_y[k] * hp_8[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pb_y, fd_s_22, fd_22, gp_7, gd_24, gd_25, \
                         hd_s_28, hd_s_29, hd_s_30, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_2 * hd_s_28[k]
                  + pb_y[k] * hp_9[k];

        t_29[k] = -f_6 * fd_s_22[k]
                  + f_3 * fd_22[k]
                  + pa_x[k] * gd_24[k]
                  + f_2 * hd_s_29[k];

        t_30[k] = f_8 * gp_7[k]
                  + pa_x[k] * gd_25[k]
                  + f_2 * hd_s_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pa_z, gp_10, gp_13, gd_18, gd_27, \
                         gd_31, gd_36, hd_s_31, hd_s_32, hd_s_33, \
                         hd_s_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pa_x[k] * gd_27[k]
                  + f_2 * hd_s_31[k];

        t_32[k] = pa_z[k] * gd_18[k]
                  + f_2 * hd_s_32[k];

        t_33[k] = f_8 * gp_10[k]
                  + pa_x[k] * gd_31[k]
                  + f_2 * hd_s_33[k];

        t_34[k] = f_8 * gp_13[k]
                  + pa_x[k] * gd_36[k]
                  + f_2 * hd_s_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pb_x, gd_40, hs_s_6, hd_s_35, hd_s_36, \
                         hd_s_37, hs_6, hp_10, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_x[k] * gd_40[k]
                  + f_2 * hd_s_35[k];

        t_36[k] = -f_1 * hs_s_6[k]
                  + f_2 * hd_s_36[k]
                  + f_3 * hs_6[k]
                  + pb_x[k] * hp_10[k];

        t_37[k] = f_2 * hd_s_37[k]
                  + pb_x[k] * hp_11[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, pb_y, pb_z, gp_8, hs_s_6, hd_s_38, hd_s_39, \
                         hd_s_40, hs_6, hp_11, hp_12, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * gp_8[k]
                  - f_1 * hs_s_6[k]
                  + f_2 * hd_s_38[k]
                  + f_3 * hs_6[k]
                  + pb_y[k] * hp_11[k];

        t_39[k] = -f_1 * hs_s_6[k]
                  + f_2 * hd_s_39[k]
                  + f_3 * hs_6[k]
                  + pb_z[k] * hp_12[k];

        t_40[k] = f_2 * hd_s_40[k]
                  + pb_x[k] * hp_13[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pa_z, pb_x, fd_s_16, fd_16, gd_27, gd_30, \
                         hs_s_8, hd_s_41, hd_s_43, hd_s_44, hs_8, \
                         hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_z[k] * gd_27[k]
                  + f_2 * hd_s_41[k];

        t_42[k] = -f_4 * fd_s_16[k]
                  + f_5 * fd_16[k]
                  + pa_y[k] * gd_30[k]
                  + f_2 * hd_s_43[k];

        t_43[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_44[k]
                  + f_3 * hs_8[k]
                  + pb_x[k] * hp_14[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_z, pb_x, fd_s_13, fd_13, gd_29, hd_s_45, \
                         hd_s_46, hd_s_47, hp_15, hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_2 * hd_s_45[k]
                  + pb_x[k] * hp_15[k];

        t_45[k] = f_2 * hd_s_46[k]
                  + pb_x[k] * hp_16[k];

        t_46[k] = -f_6 * fd_s_13[k]
                  + f_3 * fd_13[k]
                  + pa_z[k] * gd_29[k]
                  + f_2 * hd_s_47[k];
    }

#pragma omp simd aligned(t_47, t_48, pa_y, pb_y, fd_s_18, fd_18, gp_11, gd_33, hd_s_48, \
                         hd_s_49, hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * gp_11[k]
                  + f_2 * hd_s_48[k]
                  + pb_y[k] * hp_16[k];

        t_48[k] = -f_7 * fd_s_18[k]
                  + f_8 * fd_18[k]
                  + pa_y[k] * gd_33[k]
                  + f_2 * hd_s_49[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, hs_s_9, hd_s_50, hd_s_51, hd_s_52, hs_9, \
                         hp_17, hp_18, hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = -f_1 * hs_s_9[k]
                  + f_2 * hd_s_50[k]
                  + f_3 * hs_9[k]
                  + pb_x[k] * hp_17[k];

        t_50[k] = f_2 * hd_s_51[k]
                  + pb_x[k] * hp_18[k];

        t_51[k] = f_2 * hd_s_52[k]
                  + pb_x[k] * hp_19[k];
    }

#pragma omp simd aligned(t_52, t_53, pa_z, pb_y, fd_s_15, fd_15, gp_12, gd_32, hd_s_53, \
                         hd_s_54, hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -f_7 * fd_s_15[k]
                  + f_8 * fd_15[k]
                  + pa_z[k] * gd_32[k]
                  + f_2 * hd_s_53[k];

        t_53[k] = f_8 * gp_12[k]
                  + f_2 * hd_s_54[k]
                  + pb_y[k] * hp_19[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, fd_s_22, fd_22, gp_14, gd_35, gd_38, gd_40, \
                         hd_s_55, hd_s_57, hd_s_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = -f_6 * fd_s_22[k]
                  + f_3 * fd_22[k]
                  + pa_y[k] * gd_35[k]
                  + f_2 * hd_s_55[k];

        t_55[k] = f_8 * gp_14[k]
                  + pa_y[k] * gd_38[k]
                  + f_2 * hd_s_57[k];

        t_56[k] = pa_y[k] * gd_40[k]
                  + f_2 * hd_s_59[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_x, pb_y, hs_s_11, hd_s_60, hd_s_61, \
                         hd_s_62, hd_s_63, hs_10, hp_20, hp_21, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_60[k]
                  + f_3 * hs_10[k]
                  + pb_x[k] * hp_20[k];

        t_58[k] = f_2 * hd_s_61[k]
                  + pb_x[k] * hp_22[k];

        t_59[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_62[k]
                  + f_3 * hs_10[k]
                  + pb_y[k] * hp_21[k];

        t_60[k] = f_2 * hd_s_63[k]
                  + pb_y[k] * hp_22[k];
    }

#pragma omp simd aligned(t_61, pb_z, gp_15, hs_s_11, hd_s_64, hs_10, \
                         hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * gp_15[k]
                  - f_1 * hs_s_11[k]
                  + f_2 * hd_s_64[k]
                  + f_3 * hs_10[k]
                  + pb_z[k] * hp_22[k];
    }
}

auto
compute_prim_hd_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 * beta / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_19 = buffer.data(fd_s + 19);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_19 = buffer.data(fd + 19);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_15 = buffer.data(gp + 15);

    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_27 = buffer.data(gd + 27);

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_2 = buffer.data(hs_s + 2);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_6 = buffer.data(hs_s + 6);
    const auto *hs_s_8 = buffer.data(hs_s + 8);
    const auto *hs_s_9 = buffer.data(hs_s + 9);
    const auto *hs_s_11 = buffer.data(hs_s + 11);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_24 = buffer.data(hd_s + 24);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_38 = buffer.data(hd_s + 38);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_42 = buffer.data(hd_s + 42);
    const auto *hd_s_47 = buffer.data(hd_s + 47);
    const auto *hd_s_48 = buffer.data(hd_s + 48);
    const auto *hd_s_49 = buffer.data(hd_s + 49);
    const auto *hd_s_50 = buffer.data(hd_s + 50);
    const auto *hd_s_51 = buffer.data(hd_s + 51);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hs_0, hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_1[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_2[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, fd_s_4, fd_s_6, fd_4, fd_6, gd_4, gd_7, \
                         hd_s_3, hd_s_5, hd_s_6, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * fd_s_4[k]
                 + f_5 * fd_4[k]
                 + pa_x[k] * gd_4[k]
                 + f_2 * hd_s_3[k];

        t_4[k] = f_2 * hd_s_5[k]
                 + pb_y[k] * hp_3[k];

        t_5[k] = -f_4 * fd_s_6[k]
                 + f_5 * fd_6[k]
                 + pa_x[k] * gd_7[k]
                 + f_2 * hd_s_6[k];
    }

#pragma omp simd aligned(t_6, t_7, pa_x, pa_y, fd_s_0, fd_s_7, fd_0, fd_7, gd_3, gd_9, hd_s_7, \
                         hd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_6 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_y[k] * gd_3[k]
                 + f_2 * hd_s_7[k];

        t_7[k] = -f_7 * fd_s_7[k]
                 + f_8 * fd_7[k]
                 + pa_x[k] * gd_9[k]
                 + f_2 * hd_s_8[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_z, pb_z, fd_s_0, fd_0, gd_5, hs_s_2, hd_s_9, hd_s_10, \
                         hs_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_1 * hs_s_2[k]
                 + f_2 * hd_s_9[k]
                 + f_3 * hs_2[k]
                 + pb_z[k] * hp_4[k];

        t_9[k] = -f_6 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_z[k] * gd_5[k]
                 + f_2 * hd_s_10[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_y, fd_s_8, fd_8, gd_14, hs_s_3, hd_s_11, \
                         hd_s_12, hd_s_13, hs_3, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * hs_s_3[k]
                  + f_2 * hd_s_11[k]
                  + f_3 * hs_3[k]
                  + pb_y[k] * hp_5[k];

        t_11[k] = f_2 * hd_s_12[k]
                  + pb_y[k] * hp_6[k];

        t_12[k] = -f_7 * fd_s_8[k]
                  + f_8 * fd_8[k]
                  + pa_x[k] * gd_14[k]
                  + f_2 * hd_s_13[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, fd_s_3, fd_s_10, fd_3, fd_10, gd_8, gd_15, \
                         hd_s_14, hd_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_7 * fd_s_3[k]
                  + f_8 * fd_3[k]
                  + pa_y[k] * gd_8[k]
                  + f_2 * hd_s_14[k];

        t_14[k] = -f_6 * fd_s_10[k]
                  + f_3 * fd_10[k]
                  + pa_x[k] * gd_15[k]
                  + f_2 * hd_s_15[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_z, pb_z, fd_s_5, fd_5, gd_11, hs_s_4, hd_s_16, \
                         hd_s_17, hs_4, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * hs_s_4[k]
                  + f_2 * hd_s_16[k]
                  + f_3 * hs_4[k]
                  + pb_z[k] * hp_7[k];

        t_16[k] = -f_7 * fd_s_5[k]
                  + f_8 * fd_5[k]
                  + pa_z[k] * gd_11[k]
                  + f_2 * hd_s_17[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pb_y, fd_s_19, fd_19, gd_16, hs_s_5, hd_s_18, \
                         hd_s_19, hd_s_20, hs_5, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_1 * hs_s_5[k]
                  + f_2 * hd_s_18[k]
                  + f_3 * hs_5[k]
                  + pb_y[k] * hp_8[k];

        t_18[k] = f_2 * hd_s_19[k]
                  + pb_y[k] * hp_9[k];

        t_19[k] = -f_6 * fd_s_19[k]
                  + f_3 * fd_19[k]
                  + pa_x[k] * gd_16[k]
                  + f_2 * hd_s_20[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_y, gp_8, hs_s_6, hd_s_23, hd_s_24, \
                         hd_s_25, hs_6, hp_10, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * hs_s_6[k]
                  + f_2 * hd_s_23[k]
                  + f_3 * hs_6[k]
                  + pb_x[k] * hp_10[k];

        t_21[k] = f_2 * hd_s_24[k]
                  + pb_x[k] * hp_11[k];

        t_22[k] = f_0 * gp_8[k]
                  - f_1 * hs_s_6[k]
                  + f_2 * hd_s_25[k]
                  + f_3 * hs_6[k]
                  + pb_y[k] * hp_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_x, pb_z, fd_s_13, fd_13, gd_22, hs_s_6, \
                         hd_s_26, hd_s_27, hd_s_30, hs_6, hp_12, \
                         hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_1 * hs_s_6[k]
                  + f_2 * hd_s_26[k]
                  + f_3 * hs_6[k]
                  + pb_z[k] * hp_12[k];

        t_24[k] = f_2 * hd_s_27[k]
                  + pb_x[k] * hp_13[k];

        t_25[k] = -f_4 * fd_s_13[k]
                  + f_5 * fd_13[k]
                  + pa_y[k] * gd_22[k]
                  + f_2 * hd_s_30[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_x, hs_s_8, hd_s_31, hd_s_32, hd_s_33, hs_8, \
                         hp_14, hp_15, hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_31[k]
                  + f_3 * hs_8[k]
                  + pb_x[k] * hp_14[k];

        t_27[k] = f_2 * hd_s_32[k]
                  + pb_x[k] * hp_15[k];

        t_28[k] = f_2 * hd_s_33[k]
                  + pb_x[k] * hp_16[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_z, pb_y, fd_s_10, fd_10, gp_11, gd_21, hd_s_34, \
                         hd_s_35, hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_6 * fd_s_10[k]
                  + f_3 * fd_10[k]
                  + pa_z[k] * gd_21[k]
                  + f_2 * hd_s_34[k];

        t_30[k] = f_5 * gp_11[k]
                  + f_2 * hd_s_35[k]
                  + pb_y[k] * hp_16[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pb_x, fd_s_15, fd_15, gd_25, hs_s_9, hd_s_36, \
                         hd_s_37, hd_s_38, hs_9, hp_17, hp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_7 * fd_s_15[k]
                  + f_8 * fd_15[k]
                  + pa_y[k] * gd_25[k]
                  + f_2 * hd_s_36[k];

        t_32[k] = -f_1 * hs_s_9[k]
                  + f_2 * hd_s_37[k]
                  + f_3 * hs_9[k]
                  + pb_x[k] * hp_17[k];

        t_33[k] = f_2 * hd_s_38[k]
                  + pb_x[k] * hp_18[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_x, pb_y, fd_s_12, fd_12, gp_12, gd_24, \
                         hd_s_39, hd_s_40, hd_s_41, hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * hd_s_39[k]
                  + pb_x[k] * hp_19[k];

        t_35[k] = -f_7 * fd_s_12[k]
                  + f_8 * fd_12[k]
                  + pa_z[k] * gd_24[k]
                  + f_2 * hd_s_40[k];

        t_36[k] = f_8 * gp_12[k]
                  + f_2 * hd_s_41[k]
                  + pb_y[k] * hp_19[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pb_x, fd_s_19, fd_19, gd_27, hs_s_11, \
                         hd_s_42, hd_s_47, hd_s_48, hs_10, hp_20, \
                         hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_6 * fd_s_19[k]
                  + f_3 * fd_19[k]
                  + pa_y[k] * gd_27[k]
                  + f_2 * hd_s_42[k];

        t_38[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_47[k]
                  + f_3 * hs_10[k]
                  + pb_x[k] * hp_20[k];

        t_39[k] = f_2 * hd_s_48[k]
                  + pb_x[k] * hp_22[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_y, pb_z, gp_15, hs_s_11, hd_s_49, hd_s_50, \
                         hd_s_51, hs_10, hp_21, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_49[k]
                  + f_3 * hs_10[k]
                  + pb_y[k] * hp_21[k];

        t_41[k] = f_2 * hd_s_50[k]
                  + pb_y[k] * hp_22[k];

        t_42[k] = f_0 * gp_15[k]
                  - f_1 * hs_s_11[k]
                  + f_2 * hd_s_51[k]
                  + f_3 * hs_10[k]
                  + pb_z[k] * hp_22[k];
    }
}

auto
compute_prim_hd_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 * beta / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_19 = buffer.data(fd_s + 19);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_18 = buffer.data(fd + 18);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_15 = buffer.data(gp + 15);

    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_2 = buffer.data(hs_s + 2);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_6 = buffer.data(hs_s + 6);
    const auto *hs_s_8 = buffer.data(hs_s + 8);
    const auto *hs_s_9 = buffer.data(hs_s + 9);
    const auto *hs_s_11 = buffer.data(hs_s + 11);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_24 = buffer.data(hd_s + 24);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_38 = buffer.data(hd_s + 38);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_42 = buffer.data(hd_s + 42);
    const auto *hd_s_47 = buffer.data(hd_s + 47);
    const auto *hd_s_48 = buffer.data(hd_s + 48);
    const auto *hd_s_49 = buffer.data(hd_s + 49);
    const auto *hd_s_50 = buffer.data(hd_s + 50);
    const auto *hd_s_51 = buffer.data(hd_s + 51);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hs_0, hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_1[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_2[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, fd_s_4, fd_s_6, fd_4, fd_6, gd_4, gd_7, \
                         hd_s_3, hd_s_5, hd_s_6, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * fd_s_4[k]
                 + f_5 * fd_4[k]
                 + pa_x[k] * gd_4[k]
                 + f_2 * hd_s_3[k];

        t_4[k] = f_2 * hd_s_5[k]
                 + pb_y[k] * hp_3[k];

        t_5[k] = -f_4 * fd_s_6[k]
                 + f_5 * fd_6[k]
                 + pa_x[k] * gd_7[k]
                 + f_2 * hd_s_6[k];
    }

#pragma omp simd aligned(t_6, t_7, pa_x, pa_y, fd_s_0, fd_s_7, fd_0, fd_7, gd_3, gd_9, hd_s_7, \
                         hd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_6 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_y[k] * gd_3[k]
                 + f_2 * hd_s_7[k];

        t_7[k] = -f_7 * fd_s_7[k]
                 + f_8 * fd_7[k]
                 + pa_x[k] * gd_9[k]
                 + f_2 * hd_s_8[k];
    }

#pragma omp simd aligned(t_8, t_9, pa_z, pb_z, fd_s_0, fd_0, gd_5, hs_s_2, hd_s_9, hd_s_10, \
                         hs_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_1 * hs_s_2[k]
                 + f_2 * hd_s_9[k]
                 + f_3 * hs_2[k]
                 + pb_z[k] * hp_4[k];

        t_9[k] = -f_6 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_z[k] * gd_5[k]
                 + f_2 * hd_s_10[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_y, fd_s_8, fd_8, gd_14, hs_s_3, hd_s_11, \
                         hd_s_12, hd_s_13, hs_3, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * hs_s_3[k]
                  + f_2 * hd_s_11[k]
                  + f_3 * hs_3[k]
                  + pb_y[k] * hp_5[k];

        t_11[k] = f_2 * hd_s_12[k]
                  + pb_y[k] * hp_6[k];

        t_12[k] = -f_7 * fd_s_8[k]
                  + f_8 * fd_8[k]
                  + pa_x[k] * gd_14[k]
                  + f_2 * hd_s_13[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, fd_s_3, fd_s_10, fd_3, fd_10, gd_8, gd_15, \
                         hd_s_14, hd_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -f_7 * fd_s_3[k]
                  + f_8 * fd_3[k]
                  + pa_y[k] * gd_8[k]
                  + f_2 * hd_s_14[k];

        t_14[k] = -f_6 * fd_s_10[k]
                  + f_3 * fd_10[k]
                  + pa_x[k] * gd_15[k]
                  + f_2 * hd_s_15[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_z, pb_z, fd_s_5, fd_5, gd_11, hs_s_4, hd_s_16, \
                         hd_s_17, hs_4, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * hs_s_4[k]
                  + f_2 * hd_s_16[k]
                  + f_3 * hs_4[k]
                  + pb_z[k] * hp_7[k];

        t_16[k] = -f_7 * fd_s_5[k]
                  + f_8 * fd_5[k]
                  + pa_z[k] * gd_11[k]
                  + f_2 * hd_s_17[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pb_y, fd_s_19, fd_18, gd_16, hs_s_5, hd_s_18, \
                         hd_s_19, hd_s_20, hs_5, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_1 * hs_s_5[k]
                  + f_2 * hd_s_18[k]
                  + f_3 * hs_5[k]
                  + pb_y[k] * hp_8[k];

        t_18[k] = f_2 * hd_s_19[k]
                  + pb_y[k] * hp_9[k];

        t_19[k] = -f_6 * fd_s_19[k]
                  + f_3 * fd_18[k]
                  + pa_x[k] * gd_16[k]
                  + f_2 * hd_s_20[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_y, gp_8, hs_s_6, hd_s_23, hd_s_24, \
                         hd_s_25, hs_6, hp_10, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * hs_s_6[k]
                  + f_2 * hd_s_23[k]
                  + f_3 * hs_6[k]
                  + pb_x[k] * hp_10[k];

        t_21[k] = f_2 * hd_s_24[k]
                  + pb_x[k] * hp_11[k];

        t_22[k] = f_0 * gp_8[k]
                  - f_1 * hs_s_6[k]
                  + f_2 * hd_s_25[k]
                  + f_3 * hs_6[k]
                  + pb_y[k] * hp_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pb_x, pb_z, fd_s_13, fd_13, gd_22, hs_s_6, \
                         hd_s_26, hd_s_27, hd_s_30, hs_6, hp_12, \
                         hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_1 * hs_s_6[k]
                  + f_2 * hd_s_26[k]
                  + f_3 * hs_6[k]
                  + pb_z[k] * hp_12[k];

        t_24[k] = f_2 * hd_s_27[k]
                  + pb_x[k] * hp_13[k];

        t_25[k] = -f_4 * fd_s_13[k]
                  + f_5 * fd_13[k]
                  + pa_y[k] * gd_22[k]
                  + f_2 * hd_s_30[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_x, hs_s_8, hd_s_31, hd_s_32, hd_s_33, hs_8, \
                         hp_14, hp_15, hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_31[k]
                  + f_3 * hs_8[k]
                  + pb_x[k] * hp_14[k];

        t_27[k] = f_2 * hd_s_32[k]
                  + pb_x[k] * hp_15[k];

        t_28[k] = f_2 * hd_s_33[k]
                  + pb_x[k] * hp_16[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_z, pb_y, fd_s_10, fd_10, gp_11, gd_21, hd_s_34, \
                         hd_s_35, hp_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_6 * fd_s_10[k]
                  + f_3 * fd_10[k]
                  + pa_z[k] * gd_21[k]
                  + f_2 * hd_s_34[k];

        t_30[k] = f_5 * gp_11[k]
                  + f_2 * hd_s_35[k]
                  + pb_y[k] * hp_16[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_y, pb_x, fd_s_15, fd_14, gd_25, hs_s_9, hd_s_36, \
                         hd_s_37, hd_s_38, hs_9, hp_17, hp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_7 * fd_s_15[k]
                  + f_8 * fd_14[k]
                  + pa_y[k] * gd_25[k]
                  + f_2 * hd_s_36[k];

        t_32[k] = -f_1 * hs_s_9[k]
                  + f_2 * hd_s_37[k]
                  + f_3 * hs_9[k]
                  + pb_x[k] * hp_17[k];

        t_33[k] = f_2 * hd_s_38[k]
                  + pb_x[k] * hp_18[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, pb_x, pb_y, fd_s_12, fd_12, gp_12, gd_24, \
                         hd_s_39, hd_s_40, hd_s_41, hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * hd_s_39[k]
                  + pb_x[k] * hp_19[k];

        t_35[k] = -f_7 * fd_s_12[k]
                  + f_8 * fd_12[k]
                  + pa_z[k] * gd_24[k]
                  + f_2 * hd_s_40[k];

        t_36[k] = f_8 * gp_12[k]
                  + f_2 * hd_s_41[k]
                  + pb_y[k] * hp_19[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pb_x, fd_s_19, fd_18, gd_26, hs_s_11, \
                         hd_s_42, hd_s_47, hd_s_48, hs_10, hp_20, \
                         hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_6 * fd_s_19[k]
                  + f_3 * fd_18[k]
                  + pa_y[k] * gd_26[k]
                  + f_2 * hd_s_42[k];

        t_38[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_47[k]
                  + f_3 * hs_10[k]
                  + pb_x[k] * hp_20[k];

        t_39[k] = f_2 * hd_s_48[k]
                  + pb_x[k] * hp_22[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_y, pb_z, gp_15, hs_s_11, hd_s_49, hd_s_50, \
                         hd_s_51, hs_10, hp_21, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_49[k]
                  + f_3 * hs_10[k]
                  + pb_y[k] * hp_21[k];

        t_41[k] = f_2 * hd_s_50[k]
                  + pb_y[k] * hp_22[k];

        t_42[k] = f_0 * gp_15[k]
                  - f_1 * hs_s_11[k]
                  + f_2 * hd_s_51[k]
                  + f_3 * hs_10[k]
                  + pb_z[k] * hp_22[k];
    }
}

auto
compute_prim_hd_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 * beta / p;
    const auto f_6 = 1.5 / p;
    const auto f_7 = beta / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = 1.0 / p;

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

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_13 = buffer.data(fd_s + 13);

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
    const auto *fd_13 = buffer.data(fd + 13);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);

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
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_12 = buffer.data(hs_s + 12);
    const auto *hs_s_14 = buffer.data(hs_s + 14);
    const auto *hs_s_15 = buffer.data(hs_s + 15);
    const auto *hs_s_17 = buffer.data(hs_s + 17);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_4 = buffer.data(hd_s + 4);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_21 = buffer.data(hd_s + 21);
    const auto *hd_s_22 = buffer.data(hd_s + 22);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_24 = buffer.data(hd_s + 24);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_28 = buffer.data(hd_s + 28);
    const auto *hd_s_29 = buffer.data(hd_s + 29);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_38 = buffer.data(hd_s + 38);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_42 = buffer.data(hd_s + 42);
    const auto *hd_s_43 = buffer.data(hd_s + 43);
    const auto *hd_s_44 = buffer.data(hd_s + 44);
    const auto *hd_s_45 = buffer.data(hd_s + 45);
    const auto *hd_s_46 = buffer.data(hd_s + 46);
    const auto *hd_s_47 = buffer.data(hd_s + 47);
    const auto *hd_s_48 = buffer.data(hd_s + 48);
    const auto *hd_s_49 = buffer.data(hd_s + 49);
    const auto *hd_s_50 = buffer.data(hd_s + 50);
    const auto *hd_s_51 = buffer.data(hd_s + 51);
    const auto *hd_s_52 = buffer.data(hd_s + 52);
    const auto *hd_s_53 = buffer.data(hd_s + 53);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_17 = buffer.data(hs + 17);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hs_0, hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_1[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_2[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pb_x, fd_s_2, fd_2, gp_1, gd_0, gd_2, \
                         hd_s_3, hd_s_4, hd_s_5, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * gd_0[k]
                 + f_2 * hd_s_3[k];

        t_4[k] = f_4 * gp_1[k]
                 + f_2 * hd_s_4[k]
                 + pb_x[k] * hp_3[k];

        t_5[k] = -f_5 * fd_s_2[k]
                 + f_6 * fd_2[k]
                 + pa_x[k] * gd_2[k]
                 + f_2 * hd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_z, pb_x, fd_s_4, fd_4, gp_2, gd_0, gd_4, \
                         hd_s_6, hd_s_7, hd_s_8, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_z[k] * gd_0[k]
                 + f_2 * hd_s_6[k];

        t_7[k] = f_4 * gp_2[k]
                 + f_2 * hd_s_7[k]
                 + pb_x[k] * hp_4[k];

        t_8[k] = -f_5 * fd_s_4[k]
                 + f_6 * fd_4[k]
                 + pa_x[k] * gd_4[k]
                 + f_2 * hd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_y, pb_x, fd_s_0, fd_0, gp_3, gd_1, hd_s_9, hd_s_10, \
                         hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_7 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_y[k] * gd_1[k]
                 + f_2 * hd_s_9[k];

        t_10[k] = f_6 * gp_3[k]
                  + f_2 * hd_s_10[k]
                  + pb_x[k] * hp_5[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, fd_s_0, fd_s_5, fd_0, fd_5, gd_3, gd_6, \
                         hd_s_11, hd_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -f_8 * fd_s_5[k]
                  + f_9 * fd_5[k]
                  + pa_x[k] * gd_6[k]
                  + f_2 * hd_s_11[k];

        t_12[k] = -f_7 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_z[k] * gd_3[k]
                  + f_2 * hd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pb_x, fd_s_6, fd_6, gp_4, gd_8, hd_s_13, hd_s_14, \
                         hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * gp_4[k]
                  + f_2 * hd_s_13[k]
                  + pb_x[k] * hp_7[k];

        t_14[k] = -f_8 * fd_s_6[k]
                  + f_9 * fd_6[k]
                  + pa_x[k] * gd_8[k]
                  + f_2 * hd_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_y, pb_x, fd_s_1, fd_1, gp_5, gd_5, hd_s_15, hd_s_16, \
                         hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_8 * fd_s_1[k]
                  + f_9 * fd_1[k]
                  + pa_y[k] * gd_5[k]
                  + f_2 * hd_s_15[k];

        t_16[k] = f_9 * gp_5[k]
                  + f_2 * hd_s_16[k]
                  + pb_x[k] * hp_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_x, pa_y, fd_s_7, fd_s_9, fd_7, fd_9, gd_7, gd_9, \
                         gd_10, hd_s_17, hd_s_18, hd_s_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_7 * fd_s_7[k]
                  + f_3 * fd_7[k]
                  + pa_x[k] * gd_9[k]
                  + f_2 * hd_s_17[k];

        t_18[k] = -f_7 * fd_s_9[k]
                  + f_3 * fd_9[k]
                  + pa_x[k] * gd_10[k]
                  + f_2 * hd_s_18[k];

        t_19[k] = pa_y[k] * gd_7[k]
                  + f_2 * hd_s_19[k];
    }

#pragma omp simd aligned(t_20, t_21, pa_x, pa_z, fd_s_3, fd_s_10, fd_3, fd_10, gd_7, gd_11, \
                         hd_s_20, hd_s_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_7 * fd_s_10[k]
                  + f_3 * fd_10[k]
                  + pa_x[k] * gd_11[k]
                  + f_2 * hd_s_20[k];

        t_21[k] = -f_8 * fd_s_3[k]
                  + f_9 * fd_3[k]
                  + pa_z[k] * gd_7[k]
                  + f_2 * hd_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pb_x, fd_s_13, fd_13, gp_6, gp_7, gd_12, \
                         gd_13, hd_s_22, hd_s_23, hd_s_24, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_9 * gp_6[k]
                  + f_2 * hd_s_22[k]
                  + pb_x[k] * hp_12[k];

        t_23[k] = -f_7 * fd_s_13[k]
                  + f_3 * fd_13[k]
                  + pa_x[k] * gd_12[k]
                  + f_2 * hd_s_23[k];

        t_24[k] = f_9 * gp_7[k]
                  + pa_x[k] * gd_13[k]
                  + f_2 * hd_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pb_x, gp_8, gd_14, gd_16, gd_17, \
                         hd_s_25, hd_s_26, hd_s_27, hd_s_28, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * gp_8[k]
                  + f_2 * hd_s_25[k]
                  + pb_x[k] * hp_13[k];

        t_26[k] = pa_x[k] * gd_14[k]
                  + f_2 * hd_s_26[k];

        t_27[k] = pa_x[k] * gd_16[k]
                  + f_2 * hd_s_27[k];

        t_28[k] = pa_x[k] * gd_17[k]
                  + f_2 * hd_s_28[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pb_x, gp_14, gp_16, gd_18, gd_19, \
                         gd_21, hd_s_29, hd_s_30, hd_s_31, hd_s_32, \
                         hp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_x[k] * gd_18[k]
                  + f_2 * hd_s_29[k];

        t_30[k] = pa_x[k] * gd_19[k]
                  + f_2 * hd_s_30[k];

        t_31[k] = f_9 * gp_14[k]
                  + pa_x[k] * gd_21[k]
                  + f_2 * hd_s_31[k];

        t_32[k] = f_3 * gp_16[k]
                  + f_2 * hd_s_32[k]
                  + pb_x[k] * hp_18[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_x, pb_x, pb_y, gp_8, gd_23, hs_s_12, hd_s_33, \
                         hd_s_34, hd_s_35, hs_12, hp_19, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_x[k] * gd_23[k]
                  + f_2 * hd_s_33[k];

        t_34[k] = -f_1 * hs_s_12[k]
                  + f_2 * hd_s_34[k]
                  + f_3 * hs_12[k]
                  + pb_x[k] * hp_19[k];

        t_35[k] = f_0 * gp_8[k]
                  - f_1 * hs_s_12[k]
                  + f_2 * hd_s_35[k]
                  + f_3 * hs_12[k]
                  + pb_y[k] * hp_20[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_z, pb_y, pb_z, gp_9, gd_14, hs_s_12, hd_s_36, \
                         hd_s_37, hd_s_38, hs_12, hp_21, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_1 * hs_s_12[k]
                  + f_2 * hd_s_36[k]
                  + f_3 * hs_12[k]
                  + pb_z[k] * hp_21[k];

        t_37[k] = pa_z[k] * gd_14[k]
                  + f_2 * hd_s_37[k];

        t_38[k] = f_4 * gp_9[k]
                  + f_2 * hd_s_38[k]
                  + pb_y[k] * hp_22[k];
    }

#pragma omp simd aligned(t_39, t_40, pa_y, pb_x, fd_s_9, fd_9, gd_16, hs_s_14, hd_s_39, \
                         hd_s_40, hs_14, hp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_5 * fd_s_9[k]
                  + f_6 * fd_9[k]
                  + pa_y[k] * gd_16[k]
                  + f_2 * hd_s_39[k];

        t_40[k] = -f_1 * hs_s_14[k]
                  + f_2 * hd_s_40[k]
                  + f_3 * hs_14[k]
                  + pb_x[k] * hp_23[k];
    }

#pragma omp simd aligned(t_41, t_42, pa_z, pb_y, fd_s_7, fd_7, gp_11, gd_15, hd_s_41, hd_s_42, \
                         hp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -f_7 * fd_s_7[k]
                  + f_3 * fd_7[k]
                  + pa_z[k] * gd_15[k]
                  + f_2 * hd_s_41[k];

        t_42[k] = f_6 * gp_11[k]
                  + f_2 * hd_s_42[k]
                  + pb_y[k] * hp_25[k];
    }

#pragma omp simd aligned(t_43, t_44, pa_y, pb_x, fd_s_11, fd_11, gd_18, hs_s_15, hd_s_43, \
                         hd_s_44, hs_15, hp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -f_8 * fd_s_11[k]
                  + f_9 * fd_11[k]
                  + pa_y[k] * gd_18[k]
                  + f_2 * hd_s_43[k];

        t_44[k] = -f_1 * hs_s_15[k]
                  + f_2 * hd_s_44[k]
                  + f_3 * hs_15[k]
                  + pb_x[k] * hp_26[k];
    }

#pragma omp simd aligned(t_45, t_46, pa_z, pb_y, fd_s_8, fd_8, gp_13, gd_17, hd_s_45, hd_s_46, \
                         hp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -f_8 * fd_s_8[k]
                  + f_9 * fd_8[k]
                  + pa_z[k] * gd_17[k]
                  + f_2 * hd_s_45[k];

        t_46[k] = f_9 * gp_13[k]
                  + f_2 * hd_s_46[k]
                  + pb_y[k] * hp_28[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pb_y, fd_s_13, fd_13, gp_15, gp_16, gd_20, \
                         gd_22, hd_s_47, hd_s_48, hd_s_49, hp_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -f_7 * fd_s_13[k]
                  + f_3 * fd_13[k]
                  + pa_y[k] * gd_20[k]
                  + f_2 * hd_s_47[k];

        t_48[k] = f_9 * gp_15[k]
                  + pa_y[k] * gd_22[k]
                  + f_2 * hd_s_48[k];

        t_49[k] = f_3 * gp_16[k]
                  + f_2 * hd_s_49[k]
                  + pb_y[k] * hp_30[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_y, pb_x, pb_y, gd_23, hs_s_17, hd_s_50, hd_s_51, \
                         hd_s_52, hs_17, hp_31, hp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * gd_23[k]
                  + f_2 * hd_s_50[k];

        t_51[k] = -f_1 * hs_s_17[k]
                  + f_2 * hd_s_51[k]
                  + f_3 * hs_17[k]
                  + pb_x[k] * hp_31[k];

        t_52[k] = -f_1 * hs_s_17[k]
                  + f_2 * hd_s_52[k]
                  + f_3 * hs_17[k]
                  + pb_y[k] * hp_32[k];
    }

#pragma omp simd aligned(t_53, pb_z, gp_16, hs_s_17, hd_s_53, hs_17, \
                         hp_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * gp_16[k]
                  - f_1 * hs_s_17[k]
                  + f_2 * hd_s_53[k]
                  + f_3 * hs_17[k]
                  + pb_z[k] * hp_33[k];
    }
}

auto
compute_prim_hd_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 * beta / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_1 = buffer.data(fd_s + 1);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_8 = buffer.data(fd_s + 8);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_10 = buffer.data(fd_s + 10);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_15 = buffer.data(fd_s + 15);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_15 = buffer.data(fd + 15);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);

    const auto *gd_0 = buffer.data(gd + 0);
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
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_7 = buffer.data(hs_s + 7);
    const auto *hs_s_11 = buffer.data(hs_s + 11);
    const auto *hs_s_13 = buffer.data(hs_s + 13);
    const auto *hs_s_14 = buffer.data(hs_s + 14);
    const auto *hs_s_16 = buffer.data(hs_s + 16);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_4 = buffer.data(hd_s + 4);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_21 = buffer.data(hd_s + 21);
    const auto *hd_s_22 = buffer.data(hd_s + 22);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_24 = buffer.data(hd_s + 24);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_28 = buffer.data(hd_s + 28);
    const auto *hd_s_29 = buffer.data(hd_s + 29);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_38 = buffer.data(hd_s + 38);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_42 = buffer.data(hd_s + 42);
    const auto *hd_s_43 = buffer.data(hd_s + 43);
    const auto *hd_s_44 = buffer.data(hd_s + 44);
    const auto *hd_s_45 = buffer.data(hd_s + 45);
    const auto *hd_s_46 = buffer.data(hd_s + 46);
    const auto *hd_s_47 = buffer.data(hd_s + 47);
    const auto *hd_s_48 = buffer.data(hd_s + 48);
    const auto *hd_s_49 = buffer.data(hd_s + 49);
    const auto *hd_s_50 = buffer.data(hd_s + 50);
    const auto *hd_s_51 = buffer.data(hd_s + 51);
    const auto *hd_s_52 = buffer.data(hd_s + 52);
    const auto *hd_s_53 = buffer.data(hd_s + 53);
    const auto *hd_s_54 = buffer.data(hd_s + 54);
    const auto *hd_s_55 = buffer.data(hd_s + 55);
    const auto *hd_s_56 = buffer.data(hd_s + 56);
    const auto *hd_s_57 = buffer.data(hd_s + 57);
    const auto *hd_s_58 = buffer.data(hd_s + 58);
    const auto *hd_s_59 = buffer.data(hd_s + 59);
    const auto *hd_s_60 = buffer.data(hd_s + 60);
    const auto *hd_s_61 = buffer.data(hd_s + 61);
    const auto *hd_s_62 = buffer.data(hd_s + 62);
    const auto *hd_s_63 = buffer.data(hd_s + 63);
    const auto *hd_s_64 = buffer.data(hd_s + 64);
    const auto *hd_s_65 = buffer.data(hd_s + 65);
    const auto *hd_s_66 = buffer.data(hd_s + 66);
    const auto *hd_s_67 = buffer.data(hd_s + 67);
    const auto *hd_s_68 = buffer.data(hd_s + 68);
    const auto *hd_s_69 = buffer.data(hd_s + 69);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_16 = buffer.data(hs + 16);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hs_0, hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_1[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_2[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, fd_s_2, fd_2, gd_0, gd_2, gd_4, \
                         hd_s_3, hd_s_4, hd_s_5, hd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * gd_0[k]
                 + f_2 * hd_s_3[k];

        t_4[k] = -f_4 * fd_s_2[k]
                 + f_5 * fd_2[k]
                 + pa_x[k] * gd_4[k]
                 + f_2 * hd_s_4[k];

        t_5[k] = pa_y[k] * gd_2[k]
                 + f_2 * hd_s_5[k];

        t_6[k] = pa_z[k] * gd_0[k]
                 + f_2 * hd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_x, pa_y, fd_s_0, fd_s_4, fd_0, fd_4, gd_3, gd_6, hd_s_7, \
                         hd_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_4 * fd_s_4[k]
                 + f_5 * fd_4[k]
                 + pa_x[k] * gd_6[k]
                 + f_2 * hd_s_7[k];

        t_8[k] = -f_6 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_y[k] * gd_3[k]
                 + f_2 * hd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pa_z, pb_z, fd_s_5, fd_5, gd_4, gd_8, hs_s_3, \
                         hd_s_9, hd_s_10, hd_s_11, hs_3, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_7 * fd_s_5[k]
                 + f_8 * fd_5[k]
                 + pa_x[k] * gd_8[k]
                 + f_2 * hd_s_9[k];

        t_10[k] = -f_1 * hs_s_3[k]
                  + f_2 * hd_s_10[k]
                  + f_3 * hs_3[k]
                  + pb_z[k] * hp_6[k];

        t_11[k] = pa_z[k] * gd_4[k]
                  + f_2 * hd_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pa_z, pb_y, fd_s_0, fd_0, gp_3, gd_5, gd_6, \
                         hd_s_12, hd_s_13, hd_s_14, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * gp_3[k]
                  + f_2 * hd_s_12[k]
                  + pb_y[k] * hp_7[k];

        t_13[k] = pa_y[k] * gd_6[k]
                  + f_2 * hd_s_13[k];

        t_14[k] = -f_6 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_z[k] * gd_5[k]
                  + f_2 * hd_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_x, pb_y, fd_s_6, fd_6, gd_10, hs_s_4, hd_s_15, \
                         hd_s_16, hs_4, hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_1 * hs_s_4[k]
                  + f_2 * hd_s_15[k]
                  + f_3 * hs_4[k]
                  + pb_y[k] * hp_8[k];

        t_16[k] = -f_7 * fd_s_6[k]
                  + f_8 * fd_6[k]
                  + pa_x[k] * gd_10[k]
                  + f_2 * hd_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pa_y, fd_s_1, fd_s_8, fd_1, fd_8, gd_7, gd_12, \
                         hd_s_17, hd_s_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_7 * fd_s_1[k]
                  + f_8 * fd_1[k]
                  + pa_y[k] * gd_7[k]
                  + f_2 * hd_s_17[k];

        t_18[k] = -f_6 * fd_s_8[k]
                  + f_3 * fd_8[k]
                  + pa_x[k] * gd_12[k]
                  + f_2 * hd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_z, pb_z, gd_7, gd_8, hs_s_5, hd_s_19, hd_s_20, \
                         hd_s_21, hs_5, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_1 * hs_s_5[k]
                  + f_2 * hd_s_19[k]
                  + f_3 * hs_5[k]
                  + pb_z[k] * hp_11[k];

        t_20[k] = pa_z[k] * gd_7[k]
                  + f_2 * hd_s_20[k];

        t_21[k] = pa_z[k] * gd_8[k]
                  + f_2 * hd_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, pb_y, fd_s_10, fd_10, gp_4, gd_9, \
                         gd_13, hd_s_22, hd_s_23, hd_s_24, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_8 * gp_4[k]
                  + f_2 * hd_s_22[k]
                  + pb_y[k] * hp_12[k];

        t_23[k] = -f_6 * fd_s_10[k]
                  + f_3 * fd_10[k]
                  + pa_x[k] * gd_13[k]
                  + f_2 * hd_s_23[k];

        t_24[k] = pa_y[k] * gd_9[k]
                  + f_2 * hd_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_y, pb_y, fd_s_11, fd_11, gp_5, gd_10, \
                         gd_14, hd_s_25, hd_s_26, hd_s_27, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_6 * fd_s_11[k]
                  + f_3 * fd_11[k]
                  + pa_x[k] * gd_14[k]
                  + f_2 * hd_s_25[k];

        t_26[k] = f_3 * gp_5[k]
                  + f_2 * hd_s_26[k]
                  + pb_y[k] * hp_13[k];

        t_27[k] = pa_y[k] * gd_10[k]
                  + f_2 * hd_s_27[k];
    }

#pragma omp simd aligned(t_28, t_29, pa_z, pb_y, fd_s_3, fd_3, gd_9, hs_s_7, hd_s_28, hd_s_29, \
                         hs_7, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_7 * fd_s_3[k]
                  + f_8 * fd_3[k]
                  + pa_z[k] * gd_9[k]
                  + f_2 * hd_s_28[k];

        t_29[k] = -f_1 * hs_s_7[k]
                  + f_2 * hd_s_29[k]
                  + f_3 * hs_7[k]
                  + pb_y[k] * hp_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, fd_s_15, fd_15, gp_6, gd_16, gd_17, \
                         gd_18, gd_19, hd_s_30, hd_s_31, hd_s_32, \
                         hd_s_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_6 * fd_s_15[k]
                  + f_3 * fd_15[k]
                  + pa_x[k] * gd_16[k]
                  + f_2 * hd_s_30[k];

        t_31[k] = f_8 * gp_6[k]
                  + pa_x[k] * gd_17[k]
                  + f_2 * hd_s_31[k];

        t_32[k] = pa_x[k] * gd_18[k]
                  + f_2 * hd_s_32[k];

        t_33[k] = pa_x[k] * gd_19[k]
                  + f_2 * hd_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pa_z, gp_10, gd_11, gd_21, gd_22, \
                         gd_23, hd_s_34, hd_s_35, hd_s_36, hd_s_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * gd_11[k]
                  + f_2 * hd_s_34[k];

        t_35[k] = pa_x[k] * gd_21[k]
                  + f_2 * hd_s_35[k];

        t_36[k] = pa_x[k] * gd_22[k]
                  + f_2 * hd_s_36[k];

        t_37[k] = f_8 * gp_10[k]
                  + pa_x[k] * gd_23[k]
                  + f_2 * hd_s_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pa_y, gd_15, gd_24, gd_25, gd_26, \
                         hd_s_38, hd_s_39, hd_s_40, hd_s_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_x[k] * gd_24[k]
                  + f_2 * hd_s_38[k];

        t_39[k] = pa_x[k] * gd_25[k]
                  + f_2 * hd_s_39[k];

        t_40[k] = pa_x[k] * gd_26[k]
                  + f_2 * hd_s_40[k];

        t_41[k] = pa_y[k] * gd_15[k]
                  + f_2 * hd_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, gp_13, gd_27, gd_28, gd_30, gd_31, \
                         hd_s_42, hd_s_43, hd_s_44, hd_s_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_x[k] * gd_27[k]
                  + f_2 * hd_s_42[k];

        t_43[k] = pa_x[k] * gd_28[k]
                  + f_2 * hd_s_43[k];

        t_44[k] = f_8 * gp_13[k]
                  + pa_x[k] * gd_30[k]
                  + f_2 * hd_s_44[k];

        t_45[k] = pa_x[k] * gd_31[k]
                  + f_2 * hd_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_x, pb_x, gd_32, hs_s_11, hd_s_46, hd_s_47, \
                         hd_s_48, hs_11, hp_18, hp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_x[k] * gd_32[k]
                  + f_2 * hd_s_46[k];

        t_47[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_47[k]
                  + f_3 * hs_11[k]
                  + pb_x[k] * hp_18[k];

        t_48[k] = f_2 * hd_s_48[k]
                  + pb_x[k] * hp_19[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_z, pb_y, pb_z, gp_7, gd_18, hs_s_11, hd_s_49, \
                         hd_s_50, hd_s_51, hs_11, hp_19, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_0 * gp_7[k]
                  - f_1 * hs_s_11[k]
                  + f_2 * hd_s_49[k]
                  + f_3 * hs_11[k]
                  + pb_y[k] * hp_19[k];

        t_50[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_50[k]
                  + f_3 * hs_11[k]
                  + pb_z[k] * hp_20[k];

        t_51[k] = pa_z[k] * gd_18[k]
                  + f_2 * hd_s_51[k];
    }

#pragma omp simd aligned(t_52, t_53, pa_y, pb_y, fd_s_10, fd_10, gp_9, gd_22, hd_s_52, \
                         hd_s_53, hp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_9 * gp_9[k]
                  + f_2 * hd_s_52[k]
                  + pb_y[k] * hp_21[k];

        t_53[k] = -f_4 * fd_s_10[k]
                  + f_5 * fd_10[k]
                  + pa_y[k] * gd_22[k]
                  + f_2 * hd_s_53[k];
    }

#pragma omp simd aligned(t_54, t_55, pa_z, pb_x, fd_s_8, fd_8, gd_20, hs_s_13, hd_s_54, \
                         hd_s_55, hs_13, hp_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = -f_1 * hs_s_13[k]
                  + f_2 * hd_s_54[k]
                  + f_3 * hs_13[k]
                  + pb_x[k] * hp_22[k];

        t_55[k] = -f_6 * fd_s_8[k]
                  + f_3 * fd_8[k]
                  + pa_z[k] * gd_20[k]
                  + f_2 * hd_s_55[k];
    }

#pragma omp simd aligned(t_56, t_57, pa_y, pb_y, fd_s_12, fd_12, gp_11, gd_26, hd_s_56, \
                         hd_s_57, hp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_5 * gp_11[k]
                  + f_2 * hd_s_56[k]
                  + pb_y[k] * hp_23[k];

        t_57[k] = -f_7 * fd_s_12[k]
                  + f_8 * fd_12[k]
                  + pa_y[k] * gd_26[k]
                  + f_2 * hd_s_57[k];
    }

#pragma omp simd aligned(t_58, t_59, pa_z, pb_x, fd_s_9, fd_9, gd_24, hs_s_14, hd_s_58, \
                         hd_s_59, hs_14, hp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -f_1 * hs_s_14[k]
                  + f_2 * hd_s_58[k]
                  + f_3 * hs_14[k]
                  + pb_x[k] * hp_24[k];

        t_59[k] = -f_7 * fd_s_9[k]
                  + f_8 * fd_9[k]
                  + pa_z[k] * gd_24[k]
                  + f_2 * hd_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_y, pb_y, fd_s_15, fd_15, gp_12, gp_14, gd_29, \
                         gd_31, hd_s_60, hd_s_61, hd_s_62, hp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_8 * gp_12[k]
                  + f_2 * hd_s_60[k]
                  + pb_y[k] * hp_25[k];

        t_61[k] = -f_6 * fd_s_15[k]
                  + f_3 * fd_15[k]
                  + pa_y[k] * gd_29[k]
                  + f_2 * hd_s_61[k];

        t_62[k] = f_8 * gp_14[k]
                  + pa_y[k] * gd_31[k]
                  + f_2 * hd_s_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pb_x, pb_y, gp_15, gd_32, hs_s_16, hd_s_63, \
                         hd_s_64, hd_s_65, hs_16, hp_26, hp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_3 * gp_15[k]
                  + f_2 * hd_s_63[k]
                  + pb_y[k] * hp_26[k];

        t_64[k] = pa_y[k] * gd_32[k]
                  + f_2 * hd_s_64[k];

        t_65[k] = -f_1 * hs_s_16[k]
                  + f_2 * hd_s_65[k]
                  + f_3 * hs_16[k]
                  + pb_x[k] * hp_27[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_x, pb_y, pb_z, gp_15, hs_s_16, hd_s_66, \
                         hd_s_67, hd_s_68, hd_s_69, hs_16, hp_28, \
                         hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_2 * hd_s_66[k]
                  + pb_x[k] * hp_29[k];

        t_67[k] = -f_1 * hs_s_16[k]
                  + f_2 * hd_s_67[k]
                  + f_3 * hs_16[k]
                  + pb_y[k] * hp_28[k];

        t_68[k] = f_2 * hd_s_68[k]
                  + pb_y[k] * hp_29[k];

        t_69[k] = f_0 * gp_15[k]
                  - f_1 * hs_s_16[k]
                  + f_2 * hd_s_69[k]
                  + f_3 * hs_16[k]
                  + pb_z[k] * hp_29[k];
    }
}

auto
compute_prim_hd_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 * beta / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_14 = buffer.data(fd_s + 14);
    const auto *fd_s_15 = buffer.data(fd_s + 15);
    const auto *fd_s_16 = buffer.data(fd_s + 16);
    const auto *fd_s_19 = buffer.data(fd_s + 19);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_19 = buffer.data(fd + 19);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_35 = buffer.data(gd + 35);

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_2 = buffer.data(hs_s + 2);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_8 = buffer.data(hs_s + 8);
    const auto *hs_s_10 = buffer.data(hs_s + 10);
    const auto *hs_s_11 = buffer.data(hs_s + 11);
    const auto *hs_s_13 = buffer.data(hs_s + 13);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_4 = buffer.data(hd_s + 4);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_21 = buffer.data(hd_s + 21);
    const auto *hd_s_22 = buffer.data(hd_s + 22);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_24 = buffer.data(hd_s + 24);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_28 = buffer.data(hd_s + 28);
    const auto *hd_s_29 = buffer.data(hd_s + 29);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_38 = buffer.data(hd_s + 38);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_42 = buffer.data(hd_s + 42);
    const auto *hd_s_43 = buffer.data(hd_s + 43);
    const auto *hd_s_44 = buffer.data(hd_s + 44);
    const auto *hd_s_45 = buffer.data(hd_s + 45);
    const auto *hd_s_46 = buffer.data(hd_s + 46);
    const auto *hd_s_47 = buffer.data(hd_s + 47);
    const auto *hd_s_48 = buffer.data(hd_s + 48);
    const auto *hd_s_49 = buffer.data(hd_s + 49);
    const auto *hd_s_50 = buffer.data(hd_s + 50);
    const auto *hd_s_51 = buffer.data(hd_s + 51);
    const auto *hd_s_52 = buffer.data(hd_s + 52);
    const auto *hd_s_53 = buffer.data(hd_s + 53);
    const auto *hd_s_54 = buffer.data(hd_s + 54);
    const auto *hd_s_55 = buffer.data(hd_s + 55);
    const auto *hd_s_56 = buffer.data(hd_s + 56);
    const auto *hd_s_57 = buffer.data(hd_s + 57);
    const auto *hd_s_58 = buffer.data(hd_s + 58);
    const auto *hd_s_59 = buffer.data(hd_s + 59);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_13 = buffer.data(hs + 13);

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
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hs_0, hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_1[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_2[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, t_6, pa_x, pa_y, pa_z, fd_s_3, fd_3, gd_0, gd_2, gd_4, \
                         hd_s_3, hd_s_4, hd_s_5, hd_s_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * gd_0[k]
                 + f_2 * hd_s_3[k];

        t_4[k] = -f_4 * fd_s_3[k]
                 + f_5 * fd_3[k]
                 + pa_x[k] * gd_4[k]
                 + f_2 * hd_s_4[k];

        t_5[k] = pa_y[k] * gd_2[k]
                 + f_2 * hd_s_5[k];

        t_6[k] = pa_z[k] * gd_0[k]
                 + f_2 * hd_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pa_y, pb_y, fd_s_0, fd_s_5, fd_0, fd_5, gd_3, \
                         gd_7, hd_s_7, hd_s_8, hd_s_9, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * hd_s_7[k]
                 + pb_y[k] * hp_3[k];

        t_8[k] = -f_4 * fd_s_5[k]
                 + f_5 * fd_5[k]
                 + pa_x[k] * gd_7[k]
                 + f_2 * hd_s_8[k];

        t_9[k] = -f_6 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_y[k] * gd_3[k]
                 + f_2 * hd_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_z, pb_z, fd_s_7, fd_7, gd_4, gd_9, hs_s_2, \
                         hd_s_10, hd_s_11, hd_s_12, hs_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_7 * fd_s_7[k]
                  + f_8 * fd_7[k]
                  + pa_x[k] * gd_9[k]
                  + f_2 * hd_s_10[k];

        t_11[k] = -f_1 * hs_s_2[k]
                  + f_2 * hd_s_11[k]
                  + f_3 * hs_2[k]
                  + pb_z[k] * hp_4[k];

        t_12[k] = pa_z[k] * gd_4[k]
                  + f_2 * hd_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_y, fd_s_0, fd_0, gd_6, gd_7, hs_s_3, \
                         hd_s_13, hd_s_14, hd_s_15, hs_3, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * gd_7[k]
                  + f_2 * hd_s_13[k];

        t_14[k] = -f_6 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_z[k] * gd_6[k]
                  + f_2 * hd_s_14[k];

        t_15[k] = -f_1 * hs_s_3[k]
                  + f_2 * hd_s_15[k]
                  + f_3 * hs_3[k]
                  + pb_y[k] * hp_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_y, fd_s_2, fd_s_9, fd_2, fd_9, gd_8, \
                         gd_14, hd_s_16, hd_s_17, hd_s_18, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * hd_s_16[k]
                  + pb_y[k] * hp_6[k];

        t_17[k] = -f_7 * fd_s_9[k]
                  + f_8 * fd_9[k]
                  + pa_x[k] * gd_14[k]
                  + f_2 * hd_s_17[k];

        t_18[k] = -f_7 * fd_s_2[k]
                  + f_8 * fd_2[k]
                  + pa_y[k] * gd_8[k]
                  + f_2 * hd_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_z, pb_z, fd_s_11, fd_11, gd_8, gd_16, \
                         hs_s_4, hd_s_19, hd_s_20, hd_s_21, hs_4, \
                         hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_6 * fd_s_11[k]
                  + f_3 * fd_11[k]
                  + pa_x[k] * gd_16[k]
                  + f_2 * hd_s_19[k];

        t_20[k] = -f_1 * hs_s_4[k]
                  + f_2 * hd_s_20[k]
                  + f_3 * hs_4[k]
                  + pb_z[k] * hp_7[k];

        t_21[k] = pa_z[k] * gd_8[k]
                  + f_2 * hd_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_x, pa_y, pa_z, fd_s_14, fd_14, gd_9, gd_13, \
                         gd_18, hd_s_22, hd_s_23, hd_s_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_z[k] * gd_9[k]
                  + f_2 * hd_s_22[k];

        t_23[k] = -f_6 * fd_s_14[k]
                  + f_3 * fd_14[k]
                  + pa_x[k] * gd_18[k]
                  + f_2 * hd_s_23[k];

        t_24[k] = pa_y[k] * gd_13[k]
                  + f_2 * hd_s_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pa_x, pa_y, pa_z, fd_s_4, fd_s_15, fd_4, fd_15, \
                         gd_13, gd_14, gd_19, hd_s_25, hd_s_26, \
                         hd_s_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_6 * fd_s_15[k]
                  + f_3 * fd_15[k]
                  + pa_x[k] * gd_19[k]
                  + f_2 * hd_s_25[k];

        t_26[k] = pa_y[k] * gd_14[k]
                  + f_2 * hd_s_26[k];

        t_27[k] = -f_7 * fd_s_4[k]
                  + f_8 * fd_4[k]
                  + pa_z[k] * gd_13[k]
                  + f_2 * hd_s_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pb_y, fd_s_19, fd_19, gd_21, hs_s_5, hd_s_28, \
                         hd_s_29, hd_s_30, hs_5, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_1 * hs_s_5[k]
                  + f_2 * hd_s_28[k]
                  + f_3 * hs_5[k]
                  + pb_y[k] * hp_8[k];

        t_29[k] = f_2 * hd_s_29[k]
                  + pb_y[k] * hp_9[k];

        t_30[k] = -f_6 * fd_s_19[k]
                  + f_3 * fd_19[k]
                  + pa_x[k] * gd_21[k]
                  + f_2 * hd_s_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pa_z, gp_4, gd_15, gd_22, gd_23, gd_26, \
                         hd_s_31, hd_s_32, hd_s_33, hd_s_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_8 * gp_4[k]
                  + pa_x[k] * gd_22[k]
                  + f_2 * hd_s_31[k];

        t_32[k] = pa_x[k] * gd_23[k]
                  + f_2 * hd_s_32[k];

        t_33[k] = pa_z[k] * gd_15[k]
                  + f_2 * hd_s_33[k];

        t_34[k] = pa_x[k] * gd_26[k]
                  + f_2 * hd_s_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, gp_7, gd_27, gd_28, gd_29, gd_30, \
                         hd_s_35, hd_s_36, hd_s_37, hd_s_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * gp_7[k]
                  + pa_x[k] * gd_27[k]
                  + f_2 * hd_s_35[k];

        t_36[k] = pa_x[k] * gd_28[k]
                  + f_2 * hd_s_36[k];

        t_37[k] = pa_x[k] * gd_29[k]
                  + f_2 * hd_s_37[k];

        t_38[k] = pa_x[k] * gd_30[k]
                  + f_2 * hd_s_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_x, pb_x, gp_8, gd_32, gd_35, hs_s_8, hd_s_39, \
                         hd_s_40, hd_s_41, hs_8, hp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_8 * gp_8[k]
                  + pa_x[k] * gd_32[k]
                  + f_2 * hd_s_39[k];

        t_40[k] = pa_x[k] * gd_35[k]
                  + f_2 * hd_s_40[k];

        t_41[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_41[k]
                  + f_3 * hs_8[k]
                  + pb_x[k] * hp_10[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_y, pb_z, gp_5, hs_s_8, hd_s_42, hd_s_43, \
                         hd_s_44, hs_8, hp_11, hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_2 * hd_s_42[k]
                  + pb_x[k] * hp_11[k];

        t_43[k] = f_0 * gp_5[k]
                  - f_1 * hs_s_8[k]
                  + f_2 * hd_s_43[k]
                  + f_3 * hs_8[k]
                  + pb_y[k] * hp_11[k];

        t_44[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_44[k]
                  + f_3 * hs_8[k]
                  + pb_z[k] * hp_12[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pa_z, pb_x, fd_s_14, fd_14, gd_23, gd_26, \
                         hs_s_10, hd_s_45, hd_s_46, hd_s_47, hs_10, \
                         hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * gd_23[k]
                  + f_2 * hd_s_45[k];

        t_46[k] = -f_4 * fd_s_14[k]
                  + f_5 * fd_14[k]
                  + pa_y[k] * gd_26[k]
                  + f_2 * hd_s_46[k];

        t_47[k] = -f_1 * hs_s_10[k]
                  + f_2 * hd_s_47[k]
                  + f_3 * hs_10[k]
                  + pb_x[k] * hp_13[k];
    }

#pragma omp simd aligned(t_48, t_49, pa_y, pa_z, fd_s_11, fd_s_16, fd_11, fd_16, gd_25, gd_29, \
                         hd_s_48, hd_s_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -f_6 * fd_s_11[k]
                  + f_3 * fd_11[k]
                  + pa_z[k] * gd_25[k]
                  + f_2 * hd_s_48[k];

        t_49[k] = -f_7 * fd_s_16[k]
                  + f_8 * fd_16[k]
                  + pa_y[k] * gd_29[k]
                  + f_2 * hd_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, pa_z, pb_x, fd_s_13, fd_13, gd_28, hs_s_11, hd_s_50, \
                         hd_s_51, hs_11, hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_50[k]
                  + f_3 * hs_11[k]
                  + pb_x[k] * hp_15[k];

        t_51[k] = -f_7 * fd_s_13[k]
                  + f_8 * fd_13[k]
                  + pa_z[k] * gd_28[k]
                  + f_2 * hd_s_51[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_y, fd_s_19, fd_19, gp_9, gd_31, gd_33, gd_35, \
                         hd_s_52, hd_s_53, hd_s_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = -f_6 * fd_s_19[k]
                  + f_3 * fd_19[k]
                  + pa_y[k] * gd_31[k]
                  + f_2 * hd_s_52[k];

        t_53[k] = f_8 * gp_9[k]
                  + pa_y[k] * gd_33[k]
                  + f_2 * hd_s_53[k];

        t_54[k] = pa_y[k] * gd_35[k]
                  + f_2 * hd_s_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_x, pb_y, hs_s_13, hd_s_55, hd_s_56, \
                         hd_s_57, hd_s_58, hs_13, hp_18, hp_19, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -f_1 * hs_s_13[k]
                  + f_2 * hd_s_55[k]
                  + f_3 * hs_13[k]
                  + pb_x[k] * hp_18[k];

        t_56[k] = f_2 * hd_s_56[k]
                  + pb_x[k] * hp_20[k];

        t_57[k] = -f_1 * hs_s_13[k]
                  + f_2 * hd_s_57[k]
                  + f_3 * hs_13[k]
                  + pb_y[k] * hp_19[k];

        t_58[k] = f_2 * hd_s_58[k]
                  + pb_y[k] * hp_20[k];
    }

#pragma omp simd aligned(t_59, pb_z, gp_10, hs_s_13, hd_s_59, hs_13, \
                         hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * gp_10[k]
                  - f_1 * hs_s_13[k]
                  + f_2 * hd_s_59[k]
                  + f_3 * hs_13[k]
                  + pb_z[k] * hp_20[k];
    }
}

auto
compute_prim_hd_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 * beta / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_16 = buffer.data(fd_s + 16);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_16 = buffer.data(fd + 16);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);

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
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_27 = buffer.data(gd + 27);

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_2 = buffer.data(hs_s + 2);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_8 = buffer.data(hs_s + 8);
    const auto *hs_s_10 = buffer.data(hs_s + 10);
    const auto *hs_s_11 = buffer.data(hs_s + 11);
    const auto *hs_s_13 = buffer.data(hs_s + 13);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_4 = buffer.data(hd_s + 4);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_21 = buffer.data(hd_s + 21);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_28 = buffer.data(hd_s + 28);
    const auto *hd_s_29 = buffer.data(hd_s + 29);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_38 = buffer.data(hd_s + 38);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_42 = buffer.data(hd_s + 42);
    const auto *hd_s_43 = buffer.data(hd_s + 43);
    const auto *hd_s_44 = buffer.data(hd_s + 44);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_13 = buffer.data(hs + 13);

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
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hs_0, hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_1[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_2[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, fd_s_3, fd_3, gd_0, gd_4, hd_s_3, \
                         hd_s_4, hd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * gd_0[k]
                 + f_2 * hd_s_3[k];

        t_4[k] = -f_4 * fd_s_3[k]
                 + f_5 * fd_3[k]
                 + pa_x[k] * gd_4[k]
                 + f_2 * hd_s_4[k];

        t_5[k] = pa_z[k] * gd_0[k]
                 + f_2 * hd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_y, fd_s_0, fd_s_5, fd_0, fd_5, gd_3, \
                         gd_6, hd_s_6, hd_s_7, hd_s_8, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * hd_s_6[k]
                 + pb_y[k] * hp_3[k];

        t_7[k] = -f_4 * fd_s_5[k]
                 + f_5 * fd_5[k]
                 + pa_x[k] * gd_6[k]
                 + f_2 * hd_s_7[k];

        t_8[k] = -f_6 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_y[k] * gd_3[k]
                 + f_2 * hd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pb_z, fd_s_6, fd_6, gd_8, hs_s_2, hd_s_9, hd_s_10, \
                         hs_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_7 * fd_s_6[k]
                 + f_8 * fd_6[k]
                 + pa_x[k] * gd_8[k]
                 + f_2 * hd_s_9[k];

        t_10[k] = -f_1 * hs_s_2[k]
                  + f_2 * hd_s_10[k]
                  + f_3 * hs_2[k]
                  + pb_z[k] * hp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_z, pb_y, fd_s_0, fd_0, gd_5, hs_s_3, hd_s_11, \
                         hd_s_12, hd_s_13, hs_3, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -f_6 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_z[k] * gd_5[k]
                  + f_2 * hd_s_11[k];

        t_12[k] = -f_1 * hs_s_3[k]
                  + f_2 * hd_s_12[k]
                  + f_3 * hs_3[k]
                  + pb_y[k] * hp_5[k];

        t_13[k] = f_2 * hd_s_13[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, pa_x, pa_y, fd_s_2, fd_s_7, fd_2, fd_7, gd_7, gd_11, \
                         hd_s_14, hd_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -f_7 * fd_s_7[k]
                  + f_8 * fd_7[k]
                  + pa_x[k] * gd_11[k]
                  + f_2 * hd_s_14[k];

        t_15[k] = -f_7 * fd_s_2[k]
                  + f_8 * fd_2[k]
                  + pa_y[k] * gd_7[k]
                  + f_2 * hd_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, pa_x, pb_z, fd_s_9, fd_9, gd_12, hs_s_4, hd_s_16, \
                         hd_s_17, hs_4, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_6 * fd_s_9[k]
                  + f_3 * fd_9[k]
                  + pa_x[k] * gd_12[k]
                  + f_2 * hd_s_16[k];

        t_17[k] = -f_1 * hs_s_4[k]
                  + f_2 * hd_s_17[k]
                  + f_3 * hs_4[k]
                  + pb_z[k] * hp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_z, pb_y, fd_s_4, fd_4, gd_10, hs_s_5, hd_s_18, \
                         hd_s_19, hd_s_20, hs_5, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_7 * fd_s_4[k]
                  + f_8 * fd_4[k]
                  + pa_z[k] * gd_10[k]
                  + f_2 * hd_s_18[k];

        t_19[k] = -f_1 * hs_s_5[k]
                  + f_2 * hd_s_19[k]
                  + f_3 * hs_5[k]
                  + pb_y[k] * hp_8[k];

        t_20[k] = f_2 * hd_s_20[k]
                  + pb_y[k] * hp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, fd_s_16, fd_16, gd_13, gd_15, gd_27, hd_s_21, \
                         hd_s_23, hd_s_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_6 * fd_s_16[k]
                  + f_3 * fd_16[k]
                  + pa_x[k] * gd_13[k]
                  + f_2 * hd_s_21[k];

        t_22[k] = pa_x[k] * gd_15[k]
                  + f_2 * hd_s_23[k];

        t_23[k] = pa_x[k] * gd_27[k]
                  + f_2 * hd_s_25[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pb_y, gp_5, hs_s_8, hd_s_26, hd_s_27, \
                         hd_s_28, hs_8, hp_10, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_26[k]
                  + f_3 * hs_8[k]
                  + pb_x[k] * hp_10[k];

        t_25[k] = f_2 * hd_s_27[k]
                  + pb_x[k] * hp_11[k];

        t_26[k] = f_0 * gp_5[k]
                  - f_1 * hs_s_8[k]
                  + f_2 * hd_s_28[k]
                  + f_3 * hs_8[k]
                  + pb_y[k] * hp_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pa_z, pb_z, fd_s_12, fd_12, gd_15, gd_18, \
                         hs_s_8, hd_s_29, hd_s_30, hd_s_31, hs_8, \
                         hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_29[k]
                  + f_3 * hs_8[k]
                  + pb_z[k] * hp_12[k];

        t_28[k] = pa_z[k] * gd_15[k]
                  + f_2 * hd_s_30[k];

        t_29[k] = -f_4 * fd_s_12[k]
                  + f_5 * fd_12[k]
                  + pa_y[k] * gd_18[k]
                  + f_2 * hd_s_31[k];
    }

#pragma omp simd aligned(t_30, t_31, pa_z, pb_x, fd_s_9, fd_9, gd_17, hs_s_10, hd_s_32, \
                         hd_s_33, hs_10, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_1 * hs_s_10[k]
                  + f_2 * hd_s_32[k]
                  + f_3 * hs_10[k]
                  + pb_x[k] * hp_13[k];

        t_31[k] = -f_6 * fd_s_9[k]
                  + f_3 * fd_9[k]
                  + pa_z[k] * gd_17[k]
                  + f_2 * hd_s_33[k];
    }

#pragma omp simd aligned(t_32, t_33, pa_y, pb_x, fd_s_13, fd_13, gd_21, hs_s_11, hd_s_34, \
                         hd_s_35, hs_11, hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -f_7 * fd_s_13[k]
                  + f_8 * fd_13[k]
                  + pa_y[k] * gd_21[k]
                  + f_2 * hd_s_34[k];

        t_33[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_35[k]
                  + f_3 * hs_11[k]
                  + pb_x[k] * hp_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, fd_s_11, fd_s_16, fd_11, fd_16, gp_9, \
                         gd_20, gd_23, gd_25, hd_s_36, hd_s_37, \
                         hd_s_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -f_7 * fd_s_11[k]
                  + f_8 * fd_11[k]
                  + pa_z[k] * gd_20[k]
                  + f_2 * hd_s_36[k];

        t_35[k] = -f_6 * fd_s_16[k]
                  + f_3 * fd_16[k]
                  + pa_y[k] * gd_23[k]
                  + f_2 * hd_s_37[k];

        t_36[k] = f_8 * gp_9[k]
                  + pa_y[k] * gd_25[k]
                  + f_2 * hd_s_38[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pb_x, gd_27, hs_s_13, hd_s_39, hd_s_40, \
                         hd_s_41, hs_13, hp_18, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_y[k] * gd_27[k]
                  + f_2 * hd_s_39[k];

        t_38[k] = -f_1 * hs_s_13[k]
                  + f_2 * hd_s_40[k]
                  + f_3 * hs_13[k]
                  + pb_x[k] * hp_18[k];

        t_39[k] = f_2 * hd_s_41[k]
                  + pb_x[k] * hp_20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_y, pb_z, gp_10, hs_s_13, hd_s_42, hd_s_43, \
                         hd_s_44, hs_13, hp_19, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -f_1 * hs_s_13[k]
                  + f_2 * hd_s_42[k]
                  + f_3 * hs_13[k]
                  + pb_y[k] * hp_19[k];

        t_41[k] = f_2 * hd_s_43[k]
                  + pb_y[k] * hp_20[k];

        t_42[k] = f_0 * gp_10[k]
                  - f_1 * hs_s_13[k]
                  + f_2 * hd_s_44[k]
                  + f_3 * hs_13[k]
                  + pb_z[k] * hp_20[k];
    }
}

auto
compute_prim_hd_kinetic_energy_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fd_s, const size_t fd,
                                 const size_t gp, const size_t gd, const size_t hs_s,
                                 const size_t hd_s, const size_t hs, const size_t hp,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 * beta / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = beta / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 1.0 / p;

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

    const auto *fd_s_0 = buffer.data(fd_s + 0);
    const auto *fd_s_2 = buffer.data(fd_s + 2);
    const auto *fd_s_3 = buffer.data(fd_s + 3);
    const auto *fd_s_4 = buffer.data(fd_s + 4);
    const auto *fd_s_5 = buffer.data(fd_s + 5);
    const auto *fd_s_6 = buffer.data(fd_s + 6);
    const auto *fd_s_7 = buffer.data(fd_s + 7);
    const auto *fd_s_9 = buffer.data(fd_s + 9);
    const auto *fd_s_11 = buffer.data(fd_s + 11);
    const auto *fd_s_12 = buffer.data(fd_s + 12);
    const auto *fd_s_13 = buffer.data(fd_s + 13);
    const auto *fd_s_16 = buffer.data(fd_s + 16);

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_16 = buffer.data(fd + 16);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_10 = buffer.data(gp + 10);

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

    const auto *hs_s_0 = buffer.data(hs_s + 0);
    const auto *hs_s_2 = buffer.data(hs_s + 2);
    const auto *hs_s_3 = buffer.data(hs_s + 3);
    const auto *hs_s_4 = buffer.data(hs_s + 4);
    const auto *hs_s_5 = buffer.data(hs_s + 5);
    const auto *hs_s_8 = buffer.data(hs_s + 8);
    const auto *hs_s_10 = buffer.data(hs_s + 10);
    const auto *hs_s_11 = buffer.data(hs_s + 11);
    const auto *hs_s_13 = buffer.data(hs_s + 13);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_1 = buffer.data(hd_s + 1);
    const auto *hd_s_2 = buffer.data(hd_s + 2);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_4 = buffer.data(hd_s + 4);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_6 = buffer.data(hd_s + 6);
    const auto *hd_s_7 = buffer.data(hd_s + 7);
    const auto *hd_s_8 = buffer.data(hd_s + 8);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_10 = buffer.data(hd_s + 10);
    const auto *hd_s_11 = buffer.data(hd_s + 11);
    const auto *hd_s_12 = buffer.data(hd_s + 12);
    const auto *hd_s_13 = buffer.data(hd_s + 13);
    const auto *hd_s_14 = buffer.data(hd_s + 14);
    const auto *hd_s_15 = buffer.data(hd_s + 15);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_19 = buffer.data(hd_s + 19);
    const auto *hd_s_20 = buffer.data(hd_s + 20);
    const auto *hd_s_21 = buffer.data(hd_s + 21);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_25 = buffer.data(hd_s + 25);
    const auto *hd_s_26 = buffer.data(hd_s + 26);
    const auto *hd_s_27 = buffer.data(hd_s + 27);
    const auto *hd_s_28 = buffer.data(hd_s + 28);
    const auto *hd_s_29 = buffer.data(hd_s + 29);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_31 = buffer.data(hd_s + 31);
    const auto *hd_s_32 = buffer.data(hd_s + 32);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_37 = buffer.data(hd_s + 37);
    const auto *hd_s_39 = buffer.data(hd_s + 39);
    const auto *hd_s_40 = buffer.data(hd_s + 40);
    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_42 = buffer.data(hd_s + 42);
    const auto *hd_s_43 = buffer.data(hd_s + 43);
    const auto *hd_s_44 = buffer.data(hd_s + 44);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_13 = buffer.data(hs + 13);

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
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gp_0, hs_s_0, hd_s_0, hd_s_1, \
                         hd_s_2, hs_0, hp_0, hp_1, hp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gp_0[k]
                 - f_1 * hs_s_0[k]
                 + f_2 * hd_s_0[k]
                 + f_3 * hs_0[k]
                 + pb_x[k] * hp_0[k];

        t_1[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_1[k]
                 + f_3 * hs_0[k]
                 + pb_y[k] * hp_1[k];

        t_2[k] = -f_1 * hs_s_0[k]
                 + f_2 * hd_s_2[k]
                 + f_3 * hs_0[k]
                 + pb_z[k] * hp_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pa_y, pa_z, fd_s_3, fd_3, gd_0, gd_4, hd_s_3, \
                         hd_s_4, hd_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_y[k] * gd_0[k]
                 + f_2 * hd_s_3[k];

        t_4[k] = -f_4 * fd_s_3[k]
                 + f_5 * fd_3[k]
                 + pa_x[k] * gd_4[k]
                 + f_2 * hd_s_4[k];

        t_5[k] = pa_z[k] * gd_0[k]
                 + f_2 * hd_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, pb_y, fd_s_0, fd_s_5, fd_0, fd_5, gd_3, \
                         gd_6, hd_s_6, hd_s_7, hd_s_8, hp_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * hd_s_6[k]
                 + pb_y[k] * hp_3[k];

        t_7[k] = -f_4 * fd_s_5[k]
                 + f_5 * fd_5[k]
                 + pa_x[k] * gd_6[k]
                 + f_2 * hd_s_7[k];

        t_8[k] = -f_6 * fd_s_0[k]
                 + f_3 * fd_0[k]
                 + pa_y[k] * gd_3[k]
                 + f_2 * hd_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pb_z, fd_s_6, fd_6, gd_8, hs_s_2, hd_s_9, hd_s_10, \
                         hs_2, hp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_7 * fd_s_6[k]
                 + f_8 * fd_6[k]
                 + pa_x[k] * gd_8[k]
                 + f_2 * hd_s_9[k];

        t_10[k] = -f_1 * hs_s_2[k]
                  + f_2 * hd_s_10[k]
                  + f_3 * hs_2[k]
                  + pb_z[k] * hp_4[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pa_z, pb_y, fd_s_0, fd_0, gd_5, hs_s_3, hd_s_11, \
                         hd_s_12, hd_s_13, hs_3, hp_5, hp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = -f_6 * fd_s_0[k]
                  + f_3 * fd_0[k]
                  + pa_z[k] * gd_5[k]
                  + f_2 * hd_s_11[k];

        t_12[k] = -f_1 * hs_s_3[k]
                  + f_2 * hd_s_12[k]
                  + f_3 * hs_3[k]
                  + pb_y[k] * hp_5[k];

        t_13[k] = f_2 * hd_s_13[k]
                  + pb_y[k] * hp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, pa_x, pa_y, fd_s_2, fd_s_7, fd_2, fd_7, gd_7, gd_11, \
                         hd_s_14, hd_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -f_7 * fd_s_7[k]
                  + f_8 * fd_7[k]
                  + pa_x[k] * gd_11[k]
                  + f_2 * hd_s_14[k];

        t_15[k] = -f_7 * fd_s_2[k]
                  + f_8 * fd_2[k]
                  + pa_y[k] * gd_7[k]
                  + f_2 * hd_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, pa_x, pb_z, fd_s_9, fd_9, gd_12, hs_s_4, hd_s_16, \
                         hd_s_17, hs_4, hp_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_6 * fd_s_9[k]
                  + f_3 * fd_9[k]
                  + pa_x[k] * gd_12[k]
                  + f_2 * hd_s_16[k];

        t_17[k] = -f_1 * hs_s_4[k]
                  + f_2 * hd_s_17[k]
                  + f_3 * hs_4[k]
                  + pb_z[k] * hp_7[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_z, pb_y, fd_s_4, fd_4, gd_10, hs_s_5, hd_s_18, \
                         hd_s_19, hd_s_20, hs_5, hp_8, hp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_7 * fd_s_4[k]
                  + f_8 * fd_4[k]
                  + pa_z[k] * gd_10[k]
                  + f_2 * hd_s_18[k];

        t_19[k] = -f_1 * hs_s_5[k]
                  + f_2 * hd_s_19[k]
                  + f_3 * hs_5[k]
                  + pb_y[k] * hp_8[k];

        t_20[k] = f_2 * hd_s_20[k]
                  + pb_y[k] * hp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, fd_s_16, fd_16, gd_13, gd_15, gd_26, hd_s_21, \
                         hd_s_23, hd_s_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_6 * fd_s_16[k]
                  + f_3 * fd_16[k]
                  + pa_x[k] * gd_13[k]
                  + f_2 * hd_s_21[k];

        t_22[k] = pa_x[k] * gd_15[k]
                  + f_2 * hd_s_23[k];

        t_23[k] = pa_x[k] * gd_26[k]
                  + f_2 * hd_s_25[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pb_y, gp_5, hs_s_8, hd_s_26, hd_s_27, \
                         hd_s_28, hs_8, hp_10, hp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_26[k]
                  + f_3 * hs_8[k]
                  + pb_x[k] * hp_10[k];

        t_25[k] = f_2 * hd_s_27[k]
                  + pb_x[k] * hp_11[k];

        t_26[k] = f_0 * gp_5[k]
                  - f_1 * hs_s_8[k]
                  + f_2 * hd_s_28[k]
                  + f_3 * hs_8[k]
                  + pb_y[k] * hp_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pa_z, pb_z, fd_s_12, fd_12, gd_15, gd_18, \
                         hs_s_8, hd_s_29, hd_s_30, hd_s_31, hs_8, \
                         hp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * hs_s_8[k]
                  + f_2 * hd_s_29[k]
                  + f_3 * hs_8[k]
                  + pb_z[k] * hp_12[k];

        t_28[k] = pa_z[k] * gd_15[k]
                  + f_2 * hd_s_30[k];

        t_29[k] = -f_4 * fd_s_12[k]
                  + f_5 * fd_12[k]
                  + pa_y[k] * gd_18[k]
                  + f_2 * hd_s_31[k];
    }

#pragma omp simd aligned(t_30, t_31, pa_z, pb_x, fd_s_9, fd_9, gd_17, hs_s_10, hd_s_32, \
                         hd_s_33, hs_10, hp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_1 * hs_s_10[k]
                  + f_2 * hd_s_32[k]
                  + f_3 * hs_10[k]
                  + pb_x[k] * hp_13[k];

        t_31[k] = -f_6 * fd_s_9[k]
                  + f_3 * fd_9[k]
                  + pa_z[k] * gd_17[k]
                  + f_2 * hd_s_33[k];
    }

#pragma omp simd aligned(t_32, t_33, pa_y, pb_x, fd_s_13, fd_13, gd_21, hs_s_11, hd_s_34, \
                         hd_s_35, hs_11, hp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = -f_7 * fd_s_13[k]
                  + f_8 * fd_13[k]
                  + pa_y[k] * gd_21[k]
                  + f_2 * hd_s_34[k];

        t_33[k] = -f_1 * hs_s_11[k]
                  + f_2 * hd_s_35[k]
                  + f_3 * hs_11[k]
                  + pb_x[k] * hp_15[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, fd_s_11, fd_s_16, fd_11, fd_16, gd_20, \
                         gd_22, gd_26, hd_s_36, hd_s_37, hd_s_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -f_7 * fd_s_11[k]
                  + f_8 * fd_11[k]
                  + pa_z[k] * gd_20[k]
                  + f_2 * hd_s_36[k];

        t_35[k] = -f_6 * fd_s_16[k]
                  + f_3 * fd_16[k]
                  + pa_y[k] * gd_22[k]
                  + f_2 * hd_s_37[k];

        t_36[k] = pa_y[k] * gd_26[k]
                  + f_2 * hd_s_39[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_x, pb_y, hs_s_13, hd_s_40, hd_s_41, \
                         hd_s_42, hd_s_43, hs_13, hp_18, hp_19, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_1 * hs_s_13[k]
                  + f_2 * hd_s_40[k]
                  + f_3 * hs_13[k]
                  + pb_x[k] * hp_18[k];

        t_38[k] = f_2 * hd_s_41[k]
                  + pb_x[k] * hp_20[k];

        t_39[k] = -f_1 * hs_s_13[k]
                  + f_2 * hd_s_42[k]
                  + f_3 * hs_13[k]
                  + pb_y[k] * hp_19[k];

        t_40[k] = f_2 * hd_s_43[k]
                  + pb_y[k] * hp_20[k];
    }

#pragma omp simd aligned(t_41, pb_z, gp_10, hs_s_13, hd_s_44, hs_13, \
                         hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_0 * gp_10[k]
                  - f_1 * hs_s_13[k]
                  + f_2 * hd_s_44[k]
                  + f_3 * hs_13[k]
                  + pb_z[k] * hp_20[k];
    }
}

}  // namespace simdkin
