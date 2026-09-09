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


#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_0 = buffer.data(msh0 + 0);
    const auto *msh0_3 = buffer.data(msh0 + 3);
    const auto *msh0_5 = buffer.data(msh0 + 5);
    const auto *msh0_6 = buffer.data(msh0 + 6);
    const auto *msh0_9 = buffer.data(msh0 + 9);
    const auto *msh0_15 = buffer.data(msh0 + 15);
    const auto *msh0_20 = buffer.data(msh0 + 20);
    const auto *msh0_24 = buffer.data(msh0 + 24);
    const auto *msh0_27 = buffer.data(msh0 + 27);
    const auto *msh0_36 = buffer.data(msh0 + 36);
    const auto *msh0_42 = buffer.data(msh0 + 42);
    const auto *msh0_47 = buffer.data(msh0 + 47);
    const auto *msh0_51 = buffer.data(msh0 + 51);
    const auto *msh0_62 = buffer.data(msh0 + 62);

    const auto *msg_0 = buffer.data(msg + 0);
    const auto *msg_1 = buffer.data(msg + 1);
    const auto *msg_2 = buffer.data(msg + 2);
    const auto *msg_3 = buffer.data(msg + 3);
    const auto *msg_5 = buffer.data(msg + 5);
    const auto *msg_10 = buffer.data(msg + 10);
    const auto *msg_12 = buffer.data(msg + 12);
    const auto *msg_14 = buffer.data(msg + 14);
    const auto *msg_15 = buffer.data(msg + 15);
    const auto *msg_18 = buffer.data(msg + 18);
    const auto *msg_20 = buffer.data(msg + 20);
    const auto *msg_25 = buffer.data(msg + 25);
    const auto *msg_27 = buffer.data(msg + 27);
    const auto *msg_28 = buffer.data(msg + 28);
    const auto *msg_29 = buffer.data(msg + 29);
    const auto *msg_30 = buffer.data(msg + 30);
    const auto *msg_32 = buffer.data(msg + 32);
    const auto *msg_35 = buffer.data(msg + 35);
    const auto *msg_40 = buffer.data(msg + 40);
    const auto *msg_41 = buffer.data(msg + 41);
    const auto *msg_42 = buffer.data(msg + 42);
    const auto *msg_43 = buffer.data(msg + 43);
    const auto *msg_44 = buffer.data(msg + 44);
    const auto *msg_45 = buffer.data(msg + 45);
    const auto *msg_48 = buffer.data(msg + 48);
    const auto *msg_51 = buffer.data(msg + 51);
    const auto *msg_55 = buffer.data(msg + 55);
    const auto *msg_57 = buffer.data(msg + 57);
    const auto *msg_58 = buffer.data(msg + 58);
    const auto *msg_59 = buffer.data(msg + 59);
    const auto *msg_70 = buffer.data(msg + 70);
    const auto *msg_71 = buffer.data(msg + 71);
    const auto *msg_72 = buffer.data(msg + 72);
    const auto *msg_73 = buffer.data(msg + 73);
    const auto *msg_74 = buffer.data(msg + 74);
    const auto *msg_75 = buffer.data(msg + 75);
    const auto *msg_80 = buffer.data(msg + 80);
    const auto *msg_84 = buffer.data(msg + 84);
    const auto *msg_85 = buffer.data(msg + 85);
    const auto *msg_86 = buffer.data(msg + 86);
    const auto *msg_87 = buffer.data(msg + 87);
    const auto *msg_89 = buffer.data(msg + 89);
    const auto *msg_90 = buffer.data(msg + 90);
    const auto *msg_93 = buffer.data(msg + 93);

    const auto *msh1_0 = buffer.data(msh1 + 0);
    const auto *msh1_3 = buffer.data(msh1 + 3);
    const auto *msh1_5 = buffer.data(msh1 + 5);
    const auto *msh1_6 = buffer.data(msh1 + 6);
    const auto *msh1_9 = buffer.data(msh1 + 9);
    const auto *msh1_15 = buffer.data(msh1 + 15);
    const auto *msh1_20 = buffer.data(msh1 + 20);
    const auto *msh1_24 = buffer.data(msh1 + 24);
    const auto *msh1_27 = buffer.data(msh1 + 27);
    const auto *msh1_36 = buffer.data(msh1 + 36);
    const auto *msh1_42 = buffer.data(msh1 + 42);
    const auto *msh1_47 = buffer.data(msh1 + 47);
    const auto *msh1_51 = buffer.data(msh1 + 51);
    const auto *msh1_62 = buffer.data(msh1 + 62);

    const auto *nsf0_0 = buffer.data(nsf0 + 0);
    const auto *nsf0_1 = buffer.data(nsf0 + 1);
    const auto *nsf0_2 = buffer.data(nsf0 + 2);
    const auto *nsf0_6 = buffer.data(nsf0 + 6);
    const auto *nsf0_8 = buffer.data(nsf0 + 8);
    const auto *nsf0_9 = buffer.data(nsf0 + 9);
    const auto *nsf0_16 = buffer.data(nsf0 + 16);
    const auto *nsf0_17 = buffer.data(nsf0 + 17);
    const auto *nsf0_22 = buffer.data(nsf0 + 22);
    const auto *nsf0_27 = buffer.data(nsf0 + 27);
    const auto *nsf0_28 = buffer.data(nsf0 + 28);
    const auto *nsf0_29 = buffer.data(nsf0 + 29);
    const auto *nsf0_30 = buffer.data(nsf0 + 30);
    const auto *nsf0_32 = buffer.data(nsf0 + 32);
    const auto *nsf0_33 = buffer.data(nsf0 + 33);
    const auto *nsf0_36 = buffer.data(nsf0 + 36);
    const auto *nsf0_37 = buffer.data(nsf0 + 37);
    const auto *nsf0_39 = buffer.data(nsf0 + 39);
    const auto *nsf0_48 = buffer.data(nsf0 + 48);
    const auto *nsf0_49 = buffer.data(nsf0 + 49);
    const auto *nsf0_50 = buffer.data(nsf0 + 50);
    const auto *nsf0_51 = buffer.data(nsf0 + 51);
    const auto *nsf0_52 = buffer.data(nsf0 + 52);
    const auto *nsf0_55 = buffer.data(nsf0 + 55);
    const auto *nsf0_56 = buffer.data(nsf0 + 56);
    const auto *nsf0_57 = buffer.data(nsf0 + 57);
    const auto *nsf0_58 = buffer.data(nsf0 + 58);
    const auto *nsf0_59 = buffer.data(nsf0 + 59);
    const auto *nsf0_60 = buffer.data(nsf0 + 60);
    const auto *nsf0_63 = buffer.data(nsf0 + 63);

    const auto *nsf1_0 = buffer.data(nsf1 + 0);
    const auto *nsf1_1 = buffer.data(nsf1 + 1);
    const auto *nsf1_2 = buffer.data(nsf1 + 2);
    const auto *nsf1_6 = buffer.data(nsf1 + 6);
    const auto *nsf1_8 = buffer.data(nsf1 + 8);
    const auto *nsf1_9 = buffer.data(nsf1 + 9);
    const auto *nsf1_16 = buffer.data(nsf1 + 16);
    const auto *nsf1_17 = buffer.data(nsf1 + 17);
    const auto *nsf1_22 = buffer.data(nsf1 + 22);
    const auto *nsf1_27 = buffer.data(nsf1 + 27);
    const auto *nsf1_28 = buffer.data(nsf1 + 28);
    const auto *nsf1_29 = buffer.data(nsf1 + 29);
    const auto *nsf1_30 = buffer.data(nsf1 + 30);
    const auto *nsf1_32 = buffer.data(nsf1 + 32);
    const auto *nsf1_33 = buffer.data(nsf1 + 33);
    const auto *nsf1_36 = buffer.data(nsf1 + 36);
    const auto *nsf1_37 = buffer.data(nsf1 + 37);
    const auto *nsf1_39 = buffer.data(nsf1 + 39);
    const auto *nsf1_48 = buffer.data(nsf1 + 48);
    const auto *nsf1_49 = buffer.data(nsf1 + 49);
    const auto *nsf1_50 = buffer.data(nsf1 + 50);
    const auto *nsf1_51 = buffer.data(nsf1 + 51);
    const auto *nsf1_52 = buffer.data(nsf1 + 52);
    const auto *nsf1_55 = buffer.data(nsf1 + 55);
    const auto *nsf1_56 = buffer.data(nsf1 + 56);
    const auto *nsf1_57 = buffer.data(nsf1 + 57);
    const auto *nsf1_58 = buffer.data(nsf1 + 58);
    const auto *nsf1_59 = buffer.data(nsf1 + 59);
    const auto *nsf1_60 = buffer.data(nsf1 + 60);
    const auto *nsf1_63 = buffer.data(nsf1 + 63);

    const auto *nsg_0 = buffer.data(nsg + 0);
    const auto *nsg_1 = buffer.data(nsg + 1);
    const auto *nsg_2 = buffer.data(nsg + 2);
    const auto *nsg_3 = buffer.data(nsg + 3);
    const auto *nsg_5 = buffer.data(nsg + 5);
    const auto *nsg_6 = buffer.data(nsg + 6);
    const auto *nsg_9 = buffer.data(nsg + 9);
    const auto *nsg_10 = buffer.data(nsg + 10);
    const auto *nsg_12 = buffer.data(nsg + 12);
    const auto *nsg_13 = buffer.data(nsg + 13);
    const auto *nsg_14 = buffer.data(nsg + 14);
    const auto *nsg_15 = buffer.data(nsg + 15);
    const auto *nsg_16 = buffer.data(nsg + 16);
    const auto *nsg_18 = buffer.data(nsg + 18);
    const auto *nsg_20 = buffer.data(nsg + 20);
    const auto *nsg_21 = buffer.data(nsg + 21);
    const auto *nsg_25 = buffer.data(nsg + 25);
    const auto *nsg_26 = buffer.data(nsg + 26);
    const auto *nsg_27 = buffer.data(nsg + 27);
    const auto *nsg_28 = buffer.data(nsg + 28);
    const auto *nsg_29 = buffer.data(nsg + 29);
    const auto *nsg_30 = buffer.data(nsg + 30);
    const auto *nsg_32 = buffer.data(nsg + 32);
    const auto *nsg_34 = buffer.data(nsg + 34);
    const auto *nsg_35 = buffer.data(nsg + 35);
    const auto *nsg_39 = buffer.data(nsg + 39);
    const auto *nsg_40 = buffer.data(nsg + 40);
    const auto *nsg_41 = buffer.data(nsg + 41);
    const auto *nsg_42 = buffer.data(nsg + 42);
    const auto *nsg_43 = buffer.data(nsg + 43);
    const auto *nsg_44 = buffer.data(nsg + 44);
    const auto *nsg_45 = buffer.data(nsg + 45);
    const auto *nsg_46 = buffer.data(nsg + 46);
    const auto *nsg_47 = buffer.data(nsg + 47);
    const auto *nsg_48 = buffer.data(nsg + 48);
    const auto *nsg_50 = buffer.data(nsg + 50);
    const auto *nsg_51 = buffer.data(nsg + 51);
    const auto *nsg_55 = buffer.data(nsg + 55);
    const auto *nsg_56 = buffer.data(nsg + 56);
    const auto *nsg_57 = buffer.data(nsg + 57);
    const auto *nsg_58 = buffer.data(nsg + 58);
    const auto *nsg_59 = buffer.data(nsg + 59);
    const auto *nsg_60 = buffer.data(nsg + 60);
    const auto *nsg_62 = buffer.data(nsg + 62);
    const auto *nsg_63 = buffer.data(nsg + 63);
    const auto *nsg_65 = buffer.data(nsg + 65);
    const auto *nsg_70 = buffer.data(nsg + 70);
    const auto *nsg_71 = buffer.data(nsg + 71);
    const auto *nsg_72 = buffer.data(nsg + 72);
    const auto *nsg_73 = buffer.data(nsg + 73);
    const auto *nsg_74 = buffer.data(nsg + 74);
    const auto *nsg_75 = buffer.data(nsg + 75);
    const auto *nsg_76 = buffer.data(nsg + 76);
    const auto *nsg_77 = buffer.data(nsg + 77);
    const auto *nsg_78 = buffer.data(nsg + 78);
    const auto *nsg_79 = buffer.data(nsg + 79);
    const auto *nsg_80 = buffer.data(nsg + 80);
    const auto *nsg_84 = buffer.data(nsg + 84);
    const auto *nsg_85 = buffer.data(nsg + 85);
    const auto *nsg_86 = buffer.data(nsg + 86);
    const auto *nsg_87 = buffer.data(nsg + 87);
    const auto *nsg_88 = buffer.data(nsg + 88);
    const auto *nsg_89 = buffer.data(nsg + 89);
    const auto *nsg_90 = buffer.data(nsg + 90);
    const auto *nsg_93 = buffer.data(nsg + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, msg_0, nsf0_0, \
                         nsf1_0, nsg_0, nsg_1, nsg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * msg_0[k]
                 + f_1 * nsf0_0[k]
                 - f_2 * nsf1_0[k]
                 + f_3 * pc_x[k] * nsg_0[k];

        t_1[k] = f_3 * pc_y[k] * nsg_0[k];

        t_2[k] = f_3 * pc_z[k] * nsg_0[k];

        t_3[k] = f_4 * nsf0_0[k]
                 - f_5 * nsf1_0[k]
                 + f_3 * pc_y[k] * nsg_1[k];

        t_4[k] = f_3 * pc_y[k] * nsg_2[k];

        t_5[k] = f_4 * nsf0_0[k]
                 - f_5 * nsf1_0[k]
                 + f_3 * pc_z[k] * nsg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, msg_10, nsf0_1, nsf0_2, \
                         nsf1_1, nsf1_2, nsg_3, nsg_5, nsg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * nsf0_1[k]
                 - f_7 * nsf1_1[k]
                 + f_3 * pc_y[k] * nsg_3[k];

        t_7[k] = f_3 * pc_z[k] * nsg_3[k];

        t_8[k] = f_3 * pc_y[k] * nsg_5[k];

        t_9[k] = f_6 * nsf0_2[k]
                 - f_7 * nsf1_2[k]
                 + f_3 * pc_z[k] * nsg_5[k];

        t_10[k] = f_0 * msg_10[k]
                  + f_3 * pc_x[k] * nsg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, msg_12, msg_14, nsg_6, \
                         nsg_9, nsg_12, nsg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * nsg_6[k];

        t_12[k] = f_0 * msg_12[k]
                  + f_3 * pc_x[k] * nsg_12[k];

        t_13[k] = f_3 * pc_y[k] * nsg_9[k];

        t_14[k] = f_0 * msg_14[k]
                  + f_3 * pc_x[k] * nsg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, nsf0_6, nsf0_8, nsf0_9, nsf1_6, \
                         nsf1_8, nsf1_9, nsg_10, nsg_12, nsg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * nsf0_6[k]
                  - f_2 * nsf1_6[k]
                  + f_3 * pc_y[k] * nsg_10[k];

        t_16[k] = f_3 * pc_z[k] * nsg_10[k];

        t_17[k] = f_6 * nsf0_8[k]
                  - f_7 * nsf1_8[k]
                  + f_3 * pc_y[k] * nsg_12[k];

        t_18[k] = f_4 * nsf0_9[k]
                  - f_5 * nsf1_9[k]
                  + f_3 * pc_y[k] * nsg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, msh0_0, msg_0, \
                         msh1_0, nsf0_9, nsf1_9, nsg_14, nsg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * nsg_14[k];

        t_20[k] = f_1 * nsf0_9[k]
                  - f_2 * nsf1_9[k]
                  + f_3 * pc_z[k] * nsg_14[k];

        t_21[k] = pa_y[k] * msh0_0[k]
                  - f_8 * pc_y[k] * msh1_0[k];

        t_22[k] = f_9 * msg_0[k]
                  + f_3 * pc_y[k] * nsg_15[k];

        t_23[k] = f_3 * pc_z[k] * nsg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, msh0_3, msh0_5, msh0_6, \
                         msg_1, msg_3, msh1_3, msh1_5, msh1_6, nsg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * msh0_3[k]
                  + f_10 * msg_1[k]
                  - f_8 * pc_y[k] * msh1_3[k];

        t_25[k] = f_3 * pc_z[k] * nsg_16[k];

        t_26[k] = pa_y[k] * msh0_5[k]
                  - f_8 * pc_y[k] * msh1_5[k];

        t_27[k] = pa_y[k] * msh0_6[k]
                  + f_11 * msg_3[k]
                  - f_8 * pc_y[k] * msh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, msh0_9, msg_5, \
                         msg_25, msh1_9, nsg_18, nsg_20, nsg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * nsg_18[k];

        t_29[k] = f_9 * msg_5[k]
                  + f_3 * pc_y[k] * nsg_20[k];

        t_30[k] = pa_y[k] * msh0_9[k]
                  - f_8 * pc_y[k] * msh1_9[k];

        t_31[k] = f_12 * msg_25[k]
                  + f_3 * pc_x[k] * nsg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, msg_27, msg_28, msg_29, nsg_21, \
                         nsg_27, nsg_28, nsg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * nsg_21[k];

        t_33[k] = f_12 * msg_27[k]
                  + f_3 * pc_x[k] * nsg_27[k];

        t_34[k] = f_12 * msg_28[k]
                  + f_3 * pc_x[k] * nsg_28[k];

        t_35[k] = f_12 * msg_29[k]
                  + f_3 * pc_x[k] * nsg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, msg_10, nsf0_16, nsf0_17, \
                         nsf1_16, nsf1_17, nsg_25, nsg_26, nsg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * msg_10[k]
                  + f_1 * nsf0_16[k]
                  - f_2 * nsf1_16[k]
                  + f_3 * pc_y[k] * nsg_25[k];

        t_37[k] = f_3 * pc_z[k] * nsg_25[k];

        t_38[k] = f_4 * nsf0_16[k]
                  - f_5 * nsf1_16[k]
                  + f_3 * pc_z[k] * nsg_26[k];

        t_39[k] = f_6 * nsf0_17[k]
                  - f_7 * nsf1_17[k]
                  + f_3 * pc_z[k] * nsg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, msh0_0, msh0_20, \
                         msg_14, msh1_0, msh1_20, nsg_29, nsg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * msg_14[k]
                  + f_3 * pc_y[k] * nsg_29[k];

        t_41[k] = pa_y[k] * msh0_20[k]
                  - f_8 * pc_y[k] * msh1_20[k];

        t_42[k] = pa_z[k] * msh0_0[k]
                  - f_8 * pc_z[k] * msh1_0[k];

        t_43[k] = f_3 * pc_y[k] * nsg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, msh0_3, msh0_5, msg_0, \
                         msg_2, msh1_3, msh1_5, nsg_30, nsg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * msg_0[k]
                  + f_3 * pc_z[k] * nsg_30[k];

        t_45[k] = pa_z[k] * msh0_3[k]
                  - f_8 * pc_z[k] * msh1_3[k];

        t_46[k] = f_3 * pc_y[k] * nsg_32[k];

        t_47[k] = pa_z[k] * msh0_5[k]
                  + f_10 * msg_2[k]
                  - f_8 * pc_z[k] * msh1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, msh0_6, msh0_9, msg_5, \
                         msh1_6, msh1_9, nsf0_22, nsf1_22, nsg_34, \
                         nsg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * msh0_6[k]
                  - f_8 * pc_z[k] * msh1_6[k];

        t_49[k] = f_4 * nsf0_22[k]
                  - f_5 * nsf1_22[k]
                  + f_3 * pc_y[k] * nsg_34[k];

        t_50[k] = f_3 * pc_y[k] * nsg_35[k];

        t_51[k] = pa_z[k] * msh0_9[k]
                  + f_11 * msg_5[k]
                  - f_8 * pc_z[k] * msh1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, msg_40, msg_41, msg_42, \
                         msg_44, nsg_39, nsg_40, nsg_41, nsg_42, \
                         nsg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * msg_40[k]
                  + f_3 * pc_x[k] * nsg_40[k];

        t_53[k] = f_12 * msg_41[k]
                  + f_3 * pc_x[k] * nsg_41[k];

        t_54[k] = f_12 * msg_42[k]
                  + f_3 * pc_x[k] * nsg_42[k];

        t_55[k] = f_3 * pc_y[k] * nsg_39[k];

        t_56[k] = f_12 * msg_44[k]
                  + f_3 * pc_x[k] * nsg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, msh0_15, msh1_15, nsf0_27, \
                         nsf0_28, nsf1_27, nsf1_28, nsg_41, nsg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * msh0_15[k]
                  - f_8 * pc_z[k] * msh1_15[k];

        t_58[k] = f_13 * nsf0_27[k]
                  - f_14 * nsf1_27[k]
                  + f_3 * pc_y[k] * nsg_41[k];

        t_59[k] = f_6 * nsf0_28[k]
                  - f_7 * nsf1_28[k]
                  + f_3 * pc_y[k] * nsg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, msg_14, msg_45, nsf0_29, \
                         nsf0_30, nsf1_29, nsf1_30, nsg_43, nsg_44, \
                         nsg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * nsf0_29[k]
                  - f_5 * nsf1_29[k]
                  + f_3 * pc_y[k] * nsg_43[k];

        t_61[k] = f_3 * pc_y[k] * nsg_44[k];

        t_62[k] = f_9 * msg_14[k]
                  + f_1 * nsf0_29[k]
                  - f_2 * nsf1_29[k]
                  + f_3 * pc_z[k] * nsg_44[k];

        t_63[k] = f_15 * msg_45[k]
                  + f_1 * nsf0_30[k]
                  - f_2 * nsf1_30[k]
                  + f_3 * pc_x[k] * nsg_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pc_x, pc_y, pc_z, msg_15, msg_48, nsf0_33, \
                         nsf1_33, nsg_45, nsg_46, nsg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * msg_15[k]
                  + f_3 * pc_y[k] * nsg_45[k];

        t_65[k] = f_3 * pc_z[k] * nsg_45[k];

        t_66[k] = f_15 * msg_48[k]
                  + f_6 * nsf0_33[k]
                  - f_7 * nsf1_33[k]
                  + f_3 * pc_x[k] * nsg_48[k];

        t_67[k] = f_3 * pc_z[k] * nsg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, pc_z, msg_51, nsf0_30, nsf0_36, nsf1_30, \
                         nsf1_36, nsg_47, nsg_48, nsg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * nsf0_30[k]
                  - f_5 * nsf1_30[k]
                  + f_3 * pc_z[k] * nsg_47[k];

        t_69[k] = f_15 * msg_51[k]
                  + f_4 * nsf0_36[k]
                  - f_5 * nsf1_36[k]
                  + f_3 * pc_x[k] * nsg_51[k];

        t_70[k] = f_3 * pc_z[k] * nsg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, msg_20, msg_55, nsf0_32, \
                         nsf1_32, nsg_50, nsg_51, nsg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * msg_20[k]
                  + f_3 * pc_y[k] * nsg_50[k];

        t_72[k] = f_6 * nsf0_32[k]
                  - f_7 * nsf1_32[k]
                  + f_3 * pc_z[k] * nsg_50[k];

        t_73[k] = f_15 * msg_55[k]
                  + f_3 * pc_x[k] * nsg_55[k];

        t_74[k] = f_3 * pc_z[k] * nsg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, msg_25, msg_57, msg_58, msg_59, \
                         nsf0_36, nsf1_36, nsg_55, nsg_57, nsg_58, \
                         nsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_15 * msg_57[k]
                  + f_3 * pc_x[k] * nsg_57[k];

        t_76[k] = f_15 * msg_58[k]
                  + f_3 * pc_x[k] * nsg_58[k];

        t_77[k] = f_15 * msg_59[k]
                  + f_3 * pc_x[k] * nsg_59[k];

        t_78[k] = f_10 * msg_25[k]
                  + f_1 * nsf0_36[k]
                  - f_2 * nsf1_36[k]
                  + f_3 * pc_y[k] * nsg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_y, pc_z, msg_29, nsf0_36, nsf0_37, \
                         nsf1_36, nsf1_37, nsg_55, nsg_56, nsg_57, \
                         nsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * nsg_55[k];

        t_80[k] = f_4 * nsf0_36[k]
                  - f_5 * nsf1_36[k]
                  + f_3 * pc_z[k] * nsg_56[k];

        t_81[k] = f_6 * nsf0_37[k]
                  - f_7 * nsf1_37[k]
                  + f_3 * pc_z[k] * nsg_57[k];

        t_82[k] = f_10 * msg_29[k]
                  + f_3 * pc_y[k] * nsg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_y, pc_z, msh0_42, msg_15, msg_30, \
                         msh1_42, nsf0_39, nsf1_39, nsg_59, nsg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * nsf0_39[k]
                  - f_2 * nsf1_39[k]
                  + f_3 * pc_z[k] * nsg_59[k];

        t_84[k] = pa_y[k] * msh0_42[k]
                  - f_8 * pc_y[k] * msh1_42[k];

        t_85[k] = f_9 * msg_30[k]
                  + f_3 * pc_y[k] * nsg_60[k];

        t_86[k] = f_9 * msg_15[k]
                  + f_3 * pc_z[k] * nsg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, msh0_24, msh0_27, \
                         msh0_47, msg_32, msh1_24, msh1_27, msh1_47, \
                         nsg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * msh0_24[k]
                  - f_8 * pc_z[k] * msh1_24[k];

        t_88[k] = f_9 * msg_32[k]
                  + f_3 * pc_y[k] * nsg_62[k];

        t_89[k] = pa_y[k] * msh0_47[k]
                  - f_8 * pc_y[k] * msh1_47[k];

        t_90[k] = pa_z[k] * msh0_27[k]
                  - f_8 * pc_z[k] * msh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, msh0_51, msg_18, \
                         msg_35, msg_70, msh1_51, nsg_63, nsg_65, \
                         nsg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * msg_18[k]
                  + f_3 * pc_z[k] * nsg_63[k];

        t_92[k] = f_9 * msg_35[k]
                  + f_3 * pc_y[k] * nsg_65[k];

        t_93[k] = pa_y[k] * msh0_51[k]
                  - f_8 * pc_y[k] * msh1_51[k];

        t_94[k] = f_15 * msg_70[k]
                  + f_3 * pc_x[k] * nsg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, msg_71, msg_72, msg_73, msg_74, nsg_71, \
                         nsg_72, nsg_73, nsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_15 * msg_71[k]
                  + f_3 * pc_x[k] * nsg_71[k];

        t_96[k] = f_15 * msg_72[k]
                  + f_3 * pc_x[k] * nsg_72[k];

        t_97[k] = f_15 * msg_73[k]
                  + f_3 * pc_x[k] * nsg_73[k];

        t_98[k] = f_15 * msg_74[k]
                  + f_3 * pc_x[k] * nsg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, msh0_36, msg_25, msg_42, \
                         msh1_36, nsf0_48, nsf1_48, nsg_70, nsg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * msh0_36[k]
                  - f_8 * pc_z[k] * msh1_36[k];

        t_100[k] = f_9 * msg_25[k]
                   + f_3 * pc_z[k] * nsg_70[k];

        t_101[k] = f_9 * msg_42[k]
                   + f_6 * nsf0_48[k]
                   - f_7 * nsf1_48[k]
                   + f_3 * pc_y[k] * nsg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_y, msh0_62, msg_43, msg_44, msh1_62, \
                         nsf0_49, nsf1_49, nsg_73, nsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * msg_43[k]
                   + f_4 * nsf0_49[k]
                   - f_5 * nsf1_49[k]
                   + f_3 * pc_y[k] * nsg_73[k];

        t_103[k] = f_9 * msg_44[k]
                   + f_3 * pc_y[k] * nsg_74[k];

        t_104[k] = pa_y[k] * msh0_62[k]
                   - f_8 * pc_y[k] * msh1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, msg_30, msg_75, \
                         nsf0_50, nsf1_50, nsg_75, nsg_76, nsg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_15 * msg_75[k]
                   + f_1 * nsf0_50[k]
                   - f_2 * nsf1_50[k]
                   + f_3 * pc_x[k] * nsg_75[k];

        t_106[k] = f_3 * pc_y[k] * nsg_75[k];

        t_107[k] = f_10 * msg_30[k]
                   + f_3 * pc_z[k] * nsg_75[k];

        t_108[k] = f_4 * nsf0_50[k]
                   - f_5 * nsf1_50[k]
                   + f_3 * pc_y[k] * nsg_76[k];

        t_109[k] = f_3 * pc_y[k] * nsg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_y, msg_80, nsf0_51, nsf0_52, \
                         nsf0_55, nsf1_51, nsf1_52, nsf1_55, nsg_78, nsg_79, \
                         nsg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_15 * msg_80[k]
                   + f_6 * nsf0_55[k]
                   - f_7 * nsf1_55[k]
                   + f_3 * pc_x[k] * nsg_80[k];

        t_111[k] = f_6 * nsf0_51[k]
                   - f_7 * nsf1_51[k]
                   + f_3 * pc_y[k] * nsg_78[k];

        t_112[k] = f_4 * nsf0_52[k]
                   - f_5 * nsf1_52[k]
                   + f_3 * pc_y[k] * nsg_79[k];

        t_113[k] = f_3 * pc_y[k] * nsg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, msg_84, msg_85, msg_86, msg_87, \
                         nsf0_59, nsf1_59, nsg_84, nsg_85, nsg_86, \
                         nsg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_15 * msg_84[k]
                   + f_4 * nsf0_59[k]
                   - f_5 * nsf1_59[k]
                   + f_3 * pc_x[k] * nsg_84[k];

        t_115[k] = f_15 * msg_85[k]
                   + f_3 * pc_x[k] * nsg_85[k];

        t_116[k] = f_15 * msg_86[k]
                   + f_3 * pc_x[k] * nsg_86[k];

        t_117[k] = f_15 * msg_87[k]
                   + f_3 * pc_x[k] * nsg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, msg_89, nsf0_56, nsf0_57, \
                         nsf1_56, nsf1_57, nsg_84, nsg_85, nsg_86, \
                         nsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * pc_y[k] * nsg_84[k];

        t_119[k] = f_15 * msg_89[k]
                   + f_3 * pc_x[k] * nsg_89[k];

        t_120[k] = f_1 * nsf0_56[k]
                   - f_2 * nsf1_56[k]
                   + f_3 * pc_y[k] * nsg_85[k];

        t_121[k] = f_13 * nsf0_57[k]
                   - f_14 * nsf1_57[k]
                   + f_3 * pc_y[k] * nsg_86[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, msg_44, nsf0_58, nsf0_59, \
                         nsf1_58, nsf1_59, nsg_87, nsg_88, nsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_6 * nsf0_58[k]
                   - f_7 * nsf1_58[k]
                   + f_3 * pc_y[k] * nsg_87[k];

        t_123[k] = f_4 * nsf0_59[k]
                   - f_5 * nsf1_59[k]
                   + f_3 * pc_y[k] * nsg_88[k];

        t_124[k] = f_3 * pc_y[k] * nsg_89[k];

        t_125[k] = f_10 * msg_44[k]
                   + f_1 * nsf0_59[k]
                   - f_2 * nsf1_59[k]
                   + f_3 * pc_z[k] * nsg_89[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, msg_45, msg_90, msg_93, \
                         nsf0_60, nsf0_63, nsf1_60, nsf1_63, nsg_90, \
                         nsg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_16 * msg_90[k]
                   + f_1 * nsf0_60[k]
                   - f_2 * nsf1_60[k]
                   + f_3 * pc_x[k] * nsg_90[k];

        t_127[k] = f_11 * msg_45[k]
                   + f_3 * pc_y[k] * nsg_90[k];

        t_128[k] = f_3 * pc_z[k] * nsg_90[k];

        t_129[k] = f_16 * msg_93[k]
                   + f_6 * nsf0_63[k]
                   - f_7 * nsf1_63[k]
                   + f_3 * pc_x[k] * nsg_93[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_63 = buffer.data(msh0 + 63);
    const auto *msh0_66 = buffer.data(msh0 + 66);
    const auto *msh0_69 = buffer.data(msh0 + 69);
    const auto *msh0_78 = buffer.data(msh0 + 78);
    const auto *msh0_105 = buffer.data(msh0 + 105);
    const auto *msh0_108 = buffer.data(msh0 + 108);
    const auto *msh0_110 = buffer.data(msh0 + 110);
    const auto *msh0_111 = buffer.data(msh0 + 111);
    const auto *msh0_114 = buffer.data(msh0 + 114);
    const auto *msh0_125 = buffer.data(msh0 + 125);
    const auto *msh0_126 = buffer.data(msh0 + 126);
    const auto *msh0_129 = buffer.data(msh0 + 129);
    const auto *msh0_132 = buffer.data(msh0 + 132);
    const auto *msh0_141 = buffer.data(msh0 + 141);

    const auto *msg_45 = buffer.data(msg + 45);
    const auto *msg_48 = buffer.data(msg + 48);
    const auto *msg_50 = buffer.data(msg + 50);
    const auto *msg_55 = buffer.data(msg + 55);
    const auto *msg_59 = buffer.data(msg + 59);
    const auto *msg_60 = buffer.data(msg + 60);
    const auto *msg_62 = buffer.data(msg + 62);
    const auto *msg_63 = buffer.data(msg + 63);
    const auto *msg_65 = buffer.data(msg + 65);
    const auto *msg_70 = buffer.data(msg + 70);
    const auto *msg_72 = buffer.data(msg + 72);
    const auto *msg_73 = buffer.data(msg + 73);
    const auto *msg_74 = buffer.data(msg + 74);
    const auto *msg_75 = buffer.data(msg + 75);
    const auto *msg_76 = buffer.data(msg + 76);
    const auto *msg_77 = buffer.data(msg + 77);
    const auto *msg_78 = buffer.data(msg + 78);
    const auto *msg_80 = buffer.data(msg + 80);
    const auto *msg_85 = buffer.data(msg + 85);
    const auto *msg_87 = buffer.data(msg + 87);
    const auto *msg_88 = buffer.data(msg + 88);
    const auto *msg_89 = buffer.data(msg + 89);
    const auto *msg_90 = buffer.data(msg + 90);
    const auto *msg_93 = buffer.data(msg + 93);
    const auto *msg_95 = buffer.data(msg + 95);
    const auto *msg_96 = buffer.data(msg + 96);
    const auto *msg_100 = buffer.data(msg + 100);
    const auto *msg_102 = buffer.data(msg + 102);
    const auto *msg_103 = buffer.data(msg + 103);
    const auto *msg_104 = buffer.data(msg + 104);
    const auto *msg_105 = buffer.data(msg + 105);
    const auto *msg_107 = buffer.data(msg + 107);
    const auto *msg_110 = buffer.data(msg + 110);
    const auto *msg_114 = buffer.data(msg + 114);
    const auto *msg_115 = buffer.data(msg + 115);
    const auto *msg_116 = buffer.data(msg + 116);
    const auto *msg_117 = buffer.data(msg + 117);
    const auto *msg_118 = buffer.data(msg + 118);
    const auto *msg_119 = buffer.data(msg + 119);
    const auto *msg_130 = buffer.data(msg + 130);
    const auto *msg_131 = buffer.data(msg + 131);
    const auto *msg_132 = buffer.data(msg + 132);
    const auto *msg_133 = buffer.data(msg + 133);
    const auto *msg_134 = buffer.data(msg + 134);
    const auto *msg_135 = buffer.data(msg + 135);
    const auto *msg_140 = buffer.data(msg + 140);
    const auto *msg_144 = buffer.data(msg + 144);
    const auto *msg_145 = buffer.data(msg + 145);
    const auto *msg_146 = buffer.data(msg + 146);
    const auto *msg_147 = buffer.data(msg + 147);
    const auto *msg_149 = buffer.data(msg + 149);
    const auto *msg_150 = buffer.data(msg + 150);
    const auto *msg_153 = buffer.data(msg + 153);
    const auto *msg_156 = buffer.data(msg + 156);
    const auto *msg_160 = buffer.data(msg + 160);
    const auto *msg_162 = buffer.data(msg + 162);
    const auto *msg_163 = buffer.data(msg + 163);
    const auto *msg_164 = buffer.data(msg + 164);
    const auto *msg_170 = buffer.data(msg + 170);
    const auto *msg_174 = buffer.data(msg + 174);
    const auto *msg_175 = buffer.data(msg + 175);
    const auto *msg_176 = buffer.data(msg + 176);
    const auto *msg_177 = buffer.data(msg + 177);
    const auto *msg_178 = buffer.data(msg + 178);
    const auto *msg_179 = buffer.data(msg + 179);

    const auto *msh1_63 = buffer.data(msh1 + 63);
    const auto *msh1_66 = buffer.data(msh1 + 66);
    const auto *msh1_69 = buffer.data(msh1 + 69);
    const auto *msh1_78 = buffer.data(msh1 + 78);
    const auto *msh1_105 = buffer.data(msh1 + 105);
    const auto *msh1_108 = buffer.data(msh1 + 108);
    const auto *msh1_110 = buffer.data(msh1 + 110);
    const auto *msh1_111 = buffer.data(msh1 + 111);
    const auto *msh1_114 = buffer.data(msh1 + 114);
    const auto *msh1_125 = buffer.data(msh1 + 125);
    const auto *msh1_126 = buffer.data(msh1 + 126);
    const auto *msh1_129 = buffer.data(msh1 + 129);
    const auto *msh1_132 = buffer.data(msh1 + 132);
    const auto *msh1_141 = buffer.data(msh1 + 141);

    const auto *nsf0_60 = buffer.data(nsf0 + 60);
    const auto *nsf0_62 = buffer.data(nsf0 + 62);
    const auto *nsf0_66 = buffer.data(nsf0 + 66);
    const auto *nsf0_67 = buffer.data(nsf0 + 67);
    const auto *nsf0_69 = buffer.data(nsf0 + 69);
    const auto *nsf0_75 = buffer.data(nsf0 + 75);
    const auto *nsf0_78 = buffer.data(nsf0 + 78);
    const auto *nsf0_79 = buffer.data(nsf0 + 79);
    const auto *nsf0_86 = buffer.data(nsf0 + 86);
    const auto *nsf0_88 = buffer.data(nsf0 + 88);
    const auto *nsf0_89 = buffer.data(nsf0 + 89);
    const auto *nsf0_90 = buffer.data(nsf0 + 90);
    const auto *nsf0_91 = buffer.data(nsf0 + 91);
    const auto *nsf0_92 = buffer.data(nsf0 + 92);
    const auto *nsf0_95 = buffer.data(nsf0 + 95);
    const auto *nsf0_96 = buffer.data(nsf0 + 96);
    const auto *nsf0_97 = buffer.data(nsf0 + 97);
    const auto *nsf0_98 = buffer.data(nsf0 + 98);
    const auto *nsf0_99 = buffer.data(nsf0 + 99);
    const auto *nsf0_100 = buffer.data(nsf0 + 100);
    const auto *nsf0_102 = buffer.data(nsf0 + 102);
    const auto *nsf0_103 = buffer.data(nsf0 + 103);
    const auto *nsf0_106 = buffer.data(nsf0 + 106);
    const auto *nsf0_107 = buffer.data(nsf0 + 107);
    const auto *nsf0_109 = buffer.data(nsf0 + 109);
    const auto *nsf0_115 = buffer.data(nsf0 + 115);
    const auto *nsf0_118 = buffer.data(nsf0 + 118);
    const auto *nsf0_119 = buffer.data(nsf0 + 119);

    const auto *nsf1_60 = buffer.data(nsf1 + 60);
    const auto *nsf1_62 = buffer.data(nsf1 + 62);
    const auto *nsf1_66 = buffer.data(nsf1 + 66);
    const auto *nsf1_67 = buffer.data(nsf1 + 67);
    const auto *nsf1_69 = buffer.data(nsf1 + 69);
    const auto *nsf1_75 = buffer.data(nsf1 + 75);
    const auto *nsf1_78 = buffer.data(nsf1 + 78);
    const auto *nsf1_79 = buffer.data(nsf1 + 79);
    const auto *nsf1_86 = buffer.data(nsf1 + 86);
    const auto *nsf1_88 = buffer.data(nsf1 + 88);
    const auto *nsf1_89 = buffer.data(nsf1 + 89);
    const auto *nsf1_90 = buffer.data(nsf1 + 90);
    const auto *nsf1_91 = buffer.data(nsf1 + 91);
    const auto *nsf1_92 = buffer.data(nsf1 + 92);
    const auto *nsf1_95 = buffer.data(nsf1 + 95);
    const auto *nsf1_96 = buffer.data(nsf1 + 96);
    const auto *nsf1_97 = buffer.data(nsf1 + 97);
    const auto *nsf1_98 = buffer.data(nsf1 + 98);
    const auto *nsf1_99 = buffer.data(nsf1 + 99);
    const auto *nsf1_100 = buffer.data(nsf1 + 100);
    const auto *nsf1_102 = buffer.data(nsf1 + 102);
    const auto *nsf1_103 = buffer.data(nsf1 + 103);
    const auto *nsf1_106 = buffer.data(nsf1 + 106);
    const auto *nsf1_107 = buffer.data(nsf1 + 107);
    const auto *nsf1_109 = buffer.data(nsf1 + 109);
    const auto *nsf1_115 = buffer.data(nsf1 + 115);
    const auto *nsf1_118 = buffer.data(nsf1 + 118);
    const auto *nsf1_119 = buffer.data(nsf1 + 119);

    const auto *nsg_91 = buffer.data(nsg + 91);
    const auto *nsg_92 = buffer.data(nsg + 92);
    const auto *nsg_93 = buffer.data(nsg + 93);
    const auto *nsg_95 = buffer.data(nsg + 95);
    const auto *nsg_96 = buffer.data(nsg + 96);
    const auto *nsg_100 = buffer.data(nsg + 100);
    const auto *nsg_101 = buffer.data(nsg + 101);
    const auto *nsg_102 = buffer.data(nsg + 102);
    const auto *nsg_103 = buffer.data(nsg + 103);
    const auto *nsg_104 = buffer.data(nsg + 104);
    const auto *nsg_105 = buffer.data(nsg + 105);
    const auto *nsg_107 = buffer.data(nsg + 107);
    const auto *nsg_108 = buffer.data(nsg + 108);
    const auto *nsg_110 = buffer.data(nsg + 110);
    const auto *nsg_114 = buffer.data(nsg + 114);
    const auto *nsg_115 = buffer.data(nsg + 115);
    const auto *nsg_116 = buffer.data(nsg + 116);
    const auto *nsg_117 = buffer.data(nsg + 117);
    const auto *nsg_118 = buffer.data(nsg + 118);
    const auto *nsg_119 = buffer.data(nsg + 119);
    const auto *nsg_120 = buffer.data(nsg + 120);
    const auto *nsg_122 = buffer.data(nsg + 122);
    const auto *nsg_123 = buffer.data(nsg + 123);
    const auto *nsg_125 = buffer.data(nsg + 125);
    const auto *nsg_130 = buffer.data(nsg + 130);
    const auto *nsg_131 = buffer.data(nsg + 131);
    const auto *nsg_132 = buffer.data(nsg + 132);
    const auto *nsg_133 = buffer.data(nsg + 133);
    const auto *nsg_134 = buffer.data(nsg + 134);
    const auto *nsg_135 = buffer.data(nsg + 135);
    const auto *nsg_136 = buffer.data(nsg + 136);
    const auto *nsg_137 = buffer.data(nsg + 137);
    const auto *nsg_138 = buffer.data(nsg + 138);
    const auto *nsg_139 = buffer.data(nsg + 139);
    const auto *nsg_140 = buffer.data(nsg + 140);
    const auto *nsg_144 = buffer.data(nsg + 144);
    const auto *nsg_145 = buffer.data(nsg + 145);
    const auto *nsg_146 = buffer.data(nsg + 146);
    const auto *nsg_147 = buffer.data(nsg + 147);
    const auto *nsg_148 = buffer.data(nsg + 148);
    const auto *nsg_149 = buffer.data(nsg + 149);
    const auto *nsg_150 = buffer.data(nsg + 150);
    const auto *nsg_151 = buffer.data(nsg + 151);
    const auto *nsg_152 = buffer.data(nsg + 152);
    const auto *nsg_153 = buffer.data(nsg + 153);
    const auto *nsg_155 = buffer.data(nsg + 155);
    const auto *nsg_156 = buffer.data(nsg + 156);
    const auto *nsg_160 = buffer.data(nsg + 160);
    const auto *nsg_161 = buffer.data(nsg + 161);
    const auto *nsg_162 = buffer.data(nsg + 162);
    const auto *nsg_163 = buffer.data(nsg + 163);
    const auto *nsg_164 = buffer.data(nsg + 164);
    const auto *nsg_165 = buffer.data(nsg + 165);
    const auto *nsg_167 = buffer.data(nsg + 167);
    const auto *nsg_168 = buffer.data(nsg + 168);
    const auto *nsg_170 = buffer.data(nsg + 170);
    const auto *nsg_174 = buffer.data(nsg + 174);
    const auto *nsg_175 = buffer.data(nsg + 175);
    const auto *nsg_176 = buffer.data(nsg + 176);
    const auto *nsg_177 = buffer.data(nsg + 177);
    const auto *nsg_178 = buffer.data(nsg + 178);
    const auto *nsg_179 = buffer.data(nsg + 179);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_x, pc_z, msg_96, nsf0_60, nsf0_66, \
                         nsf1_60, nsf1_66, nsg_91, nsg_92, nsg_93, \
                         nsg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_z[k] * nsg_91[k];

        t_131[k] = f_4 * nsf0_60[k]
                   - f_5 * nsf1_60[k]
                   + f_3 * pc_z[k] * nsg_92[k];

        t_132[k] = f_16 * msg_96[k]
                   + f_4 * nsf0_66[k]
                   - f_5 * nsf1_66[k]
                   + f_3 * pc_x[k] * nsg_96[k];

        t_133[k] = f_3 * pc_z[k] * nsg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, msg_50, msg_100, \
                         nsf0_62, nsf1_62, nsg_95, nsg_96, nsg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * msg_50[k]
                   + f_3 * pc_y[k] * nsg_95[k];

        t_135[k] = f_6 * nsf0_62[k]
                   - f_7 * nsf1_62[k]
                   + f_3 * pc_z[k] * nsg_95[k];

        t_136[k] = f_16 * msg_100[k]
                   + f_3 * pc_x[k] * nsg_100[k];

        t_137[k] = f_3 * pc_z[k] * nsg_96[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pc_x, pc_y, msg_55, msg_102, msg_103, \
                         msg_104, nsf0_66, nsf1_66, nsg_100, nsg_102, nsg_103, \
                         nsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_16 * msg_102[k]
                   + f_3 * pc_x[k] * nsg_102[k];

        t_139[k] = f_16 * msg_103[k]
                   + f_3 * pc_x[k] * nsg_103[k];

        t_140[k] = f_16 * msg_104[k]
                   + f_3 * pc_x[k] * nsg_104[k];

        t_141[k] = f_11 * msg_55[k]
                   + f_1 * nsf0_66[k]
                   - f_2 * nsf1_66[k]
                   + f_3 * pc_y[k] * nsg_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, pc_z, msg_59, nsf0_66, nsf0_67, \
                         nsf1_66, nsf1_67, nsg_100, nsg_101, nsg_102, \
                         nsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * nsg_100[k];

        t_143[k] = f_4 * nsf0_66[k]
                   - f_5 * nsf1_66[k]
                   + f_3 * pc_z[k] * nsg_101[k];

        t_144[k] = f_6 * nsf0_67[k]
                   - f_7 * nsf1_67[k]
                   + f_3 * pc_z[k] * nsg_102[k];

        t_145[k] = f_11 * msg_59[k]
                   + f_3 * pc_y[k] * nsg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, msh0_63, msg_45, \
                         msg_60, msh1_63, nsf0_69, nsf1_69, nsg_104, \
                         nsg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * nsf0_69[k]
                   - f_2 * nsf1_69[k]
                   + f_3 * pc_z[k] * nsg_104[k];

        t_147[k] = pa_z[k] * msh0_63[k]
                   - f_8 * pc_z[k] * msh1_63[k];

        t_148[k] = f_10 * msg_60[k]
                   + f_3 * pc_y[k] * nsg_105[k];

        t_149[k] = f_9 * msg_45[k]
                   + f_3 * pc_z[k] * nsg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_z, pc_x, pc_y, pc_z, msh0_66, msg_62, \
                         msg_110, msh1_66, nsf0_75, nsf1_75, nsg_107, \
                         nsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * msh0_66[k]
                   - f_8 * pc_z[k] * msh1_66[k];

        t_151[k] = f_10 * msg_62[k]
                   + f_3 * pc_y[k] * nsg_107[k];

        t_152[k] = f_16 * msg_110[k]
                   + f_6 * nsf0_75[k]
                   - f_7 * nsf1_75[k]
                   + f_3 * pc_x[k] * nsg_110[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_y, pc_z, msh0_69, msg_48, msg_65, \
                         msh1_69, nsg_108, nsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * msh0_69[k]
                   - f_8 * pc_z[k] * msh1_69[k];

        t_154[k] = f_9 * msg_48[k]
                   + f_3 * pc_z[k] * nsg_108[k];

        t_155[k] = f_10 * msg_65[k]
                   + f_3 * pc_y[k] * nsg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, msg_114, msg_115, msg_116, msg_117, \
                         nsf0_79, nsf1_79, nsg_114, nsg_115, nsg_116, \
                         nsg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * msg_114[k]
                   + f_4 * nsf0_79[k]
                   - f_5 * nsf1_79[k]
                   + f_3 * pc_x[k] * nsg_114[k];

        t_157[k] = f_16 * msg_115[k]
                   + f_3 * pc_x[k] * nsg_115[k];

        t_158[k] = f_16 * msg_116[k]
                   + f_3 * pc_x[k] * nsg_116[k];

        t_159[k] = f_16 * msg_117[k]
                   + f_3 * pc_x[k] * nsg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pc_x, pc_z, msh0_78, msg_55, \
                         msg_118, msg_119, msh1_78, nsg_115, nsg_118, \
                         nsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * msg_118[k]
                   + f_3 * pc_x[k] * nsg_118[k];

        t_161[k] = f_16 * msg_119[k]
                   + f_3 * pc_x[k] * nsg_119[k];

        t_162[k] = pa_z[k] * msh0_78[k]
                   - f_8 * pc_z[k] * msh1_78[k];

        t_163[k] = f_9 * msg_55[k]
                   + f_3 * pc_z[k] * nsg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pc_y, msg_72, msg_73, msg_74, nsf0_78, nsf0_79, \
                         nsf1_78, nsf1_79, nsg_117, nsg_118, nsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_10 * msg_72[k]
                   + f_6 * nsf0_78[k]
                   - f_7 * nsf1_78[k]
                   + f_3 * pc_y[k] * nsg_117[k];

        t_165[k] = f_10 * msg_73[k]
                   + f_4 * nsf0_79[k]
                   - f_5 * nsf1_79[k]
                   + f_3 * pc_y[k] * nsg_118[k];

        t_166[k] = f_10 * msg_74[k]
                   + f_3 * pc_y[k] * nsg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pc_y, pc_z, msh0_105, msg_59, \
                         msg_60, msg_75, msh1_105, nsf0_79, nsf1_79, nsg_119, \
                         nsg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_9 * msg_59[k]
                   + f_1 * nsf0_79[k]
                   - f_2 * nsf1_79[k]
                   + f_3 * pc_z[k] * nsg_119[k];

        t_168[k] = pa_y[k] * msh0_105[k]
                   - f_8 * pc_y[k] * msh1_105[k];

        t_169[k] = f_9 * msg_75[k]
                   + f_3 * pc_y[k] * nsg_120[k];

        t_170[k] = f_10 * msg_60[k]
                   + f_3 * pc_z[k] * nsg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_y, pc_y, msh0_108, msh0_110, msh0_111, \
                         msg_76, msg_77, msg_78, msh1_108, msh1_110, msh1_111, \
                         nsg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_y[k] * msh0_108[k]
                   + f_10 * msg_76[k]
                   - f_8 * pc_y[k] * msh1_108[k];

        t_172[k] = f_9 * msg_77[k]
                   + f_3 * pc_y[k] * nsg_122[k];

        t_173[k] = pa_y[k] * msh0_110[k]
                   - f_8 * pc_y[k] * msh1_110[k];

        t_174[k] = pa_y[k] * msh0_111[k]
                   + f_11 * msg_78[k]
                   - f_8 * pc_y[k] * msh1_111[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pc_x, pc_y, pc_z, msh0_114, msg_63, \
                         msg_80, msg_130, msh1_114, nsg_123, nsg_125, \
                         nsg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * msg_63[k]
                   + f_3 * pc_z[k] * nsg_123[k];

        t_176[k] = f_9 * msg_80[k]
                   + f_3 * pc_y[k] * nsg_125[k];

        t_177[k] = pa_y[k] * msh0_114[k]
                   - f_8 * pc_y[k] * msh1_114[k];

        t_178[k] = f_16 * msg_130[k]
                   + f_3 * pc_x[k] * nsg_130[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, msg_131, msg_132, msg_133, msg_134, \
                         nsg_131, nsg_132, nsg_133, nsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_16 * msg_131[k]
                   + f_3 * pc_x[k] * nsg_131[k];

        t_180[k] = f_16 * msg_132[k]
                   + f_3 * pc_x[k] * nsg_132[k];

        t_181[k] = f_16 * msg_133[k]
                   + f_3 * pc_x[k] * nsg_133[k];

        t_182[k] = f_16 * msg_134[k]
                   + f_3 * pc_x[k] * nsg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, msg_70, msg_85, msg_87, nsf0_86, \
                         nsf0_88, nsf1_86, nsf1_88, nsg_130, nsg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * msg_85[k]
                   + f_1 * nsf0_86[k]
                   - f_2 * nsf1_86[k]
                   + f_3 * pc_y[k] * nsg_130[k];

        t_184[k] = f_10 * msg_70[k]
                   + f_3 * pc_z[k] * nsg_130[k];

        t_185[k] = f_9 * msg_87[k]
                   + f_6 * nsf0_88[k]
                   - f_7 * nsf1_88[k]
                   + f_3 * pc_y[k] * nsg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pc_y, msh0_125, msg_88, msg_89, msh1_125, \
                         nsf0_89, nsf1_89, nsg_133, nsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * msg_88[k]
                   + f_4 * nsf0_89[k]
                   - f_5 * nsf1_89[k]
                   + f_3 * pc_y[k] * nsg_133[k];

        t_187[k] = f_9 * msg_89[k]
                   + f_3 * pc_y[k] * nsg_134[k];

        t_188[k] = pa_y[k] * msh0_125[k]
                   - f_8 * pc_y[k] * msh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, msg_75, msg_135, \
                         nsf0_90, nsf1_90, nsg_135, nsg_136, nsg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_16 * msg_135[k]
                   + f_1 * nsf0_90[k]
                   - f_2 * nsf1_90[k]
                   + f_3 * pc_x[k] * nsg_135[k];

        t_190[k] = f_3 * pc_y[k] * nsg_135[k];

        t_191[k] = f_11 * msg_75[k]
                   + f_3 * pc_z[k] * nsg_135[k];

        t_192[k] = f_4 * nsf0_90[k]
                   - f_5 * nsf1_90[k]
                   + f_3 * pc_y[k] * nsg_136[k];

        t_193[k] = f_3 * pc_y[k] * nsg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, msg_140, nsf0_91, nsf0_92, \
                         nsf0_95, nsf1_91, nsf1_92, nsf1_95, nsg_138, nsg_139, \
                         nsg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_16 * msg_140[k]
                   + f_6 * nsf0_95[k]
                   - f_7 * nsf1_95[k]
                   + f_3 * pc_x[k] * nsg_140[k];

        t_195[k] = f_6 * nsf0_91[k]
                   - f_7 * nsf1_91[k]
                   + f_3 * pc_y[k] * nsg_138[k];

        t_196[k] = f_4 * nsf0_92[k]
                   - f_5 * nsf1_92[k]
                   + f_3 * pc_y[k] * nsg_139[k];

        t_197[k] = f_3 * pc_y[k] * nsg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, msg_144, msg_145, msg_146, msg_147, \
                         nsf0_99, nsf1_99, nsg_144, nsg_145, nsg_146, \
                         nsg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_16 * msg_144[k]
                   + f_4 * nsf0_99[k]
                   - f_5 * nsf1_99[k]
                   + f_3 * pc_x[k] * nsg_144[k];

        t_199[k] = f_16 * msg_145[k]
                   + f_3 * pc_x[k] * nsg_145[k];

        t_200[k] = f_16 * msg_146[k]
                   + f_3 * pc_x[k] * nsg_146[k];

        t_201[k] = f_16 * msg_147[k]
                   + f_3 * pc_x[k] * nsg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pc_x, pc_y, msg_149, nsf0_96, nsf0_97, \
                         nsf1_96, nsf1_97, nsg_144, nsg_145, nsg_146, \
                         nsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_3 * pc_y[k] * nsg_144[k];

        t_203[k] = f_16 * msg_149[k]
                   + f_3 * pc_x[k] * nsg_149[k];

        t_204[k] = f_1 * nsf0_96[k]
                   - f_2 * nsf1_96[k]
                   + f_3 * pc_y[k] * nsg_145[k];

        t_205[k] = f_13 * nsf0_97[k]
                   - f_14 * nsf1_97[k]
                   + f_3 * pc_y[k] * nsg_146[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_y, pc_z, msg_89, nsf0_98, nsf0_99, \
                         nsf1_98, nsf1_99, nsg_147, nsg_148, nsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_6 * nsf0_98[k]
                   - f_7 * nsf1_98[k]
                   + f_3 * pc_y[k] * nsg_147[k];

        t_207[k] = f_4 * nsf0_99[k]
                   - f_5 * nsf1_99[k]
                   + f_3 * pc_y[k] * nsg_148[k];

        t_208[k] = f_3 * pc_y[k] * nsg_149[k];

        t_209[k] = f_11 * msg_89[k]
                   + f_1 * nsf0_99[k]
                   - f_2 * nsf1_99[k]
                   + f_3 * pc_z[k] * nsg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, pc_y, pc_z, msg_90, msg_150, \
                         msg_153, nsf0_100, nsf0_103, nsf1_100, nsf1_103, nsg_150, \
                         nsg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * msg_150[k]
                   + f_1 * nsf0_100[k]
                   - f_2 * nsf1_100[k]
                   + f_3 * pc_x[k] * nsg_150[k];

        t_211[k] = f_18 * msg_90[k]
                   + f_3 * pc_y[k] * nsg_150[k];

        t_212[k] = f_3 * pc_z[k] * nsg_150[k];

        t_213[k] = f_17 * msg_153[k]
                   + f_6 * nsf0_103[k]
                   - f_7 * nsf1_103[k]
                   + f_3 * pc_x[k] * nsg_153[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pc_x, pc_z, msg_156, nsf0_100, nsf0_106, \
                         nsf1_100, nsf1_106, nsg_151, nsg_152, nsg_153, \
                         nsg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_3 * pc_z[k] * nsg_151[k];

        t_215[k] = f_4 * nsf0_100[k]
                   - f_5 * nsf1_100[k]
                   + f_3 * pc_z[k] * nsg_152[k];

        t_216[k] = f_17 * msg_156[k]
                   + f_4 * nsf0_106[k]
                   - f_5 * nsf1_106[k]
                   + f_3 * pc_x[k] * nsg_156[k];

        t_217[k] = f_3 * pc_z[k] * nsg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, msg_95, msg_160, \
                         nsf0_102, nsf1_102, nsg_155, nsg_156, \
                         nsg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_18 * msg_95[k]
                   + f_3 * pc_y[k] * nsg_155[k];

        t_219[k] = f_6 * nsf0_102[k]
                   - f_7 * nsf1_102[k]
                   + f_3 * pc_z[k] * nsg_155[k];

        t_220[k] = f_17 * msg_160[k]
                   + f_3 * pc_x[k] * nsg_160[k];

        t_221[k] = f_3 * pc_z[k] * nsg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, msg_100, msg_162, msg_163, \
                         msg_164, nsf0_106, nsf1_106, nsg_160, nsg_162, nsg_163, \
                         nsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_17 * msg_162[k]
                   + f_3 * pc_x[k] * nsg_162[k];

        t_223[k] = f_17 * msg_163[k]
                   + f_3 * pc_x[k] * nsg_163[k];

        t_224[k] = f_17 * msg_164[k]
                   + f_3 * pc_x[k] * nsg_164[k];

        t_225[k] = f_18 * msg_100[k]
                   + f_1 * nsf0_106[k]
                   - f_2 * nsf1_106[k]
                   + f_3 * pc_y[k] * nsg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_y, pc_z, msg_104, nsf0_106, nsf0_107, \
                         nsf1_106, nsf1_107, nsg_160, nsg_161, nsg_162, \
                         nsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * nsg_160[k];

        t_227[k] = f_4 * nsf0_106[k]
                   - f_5 * nsf1_106[k]
                   + f_3 * pc_z[k] * nsg_161[k];

        t_228[k] = f_6 * nsf0_107[k]
                   - f_7 * nsf1_107[k]
                   + f_3 * pc_z[k] * nsg_162[k];

        t_229[k] = f_18 * msg_104[k]
                   + f_3 * pc_y[k] * nsg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pa_z, pc_y, pc_z, msh0_126, msg_90, \
                         msg_105, msh1_126, nsf0_109, nsf1_109, nsg_164, \
                         nsg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * nsf0_109[k]
                   - f_2 * nsf1_109[k]
                   + f_3 * pc_z[k] * nsg_164[k];

        t_231[k] = pa_z[k] * msh0_126[k]
                   - f_8 * pc_z[k] * msh1_126[k];

        t_232[k] = f_11 * msg_105[k]
                   + f_3 * pc_y[k] * nsg_165[k];

        t_233[k] = f_9 * msg_90[k]
                   + f_3 * pc_z[k] * nsg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_z, pc_x, pc_y, pc_z, msh0_129, msg_107, \
                         msg_170, msh1_129, nsf0_115, nsf1_115, nsg_167, \
                         nsg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * msh0_129[k]
                   - f_8 * pc_z[k] * msh1_129[k];

        t_235[k] = f_11 * msg_107[k]
                   + f_3 * pc_y[k] * nsg_167[k];

        t_236[k] = f_17 * msg_170[k]
                   + f_6 * nsf0_115[k]
                   - f_7 * nsf1_115[k]
                   + f_3 * pc_x[k] * nsg_170[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_z, pc_y, pc_z, msh0_132, msg_93, msg_110, \
                         msh1_132, nsg_168, nsg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * msh0_132[k]
                   - f_8 * pc_z[k] * msh1_132[k];

        t_238[k] = f_9 * msg_93[k]
                   + f_3 * pc_z[k] * nsg_168[k];

        t_239[k] = f_11 * msg_110[k]
                   + f_3 * pc_y[k] * nsg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, msg_174, msg_175, msg_176, msg_177, \
                         nsf0_119, nsf1_119, nsg_174, nsg_175, nsg_176, \
                         nsg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * msg_174[k]
                   + f_4 * nsf0_119[k]
                   - f_5 * nsf1_119[k]
                   + f_3 * pc_x[k] * nsg_174[k];

        t_241[k] = f_17 * msg_175[k]
                   + f_3 * pc_x[k] * nsg_175[k];

        t_242[k] = f_17 * msg_176[k]
                   + f_3 * pc_x[k] * nsg_176[k];

        t_243[k] = f_17 * msg_177[k]
                   + f_3 * pc_x[k] * nsg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pc_x, pc_z, msh0_141, msg_100, \
                         msg_178, msg_179, msh1_141, nsg_175, nsg_178, \
                         nsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_17 * msg_178[k]
                   + f_3 * pc_x[k] * nsg_178[k];

        t_245[k] = f_17 * msg_179[k]
                   + f_3 * pc_x[k] * nsg_179[k];

        t_246[k] = pa_z[k] * msh0_141[k]
                   - f_8 * pc_z[k] * msh1_141[k];

        t_247[k] = f_9 * msg_100[k]
                   + f_3 * pc_z[k] * nsg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, msg_117, msg_118, msg_119, nsf0_118, \
                         nsf0_119, nsf1_118, nsf1_119, nsg_177, nsg_178, \
                         nsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * msg_117[k]
                   + f_6 * nsf0_118[k]
                   - f_7 * nsf1_118[k]
                   + f_3 * pc_y[k] * nsg_177[k];

        t_249[k] = f_11 * msg_118[k]
                   + f_4 * nsf0_119[k]
                   - f_5 * nsf1_119[k]
                   + f_3 * pc_y[k] * nsg_178[k];

        t_250[k] = f_11 * msg_119[k]
                   + f_3 * pc_y[k] * nsg_179[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

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
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);
    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);
    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_189 = buffer.data(msh0 + 189);
    const auto *msh0_192 = buffer.data(msh0 + 192);
    const auto *msh0_194 = buffer.data(msh0 + 194);
    const auto *msh0_195 = buffer.data(msh0 + 195);
    const auto *msh0_198 = buffer.data(msh0 + 198);
    const auto *msh0_209 = buffer.data(msh0 + 209);
    const auto *msh0_210 = buffer.data(msh0 + 210);
    const auto *msh0_213 = buffer.data(msh0 + 213);
    const auto *msh0_216 = buffer.data(msh0 + 216);
    const auto *msh0_225 = buffer.data(msh0 + 225);

    const auto *msg_104 = buffer.data(msg + 104);
    const auto *msg_105 = buffer.data(msg + 105);
    const auto *msg_108 = buffer.data(msg + 108);
    const auto *msg_115 = buffer.data(msg + 115);
    const auto *msg_119 = buffer.data(msg + 119);
    const auto *msg_120 = buffer.data(msg + 120);
    const auto *msg_122 = buffer.data(msg + 122);
    const auto *msg_123 = buffer.data(msg + 123);
    const auto *msg_125 = buffer.data(msg + 125);
    const auto *msg_130 = buffer.data(msg + 130);
    const auto *msg_132 = buffer.data(msg + 132);
    const auto *msg_133 = buffer.data(msg + 133);
    const auto *msg_134 = buffer.data(msg + 134);
    const auto *msg_135 = buffer.data(msg + 135);
    const auto *msg_136 = buffer.data(msg + 136);
    const auto *msg_137 = buffer.data(msg + 137);
    const auto *msg_138 = buffer.data(msg + 138);
    const auto *msg_140 = buffer.data(msg + 140);
    const auto *msg_145 = buffer.data(msg + 145);
    const auto *msg_147 = buffer.data(msg + 147);
    const auto *msg_148 = buffer.data(msg + 148);
    const auto *msg_149 = buffer.data(msg + 149);
    const auto *msg_150 = buffer.data(msg + 150);
    const auto *msg_153 = buffer.data(msg + 153);
    const auto *msg_155 = buffer.data(msg + 155);
    const auto *msg_160 = buffer.data(msg + 160);
    const auto *msg_164 = buffer.data(msg + 164);
    const auto *msg_165 = buffer.data(msg + 165);
    const auto *msg_167 = buffer.data(msg + 167);
    const auto *msg_168 = buffer.data(msg + 168);
    const auto *msg_170 = buffer.data(msg + 170);
    const auto *msg_177 = buffer.data(msg + 177);
    const auto *msg_178 = buffer.data(msg + 178);
    const auto *msg_179 = buffer.data(msg + 179);
    const auto *msg_180 = buffer.data(msg + 180);
    const auto *msg_182 = buffer.data(msg + 182);
    const auto *msg_183 = buffer.data(msg + 183);
    const auto *msg_185 = buffer.data(msg + 185);
    const auto *msg_186 = buffer.data(msg + 186);
    const auto *msg_189 = buffer.data(msg + 189);
    const auto *msg_190 = buffer.data(msg + 190);
    const auto *msg_191 = buffer.data(msg + 191);
    const auto *msg_192 = buffer.data(msg + 192);
    const auto *msg_193 = buffer.data(msg + 193);
    const auto *msg_194 = buffer.data(msg + 194);
    const auto *msg_205 = buffer.data(msg + 205);
    const auto *msg_206 = buffer.data(msg + 206);
    const auto *msg_207 = buffer.data(msg + 207);
    const auto *msg_208 = buffer.data(msg + 208);
    const auto *msg_209 = buffer.data(msg + 209);
    const auto *msg_210 = buffer.data(msg + 210);
    const auto *msg_215 = buffer.data(msg + 215);
    const auto *msg_219 = buffer.data(msg + 219);
    const auto *msg_220 = buffer.data(msg + 220);
    const auto *msg_221 = buffer.data(msg + 221);
    const auto *msg_222 = buffer.data(msg + 222);
    const auto *msg_224 = buffer.data(msg + 224);
    const auto *msg_225 = buffer.data(msg + 225);
    const auto *msg_228 = buffer.data(msg + 228);
    const auto *msg_231 = buffer.data(msg + 231);
    const auto *msg_235 = buffer.data(msg + 235);
    const auto *msg_237 = buffer.data(msg + 237);
    const auto *msg_238 = buffer.data(msg + 238);
    const auto *msg_239 = buffer.data(msg + 239);
    const auto *msg_245 = buffer.data(msg + 245);
    const auto *msg_249 = buffer.data(msg + 249);
    const auto *msg_250 = buffer.data(msg + 250);
    const auto *msg_251 = buffer.data(msg + 251);
    const auto *msg_252 = buffer.data(msg + 252);
    const auto *msg_253 = buffer.data(msg + 253);
    const auto *msg_254 = buffer.data(msg + 254);
    const auto *msg_255 = buffer.data(msg + 255);
    const auto *msg_258 = buffer.data(msg + 258);
    const auto *msg_260 = buffer.data(msg + 260);
    const auto *msg_261 = buffer.data(msg + 261);
    const auto *msg_264 = buffer.data(msg + 264);
    const auto *msg_265 = buffer.data(msg + 265);
    const auto *msg_266 = buffer.data(msg + 266);

    const auto *msh1_189 = buffer.data(msh1 + 189);
    const auto *msh1_192 = buffer.data(msh1 + 192);
    const auto *msh1_194 = buffer.data(msh1 + 194);
    const auto *msh1_195 = buffer.data(msh1 + 195);
    const auto *msh1_198 = buffer.data(msh1 + 198);
    const auto *msh1_209 = buffer.data(msh1 + 209);
    const auto *msh1_210 = buffer.data(msh1 + 210);
    const auto *msh1_213 = buffer.data(msh1 + 213);
    const auto *msh1_216 = buffer.data(msh1 + 216);
    const auto *msh1_225 = buffer.data(msh1 + 225);

    const auto *nsf0_119 = buffer.data(nsf0 + 119);
    const auto *nsf0_120 = buffer.data(nsf0 + 120);
    const auto *nsf0_123 = buffer.data(nsf0 + 123);
    const auto *nsf0_125 = buffer.data(nsf0 + 125);
    const auto *nsf0_126 = buffer.data(nsf0 + 126);
    const auto *nsf0_128 = buffer.data(nsf0 + 128);
    const auto *nsf0_129 = buffer.data(nsf0 + 129);
    const auto *nsf0_136 = buffer.data(nsf0 + 136);
    const auto *nsf0_138 = buffer.data(nsf0 + 138);
    const auto *nsf0_139 = buffer.data(nsf0 + 139);
    const auto *nsf0_140 = buffer.data(nsf0 + 140);
    const auto *nsf0_141 = buffer.data(nsf0 + 141);
    const auto *nsf0_142 = buffer.data(nsf0 + 142);
    const auto *nsf0_145 = buffer.data(nsf0 + 145);
    const auto *nsf0_146 = buffer.data(nsf0 + 146);
    const auto *nsf0_147 = buffer.data(nsf0 + 147);
    const auto *nsf0_148 = buffer.data(nsf0 + 148);
    const auto *nsf0_149 = buffer.data(nsf0 + 149);
    const auto *nsf0_150 = buffer.data(nsf0 + 150);
    const auto *nsf0_152 = buffer.data(nsf0 + 152);
    const auto *nsf0_153 = buffer.data(nsf0 + 153);
    const auto *nsf0_156 = buffer.data(nsf0 + 156);
    const auto *nsf0_157 = buffer.data(nsf0 + 157);
    const auto *nsf0_159 = buffer.data(nsf0 + 159);
    const auto *nsf0_165 = buffer.data(nsf0 + 165);
    const auto *nsf0_168 = buffer.data(nsf0 + 168);
    const auto *nsf0_169 = buffer.data(nsf0 + 169);
    const auto *nsf0_170 = buffer.data(nsf0 + 170);
    const auto *nsf0_173 = buffer.data(nsf0 + 173);
    const auto *nsf0_175 = buffer.data(nsf0 + 175);
    const auto *nsf0_176 = buffer.data(nsf0 + 176);
    const auto *nsf0_179 = buffer.data(nsf0 + 179);

    const auto *nsf1_119 = buffer.data(nsf1 + 119);
    const auto *nsf1_120 = buffer.data(nsf1 + 120);
    const auto *nsf1_123 = buffer.data(nsf1 + 123);
    const auto *nsf1_125 = buffer.data(nsf1 + 125);
    const auto *nsf1_126 = buffer.data(nsf1 + 126);
    const auto *nsf1_128 = buffer.data(nsf1 + 128);
    const auto *nsf1_129 = buffer.data(nsf1 + 129);
    const auto *nsf1_136 = buffer.data(nsf1 + 136);
    const auto *nsf1_138 = buffer.data(nsf1 + 138);
    const auto *nsf1_139 = buffer.data(nsf1 + 139);
    const auto *nsf1_140 = buffer.data(nsf1 + 140);
    const auto *nsf1_141 = buffer.data(nsf1 + 141);
    const auto *nsf1_142 = buffer.data(nsf1 + 142);
    const auto *nsf1_145 = buffer.data(nsf1 + 145);
    const auto *nsf1_146 = buffer.data(nsf1 + 146);
    const auto *nsf1_147 = buffer.data(nsf1 + 147);
    const auto *nsf1_148 = buffer.data(nsf1 + 148);
    const auto *nsf1_149 = buffer.data(nsf1 + 149);
    const auto *nsf1_150 = buffer.data(nsf1 + 150);
    const auto *nsf1_152 = buffer.data(nsf1 + 152);
    const auto *nsf1_153 = buffer.data(nsf1 + 153);
    const auto *nsf1_156 = buffer.data(nsf1 + 156);
    const auto *nsf1_157 = buffer.data(nsf1 + 157);
    const auto *nsf1_159 = buffer.data(nsf1 + 159);
    const auto *nsf1_165 = buffer.data(nsf1 + 165);
    const auto *nsf1_168 = buffer.data(nsf1 + 168);
    const auto *nsf1_169 = buffer.data(nsf1 + 169);
    const auto *nsf1_170 = buffer.data(nsf1 + 170);
    const auto *nsf1_173 = buffer.data(nsf1 + 173);
    const auto *nsf1_175 = buffer.data(nsf1 + 175);
    const auto *nsf1_176 = buffer.data(nsf1 + 176);
    const auto *nsf1_179 = buffer.data(nsf1 + 179);

    const auto *nsg_179 = buffer.data(nsg + 179);
    const auto *nsg_180 = buffer.data(nsg + 180);
    const auto *nsg_182 = buffer.data(nsg + 182);
    const auto *nsg_183 = buffer.data(nsg + 183);
    const auto *nsg_185 = buffer.data(nsg + 185);
    const auto *nsg_186 = buffer.data(nsg + 186);
    const auto *nsg_189 = buffer.data(nsg + 189);
    const auto *nsg_190 = buffer.data(nsg + 190);
    const auto *nsg_191 = buffer.data(nsg + 191);
    const auto *nsg_192 = buffer.data(nsg + 192);
    const auto *nsg_193 = buffer.data(nsg + 193);
    const auto *nsg_194 = buffer.data(nsg + 194);
    const auto *nsg_195 = buffer.data(nsg + 195);
    const auto *nsg_197 = buffer.data(nsg + 197);
    const auto *nsg_198 = buffer.data(nsg + 198);
    const auto *nsg_200 = buffer.data(nsg + 200);
    const auto *nsg_205 = buffer.data(nsg + 205);
    const auto *nsg_206 = buffer.data(nsg + 206);
    const auto *nsg_207 = buffer.data(nsg + 207);
    const auto *nsg_208 = buffer.data(nsg + 208);
    const auto *nsg_209 = buffer.data(nsg + 209);
    const auto *nsg_210 = buffer.data(nsg + 210);
    const auto *nsg_211 = buffer.data(nsg + 211);
    const auto *nsg_212 = buffer.data(nsg + 212);
    const auto *nsg_213 = buffer.data(nsg + 213);
    const auto *nsg_214 = buffer.data(nsg + 214);
    const auto *nsg_215 = buffer.data(nsg + 215);
    const auto *nsg_219 = buffer.data(nsg + 219);
    const auto *nsg_220 = buffer.data(nsg + 220);
    const auto *nsg_221 = buffer.data(nsg + 221);
    const auto *nsg_222 = buffer.data(nsg + 222);
    const auto *nsg_223 = buffer.data(nsg + 223);
    const auto *nsg_224 = buffer.data(nsg + 224);
    const auto *nsg_225 = buffer.data(nsg + 225);
    const auto *nsg_226 = buffer.data(nsg + 226);
    const auto *nsg_227 = buffer.data(nsg + 227);
    const auto *nsg_228 = buffer.data(nsg + 228);
    const auto *nsg_230 = buffer.data(nsg + 230);
    const auto *nsg_231 = buffer.data(nsg + 231);
    const auto *nsg_235 = buffer.data(nsg + 235);
    const auto *nsg_236 = buffer.data(nsg + 236);
    const auto *nsg_237 = buffer.data(nsg + 237);
    const auto *nsg_238 = buffer.data(nsg + 238);
    const auto *nsg_239 = buffer.data(nsg + 239);
    const auto *nsg_240 = buffer.data(nsg + 240);
    const auto *nsg_242 = buffer.data(nsg + 242);
    const auto *nsg_243 = buffer.data(nsg + 243);
    const auto *nsg_245 = buffer.data(nsg + 245);
    const auto *nsg_249 = buffer.data(nsg + 249);
    const auto *nsg_250 = buffer.data(nsg + 250);
    const auto *nsg_251 = buffer.data(nsg + 251);
    const auto *nsg_252 = buffer.data(nsg + 252);
    const auto *nsg_253 = buffer.data(nsg + 253);
    const auto *nsg_254 = buffer.data(nsg + 254);
    const auto *nsg_255 = buffer.data(nsg + 255);
    const auto *nsg_257 = buffer.data(nsg + 257);
    const auto *nsg_258 = buffer.data(nsg + 258);
    const auto *nsg_260 = buffer.data(nsg + 260);
    const auto *nsg_261 = buffer.data(nsg + 261);
    const auto *nsg_264 = buffer.data(nsg + 264);
    const auto *nsg_265 = buffer.data(nsg + 265);
    const auto *nsg_266 = buffer.data(nsg + 266);

#pragma omp simd aligned(t_251, t_252, t_253, pc_x, pc_y, pc_z, msg_104, msg_120, msg_180, \
                         nsf0_119, nsf0_120, nsf1_119, nsf1_120, nsg_179, \
                         nsg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * msg_104[k]
                   + f_1 * nsf0_119[k]
                   - f_2 * nsf1_119[k]
                   + f_3 * pc_z[k] * nsg_179[k];

        t_252[k] = f_17 * msg_180[k]
                   + f_1 * nsf0_120[k]
                   - f_2 * nsf1_120[k]
                   + f_3 * pc_x[k] * nsg_180[k];

        t_253[k] = f_10 * msg_120[k]
                   + f_3 * pc_y[k] * nsg_180[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, pc_y, pc_z, msg_105, msg_122, msg_183, \
                         nsf0_123, nsf1_123, nsg_180, nsg_182, \
                         nsg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_10 * msg_105[k]
                   + f_3 * pc_z[k] * nsg_180[k];

        t_255[k] = f_17 * msg_183[k]
                   + f_6 * nsf0_123[k]
                   - f_7 * nsf1_123[k]
                   + f_3 * pc_x[k] * nsg_183[k];

        t_256[k] = f_10 * msg_122[k]
                   + f_3 * pc_y[k] * nsg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, msg_108, msg_185, msg_186, nsf0_125, \
                         nsf0_126, nsf1_125, nsf1_126, nsg_183, nsg_185, \
                         nsg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_17 * msg_185[k]
                   + f_6 * nsf0_125[k]
                   - f_7 * nsf1_125[k]
                   + f_3 * pc_x[k] * nsg_185[k];

        t_258[k] = f_17 * msg_186[k]
                   + f_4 * nsf0_126[k]
                   - f_5 * nsf1_126[k]
                   + f_3 * pc_x[k] * nsg_186[k];

        t_259[k] = f_10 * msg_108[k]
                   + f_3 * pc_z[k] * nsg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, msg_125, msg_189, msg_190, \
                         msg_191, nsf0_129, nsf1_129, nsg_185, nsg_189, nsg_190, \
                         nsg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * msg_125[k]
                   + f_3 * pc_y[k] * nsg_185[k];

        t_261[k] = f_17 * msg_189[k]
                   + f_4 * nsf0_129[k]
                   - f_5 * nsf1_129[k]
                   + f_3 * pc_x[k] * nsg_189[k];

        t_262[k] = f_17 * msg_190[k]
                   + f_3 * pc_x[k] * nsg_190[k];

        t_263[k] = f_17 * msg_191[k]
                   + f_3 * pc_x[k] * nsg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, msg_130, msg_192, msg_193, \
                         msg_194, nsf0_126, nsf1_126, nsg_190, nsg_192, nsg_193, \
                         nsg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * msg_192[k]
                   + f_3 * pc_x[k] * nsg_192[k];

        t_265[k] = f_17 * msg_193[k]
                   + f_3 * pc_x[k] * nsg_193[k];

        t_266[k] = f_17 * msg_194[k]
                   + f_3 * pc_x[k] * nsg_194[k];

        t_267[k] = f_10 * msg_130[k]
                   + f_1 * nsf0_126[k]
                   - f_2 * nsf1_126[k]
                   + f_3 * pc_y[k] * nsg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, msg_115, msg_132, msg_133, nsf0_128, \
                         nsf0_129, nsf1_128, nsf1_129, nsg_190, nsg_192, \
                         nsg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * msg_115[k]
                   + f_3 * pc_z[k] * nsg_190[k];

        t_269[k] = f_10 * msg_132[k]
                   + f_6 * nsf0_128[k]
                   - f_7 * nsf1_128[k]
                   + f_3 * pc_y[k] * nsg_192[k];

        t_270[k] = f_10 * msg_133[k]
                   + f_4 * nsf0_129[k]
                   - f_5 * nsf1_129[k]
                   + f_3 * pc_y[k] * nsg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_y, pc_y, pc_z, msh0_189, msg_119, \
                         msg_134, msg_135, msh1_189, nsf0_129, nsf1_129, nsg_194, \
                         nsg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * msg_134[k]
                   + f_3 * pc_y[k] * nsg_194[k];

        t_272[k] = f_10 * msg_119[k]
                   + f_1 * nsf0_129[k]
                   - f_2 * nsf1_129[k]
                   + f_3 * pc_z[k] * nsg_194[k];

        t_273[k] = pa_y[k] * msh0_189[k]
                   - f_8 * pc_y[k] * msh1_189[k];

        t_274[k] = f_9 * msg_135[k]
                   + f_3 * pc_y[k] * nsg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pa_y, pc_y, pc_z, msh0_192, msh0_194, \
                         msg_120, msg_136, msg_137, msh1_192, msh1_194, nsg_195, \
                         nsg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * msg_120[k]
                   + f_3 * pc_z[k] * nsg_195[k];

        t_276[k] = pa_y[k] * msh0_192[k]
                   + f_10 * msg_136[k]
                   - f_8 * pc_y[k] * msh1_192[k];

        t_277[k] = f_9 * msg_137[k]
                   + f_3 * pc_y[k] * nsg_197[k];

        t_278[k] = pa_y[k] * msh0_194[k]
                   - f_8 * pc_y[k] * msh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_y, pc_y, pc_z, msh0_195, msh0_198, \
                         msg_123, msg_138, msg_140, msh1_195, msh1_198, nsg_198, \
                         nsg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pa_y[k] * msh0_195[k]
                   + f_11 * msg_138[k]
                   - f_8 * pc_y[k] * msh1_195[k];

        t_280[k] = f_11 * msg_123[k]
                   + f_3 * pc_z[k] * nsg_198[k];

        t_281[k] = f_9 * msg_140[k]
                   + f_3 * pc_y[k] * nsg_200[k];

        t_282[k] = pa_y[k] * msh0_198[k]
                   - f_8 * pc_y[k] * msh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, msg_205, msg_206, msg_207, \
                         msg_208, msg_209, nsg_205, nsg_206, nsg_207, nsg_208, \
                         nsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_17 * msg_205[k]
                   + f_3 * pc_x[k] * nsg_205[k];

        t_284[k] = f_17 * msg_206[k]
                   + f_3 * pc_x[k] * nsg_206[k];

        t_285[k] = f_17 * msg_207[k]
                   + f_3 * pc_x[k] * nsg_207[k];

        t_286[k] = f_17 * msg_208[k]
                   + f_3 * pc_x[k] * nsg_208[k];

        t_287[k] = f_17 * msg_209[k]
                   + f_3 * pc_x[k] * nsg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, msg_130, msg_145, msg_147, nsf0_136, \
                         nsf0_138, nsf1_136, nsf1_138, nsg_205, \
                         nsg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * msg_145[k]
                   + f_1 * nsf0_136[k]
                   - f_2 * nsf1_136[k]
                   + f_3 * pc_y[k] * nsg_205[k];

        t_289[k] = f_11 * msg_130[k]
                   + f_3 * pc_z[k] * nsg_205[k];

        t_290[k] = f_9 * msg_147[k]
                   + f_6 * nsf0_138[k]
                   - f_7 * nsf1_138[k]
                   + f_3 * pc_y[k] * nsg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pc_y, msh0_209, msg_148, msg_149, \
                         msh1_209, nsf0_139, nsf1_139, nsg_208, \
                         nsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * msg_148[k]
                   + f_4 * nsf0_139[k]
                   - f_5 * nsf1_139[k]
                   + f_3 * pc_y[k] * nsg_208[k];

        t_292[k] = f_9 * msg_149[k]
                   + f_3 * pc_y[k] * nsg_209[k];

        t_293[k] = pa_y[k] * msh0_209[k]
                   - f_8 * pc_y[k] * msh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pc_x, pc_y, pc_z, msg_135, \
                         msg_210, nsf0_140, nsf1_140, nsg_210, nsg_211, \
                         nsg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_17 * msg_210[k]
                   + f_1 * nsf0_140[k]
                   - f_2 * nsf1_140[k]
                   + f_3 * pc_x[k] * nsg_210[k];

        t_295[k] = f_3 * pc_y[k] * nsg_210[k];

        t_296[k] = f_18 * msg_135[k]
                   + f_3 * pc_z[k] * nsg_210[k];

        t_297[k] = f_4 * nsf0_140[k]
                   - f_5 * nsf1_140[k]
                   + f_3 * pc_y[k] * nsg_211[k];

        t_298[k] = f_3 * pc_y[k] * nsg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, msg_215, nsf0_141, nsf0_142, \
                         nsf0_145, nsf1_141, nsf1_142, nsf1_145, nsg_213, nsg_214, \
                         nsg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_17 * msg_215[k]
                   + f_6 * nsf0_145[k]
                   - f_7 * nsf1_145[k]
                   + f_3 * pc_x[k] * nsg_215[k];

        t_300[k] = f_6 * nsf0_141[k]
                   - f_7 * nsf1_141[k]
                   + f_3 * pc_y[k] * nsg_213[k];

        t_301[k] = f_4 * nsf0_142[k]
                   - f_5 * nsf1_142[k]
                   + f_3 * pc_y[k] * nsg_214[k];

        t_302[k] = f_3 * pc_y[k] * nsg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, msg_219, msg_220, msg_221, msg_222, \
                         nsf0_149, nsf1_149, nsg_219, nsg_220, nsg_221, \
                         nsg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_17 * msg_219[k]
                   + f_4 * nsf0_149[k]
                   - f_5 * nsf1_149[k]
                   + f_3 * pc_x[k] * nsg_219[k];

        t_304[k] = f_17 * msg_220[k]
                   + f_3 * pc_x[k] * nsg_220[k];

        t_305[k] = f_17 * msg_221[k]
                   + f_3 * pc_x[k] * nsg_221[k];

        t_306[k] = f_17 * msg_222[k]
                   + f_3 * pc_x[k] * nsg_222[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_x, pc_y, msg_224, nsf0_146, nsf0_147, \
                         nsf1_146, nsf1_147, nsg_219, nsg_220, nsg_221, \
                         nsg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_y[k] * nsg_219[k];

        t_308[k] = f_17 * msg_224[k]
                   + f_3 * pc_x[k] * nsg_224[k];

        t_309[k] = f_1 * nsf0_146[k]
                   - f_2 * nsf1_146[k]
                   + f_3 * pc_y[k] * nsg_220[k];

        t_310[k] = f_13 * nsf0_147[k]
                   - f_14 * nsf1_147[k]
                   + f_3 * pc_y[k] * nsg_221[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, msg_149, nsf0_148, nsf0_149, \
                         nsf1_148, nsf1_149, nsg_222, nsg_223, \
                         nsg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_6 * nsf0_148[k]
                   - f_7 * nsf1_148[k]
                   + f_3 * pc_y[k] * nsg_222[k];

        t_312[k] = f_4 * nsf0_149[k]
                   - f_5 * nsf1_149[k]
                   + f_3 * pc_y[k] * nsg_223[k];

        t_313[k] = f_3 * pc_y[k] * nsg_224[k];

        t_314[k] = f_18 * msg_149[k]
                   + f_1 * nsf0_149[k]
                   - f_2 * nsf1_149[k]
                   + f_3 * pc_z[k] * nsg_224[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, msg_150, msg_225, \
                         msg_228, nsf0_150, nsf0_153, nsf1_150, nsf1_153, nsg_225, \
                         nsg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_19 * msg_225[k]
                   + f_1 * nsf0_150[k]
                   - f_2 * nsf1_150[k]
                   + f_3 * pc_x[k] * nsg_225[k];

        t_316[k] = f_19 * msg_150[k]
                   + f_3 * pc_y[k] * nsg_225[k];

        t_317[k] = f_3 * pc_z[k] * nsg_225[k];

        t_318[k] = f_19 * msg_228[k]
                   + f_6 * nsf0_153[k]
                   - f_7 * nsf1_153[k]
                   + f_3 * pc_x[k] * nsg_228[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pc_x, pc_z, msg_231, nsf0_150, nsf0_156, \
                         nsf1_150, nsf1_156, nsg_226, nsg_227, nsg_228, \
                         nsg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_3 * pc_z[k] * nsg_226[k];

        t_320[k] = f_4 * nsf0_150[k]
                   - f_5 * nsf1_150[k]
                   + f_3 * pc_z[k] * nsg_227[k];

        t_321[k] = f_19 * msg_231[k]
                   + f_4 * nsf0_156[k]
                   - f_5 * nsf1_156[k]
                   + f_3 * pc_x[k] * nsg_231[k];

        t_322[k] = f_3 * pc_z[k] * nsg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, pc_z, msg_155, msg_235, \
                         nsf0_152, nsf1_152, nsg_230, nsg_231, \
                         nsg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_19 * msg_155[k]
                   + f_3 * pc_y[k] * nsg_230[k];

        t_324[k] = f_6 * nsf0_152[k]
                   - f_7 * nsf1_152[k]
                   + f_3 * pc_z[k] * nsg_230[k];

        t_325[k] = f_19 * msg_235[k]
                   + f_3 * pc_x[k] * nsg_235[k];

        t_326[k] = f_3 * pc_z[k] * nsg_231[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, msg_160, msg_237, msg_238, \
                         msg_239, nsf0_156, nsf1_156, nsg_235, nsg_237, nsg_238, \
                         nsg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_19 * msg_237[k]
                   + f_3 * pc_x[k] * nsg_237[k];

        t_328[k] = f_19 * msg_238[k]
                   + f_3 * pc_x[k] * nsg_238[k];

        t_329[k] = f_19 * msg_239[k]
                   + f_3 * pc_x[k] * nsg_239[k];

        t_330[k] = f_19 * msg_160[k]
                   + f_1 * nsf0_156[k]
                   - f_2 * nsf1_156[k]
                   + f_3 * pc_y[k] * nsg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pc_y, pc_z, msg_164, nsf0_156, nsf0_157, \
                         nsf1_156, nsf1_157, nsg_235, nsg_236, nsg_237, \
                         nsg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * nsg_235[k];

        t_332[k] = f_4 * nsf0_156[k]
                   - f_5 * nsf1_156[k]
                   + f_3 * pc_z[k] * nsg_236[k];

        t_333[k] = f_6 * nsf0_157[k]
                   - f_7 * nsf1_157[k]
                   + f_3 * pc_z[k] * nsg_237[k];

        t_334[k] = f_19 * msg_164[k]
                   + f_3 * pc_y[k] * nsg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_z, pc_y, pc_z, msh0_210, msg_150, \
                         msg_165, msh1_210, nsf0_159, nsf1_159, nsg_239, \
                         nsg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_1 * nsf0_159[k]
                   - f_2 * nsf1_159[k]
                   + f_3 * pc_z[k] * nsg_239[k];

        t_336[k] = pa_z[k] * msh0_210[k]
                   - f_8 * pc_z[k] * msh1_210[k];

        t_337[k] = f_18 * msg_165[k]
                   + f_3 * pc_y[k] * nsg_240[k];

        t_338[k] = f_9 * msg_150[k]
                   + f_3 * pc_z[k] * nsg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_z, pc_x, pc_y, pc_z, msh0_213, msg_167, \
                         msg_245, msh1_213, nsf0_165, nsf1_165, nsg_242, \
                         nsg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * msh0_213[k]
                   - f_8 * pc_z[k] * msh1_213[k];

        t_340[k] = f_18 * msg_167[k]
                   + f_3 * pc_y[k] * nsg_242[k];

        t_341[k] = f_19 * msg_245[k]
                   + f_6 * nsf0_165[k]
                   - f_7 * nsf1_165[k]
                   + f_3 * pc_x[k] * nsg_245[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_z, pc_y, pc_z, msh0_216, msg_153, msg_170, \
                         msh1_216, nsg_243, nsg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * msh0_216[k]
                   - f_8 * pc_z[k] * msh1_216[k];

        t_343[k] = f_9 * msg_153[k]
                   + f_3 * pc_z[k] * nsg_243[k];

        t_344[k] = f_18 * msg_170[k]
                   + f_3 * pc_y[k] * nsg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pc_x, msg_249, msg_250, msg_251, msg_252, \
                         nsf0_169, nsf1_169, nsg_249, nsg_250, nsg_251, \
                         nsg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_19 * msg_249[k]
                   + f_4 * nsf0_169[k]
                   - f_5 * nsf1_169[k]
                   + f_3 * pc_x[k] * nsg_249[k];

        t_346[k] = f_19 * msg_250[k]
                   + f_3 * pc_x[k] * nsg_250[k];

        t_347[k] = f_19 * msg_251[k]
                   + f_3 * pc_x[k] * nsg_251[k];

        t_348[k] = f_19 * msg_252[k]
                   + f_3 * pc_x[k] * nsg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_z, pc_x, pc_z, msh0_225, msg_160, \
                         msg_253, msg_254, msh1_225, nsg_250, nsg_253, \
                         nsg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_19 * msg_253[k]
                   + f_3 * pc_x[k] * nsg_253[k];

        t_350[k] = f_19 * msg_254[k]
                   + f_3 * pc_x[k] * nsg_254[k];

        t_351[k] = pa_z[k] * msh0_225[k]
                   - f_8 * pc_z[k] * msh1_225[k];

        t_352[k] = f_9 * msg_160[k]
                   + f_3 * pc_z[k] * nsg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, msg_177, msg_178, msg_179, nsf0_168, \
                         nsf0_169, nsf1_168, nsf1_169, nsg_252, nsg_253, \
                         nsg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_18 * msg_177[k]
                   + f_6 * nsf0_168[k]
                   - f_7 * nsf1_168[k]
                   + f_3 * pc_y[k] * nsg_252[k];

        t_354[k] = f_18 * msg_178[k]
                   + f_4 * nsf0_169[k]
                   - f_5 * nsf1_169[k]
                   + f_3 * pc_y[k] * nsg_253[k];

        t_355[k] = f_18 * msg_179[k]
                   + f_3 * pc_y[k] * nsg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, msg_164, msg_180, msg_255, \
                         nsf0_169, nsf0_170, nsf1_169, nsf1_170, nsg_254, \
                         nsg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_9 * msg_164[k]
                   + f_1 * nsf0_169[k]
                   - f_2 * nsf1_169[k]
                   + f_3 * pc_z[k] * nsg_254[k];

        t_357[k] = f_19 * msg_255[k]
                   + f_1 * nsf0_170[k]
                   - f_2 * nsf1_170[k]
                   + f_3 * pc_x[k] * nsg_255[k];

        t_358[k] = f_11 * msg_180[k]
                   + f_3 * pc_y[k] * nsg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_x, pc_y, pc_z, msg_165, msg_182, msg_258, \
                         nsf0_173, nsf1_173, nsg_255, nsg_257, \
                         nsg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_10 * msg_165[k]
                   + f_3 * pc_z[k] * nsg_255[k];

        t_360[k] = f_19 * msg_258[k]
                   + f_6 * nsf0_173[k]
                   - f_7 * nsf1_173[k]
                   + f_3 * pc_x[k] * nsg_258[k];

        t_361[k] = f_11 * msg_182[k]
                   + f_3 * pc_y[k] * nsg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, msg_168, msg_260, msg_261, nsf0_175, \
                         nsf0_176, nsf1_175, nsf1_176, nsg_258, nsg_260, \
                         nsg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_19 * msg_260[k]
                   + f_6 * nsf0_175[k]
                   - f_7 * nsf1_175[k]
                   + f_3 * pc_x[k] * nsg_260[k];

        t_363[k] = f_19 * msg_261[k]
                   + f_4 * nsf0_176[k]
                   - f_5 * nsf1_176[k]
                   + f_3 * pc_x[k] * nsg_261[k];

        t_364[k] = f_10 * msg_168[k]
                   + f_3 * pc_z[k] * nsg_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, msg_185, msg_264, msg_265, \
                         msg_266, nsf0_179, nsf1_179, nsg_260, nsg_264, nsg_265, \
                         nsg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * msg_185[k]
                   + f_3 * pc_y[k] * nsg_260[k];

        t_366[k] = f_19 * msg_264[k]
                   + f_4 * nsf0_179[k]
                   - f_5 * nsf1_179[k]
                   + f_3 * pc_x[k] * nsg_264[k];

        t_367[k] = f_19 * msg_265[k]
                   + f_3 * pc_x[k] * nsg_265[k];

        t_368[k] = f_19 * msg_266[k]
                   + f_3 * pc_x[k] * nsg_266[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);
    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);
    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_294 = buffer.data(msh0 + 294);
    const auto *msh0_297 = buffer.data(msh0 + 297);
    const auto *msh0_299 = buffer.data(msh0 + 299);
    const auto *msh0_300 = buffer.data(msh0 + 300);
    const auto *msh0_303 = buffer.data(msh0 + 303);
    const auto *msh0_314 = buffer.data(msh0 + 314);
    const auto *msh0_315 = buffer.data(msh0 + 315);
    const auto *msh0_318 = buffer.data(msh0 + 318);
    const auto *msh0_321 = buffer.data(msh0 + 321);
    const auto *msh0_330 = buffer.data(msh0 + 330);

    const auto *msg_175 = buffer.data(msg + 175);
    const auto *msg_179 = buffer.data(msg + 179);
    const auto *msg_180 = buffer.data(msg + 180);
    const auto *msg_183 = buffer.data(msg + 183);
    const auto *msg_190 = buffer.data(msg + 190);
    const auto *msg_192 = buffer.data(msg + 192);
    const auto *msg_193 = buffer.data(msg + 193);
    const auto *msg_194 = buffer.data(msg + 194);
    const auto *msg_195 = buffer.data(msg + 195);
    const auto *msg_197 = buffer.data(msg + 197);
    const auto *msg_198 = buffer.data(msg + 198);
    const auto *msg_200 = buffer.data(msg + 200);
    const auto *msg_205 = buffer.data(msg + 205);
    const auto *msg_207 = buffer.data(msg + 207);
    const auto *msg_208 = buffer.data(msg + 208);
    const auto *msg_209 = buffer.data(msg + 209);
    const auto *msg_210 = buffer.data(msg + 210);
    const auto *msg_211 = buffer.data(msg + 211);
    const auto *msg_212 = buffer.data(msg + 212);
    const auto *msg_213 = buffer.data(msg + 213);
    const auto *msg_215 = buffer.data(msg + 215);
    const auto *msg_220 = buffer.data(msg + 220);
    const auto *msg_222 = buffer.data(msg + 222);
    const auto *msg_223 = buffer.data(msg + 223);
    const auto *msg_224 = buffer.data(msg + 224);
    const auto *msg_225 = buffer.data(msg + 225);
    const auto *msg_228 = buffer.data(msg + 228);
    const auto *msg_230 = buffer.data(msg + 230);
    const auto *msg_235 = buffer.data(msg + 235);
    const auto *msg_239 = buffer.data(msg + 239);
    const auto *msg_240 = buffer.data(msg + 240);
    const auto *msg_242 = buffer.data(msg + 242);
    const auto *msg_245 = buffer.data(msg + 245);
    const auto *msg_252 = buffer.data(msg + 252);
    const auto *msg_253 = buffer.data(msg + 253);
    const auto *msg_254 = buffer.data(msg + 254);
    const auto *msg_255 = buffer.data(msg + 255);
    const auto *msg_257 = buffer.data(msg + 257);
    const auto *msg_267 = buffer.data(msg + 267);
    const auto *msg_268 = buffer.data(msg + 268);
    const auto *msg_269 = buffer.data(msg + 269);
    const auto *msg_270 = buffer.data(msg + 270);
    const auto *msg_273 = buffer.data(msg + 273);
    const auto *msg_275 = buffer.data(msg + 275);
    const auto *msg_276 = buffer.data(msg + 276);
    const auto *msg_279 = buffer.data(msg + 279);
    const auto *msg_280 = buffer.data(msg + 280);
    const auto *msg_281 = buffer.data(msg + 281);
    const auto *msg_282 = buffer.data(msg + 282);
    const auto *msg_283 = buffer.data(msg + 283);
    const auto *msg_284 = buffer.data(msg + 284);
    const auto *msg_295 = buffer.data(msg + 295);
    const auto *msg_296 = buffer.data(msg + 296);
    const auto *msg_297 = buffer.data(msg + 297);
    const auto *msg_298 = buffer.data(msg + 298);
    const auto *msg_299 = buffer.data(msg + 299);
    const auto *msg_300 = buffer.data(msg + 300);
    const auto *msg_305 = buffer.data(msg + 305);
    const auto *msg_309 = buffer.data(msg + 309);
    const auto *msg_310 = buffer.data(msg + 310);
    const auto *msg_311 = buffer.data(msg + 311);
    const auto *msg_312 = buffer.data(msg + 312);
    const auto *msg_314 = buffer.data(msg + 314);
    const auto *msg_315 = buffer.data(msg + 315);
    const auto *msg_318 = buffer.data(msg + 318);
    const auto *msg_321 = buffer.data(msg + 321);
    const auto *msg_325 = buffer.data(msg + 325);
    const auto *msg_327 = buffer.data(msg + 327);
    const auto *msg_328 = buffer.data(msg + 328);
    const auto *msg_329 = buffer.data(msg + 329);
    const auto *msg_335 = buffer.data(msg + 335);
    const auto *msg_339 = buffer.data(msg + 339);
    const auto *msg_340 = buffer.data(msg + 340);
    const auto *msg_341 = buffer.data(msg + 341);
    const auto *msg_342 = buffer.data(msg + 342);
    const auto *msg_343 = buffer.data(msg + 343);
    const auto *msg_344 = buffer.data(msg + 344);
    const auto *msg_345 = buffer.data(msg + 345);
    const auto *msg_348 = buffer.data(msg + 348);

    const auto *msh1_294 = buffer.data(msh1 + 294);
    const auto *msh1_297 = buffer.data(msh1 + 297);
    const auto *msh1_299 = buffer.data(msh1 + 299);
    const auto *msh1_300 = buffer.data(msh1 + 300);
    const auto *msh1_303 = buffer.data(msh1 + 303);
    const auto *msh1_314 = buffer.data(msh1 + 314);
    const auto *msh1_315 = buffer.data(msh1 + 315);
    const auto *msh1_318 = buffer.data(msh1 + 318);
    const auto *msh1_321 = buffer.data(msh1 + 321);
    const auto *msh1_330 = buffer.data(msh1 + 330);

    const auto *nsf0_176 = buffer.data(nsf0 + 176);
    const auto *nsf0_178 = buffer.data(nsf0 + 178);
    const auto *nsf0_179 = buffer.data(nsf0 + 179);
    const auto *nsf0_180 = buffer.data(nsf0 + 180);
    const auto *nsf0_183 = buffer.data(nsf0 + 183);
    const auto *nsf0_185 = buffer.data(nsf0 + 185);
    const auto *nsf0_186 = buffer.data(nsf0 + 186);
    const auto *nsf0_188 = buffer.data(nsf0 + 188);
    const auto *nsf0_189 = buffer.data(nsf0 + 189);
    const auto *nsf0_196 = buffer.data(nsf0 + 196);
    const auto *nsf0_198 = buffer.data(nsf0 + 198);
    const auto *nsf0_199 = buffer.data(nsf0 + 199);
    const auto *nsf0_200 = buffer.data(nsf0 + 200);
    const auto *nsf0_201 = buffer.data(nsf0 + 201);
    const auto *nsf0_202 = buffer.data(nsf0 + 202);
    const auto *nsf0_205 = buffer.data(nsf0 + 205);
    const auto *nsf0_206 = buffer.data(nsf0 + 206);
    const auto *nsf0_207 = buffer.data(nsf0 + 207);
    const auto *nsf0_208 = buffer.data(nsf0 + 208);
    const auto *nsf0_209 = buffer.data(nsf0 + 209);
    const auto *nsf0_210 = buffer.data(nsf0 + 210);
    const auto *nsf0_212 = buffer.data(nsf0 + 212);
    const auto *nsf0_213 = buffer.data(nsf0 + 213);
    const auto *nsf0_216 = buffer.data(nsf0 + 216);
    const auto *nsf0_217 = buffer.data(nsf0 + 217);
    const auto *nsf0_219 = buffer.data(nsf0 + 219);
    const auto *nsf0_225 = buffer.data(nsf0 + 225);
    const auto *nsf0_228 = buffer.data(nsf0 + 228);
    const auto *nsf0_229 = buffer.data(nsf0 + 229);
    const auto *nsf0_230 = buffer.data(nsf0 + 230);
    const auto *nsf0_233 = buffer.data(nsf0 + 233);

    const auto *nsf1_176 = buffer.data(nsf1 + 176);
    const auto *nsf1_178 = buffer.data(nsf1 + 178);
    const auto *nsf1_179 = buffer.data(nsf1 + 179);
    const auto *nsf1_180 = buffer.data(nsf1 + 180);
    const auto *nsf1_183 = buffer.data(nsf1 + 183);
    const auto *nsf1_185 = buffer.data(nsf1 + 185);
    const auto *nsf1_186 = buffer.data(nsf1 + 186);
    const auto *nsf1_188 = buffer.data(nsf1 + 188);
    const auto *nsf1_189 = buffer.data(nsf1 + 189);
    const auto *nsf1_196 = buffer.data(nsf1 + 196);
    const auto *nsf1_198 = buffer.data(nsf1 + 198);
    const auto *nsf1_199 = buffer.data(nsf1 + 199);
    const auto *nsf1_200 = buffer.data(nsf1 + 200);
    const auto *nsf1_201 = buffer.data(nsf1 + 201);
    const auto *nsf1_202 = buffer.data(nsf1 + 202);
    const auto *nsf1_205 = buffer.data(nsf1 + 205);
    const auto *nsf1_206 = buffer.data(nsf1 + 206);
    const auto *nsf1_207 = buffer.data(nsf1 + 207);
    const auto *nsf1_208 = buffer.data(nsf1 + 208);
    const auto *nsf1_209 = buffer.data(nsf1 + 209);
    const auto *nsf1_210 = buffer.data(nsf1 + 210);
    const auto *nsf1_212 = buffer.data(nsf1 + 212);
    const auto *nsf1_213 = buffer.data(nsf1 + 213);
    const auto *nsf1_216 = buffer.data(nsf1 + 216);
    const auto *nsf1_217 = buffer.data(nsf1 + 217);
    const auto *nsf1_219 = buffer.data(nsf1 + 219);
    const auto *nsf1_225 = buffer.data(nsf1 + 225);
    const auto *nsf1_228 = buffer.data(nsf1 + 228);
    const auto *nsf1_229 = buffer.data(nsf1 + 229);
    const auto *nsf1_230 = buffer.data(nsf1 + 230);
    const auto *nsf1_233 = buffer.data(nsf1 + 233);

    const auto *nsg_265 = buffer.data(nsg + 265);
    const auto *nsg_267 = buffer.data(nsg + 267);
    const auto *nsg_268 = buffer.data(nsg + 268);
    const auto *nsg_269 = buffer.data(nsg + 269);
    const auto *nsg_270 = buffer.data(nsg + 270);
    const auto *nsg_272 = buffer.data(nsg + 272);
    const auto *nsg_273 = buffer.data(nsg + 273);
    const auto *nsg_275 = buffer.data(nsg + 275);
    const auto *nsg_276 = buffer.data(nsg + 276);
    const auto *nsg_279 = buffer.data(nsg + 279);
    const auto *nsg_280 = buffer.data(nsg + 280);
    const auto *nsg_281 = buffer.data(nsg + 281);
    const auto *nsg_282 = buffer.data(nsg + 282);
    const auto *nsg_283 = buffer.data(nsg + 283);
    const auto *nsg_284 = buffer.data(nsg + 284);
    const auto *nsg_285 = buffer.data(nsg + 285);
    const auto *nsg_287 = buffer.data(nsg + 287);
    const auto *nsg_288 = buffer.data(nsg + 288);
    const auto *nsg_290 = buffer.data(nsg + 290);
    const auto *nsg_295 = buffer.data(nsg + 295);
    const auto *nsg_296 = buffer.data(nsg + 296);
    const auto *nsg_297 = buffer.data(nsg + 297);
    const auto *nsg_298 = buffer.data(nsg + 298);
    const auto *nsg_299 = buffer.data(nsg + 299);
    const auto *nsg_300 = buffer.data(nsg + 300);
    const auto *nsg_301 = buffer.data(nsg + 301);
    const auto *nsg_302 = buffer.data(nsg + 302);
    const auto *nsg_303 = buffer.data(nsg + 303);
    const auto *nsg_304 = buffer.data(nsg + 304);
    const auto *nsg_305 = buffer.data(nsg + 305);
    const auto *nsg_309 = buffer.data(nsg + 309);
    const auto *nsg_310 = buffer.data(nsg + 310);
    const auto *nsg_311 = buffer.data(nsg + 311);
    const auto *nsg_312 = buffer.data(nsg + 312);
    const auto *nsg_313 = buffer.data(nsg + 313);
    const auto *nsg_314 = buffer.data(nsg + 314);
    const auto *nsg_315 = buffer.data(nsg + 315);
    const auto *nsg_316 = buffer.data(nsg + 316);
    const auto *nsg_317 = buffer.data(nsg + 317);
    const auto *nsg_318 = buffer.data(nsg + 318);
    const auto *nsg_320 = buffer.data(nsg + 320);
    const auto *nsg_321 = buffer.data(nsg + 321);
    const auto *nsg_325 = buffer.data(nsg + 325);
    const auto *nsg_326 = buffer.data(nsg + 326);
    const auto *nsg_327 = buffer.data(nsg + 327);
    const auto *nsg_328 = buffer.data(nsg + 328);
    const auto *nsg_329 = buffer.data(nsg + 329);
    const auto *nsg_330 = buffer.data(nsg + 330);
    const auto *nsg_332 = buffer.data(nsg + 332);
    const auto *nsg_333 = buffer.data(nsg + 333);
    const auto *nsg_335 = buffer.data(nsg + 335);
    const auto *nsg_339 = buffer.data(nsg + 339);
    const auto *nsg_340 = buffer.data(nsg + 340);
    const auto *nsg_341 = buffer.data(nsg + 341);
    const auto *nsg_342 = buffer.data(nsg + 342);
    const auto *nsg_343 = buffer.data(nsg + 343);
    const auto *nsg_344 = buffer.data(nsg + 344);
    const auto *nsg_345 = buffer.data(nsg + 345);
    const auto *nsg_347 = buffer.data(nsg + 347);
    const auto *nsg_348 = buffer.data(nsg + 348);

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, msg_190, msg_267, msg_268, \
                         msg_269, nsf0_176, nsf1_176, nsg_265, nsg_267, nsg_268, \
                         nsg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_19 * msg_267[k]
                   + f_3 * pc_x[k] * nsg_267[k];

        t_370[k] = f_19 * msg_268[k]
                   + f_3 * pc_x[k] * nsg_268[k];

        t_371[k] = f_19 * msg_269[k]
                   + f_3 * pc_x[k] * nsg_269[k];

        t_372[k] = f_11 * msg_190[k]
                   + f_1 * nsf0_176[k]
                   - f_2 * nsf1_176[k]
                   + f_3 * pc_y[k] * nsg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, msg_175, msg_192, msg_193, nsf0_178, \
                         nsf0_179, nsf1_178, nsf1_179, nsg_265, nsg_267, \
                         nsg_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * msg_175[k]
                   + f_3 * pc_z[k] * nsg_265[k];

        t_374[k] = f_11 * msg_192[k]
                   + f_6 * nsf0_178[k]
                   - f_7 * nsf1_178[k]
                   + f_3 * pc_y[k] * nsg_267[k];

        t_375[k] = f_11 * msg_193[k]
                   + f_4 * nsf0_179[k]
                   - f_5 * nsf1_179[k]
                   + f_3 * pc_y[k] * nsg_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, msg_179, msg_194, msg_270, \
                         nsf0_179, nsf0_180, nsf1_179, nsf1_180, nsg_269, \
                         nsg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * msg_194[k]
                   + f_3 * pc_y[k] * nsg_269[k];

        t_377[k] = f_10 * msg_179[k]
                   + f_1 * nsf0_179[k]
                   - f_2 * nsf1_179[k]
                   + f_3 * pc_z[k] * nsg_269[k];

        t_378[k] = f_19 * msg_270[k]
                   + f_1 * nsf0_180[k]
                   - f_2 * nsf1_180[k]
                   + f_3 * pc_x[k] * nsg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, msg_180, msg_195, \
                         msg_197, msg_273, nsf0_183, nsf1_183, nsg_270, nsg_272, \
                         nsg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * msg_195[k]
                   + f_3 * pc_y[k] * nsg_270[k];

        t_380[k] = f_11 * msg_180[k]
                   + f_3 * pc_z[k] * nsg_270[k];

        t_381[k] = f_19 * msg_273[k]
                   + f_6 * nsf0_183[k]
                   - f_7 * nsf1_183[k]
                   + f_3 * pc_x[k] * nsg_273[k];

        t_382[k] = f_10 * msg_197[k]
                   + f_3 * pc_y[k] * nsg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, msg_183, msg_275, msg_276, nsf0_185, \
                         nsf0_186, nsf1_185, nsf1_186, nsg_273, nsg_275, \
                         nsg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_19 * msg_275[k]
                   + f_6 * nsf0_185[k]
                   - f_7 * nsf1_185[k]
                   + f_3 * pc_x[k] * nsg_275[k];

        t_384[k] = f_19 * msg_276[k]
                   + f_4 * nsf0_186[k]
                   - f_5 * nsf1_186[k]
                   + f_3 * pc_x[k] * nsg_276[k];

        t_385[k] = f_11 * msg_183[k]
                   + f_3 * pc_z[k] * nsg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, msg_200, msg_279, msg_280, \
                         msg_281, nsf0_189, nsf1_189, nsg_275, nsg_279, nsg_280, \
                         nsg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * msg_200[k]
                   + f_3 * pc_y[k] * nsg_275[k];

        t_387[k] = f_19 * msg_279[k]
                   + f_4 * nsf0_189[k]
                   - f_5 * nsf1_189[k]
                   + f_3 * pc_x[k] * nsg_279[k];

        t_388[k] = f_19 * msg_280[k]
                   + f_3 * pc_x[k] * nsg_280[k];

        t_389[k] = f_19 * msg_281[k]
                   + f_3 * pc_x[k] * nsg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, msg_205, msg_282, msg_283, \
                         msg_284, nsf0_186, nsf1_186, nsg_280, nsg_282, nsg_283, \
                         nsg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_19 * msg_282[k]
                   + f_3 * pc_x[k] * nsg_282[k];

        t_391[k] = f_19 * msg_283[k]
                   + f_3 * pc_x[k] * nsg_283[k];

        t_392[k] = f_19 * msg_284[k]
                   + f_3 * pc_x[k] * nsg_284[k];

        t_393[k] = f_10 * msg_205[k]
                   + f_1 * nsf0_186[k]
                   - f_2 * nsf1_186[k]
                   + f_3 * pc_y[k] * nsg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, msg_190, msg_207, msg_208, nsf0_188, \
                         nsf0_189, nsf1_188, nsf1_189, nsg_280, nsg_282, \
                         nsg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * msg_190[k]
                   + f_3 * pc_z[k] * nsg_280[k];

        t_395[k] = f_10 * msg_207[k]
                   + f_6 * nsf0_188[k]
                   - f_7 * nsf1_188[k]
                   + f_3 * pc_y[k] * nsg_282[k];

        t_396[k] = f_10 * msg_208[k]
                   + f_4 * nsf0_189[k]
                   - f_5 * nsf1_189[k]
                   + f_3 * pc_y[k] * nsg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pc_y, pc_z, msh0_294, msg_194, \
                         msg_209, msg_210, msh1_294, nsf0_189, nsf1_189, nsg_284, \
                         nsg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * msg_209[k]
                   + f_3 * pc_y[k] * nsg_284[k];

        t_398[k] = f_11 * msg_194[k]
                   + f_1 * nsf0_189[k]
                   - f_2 * nsf1_189[k]
                   + f_3 * pc_z[k] * nsg_284[k];

        t_399[k] = pa_y[k] * msh0_294[k]
                   - f_8 * pc_y[k] * msh1_294[k];

        t_400[k] = f_9 * msg_210[k]
                   + f_3 * pc_y[k] * nsg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pc_y, pc_z, msh0_297, msh0_299, \
                         msg_195, msg_211, msg_212, msh1_297, msh1_299, nsg_285, \
                         nsg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_18 * msg_195[k]
                   + f_3 * pc_z[k] * nsg_285[k];

        t_402[k] = pa_y[k] * msh0_297[k]
                   + f_10 * msg_211[k]
                   - f_8 * pc_y[k] * msh1_297[k];

        t_403[k] = f_9 * msg_212[k]
                   + f_3 * pc_y[k] * nsg_287[k];

        t_404[k] = pa_y[k] * msh0_299[k]
                   - f_8 * pc_y[k] * msh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_y, pc_y, pc_z, msh0_300, msh0_303, \
                         msg_198, msg_213, msg_215, msh1_300, msh1_303, nsg_288, \
                         nsg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_y[k] * msh0_300[k]
                   + f_11 * msg_213[k]
                   - f_8 * pc_y[k] * msh1_300[k];

        t_406[k] = f_18 * msg_198[k]
                   + f_3 * pc_z[k] * nsg_288[k];

        t_407[k] = f_9 * msg_215[k]
                   + f_3 * pc_y[k] * nsg_290[k];

        t_408[k] = pa_y[k] * msh0_303[k]
                   - f_8 * pc_y[k] * msh1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, msg_295, msg_296, msg_297, \
                         msg_298, msg_299, nsg_295, nsg_296, nsg_297, nsg_298, \
                         nsg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_19 * msg_295[k]
                   + f_3 * pc_x[k] * nsg_295[k];

        t_410[k] = f_19 * msg_296[k]
                   + f_3 * pc_x[k] * nsg_296[k];

        t_411[k] = f_19 * msg_297[k]
                   + f_3 * pc_x[k] * nsg_297[k];

        t_412[k] = f_19 * msg_298[k]
                   + f_3 * pc_x[k] * nsg_298[k];

        t_413[k] = f_19 * msg_299[k]
                   + f_3 * pc_x[k] * nsg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, msg_205, msg_220, msg_222, nsf0_196, \
                         nsf0_198, nsf1_196, nsf1_198, nsg_295, \
                         nsg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * msg_220[k]
                   + f_1 * nsf0_196[k]
                   - f_2 * nsf1_196[k]
                   + f_3 * pc_y[k] * nsg_295[k];

        t_415[k] = f_18 * msg_205[k]
                   + f_3 * pc_z[k] * nsg_295[k];

        t_416[k] = f_9 * msg_222[k]
                   + f_6 * nsf0_198[k]
                   - f_7 * nsf1_198[k]
                   + f_3 * pc_y[k] * nsg_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_y, pc_y, msh0_314, msg_223, msg_224, \
                         msh1_314, nsf0_199, nsf1_199, nsg_298, \
                         nsg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * msg_223[k]
                   + f_4 * nsf0_199[k]
                   - f_5 * nsf1_199[k]
                   + f_3 * pc_y[k] * nsg_298[k];

        t_418[k] = f_9 * msg_224[k]
                   + f_3 * pc_y[k] * nsg_299[k];

        t_419[k] = pa_y[k] * msh0_314[k]
                   - f_8 * pc_y[k] * msh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, msg_210, \
                         msg_300, nsf0_200, nsf1_200, nsg_300, nsg_301, \
                         nsg_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_19 * msg_300[k]
                   + f_1 * nsf0_200[k]
                   - f_2 * nsf1_200[k]
                   + f_3 * pc_x[k] * nsg_300[k];

        t_421[k] = f_3 * pc_y[k] * nsg_300[k];

        t_422[k] = f_19 * msg_210[k]
                   + f_3 * pc_z[k] * nsg_300[k];

        t_423[k] = f_4 * nsf0_200[k]
                   - f_5 * nsf1_200[k]
                   + f_3 * pc_y[k] * nsg_301[k];

        t_424[k] = f_3 * pc_y[k] * nsg_302[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, msg_305, nsf0_201, nsf0_202, \
                         nsf0_205, nsf1_201, nsf1_202, nsf1_205, nsg_303, nsg_304, \
                         nsg_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_19 * msg_305[k]
                   + f_6 * nsf0_205[k]
                   - f_7 * nsf1_205[k]
                   + f_3 * pc_x[k] * nsg_305[k];

        t_426[k] = f_6 * nsf0_201[k]
                   - f_7 * nsf1_201[k]
                   + f_3 * pc_y[k] * nsg_303[k];

        t_427[k] = f_4 * nsf0_202[k]
                   - f_5 * nsf1_202[k]
                   + f_3 * pc_y[k] * nsg_304[k];

        t_428[k] = f_3 * pc_y[k] * nsg_305[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pc_x, msg_309, msg_310, msg_311, msg_312, \
                         nsf0_209, nsf1_209, nsg_309, nsg_310, nsg_311, \
                         nsg_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_19 * msg_309[k]
                   + f_4 * nsf0_209[k]
                   - f_5 * nsf1_209[k]
                   + f_3 * pc_x[k] * nsg_309[k];

        t_430[k] = f_19 * msg_310[k]
                   + f_3 * pc_x[k] * nsg_310[k];

        t_431[k] = f_19 * msg_311[k]
                   + f_3 * pc_x[k] * nsg_311[k];

        t_432[k] = f_19 * msg_312[k]
                   + f_3 * pc_x[k] * nsg_312[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, msg_314, nsf0_206, nsf0_207, \
                         nsf1_206, nsf1_207, nsg_309, nsg_310, nsg_311, \
                         nsg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_3 * pc_y[k] * nsg_309[k];

        t_434[k] = f_19 * msg_314[k]
                   + f_3 * pc_x[k] * nsg_314[k];

        t_435[k] = f_1 * nsf0_206[k]
                   - f_2 * nsf1_206[k]
                   + f_3 * pc_y[k] * nsg_310[k];

        t_436[k] = f_13 * nsf0_207[k]
                   - f_14 * nsf1_207[k]
                   + f_3 * pc_y[k] * nsg_311[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_y, pc_z, msg_224, nsf0_208, nsf0_209, \
                         nsf1_208, nsf1_209, nsg_312, nsg_313, \
                         nsg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_6 * nsf0_208[k]
                   - f_7 * nsf1_208[k]
                   + f_3 * pc_y[k] * nsg_312[k];

        t_438[k] = f_4 * nsf0_209[k]
                   - f_5 * nsf1_209[k]
                   + f_3 * pc_y[k] * nsg_313[k];

        t_439[k] = f_3 * pc_y[k] * nsg_314[k];

        t_440[k] = f_19 * msg_224[k]
                   + f_1 * nsf0_209[k]
                   - f_2 * nsf1_209[k]
                   + f_3 * pc_z[k] * nsg_314[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_x, pc_y, pc_z, msg_225, msg_315, \
                         msg_318, nsf0_210, nsf0_213, nsf1_210, nsf1_213, nsg_315, \
                         nsg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_18 * msg_315[k]
                   + f_1 * nsf0_210[k]
                   - f_2 * nsf1_210[k]
                   + f_3 * pc_x[k] * nsg_315[k];

        t_442[k] = f_17 * msg_225[k]
                   + f_3 * pc_y[k] * nsg_315[k];

        t_443[k] = f_3 * pc_z[k] * nsg_315[k];

        t_444[k] = f_18 * msg_318[k]
                   + f_6 * nsf0_213[k]
                   - f_7 * nsf1_213[k]
                   + f_3 * pc_x[k] * nsg_318[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_x, pc_z, msg_321, nsf0_210, nsf0_216, \
                         nsf1_210, nsf1_216, nsg_316, nsg_317, nsg_318, \
                         nsg_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_3 * pc_z[k] * nsg_316[k];

        t_446[k] = f_4 * nsf0_210[k]
                   - f_5 * nsf1_210[k]
                   + f_3 * pc_z[k] * nsg_317[k];

        t_447[k] = f_18 * msg_321[k]
                   + f_4 * nsf0_216[k]
                   - f_5 * nsf1_216[k]
                   + f_3 * pc_x[k] * nsg_321[k];

        t_448[k] = f_3 * pc_z[k] * nsg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, msg_230, msg_325, \
                         nsf0_212, nsf1_212, nsg_320, nsg_321, \
                         nsg_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_17 * msg_230[k]
                   + f_3 * pc_y[k] * nsg_320[k];

        t_450[k] = f_6 * nsf0_212[k]
                   - f_7 * nsf1_212[k]
                   + f_3 * pc_z[k] * nsg_320[k];

        t_451[k] = f_18 * msg_325[k]
                   + f_3 * pc_x[k] * nsg_325[k];

        t_452[k] = f_3 * pc_z[k] * nsg_321[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, msg_235, msg_327, msg_328, \
                         msg_329, nsf0_216, nsf1_216, nsg_325, nsg_327, nsg_328, \
                         nsg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_18 * msg_327[k]
                   + f_3 * pc_x[k] * nsg_327[k];

        t_454[k] = f_18 * msg_328[k]
                   + f_3 * pc_x[k] * nsg_328[k];

        t_455[k] = f_18 * msg_329[k]
                   + f_3 * pc_x[k] * nsg_329[k];

        t_456[k] = f_17 * msg_235[k]
                   + f_1 * nsf0_216[k]
                   - f_2 * nsf1_216[k]
                   + f_3 * pc_y[k] * nsg_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pc_y, pc_z, msg_239, nsf0_216, nsf0_217, \
                         nsf1_216, nsf1_217, nsg_325, nsg_326, nsg_327, \
                         nsg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * nsg_325[k];

        t_458[k] = f_4 * nsf0_216[k]
                   - f_5 * nsf1_216[k]
                   + f_3 * pc_z[k] * nsg_326[k];

        t_459[k] = f_6 * nsf0_217[k]
                   - f_7 * nsf1_217[k]
                   + f_3 * pc_z[k] * nsg_327[k];

        t_460[k] = f_17 * msg_239[k]
                   + f_3 * pc_y[k] * nsg_329[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pa_z, pc_y, pc_z, msh0_315, msg_225, \
                         msg_240, msh1_315, nsf0_219, nsf1_219, nsg_329, \
                         nsg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_1 * nsf0_219[k]
                   - f_2 * nsf1_219[k]
                   + f_3 * pc_z[k] * nsg_329[k];

        t_462[k] = pa_z[k] * msh0_315[k]
                   - f_8 * pc_z[k] * msh1_315[k];

        t_463[k] = f_19 * msg_240[k]
                   + f_3 * pc_y[k] * nsg_330[k];

        t_464[k] = f_9 * msg_225[k]
                   + f_3 * pc_z[k] * nsg_330[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pa_z, pc_x, pc_y, pc_z, msh0_318, msg_242, \
                         msg_335, msh1_318, nsf0_225, nsf1_225, nsg_332, \
                         nsg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_z[k] * msh0_318[k]
                   - f_8 * pc_z[k] * msh1_318[k];

        t_466[k] = f_19 * msg_242[k]
                   + f_3 * pc_y[k] * nsg_332[k];

        t_467[k] = f_18 * msg_335[k]
                   + f_6 * nsf0_225[k]
                   - f_7 * nsf1_225[k]
                   + f_3 * pc_x[k] * nsg_335[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pa_z, pc_y, pc_z, msh0_321, msg_228, msg_245, \
                         msh1_321, nsg_333, nsg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = pa_z[k] * msh0_321[k]
                   - f_8 * pc_z[k] * msh1_321[k];

        t_469[k] = f_9 * msg_228[k]
                   + f_3 * pc_z[k] * nsg_333[k];

        t_470[k] = f_19 * msg_245[k]
                   + f_3 * pc_y[k] * nsg_335[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, msg_339, msg_340, msg_341, msg_342, \
                         nsf0_229, nsf1_229, nsg_339, nsg_340, nsg_341, \
                         nsg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_18 * msg_339[k]
                   + f_4 * nsf0_229[k]
                   - f_5 * nsf1_229[k]
                   + f_3 * pc_x[k] * nsg_339[k];

        t_472[k] = f_18 * msg_340[k]
                   + f_3 * pc_x[k] * nsg_340[k];

        t_473[k] = f_18 * msg_341[k]
                   + f_3 * pc_x[k] * nsg_341[k];

        t_474[k] = f_18 * msg_342[k]
                   + f_3 * pc_x[k] * nsg_342[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_z, pc_x, pc_z, msh0_330, msg_235, \
                         msg_343, msg_344, msh1_330, nsg_340, nsg_343, \
                         nsg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_18 * msg_343[k]
                   + f_3 * pc_x[k] * nsg_343[k];

        t_476[k] = f_18 * msg_344[k]
                   + f_3 * pc_x[k] * nsg_344[k];

        t_477[k] = pa_z[k] * msh0_330[k]
                   - f_8 * pc_z[k] * msh1_330[k];

        t_478[k] = f_9 * msg_235[k]
                   + f_3 * pc_z[k] * nsg_340[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_y, msg_252, msg_253, msg_254, nsf0_228, \
                         nsf0_229, nsf1_228, nsf1_229, nsg_342, nsg_343, \
                         nsg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_19 * msg_252[k]
                   + f_6 * nsf0_228[k]
                   - f_7 * nsf1_228[k]
                   + f_3 * pc_y[k] * nsg_342[k];

        t_480[k] = f_19 * msg_253[k]
                   + f_4 * nsf0_229[k]
                   - f_5 * nsf1_229[k]
                   + f_3 * pc_y[k] * nsg_343[k];

        t_481[k] = f_19 * msg_254[k]
                   + f_3 * pc_y[k] * nsg_344[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, msg_239, msg_255, msg_345, \
                         nsf0_229, nsf0_230, nsf1_229, nsf1_230, nsg_344, \
                         nsg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_9 * msg_239[k]
                   + f_1 * nsf0_229[k]
                   - f_2 * nsf1_229[k]
                   + f_3 * pc_z[k] * nsg_344[k];

        t_483[k] = f_18 * msg_345[k]
                   + f_1 * nsf0_230[k]
                   - f_2 * nsf1_230[k]
                   + f_3 * pc_x[k] * nsg_345[k];

        t_484[k] = f_18 * msg_255[k]
                   + f_3 * pc_y[k] * nsg_345[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_y, pc_z, msg_240, msg_257, msg_348, \
                         nsf0_233, nsf1_233, nsg_345, nsg_347, \
                         nsg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_10 * msg_240[k]
                   + f_3 * pc_z[k] * nsg_345[k];

        t_486[k] = f_18 * msg_348[k]
                   + f_6 * nsf0_233[k]
                   - f_7 * nsf1_233[k]
                   + f_3 * pc_x[k] * nsg_348[k];

        t_487[k] = f_18 * msg_257[k]
                   + f_3 * pc_y[k] * nsg_347[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);
    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_420 = buffer.data(msh0 + 420);
    const auto *msh0_423 = buffer.data(msh0 + 423);
    const auto *msh0_425 = buffer.data(msh0 + 425);
    const auto *msh0_426 = buffer.data(msh0 + 426);
    const auto *msh0_429 = buffer.data(msh0 + 429);
    const auto *msh0_440 = buffer.data(msh0 + 440);

    const auto *msg_243 = buffer.data(msg + 243);
    const auto *msg_250 = buffer.data(msg + 250);
    const auto *msg_254 = buffer.data(msg + 254);
    const auto *msg_255 = buffer.data(msg + 255);
    const auto *msg_258 = buffer.data(msg + 258);
    const auto *msg_260 = buffer.data(msg + 260);
    const auto *msg_265 = buffer.data(msg + 265);
    const auto *msg_267 = buffer.data(msg + 267);
    const auto *msg_268 = buffer.data(msg + 268);
    const auto *msg_269 = buffer.data(msg + 269);
    const auto *msg_270 = buffer.data(msg + 270);
    const auto *msg_272 = buffer.data(msg + 272);
    const auto *msg_273 = buffer.data(msg + 273);
    const auto *msg_275 = buffer.data(msg + 275);
    const auto *msg_280 = buffer.data(msg + 280);
    const auto *msg_282 = buffer.data(msg + 282);
    const auto *msg_283 = buffer.data(msg + 283);
    const auto *msg_284 = buffer.data(msg + 284);
    const auto *msg_285 = buffer.data(msg + 285);
    const auto *msg_287 = buffer.data(msg + 287);
    const auto *msg_288 = buffer.data(msg + 288);
    const auto *msg_290 = buffer.data(msg + 290);
    const auto *msg_295 = buffer.data(msg + 295);
    const auto *msg_297 = buffer.data(msg + 297);
    const auto *msg_298 = buffer.data(msg + 298);
    const auto *msg_299 = buffer.data(msg + 299);
    const auto *msg_300 = buffer.data(msg + 300);
    const auto *msg_301 = buffer.data(msg + 301);
    const auto *msg_302 = buffer.data(msg + 302);
    const auto *msg_303 = buffer.data(msg + 303);
    const auto *msg_305 = buffer.data(msg + 305);
    const auto *msg_310 = buffer.data(msg + 310);
    const auto *msg_312 = buffer.data(msg + 312);
    const auto *msg_313 = buffer.data(msg + 313);
    const auto *msg_314 = buffer.data(msg + 314);
    const auto *msg_315 = buffer.data(msg + 315);
    const auto *msg_320 = buffer.data(msg + 320);
    const auto *msg_325 = buffer.data(msg + 325);
    const auto *msg_350 = buffer.data(msg + 350);
    const auto *msg_351 = buffer.data(msg + 351);
    const auto *msg_354 = buffer.data(msg + 354);
    const auto *msg_355 = buffer.data(msg + 355);
    const auto *msg_356 = buffer.data(msg + 356);
    const auto *msg_357 = buffer.data(msg + 357);
    const auto *msg_358 = buffer.data(msg + 358);
    const auto *msg_359 = buffer.data(msg + 359);
    const auto *msg_360 = buffer.data(msg + 360);
    const auto *msg_363 = buffer.data(msg + 363);
    const auto *msg_365 = buffer.data(msg + 365);
    const auto *msg_366 = buffer.data(msg + 366);
    const auto *msg_369 = buffer.data(msg + 369);
    const auto *msg_370 = buffer.data(msg + 370);
    const auto *msg_371 = buffer.data(msg + 371);
    const auto *msg_372 = buffer.data(msg + 372);
    const auto *msg_373 = buffer.data(msg + 373);
    const auto *msg_374 = buffer.data(msg + 374);
    const auto *msg_375 = buffer.data(msg + 375);
    const auto *msg_378 = buffer.data(msg + 378);
    const auto *msg_380 = buffer.data(msg + 380);
    const auto *msg_381 = buffer.data(msg + 381);
    const auto *msg_384 = buffer.data(msg + 384);
    const auto *msg_385 = buffer.data(msg + 385);
    const auto *msg_386 = buffer.data(msg + 386);
    const auto *msg_387 = buffer.data(msg + 387);
    const auto *msg_388 = buffer.data(msg + 388);
    const auto *msg_389 = buffer.data(msg + 389);
    const auto *msg_400 = buffer.data(msg + 400);
    const auto *msg_401 = buffer.data(msg + 401);
    const auto *msg_402 = buffer.data(msg + 402);
    const auto *msg_403 = buffer.data(msg + 403);
    const auto *msg_404 = buffer.data(msg + 404);
    const auto *msg_405 = buffer.data(msg + 405);
    const auto *msg_410 = buffer.data(msg + 410);
    const auto *msg_414 = buffer.data(msg + 414);
    const auto *msg_415 = buffer.data(msg + 415);
    const auto *msg_416 = buffer.data(msg + 416);
    const auto *msg_417 = buffer.data(msg + 417);
    const auto *msg_419 = buffer.data(msg + 419);
    const auto *msg_420 = buffer.data(msg + 420);
    const auto *msg_423 = buffer.data(msg + 423);
    const auto *msg_426 = buffer.data(msg + 426);
    const auto *msg_430 = buffer.data(msg + 430);
    const auto *msg_432 = buffer.data(msg + 432);
    const auto *msg_433 = buffer.data(msg + 433);
    const auto *msg_434 = buffer.data(msg + 434);

    const auto *msh1_420 = buffer.data(msh1 + 420);
    const auto *msh1_423 = buffer.data(msh1 + 423);
    const auto *msh1_425 = buffer.data(msh1 + 425);
    const auto *msh1_426 = buffer.data(msh1 + 426);
    const auto *msh1_429 = buffer.data(msh1 + 429);
    const auto *msh1_440 = buffer.data(msh1 + 440);

    const auto *nsf0_235 = buffer.data(nsf0 + 235);
    const auto *nsf0_236 = buffer.data(nsf0 + 236);
    const auto *nsf0_238 = buffer.data(nsf0 + 238);
    const auto *nsf0_239 = buffer.data(nsf0 + 239);
    const auto *nsf0_240 = buffer.data(nsf0 + 240);
    const auto *nsf0_243 = buffer.data(nsf0 + 243);
    const auto *nsf0_245 = buffer.data(nsf0 + 245);
    const auto *nsf0_246 = buffer.data(nsf0 + 246);
    const auto *nsf0_248 = buffer.data(nsf0 + 248);
    const auto *nsf0_249 = buffer.data(nsf0 + 249);
    const auto *nsf0_250 = buffer.data(nsf0 + 250);
    const auto *nsf0_253 = buffer.data(nsf0 + 253);
    const auto *nsf0_255 = buffer.data(nsf0 + 255);
    const auto *nsf0_256 = buffer.data(nsf0 + 256);
    const auto *nsf0_258 = buffer.data(nsf0 + 258);
    const auto *nsf0_259 = buffer.data(nsf0 + 259);
    const auto *nsf0_266 = buffer.data(nsf0 + 266);
    const auto *nsf0_268 = buffer.data(nsf0 + 268);
    const auto *nsf0_269 = buffer.data(nsf0 + 269);
    const auto *nsf0_270 = buffer.data(nsf0 + 270);
    const auto *nsf0_271 = buffer.data(nsf0 + 271);
    const auto *nsf0_272 = buffer.data(nsf0 + 272);
    const auto *nsf0_275 = buffer.data(nsf0 + 275);
    const auto *nsf0_276 = buffer.data(nsf0 + 276);
    const auto *nsf0_277 = buffer.data(nsf0 + 277);
    const auto *nsf0_278 = buffer.data(nsf0 + 278);
    const auto *nsf0_279 = buffer.data(nsf0 + 279);
    const auto *nsf0_280 = buffer.data(nsf0 + 280);
    const auto *nsf0_282 = buffer.data(nsf0 + 282);
    const auto *nsf0_283 = buffer.data(nsf0 + 283);
    const auto *nsf0_286 = buffer.data(nsf0 + 286);

    const auto *nsf1_235 = buffer.data(nsf1 + 235);
    const auto *nsf1_236 = buffer.data(nsf1 + 236);
    const auto *nsf1_238 = buffer.data(nsf1 + 238);
    const auto *nsf1_239 = buffer.data(nsf1 + 239);
    const auto *nsf1_240 = buffer.data(nsf1 + 240);
    const auto *nsf1_243 = buffer.data(nsf1 + 243);
    const auto *nsf1_245 = buffer.data(nsf1 + 245);
    const auto *nsf1_246 = buffer.data(nsf1 + 246);
    const auto *nsf1_248 = buffer.data(nsf1 + 248);
    const auto *nsf1_249 = buffer.data(nsf1 + 249);
    const auto *nsf1_250 = buffer.data(nsf1 + 250);
    const auto *nsf1_253 = buffer.data(nsf1 + 253);
    const auto *nsf1_255 = buffer.data(nsf1 + 255);
    const auto *nsf1_256 = buffer.data(nsf1 + 256);
    const auto *nsf1_258 = buffer.data(nsf1 + 258);
    const auto *nsf1_259 = buffer.data(nsf1 + 259);
    const auto *nsf1_266 = buffer.data(nsf1 + 266);
    const auto *nsf1_268 = buffer.data(nsf1 + 268);
    const auto *nsf1_269 = buffer.data(nsf1 + 269);
    const auto *nsf1_270 = buffer.data(nsf1 + 270);
    const auto *nsf1_271 = buffer.data(nsf1 + 271);
    const auto *nsf1_272 = buffer.data(nsf1 + 272);
    const auto *nsf1_275 = buffer.data(nsf1 + 275);
    const auto *nsf1_276 = buffer.data(nsf1 + 276);
    const auto *nsf1_277 = buffer.data(nsf1 + 277);
    const auto *nsf1_278 = buffer.data(nsf1 + 278);
    const auto *nsf1_279 = buffer.data(nsf1 + 279);
    const auto *nsf1_280 = buffer.data(nsf1 + 280);
    const auto *nsf1_282 = buffer.data(nsf1 + 282);
    const auto *nsf1_283 = buffer.data(nsf1 + 283);
    const auto *nsf1_286 = buffer.data(nsf1 + 286);

    const auto *nsg_348 = buffer.data(nsg + 348);
    const auto *nsg_350 = buffer.data(nsg + 350);
    const auto *nsg_351 = buffer.data(nsg + 351);
    const auto *nsg_354 = buffer.data(nsg + 354);
    const auto *nsg_355 = buffer.data(nsg + 355);
    const auto *nsg_356 = buffer.data(nsg + 356);
    const auto *nsg_357 = buffer.data(nsg + 357);
    const auto *nsg_358 = buffer.data(nsg + 358);
    const auto *nsg_359 = buffer.data(nsg + 359);
    const auto *nsg_360 = buffer.data(nsg + 360);
    const auto *nsg_362 = buffer.data(nsg + 362);
    const auto *nsg_363 = buffer.data(nsg + 363);
    const auto *nsg_365 = buffer.data(nsg + 365);
    const auto *nsg_366 = buffer.data(nsg + 366);
    const auto *nsg_369 = buffer.data(nsg + 369);
    const auto *nsg_370 = buffer.data(nsg + 370);
    const auto *nsg_371 = buffer.data(nsg + 371);
    const auto *nsg_372 = buffer.data(nsg + 372);
    const auto *nsg_373 = buffer.data(nsg + 373);
    const auto *nsg_374 = buffer.data(nsg + 374);
    const auto *nsg_375 = buffer.data(nsg + 375);
    const auto *nsg_377 = buffer.data(nsg + 377);
    const auto *nsg_378 = buffer.data(nsg + 378);
    const auto *nsg_380 = buffer.data(nsg + 380);
    const auto *nsg_381 = buffer.data(nsg + 381);
    const auto *nsg_384 = buffer.data(nsg + 384);
    const auto *nsg_385 = buffer.data(nsg + 385);
    const auto *nsg_386 = buffer.data(nsg + 386);
    const auto *nsg_387 = buffer.data(nsg + 387);
    const auto *nsg_388 = buffer.data(nsg + 388);
    const auto *nsg_389 = buffer.data(nsg + 389);
    const auto *nsg_390 = buffer.data(nsg + 390);
    const auto *nsg_392 = buffer.data(nsg + 392);
    const auto *nsg_393 = buffer.data(nsg + 393);
    const auto *nsg_395 = buffer.data(nsg + 395);
    const auto *nsg_400 = buffer.data(nsg + 400);
    const auto *nsg_401 = buffer.data(nsg + 401);
    const auto *nsg_402 = buffer.data(nsg + 402);
    const auto *nsg_403 = buffer.data(nsg + 403);
    const auto *nsg_404 = buffer.data(nsg + 404);
    const auto *nsg_405 = buffer.data(nsg + 405);
    const auto *nsg_406 = buffer.data(nsg + 406);
    const auto *nsg_407 = buffer.data(nsg + 407);
    const auto *nsg_408 = buffer.data(nsg + 408);
    const auto *nsg_409 = buffer.data(nsg + 409);
    const auto *nsg_410 = buffer.data(nsg + 410);
    const auto *nsg_414 = buffer.data(nsg + 414);
    const auto *nsg_415 = buffer.data(nsg + 415);
    const auto *nsg_416 = buffer.data(nsg + 416);
    const auto *nsg_417 = buffer.data(nsg + 417);
    const auto *nsg_418 = buffer.data(nsg + 418);
    const auto *nsg_419 = buffer.data(nsg + 419);
    const auto *nsg_420 = buffer.data(nsg + 420);
    const auto *nsg_421 = buffer.data(nsg + 421);
    const auto *nsg_422 = buffer.data(nsg + 422);
    const auto *nsg_423 = buffer.data(nsg + 423);
    const auto *nsg_425 = buffer.data(nsg + 425);
    const auto *nsg_426 = buffer.data(nsg + 426);
    const auto *nsg_430 = buffer.data(nsg + 430);
    const auto *nsg_432 = buffer.data(nsg + 432);
    const auto *nsg_433 = buffer.data(nsg + 433);
    const auto *nsg_434 = buffer.data(nsg + 434);

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_z, msg_243, msg_350, msg_351, nsf0_235, \
                         nsf0_236, nsf1_235, nsf1_236, nsg_348, nsg_350, \
                         nsg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_18 * msg_350[k]
                   + f_6 * nsf0_235[k]
                   - f_7 * nsf1_235[k]
                   + f_3 * pc_x[k] * nsg_350[k];

        t_489[k] = f_18 * msg_351[k]
                   + f_4 * nsf0_236[k]
                   - f_5 * nsf1_236[k]
                   + f_3 * pc_x[k] * nsg_351[k];

        t_490[k] = f_10 * msg_243[k]
                   + f_3 * pc_z[k] * nsg_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, msg_260, msg_354, msg_355, \
                         msg_356, nsf0_239, nsf1_239, nsg_350, nsg_354, nsg_355, \
                         nsg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_18 * msg_260[k]
                   + f_3 * pc_y[k] * nsg_350[k];

        t_492[k] = f_18 * msg_354[k]
                   + f_4 * nsf0_239[k]
                   - f_5 * nsf1_239[k]
                   + f_3 * pc_x[k] * nsg_354[k];

        t_493[k] = f_18 * msg_355[k]
                   + f_3 * pc_x[k] * nsg_355[k];

        t_494[k] = f_18 * msg_356[k]
                   + f_3 * pc_x[k] * nsg_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, msg_265, msg_357, msg_358, \
                         msg_359, nsf0_236, nsf1_236, nsg_355, nsg_357, nsg_358, \
                         nsg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_18 * msg_357[k]
                   + f_3 * pc_x[k] * nsg_357[k];

        t_496[k] = f_18 * msg_358[k]
                   + f_3 * pc_x[k] * nsg_358[k];

        t_497[k] = f_18 * msg_359[k]
                   + f_3 * pc_x[k] * nsg_359[k];

        t_498[k] = f_18 * msg_265[k]
                   + f_1 * nsf0_236[k]
                   - f_2 * nsf1_236[k]
                   + f_3 * pc_y[k] * nsg_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, msg_250, msg_267, msg_268, nsf0_238, \
                         nsf0_239, nsf1_238, nsf1_239, nsg_355, nsg_357, \
                         nsg_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_10 * msg_250[k]
                   + f_3 * pc_z[k] * nsg_355[k];

        t_500[k] = f_18 * msg_267[k]
                   + f_6 * nsf0_238[k]
                   - f_7 * nsf1_238[k]
                   + f_3 * pc_y[k] * nsg_357[k];

        t_501[k] = f_18 * msg_268[k]
                   + f_4 * nsf0_239[k]
                   - f_5 * nsf1_239[k]
                   + f_3 * pc_y[k] * nsg_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, msg_254, msg_269, msg_360, \
                         nsf0_239, nsf0_240, nsf1_239, nsf1_240, nsg_359, \
                         nsg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_18 * msg_269[k]
                   + f_3 * pc_y[k] * nsg_359[k];

        t_503[k] = f_10 * msg_254[k]
                   + f_1 * nsf0_239[k]
                   - f_2 * nsf1_239[k]
                   + f_3 * pc_z[k] * nsg_359[k];

        t_504[k] = f_18 * msg_360[k]
                   + f_1 * nsf0_240[k]
                   - f_2 * nsf1_240[k]
                   + f_3 * pc_x[k] * nsg_360[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, msg_255, msg_270, \
                         msg_272, msg_363, nsf0_243, nsf1_243, nsg_360, nsg_362, \
                         nsg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_11 * msg_270[k]
                   + f_3 * pc_y[k] * nsg_360[k];

        t_506[k] = f_11 * msg_255[k]
                   + f_3 * pc_z[k] * nsg_360[k];

        t_507[k] = f_18 * msg_363[k]
                   + f_6 * nsf0_243[k]
                   - f_7 * nsf1_243[k]
                   + f_3 * pc_x[k] * nsg_363[k];

        t_508[k] = f_11 * msg_272[k]
                   + f_3 * pc_y[k] * nsg_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, msg_258, msg_365, msg_366, nsf0_245, \
                         nsf0_246, nsf1_245, nsf1_246, nsg_363, nsg_365, \
                         nsg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_18 * msg_365[k]
                   + f_6 * nsf0_245[k]
                   - f_7 * nsf1_245[k]
                   + f_3 * pc_x[k] * nsg_365[k];

        t_510[k] = f_18 * msg_366[k]
                   + f_4 * nsf0_246[k]
                   - f_5 * nsf1_246[k]
                   + f_3 * pc_x[k] * nsg_366[k];

        t_511[k] = f_11 * msg_258[k]
                   + f_3 * pc_z[k] * nsg_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_y, msg_275, msg_369, msg_370, \
                         msg_371, nsf0_249, nsf1_249, nsg_365, nsg_369, nsg_370, \
                         nsg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * msg_275[k]
                   + f_3 * pc_y[k] * nsg_365[k];

        t_513[k] = f_18 * msg_369[k]
                   + f_4 * nsf0_249[k]
                   - f_5 * nsf1_249[k]
                   + f_3 * pc_x[k] * nsg_369[k];

        t_514[k] = f_18 * msg_370[k]
                   + f_3 * pc_x[k] * nsg_370[k];

        t_515[k] = f_18 * msg_371[k]
                   + f_3 * pc_x[k] * nsg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, msg_280, msg_372, msg_373, \
                         msg_374, nsf0_246, nsf1_246, nsg_370, nsg_372, nsg_373, \
                         nsg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_18 * msg_372[k]
                   + f_3 * pc_x[k] * nsg_372[k];

        t_517[k] = f_18 * msg_373[k]
                   + f_3 * pc_x[k] * nsg_373[k];

        t_518[k] = f_18 * msg_374[k]
                   + f_3 * pc_x[k] * nsg_374[k];

        t_519[k] = f_11 * msg_280[k]
                   + f_1 * nsf0_246[k]
                   - f_2 * nsf1_246[k]
                   + f_3 * pc_y[k] * nsg_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, msg_265, msg_282, msg_283, nsf0_248, \
                         nsf0_249, nsf1_248, nsf1_249, nsg_370, nsg_372, \
                         nsg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * msg_265[k]
                   + f_3 * pc_z[k] * nsg_370[k];

        t_521[k] = f_11 * msg_282[k]
                   + f_6 * nsf0_248[k]
                   - f_7 * nsf1_248[k]
                   + f_3 * pc_y[k] * nsg_372[k];

        t_522[k] = f_11 * msg_283[k]
                   + f_4 * nsf0_249[k]
                   - f_5 * nsf1_249[k]
                   + f_3 * pc_y[k] * nsg_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, msg_269, msg_284, msg_375, \
                         nsf0_249, nsf0_250, nsf1_249, nsf1_250, nsg_374, \
                         nsg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * msg_284[k]
                   + f_3 * pc_y[k] * nsg_374[k];

        t_524[k] = f_11 * msg_269[k]
                   + f_1 * nsf0_249[k]
                   - f_2 * nsf1_249[k]
                   + f_3 * pc_z[k] * nsg_374[k];

        t_525[k] = f_18 * msg_375[k]
                   + f_1 * nsf0_250[k]
                   - f_2 * nsf1_250[k]
                   + f_3 * pc_x[k] * nsg_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_x, pc_y, pc_z, msg_270, msg_285, \
                         msg_287, msg_378, nsf0_253, nsf1_253, nsg_375, nsg_377, \
                         nsg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * msg_285[k]
                   + f_3 * pc_y[k] * nsg_375[k];

        t_527[k] = f_18 * msg_270[k]
                   + f_3 * pc_z[k] * nsg_375[k];

        t_528[k] = f_18 * msg_378[k]
                   + f_6 * nsf0_253[k]
                   - f_7 * nsf1_253[k]
                   + f_3 * pc_x[k] * nsg_378[k];

        t_529[k] = f_10 * msg_287[k]
                   + f_3 * pc_y[k] * nsg_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_z, msg_273, msg_380, msg_381, nsf0_255, \
                         nsf0_256, nsf1_255, nsf1_256, nsg_378, nsg_380, \
                         nsg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_18 * msg_380[k]
                   + f_6 * nsf0_255[k]
                   - f_7 * nsf1_255[k]
                   + f_3 * pc_x[k] * nsg_380[k];

        t_531[k] = f_18 * msg_381[k]
                   + f_4 * nsf0_256[k]
                   - f_5 * nsf1_256[k]
                   + f_3 * pc_x[k] * nsg_381[k];

        t_532[k] = f_18 * msg_273[k]
                   + f_3 * pc_z[k] * nsg_378[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, msg_290, msg_384, msg_385, \
                         msg_386, nsf0_259, nsf1_259, nsg_380, nsg_384, nsg_385, \
                         nsg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_10 * msg_290[k]
                   + f_3 * pc_y[k] * nsg_380[k];

        t_534[k] = f_18 * msg_384[k]
                   + f_4 * nsf0_259[k]
                   - f_5 * nsf1_259[k]
                   + f_3 * pc_x[k] * nsg_384[k];

        t_535[k] = f_18 * msg_385[k]
                   + f_3 * pc_x[k] * nsg_385[k];

        t_536[k] = f_18 * msg_386[k]
                   + f_3 * pc_x[k] * nsg_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, msg_295, msg_387, msg_388, \
                         msg_389, nsf0_256, nsf1_256, nsg_385, nsg_387, nsg_388, \
                         nsg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_18 * msg_387[k]
                   + f_3 * pc_x[k] * nsg_387[k];

        t_538[k] = f_18 * msg_388[k]
                   + f_3 * pc_x[k] * nsg_388[k];

        t_539[k] = f_18 * msg_389[k]
                   + f_3 * pc_x[k] * nsg_389[k];

        t_540[k] = f_10 * msg_295[k]
                   + f_1 * nsf0_256[k]
                   - f_2 * nsf1_256[k]
                   + f_3 * pc_y[k] * nsg_385[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_y, pc_z, msg_280, msg_297, msg_298, nsf0_258, \
                         nsf0_259, nsf1_258, nsf1_259, nsg_385, nsg_387, \
                         nsg_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_18 * msg_280[k]
                   + f_3 * pc_z[k] * nsg_385[k];

        t_542[k] = f_10 * msg_297[k]
                   + f_6 * nsf0_258[k]
                   - f_7 * nsf1_258[k]
                   + f_3 * pc_y[k] * nsg_387[k];

        t_543[k] = f_10 * msg_298[k]
                   + f_4 * nsf0_259[k]
                   - f_5 * nsf1_259[k]
                   + f_3 * pc_y[k] * nsg_388[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_y, pc_z, msh0_420, msg_284, \
                         msg_299, msg_300, msh1_420, nsf0_259, nsf1_259, nsg_389, \
                         nsg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_10 * msg_299[k]
                   + f_3 * pc_y[k] * nsg_389[k];

        t_545[k] = f_18 * msg_284[k]
                   + f_1 * nsf0_259[k]
                   - f_2 * nsf1_259[k]
                   + f_3 * pc_z[k] * nsg_389[k];

        t_546[k] = pa_y[k] * msh0_420[k]
                   - f_8 * pc_y[k] * msh1_420[k];

        t_547[k] = f_9 * msg_300[k]
                   + f_3 * pc_y[k] * nsg_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pa_y, pc_y, pc_z, msh0_423, msh0_425, \
                         msg_285, msg_301, msg_302, msh1_423, msh1_425, nsg_390, \
                         nsg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_19 * msg_285[k]
                   + f_3 * pc_z[k] * nsg_390[k];

        t_549[k] = pa_y[k] * msh0_423[k]
                   + f_10 * msg_301[k]
                   - f_8 * pc_y[k] * msh1_423[k];

        t_550[k] = f_9 * msg_302[k]
                   + f_3 * pc_y[k] * nsg_392[k];

        t_551[k] = pa_y[k] * msh0_425[k]
                   - f_8 * pc_y[k] * msh1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pa_y, pc_y, pc_z, msh0_426, msh0_429, \
                         msg_288, msg_303, msg_305, msh1_426, msh1_429, nsg_393, \
                         nsg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pa_y[k] * msh0_426[k]
                   + f_11 * msg_303[k]
                   - f_8 * pc_y[k] * msh1_426[k];

        t_553[k] = f_19 * msg_288[k]
                   + f_3 * pc_z[k] * nsg_393[k];

        t_554[k] = f_9 * msg_305[k]
                   + f_3 * pc_y[k] * nsg_395[k];

        t_555[k] = pa_y[k] * msh0_429[k]
                   - f_8 * pc_y[k] * msh1_429[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, pc_x, msg_400, msg_401, msg_402, \
                         msg_403, msg_404, nsg_400, nsg_401, nsg_402, nsg_403, \
                         nsg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_18 * msg_400[k]
                   + f_3 * pc_x[k] * nsg_400[k];

        t_557[k] = f_18 * msg_401[k]
                   + f_3 * pc_x[k] * nsg_401[k];

        t_558[k] = f_18 * msg_402[k]
                   + f_3 * pc_x[k] * nsg_402[k];

        t_559[k] = f_18 * msg_403[k]
                   + f_3 * pc_x[k] * nsg_403[k];

        t_560[k] = f_18 * msg_404[k]
                   + f_3 * pc_x[k] * nsg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_y, pc_z, msg_295, msg_310, msg_312, nsf0_266, \
                         nsf0_268, nsf1_266, nsf1_268, nsg_400, \
                         nsg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_9 * msg_310[k]
                   + f_1 * nsf0_266[k]
                   - f_2 * nsf1_266[k]
                   + f_3 * pc_y[k] * nsg_400[k];

        t_562[k] = f_19 * msg_295[k]
                   + f_3 * pc_z[k] * nsg_400[k];

        t_563[k] = f_9 * msg_312[k]
                   + f_6 * nsf0_268[k]
                   - f_7 * nsf1_268[k]
                   + f_3 * pc_y[k] * nsg_402[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_y, pc_y, msh0_440, msg_313, msg_314, \
                         msh1_440, nsf0_269, nsf1_269, nsg_403, \
                         nsg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_9 * msg_313[k]
                   + f_4 * nsf0_269[k]
                   - f_5 * nsf1_269[k]
                   + f_3 * pc_y[k] * nsg_403[k];

        t_565[k] = f_9 * msg_314[k]
                   + f_3 * pc_y[k] * nsg_404[k];

        t_566[k] = pa_y[k] * msh0_440[k]
                   - f_8 * pc_y[k] * msh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, pc_x, pc_y, pc_z, msg_300, \
                         msg_405, nsf0_270, nsf1_270, nsg_405, nsg_406, \
                         nsg_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_18 * msg_405[k]
                   + f_1 * nsf0_270[k]
                   - f_2 * nsf1_270[k]
                   + f_3 * pc_x[k] * nsg_405[k];

        t_568[k] = f_3 * pc_y[k] * nsg_405[k];

        t_569[k] = f_17 * msg_300[k]
                   + f_3 * pc_z[k] * nsg_405[k];

        t_570[k] = f_4 * nsf0_270[k]
                   - f_5 * nsf1_270[k]
                   + f_3 * pc_y[k] * nsg_406[k];

        t_571[k] = f_3 * pc_y[k] * nsg_407[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, msg_410, nsf0_271, nsf0_272, \
                         nsf0_275, nsf1_271, nsf1_272, nsf1_275, nsg_408, nsg_409, \
                         nsg_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_18 * msg_410[k]
                   + f_6 * nsf0_275[k]
                   - f_7 * nsf1_275[k]
                   + f_3 * pc_x[k] * nsg_410[k];

        t_573[k] = f_6 * nsf0_271[k]
                   - f_7 * nsf1_271[k]
                   + f_3 * pc_y[k] * nsg_408[k];

        t_574[k] = f_4 * nsf0_272[k]
                   - f_5 * nsf1_272[k]
                   + f_3 * pc_y[k] * nsg_409[k];

        t_575[k] = f_3 * pc_y[k] * nsg_410[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pc_x, msg_414, msg_415, msg_416, msg_417, \
                         nsf0_279, nsf1_279, nsg_414, nsg_415, nsg_416, \
                         nsg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_18 * msg_414[k]
                   + f_4 * nsf0_279[k]
                   - f_5 * nsf1_279[k]
                   + f_3 * pc_x[k] * nsg_414[k];

        t_577[k] = f_18 * msg_415[k]
                   + f_3 * pc_x[k] * nsg_415[k];

        t_578[k] = f_18 * msg_416[k]
                   + f_3 * pc_x[k] * nsg_416[k];

        t_579[k] = f_18 * msg_417[k]
                   + f_3 * pc_x[k] * nsg_417[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, msg_419, nsf0_276, nsf0_277, \
                         nsf1_276, nsf1_277, nsg_414, nsg_415, nsg_416, \
                         nsg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_3 * pc_y[k] * nsg_414[k];

        t_581[k] = f_18 * msg_419[k]
                   + f_3 * pc_x[k] * nsg_419[k];

        t_582[k] = f_1 * nsf0_276[k]
                   - f_2 * nsf1_276[k]
                   + f_3 * pc_y[k] * nsg_415[k];

        t_583[k] = f_13 * nsf0_277[k]
                   - f_14 * nsf1_277[k]
                   + f_3 * pc_y[k] * nsg_416[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, msg_314, nsf0_278, nsf0_279, \
                         nsf1_278, nsf1_279, nsg_417, nsg_418, \
                         nsg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * nsf0_278[k]
                   - f_7 * nsf1_278[k]
                   + f_3 * pc_y[k] * nsg_417[k];

        t_585[k] = f_4 * nsf0_279[k]
                   - f_5 * nsf1_279[k]
                   + f_3 * pc_y[k] * nsg_418[k];

        t_586[k] = f_3 * pc_y[k] * nsg_419[k];

        t_587[k] = f_17 * msg_314[k]
                   + f_1 * nsf0_279[k]
                   - f_2 * nsf1_279[k]
                   + f_3 * pc_z[k] * nsg_419[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, msg_315, msg_420, \
                         msg_423, nsf0_280, nsf0_283, nsf1_280, nsf1_283, nsg_420, \
                         nsg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_11 * msg_420[k]
                   + f_1 * nsf0_280[k]
                   - f_2 * nsf1_280[k]
                   + f_3 * pc_x[k] * nsg_420[k];

        t_589[k] = f_16 * msg_315[k]
                   + f_3 * pc_y[k] * nsg_420[k];

        t_590[k] = f_3 * pc_z[k] * nsg_420[k];

        t_591[k] = f_11 * msg_423[k]
                   + f_6 * nsf0_283[k]
                   - f_7 * nsf1_283[k]
                   + f_3 * pc_x[k] * nsg_423[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, msg_426, nsf0_280, nsf0_286, \
                         nsf1_280, nsf1_286, nsg_421, nsg_422, nsg_423, \
                         nsg_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * nsg_421[k];

        t_593[k] = f_4 * nsf0_280[k]
                   - f_5 * nsf1_280[k]
                   + f_3 * pc_z[k] * nsg_422[k];

        t_594[k] = f_11 * msg_426[k]
                   + f_4 * nsf0_286[k]
                   - f_5 * nsf1_286[k]
                   + f_3 * pc_x[k] * nsg_426[k];

        t_595[k] = f_3 * pc_z[k] * nsg_423[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, msg_320, msg_430, \
                         nsf0_282, nsf1_282, nsg_425, nsg_426, \
                         nsg_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_16 * msg_320[k]
                   + f_3 * pc_y[k] * nsg_425[k];

        t_597[k] = f_6 * nsf0_282[k]
                   - f_7 * nsf1_282[k]
                   + f_3 * pc_z[k] * nsg_425[k];

        t_598[k] = f_11 * msg_430[k]
                   + f_3 * pc_x[k] * nsg_430[k];

        t_599[k] = f_3 * pc_z[k] * nsg_426[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, msg_325, msg_432, msg_433, \
                         msg_434, nsf0_286, nsf1_286, nsg_430, nsg_432, nsg_433, \
                         nsg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_11 * msg_432[k]
                   + f_3 * pc_x[k] * nsg_432[k];

        t_601[k] = f_11 * msg_433[k]
                   + f_3 * pc_x[k] * nsg_433[k];

        t_602[k] = f_11 * msg_434[k]
                   + f_3 * pc_x[k] * nsg_434[k];

        t_603[k] = f_16 * msg_325[k]
                   + f_1 * nsf0_286[k]
                   - f_2 * nsf1_286[k]
                   + f_3 * pc_y[k] * nsg_430[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_441 = buffer.data(msh0 + 441);
    const auto *msh0_444 = buffer.data(msh0 + 444);
    const auto *msh0_447 = buffer.data(msh0 + 447);
    const auto *msh0_456 = buffer.data(msh0 + 456);

    const auto *msg_315 = buffer.data(msg + 315);
    const auto *msg_318 = buffer.data(msg + 318);
    const auto *msg_325 = buffer.data(msg + 325);
    const auto *msg_329 = buffer.data(msg + 329);
    const auto *msg_330 = buffer.data(msg + 330);
    const auto *msg_332 = buffer.data(msg + 332);
    const auto *msg_333 = buffer.data(msg + 333);
    const auto *msg_335 = buffer.data(msg + 335);
    const auto *msg_340 = buffer.data(msg + 340);
    const auto *msg_342 = buffer.data(msg + 342);
    const auto *msg_343 = buffer.data(msg + 343);
    const auto *msg_344 = buffer.data(msg + 344);
    const auto *msg_345 = buffer.data(msg + 345);
    const auto *msg_347 = buffer.data(msg + 347);
    const auto *msg_348 = buffer.data(msg + 348);
    const auto *msg_350 = buffer.data(msg + 350);
    const auto *msg_355 = buffer.data(msg + 355);
    const auto *msg_357 = buffer.data(msg + 357);
    const auto *msg_358 = buffer.data(msg + 358);
    const auto *msg_359 = buffer.data(msg + 359);
    const auto *msg_360 = buffer.data(msg + 360);
    const auto *msg_362 = buffer.data(msg + 362);
    const auto *msg_363 = buffer.data(msg + 363);
    const auto *msg_365 = buffer.data(msg + 365);
    const auto *msg_370 = buffer.data(msg + 370);
    const auto *msg_372 = buffer.data(msg + 372);
    const auto *msg_373 = buffer.data(msg + 373);
    const auto *msg_374 = buffer.data(msg + 374);
    const auto *msg_375 = buffer.data(msg + 375);
    const auto *msg_377 = buffer.data(msg + 377);
    const auto *msg_378 = buffer.data(msg + 378);
    const auto *msg_380 = buffer.data(msg + 380);
    const auto *msg_385 = buffer.data(msg + 385);
    const auto *msg_387 = buffer.data(msg + 387);
    const auto *msg_388 = buffer.data(msg + 388);
    const auto *msg_389 = buffer.data(msg + 389);
    const auto *msg_390 = buffer.data(msg + 390);
    const auto *msg_392 = buffer.data(msg + 392);
    const auto *msg_395 = buffer.data(msg + 395);
    const auto *msg_400 = buffer.data(msg + 400);
    const auto *msg_402 = buffer.data(msg + 402);
    const auto *msg_403 = buffer.data(msg + 403);
    const auto *msg_440 = buffer.data(msg + 440);
    const auto *msg_444 = buffer.data(msg + 444);
    const auto *msg_445 = buffer.data(msg + 445);
    const auto *msg_446 = buffer.data(msg + 446);
    const auto *msg_447 = buffer.data(msg + 447);
    const auto *msg_448 = buffer.data(msg + 448);
    const auto *msg_449 = buffer.data(msg + 449);
    const auto *msg_450 = buffer.data(msg + 450);
    const auto *msg_453 = buffer.data(msg + 453);
    const auto *msg_455 = buffer.data(msg + 455);
    const auto *msg_456 = buffer.data(msg + 456);
    const auto *msg_459 = buffer.data(msg + 459);
    const auto *msg_460 = buffer.data(msg + 460);
    const auto *msg_461 = buffer.data(msg + 461);
    const auto *msg_462 = buffer.data(msg + 462);
    const auto *msg_463 = buffer.data(msg + 463);
    const auto *msg_464 = buffer.data(msg + 464);
    const auto *msg_465 = buffer.data(msg + 465);
    const auto *msg_468 = buffer.data(msg + 468);
    const auto *msg_470 = buffer.data(msg + 470);
    const auto *msg_471 = buffer.data(msg + 471);
    const auto *msg_474 = buffer.data(msg + 474);
    const auto *msg_475 = buffer.data(msg + 475);
    const auto *msg_476 = buffer.data(msg + 476);
    const auto *msg_477 = buffer.data(msg + 477);
    const auto *msg_478 = buffer.data(msg + 478);
    const auto *msg_479 = buffer.data(msg + 479);
    const auto *msg_480 = buffer.data(msg + 480);
    const auto *msg_483 = buffer.data(msg + 483);
    const auto *msg_485 = buffer.data(msg + 485);
    const auto *msg_486 = buffer.data(msg + 486);
    const auto *msg_489 = buffer.data(msg + 489);
    const auto *msg_490 = buffer.data(msg + 490);
    const auto *msg_491 = buffer.data(msg + 491);
    const auto *msg_492 = buffer.data(msg + 492);
    const auto *msg_493 = buffer.data(msg + 493);
    const auto *msg_494 = buffer.data(msg + 494);
    const auto *msg_495 = buffer.data(msg + 495);
    const auto *msg_498 = buffer.data(msg + 498);
    const auto *msg_500 = buffer.data(msg + 500);
    const auto *msg_501 = buffer.data(msg + 501);
    const auto *msg_504 = buffer.data(msg + 504);
    const auto *msg_505 = buffer.data(msg + 505);
    const auto *msg_506 = buffer.data(msg + 506);
    const auto *msg_507 = buffer.data(msg + 507);
    const auto *msg_508 = buffer.data(msg + 508);
    const auto *msg_509 = buffer.data(msg + 509);

    const auto *msh1_441 = buffer.data(msh1 + 441);
    const auto *msh1_444 = buffer.data(msh1 + 444);
    const auto *msh1_447 = buffer.data(msh1 + 447);
    const auto *msh1_456 = buffer.data(msh1 + 456);

    const auto *nsf0_286 = buffer.data(nsf0 + 286);
    const auto *nsf0_287 = buffer.data(nsf0 + 287);
    const auto *nsf0_289 = buffer.data(nsf0 + 289);
    const auto *nsf0_295 = buffer.data(nsf0 + 295);
    const auto *nsf0_298 = buffer.data(nsf0 + 298);
    const auto *nsf0_299 = buffer.data(nsf0 + 299);
    const auto *nsf0_300 = buffer.data(nsf0 + 300);
    const auto *nsf0_303 = buffer.data(nsf0 + 303);
    const auto *nsf0_305 = buffer.data(nsf0 + 305);
    const auto *nsf0_306 = buffer.data(nsf0 + 306);
    const auto *nsf0_308 = buffer.data(nsf0 + 308);
    const auto *nsf0_309 = buffer.data(nsf0 + 309);
    const auto *nsf0_310 = buffer.data(nsf0 + 310);
    const auto *nsf0_313 = buffer.data(nsf0 + 313);
    const auto *nsf0_315 = buffer.data(nsf0 + 315);
    const auto *nsf0_316 = buffer.data(nsf0 + 316);
    const auto *nsf0_318 = buffer.data(nsf0 + 318);
    const auto *nsf0_319 = buffer.data(nsf0 + 319);
    const auto *nsf0_320 = buffer.data(nsf0 + 320);
    const auto *nsf0_323 = buffer.data(nsf0 + 323);
    const auto *nsf0_325 = buffer.data(nsf0 + 325);
    const auto *nsf0_326 = buffer.data(nsf0 + 326);
    const auto *nsf0_328 = buffer.data(nsf0 + 328);
    const auto *nsf0_329 = buffer.data(nsf0 + 329);
    const auto *nsf0_330 = buffer.data(nsf0 + 330);
    const auto *nsf0_333 = buffer.data(nsf0 + 333);
    const auto *nsf0_335 = buffer.data(nsf0 + 335);
    const auto *nsf0_336 = buffer.data(nsf0 + 336);
    const auto *nsf0_338 = buffer.data(nsf0 + 338);
    const auto *nsf0_339 = buffer.data(nsf0 + 339);

    const auto *nsf1_286 = buffer.data(nsf1 + 286);
    const auto *nsf1_287 = buffer.data(nsf1 + 287);
    const auto *nsf1_289 = buffer.data(nsf1 + 289);
    const auto *nsf1_295 = buffer.data(nsf1 + 295);
    const auto *nsf1_298 = buffer.data(nsf1 + 298);
    const auto *nsf1_299 = buffer.data(nsf1 + 299);
    const auto *nsf1_300 = buffer.data(nsf1 + 300);
    const auto *nsf1_303 = buffer.data(nsf1 + 303);
    const auto *nsf1_305 = buffer.data(nsf1 + 305);
    const auto *nsf1_306 = buffer.data(nsf1 + 306);
    const auto *nsf1_308 = buffer.data(nsf1 + 308);
    const auto *nsf1_309 = buffer.data(nsf1 + 309);
    const auto *nsf1_310 = buffer.data(nsf1 + 310);
    const auto *nsf1_313 = buffer.data(nsf1 + 313);
    const auto *nsf1_315 = buffer.data(nsf1 + 315);
    const auto *nsf1_316 = buffer.data(nsf1 + 316);
    const auto *nsf1_318 = buffer.data(nsf1 + 318);
    const auto *nsf1_319 = buffer.data(nsf1 + 319);
    const auto *nsf1_320 = buffer.data(nsf1 + 320);
    const auto *nsf1_323 = buffer.data(nsf1 + 323);
    const auto *nsf1_325 = buffer.data(nsf1 + 325);
    const auto *nsf1_326 = buffer.data(nsf1 + 326);
    const auto *nsf1_328 = buffer.data(nsf1 + 328);
    const auto *nsf1_329 = buffer.data(nsf1 + 329);
    const auto *nsf1_330 = buffer.data(nsf1 + 330);
    const auto *nsf1_333 = buffer.data(nsf1 + 333);
    const auto *nsf1_335 = buffer.data(nsf1 + 335);
    const auto *nsf1_336 = buffer.data(nsf1 + 336);
    const auto *nsf1_338 = buffer.data(nsf1 + 338);
    const auto *nsf1_339 = buffer.data(nsf1 + 339);

    const auto *nsg_430 = buffer.data(nsg + 430);
    const auto *nsg_431 = buffer.data(nsg + 431);
    const auto *nsg_432 = buffer.data(nsg + 432);
    const auto *nsg_434 = buffer.data(nsg + 434);
    const auto *nsg_435 = buffer.data(nsg + 435);
    const auto *nsg_437 = buffer.data(nsg + 437);
    const auto *nsg_438 = buffer.data(nsg + 438);
    const auto *nsg_440 = buffer.data(nsg + 440);
    const auto *nsg_444 = buffer.data(nsg + 444);
    const auto *nsg_445 = buffer.data(nsg + 445);
    const auto *nsg_446 = buffer.data(nsg + 446);
    const auto *nsg_447 = buffer.data(nsg + 447);
    const auto *nsg_448 = buffer.data(nsg + 448);
    const auto *nsg_449 = buffer.data(nsg + 449);
    const auto *nsg_450 = buffer.data(nsg + 450);
    const auto *nsg_452 = buffer.data(nsg + 452);
    const auto *nsg_453 = buffer.data(nsg + 453);
    const auto *nsg_455 = buffer.data(nsg + 455);
    const auto *nsg_456 = buffer.data(nsg + 456);
    const auto *nsg_459 = buffer.data(nsg + 459);
    const auto *nsg_460 = buffer.data(nsg + 460);
    const auto *nsg_461 = buffer.data(nsg + 461);
    const auto *nsg_462 = buffer.data(nsg + 462);
    const auto *nsg_463 = buffer.data(nsg + 463);
    const auto *nsg_464 = buffer.data(nsg + 464);
    const auto *nsg_465 = buffer.data(nsg + 465);
    const auto *nsg_467 = buffer.data(nsg + 467);
    const auto *nsg_468 = buffer.data(nsg + 468);
    const auto *nsg_470 = buffer.data(nsg + 470);
    const auto *nsg_471 = buffer.data(nsg + 471);
    const auto *nsg_474 = buffer.data(nsg + 474);
    const auto *nsg_475 = buffer.data(nsg + 475);
    const auto *nsg_476 = buffer.data(nsg + 476);
    const auto *nsg_477 = buffer.data(nsg + 477);
    const auto *nsg_478 = buffer.data(nsg + 478);
    const auto *nsg_479 = buffer.data(nsg + 479);
    const auto *nsg_480 = buffer.data(nsg + 480);
    const auto *nsg_482 = buffer.data(nsg + 482);
    const auto *nsg_483 = buffer.data(nsg + 483);
    const auto *nsg_485 = buffer.data(nsg + 485);
    const auto *nsg_486 = buffer.data(nsg + 486);
    const auto *nsg_489 = buffer.data(nsg + 489);
    const auto *nsg_490 = buffer.data(nsg + 490);
    const auto *nsg_491 = buffer.data(nsg + 491);
    const auto *nsg_492 = buffer.data(nsg + 492);
    const auto *nsg_493 = buffer.data(nsg + 493);
    const auto *nsg_494 = buffer.data(nsg + 494);
    const auto *nsg_495 = buffer.data(nsg + 495);
    const auto *nsg_497 = buffer.data(nsg + 497);
    const auto *nsg_498 = buffer.data(nsg + 498);
    const auto *nsg_500 = buffer.data(nsg + 500);
    const auto *nsg_501 = buffer.data(nsg + 501);
    const auto *nsg_504 = buffer.data(nsg + 504);
    const auto *nsg_505 = buffer.data(nsg + 505);
    const auto *nsg_506 = buffer.data(nsg + 506);
    const auto *nsg_507 = buffer.data(nsg + 507);
    const auto *nsg_508 = buffer.data(nsg + 508);
    const auto *nsg_509 = buffer.data(nsg + 509);

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pc_y, pc_z, msg_329, nsf0_286, nsf0_287, \
                         nsf1_286, nsf1_287, nsg_430, nsg_431, nsg_432, \
                         nsg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * nsg_430[k];

        t_605[k] = f_4 * nsf0_286[k]
                   - f_5 * nsf1_286[k]
                   + f_3 * pc_z[k] * nsg_431[k];

        t_606[k] = f_6 * nsf0_287[k]
                   - f_7 * nsf1_287[k]
                   + f_3 * pc_z[k] * nsg_432[k];

        t_607[k] = f_16 * msg_329[k]
                   + f_3 * pc_y[k] * nsg_434[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pa_z, pc_y, pc_z, msh0_441, msg_315, \
                         msg_330, msh1_441, nsf0_289, nsf1_289, nsg_434, \
                         nsg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_1 * nsf0_289[k]
                   - f_2 * nsf1_289[k]
                   + f_3 * pc_z[k] * nsg_434[k];

        t_609[k] = pa_z[k] * msh0_441[k]
                   - f_8 * pc_z[k] * msh1_441[k];

        t_610[k] = f_17 * msg_330[k]
                   + f_3 * pc_y[k] * nsg_435[k];

        t_611[k] = f_9 * msg_315[k]
                   + f_3 * pc_z[k] * nsg_435[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pa_z, pc_x, pc_y, pc_z, msh0_444, msg_332, \
                         msg_440, msh1_444, nsf0_295, nsf1_295, nsg_437, \
                         nsg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pa_z[k] * msh0_444[k]
                   - f_8 * pc_z[k] * msh1_444[k];

        t_613[k] = f_17 * msg_332[k]
                   + f_3 * pc_y[k] * nsg_437[k];

        t_614[k] = f_11 * msg_440[k]
                   + f_6 * nsf0_295[k]
                   - f_7 * nsf1_295[k]
                   + f_3 * pc_x[k] * nsg_440[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pa_z, pc_y, pc_z, msh0_447, msg_318, msg_335, \
                         msh1_447, nsg_438, nsg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = pa_z[k] * msh0_447[k]
                   - f_8 * pc_z[k] * msh1_447[k];

        t_616[k] = f_9 * msg_318[k]
                   + f_3 * pc_z[k] * nsg_438[k];

        t_617[k] = f_17 * msg_335[k]
                   + f_3 * pc_y[k] * nsg_440[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pc_x, msg_444, msg_445, msg_446, msg_447, \
                         nsf0_299, nsf1_299, nsg_444, nsg_445, nsg_446, \
                         nsg_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_11 * msg_444[k]
                   + f_4 * nsf0_299[k]
                   - f_5 * nsf1_299[k]
                   + f_3 * pc_x[k] * nsg_444[k];

        t_619[k] = f_11 * msg_445[k]
                   + f_3 * pc_x[k] * nsg_445[k];

        t_620[k] = f_11 * msg_446[k]
                   + f_3 * pc_x[k] * nsg_446[k];

        t_621[k] = f_11 * msg_447[k]
                   + f_3 * pc_x[k] * nsg_447[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pa_z, pc_x, pc_z, msh0_456, msg_325, \
                         msg_448, msg_449, msh1_456, nsg_445, nsg_448, \
                         nsg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_11 * msg_448[k]
                   + f_3 * pc_x[k] * nsg_448[k];

        t_623[k] = f_11 * msg_449[k]
                   + f_3 * pc_x[k] * nsg_449[k];

        t_624[k] = pa_z[k] * msh0_456[k]
                   - f_8 * pc_z[k] * msh1_456[k];

        t_625[k] = f_9 * msg_325[k]
                   + f_3 * pc_z[k] * nsg_445[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, msg_342, msg_343, msg_344, nsf0_298, \
                         nsf0_299, nsf1_298, nsf1_299, nsg_447, nsg_448, \
                         nsg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_17 * msg_342[k]
                   + f_6 * nsf0_298[k]
                   - f_7 * nsf1_298[k]
                   + f_3 * pc_y[k] * nsg_447[k];

        t_627[k] = f_17 * msg_343[k]
                   + f_4 * nsf0_299[k]
                   - f_5 * nsf1_299[k]
                   + f_3 * pc_y[k] * nsg_448[k];

        t_628[k] = f_17 * msg_344[k]
                   + f_3 * pc_y[k] * nsg_449[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, pc_z, msg_329, msg_345, msg_450, \
                         nsf0_299, nsf0_300, nsf1_299, nsf1_300, nsg_449, \
                         nsg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_9 * msg_329[k]
                   + f_1 * nsf0_299[k]
                   - f_2 * nsf1_299[k]
                   + f_3 * pc_z[k] * nsg_449[k];

        t_630[k] = f_11 * msg_450[k]
                   + f_1 * nsf0_300[k]
                   - f_2 * nsf1_300[k]
                   + f_3 * pc_x[k] * nsg_450[k];

        t_631[k] = f_19 * msg_345[k]
                   + f_3 * pc_y[k] * nsg_450[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, pc_x, pc_y, pc_z, msg_330, msg_347, msg_453, \
                         nsf0_303, nsf1_303, nsg_450, nsg_452, \
                         nsg_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_10 * msg_330[k]
                   + f_3 * pc_z[k] * nsg_450[k];

        t_633[k] = f_11 * msg_453[k]
                   + f_6 * nsf0_303[k]
                   - f_7 * nsf1_303[k]
                   + f_3 * pc_x[k] * nsg_453[k];

        t_634[k] = f_19 * msg_347[k]
                   + f_3 * pc_y[k] * nsg_452[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, pc_x, pc_z, msg_333, msg_455, msg_456, nsf0_305, \
                         nsf0_306, nsf1_305, nsf1_306, nsg_453, nsg_455, \
                         nsg_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_11 * msg_455[k]
                   + f_6 * nsf0_305[k]
                   - f_7 * nsf1_305[k]
                   + f_3 * pc_x[k] * nsg_455[k];

        t_636[k] = f_11 * msg_456[k]
                   + f_4 * nsf0_306[k]
                   - f_5 * nsf1_306[k]
                   + f_3 * pc_x[k] * nsg_456[k];

        t_637[k] = f_10 * msg_333[k]
                   + f_3 * pc_z[k] * nsg_453[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pc_x, pc_y, msg_350, msg_459, msg_460, \
                         msg_461, nsf0_309, nsf1_309, nsg_455, nsg_459, nsg_460, \
                         nsg_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_19 * msg_350[k]
                   + f_3 * pc_y[k] * nsg_455[k];

        t_639[k] = f_11 * msg_459[k]
                   + f_4 * nsf0_309[k]
                   - f_5 * nsf1_309[k]
                   + f_3 * pc_x[k] * nsg_459[k];

        t_640[k] = f_11 * msg_460[k]
                   + f_3 * pc_x[k] * nsg_460[k];

        t_641[k] = f_11 * msg_461[k]
                   + f_3 * pc_x[k] * nsg_461[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_y, msg_355, msg_462, msg_463, \
                         msg_464, nsf0_306, nsf1_306, nsg_460, nsg_462, nsg_463, \
                         nsg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_11 * msg_462[k]
                   + f_3 * pc_x[k] * nsg_462[k];

        t_643[k] = f_11 * msg_463[k]
                   + f_3 * pc_x[k] * nsg_463[k];

        t_644[k] = f_11 * msg_464[k]
                   + f_3 * pc_x[k] * nsg_464[k];

        t_645[k] = f_19 * msg_355[k]
                   + f_1 * nsf0_306[k]
                   - f_2 * nsf1_306[k]
                   + f_3 * pc_y[k] * nsg_460[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_y, pc_z, msg_340, msg_357, msg_358, nsf0_308, \
                         nsf0_309, nsf1_308, nsf1_309, nsg_460, nsg_462, \
                         nsg_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_10 * msg_340[k]
                   + f_3 * pc_z[k] * nsg_460[k];

        t_647[k] = f_19 * msg_357[k]
                   + f_6 * nsf0_308[k]
                   - f_7 * nsf1_308[k]
                   + f_3 * pc_y[k] * nsg_462[k];

        t_648[k] = f_19 * msg_358[k]
                   + f_4 * nsf0_309[k]
                   - f_5 * nsf1_309[k]
                   + f_3 * pc_y[k] * nsg_463[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pc_x, pc_y, pc_z, msg_344, msg_359, msg_465, \
                         nsf0_309, nsf0_310, nsf1_309, nsf1_310, nsg_464, \
                         nsg_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_19 * msg_359[k]
                   + f_3 * pc_y[k] * nsg_464[k];

        t_650[k] = f_10 * msg_344[k]
                   + f_1 * nsf0_309[k]
                   - f_2 * nsf1_309[k]
                   + f_3 * pc_z[k] * nsg_464[k];

        t_651[k] = f_11 * msg_465[k]
                   + f_1 * nsf0_310[k]
                   - f_2 * nsf1_310[k]
                   + f_3 * pc_x[k] * nsg_465[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, pc_x, pc_y, pc_z, msg_345, msg_360, \
                         msg_362, msg_468, nsf0_313, nsf1_313, nsg_465, nsg_467, \
                         nsg_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_18 * msg_360[k]
                   + f_3 * pc_y[k] * nsg_465[k];

        t_653[k] = f_11 * msg_345[k]
                   + f_3 * pc_z[k] * nsg_465[k];

        t_654[k] = f_11 * msg_468[k]
                   + f_6 * nsf0_313[k]
                   - f_7 * nsf1_313[k]
                   + f_3 * pc_x[k] * nsg_468[k];

        t_655[k] = f_18 * msg_362[k]
                   + f_3 * pc_y[k] * nsg_467[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_z, msg_348, msg_470, msg_471, nsf0_315, \
                         nsf0_316, nsf1_315, nsf1_316, nsg_468, nsg_470, \
                         nsg_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_11 * msg_470[k]
                   + f_6 * nsf0_315[k]
                   - f_7 * nsf1_315[k]
                   + f_3 * pc_x[k] * nsg_470[k];

        t_657[k] = f_11 * msg_471[k]
                   + f_4 * nsf0_316[k]
                   - f_5 * nsf1_316[k]
                   + f_3 * pc_x[k] * nsg_471[k];

        t_658[k] = f_11 * msg_348[k]
                   + f_3 * pc_z[k] * nsg_468[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, pc_x, pc_y, msg_365, msg_474, msg_475, \
                         msg_476, nsf0_319, nsf1_319, nsg_470, nsg_474, nsg_475, \
                         nsg_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_18 * msg_365[k]
                   + f_3 * pc_y[k] * nsg_470[k];

        t_660[k] = f_11 * msg_474[k]
                   + f_4 * nsf0_319[k]
                   - f_5 * nsf1_319[k]
                   + f_3 * pc_x[k] * nsg_474[k];

        t_661[k] = f_11 * msg_475[k]
                   + f_3 * pc_x[k] * nsg_475[k];

        t_662[k] = f_11 * msg_476[k]
                   + f_3 * pc_x[k] * nsg_476[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, pc_x, pc_y, msg_370, msg_477, msg_478, \
                         msg_479, nsf0_316, nsf1_316, nsg_475, nsg_477, nsg_478, \
                         nsg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_11 * msg_477[k]
                   + f_3 * pc_x[k] * nsg_477[k];

        t_664[k] = f_11 * msg_478[k]
                   + f_3 * pc_x[k] * nsg_478[k];

        t_665[k] = f_11 * msg_479[k]
                   + f_3 * pc_x[k] * nsg_479[k];

        t_666[k] = f_18 * msg_370[k]
                   + f_1 * nsf0_316[k]
                   - f_2 * nsf1_316[k]
                   + f_3 * pc_y[k] * nsg_475[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, pc_z, msg_355, msg_372, msg_373, nsf0_318, \
                         nsf0_319, nsf1_318, nsf1_319, nsg_475, nsg_477, \
                         nsg_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_11 * msg_355[k]
                   + f_3 * pc_z[k] * nsg_475[k];

        t_668[k] = f_18 * msg_372[k]
                   + f_6 * nsf0_318[k]
                   - f_7 * nsf1_318[k]
                   + f_3 * pc_y[k] * nsg_477[k];

        t_669[k] = f_18 * msg_373[k]
                   + f_4 * nsf0_319[k]
                   - f_5 * nsf1_319[k]
                   + f_3 * pc_y[k] * nsg_478[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, msg_359, msg_374, msg_480, \
                         nsf0_319, nsf0_320, nsf1_319, nsf1_320, nsg_479, \
                         nsg_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_18 * msg_374[k]
                   + f_3 * pc_y[k] * nsg_479[k];

        t_671[k] = f_11 * msg_359[k]
                   + f_1 * nsf0_319[k]
                   - f_2 * nsf1_319[k]
                   + f_3 * pc_z[k] * nsg_479[k];

        t_672[k] = f_11 * msg_480[k]
                   + f_1 * nsf0_320[k]
                   - f_2 * nsf1_320[k]
                   + f_3 * pc_x[k] * nsg_480[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, msg_360, msg_375, \
                         msg_377, msg_483, nsf0_323, nsf1_323, nsg_480, nsg_482, \
                         nsg_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_11 * msg_375[k]
                   + f_3 * pc_y[k] * nsg_480[k];

        t_674[k] = f_18 * msg_360[k]
                   + f_3 * pc_z[k] * nsg_480[k];

        t_675[k] = f_11 * msg_483[k]
                   + f_6 * nsf0_323[k]
                   - f_7 * nsf1_323[k]
                   + f_3 * pc_x[k] * nsg_483[k];

        t_676[k] = f_11 * msg_377[k]
                   + f_3 * pc_y[k] * nsg_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, msg_363, msg_485, msg_486, nsf0_325, \
                         nsf0_326, nsf1_325, nsf1_326, nsg_483, nsg_485, \
                         nsg_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_11 * msg_485[k]
                   + f_6 * nsf0_325[k]
                   - f_7 * nsf1_325[k]
                   + f_3 * pc_x[k] * nsg_485[k];

        t_678[k] = f_11 * msg_486[k]
                   + f_4 * nsf0_326[k]
                   - f_5 * nsf1_326[k]
                   + f_3 * pc_x[k] * nsg_486[k];

        t_679[k] = f_18 * msg_363[k]
                   + f_3 * pc_z[k] * nsg_483[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_y, msg_380, msg_489, msg_490, \
                         msg_491, nsf0_329, nsf1_329, nsg_485, nsg_489, nsg_490, \
                         nsg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_11 * msg_380[k]
                   + f_3 * pc_y[k] * nsg_485[k];

        t_681[k] = f_11 * msg_489[k]
                   + f_4 * nsf0_329[k]
                   - f_5 * nsf1_329[k]
                   + f_3 * pc_x[k] * nsg_489[k];

        t_682[k] = f_11 * msg_490[k]
                   + f_3 * pc_x[k] * nsg_490[k];

        t_683[k] = f_11 * msg_491[k]
                   + f_3 * pc_x[k] * nsg_491[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pc_x, pc_y, msg_385, msg_492, msg_493, \
                         msg_494, nsf0_326, nsf1_326, nsg_490, nsg_492, nsg_493, \
                         nsg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_11 * msg_492[k]
                   + f_3 * pc_x[k] * nsg_492[k];

        t_685[k] = f_11 * msg_493[k]
                   + f_3 * pc_x[k] * nsg_493[k];

        t_686[k] = f_11 * msg_494[k]
                   + f_3 * pc_x[k] * nsg_494[k];

        t_687[k] = f_11 * msg_385[k]
                   + f_1 * nsf0_326[k]
                   - f_2 * nsf1_326[k]
                   + f_3 * pc_y[k] * nsg_490[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pc_y, pc_z, msg_370, msg_387, msg_388, nsf0_328, \
                         nsf0_329, nsf1_328, nsf1_329, nsg_490, nsg_492, \
                         nsg_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_18 * msg_370[k]
                   + f_3 * pc_z[k] * nsg_490[k];

        t_689[k] = f_11 * msg_387[k]
                   + f_6 * nsf0_328[k]
                   - f_7 * nsf1_328[k]
                   + f_3 * pc_y[k] * nsg_492[k];

        t_690[k] = f_11 * msg_388[k]
                   + f_4 * nsf0_329[k]
                   - f_5 * nsf1_329[k]
                   + f_3 * pc_y[k] * nsg_493[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, pc_x, pc_y, pc_z, msg_374, msg_389, msg_495, \
                         nsf0_329, nsf0_330, nsf1_329, nsf1_330, nsg_494, \
                         nsg_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_11 * msg_389[k]
                   + f_3 * pc_y[k] * nsg_494[k];

        t_692[k] = f_18 * msg_374[k]
                   + f_1 * nsf0_329[k]
                   - f_2 * nsf1_329[k]
                   + f_3 * pc_z[k] * nsg_494[k];

        t_693[k] = f_11 * msg_495[k]
                   + f_1 * nsf0_330[k]
                   - f_2 * nsf1_330[k]
                   + f_3 * pc_x[k] * nsg_495[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pc_x, pc_y, pc_z, msg_375, msg_390, \
                         msg_392, msg_498, nsf0_333, nsf1_333, nsg_495, nsg_497, \
                         nsg_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_10 * msg_390[k]
                   + f_3 * pc_y[k] * nsg_495[k];

        t_695[k] = f_19 * msg_375[k]
                   + f_3 * pc_z[k] * nsg_495[k];

        t_696[k] = f_11 * msg_498[k]
                   + f_6 * nsf0_333[k]
                   - f_7 * nsf1_333[k]
                   + f_3 * pc_x[k] * nsg_498[k];

        t_697[k] = f_10 * msg_392[k]
                   + f_3 * pc_y[k] * nsg_497[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, pc_z, msg_378, msg_500, msg_501, nsf0_335, \
                         nsf0_336, nsf1_335, nsf1_336, nsg_498, nsg_500, \
                         nsg_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_11 * msg_500[k]
                   + f_6 * nsf0_335[k]
                   - f_7 * nsf1_335[k]
                   + f_3 * pc_x[k] * nsg_500[k];

        t_699[k] = f_11 * msg_501[k]
                   + f_4 * nsf0_336[k]
                   - f_5 * nsf1_336[k]
                   + f_3 * pc_x[k] * nsg_501[k];

        t_700[k] = f_19 * msg_378[k]
                   + f_3 * pc_z[k] * nsg_498[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pc_x, pc_y, msg_395, msg_504, msg_505, \
                         msg_506, nsf0_339, nsf1_339, nsg_500, nsg_504, nsg_505, \
                         nsg_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_10 * msg_395[k]
                   + f_3 * pc_y[k] * nsg_500[k];

        t_702[k] = f_11 * msg_504[k]
                   + f_4 * nsf0_339[k]
                   - f_5 * nsf1_339[k]
                   + f_3 * pc_x[k] * nsg_504[k];

        t_703[k] = f_11 * msg_505[k]
                   + f_3 * pc_x[k] * nsg_505[k];

        t_704[k] = f_11 * msg_506[k]
                   + f_3 * pc_x[k] * nsg_506[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pc_x, pc_y, msg_400, msg_507, msg_508, \
                         msg_509, nsf0_336, nsf1_336, nsg_505, nsg_507, nsg_508, \
                         nsg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_11 * msg_507[k]
                   + f_3 * pc_x[k] * nsg_507[k];

        t_706[k] = f_11 * msg_508[k]
                   + f_3 * pc_x[k] * nsg_508[k];

        t_707[k] = f_11 * msg_509[k]
                   + f_3 * pc_x[k] * nsg_509[k];

        t_708[k] = f_10 * msg_400[k]
                   + f_1 * nsf0_336[k]
                   - f_2 * nsf1_336[k]
                   + f_3 * pc_y[k] * nsg_505[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_y, pc_z, msg_385, msg_402, msg_403, nsf0_338, \
                         nsf0_339, nsf1_338, nsf1_339, nsg_505, nsg_507, \
                         nsg_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_19 * msg_385[k]
                   + f_3 * pc_z[k] * nsg_505[k];

        t_710[k] = f_10 * msg_402[k]
                   + f_6 * nsf0_338[k]
                   - f_7 * nsf1_338[k]
                   + f_3 * pc_y[k] * nsg_507[k];

        t_711[k] = f_10 * msg_403[k]
                   + f_4 * nsf0_339[k]
                   - f_5 * nsf1_339[k]
                   + f_3 * pc_y[k] * nsg_508[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_567 = buffer.data(msh0 + 567);
    const auto *msh0_570 = buffer.data(msh0 + 570);
    const auto *msh0_572 = buffer.data(msh0 + 572);
    const auto *msh0_573 = buffer.data(msh0 + 573);
    const auto *msh0_576 = buffer.data(msh0 + 576);
    const auto *msh0_587 = buffer.data(msh0 + 587);
    const auto *msh0_588 = buffer.data(msh0 + 588);
    const auto *msh0_591 = buffer.data(msh0 + 591);
    const auto *msh0_594 = buffer.data(msh0 + 594);
    const auto *msh0_603 = buffer.data(msh0 + 603);

    const auto *msg_389 = buffer.data(msg + 389);
    const auto *msg_390 = buffer.data(msg + 390);
    const auto *msg_393 = buffer.data(msg + 393);
    const auto *msg_400 = buffer.data(msg + 400);
    const auto *msg_404 = buffer.data(msg + 404);
    const auto *msg_405 = buffer.data(msg + 405);
    const auto *msg_406 = buffer.data(msg + 406);
    const auto *msg_407 = buffer.data(msg + 407);
    const auto *msg_408 = buffer.data(msg + 408);
    const auto *msg_410 = buffer.data(msg + 410);
    const auto *msg_415 = buffer.data(msg + 415);
    const auto *msg_417 = buffer.data(msg + 417);
    const auto *msg_418 = buffer.data(msg + 418);
    const auto *msg_419 = buffer.data(msg + 419);
    const auto *msg_420 = buffer.data(msg + 420);
    const auto *msg_423 = buffer.data(msg + 423);
    const auto *msg_425 = buffer.data(msg + 425);
    const auto *msg_430 = buffer.data(msg + 430);
    const auto *msg_434 = buffer.data(msg + 434);
    const auto *msg_435 = buffer.data(msg + 435);
    const auto *msg_437 = buffer.data(msg + 437);
    const auto *msg_438 = buffer.data(msg + 438);
    const auto *msg_440 = buffer.data(msg + 440);
    const auto *msg_445 = buffer.data(msg + 445);
    const auto *msg_447 = buffer.data(msg + 447);
    const auto *msg_448 = buffer.data(msg + 448);
    const auto *msg_449 = buffer.data(msg + 449);
    const auto *msg_450 = buffer.data(msg + 450);
    const auto *msg_452 = buffer.data(msg + 452);
    const auto *msg_453 = buffer.data(msg + 453);
    const auto *msg_455 = buffer.data(msg + 455);
    const auto *msg_460 = buffer.data(msg + 460);
    const auto *msg_462 = buffer.data(msg + 462);
    const auto *msg_463 = buffer.data(msg + 463);
    const auto *msg_464 = buffer.data(msg + 464);
    const auto *msg_465 = buffer.data(msg + 465);
    const auto *msg_467 = buffer.data(msg + 467);
    const auto *msg_470 = buffer.data(msg + 470);
    const auto *msg_520 = buffer.data(msg + 520);
    const auto *msg_521 = buffer.data(msg + 521);
    const auto *msg_522 = buffer.data(msg + 522);
    const auto *msg_523 = buffer.data(msg + 523);
    const auto *msg_524 = buffer.data(msg + 524);
    const auto *msg_525 = buffer.data(msg + 525);
    const auto *msg_530 = buffer.data(msg + 530);
    const auto *msg_534 = buffer.data(msg + 534);
    const auto *msg_535 = buffer.data(msg + 535);
    const auto *msg_536 = buffer.data(msg + 536);
    const auto *msg_537 = buffer.data(msg + 537);
    const auto *msg_539 = buffer.data(msg + 539);
    const auto *msg_540 = buffer.data(msg + 540);
    const auto *msg_543 = buffer.data(msg + 543);
    const auto *msg_546 = buffer.data(msg + 546);
    const auto *msg_550 = buffer.data(msg + 550);
    const auto *msg_552 = buffer.data(msg + 552);
    const auto *msg_553 = buffer.data(msg + 553);
    const auto *msg_554 = buffer.data(msg + 554);
    const auto *msg_560 = buffer.data(msg + 560);
    const auto *msg_564 = buffer.data(msg + 564);
    const auto *msg_565 = buffer.data(msg + 565);
    const auto *msg_566 = buffer.data(msg + 566);
    const auto *msg_567 = buffer.data(msg + 567);
    const auto *msg_568 = buffer.data(msg + 568);
    const auto *msg_569 = buffer.data(msg + 569);
    const auto *msg_570 = buffer.data(msg + 570);
    const auto *msg_573 = buffer.data(msg + 573);
    const auto *msg_575 = buffer.data(msg + 575);
    const auto *msg_576 = buffer.data(msg + 576);
    const auto *msg_579 = buffer.data(msg + 579);
    const auto *msg_580 = buffer.data(msg + 580);
    const auto *msg_581 = buffer.data(msg + 581);
    const auto *msg_582 = buffer.data(msg + 582);
    const auto *msg_583 = buffer.data(msg + 583);
    const auto *msg_584 = buffer.data(msg + 584);
    const auto *msg_585 = buffer.data(msg + 585);
    const auto *msg_588 = buffer.data(msg + 588);
    const auto *msg_590 = buffer.data(msg + 590);
    const auto *msg_591 = buffer.data(msg + 591);
    const auto *msg_594 = buffer.data(msg + 594);
    const auto *msg_595 = buffer.data(msg + 595);
    const auto *msg_596 = buffer.data(msg + 596);

    const auto *msh1_567 = buffer.data(msh1 + 567);
    const auto *msh1_570 = buffer.data(msh1 + 570);
    const auto *msh1_572 = buffer.data(msh1 + 572);
    const auto *msh1_573 = buffer.data(msh1 + 573);
    const auto *msh1_576 = buffer.data(msh1 + 576);
    const auto *msh1_587 = buffer.data(msh1 + 587);
    const auto *msh1_588 = buffer.data(msh1 + 588);
    const auto *msh1_591 = buffer.data(msh1 + 591);
    const auto *msh1_594 = buffer.data(msh1 + 594);
    const auto *msh1_603 = buffer.data(msh1 + 603);

    const auto *nsf0_339 = buffer.data(nsf0 + 339);
    const auto *nsf0_346 = buffer.data(nsf0 + 346);
    const auto *nsf0_348 = buffer.data(nsf0 + 348);
    const auto *nsf0_349 = buffer.data(nsf0 + 349);
    const auto *nsf0_350 = buffer.data(nsf0 + 350);
    const auto *nsf0_351 = buffer.data(nsf0 + 351);
    const auto *nsf0_352 = buffer.data(nsf0 + 352);
    const auto *nsf0_355 = buffer.data(nsf0 + 355);
    const auto *nsf0_356 = buffer.data(nsf0 + 356);
    const auto *nsf0_357 = buffer.data(nsf0 + 357);
    const auto *nsf0_358 = buffer.data(nsf0 + 358);
    const auto *nsf0_359 = buffer.data(nsf0 + 359);
    const auto *nsf0_360 = buffer.data(nsf0 + 360);
    const auto *nsf0_362 = buffer.data(nsf0 + 362);
    const auto *nsf0_363 = buffer.data(nsf0 + 363);
    const auto *nsf0_366 = buffer.data(nsf0 + 366);
    const auto *nsf0_367 = buffer.data(nsf0 + 367);
    const auto *nsf0_369 = buffer.data(nsf0 + 369);
    const auto *nsf0_375 = buffer.data(nsf0 + 375);
    const auto *nsf0_378 = buffer.data(nsf0 + 378);
    const auto *nsf0_379 = buffer.data(nsf0 + 379);
    const auto *nsf0_380 = buffer.data(nsf0 + 380);
    const auto *nsf0_383 = buffer.data(nsf0 + 383);
    const auto *nsf0_385 = buffer.data(nsf0 + 385);
    const auto *nsf0_386 = buffer.data(nsf0 + 386);
    const auto *nsf0_388 = buffer.data(nsf0 + 388);
    const auto *nsf0_389 = buffer.data(nsf0 + 389);
    const auto *nsf0_390 = buffer.data(nsf0 + 390);
    const auto *nsf0_393 = buffer.data(nsf0 + 393);
    const auto *nsf0_395 = buffer.data(nsf0 + 395);
    const auto *nsf0_396 = buffer.data(nsf0 + 396);
    const auto *nsf0_399 = buffer.data(nsf0 + 399);

    const auto *nsf1_339 = buffer.data(nsf1 + 339);
    const auto *nsf1_346 = buffer.data(nsf1 + 346);
    const auto *nsf1_348 = buffer.data(nsf1 + 348);
    const auto *nsf1_349 = buffer.data(nsf1 + 349);
    const auto *nsf1_350 = buffer.data(nsf1 + 350);
    const auto *nsf1_351 = buffer.data(nsf1 + 351);
    const auto *nsf1_352 = buffer.data(nsf1 + 352);
    const auto *nsf1_355 = buffer.data(nsf1 + 355);
    const auto *nsf1_356 = buffer.data(nsf1 + 356);
    const auto *nsf1_357 = buffer.data(nsf1 + 357);
    const auto *nsf1_358 = buffer.data(nsf1 + 358);
    const auto *nsf1_359 = buffer.data(nsf1 + 359);
    const auto *nsf1_360 = buffer.data(nsf1 + 360);
    const auto *nsf1_362 = buffer.data(nsf1 + 362);
    const auto *nsf1_363 = buffer.data(nsf1 + 363);
    const auto *nsf1_366 = buffer.data(nsf1 + 366);
    const auto *nsf1_367 = buffer.data(nsf1 + 367);
    const auto *nsf1_369 = buffer.data(nsf1 + 369);
    const auto *nsf1_375 = buffer.data(nsf1 + 375);
    const auto *nsf1_378 = buffer.data(nsf1 + 378);
    const auto *nsf1_379 = buffer.data(nsf1 + 379);
    const auto *nsf1_380 = buffer.data(nsf1 + 380);
    const auto *nsf1_383 = buffer.data(nsf1 + 383);
    const auto *nsf1_385 = buffer.data(nsf1 + 385);
    const auto *nsf1_386 = buffer.data(nsf1 + 386);
    const auto *nsf1_388 = buffer.data(nsf1 + 388);
    const auto *nsf1_389 = buffer.data(nsf1 + 389);
    const auto *nsf1_390 = buffer.data(nsf1 + 390);
    const auto *nsf1_393 = buffer.data(nsf1 + 393);
    const auto *nsf1_395 = buffer.data(nsf1 + 395);
    const auto *nsf1_396 = buffer.data(nsf1 + 396);
    const auto *nsf1_399 = buffer.data(nsf1 + 399);

    const auto *nsg_509 = buffer.data(nsg + 509);
    const auto *nsg_510 = buffer.data(nsg + 510);
    const auto *nsg_512 = buffer.data(nsg + 512);
    const auto *nsg_513 = buffer.data(nsg + 513);
    const auto *nsg_515 = buffer.data(nsg + 515);
    const auto *nsg_520 = buffer.data(nsg + 520);
    const auto *nsg_521 = buffer.data(nsg + 521);
    const auto *nsg_522 = buffer.data(nsg + 522);
    const auto *nsg_523 = buffer.data(nsg + 523);
    const auto *nsg_524 = buffer.data(nsg + 524);
    const auto *nsg_525 = buffer.data(nsg + 525);
    const auto *nsg_526 = buffer.data(nsg + 526);
    const auto *nsg_527 = buffer.data(nsg + 527);
    const auto *nsg_528 = buffer.data(nsg + 528);
    const auto *nsg_529 = buffer.data(nsg + 529);
    const auto *nsg_530 = buffer.data(nsg + 530);
    const auto *nsg_534 = buffer.data(nsg + 534);
    const auto *nsg_535 = buffer.data(nsg + 535);
    const auto *nsg_536 = buffer.data(nsg + 536);
    const auto *nsg_537 = buffer.data(nsg + 537);
    const auto *nsg_538 = buffer.data(nsg + 538);
    const auto *nsg_539 = buffer.data(nsg + 539);
    const auto *nsg_540 = buffer.data(nsg + 540);
    const auto *nsg_541 = buffer.data(nsg + 541);
    const auto *nsg_542 = buffer.data(nsg + 542);
    const auto *nsg_543 = buffer.data(nsg + 543);
    const auto *nsg_545 = buffer.data(nsg + 545);
    const auto *nsg_546 = buffer.data(nsg + 546);
    const auto *nsg_550 = buffer.data(nsg + 550);
    const auto *nsg_551 = buffer.data(nsg + 551);
    const auto *nsg_552 = buffer.data(nsg + 552);
    const auto *nsg_553 = buffer.data(nsg + 553);
    const auto *nsg_554 = buffer.data(nsg + 554);
    const auto *nsg_555 = buffer.data(nsg + 555);
    const auto *nsg_557 = buffer.data(nsg + 557);
    const auto *nsg_558 = buffer.data(nsg + 558);
    const auto *nsg_560 = buffer.data(nsg + 560);
    const auto *nsg_564 = buffer.data(nsg + 564);
    const auto *nsg_565 = buffer.data(nsg + 565);
    const auto *nsg_566 = buffer.data(nsg + 566);
    const auto *nsg_567 = buffer.data(nsg + 567);
    const auto *nsg_568 = buffer.data(nsg + 568);
    const auto *nsg_569 = buffer.data(nsg + 569);
    const auto *nsg_570 = buffer.data(nsg + 570);
    const auto *nsg_572 = buffer.data(nsg + 572);
    const auto *nsg_573 = buffer.data(nsg + 573);
    const auto *nsg_575 = buffer.data(nsg + 575);
    const auto *nsg_576 = buffer.data(nsg + 576);
    const auto *nsg_579 = buffer.data(nsg + 579);
    const auto *nsg_580 = buffer.data(nsg + 580);
    const auto *nsg_581 = buffer.data(nsg + 581);
    const auto *nsg_582 = buffer.data(nsg + 582);
    const auto *nsg_583 = buffer.data(nsg + 583);
    const auto *nsg_584 = buffer.data(nsg + 584);
    const auto *nsg_585 = buffer.data(nsg + 585);
    const auto *nsg_587 = buffer.data(nsg + 587);
    const auto *nsg_588 = buffer.data(nsg + 588);
    const auto *nsg_590 = buffer.data(nsg + 590);
    const auto *nsg_591 = buffer.data(nsg + 591);
    const auto *nsg_594 = buffer.data(nsg + 594);
    const auto *nsg_595 = buffer.data(nsg + 595);
    const auto *nsg_596 = buffer.data(nsg + 596);

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pa_y, pc_y, pc_z, msh0_567, msg_389, \
                         msg_404, msg_405, msh1_567, nsf0_339, nsf1_339, nsg_509, \
                         nsg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_10 * msg_404[k]
                   + f_3 * pc_y[k] * nsg_509[k];

        t_713[k] = f_19 * msg_389[k]
                   + f_1 * nsf0_339[k]
                   - f_2 * nsf1_339[k]
                   + f_3 * pc_z[k] * nsg_509[k];

        t_714[k] = pa_y[k] * msh0_567[k]
                   - f_8 * pc_y[k] * msh1_567[k];

        t_715[k] = f_9 * msg_405[k]
                   + f_3 * pc_y[k] * nsg_510[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pa_y, pc_y, pc_z, msh0_570, msh0_572, \
                         msg_390, msg_406, msg_407, msh1_570, msh1_572, nsg_510, \
                         nsg_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_17 * msg_390[k]
                   + f_3 * pc_z[k] * nsg_510[k];

        t_717[k] = pa_y[k] * msh0_570[k]
                   + f_10 * msg_406[k]
                   - f_8 * pc_y[k] * msh1_570[k];

        t_718[k] = f_9 * msg_407[k]
                   + f_3 * pc_y[k] * nsg_512[k];

        t_719[k] = pa_y[k] * msh0_572[k]
                   - f_8 * pc_y[k] * msh1_572[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_y, pc_y, pc_z, msh0_573, msh0_576, \
                         msg_393, msg_408, msg_410, msh1_573, msh1_576, nsg_513, \
                         nsg_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_y[k] * msh0_573[k]
                   + f_11 * msg_408[k]
                   - f_8 * pc_y[k] * msh1_573[k];

        t_721[k] = f_17 * msg_393[k]
                   + f_3 * pc_z[k] * nsg_513[k];

        t_722[k] = f_9 * msg_410[k]
                   + f_3 * pc_y[k] * nsg_515[k];

        t_723[k] = pa_y[k] * msh0_576[k]
                   - f_8 * pc_y[k] * msh1_576[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, msg_520, msg_521, msg_522, \
                         msg_523, msg_524, nsg_520, nsg_521, nsg_522, nsg_523, \
                         nsg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_11 * msg_520[k]
                   + f_3 * pc_x[k] * nsg_520[k];

        t_725[k] = f_11 * msg_521[k]
                   + f_3 * pc_x[k] * nsg_521[k];

        t_726[k] = f_11 * msg_522[k]
                   + f_3 * pc_x[k] * nsg_522[k];

        t_727[k] = f_11 * msg_523[k]
                   + f_3 * pc_x[k] * nsg_523[k];

        t_728[k] = f_11 * msg_524[k]
                   + f_3 * pc_x[k] * nsg_524[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_y, pc_z, msg_400, msg_415, msg_417, nsf0_346, \
                         nsf0_348, nsf1_346, nsf1_348, nsg_520, \
                         nsg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_9 * msg_415[k]
                   + f_1 * nsf0_346[k]
                   - f_2 * nsf1_346[k]
                   + f_3 * pc_y[k] * nsg_520[k];

        t_730[k] = f_17 * msg_400[k]
                   + f_3 * pc_z[k] * nsg_520[k];

        t_731[k] = f_9 * msg_417[k]
                   + f_6 * nsf0_348[k]
                   - f_7 * nsf1_348[k]
                   + f_3 * pc_y[k] * nsg_522[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pa_y, pc_y, msh0_587, msg_418, msg_419, \
                         msh1_587, nsf0_349, nsf1_349, nsg_523, \
                         nsg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_9 * msg_418[k]
                   + f_4 * nsf0_349[k]
                   - f_5 * nsf1_349[k]
                   + f_3 * pc_y[k] * nsg_523[k];

        t_733[k] = f_9 * msg_419[k]
                   + f_3 * pc_y[k] * nsg_524[k];

        t_734[k] = pa_y[k] * msh0_587[k]
                   - f_8 * pc_y[k] * msh1_587[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, pc_x, pc_y, pc_z, msg_405, \
                         msg_525, nsf0_350, nsf1_350, nsg_525, nsg_526, \
                         nsg_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_11 * msg_525[k]
                   + f_1 * nsf0_350[k]
                   - f_2 * nsf1_350[k]
                   + f_3 * pc_x[k] * nsg_525[k];

        t_736[k] = f_3 * pc_y[k] * nsg_525[k];

        t_737[k] = f_16 * msg_405[k]
                   + f_3 * pc_z[k] * nsg_525[k];

        t_738[k] = f_4 * nsf0_350[k]
                   - f_5 * nsf1_350[k]
                   + f_3 * pc_y[k] * nsg_526[k];

        t_739[k] = f_3 * pc_y[k] * nsg_527[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, pc_y, msg_530, nsf0_351, nsf0_352, \
                         nsf0_355, nsf1_351, nsf1_352, nsf1_355, nsg_528, nsg_529, \
                         nsg_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_11 * msg_530[k]
                   + f_6 * nsf0_355[k]
                   - f_7 * nsf1_355[k]
                   + f_3 * pc_x[k] * nsg_530[k];

        t_741[k] = f_6 * nsf0_351[k]
                   - f_7 * nsf1_351[k]
                   + f_3 * pc_y[k] * nsg_528[k];

        t_742[k] = f_4 * nsf0_352[k]
                   - f_5 * nsf1_352[k]
                   + f_3 * pc_y[k] * nsg_529[k];

        t_743[k] = f_3 * pc_y[k] * nsg_530[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, pc_x, msg_534, msg_535, msg_536, msg_537, \
                         nsf0_359, nsf1_359, nsg_534, nsg_535, nsg_536, \
                         nsg_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_11 * msg_534[k]
                   + f_4 * nsf0_359[k]
                   - f_5 * nsf1_359[k]
                   + f_3 * pc_x[k] * nsg_534[k];

        t_745[k] = f_11 * msg_535[k]
                   + f_3 * pc_x[k] * nsg_535[k];

        t_746[k] = f_11 * msg_536[k]
                   + f_3 * pc_x[k] * nsg_536[k];

        t_747[k] = f_11 * msg_537[k]
                   + f_3 * pc_x[k] * nsg_537[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pc_x, pc_y, msg_539, nsf0_356, nsf0_357, \
                         nsf1_356, nsf1_357, nsg_534, nsg_535, nsg_536, \
                         nsg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_3 * pc_y[k] * nsg_534[k];

        t_749[k] = f_11 * msg_539[k]
                   + f_3 * pc_x[k] * nsg_539[k];

        t_750[k] = f_1 * nsf0_356[k]
                   - f_2 * nsf1_356[k]
                   + f_3 * pc_y[k] * nsg_535[k];

        t_751[k] = f_13 * nsf0_357[k]
                   - f_14 * nsf1_357[k]
                   + f_3 * pc_y[k] * nsg_536[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pc_y, pc_z, msg_419, nsf0_358, nsf0_359, \
                         nsf1_358, nsf1_359, nsg_537, nsg_538, \
                         nsg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_6 * nsf0_358[k]
                   - f_7 * nsf1_358[k]
                   + f_3 * pc_y[k] * nsg_537[k];

        t_753[k] = f_4 * nsf0_359[k]
                   - f_5 * nsf1_359[k]
                   + f_3 * pc_y[k] * nsg_538[k];

        t_754[k] = f_3 * pc_y[k] * nsg_539[k];

        t_755[k] = f_16 * msg_419[k]
                   + f_1 * nsf0_359[k]
                   - f_2 * nsf1_359[k]
                   + f_3 * pc_z[k] * nsg_539[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, pc_y, pc_z, msg_420, msg_540, \
                         msg_543, nsf0_360, nsf0_363, nsf1_360, nsf1_363, nsg_540, \
                         nsg_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_10 * msg_540[k]
                   + f_1 * nsf0_360[k]
                   - f_2 * nsf1_360[k]
                   + f_3 * pc_x[k] * nsg_540[k];

        t_757[k] = f_15 * msg_420[k]
                   + f_3 * pc_y[k] * nsg_540[k];

        t_758[k] = f_3 * pc_z[k] * nsg_540[k];

        t_759[k] = f_10 * msg_543[k]
                   + f_6 * nsf0_363[k]
                   - f_7 * nsf1_363[k]
                   + f_3 * pc_x[k] * nsg_543[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, pc_x, pc_z, msg_546, nsf0_360, nsf0_366, \
                         nsf1_360, nsf1_366, nsg_541, nsg_542, nsg_543, \
                         nsg_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_3 * pc_z[k] * nsg_541[k];

        t_761[k] = f_4 * nsf0_360[k]
                   - f_5 * nsf1_360[k]
                   + f_3 * pc_z[k] * nsg_542[k];

        t_762[k] = f_10 * msg_546[k]
                   + f_4 * nsf0_366[k]
                   - f_5 * nsf1_366[k]
                   + f_3 * pc_x[k] * nsg_546[k];

        t_763[k] = f_3 * pc_z[k] * nsg_543[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, t_767, pc_x, pc_y, pc_z, msg_425, msg_550, \
                         nsf0_362, nsf1_362, nsg_545, nsg_546, \
                         nsg_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_15 * msg_425[k]
                   + f_3 * pc_y[k] * nsg_545[k];

        t_765[k] = f_6 * nsf0_362[k]
                   - f_7 * nsf1_362[k]
                   + f_3 * pc_z[k] * nsg_545[k];

        t_766[k] = f_10 * msg_550[k]
                   + f_3 * pc_x[k] * nsg_550[k];

        t_767[k] = f_3 * pc_z[k] * nsg_546[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, msg_430, msg_552, msg_553, \
                         msg_554, nsf0_366, nsf1_366, nsg_550, nsg_552, nsg_553, \
                         nsg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_10 * msg_552[k]
                   + f_3 * pc_x[k] * nsg_552[k];

        t_769[k] = f_10 * msg_553[k]
                   + f_3 * pc_x[k] * nsg_553[k];

        t_770[k] = f_10 * msg_554[k]
                   + f_3 * pc_x[k] * nsg_554[k];

        t_771[k] = f_15 * msg_430[k]
                   + f_1 * nsf0_366[k]
                   - f_2 * nsf1_366[k]
                   + f_3 * pc_y[k] * nsg_550[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, pc_y, pc_z, msg_434, nsf0_366, nsf0_367, \
                         nsf1_366, nsf1_367, nsg_550, nsg_551, nsg_552, \
                         nsg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * nsg_550[k];

        t_773[k] = f_4 * nsf0_366[k]
                   - f_5 * nsf1_366[k]
                   + f_3 * pc_z[k] * nsg_551[k];

        t_774[k] = f_6 * nsf0_367[k]
                   - f_7 * nsf1_367[k]
                   + f_3 * pc_z[k] * nsg_552[k];

        t_775[k] = f_15 * msg_434[k]
                   + f_3 * pc_y[k] * nsg_554[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pa_z, pc_y, pc_z, msh0_588, msg_420, \
                         msg_435, msh1_588, nsf0_369, nsf1_369, nsg_554, \
                         nsg_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_1 * nsf0_369[k]
                   - f_2 * nsf1_369[k]
                   + f_3 * pc_z[k] * nsg_554[k];

        t_777[k] = pa_z[k] * msh0_588[k]
                   - f_8 * pc_z[k] * msh1_588[k];

        t_778[k] = f_16 * msg_435[k]
                   + f_3 * pc_y[k] * nsg_555[k];

        t_779[k] = f_9 * msg_420[k]
                   + f_3 * pc_z[k] * nsg_555[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pa_z, pc_x, pc_y, pc_z, msh0_591, msg_437, \
                         msg_560, msh1_591, nsf0_375, nsf1_375, nsg_557, \
                         nsg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = pa_z[k] * msh0_591[k]
                   - f_8 * pc_z[k] * msh1_591[k];

        t_781[k] = f_16 * msg_437[k]
                   + f_3 * pc_y[k] * nsg_557[k];

        t_782[k] = f_10 * msg_560[k]
                   + f_6 * nsf0_375[k]
                   - f_7 * nsf1_375[k]
                   + f_3 * pc_x[k] * nsg_560[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, pa_z, pc_y, pc_z, msh0_594, msg_423, msg_440, \
                         msh1_594, nsg_558, nsg_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = pa_z[k] * msh0_594[k]
                   - f_8 * pc_z[k] * msh1_594[k];

        t_784[k] = f_9 * msg_423[k]
                   + f_3 * pc_z[k] * nsg_558[k];

        t_785[k] = f_16 * msg_440[k]
                   + f_3 * pc_y[k] * nsg_560[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, pc_x, msg_564, msg_565, msg_566, msg_567, \
                         nsf0_379, nsf1_379, nsg_564, nsg_565, nsg_566, \
                         nsg_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_10 * msg_564[k]
                   + f_4 * nsf0_379[k]
                   - f_5 * nsf1_379[k]
                   + f_3 * pc_x[k] * nsg_564[k];

        t_787[k] = f_10 * msg_565[k]
                   + f_3 * pc_x[k] * nsg_565[k];

        t_788[k] = f_10 * msg_566[k]
                   + f_3 * pc_x[k] * nsg_566[k];

        t_789[k] = f_10 * msg_567[k]
                   + f_3 * pc_x[k] * nsg_567[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_z, pc_x, pc_z, msh0_603, msg_430, \
                         msg_568, msg_569, msh1_603, nsg_565, nsg_568, \
                         nsg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_10 * msg_568[k]
                   + f_3 * pc_x[k] * nsg_568[k];

        t_791[k] = f_10 * msg_569[k]
                   + f_3 * pc_x[k] * nsg_569[k];

        t_792[k] = pa_z[k] * msh0_603[k]
                   - f_8 * pc_z[k] * msh1_603[k];

        t_793[k] = f_9 * msg_430[k]
                   + f_3 * pc_z[k] * nsg_565[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pc_y, msg_447, msg_448, msg_449, nsf0_378, \
                         nsf0_379, nsf1_378, nsf1_379, nsg_567, nsg_568, \
                         nsg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_16 * msg_447[k]
                   + f_6 * nsf0_378[k]
                   - f_7 * nsf1_378[k]
                   + f_3 * pc_y[k] * nsg_567[k];

        t_795[k] = f_16 * msg_448[k]
                   + f_4 * nsf0_379[k]
                   - f_5 * nsf1_379[k]
                   + f_3 * pc_y[k] * nsg_568[k];

        t_796[k] = f_16 * msg_449[k]
                   + f_3 * pc_y[k] * nsg_569[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pc_x, pc_y, pc_z, msg_434, msg_450, msg_570, \
                         nsf0_379, nsf0_380, nsf1_379, nsf1_380, nsg_569, \
                         nsg_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_9 * msg_434[k]
                   + f_1 * nsf0_379[k]
                   - f_2 * nsf1_379[k]
                   + f_3 * pc_z[k] * nsg_569[k];

        t_798[k] = f_10 * msg_570[k]
                   + f_1 * nsf0_380[k]
                   - f_2 * nsf1_380[k]
                   + f_3 * pc_x[k] * nsg_570[k];

        t_799[k] = f_17 * msg_450[k]
                   + f_3 * pc_y[k] * nsg_570[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, msg_435, msg_452, msg_573, \
                         nsf0_383, nsf1_383, nsg_570, nsg_572, \
                         nsg_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_10 * msg_435[k]
                   + f_3 * pc_z[k] * nsg_570[k];

        t_801[k] = f_10 * msg_573[k]
                   + f_6 * nsf0_383[k]
                   - f_7 * nsf1_383[k]
                   + f_3 * pc_x[k] * nsg_573[k];

        t_802[k] = f_17 * msg_452[k]
                   + f_3 * pc_y[k] * nsg_572[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_x, pc_z, msg_438, msg_575, msg_576, nsf0_385, \
                         nsf0_386, nsf1_385, nsf1_386, nsg_573, nsg_575, \
                         nsg_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_10 * msg_575[k]
                   + f_6 * nsf0_385[k]
                   - f_7 * nsf1_385[k]
                   + f_3 * pc_x[k] * nsg_575[k];

        t_804[k] = f_10 * msg_576[k]
                   + f_4 * nsf0_386[k]
                   - f_5 * nsf1_386[k]
                   + f_3 * pc_x[k] * nsg_576[k];

        t_805[k] = f_10 * msg_438[k]
                   + f_3 * pc_z[k] * nsg_573[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, t_809, pc_x, pc_y, msg_455, msg_579, msg_580, \
                         msg_581, nsf0_389, nsf1_389, nsg_575, nsg_579, nsg_580, \
                         nsg_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_17 * msg_455[k]
                   + f_3 * pc_y[k] * nsg_575[k];

        t_807[k] = f_10 * msg_579[k]
                   + f_4 * nsf0_389[k]
                   - f_5 * nsf1_389[k]
                   + f_3 * pc_x[k] * nsg_579[k];

        t_808[k] = f_10 * msg_580[k]
                   + f_3 * pc_x[k] * nsg_580[k];

        t_809[k] = f_10 * msg_581[k]
                   + f_3 * pc_x[k] * nsg_581[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, pc_x, pc_y, msg_460, msg_582, msg_583, \
                         msg_584, nsf0_386, nsf1_386, nsg_580, nsg_582, nsg_583, \
                         nsg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_10 * msg_582[k]
                   + f_3 * pc_x[k] * nsg_582[k];

        t_811[k] = f_10 * msg_583[k]
                   + f_3 * pc_x[k] * nsg_583[k];

        t_812[k] = f_10 * msg_584[k]
                   + f_3 * pc_x[k] * nsg_584[k];

        t_813[k] = f_17 * msg_460[k]
                   + f_1 * nsf0_386[k]
                   - f_2 * nsf1_386[k]
                   + f_3 * pc_y[k] * nsg_580[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, pc_y, pc_z, msg_445, msg_462, msg_463, nsf0_388, \
                         nsf0_389, nsf1_388, nsf1_389, nsg_580, nsg_582, \
                         nsg_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_10 * msg_445[k]
                   + f_3 * pc_z[k] * nsg_580[k];

        t_815[k] = f_17 * msg_462[k]
                   + f_6 * nsf0_388[k]
                   - f_7 * nsf1_388[k]
                   + f_3 * pc_y[k] * nsg_582[k];

        t_816[k] = f_17 * msg_463[k]
                   + f_4 * nsf0_389[k]
                   - f_5 * nsf1_389[k]
                   + f_3 * pc_y[k] * nsg_583[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pc_x, pc_y, pc_z, msg_449, msg_464, msg_585, \
                         nsf0_389, nsf0_390, nsf1_389, nsf1_390, nsg_584, \
                         nsg_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_17 * msg_464[k]
                   + f_3 * pc_y[k] * nsg_584[k];

        t_818[k] = f_10 * msg_449[k]
                   + f_1 * nsf0_389[k]
                   - f_2 * nsf1_389[k]
                   + f_3 * pc_z[k] * nsg_584[k];

        t_819[k] = f_10 * msg_585[k]
                   + f_1 * nsf0_390[k]
                   - f_2 * nsf1_390[k]
                   + f_3 * pc_x[k] * nsg_585[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, t_823, pc_x, pc_y, pc_z, msg_450, msg_465, \
                         msg_467, msg_588, nsf0_393, nsf1_393, nsg_585, nsg_587, \
                         nsg_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_19 * msg_465[k]
                   + f_3 * pc_y[k] * nsg_585[k];

        t_821[k] = f_11 * msg_450[k]
                   + f_3 * pc_z[k] * nsg_585[k];

        t_822[k] = f_10 * msg_588[k]
                   + f_6 * nsf0_393[k]
                   - f_7 * nsf1_393[k]
                   + f_3 * pc_x[k] * nsg_588[k];

        t_823[k] = f_19 * msg_467[k]
                   + f_3 * pc_y[k] * nsg_587[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, msg_453, msg_590, msg_591, nsf0_395, \
                         nsf0_396, nsf1_395, nsf1_396, nsg_588, nsg_590, \
                         nsg_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_10 * msg_590[k]
                   + f_6 * nsf0_395[k]
                   - f_7 * nsf1_395[k]
                   + f_3 * pc_x[k] * nsg_590[k];

        t_825[k] = f_10 * msg_591[k]
                   + f_4 * nsf0_396[k]
                   - f_5 * nsf1_396[k]
                   + f_3 * pc_x[k] * nsg_591[k];

        t_826[k] = f_11 * msg_453[k]
                   + f_3 * pc_z[k] * nsg_588[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, pc_x, pc_y, msg_470, msg_594, msg_595, \
                         msg_596, nsf0_399, nsf1_399, nsg_590, nsg_594, nsg_595, \
                         nsg_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_19 * msg_470[k]
                   + f_3 * pc_y[k] * nsg_590[k];

        t_828[k] = f_10 * msg_594[k]
                   + f_4 * nsf0_399[k]
                   - f_5 * nsf1_399[k]
                   + f_3 * pc_x[k] * nsg_594[k];

        t_829[k] = f_10 * msg_595[k]
                   + f_3 * pc_x[k] * nsg_595[k];

        t_830[k] = f_10 * msg_596[k]
                   + f_3 * pc_x[k] * nsg_596[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);
    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_735 = buffer.data(msh0 + 735);
    const auto *msh0_738 = buffer.data(msh0 + 738);
    const auto *msh0_740 = buffer.data(msh0 + 740);
    const auto *msh0_741 = buffer.data(msh0 + 741);
    const auto *msh0_744 = buffer.data(msh0 + 744);
    const auto *msh0_755 = buffer.data(msh0 + 755);

    const auto *msg_460 = buffer.data(msg + 460);
    const auto *msg_464 = buffer.data(msg + 464);
    const auto *msg_465 = buffer.data(msg + 465);
    const auto *msg_468 = buffer.data(msg + 468);
    const auto *msg_475 = buffer.data(msg + 475);
    const auto *msg_477 = buffer.data(msg + 477);
    const auto *msg_478 = buffer.data(msg + 478);
    const auto *msg_479 = buffer.data(msg + 479);
    const auto *msg_480 = buffer.data(msg + 480);
    const auto *msg_482 = buffer.data(msg + 482);
    const auto *msg_483 = buffer.data(msg + 483);
    const auto *msg_485 = buffer.data(msg + 485);
    const auto *msg_490 = buffer.data(msg + 490);
    const auto *msg_492 = buffer.data(msg + 492);
    const auto *msg_493 = buffer.data(msg + 493);
    const auto *msg_494 = buffer.data(msg + 494);
    const auto *msg_495 = buffer.data(msg + 495);
    const auto *msg_497 = buffer.data(msg + 497);
    const auto *msg_498 = buffer.data(msg + 498);
    const auto *msg_500 = buffer.data(msg + 500);
    const auto *msg_505 = buffer.data(msg + 505);
    const auto *msg_507 = buffer.data(msg + 507);
    const auto *msg_508 = buffer.data(msg + 508);
    const auto *msg_509 = buffer.data(msg + 509);
    const auto *msg_510 = buffer.data(msg + 510);
    const auto *msg_512 = buffer.data(msg + 512);
    const auto *msg_513 = buffer.data(msg + 513);
    const auto *msg_515 = buffer.data(msg + 515);
    const auto *msg_520 = buffer.data(msg + 520);
    const auto *msg_522 = buffer.data(msg + 522);
    const auto *msg_523 = buffer.data(msg + 523);
    const auto *msg_524 = buffer.data(msg + 524);
    const auto *msg_525 = buffer.data(msg + 525);
    const auto *msg_526 = buffer.data(msg + 526);
    const auto *msg_527 = buffer.data(msg + 527);
    const auto *msg_528 = buffer.data(msg + 528);
    const auto *msg_530 = buffer.data(msg + 530);
    const auto *msg_535 = buffer.data(msg + 535);
    const auto *msg_537 = buffer.data(msg + 537);
    const auto *msg_538 = buffer.data(msg + 538);
    const auto *msg_539 = buffer.data(msg + 539);
    const auto *msg_597 = buffer.data(msg + 597);
    const auto *msg_598 = buffer.data(msg + 598);
    const auto *msg_599 = buffer.data(msg + 599);
    const auto *msg_600 = buffer.data(msg + 600);
    const auto *msg_603 = buffer.data(msg + 603);
    const auto *msg_605 = buffer.data(msg + 605);
    const auto *msg_606 = buffer.data(msg + 606);
    const auto *msg_609 = buffer.data(msg + 609);
    const auto *msg_610 = buffer.data(msg + 610);
    const auto *msg_611 = buffer.data(msg + 611);
    const auto *msg_612 = buffer.data(msg + 612);
    const auto *msg_613 = buffer.data(msg + 613);
    const auto *msg_614 = buffer.data(msg + 614);
    const auto *msg_615 = buffer.data(msg + 615);
    const auto *msg_618 = buffer.data(msg + 618);
    const auto *msg_620 = buffer.data(msg + 620);
    const auto *msg_621 = buffer.data(msg + 621);
    const auto *msg_624 = buffer.data(msg + 624);
    const auto *msg_625 = buffer.data(msg + 625);
    const auto *msg_626 = buffer.data(msg + 626);
    const auto *msg_627 = buffer.data(msg + 627);
    const auto *msg_628 = buffer.data(msg + 628);
    const auto *msg_629 = buffer.data(msg + 629);
    const auto *msg_630 = buffer.data(msg + 630);
    const auto *msg_633 = buffer.data(msg + 633);
    const auto *msg_635 = buffer.data(msg + 635);
    const auto *msg_636 = buffer.data(msg + 636);
    const auto *msg_639 = buffer.data(msg + 639);
    const auto *msg_640 = buffer.data(msg + 640);
    const auto *msg_641 = buffer.data(msg + 641);
    const auto *msg_642 = buffer.data(msg + 642);
    const auto *msg_643 = buffer.data(msg + 643);
    const auto *msg_644 = buffer.data(msg + 644);
    const auto *msg_655 = buffer.data(msg + 655);
    const auto *msg_656 = buffer.data(msg + 656);
    const auto *msg_657 = buffer.data(msg + 657);
    const auto *msg_658 = buffer.data(msg + 658);
    const auto *msg_659 = buffer.data(msg + 659);
    const auto *msg_660 = buffer.data(msg + 660);
    const auto *msg_665 = buffer.data(msg + 665);
    const auto *msg_669 = buffer.data(msg + 669);
    const auto *msg_670 = buffer.data(msg + 670);
    const auto *msg_671 = buffer.data(msg + 671);
    const auto *msg_672 = buffer.data(msg + 672);
    const auto *msg_674 = buffer.data(msg + 674);

    const auto *msh1_735 = buffer.data(msh1 + 735);
    const auto *msh1_738 = buffer.data(msh1 + 738);
    const auto *msh1_740 = buffer.data(msh1 + 740);
    const auto *msh1_741 = buffer.data(msh1 + 741);
    const auto *msh1_744 = buffer.data(msh1 + 744);
    const auto *msh1_755 = buffer.data(msh1 + 755);

    const auto *nsf0_396 = buffer.data(nsf0 + 396);
    const auto *nsf0_398 = buffer.data(nsf0 + 398);
    const auto *nsf0_399 = buffer.data(nsf0 + 399);
    const auto *nsf0_400 = buffer.data(nsf0 + 400);
    const auto *nsf0_403 = buffer.data(nsf0 + 403);
    const auto *nsf0_405 = buffer.data(nsf0 + 405);
    const auto *nsf0_406 = buffer.data(nsf0 + 406);
    const auto *nsf0_408 = buffer.data(nsf0 + 408);
    const auto *nsf0_409 = buffer.data(nsf0 + 409);
    const auto *nsf0_410 = buffer.data(nsf0 + 410);
    const auto *nsf0_413 = buffer.data(nsf0 + 413);
    const auto *nsf0_415 = buffer.data(nsf0 + 415);
    const auto *nsf0_416 = buffer.data(nsf0 + 416);
    const auto *nsf0_418 = buffer.data(nsf0 + 418);
    const auto *nsf0_419 = buffer.data(nsf0 + 419);
    const auto *nsf0_420 = buffer.data(nsf0 + 420);
    const auto *nsf0_423 = buffer.data(nsf0 + 423);
    const auto *nsf0_425 = buffer.data(nsf0 + 425);
    const auto *nsf0_426 = buffer.data(nsf0 + 426);
    const auto *nsf0_428 = buffer.data(nsf0 + 428);
    const auto *nsf0_429 = buffer.data(nsf0 + 429);
    const auto *nsf0_436 = buffer.data(nsf0 + 436);
    const auto *nsf0_438 = buffer.data(nsf0 + 438);
    const auto *nsf0_439 = buffer.data(nsf0 + 439);
    const auto *nsf0_440 = buffer.data(nsf0 + 440);
    const auto *nsf0_441 = buffer.data(nsf0 + 441);
    const auto *nsf0_442 = buffer.data(nsf0 + 442);
    const auto *nsf0_445 = buffer.data(nsf0 + 445);
    const auto *nsf0_446 = buffer.data(nsf0 + 446);
    const auto *nsf0_447 = buffer.data(nsf0 + 447);
    const auto *nsf0_448 = buffer.data(nsf0 + 448);
    const auto *nsf0_449 = buffer.data(nsf0 + 449);

    const auto *nsf1_396 = buffer.data(nsf1 + 396);
    const auto *nsf1_398 = buffer.data(nsf1 + 398);
    const auto *nsf1_399 = buffer.data(nsf1 + 399);
    const auto *nsf1_400 = buffer.data(nsf1 + 400);
    const auto *nsf1_403 = buffer.data(nsf1 + 403);
    const auto *nsf1_405 = buffer.data(nsf1 + 405);
    const auto *nsf1_406 = buffer.data(nsf1 + 406);
    const auto *nsf1_408 = buffer.data(nsf1 + 408);
    const auto *nsf1_409 = buffer.data(nsf1 + 409);
    const auto *nsf1_410 = buffer.data(nsf1 + 410);
    const auto *nsf1_413 = buffer.data(nsf1 + 413);
    const auto *nsf1_415 = buffer.data(nsf1 + 415);
    const auto *nsf1_416 = buffer.data(nsf1 + 416);
    const auto *nsf1_418 = buffer.data(nsf1 + 418);
    const auto *nsf1_419 = buffer.data(nsf1 + 419);
    const auto *nsf1_420 = buffer.data(nsf1 + 420);
    const auto *nsf1_423 = buffer.data(nsf1 + 423);
    const auto *nsf1_425 = buffer.data(nsf1 + 425);
    const auto *nsf1_426 = buffer.data(nsf1 + 426);
    const auto *nsf1_428 = buffer.data(nsf1 + 428);
    const auto *nsf1_429 = buffer.data(nsf1 + 429);
    const auto *nsf1_436 = buffer.data(nsf1 + 436);
    const auto *nsf1_438 = buffer.data(nsf1 + 438);
    const auto *nsf1_439 = buffer.data(nsf1 + 439);
    const auto *nsf1_440 = buffer.data(nsf1 + 440);
    const auto *nsf1_441 = buffer.data(nsf1 + 441);
    const auto *nsf1_442 = buffer.data(nsf1 + 442);
    const auto *nsf1_445 = buffer.data(nsf1 + 445);
    const auto *nsf1_446 = buffer.data(nsf1 + 446);
    const auto *nsf1_447 = buffer.data(nsf1 + 447);
    const auto *nsf1_448 = buffer.data(nsf1 + 448);
    const auto *nsf1_449 = buffer.data(nsf1 + 449);

    const auto *nsg_595 = buffer.data(nsg + 595);
    const auto *nsg_597 = buffer.data(nsg + 597);
    const auto *nsg_598 = buffer.data(nsg + 598);
    const auto *nsg_599 = buffer.data(nsg + 599);
    const auto *nsg_600 = buffer.data(nsg + 600);
    const auto *nsg_602 = buffer.data(nsg + 602);
    const auto *nsg_603 = buffer.data(nsg + 603);
    const auto *nsg_605 = buffer.data(nsg + 605);
    const auto *nsg_606 = buffer.data(nsg + 606);
    const auto *nsg_609 = buffer.data(nsg + 609);
    const auto *nsg_610 = buffer.data(nsg + 610);
    const auto *nsg_611 = buffer.data(nsg + 611);
    const auto *nsg_612 = buffer.data(nsg + 612);
    const auto *nsg_613 = buffer.data(nsg + 613);
    const auto *nsg_614 = buffer.data(nsg + 614);
    const auto *nsg_615 = buffer.data(nsg + 615);
    const auto *nsg_617 = buffer.data(nsg + 617);
    const auto *nsg_618 = buffer.data(nsg + 618);
    const auto *nsg_620 = buffer.data(nsg + 620);
    const auto *nsg_621 = buffer.data(nsg + 621);
    const auto *nsg_624 = buffer.data(nsg + 624);
    const auto *nsg_625 = buffer.data(nsg + 625);
    const auto *nsg_626 = buffer.data(nsg + 626);
    const auto *nsg_627 = buffer.data(nsg + 627);
    const auto *nsg_628 = buffer.data(nsg + 628);
    const auto *nsg_629 = buffer.data(nsg + 629);
    const auto *nsg_630 = buffer.data(nsg + 630);
    const auto *nsg_632 = buffer.data(nsg + 632);
    const auto *nsg_633 = buffer.data(nsg + 633);
    const auto *nsg_635 = buffer.data(nsg + 635);
    const auto *nsg_636 = buffer.data(nsg + 636);
    const auto *nsg_639 = buffer.data(nsg + 639);
    const auto *nsg_640 = buffer.data(nsg + 640);
    const auto *nsg_641 = buffer.data(nsg + 641);
    const auto *nsg_642 = buffer.data(nsg + 642);
    const auto *nsg_643 = buffer.data(nsg + 643);
    const auto *nsg_644 = buffer.data(nsg + 644);
    const auto *nsg_645 = buffer.data(nsg + 645);
    const auto *nsg_647 = buffer.data(nsg + 647);
    const auto *nsg_648 = buffer.data(nsg + 648);
    const auto *nsg_650 = buffer.data(nsg + 650);
    const auto *nsg_655 = buffer.data(nsg + 655);
    const auto *nsg_656 = buffer.data(nsg + 656);
    const auto *nsg_657 = buffer.data(nsg + 657);
    const auto *nsg_658 = buffer.data(nsg + 658);
    const auto *nsg_659 = buffer.data(nsg + 659);
    const auto *nsg_660 = buffer.data(nsg + 660);
    const auto *nsg_661 = buffer.data(nsg + 661);
    const auto *nsg_662 = buffer.data(nsg + 662);
    const auto *nsg_663 = buffer.data(nsg + 663);
    const auto *nsg_664 = buffer.data(nsg + 664);
    const auto *nsg_665 = buffer.data(nsg + 665);
    const auto *nsg_669 = buffer.data(nsg + 669);
    const auto *nsg_670 = buffer.data(nsg + 670);
    const auto *nsg_671 = buffer.data(nsg + 671);
    const auto *nsg_672 = buffer.data(nsg + 672);
    const auto *nsg_673 = buffer.data(nsg + 673);
    const auto *nsg_674 = buffer.data(nsg + 674);

#pragma omp simd aligned(t_831, t_832, t_833, t_834, pc_x, pc_y, msg_475, msg_597, msg_598, \
                         msg_599, nsf0_396, nsf1_396, nsg_595, nsg_597, nsg_598, \
                         nsg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_10 * msg_597[k]
                   + f_3 * pc_x[k] * nsg_597[k];

        t_832[k] = f_10 * msg_598[k]
                   + f_3 * pc_x[k] * nsg_598[k];

        t_833[k] = f_10 * msg_599[k]
                   + f_3 * pc_x[k] * nsg_599[k];

        t_834[k] = f_19 * msg_475[k]
                   + f_1 * nsf0_396[k]
                   - f_2 * nsf1_396[k]
                   + f_3 * pc_y[k] * nsg_595[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, pc_y, pc_z, msg_460, msg_477, msg_478, nsf0_398, \
                         nsf0_399, nsf1_398, nsf1_399, nsg_595, nsg_597, \
                         nsg_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_11 * msg_460[k]
                   + f_3 * pc_z[k] * nsg_595[k];

        t_836[k] = f_19 * msg_477[k]
                   + f_6 * nsf0_398[k]
                   - f_7 * nsf1_398[k]
                   + f_3 * pc_y[k] * nsg_597[k];

        t_837[k] = f_19 * msg_478[k]
                   + f_4 * nsf0_399[k]
                   - f_5 * nsf1_399[k]
                   + f_3 * pc_y[k] * nsg_598[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, pc_x, pc_y, pc_z, msg_464, msg_479, msg_600, \
                         nsf0_399, nsf0_400, nsf1_399, nsf1_400, nsg_599, \
                         nsg_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_19 * msg_479[k]
                   + f_3 * pc_y[k] * nsg_599[k];

        t_839[k] = f_11 * msg_464[k]
                   + f_1 * nsf0_399[k]
                   - f_2 * nsf1_399[k]
                   + f_3 * pc_z[k] * nsg_599[k];

        t_840[k] = f_10 * msg_600[k]
                   + f_1 * nsf0_400[k]
                   - f_2 * nsf1_400[k]
                   + f_3 * pc_x[k] * nsg_600[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pc_x, pc_y, pc_z, msg_465, msg_480, \
                         msg_482, msg_603, nsf0_403, nsf1_403, nsg_600, nsg_602, \
                         nsg_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_18 * msg_480[k]
                   + f_3 * pc_y[k] * nsg_600[k];

        t_842[k] = f_18 * msg_465[k]
                   + f_3 * pc_z[k] * nsg_600[k];

        t_843[k] = f_10 * msg_603[k]
                   + f_6 * nsf0_403[k]
                   - f_7 * nsf1_403[k]
                   + f_3 * pc_x[k] * nsg_603[k];

        t_844[k] = f_18 * msg_482[k]
                   + f_3 * pc_y[k] * nsg_602[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_z, msg_468, msg_605, msg_606, nsf0_405, \
                         nsf0_406, nsf1_405, nsf1_406, nsg_603, nsg_605, \
                         nsg_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_10 * msg_605[k]
                   + f_6 * nsf0_405[k]
                   - f_7 * nsf1_405[k]
                   + f_3 * pc_x[k] * nsg_605[k];

        t_846[k] = f_10 * msg_606[k]
                   + f_4 * nsf0_406[k]
                   - f_5 * nsf1_406[k]
                   + f_3 * pc_x[k] * nsg_606[k];

        t_847[k] = f_18 * msg_468[k]
                   + f_3 * pc_z[k] * nsg_603[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, pc_y, msg_485, msg_609, msg_610, \
                         msg_611, nsf0_409, nsf1_409, nsg_605, nsg_609, nsg_610, \
                         nsg_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_18 * msg_485[k]
                   + f_3 * pc_y[k] * nsg_605[k];

        t_849[k] = f_10 * msg_609[k]
                   + f_4 * nsf0_409[k]
                   - f_5 * nsf1_409[k]
                   + f_3 * pc_x[k] * nsg_609[k];

        t_850[k] = f_10 * msg_610[k]
                   + f_3 * pc_x[k] * nsg_610[k];

        t_851[k] = f_10 * msg_611[k]
                   + f_3 * pc_x[k] * nsg_611[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, pc_y, msg_490, msg_612, msg_613, \
                         msg_614, nsf0_406, nsf1_406, nsg_610, nsg_612, nsg_613, \
                         nsg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_10 * msg_612[k]
                   + f_3 * pc_x[k] * nsg_612[k];

        t_853[k] = f_10 * msg_613[k]
                   + f_3 * pc_x[k] * nsg_613[k];

        t_854[k] = f_10 * msg_614[k]
                   + f_3 * pc_x[k] * nsg_614[k];

        t_855[k] = f_18 * msg_490[k]
                   + f_1 * nsf0_406[k]
                   - f_2 * nsf1_406[k]
                   + f_3 * pc_y[k] * nsg_610[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, msg_475, msg_492, msg_493, nsf0_408, \
                         nsf0_409, nsf1_408, nsf1_409, nsg_610, nsg_612, \
                         nsg_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_18 * msg_475[k]
                   + f_3 * pc_z[k] * nsg_610[k];

        t_857[k] = f_18 * msg_492[k]
                   + f_6 * nsf0_408[k]
                   - f_7 * nsf1_408[k]
                   + f_3 * pc_y[k] * nsg_612[k];

        t_858[k] = f_18 * msg_493[k]
                   + f_4 * nsf0_409[k]
                   - f_5 * nsf1_409[k]
                   + f_3 * pc_y[k] * nsg_613[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_x, pc_y, pc_z, msg_479, msg_494, msg_615, \
                         nsf0_409, nsf0_410, nsf1_409, nsf1_410, nsg_614, \
                         nsg_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_18 * msg_494[k]
                   + f_3 * pc_y[k] * nsg_614[k];

        t_860[k] = f_18 * msg_479[k]
                   + f_1 * nsf0_409[k]
                   - f_2 * nsf1_409[k]
                   + f_3 * pc_z[k] * nsg_614[k];

        t_861[k] = f_10 * msg_615[k]
                   + f_1 * nsf0_410[k]
                   - f_2 * nsf1_410[k]
                   + f_3 * pc_x[k] * nsg_615[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pc_x, pc_y, pc_z, msg_480, msg_495, \
                         msg_497, msg_618, nsf0_413, nsf1_413, nsg_615, nsg_617, \
                         nsg_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_11 * msg_495[k]
                   + f_3 * pc_y[k] * nsg_615[k];

        t_863[k] = f_19 * msg_480[k]
                   + f_3 * pc_z[k] * nsg_615[k];

        t_864[k] = f_10 * msg_618[k]
                   + f_6 * nsf0_413[k]
                   - f_7 * nsf1_413[k]
                   + f_3 * pc_x[k] * nsg_618[k];

        t_865[k] = f_11 * msg_497[k]
                   + f_3 * pc_y[k] * nsg_617[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_z, msg_483, msg_620, msg_621, nsf0_415, \
                         nsf0_416, nsf1_415, nsf1_416, nsg_618, nsg_620, \
                         nsg_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_10 * msg_620[k]
                   + f_6 * nsf0_415[k]
                   - f_7 * nsf1_415[k]
                   + f_3 * pc_x[k] * nsg_620[k];

        t_867[k] = f_10 * msg_621[k]
                   + f_4 * nsf0_416[k]
                   - f_5 * nsf1_416[k]
                   + f_3 * pc_x[k] * nsg_621[k];

        t_868[k] = f_19 * msg_483[k]
                   + f_3 * pc_z[k] * nsg_618[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, msg_500, msg_624, msg_625, \
                         msg_626, nsf0_419, nsf1_419, nsg_620, nsg_624, nsg_625, \
                         nsg_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_11 * msg_500[k]
                   + f_3 * pc_y[k] * nsg_620[k];

        t_870[k] = f_10 * msg_624[k]
                   + f_4 * nsf0_419[k]
                   - f_5 * nsf1_419[k]
                   + f_3 * pc_x[k] * nsg_624[k];

        t_871[k] = f_10 * msg_625[k]
                   + f_3 * pc_x[k] * nsg_625[k];

        t_872[k] = f_10 * msg_626[k]
                   + f_3 * pc_x[k] * nsg_626[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pc_x, pc_y, msg_505, msg_627, msg_628, \
                         msg_629, nsf0_416, nsf1_416, nsg_625, nsg_627, nsg_628, \
                         nsg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_10 * msg_627[k]
                   + f_3 * pc_x[k] * nsg_627[k];

        t_874[k] = f_10 * msg_628[k]
                   + f_3 * pc_x[k] * nsg_628[k];

        t_875[k] = f_10 * msg_629[k]
                   + f_3 * pc_x[k] * nsg_629[k];

        t_876[k] = f_11 * msg_505[k]
                   + f_1 * nsf0_416[k]
                   - f_2 * nsf1_416[k]
                   + f_3 * pc_y[k] * nsg_625[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_y, pc_z, msg_490, msg_507, msg_508, nsf0_418, \
                         nsf0_419, nsf1_418, nsf1_419, nsg_625, nsg_627, \
                         nsg_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_19 * msg_490[k]
                   + f_3 * pc_z[k] * nsg_625[k];

        t_878[k] = f_11 * msg_507[k]
                   + f_6 * nsf0_418[k]
                   - f_7 * nsf1_418[k]
                   + f_3 * pc_y[k] * nsg_627[k];

        t_879[k] = f_11 * msg_508[k]
                   + f_4 * nsf0_419[k]
                   - f_5 * nsf1_419[k]
                   + f_3 * pc_y[k] * nsg_628[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, pc_x, pc_y, pc_z, msg_494, msg_509, msg_630, \
                         nsf0_419, nsf0_420, nsf1_419, nsf1_420, nsg_629, \
                         nsg_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_11 * msg_509[k]
                   + f_3 * pc_y[k] * nsg_629[k];

        t_881[k] = f_19 * msg_494[k]
                   + f_1 * nsf0_419[k]
                   - f_2 * nsf1_419[k]
                   + f_3 * pc_z[k] * nsg_629[k];

        t_882[k] = f_10 * msg_630[k]
                   + f_1 * nsf0_420[k]
                   - f_2 * nsf1_420[k]
                   + f_3 * pc_x[k] * nsg_630[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, pc_x, pc_y, pc_z, msg_495, msg_510, \
                         msg_512, msg_633, nsf0_423, nsf1_423, nsg_630, nsg_632, \
                         nsg_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_10 * msg_510[k]
                   + f_3 * pc_y[k] * nsg_630[k];

        t_884[k] = f_17 * msg_495[k]
                   + f_3 * pc_z[k] * nsg_630[k];

        t_885[k] = f_10 * msg_633[k]
                   + f_6 * nsf0_423[k]
                   - f_7 * nsf1_423[k]
                   + f_3 * pc_x[k] * nsg_633[k];

        t_886[k] = f_10 * msg_512[k]
                   + f_3 * pc_y[k] * nsg_632[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, pc_x, pc_z, msg_498, msg_635, msg_636, nsf0_425, \
                         nsf0_426, nsf1_425, nsf1_426, nsg_633, nsg_635, \
                         nsg_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_10 * msg_635[k]
                   + f_6 * nsf0_425[k]
                   - f_7 * nsf1_425[k]
                   + f_3 * pc_x[k] * nsg_635[k];

        t_888[k] = f_10 * msg_636[k]
                   + f_4 * nsf0_426[k]
                   - f_5 * nsf1_426[k]
                   + f_3 * pc_x[k] * nsg_636[k];

        t_889[k] = f_17 * msg_498[k]
                   + f_3 * pc_z[k] * nsg_633[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pc_x, pc_y, msg_515, msg_639, msg_640, \
                         msg_641, nsf0_429, nsf1_429, nsg_635, nsg_639, nsg_640, \
                         nsg_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_10 * msg_515[k]
                   + f_3 * pc_y[k] * nsg_635[k];

        t_891[k] = f_10 * msg_639[k]
                   + f_4 * nsf0_429[k]
                   - f_5 * nsf1_429[k]
                   + f_3 * pc_x[k] * nsg_639[k];

        t_892[k] = f_10 * msg_640[k]
                   + f_3 * pc_x[k] * nsg_640[k];

        t_893[k] = f_10 * msg_641[k]
                   + f_3 * pc_x[k] * nsg_641[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pc_x, pc_y, msg_520, msg_642, msg_643, \
                         msg_644, nsf0_426, nsf1_426, nsg_640, nsg_642, nsg_643, \
                         nsg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_10 * msg_642[k]
                   + f_3 * pc_x[k] * nsg_642[k];

        t_895[k] = f_10 * msg_643[k]
                   + f_3 * pc_x[k] * nsg_643[k];

        t_896[k] = f_10 * msg_644[k]
                   + f_3 * pc_x[k] * nsg_644[k];

        t_897[k] = f_10 * msg_520[k]
                   + f_1 * nsf0_426[k]
                   - f_2 * nsf1_426[k]
                   + f_3 * pc_y[k] * nsg_640[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pc_y, pc_z, msg_505, msg_522, msg_523, nsf0_428, \
                         nsf0_429, nsf1_428, nsf1_429, nsg_640, nsg_642, \
                         nsg_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_17 * msg_505[k]
                   + f_3 * pc_z[k] * nsg_640[k];

        t_899[k] = f_10 * msg_522[k]
                   + f_6 * nsf0_428[k]
                   - f_7 * nsf1_428[k]
                   + f_3 * pc_y[k] * nsg_642[k];

        t_900[k] = f_10 * msg_523[k]
                   + f_4 * nsf0_429[k]
                   - f_5 * nsf1_429[k]
                   + f_3 * pc_y[k] * nsg_643[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pa_y, pc_y, pc_z, msh0_735, msg_509, \
                         msg_524, msg_525, msh1_735, nsf0_429, nsf1_429, nsg_644, \
                         nsg_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_10 * msg_524[k]
                   + f_3 * pc_y[k] * nsg_644[k];

        t_902[k] = f_17 * msg_509[k]
                   + f_1 * nsf0_429[k]
                   - f_2 * nsf1_429[k]
                   + f_3 * pc_z[k] * nsg_644[k];

        t_903[k] = pa_y[k] * msh0_735[k]
                   - f_8 * pc_y[k] * msh1_735[k];

        t_904[k] = f_9 * msg_525[k]
                   + f_3 * pc_y[k] * nsg_645[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, pa_y, pc_y, pc_z, msh0_738, msh0_740, \
                         msg_510, msg_526, msg_527, msh1_738, msh1_740, nsg_645, \
                         nsg_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_16 * msg_510[k]
                   + f_3 * pc_z[k] * nsg_645[k];

        t_906[k] = pa_y[k] * msh0_738[k]
                   + f_10 * msg_526[k]
                   - f_8 * pc_y[k] * msh1_738[k];

        t_907[k] = f_9 * msg_527[k]
                   + f_3 * pc_y[k] * nsg_647[k];

        t_908[k] = pa_y[k] * msh0_740[k]
                   - f_8 * pc_y[k] * msh1_740[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, pa_y, pc_y, pc_z, msh0_741, msh0_744, \
                         msg_513, msg_528, msg_530, msh1_741, msh1_744, nsg_648, \
                         nsg_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = pa_y[k] * msh0_741[k]
                   + f_11 * msg_528[k]
                   - f_8 * pc_y[k] * msh1_741[k];

        t_910[k] = f_16 * msg_513[k]
                   + f_3 * pc_z[k] * nsg_648[k];

        t_911[k] = f_9 * msg_530[k]
                   + f_3 * pc_y[k] * nsg_650[k];

        t_912[k] = pa_y[k] * msh0_744[k]
                   - f_8 * pc_y[k] * msh1_744[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, t_917, pc_x, msg_655, msg_656, msg_657, \
                         msg_658, msg_659, nsg_655, nsg_656, nsg_657, nsg_658, \
                         nsg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = f_10 * msg_655[k]
                   + f_3 * pc_x[k] * nsg_655[k];

        t_914[k] = f_10 * msg_656[k]
                   + f_3 * pc_x[k] * nsg_656[k];

        t_915[k] = f_10 * msg_657[k]
                   + f_3 * pc_x[k] * nsg_657[k];

        t_916[k] = f_10 * msg_658[k]
                   + f_3 * pc_x[k] * nsg_658[k];

        t_917[k] = f_10 * msg_659[k]
                   + f_3 * pc_x[k] * nsg_659[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_y, pc_z, msg_520, msg_535, msg_537, nsf0_436, \
                         nsf0_438, nsf1_436, nsf1_438, nsg_655, \
                         nsg_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_9 * msg_535[k]
                   + f_1 * nsf0_436[k]
                   - f_2 * nsf1_436[k]
                   + f_3 * pc_y[k] * nsg_655[k];

        t_919[k] = f_16 * msg_520[k]
                   + f_3 * pc_z[k] * nsg_655[k];

        t_920[k] = f_9 * msg_537[k]
                   + f_6 * nsf0_438[k]
                   - f_7 * nsf1_438[k]
                   + f_3 * pc_y[k] * nsg_657[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, pa_y, pc_y, msh0_755, msg_538, msg_539, \
                         msh1_755, nsf0_439, nsf1_439, nsg_658, \
                         nsg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_9 * msg_538[k]
                   + f_4 * nsf0_439[k]
                   - f_5 * nsf1_439[k]
                   + f_3 * pc_y[k] * nsg_658[k];

        t_922[k] = f_9 * msg_539[k]
                   + f_3 * pc_y[k] * nsg_659[k];

        t_923[k] = pa_y[k] * msh0_755[k]
                   - f_8 * pc_y[k] * msh1_755[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, msg_525, \
                         msg_660, nsf0_440, nsf1_440, nsg_660, nsg_661, \
                         nsg_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_10 * msg_660[k]
                   + f_1 * nsf0_440[k]
                   - f_2 * nsf1_440[k]
                   + f_3 * pc_x[k] * nsg_660[k];

        t_925[k] = f_3 * pc_y[k] * nsg_660[k];

        t_926[k] = f_15 * msg_525[k]
                   + f_3 * pc_z[k] * nsg_660[k];

        t_927[k] = f_4 * nsf0_440[k]
                   - f_5 * nsf1_440[k]
                   + f_3 * pc_y[k] * nsg_661[k];

        t_928[k] = f_3 * pc_y[k] * nsg_662[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, t_932, pc_x, pc_y, msg_665, nsf0_441, nsf0_442, \
                         nsf0_445, nsf1_441, nsf1_442, nsf1_445, nsg_663, nsg_664, \
                         nsg_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_10 * msg_665[k]
                   + f_6 * nsf0_445[k]
                   - f_7 * nsf1_445[k]
                   + f_3 * pc_x[k] * nsg_665[k];

        t_930[k] = f_6 * nsf0_441[k]
                   - f_7 * nsf1_441[k]
                   + f_3 * pc_y[k] * nsg_663[k];

        t_931[k] = f_4 * nsf0_442[k]
                   - f_5 * nsf1_442[k]
                   + f_3 * pc_y[k] * nsg_664[k];

        t_932[k] = f_3 * pc_y[k] * nsg_665[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, pc_x, msg_669, msg_670, msg_671, msg_672, \
                         nsf0_449, nsf1_449, nsg_669, nsg_670, nsg_671, \
                         nsg_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_10 * msg_669[k]
                   + f_4 * nsf0_449[k]
                   - f_5 * nsf1_449[k]
                   + f_3 * pc_x[k] * nsg_669[k];

        t_934[k] = f_10 * msg_670[k]
                   + f_3 * pc_x[k] * nsg_670[k];

        t_935[k] = f_10 * msg_671[k]
                   + f_3 * pc_x[k] * nsg_671[k];

        t_936[k] = f_10 * msg_672[k]
                   + f_3 * pc_x[k] * nsg_672[k];
    }

#pragma omp simd aligned(t_937, t_938, t_939, t_940, pc_x, pc_y, msg_674, nsf0_446, nsf0_447, \
                         nsf1_446, nsf1_447, nsg_669, nsg_670, nsg_671, \
                         nsg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_937[k] = f_3 * pc_y[k] * nsg_669[k];

        t_938[k] = f_10 * msg_674[k]
                   + f_3 * pc_x[k] * nsg_674[k];

        t_939[k] = f_1 * nsf0_446[k]
                   - f_2 * nsf1_446[k]
                   + f_3 * pc_y[k] * nsg_670[k];

        t_940[k] = f_13 * nsf0_447[k]
                   - f_14 * nsf1_447[k]
                   + f_3 * pc_y[k] * nsg_671[k];
    }

#pragma omp simd aligned(t_941, t_942, t_943, t_944, pc_y, pc_z, msg_539, nsf0_448, nsf0_449, \
                         nsf1_448, nsf1_449, nsg_672, nsg_673, \
                         nsg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_941[k] = f_6 * nsf0_448[k]
                   - f_7 * nsf1_448[k]
                   + f_3 * pc_y[k] * nsg_672[k];

        t_942[k] = f_4 * nsf0_449[k]
                   - f_5 * nsf1_449[k]
                   + f_3 * pc_y[k] * nsg_673[k];

        t_943[k] = f_3 * pc_y[k] * nsg_674[k];

        t_944[k] = f_15 * msg_539[k]
                   + f_1 * nsf0_449[k]
                   - f_2 * nsf1_449[k]
                   + f_3 * pc_z[k] * nsg_674[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.5 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_756 = buffer.data(msh0 + 756);
    const auto *msh0_759 = buffer.data(msh0 + 759);
    const auto *msh0_762 = buffer.data(msh0 + 762);
    const auto *msh0_945 = buffer.data(msh0 + 945);
    const auto *msh0_948 = buffer.data(msh0 + 948);
    const auto *msh0_951 = buffer.data(msh0 + 951);
    const auto *msh0_960 = buffer.data(msh0 + 960);
    const auto *msh0_962 = buffer.data(msh0 + 962);
    const auto *msh0_963 = buffer.data(msh0 + 963);
    const auto *msh0_965 = buffer.data(msh0 + 965);
    const auto *msh0_971 = buffer.data(msh0 + 971);
    const auto *msh0_975 = buffer.data(msh0 + 975);
    const auto *msh0_981 = buffer.data(msh0 + 981);
    const auto *msh0_983 = buffer.data(msh0 + 983);
    const auto *msh0_984 = buffer.data(msh0 + 984);
    const auto *msh0_986 = buffer.data(msh0 + 986);
    const auto *msh0_987 = buffer.data(msh0 + 987);
    const auto *msh0_990 = buffer.data(msh0 + 990);
    const auto *msh0_992 = buffer.data(msh0 + 992);
    const auto *msh0_993 = buffer.data(msh0 + 993);
    const auto *msh0_996 = buffer.data(msh0 + 996);
    const auto *msh0_1002 = buffer.data(msh0 + 1002);
    const auto *msh0_1004 = buffer.data(msh0 + 1004);
    const auto *msh0_1005 = buffer.data(msh0 + 1005);
    const auto *msh0_1007 = buffer.data(msh0 + 1007);
    const auto *msh0_1008 = buffer.data(msh0 + 1008);
    const auto *msh0_1011 = buffer.data(msh0 + 1011);
    const auto *msh0_1013 = buffer.data(msh0 + 1013);
    const auto *msh0_1014 = buffer.data(msh0 + 1014);
    const auto *msh0_1017 = buffer.data(msh0 + 1017);
    const auto *msh0_1023 = buffer.data(msh0 + 1023);
    const auto *msh0_1025 = buffer.data(msh0 + 1025);
    const auto *msh0_1026 = buffer.data(msh0 + 1026);
    const auto *msh0_1028 = buffer.data(msh0 + 1028);
    const auto *msh0_1029 = buffer.data(msh0 + 1029);
    const auto *msh0_1032 = buffer.data(msh0 + 1032);
    const auto *msh0_1034 = buffer.data(msh0 + 1034);
    const auto *msh0_1035 = buffer.data(msh0 + 1035);
    const auto *msh0_1038 = buffer.data(msh0 + 1038);
    const auto *msh0_1044 = buffer.data(msh0 + 1044);
    const auto *msh0_1046 = buffer.data(msh0 + 1046);
    const auto *msh0_1047 = buffer.data(msh0 + 1047);
    const auto *msh0_1049 = buffer.data(msh0 + 1049);
    const auto *msh0_1050 = buffer.data(msh0 + 1050);
    const auto *msh0_1053 = buffer.data(msh0 + 1053);
    const auto *msh0_1055 = buffer.data(msh0 + 1055);
    const auto *msh0_1056 = buffer.data(msh0 + 1056);
    const auto *msh0_1059 = buffer.data(msh0 + 1059);
    const auto *msh0_1065 = buffer.data(msh0 + 1065);

    const auto *msg_540 = buffer.data(msg + 540);
    const auto *msg_543 = buffer.data(msg + 543);
    const auto *msg_545 = buffer.data(msg + 545);
    const auto *msg_550 = buffer.data(msg + 550);
    const auto *msg_554 = buffer.data(msg + 554);
    const auto *msg_555 = buffer.data(msg + 555);
    const auto *msg_557 = buffer.data(msg + 557);
    const auto *msg_558 = buffer.data(msg + 558);
    const auto *msg_560 = buffer.data(msg + 560);
    const auto *msg_565 = buffer.data(msg + 565);
    const auto *msg_569 = buffer.data(msg + 569);
    const auto *msg_570 = buffer.data(msg + 570);
    const auto *msg_572 = buffer.data(msg + 572);
    const auto *msg_573 = buffer.data(msg + 573);
    const auto *msg_575 = buffer.data(msg + 575);
    const auto *msg_580 = buffer.data(msg + 580);
    const auto *msg_584 = buffer.data(msg + 584);
    const auto *msg_585 = buffer.data(msg + 585);
    const auto *msg_587 = buffer.data(msg + 587);
    const auto *msg_588 = buffer.data(msg + 588);
    const auto *msg_590 = buffer.data(msg + 590);
    const auto *msg_595 = buffer.data(msg + 595);
    const auto *msg_599 = buffer.data(msg + 599);
    const auto *msg_600 = buffer.data(msg + 600);
    const auto *msg_602 = buffer.data(msg + 602);
    const auto *msg_603 = buffer.data(msg + 603);
    const auto *msg_605 = buffer.data(msg + 605);
    const auto *msg_610 = buffer.data(msg + 610);
    const auto *msg_614 = buffer.data(msg + 614);
    const auto *msg_615 = buffer.data(msg + 615);
    const auto *msg_617 = buffer.data(msg + 617);
    const auto *msg_620 = buffer.data(msg + 620);
    const auto *msg_675 = buffer.data(msg + 675);
    const auto *msg_678 = buffer.data(msg + 678);
    const auto *msg_681 = buffer.data(msg + 681);
    const auto *msg_685 = buffer.data(msg + 685);
    const auto *msg_687 = buffer.data(msg + 687);
    const auto *msg_688 = buffer.data(msg + 688);
    const auto *msg_689 = buffer.data(msg + 689);
    const auto *msg_695 = buffer.data(msg + 695);
    const auto *msg_699 = buffer.data(msg + 699);
    const auto *msg_700 = buffer.data(msg + 700);
    const auto *msg_701 = buffer.data(msg + 701);
    const auto *msg_702 = buffer.data(msg + 702);
    const auto *msg_703 = buffer.data(msg + 703);
    const auto *msg_704 = buffer.data(msg + 704);
    const auto *msg_705 = buffer.data(msg + 705);
    const auto *msg_708 = buffer.data(msg + 708);
    const auto *msg_710 = buffer.data(msg + 710);
    const auto *msg_711 = buffer.data(msg + 711);
    const auto *msg_714 = buffer.data(msg + 714);
    const auto *msg_715 = buffer.data(msg + 715);
    const auto *msg_716 = buffer.data(msg + 716);
    const auto *msg_717 = buffer.data(msg + 717);
    const auto *msg_718 = buffer.data(msg + 718);
    const auto *msg_719 = buffer.data(msg + 719);
    const auto *msg_720 = buffer.data(msg + 720);
    const auto *msg_723 = buffer.data(msg + 723);
    const auto *msg_725 = buffer.data(msg + 725);
    const auto *msg_726 = buffer.data(msg + 726);
    const auto *msg_729 = buffer.data(msg + 729);
    const auto *msg_730 = buffer.data(msg + 730);
    const auto *msg_731 = buffer.data(msg + 731);
    const auto *msg_732 = buffer.data(msg + 732);
    const auto *msg_733 = buffer.data(msg + 733);
    const auto *msg_734 = buffer.data(msg + 734);
    const auto *msg_735 = buffer.data(msg + 735);
    const auto *msg_738 = buffer.data(msg + 738);
    const auto *msg_740 = buffer.data(msg + 740);
    const auto *msg_741 = buffer.data(msg + 741);
    const auto *msg_744 = buffer.data(msg + 744);
    const auto *msg_745 = buffer.data(msg + 745);
    const auto *msg_746 = buffer.data(msg + 746);
    const auto *msg_747 = buffer.data(msg + 747);
    const auto *msg_748 = buffer.data(msg + 748);
    const auto *msg_749 = buffer.data(msg + 749);
    const auto *msg_750 = buffer.data(msg + 750);
    const auto *msg_753 = buffer.data(msg + 753);
    const auto *msg_755 = buffer.data(msg + 755);
    const auto *msg_756 = buffer.data(msg + 756);
    const auto *msg_759 = buffer.data(msg + 759);
    const auto *msg_760 = buffer.data(msg + 760);
    const auto *msg_761 = buffer.data(msg + 761);
    const auto *msg_762 = buffer.data(msg + 762);
    const auto *msg_763 = buffer.data(msg + 763);
    const auto *msg_764 = buffer.data(msg + 764);

    const auto *msh1_756 = buffer.data(msh1 + 756);
    const auto *msh1_759 = buffer.data(msh1 + 759);
    const auto *msh1_762 = buffer.data(msh1 + 762);
    const auto *msh1_945 = buffer.data(msh1 + 945);
    const auto *msh1_948 = buffer.data(msh1 + 948);
    const auto *msh1_951 = buffer.data(msh1 + 951);
    const auto *msh1_960 = buffer.data(msh1 + 960);
    const auto *msh1_962 = buffer.data(msh1 + 962);
    const auto *msh1_963 = buffer.data(msh1 + 963);
    const auto *msh1_965 = buffer.data(msh1 + 965);
    const auto *msh1_971 = buffer.data(msh1 + 971);
    const auto *msh1_975 = buffer.data(msh1 + 975);
    const auto *msh1_981 = buffer.data(msh1 + 981);
    const auto *msh1_983 = buffer.data(msh1 + 983);
    const auto *msh1_984 = buffer.data(msh1 + 984);
    const auto *msh1_986 = buffer.data(msh1 + 986);
    const auto *msh1_987 = buffer.data(msh1 + 987);
    const auto *msh1_990 = buffer.data(msh1 + 990);
    const auto *msh1_992 = buffer.data(msh1 + 992);
    const auto *msh1_993 = buffer.data(msh1 + 993);
    const auto *msh1_996 = buffer.data(msh1 + 996);
    const auto *msh1_1002 = buffer.data(msh1 + 1002);
    const auto *msh1_1004 = buffer.data(msh1 + 1004);
    const auto *msh1_1005 = buffer.data(msh1 + 1005);
    const auto *msh1_1007 = buffer.data(msh1 + 1007);
    const auto *msh1_1008 = buffer.data(msh1 + 1008);
    const auto *msh1_1011 = buffer.data(msh1 + 1011);
    const auto *msh1_1013 = buffer.data(msh1 + 1013);
    const auto *msh1_1014 = buffer.data(msh1 + 1014);
    const auto *msh1_1017 = buffer.data(msh1 + 1017);
    const auto *msh1_1023 = buffer.data(msh1 + 1023);
    const auto *msh1_1025 = buffer.data(msh1 + 1025);
    const auto *msh1_1026 = buffer.data(msh1 + 1026);
    const auto *msh1_1028 = buffer.data(msh1 + 1028);
    const auto *msh1_1029 = buffer.data(msh1 + 1029);
    const auto *msh1_1032 = buffer.data(msh1 + 1032);
    const auto *msh1_1034 = buffer.data(msh1 + 1034);
    const auto *msh1_1035 = buffer.data(msh1 + 1035);
    const auto *msh1_1038 = buffer.data(msh1 + 1038);
    const auto *msh1_1044 = buffer.data(msh1 + 1044);
    const auto *msh1_1046 = buffer.data(msh1 + 1046);
    const auto *msh1_1047 = buffer.data(msh1 + 1047);
    const auto *msh1_1049 = buffer.data(msh1 + 1049);
    const auto *msh1_1050 = buffer.data(msh1 + 1050);
    const auto *msh1_1053 = buffer.data(msh1 + 1053);
    const auto *msh1_1055 = buffer.data(msh1 + 1055);
    const auto *msh1_1056 = buffer.data(msh1 + 1056);
    const auto *msh1_1059 = buffer.data(msh1 + 1059);
    const auto *msh1_1065 = buffer.data(msh1 + 1065);

    const auto *nsf0_450 = buffer.data(nsf0 + 450);
    const auto *nsf0_452 = buffer.data(nsf0 + 452);

    const auto *nsf1_450 = buffer.data(nsf1 + 450);
    const auto *nsf1_452 = buffer.data(nsf1 + 452);

    const auto *nsg_675 = buffer.data(nsg + 675);
    const auto *nsg_676 = buffer.data(nsg + 676);
    const auto *nsg_677 = buffer.data(nsg + 677);
    const auto *nsg_678 = buffer.data(nsg + 678);
    const auto *nsg_680 = buffer.data(nsg + 680);
    const auto *nsg_681 = buffer.data(nsg + 681);
    const auto *nsg_685 = buffer.data(nsg + 685);
    const auto *nsg_687 = buffer.data(nsg + 687);
    const auto *nsg_688 = buffer.data(nsg + 688);
    const auto *nsg_689 = buffer.data(nsg + 689);
    const auto *nsg_690 = buffer.data(nsg + 690);
    const auto *nsg_692 = buffer.data(nsg + 692);
    const auto *nsg_693 = buffer.data(nsg + 693);
    const auto *nsg_695 = buffer.data(nsg + 695);
    const auto *nsg_700 = buffer.data(nsg + 700);
    const auto *nsg_701 = buffer.data(nsg + 701);
    const auto *nsg_702 = buffer.data(nsg + 702);
    const auto *nsg_703 = buffer.data(nsg + 703);
    const auto *nsg_704 = buffer.data(nsg + 704);
    const auto *nsg_705 = buffer.data(nsg + 705);
    const auto *nsg_707 = buffer.data(nsg + 707);
    const auto *nsg_708 = buffer.data(nsg + 708);
    const auto *nsg_710 = buffer.data(nsg + 710);
    const auto *nsg_715 = buffer.data(nsg + 715);
    const auto *nsg_716 = buffer.data(nsg + 716);
    const auto *nsg_717 = buffer.data(nsg + 717);
    const auto *nsg_718 = buffer.data(nsg + 718);
    const auto *nsg_719 = buffer.data(nsg + 719);
    const auto *nsg_720 = buffer.data(nsg + 720);
    const auto *nsg_722 = buffer.data(nsg + 722);
    const auto *nsg_723 = buffer.data(nsg + 723);
    const auto *nsg_725 = buffer.data(nsg + 725);
    const auto *nsg_730 = buffer.data(nsg + 730);
    const auto *nsg_731 = buffer.data(nsg + 731);
    const auto *nsg_732 = buffer.data(nsg + 732);
    const auto *nsg_733 = buffer.data(nsg + 733);
    const auto *nsg_734 = buffer.data(nsg + 734);
    const auto *nsg_735 = buffer.data(nsg + 735);
    const auto *nsg_737 = buffer.data(nsg + 737);
    const auto *nsg_738 = buffer.data(nsg + 738);
    const auto *nsg_740 = buffer.data(nsg + 740);
    const auto *nsg_745 = buffer.data(nsg + 745);
    const auto *nsg_746 = buffer.data(nsg + 746);
    const auto *nsg_747 = buffer.data(nsg + 747);
    const auto *nsg_748 = buffer.data(nsg + 748);
    const auto *nsg_749 = buffer.data(nsg + 749);
    const auto *nsg_750 = buffer.data(nsg + 750);
    const auto *nsg_752 = buffer.data(nsg + 752);
    const auto *nsg_753 = buffer.data(nsg + 753);
    const auto *nsg_755 = buffer.data(nsg + 755);
    const auto *nsg_760 = buffer.data(nsg + 760);
    const auto *nsg_761 = buffer.data(nsg + 761);
    const auto *nsg_762 = buffer.data(nsg + 762);
    const auto *nsg_763 = buffer.data(nsg + 763);
    const auto *nsg_764 = buffer.data(nsg + 764);

#pragma omp simd aligned(t_945, t_946, t_947, t_948, pa_x, pc_x, pc_y, pc_z, msh0_945, \
                         msh0_948, msg_540, msg_675, msg_678, msh1_945, msh1_948, \
                         nsg_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = pa_x[k] * msh0_945[k]
                   + f_19 * msg_675[k]
                   - f_8 * pc_x[k] * msh1_945[k];

        t_946[k] = f_12 * msg_540[k]
                   + f_3 * pc_y[k] * nsg_675[k];

        t_947[k] = f_3 * pc_z[k] * nsg_675[k];

        t_948[k] = pa_x[k] * msh0_948[k]
                   + f_11 * msg_678[k]
                   - f_8 * pc_x[k] * msh1_948[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pa_x, pc_x, pc_z, msh0_951, msg_681, \
                         msh1_951, nsf0_450, nsf1_450, nsg_676, nsg_677, \
                         nsg_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_3 * pc_z[k] * nsg_676[k];

        t_950[k] = f_4 * nsf0_450[k]
                   - f_5 * nsf1_450[k]
                   + f_3 * pc_z[k] * nsg_677[k];

        t_951[k] = pa_x[k] * msh0_951[k]
                   + f_10 * msg_681[k]
                   - f_8 * pc_x[k] * msh1_951[k];

        t_952[k] = f_3 * pc_z[k] * nsg_678[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pc_x, pc_y, pc_z, msg_545, msg_685, \
                         nsf0_452, nsf1_452, nsg_680, nsg_681, \
                         nsg_685 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_12 * msg_545[k]
                   + f_3 * pc_y[k] * nsg_680[k];

        t_954[k] = f_6 * nsf0_452[k]
                   - f_7 * nsf1_452[k]
                   + f_3 * pc_z[k] * nsg_680[k];

        t_955[k] = f_9 * msg_685[k]
                   + f_3 * pc_x[k] * nsg_685[k];

        t_956[k] = f_3 * pc_z[k] * nsg_681[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pa_x, pc_x, msh0_960, msg_687, msg_688, \
                         msg_689, msh1_960, nsg_687, nsg_688, nsg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_9 * msg_687[k]
                   + f_3 * pc_x[k] * nsg_687[k];

        t_958[k] = f_9 * msg_688[k]
                   + f_3 * pc_x[k] * nsg_688[k];

        t_959[k] = f_9 * msg_689[k]
                   + f_3 * pc_x[k] * nsg_689[k];

        t_960[k] = pa_x[k] * msh0_960[k]
                   - f_8 * pc_x[k] * msh1_960[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, t_964, pa_x, pc_x, pc_y, pc_z, msh0_962, \
                         msh0_963, msg_554, msh1_962, msh1_963, nsg_685, \
                         nsg_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_3 * pc_z[k] * nsg_685[k];

        t_962[k] = pa_x[k] * msh0_962[k]
                   - f_8 * pc_x[k] * msh1_962[k];

        t_963[k] = pa_x[k] * msh0_963[k]
                   - f_8 * pc_x[k] * msh1_963[k];

        t_964[k] = f_12 * msg_554[k]
                   + f_3 * pc_y[k] * nsg_689[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, pa_x, pa_z, pc_x, pc_y, pc_z, msh0_756, \
                         msh0_965, msg_540, msg_555, msh1_756, msh1_965, \
                         nsg_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = pa_x[k] * msh0_965[k]
                   - f_8 * pc_x[k] * msh1_965[k];

        t_966[k] = pa_z[k] * msh0_756[k]
                   - f_8 * pc_z[k] * msh1_756[k];

        t_967[k] = f_15 * msg_555[k]
                   + f_3 * pc_y[k] * nsg_690[k];

        t_968[k] = f_9 * msg_540[k]
                   + f_3 * pc_z[k] * nsg_690[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pa_x, pa_z, pc_x, pc_y, pc_z, msh0_759, \
                         msh0_971, msg_557, msg_695, msh1_759, msh1_971, \
                         nsg_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = pa_z[k] * msh0_759[k]
                   - f_8 * pc_z[k] * msh1_759[k];

        t_970[k] = f_15 * msg_557[k]
                   + f_3 * pc_y[k] * nsg_692[k];

        t_971[k] = pa_x[k] * msh0_971[k]
                   + f_11 * msg_695[k]
                   - f_8 * pc_x[k] * msh1_971[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, pa_z, pc_y, pc_z, msh0_762, msg_543, msg_560, \
                         msh1_762, nsg_693, nsg_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = pa_z[k] * msh0_762[k]
                   - f_8 * pc_z[k] * msh1_762[k];

        t_973[k] = f_9 * msg_543[k]
                   + f_3 * pc_z[k] * nsg_693[k];

        t_974[k] = f_15 * msg_560[k]
                   + f_3 * pc_y[k] * nsg_695[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, pa_x, pc_x, msh0_975, msg_699, msg_700, \
                         msg_701, msg_702, msh1_975, nsg_700, nsg_701, \
                         nsg_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = pa_x[k] * msh0_975[k]
                   + f_10 * msg_699[k]
                   - f_8 * pc_x[k] * msh1_975[k];

        t_976[k] = f_9 * msg_700[k]
                   + f_3 * pc_x[k] * nsg_700[k];

        t_977[k] = f_9 * msg_701[k]
                   + f_3 * pc_x[k] * nsg_701[k];

        t_978[k] = f_9 * msg_702[k]
                   + f_3 * pc_x[k] * nsg_702[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pa_x, pc_x, pc_z, msh0_981, msg_550, \
                         msg_703, msg_704, msh1_981, nsg_700, nsg_703, \
                         nsg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = f_9 * msg_703[k]
                   + f_3 * pc_x[k] * nsg_703[k];

        t_980[k] = f_9 * msg_704[k]
                   + f_3 * pc_x[k] * nsg_704[k];

        t_981[k] = pa_x[k] * msh0_981[k]
                   - f_8 * pc_x[k] * msh1_981[k];

        t_982[k] = f_9 * msg_550[k]
                   + f_3 * pc_z[k] * nsg_700[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, t_986, pa_x, pc_x, pc_y, msh0_983, msh0_984, \
                         msh0_986, msg_569, msh1_983, msh1_984, msh1_986, \
                         nsg_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = pa_x[k] * msh0_983[k]
                   - f_8 * pc_x[k] * msh1_983[k];

        t_984[k] = pa_x[k] * msh0_984[k]
                   - f_8 * pc_x[k] * msh1_984[k];

        t_985[k] = f_15 * msg_569[k]
                   + f_3 * pc_y[k] * nsg_704[k];

        t_986[k] = pa_x[k] * msh0_986[k]
                   - f_8 * pc_x[k] * msh1_986[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, pa_x, pc_x, pc_y, pc_z, msh0_987, msg_555, \
                         msg_570, msg_705, msh1_987, nsg_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = pa_x[k] * msh0_987[k]
                   + f_19 * msg_705[k]
                   - f_8 * pc_x[k] * msh1_987[k];

        t_988[k] = f_16 * msg_570[k]
                   + f_3 * pc_y[k] * nsg_705[k];

        t_989[k] = f_10 * msg_555[k]
                   + f_3 * pc_z[k] * nsg_705[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, pa_x, pc_x, pc_y, msh0_990, msh0_992, msg_572, \
                         msg_708, msg_710, msh1_990, msh1_992, \
                         nsg_707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = pa_x[k] * msh0_990[k]
                   + f_11 * msg_708[k]
                   - f_8 * pc_x[k] * msh1_990[k];

        t_991[k] = f_16 * msg_572[k]
                   + f_3 * pc_y[k] * nsg_707[k];

        t_992[k] = pa_x[k] * msh0_992[k]
                   + f_11 * msg_710[k]
                   - f_8 * pc_x[k] * msh1_992[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, pa_x, pc_x, pc_y, pc_z, msh0_993, msg_558, \
                         msg_575, msg_711, msh1_993, nsg_708, nsg_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = pa_x[k] * msh0_993[k]
                   + f_10 * msg_711[k]
                   - f_8 * pc_x[k] * msh1_993[k];

        t_994[k] = f_10 * msg_558[k]
                   + f_3 * pc_z[k] * nsg_708[k];

        t_995[k] = f_16 * msg_575[k]
                   + f_3 * pc_y[k] * nsg_710[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, pa_x, pc_x, msh0_996, msg_714, msg_715, \
                         msg_716, msg_717, msh1_996, nsg_715, nsg_716, \
                         nsg_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = pa_x[k] * msh0_996[k]
                   + f_10 * msg_714[k]
                   - f_8 * pc_x[k] * msh1_996[k];

        t_997[k] = f_9 * msg_715[k]
                   + f_3 * pc_x[k] * nsg_715[k];

        t_998[k] = f_9 * msg_716[k]
                   + f_3 * pc_x[k] * nsg_716[k];

        t_999[k] = f_9 * msg_717[k]
                   + f_3 * pc_x[k] * nsg_717[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pa_x, pc_x, pc_z, msh0_1002, msg_565, \
                         msg_718, msg_719, msh1_1002, nsg_715, nsg_718, \
                         nsg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_9 * msg_718[k]
                    + f_3 * pc_x[k] * nsg_718[k];

        t_1001[k] = f_9 * msg_719[k]
                    + f_3 * pc_x[k] * nsg_719[k];

        t_1002[k] = pa_x[k] * msh0_1002[k]
                    - f_8 * pc_x[k] * msh1_1002[k];

        t_1003[k] = f_10 * msg_565[k]
                    + f_3 * pc_z[k] * nsg_715[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pa_x, pc_x, pc_y, msh0_1004, \
                         msh0_1005, msh0_1007, msg_584, msh1_1004, msh1_1005, msh1_1007, \
                         nsg_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = pa_x[k] * msh0_1004[k]
                    - f_8 * pc_x[k] * msh1_1004[k];

        t_1005[k] = pa_x[k] * msh0_1005[k]
                    - f_8 * pc_x[k] * msh1_1005[k];

        t_1006[k] = f_16 * msg_584[k]
                    + f_3 * pc_y[k] * nsg_719[k];

        t_1007[k] = pa_x[k] * msh0_1007[k]
                    - f_8 * pc_x[k] * msh1_1007[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, pa_x, pc_x, pc_y, pc_z, msh0_1008, msg_570, \
                         msg_585, msg_720, msh1_1008, nsg_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = pa_x[k] * msh0_1008[k]
                    + f_19 * msg_720[k]
                    - f_8 * pc_x[k] * msh1_1008[k];

        t_1009[k] = f_17 * msg_585[k]
                    + f_3 * pc_y[k] * nsg_720[k];

        t_1010[k] = f_11 * msg_570[k]
                    + f_3 * pc_z[k] * nsg_720[k];
    }

#pragma omp simd aligned(t_1011, t_1012, t_1013, pa_x, pc_x, pc_y, msh0_1011, msh0_1013, \
                         msg_587, msg_723, msg_725, msh1_1011, msh1_1013, \
                         nsg_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1011[k] = pa_x[k] * msh0_1011[k]
                    + f_11 * msg_723[k]
                    - f_8 * pc_x[k] * msh1_1011[k];

        t_1012[k] = f_17 * msg_587[k]
                    + f_3 * pc_y[k] * nsg_722[k];

        t_1013[k] = pa_x[k] * msh0_1013[k]
                    + f_11 * msg_725[k]
                    - f_8 * pc_x[k] * msh1_1013[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, pa_x, pc_x, pc_y, pc_z, msh0_1014, msg_573, \
                         msg_590, msg_726, msh1_1014, nsg_723, \
                         nsg_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = pa_x[k] * msh0_1014[k]
                    + f_10 * msg_726[k]
                    - f_8 * pc_x[k] * msh1_1014[k];

        t_1015[k] = f_11 * msg_573[k]
                    + f_3 * pc_z[k] * nsg_723[k];

        t_1016[k] = f_17 * msg_590[k]
                    + f_3 * pc_y[k] * nsg_725[k];
    }

#pragma omp simd aligned(t_1017, t_1018, t_1019, t_1020, pa_x, pc_x, msh0_1017, msg_729, \
                         msg_730, msg_731, msg_732, msh1_1017, nsg_730, nsg_731, \
                         nsg_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1017[k] = pa_x[k] * msh0_1017[k]
                    + f_10 * msg_729[k]
                    - f_8 * pc_x[k] * msh1_1017[k];

        t_1018[k] = f_9 * msg_730[k]
                    + f_3 * pc_x[k] * nsg_730[k];

        t_1019[k] = f_9 * msg_731[k]
                    + f_3 * pc_x[k] * nsg_731[k];

        t_1020[k] = f_9 * msg_732[k]
                    + f_3 * pc_x[k] * nsg_732[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, pa_x, pc_x, pc_z, msh0_1023, msg_580, \
                         msg_733, msg_734, msh1_1023, nsg_730, nsg_733, \
                         nsg_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_9 * msg_733[k]
                    + f_3 * pc_x[k] * nsg_733[k];

        t_1022[k] = f_9 * msg_734[k]
                    + f_3 * pc_x[k] * nsg_734[k];

        t_1023[k] = pa_x[k] * msh0_1023[k]
                    - f_8 * pc_x[k] * msh1_1023[k];

        t_1024[k] = f_11 * msg_580[k]
                    + f_3 * pc_z[k] * nsg_730[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, pa_x, pc_x, pc_y, msh0_1025, \
                         msh0_1026, msh0_1028, msg_599, msh1_1025, msh1_1026, msh1_1028, \
                         nsg_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = pa_x[k] * msh0_1025[k]
                    - f_8 * pc_x[k] * msh1_1025[k];

        t_1026[k] = pa_x[k] * msh0_1026[k]
                    - f_8 * pc_x[k] * msh1_1026[k];

        t_1027[k] = f_17 * msg_599[k]
                    + f_3 * pc_y[k] * nsg_734[k];

        t_1028[k] = pa_x[k] * msh0_1028[k]
                    - f_8 * pc_x[k] * msh1_1028[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, pa_x, pc_x, pc_y, pc_z, msh0_1029, msg_585, \
                         msg_600, msg_735, msh1_1029, nsg_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = pa_x[k] * msh0_1029[k]
                    + f_19 * msg_735[k]
                    - f_8 * pc_x[k] * msh1_1029[k];

        t_1030[k] = f_19 * msg_600[k]
                    + f_3 * pc_y[k] * nsg_735[k];

        t_1031[k] = f_18 * msg_585[k]
                    + f_3 * pc_z[k] * nsg_735[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, pa_x, pc_x, pc_y, msh0_1032, msh0_1034, \
                         msg_602, msg_738, msg_740, msh1_1032, msh1_1034, \
                         nsg_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = pa_x[k] * msh0_1032[k]
                    + f_11 * msg_738[k]
                    - f_8 * pc_x[k] * msh1_1032[k];

        t_1033[k] = f_19 * msg_602[k]
                    + f_3 * pc_y[k] * nsg_737[k];

        t_1034[k] = pa_x[k] * msh0_1034[k]
                    + f_11 * msg_740[k]
                    - f_8 * pc_x[k] * msh1_1034[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, pa_x, pc_x, pc_y, pc_z, msh0_1035, msg_588, \
                         msg_605, msg_741, msh1_1035, nsg_738, \
                         nsg_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = pa_x[k] * msh0_1035[k]
                    + f_10 * msg_741[k]
                    - f_8 * pc_x[k] * msh1_1035[k];

        t_1036[k] = f_18 * msg_588[k]
                    + f_3 * pc_z[k] * nsg_738[k];

        t_1037[k] = f_19 * msg_605[k]
                    + f_3 * pc_y[k] * nsg_740[k];
    }

#pragma omp simd aligned(t_1038, t_1039, t_1040, t_1041, pa_x, pc_x, msh0_1038, msg_744, \
                         msg_745, msg_746, msg_747, msh1_1038, nsg_745, nsg_746, \
                         nsg_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1038[k] = pa_x[k] * msh0_1038[k]
                    + f_10 * msg_744[k]
                    - f_8 * pc_x[k] * msh1_1038[k];

        t_1039[k] = f_9 * msg_745[k]
                    + f_3 * pc_x[k] * nsg_745[k];

        t_1040[k] = f_9 * msg_746[k]
                    + f_3 * pc_x[k] * nsg_746[k];

        t_1041[k] = f_9 * msg_747[k]
                    + f_3 * pc_x[k] * nsg_747[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, pa_x, pc_x, pc_z, msh0_1044, msg_595, \
                         msg_748, msg_749, msh1_1044, nsg_745, nsg_748, \
                         nsg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_9 * msg_748[k]
                    + f_3 * pc_x[k] * nsg_748[k];

        t_1043[k] = f_9 * msg_749[k]
                    + f_3 * pc_x[k] * nsg_749[k];

        t_1044[k] = pa_x[k] * msh0_1044[k]
                    - f_8 * pc_x[k] * msh1_1044[k];

        t_1045[k] = f_18 * msg_595[k]
                    + f_3 * pc_z[k] * nsg_745[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, t_1049, pa_x, pc_x, pc_y, msh0_1046, \
                         msh0_1047, msh0_1049, msg_614, msh1_1046, msh1_1047, msh1_1049, \
                         nsg_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = pa_x[k] * msh0_1046[k]
                    - f_8 * pc_x[k] * msh1_1046[k];

        t_1047[k] = pa_x[k] * msh0_1047[k]
                    - f_8 * pc_x[k] * msh1_1047[k];

        t_1048[k] = f_19 * msg_614[k]
                    + f_3 * pc_y[k] * nsg_749[k];

        t_1049[k] = pa_x[k] * msh0_1049[k]
                    - f_8 * pc_x[k] * msh1_1049[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, pa_x, pc_x, pc_y, pc_z, msh0_1050, msg_600, \
                         msg_615, msg_750, msh1_1050, nsg_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = pa_x[k] * msh0_1050[k]
                    + f_19 * msg_750[k]
                    - f_8 * pc_x[k] * msh1_1050[k];

        t_1051[k] = f_18 * msg_615[k]
                    + f_3 * pc_y[k] * nsg_750[k];

        t_1052[k] = f_19 * msg_600[k]
                    + f_3 * pc_z[k] * nsg_750[k];
    }

#pragma omp simd aligned(t_1053, t_1054, t_1055, pa_x, pc_x, pc_y, msh0_1053, msh0_1055, \
                         msg_617, msg_753, msg_755, msh1_1053, msh1_1055, \
                         nsg_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1053[k] = pa_x[k] * msh0_1053[k]
                    + f_11 * msg_753[k]
                    - f_8 * pc_x[k] * msh1_1053[k];

        t_1054[k] = f_18 * msg_617[k]
                    + f_3 * pc_y[k] * nsg_752[k];

        t_1055[k] = pa_x[k] * msh0_1055[k]
                    + f_11 * msg_755[k]
                    - f_8 * pc_x[k] * msh1_1055[k];
    }

#pragma omp simd aligned(t_1056, t_1057, t_1058, pa_x, pc_x, pc_y, pc_z, msh0_1056, msg_603, \
                         msg_620, msg_756, msh1_1056, nsg_753, \
                         nsg_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1056[k] = pa_x[k] * msh0_1056[k]
                    + f_10 * msg_756[k]
                    - f_8 * pc_x[k] * msh1_1056[k];

        t_1057[k] = f_19 * msg_603[k]
                    + f_3 * pc_z[k] * nsg_753[k];

        t_1058[k] = f_18 * msg_620[k]
                    + f_3 * pc_y[k] * nsg_755[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, t_1062, pa_x, pc_x, msh0_1059, msg_759, \
                         msg_760, msg_761, msg_762, msh1_1059, nsg_760, nsg_761, \
                         nsg_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = pa_x[k] * msh0_1059[k]
                    + f_10 * msg_759[k]
                    - f_8 * pc_x[k] * msh1_1059[k];

        t_1060[k] = f_9 * msg_760[k]
                    + f_3 * pc_x[k] * nsg_760[k];

        t_1061[k] = f_9 * msg_761[k]
                    + f_3 * pc_x[k] * nsg_761[k];

        t_1062[k] = f_9 * msg_762[k]
                    + f_3 * pc_x[k] * nsg_762[k];
    }

#pragma omp simd aligned(t_1063, t_1064, t_1065, t_1066, pa_x, pc_x, pc_z, msh0_1065, msg_610, \
                         msg_763, msg_764, msh1_1065, nsg_760, nsg_763, \
                         nsg_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1063[k] = f_9 * msg_763[k]
                    + f_3 * pc_x[k] * nsg_763[k];

        t_1064[k] = f_9 * msg_764[k]
                    + f_3 * pc_x[k] * nsg_764[k];

        t_1065[k] = pa_x[k] * msh0_1065[k]
                    - f_8 * pc_x[k] * msh1_1065[k];

        t_1066[k] = f_19 * msg_610[k]
                    + f_3 * pc_z[k] * nsg_760[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msh0,
                                                          const size_t msg, const size_t msh1,
                                                          const size_t nsf0, const size_t nsf1,
                                                          const size_t nsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
    auto *t_1182 = buffer.data(target + 1182);
    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_924 = buffer.data(msh0 + 924);
    const auto *msh0_929 = buffer.data(msh0 + 929);
    const auto *msh0_933 = buffer.data(msh0 + 933);
    const auto *msh0_945 = buffer.data(msh0 + 945);
    const auto *msh0_946 = buffer.data(msh0 + 946);
    const auto *msh0_948 = buffer.data(msh0 + 948);
    const auto *msh0_951 = buffer.data(msh0 + 951);
    const auto *msh0_960 = buffer.data(msh0 + 960);
    const auto *msh0_1067 = buffer.data(msh0 + 1067);
    const auto *msh0_1068 = buffer.data(msh0 + 1068);
    const auto *msh0_1070 = buffer.data(msh0 + 1070);
    const auto *msh0_1071 = buffer.data(msh0 + 1071);
    const auto *msh0_1074 = buffer.data(msh0 + 1074);
    const auto *msh0_1076 = buffer.data(msh0 + 1076);
    const auto *msh0_1077 = buffer.data(msh0 + 1077);
    const auto *msh0_1080 = buffer.data(msh0 + 1080);
    const auto *msh0_1086 = buffer.data(msh0 + 1086);
    const auto *msh0_1088 = buffer.data(msh0 + 1088);
    const auto *msh0_1089 = buffer.data(msh0 + 1089);
    const auto *msh0_1091 = buffer.data(msh0 + 1091);
    const auto *msh0_1092 = buffer.data(msh0 + 1092);
    const auto *msh0_1095 = buffer.data(msh0 + 1095);
    const auto *msh0_1097 = buffer.data(msh0 + 1097);
    const auto *msh0_1098 = buffer.data(msh0 + 1098);
    const auto *msh0_1101 = buffer.data(msh0 + 1101);
    const auto *msh0_1107 = buffer.data(msh0 + 1107);
    const auto *msh0_1109 = buffer.data(msh0 + 1109);
    const auto *msh0_1110 = buffer.data(msh0 + 1110);
    const auto *msh0_1112 = buffer.data(msh0 + 1112);
    const auto *msh0_1116 = buffer.data(msh0 + 1116);
    const auto *msh0_1119 = buffer.data(msh0 + 1119);
    const auto *msh0_1128 = buffer.data(msh0 + 1128);
    const auto *msh0_1130 = buffer.data(msh0 + 1130);
    const auto *msh0_1131 = buffer.data(msh0 + 1131);
    const auto *msh0_1133 = buffer.data(msh0 + 1133);
    const auto *msh0_1134 = buffer.data(msh0 + 1134);
    const auto *msh0_1139 = buffer.data(msh0 + 1139);
    const auto *msh0_1143 = buffer.data(msh0 + 1143);
    const auto *msh0_1149 = buffer.data(msh0 + 1149);
    const auto *msh0_1150 = buffer.data(msh0 + 1150);
    const auto *msh0_1151 = buffer.data(msh0 + 1151);
    const auto *msh0_1152 = buffer.data(msh0 + 1152);
    const auto *msh0_1154 = buffer.data(msh0 + 1154);

    const auto *msg_615 = buffer.data(msg + 615);
    const auto *msg_618 = buffer.data(msg + 618);
    const auto *msg_625 = buffer.data(msg + 625);
    const auto *msg_629 = buffer.data(msg + 629);
    const auto *msg_630 = buffer.data(msg + 630);
    const auto *msg_632 = buffer.data(msg + 632);
    const auto *msg_633 = buffer.data(msg + 633);
    const auto *msg_635 = buffer.data(msg + 635);
    const auto *msg_640 = buffer.data(msg + 640);
    const auto *msg_644 = buffer.data(msg + 644);
    const auto *msg_645 = buffer.data(msg + 645);
    const auto *msg_647 = buffer.data(msg + 647);
    const auto *msg_648 = buffer.data(msg + 648);
    const auto *msg_650 = buffer.data(msg + 650);
    const auto *msg_655 = buffer.data(msg + 655);
    const auto *msg_659 = buffer.data(msg + 659);
    const auto *msg_660 = buffer.data(msg + 660);
    const auto *msg_662 = buffer.data(msg + 662);
    const auto *msg_665 = buffer.data(msg + 665);
    const auto *msg_674 = buffer.data(msg + 674);
    const auto *msg_685 = buffer.data(msg + 685);
    const auto *msg_689 = buffer.data(msg + 689);
    const auto *msg_765 = buffer.data(msg + 765);
    const auto *msg_768 = buffer.data(msg + 768);
    const auto *msg_770 = buffer.data(msg + 770);
    const auto *msg_771 = buffer.data(msg + 771);
    const auto *msg_774 = buffer.data(msg + 774);
    const auto *msg_775 = buffer.data(msg + 775);
    const auto *msg_776 = buffer.data(msg + 776);
    const auto *msg_777 = buffer.data(msg + 777);
    const auto *msg_778 = buffer.data(msg + 778);
    const auto *msg_779 = buffer.data(msg + 779);
    const auto *msg_780 = buffer.data(msg + 780);
    const auto *msg_783 = buffer.data(msg + 783);
    const auto *msg_785 = buffer.data(msg + 785);
    const auto *msg_786 = buffer.data(msg + 786);
    const auto *msg_789 = buffer.data(msg + 789);
    const auto *msg_790 = buffer.data(msg + 790);
    const auto *msg_791 = buffer.data(msg + 791);
    const auto *msg_792 = buffer.data(msg + 792);
    const auto *msg_793 = buffer.data(msg + 793);
    const auto *msg_794 = buffer.data(msg + 794);
    const auto *msg_798 = buffer.data(msg + 798);
    const auto *msg_801 = buffer.data(msg + 801);
    const auto *msg_805 = buffer.data(msg + 805);
    const auto *msg_806 = buffer.data(msg + 806);
    const auto *msg_807 = buffer.data(msg + 807);
    const auto *msg_808 = buffer.data(msg + 808);
    const auto *msg_809 = buffer.data(msg + 809);
    const auto *msg_810 = buffer.data(msg + 810);
    const auto *msg_815 = buffer.data(msg + 815);
    const auto *msg_819 = buffer.data(msg + 819);
    const auto *msg_820 = buffer.data(msg + 820);
    const auto *msg_821 = buffer.data(msg + 821);
    const auto *msg_822 = buffer.data(msg + 822);
    const auto *msg_824 = buffer.data(msg + 824);

    const auto *msh1_924 = buffer.data(msh1 + 924);
    const auto *msh1_929 = buffer.data(msh1 + 929);
    const auto *msh1_933 = buffer.data(msh1 + 933);
    const auto *msh1_945 = buffer.data(msh1 + 945);
    const auto *msh1_946 = buffer.data(msh1 + 946);
    const auto *msh1_948 = buffer.data(msh1 + 948);
    const auto *msh1_951 = buffer.data(msh1 + 951);
    const auto *msh1_960 = buffer.data(msh1 + 960);
    const auto *msh1_1067 = buffer.data(msh1 + 1067);
    const auto *msh1_1068 = buffer.data(msh1 + 1068);
    const auto *msh1_1070 = buffer.data(msh1 + 1070);
    const auto *msh1_1071 = buffer.data(msh1 + 1071);
    const auto *msh1_1074 = buffer.data(msh1 + 1074);
    const auto *msh1_1076 = buffer.data(msh1 + 1076);
    const auto *msh1_1077 = buffer.data(msh1 + 1077);
    const auto *msh1_1080 = buffer.data(msh1 + 1080);
    const auto *msh1_1086 = buffer.data(msh1 + 1086);
    const auto *msh1_1088 = buffer.data(msh1 + 1088);
    const auto *msh1_1089 = buffer.data(msh1 + 1089);
    const auto *msh1_1091 = buffer.data(msh1 + 1091);
    const auto *msh1_1092 = buffer.data(msh1 + 1092);
    const auto *msh1_1095 = buffer.data(msh1 + 1095);
    const auto *msh1_1097 = buffer.data(msh1 + 1097);
    const auto *msh1_1098 = buffer.data(msh1 + 1098);
    const auto *msh1_1101 = buffer.data(msh1 + 1101);
    const auto *msh1_1107 = buffer.data(msh1 + 1107);
    const auto *msh1_1109 = buffer.data(msh1 + 1109);
    const auto *msh1_1110 = buffer.data(msh1 + 1110);
    const auto *msh1_1112 = buffer.data(msh1 + 1112);
    const auto *msh1_1116 = buffer.data(msh1 + 1116);
    const auto *msh1_1119 = buffer.data(msh1 + 1119);
    const auto *msh1_1128 = buffer.data(msh1 + 1128);
    const auto *msh1_1130 = buffer.data(msh1 + 1130);
    const auto *msh1_1131 = buffer.data(msh1 + 1131);
    const auto *msh1_1133 = buffer.data(msh1 + 1133);
    const auto *msh1_1134 = buffer.data(msh1 + 1134);
    const auto *msh1_1139 = buffer.data(msh1 + 1139);
    const auto *msh1_1143 = buffer.data(msh1 + 1143);
    const auto *msh1_1149 = buffer.data(msh1 + 1149);
    const auto *msh1_1150 = buffer.data(msh1 + 1150);
    const auto *msh1_1151 = buffer.data(msh1 + 1151);
    const auto *msh1_1152 = buffer.data(msh1 + 1152);
    const auto *msh1_1154 = buffer.data(msh1 + 1154);

    const auto *nsf0_540 = buffer.data(nsf0 + 540);
    const auto *nsf0_541 = buffer.data(nsf0 + 541);
    const auto *nsf0_542 = buffer.data(nsf0 + 542);
    const auto *nsf0_550 = buffer.data(nsf0 + 550);
    const auto *nsf0_551 = buffer.data(nsf0 + 551);
    const auto *nsf0_553 = buffer.data(nsf0 + 553);
    const auto *nsf0_555 = buffer.data(nsf0 + 555);
    const auto *nsf0_556 = buffer.data(nsf0 + 556);
    const auto *nsf0_557 = buffer.data(nsf0 + 557);
    const auto *nsf0_558 = buffer.data(nsf0 + 558);
    const auto *nsf0_559 = buffer.data(nsf0 + 559);
    const auto *nsf0_562 = buffer.data(nsf0 + 562);
    const auto *nsf0_564 = buffer.data(nsf0 + 564);
    const auto *nsf0_565 = buffer.data(nsf0 + 565);
    const auto *nsf0_567 = buffer.data(nsf0 + 567);
    const auto *nsf0_568 = buffer.data(nsf0 + 568);
    const auto *nsf0_569 = buffer.data(nsf0 + 569);

    const auto *nsf1_540 = buffer.data(nsf1 + 540);
    const auto *nsf1_541 = buffer.data(nsf1 + 541);
    const auto *nsf1_542 = buffer.data(nsf1 + 542);
    const auto *nsf1_550 = buffer.data(nsf1 + 550);
    const auto *nsf1_551 = buffer.data(nsf1 + 551);
    const auto *nsf1_553 = buffer.data(nsf1 + 553);
    const auto *nsf1_555 = buffer.data(nsf1 + 555);
    const auto *nsf1_556 = buffer.data(nsf1 + 556);
    const auto *nsf1_557 = buffer.data(nsf1 + 557);
    const auto *nsf1_558 = buffer.data(nsf1 + 558);
    const auto *nsf1_559 = buffer.data(nsf1 + 559);
    const auto *nsf1_562 = buffer.data(nsf1 + 562);
    const auto *nsf1_564 = buffer.data(nsf1 + 564);
    const auto *nsf1_565 = buffer.data(nsf1 + 565);
    const auto *nsf1_567 = buffer.data(nsf1 + 567);
    const auto *nsf1_568 = buffer.data(nsf1 + 568);
    const auto *nsf1_569 = buffer.data(nsf1 + 569);

    const auto *nsg_764 = buffer.data(nsg + 764);
    const auto *nsg_765 = buffer.data(nsg + 765);
    const auto *nsg_767 = buffer.data(nsg + 767);
    const auto *nsg_768 = buffer.data(nsg + 768);
    const auto *nsg_770 = buffer.data(nsg + 770);
    const auto *nsg_775 = buffer.data(nsg + 775);
    const auto *nsg_776 = buffer.data(nsg + 776);
    const auto *nsg_777 = buffer.data(nsg + 777);
    const auto *nsg_778 = buffer.data(nsg + 778);
    const auto *nsg_779 = buffer.data(nsg + 779);
    const auto *nsg_780 = buffer.data(nsg + 780);
    const auto *nsg_782 = buffer.data(nsg + 782);
    const auto *nsg_783 = buffer.data(nsg + 783);
    const auto *nsg_785 = buffer.data(nsg + 785);
    const auto *nsg_790 = buffer.data(nsg + 790);
    const auto *nsg_791 = buffer.data(nsg + 791);
    const auto *nsg_792 = buffer.data(nsg + 792);
    const auto *nsg_793 = buffer.data(nsg + 793);
    const auto *nsg_794 = buffer.data(nsg + 794);
    const auto *nsg_795 = buffer.data(nsg + 795);
    const auto *nsg_797 = buffer.data(nsg + 797);
    const auto *nsg_798 = buffer.data(nsg + 798);
    const auto *nsg_800 = buffer.data(nsg + 800);
    const auto *nsg_805 = buffer.data(nsg + 805);
    const auto *nsg_806 = buffer.data(nsg + 806);
    const auto *nsg_807 = buffer.data(nsg + 807);
    const auto *nsg_808 = buffer.data(nsg + 808);
    const auto *nsg_809 = buffer.data(nsg + 809);
    const auto *nsg_810 = buffer.data(nsg + 810);
    const auto *nsg_811 = buffer.data(nsg + 811);
    const auto *nsg_812 = buffer.data(nsg + 812);
    const auto *nsg_813 = buffer.data(nsg + 813);
    const auto *nsg_814 = buffer.data(nsg + 814);
    const auto *nsg_815 = buffer.data(nsg + 815);
    const auto *nsg_819 = buffer.data(nsg + 819);
    const auto *nsg_820 = buffer.data(nsg + 820);
    const auto *nsg_821 = buffer.data(nsg + 821);
    const auto *nsg_822 = buffer.data(nsg + 822);
    const auto *nsg_824 = buffer.data(nsg + 824);
    const auto *nsg_825 = buffer.data(nsg + 825);
    const auto *nsg_826 = buffer.data(nsg + 826);
    const auto *nsg_828 = buffer.data(nsg + 828);
    const auto *nsg_830 = buffer.data(nsg + 830);
    const auto *nsg_831 = buffer.data(nsg + 831);
    const auto *nsg_833 = buffer.data(nsg + 833);
    const auto *nsg_834 = buffer.data(nsg + 834);
    const auto *nsg_835 = buffer.data(nsg + 835);
    const auto *nsg_836 = buffer.data(nsg + 836);
    const auto *nsg_837 = buffer.data(nsg + 837);
    const auto *nsg_838 = buffer.data(nsg + 838);
    const auto *nsg_839 = buffer.data(nsg + 839);
    const auto *nsg_842 = buffer.data(nsg + 842);
    const auto *nsg_844 = buffer.data(nsg + 844);
    const auto *nsg_845 = buffer.data(nsg + 845);
    const auto *nsg_847 = buffer.data(nsg + 847);
    const auto *nsg_848 = buffer.data(nsg + 848);
    const auto *nsg_849 = buffer.data(nsg + 849);
    const auto *nsg_850 = buffer.data(nsg + 850);
    const auto *nsg_851 = buffer.data(nsg + 851);
    const auto *nsg_852 = buffer.data(nsg + 852);
    const auto *nsg_853 = buffer.data(nsg + 853);
    const auto *nsg_854 = buffer.data(nsg + 854);

#pragma omp simd aligned(t_1067, t_1068, t_1069, t_1070, pa_x, pc_x, pc_y, msh0_1067, \
                         msh0_1068, msh0_1070, msg_629, msh1_1067, msh1_1068, msh1_1070, \
                         nsg_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = pa_x[k] * msh0_1067[k]
                    - f_8 * pc_x[k] * msh1_1067[k];

        t_1068[k] = pa_x[k] * msh0_1068[k]
                    - f_8 * pc_x[k] * msh1_1068[k];

        t_1069[k] = f_18 * msg_629[k]
                    + f_3 * pc_y[k] * nsg_764[k];

        t_1070[k] = pa_x[k] * msh0_1070[k]
                    - f_8 * pc_x[k] * msh1_1070[k];
    }

#pragma omp simd aligned(t_1071, t_1072, t_1073, pa_x, pc_x, pc_y, pc_z, msh0_1071, msg_615, \
                         msg_630, msg_765, msh1_1071, nsg_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1071[k] = pa_x[k] * msh0_1071[k]
                    + f_19 * msg_765[k]
                    - f_8 * pc_x[k] * msh1_1071[k];

        t_1072[k] = f_11 * msg_630[k]
                    + f_3 * pc_y[k] * nsg_765[k];

        t_1073[k] = f_17 * msg_615[k]
                    + f_3 * pc_z[k] * nsg_765[k];
    }

#pragma omp simd aligned(t_1074, t_1075, t_1076, pa_x, pc_x, pc_y, msh0_1074, msh0_1076, \
                         msg_632, msg_768, msg_770, msh1_1074, msh1_1076, \
                         nsg_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = pa_x[k] * msh0_1074[k]
                    + f_11 * msg_768[k]
                    - f_8 * pc_x[k] * msh1_1074[k];

        t_1075[k] = f_11 * msg_632[k]
                    + f_3 * pc_y[k] * nsg_767[k];

        t_1076[k] = pa_x[k] * msh0_1076[k]
                    + f_11 * msg_770[k]
                    - f_8 * pc_x[k] * msh1_1076[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, pa_x, pc_x, pc_y, pc_z, msh0_1077, msg_618, \
                         msg_635, msg_771, msh1_1077, nsg_768, \
                         nsg_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = pa_x[k] * msh0_1077[k]
                    + f_10 * msg_771[k]
                    - f_8 * pc_x[k] * msh1_1077[k];

        t_1078[k] = f_17 * msg_618[k]
                    + f_3 * pc_z[k] * nsg_768[k];

        t_1079[k] = f_11 * msg_635[k]
                    + f_3 * pc_y[k] * nsg_770[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, pa_x, pc_x, msh0_1080, msg_774, \
                         msg_775, msg_776, msg_777, msh1_1080, nsg_775, nsg_776, \
                         nsg_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = pa_x[k] * msh0_1080[k]
                    + f_10 * msg_774[k]
                    - f_8 * pc_x[k] * msh1_1080[k];

        t_1081[k] = f_9 * msg_775[k]
                    + f_3 * pc_x[k] * nsg_775[k];

        t_1082[k] = f_9 * msg_776[k]
                    + f_3 * pc_x[k] * nsg_776[k];

        t_1083[k] = f_9 * msg_777[k]
                    + f_3 * pc_x[k] * nsg_777[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, t_1087, pa_x, pc_x, pc_z, msh0_1086, msg_625, \
                         msg_778, msg_779, msh1_1086, nsg_775, nsg_778, \
                         nsg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_9 * msg_778[k]
                    + f_3 * pc_x[k] * nsg_778[k];

        t_1085[k] = f_9 * msg_779[k]
                    + f_3 * pc_x[k] * nsg_779[k];

        t_1086[k] = pa_x[k] * msh0_1086[k]
                    - f_8 * pc_x[k] * msh1_1086[k];

        t_1087[k] = f_17 * msg_625[k]
                    + f_3 * pc_z[k] * nsg_775[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, t_1091, pa_x, pc_x, pc_y, msh0_1088, \
                         msh0_1089, msh0_1091, msg_644, msh1_1088, msh1_1089, msh1_1091, \
                         nsg_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = pa_x[k] * msh0_1088[k]
                    - f_8 * pc_x[k] * msh1_1088[k];

        t_1089[k] = pa_x[k] * msh0_1089[k]
                    - f_8 * pc_x[k] * msh1_1089[k];

        t_1090[k] = f_11 * msg_644[k]
                    + f_3 * pc_y[k] * nsg_779[k];

        t_1091[k] = pa_x[k] * msh0_1091[k]
                    - f_8 * pc_x[k] * msh1_1091[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pa_x, pc_x, pc_y, pc_z, msh0_1092, msg_630, \
                         msg_645, msg_780, msh1_1092, nsg_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = pa_x[k] * msh0_1092[k]
                    + f_19 * msg_780[k]
                    - f_8 * pc_x[k] * msh1_1092[k];

        t_1093[k] = f_10 * msg_645[k]
                    + f_3 * pc_y[k] * nsg_780[k];

        t_1094[k] = f_16 * msg_630[k]
                    + f_3 * pc_z[k] * nsg_780[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pa_x, pc_x, pc_y, msh0_1095, msh0_1097, \
                         msg_647, msg_783, msg_785, msh1_1095, msh1_1097, \
                         nsg_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = pa_x[k] * msh0_1095[k]
                    + f_11 * msg_783[k]
                    - f_8 * pc_x[k] * msh1_1095[k];

        t_1096[k] = f_10 * msg_647[k]
                    + f_3 * pc_y[k] * nsg_782[k];

        t_1097[k] = pa_x[k] * msh0_1097[k]
                    + f_11 * msg_785[k]
                    - f_8 * pc_x[k] * msh1_1097[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, pa_x, pc_x, pc_y, pc_z, msh0_1098, msg_633, \
                         msg_650, msg_786, msh1_1098, nsg_783, \
                         nsg_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = pa_x[k] * msh0_1098[k]
                    + f_10 * msg_786[k]
                    - f_8 * pc_x[k] * msh1_1098[k];

        t_1099[k] = f_16 * msg_633[k]
                    + f_3 * pc_z[k] * nsg_783[k];

        t_1100[k] = f_10 * msg_650[k]
                    + f_3 * pc_y[k] * nsg_785[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, t_1104, pa_x, pc_x, msh0_1101, msg_789, \
                         msg_790, msg_791, msg_792, msh1_1101, nsg_790, nsg_791, \
                         nsg_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = pa_x[k] * msh0_1101[k]
                    + f_10 * msg_789[k]
                    - f_8 * pc_x[k] * msh1_1101[k];

        t_1102[k] = f_9 * msg_790[k]
                    + f_3 * pc_x[k] * nsg_790[k];

        t_1103[k] = f_9 * msg_791[k]
                    + f_3 * pc_x[k] * nsg_791[k];

        t_1104[k] = f_9 * msg_792[k]
                    + f_3 * pc_x[k] * nsg_792[k];
    }

#pragma omp simd aligned(t_1105, t_1106, t_1107, t_1108, pa_x, pc_x, pc_z, msh0_1107, msg_640, \
                         msg_793, msg_794, msh1_1107, nsg_790, nsg_793, \
                         nsg_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1105[k] = f_9 * msg_793[k]
                    + f_3 * pc_x[k] * nsg_793[k];

        t_1106[k] = f_9 * msg_794[k]
                    + f_3 * pc_x[k] * nsg_794[k];

        t_1107[k] = pa_x[k] * msh0_1107[k]
                    - f_8 * pc_x[k] * msh1_1107[k];

        t_1108[k] = f_16 * msg_640[k]
                    + f_3 * pc_z[k] * nsg_790[k];
    }

#pragma omp simd aligned(t_1109, t_1110, t_1111, t_1112, pa_x, pc_x, pc_y, msh0_1109, \
                         msh0_1110, msh0_1112, msg_659, msh1_1109, msh1_1110, msh1_1112, \
                         nsg_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1109[k] = pa_x[k] * msh0_1109[k]
                    - f_8 * pc_x[k] * msh1_1109[k];

        t_1110[k] = pa_x[k] * msh0_1110[k]
                    - f_8 * pc_x[k] * msh1_1110[k];

        t_1111[k] = f_10 * msg_659[k]
                    + f_3 * pc_y[k] * nsg_794[k];

        t_1112[k] = pa_x[k] * msh0_1112[k]
                    - f_8 * pc_x[k] * msh1_1112[k];
    }

#pragma omp simd aligned(t_1113, t_1114, t_1115, pa_y, pc_y, pc_z, msh0_924, msg_645, msg_660, \
                         msh1_924, nsg_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1113[k] = pa_y[k] * msh0_924[k]
                    - f_8 * pc_y[k] * msh1_924[k];

        t_1114[k] = f_9 * msg_660[k]
                    + f_3 * pc_y[k] * nsg_795[k];

        t_1115[k] = f_15 * msg_645[k]
                    + f_3 * pc_z[k] * nsg_795[k];
    }

#pragma omp simd aligned(t_1116, t_1117, t_1118, pa_x, pa_y, pc_x, pc_y, msh0_929, msh0_1116, \
                         msg_662, msg_798, msh1_929, msh1_1116, \
                         nsg_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1116[k] = pa_x[k] * msh0_1116[k]
                    + f_11 * msg_798[k]
                    - f_8 * pc_x[k] * msh1_1116[k];

        t_1117[k] = f_9 * msg_662[k]
                    + f_3 * pc_y[k] * nsg_797[k];

        t_1118[k] = pa_y[k] * msh0_929[k]
                    - f_8 * pc_y[k] * msh1_929[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, pa_x, pc_x, pc_y, pc_z, msh0_1119, msg_648, \
                         msg_665, msg_801, msh1_1119, nsg_798, \
                         nsg_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = pa_x[k] * msh0_1119[k]
                    + f_10 * msg_801[k]
                    - f_8 * pc_x[k] * msh1_1119[k];

        t_1120[k] = f_15 * msg_648[k]
                    + f_3 * pc_z[k] * nsg_798[k];

        t_1121[k] = f_9 * msg_665[k]
                    + f_3 * pc_y[k] * nsg_800[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, t_1125, pa_y, pc_x, pc_y, msh0_933, msg_805, \
                         msg_806, msg_807, msh1_933, nsg_805, nsg_806, \
                         nsg_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = pa_y[k] * msh0_933[k]
                    - f_8 * pc_y[k] * msh1_933[k];

        t_1123[k] = f_9 * msg_805[k]
                    + f_3 * pc_x[k] * nsg_805[k];

        t_1124[k] = f_9 * msg_806[k]
                    + f_3 * pc_x[k] * nsg_806[k];

        t_1125[k] = f_9 * msg_807[k]
                    + f_3 * pc_x[k] * nsg_807[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, t_1129, pa_x, pc_x, pc_z, msh0_1128, msg_655, \
                         msg_808, msg_809, msh1_1128, nsg_805, nsg_808, \
                         nsg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_9 * msg_808[k]
                    + f_3 * pc_x[k] * nsg_808[k];

        t_1127[k] = f_9 * msg_809[k]
                    + f_3 * pc_x[k] * nsg_809[k];

        t_1128[k] = pa_x[k] * msh0_1128[k]
                    - f_8 * pc_x[k] * msh1_1128[k];

        t_1129[k] = f_15 * msg_655[k]
                    + f_3 * pc_z[k] * nsg_805[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, t_1133, pa_x, pc_x, pc_y, msh0_1130, \
                         msh0_1131, msh0_1133, msg_674, msh1_1130, msh1_1131, msh1_1133, \
                         nsg_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = pa_x[k] * msh0_1130[k]
                    - f_8 * pc_x[k] * msh1_1130[k];

        t_1131[k] = pa_x[k] * msh0_1131[k]
                    - f_8 * pc_x[k] * msh1_1131[k];

        t_1132[k] = f_9 * msg_674[k]
                    + f_3 * pc_y[k] * nsg_809[k];

        t_1133[k] = pa_x[k] * msh0_1133[k]
                    - f_8 * pc_x[k] * msh1_1133[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, pa_x, pc_x, pc_y, pc_z, msh0_1134, \
                         msg_660, msg_810, msh1_1134, nsf0_540, nsf1_540, nsg_810, \
                         nsg_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = pa_x[k] * msh0_1134[k]
                    + f_19 * msg_810[k]
                    - f_8 * pc_x[k] * msh1_1134[k];

        t_1135[k] = f_3 * pc_y[k] * nsg_810[k];

        t_1136[k] = f_12 * msg_660[k]
                    + f_3 * pc_z[k] * nsg_810[k];

        t_1137[k] = f_4 * nsf0_540[k]
                    - f_5 * nsf1_540[k]
                    + f_3 * pc_y[k] * nsg_811[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, pa_x, pc_x, pc_y, msh0_1139, msg_815, \
                         msh1_1139, nsf0_541, nsf1_541, nsg_812, \
                         nsg_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = f_3 * pc_y[k] * nsg_812[k];

        t_1139[k] = pa_x[k] * msh0_1139[k]
                    + f_11 * msg_815[k]
                    - f_8 * pc_x[k] * msh1_1139[k];

        t_1140[k] = f_6 * nsf0_541[k]
                    - f_7 * nsf1_541[k]
                    + f_3 * pc_y[k] * nsg_813[k];
    }

#pragma omp simd aligned(t_1141, t_1142, t_1143, t_1144, pa_x, pc_x, pc_y, msh0_1143, msg_819, \
                         msg_820, msh1_1143, nsf0_542, nsf1_542, nsg_814, nsg_815, \
                         nsg_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1141[k] = f_4 * nsf0_542[k]
                    - f_5 * nsf1_542[k]
                    + f_3 * pc_y[k] * nsg_814[k];

        t_1142[k] = f_3 * pc_y[k] * nsg_815[k];

        t_1143[k] = pa_x[k] * msh0_1143[k]
                    + f_10 * msg_819[k]
                    - f_8 * pc_x[k] * msh1_1143[k];

        t_1144[k] = f_9 * msg_820[k]
                    + f_3 * pc_x[k] * nsg_820[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, t_1148, pc_x, pc_y, msg_821, msg_822, \
                         msg_824, nsg_819, nsg_821, nsg_822, nsg_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_9 * msg_821[k]
                    + f_3 * pc_x[k] * nsg_821[k];

        t_1146[k] = f_9 * msg_822[k]
                    + f_3 * pc_x[k] * nsg_822[k];

        t_1147[k] = f_3 * pc_y[k] * nsg_819[k];

        t_1148[k] = f_9 * msg_824[k]
                    + f_3 * pc_x[k] * nsg_824[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pa_x, pc_x, msh0_1149, msh0_1150, \
                         msh0_1151, msh0_1152, msh1_1149, msh1_1150, msh1_1151, \
                         msh1_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = pa_x[k] * msh0_1149[k]
                    - f_8 * pc_x[k] * msh1_1149[k];

        t_1150[k] = pa_x[k] * msh0_1150[k]
                    - f_8 * pc_x[k] * msh1_1150[k];

        t_1151[k] = pa_x[k] * msh0_1151[k]
                    - f_8 * pc_x[k] * msh1_1151[k];

        t_1152[k] = pa_x[k] * msh0_1152[k]
                    - f_8 * pc_x[k] * msh1_1152[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, pa_x, pc_x, pc_y, msh0_1154, \
                         msh1_1154, nsf0_550, nsf0_551, nsf1_550, nsf1_551, nsg_824, nsg_825, \
                         nsg_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_3 * pc_y[k] * nsg_824[k];

        t_1154[k] = pa_x[k] * msh0_1154[k]
                    - f_8 * pc_x[k] * msh1_1154[k];

        t_1155[k] = f_1 * nsf0_550[k]
                    - f_2 * nsf1_550[k]
                    + f_3 * pc_x[k] * nsg_825[k];

        t_1156[k] = f_13 * nsf0_551[k]
                    - f_14 * nsf1_551[k]
                    + f_3 * pc_x[k] * nsg_826[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, t_1160, pc_x, pc_z, nsf0_553, nsf0_555, \
                         nsf1_553, nsf1_555, nsg_825, nsg_826, nsg_828, \
                         nsg_830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_3 * pc_z[k] * nsg_825[k];

        t_1158[k] = f_6 * nsf0_553[k]
                    - f_7 * nsf1_553[k]
                    + f_3 * pc_x[k] * nsg_828[k];

        t_1159[k] = f_3 * pc_z[k] * nsg_826[k];

        t_1160[k] = f_6 * nsf0_555[k]
                    - f_7 * nsf1_555[k]
                    + f_3 * pc_x[k] * nsg_830[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, t_1164, pc_x, pc_z, nsf0_556, nsf0_558, \
                         nsf0_559, nsf1_556, nsf1_558, nsf1_559, nsg_828, nsg_831, nsg_833, \
                         nsg_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_4 * nsf0_556[k]
                    - f_5 * nsf1_556[k]
                    + f_3 * pc_x[k] * nsg_831[k];

        t_1162[k] = f_3 * pc_z[k] * nsg_828[k];

        t_1163[k] = f_4 * nsf0_558[k]
                    - f_5 * nsf1_558[k]
                    + f_3 * pc_x[k] * nsg_833[k];

        t_1164[k] = f_4 * nsf0_559[k]
                    - f_5 * nsf1_559[k]
                    + f_3 * pc_x[k] * nsg_834[k];
    }

#pragma omp simd aligned(t_1165, t_1166, t_1167, t_1168, t_1169, t_1170, pc_x, pc_y, msg_685, \
                         nsf0_556, nsf1_556, nsg_835, nsg_836, nsg_837, nsg_838, \
                         nsg_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1165[k] = f_3 * pc_x[k] * nsg_835[k];

        t_1166[k] = f_3 * pc_x[k] * nsg_836[k];

        t_1167[k] = f_3 * pc_x[k] * nsg_837[k];

        t_1168[k] = f_3 * pc_x[k] * nsg_838[k];

        t_1169[k] = f_3 * pc_x[k] * nsg_839[k];

        t_1170[k] = f_0 * msg_685[k]
                    + f_1 * nsf0_556[k]
                    - f_2 * nsf1_556[k]
                    + f_3 * pc_y[k] * nsg_835[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, t_1174, pc_y, pc_z, msg_689, nsf0_556, \
                         nsf0_557, nsf1_556, nsf1_557, nsg_835, nsg_836, nsg_837, \
                         nsg_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_3 * pc_z[k] * nsg_835[k];

        t_1172[k] = f_4 * nsf0_556[k]
                    - f_5 * nsf1_556[k]
                    + f_3 * pc_z[k] * nsg_836[k];

        t_1173[k] = f_6 * nsf0_557[k]
                    - f_7 * nsf1_557[k]
                    + f_3 * pc_z[k] * nsg_837[k];

        t_1174[k] = f_0 * msg_689[k]
                    + f_3 * pc_y[k] * nsg_839[k];
    }

#pragma omp simd aligned(t_1175, t_1176, t_1177, pa_z, pc_z, msh0_945, msh0_946, msh1_945, \
                         msh1_946, nsf0_559, nsf1_559, nsg_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1175[k] = f_1 * nsf0_559[k]
                    - f_2 * nsf1_559[k]
                    + f_3 * pc_z[k] * nsg_839[k];

        t_1176[k] = pa_z[k] * msh0_945[k]
                    - f_8 * pc_z[k] * msh1_945[k];

        t_1177[k] = pa_z[k] * msh0_946[k]
                    - f_8 * pc_z[k] * msh1_946[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, pa_z, pc_x, pc_z, msh0_948, msh1_948, \
                         nsf0_562, nsf0_564, nsf1_562, nsf1_564, nsg_842, \
                         nsg_844 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = f_13 * nsf0_562[k]
                    - f_14 * nsf1_562[k]
                    + f_3 * pc_x[k] * nsg_842[k];

        t_1179[k] = pa_z[k] * msh0_948[k]
                    - f_8 * pc_z[k] * msh1_948[k];

        t_1180[k] = f_6 * nsf0_564[k]
                    - f_7 * nsf1_564[k]
                    + f_3 * pc_x[k] * nsg_844[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, pa_z, pc_x, pc_z, msh0_951, msh1_951, \
                         nsf0_565, nsf0_567, nsf1_565, nsf1_567, nsg_845, \
                         nsg_847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = f_6 * nsf0_565[k]
                    - f_7 * nsf1_565[k]
                    + f_3 * pc_x[k] * nsg_845[k];

        t_1182[k] = pa_z[k] * msh0_951[k]
                    - f_8 * pc_z[k] * msh1_951[k];

        t_1183[k] = f_4 * nsf0_567[k]
                    - f_5 * nsf1_567[k]
                    + f_3 * pc_x[k] * nsg_847[k];
    }

#pragma omp simd aligned(t_1184, t_1185, t_1186, t_1187, t_1188, pc_x, nsf0_568, nsf0_569, \
                         nsf1_568, nsf1_569, nsg_848, nsg_849, nsg_850, nsg_851, \
                         nsg_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1184[k] = f_4 * nsf0_568[k]
                    - f_5 * nsf1_568[k]
                    + f_3 * pc_x[k] * nsg_848[k];

        t_1185[k] = f_4 * nsf0_569[k]
                    - f_5 * nsf1_569[k]
                    + f_3 * pc_x[k] * nsg_849[k];

        t_1186[k] = f_3 * pc_x[k] * nsg_850[k];

        t_1187[k] = f_3 * pc_x[k] * nsg_851[k];

        t_1188[k] = f_3 * pc_x[k] * nsg_852[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pa_z, pc_x, pc_z, msh0_960, msg_685, \
                         msh1_960, nsg_850, nsg_853, nsg_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_3 * pc_x[k] * nsg_853[k];

        t_1190[k] = f_3 * pc_x[k] * nsg_854[k];

        t_1191[k] = pa_z[k] * msh0_960[k]
                    - f_8 * pc_z[k] * msh1_960[k];

        t_1192[k] = f_9 * msg_685[k]
                    + f_3 * pc_z[k] * nsg_850[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msh0,
                                                           const size_t msg, const size_t msh1,
                                                           const size_t nsf0, const size_t nsf1,
                                                           const size_t nsg, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);
    auto *t_1260 = buffer.data(target + 1260);
    auto *t_1261 = buffer.data(target + 1261);
    auto *t_1262 = buffer.data(target + 1262);
    auto *t_1263 = buffer.data(target + 1263);
    auto *t_1264 = buffer.data(target + 1264);
    auto *t_1265 = buffer.data(target + 1265);
    auto *t_1266 = buffer.data(target + 1266);
    auto *t_1267 = buffer.data(target + 1267);
    auto *t_1268 = buffer.data(target + 1268);
    auto *t_1269 = buffer.data(target + 1269);
    auto *t_1270 = buffer.data(target + 1270);
    auto *t_1271 = buffer.data(target + 1271);
    auto *t_1272 = buffer.data(target + 1272);
    auto *t_1273 = buffer.data(target + 1273);
    auto *t_1274 = buffer.data(target + 1274);
    auto *t_1275 = buffer.data(target + 1275);
    auto *t_1276 = buffer.data(target + 1276);
    auto *t_1277 = buffer.data(target + 1277);
    auto *t_1278 = buffer.data(target + 1278);
    auto *t_1279 = buffer.data(target + 1279);
    auto *t_1280 = buffer.data(target + 1280);
    auto *t_1281 = buffer.data(target + 1281);
    auto *t_1282 = buffer.data(target + 1282);
    auto *t_1283 = buffer.data(target + 1283);
    auto *t_1284 = buffer.data(target + 1284);
    auto *t_1285 = buffer.data(target + 1285);
    auto *t_1286 = buffer.data(target + 1286);
    auto *t_1287 = buffer.data(target + 1287);
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);
    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);
    auto *t_1296 = buffer.data(target + 1296);
    auto *t_1297 = buffer.data(target + 1297);
    auto *t_1298 = buffer.data(target + 1298);
    auto *t_1299 = buffer.data(target + 1299);
    auto *t_1300 = buffer.data(target + 1300);
    auto *t_1301 = buffer.data(target + 1301);
    auto *t_1302 = buffer.data(target + 1302);
    auto *t_1303 = buffer.data(target + 1303);
    auto *t_1304 = buffer.data(target + 1304);
    auto *t_1305 = buffer.data(target + 1305);
    auto *t_1306 = buffer.data(target + 1306);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_962 = buffer.data(msh0 + 962);
    const auto *msh0_963 = buffer.data(msh0 + 963);

    const auto *msg_686 = buffer.data(msg + 686);
    const auto *msg_687 = buffer.data(msg + 687);
    const auto *msg_689 = buffer.data(msg + 689);
    const auto *msg_700 = buffer.data(msg + 700);
    const auto *msg_704 = buffer.data(msg + 704);
    const auto *msg_715 = buffer.data(msg + 715);
    const auto *msg_717 = buffer.data(msg + 717);
    const auto *msg_718 = buffer.data(msg + 718);
    const auto *msg_719 = buffer.data(msg + 719);
    const auto *msg_730 = buffer.data(msg + 730);
    const auto *msg_732 = buffer.data(msg + 732);
    const auto *msg_733 = buffer.data(msg + 733);
    const auto *msg_734 = buffer.data(msg + 734);
    const auto *msg_745 = buffer.data(msg + 745);
    const auto *msg_747 = buffer.data(msg + 747);
    const auto *msg_748 = buffer.data(msg + 748);
    const auto *msg_749 = buffer.data(msg + 749);
    const auto *msg_760 = buffer.data(msg + 760);
    const auto *msg_762 = buffer.data(msg + 762);
    const auto *msg_763 = buffer.data(msg + 763);
    const auto *msg_764 = buffer.data(msg + 764);
    const auto *msg_775 = buffer.data(msg + 775);
    const auto *msg_777 = buffer.data(msg + 777);
    const auto *msg_778 = buffer.data(msg + 778);
    const auto *msg_779 = buffer.data(msg + 779);

    const auto *msh1_962 = buffer.data(msh1 + 962);
    const auto *msh1_963 = buffer.data(msh1 + 963);

    const auto *nsf0_569 = buffer.data(nsf0 + 569);
    const auto *nsf0_570 = buffer.data(nsf0 + 570);
    const auto *nsf0_571 = buffer.data(nsf0 + 571);
    const auto *nsf0_572 = buffer.data(nsf0 + 572);
    const auto *nsf0_573 = buffer.data(nsf0 + 573);
    const auto *nsf0_574 = buffer.data(nsf0 + 574);
    const auto *nsf0_575 = buffer.data(nsf0 + 575);
    const auto *nsf0_576 = buffer.data(nsf0 + 576);
    const auto *nsf0_577 = buffer.data(nsf0 + 577);
    const auto *nsf0_578 = buffer.data(nsf0 + 578);
    const auto *nsf0_579 = buffer.data(nsf0 + 579);
    const auto *nsf0_580 = buffer.data(nsf0 + 580);
    const auto *nsf0_581 = buffer.data(nsf0 + 581);
    const auto *nsf0_582 = buffer.data(nsf0 + 582);
    const auto *nsf0_583 = buffer.data(nsf0 + 583);
    const auto *nsf0_584 = buffer.data(nsf0 + 584);
    const auto *nsf0_585 = buffer.data(nsf0 + 585);
    const auto *nsf0_586 = buffer.data(nsf0 + 586);
    const auto *nsf0_587 = buffer.data(nsf0 + 587);
    const auto *nsf0_588 = buffer.data(nsf0 + 588);
    const auto *nsf0_589 = buffer.data(nsf0 + 589);
    const auto *nsf0_590 = buffer.data(nsf0 + 590);
    const auto *nsf0_591 = buffer.data(nsf0 + 591);
    const auto *nsf0_592 = buffer.data(nsf0 + 592);
    const auto *nsf0_593 = buffer.data(nsf0 + 593);
    const auto *nsf0_594 = buffer.data(nsf0 + 594);
    const auto *nsf0_595 = buffer.data(nsf0 + 595);
    const auto *nsf0_596 = buffer.data(nsf0 + 596);
    const auto *nsf0_597 = buffer.data(nsf0 + 597);
    const auto *nsf0_598 = buffer.data(nsf0 + 598);
    const auto *nsf0_599 = buffer.data(nsf0 + 599);
    const auto *nsf0_600 = buffer.data(nsf0 + 600);
    const auto *nsf0_601 = buffer.data(nsf0 + 601);
    const auto *nsf0_602 = buffer.data(nsf0 + 602);
    const auto *nsf0_603 = buffer.data(nsf0 + 603);
    const auto *nsf0_604 = buffer.data(nsf0 + 604);
    const auto *nsf0_605 = buffer.data(nsf0 + 605);
    const auto *nsf0_606 = buffer.data(nsf0 + 606);
    const auto *nsf0_607 = buffer.data(nsf0 + 607);
    const auto *nsf0_608 = buffer.data(nsf0 + 608);
    const auto *nsf0_609 = buffer.data(nsf0 + 609);
    const auto *nsf0_610 = buffer.data(nsf0 + 610);
    const auto *nsf0_611 = buffer.data(nsf0 + 611);
    const auto *nsf0_612 = buffer.data(nsf0 + 612);
    const auto *nsf0_613 = buffer.data(nsf0 + 613);
    const auto *nsf0_614 = buffer.data(nsf0 + 614);
    const auto *nsf0_615 = buffer.data(nsf0 + 615);
    const auto *nsf0_616 = buffer.data(nsf0 + 616);
    const auto *nsf0_617 = buffer.data(nsf0 + 617);
    const auto *nsf0_618 = buffer.data(nsf0 + 618);
    const auto *nsf0_619 = buffer.data(nsf0 + 619);
    const auto *nsf0_620 = buffer.data(nsf0 + 620);
    const auto *nsf0_621 = buffer.data(nsf0 + 621);
    const auto *nsf0_622 = buffer.data(nsf0 + 622);
    const auto *nsf0_623 = buffer.data(nsf0 + 623);
    const auto *nsf0_624 = buffer.data(nsf0 + 624);

    const auto *nsf1_569 = buffer.data(nsf1 + 569);
    const auto *nsf1_570 = buffer.data(nsf1 + 570);
    const auto *nsf1_571 = buffer.data(nsf1 + 571);
    const auto *nsf1_572 = buffer.data(nsf1 + 572);
    const auto *nsf1_573 = buffer.data(nsf1 + 573);
    const auto *nsf1_574 = buffer.data(nsf1 + 574);
    const auto *nsf1_575 = buffer.data(nsf1 + 575);
    const auto *nsf1_576 = buffer.data(nsf1 + 576);
    const auto *nsf1_577 = buffer.data(nsf1 + 577);
    const auto *nsf1_578 = buffer.data(nsf1 + 578);
    const auto *nsf1_579 = buffer.data(nsf1 + 579);
    const auto *nsf1_580 = buffer.data(nsf1 + 580);
    const auto *nsf1_581 = buffer.data(nsf1 + 581);
    const auto *nsf1_582 = buffer.data(nsf1 + 582);
    const auto *nsf1_583 = buffer.data(nsf1 + 583);
    const auto *nsf1_584 = buffer.data(nsf1 + 584);
    const auto *nsf1_585 = buffer.data(nsf1 + 585);
    const auto *nsf1_586 = buffer.data(nsf1 + 586);
    const auto *nsf1_587 = buffer.data(nsf1 + 587);
    const auto *nsf1_588 = buffer.data(nsf1 + 588);
    const auto *nsf1_589 = buffer.data(nsf1 + 589);
    const auto *nsf1_590 = buffer.data(nsf1 + 590);
    const auto *nsf1_591 = buffer.data(nsf1 + 591);
    const auto *nsf1_592 = buffer.data(nsf1 + 592);
    const auto *nsf1_593 = buffer.data(nsf1 + 593);
    const auto *nsf1_594 = buffer.data(nsf1 + 594);
    const auto *nsf1_595 = buffer.data(nsf1 + 595);
    const auto *nsf1_596 = buffer.data(nsf1 + 596);
    const auto *nsf1_597 = buffer.data(nsf1 + 597);
    const auto *nsf1_598 = buffer.data(nsf1 + 598);
    const auto *nsf1_599 = buffer.data(nsf1 + 599);
    const auto *nsf1_600 = buffer.data(nsf1 + 600);
    const auto *nsf1_601 = buffer.data(nsf1 + 601);
    const auto *nsf1_602 = buffer.data(nsf1 + 602);
    const auto *nsf1_603 = buffer.data(nsf1 + 603);
    const auto *nsf1_604 = buffer.data(nsf1 + 604);
    const auto *nsf1_605 = buffer.data(nsf1 + 605);
    const auto *nsf1_606 = buffer.data(nsf1 + 606);
    const auto *nsf1_607 = buffer.data(nsf1 + 607);
    const auto *nsf1_608 = buffer.data(nsf1 + 608);
    const auto *nsf1_609 = buffer.data(nsf1 + 609);
    const auto *nsf1_610 = buffer.data(nsf1 + 610);
    const auto *nsf1_611 = buffer.data(nsf1 + 611);
    const auto *nsf1_612 = buffer.data(nsf1 + 612);
    const auto *nsf1_613 = buffer.data(nsf1 + 613);
    const auto *nsf1_614 = buffer.data(nsf1 + 614);
    const auto *nsf1_615 = buffer.data(nsf1 + 615);
    const auto *nsf1_616 = buffer.data(nsf1 + 616);
    const auto *nsf1_617 = buffer.data(nsf1 + 617);
    const auto *nsf1_618 = buffer.data(nsf1 + 618);
    const auto *nsf1_619 = buffer.data(nsf1 + 619);
    const auto *nsf1_620 = buffer.data(nsf1 + 620);
    const auto *nsf1_621 = buffer.data(nsf1 + 621);
    const auto *nsf1_622 = buffer.data(nsf1 + 622);
    const auto *nsf1_623 = buffer.data(nsf1 + 623);
    const auto *nsf1_624 = buffer.data(nsf1 + 624);

    const auto *nsg_854 = buffer.data(nsg + 854);
    const auto *nsg_855 = buffer.data(nsg + 855);
    const auto *nsg_856 = buffer.data(nsg + 856);
    const auto *nsg_857 = buffer.data(nsg + 857);
    const auto *nsg_858 = buffer.data(nsg + 858);
    const auto *nsg_859 = buffer.data(nsg + 859);
    const auto *nsg_860 = buffer.data(nsg + 860);
    const auto *nsg_861 = buffer.data(nsg + 861);
    const auto *nsg_862 = buffer.data(nsg + 862);
    const auto *nsg_863 = buffer.data(nsg + 863);
    const auto *nsg_864 = buffer.data(nsg + 864);
    const auto *nsg_865 = buffer.data(nsg + 865);
    const auto *nsg_866 = buffer.data(nsg + 866);
    const auto *nsg_867 = buffer.data(nsg + 867);
    const auto *nsg_868 = buffer.data(nsg + 868);
    const auto *nsg_869 = buffer.data(nsg + 869);
    const auto *nsg_870 = buffer.data(nsg + 870);
    const auto *nsg_871 = buffer.data(nsg + 871);
    const auto *nsg_872 = buffer.data(nsg + 872);
    const auto *nsg_873 = buffer.data(nsg + 873);
    const auto *nsg_874 = buffer.data(nsg + 874);
    const auto *nsg_875 = buffer.data(nsg + 875);
    const auto *nsg_876 = buffer.data(nsg + 876);
    const auto *nsg_877 = buffer.data(nsg + 877);
    const auto *nsg_878 = buffer.data(nsg + 878);
    const auto *nsg_879 = buffer.data(nsg + 879);
    const auto *nsg_880 = buffer.data(nsg + 880);
    const auto *nsg_881 = buffer.data(nsg + 881);
    const auto *nsg_882 = buffer.data(nsg + 882);
    const auto *nsg_883 = buffer.data(nsg + 883);
    const auto *nsg_884 = buffer.data(nsg + 884);
    const auto *nsg_885 = buffer.data(nsg + 885);
    const auto *nsg_886 = buffer.data(nsg + 886);
    const auto *nsg_887 = buffer.data(nsg + 887);
    const auto *nsg_888 = buffer.data(nsg + 888);
    const auto *nsg_889 = buffer.data(nsg + 889);
    const auto *nsg_890 = buffer.data(nsg + 890);
    const auto *nsg_891 = buffer.data(nsg + 891);
    const auto *nsg_892 = buffer.data(nsg + 892);
    const auto *nsg_893 = buffer.data(nsg + 893);
    const auto *nsg_894 = buffer.data(nsg + 894);
    const auto *nsg_895 = buffer.data(nsg + 895);
    const auto *nsg_896 = buffer.data(nsg + 896);
    const auto *nsg_897 = buffer.data(nsg + 897);
    const auto *nsg_898 = buffer.data(nsg + 898);
    const auto *nsg_899 = buffer.data(nsg + 899);
    const auto *nsg_900 = buffer.data(nsg + 900);
    const auto *nsg_901 = buffer.data(nsg + 901);
    const auto *nsg_902 = buffer.data(nsg + 902);
    const auto *nsg_903 = buffer.data(nsg + 903);
    const auto *nsg_904 = buffer.data(nsg + 904);
    const auto *nsg_905 = buffer.data(nsg + 905);
    const auto *nsg_906 = buffer.data(nsg + 906);
    const auto *nsg_907 = buffer.data(nsg + 907);
    const auto *nsg_908 = buffer.data(nsg + 908);
    const auto *nsg_909 = buffer.data(nsg + 909);
    const auto *nsg_910 = buffer.data(nsg + 910);
    const auto *nsg_911 = buffer.data(nsg + 911);
    const auto *nsg_912 = buffer.data(nsg + 912);
    const auto *nsg_913 = buffer.data(nsg + 913);
    const auto *nsg_914 = buffer.data(nsg + 914);
    const auto *nsg_915 = buffer.data(nsg + 915);
    const auto *nsg_916 = buffer.data(nsg + 916);
    const auto *nsg_917 = buffer.data(nsg + 917);
    const auto *nsg_918 = buffer.data(nsg + 918);
    const auto *nsg_919 = buffer.data(nsg + 919);
    const auto *nsg_920 = buffer.data(nsg + 920);
    const auto *nsg_921 = buffer.data(nsg + 921);
    const auto *nsg_922 = buffer.data(nsg + 922);
    const auto *nsg_923 = buffer.data(nsg + 923);
    const auto *nsg_924 = buffer.data(nsg + 924);
    const auto *nsg_925 = buffer.data(nsg + 925);
    const auto *nsg_926 = buffer.data(nsg + 926);
    const auto *nsg_927 = buffer.data(nsg + 927);
    const auto *nsg_928 = buffer.data(nsg + 928);
    const auto *nsg_929 = buffer.data(nsg + 929);
    const auto *nsg_930 = buffer.data(nsg + 930);
    const auto *nsg_931 = buffer.data(nsg + 931);
    const auto *nsg_932 = buffer.data(nsg + 932);
    const auto *nsg_933 = buffer.data(nsg + 933);
    const auto *nsg_934 = buffer.data(nsg + 934);

#pragma omp simd aligned(t_1193, t_1194, t_1195, pa_z, pc_y, pc_z, msh0_962, msh0_963, \
                         msg_686, msg_687, msg_704, msh1_962, msh1_963, \
                         nsg_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = pa_z[k] * msh0_962[k]
                    + f_10 * msg_686[k]
                    - f_8 * pc_z[k] * msh1_962[k];

        t_1194[k] = pa_z[k] * msh0_963[k]
                    + f_11 * msg_687[k]
                    - f_8 * pc_z[k] * msh1_963[k];

        t_1195[k] = f_12 * msg_704[k]
                    + f_3 * pc_y[k] * nsg_854[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_z, msg_689, nsf0_569, nsf0_570, \
                         nsf0_571, nsf1_569, nsf1_570, nsf1_571, nsg_854, nsg_855, \
                         nsg_856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_9 * msg_689[k]
                    + f_1 * nsf0_569[k]
                    - f_2 * nsf1_569[k]
                    + f_3 * pc_z[k] * nsg_854[k];

        t_1197[k] = f_1 * nsf0_570[k]
                    - f_2 * nsf1_570[k]
                    + f_3 * pc_x[k] * nsg_855[k];

        t_1198[k] = f_13 * nsf0_571[k]
                    - f_14 * nsf1_571[k]
                    + f_3 * pc_x[k] * nsg_856[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_x, nsf0_572, nsf0_573, nsf0_574, nsf1_572, \
                         nsf1_573, nsf1_574, nsg_857, nsg_858, \
                         nsg_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_13 * nsf0_572[k]
                    - f_14 * nsf1_572[k]
                    + f_3 * pc_x[k] * nsg_857[k];

        t_1200[k] = f_6 * nsf0_573[k]
                    - f_7 * nsf1_573[k]
                    + f_3 * pc_x[k] * nsg_858[k];

        t_1201[k] = f_6 * nsf0_574[k]
                    - f_7 * nsf1_574[k]
                    + f_3 * pc_x[k] * nsg_859[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pc_x, nsf0_575, nsf0_576, nsf0_577, nsf1_575, \
                         nsf1_576, nsf1_577, nsg_860, nsg_861, \
                         nsg_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_6 * nsf0_575[k]
                    - f_7 * nsf1_575[k]
                    + f_3 * pc_x[k] * nsg_860[k];

        t_1203[k] = f_4 * nsf0_576[k]
                    - f_5 * nsf1_576[k]
                    + f_3 * pc_x[k] * nsg_861[k];

        t_1204[k] = f_4 * nsf0_577[k]
                    - f_5 * nsf1_577[k]
                    + f_3 * pc_x[k] * nsg_862[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, t_1208, t_1209, pc_x, nsf0_578, nsf0_579, \
                         nsf1_578, nsf1_579, nsg_863, nsg_864, nsg_865, nsg_866, \
                         nsg_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_4 * nsf0_578[k]
                    - f_5 * nsf1_578[k]
                    + f_3 * pc_x[k] * nsg_863[k];

        t_1206[k] = f_4 * nsf0_579[k]
                    - f_5 * nsf1_579[k]
                    + f_3 * pc_x[k] * nsg_864[k];

        t_1207[k] = f_3 * pc_x[k] * nsg_865[k];

        t_1208[k] = f_3 * pc_x[k] * nsg_866[k];

        t_1209[k] = f_3 * pc_x[k] * nsg_867[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pc_x, pc_y, pc_z, msg_700, msg_715, \
                         nsf0_576, nsf1_576, nsg_865, nsg_868, \
                         nsg_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = f_3 * pc_x[k] * nsg_868[k];

        t_1211[k] = f_3 * pc_x[k] * nsg_869[k];

        t_1212[k] = f_15 * msg_715[k]
                    + f_1 * nsf0_576[k]
                    - f_2 * nsf1_576[k]
                    + f_3 * pc_y[k] * nsg_865[k];

        t_1213[k] = f_10 * msg_700[k]
                    + f_3 * pc_z[k] * nsg_865[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, pc_y, msg_717, msg_718, msg_719, nsf0_578, \
                         nsf0_579, nsf1_578, nsf1_579, nsg_867, nsg_868, \
                         nsg_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = f_15 * msg_717[k]
                    + f_6 * nsf0_578[k]
                    - f_7 * nsf1_578[k]
                    + f_3 * pc_y[k] * nsg_867[k];

        t_1215[k] = f_15 * msg_718[k]
                    + f_4 * nsf0_579[k]
                    - f_5 * nsf1_579[k]
                    + f_3 * pc_y[k] * nsg_868[k];

        t_1216[k] = f_15 * msg_719[k]
                    + f_3 * pc_y[k] * nsg_869[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, pc_x, pc_z, msg_704, nsf0_579, nsf0_580, \
                         nsf0_581, nsf1_579, nsf1_580, nsf1_581, nsg_869, nsg_870, \
                         nsg_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_10 * msg_704[k]
                    + f_1 * nsf0_579[k]
                    - f_2 * nsf1_579[k]
                    + f_3 * pc_z[k] * nsg_869[k];

        t_1218[k] = f_1 * nsf0_580[k]
                    - f_2 * nsf1_580[k]
                    + f_3 * pc_x[k] * nsg_870[k];

        t_1219[k] = f_13 * nsf0_581[k]
                    - f_14 * nsf1_581[k]
                    + f_3 * pc_x[k] * nsg_871[k];
    }

#pragma omp simd aligned(t_1220, t_1221, t_1222, pc_x, nsf0_582, nsf0_583, nsf0_584, nsf1_582, \
                         nsf1_583, nsf1_584, nsg_872, nsg_873, \
                         nsg_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1220[k] = f_13 * nsf0_582[k]
                    - f_14 * nsf1_582[k]
                    + f_3 * pc_x[k] * nsg_872[k];

        t_1221[k] = f_6 * nsf0_583[k]
                    - f_7 * nsf1_583[k]
                    + f_3 * pc_x[k] * nsg_873[k];

        t_1222[k] = f_6 * nsf0_584[k]
                    - f_7 * nsf1_584[k]
                    + f_3 * pc_x[k] * nsg_874[k];
    }

#pragma omp simd aligned(t_1223, t_1224, t_1225, pc_x, nsf0_585, nsf0_586, nsf0_587, nsf1_585, \
                         nsf1_586, nsf1_587, nsg_875, nsg_876, \
                         nsg_877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1223[k] = f_6 * nsf0_585[k]
                    - f_7 * nsf1_585[k]
                    + f_3 * pc_x[k] * nsg_875[k];

        t_1224[k] = f_4 * nsf0_586[k]
                    - f_5 * nsf1_586[k]
                    + f_3 * pc_x[k] * nsg_876[k];

        t_1225[k] = f_4 * nsf0_587[k]
                    - f_5 * nsf1_587[k]
                    + f_3 * pc_x[k] * nsg_877[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, t_1230, pc_x, nsf0_588, nsf0_589, \
                         nsf1_588, nsf1_589, nsg_878, nsg_879, nsg_880, nsg_881, \
                         nsg_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_4 * nsf0_588[k]
                    - f_5 * nsf1_588[k]
                    + f_3 * pc_x[k] * nsg_878[k];

        t_1227[k] = f_4 * nsf0_589[k]
                    - f_5 * nsf1_589[k]
                    + f_3 * pc_x[k] * nsg_879[k];

        t_1228[k] = f_3 * pc_x[k] * nsg_880[k];

        t_1229[k] = f_3 * pc_x[k] * nsg_881[k];

        t_1230[k] = f_3 * pc_x[k] * nsg_882[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pc_x, pc_y, pc_z, msg_715, msg_730, \
                         nsf0_586, nsf1_586, nsg_880, nsg_883, \
                         nsg_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = f_3 * pc_x[k] * nsg_883[k];

        t_1232[k] = f_3 * pc_x[k] * nsg_884[k];

        t_1233[k] = f_16 * msg_730[k]
                    + f_1 * nsf0_586[k]
                    - f_2 * nsf1_586[k]
                    + f_3 * pc_y[k] * nsg_880[k];

        t_1234[k] = f_11 * msg_715[k]
                    + f_3 * pc_z[k] * nsg_880[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pc_y, msg_732, msg_733, msg_734, nsf0_588, \
                         nsf0_589, nsf1_588, nsf1_589, nsg_882, nsg_883, \
                         nsg_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_16 * msg_732[k]
                    + f_6 * nsf0_588[k]
                    - f_7 * nsf1_588[k]
                    + f_3 * pc_y[k] * nsg_882[k];

        t_1236[k] = f_16 * msg_733[k]
                    + f_4 * nsf0_589[k]
                    - f_5 * nsf1_589[k]
                    + f_3 * pc_y[k] * nsg_883[k];

        t_1237[k] = f_16 * msg_734[k]
                    + f_3 * pc_y[k] * nsg_884[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, pc_x, pc_z, msg_719, nsf0_589, nsf0_590, \
                         nsf0_591, nsf1_589, nsf1_590, nsf1_591, nsg_884, nsg_885, \
                         nsg_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = f_11 * msg_719[k]
                    + f_1 * nsf0_589[k]
                    - f_2 * nsf1_589[k]
                    + f_3 * pc_z[k] * nsg_884[k];

        t_1239[k] = f_1 * nsf0_590[k]
                    - f_2 * nsf1_590[k]
                    + f_3 * pc_x[k] * nsg_885[k];

        t_1240[k] = f_13 * nsf0_591[k]
                    - f_14 * nsf1_591[k]
                    + f_3 * pc_x[k] * nsg_886[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, nsf0_592, nsf0_593, nsf0_594, nsf1_592, \
                         nsf1_593, nsf1_594, nsg_887, nsg_888, \
                         nsg_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_13 * nsf0_592[k]
                    - f_14 * nsf1_592[k]
                    + f_3 * pc_x[k] * nsg_887[k];

        t_1242[k] = f_6 * nsf0_593[k]
                    - f_7 * nsf1_593[k]
                    + f_3 * pc_x[k] * nsg_888[k];

        t_1243[k] = f_6 * nsf0_594[k]
                    - f_7 * nsf1_594[k]
                    + f_3 * pc_x[k] * nsg_889[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, pc_x, nsf0_595, nsf0_596, nsf0_597, nsf1_595, \
                         nsf1_596, nsf1_597, nsg_890, nsg_891, \
                         nsg_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_6 * nsf0_595[k]
                    - f_7 * nsf1_595[k]
                    + f_3 * pc_x[k] * nsg_890[k];

        t_1245[k] = f_4 * nsf0_596[k]
                    - f_5 * nsf1_596[k]
                    + f_3 * pc_x[k] * nsg_891[k];

        t_1246[k] = f_4 * nsf0_597[k]
                    - f_5 * nsf1_597[k]
                    + f_3 * pc_x[k] * nsg_892[k];
    }

#pragma omp simd aligned(t_1247, t_1248, t_1249, t_1250, t_1251, pc_x, nsf0_598, nsf0_599, \
                         nsf1_598, nsf1_599, nsg_893, nsg_894, nsg_895, nsg_896, \
                         nsg_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1247[k] = f_4 * nsf0_598[k]
                    - f_5 * nsf1_598[k]
                    + f_3 * pc_x[k] * nsg_893[k];

        t_1248[k] = f_4 * nsf0_599[k]
                    - f_5 * nsf1_599[k]
                    + f_3 * pc_x[k] * nsg_894[k];

        t_1249[k] = f_3 * pc_x[k] * nsg_895[k];

        t_1250[k] = f_3 * pc_x[k] * nsg_896[k];

        t_1251[k] = f_3 * pc_x[k] * nsg_897[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, pc_x, pc_y, pc_z, msg_730, msg_745, \
                         nsf0_596, nsf1_596, nsg_895, nsg_898, \
                         nsg_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = f_3 * pc_x[k] * nsg_898[k];

        t_1253[k] = f_3 * pc_x[k] * nsg_899[k];

        t_1254[k] = f_17 * msg_745[k]
                    + f_1 * nsf0_596[k]
                    - f_2 * nsf1_596[k]
                    + f_3 * pc_y[k] * nsg_895[k];

        t_1255[k] = f_18 * msg_730[k]
                    + f_3 * pc_z[k] * nsg_895[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, pc_y, msg_747, msg_748, msg_749, nsf0_598, \
                         nsf0_599, nsf1_598, nsf1_599, nsg_897, nsg_898, \
                         nsg_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = f_17 * msg_747[k]
                    + f_6 * nsf0_598[k]
                    - f_7 * nsf1_598[k]
                    + f_3 * pc_y[k] * nsg_897[k];

        t_1257[k] = f_17 * msg_748[k]
                    + f_4 * nsf0_599[k]
                    - f_5 * nsf1_599[k]
                    + f_3 * pc_y[k] * nsg_898[k];

        t_1258[k] = f_17 * msg_749[k]
                    + f_3 * pc_y[k] * nsg_899[k];
    }

#pragma omp simd aligned(t_1259, t_1260, t_1261, pc_x, pc_z, msg_734, nsf0_599, nsf0_600, \
                         nsf0_601, nsf1_599, nsf1_600, nsf1_601, nsg_899, nsg_900, \
                         nsg_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1259[k] = f_18 * msg_734[k]
                    + f_1 * nsf0_599[k]
                    - f_2 * nsf1_599[k]
                    + f_3 * pc_z[k] * nsg_899[k];

        t_1260[k] = f_1 * nsf0_600[k]
                    - f_2 * nsf1_600[k]
                    + f_3 * pc_x[k] * nsg_900[k];

        t_1261[k] = f_13 * nsf0_601[k]
                    - f_14 * nsf1_601[k]
                    + f_3 * pc_x[k] * nsg_901[k];
    }

#pragma omp simd aligned(t_1262, t_1263, t_1264, pc_x, nsf0_602, nsf0_603, nsf0_604, nsf1_602, \
                         nsf1_603, nsf1_604, nsg_902, nsg_903, \
                         nsg_904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1262[k] = f_13 * nsf0_602[k]
                    - f_14 * nsf1_602[k]
                    + f_3 * pc_x[k] * nsg_902[k];

        t_1263[k] = f_6 * nsf0_603[k]
                    - f_7 * nsf1_603[k]
                    + f_3 * pc_x[k] * nsg_903[k];

        t_1264[k] = f_6 * nsf0_604[k]
                    - f_7 * nsf1_604[k]
                    + f_3 * pc_x[k] * nsg_904[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, pc_x, nsf0_605, nsf0_606, nsf0_607, nsf1_605, \
                         nsf1_606, nsf1_607, nsg_905, nsg_906, \
                         nsg_907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = f_6 * nsf0_605[k]
                    - f_7 * nsf1_605[k]
                    + f_3 * pc_x[k] * nsg_905[k];

        t_1266[k] = f_4 * nsf0_606[k]
                    - f_5 * nsf1_606[k]
                    + f_3 * pc_x[k] * nsg_906[k];

        t_1267[k] = f_4 * nsf0_607[k]
                    - f_5 * nsf1_607[k]
                    + f_3 * pc_x[k] * nsg_907[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, t_1271, t_1272, pc_x, nsf0_608, nsf0_609, \
                         nsf1_608, nsf1_609, nsg_908, nsg_909, nsg_910, nsg_911, \
                         nsg_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = f_4 * nsf0_608[k]
                    - f_5 * nsf1_608[k]
                    + f_3 * pc_x[k] * nsg_908[k];

        t_1269[k] = f_4 * nsf0_609[k]
                    - f_5 * nsf1_609[k]
                    + f_3 * pc_x[k] * nsg_909[k];

        t_1270[k] = f_3 * pc_x[k] * nsg_910[k];

        t_1271[k] = f_3 * pc_x[k] * nsg_911[k];

        t_1272[k] = f_3 * pc_x[k] * nsg_912[k];
    }

#pragma omp simd aligned(t_1273, t_1274, t_1275, t_1276, pc_x, pc_y, pc_z, msg_745, msg_760, \
                         nsf0_606, nsf1_606, nsg_910, nsg_913, \
                         nsg_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1273[k] = f_3 * pc_x[k] * nsg_913[k];

        t_1274[k] = f_3 * pc_x[k] * nsg_914[k];

        t_1275[k] = f_19 * msg_760[k]
                    + f_1 * nsf0_606[k]
                    - f_2 * nsf1_606[k]
                    + f_3 * pc_y[k] * nsg_910[k];

        t_1276[k] = f_19 * msg_745[k]
                    + f_3 * pc_z[k] * nsg_910[k];
    }

#pragma omp simd aligned(t_1277, t_1278, t_1279, pc_y, msg_762, msg_763, msg_764, nsf0_608, \
                         nsf0_609, nsf1_608, nsf1_609, nsg_912, nsg_913, \
                         nsg_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1277[k] = f_19 * msg_762[k]
                    + f_6 * nsf0_608[k]
                    - f_7 * nsf1_608[k]
                    + f_3 * pc_y[k] * nsg_912[k];

        t_1278[k] = f_19 * msg_763[k]
                    + f_4 * nsf0_609[k]
                    - f_5 * nsf1_609[k]
                    + f_3 * pc_y[k] * nsg_913[k];

        t_1279[k] = f_19 * msg_764[k]
                    + f_3 * pc_y[k] * nsg_914[k];
    }

#pragma omp simd aligned(t_1280, t_1281, t_1282, pc_x, pc_z, msg_749, nsf0_609, nsf0_610, \
                         nsf0_611, nsf1_609, nsf1_610, nsf1_611, nsg_914, nsg_915, \
                         nsg_916 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1280[k] = f_19 * msg_749[k]
                    + f_1 * nsf0_609[k]
                    - f_2 * nsf1_609[k]
                    + f_3 * pc_z[k] * nsg_914[k];

        t_1281[k] = f_1 * nsf0_610[k]
                    - f_2 * nsf1_610[k]
                    + f_3 * pc_x[k] * nsg_915[k];

        t_1282[k] = f_13 * nsf0_611[k]
                    - f_14 * nsf1_611[k]
                    + f_3 * pc_x[k] * nsg_916[k];
    }

#pragma omp simd aligned(t_1283, t_1284, t_1285, pc_x, nsf0_612, nsf0_613, nsf0_614, nsf1_612, \
                         nsf1_613, nsf1_614, nsg_917, nsg_918, \
                         nsg_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1283[k] = f_13 * nsf0_612[k]
                    - f_14 * nsf1_612[k]
                    + f_3 * pc_x[k] * nsg_917[k];

        t_1284[k] = f_6 * nsf0_613[k]
                    - f_7 * nsf1_613[k]
                    + f_3 * pc_x[k] * nsg_918[k];

        t_1285[k] = f_6 * nsf0_614[k]
                    - f_7 * nsf1_614[k]
                    + f_3 * pc_x[k] * nsg_919[k];
    }

#pragma omp simd aligned(t_1286, t_1287, t_1288, pc_x, nsf0_615, nsf0_616, nsf0_617, nsf1_615, \
                         nsf1_616, nsf1_617, nsg_920, nsg_921, \
                         nsg_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1286[k] = f_6 * nsf0_615[k]
                    - f_7 * nsf1_615[k]
                    + f_3 * pc_x[k] * nsg_920[k];

        t_1287[k] = f_4 * nsf0_616[k]
                    - f_5 * nsf1_616[k]
                    + f_3 * pc_x[k] * nsg_921[k];

        t_1288[k] = f_4 * nsf0_617[k]
                    - f_5 * nsf1_617[k]
                    + f_3 * pc_x[k] * nsg_922[k];
    }

#pragma omp simd aligned(t_1289, t_1290, t_1291, t_1292, t_1293, pc_x, nsf0_618, nsf0_619, \
                         nsf1_618, nsf1_619, nsg_923, nsg_924, nsg_925, nsg_926, \
                         nsg_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1289[k] = f_4 * nsf0_618[k]
                    - f_5 * nsf1_618[k]
                    + f_3 * pc_x[k] * nsg_923[k];

        t_1290[k] = f_4 * nsf0_619[k]
                    - f_5 * nsf1_619[k]
                    + f_3 * pc_x[k] * nsg_924[k];

        t_1291[k] = f_3 * pc_x[k] * nsg_925[k];

        t_1292[k] = f_3 * pc_x[k] * nsg_926[k];

        t_1293[k] = f_3 * pc_x[k] * nsg_927[k];
    }

#pragma omp simd aligned(t_1294, t_1295, t_1296, t_1297, pc_x, pc_y, pc_z, msg_760, msg_775, \
                         nsf0_616, nsf1_616, nsg_925, nsg_928, \
                         nsg_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1294[k] = f_3 * pc_x[k] * nsg_928[k];

        t_1295[k] = f_3 * pc_x[k] * nsg_929[k];

        t_1296[k] = f_18 * msg_775[k]
                    + f_1 * nsf0_616[k]
                    - f_2 * nsf1_616[k]
                    + f_3 * pc_y[k] * nsg_925[k];

        t_1297[k] = f_17 * msg_760[k]
                    + f_3 * pc_z[k] * nsg_925[k];
    }

#pragma omp simd aligned(t_1298, t_1299, t_1300, pc_y, msg_777, msg_778, msg_779, nsf0_618, \
                         nsf0_619, nsf1_618, nsf1_619, nsg_927, nsg_928, \
                         nsg_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1298[k] = f_18 * msg_777[k]
                    + f_6 * nsf0_618[k]
                    - f_7 * nsf1_618[k]
                    + f_3 * pc_y[k] * nsg_927[k];

        t_1299[k] = f_18 * msg_778[k]
                    + f_4 * nsf0_619[k]
                    - f_5 * nsf1_619[k]
                    + f_3 * pc_y[k] * nsg_928[k];

        t_1300[k] = f_18 * msg_779[k]
                    + f_3 * pc_y[k] * nsg_929[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, pc_x, pc_z, msg_764, nsf0_619, nsf0_620, \
                         nsf0_621, nsf1_619, nsf1_620, nsf1_621, nsg_929, nsg_930, \
                         nsg_931 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_17 * msg_764[k]
                    + f_1 * nsf0_619[k]
                    - f_2 * nsf1_619[k]
                    + f_3 * pc_z[k] * nsg_929[k];

        t_1302[k] = f_1 * nsf0_620[k]
                    - f_2 * nsf1_620[k]
                    + f_3 * pc_x[k] * nsg_930[k];

        t_1303[k] = f_13 * nsf0_621[k]
                    - f_14 * nsf1_621[k]
                    + f_3 * pc_x[k] * nsg_931[k];
    }

#pragma omp simd aligned(t_1304, t_1305, t_1306, pc_x, nsf0_622, nsf0_623, nsf0_624, nsf1_622, \
                         nsf1_623, nsf1_624, nsg_932, nsg_933, \
                         nsg_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1304[k] = f_13 * nsf0_622[k]
                    - f_14 * nsf1_622[k]
                    + f_3 * pc_x[k] * nsg_932[k];

        t_1305[k] = f_6 * nsf0_623[k]
                    - f_7 * nsf1_623[k]
                    + f_3 * pc_x[k] * nsg_933[k];

        t_1306[k] = f_6 * nsf0_624[k]
                    - f_7 * nsf1_624[k]
                    + f_3 * pc_x[k] * nsg_934[k];
    }
}

static auto
compute_prim_nsh_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msh0,
                                                           const size_t msg, const size_t msh1,
                                                           const size_t nsf0, const size_t nsf1,
                                                           const size_t nsg, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 4.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_19 = 2.5 / q;

    auto *t_1307 = buffer.data(target + 1307);
    auto *t_1308 = buffer.data(target + 1308);
    auto *t_1309 = buffer.data(target + 1309);
    auto *t_1310 = buffer.data(target + 1310);
    auto *t_1311 = buffer.data(target + 1311);
    auto *t_1312 = buffer.data(target + 1312);
    auto *t_1313 = buffer.data(target + 1313);
    auto *t_1314 = buffer.data(target + 1314);
    auto *t_1315 = buffer.data(target + 1315);
    auto *t_1316 = buffer.data(target + 1316);
    auto *t_1317 = buffer.data(target + 1317);
    auto *t_1318 = buffer.data(target + 1318);
    auto *t_1319 = buffer.data(target + 1319);
    auto *t_1320 = buffer.data(target + 1320);
    auto *t_1321 = buffer.data(target + 1321);
    auto *t_1322 = buffer.data(target + 1322);
    auto *t_1323 = buffer.data(target + 1323);
    auto *t_1324 = buffer.data(target + 1324);
    auto *t_1325 = buffer.data(target + 1325);
    auto *t_1326 = buffer.data(target + 1326);
    auto *t_1327 = buffer.data(target + 1327);
    auto *t_1328 = buffer.data(target + 1328);
    auto *t_1329 = buffer.data(target + 1329);
    auto *t_1330 = buffer.data(target + 1330);
    auto *t_1331 = buffer.data(target + 1331);
    auto *t_1332 = buffer.data(target + 1332);
    auto *t_1333 = buffer.data(target + 1333);
    auto *t_1334 = buffer.data(target + 1334);
    auto *t_1335 = buffer.data(target + 1335);
    auto *t_1336 = buffer.data(target + 1336);
    auto *t_1337 = buffer.data(target + 1337);
    auto *t_1338 = buffer.data(target + 1338);
    auto *t_1339 = buffer.data(target + 1339);
    auto *t_1340 = buffer.data(target + 1340);
    auto *t_1341 = buffer.data(target + 1341);
    auto *t_1342 = buffer.data(target + 1342);
    auto *t_1343 = buffer.data(target + 1343);
    auto *t_1344 = buffer.data(target + 1344);
    auto *t_1345 = buffer.data(target + 1345);
    auto *t_1346 = buffer.data(target + 1346);
    auto *t_1347 = buffer.data(target + 1347);
    auto *t_1348 = buffer.data(target + 1348);
    auto *t_1349 = buffer.data(target + 1349);
    auto *t_1350 = buffer.data(target + 1350);
    auto *t_1351 = buffer.data(target + 1351);
    auto *t_1352 = buffer.data(target + 1352);
    auto *t_1353 = buffer.data(target + 1353);
    auto *t_1354 = buffer.data(target + 1354);
    auto *t_1355 = buffer.data(target + 1355);
    auto *t_1356 = buffer.data(target + 1356);
    auto *t_1357 = buffer.data(target + 1357);
    auto *t_1358 = buffer.data(target + 1358);
    auto *t_1359 = buffer.data(target + 1359);
    auto *t_1360 = buffer.data(target + 1360);
    auto *t_1361 = buffer.data(target + 1361);
    auto *t_1362 = buffer.data(target + 1362);
    auto *t_1363 = buffer.data(target + 1363);
    auto *t_1364 = buffer.data(target + 1364);
    auto *t_1365 = buffer.data(target + 1365);
    auto *t_1366 = buffer.data(target + 1366);
    auto *t_1367 = buffer.data(target + 1367);
    auto *t_1368 = buffer.data(target + 1368);
    auto *t_1369 = buffer.data(target + 1369);
    auto *t_1370 = buffer.data(target + 1370);
    auto *t_1371 = buffer.data(target + 1371);
    auto *t_1372 = buffer.data(target + 1372);
    auto *t_1373 = buffer.data(target + 1373);
    auto *t_1374 = buffer.data(target + 1374);
    auto *t_1375 = buffer.data(target + 1375);
    auto *t_1376 = buffer.data(target + 1376);
    auto *t_1377 = buffer.data(target + 1377);
    auto *t_1378 = buffer.data(target + 1378);
    auto *t_1379 = buffer.data(target + 1379);
    auto *t_1380 = buffer.data(target + 1380);
    auto *t_1381 = buffer.data(target + 1381);
    auto *t_1382 = buffer.data(target + 1382);
    auto *t_1383 = buffer.data(target + 1383);
    auto *t_1384 = buffer.data(target + 1384);
    auto *t_1385 = buffer.data(target + 1385);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh0_1134 = buffer.data(msh0 + 1134);
    const auto *msh0_1136 = buffer.data(msh0 + 1136);
    const auto *msh0_1139 = buffer.data(msh0 + 1139);
    const auto *msh0_1143 = buffer.data(msh0 + 1143);
    const auto *msh0_1149 = buffer.data(msh0 + 1149);
    const auto *msh0_1151 = buffer.data(msh0 + 1151);
    const auto *msh0_1152 = buffer.data(msh0 + 1152);
    const auto *msh0_1154 = buffer.data(msh0 + 1154);

    const auto *msg_775 = buffer.data(msg + 775);
    const auto *msg_779 = buffer.data(msg + 779);
    const auto *msg_790 = buffer.data(msg + 790);
    const auto *msg_792 = buffer.data(msg + 792);
    const auto *msg_793 = buffer.data(msg + 793);
    const auto *msg_794 = buffer.data(msg + 794);
    const auto *msg_805 = buffer.data(msg + 805);
    const auto *msg_807 = buffer.data(msg + 807);
    const auto *msg_808 = buffer.data(msg + 808);
    const auto *msg_809 = buffer.data(msg + 809);
    const auto *msg_820 = buffer.data(msg + 820);
    const auto *msg_822 = buffer.data(msg + 822);
    const auto *msg_823 = buffer.data(msg + 823);
    const auto *msg_824 = buffer.data(msg + 824);

    const auto *msh1_1134 = buffer.data(msh1 + 1134);
    const auto *msh1_1136 = buffer.data(msh1 + 1136);
    const auto *msh1_1139 = buffer.data(msh1 + 1139);
    const auto *msh1_1143 = buffer.data(msh1 + 1143);
    const auto *msh1_1149 = buffer.data(msh1 + 1149);
    const auto *msh1_1151 = buffer.data(msh1 + 1151);
    const auto *msh1_1152 = buffer.data(msh1 + 1152);
    const auto *msh1_1154 = buffer.data(msh1 + 1154);

    const auto *nsf0_625 = buffer.data(nsf0 + 625);
    const auto *nsf0_626 = buffer.data(nsf0 + 626);
    const auto *nsf0_627 = buffer.data(nsf0 + 627);
    const auto *nsf0_628 = buffer.data(nsf0 + 628);
    const auto *nsf0_629 = buffer.data(nsf0 + 629);
    const auto *nsf0_630 = buffer.data(nsf0 + 630);
    const auto *nsf0_631 = buffer.data(nsf0 + 631);
    const auto *nsf0_632 = buffer.data(nsf0 + 632);
    const auto *nsf0_633 = buffer.data(nsf0 + 633);
    const auto *nsf0_634 = buffer.data(nsf0 + 634);
    const auto *nsf0_635 = buffer.data(nsf0 + 635);
    const auto *nsf0_636 = buffer.data(nsf0 + 636);
    const auto *nsf0_637 = buffer.data(nsf0 + 637);
    const auto *nsf0_638 = buffer.data(nsf0 + 638);
    const auto *nsf0_639 = buffer.data(nsf0 + 639);
    const auto *nsf0_641 = buffer.data(nsf0 + 641);
    const auto *nsf0_643 = buffer.data(nsf0 + 643);
    const auto *nsf0_644 = buffer.data(nsf0 + 644);
    const auto *nsf0_646 = buffer.data(nsf0 + 646);
    const auto *nsf0_647 = buffer.data(nsf0 + 647);
    const auto *nsf0_648 = buffer.data(nsf0 + 648);
    const auto *nsf0_650 = buffer.data(nsf0 + 650);
    const auto *nsf0_652 = buffer.data(nsf0 + 652);
    const auto *nsf0_653 = buffer.data(nsf0 + 653);
    const auto *nsf0_655 = buffer.data(nsf0 + 655);
    const auto *nsf0_656 = buffer.data(nsf0 + 656);
    const auto *nsf0_657 = buffer.data(nsf0 + 657);
    const auto *nsf0_658 = buffer.data(nsf0 + 658);
    const auto *nsf0_659 = buffer.data(nsf0 + 659);

    const auto *nsf1_625 = buffer.data(nsf1 + 625);
    const auto *nsf1_626 = buffer.data(nsf1 + 626);
    const auto *nsf1_627 = buffer.data(nsf1 + 627);
    const auto *nsf1_628 = buffer.data(nsf1 + 628);
    const auto *nsf1_629 = buffer.data(nsf1 + 629);
    const auto *nsf1_630 = buffer.data(nsf1 + 630);
    const auto *nsf1_631 = buffer.data(nsf1 + 631);
    const auto *nsf1_632 = buffer.data(nsf1 + 632);
    const auto *nsf1_633 = buffer.data(nsf1 + 633);
    const auto *nsf1_634 = buffer.data(nsf1 + 634);
    const auto *nsf1_635 = buffer.data(nsf1 + 635);
    const auto *nsf1_636 = buffer.data(nsf1 + 636);
    const auto *nsf1_637 = buffer.data(nsf1 + 637);
    const auto *nsf1_638 = buffer.data(nsf1 + 638);
    const auto *nsf1_639 = buffer.data(nsf1 + 639);
    const auto *nsf1_641 = buffer.data(nsf1 + 641);
    const auto *nsf1_643 = buffer.data(nsf1 + 643);
    const auto *nsf1_644 = buffer.data(nsf1 + 644);
    const auto *nsf1_646 = buffer.data(nsf1 + 646);
    const auto *nsf1_647 = buffer.data(nsf1 + 647);
    const auto *nsf1_648 = buffer.data(nsf1 + 648);
    const auto *nsf1_650 = buffer.data(nsf1 + 650);
    const auto *nsf1_652 = buffer.data(nsf1 + 652);
    const auto *nsf1_653 = buffer.data(nsf1 + 653);
    const auto *nsf1_655 = buffer.data(nsf1 + 655);
    const auto *nsf1_656 = buffer.data(nsf1 + 656);
    const auto *nsf1_657 = buffer.data(nsf1 + 657);
    const auto *nsf1_658 = buffer.data(nsf1 + 658);
    const auto *nsf1_659 = buffer.data(nsf1 + 659);

    const auto *nsg_935 = buffer.data(nsg + 935);
    const auto *nsg_936 = buffer.data(nsg + 936);
    const auto *nsg_937 = buffer.data(nsg + 937);
    const auto *nsg_938 = buffer.data(nsg + 938);
    const auto *nsg_939 = buffer.data(nsg + 939);
    const auto *nsg_940 = buffer.data(nsg + 940);
    const auto *nsg_941 = buffer.data(nsg + 941);
    const auto *nsg_942 = buffer.data(nsg + 942);
    const auto *nsg_943 = buffer.data(nsg + 943);
    const auto *nsg_944 = buffer.data(nsg + 944);
    const auto *nsg_945 = buffer.data(nsg + 945);
    const auto *nsg_946 = buffer.data(nsg + 946);
    const auto *nsg_947 = buffer.data(nsg + 947);
    const auto *nsg_948 = buffer.data(nsg + 948);
    const auto *nsg_949 = buffer.data(nsg + 949);
    const auto *nsg_950 = buffer.data(nsg + 950);
    const auto *nsg_951 = buffer.data(nsg + 951);
    const auto *nsg_952 = buffer.data(nsg + 952);
    const auto *nsg_953 = buffer.data(nsg + 953);
    const auto *nsg_954 = buffer.data(nsg + 954);
    const auto *nsg_955 = buffer.data(nsg + 955);
    const auto *nsg_956 = buffer.data(nsg + 956);
    const auto *nsg_957 = buffer.data(nsg + 957);
    const auto *nsg_958 = buffer.data(nsg + 958);
    const auto *nsg_959 = buffer.data(nsg + 959);
    const auto *nsg_961 = buffer.data(nsg + 961);
    const auto *nsg_963 = buffer.data(nsg + 963);
    const auto *nsg_964 = buffer.data(nsg + 964);
    const auto *nsg_966 = buffer.data(nsg + 966);
    const auto *nsg_967 = buffer.data(nsg + 967);
    const auto *nsg_968 = buffer.data(nsg + 968);
    const auto *nsg_970 = buffer.data(nsg + 970);
    const auto *nsg_971 = buffer.data(nsg + 971);
    const auto *nsg_972 = buffer.data(nsg + 972);
    const auto *nsg_973 = buffer.data(nsg + 973);
    const auto *nsg_974 = buffer.data(nsg + 974);
    const auto *nsg_975 = buffer.data(nsg + 975);
    const auto *nsg_977 = buffer.data(nsg + 977);
    const auto *nsg_978 = buffer.data(nsg + 978);
    const auto *nsg_980 = buffer.data(nsg + 980);
    const auto *nsg_981 = buffer.data(nsg + 981);
    const auto *nsg_982 = buffer.data(nsg + 982);
    const auto *nsg_984 = buffer.data(nsg + 984);
    const auto *nsg_985 = buffer.data(nsg + 985);
    const auto *nsg_986 = buffer.data(nsg + 986);
    const auto *nsg_987 = buffer.data(nsg + 987);
    const auto *nsg_988 = buffer.data(nsg + 988);
    const auto *nsg_989 = buffer.data(nsg + 989);

#pragma omp simd aligned(t_1307, t_1308, t_1309, pc_x, nsf0_625, nsf0_626, nsf0_627, nsf1_625, \
                         nsf1_626, nsf1_627, nsg_935, nsg_936, \
                         nsg_937 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1307[k] = f_6 * nsf0_625[k]
                    - f_7 * nsf1_625[k]
                    + f_3 * pc_x[k] * nsg_935[k];

        t_1308[k] = f_4 * nsf0_626[k]
                    - f_5 * nsf1_626[k]
                    + f_3 * pc_x[k] * nsg_936[k];

        t_1309[k] = f_4 * nsf0_627[k]
                    - f_5 * nsf1_627[k]
                    + f_3 * pc_x[k] * nsg_937[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, t_1313, t_1314, pc_x, nsf0_628, nsf0_629, \
                         nsf1_628, nsf1_629, nsg_938, nsg_939, nsg_940, nsg_941, \
                         nsg_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_4 * nsf0_628[k]
                    - f_5 * nsf1_628[k]
                    + f_3 * pc_x[k] * nsg_938[k];

        t_1311[k] = f_4 * nsf0_629[k]
                    - f_5 * nsf1_629[k]
                    + f_3 * pc_x[k] * nsg_939[k];

        t_1312[k] = f_3 * pc_x[k] * nsg_940[k];

        t_1313[k] = f_3 * pc_x[k] * nsg_941[k];

        t_1314[k] = f_3 * pc_x[k] * nsg_942[k];
    }

#pragma omp simd aligned(t_1315, t_1316, t_1317, t_1318, pc_x, pc_y, pc_z, msg_775, msg_790, \
                         nsf0_626, nsf1_626, nsg_940, nsg_943, \
                         nsg_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1315[k] = f_3 * pc_x[k] * nsg_943[k];

        t_1316[k] = f_3 * pc_x[k] * nsg_944[k];

        t_1317[k] = f_11 * msg_790[k]
                    + f_1 * nsf0_626[k]
                    - f_2 * nsf1_626[k]
                    + f_3 * pc_y[k] * nsg_940[k];

        t_1318[k] = f_16 * msg_775[k]
                    + f_3 * pc_z[k] * nsg_940[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, pc_y, msg_792, msg_793, msg_794, nsf0_628, \
                         nsf0_629, nsf1_628, nsf1_629, nsg_942, nsg_943, \
                         nsg_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = f_11 * msg_792[k]
                    + f_6 * nsf0_628[k]
                    - f_7 * nsf1_628[k]
                    + f_3 * pc_y[k] * nsg_942[k];

        t_1320[k] = f_11 * msg_793[k]
                    + f_4 * nsf0_629[k]
                    - f_5 * nsf1_629[k]
                    + f_3 * pc_y[k] * nsg_943[k];

        t_1321[k] = f_11 * msg_794[k]
                    + f_3 * pc_y[k] * nsg_944[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, pc_x, pc_z, msg_779, nsf0_629, nsf0_630, \
                         nsf0_631, nsf1_629, nsf1_630, nsf1_631, nsg_944, nsg_945, \
                         nsg_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = f_16 * msg_779[k]
                    + f_1 * nsf0_629[k]
                    - f_2 * nsf1_629[k]
                    + f_3 * pc_z[k] * nsg_944[k];

        t_1323[k] = f_1 * nsf0_630[k]
                    - f_2 * nsf1_630[k]
                    + f_3 * pc_x[k] * nsg_945[k];

        t_1324[k] = f_13 * nsf0_631[k]
                    - f_14 * nsf1_631[k]
                    + f_3 * pc_x[k] * nsg_946[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pc_x, nsf0_632, nsf0_633, nsf0_634, nsf1_632, \
                         nsf1_633, nsf1_634, nsg_947, nsg_948, \
                         nsg_949 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_13 * nsf0_632[k]
                    - f_14 * nsf1_632[k]
                    + f_3 * pc_x[k] * nsg_947[k];

        t_1326[k] = f_6 * nsf0_633[k]
                    - f_7 * nsf1_633[k]
                    + f_3 * pc_x[k] * nsg_948[k];

        t_1327[k] = f_6 * nsf0_634[k]
                    - f_7 * nsf1_634[k]
                    + f_3 * pc_x[k] * nsg_949[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, pc_x, nsf0_635, nsf0_636, nsf0_637, nsf1_635, \
                         nsf1_636, nsf1_637, nsg_950, nsg_951, \
                         nsg_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_6 * nsf0_635[k]
                    - f_7 * nsf1_635[k]
                    + f_3 * pc_x[k] * nsg_950[k];

        t_1329[k] = f_4 * nsf0_636[k]
                    - f_5 * nsf1_636[k]
                    + f_3 * pc_x[k] * nsg_951[k];

        t_1330[k] = f_4 * nsf0_637[k]
                    - f_5 * nsf1_637[k]
                    + f_3 * pc_x[k] * nsg_952[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, t_1335, pc_x, nsf0_638, nsf0_639, \
                         nsf1_638, nsf1_639, nsg_953, nsg_954, nsg_955, nsg_956, \
                         nsg_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_4 * nsf0_638[k]
                    - f_5 * nsf1_638[k]
                    + f_3 * pc_x[k] * nsg_953[k];

        t_1332[k] = f_4 * nsf0_639[k]
                    - f_5 * nsf1_639[k]
                    + f_3 * pc_x[k] * nsg_954[k];

        t_1333[k] = f_3 * pc_x[k] * nsg_955[k];

        t_1334[k] = f_3 * pc_x[k] * nsg_956[k];

        t_1335[k] = f_3 * pc_x[k] * nsg_957[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, t_1339, pc_x, pc_y, pc_z, msg_790, msg_805, \
                         nsf0_636, nsf1_636, nsg_955, nsg_958, \
                         nsg_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_3 * pc_x[k] * nsg_958[k];

        t_1337[k] = f_3 * pc_x[k] * nsg_959[k];

        t_1338[k] = f_10 * msg_805[k]
                    + f_1 * nsf0_636[k]
                    - f_2 * nsf1_636[k]
                    + f_3 * pc_y[k] * nsg_955[k];

        t_1339[k] = f_15 * msg_790[k]
                    + f_3 * pc_z[k] * nsg_955[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, pc_y, msg_807, msg_808, msg_809, nsf0_638, \
                         nsf0_639, nsf1_638, nsf1_639, nsg_957, nsg_958, \
                         nsg_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = f_10 * msg_807[k]
                    + f_6 * nsf0_638[k]
                    - f_7 * nsf1_638[k]
                    + f_3 * pc_y[k] * nsg_957[k];

        t_1341[k] = f_10 * msg_808[k]
                    + f_4 * nsf0_639[k]
                    - f_5 * nsf1_639[k]
                    + f_3 * pc_y[k] * nsg_958[k];

        t_1342[k] = f_10 * msg_809[k]
                    + f_3 * pc_y[k] * nsg_959[k];
    }

#pragma omp simd aligned(t_1343, t_1344, t_1345, pa_y, pc_x, pc_y, pc_z, msh0_1134, msg_794, \
                         msh1_1134, nsf0_639, nsf0_641, nsf1_639, nsf1_641, nsg_959, \
                         nsg_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1343[k] = f_15 * msg_794[k]
                    + f_1 * nsf0_639[k]
                    - f_2 * nsf1_639[k]
                    + f_3 * pc_z[k] * nsg_959[k];

        t_1344[k] = pa_y[k] * msh0_1134[k]
                    - f_8 * pc_y[k] * msh1_1134[k];

        t_1345[k] = f_13 * nsf0_641[k]
                    - f_14 * nsf1_641[k]
                    + f_3 * pc_x[k] * nsg_961[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, pa_y, pc_x, pc_y, msh0_1136, msh1_1136, \
                         nsf0_643, nsf0_644, nsf1_643, nsf1_644, nsg_963, \
                         nsg_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = pa_y[k] * msh0_1136[k]
                    - f_8 * pc_y[k] * msh1_1136[k];

        t_1347[k] = f_6 * nsf0_643[k]
                    - f_7 * nsf1_643[k]
                    + f_3 * pc_x[k] * nsg_963[k];

        t_1348[k] = f_6 * nsf0_644[k]
                    - f_7 * nsf1_644[k]
                    + f_3 * pc_x[k] * nsg_964[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pa_y, pc_x, pc_y, msh0_1139, msh1_1139, \
                         nsf0_646, nsf0_647, nsf1_646, nsf1_647, nsg_966, \
                         nsg_967 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = pa_y[k] * msh0_1139[k]
                    - f_8 * pc_y[k] * msh1_1139[k];

        t_1350[k] = f_4 * nsf0_646[k]
                    - f_5 * nsf1_646[k]
                    + f_3 * pc_x[k] * nsg_966[k];

        t_1351[k] = f_4 * nsf0_647[k]
                    - f_5 * nsf1_647[k]
                    + f_3 * pc_x[k] * nsg_967[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, t_1355, t_1356, pa_y, pc_x, pc_y, msh0_1143, \
                         msh1_1143, nsf0_648, nsf1_648, nsg_968, nsg_970, nsg_971, \
                         nsg_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_4 * nsf0_648[k]
                    - f_5 * nsf1_648[k]
                    + f_3 * pc_x[k] * nsg_968[k];

        t_1353[k] = pa_y[k] * msh0_1143[k]
                    - f_8 * pc_y[k] * msh1_1143[k];

        t_1354[k] = f_3 * pc_x[k] * nsg_970[k];

        t_1355[k] = f_3 * pc_x[k] * nsg_971[k];

        t_1356[k] = f_3 * pc_x[k] * nsg_972[k];
    }

#pragma omp simd aligned(t_1357, t_1358, t_1359, t_1360, pa_y, pc_x, pc_y, pc_z, msh0_1149, \
                         msg_805, msg_820, msh1_1149, nsg_970, nsg_973, \
                         nsg_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1357[k] = f_3 * pc_x[k] * nsg_973[k];

        t_1358[k] = f_3 * pc_x[k] * nsg_974[k];

        t_1359[k] = pa_y[k] * msh0_1149[k]
                    + f_19 * msg_820[k]
                    - f_8 * pc_y[k] * msh1_1149[k];

        t_1360[k] = f_12 * msg_805[k]
                    + f_3 * pc_z[k] * nsg_970[k];
    }

#pragma omp simd aligned(t_1361, t_1362, t_1363, t_1364, pa_y, pc_y, msh0_1151, msh0_1152, \
                         msh0_1154, msg_822, msg_823, msg_824, msh1_1151, msh1_1152, \
                         msh1_1154, nsg_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1361[k] = pa_y[k] * msh0_1151[k]
                    + f_11 * msg_822[k]
                    - f_8 * pc_y[k] * msh1_1151[k];

        t_1362[k] = pa_y[k] * msh0_1152[k]
                    + f_10 * msg_823[k]
                    - f_8 * pc_y[k] * msh1_1152[k];

        t_1363[k] = f_9 * msg_824[k]
                    + f_3 * pc_y[k] * nsg_974[k];

        t_1364[k] = pa_y[k] * msh0_1154[k]
                    - f_8 * pc_y[k] * msh1_1154[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, t_1368, t_1369, pc_x, pc_y, nsf0_650, \
                         nsf0_652, nsf0_653, nsf1_650, nsf1_652, nsf1_653, nsg_975, nsg_977, \
                         nsg_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = f_1 * nsf0_650[k]
                    - f_2 * nsf1_650[k]
                    + f_3 * pc_x[k] * nsg_975[k];

        t_1366[k] = f_3 * pc_y[k] * nsg_975[k];

        t_1367[k] = f_13 * nsf0_652[k]
                    - f_14 * nsf1_652[k]
                    + f_3 * pc_x[k] * nsg_977[k];

        t_1368[k] = f_6 * nsf0_653[k]
                    - f_7 * nsf1_653[k]
                    + f_3 * pc_x[k] * nsg_978[k];

        t_1369[k] = f_3 * pc_y[k] * nsg_977[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, t_1373, pc_x, pc_y, nsf0_655, nsf0_656, \
                         nsf0_657, nsf1_655, nsf1_656, nsf1_657, nsg_980, nsg_981, \
                         nsg_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_6 * nsf0_655[k]
                    - f_7 * nsf1_655[k]
                    + f_3 * pc_x[k] * nsg_980[k];

        t_1371[k] = f_4 * nsf0_656[k]
                    - f_5 * nsf1_656[k]
                    + f_3 * pc_x[k] * nsg_981[k];

        t_1372[k] = f_4 * nsf0_657[k]
                    - f_5 * nsf1_657[k]
                    + f_3 * pc_x[k] * nsg_982[k];

        t_1373[k] = f_3 * pc_y[k] * nsg_980[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, t_1377, t_1378, t_1379, pc_x, nsf0_659, \
                         nsf1_659, nsg_984, nsg_985, nsg_986, nsg_987, nsg_988, \
                         nsg_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_4 * nsf0_659[k]
                    - f_5 * nsf1_659[k]
                    + f_3 * pc_x[k] * nsg_984[k];

        t_1375[k] = f_3 * pc_x[k] * nsg_985[k];

        t_1376[k] = f_3 * pc_x[k] * nsg_986[k];

        t_1377[k] = f_3 * pc_x[k] * nsg_987[k];

        t_1378[k] = f_3 * pc_x[k] * nsg_988[k];

        t_1379[k] = f_3 * pc_x[k] * nsg_989[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, pc_y, nsf0_656, nsf0_657, nsf0_658, nsf1_656, \
                         nsf1_657, nsf1_658, nsg_985, nsg_986, \
                         nsg_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_1 * nsf0_656[k]
                    - f_2 * nsf1_656[k]
                    + f_3 * pc_y[k] * nsg_985[k];

        t_1381[k] = f_13 * nsf0_657[k]
                    - f_14 * nsf1_657[k]
                    + f_3 * pc_y[k] * nsg_986[k];

        t_1382[k] = f_6 * nsf0_658[k]
                    - f_7 * nsf1_658[k]
                    + f_3 * pc_y[k] * nsg_987[k];
    }

#pragma omp simd aligned(t_1383, t_1384, t_1385, pc_y, pc_z, msg_824, nsf0_659, nsf1_659, \
                         nsg_988, nsg_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1383[k] = f_4 * nsf0_659[k]
                    - f_5 * nsf1_659[k]
                    + f_3 * pc_y[k] * nsg_988[k];

        t_1384[k] = f_3 * pc_y[k] * nsg_989[k];

        t_1385[k] = f_0 * msg_824[k]
                    + f_1 * nsf0_659[k]
                    - f_2 * nsf1_659[k]
                    + f_3 * pc_z[k] * nsg_989[k];
    }
}

auto
compute_prim_nsh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t msh0, const size_t msg,
                                                   const size_t msh1, const size_t nsf0,
                                                   const size_t nsf1, const size_t nsg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_nsh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, msh0, msg,
                                                              msh1, nsf0, nsf1, nsg, ncols,
                                                              gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, msh0,
                                                               msg, msh1, nsf0, nsf1, nsg,
                                                               ncols, gamma, p, q);

    compute_prim_nsh_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, msh0,
                                                               msg, msh1, nsf0, nsf1, nsg,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
