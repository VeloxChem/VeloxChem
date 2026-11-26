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


#include "SimdThreeCenterElectronRepulsionVrrRecPPH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pph_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sph0, const size_t spg,
                                                          const size_t sph1, const size_t psh0,
                                                          const size_t psg, const size_t psh1,
                                                          const size_t ppf0, const size_t ppf1,
                                                          const size_t ppg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / gamma;
    const auto f_12 = 1.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sph0_0 = buffer.data(sph0 + 0);
    const auto *sph0_1 = buffer.data(sph0 + 1);
    const auto *sph0_2 = buffer.data(sph0 + 2);
    const auto *sph0_3 = buffer.data(sph0 + 3);
    const auto *sph0_5 = buffer.data(sph0 + 5);
    const auto *sph0_6 = buffer.data(sph0 + 6);
    const auto *sph0_8 = buffer.data(sph0 + 8);
    const auto *sph0_9 = buffer.data(sph0 + 9);
    const auto *sph0_20 = buffer.data(sph0 + 20);
    const auto *sph0_24 = buffer.data(sph0 + 24);
    const auto *sph0_27 = buffer.data(sph0 + 27);
    const auto *sph0_36 = buffer.data(sph0 + 36);
    const auto *sph0_38 = buffer.data(sph0 + 38);
    const auto *sph0_39 = buffer.data(sph0 + 39);
    const auto *sph0_41 = buffer.data(sph0 + 41);
    const auto *sph0_42 = buffer.data(sph0 + 42);
    const auto *sph0_47 = buffer.data(sph0 + 47);
    const auto *sph0_51 = buffer.data(sph0 + 51);
    const auto *sph0_57 = buffer.data(sph0 + 57);
    const auto *sph0_59 = buffer.data(sph0 + 59);
    const auto *sph0_60 = buffer.data(sph0 + 60);
    const auto *sph0_62 = buffer.data(sph0 + 62);

    const auto *spg_0 = buffer.data(spg + 0);
    const auto *spg_1 = buffer.data(spg + 1);
    const auto *spg_3 = buffer.data(spg + 3);
    const auto *spg_5 = buffer.data(spg + 5);
    const auto *spg_10 = buffer.data(spg + 10);
    const auto *spg_12 = buffer.data(spg + 12);
    const auto *spg_14 = buffer.data(spg + 14);
    const auto *spg_18 = buffer.data(spg + 18);
    const auto *spg_21 = buffer.data(spg + 21);
    const auto *spg_25 = buffer.data(spg + 25);
    const auto *spg_27 = buffer.data(spg + 27);
    const auto *spg_29 = buffer.data(spg + 29);
    const auto *spg_35 = buffer.data(spg + 35);
    const auto *spg_39 = buffer.data(spg + 39);
    const auto *spg_40 = buffer.data(spg + 40);
    const auto *spg_42 = buffer.data(spg + 42);
    const auto *spg_43 = buffer.data(spg + 43);
    const auto *spg_44 = buffer.data(spg + 44);

    const auto *sph1_0 = buffer.data(sph1 + 0);
    const auto *sph1_1 = buffer.data(sph1 + 1);
    const auto *sph1_2 = buffer.data(sph1 + 2);
    const auto *sph1_3 = buffer.data(sph1 + 3);
    const auto *sph1_5 = buffer.data(sph1 + 5);
    const auto *sph1_6 = buffer.data(sph1 + 6);
    const auto *sph1_8 = buffer.data(sph1 + 8);
    const auto *sph1_9 = buffer.data(sph1 + 9);
    const auto *sph1_20 = buffer.data(sph1 + 20);
    const auto *sph1_24 = buffer.data(sph1 + 24);
    const auto *sph1_27 = buffer.data(sph1 + 27);
    const auto *sph1_36 = buffer.data(sph1 + 36);
    const auto *sph1_38 = buffer.data(sph1 + 38);
    const auto *sph1_39 = buffer.data(sph1 + 39);
    const auto *sph1_41 = buffer.data(sph1 + 41);
    const auto *sph1_42 = buffer.data(sph1 + 42);
    const auto *sph1_47 = buffer.data(sph1 + 47);
    const auto *sph1_51 = buffer.data(sph1 + 51);
    const auto *sph1_57 = buffer.data(sph1 + 57);
    const auto *sph1_59 = buffer.data(sph1 + 59);
    const auto *sph1_60 = buffer.data(sph1 + 60);
    const auto *sph1_62 = buffer.data(sph1 + 62);

    const auto *psh0_0 = buffer.data(psh0 + 0);
    const auto *psh0_3 = buffer.data(psh0 + 3);
    const auto *psh0_5 = buffer.data(psh0 + 5);
    const auto *psh0_6 = buffer.data(psh0 + 6);
    const auto *psh0_9 = buffer.data(psh0 + 9);
    const auto *psh0_22 = buffer.data(psh0 + 22);
    const auto *psh0_24 = buffer.data(psh0 + 24);
    const auto *psh0_27 = buffer.data(psh0 + 27);
    const auto *psh0_36 = buffer.data(psh0 + 36);
    const auto *psh0_38 = buffer.data(psh0 + 38);
    const auto *psh0_39 = buffer.data(psh0 + 39);

    const auto *psg_0 = buffer.data(psg + 0);
    const auto *psg_2 = buffer.data(psg + 2);
    const auto *psg_3 = buffer.data(psg + 3);
    const auto *psg_5 = buffer.data(psg + 5);
    const auto *psg_6 = buffer.data(psg + 6);
    const auto *psg_9 = buffer.data(psg + 9);
    const auto *psg_10 = buffer.data(psg + 10);
    const auto *psg_12 = buffer.data(psg + 12);
    const auto *psg_14 = buffer.data(psg + 14);
    const auto *psg_15 = buffer.data(psg + 15);
    const auto *psg_16 = buffer.data(psg + 16);
    const auto *psg_18 = buffer.data(psg + 18);
    const auto *psg_25 = buffer.data(psg + 25);
    const auto *psg_26 = buffer.data(psg + 26);
    const auto *psg_27 = buffer.data(psg + 27);
    const auto *psg_28 = buffer.data(psg + 28);
    const auto *psg_29 = buffer.data(psg + 29);

    const auto *psh1_0 = buffer.data(psh1 + 0);
    const auto *psh1_3 = buffer.data(psh1 + 3);
    const auto *psh1_5 = buffer.data(psh1 + 5);
    const auto *psh1_6 = buffer.data(psh1 + 6);
    const auto *psh1_9 = buffer.data(psh1 + 9);
    const auto *psh1_22 = buffer.data(psh1 + 22);
    const auto *psh1_24 = buffer.data(psh1 + 24);
    const auto *psh1_27 = buffer.data(psh1 + 27);
    const auto *psh1_36 = buffer.data(psh1 + 36);
    const auto *psh1_38 = buffer.data(psh1 + 38);
    const auto *psh1_39 = buffer.data(psh1 + 39);

    const auto *ppf0_0 = buffer.data(ppf0 + 0);
    const auto *ppf0_1 = buffer.data(ppf0 + 1);
    const auto *ppf0_2 = buffer.data(ppf0 + 2);
    const auto *ppf0_6 = buffer.data(ppf0 + 6);
    const auto *ppf0_8 = buffer.data(ppf0 + 8);
    const auto *ppf0_9 = buffer.data(ppf0 + 9);
    const auto *ppf0_40 = buffer.data(ppf0 + 40);
    const auto *ppf0_41 = buffer.data(ppf0 + 41);
    const auto *ppf0_43 = buffer.data(ppf0 + 43);
    const auto *ppf0_45 = buffer.data(ppf0 + 45);
    const auto *ppf0_46 = buffer.data(ppf0 + 46);
    const auto *ppf0_47 = buffer.data(ppf0 + 47);
    const auto *ppf0_48 = buffer.data(ppf0 + 48);
    const auto *ppf0_49 = buffer.data(ppf0 + 49);
    const auto *ppf0_58 = buffer.data(ppf0 + 58);

    const auto *ppf1_0 = buffer.data(ppf1 + 0);
    const auto *ppf1_1 = buffer.data(ppf1 + 1);
    const auto *ppf1_2 = buffer.data(ppf1 + 2);
    const auto *ppf1_6 = buffer.data(ppf1 + 6);
    const auto *ppf1_8 = buffer.data(ppf1 + 8);
    const auto *ppf1_9 = buffer.data(ppf1 + 9);
    const auto *ppf1_40 = buffer.data(ppf1 + 40);
    const auto *ppf1_41 = buffer.data(ppf1 + 41);
    const auto *ppf1_43 = buffer.data(ppf1 + 43);
    const auto *ppf1_45 = buffer.data(ppf1 + 45);
    const auto *ppf1_46 = buffer.data(ppf1 + 46);
    const auto *ppf1_47 = buffer.data(ppf1 + 47);
    const auto *ppf1_48 = buffer.data(ppf1 + 48);
    const auto *ppf1_49 = buffer.data(ppf1 + 49);
    const auto *ppf1_58 = buffer.data(ppf1 + 58);

    const auto *ppg_0 = buffer.data(ppg + 0);
    const auto *ppg_1 = buffer.data(ppg + 1);
    const auto *ppg_2 = buffer.data(ppg + 2);
    const auto *ppg_3 = buffer.data(ppg + 3);
    const auto *ppg_5 = buffer.data(ppg + 5);
    const auto *ppg_6 = buffer.data(ppg + 6);
    const auto *ppg_9 = buffer.data(ppg + 9);
    const auto *ppg_10 = buffer.data(ppg + 10);
    const auto *ppg_12 = buffer.data(ppg + 12);
    const auto *ppg_13 = buffer.data(ppg + 13);
    const auto *ppg_14 = buffer.data(ppg + 14);
    const auto *ppg_15 = buffer.data(ppg + 15);
    const auto *ppg_17 = buffer.data(ppg + 17);
    const auto *ppg_18 = buffer.data(ppg + 18);
    const auto *ppg_20 = buffer.data(ppg + 20);
    const auto *ppg_21 = buffer.data(ppg + 21);
    const auto *ppg_24 = buffer.data(ppg + 24);
    const auto *ppg_25 = buffer.data(ppg + 25);
    const auto *ppg_27 = buffer.data(ppg + 27);
    const auto *ppg_29 = buffer.data(ppg + 29);
    const auto *ppg_30 = buffer.data(ppg + 30);
    const auto *ppg_32 = buffer.data(ppg + 32);
    const auto *ppg_33 = buffer.data(ppg + 33);
    const auto *ppg_35 = buffer.data(ppg + 35);
    const auto *ppg_36 = buffer.data(ppg + 36);
    const auto *ppg_39 = buffer.data(ppg + 39);
    const auto *ppg_40 = buffer.data(ppg + 40);
    const auto *ppg_42 = buffer.data(ppg + 42);
    const auto *ppg_44 = buffer.data(ppg + 44);
    const auto *ppg_45 = buffer.data(ppg + 45);
    const auto *ppg_46 = buffer.data(ppg + 46);
    const auto *ppg_48 = buffer.data(ppg + 48);
    const auto *ppg_55 = buffer.data(ppg + 55);
    const auto *ppg_56 = buffer.data(ppg + 56);
    const auto *ppg_57 = buffer.data(ppg + 57);
    const auto *ppg_58 = buffer.data(ppg + 58);
    const auto *ppg_59 = buffer.data(ppg + 59);
    const auto *ppg_60 = buffer.data(ppg + 60);
    const auto *ppg_61 = buffer.data(ppg + 61);
    const auto *ppg_63 = buffer.data(ppg + 63);
    const auto *ppg_65 = buffer.data(ppg + 65);
    const auto *ppg_66 = buffer.data(ppg + 66);
    const auto *ppg_68 = buffer.data(ppg + 68);
    const auto *ppg_69 = buffer.data(ppg + 69);
    const auto *ppg_70 = buffer.data(ppg + 70);
    const auto *ppg_71 = buffer.data(ppg + 71);
    const auto *ppg_72 = buffer.data(ppg + 72);
    const auto *ppg_73 = buffer.data(ppg + 73);
    const auto *ppg_74 = buffer.data(ppg + 74);
    const auto *ppg_75 = buffer.data(ppg + 75);
    const auto *ppg_76 = buffer.data(ppg + 76);
    const auto *ppg_78 = buffer.data(ppg + 78);
    const auto *ppg_83 = buffer.data(ppg + 83);
    const auto *ppg_85 = buffer.data(ppg + 85);
    const auto *ppg_86 = buffer.data(ppg + 86);
    const auto *ppg_87 = buffer.data(ppg + 87);
    const auto *ppg_88 = buffer.data(ppg + 88);
    const auto *ppg_89 = buffer.data(ppg + 89);
    const auto *ppg_90 = buffer.data(ppg + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, spg_0, psg_0, ppf0_0, \
                         ppf1_0, ppg_0, ppg_1, ppg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spg_0[k]
                 + f_0 * psg_0[k]
                 + f_1 * ppf0_0[k]
                 - f_2 * ppf1_0[k]
                 + f_3 * pc_x[k] * ppg_0[k];

        t_1[k] = f_3 * pc_y[k] * ppg_0[k];

        t_2[k] = f_3 * pc_z[k] * ppg_0[k];

        t_3[k] = f_4 * ppf0_0[k]
                 - f_5 * ppf1_0[k]
                 + f_3 * pc_y[k] * ppg_1[k];

        t_4[k] = f_3 * pc_y[k] * ppg_2[k];

        t_5[k] = f_4 * ppf0_0[k]
                 - f_5 * ppf1_0[k]
                 + f_3 * pc_z[k] * ppg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_y, pc_z, ppf0_1, ppf0_2, ppf1_1, ppf1_2, \
                         ppg_3, ppg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * ppf0_1[k]
                 - f_7 * ppf1_1[k]
                 + f_3 * pc_y[k] * ppg_3[k];

        t_7[k] = f_3 * pc_z[k] * ppg_3[k];

        t_8[k] = f_3 * pc_y[k] * ppg_5[k];

        t_9[k] = f_6 * ppf0_2[k]
                 - f_7 * ppf1_2[k]
                 + f_3 * pc_z[k] * ppg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pc_x, pc_y, pc_z, spg_10, spg_12, psg_10, \
                         psg_12, ppg_6, ppg_9, ppg_10, ppg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * spg_10[k]
                  + f_0 * psg_10[k]
                  + f_3 * pc_x[k] * ppg_10[k];

        t_11[k] = f_3 * pc_z[k] * ppg_6[k];

        t_12[k] = f_0 * spg_12[k]
                  + f_0 * psg_12[k]
                  + f_3 * pc_x[k] * ppg_12[k];

        t_13[k] = f_3 * pc_y[k] * ppg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, spg_14, psg_14, ppf0_6, \
                         ppf0_8, ppf1_6, ppf1_8, ppg_10, ppg_12, \
                         ppg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * spg_14[k]
                  + f_0 * psg_14[k]
                  + f_3 * pc_x[k] * ppg_14[k];

        t_15[k] = f_1 * ppf0_6[k]
                  - f_2 * ppf1_6[k]
                  + f_3 * pc_y[k] * ppg_10[k];

        t_16[k] = f_3 * pc_z[k] * ppg_10[k];

        t_17[k] = f_6 * ppf0_8[k]
                  - f_7 * ppf1_8[k]
                  + f_3 * pc_y[k] * ppg_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pb_y, pc_y, pc_z, psh0_0, psg_0, \
                         psh1_0, ppf0_9, ppf1_9, ppg_13, ppg_14, \
                         ppg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_4 * ppf0_9[k]
                  - f_5 * ppf1_9[k]
                  + f_3 * pc_y[k] * ppg_13[k];

        t_19[k] = f_3 * pc_y[k] * ppg_14[k];

        t_20[k] = f_1 * ppf0_9[k]
                  - f_2 * ppf1_9[k]
                  + f_3 * pc_z[k] * ppg_14[k];

        t_21[k] = pb_y[k] * psh0_0[k]
                  - f_8 * pc_y[k] * psh1_0[k];

        t_22[k] = f_0 * psg_0[k]
                  + f_3 * pc_y[k] * ppg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pc_x, pc_y, pc_z, sph0_24, spg_18, sph1_24, \
                         psg_2, ppg_15, ppg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * pc_z[k] * ppg_15[k];

        t_24[k] = pa_x[k] * sph0_24[k]
                  + f_9 * spg_18[k]
                  - f_8 * pc_x[k] * sph1_24[k];

        t_25[k] = f_0 * psg_2[k]
                  + f_3 * pc_y[k] * ppg_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_x, pb_y, pc_x, pc_y, pc_z, sph0_27, spg_21, \
                         sph1_27, psh0_5, psh1_5, ppg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * psh0_5[k]
                  - f_8 * pc_y[k] * psh1_5[k];

        t_27[k] = pa_x[k] * sph0_27[k]
                  + f_10 * spg_21[k]
                  - f_8 * pc_x[k] * sph1_27[k];

        t_28[k] = f_3 * pc_z[k] * ppg_18[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pc_x, pc_y, pc_z, spg_25, psh0_9, \
                         psg_5, psh1_9, ppg_20, ppg_21, ppg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * psg_5[k]
                  + f_3 * pc_y[k] * ppg_20[k];

        t_30[k] = pb_y[k] * psh0_9[k]
                  - f_8 * pc_y[k] * psh1_9[k];

        t_31[k] = f_0 * spg_25[k]
                  + f_3 * pc_x[k] * ppg_25[k];

        t_32[k] = f_3 * pc_z[k] * ppg_21[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_x, pc_x, pc_y, sph0_36, spg_27, spg_29, \
                         sph1_36, psg_9, ppg_24, ppg_27, ppg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * spg_27[k]
                  + f_3 * pc_x[k] * ppg_27[k];

        t_34[k] = f_0 * psg_9[k]
                  + f_3 * pc_y[k] * ppg_24[k];

        t_35[k] = f_0 * spg_29[k]
                  + f_3 * pc_x[k] * ppg_29[k];

        t_36[k] = pa_x[k] * sph0_36[k]
                  - f_8 * pc_x[k] * sph1_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pc_x, pc_y, pc_z, sph0_38, sph0_39, \
                         sph1_38, sph1_39, psg_14, ppg_25, ppg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * pc_z[k] * ppg_25[k];

        t_38[k] = pa_x[k] * sph0_38[k]
                  - f_8 * pc_x[k] * sph1_38[k];

        t_39[k] = pa_x[k] * sph0_39[k]
                  - f_8 * pc_x[k] * sph1_39[k];

        t_40[k] = f_0 * psg_14[k]
                  + f_3 * pc_y[k] * ppg_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pb_z, pc_x, pc_y, pc_z, sph0_41, \
                         sph1_41, psh0_0, psg_0, psh1_0, ppg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_x[k] * sph0_41[k]
                  - f_8 * pc_x[k] * sph1_41[k];

        t_42[k] = pb_z[k] * psh0_0[k]
                  - f_8 * pc_z[k] * psh1_0[k];

        t_43[k] = f_3 * pc_y[k] * ppg_30[k];

        t_44[k] = f_0 * psg_0[k]
                  + f_3 * pc_z[k] * ppg_30[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_x, pb_z, pc_x, pc_y, pc_z, sph0_47, spg_35, \
                         sph1_47, psh0_3, psh1_3, ppg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_z[k] * psh0_3[k]
                  - f_8 * pc_z[k] * psh1_3[k];

        t_46[k] = f_3 * pc_y[k] * ppg_32[k];

        t_47[k] = pa_x[k] * sph0_47[k]
                  + f_9 * spg_35[k]
                  - f_8 * pc_x[k] * sph1_47[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_z, pc_y, pc_z, psh0_6, psg_3, psh1_6, ppg_33, \
                         ppg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_z[k] * psh0_6[k]
                  - f_8 * pc_z[k] * psh1_6[k];

        t_49[k] = f_0 * psg_3[k]
                  + f_3 * pc_z[k] * ppg_33[k];

        t_50[k] = f_3 * pc_y[k] * ppg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_x, pc_x, pc_z, sph0_51, spg_39, spg_40, \
                         spg_42, sph1_51, psg_6, ppg_36, ppg_40, \
                         ppg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_x[k] * sph0_51[k]
                  + f_10 * spg_39[k]
                  - f_8 * pc_x[k] * sph1_51[k];

        t_52[k] = f_0 * spg_40[k]
                  + f_3 * pc_x[k] * ppg_40[k];

        t_53[k] = f_0 * psg_6[k]
                  + f_3 * pc_z[k] * ppg_36[k];

        t_54[k] = f_0 * spg_42[k]
                  + f_3 * pc_x[k] * ppg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pc_x, pc_y, pc_z, sph0_57, spg_44, \
                         sph1_57, psg_10, ppg_39, ppg_40, ppg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_3 * pc_y[k] * ppg_39[k];

        t_56[k] = f_0 * spg_44[k]
                  + f_3 * pc_x[k] * ppg_44[k];

        t_57[k] = pa_x[k] * sph0_57[k]
                  - f_8 * pc_x[k] * sph1_57[k];

        t_58[k] = f_0 * psg_10[k]
                  + f_3 * pc_z[k] * ppg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_x, pc_x, pc_y, sph0_59, sph0_60, sph0_62, \
                         sph1_59, sph1_60, sph1_62, ppg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_x[k] * sph0_59[k]
                  - f_8 * pc_x[k] * sph1_59[k];

        t_60[k] = pa_x[k] * sph0_60[k]
                  - f_8 * pc_x[k] * sph1_60[k];

        t_61[k] = f_3 * pc_y[k] * ppg_44[k];

        t_62[k] = pa_x[k] * sph0_62[k]
                  - f_8 * pc_x[k] * sph1_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_y, pc_y, pc_z, sph0_0, sph0_1, sph0_3, \
                         spg_0, spg_1, sph1_0, sph1_1, sph1_3, ppg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * sph0_0[k]
                  - f_8 * pc_y[k] * sph1_0[k];

        t_64[k] = pa_y[k] * sph0_1[k]
                  + f_0 * spg_0[k]
                  - f_8 * pc_y[k] * sph1_1[k];

        t_65[k] = f_3 * pc_z[k] * ppg_45[k];

        t_66[k] = pa_y[k] * sph0_3[k]
                  + f_10 * spg_1[k]
                  - f_8 * pc_y[k] * sph1_3[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pc_y, pc_z, sph0_5, sph0_6, spg_3, \
                         sph1_5, sph1_6, ppg_46, ppg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_z[k] * ppg_46[k];

        t_68[k] = pa_y[k] * sph0_5[k]
                  - f_8 * pc_y[k] * sph1_5[k];

        t_69[k] = pa_y[k] * sph0_6[k]
                  + f_9 * spg_3[k]
                  - f_8 * pc_y[k] * sph1_6[k];

        t_70[k] = f_3 * pc_z[k] * ppg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pc_x, pc_y, sph0_8, sph0_9, spg_5, \
                         sph1_8, sph1_9, psg_25, psg_26, ppg_55, \
                         ppg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_y[k] * sph0_8[k]
                  + f_0 * spg_5[k]
                  - f_8 * pc_y[k] * sph1_8[k];

        t_72[k] = pa_y[k] * sph0_9[k]
                  - f_8 * pc_y[k] * sph1_9[k];

        t_73[k] = f_0 * psg_25[k]
                  + f_3 * pc_x[k] * ppg_55[k];

        t_74[k] = f_0 * psg_26[k]
                  + f_3 * pc_x[k] * ppg_56[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_x, pc_x, psh0_36, psg_27, psg_28, psg_29, \
                         psh1_36, ppg_57, ppg_58, ppg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * psg_27[k]
                  + f_3 * pc_x[k] * ppg_57[k];

        t_76[k] = f_0 * psg_28[k]
                  + f_3 * pc_x[k] * ppg_58[k];

        t_77[k] = f_0 * psg_29[k]
                  + f_3 * pc_x[k] * ppg_59[k];

        t_78[k] = pb_x[k] * psh0_36[k]
                  - f_8 * pc_x[k] * psh1_36[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_x, pc_x, pc_y, pc_z, spg_14, psh0_38, \
                         psh0_39, psh1_38, psh1_39, ppg_55, ppg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * ppg_55[k];

        t_80[k] = pb_x[k] * psh0_38[k]
                  - f_8 * pc_x[k] * psh1_38[k];

        t_81[k] = pb_x[k] * psh0_39[k]
                  - f_8 * pc_x[k] * psh1_39[k];

        t_82[k] = f_0 * spg_14[k]
                  + f_3 * pc_y[k] * ppg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_x, pc_y, pc_z, sph0_20, sph1_20, \
                         ppf0_40, ppf0_41, ppf1_40, ppf1_41, ppg_60, \
                         ppg_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pa_y[k] * sph0_20[k]
                  - f_8 * pc_y[k] * sph1_20[k];

        t_84[k] = f_1 * ppf0_40[k]
                  - f_2 * ppf1_40[k]
                  + f_3 * pc_x[k] * ppg_60[k];

        t_85[k] = f_11 * ppf0_41[k]
                  - f_12 * ppf1_41[k]
                  + f_3 * pc_x[k] * ppg_61[k];

        t_86[k] = f_3 * pc_z[k] * ppg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_z, ppf0_43, ppf0_45, ppf0_46, \
                         ppf1_43, ppf1_45, ppf1_46, ppg_61, ppg_63, ppg_65, \
                         ppg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_6 * ppf0_43[k]
                  - f_7 * ppf1_43[k]
                  + f_3 * pc_x[k] * ppg_63[k];

        t_88[k] = f_3 * pc_z[k] * ppg_61[k];

        t_89[k] = f_6 * ppf0_45[k]
                  - f_7 * ppf1_45[k]
                  + f_3 * pc_x[k] * ppg_65[k];

        t_90[k] = f_4 * ppf0_46[k]
                  - f_5 * ppf1_46[k]
                  + f_3 * pc_x[k] * ppg_66[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pc_x, pc_z, ppf0_48, ppf0_49, ppf1_48, \
                         ppf1_49, ppg_63, ppg_68, ppg_69, ppg_70, \
                         ppg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * ppg_63[k];

        t_92[k] = f_4 * ppf0_48[k]
                  - f_5 * ppf1_48[k]
                  + f_3 * pc_x[k] * ppg_68[k];

        t_93[k] = f_4 * ppf0_49[k]
                  - f_5 * ppf1_49[k]
                  + f_3 * pc_x[k] * ppg_69[k];

        t_94[k] = f_3 * pc_x[k] * ppg_70[k];

        t_95[k] = f_3 * pc_x[k] * ppg_71[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pc_x, pc_y, pc_z, spg_25, psg_25, \
                         ppf0_46, ppf1_46, ppg_70, ppg_72, ppg_73, \
                         ppg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_3 * pc_x[k] * ppg_72[k];

        t_97[k] = f_3 * pc_x[k] * ppg_73[k];

        t_98[k] = f_3 * pc_x[k] * ppg_74[k];

        t_99[k] = f_0 * spg_25[k]
                  + f_0 * psg_25[k]
                  + f_1 * ppf0_46[k]
                  - f_2 * ppf1_46[k]
                  + f_3 * pc_y[k] * ppg_70[k];

        t_100[k] = f_3 * pc_z[k] * ppg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, pc_z, spg_29, psg_29, ppf0_46, ppf0_47, \
                         ppf1_46, ppf1_47, ppg_71, ppg_72, ppg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_4 * ppf0_46[k]
                   - f_5 * ppf1_46[k]
                   + f_3 * pc_z[k] * ppg_71[k];

        t_102[k] = f_6 * ppf0_47[k]
                   - f_7 * ppf1_47[k]
                   + f_3 * pc_z[k] * ppg_72[k];

        t_103[k] = f_0 * spg_29[k]
                   + f_0 * psg_29[k]
                   + f_3 * pc_y[k] * ppg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pb_z, pc_y, pc_z, sph0_42, sph1_42, \
                         psh0_22, psh1_22, ppf0_49, ppf1_49, ppg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * ppf0_49[k]
                   - f_2 * ppf1_49[k]
                   + f_3 * pc_z[k] * ppg_74[k];

        t_105[k] = pa_y[k] * sph0_42[k]
                   - f_8 * pc_y[k] * sph1_42[k];

        t_106[k] = pb_z[k] * psh0_22[k]
                   - f_8 * pc_z[k] * psh1_22[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_y, pb_z, pc_y, pc_z, sph0_47, sph1_47, \
                         psh0_24, psg_15, psg_16, psh1_24, ppg_75, \
                         ppg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * psg_15[k]
                   + f_3 * pc_z[k] * ppg_75[k];

        t_108[k] = pb_z[k] * psh0_24[k]
                   - f_8 * pc_z[k] * psh1_24[k];

        t_109[k] = f_0 * psg_16[k]
                   + f_3 * pc_z[k] * ppg_76[k];

        t_110[k] = pa_y[k] * sph0_47[k]
                   - f_8 * pc_y[k] * sph1_47[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_z, pc_x, pc_z, psh0_27, psg_18, psh1_27, \
                         ppf0_58, ppf1_58, ppg_78, ppg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pb_z[k] * psh0_27[k]
                   - f_8 * pc_z[k] * psh1_27[k];

        t_112[k] = f_0 * psg_18[k]
                   + f_3 * pc_z[k] * ppg_78[k];

        t_113[k] = f_4 * ppf0_58[k]
                   - f_5 * ppf1_58[k]
                   + f_3 * pc_x[k] * ppg_83[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, pa_y, pc_x, pc_y, sph0_51, \
                         sph1_51, ppg_85, ppg_86, ppg_87, ppg_88, \
                         ppg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pa_y[k] * sph0_51[k]
                   - f_8 * pc_y[k] * sph1_51[k];

        t_115[k] = f_3 * pc_x[k] * ppg_85[k];

        t_116[k] = f_3 * pc_x[k] * ppg_86[k];

        t_117[k] = f_3 * pc_x[k] * ppg_87[k];

        t_118[k] = f_3 * pc_x[k] * ppg_88[k];

        t_119[k] = f_3 * pc_x[k] * ppg_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pa_y, pb_z, pc_y, pc_z, sph0_59, spg_42, \
                         sph1_59, psh0_36, psg_25, psh1_36, ppg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pb_z[k] * psh0_36[k]
                   - f_8 * pc_z[k] * psh1_36[k];

        t_121[k] = f_0 * psg_25[k]
                   + f_3 * pc_z[k] * ppg_85[k];

        t_122[k] = pa_y[k] * sph0_59[k]
                   + f_9 * spg_42[k]
                   - f_8 * pc_y[k] * sph1_59[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pa_y, pc_y, sph0_60, sph0_62, spg_43, spg_44, \
                         sph1_60, sph1_62, ppg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = pa_y[k] * sph0_60[k]
                   + f_10 * spg_43[k]
                   - f_8 * pc_y[k] * sph1_60[k];

        t_124[k] = f_0 * spg_44[k]
                   + f_3 * pc_y[k] * ppg_89[k];

        t_125[k] = pa_y[k] * sph0_62[k]
                   - f_8 * pc_y[k] * sph1_62[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_z, pc_y, pc_z, sph0_0, sph0_2, sph0_3, \
                         spg_0, sph1_0, sph1_2, sph1_3, ppg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_z[k] * sph0_0[k]
                   - f_8 * pc_z[k] * sph1_0[k];

        t_127[k] = f_3 * pc_y[k] * ppg_90[k];

        t_128[k] = pa_z[k] * sph0_2[k]
                   + f_0 * spg_0[k]
                   - f_8 * pc_z[k] * sph1_2[k];

        t_129[k] = pa_z[k] * sph0_3[k]
                   - f_8 * pc_z[k] * sph1_3[k];
    }
}

static auto
compute_prim_pph_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sph0, const size_t spg,
                                                          const size_t sph1, const size_t psh0,
                                                          const size_t psg, const size_t psh1,
                                                          const size_t ppf0, const size_t ppf1,
                                                          const size_t ppg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / gamma;
    const auto f_12 = 1.5 * p / (gamma * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sph0_5 = buffer.data(sph0 + 5);
    const auto *sph0_6 = buffer.data(sph0 + 6);
    const auto *sph0_7 = buffer.data(sph0 + 7);
    const auto *sph0_9 = buffer.data(sph0 + 9);
    const auto *sph0_15 = buffer.data(sph0 + 15);
    const auto *sph0_21 = buffer.data(sph0 + 21);
    const auto *sph0_24 = buffer.data(sph0 + 24);
    const auto *sph0_27 = buffer.data(sph0 + 27);
    const auto *sph0_36 = buffer.data(sph0 + 36);
    const auto *sph0_37 = buffer.data(sph0 + 37);
    const auto *sph0_38 = buffer.data(sph0 + 38);
    const auto *sph0_39 = buffer.data(sph0 + 39);

    const auto *spg_2 = buffer.data(spg + 2);
    const auto *spg_3 = buffer.data(spg + 3);
    const auto *spg_5 = buffer.data(spg + 5);
    const auto *spg_25 = buffer.data(spg + 25);
    const auto *spg_26 = buffer.data(spg + 26);
    const auto *spg_27 = buffer.data(spg + 27);
    const auto *spg_44 = buffer.data(spg + 44);

    const auto *sph1_5 = buffer.data(sph1 + 5);
    const auto *sph1_6 = buffer.data(sph1 + 6);
    const auto *sph1_7 = buffer.data(sph1 + 7);
    const auto *sph1_9 = buffer.data(sph1 + 9);
    const auto *sph1_15 = buffer.data(sph1 + 15);
    const auto *sph1_21 = buffer.data(sph1 + 21);
    const auto *sph1_24 = buffer.data(sph1 + 24);
    const auto *sph1_27 = buffer.data(sph1 + 27);
    const auto *sph1_36 = buffer.data(sph1 + 36);
    const auto *sph1_37 = buffer.data(sph1 + 37);
    const auto *sph1_38 = buffer.data(sph1 + 38);
    const auto *sph1_39 = buffer.data(sph1 + 39);

    const auto *psh0_44 = buffer.data(psh0 + 44);
    const auto *psh0_47 = buffer.data(psh0 + 47);
    const auto *psh0_51 = buffer.data(psh0 + 51);
    const auto *psh0_58 = buffer.data(psh0 + 58);
    const auto *psh0_59 = buffer.data(psh0 + 59);
    const auto *psh0_60 = buffer.data(psh0 + 60);
    const auto *psh0_62 = buffer.data(psh0 + 62);

    const auto *psg_30 = buffer.data(psg + 30);
    const auto *psg_32 = buffer.data(psg + 32);
    const auto *psg_35 = buffer.data(psg + 35);
    const auto *psg_40 = buffer.data(psg + 40);
    const auto *psg_41 = buffer.data(psg + 41);
    const auto *psg_42 = buffer.data(psg + 42);
    const auto *psg_43 = buffer.data(psg + 43);
    const auto *psg_44 = buffer.data(psg + 44);

    const auto *psh1_44 = buffer.data(psh1 + 44);
    const auto *psh1_47 = buffer.data(psh1 + 47);
    const auto *psh1_51 = buffer.data(psh1 + 51);
    const auto *psh1_58 = buffer.data(psh1 + 58);
    const auto *psh1_59 = buffer.data(psh1 + 59);
    const auto *psh1_60 = buffer.data(psh1 + 60);
    const auto *psh1_62 = buffer.data(psh1 + 62);

    const auto *ppf0_77 = buffer.data(ppf0 + 77);
    const auto *ppf0_80 = buffer.data(ppf0 + 80);
    const auto *ppf0_82 = buffer.data(ppf0 + 82);
    const auto *ppf0_83 = buffer.data(ppf0 + 83);
    const auto *ppf0_85 = buffer.data(ppf0 + 85);
    const auto *ppf0_86 = buffer.data(ppf0 + 86);
    const auto *ppf0_87 = buffer.data(ppf0 + 87);
    const auto *ppf0_88 = buffer.data(ppf0 + 88);
    const auto *ppf0_89 = buffer.data(ppf0 + 89);

    const auto *ppf1_77 = buffer.data(ppf1 + 77);
    const auto *ppf1_80 = buffer.data(ppf1 + 80);
    const auto *ppf1_82 = buffer.data(ppf1 + 82);
    const auto *ppf1_83 = buffer.data(ppf1 + 83);
    const auto *ppf1_85 = buffer.data(ppf1 + 85);
    const auto *ppf1_86 = buffer.data(ppf1 + 86);
    const auto *ppf1_87 = buffer.data(ppf1 + 87);
    const auto *ppf1_88 = buffer.data(ppf1 + 88);
    const auto *ppf1_89 = buffer.data(ppf1 + 89);

    const auto *ppg_92 = buffer.data(ppg + 92);
    const auto *ppg_95 = buffer.data(ppg + 95);
    const auto *ppg_100 = buffer.data(ppg + 100);
    const auto *ppg_101 = buffer.data(ppg + 101);
    const auto *ppg_102 = buffer.data(ppg + 102);
    const auto *ppg_103 = buffer.data(ppg + 103);
    const auto *ppg_104 = buffer.data(ppg + 104);
    const auto *ppg_105 = buffer.data(ppg + 105);
    const auto *ppg_107 = buffer.data(ppg + 107);
    const auto *ppg_110 = buffer.data(ppg + 110);
    const auto *ppg_112 = buffer.data(ppg + 112);
    const auto *ppg_115 = buffer.data(ppg + 115);
    const auto *ppg_116 = buffer.data(ppg + 116);
    const auto *ppg_117 = buffer.data(ppg + 117);
    const auto *ppg_118 = buffer.data(ppg + 118);
    const auto *ppg_119 = buffer.data(ppg + 119);
    const auto *ppg_120 = buffer.data(ppg + 120);
    const auto *ppg_122 = buffer.data(ppg + 122);
    const auto *ppg_123 = buffer.data(ppg + 123);
    const auto *ppg_125 = buffer.data(ppg + 125);
    const auto *ppg_126 = buffer.data(ppg + 126);
    const auto *ppg_127 = buffer.data(ppg + 127);
    const auto *ppg_129 = buffer.data(ppg + 129);
    const auto *ppg_130 = buffer.data(ppg + 130);
    const auto *ppg_131 = buffer.data(ppg + 131);
    const auto *ppg_132 = buffer.data(ppg + 132);
    const auto *ppg_133 = buffer.data(ppg + 133);
    const auto *ppg_134 = buffer.data(ppg + 134);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_z, pc_y, pc_z, sph0_5, sph0_6, sph0_7, \
                         spg_2, spg_3, sph1_5, sph1_6, sph1_7, ppg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_y[k] * ppg_92[k];

        t_131[k] = pa_z[k] * sph0_5[k]
                   + f_10 * spg_2[k]
                   - f_8 * pc_z[k] * sph1_5[k];

        t_132[k] = pa_z[k] * sph0_6[k]
                   - f_8 * pc_z[k] * sph1_6[k];

        t_133[k] = pa_z[k] * sph0_7[k]
                   + f_0 * spg_3[k]
                   - f_8 * pc_z[k] * sph1_7[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_z, pc_x, pc_y, pc_z, sph0_9, spg_5, \
                         sph1_9, psg_40, psg_41, ppg_95, ppg_100, \
                         ppg_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_3 * pc_y[k] * ppg_95[k];

        t_135[k] = pa_z[k] * sph0_9[k]
                   + f_9 * spg_5[k]
                   - f_8 * pc_z[k] * sph1_9[k];

        t_136[k] = f_0 * psg_40[k]
                   + f_3 * pc_x[k] * ppg_100[k];

        t_137[k] = f_0 * psg_41[k]
                   + f_3 * pc_x[k] * ppg_101[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_z, pc_x, pc_z, sph0_15, sph1_15, \
                         psg_42, psg_43, psg_44, ppg_102, ppg_103, \
                         ppg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * psg_42[k]
                   + f_3 * pc_x[k] * ppg_102[k];

        t_139[k] = f_0 * psg_43[k]
                   + f_3 * pc_x[k] * ppg_103[k];

        t_140[k] = f_0 * psg_44[k]
                   + f_3 * pc_x[k] * ppg_104[k];

        t_141[k] = pa_z[k] * sph0_15[k]
                   - f_8 * pc_z[k] * sph1_15[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_x, pc_x, pc_y, psh0_58, psh0_59, \
                         psh0_60, psh1_58, psh1_59, psh1_60, ppg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = pb_x[k] * psh0_58[k]
                   - f_8 * pc_x[k] * psh1_58[k];

        t_143[k] = pb_x[k] * psh0_59[k]
                   - f_8 * pc_x[k] * psh1_59[k];

        t_144[k] = pb_x[k] * psh0_60[k]
                   - f_8 * pc_x[k] * psh1_60[k];

        t_145[k] = f_3 * pc_y[k] * ppg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_z, pb_x, pc_x, pc_y, pc_z, sph0_21, sph1_21, \
                         psh0_62, psg_30, psh1_62, ppg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pb_x[k] * psh0_62[k]
                   - f_8 * pc_x[k] * psh1_62[k];

        t_147[k] = pa_z[k] * sph0_21[k]
                   - f_8 * pc_z[k] * sph1_21[k];

        t_148[k] = f_0 * psg_30[k]
                   + f_3 * pc_y[k] * ppg_105[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_z, pb_y, pc_y, pc_z, sph0_24, sph1_24, \
                         psh0_44, psh0_47, psg_32, psh1_44, psh1_47, \
                         ppg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pb_y[k] * psh0_44[k]
                   - f_8 * pc_y[k] * psh1_44[k];

        t_150[k] = pa_z[k] * sph0_24[k]
                   - f_8 * pc_z[k] * sph1_24[k];

        t_151[k] = f_0 * psg_32[k]
                   + f_3 * pc_y[k] * ppg_107[k];

        t_152[k] = pb_y[k] * psh0_47[k]
                   - f_8 * pc_y[k] * psh1_47[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_x, pc_y, pc_z, sph0_27, sph1_27, \
                         psg_35, ppf0_77, ppf1_77, ppg_110, ppg_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * sph0_27[k]
                   - f_8 * pc_z[k] * sph1_27[k];

        t_154[k] = f_4 * ppf0_77[k]
                   - f_5 * ppf1_77[k]
                   + f_3 * pc_x[k] * ppg_112[k];

        t_155[k] = f_0 * psg_35[k]
                   + f_3 * pc_y[k] * ppg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, pb_y, pc_x, pc_y, psh0_51, \
                         psh1_51, ppg_115, ppg_116, ppg_117, ppg_118, \
                         ppg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pb_y[k] * psh0_51[k]
                   - f_8 * pc_y[k] * psh1_51[k];

        t_157[k] = f_3 * pc_x[k] * ppg_115[k];

        t_158[k] = f_3 * pc_x[k] * ppg_116[k];

        t_159[k] = f_3 * pc_x[k] * ppg_117[k];

        t_160[k] = f_3 * pc_x[k] * ppg_118[k];

        t_161[k] = f_3 * pc_x[k] * ppg_119[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pa_z, pc_z, sph0_36, sph0_37, sph0_38, spg_25, \
                         spg_26, sph1_36, sph1_37, sph1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_z[k] * sph0_36[k]
                   - f_8 * pc_z[k] * sph1_36[k];

        t_163[k] = pa_z[k] * sph0_37[k]
                   + f_0 * spg_25[k]
                   - f_8 * pc_z[k] * sph1_37[k];

        t_164[k] = pa_z[k] * sph0_38[k]
                   + f_10 * spg_26[k]
                   - f_8 * pc_z[k] * sph1_38[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_z, pb_y, pc_y, pc_z, sph0_39, spg_27, \
                         sph1_39, psh0_62, psg_44, psh1_62, ppg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pa_z[k] * sph0_39[k]
                   + f_9 * spg_27[k]
                   - f_8 * pc_z[k] * sph1_39[k];

        t_166[k] = f_0 * psg_44[k]
                   + f_3 * pc_y[k] * ppg_119[k];

        t_167[k] = pb_y[k] * psh0_62[k]
                   - f_8 * pc_y[k] * psh1_62[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pc_x, pc_y, ppf0_80, ppf0_82, \
                         ppf0_83, ppf1_80, ppf1_82, ppf1_83, ppg_120, ppg_122, \
                         ppg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_1 * ppf0_80[k]
                   - f_2 * ppf1_80[k]
                   + f_3 * pc_x[k] * ppg_120[k];

        t_169[k] = f_3 * pc_y[k] * ppg_120[k];

        t_170[k] = f_11 * ppf0_82[k]
                   - f_12 * ppf1_82[k]
                   + f_3 * pc_x[k] * ppg_122[k];

        t_171[k] = f_6 * ppf0_83[k]
                   - f_7 * ppf1_83[k]
                   + f_3 * pc_x[k] * ppg_123[k];

        t_172[k] = f_3 * pc_y[k] * ppg_122[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pc_x, pc_y, ppf0_85, ppf0_86, ppf0_87, \
                         ppf1_85, ppf1_86, ppf1_87, ppg_125, ppg_126, \
                         ppg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_6 * ppf0_85[k]
                   - f_7 * ppf1_85[k]
                   + f_3 * pc_x[k] * ppg_125[k];

        t_174[k] = f_4 * ppf0_86[k]
                   - f_5 * ppf1_86[k]
                   + f_3 * pc_x[k] * ppg_126[k];

        t_175[k] = f_4 * ppf0_87[k]
                   - f_5 * ppf1_87[k]
                   + f_3 * pc_x[k] * ppg_127[k];

        t_176[k] = f_3 * pc_y[k] * ppg_125[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, t_181, t_182, pc_x, ppf0_89, ppf1_89, \
                         ppg_129, ppg_130, ppg_131, ppg_132, ppg_133, \
                         ppg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_4 * ppf0_89[k]
                   - f_5 * ppf1_89[k]
                   + f_3 * pc_x[k] * ppg_129[k];

        t_178[k] = f_3 * pc_x[k] * ppg_130[k];

        t_179[k] = f_3 * pc_x[k] * ppg_131[k];

        t_180[k] = f_3 * pc_x[k] * ppg_132[k];

        t_181[k] = f_3 * pc_x[k] * ppg_133[k];

        t_182[k] = f_3 * pc_x[k] * ppg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, ppf0_86, ppf0_87, ppf0_88, ppf1_86, \
                         ppf1_87, ppf1_88, ppg_130, ppg_131, ppg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_1 * ppf0_86[k]
                   - f_2 * ppf1_86[k]
                   + f_3 * pc_y[k] * ppg_130[k];

        t_184[k] = f_11 * ppf0_87[k]
                   - f_12 * ppf1_87[k]
                   + f_3 * pc_y[k] * ppg_131[k];

        t_185[k] = f_6 * ppf0_88[k]
                   - f_7 * ppf1_88[k]
                   + f_3 * pc_y[k] * ppg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pc_y, pc_z, spg_44, psg_44, ppf0_89, ppf1_89, \
                         ppg_133, ppg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_4 * ppf0_89[k]
                   - f_5 * ppf1_89[k]
                   + f_3 * pc_y[k] * ppg_133[k];

        t_187[k] = f_3 * pc_y[k] * ppg_134[k];

        t_188[k] = f_0 * spg_44[k]
                   + f_0 * psg_44[k]
                   + f_1 * ppf0_89[k]
                   - f_2 * ppf1_89[k]
                   + f_3 * pc_z[k] * ppg_134[k];
    }
}

auto
compute_prim_pph_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sph0,
                                                   const size_t spg, const size_t sph1,
                                                   const size_t psh0, const size_t psg,
                                                   const size_t psh1, const size_t ppf0,
                                                   const size_t ppf1, const size_t ppg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pph_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sph0,
                                                              spg, sph1, psh0, psg, psh1, ppf0,
                                                              ppf1, ppg, ncols, gamma, p, q);

    compute_prim_pph_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sph0,
                                                              spg, sph1, psh0, psg, psh1, ppf0,
                                                              ppf1, ppg, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
