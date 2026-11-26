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


#include "SimdThreeCenterElectronRepulsionVrrRecDPH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_dph_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t pph0, const size_t ppg,
                                                          const size_t pph1, const size_t dsh0,
                                                          const size_t dsg, const size_t dsh1,
                                                          const size_t dpf0, const size_t dpf1,
                                                          const size_t dpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *pph0_0 = buffer.data(pph0 + 0);
    const auto *pph0_5 = buffer.data(pph0 + 5);
    const auto *pph0_9 = buffer.data(pph0 + 9);
    const auto *pph0_14 = buffer.data(pph0 + 14);
    const auto *pph0_42 = buffer.data(pph0 + 42);
    const auto *pph0_47 = buffer.data(pph0 + 47);
    const auto *pph0_51 = buffer.data(pph0 + 51);
    const auto *pph0_87 = buffer.data(pph0 + 87);
    const auto *pph0_90 = buffer.data(pph0 + 90);
    const auto *pph0_99 = buffer.data(pph0 + 99);
    const auto *pph0_101 = buffer.data(pph0 + 101);
    const auto *pph0_102 = buffer.data(pph0 + 102);
    const auto *pph0_103 = buffer.data(pph0 + 103);
    const auto *pph0_104 = buffer.data(pph0 + 104);
    const auto *pph0_120 = buffer.data(pph0 + 120);
    const auto *pph0_122 = buffer.data(pph0 + 122);
    const auto *pph0_123 = buffer.data(pph0 + 123);
    const auto *pph0_125 = buffer.data(pph0 + 125);

    const auto *ppg_0 = buffer.data(ppg + 0);
    const auto *ppg_5 = buffer.data(ppg + 5);
    const auto *ppg_10 = buffer.data(ppg + 10);
    const auto *ppg_12 = buffer.data(ppg + 12);
    const auto *ppg_14 = buffer.data(ppg + 14);
    const auto *ppg_15 = buffer.data(ppg + 15);
    const auto *ppg_20 = buffer.data(ppg + 20);
    const auto *ppg_25 = buffer.data(ppg + 25);
    const auto *ppg_27 = buffer.data(ppg + 27);
    const auto *ppg_29 = buffer.data(ppg + 29);
    const auto *ppg_30 = buffer.data(ppg + 30);
    const auto *ppg_35 = buffer.data(ppg + 35);
    const auto *ppg_40 = buffer.data(ppg + 40);
    const auto *ppg_42 = buffer.data(ppg + 42);
    const auto *ppg_44 = buffer.data(ppg + 44);
    const auto *ppg_48 = buffer.data(ppg + 48);
    const auto *ppg_51 = buffer.data(ppg + 51);
    const auto *ppg_55 = buffer.data(ppg + 55);
    const auto *ppg_57 = buffer.data(ppg + 57);
    const auto *ppg_58 = buffer.data(ppg + 58);
    const auto *ppg_60 = buffer.data(ppg + 60);
    const auto *ppg_63 = buffer.data(ppg + 63);
    const auto *ppg_66 = buffer.data(ppg + 66);
    const auto *ppg_70 = buffer.data(ppg + 70);
    const auto *ppg_72 = buffer.data(ppg + 72);
    const auto *ppg_73 = buffer.data(ppg + 73);
    const auto *ppg_74 = buffer.data(ppg + 74);
    const auto *ppg_85 = buffer.data(ppg + 85);
    const auto *ppg_87 = buffer.data(ppg + 87);
    const auto *ppg_88 = buffer.data(ppg + 88);
    const auto *ppg_89 = buffer.data(ppg + 89);

    const auto *pph1_0 = buffer.data(pph1 + 0);
    const auto *pph1_5 = buffer.data(pph1 + 5);
    const auto *pph1_9 = buffer.data(pph1 + 9);
    const auto *pph1_14 = buffer.data(pph1 + 14);
    const auto *pph1_42 = buffer.data(pph1 + 42);
    const auto *pph1_47 = buffer.data(pph1 + 47);
    const auto *pph1_51 = buffer.data(pph1 + 51);
    const auto *pph1_87 = buffer.data(pph1 + 87);
    const auto *pph1_90 = buffer.data(pph1 + 90);
    const auto *pph1_99 = buffer.data(pph1 + 99);
    const auto *pph1_101 = buffer.data(pph1 + 101);
    const auto *pph1_102 = buffer.data(pph1 + 102);
    const auto *pph1_103 = buffer.data(pph1 + 103);
    const auto *pph1_104 = buffer.data(pph1 + 104);
    const auto *pph1_120 = buffer.data(pph1 + 120);
    const auto *pph1_122 = buffer.data(pph1 + 122);
    const auto *pph1_123 = buffer.data(pph1 + 123);
    const auto *pph1_125 = buffer.data(pph1 + 125);

    const auto *dsh0_0 = buffer.data(dsh0 + 0);
    const auto *dsh0_3 = buffer.data(dsh0 + 3);
    const auto *dsh0_5 = buffer.data(dsh0 + 5);
    const auto *dsh0_6 = buffer.data(dsh0 + 6);
    const auto *dsh0_9 = buffer.data(dsh0 + 9);
    const auto *dsh0_15 = buffer.data(dsh0 + 15);
    const auto *dsh0_17 = buffer.data(dsh0 + 17);
    const auto *dsh0_18 = buffer.data(dsh0 + 18);
    const auto *dsh0_20 = buffer.data(dsh0 + 20);
    const auto *dsh0_24 = buffer.data(dsh0 + 24);
    const auto *dsh0_27 = buffer.data(dsh0 + 27);

    const auto *dsg_0 = buffer.data(dsg + 0);
    const auto *dsg_1 = buffer.data(dsg + 1);
    const auto *dsg_2 = buffer.data(dsg + 2);
    const auto *dsg_3 = buffer.data(dsg + 3);
    const auto *dsg_5 = buffer.data(dsg + 5);
    const auto *dsg_6 = buffer.data(dsg + 6);
    const auto *dsg_9 = buffer.data(dsg + 9);
    const auto *dsg_10 = buffer.data(dsg + 10);
    const auto *dsg_12 = buffer.data(dsg + 12);
    const auto *dsg_13 = buffer.data(dsg + 13);
    const auto *dsg_14 = buffer.data(dsg + 14);
    const auto *dsg_15 = buffer.data(dsg + 15);
    const auto *dsg_16 = buffer.data(dsg + 16);
    const auto *dsg_18 = buffer.data(dsg + 18);
    const auto *dsg_20 = buffer.data(dsg + 20);
    const auto *dsg_21 = buffer.data(dsg + 21);
    const auto *dsg_25 = buffer.data(dsg + 25);
    const auto *dsg_27 = buffer.data(dsg + 27);
    const auto *dsg_28 = buffer.data(dsg + 28);

    const auto *dsh1_0 = buffer.data(dsh1 + 0);
    const auto *dsh1_3 = buffer.data(dsh1 + 3);
    const auto *dsh1_5 = buffer.data(dsh1 + 5);
    const auto *dsh1_6 = buffer.data(dsh1 + 6);
    const auto *dsh1_9 = buffer.data(dsh1 + 9);
    const auto *dsh1_15 = buffer.data(dsh1 + 15);
    const auto *dsh1_17 = buffer.data(dsh1 + 17);
    const auto *dsh1_18 = buffer.data(dsh1 + 18);
    const auto *dsh1_20 = buffer.data(dsh1 + 20);
    const auto *dsh1_24 = buffer.data(dsh1 + 24);
    const auto *dsh1_27 = buffer.data(dsh1 + 27);

    const auto *dpf0_0 = buffer.data(dpf0 + 0);
    const auto *dpf0_1 = buffer.data(dpf0 + 1);
    const auto *dpf0_2 = buffer.data(dpf0 + 2);
    const auto *dpf0_6 = buffer.data(dpf0 + 6);
    const auto *dpf0_8 = buffer.data(dpf0 + 8);
    const auto *dpf0_9 = buffer.data(dpf0 + 9);
    const auto *dpf0_28 = buffer.data(dpf0 + 28);
    const auto *dpf0_29 = buffer.data(dpf0 + 29);
    const auto *dpf0_33 = buffer.data(dpf0 + 33);
    const auto *dpf0_36 = buffer.data(dpf0 + 36);
    const auto *dpf0_37 = buffer.data(dpf0 + 37);
    const auto *dpf0_39 = buffer.data(dpf0 + 39);
    const auto *dpf0_40 = buffer.data(dpf0 + 40);
    const auto *dpf0_42 = buffer.data(dpf0 + 42);

    const auto *dpf1_0 = buffer.data(dpf1 + 0);
    const auto *dpf1_1 = buffer.data(dpf1 + 1);
    const auto *dpf1_2 = buffer.data(dpf1 + 2);
    const auto *dpf1_6 = buffer.data(dpf1 + 6);
    const auto *dpf1_8 = buffer.data(dpf1 + 8);
    const auto *dpf1_9 = buffer.data(dpf1 + 9);
    const auto *dpf1_28 = buffer.data(dpf1 + 28);
    const auto *dpf1_29 = buffer.data(dpf1 + 29);
    const auto *dpf1_33 = buffer.data(dpf1 + 33);
    const auto *dpf1_36 = buffer.data(dpf1 + 36);
    const auto *dpf1_37 = buffer.data(dpf1 + 37);
    const auto *dpf1_39 = buffer.data(dpf1 + 39);
    const auto *dpf1_40 = buffer.data(dpf1 + 40);
    const auto *dpf1_42 = buffer.data(dpf1 + 42);

    const auto *dpg_0 = buffer.data(dpg + 0);
    const auto *dpg_1 = buffer.data(dpg + 1);
    const auto *dpg_2 = buffer.data(dpg + 2);
    const auto *dpg_3 = buffer.data(dpg + 3);
    const auto *dpg_5 = buffer.data(dpg + 5);
    const auto *dpg_6 = buffer.data(dpg + 6);
    const auto *dpg_9 = buffer.data(dpg + 9);
    const auto *dpg_10 = buffer.data(dpg + 10);
    const auto *dpg_12 = buffer.data(dpg + 12);
    const auto *dpg_13 = buffer.data(dpg + 13);
    const auto *dpg_14 = buffer.data(dpg + 14);
    const auto *dpg_15 = buffer.data(dpg + 15);
    const auto *dpg_17 = buffer.data(dpg + 17);
    const auto *dpg_18 = buffer.data(dpg + 18);
    const auto *dpg_20 = buffer.data(dpg + 20);
    const auto *dpg_21 = buffer.data(dpg + 21);
    const auto *dpg_24 = buffer.data(dpg + 24);
    const auto *dpg_25 = buffer.data(dpg + 25);
    const auto *dpg_27 = buffer.data(dpg + 27);
    const auto *dpg_29 = buffer.data(dpg + 29);
    const auto *dpg_30 = buffer.data(dpg + 30);
    const auto *dpg_32 = buffer.data(dpg + 32);
    const auto *dpg_33 = buffer.data(dpg + 33);
    const auto *dpg_35 = buffer.data(dpg + 35);
    const auto *dpg_36 = buffer.data(dpg + 36);
    const auto *dpg_39 = buffer.data(dpg + 39);
    const auto *dpg_40 = buffer.data(dpg + 40);
    const auto *dpg_42 = buffer.data(dpg + 42);
    const auto *dpg_43 = buffer.data(dpg + 43);
    const auto *dpg_44 = buffer.data(dpg + 44);
    const auto *dpg_45 = buffer.data(dpg + 45);
    const auto *dpg_46 = buffer.data(dpg + 46);
    const auto *dpg_48 = buffer.data(dpg + 48);
    const auto *dpg_50 = buffer.data(dpg + 50);
    const auto *dpg_51 = buffer.data(dpg + 51);
    const auto *dpg_55 = buffer.data(dpg + 55);
    const auto *dpg_56 = buffer.data(dpg + 56);
    const auto *dpg_57 = buffer.data(dpg + 57);
    const auto *dpg_58 = buffer.data(dpg + 58);
    const auto *dpg_59 = buffer.data(dpg + 59);
    const auto *dpg_60 = buffer.data(dpg + 60);
    const auto *dpg_61 = buffer.data(dpg + 61);
    const auto *dpg_62 = buffer.data(dpg + 62);
    const auto *dpg_63 = buffer.data(dpg + 63);
    const auto *dpg_65 = buffer.data(dpg + 65);
    const auto *dpg_66 = buffer.data(dpg + 66);
    const auto *dpg_70 = buffer.data(dpg + 70);
    const auto *dpg_72 = buffer.data(dpg + 72);
    const auto *dpg_73 = buffer.data(dpg + 73);
    const auto *dpg_74 = buffer.data(dpg + 74);
    const auto *dpg_75 = buffer.data(dpg + 75);
    const auto *dpg_76 = buffer.data(dpg + 76);
    const auto *dpg_78 = buffer.data(dpg + 78);
    const auto *dpg_80 = buffer.data(dpg + 80);
    const auto *dpg_81 = buffer.data(dpg + 81);
    const auto *dpg_85 = buffer.data(dpg + 85);
    const auto *dpg_87 = buffer.data(dpg + 87);
    const auto *dpg_88 = buffer.data(dpg + 88);
    const auto *dpg_89 = buffer.data(dpg + 89);
    const auto *dpg_90 = buffer.data(dpg + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ppg_0, dsg_0, dpf0_0, \
                         dpf1_0, dpg_0, dpg_1, dpg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ppg_0[k]
                 + f_1 * dsg_0[k]
                 + f_2 * dpf0_0[k]
                 - f_3 * dpf1_0[k]
                 + f_4 * pc_x[k] * dpg_0[k];

        t_1[k] = f_4 * pc_y[k] * dpg_0[k];

        t_2[k] = f_4 * pc_z[k] * dpg_0[k];

        t_3[k] = f_5 * dpf0_0[k]
                 - f_6 * dpf1_0[k]
                 + f_4 * pc_y[k] * dpg_1[k];

        t_4[k] = f_4 * pc_y[k] * dpg_2[k];

        t_5[k] = f_5 * dpf0_0[k]
                 - f_6 * dpf1_0[k]
                 + f_4 * pc_z[k] * dpg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_y, pc_z, dpf0_1, dpf0_2, dpf1_1, dpf1_2, \
                         dpg_3, dpg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * dpf0_1[k]
                 - f_8 * dpf1_1[k]
                 + f_4 * pc_y[k] * dpg_3[k];

        t_7[k] = f_4 * pc_z[k] * dpg_3[k];

        t_8[k] = f_4 * pc_y[k] * dpg_5[k];

        t_9[k] = f_7 * dpf0_2[k]
                 - f_8 * dpf1_2[k]
                 + f_4 * pc_z[k] * dpg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pc_x, pc_y, pc_z, ppg_10, ppg_12, dsg_10, \
                         dsg_12, dpg_6, dpg_9, dpg_10, dpg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * ppg_10[k]
                  + f_1 * dsg_10[k]
                  + f_4 * pc_x[k] * dpg_10[k];

        t_11[k] = f_4 * pc_z[k] * dpg_6[k];

        t_12[k] = f_0 * ppg_12[k]
                  + f_1 * dsg_12[k]
                  + f_4 * pc_x[k] * dpg_12[k];

        t_13[k] = f_4 * pc_y[k] * dpg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, ppg_14, dsg_14, dpf0_6, \
                         dpf0_8, dpf1_6, dpf1_8, dpg_10, dpg_12, \
                         dpg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * ppg_14[k]
                  + f_1 * dsg_14[k]
                  + f_4 * pc_x[k] * dpg_14[k];

        t_15[k] = f_2 * dpf0_6[k]
                  - f_3 * dpf1_6[k]
                  + f_4 * pc_y[k] * dpg_10[k];

        t_16[k] = f_4 * pc_z[k] * dpg_10[k];

        t_17[k] = f_7 * dpf0_8[k]
                  - f_8 * dpf1_8[k]
                  + f_4 * pc_y[k] * dpg_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pb_y, pc_y, pc_z, dsh0_0, dsg_0, \
                         dsh1_0, dpf0_9, dpf1_9, dpg_13, dpg_14, \
                         dpg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * dpf0_9[k]
                  - f_6 * dpf1_9[k]
                  + f_4 * pc_y[k] * dpg_13[k];

        t_19[k] = f_4 * pc_y[k] * dpg_14[k];

        t_20[k] = f_2 * dpf0_9[k]
                  - f_3 * dpf1_9[k]
                  + f_4 * pc_z[k] * dpg_14[k];

        t_21[k] = pb_y[k] * dsh0_0[k]
                  - f_9 * pc_y[k] * dsh1_0[k];

        t_22[k] = f_1 * dsg_0[k]
                  + f_4 * pc_y[k] * dpg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_y, pc_z, dsh0_3, dsh0_5, dsg_1, \
                         dsg_2, dsh1_3, dsh1_5, dpg_15, dpg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * pc_z[k] * dpg_15[k];

        t_24[k] = pb_y[k] * dsh0_3[k]
                  + f_0 * dsg_1[k]
                  - f_9 * pc_y[k] * dsh1_3[k];

        t_25[k] = f_1 * dsg_2[k]
                  + f_4 * pc_y[k] * dpg_17[k];

        t_26[k] = pb_y[k] * dsh0_5[k]
                  - f_9 * pc_y[k] * dsh1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_y, pc_y, pc_z, dsh0_6, dsh0_9, dsg_3, \
                         dsg_5, dsh1_6, dsh1_9, dpg_18, dpg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * dsh0_6[k]
                  + f_10 * dsg_3[k]
                  - f_9 * pc_y[k] * dsh1_6[k];

        t_28[k] = f_4 * pc_z[k] * dpg_18[k];

        t_29[k] = f_1 * dsg_5[k]
                  + f_4 * pc_y[k] * dpg_20[k];

        t_30[k] = pb_y[k] * dsh0_9[k]
                  - f_9 * pc_y[k] * dsh1_9[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pc_x, pc_y, pc_z, ppg_25, ppg_27, dsg_9, \
                         dpg_21, dpg_24, dpg_25, dpg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * ppg_25[k]
                  + f_4 * pc_x[k] * dpg_25[k];

        t_32[k] = f_4 * pc_z[k] * dpg_21[k];

        t_33[k] = f_0 * ppg_27[k]
                  + f_4 * pc_x[k] * dpg_27[k];

        t_34[k] = f_1 * dsg_9[k]
                  + f_4 * pc_y[k] * dpg_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pc_x, pc_y, pc_z, ppg_29, dsh0_15, dsg_10, \
                         dsh1_15, dpg_25, dpg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * ppg_29[k]
                  + f_4 * pc_x[k] * dpg_29[k];

        t_36[k] = pb_y[k] * dsh0_15[k]
                  + f_11 * dsg_10[k]
                  - f_9 * pc_y[k] * dsh1_15[k];

        t_37[k] = f_4 * pc_z[k] * dpg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, dsh0_17, dsh0_18, dsh0_20, \
                         dsg_12, dsg_13, dsg_14, dsh1_17, dsh1_18, dsh1_20, \
                         dpg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * dsh0_17[k]
                  + f_10 * dsg_12[k]
                  - f_9 * pc_y[k] * dsh1_17[k];

        t_39[k] = pb_y[k] * dsh0_18[k]
                  + f_0 * dsg_13[k]
                  - f_9 * pc_y[k] * dsh1_18[k];

        t_40[k] = f_1 * dsg_14[k]
                  + f_4 * pc_y[k] * dpg_29[k];

        t_41[k] = pb_y[k] * dsh0_20[k]
                  - f_9 * pc_y[k] * dsh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, dsh0_0, dsh0_3, \
                         dsg_0, dsh1_0, dsh1_3, dpg_30, dpg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * dsh0_0[k]
                  - f_9 * pc_z[k] * dsh1_0[k];

        t_43[k] = f_4 * pc_y[k] * dpg_30[k];

        t_44[k] = f_1 * dsg_0[k]
                  + f_4 * pc_z[k] * dpg_30[k];

        t_45[k] = pb_z[k] * dsh0_3[k]
                  - f_9 * pc_z[k] * dsh1_3[k];

        t_46[k] = f_4 * pc_y[k] * dpg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, dsh0_5, dsh0_6, dsg_2, \
                         dsg_3, dsh1_5, dsh1_6, dpg_33, dpg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * dsh0_5[k]
                  + f_0 * dsg_2[k]
                  - f_9 * pc_z[k] * dsh1_5[k];

        t_48[k] = pb_z[k] * dsh0_6[k]
                  - f_9 * pc_z[k] * dsh1_6[k];

        t_49[k] = f_1 * dsg_3[k]
                  + f_4 * pc_z[k] * dpg_33[k];

        t_50[k] = f_4 * pc_y[k] * dpg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, ppg_40, ppg_42, dsh0_9, \
                         dsg_5, dsg_6, dsh1_9, dpg_36, dpg_40, dpg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * dsh0_9[k]
                  + f_10 * dsg_5[k]
                  - f_9 * pc_z[k] * dsh1_9[k];

        t_52[k] = f_0 * ppg_40[k]
                  + f_4 * pc_x[k] * dpg_40[k];

        t_53[k] = f_1 * dsg_6[k]
                  + f_4 * pc_z[k] * dpg_36[k];

        t_54[k] = f_0 * ppg_42[k]
                  + f_4 * pc_x[k] * dpg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_y, pc_z, ppg_44, dsh0_15, \
                         dsg_10, dsh1_15, dpg_39, dpg_40, dpg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_4 * pc_y[k] * dpg_39[k];

        t_56[k] = f_0 * ppg_44[k]
                  + f_4 * pc_x[k] * dpg_44[k];

        t_57[k] = pb_z[k] * dsh0_15[k]
                  - f_9 * pc_z[k] * dsh1_15[k];

        t_58[k] = f_1 * dsg_10[k]
                  + f_4 * pc_z[k] * dpg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pc_y, dpf0_28, dpf0_29, dpf1_28, dpf1_29, dpg_42, \
                         dpg_43, dpg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_7 * dpf0_28[k]
                  - f_8 * dpf1_28[k]
                  + f_4 * pc_y[k] * dpg_42[k];

        t_60[k] = f_5 * dpf0_29[k]
                  - f_6 * dpf1_29[k]
                  + f_4 * pc_y[k] * dpg_43[k];

        t_61[k] = f_4 * pc_y[k] * dpg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pb_z, pc_y, pc_z, pph0_0, ppg_0, \
                         pph1_0, dsh0_20, dsg_14, dsh1_20, dpg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_z[k] * dsh0_20[k]
                  + f_11 * dsg_14[k]
                  - f_9 * pc_z[k] * dsh1_20[k];

        t_63[k] = pa_y[k] * pph0_0[k]
                  - f_9 * pc_y[k] * pph1_0[k];

        t_64[k] = f_1 * ppg_0[k]
                  + f_4 * pc_y[k] * dpg_45[k];

        t_65[k] = f_4 * pc_z[k] * dpg_45[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_y, pc_x, pc_y, pc_z, pph0_5, ppg_48, pph1_5, \
                         dsg_18, dpf0_33, dpf1_33, dpg_46, dpg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_1 * ppg_48[k]
                  + f_1 * dsg_18[k]
                  + f_7 * dpf0_33[k]
                  - f_8 * dpf1_33[k]
                  + f_4 * pc_x[k] * dpg_48[k];

        t_67[k] = f_4 * pc_z[k] * dpg_46[k];

        t_68[k] = pa_y[k] * pph0_5[k]
                  - f_9 * pc_y[k] * pph1_5[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pc_x, pc_y, pc_z, ppg_5, ppg_51, dsg_21, dpf0_36, \
                         dpf1_36, dpg_48, dpg_50, dpg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * ppg_51[k]
                  + f_1 * dsg_21[k]
                  + f_5 * dpf0_36[k]
                  - f_6 * dpf1_36[k]
                  + f_4 * pc_x[k] * dpg_51[k];

        t_70[k] = f_4 * pc_z[k] * dpg_48[k];

        t_71[k] = f_1 * ppg_5[k]
                  + f_4 * pc_y[k] * dpg_50[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_y, pc_x, pc_y, pc_z, pph0_9, ppg_55, pph1_9, \
                         dsg_25, dpg_51, dpg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * pph0_9[k]
                  - f_9 * pc_y[k] * pph1_9[k];

        t_73[k] = f_1 * ppg_55[k]
                  + f_1 * dsg_25[k]
                  + f_4 * pc_x[k] * dpg_55[k];

        t_74[k] = f_4 * pc_z[k] * dpg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_y, pc_x, pc_y, pph0_14, ppg_57, ppg_58, pph1_14, \
                         dsg_27, dsg_28, dpg_57, dpg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * ppg_57[k]
                  + f_1 * dsg_27[k]
                  + f_4 * pc_x[k] * dpg_57[k];

        t_76[k] = f_1 * ppg_58[k]
                  + f_1 * dsg_28[k]
                  + f_4 * pc_x[k] * dpg_58[k];

        t_77[k] = pa_y[k] * pph0_14[k]
                  - f_9 * pc_y[k] * pph1_14[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pc_y, pc_z, ppg_10, dpf0_36, dpf0_37, \
                         dpf1_36, dpf1_37, dpg_55, dpg_56, dpg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * ppg_10[k]
                  + f_2 * dpf0_36[k]
                  - f_3 * dpf1_36[k]
                  + f_4 * pc_y[k] * dpg_55[k];

        t_79[k] = f_4 * pc_z[k] * dpg_55[k];

        t_80[k] = f_5 * dpf0_36[k]
                  - f_6 * dpf1_36[k]
                  + f_4 * pc_z[k] * dpg_56[k];

        t_81[k] = f_7 * dpf0_37[k]
                  - f_8 * dpf1_37[k]
                  + f_4 * pc_z[k] * dpg_57[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pc_x, pc_y, pc_z, ppg_14, ppg_60, dpf0_39, dpf0_40, \
                         dpf1_39, dpf1_40, dpg_59, dpg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_1 * ppg_14[k]
                  + f_4 * pc_y[k] * dpg_59[k];

        t_83[k] = f_2 * dpf0_39[k]
                  - f_3 * dpf1_39[k]
                  + f_4 * pc_z[k] * dpg_59[k];

        t_84[k] = f_1 * ppg_60[k]
                  + f_2 * dpf0_40[k]
                  - f_3 * dpf1_40[k]
                  + f_4 * pc_x[k] * dpg_60[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_x, pc_x, pc_y, pc_z, pph0_87, ppg_15, \
                         ppg_63, pph1_87, dsg_15, dpg_60, dpg_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_1 * ppg_15[k]
                  + f_1 * dsg_15[k]
                  + f_4 * pc_y[k] * dpg_60[k];

        t_86[k] = f_4 * pc_z[k] * dpg_60[k];

        t_87[k] = pa_x[k] * pph0_87[k]
                  + f_10 * ppg_63[k]
                  - f_9 * pc_x[k] * pph1_87[k];

        t_88[k] = f_4 * pc_z[k] * dpg_61[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pc_x, pc_z, pph0_90, ppg_66, pph1_90, \
                         dpf0_40, dpf1_40, dpg_62, dpg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_5 * dpf0_40[k]
                  - f_6 * dpf1_40[k]
                  + f_4 * pc_z[k] * dpg_62[k];

        t_90[k] = pa_x[k] * pph0_90[k]
                  + f_0 * ppg_66[k]
                  - f_9 * pc_x[k] * pph1_90[k];

        t_91[k] = f_4 * pc_z[k] * dpg_63[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, ppg_20, ppg_70, dsg_20, \
                         dpf0_42, dpf1_42, dpg_65, dpg_66, dpg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_1 * ppg_20[k]
                  + f_1 * dsg_20[k]
                  + f_4 * pc_y[k] * dpg_65[k];

        t_93[k] = f_7 * dpf0_42[k]
                  - f_8 * dpf1_42[k]
                  + f_4 * pc_z[k] * dpg_65[k];

        t_94[k] = f_1 * ppg_70[k]
                  + f_4 * pc_x[k] * dpg_70[k];

        t_95[k] = f_4 * pc_z[k] * dpg_66[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pc_x, pph0_99, ppg_72, ppg_73, ppg_74, \
                         pph1_99, dpg_72, dpg_73, dpg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * ppg_72[k]
                  + f_4 * pc_x[k] * dpg_72[k];

        t_97[k] = f_1 * ppg_73[k]
                  + f_4 * pc_x[k] * dpg_73[k];

        t_98[k] = f_1 * ppg_74[k]
                  + f_4 * pc_x[k] * dpg_74[k];

        t_99[k] = pa_x[k] * pph0_99[k]
                  - f_9 * pc_x[k] * pph1_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pc_x, pc_z, pph0_101, pph0_102, \
                         pph0_103, pph1_101, pph1_102, pph1_103, \
                         dpg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_4 * pc_z[k] * dpg_70[k];

        t_101[k] = pa_x[k] * pph0_101[k]
                   - f_9 * pc_x[k] * pph1_101[k];

        t_102[k] = pa_x[k] * pph0_102[k]
                   - f_9 * pc_x[k] * pph1_102[k];

        t_103[k] = pa_x[k] * pph0_103[k]
                   - f_9 * pc_x[k] * pph1_103[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pa_y, pc_x, pc_y, pc_z, pph0_42, \
                         pph0_104, ppg_30, pph1_42, pph1_104, dsg_15, \
                         dpg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_x[k] * pph0_104[k]
                   - f_9 * pc_x[k] * pph1_104[k];

        t_105[k] = pa_y[k] * pph0_42[k]
                   - f_9 * pc_y[k] * pph1_42[k];

        t_106[k] = f_1 * ppg_30[k]
                   + f_4 * pc_y[k] * dpg_75[k];

        t_107[k] = f_1 * dsg_15[k]
                   + f_4 * pc_z[k] * dpg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pa_y, pb_z, pc_y, pc_z, pph0_47, pph1_47, \
                         dsh0_24, dsh0_27, dsg_16, dsh1_24, dsh1_27, \
                         dpg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * dsh0_24[k]
                   - f_9 * pc_z[k] * dsh1_24[k];

        t_109[k] = f_1 * dsg_16[k]
                   + f_4 * pc_z[k] * dpg_76[k];

        t_110[k] = pa_y[k] * pph0_47[k]
                   - f_9 * pc_y[k] * pph1_47[k];

        t_111[k] = pb_z[k] * dsh0_27[k]
                   - f_9 * pc_z[k] * dsh1_27[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_y, pc_x, pc_y, pc_z, pph0_51, ppg_35, \
                         ppg_85, pph1_51, dsg_18, dpg_78, dpg_80, \
                         dpg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_1 * dsg_18[k]
                   + f_4 * pc_z[k] * dpg_78[k];

        t_113[k] = f_1 * ppg_35[k]
                   + f_4 * pc_y[k] * dpg_80[k];

        t_114[k] = pa_y[k] * pph0_51[k]
                   - f_9 * pc_y[k] * pph1_51[k];

        t_115[k] = f_1 * ppg_85[k]
                   + f_4 * pc_x[k] * dpg_85[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pc_x, pc_z, ppg_87, ppg_88, ppg_89, \
                         dsg_21, dpg_81, dpg_87, dpg_88, dpg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_1 * dsg_21[k]
                   + f_4 * pc_z[k] * dpg_81[k];

        t_117[k] = f_1 * ppg_87[k]
                   + f_4 * pc_x[k] * dpg_87[k];

        t_118[k] = f_1 * ppg_88[k]
                   + f_4 * pc_x[k] * dpg_88[k];

        t_119[k] = f_1 * ppg_89[k]
                   + f_4 * pc_x[k] * dpg_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_x, pc_x, pc_z, pph0_120, pph0_122, \
                         pph0_123, pph1_120, pph1_122, pph1_123, dsg_25, \
                         dpg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_x[k] * pph0_120[k]
                   - f_9 * pc_x[k] * pph1_120[k];

        t_121[k] = f_1 * dsg_25[k]
                   + f_4 * pc_z[k] * dpg_85[k];

        t_122[k] = pa_x[k] * pph0_122[k]
                   - f_9 * pc_x[k] * pph1_122[k];

        t_123[k] = pa_x[k] * pph0_123[k]
                   - f_9 * pc_x[k] * pph1_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_x, pa_z, pc_x, pc_y, pc_z, pph0_0, \
                         pph0_125, ppg_44, pph1_0, pph1_125, dpg_89, \
                         dpg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * ppg_44[k]
                   + f_4 * pc_y[k] * dpg_89[k];

        t_125[k] = pa_x[k] * pph0_125[k]
                   - f_9 * pc_x[k] * pph1_125[k];

        t_126[k] = pa_z[k] * pph0_0[k]
                   - f_9 * pc_z[k] * pph1_0[k];

        t_127[k] = f_4 * pc_y[k] * dpg_90[k];
    }
}

static auto
compute_prim_dph_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t pph0, const size_t ppg,
                                                          const size_t pph1, const size_t dsh0,
                                                          const size_t dsg, const size_t dsh1,
                                                          const size_t dpf0, const size_t dpf1,
                                                          const size_t dpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);
    const auto f_14 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *pph0_3 = buffer.data(pph0 + 3);
    const auto *pph0_6 = buffer.data(pph0 + 6);
    const auto *pph0_10 = buffer.data(pph0 + 10);
    const auto *pph0_21 = buffer.data(pph0 + 21);
    const auto *pph0_24 = buffer.data(pph0 + 24);
    const auto *pph0_27 = buffer.data(pph0 + 27);
    const auto *pph0_64 = buffer.data(pph0 + 64);
    const auto *pph0_66 = buffer.data(pph0 + 66);
    const auto *pph0_126 = buffer.data(pph0 + 126);
    const auto *pph0_128 = buffer.data(pph0 + 128);
    const auto *pph0_154 = buffer.data(pph0 + 154);
    const auto *pph0_162 = buffer.data(pph0 + 162);
    const auto *pph0_163 = buffer.data(pph0 + 163);
    const auto *pph0_164 = buffer.data(pph0 + 164);
    const auto *pph0_165 = buffer.data(pph0 + 165);
    const auto *pph0_167 = buffer.data(pph0 + 167);
    const auto *pph0_173 = buffer.data(pph0 + 173);
    const auto *pph0_177 = buffer.data(pph0 + 177);
    const auto *pph0_183 = buffer.data(pph0 + 183);
    const auto *pph0_184 = buffer.data(pph0 + 184);
    const auto *pph0_185 = buffer.data(pph0 + 185);
    const auto *pph0_186 = buffer.data(pph0 + 186);
    const auto *pph0_188 = buffer.data(pph0 + 188);

    const auto *ppg_0 = buffer.data(ppg + 0);
    const auto *ppg_14 = buffer.data(ppg + 14);
    const auto *ppg_15 = buffer.data(ppg + 15);
    const auto *ppg_30 = buffer.data(ppg + 30);
    const auto *ppg_59 = buffer.data(ppg + 59);
    const auto *ppg_70 = buffer.data(ppg + 70);
    const auto *ppg_74 = buffer.data(ppg + 74);
    const auto *ppg_89 = buffer.data(ppg + 89);
    const auto *ppg_95 = buffer.data(ppg + 95);
    const auto *ppg_99 = buffer.data(ppg + 99);
    const auto *ppg_101 = buffer.data(ppg + 101);
    const auto *ppg_102 = buffer.data(ppg + 102);
    const auto *ppg_104 = buffer.data(ppg + 104);
    const auto *ppg_112 = buffer.data(ppg + 112);
    const auto *ppg_115 = buffer.data(ppg + 115);
    const auto *ppg_116 = buffer.data(ppg + 116);
    const auto *ppg_117 = buffer.data(ppg + 117);
    const auto *ppg_119 = buffer.data(ppg + 119);
    const auto *ppg_120 = buffer.data(ppg + 120);
    const auto *ppg_125 = buffer.data(ppg + 125);
    const auto *ppg_129 = buffer.data(ppg + 129);
    const auto *ppg_130 = buffer.data(ppg + 130);
    const auto *ppg_131 = buffer.data(ppg + 131);
    const auto *ppg_132 = buffer.data(ppg + 132);
    const auto *ppg_134 = buffer.data(ppg + 134);

    const auto *pph1_3 = buffer.data(pph1 + 3);
    const auto *pph1_6 = buffer.data(pph1 + 6);
    const auto *pph1_10 = buffer.data(pph1 + 10);
    const auto *pph1_21 = buffer.data(pph1 + 21);
    const auto *pph1_24 = buffer.data(pph1 + 24);
    const auto *pph1_27 = buffer.data(pph1 + 27);
    const auto *pph1_64 = buffer.data(pph1 + 64);
    const auto *pph1_66 = buffer.data(pph1 + 66);
    const auto *pph1_126 = buffer.data(pph1 + 126);
    const auto *pph1_128 = buffer.data(pph1 + 128);
    const auto *pph1_154 = buffer.data(pph1 + 154);
    const auto *pph1_162 = buffer.data(pph1 + 162);
    const auto *pph1_163 = buffer.data(pph1 + 163);
    const auto *pph1_164 = buffer.data(pph1 + 164);
    const auto *pph1_165 = buffer.data(pph1 + 165);
    const auto *pph1_167 = buffer.data(pph1 + 167);
    const auto *pph1_173 = buffer.data(pph1 + 173);
    const auto *pph1_177 = buffer.data(pph1 + 177);
    const auto *pph1_183 = buffer.data(pph1 + 183);
    const auto *pph1_184 = buffer.data(pph1 + 184);
    const auto *pph1_185 = buffer.data(pph1 + 185);
    const auto *pph1_186 = buffer.data(pph1 + 186);
    const auto *pph1_188 = buffer.data(pph1 + 188);

    const auto *dsh0_47 = buffer.data(dsh0 + 47);
    const auto *dsh0_51 = buffer.data(dsh0 + 51);
    const auto *dsh0_63 = buffer.data(dsh0 + 63);
    const auto *dsh0_64 = buffer.data(dsh0 + 64);
    const auto *dsh0_66 = buffer.data(dsh0 + 66);
    const auto *dsh0_68 = buffer.data(dsh0 + 68);
    const auto *dsh0_69 = buffer.data(dsh0 + 69);
    const auto *dsh0_71 = buffer.data(dsh0 + 71);
    const auto *dsh0_72 = buffer.data(dsh0 + 72);
    const auto *dsh0_78 = buffer.data(dsh0 + 78);
    const auto *dsh0_80 = buffer.data(dsh0 + 80);
    const auto *dsh0_81 = buffer.data(dsh0 + 81);
    const auto *dsh0_83 = buffer.data(dsh0 + 83);

    const auto *dsg_30 = buffer.data(dsg + 30);
    const auto *dsg_32 = buffer.data(dsg + 32);
    const auto *dsg_35 = buffer.data(dsg + 35);
    const auto *dsg_39 = buffer.data(dsg + 39);
    const auto *dsg_41 = buffer.data(dsg + 41);
    const auto *dsg_42 = buffer.data(dsg + 42);
    const auto *dsg_44 = buffer.data(dsg + 44);
    const auto *dsg_45 = buffer.data(dsg + 45);
    const auto *dsg_46 = buffer.data(dsg + 46);
    const auto *dsg_48 = buffer.data(dsg + 48);
    const auto *dsg_50 = buffer.data(dsg + 50);
    const auto *dsg_51 = buffer.data(dsg + 51);
    const auto *dsg_53 = buffer.data(dsg + 53);
    const auto *dsg_54 = buffer.data(dsg + 54);
    const auto *dsg_55 = buffer.data(dsg + 55);
    const auto *dsg_56 = buffer.data(dsg + 56);
    const auto *dsg_57 = buffer.data(dsg + 57);
    const auto *dsg_58 = buffer.data(dsg + 58);
    const auto *dsg_59 = buffer.data(dsg + 59);

    const auto *dsh1_47 = buffer.data(dsh1 + 47);
    const auto *dsh1_51 = buffer.data(dsh1 + 51);
    const auto *dsh1_63 = buffer.data(dsh1 + 63);
    const auto *dsh1_64 = buffer.data(dsh1 + 64);
    const auto *dsh1_66 = buffer.data(dsh1 + 66);
    const auto *dsh1_68 = buffer.data(dsh1 + 68);
    const auto *dsh1_69 = buffer.data(dsh1 + 69);
    const auto *dsh1_71 = buffer.data(dsh1 + 71);
    const auto *dsh1_72 = buffer.data(dsh1 + 72);
    const auto *dsh1_78 = buffer.data(dsh1 + 78);
    const auto *dsh1_80 = buffer.data(dsh1 + 80);
    const auto *dsh1_81 = buffer.data(dsh1 + 81);
    const auto *dsh1_83 = buffer.data(dsh1 + 83);

    const auto *dpf0_62 = buffer.data(dpf0 + 62);
    const auto *dpf0_65 = buffer.data(dpf0 + 65);
    const auto *dpf0_66 = buffer.data(dpf0 + 66);
    const auto *dpf0_67 = buffer.data(dpf0 + 67);
    const auto *dpf0_68 = buffer.data(dpf0 + 68);
    const auto *dpf0_69 = buffer.data(dpf0 + 69);
    const auto *dpf0_80 = buffer.data(dpf0 + 80);
    const auto *dpf0_81 = buffer.data(dpf0 + 81);
    const auto *dpf0_82 = buffer.data(dpf0 + 82);
    const auto *dpf0_100 = buffer.data(dpf0 + 100);
    const auto *dpf0_101 = buffer.data(dpf0 + 101);
    const auto *dpf0_103 = buffer.data(dpf0 + 103);
    const auto *dpf0_105 = buffer.data(dpf0 + 105);
    const auto *dpf0_106 = buffer.data(dpf0 + 106);
    const auto *dpf0_107 = buffer.data(dpf0 + 107);
    const auto *dpf0_108 = buffer.data(dpf0 + 108);
    const auto *dpf0_109 = buffer.data(dpf0 + 109);
    const auto *dpf0_115 = buffer.data(dpf0 + 115);
    const auto *dpf0_118 = buffer.data(dpf0 + 118);
    const auto *dpf0_119 = buffer.data(dpf0 + 119);

    const auto *dpf1_62 = buffer.data(dpf1 + 62);
    const auto *dpf1_65 = buffer.data(dpf1 + 65);
    const auto *dpf1_66 = buffer.data(dpf1 + 66);
    const auto *dpf1_67 = buffer.data(dpf1 + 67);
    const auto *dpf1_68 = buffer.data(dpf1 + 68);
    const auto *dpf1_69 = buffer.data(dpf1 + 69);
    const auto *dpf1_80 = buffer.data(dpf1 + 80);
    const auto *dpf1_81 = buffer.data(dpf1 + 81);
    const auto *dpf1_82 = buffer.data(dpf1 + 82);
    const auto *dpf1_100 = buffer.data(dpf1 + 100);
    const auto *dpf1_101 = buffer.data(dpf1 + 101);
    const auto *dpf1_103 = buffer.data(dpf1 + 103);
    const auto *dpf1_105 = buffer.data(dpf1 + 105);
    const auto *dpf1_106 = buffer.data(dpf1 + 106);
    const auto *dpf1_107 = buffer.data(dpf1 + 107);
    const auto *dpf1_108 = buffer.data(dpf1 + 108);
    const auto *dpf1_109 = buffer.data(dpf1 + 109);
    const auto *dpf1_115 = buffer.data(dpf1 + 115);
    const auto *dpf1_118 = buffer.data(dpf1 + 118);
    const auto *dpf1_119 = buffer.data(dpf1 + 119);

    const auto *dpg_90 = buffer.data(dpg + 90);
    const auto *dpg_92 = buffer.data(dpg + 92);
    const auto *dpg_94 = buffer.data(dpg + 94);
    const auto *dpg_95 = buffer.data(dpg + 95);
    const auto *dpg_99 = buffer.data(dpg + 99);
    const auto *dpg_100 = buffer.data(dpg + 100);
    const auto *dpg_101 = buffer.data(dpg + 101);
    const auto *dpg_102 = buffer.data(dpg + 102);
    const auto *dpg_103 = buffer.data(dpg + 103);
    const auto *dpg_104 = buffer.data(dpg + 104);
    const auto *dpg_105 = buffer.data(dpg + 105);
    const auto *dpg_107 = buffer.data(dpg + 107);
    const auto *dpg_110 = buffer.data(dpg + 110);
    const auto *dpg_114 = buffer.data(dpg + 114);
    const auto *dpg_115 = buffer.data(dpg + 115);
    const auto *dpg_116 = buffer.data(dpg + 116);
    const auto *dpg_117 = buffer.data(dpg + 117);
    const auto *dpg_119 = buffer.data(dpg + 119);
    const auto *dpg_120 = buffer.data(dpg + 120);
    const auto *dpg_121 = buffer.data(dpg + 121);
    const auto *dpg_122 = buffer.data(dpg + 122);
    const auto *dpg_123 = buffer.data(dpg + 123);
    const auto *dpg_124 = buffer.data(dpg + 124);
    const auto *dpg_125 = buffer.data(dpg + 125);
    const auto *dpg_129 = buffer.data(dpg + 129);
    const auto *dpg_130 = buffer.data(dpg + 130);
    const auto *dpg_131 = buffer.data(dpg + 131);
    const auto *dpg_132 = buffer.data(dpg + 132);
    const auto *dpg_134 = buffer.data(dpg + 134);
    const auto *dpg_135 = buffer.data(dpg + 135);
    const auto *dpg_136 = buffer.data(dpg + 136);
    const auto *dpg_138 = buffer.data(dpg + 138);
    const auto *dpg_145 = buffer.data(dpg + 145);
    const auto *dpg_146 = buffer.data(dpg + 146);
    const auto *dpg_147 = buffer.data(dpg + 147);
    const auto *dpg_148 = buffer.data(dpg + 148);
    const auto *dpg_149 = buffer.data(dpg + 149);
    const auto *dpg_150 = buffer.data(dpg + 150);
    const auto *dpg_151 = buffer.data(dpg + 151);
    const auto *dpg_153 = buffer.data(dpg + 153);
    const auto *dpg_155 = buffer.data(dpg + 155);
    const auto *dpg_156 = buffer.data(dpg + 156);
    const auto *dpg_158 = buffer.data(dpg + 158);
    const auto *dpg_159 = buffer.data(dpg + 159);
    const auto *dpg_160 = buffer.data(dpg + 160);
    const auto *dpg_161 = buffer.data(dpg + 161);
    const auto *dpg_162 = buffer.data(dpg + 162);
    const auto *dpg_163 = buffer.data(dpg + 163);
    const auto *dpg_164 = buffer.data(dpg + 164);
    const auto *dpg_165 = buffer.data(dpg + 165);
    const auto *dpg_166 = buffer.data(dpg + 166);
    const auto *dpg_168 = buffer.data(dpg + 168);
    const auto *dpg_170 = buffer.data(dpg + 170);
    const auto *dpg_173 = buffer.data(dpg + 173);
    const auto *dpg_174 = buffer.data(dpg + 174);
    const auto *dpg_175 = buffer.data(dpg + 175);
    const auto *dpg_176 = buffer.data(dpg + 176);
    const auto *dpg_177 = buffer.data(dpg + 177);
    const auto *dpg_178 = buffer.data(dpg + 178);
    const auto *dpg_179 = buffer.data(dpg + 179);

#pragma omp simd aligned(t_128, t_129, t_130, pa_z, pc_y, pc_z, pph0_3, ppg_0, pph1_3, dpg_90, \
                         dpg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_1 * ppg_0[k]
                   + f_4 * pc_z[k] * dpg_90[k];

        t_129[k] = pa_z[k] * pph0_3[k]
                   - f_9 * pc_z[k] * pph1_3[k];

        t_130[k] = f_4 * pc_y[k] * dpg_92[k];
    }

#pragma omp simd aligned(t_131, t_132, pa_z, pc_x, pc_z, pph0_6, ppg_95, pph1_6, dsg_35, \
                         dpf0_65, dpf1_65, dpg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_1 * ppg_95[k]
                   + f_1 * dsg_35[k]
                   + f_7 * dpf0_65[k]
                   - f_8 * dpf1_65[k]
                   + f_4 * pc_x[k] * dpg_95[k];

        t_132[k] = pa_z[k] * pph0_6[k]
                   - f_9 * pc_z[k] * pph1_6[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pc_x, pc_y, ppg_99, dsg_39, dpf0_62, dpf0_69, \
                         dpf1_62, dpf1_69, dpg_94, dpg_95, dpg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_5 * dpf0_62[k]
                   - f_6 * dpf1_62[k]
                   + f_4 * pc_y[k] * dpg_94[k];

        t_134[k] = f_4 * pc_y[k] * dpg_95[k];

        t_135[k] = f_1 * ppg_99[k]
                   + f_1 * dsg_39[k]
                   + f_5 * dpf0_69[k]
                   - f_6 * dpf1_69[k]
                   + f_4 * pc_x[k] * dpg_99[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_z, pc_x, pc_z, pph0_10, ppg_101, ppg_102, \
                         pph1_10, dsg_41, dsg_42, dpg_101, dpg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_z[k] * pph0_10[k]
                   - f_9 * pc_z[k] * pph1_10[k];

        t_137[k] = f_1 * ppg_101[k]
                   + f_1 * dsg_41[k]
                   + f_4 * pc_x[k] * dpg_101[k];

        t_138[k] = f_1 * ppg_102[k]
                   + f_1 * dsg_42[k]
                   + f_4 * pc_x[k] * dpg_102[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, ppg_104, dsg_44, dpf0_66, \
                         dpf0_67, dpf1_66, dpf1_67, dpg_99, dpg_100, dpg_101, \
                         dpg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_4 * pc_y[k] * dpg_99[k];

        t_140[k] = f_1 * ppg_104[k]
                   + f_1 * dsg_44[k]
                   + f_4 * pc_x[k] * dpg_104[k];

        t_141[k] = f_2 * dpf0_66[k]
                   - f_3 * dpf1_66[k]
                   + f_4 * pc_y[k] * dpg_100[k];

        t_142[k] = f_12 * dpf0_67[k]
                   - f_13 * dpf1_67[k]
                   + f_4 * pc_y[k] * dpg_101[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_y, pc_z, ppg_14, dpf0_68, dpf0_69, \
                         dpf1_68, dpf1_69, dpg_102, dpg_103, dpg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_7 * dpf0_68[k]
                   - f_8 * dpf1_68[k]
                   + f_4 * pc_y[k] * dpg_102[k];

        t_144[k] = f_5 * dpf0_69[k]
                   - f_6 * dpf1_69[k]
                   + f_4 * pc_y[k] * dpg_103[k];

        t_145[k] = f_4 * pc_y[k] * dpg_104[k];

        t_146[k] = f_1 * ppg_14[k]
                   + f_2 * dpf0_69[k]
                   - f_3 * dpf1_69[k]
                   + f_4 * pc_z[k] * dpg_104[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_z, pc_y, pc_z, pph0_21, pph0_24, \
                         ppg_15, pph1_21, pph1_24, dsg_30, dpg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pa_z[k] * pph0_21[k]
                   - f_9 * pc_z[k] * pph1_21[k];

        t_148[k] = f_1 * dsg_30[k]
                   + f_4 * pc_y[k] * dpg_105[k];

        t_149[k] = f_1 * ppg_15[k]
                   + f_4 * pc_z[k] * dpg_105[k];

        t_150[k] = pa_z[k] * pph0_24[k]
                   - f_9 * pc_z[k] * pph1_24[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_z, pb_y, pc_y, pc_z, pph0_27, pph1_27, \
                         dsh0_47, dsg_32, dsh1_47, dpg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_1 * dsg_32[k]
                   + f_4 * pc_y[k] * dpg_107[k];

        t_152[k] = pb_y[k] * dsh0_47[k]
                   - f_9 * pc_y[k] * dsh1_47[k];

        t_153[k] = pa_z[k] * pph0_27[k]
                   - f_9 * pc_z[k] * pph1_27[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_x, pb_y, pc_x, pc_y, pph0_154, ppg_112, \
                         pph1_154, dsh0_51, dsg_35, dsh1_51, dpg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pa_x[k] * pph0_154[k]
                   + f_0 * ppg_112[k]
                   - f_9 * pc_x[k] * pph1_154[k];

        t_155[k] = f_1 * dsg_35[k]
                   + f_4 * pc_y[k] * dpg_110[k];

        t_156[k] = pb_y[k] * dsh0_51[k]
                   - f_9 * pc_y[k] * dsh1_51[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, pc_y, ppg_115, ppg_116, ppg_117, \
                         dsg_39, dpg_114, dpg_115, dpg_116, dpg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_1 * ppg_115[k]
                   + f_4 * pc_x[k] * dpg_115[k];

        t_158[k] = f_1 * ppg_116[k]
                   + f_4 * pc_x[k] * dpg_116[k];

        t_159[k] = f_1 * ppg_117[k]
                   + f_4 * pc_x[k] * dpg_117[k];

        t_160[k] = f_1 * dsg_39[k]
                   + f_4 * pc_y[k] * dpg_114[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pa_x, pc_x, pph0_162, pph0_163, pph0_164, \
                         ppg_119, pph1_162, pph1_163, pph1_164, \
                         dpg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * ppg_119[k]
                   + f_4 * pc_x[k] * dpg_119[k];

        t_162[k] = pa_x[k] * pph0_162[k]
                   - f_9 * pc_x[k] * pph1_162[k];

        t_163[k] = pa_x[k] * pph0_163[k]
                   - f_9 * pc_x[k] * pph1_163[k];

        t_164[k] = pa_x[k] * pph0_164[k]
                   - f_9 * pc_x[k] * pph1_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_x, pc_x, pc_y, pph0_165, pph0_167, pph1_165, \
                         pph1_167, dsg_44, dpg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pa_x[k] * pph0_165[k]
                   - f_9 * pc_x[k] * pph1_165[k];

        t_166[k] = f_1 * dsg_44[k]
                   + f_4 * pc_y[k] * dpg_119[k];

        t_167[k] = pa_x[k] * pph0_167[k]
                   - f_9 * pc_x[k] * pph1_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pc_x, pc_y, pc_z, ppg_30, ppg_120, \
                         dsg_30, dpf0_80, dpf1_80, dpg_120, dpg_121, \
                         dpg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_1 * ppg_120[k]
                   + f_2 * dpf0_80[k]
                   - f_3 * dpf1_80[k]
                   + f_4 * pc_x[k] * dpg_120[k];

        t_169[k] = f_4 * pc_y[k] * dpg_120[k];

        t_170[k] = f_1 * ppg_30[k]
                   + f_1 * dsg_30[k]
                   + f_4 * pc_z[k] * dpg_120[k];

        t_171[k] = f_5 * dpf0_80[k]
                   - f_6 * dpf1_80[k]
                   + f_4 * pc_y[k] * dpg_121[k];

        t_172[k] = f_4 * pc_y[k] * dpg_122[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, pa_x, pc_x, pc_y, pph0_173, ppg_125, pph1_173, \
                         dpf0_81, dpf0_82, dpf1_81, dpf1_82, dpg_123, \
                         dpg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = pa_x[k] * pph0_173[k]
                   + f_10 * ppg_125[k]
                   - f_9 * pc_x[k] * pph1_173[k];

        t_174[k] = f_7 * dpf0_81[k]
                   - f_8 * dpf1_81[k]
                   + f_4 * pc_y[k] * dpg_123[k];

        t_175[k] = f_5 * dpf0_82[k]
                   - f_6 * dpf1_82[k]
                   + f_4 * pc_y[k] * dpg_124[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_x, pc_x, pc_y, pph0_177, ppg_129, \
                         ppg_130, ppg_131, pph1_177, dpg_125, dpg_130, \
                         dpg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_4 * pc_y[k] * dpg_125[k];

        t_177[k] = pa_x[k] * pph0_177[k]
                   + f_0 * ppg_129[k]
                   - f_9 * pc_x[k] * pph1_177[k];

        t_178[k] = f_1 * ppg_130[k]
                   + f_4 * pc_x[k] * dpg_130[k];

        t_179[k] = f_1 * ppg_131[k]
                   + f_4 * pc_x[k] * dpg_131[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_x, pc_x, pc_y, pph0_183, ppg_132, \
                         ppg_134, pph1_183, dpg_129, dpg_132, dpg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * ppg_132[k]
                   + f_4 * pc_x[k] * dpg_132[k];

        t_181[k] = f_4 * pc_y[k] * dpg_129[k];

        t_182[k] = f_1 * ppg_134[k]
                   + f_4 * pc_x[k] * dpg_134[k];

        t_183[k] = pa_x[k] * pph0_183[k]
                   - f_9 * pc_x[k] * pph1_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_x, pc_x, pc_y, pph0_184, pph0_185, \
                         pph0_186, pph1_184, pph1_185, pph1_186, \
                         dpg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pa_x[k] * pph0_184[k]
                   - f_9 * pc_x[k] * pph1_184[k];

        t_185[k] = pa_x[k] * pph0_185[k]
                   - f_9 * pc_x[k] * pph1_185[k];

        t_186[k] = pa_x[k] * pph0_186[k]
                   - f_9 * pc_x[k] * pph1_186[k];

        t_187[k] = f_4 * pc_y[k] * dpg_134[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_x, pb_x, pc_x, pph0_188, pph1_188, dsh0_63, \
                         dsh0_64, dsg_45, dsg_46, dsh1_63, dsh1_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_x[k] * pph0_188[k]
                   - f_9 * pc_x[k] * pph1_188[k];

        t_189[k] = pb_x[k] * dsh0_63[k]
                   + f_11 * dsg_45[k]
                   - f_9 * pc_x[k] * dsh1_63[k];

        t_190[k] = pb_x[k] * dsh0_64[k]
                   + f_14 * dsg_46[k]
                   - f_9 * pc_x[k] * dsh1_64[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pb_x, pc_x, pc_z, dsh0_66, dsh0_68, \
                         dsg_48, dsg_50, dsh1_66, dsh1_68, dpg_135, \
                         dpg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_4 * pc_z[k] * dpg_135[k];

        t_192[k] = pb_x[k] * dsh0_66[k]
                   + f_10 * dsg_48[k]
                   - f_9 * pc_x[k] * dsh1_66[k];

        t_193[k] = f_4 * pc_z[k] * dpg_136[k];

        t_194[k] = pb_x[k] * dsh0_68[k]
                   + f_10 * dsg_50[k]
                   - f_9 * pc_x[k] * dsh1_68[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_x, pc_x, pc_z, dsh0_69, dsh0_71, dsg_51, \
                         dsg_53, dsh1_69, dsh1_71, dpg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_x[k] * dsh0_69[k]
                   + f_0 * dsg_51[k]
                   - f_9 * pc_x[k] * dsh1_69[k];

        t_196[k] = f_4 * pc_z[k] * dpg_138[k];

        t_197[k] = pb_x[k] * dsh0_71[k]
                   + f_0 * dsg_53[k]
                   - f_9 * pc_x[k] * dsh1_71[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pb_x, pc_x, dsh0_72, dsg_54, dsg_55, \
                         dsg_56, dsg_57, dsh1_72, dpg_145, dpg_146, \
                         dpg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_x[k] * dsh0_72[k]
                   + f_0 * dsg_54[k]
                   - f_9 * pc_x[k] * dsh1_72[k];

        t_199[k] = f_1 * dsg_55[k]
                   + f_4 * pc_x[k] * dpg_145[k];

        t_200[k] = f_1 * dsg_56[k]
                   + f_4 * pc_x[k] * dpg_146[k];

        t_201[k] = f_1 * dsg_57[k]
                   + f_4 * pc_x[k] * dpg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_x, pc_x, pc_z, dsh0_78, dsg_58, \
                         dsg_59, dsh1_78, dpg_145, dpg_148, dpg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_1 * dsg_58[k]
                   + f_4 * pc_x[k] * dpg_148[k];

        t_203[k] = f_1 * dsg_59[k]
                   + f_4 * pc_x[k] * dpg_149[k];

        t_204[k] = pb_x[k] * dsh0_78[k]
                   - f_9 * pc_x[k] * dsh1_78[k];

        t_205[k] = f_4 * pc_z[k] * dpg_145[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_x, pc_x, pc_y, ppg_59, dsh0_80, \
                         dsh0_81, dsh0_83, dsh1_80, dsh1_81, dsh1_83, \
                         dpg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_x[k] * dsh0_80[k]
                   - f_9 * pc_x[k] * dsh1_80[k];

        t_207[k] = pb_x[k] * dsh0_81[k]
                   - f_9 * pc_x[k] * dsh1_81[k];

        t_208[k] = f_0 * ppg_59[k]
                   + f_4 * pc_y[k] * dpg_149[k];

        t_209[k] = pb_x[k] * dsh0_83[k]
                   - f_9 * pc_x[k] * dsh1_83[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, pc_z, dpf0_100, dpf0_101, \
                         dpf0_103, dpf1_100, dpf1_101, dpf1_103, dpg_150, dpg_151, \
                         dpg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_2 * dpf0_100[k]
                   - f_3 * dpf1_100[k]
                   + f_4 * pc_x[k] * dpg_150[k];

        t_211[k] = f_12 * dpf0_101[k]
                   - f_13 * dpf1_101[k]
                   + f_4 * pc_x[k] * dpg_151[k];

        t_212[k] = f_4 * pc_z[k] * dpg_150[k];

        t_213[k] = f_7 * dpf0_103[k]
                   - f_8 * dpf1_103[k]
                   + f_4 * pc_x[k] * dpg_153[k];

        t_214[k] = f_4 * pc_z[k] * dpg_151[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, pc_z, dpf0_105, dpf0_106, dpf0_108, \
                         dpf1_105, dpf1_106, dpf1_108, dpg_153, dpg_155, dpg_156, \
                         dpg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_7 * dpf0_105[k]
                   - f_8 * dpf1_105[k]
                   + f_4 * pc_x[k] * dpg_155[k];

        t_216[k] = f_5 * dpf0_106[k]
                   - f_6 * dpf1_106[k]
                   + f_4 * pc_x[k] * dpg_156[k];

        t_217[k] = f_4 * pc_z[k] * dpg_153[k];

        t_218[k] = f_5 * dpf0_108[k]
                   - f_6 * dpf1_108[k]
                   + f_4 * pc_x[k] * dpg_158[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, t_224, pc_x, dpf0_109, dpf1_109, \
                         dpg_159, dpg_160, dpg_161, dpg_162, dpg_163, \
                         dpg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_5 * dpf0_109[k]
                   - f_6 * dpf1_109[k]
                   + f_4 * pc_x[k] * dpg_159[k];

        t_220[k] = f_4 * pc_x[k] * dpg_160[k];

        t_221[k] = f_4 * pc_x[k] * dpg_161[k];

        t_222[k] = f_4 * pc_x[k] * dpg_162[k];

        t_223[k] = f_4 * pc_x[k] * dpg_163[k];

        t_224[k] = f_4 * pc_x[k] * dpg_164[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pc_y, pc_z, ppg_70, dsg_55, dpf0_106, \
                         dpf0_107, dpf1_106, dpf1_107, dpg_160, dpg_161, \
                         dpg_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_0 * ppg_70[k]
                   + f_1 * dsg_55[k]
                   + f_2 * dpf0_106[k]
                   - f_3 * dpf1_106[k]
                   + f_4 * pc_y[k] * dpg_160[k];

        t_226[k] = f_4 * pc_z[k] * dpg_160[k];

        t_227[k] = f_5 * dpf0_106[k]
                   - f_6 * dpf1_106[k]
                   + f_4 * pc_z[k] * dpg_161[k];

        t_228[k] = f_7 * dpf0_107[k]
                   - f_8 * dpf1_107[k]
                   + f_4 * pc_z[k] * dpg_162[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_z, pc_y, pc_z, ppg_74, dsh0_63, \
                         dsh0_64, dsg_59, dsh1_63, dsh1_64, dpf0_109, dpf1_109, \
                         dpg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * ppg_74[k]
                   + f_1 * dsg_59[k]
                   + f_4 * pc_y[k] * dpg_164[k];

        t_230[k] = f_2 * dpf0_109[k]
                   - f_3 * dpf1_109[k]
                   + f_4 * pc_z[k] * dpg_164[k];

        t_231[k] = pb_z[k] * dsh0_63[k]
                   - f_9 * pc_z[k] * dsh1_63[k];

        t_232[k] = pb_z[k] * dsh0_64[k]
                   - f_9 * pc_z[k] * dsh1_64[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_z, pc_x, pc_z, dsh0_66, dsg_45, \
                         dsg_46, dsh1_66, dpf0_115, dpf1_115, dpg_165, dpg_166, \
                         dpg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_1 * dsg_45[k]
                   + f_4 * pc_z[k] * dpg_165[k];

        t_234[k] = pb_z[k] * dsh0_66[k]
                   - f_9 * pc_z[k] * dsh1_66[k];

        t_235[k] = f_1 * dsg_46[k]
                   + f_4 * pc_z[k] * dpg_166[k];

        t_236[k] = f_7 * dpf0_115[k]
                   - f_8 * dpf1_115[k]
                   + f_4 * pc_x[k] * dpg_170[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pb_z, pc_x, pc_z, dsh0_69, dsg_48, dsh1_69, \
                         dpf0_118, dpf1_118, dpg_168, dpg_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pb_z[k] * dsh0_69[k]
                   - f_9 * pc_z[k] * dsh1_69[k];

        t_238[k] = f_1 * dsg_48[k]
                   + f_4 * pc_z[k] * dpg_168[k];

        t_239[k] = f_5 * dpf0_118[k]
                   - f_6 * dpf1_118[k]
                   + f_4 * pc_x[k] * dpg_173[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, t_245, pc_x, dpf0_119, dpf1_119, \
                         dpg_174, dpg_175, dpg_176, dpg_177, dpg_178, \
                         dpg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_5 * dpf0_119[k]
                   - f_6 * dpf1_119[k]
                   + f_4 * pc_x[k] * dpg_174[k];

        t_241[k] = f_4 * pc_x[k] * dpg_175[k];

        t_242[k] = f_4 * pc_x[k] * dpg_176[k];

        t_243[k] = f_4 * pc_x[k] * dpg_177[k];

        t_244[k] = f_4 * pc_x[k] * dpg_178[k];

        t_245[k] = f_4 * pc_x[k] * dpg_179[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pb_z, pc_z, dsh0_78, dsh0_80, dsh0_81, \
                         dsg_55, dsg_56, dsg_57, dsh1_78, dsh1_80, dsh1_81, \
                         dpg_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = pb_z[k] * dsh0_78[k]
                   - f_9 * pc_z[k] * dsh1_78[k];

        t_247[k] = f_1 * dsg_55[k]
                   + f_4 * pc_z[k] * dpg_175[k];

        t_248[k] = pb_z[k] * dsh0_80[k]
                   + f_0 * dsg_56[k]
                   - f_9 * pc_z[k] * dsh1_80[k];

        t_249[k] = pb_z[k] * dsh0_81[k]
                   + f_10 * dsg_57[k]
                   - f_9 * pc_z[k] * dsh1_81[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pa_y, pb_z, pc_y, pc_z, pph0_126, ppg_89, \
                         pph1_126, dsh0_83, dsg_59, dsh1_83, dpg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_0 * ppg_89[k]
                   + f_4 * pc_y[k] * dpg_179[k];

        t_251[k] = pb_z[k] * dsh0_83[k]
                   + f_11 * dsg_59[k]
                   - f_9 * pc_z[k] * dsh1_83[k];

        t_252[k] = pa_y[k] * pph0_126[k]
                   - f_9 * pc_y[k] * pph1_126[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, pa_y, pa_z, pc_y, pc_z, pph0_64, pph0_66, \
                         pph0_128, pph1_64, pph1_66, pph1_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = pa_z[k] * pph0_64[k]
                   - f_9 * pc_z[k] * pph1_64[k];

        t_254[k] = pa_y[k] * pph0_128[k]
                   - f_9 * pc_y[k] * pph1_128[k];

        t_255[k] = pa_z[k] * pph0_66[k]
                   - f_9 * pc_z[k] * pph1_66[k];
    }
}

static auto
compute_prim_dph_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t pph0, const size_t ppg,
                                                          const size_t pph1, const size_t dsh0,
                                                          const size_t dsg, const size_t dsh1,
                                                          const size_t dpf0, const size_t dpf1,
                                                          const size_t dpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);
    const auto f_14 = 2.0 / q;

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
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *pph0_69 = buffer.data(pph0 + 69);
    const auto *pph0_78 = buffer.data(pph0 + 78);
    const auto *pph0_85 = buffer.data(pph0 + 85);
    const auto *pph0_87 = buffer.data(pph0 + 87);
    const auto *pph0_90 = buffer.data(pph0 + 90);
    const auto *pph0_99 = buffer.data(pph0 + 99);
    const auto *pph0_130 = buffer.data(pph0 + 130);
    const auto *pph0_131 = buffer.data(pph0 + 131);
    const auto *pph0_133 = buffer.data(pph0 + 133);
    const auto *pph0_134 = buffer.data(pph0 + 134);
    const auto *pph0_135 = buffer.data(pph0 + 135);
    const auto *pph0_146 = buffer.data(pph0 + 146);
    const auto *pph0_168 = buffer.data(pph0 + 168);
    const auto *pph0_170 = buffer.data(pph0 + 170);
    const auto *pph0_173 = buffer.data(pph0 + 173);
    const auto *pph0_177 = buffer.data(pph0 + 177);
    const auto *pph0_185 = buffer.data(pph0 + 185);
    const auto *pph0_186 = buffer.data(pph0 + 186);
    const auto *pph0_188 = buffer.data(pph0 + 188);

    const auto *ppg_55 = buffer.data(ppg + 55);
    const auto *ppg_70 = buffer.data(ppg + 70);
    const auto *ppg_74 = buffer.data(ppg + 74);
    const auto *ppg_85 = buffer.data(ppg + 85);
    const auto *ppg_92 = buffer.data(ppg + 92);
    const auto *ppg_94 = buffer.data(ppg + 94);
    const auto *ppg_95 = buffer.data(ppg + 95);
    const auto *ppg_104 = buffer.data(ppg + 104);
    const auto *ppg_117 = buffer.data(ppg + 117);
    const auto *ppg_118 = buffer.data(ppg + 118);
    const auto *ppg_119 = buffer.data(ppg + 119);
    const auto *ppg_130 = buffer.data(ppg + 130);
    const auto *ppg_132 = buffer.data(ppg + 132);
    const auto *ppg_133 = buffer.data(ppg + 133);
    const auto *ppg_134 = buffer.data(ppg + 134);

    const auto *pph1_69 = buffer.data(pph1 + 69);
    const auto *pph1_78 = buffer.data(pph1 + 78);
    const auto *pph1_85 = buffer.data(pph1 + 85);
    const auto *pph1_87 = buffer.data(pph1 + 87);
    const auto *pph1_90 = buffer.data(pph1 + 90);
    const auto *pph1_99 = buffer.data(pph1 + 99);
    const auto *pph1_130 = buffer.data(pph1 + 130);
    const auto *pph1_131 = buffer.data(pph1 + 131);
    const auto *pph1_133 = buffer.data(pph1 + 133);
    const auto *pph1_134 = buffer.data(pph1 + 134);
    const auto *pph1_135 = buffer.data(pph1 + 135);
    const auto *pph1_146 = buffer.data(pph1 + 146);
    const auto *pph1_168 = buffer.data(pph1 + 168);
    const auto *pph1_170 = buffer.data(pph1 + 170);
    const auto *pph1_173 = buffer.data(pph1 + 173);
    const auto *pph1_177 = buffer.data(pph1 + 177);
    const auto *pph1_185 = buffer.data(pph1 + 185);
    const auto *pph1_186 = buffer.data(pph1 + 186);
    const auto *pph1_188 = buffer.data(pph1 + 188);

    const auto *dsh0_101 = buffer.data(dsh0 + 101);
    const auto *dsh0_102 = buffer.data(dsh0 + 102);
    const auto *dsh0_105 = buffer.data(dsh0 + 105);
    const auto *dsh0_107 = buffer.data(dsh0 + 107);
    const auto *dsh0_108 = buffer.data(dsh0 + 108);
    const auto *dsh0_110 = buffer.data(dsh0 + 110);
    const auto *dsh0_111 = buffer.data(dsh0 + 111);
    const auto *dsh0_112 = buffer.data(dsh0 + 112);
    const auto *dsh0_114 = buffer.data(dsh0 + 114);
    const auto *dsh0_120 = buffer.data(dsh0 + 120);
    const auto *dsh0_121 = buffer.data(dsh0 + 121);
    const auto *dsh0_122 = buffer.data(dsh0 + 122);
    const auto *dsh0_123 = buffer.data(dsh0 + 123);
    const auto *dsh0_125 = buffer.data(dsh0 + 125);

    const auto *dsg_70 = buffer.data(dsg + 70);
    const auto *dsg_71 = buffer.data(dsg + 71);
    const auto *dsg_72 = buffer.data(dsg + 72);
    const auto *dsg_73 = buffer.data(dsg + 73);
    const auto *dsg_74 = buffer.data(dsg + 74);
    const auto *dsg_75 = buffer.data(dsg + 75);
    const auto *dsg_77 = buffer.data(dsg + 77);
    const auto *dsg_78 = buffer.data(dsg + 78);
    const auto *dsg_80 = buffer.data(dsg + 80);
    const auto *dsg_81 = buffer.data(dsg + 81);
    const auto *dsg_82 = buffer.data(dsg + 82);
    const auto *dsg_84 = buffer.data(dsg + 84);
    const auto *dsg_85 = buffer.data(dsg + 85);
    const auto *dsg_86 = buffer.data(dsg + 86);
    const auto *dsg_87 = buffer.data(dsg + 87);
    const auto *dsg_88 = buffer.data(dsg + 88);
    const auto *dsg_89 = buffer.data(dsg + 89);

    const auto *dsh1_101 = buffer.data(dsh1 + 101);
    const auto *dsh1_102 = buffer.data(dsh1 + 102);
    const auto *dsh1_105 = buffer.data(dsh1 + 105);
    const auto *dsh1_107 = buffer.data(dsh1 + 107);
    const auto *dsh1_108 = buffer.data(dsh1 + 108);
    const auto *dsh1_110 = buffer.data(dsh1 + 110);
    const auto *dsh1_111 = buffer.data(dsh1 + 111);
    const auto *dsh1_112 = buffer.data(dsh1 + 112);
    const auto *dsh1_114 = buffer.data(dsh1 + 114);
    const auto *dsh1_120 = buffer.data(dsh1 + 120);
    const auto *dsh1_121 = buffer.data(dsh1 + 121);
    const auto *dsh1_122 = buffer.data(dsh1 + 122);
    const auto *dsh1_123 = buffer.data(dsh1 + 123);
    const auto *dsh1_125 = buffer.data(dsh1 + 125);

    const auto *dpf0_130 = buffer.data(dpf0 + 130);
    const auto *dpf0_132 = buffer.data(dpf0 + 132);
    const auto *dpf0_134 = buffer.data(dpf0 + 134);
    const auto *dpf0_135 = buffer.data(dpf0 + 135);
    const auto *dpf0_137 = buffer.data(dpf0 + 137);
    const auto *dpf0_138 = buffer.data(dpf0 + 138);
    const auto *dpf0_139 = buffer.data(dpf0 + 139);
    const auto *dpf0_141 = buffer.data(dpf0 + 141);
    const auto *dpf0_143 = buffer.data(dpf0 + 143);
    const auto *dpf0_144 = buffer.data(dpf0 + 144);
    const auto *dpf0_146 = buffer.data(dpf0 + 146);
    const auto *dpf0_147 = buffer.data(dpf0 + 147);
    const auto *dpf0_148 = buffer.data(dpf0 + 148);
    const auto *dpf0_163 = buffer.data(dpf0 + 163);
    const auto *dpf0_166 = buffer.data(dpf0 + 166);
    const auto *dpf0_167 = buffer.data(dpf0 + 167);
    const auto *dpf0_170 = buffer.data(dpf0 + 170);
    const auto *dpf0_172 = buffer.data(dpf0 + 172);
    const auto *dpf0_173 = buffer.data(dpf0 + 173);
    const auto *dpf0_175 = buffer.data(dpf0 + 175);
    const auto *dpf0_176 = buffer.data(dpf0 + 176);
    const auto *dpf0_177 = buffer.data(dpf0 + 177);
    const auto *dpf0_178 = buffer.data(dpf0 + 178);
    const auto *dpf0_179 = buffer.data(dpf0 + 179);

    const auto *dpf1_130 = buffer.data(dpf1 + 130);
    const auto *dpf1_132 = buffer.data(dpf1 + 132);
    const auto *dpf1_134 = buffer.data(dpf1 + 134);
    const auto *dpf1_135 = buffer.data(dpf1 + 135);
    const auto *dpf1_137 = buffer.data(dpf1 + 137);
    const auto *dpf1_138 = buffer.data(dpf1 + 138);
    const auto *dpf1_139 = buffer.data(dpf1 + 139);
    const auto *dpf1_141 = buffer.data(dpf1 + 141);
    const auto *dpf1_143 = buffer.data(dpf1 + 143);
    const auto *dpf1_144 = buffer.data(dpf1 + 144);
    const auto *dpf1_146 = buffer.data(dpf1 + 146);
    const auto *dpf1_147 = buffer.data(dpf1 + 147);
    const auto *dpf1_148 = buffer.data(dpf1 + 148);
    const auto *dpf1_163 = buffer.data(dpf1 + 163);
    const auto *dpf1_166 = buffer.data(dpf1 + 166);
    const auto *dpf1_167 = buffer.data(dpf1 + 167);
    const auto *dpf1_170 = buffer.data(dpf1 + 170);
    const auto *dpf1_172 = buffer.data(dpf1 + 172);
    const auto *dpf1_173 = buffer.data(dpf1 + 173);
    const auto *dpf1_175 = buffer.data(dpf1 + 175);
    const auto *dpf1_176 = buffer.data(dpf1 + 176);
    const auto *dpf1_177 = buffer.data(dpf1 + 177);
    const auto *dpf1_178 = buffer.data(dpf1 + 178);
    const auto *dpf1_179 = buffer.data(dpf1 + 179);

    const auto *dpg_190 = buffer.data(dpg + 190);
    const auto *dpg_191 = buffer.data(dpg + 191);
    const auto *dpg_192 = buffer.data(dpg + 192);
    const auto *dpg_193 = buffer.data(dpg + 193);
    const auto *dpg_194 = buffer.data(dpg + 194);
    const auto *dpg_195 = buffer.data(dpg + 195);
    const auto *dpg_197 = buffer.data(dpg + 197);
    const auto *dpg_199 = buffer.data(dpg + 199);
    const auto *dpg_200 = buffer.data(dpg + 200);
    const auto *dpg_202 = buffer.data(dpg + 202);
    const auto *dpg_203 = buffer.data(dpg + 203);
    const auto *dpg_204 = buffer.data(dpg + 204);
    const auto *dpg_205 = buffer.data(dpg + 205);
    const auto *dpg_206 = buffer.data(dpg + 206);
    const auto *dpg_207 = buffer.data(dpg + 207);
    const auto *dpg_208 = buffer.data(dpg + 208);
    const auto *dpg_209 = buffer.data(dpg + 209);
    const auto *dpg_211 = buffer.data(dpg + 211);
    const auto *dpg_213 = buffer.data(dpg + 213);
    const auto *dpg_214 = buffer.data(dpg + 214);
    const auto *dpg_216 = buffer.data(dpg + 216);
    const auto *dpg_217 = buffer.data(dpg + 217);
    const auto *dpg_218 = buffer.data(dpg + 218);
    const auto *dpg_220 = buffer.data(dpg + 220);
    const auto *dpg_221 = buffer.data(dpg + 221);
    const auto *dpg_222 = buffer.data(dpg + 222);
    const auto *dpg_223 = buffer.data(dpg + 223);
    const auto *dpg_224 = buffer.data(dpg + 224);
    const auto *dpg_225 = buffer.data(dpg + 225);
    const auto *dpg_227 = buffer.data(dpg + 227);
    const auto *dpg_230 = buffer.data(dpg + 230);
    const auto *dpg_235 = buffer.data(dpg + 235);
    const auto *dpg_236 = buffer.data(dpg + 236);
    const auto *dpg_237 = buffer.data(dpg + 237);
    const auto *dpg_238 = buffer.data(dpg + 238);
    const auto *dpg_239 = buffer.data(dpg + 239);
    const auto *dpg_240 = buffer.data(dpg + 240);
    const auto *dpg_242 = buffer.data(dpg + 242);
    const auto *dpg_243 = buffer.data(dpg + 243);
    const auto *dpg_245 = buffer.data(dpg + 245);
    const auto *dpg_246 = buffer.data(dpg + 246);
    const auto *dpg_247 = buffer.data(dpg + 247);
    const auto *dpg_250 = buffer.data(dpg + 250);
    const auto *dpg_251 = buffer.data(dpg + 251);
    const auto *dpg_252 = buffer.data(dpg + 252);
    const auto *dpg_253 = buffer.data(dpg + 253);
    const auto *dpg_254 = buffer.data(dpg + 254);
    const auto *dpg_255 = buffer.data(dpg + 255);
    const auto *dpg_257 = buffer.data(dpg + 257);
    const auto *dpg_258 = buffer.data(dpg + 258);
    const auto *dpg_260 = buffer.data(dpg + 260);
    const auto *dpg_261 = buffer.data(dpg + 261);
    const auto *dpg_262 = buffer.data(dpg + 262);
    const auto *dpg_264 = buffer.data(dpg + 264);
    const auto *dpg_265 = buffer.data(dpg + 265);
    const auto *dpg_266 = buffer.data(dpg + 266);
    const auto *dpg_267 = buffer.data(dpg + 267);
    const auto *dpg_268 = buffer.data(dpg + 268);
    const auto *dpg_269 = buffer.data(dpg + 269);

#pragma omp simd aligned(t_256, t_257, t_258, pa_y, pa_z, pc_y, pc_z, pph0_69, pph0_130, \
                         pph0_131, ppg_92, pph1_69, pph1_130, \
                         pph1_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = pa_y[k] * pph0_130[k]
                   + f_1 * ppg_92[k]
                   - f_9 * pc_y[k] * pph1_130[k];

        t_257[k] = pa_y[k] * pph0_131[k]
                   - f_9 * pc_y[k] * pph1_131[k];

        t_258[k] = pa_z[k] * pph0_69[k]
                   - f_9 * pc_z[k] * pph1_69[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, pa_y, pc_y, pph0_133, pph0_134, pph0_135, \
                         ppg_94, ppg_95, pph1_133, pph1_134, pph1_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = pa_y[k] * pph0_133[k]
                   + f_0 * ppg_94[k]
                   - f_9 * pc_y[k] * pph1_133[k];

        t_260[k] = pa_y[k] * pph0_134[k]
                   + f_1 * ppg_95[k]
                   - f_9 * pc_y[k] * pph1_134[k];

        t_261[k] = pa_y[k] * pph0_135[k]
                   - f_9 * pc_y[k] * pph1_135[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, pc_x, dsg_70, dsg_71, dsg_72, \
                         dsg_73, dsg_74, dpg_190, dpg_191, dpg_192, dpg_193, \
                         dpg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_1 * dsg_70[k]
                   + f_4 * pc_x[k] * dpg_190[k];

        t_263[k] = f_1 * dsg_71[k]
                   + f_4 * pc_x[k] * dpg_191[k];

        t_264[k] = f_1 * dsg_72[k]
                   + f_4 * pc_x[k] * dpg_192[k];

        t_265[k] = f_1 * dsg_73[k]
                   + f_4 * pc_x[k] * dpg_193[k];

        t_266[k] = f_1 * dsg_74[k]
                   + f_4 * pc_x[k] * dpg_194[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pa_z, pb_x, pc_x, pc_z, pph0_78, ppg_55, \
                         pph1_78, dsh0_101, dsh0_102, dsh1_101, dsh1_102, \
                         dpg_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pa_z[k] * pph0_78[k]
                   - f_9 * pc_z[k] * pph1_78[k];

        t_268[k] = f_1 * ppg_55[k]
                   + f_4 * pc_z[k] * dpg_190[k];

        t_269[k] = pb_x[k] * dsh0_101[k]
                   - f_9 * pc_x[k] * dsh1_101[k];

        t_270[k] = pb_x[k] * dsh0_102[k]
                   - f_9 * pc_x[k] * dsh1_102[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, pa_y, pc_x, pc_y, pph0_146, ppg_104, pph1_146, \
                         dpf0_130, dpf1_130, dpg_194, dpg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_1 * ppg_104[k]
                   + f_4 * pc_y[k] * dpg_194[k];

        t_272[k] = pa_y[k] * pph0_146[k]
                   - f_9 * pc_y[k] * pph1_146[k];

        t_273[k] = f_2 * dpf0_130[k]
                   - f_3 * dpf1_130[k]
                   + f_4 * pc_x[k] * dpg_195[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, pa_z, pc_x, pc_z, pph0_85, pph0_87, pph1_85, \
                         pph1_87, dpf0_132, dpf1_132, dpg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = pa_z[k] * pph0_85[k]
                   - f_9 * pc_z[k] * pph1_85[k];

        t_275[k] = f_12 * dpf0_132[k]
                   - f_13 * dpf1_132[k]
                   + f_4 * pc_x[k] * dpg_197[k];

        t_276[k] = pa_z[k] * pph0_87[k]
                   - f_9 * pc_z[k] * pph1_87[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, pa_z, pc_x, pc_z, pph0_90, pph1_90, dpf0_134, \
                         dpf0_135, dpf1_134, dpf1_135, dpg_199, \
                         dpg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_7 * dpf0_134[k]
                   - f_8 * dpf1_134[k]
                   + f_4 * pc_x[k] * dpg_199[k];

        t_278[k] = f_7 * dpf0_135[k]
                   - f_8 * dpf1_135[k]
                   + f_4 * pc_x[k] * dpg_200[k];

        t_279[k] = pa_z[k] * pph0_90[k]
                   - f_9 * pc_z[k] * pph1_90[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, dpf0_137, dpf0_138, dpf0_139, \
                         dpf1_137, dpf1_138, dpf1_139, dpg_202, dpg_203, dpg_204, \
                         dpg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_5 * dpf0_137[k]
                   - f_6 * dpf1_137[k]
                   + f_4 * pc_x[k] * dpg_202[k];

        t_281[k] = f_5 * dpf0_138[k]
                   - f_6 * dpf1_138[k]
                   + f_4 * pc_x[k] * dpg_203[k];

        t_282[k] = f_5 * dpf0_139[k]
                   - f_6 * dpf1_139[k]
                   + f_4 * pc_x[k] * dpg_204[k];

        t_283[k] = f_4 * pc_x[k] * dpg_205[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_z, pc_x, pc_z, pph0_99, \
                         pph1_99, dpg_206, dpg_207, dpg_208, dpg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_4 * pc_x[k] * dpg_206[k];

        t_285[k] = f_4 * pc_x[k] * dpg_207[k];

        t_286[k] = f_4 * pc_x[k] * dpg_208[k];

        t_287[k] = f_4 * pc_x[k] * dpg_209[k];

        t_288[k] = pa_z[k] * pph0_99[k]
                   - f_9 * pc_z[k] * pph1_99[k];
    }

#pragma omp simd aligned(t_289, t_290, pc_y, pc_z, ppg_70, ppg_117, dsg_72, dpf0_138, \
                         dpf1_138, dpg_205, dpg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_1 * ppg_70[k]
                   + f_4 * pc_z[k] * dpg_205[k];

        t_290[k] = f_1 * ppg_117[k]
                   + f_1 * dsg_72[k]
                   + f_7 * dpf0_138[k]
                   - f_8 * dpf1_138[k]
                   + f_4 * pc_y[k] * dpg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pc_y, pc_z, ppg_74, ppg_118, ppg_119, dsg_73, \
                         dsg_74, dpf0_139, dpf1_139, dpg_208, dpg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_1 * ppg_118[k]
                   + f_1 * dsg_73[k]
                   + f_5 * dpf0_139[k]
                   - f_6 * dpf1_139[k]
                   + f_4 * pc_y[k] * dpg_208[k];

        t_292[k] = f_1 * ppg_119[k]
                   + f_1 * dsg_74[k]
                   + f_4 * pc_y[k] * dpg_209[k];

        t_293[k] = f_1 * ppg_74[k]
                   + f_2 * dpf0_139[k]
                   - f_3 * dpf1_139[k]
                   + f_4 * pc_z[k] * dpg_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pa_y, pc_x, pc_y, pph0_168, pph0_170, pph1_168, \
                         pph1_170, dpf0_141, dpf1_141, dpg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pa_y[k] * pph0_168[k]
                   - f_9 * pc_y[k] * pph1_168[k];

        t_295[k] = f_12 * dpf0_141[k]
                   - f_13 * dpf1_141[k]
                   + f_4 * pc_x[k] * dpg_211[k];

        t_296[k] = pa_y[k] * pph0_170[k]
                   - f_9 * pc_y[k] * pph1_170[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pa_y, pc_x, pc_y, pph0_173, pph1_173, dpf0_143, \
                         dpf0_144, dpf1_143, dpf1_144, dpg_213, \
                         dpg_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_7 * dpf0_143[k]
                   - f_8 * dpf1_143[k]
                   + f_4 * pc_x[k] * dpg_213[k];

        t_298[k] = f_7 * dpf0_144[k]
                   - f_8 * dpf1_144[k]
                   + f_4 * pc_x[k] * dpg_214[k];

        t_299[k] = pa_y[k] * pph0_173[k]
                   - f_9 * pc_y[k] * pph1_173[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pc_x, dpf0_146, dpf0_147, dpf0_148, dpf1_146, \
                         dpf1_147, dpf1_148, dpg_216, dpg_217, \
                         dpg_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_5 * dpf0_146[k]
                   - f_6 * dpf1_146[k]
                   + f_4 * pc_x[k] * dpg_216[k];

        t_301[k] = f_5 * dpf0_147[k]
                   - f_6 * dpf1_147[k]
                   + f_4 * pc_x[k] * dpg_217[k];

        t_302[k] = f_5 * dpf0_148[k]
                   - f_6 * dpf1_148[k]
                   + f_4 * pc_x[k] * dpg_218[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, t_308, pa_y, pc_x, pc_y, pph0_177, \
                         pph1_177, dpg_220, dpg_221, dpg_222, dpg_223, \
                         dpg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pa_y[k] * pph0_177[k]
                   - f_9 * pc_y[k] * pph1_177[k];

        t_304[k] = f_4 * pc_x[k] * dpg_220[k];

        t_305[k] = f_4 * pc_x[k] * dpg_221[k];

        t_306[k] = f_4 * pc_x[k] * dpg_222[k];

        t_307[k] = f_4 * pc_x[k] * dpg_223[k];

        t_308[k] = f_4 * pc_x[k] * dpg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pa_y, pc_y, pc_z, pph0_185, ppg_85, ppg_130, \
                         ppg_132, pph1_185, dsg_70, dpf0_146, dpf1_146, \
                         dpg_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * ppg_130[k]
                   + f_2 * dpf0_146[k]
                   - f_3 * dpf1_146[k]
                   + f_4 * pc_y[k] * dpg_220[k];

        t_310[k] = f_1 * ppg_85[k]
                   + f_1 * dsg_70[k]
                   + f_4 * pc_z[k] * dpg_220[k];

        t_311[k] = pa_y[k] * pph0_185[k]
                   + f_10 * ppg_132[k]
                   - f_9 * pc_y[k] * pph1_185[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pa_y, pc_y, pph0_186, pph0_188, ppg_133, \
                         ppg_134, pph1_186, pph1_188, dpg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = pa_y[k] * pph0_186[k]
                   + f_0 * ppg_133[k]
                   - f_9 * pc_y[k] * pph1_186[k];

        t_313[k] = f_1 * ppg_134[k]
                   + f_4 * pc_y[k] * dpg_224[k];

        t_314[k] = pa_y[k] * pph0_188[k]
                   - f_9 * pc_y[k] * pph1_188[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pb_x, pc_x, pc_y, dsh0_105, dsh0_107, dsg_75, \
                         dsg_77, dsh1_105, dsh1_107, dpg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pb_x[k] * dsh0_105[k]
                   + f_11 * dsg_75[k]
                   - f_9 * pc_x[k] * dsh1_105[k];

        t_316[k] = f_4 * pc_y[k] * dpg_225[k];

        t_317[k] = pb_x[k] * dsh0_107[k]
                   + f_14 * dsg_77[k]
                   - f_9 * pc_x[k] * dsh1_107[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pb_x, pc_x, pc_y, dsh0_108, dsh0_110, dsg_78, \
                         dsg_80, dsh1_108, dsh1_110, dpg_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_x[k] * dsh0_108[k]
                   + f_10 * dsg_78[k]
                   - f_9 * pc_x[k] * dsh1_108[k];

        t_319[k] = f_4 * pc_y[k] * dpg_227[k];

        t_320[k] = pb_x[k] * dsh0_110[k]
                   + f_10 * dsg_80[k]
                   - f_9 * pc_x[k] * dsh1_110[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pb_x, pc_x, pc_y, dsh0_111, dsh0_112, dsg_81, \
                         dsg_82, dsh1_111, dsh1_112, dpg_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = pb_x[k] * dsh0_111[k]
                   + f_0 * dsg_81[k]
                   - f_9 * pc_x[k] * dsh1_111[k];

        t_322[k] = pb_x[k] * dsh0_112[k]
                   + f_0 * dsg_82[k]
                   - f_9 * pc_x[k] * dsh1_112[k];

        t_323[k] = f_4 * pc_y[k] * dpg_230[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pb_x, pc_x, dsh0_114, dsg_84, dsg_85, \
                         dsg_86, dsg_87, dsh1_114, dpg_235, dpg_236, \
                         dpg_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = pb_x[k] * dsh0_114[k]
                   + f_0 * dsg_84[k]
                   - f_9 * pc_x[k] * dsh1_114[k];

        t_325[k] = f_1 * dsg_85[k]
                   + f_4 * pc_x[k] * dpg_235[k];

        t_326[k] = f_1 * dsg_86[k]
                   + f_4 * pc_x[k] * dpg_236[k];

        t_327[k] = f_1 * dsg_87[k]
                   + f_4 * pc_x[k] * dpg_237[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pb_x, pc_x, dsh0_120, dsh0_121, dsg_88, \
                         dsg_89, dsh1_120, dsh1_121, dpg_238, dpg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_1 * dsg_88[k]
                   + f_4 * pc_x[k] * dpg_238[k];

        t_329[k] = f_1 * dsg_89[k]
                   + f_4 * pc_x[k] * dpg_239[k];

        t_330[k] = pb_x[k] * dsh0_120[k]
                   - f_9 * pc_x[k] * dsh1_120[k];

        t_331[k] = pb_x[k] * dsh0_121[k]
                   - f_9 * pc_x[k] * dsh1_121[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pb_x, pc_x, pc_y, dsh0_122, dsh0_123, \
                         dsh0_125, dsh1_122, dsh1_123, dsh1_125, \
                         dpg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_x[k] * dsh0_122[k]
                   - f_9 * pc_x[k] * dsh1_122[k];

        t_333[k] = pb_x[k] * dsh0_123[k]
                   - f_9 * pc_x[k] * dsh1_123[k];

        t_334[k] = f_4 * pc_y[k] * dpg_239[k];

        t_335[k] = pb_x[k] * dsh0_125[k]
                   - f_9 * pc_x[k] * dsh1_125[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pb_y, pc_x, pc_y, dsh0_105, dsh0_107, \
                         dsg_75, dsh1_105, dsh1_107, dpf0_163, dpf1_163, dpg_240, \
                         dpg_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pb_y[k] * dsh0_105[k]
                   - f_9 * pc_y[k] * dsh1_105[k];

        t_337[k] = f_1 * dsg_75[k]
                   + f_4 * pc_y[k] * dpg_240[k];

        t_338[k] = pb_y[k] * dsh0_107[k]
                   - f_9 * pc_y[k] * dsh1_107[k];

        t_339[k] = f_7 * dpf0_163[k]
                   - f_8 * dpf1_163[k]
                   + f_4 * pc_x[k] * dpg_243[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pb_y, pc_x, pc_y, dsh0_110, dsg_77, dsh1_110, \
                         dpf0_166, dpf1_166, dpg_242, dpg_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_1 * dsg_77[k]
                   + f_4 * pc_y[k] * dpg_242[k];

        t_341[k] = pb_y[k] * dsh0_110[k]
                   - f_9 * pc_y[k] * dsh1_110[k];

        t_342[k] = f_5 * dpf0_166[k]
                   - f_6 * dpf1_166[k]
                   + f_4 * pc_x[k] * dpg_246[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pb_y, pc_x, pc_y, dsh0_114, dsg_80, \
                         dsh1_114, dpf0_167, dpf1_167, dpg_245, dpg_247, \
                         dpg_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_5 * dpf0_167[k]
                   - f_6 * dpf1_167[k]
                   + f_4 * pc_x[k] * dpg_247[k];

        t_344[k] = f_1 * dsg_80[k]
                   + f_4 * pc_y[k] * dpg_245[k];

        t_345[k] = pb_y[k] * dsh0_114[k]
                   - f_9 * pc_y[k] * dsh1_114[k];

        t_346[k] = f_4 * pc_x[k] * dpg_250[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pb_y, pc_x, pc_y, dsh0_120, \
                         dsg_85, dsh1_120, dpg_251, dpg_252, dpg_253, \
                         dpg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_4 * pc_x[k] * dpg_251[k];

        t_348[k] = f_4 * pc_x[k] * dpg_252[k];

        t_349[k] = f_4 * pc_x[k] * dpg_253[k];

        t_350[k] = f_4 * pc_x[k] * dpg_254[k];

        t_351[k] = pb_y[k] * dsh0_120[k]
                   + f_11 * dsg_85[k]
                   - f_9 * pc_y[k] * dsh1_120[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pb_y, pc_y, dsh0_121, dsh0_122, dsh0_123, \
                         dsg_86, dsg_87, dsg_88, dsh1_121, dsh1_122, \
                         dsh1_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = pb_y[k] * dsh0_121[k]
                   + f_14 * dsg_86[k]
                   - f_9 * pc_y[k] * dsh1_121[k];

        t_353[k] = pb_y[k] * dsh0_122[k]
                   + f_10 * dsg_87[k]
                   - f_9 * pc_y[k] * dsh1_122[k];

        t_354[k] = pb_y[k] * dsh0_123[k]
                   + f_0 * dsg_88[k]
                   - f_9 * pc_y[k] * dsh1_123[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pb_y, pc_x, pc_y, dsh0_125, dsg_89, \
                         dsh1_125, dpf0_170, dpf1_170, dpg_254, \
                         dpg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_1 * dsg_89[k]
                   + f_4 * pc_y[k] * dpg_254[k];

        t_356[k] = pb_y[k] * dsh0_125[k]
                   - f_9 * pc_y[k] * dsh1_125[k];

        t_357[k] = f_2 * dpf0_170[k]
                   - f_3 * dpf1_170[k]
                   + f_4 * pc_x[k] * dpg_255[k];

        t_358[k] = f_4 * pc_y[k] * dpg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pc_x, pc_y, dpf0_172, dpf0_173, dpf0_175, \
                         dpf1_172, dpf1_173, dpf1_175, dpg_257, dpg_258, \
                         dpg_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * dpf0_172[k]
                   - f_13 * dpf1_172[k]
                   + f_4 * pc_x[k] * dpg_257[k];

        t_360[k] = f_7 * dpf0_173[k]
                   - f_8 * dpf1_173[k]
                   + f_4 * pc_x[k] * dpg_258[k];

        t_361[k] = f_4 * pc_y[k] * dpg_257[k];

        t_362[k] = f_7 * dpf0_175[k]
                   - f_8 * dpf1_175[k]
                   + f_4 * pc_x[k] * dpg_260[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pc_x, pc_y, dpf0_176, dpf0_177, dpf0_179, \
                         dpf1_176, dpf1_177, dpf1_179, dpg_260, dpg_261, dpg_262, \
                         dpg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_5 * dpf0_176[k]
                   - f_6 * dpf1_176[k]
                   + f_4 * pc_x[k] * dpg_261[k];

        t_364[k] = f_5 * dpf0_177[k]
                   - f_6 * dpf1_177[k]
                   + f_4 * pc_x[k] * dpg_262[k];

        t_365[k] = f_4 * pc_y[k] * dpg_260[k];

        t_366[k] = f_5 * dpf0_179[k]
                   - f_6 * dpf1_179[k]
                   + f_4 * pc_x[k] * dpg_264[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, t_372, pc_x, pc_y, dpf0_176, \
                         dpf1_176, dpg_265, dpg_266, dpg_267, dpg_268, \
                         dpg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_4 * pc_x[k] * dpg_265[k];

        t_368[k] = f_4 * pc_x[k] * dpg_266[k];

        t_369[k] = f_4 * pc_x[k] * dpg_267[k];

        t_370[k] = f_4 * pc_x[k] * dpg_268[k];

        t_371[k] = f_4 * pc_x[k] * dpg_269[k];

        t_372[k] = f_2 * dpf0_176[k]
                   - f_3 * dpf1_176[k]
                   + f_4 * pc_y[k] * dpg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pc_y, dpf0_177, dpf0_178, dpf0_179, \
                         dpf1_177, dpf1_178, dpf1_179, dpg_266, dpg_267, dpg_268, \
                         dpg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_12 * dpf0_177[k]
                   - f_13 * dpf1_177[k]
                   + f_4 * pc_y[k] * dpg_266[k];

        t_374[k] = f_7 * dpf0_178[k]
                   - f_8 * dpf1_178[k]
                   + f_4 * pc_y[k] * dpg_267[k];

        t_375[k] = f_5 * dpf0_179[k]
                   - f_6 * dpf1_179[k]
                   + f_4 * pc_y[k] * dpg_268[k];

        t_376[k] = f_4 * pc_y[k] * dpg_269[k];
    }

#pragma omp simd aligned(t_377, pc_z, ppg_134, dsg_89, dpf0_179, dpf1_179, \
                         dpg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_0 * ppg_134[k]
                   + f_1 * dsg_89[k]
                   + f_2 * dpf0_179[k]
                   - f_3 * dpf1_179[k]
                   + f_4 * pc_z[k] * dpg_269[k];
    }
}

auto
compute_prim_dph_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t pph0,
                                                   const size_t ppg, const size_t pph1,
                                                   const size_t dsh0, const size_t dsg,
                                                   const size_t dsh1, const size_t dpf0,
                                                   const size_t dpf1, const size_t dpg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_dph_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, pph0,
                                                              ppg, pph1, dsh0, dsg, dsh1, dpf0,
                                                              dpf1, dpg, ncols, gamma, p, q);

    compute_prim_dph_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, pph0,
                                                              ppg, pph1, dsh0, dsg, dsh1, dpf0,
                                                              dpf1, dpg, ncols, gamma, p, q);

    compute_prim_dph_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, pph0,
                                                              ppg, pph1, dsh0, dsg, dsh1, dpf0,
                                                              dpf1, dpg, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
