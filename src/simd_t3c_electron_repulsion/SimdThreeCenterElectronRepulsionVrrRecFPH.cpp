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


#include "SimdThreeCenterElectronRepulsionVrrRecFPH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fph_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t pph0, const size_t pph1,
                                                          const size_t dph0, const size_t dpg,
                                                          const size_t dph1, const size_t fsh0,
                                                          const size_t fsg, const size_t fsh1,
                                                          const size_t fpf0, const size_t fpf1,
                                                          const size_t fpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);

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

    const auto *pph0_99 = buffer.data(pph0 + 99);

    const auto *pph1_99 = buffer.data(pph1 + 99);

    const auto *dph0_0 = buffer.data(dph0 + 0);
    const auto *dph0_3 = buffer.data(dph0 + 3);
    const auto *dph0_5 = buffer.data(dph0 + 5);
    const auto *dph0_6 = buffer.data(dph0 + 6);
    const auto *dph0_9 = buffer.data(dph0 + 9);
    const auto *dph0_14 = buffer.data(dph0 + 14);
    const auto *dph0_20 = buffer.data(dph0 + 20);
    const auto *dph0_42 = buffer.data(dph0 + 42);
    const auto *dph0_47 = buffer.data(dph0 + 47);
    const auto *dph0_51 = buffer.data(dph0 + 51);
    const auto *dph0_62 = buffer.data(dph0 + 62);
    const auto *dph0_99 = buffer.data(dph0 + 99);

    const auto *dpg_0 = buffer.data(dpg + 0);
    const auto *dpg_1 = buffer.data(dpg + 1);
    const auto *dpg_3 = buffer.data(dpg + 3);
    const auto *dpg_5 = buffer.data(dpg + 5);
    const auto *dpg_10 = buffer.data(dpg + 10);
    const auto *dpg_12 = buffer.data(dpg + 12);
    const auto *dpg_14 = buffer.data(dpg + 14);
    const auto *dpg_15 = buffer.data(dpg + 15);
    const auto *dpg_20 = buffer.data(dpg + 20);
    const auto *dpg_25 = buffer.data(dpg + 25);
    const auto *dpg_27 = buffer.data(dpg + 27);
    const auto *dpg_29 = buffer.data(dpg + 29);
    const auto *dpg_30 = buffer.data(dpg + 30);
    const auto *dpg_35 = buffer.data(dpg + 35);
    const auto *dpg_40 = buffer.data(dpg + 40);
    const auto *dpg_42 = buffer.data(dpg + 42);
    const auto *dpg_44 = buffer.data(dpg + 44);
    const auto *dpg_55 = buffer.data(dpg + 55);
    const auto *dpg_57 = buffer.data(dpg + 57);
    const auto *dpg_58 = buffer.data(dpg + 58);
    const auto *dpg_60 = buffer.data(dpg + 60);
    const auto *dpg_63 = buffer.data(dpg + 63);
    const auto *dpg_66 = buffer.data(dpg + 66);
    const auto *dpg_70 = buffer.data(dpg + 70);
    const auto *dpg_72 = buffer.data(dpg + 72);
    const auto *dpg_73 = buffer.data(dpg + 73);
    const auto *dpg_74 = buffer.data(dpg + 74);
    const auto *dpg_85 = buffer.data(dpg + 85);
    const auto *dpg_87 = buffer.data(dpg + 87);
    const auto *dpg_88 = buffer.data(dpg + 88);
    const auto *dpg_89 = buffer.data(dpg + 89);

    const auto *dph1_0 = buffer.data(dph1 + 0);
    const auto *dph1_3 = buffer.data(dph1 + 3);
    const auto *dph1_5 = buffer.data(dph1 + 5);
    const auto *dph1_6 = buffer.data(dph1 + 6);
    const auto *dph1_9 = buffer.data(dph1 + 9);
    const auto *dph1_14 = buffer.data(dph1 + 14);
    const auto *dph1_20 = buffer.data(dph1 + 20);
    const auto *dph1_42 = buffer.data(dph1 + 42);
    const auto *dph1_47 = buffer.data(dph1 + 47);
    const auto *dph1_51 = buffer.data(dph1 + 51);
    const auto *dph1_62 = buffer.data(dph1 + 62);
    const auto *dph1_99 = buffer.data(dph1 + 99);

    const auto *fsh0_0 = buffer.data(fsh0 + 0);
    const auto *fsh0_3 = buffer.data(fsh0 + 3);
    const auto *fsh0_5 = buffer.data(fsh0 + 5);
    const auto *fsh0_6 = buffer.data(fsh0 + 6);
    const auto *fsh0_9 = buffer.data(fsh0 + 9);
    const auto *fsh0_15 = buffer.data(fsh0 + 15);
    const auto *fsh0_17 = buffer.data(fsh0 + 17);
    const auto *fsh0_18 = buffer.data(fsh0 + 18);
    const auto *fsh0_20 = buffer.data(fsh0 + 20);
    const auto *fsh0_24 = buffer.data(fsh0 + 24);
    const auto *fsh0_27 = buffer.data(fsh0 + 27);
    const auto *fsh0_36 = buffer.data(fsh0 + 36);
    const auto *fsh0_38 = buffer.data(fsh0 + 38);
    const auto *fsh0_39 = buffer.data(fsh0 + 39);

    const auto *fsg_0 = buffer.data(fsg + 0);
    const auto *fsg_1 = buffer.data(fsg + 1);
    const auto *fsg_2 = buffer.data(fsg + 2);
    const auto *fsg_3 = buffer.data(fsg + 3);
    const auto *fsg_5 = buffer.data(fsg + 5);
    const auto *fsg_6 = buffer.data(fsg + 6);
    const auto *fsg_9 = buffer.data(fsg + 9);
    const auto *fsg_10 = buffer.data(fsg + 10);
    const auto *fsg_12 = buffer.data(fsg + 12);
    const auto *fsg_13 = buffer.data(fsg + 13);
    const auto *fsg_14 = buffer.data(fsg + 14);
    const auto *fsg_15 = buffer.data(fsg + 15);
    const auto *fsg_16 = buffer.data(fsg + 16);
    const auto *fsg_18 = buffer.data(fsg + 18);
    const auto *fsg_20 = buffer.data(fsg + 20);
    const auto *fsg_21 = buffer.data(fsg + 21);
    const auto *fsg_25 = buffer.data(fsg + 25);
    const auto *fsg_26 = buffer.data(fsg + 26);
    const auto *fsg_27 = buffer.data(fsg + 27);
    const auto *fsg_28 = buffer.data(fsg + 28);
    const auto *fsg_29 = buffer.data(fsg + 29);

    const auto *fsh1_0 = buffer.data(fsh1 + 0);
    const auto *fsh1_3 = buffer.data(fsh1 + 3);
    const auto *fsh1_5 = buffer.data(fsh1 + 5);
    const auto *fsh1_6 = buffer.data(fsh1 + 6);
    const auto *fsh1_9 = buffer.data(fsh1 + 9);
    const auto *fsh1_15 = buffer.data(fsh1 + 15);
    const auto *fsh1_17 = buffer.data(fsh1 + 17);
    const auto *fsh1_18 = buffer.data(fsh1 + 18);
    const auto *fsh1_20 = buffer.data(fsh1 + 20);
    const auto *fsh1_24 = buffer.data(fsh1 + 24);
    const auto *fsh1_27 = buffer.data(fsh1 + 27);
    const auto *fsh1_36 = buffer.data(fsh1 + 36);
    const auto *fsh1_38 = buffer.data(fsh1 + 38);
    const auto *fsh1_39 = buffer.data(fsh1 + 39);

    const auto *fpf0_0 = buffer.data(fpf0 + 0);
    const auto *fpf0_1 = buffer.data(fpf0 + 1);
    const auto *fpf0_2 = buffer.data(fpf0 + 2);
    const auto *fpf0_6 = buffer.data(fpf0 + 6);
    const auto *fpf0_8 = buffer.data(fpf0 + 8);
    const auto *fpf0_9 = buffer.data(fpf0 + 9);
    const auto *fpf0_28 = buffer.data(fpf0 + 28);
    const auto *fpf0_29 = buffer.data(fpf0 + 29);
    const auto *fpf0_36 = buffer.data(fpf0 + 36);
    const auto *fpf0_37 = buffer.data(fpf0 + 37);
    const auto *fpf0_40 = buffer.data(fpf0 + 40);
    const auto *fpf0_42 = buffer.data(fpf0 + 42);
    const auto *fpf0_43 = buffer.data(fpf0 + 43);
    const auto *fpf0_46 = buffer.data(fpf0 + 46);
    const auto *fpf0_47 = buffer.data(fpf0 + 47);
    const auto *fpf0_49 = buffer.data(fpf0 + 49);

    const auto *fpf1_0 = buffer.data(fpf1 + 0);
    const auto *fpf1_1 = buffer.data(fpf1 + 1);
    const auto *fpf1_2 = buffer.data(fpf1 + 2);
    const auto *fpf1_6 = buffer.data(fpf1 + 6);
    const auto *fpf1_8 = buffer.data(fpf1 + 8);
    const auto *fpf1_9 = buffer.data(fpf1 + 9);
    const auto *fpf1_28 = buffer.data(fpf1 + 28);
    const auto *fpf1_29 = buffer.data(fpf1 + 29);
    const auto *fpf1_36 = buffer.data(fpf1 + 36);
    const auto *fpf1_37 = buffer.data(fpf1 + 37);
    const auto *fpf1_40 = buffer.data(fpf1 + 40);
    const auto *fpf1_42 = buffer.data(fpf1 + 42);
    const auto *fpf1_43 = buffer.data(fpf1 + 43);
    const auto *fpf1_46 = buffer.data(fpf1 + 46);
    const auto *fpf1_47 = buffer.data(fpf1 + 47);
    const auto *fpf1_49 = buffer.data(fpf1 + 49);

    const auto *fpg_0 = buffer.data(fpg + 0);
    const auto *fpg_1 = buffer.data(fpg + 1);
    const auto *fpg_2 = buffer.data(fpg + 2);
    const auto *fpg_3 = buffer.data(fpg + 3);
    const auto *fpg_5 = buffer.data(fpg + 5);
    const auto *fpg_6 = buffer.data(fpg + 6);
    const auto *fpg_9 = buffer.data(fpg + 9);
    const auto *fpg_10 = buffer.data(fpg + 10);
    const auto *fpg_12 = buffer.data(fpg + 12);
    const auto *fpg_13 = buffer.data(fpg + 13);
    const auto *fpg_14 = buffer.data(fpg + 14);
    const auto *fpg_15 = buffer.data(fpg + 15);
    const auto *fpg_17 = buffer.data(fpg + 17);
    const auto *fpg_18 = buffer.data(fpg + 18);
    const auto *fpg_20 = buffer.data(fpg + 20);
    const auto *fpg_21 = buffer.data(fpg + 21);
    const auto *fpg_24 = buffer.data(fpg + 24);
    const auto *fpg_25 = buffer.data(fpg + 25);
    const auto *fpg_27 = buffer.data(fpg + 27);
    const auto *fpg_29 = buffer.data(fpg + 29);
    const auto *fpg_30 = buffer.data(fpg + 30);
    const auto *fpg_32 = buffer.data(fpg + 32);
    const auto *fpg_33 = buffer.data(fpg + 33);
    const auto *fpg_35 = buffer.data(fpg + 35);
    const auto *fpg_36 = buffer.data(fpg + 36);
    const auto *fpg_39 = buffer.data(fpg + 39);
    const auto *fpg_40 = buffer.data(fpg + 40);
    const auto *fpg_42 = buffer.data(fpg + 42);
    const auto *fpg_43 = buffer.data(fpg + 43);
    const auto *fpg_44 = buffer.data(fpg + 44);
    const auto *fpg_45 = buffer.data(fpg + 45);
    const auto *fpg_46 = buffer.data(fpg + 46);
    const auto *fpg_48 = buffer.data(fpg + 48);
    const auto *fpg_50 = buffer.data(fpg + 50);
    const auto *fpg_51 = buffer.data(fpg + 51);
    const auto *fpg_55 = buffer.data(fpg + 55);
    const auto *fpg_56 = buffer.data(fpg + 56);
    const auto *fpg_57 = buffer.data(fpg + 57);
    const auto *fpg_58 = buffer.data(fpg + 58);
    const auto *fpg_59 = buffer.data(fpg + 59);
    const auto *fpg_60 = buffer.data(fpg + 60);
    const auto *fpg_61 = buffer.data(fpg + 61);
    const auto *fpg_62 = buffer.data(fpg + 62);
    const auto *fpg_63 = buffer.data(fpg + 63);
    const auto *fpg_65 = buffer.data(fpg + 65);
    const auto *fpg_66 = buffer.data(fpg + 66);
    const auto *fpg_70 = buffer.data(fpg + 70);
    const auto *fpg_71 = buffer.data(fpg + 71);
    const auto *fpg_72 = buffer.data(fpg + 72);
    const auto *fpg_73 = buffer.data(fpg + 73);
    const auto *fpg_74 = buffer.data(fpg + 74);
    const auto *fpg_75 = buffer.data(fpg + 75);
    const auto *fpg_76 = buffer.data(fpg + 76);
    const auto *fpg_78 = buffer.data(fpg + 78);
    const auto *fpg_80 = buffer.data(fpg + 80);
    const auto *fpg_81 = buffer.data(fpg + 81);
    const auto *fpg_85 = buffer.data(fpg + 85);
    const auto *fpg_87 = buffer.data(fpg + 87);
    const auto *fpg_88 = buffer.data(fpg + 88);
    const auto *fpg_89 = buffer.data(fpg + 89);
    const auto *fpg_90 = buffer.data(fpg + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dpg_0, fsg_0, fpf0_0, \
                         fpf1_0, fpg_0, fpg_1, fpg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dpg_0[k]
                 + f_1 * fsg_0[k]
                 + f_2 * fpf0_0[k]
                 - f_3 * fpf1_0[k]
                 + f_4 * pc_x[k] * fpg_0[k];

        t_1[k] = f_4 * pc_y[k] * fpg_0[k];

        t_2[k] = f_4 * pc_z[k] * fpg_0[k];

        t_3[k] = f_5 * fpf0_0[k]
                 - f_6 * fpf1_0[k]
                 + f_4 * pc_y[k] * fpg_1[k];

        t_4[k] = f_4 * pc_y[k] * fpg_2[k];

        t_5[k] = f_5 * fpf0_0[k]
                 - f_6 * fpf1_0[k]
                 + f_4 * pc_z[k] * fpg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_y, pc_z, fpf0_1, fpf0_2, fpf1_1, fpf1_2, \
                         fpg_3, fpg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * fpf0_1[k]
                 - f_8 * fpf1_1[k]
                 + f_4 * pc_y[k] * fpg_3[k];

        t_7[k] = f_4 * pc_z[k] * fpg_3[k];

        t_8[k] = f_4 * pc_y[k] * fpg_5[k];

        t_9[k] = f_7 * fpf0_2[k]
                 - f_8 * fpf1_2[k]
                 + f_4 * pc_z[k] * fpg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pc_x, pc_y, pc_z, dpg_10, dpg_12, fsg_10, \
                         fsg_12, fpg_6, fpg_9, fpg_10, fpg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * dpg_10[k]
                  + f_1 * fsg_10[k]
                  + f_4 * pc_x[k] * fpg_10[k];

        t_11[k] = f_4 * pc_z[k] * fpg_6[k];

        t_12[k] = f_0 * dpg_12[k]
                  + f_1 * fsg_12[k]
                  + f_4 * pc_x[k] * fpg_12[k];

        t_13[k] = f_4 * pc_y[k] * fpg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, dpg_14, fsg_14, fpf0_6, \
                         fpf0_8, fpf1_6, fpf1_8, fpg_10, fpg_12, \
                         fpg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * dpg_14[k]
                  + f_1 * fsg_14[k]
                  + f_4 * pc_x[k] * fpg_14[k];

        t_15[k] = f_2 * fpf0_6[k]
                  - f_3 * fpf1_6[k]
                  + f_4 * pc_y[k] * fpg_10[k];

        t_16[k] = f_4 * pc_z[k] * fpg_10[k];

        t_17[k] = f_7 * fpf0_8[k]
                  - f_8 * fpf1_8[k]
                  + f_4 * pc_y[k] * fpg_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pb_y, pc_y, pc_z, fsh0_0, fsg_0, \
                         fsh1_0, fpf0_9, fpf1_9, fpg_13, fpg_14, \
                         fpg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * fpf0_9[k]
                  - f_6 * fpf1_9[k]
                  + f_4 * pc_y[k] * fpg_13[k];

        t_19[k] = f_4 * pc_y[k] * fpg_14[k];

        t_20[k] = f_2 * fpf0_9[k]
                  - f_3 * fpf1_9[k]
                  + f_4 * pc_z[k] * fpg_14[k];

        t_21[k] = pb_y[k] * fsh0_0[k]
                  - f_9 * pc_y[k] * fsh1_0[k];

        t_22[k] = f_1 * fsg_0[k]
                  + f_4 * pc_y[k] * fpg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_y, pc_z, fsh0_3, fsh0_5, fsg_1, \
                         fsg_2, fsh1_3, fsh1_5, fpg_15, fpg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * pc_z[k] * fpg_15[k];

        t_24[k] = pb_y[k] * fsh0_3[k]
                  + f_10 * fsg_1[k]
                  - f_9 * pc_y[k] * fsh1_3[k];

        t_25[k] = f_1 * fsg_2[k]
                  + f_4 * pc_y[k] * fpg_17[k];

        t_26[k] = pb_y[k] * fsh0_5[k]
                  - f_9 * pc_y[k] * fsh1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_y, pc_y, pc_z, fsh0_6, fsh0_9, fsg_3, \
                         fsg_5, fsh1_6, fsh1_9, fpg_18, fpg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * fsh0_6[k]
                  + f_0 * fsg_3[k]
                  - f_9 * pc_y[k] * fsh1_6[k];

        t_28[k] = f_4 * pc_z[k] * fpg_18[k];

        t_29[k] = f_1 * fsg_5[k]
                  + f_4 * pc_y[k] * fpg_20[k];

        t_30[k] = pb_y[k] * fsh0_9[k]
                  - f_9 * pc_y[k] * fsh1_9[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pc_x, pc_y, pc_z, dpg_25, dpg_27, fsg_9, \
                         fpg_21, fpg_24, fpg_25, fpg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * dpg_25[k]
                  + f_4 * pc_x[k] * fpg_25[k];

        t_32[k] = f_4 * pc_z[k] * fpg_21[k];

        t_33[k] = f_0 * dpg_27[k]
                  + f_4 * pc_x[k] * fpg_27[k];

        t_34[k] = f_1 * fsg_9[k]
                  + f_4 * pc_y[k] * fpg_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pc_x, pc_y, pc_z, dpg_29, fsh0_15, fsg_10, \
                         fsh1_15, fpg_25, fpg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * dpg_29[k]
                  + f_4 * pc_x[k] * fpg_29[k];

        t_36[k] = pb_y[k] * fsh0_15[k]
                  + f_11 * fsg_10[k]
                  - f_9 * pc_y[k] * fsh1_15[k];

        t_37[k] = f_4 * pc_z[k] * fpg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, fsh0_17, fsh0_18, fsh0_20, \
                         fsg_12, fsg_13, fsg_14, fsh1_17, fsh1_18, fsh1_20, \
                         fpg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * fsh0_17[k]
                  + f_0 * fsg_12[k]
                  - f_9 * pc_y[k] * fsh1_17[k];

        t_39[k] = pb_y[k] * fsh0_18[k]
                  + f_10 * fsg_13[k]
                  - f_9 * pc_y[k] * fsh1_18[k];

        t_40[k] = f_1 * fsg_14[k]
                  + f_4 * pc_y[k] * fpg_29[k];

        t_41[k] = pb_y[k] * fsh0_20[k]
                  - f_9 * pc_y[k] * fsh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, fsh0_0, fsh0_3, \
                         fsg_0, fsh1_0, fsh1_3, fpg_30, fpg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * fsh0_0[k]
                  - f_9 * pc_z[k] * fsh1_0[k];

        t_43[k] = f_4 * pc_y[k] * fpg_30[k];

        t_44[k] = f_1 * fsg_0[k]
                  + f_4 * pc_z[k] * fpg_30[k];

        t_45[k] = pb_z[k] * fsh0_3[k]
                  - f_9 * pc_z[k] * fsh1_3[k];

        t_46[k] = f_4 * pc_y[k] * fpg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, fsh0_5, fsh0_6, fsg_2, \
                         fsg_3, fsh1_5, fsh1_6, fpg_33, fpg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * fsh0_5[k]
                  + f_10 * fsg_2[k]
                  - f_9 * pc_z[k] * fsh1_5[k];

        t_48[k] = pb_z[k] * fsh0_6[k]
                  - f_9 * pc_z[k] * fsh1_6[k];

        t_49[k] = f_1 * fsg_3[k]
                  + f_4 * pc_z[k] * fpg_33[k];

        t_50[k] = f_4 * pc_y[k] * fpg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, dpg_40, dpg_42, fsh0_9, \
                         fsg_5, fsg_6, fsh1_9, fpg_36, fpg_40, fpg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * fsh0_9[k]
                  + f_0 * fsg_5[k]
                  - f_9 * pc_z[k] * fsh1_9[k];

        t_52[k] = f_0 * dpg_40[k]
                  + f_4 * pc_x[k] * fpg_40[k];

        t_53[k] = f_1 * fsg_6[k]
                  + f_4 * pc_z[k] * fpg_36[k];

        t_54[k] = f_0 * dpg_42[k]
                  + f_4 * pc_x[k] * fpg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_y, pc_z, dpg_44, fsh0_15, \
                         fsg_10, fsh1_15, fpg_39, fpg_40, fpg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_4 * pc_y[k] * fpg_39[k];

        t_56[k] = f_0 * dpg_44[k]
                  + f_4 * pc_x[k] * fpg_44[k];

        t_57[k] = pb_z[k] * fsh0_15[k]
                  - f_9 * pc_z[k] * fsh1_15[k];

        t_58[k] = f_1 * fsg_10[k]
                  + f_4 * pc_z[k] * fpg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pc_y, fpf0_28, fpf0_29, fpf1_28, fpf1_29, fpg_42, \
                         fpg_43, fpg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_7 * fpf0_28[k]
                  - f_8 * fpf1_28[k]
                  + f_4 * pc_y[k] * fpg_42[k];

        t_60[k] = f_5 * fpf0_29[k]
                  - f_6 * fpf1_29[k]
                  + f_4 * pc_y[k] * fpg_43[k];

        t_61[k] = f_4 * pc_y[k] * fpg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pb_z, pc_y, pc_z, dph0_0, dpg_0, \
                         dph1_0, fsh0_20, fsg_14, fsh1_20, fpg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_z[k] * fsh0_20[k]
                  + f_11 * fsg_14[k]
                  - f_9 * pc_z[k] * fsh1_20[k];

        t_63[k] = pa_y[k] * dph0_0[k]
                  - f_9 * pc_y[k] * dph1_0[k];

        t_64[k] = f_1 * dpg_0[k]
                  + f_4 * pc_y[k] * fpg_45[k];

        t_65[k] = f_4 * pc_z[k] * fpg_45[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_y, pc_y, pc_z, dph0_3, dph0_5, dph0_6, \
                         dpg_1, dpg_3, dph1_3, dph1_5, dph1_6, fpg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_y[k] * dph0_3[k]
                  + f_10 * dpg_1[k]
                  - f_9 * pc_y[k] * dph1_3[k];

        t_67[k] = f_4 * pc_z[k] * fpg_46[k];

        t_68[k] = pa_y[k] * dph0_5[k]
                  - f_9 * pc_y[k] * dph1_5[k];

        t_69[k] = pa_y[k] * dph0_6[k]
                  + f_0 * dpg_3[k]
                  - f_9 * pc_y[k] * dph1_6[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_y, pc_x, pc_y, pc_z, dph0_9, dpg_5, \
                         dpg_55, dph1_9, fsg_25, fpg_48, fpg_50, \
                         fpg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_4 * pc_z[k] * fpg_48[k];

        t_71[k] = f_1 * dpg_5[k]
                  + f_4 * pc_y[k] * fpg_50[k];

        t_72[k] = pa_y[k] * dph0_9[k]
                  - f_9 * pc_y[k] * dph1_9[k];

        t_73[k] = f_10 * dpg_55[k]
                  + f_1 * fsg_25[k]
                  + f_4 * pc_x[k] * fpg_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pc_x, pc_z, dpg_57, dpg_58, fsg_27, fsg_28, fpg_51, \
                         fpg_57, fpg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_4 * pc_z[k] * fpg_51[k];

        t_75[k] = f_10 * dpg_57[k]
                  + f_1 * fsg_27[k]
                  + f_4 * pc_x[k] * fpg_57[k];

        t_76[k] = f_10 * dpg_58[k]
                  + f_1 * fsg_28[k]
                  + f_4 * pc_x[k] * fpg_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, dph0_14, dpg_10, dph1_14, \
                         fpf0_36, fpf1_36, fpg_55, fpg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_y[k] * dph0_14[k]
                  - f_9 * pc_y[k] * dph1_14[k];

        t_78[k] = f_1 * dpg_10[k]
                  + f_2 * fpf0_36[k]
                  - f_3 * fpf1_36[k]
                  + f_4 * pc_y[k] * fpg_55[k];

        t_79[k] = f_4 * pc_z[k] * fpg_55[k];

        t_80[k] = f_5 * fpf0_36[k]
                  - f_6 * fpf1_36[k]
                  + f_4 * pc_z[k] * fpg_56[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_y, pc_y, pc_z, dph0_20, dpg_14, dph1_20, \
                         fpf0_37, fpf1_37, fpg_57, fpg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_7 * fpf0_37[k]
                  - f_8 * fpf1_37[k]
                  + f_4 * pc_z[k] * fpg_57[k];

        t_82[k] = f_1 * dpg_14[k]
                  + f_4 * pc_y[k] * fpg_59[k];

        t_83[k] = pa_y[k] * dph0_20[k]
                  - f_9 * pc_y[k] * dph1_20[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, pc_z, dpg_15, dpg_60, fsg_15, fpf0_40, \
                         fpf1_40, fpg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_10 * dpg_60[k]
                  + f_2 * fpf0_40[k]
                  - f_3 * fpf1_40[k]
                  + f_4 * pc_x[k] * fpg_60[k];

        t_85[k] = f_1 * dpg_15[k]
                  + f_1 * fsg_15[k]
                  + f_4 * pc_y[k] * fpg_60[k];

        t_86[k] = f_4 * pc_z[k] * fpg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pc_x, pc_z, dpg_63, fpf0_40, fpf0_43, fpf1_40, \
                         fpf1_43, fpg_61, fpg_62, fpg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_10 * dpg_63[k]
                  + f_7 * fpf0_43[k]
                  - f_8 * fpf1_43[k]
                  + f_4 * pc_x[k] * fpg_63[k];

        t_88[k] = f_4 * pc_z[k] * fpg_61[k];

        t_89[k] = f_5 * fpf0_40[k]
                  - f_6 * fpf1_40[k]
                  + f_4 * pc_z[k] * fpg_62[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pc_x, pc_y, pc_z, dpg_20, dpg_66, fsg_20, fpf0_46, \
                         fpf1_46, fpg_63, fpg_65, fpg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_10 * dpg_66[k]
                  + f_5 * fpf0_46[k]
                  - f_6 * fpf1_46[k]
                  + f_4 * pc_x[k] * fpg_66[k];

        t_91[k] = f_4 * pc_z[k] * fpg_63[k];

        t_92[k] = f_1 * dpg_20[k]
                  + f_1 * fsg_20[k]
                  + f_4 * pc_y[k] * fpg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_z, dpg_70, dpg_72, fpf0_42, fpf1_42, \
                         fpg_65, fpg_66, fpg_70, fpg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_7 * fpf0_42[k]
                  - f_8 * fpf1_42[k]
                  + f_4 * pc_z[k] * fpg_65[k];

        t_94[k] = f_10 * dpg_70[k]
                  + f_4 * pc_x[k] * fpg_70[k];

        t_95[k] = f_4 * pc_z[k] * fpg_66[k];

        t_96[k] = f_10 * dpg_72[k]
                  + f_4 * pc_x[k] * fpg_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_x, pc_x, pc_z, pph0_99, pph1_99, dph0_99, \
                         dpg_73, dpg_74, dph1_99, fpg_70, fpg_73, \
                         fpg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_10 * dpg_73[k]
                  + f_4 * pc_x[k] * fpg_73[k];

        t_98[k] = f_10 * dpg_74[k]
                  + f_4 * pc_x[k] * fpg_74[k];

        t_99[k] = f_12 * pph0_99[k]
                  - f_13 * pph1_99[k]
                  + pa_x[k] * dph0_99[k]
                  - f_9 * pc_x[k] * dph1_99[k];

        t_100[k] = f_4 * pc_z[k] * fpg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, pc_z, dpg_29, fsg_29, fpf0_46, fpf0_47, \
                         fpf1_46, fpf1_47, fpg_71, fpg_72, fpg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_5 * fpf0_46[k]
                   - f_6 * fpf1_46[k]
                   + f_4 * pc_z[k] * fpg_71[k];

        t_102[k] = f_7 * fpf0_47[k]
                   - f_8 * fpf1_47[k]
                   + f_4 * pc_z[k] * fpg_72[k];

        t_103[k] = f_1 * dpg_29[k]
                   + f_1 * fsg_29[k]
                   + f_4 * pc_y[k] * fpg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pc_y, pc_z, dph0_42, dpg_30, \
                         dph1_42, fsg_15, fpf0_49, fpf1_49, fpg_74, \
                         fpg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_2 * fpf0_49[k]
                   - f_3 * fpf1_49[k]
                   + f_4 * pc_z[k] * fpg_74[k];

        t_105[k] = pa_y[k] * dph0_42[k]
                   - f_9 * pc_y[k] * dph1_42[k];

        t_106[k] = f_1 * dpg_30[k]
                   + f_4 * pc_y[k] * fpg_75[k];

        t_107[k] = f_1 * fsg_15[k]
                   + f_4 * pc_z[k] * fpg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pa_y, pb_z, pc_y, pc_z, dph0_47, dph1_47, \
                         fsh0_24, fsh0_27, fsg_16, fsh1_24, fsh1_27, \
                         fpg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * fsh0_24[k]
                   - f_9 * pc_z[k] * fsh1_24[k];

        t_109[k] = f_1 * fsg_16[k]
                   + f_4 * pc_z[k] * fpg_76[k];

        t_110[k] = pa_y[k] * dph0_47[k]
                   - f_9 * pc_y[k] * dph1_47[k];

        t_111[k] = pb_z[k] * fsh0_27[k]
                   - f_9 * pc_z[k] * fsh1_27[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_y, pc_x, pc_y, pc_z, dph0_51, dpg_35, \
                         dpg_85, dph1_51, fsg_18, fpg_78, fpg_80, \
                         fpg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_1 * fsg_18[k]
                   + f_4 * pc_z[k] * fpg_78[k];

        t_113[k] = f_1 * dpg_35[k]
                   + f_4 * pc_y[k] * fpg_80[k];

        t_114[k] = pa_y[k] * dph0_51[k]
                   - f_9 * pc_y[k] * dph1_51[k];

        t_115[k] = f_10 * dpg_85[k]
                   + f_4 * pc_x[k] * fpg_85[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pc_x, pc_z, dpg_87, dpg_88, dpg_89, \
                         fsg_21, fpg_81, fpg_87, fpg_88, fpg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_1 * fsg_21[k]
                   + f_4 * pc_z[k] * fpg_81[k];

        t_117[k] = f_10 * dpg_87[k]
                   + f_4 * pc_x[k] * fpg_87[k];

        t_118[k] = f_10 * dpg_88[k]
                   + f_4 * pc_x[k] * fpg_88[k];

        t_119[k] = f_10 * dpg_89[k]
                   + f_4 * pc_x[k] * fpg_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_z, pc_z, fsh0_36, fsh0_38, fsh0_39, \
                         fsg_25, fsg_26, fsg_27, fsh1_36, fsh1_38, fsh1_39, \
                         fpg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pb_z[k] * fsh0_36[k]
                   - f_9 * pc_z[k] * fsh1_36[k];

        t_121[k] = f_1 * fsg_25[k]
                   + f_4 * pc_z[k] * fpg_85[k];

        t_122[k] = pb_z[k] * fsh0_38[k]
                   + f_10 * fsg_26[k]
                   - f_9 * pc_z[k] * fsh1_38[k];

        t_123[k] = pb_z[k] * fsh0_39[k]
                   + f_0 * fsg_27[k]
                   - f_9 * pc_z[k] * fsh1_39[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pa_z, pc_y, pc_z, dph0_0, dph0_62, \
                         dpg_44, dph1_0, dph1_62, fpg_89, fpg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * dpg_44[k]
                   + f_4 * pc_y[k] * fpg_89[k];

        t_125[k] = pa_y[k] * dph0_62[k]
                   - f_9 * pc_y[k] * dph1_62[k];

        t_126[k] = pa_z[k] * dph0_0[k]
                   - f_9 * pc_z[k] * dph1_0[k];

        t_127[k] = f_4 * pc_y[k] * fpg_90[k];
    }
}

static auto
compute_prim_fph_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t pph0, const size_t pph1,
                                                          const size_t dph0, const size_t dpg,
                                                          const size_t dph1, const size_t fsh0,
                                                          const size_t fsg, const size_t fsh1,
                                                          const size_t fpf0, const size_t fpf1,
                                                          const size_t fpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);
    const auto f_14 = 1.5 / gamma;
    const auto f_15 = 1.5 * p / (gamma * q);
    const auto f_16 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *pph0_188 = buffer.data(pph0 + 188);

    const auto *pph1_188 = buffer.data(pph1 + 188);

    const auto *dph0_3 = buffer.data(dph0 + 3);
    const auto *dph0_5 = buffer.data(dph0 + 5);
    const auto *dph0_6 = buffer.data(dph0 + 6);
    const auto *dph0_9 = buffer.data(dph0 + 9);
    const auto *dph0_10 = buffer.data(dph0 + 10);
    const auto *dph0_15 = buffer.data(dph0 + 15);
    const auto *dph0_21 = buffer.data(dph0 + 21);
    const auto *dph0_24 = buffer.data(dph0 + 24);
    const auto *dph0_27 = buffer.data(dph0 + 27);
    const auto *dph0_28 = buffer.data(dph0 + 28);
    const auto *dph0_36 = buffer.data(dph0 + 36);
    const auto *dph0_188 = buffer.data(dph0 + 188);
    const auto *dph0_210 = buffer.data(dph0 + 210);
    const auto *dph0_213 = buffer.data(dph0 + 213);
    const auto *dph0_216 = buffer.data(dph0 + 216);
    const auto *dph0_225 = buffer.data(dph0 + 225);
    const auto *dph0_227 = buffer.data(dph0 + 227);
    const auto *dph0_228 = buffer.data(dph0 + 228);
    const auto *dph0_229 = buffer.data(dph0 + 229);
    const auto *dph0_230 = buffer.data(dph0 + 230);
    const auto *dph0_236 = buffer.data(dph0 + 236);
    const auto *dph0_240 = buffer.data(dph0 + 240);
    const auto *dph0_246 = buffer.data(dph0 + 246);
    const auto *dph0_248 = buffer.data(dph0 + 248);

    const auto *dpg_0 = buffer.data(dpg + 0);
    const auto *dpg_2 = buffer.data(dpg + 2);
    const auto *dpg_5 = buffer.data(dpg + 5);
    const auto *dpg_14 = buffer.data(dpg + 14);
    const auto *dpg_15 = buffer.data(dpg + 15);
    const auto *dpg_18 = buffer.data(dpg + 18);
    const auto *dpg_30 = buffer.data(dpg + 30);
    const auto *dpg_45 = buffer.data(dpg + 45);
    const auto *dpg_50 = buffer.data(dpg + 50);
    const auto *dpg_55 = buffer.data(dpg + 55);
    const auto *dpg_59 = buffer.data(dpg + 59);
    const auto *dpg_60 = buffer.data(dpg + 60);
    const auto *dpg_65 = buffer.data(dpg + 65);
    const auto *dpg_75 = buffer.data(dpg + 75);
    const auto *dpg_80 = buffer.data(dpg + 80);
    const auto *dpg_101 = buffer.data(dpg + 101);
    const auto *dpg_102 = buffer.data(dpg + 102);
    const auto *dpg_104 = buffer.data(dpg + 104);
    const auto *dpg_115 = buffer.data(dpg + 115);
    const auto *dpg_116 = buffer.data(dpg + 116);
    const auto *dpg_117 = buffer.data(dpg + 117);
    const auto *dpg_119 = buffer.data(dpg + 119);
    const auto *dpg_120 = buffer.data(dpg + 120);
    const auto *dpg_125 = buffer.data(dpg + 125);
    const auto *dpg_129 = buffer.data(dpg + 129);
    const auto *dpg_130 = buffer.data(dpg + 130);
    const auto *dpg_131 = buffer.data(dpg + 131);
    const auto *dpg_132 = buffer.data(dpg + 132);
    const auto *dpg_134 = buffer.data(dpg + 134);
    const auto *dpg_135 = buffer.data(dpg + 135);
    const auto *dpg_138 = buffer.data(dpg + 138);
    const auto *dpg_141 = buffer.data(dpg + 141);
    const auto *dpg_145 = buffer.data(dpg + 145);
    const auto *dpg_147 = buffer.data(dpg + 147);
    const auto *dpg_148 = buffer.data(dpg + 148);
    const auto *dpg_149 = buffer.data(dpg + 149);
    const auto *dpg_150 = buffer.data(dpg + 150);
    const auto *dpg_153 = buffer.data(dpg + 153);
    const auto *dpg_156 = buffer.data(dpg + 156);
    const auto *dpg_160 = buffer.data(dpg + 160);
    const auto *dpg_162 = buffer.data(dpg + 162);
    const auto *dpg_163 = buffer.data(dpg + 163);
    const auto *dpg_164 = buffer.data(dpg + 164);
    const auto *dpg_170 = buffer.data(dpg + 170);
    const auto *dpg_174 = buffer.data(dpg + 174);
    const auto *dpg_175 = buffer.data(dpg + 175);
    const auto *dpg_177 = buffer.data(dpg + 177);
    const auto *dpg_178 = buffer.data(dpg + 178);
    const auto *dpg_179 = buffer.data(dpg + 179);

    const auto *dph1_3 = buffer.data(dph1 + 3);
    const auto *dph1_5 = buffer.data(dph1 + 5);
    const auto *dph1_6 = buffer.data(dph1 + 6);
    const auto *dph1_9 = buffer.data(dph1 + 9);
    const auto *dph1_10 = buffer.data(dph1 + 10);
    const auto *dph1_15 = buffer.data(dph1 + 15);
    const auto *dph1_21 = buffer.data(dph1 + 21);
    const auto *dph1_24 = buffer.data(dph1 + 24);
    const auto *dph1_27 = buffer.data(dph1 + 27);
    const auto *dph1_28 = buffer.data(dph1 + 28);
    const auto *dph1_36 = buffer.data(dph1 + 36);
    const auto *dph1_188 = buffer.data(dph1 + 188);
    const auto *dph1_210 = buffer.data(dph1 + 210);
    const auto *dph1_213 = buffer.data(dph1 + 213);
    const auto *dph1_216 = buffer.data(dph1 + 216);
    const auto *dph1_225 = buffer.data(dph1 + 225);
    const auto *dph1_227 = buffer.data(dph1 + 227);
    const auto *dph1_228 = buffer.data(dph1 + 228);
    const auto *dph1_229 = buffer.data(dph1 + 229);
    const auto *dph1_230 = buffer.data(dph1 + 230);
    const auto *dph1_236 = buffer.data(dph1 + 236);
    const auto *dph1_240 = buffer.data(dph1 + 240);
    const auto *dph1_246 = buffer.data(dph1 + 246);
    const auto *dph1_248 = buffer.data(dph1 + 248);

    const auto *fsh0_47 = buffer.data(fsh0 + 47);
    const auto *fsh0_51 = buffer.data(fsh0 + 51);
    const auto *fsh0_58 = buffer.data(fsh0 + 58);
    const auto *fsh0_59 = buffer.data(fsh0 + 59);
    const auto *fsh0_60 = buffer.data(fsh0 + 60);
    const auto *fsh0_62 = buffer.data(fsh0 + 62);
    const auto *fsh0_63 = buffer.data(fsh0 + 63);
    const auto *fsh0_66 = buffer.data(fsh0 + 66);
    const auto *fsh0_69 = buffer.data(fsh0 + 69);

    const auto *fsg_30 = buffer.data(fsg + 30);
    const auto *fsg_32 = buffer.data(fsg + 32);
    const auto *fsg_35 = buffer.data(fsg + 35);
    const auto *fsg_39 = buffer.data(fsg + 39);
    const auto *fsg_41 = buffer.data(fsg + 41);
    const auto *fsg_42 = buffer.data(fsg + 42);
    const auto *fsg_43 = buffer.data(fsg + 43);
    const auto *fsg_44 = buffer.data(fsg + 44);
    const auto *fsg_45 = buffer.data(fsg + 45);
    const auto *fsg_46 = buffer.data(fsg + 46);
    const auto *fsg_48 = buffer.data(fsg + 48);
    const auto *fsg_50 = buffer.data(fsg + 50);
    const auto *fsg_51 = buffer.data(fsg + 51);
    const auto *fsg_55 = buffer.data(fsg + 55);
    const auto *fsg_57 = buffer.data(fsg + 57);
    const auto *fsg_58 = buffer.data(fsg + 58);
    const auto *fsg_59 = buffer.data(fsg + 59);

    const auto *fsh1_47 = buffer.data(fsh1 + 47);
    const auto *fsh1_51 = buffer.data(fsh1 + 51);
    const auto *fsh1_58 = buffer.data(fsh1 + 58);
    const auto *fsh1_59 = buffer.data(fsh1 + 59);
    const auto *fsh1_60 = buffer.data(fsh1 + 60);
    const auto *fsh1_62 = buffer.data(fsh1 + 62);
    const auto *fsh1_63 = buffer.data(fsh1 + 63);
    const auto *fsh1_66 = buffer.data(fsh1 + 66);
    const auto *fsh1_69 = buffer.data(fsh1 + 69);

    const auto *fpf0_62 = buffer.data(fpf0 + 62);
    const auto *fpf0_67 = buffer.data(fpf0 + 67);
    const auto *fpf0_68 = buffer.data(fpf0 + 68);
    const auto *fpf0_69 = buffer.data(fpf0 + 69);
    const auto *fpf0_80 = buffer.data(fpf0 + 80);
    const auto *fpf0_81 = buffer.data(fpf0 + 81);
    const auto *fpf0_82 = buffer.data(fpf0 + 82);
    const auto *fpf0_85 = buffer.data(fpf0 + 85);
    const auto *fpf0_86 = buffer.data(fpf0 + 86);
    const auto *fpf0_87 = buffer.data(fpf0 + 87);
    const auto *fpf0_88 = buffer.data(fpf0 + 88);
    const auto *fpf0_89 = buffer.data(fpf0 + 89);
    const auto *fpf0_90 = buffer.data(fpf0 + 90);
    const auto *fpf0_92 = buffer.data(fpf0 + 92);
    const auto *fpf0_93 = buffer.data(fpf0 + 93);
    const auto *fpf0_96 = buffer.data(fpf0 + 96);
    const auto *fpf0_97 = buffer.data(fpf0 + 97);
    const auto *fpf0_99 = buffer.data(fpf0 + 99);
    const auto *fpf0_100 = buffer.data(fpf0 + 100);
    const auto *fpf0_102 = buffer.data(fpf0 + 102);

    const auto *fpf1_62 = buffer.data(fpf1 + 62);
    const auto *fpf1_67 = buffer.data(fpf1 + 67);
    const auto *fpf1_68 = buffer.data(fpf1 + 68);
    const auto *fpf1_69 = buffer.data(fpf1 + 69);
    const auto *fpf1_80 = buffer.data(fpf1 + 80);
    const auto *fpf1_81 = buffer.data(fpf1 + 81);
    const auto *fpf1_82 = buffer.data(fpf1 + 82);
    const auto *fpf1_85 = buffer.data(fpf1 + 85);
    const auto *fpf1_86 = buffer.data(fpf1 + 86);
    const auto *fpf1_87 = buffer.data(fpf1 + 87);
    const auto *fpf1_88 = buffer.data(fpf1 + 88);
    const auto *fpf1_89 = buffer.data(fpf1 + 89);
    const auto *fpf1_90 = buffer.data(fpf1 + 90);
    const auto *fpf1_92 = buffer.data(fpf1 + 92);
    const auto *fpf1_93 = buffer.data(fpf1 + 93);
    const auto *fpf1_96 = buffer.data(fpf1 + 96);
    const auto *fpf1_97 = buffer.data(fpf1 + 97);
    const auto *fpf1_99 = buffer.data(fpf1 + 99);
    const auto *fpf1_100 = buffer.data(fpf1 + 100);
    const auto *fpf1_102 = buffer.data(fpf1 + 102);

    const auto *fpg_90 = buffer.data(fpg + 90);
    const auto *fpg_92 = buffer.data(fpg + 92);
    const auto *fpg_94 = buffer.data(fpg + 94);
    const auto *fpg_95 = buffer.data(fpg + 95);
    const auto *fpg_99 = buffer.data(fpg + 99);
    const auto *fpg_101 = buffer.data(fpg + 101);
    const auto *fpg_102 = buffer.data(fpg + 102);
    const auto *fpg_103 = buffer.data(fpg + 103);
    const auto *fpg_104 = buffer.data(fpg + 104);
    const auto *fpg_105 = buffer.data(fpg + 105);
    const auto *fpg_107 = buffer.data(fpg + 107);
    const auto *fpg_110 = buffer.data(fpg + 110);
    const auto *fpg_114 = buffer.data(fpg + 114);
    const auto *fpg_115 = buffer.data(fpg + 115);
    const auto *fpg_116 = buffer.data(fpg + 116);
    const auto *fpg_117 = buffer.data(fpg + 117);
    const auto *fpg_119 = buffer.data(fpg + 119);
    const auto *fpg_120 = buffer.data(fpg + 120);
    const auto *fpg_121 = buffer.data(fpg + 121);
    const auto *fpg_122 = buffer.data(fpg + 122);
    const auto *fpg_123 = buffer.data(fpg + 123);
    const auto *fpg_124 = buffer.data(fpg + 124);
    const auto *fpg_125 = buffer.data(fpg + 125);
    const auto *fpg_129 = buffer.data(fpg + 129);
    const auto *fpg_130 = buffer.data(fpg + 130);
    const auto *fpg_131 = buffer.data(fpg + 131);
    const auto *fpg_132 = buffer.data(fpg + 132);
    const auto *fpg_133 = buffer.data(fpg + 133);
    const auto *fpg_134 = buffer.data(fpg + 134);
    const auto *fpg_135 = buffer.data(fpg + 135);
    const auto *fpg_136 = buffer.data(fpg + 136);
    const auto *fpg_137 = buffer.data(fpg + 137);
    const auto *fpg_138 = buffer.data(fpg + 138);
    const auto *fpg_140 = buffer.data(fpg + 140);
    const auto *fpg_141 = buffer.data(fpg + 141);
    const auto *fpg_145 = buffer.data(fpg + 145);
    const auto *fpg_146 = buffer.data(fpg + 146);
    const auto *fpg_147 = buffer.data(fpg + 147);
    const auto *fpg_148 = buffer.data(fpg + 148);
    const auto *fpg_149 = buffer.data(fpg + 149);
    const auto *fpg_150 = buffer.data(fpg + 150);
    const auto *fpg_151 = buffer.data(fpg + 151);
    const auto *fpg_152 = buffer.data(fpg + 152);
    const auto *fpg_153 = buffer.data(fpg + 153);
    const auto *fpg_155 = buffer.data(fpg + 155);
    const auto *fpg_156 = buffer.data(fpg + 156);
    const auto *fpg_160 = buffer.data(fpg + 160);
    const auto *fpg_162 = buffer.data(fpg + 162);
    const auto *fpg_163 = buffer.data(fpg + 163);
    const auto *fpg_164 = buffer.data(fpg + 164);
    const auto *fpg_165 = buffer.data(fpg + 165);
    const auto *fpg_166 = buffer.data(fpg + 166);
    const auto *fpg_168 = buffer.data(fpg + 168);
    const auto *fpg_170 = buffer.data(fpg + 170);
    const auto *fpg_171 = buffer.data(fpg + 171);
    const auto *fpg_175 = buffer.data(fpg + 175);
    const auto *fpg_177 = buffer.data(fpg + 177);
    const auto *fpg_178 = buffer.data(fpg + 178);
    const auto *fpg_179 = buffer.data(fpg + 179);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pa_z, pc_y, pc_z, dph0_3, dph0_5, dpg_0, \
                         dpg_2, dph1_3, dph1_5, fpg_90, fpg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_1 * dpg_0[k]
                   + f_4 * pc_z[k] * fpg_90[k];

        t_129[k] = pa_z[k] * dph0_3[k]
                   - f_9 * pc_z[k] * dph1_3[k];

        t_130[k] = f_4 * pc_y[k] * fpg_92[k];

        t_131[k] = pa_z[k] * dph0_5[k]
                   + f_10 * dpg_2[k]
                   - f_9 * pc_z[k] * dph1_5[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_z, pc_y, pc_z, dph0_6, dph0_9, dpg_5, \
                         dph1_6, dph1_9, fpf0_62, fpf1_62, fpg_94, \
                         fpg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_z[k] * dph0_6[k]
                   - f_9 * pc_z[k] * dph1_6[k];

        t_133[k] = f_5 * fpf0_62[k]
                   - f_6 * fpf1_62[k]
                   + f_4 * pc_y[k] * fpg_94[k];

        t_134[k] = f_4 * pc_y[k] * fpg_95[k];

        t_135[k] = pa_z[k] * dph0_9[k]
                   + f_0 * dpg_5[k]
                   - f_9 * pc_z[k] * dph1_9[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_z, pc_x, pc_z, dph0_10, dpg_101, dpg_102, \
                         dph1_10, fsg_41, fsg_42, fpg_101, fpg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_z[k] * dph0_10[k]
                   - f_9 * pc_z[k] * dph1_10[k];

        t_137[k] = f_10 * dpg_101[k]
                   + f_1 * fsg_41[k]
                   + f_4 * pc_x[k] * fpg_101[k];

        t_138[k] = f_10 * dpg_102[k]
                   + f_1 * fsg_42[k]
                   + f_4 * pc_x[k] * fpg_102[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pa_z, pc_x, pc_y, pc_z, dph0_15, dpg_104, \
                         dph1_15, fsg_44, fpg_99, fpg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_4 * pc_y[k] * fpg_99[k];

        t_140[k] = f_10 * dpg_104[k]
                   + f_1 * fsg_44[k]
                   + f_4 * pc_x[k] * fpg_104[k];

        t_141[k] = pa_z[k] * dph0_15[k]
                   - f_9 * pc_z[k] * dph1_15[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, fpf0_67, fpf0_68, fpf0_69, fpf1_67, \
                         fpf1_68, fpf1_69, fpg_101, fpg_102, fpg_103, \
                         fpg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * fpf0_67[k]
                   - f_15 * fpf1_67[k]
                   + f_4 * pc_y[k] * fpg_101[k];

        t_143[k] = f_7 * fpf0_68[k]
                   - f_8 * fpf1_68[k]
                   + f_4 * pc_y[k] * fpg_102[k];

        t_144[k] = f_5 * fpf0_69[k]
                   - f_6 * fpf1_69[k]
                   + f_4 * pc_y[k] * fpg_103[k];

        t_145[k] = f_4 * pc_y[k] * fpg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, dph0_21, dpg_14, \
                         dpg_15, dph1_21, fsg_30, fpf0_69, fpf1_69, fpg_104, \
                         fpg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * dpg_14[k]
                   + f_2 * fpf0_69[k]
                   - f_3 * fpf1_69[k]
                   + f_4 * pc_z[k] * fpg_104[k];

        t_147[k] = pa_z[k] * dph0_21[k]
                   - f_9 * pc_z[k] * dph1_21[k];

        t_148[k] = f_1 * fsg_30[k]
                   + f_4 * pc_y[k] * fpg_105[k];

        t_149[k] = f_1 * dpg_15[k]
                   + f_4 * pc_z[k] * fpg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_z, pb_y, pc_y, pc_z, dph0_24, dph0_27, \
                         dph1_24, dph1_27, fsh0_47, fsg_32, fsh1_47, \
                         fpg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * dph0_24[k]
                   - f_9 * pc_z[k] * dph1_24[k];

        t_151[k] = f_1 * fsg_32[k]
                   + f_4 * pc_y[k] * fpg_107[k];

        t_152[k] = pb_y[k] * fsh0_47[k]
                   - f_9 * pc_y[k] * fsh1_47[k];

        t_153[k] = pa_z[k] * dph0_27[k]
                   - f_9 * pc_z[k] * dph1_27[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_z, pb_y, pc_y, pc_z, dph0_28, dpg_18, \
                         dph1_28, fsh0_51, fsg_35, fsh1_51, fpg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pa_z[k] * dph0_28[k]
                   + f_1 * dpg_18[k]
                   - f_9 * pc_z[k] * dph1_28[k];

        t_155[k] = f_1 * fsg_35[k]
                   + f_4 * pc_y[k] * fpg_110[k];

        t_156[k] = pb_y[k] * fsh0_51[k]
                   - f_9 * pc_y[k] * fsh1_51[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, pc_y, dpg_115, dpg_116, dpg_117, \
                         fsg_39, fpg_114, fpg_115, fpg_116, fpg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_10 * dpg_115[k]
                   + f_4 * pc_x[k] * fpg_115[k];

        t_158[k] = f_10 * dpg_116[k]
                   + f_4 * pc_x[k] * fpg_116[k];

        t_159[k] = f_10 * dpg_117[k]
                   + f_4 * pc_x[k] * fpg_117[k];

        t_160[k] = f_1 * fsg_39[k]
                   + f_4 * pc_y[k] * fpg_114[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pa_z, pb_y, pc_x, pc_y, pc_z, dph0_36, dpg_119, \
                         dph1_36, fsh0_58, fsg_41, fsh1_58, fpg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_10 * dpg_119[k]
                   + f_4 * pc_x[k] * fpg_119[k];

        t_162[k] = pa_z[k] * dph0_36[k]
                   - f_9 * pc_z[k] * dph1_36[k];

        t_163[k] = pb_y[k] * fsh0_58[k]
                   + f_16 * fsg_41[k]
                   - f_9 * pc_y[k] * fsh1_58[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pb_y, pc_y, fsh0_59, fsh0_60, fsh0_62, \
                         fsg_42, fsg_43, fsg_44, fsh1_59, fsh1_60, fsh1_62, \
                         fpg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pb_y[k] * fsh0_59[k]
                   + f_0 * fsg_42[k]
                   - f_9 * pc_y[k] * fsh1_59[k];

        t_165[k] = pb_y[k] * fsh0_60[k]
                   + f_10 * fsg_43[k]
                   - f_9 * pc_y[k] * fsh1_60[k];

        t_166[k] = f_1 * fsg_44[k]
                   + f_4 * pc_y[k] * fpg_119[k];

        t_167[k] = pb_y[k] * fsh0_62[k]
                   - f_9 * pc_y[k] * fsh1_62[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pc_x, pc_y, pc_z, dpg_30, dpg_120, \
                         fsg_30, fpf0_80, fpf1_80, fpg_120, fpg_121, \
                         fpg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_10 * dpg_120[k]
                   + f_2 * fpf0_80[k]
                   - f_3 * fpf1_80[k]
                   + f_4 * pc_x[k] * fpg_120[k];

        t_169[k] = f_4 * pc_y[k] * fpg_120[k];

        t_170[k] = f_1 * dpg_30[k]
                   + f_1 * fsg_30[k]
                   + f_4 * pc_z[k] * fpg_120[k];

        t_171[k] = f_5 * fpf0_80[k]
                   - f_6 * fpf1_80[k]
                   + f_4 * pc_y[k] * fpg_121[k];

        t_172[k] = f_4 * pc_y[k] * fpg_122[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pc_x, pc_y, dpg_125, fpf0_81, fpf0_82, \
                         fpf0_85, fpf1_81, fpf1_82, fpf1_85, fpg_123, fpg_124, \
                         fpg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_10 * dpg_125[k]
                   + f_7 * fpf0_85[k]
                   - f_8 * fpf1_85[k]
                   + f_4 * pc_x[k] * fpg_125[k];

        t_174[k] = f_7 * fpf0_81[k]
                   - f_8 * fpf1_81[k]
                   + f_4 * pc_y[k] * fpg_123[k];

        t_175[k] = f_5 * fpf0_82[k]
                   - f_6 * fpf1_82[k]
                   + f_4 * pc_y[k] * fpg_124[k];

        t_176[k] = f_4 * pc_y[k] * fpg_125[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pc_x, dpg_129, dpg_130, dpg_131, dpg_132, \
                         fpf0_89, fpf1_89, fpg_129, fpg_130, fpg_131, \
                         fpg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_10 * dpg_129[k]
                   + f_5 * fpf0_89[k]
                   - f_6 * fpf1_89[k]
                   + f_4 * pc_x[k] * fpg_129[k];

        t_178[k] = f_10 * dpg_130[k]
                   + f_4 * pc_x[k] * fpg_130[k];

        t_179[k] = f_10 * dpg_131[k]
                   + f_4 * pc_x[k] * fpg_131[k];

        t_180[k] = f_10 * dpg_132[k]
                   + f_4 * pc_x[k] * fpg_132[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, dpg_134, fpf0_86, fpf0_87, \
                         fpf1_86, fpf1_87, fpg_129, fpg_130, fpg_131, \
                         fpg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_4 * pc_y[k] * fpg_129[k];

        t_182[k] = f_10 * dpg_134[k]
                   + f_4 * pc_x[k] * fpg_134[k];

        t_183[k] = f_2 * fpf0_86[k]
                   - f_3 * fpf1_86[k]
                   + f_4 * pc_y[k] * fpg_130[k];

        t_184[k] = f_14 * fpf0_87[k]
                   - f_15 * fpf1_87[k]
                   + f_4 * pc_y[k] * fpg_131[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, pc_y, fpf0_88, fpf0_89, fpf1_88, fpf1_89, \
                         fpg_132, fpg_133, fpg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_7 * fpf0_88[k]
                   - f_8 * fpf1_88[k]
                   + f_4 * pc_y[k] * fpg_132[k];

        t_186[k] = f_5 * fpf0_89[k]
                   - f_6 * fpf1_89[k]
                   + f_4 * pc_y[k] * fpg_133[k];

        t_187[k] = f_4 * pc_y[k] * fpg_134[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_x, pc_x, pc_y, pph0_188, pph1_188, dph0_188, \
                         dpg_45, dpg_135, dph1_188, fsg_45, fpf0_90, fpf1_90, \
                         fpg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_12 * pph0_188[k]
                   - f_13 * pph1_188[k]
                   + pa_x[k] * dph0_188[k]
                   - f_9 * pc_x[k] * dph1_188[k];

        t_189[k] = f_1 * dpg_135[k]
                   + f_1 * fsg_45[k]
                   + f_2 * fpf0_90[k]
                   - f_3 * fpf1_90[k]
                   + f_4 * pc_x[k] * fpg_135[k];

        t_190[k] = f_10 * dpg_45[k]
                   + f_4 * pc_y[k] * fpg_135[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pc_x, pc_z, dpg_138, fsg_48, fpf0_90, \
                         fpf0_93, fpf1_90, fpf1_93, fpg_135, fpg_136, fpg_137, \
                         fpg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_4 * pc_z[k] * fpg_135[k];

        t_192[k] = f_1 * dpg_138[k]
                   + f_1 * fsg_48[k]
                   + f_7 * fpf0_93[k]
                   - f_8 * fpf1_93[k]
                   + f_4 * pc_x[k] * fpg_138[k];

        t_193[k] = f_4 * pc_z[k] * fpg_136[k];

        t_194[k] = f_5 * fpf0_90[k]
                   - f_6 * fpf1_90[k]
                   + f_4 * pc_z[k] * fpg_137[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pc_x, pc_y, pc_z, dpg_50, dpg_141, fsg_51, \
                         fpf0_96, fpf1_96, fpg_138, fpg_140, fpg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_1 * dpg_141[k]
                   + f_1 * fsg_51[k]
                   + f_5 * fpf0_96[k]
                   - f_6 * fpf1_96[k]
                   + f_4 * pc_x[k] * fpg_141[k];

        t_196[k] = f_4 * pc_z[k] * fpg_138[k];

        t_197[k] = f_10 * dpg_50[k]
                   + f_4 * pc_y[k] * fpg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_z, dpg_145, dpg_147, fsg_55, \
                         fsg_57, fpf0_92, fpf1_92, fpg_140, fpg_141, fpg_145, \
                         fpg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_7 * fpf0_92[k]
                   - f_8 * fpf1_92[k]
                   + f_4 * pc_z[k] * fpg_140[k];

        t_199[k] = f_1 * dpg_145[k]
                   + f_1 * fsg_55[k]
                   + f_4 * pc_x[k] * fpg_145[k];

        t_200[k] = f_4 * pc_z[k] * fpg_141[k];

        t_201[k] = f_1 * dpg_147[k]
                   + f_1 * fsg_57[k]
                   + f_4 * pc_x[k] * fpg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, pc_x, pc_y, dpg_55, dpg_148, dpg_149, fsg_58, \
                         fsg_59, fpf0_96, fpf1_96, fpg_145, fpg_148, \
                         fpg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_1 * dpg_148[k]
                   + f_1 * fsg_58[k]
                   + f_4 * pc_x[k] * fpg_148[k];

        t_203[k] = f_1 * dpg_149[k]
                   + f_1 * fsg_59[k]
                   + f_4 * pc_x[k] * fpg_149[k];

        t_204[k] = f_10 * dpg_55[k]
                   + f_2 * fpf0_96[k]
                   - f_3 * fpf1_96[k]
                   + f_4 * pc_y[k] * fpg_145[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_y, pc_z, dpg_59, fpf0_96, fpf0_97, \
                         fpf1_96, fpf1_97, fpg_145, fpg_146, fpg_147, \
                         fpg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_4 * pc_z[k] * fpg_145[k];

        t_206[k] = f_5 * fpf0_96[k]
                   - f_6 * fpf1_96[k]
                   + f_4 * pc_z[k] * fpg_146[k];

        t_207[k] = f_7 * fpf0_97[k]
                   - f_8 * fpf1_97[k]
                   + f_4 * pc_z[k] * fpg_147[k];

        t_208[k] = f_10 * dpg_59[k]
                   + f_4 * pc_y[k] * fpg_149[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pa_x, pc_x, pc_y, pc_z, dph0_210, dpg_60, \
                         dpg_150, dph1_210, fsg_45, fpf0_99, fpf1_99, fpg_149, \
                         fpg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_2 * fpf0_99[k]
                   - f_3 * fpf1_99[k]
                   + f_4 * pc_z[k] * fpg_149[k];

        t_210[k] = pa_x[k] * dph0_210[k]
                   + f_11 * dpg_150[k]
                   - f_9 * pc_x[k] * dph1_210[k];

        t_211[k] = f_10 * dpg_60[k]
                   + f_1 * fsg_45[k]
                   + f_4 * pc_y[k] * fpg_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_x, pc_x, pc_z, dph0_213, dpg_153, \
                         dph1_213, fpf0_100, fpf1_100, fpg_150, fpg_151, \
                         fpg_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_4 * pc_z[k] * fpg_150[k];

        t_213[k] = pa_x[k] * dph0_213[k]
                   + f_0 * dpg_153[k]
                   - f_9 * pc_x[k] * dph1_213[k];

        t_214[k] = f_4 * pc_z[k] * fpg_151[k];

        t_215[k] = f_5 * fpf0_100[k]
                   - f_6 * fpf1_100[k]
                   + f_4 * pc_z[k] * fpg_152[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pa_x, pc_x, pc_y, pc_z, dph0_216, dpg_65, \
                         dpg_156, dph1_216, fsg_50, fpg_153, fpg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pa_x[k] * dph0_216[k]
                   + f_10 * dpg_156[k]
                   - f_9 * pc_x[k] * dph1_216[k];

        t_217[k] = f_4 * pc_z[k] * fpg_153[k];

        t_218[k] = f_10 * dpg_65[k]
                   + f_1 * fsg_50[k]
                   + f_4 * pc_y[k] * fpg_155[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pc_x, pc_z, dpg_160, dpg_162, fpf0_102, \
                         fpf1_102, fpg_155, fpg_156, fpg_160, fpg_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_7 * fpf0_102[k]
                   - f_8 * fpf1_102[k]
                   + f_4 * pc_z[k] * fpg_155[k];

        t_220[k] = f_1 * dpg_160[k]
                   + f_4 * pc_x[k] * fpg_160[k];

        t_221[k] = f_4 * pc_z[k] * fpg_156[k];

        t_222[k] = f_1 * dpg_162[k]
                   + f_4 * pc_x[k] * fpg_162[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_x, pc_x, pc_z, dph0_225, dpg_163, \
                         dpg_164, dph1_225, fpg_160, fpg_163, fpg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_1 * dpg_163[k]
                   + f_4 * pc_x[k] * fpg_163[k];

        t_224[k] = f_1 * dpg_164[k]
                   + f_4 * pc_x[k] * fpg_164[k];

        t_225[k] = pa_x[k] * dph0_225[k]
                   - f_9 * pc_x[k] * dph1_225[k];

        t_226[k] = f_4 * pc_z[k] * fpg_160[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, pa_x, pc_x, dph0_227, dph0_228, dph0_229, \
                         dph0_230, dph1_227, dph1_228, dph1_229, \
                         dph1_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pa_x[k] * dph0_227[k]
                   - f_9 * pc_x[k] * dph1_227[k];

        t_228[k] = pa_x[k] * dph0_228[k]
                   - f_9 * pc_x[k] * dph1_228[k];

        t_229[k] = pa_x[k] * dph0_229[k]
                   - f_9 * pc_x[k] * dph1_229[k];

        t_230[k] = pa_x[k] * dph0_230[k]
                   - f_9 * pc_x[k] * dph1_230[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_z, pc_y, pc_z, dpg_75, fsh0_63, \
                         fsh0_66, fsg_45, fsh1_63, fsh1_66, fpg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = pb_z[k] * fsh0_63[k]
                   - f_9 * pc_z[k] * fsh1_63[k];

        t_232[k] = f_10 * dpg_75[k]
                   + f_4 * pc_y[k] * fpg_165[k];

        t_233[k] = f_1 * fsg_45[k]
                   + f_4 * pc_z[k] * fpg_165[k];

        t_234[k] = pb_z[k] * fsh0_66[k]
                   - f_9 * pc_z[k] * fsh1_66[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pa_x, pb_z, pc_x, pc_z, dph0_236, dpg_170, \
                         dph1_236, fsh0_69, fsg_46, fsh1_69, fpg_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_1 * fsg_46[k]
                   + f_4 * pc_z[k] * fpg_166[k];

        t_236[k] = pa_x[k] * dph0_236[k]
                   + f_0 * dpg_170[k]
                   - f_9 * pc_x[k] * dph1_236[k];

        t_237[k] = pb_z[k] * fsh0_69[k]
                   - f_9 * pc_z[k] * fsh1_69[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pa_x, pc_x, pc_y, pc_z, dph0_240, dpg_80, \
                         dpg_174, dph1_240, fsg_48, fpg_168, fpg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_1 * fsg_48[k]
                   + f_4 * pc_z[k] * fpg_168[k];

        t_239[k] = f_10 * dpg_80[k]
                   + f_4 * pc_y[k] * fpg_170[k];

        t_240[k] = pa_x[k] * dph0_240[k]
                   + f_10 * dpg_174[k]
                   - f_9 * pc_x[k] * dph1_240[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pc_x, pc_z, dpg_175, dpg_177, dpg_178, \
                         fsg_51, fpg_171, fpg_175, fpg_177, fpg_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_1 * dpg_175[k]
                   + f_4 * pc_x[k] * fpg_175[k];

        t_242[k] = f_1 * fsg_51[k]
                   + f_4 * pc_z[k] * fpg_171[k];

        t_243[k] = f_1 * dpg_177[k]
                   + f_4 * pc_x[k] * fpg_177[k];

        t_244[k] = f_1 * dpg_178[k]
                   + f_4 * pc_x[k] * fpg_178[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pa_x, pc_x, pc_z, dph0_246, dph0_248, \
                         dpg_179, dph1_246, dph1_248, fsg_55, fpg_175, \
                         fpg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_1 * dpg_179[k]
                   + f_4 * pc_x[k] * fpg_179[k];

        t_246[k] = pa_x[k] * dph0_246[k]
                   - f_9 * pc_x[k] * dph1_246[k];

        t_247[k] = f_1 * fsg_55[k]
                   + f_4 * pc_z[k] * fpg_175[k];

        t_248[k] = pa_x[k] * dph0_248[k]
                   - f_9 * pc_x[k] * dph1_248[k];
    }
}

static auto
compute_prim_fph_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dph0, const size_t dpg,
                                                          const size_t dph1, const size_t fsh0,
                                                          const size_t fsg, const size_t fsh1,
                                                          const size_t fpf0, const size_t fpf1,
                                                          const size_t fpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_14 = 1.5 / gamma;
    const auto f_15 = 1.5 * p / (gamma * q);

    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dph0_66 = buffer.data(dph0 + 66);
    const auto *dph0_69 = buffer.data(dph0 + 69);
    const auto *dph0_73 = buffer.data(dph0 + 73);
    const auto *dph0_85 = buffer.data(dph0 + 85);
    const auto *dph0_87 = buffer.data(dph0 + 87);
    const auto *dph0_90 = buffer.data(dph0 + 90);
    const auto *dph0_126 = buffer.data(dph0 + 126);
    const auto *dph0_131 = buffer.data(dph0 + 131);
    const auto *dph0_135 = buffer.data(dph0 + 135);
    const auto *dph0_140 = buffer.data(dph0 + 140);
    const auto *dph0_168 = buffer.data(dph0 + 168);
    const auto *dph0_170 = buffer.data(dph0 + 170);
    const auto *dph0_173 = buffer.data(dph0 + 173);
    const auto *dph0_177 = buffer.data(dph0 + 177);
    const auto *dph0_249 = buffer.data(dph0 + 249);
    const auto *dph0_251 = buffer.data(dph0 + 251);
    const auto *dph0_288 = buffer.data(dph0 + 288);
    const auto *dph0_290 = buffer.data(dph0 + 290);
    const auto *dph0_291 = buffer.data(dph0 + 291);
    const auto *dph0_292 = buffer.data(dph0 + 292);
    const auto *dph0_293 = buffer.data(dph0 + 293);
    const auto *dph0_309 = buffer.data(dph0 + 309);
    const auto *dph0_310 = buffer.data(dph0 + 310);
    const auto *dph0_311 = buffer.data(dph0 + 311);
    const auto *dph0_312 = buffer.data(dph0 + 312);
    const auto *dph0_314 = buffer.data(dph0 + 314);
    const auto *dph0_339 = buffer.data(dph0 + 339);
    const auto *dph0_342 = buffer.data(dph0 + 342);
    const auto *dph0_343 = buffer.data(dph0 + 343);
    const auto *dph0_351 = buffer.data(dph0 + 351);
    const auto *dph0_352 = buffer.data(dph0 + 352);
    const auto *dph0_353 = buffer.data(dph0 + 353);
    const auto *dph0_354 = buffer.data(dph0 + 354);
    const auto *dph0_356 = buffer.data(dph0 + 356);
    const auto *dph0_357 = buffer.data(dph0 + 357);
    const auto *dph0_362 = buffer.data(dph0 + 362);

    const auto *dpg_45 = buffer.data(dpg + 45);
    const auto *dpg_48 = buffer.data(dpg + 48);
    const auto *dpg_55 = buffer.data(dpg + 55);
    const auto *dpg_59 = buffer.data(dpg + 59);
    const auto *dpg_60 = buffer.data(dpg + 60);
    const auto *dpg_63 = buffer.data(dpg + 63);
    const auto *dpg_70 = buffer.data(dpg + 70);
    const auto *dpg_78 = buffer.data(dpg + 78);
    const auto *dpg_89 = buffer.data(dpg + 89);
    const auto *dpg_90 = buffer.data(dpg + 90);
    const auto *dpg_92 = buffer.data(dpg + 92);
    const auto *dpg_95 = buffer.data(dpg + 95);
    const auto *dpg_100 = buffer.data(dpg + 100);
    const auto *dpg_102 = buffer.data(dpg + 102);
    const auto *dpg_103 = buffer.data(dpg + 103);
    const auto *dpg_104 = buffer.data(dpg + 104);
    const auto *dpg_105 = buffer.data(dpg + 105);
    const auto *dpg_107 = buffer.data(dpg + 107);
    const auto *dpg_110 = buffer.data(dpg + 110);
    const auto *dpg_120 = buffer.data(dpg + 120);
    const auto *dpg_122 = buffer.data(dpg + 122);
    const auto *dpg_125 = buffer.data(dpg + 125);
    const auto *dpg_134 = buffer.data(dpg + 134);
    const auto *dpg_191 = buffer.data(dpg + 191);
    const auto *dpg_192 = buffer.data(dpg + 192);
    const auto *dpg_193 = buffer.data(dpg + 193);
    const auto *dpg_195 = buffer.data(dpg + 195);
    const auto *dpg_200 = buffer.data(dpg + 200);
    const auto *dpg_204 = buffer.data(dpg + 204);
    const auto *dpg_205 = buffer.data(dpg + 205);
    const auto *dpg_206 = buffer.data(dpg + 206);
    const auto *dpg_207 = buffer.data(dpg + 207);
    const auto *dpg_208 = buffer.data(dpg + 208);
    const auto *dpg_209 = buffer.data(dpg + 209);
    const auto *dpg_213 = buffer.data(dpg + 213);
    const auto *dpg_216 = buffer.data(dpg + 216);
    const auto *dpg_220 = buffer.data(dpg + 220);
    const auto *dpg_221 = buffer.data(dpg + 221);
    const auto *dpg_222 = buffer.data(dpg + 222);
    const auto *dpg_223 = buffer.data(dpg + 223);
    const auto *dpg_224 = buffer.data(dpg + 224);
    const auto *dpg_225 = buffer.data(dpg + 225);
    const auto *dpg_230 = buffer.data(dpg + 230);
    const auto *dpg_234 = buffer.data(dpg + 234);
    const auto *dpg_235 = buffer.data(dpg + 235);
    const auto *dpg_236 = buffer.data(dpg + 236);
    const auto *dpg_237 = buffer.data(dpg + 237);
    const auto *dpg_239 = buffer.data(dpg + 239);
    const auto *dpg_243 = buffer.data(dpg + 243);
    const auto *dpg_246 = buffer.data(dpg + 246);
    const auto *dpg_247 = buffer.data(dpg + 247);
    const auto *dpg_250 = buffer.data(dpg + 250);
    const auto *dpg_251 = buffer.data(dpg + 251);
    const auto *dpg_252 = buffer.data(dpg + 252);
    const auto *dpg_254 = buffer.data(dpg + 254);
    const auto *dpg_255 = buffer.data(dpg + 255);
    const auto *dpg_260 = buffer.data(dpg + 260);

    const auto *dph1_66 = buffer.data(dph1 + 66);
    const auto *dph1_69 = buffer.data(dph1 + 69);
    const auto *dph1_73 = buffer.data(dph1 + 73);
    const auto *dph1_85 = buffer.data(dph1 + 85);
    const auto *dph1_87 = buffer.data(dph1 + 87);
    const auto *dph1_90 = buffer.data(dph1 + 90);
    const auto *dph1_126 = buffer.data(dph1 + 126);
    const auto *dph1_131 = buffer.data(dph1 + 131);
    const auto *dph1_135 = buffer.data(dph1 + 135);
    const auto *dph1_140 = buffer.data(dph1 + 140);
    const auto *dph1_168 = buffer.data(dph1 + 168);
    const auto *dph1_170 = buffer.data(dph1 + 170);
    const auto *dph1_173 = buffer.data(dph1 + 173);
    const auto *dph1_177 = buffer.data(dph1 + 177);
    const auto *dph1_249 = buffer.data(dph1 + 249);
    const auto *dph1_251 = buffer.data(dph1 + 251);
    const auto *dph1_288 = buffer.data(dph1 + 288);
    const auto *dph1_290 = buffer.data(dph1 + 290);
    const auto *dph1_291 = buffer.data(dph1 + 291);
    const auto *dph1_292 = buffer.data(dph1 + 292);
    const auto *dph1_293 = buffer.data(dph1 + 293);
    const auto *dph1_309 = buffer.data(dph1 + 309);
    const auto *dph1_310 = buffer.data(dph1 + 310);
    const auto *dph1_311 = buffer.data(dph1 + 311);
    const auto *dph1_312 = buffer.data(dph1 + 312);
    const auto *dph1_314 = buffer.data(dph1 + 314);
    const auto *dph1_339 = buffer.data(dph1 + 339);
    const auto *dph1_342 = buffer.data(dph1 + 342);
    const auto *dph1_343 = buffer.data(dph1 + 343);
    const auto *dph1_351 = buffer.data(dph1 + 351);
    const auto *dph1_352 = buffer.data(dph1 + 352);
    const auto *dph1_353 = buffer.data(dph1 + 353);
    const auto *dph1_354 = buffer.data(dph1 + 354);
    const auto *dph1_356 = buffer.data(dph1 + 356);
    const auto *dph1_357 = buffer.data(dph1 + 357);
    const auto *dph1_362 = buffer.data(dph1 + 362);

    const auto *fsh0_105 = buffer.data(fsh0 + 105);
    const auto *fsh0_110 = buffer.data(fsh0 + 110);
    const auto *fsh0_114 = buffer.data(fsh0 + 114);

    const auto *fsg_62 = buffer.data(fsg + 62);
    const auto *fsg_63 = buffer.data(fsg + 63);
    const auto *fsg_65 = buffer.data(fsg + 65);
    const auto *fsg_71 = buffer.data(fsg + 71);
    const auto *fsg_72 = buffer.data(fsg + 72);
    const auto *fsg_73 = buffer.data(fsg + 73);
    const auto *fsg_75 = buffer.data(fsg + 75);
    const auto *fsg_77 = buffer.data(fsg + 77);
    const auto *fsg_80 = buffer.data(fsg + 80);
    const auto *fsg_84 = buffer.data(fsg + 84);
    const auto *fsg_85 = buffer.data(fsg + 85);
    const auto *fsg_86 = buffer.data(fsg + 86);
    const auto *fsg_87 = buffer.data(fsg + 87);
    const auto *fsg_89 = buffer.data(fsg + 89);

    const auto *fsh1_105 = buffer.data(fsh1 + 105);
    const auto *fsh1_110 = buffer.data(fsh1 + 110);
    const auto *fsh1_114 = buffer.data(fsh1 + 114);

    const auto *fpf0_126 = buffer.data(fpf0 + 126);
    const auto *fpf0_128 = buffer.data(fpf0 + 128);
    const auto *fpf0_129 = buffer.data(fpf0 + 129);
    const auto *fpf0_130 = buffer.data(fpf0 + 130);
    const auto *fpf0_135 = buffer.data(fpf0 + 135);
    const auto *fpf0_139 = buffer.data(fpf0 + 139);
    const auto *fpf0_143 = buffer.data(fpf0 + 143);
    const auto *fpf0_146 = buffer.data(fpf0 + 146);
    const auto *fpf0_150 = buffer.data(fpf0 + 150);
    const auto *fpf0_151 = buffer.data(fpf0 + 151);
    const auto *fpf0_152 = buffer.data(fpf0 + 152);
    const auto *fpf0_155 = buffer.data(fpf0 + 155);
    const auto *fpf0_156 = buffer.data(fpf0 + 156);
    const auto *fpf0_157 = buffer.data(fpf0 + 157);
    const auto *fpf0_158 = buffer.data(fpf0 + 158);
    const auto *fpf0_159 = buffer.data(fpf0 + 159);
    const auto *fpf0_170 = buffer.data(fpf0 + 170);
    const auto *fpf0_171 = buffer.data(fpf0 + 171);
    const auto *fpf0_172 = buffer.data(fpf0 + 172);

    const auto *fpf1_126 = buffer.data(fpf1 + 126);
    const auto *fpf1_128 = buffer.data(fpf1 + 128);
    const auto *fpf1_129 = buffer.data(fpf1 + 129);
    const auto *fpf1_130 = buffer.data(fpf1 + 130);
    const auto *fpf1_135 = buffer.data(fpf1 + 135);
    const auto *fpf1_139 = buffer.data(fpf1 + 139);
    const auto *fpf1_143 = buffer.data(fpf1 + 143);
    const auto *fpf1_146 = buffer.data(fpf1 + 146);
    const auto *fpf1_150 = buffer.data(fpf1 + 150);
    const auto *fpf1_151 = buffer.data(fpf1 + 151);
    const auto *fpf1_152 = buffer.data(fpf1 + 152);
    const auto *fpf1_155 = buffer.data(fpf1 + 155);
    const auto *fpf1_156 = buffer.data(fpf1 + 156);
    const auto *fpf1_157 = buffer.data(fpf1 + 157);
    const auto *fpf1_158 = buffer.data(fpf1 + 158);
    const auto *fpf1_159 = buffer.data(fpf1 + 159);
    const auto *fpf1_170 = buffer.data(fpf1 + 170);
    const auto *fpf1_171 = buffer.data(fpf1 + 171);
    const auto *fpf1_172 = buffer.data(fpf1 + 172);

    const auto *fpg_179 = buffer.data(fpg + 179);
    const auto *fpg_180 = buffer.data(fpg + 180);
    const auto *fpg_182 = buffer.data(fpg + 182);
    const auto *fpg_183 = buffer.data(fpg + 183);
    const auto *fpg_185 = buffer.data(fpg + 185);
    const auto *fpg_190 = buffer.data(fpg + 190);
    const auto *fpg_191 = buffer.data(fpg + 191);
    const auto *fpg_192 = buffer.data(fpg + 192);
    const auto *fpg_193 = buffer.data(fpg + 193);
    const auto *fpg_194 = buffer.data(fpg + 194);
    const auto *fpg_195 = buffer.data(fpg + 195);
    const auto *fpg_197 = buffer.data(fpg + 197);
    const auto *fpg_198 = buffer.data(fpg + 198);
    const auto *fpg_200 = buffer.data(fpg + 200);
    const auto *fpg_204 = buffer.data(fpg + 204);
    const auto *fpg_205 = buffer.data(fpg + 205);
    const auto *fpg_206 = buffer.data(fpg + 206);
    const auto *fpg_207 = buffer.data(fpg + 207);
    const auto *fpg_208 = buffer.data(fpg + 208);
    const auto *fpg_209 = buffer.data(fpg + 209);
    const auto *fpg_210 = buffer.data(fpg + 210);
    const auto *fpg_212 = buffer.data(fpg + 212);
    const auto *fpg_213 = buffer.data(fpg + 213);
    const auto *fpg_215 = buffer.data(fpg + 215);
    const auto *fpg_216 = buffer.data(fpg + 216);
    const auto *fpg_220 = buffer.data(fpg + 220);
    const auto *fpg_221 = buffer.data(fpg + 221);
    const auto *fpg_222 = buffer.data(fpg + 222);
    const auto *fpg_223 = buffer.data(fpg + 223);
    const auto *fpg_224 = buffer.data(fpg + 224);
    const auto *fpg_225 = buffer.data(fpg + 225);
    const auto *fpg_226 = buffer.data(fpg + 226);
    const auto *fpg_227 = buffer.data(fpg + 227);
    const auto *fpg_228 = buffer.data(fpg + 228);
    const auto *fpg_229 = buffer.data(fpg + 229);
    const auto *fpg_230 = buffer.data(fpg + 230);
    const auto *fpg_234 = buffer.data(fpg + 234);
    const auto *fpg_235 = buffer.data(fpg + 235);
    const auto *fpg_236 = buffer.data(fpg + 236);
    const auto *fpg_237 = buffer.data(fpg + 237);
    const auto *fpg_238 = buffer.data(fpg + 238);
    const auto *fpg_239 = buffer.data(fpg + 239);
    const auto *fpg_240 = buffer.data(fpg + 240);
    const auto *fpg_242 = buffer.data(fpg + 242);
    const auto *fpg_245 = buffer.data(fpg + 245);
    const auto *fpg_249 = buffer.data(fpg + 249);
    const auto *fpg_250 = buffer.data(fpg + 250);
    const auto *fpg_251 = buffer.data(fpg + 251);
    const auto *fpg_252 = buffer.data(fpg + 252);
    const auto *fpg_254 = buffer.data(fpg + 254);
    const auto *fpg_255 = buffer.data(fpg + 255);
    const auto *fpg_256 = buffer.data(fpg + 256);
    const auto *fpg_257 = buffer.data(fpg + 257);
    const auto *fpg_258 = buffer.data(fpg + 258);
    const auto *fpg_259 = buffer.data(fpg + 259);
    const auto *fpg_260 = buffer.data(fpg + 260);

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pa_x, pa_y, pc_x, pc_y, dph0_126, \
                         dph0_249, dph0_251, dpg_89, dph1_126, dph1_249, dph1_251, \
                         fpg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = pa_x[k] * dph0_249[k]
                   - f_9 * pc_x[k] * dph1_249[k];

        t_250[k] = f_10 * dpg_89[k]
                   + f_4 * pc_y[k] * fpg_179[k];

        t_251[k] = pa_x[k] * dph0_251[k]
                   - f_9 * pc_x[k] * dph1_251[k];

        t_252[k] = pa_y[k] * dph0_126[k]
                   - f_9 * pc_y[k] * dph1_126[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pa_z, pc_y, pc_z, dph0_66, dpg_45, \
                         dpg_90, dpg_92, dph1_66, fpg_180, fpg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_1 * dpg_90[k]
                   + f_4 * pc_y[k] * fpg_180[k];

        t_254[k] = f_1 * dpg_45[k]
                   + f_4 * pc_z[k] * fpg_180[k];

        t_255[k] = pa_z[k] * dph0_66[k]
                   - f_9 * pc_z[k] * dph1_66[k];

        t_256[k] = f_1 * dpg_92[k]
                   + f_4 * pc_y[k] * fpg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pa_y, pa_z, pc_y, pc_z, dph0_69, \
                         dph0_131, dpg_48, dpg_95, dph1_69, dph1_131, fpg_183, \
                         fpg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = pa_y[k] * dph0_131[k]
                   - f_9 * pc_y[k] * dph1_131[k];

        t_258[k] = pa_z[k] * dph0_69[k]
                   - f_9 * pc_z[k] * dph1_69[k];

        t_259[k] = f_1 * dpg_48[k]
                   + f_4 * pc_z[k] * fpg_183[k];

        t_260[k] = f_1 * dpg_95[k]
                   + f_4 * pc_y[k] * fpg_185[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pa_y, pa_z, pc_x, pc_y, pc_z, dph0_73, dph0_135, \
                         dpg_191, dph1_73, dph1_135, fsg_71, fpg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pa_y[k] * dph0_135[k]
                   - f_9 * pc_y[k] * dph1_135[k];

        t_262[k] = pa_z[k] * dph0_73[k]
                   - f_9 * pc_z[k] * dph1_73[k];

        t_263[k] = f_1 * dpg_191[k]
                   + f_1 * fsg_71[k]
                   + f_4 * pc_x[k] * fpg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pa_y, pc_x, pc_y, dph0_140, dpg_192, dpg_193, \
                         dph1_140, fsg_72, fsg_73, fpg_192, fpg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_1 * dpg_192[k]
                   + f_1 * fsg_72[k]
                   + f_4 * pc_x[k] * fpg_192[k];

        t_265[k] = f_1 * dpg_193[k]
                   + f_1 * fsg_73[k]
                   + f_4 * pc_x[k] * fpg_193[k];

        t_266[k] = pa_y[k] * dph0_140[k]
                   - f_9 * pc_y[k] * dph1_140[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, dpg_55, dpg_100, dpg_102, fpf0_126, \
                         fpf0_128, fpf1_126, fpf1_128, fpg_190, \
                         fpg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_1 * dpg_100[k]
                   + f_2 * fpf0_126[k]
                   - f_3 * fpf1_126[k]
                   + f_4 * pc_y[k] * fpg_190[k];

        t_268[k] = f_1 * dpg_55[k]
                   + f_4 * pc_z[k] * fpg_190[k];

        t_269[k] = f_1 * dpg_102[k]
                   + f_7 * fpf0_128[k]
                   - f_8 * fpf1_128[k]
                   + f_4 * pc_y[k] * fpg_192[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pc_y, pc_z, dpg_59, dpg_103, dpg_104, fpf0_129, \
                         fpf1_129, fpg_193, fpg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * dpg_103[k]
                   + f_5 * fpf0_129[k]
                   - f_6 * fpf1_129[k]
                   + f_4 * pc_y[k] * fpg_193[k];

        t_271[k] = f_1 * dpg_104[k]
                   + f_4 * pc_y[k] * fpg_194[k];

        t_272[k] = f_1 * dpg_59[k]
                   + f_2 * fpf0_129[k]
                   - f_3 * fpf1_129[k]
                   + f_4 * pc_z[k] * fpg_194[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, pa_z, pc_x, pc_z, dph0_85, dph0_87, \
                         dpg_60, dpg_195, dph1_85, dph1_87, fpf0_130, fpf1_130, \
                         fpg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * dpg_195[k]
                   + f_2 * fpf0_130[k]
                   - f_3 * fpf1_130[k]
                   + f_4 * pc_x[k] * fpg_195[k];

        t_274[k] = pa_z[k] * dph0_85[k]
                   - f_9 * pc_z[k] * dph1_85[k];

        t_275[k] = f_1 * dpg_60[k]
                   + f_4 * pc_z[k] * fpg_195[k];

        t_276[k] = pa_z[k] * dph0_87[k]
                   - f_9 * pc_z[k] * dph1_87[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, pa_z, pc_x, pc_y, pc_z, dph0_90, dpg_107, \
                         dpg_200, dph1_90, fsg_62, fpf0_135, fpf1_135, fpg_197, \
                         fpg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_1 * dpg_107[k]
                   + f_1 * fsg_62[k]
                   + f_4 * pc_y[k] * fpg_197[k];

        t_278[k] = f_1 * dpg_200[k]
                   + f_7 * fpf0_135[k]
                   - f_8 * fpf1_135[k]
                   + f_4 * pc_x[k] * fpg_200[k];

        t_279[k] = pa_z[k] * dph0_90[k]
                   - f_9 * pc_z[k] * dph1_90[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pc_x, pc_y, pc_z, dpg_63, dpg_110, dpg_204, \
                         fsg_65, fpf0_139, fpf1_139, fpg_198, fpg_200, \
                         fpg_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_1 * dpg_63[k]
                   + f_4 * pc_z[k] * fpg_198[k];

        t_281[k] = f_1 * dpg_110[k]
                   + f_1 * fsg_65[k]
                   + f_4 * pc_y[k] * fpg_200[k];

        t_282[k] = f_1 * dpg_204[k]
                   + f_5 * fpf0_139[k]
                   - f_6 * fpf1_139[k]
                   + f_4 * pc_x[k] * fpg_204[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, dpg_205, dpg_206, dpg_207, \
                         dpg_208, dpg_209, fpg_205, fpg_206, fpg_207, fpg_208, \
                         fpg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_1 * dpg_205[k]
                   + f_4 * pc_x[k] * fpg_205[k];

        t_284[k] = f_1 * dpg_206[k]
                   + f_4 * pc_x[k] * fpg_206[k];

        t_285[k] = f_1 * dpg_207[k]
                   + f_4 * pc_x[k] * fpg_207[k];

        t_286[k] = f_1 * dpg_208[k]
                   + f_4 * pc_x[k] * fpg_208[k];

        t_287[k] = f_1 * dpg_209[k]
                   + f_4 * pc_x[k] * fpg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pa_x, pc_x, pc_z, dph0_288, dph0_290, \
                         dph0_291, dpg_70, dph1_288, dph1_290, dph1_291, \
                         fpg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = pa_x[k] * dph0_288[k]
                   - f_9 * pc_x[k] * dph1_288[k];

        t_289[k] = f_1 * dpg_70[k]
                   + f_4 * pc_z[k] * fpg_205[k];

        t_290[k] = pa_x[k] * dph0_290[k]
                   - f_9 * pc_x[k] * dph1_290[k];

        t_291[k] = pa_x[k] * dph0_291[k]
                   - f_9 * pc_x[k] * dph1_291[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_x, pa_y, pc_x, pc_y, dph0_168, \
                         dph0_292, dph0_293, dpg_120, dph1_168, dph1_292, dph1_293, \
                         fpg_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = pa_x[k] * dph0_292[k]
                   - f_9 * pc_x[k] * dph1_292[k];

        t_293[k] = pa_x[k] * dph0_293[k]
                   - f_9 * pc_x[k] * dph1_293[k];

        t_294[k] = pa_y[k] * dph0_168[k]
                   - f_9 * pc_y[k] * dph1_168[k];

        t_295[k] = f_1 * dpg_120[k]
                   + f_4 * pc_y[k] * fpg_210[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, pa_y, pc_x, pc_y, dph0_170, dpg_122, dpg_213, \
                         dph1_170, fpf0_143, fpf1_143, fpg_212, \
                         fpg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = pa_y[k] * dph0_170[k]
                   - f_9 * pc_y[k] * dph1_170[k];

        t_297[k] = f_1 * dpg_213[k]
                   + f_7 * fpf0_143[k]
                   - f_8 * fpf1_143[k]
                   + f_4 * pc_x[k] * fpg_213[k];

        t_298[k] = f_1 * dpg_122[k]
                   + f_4 * pc_y[k] * fpg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pa_y, pc_x, pc_y, pc_z, dph0_173, dpg_78, \
                         dpg_216, dph1_173, fsg_63, fpf0_146, fpf1_146, fpg_213, \
                         fpg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pa_y[k] * dph0_173[k]
                   - f_9 * pc_y[k] * dph1_173[k];

        t_300[k] = f_1 * dpg_216[k]
                   + f_5 * fpf0_146[k]
                   - f_6 * fpf1_146[k]
                   + f_4 * pc_x[k] * fpg_216[k];

        t_301[k] = f_1 * dpg_78[k]
                   + f_1 * fsg_63[k]
                   + f_4 * pc_z[k] * fpg_213[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pa_y, pc_x, pc_y, dph0_177, dpg_125, \
                         dpg_220, dpg_221, dph1_177, fpg_215, fpg_220, \
                         fpg_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_1 * dpg_125[k]
                   + f_4 * pc_y[k] * fpg_215[k];

        t_303[k] = pa_y[k] * dph0_177[k]
                   - f_9 * pc_y[k] * dph1_177[k];

        t_304[k] = f_1 * dpg_220[k]
                   + f_4 * pc_x[k] * fpg_220[k];

        t_305[k] = f_1 * dpg_221[k]
                   + f_4 * pc_x[k] * fpg_221[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pa_x, pc_x, dph0_309, dpg_222, dpg_223, \
                         dpg_224, dph1_309, fpg_222, fpg_223, fpg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_1 * dpg_222[k]
                   + f_4 * pc_x[k] * fpg_222[k];

        t_307[k] = f_1 * dpg_223[k]
                   + f_4 * pc_x[k] * fpg_223[k];

        t_308[k] = f_1 * dpg_224[k]
                   + f_4 * pc_x[k] * fpg_224[k];

        t_309[k] = pa_x[k] * dph0_309[k]
                   - f_9 * pc_x[k] * dph1_309[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pa_x, pc_x, pc_y, dph0_310, dph0_311, \
                         dph0_312, dpg_134, dph1_310, dph1_311, dph1_312, \
                         fpg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = pa_x[k] * dph0_310[k]
                   - f_9 * pc_x[k] * dph1_310[k];

        t_311[k] = pa_x[k] * dph0_311[k]
                   - f_9 * pc_x[k] * dph1_311[k];

        t_312[k] = pa_x[k] * dph0_312[k]
                   - f_9 * pc_x[k] * dph1_312[k];

        t_313[k] = f_1 * dpg_134[k]
                   + f_4 * pc_y[k] * fpg_224[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pa_x, pc_x, pc_y, pc_z, dph0_314, dpg_90, \
                         dpg_225, dph1_314, fsg_75, fpf0_150, fpf1_150, \
                         fpg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pa_x[k] * dph0_314[k]
                   - f_9 * pc_x[k] * dph1_314[k];

        t_315[k] = f_1 * dpg_225[k]
                   + f_1 * fsg_75[k]
                   + f_2 * fpf0_150[k]
                   - f_3 * fpf1_150[k]
                   + f_4 * pc_x[k] * fpg_225[k];

        t_316[k] = f_4 * pc_y[k] * fpg_225[k];

        t_317[k] = f_10 * dpg_90[k]
                   + f_4 * pc_z[k] * fpg_225[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_x, pc_y, dpg_230, fsg_80, fpf0_150, fpf0_155, \
                         fpf1_150, fpf1_155, fpg_226, fpg_227, \
                         fpg_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_5 * fpf0_150[k]
                   - f_6 * fpf1_150[k]
                   + f_4 * pc_y[k] * fpg_226[k];

        t_319[k] = f_4 * pc_y[k] * fpg_227[k];

        t_320[k] = f_1 * dpg_230[k]
                   + f_1 * fsg_80[k]
                   + f_7 * fpf0_155[k]
                   - f_8 * fpf1_155[k]
                   + f_4 * pc_x[k] * fpg_230[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pc_y, fpf0_151, fpf0_152, fpf1_151, fpf1_152, \
                         fpg_228, fpg_229, fpg_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_7 * fpf0_151[k]
                   - f_8 * fpf1_151[k]
                   + f_4 * pc_y[k] * fpg_228[k];

        t_322[k] = f_5 * fpf0_152[k]
                   - f_6 * fpf1_152[k]
                   + f_4 * pc_y[k] * fpg_229[k];

        t_323[k] = f_4 * pc_y[k] * fpg_230[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, pc_x, dpg_234, dpg_235, dpg_236, fsg_84, fsg_85, \
                         fsg_86, fpf0_159, fpf1_159, fpg_234, fpg_235, \
                         fpg_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_1 * dpg_234[k]
                   + f_1 * fsg_84[k]
                   + f_5 * fpf0_159[k]
                   - f_6 * fpf1_159[k]
                   + f_4 * pc_x[k] * fpg_234[k];

        t_325[k] = f_1 * dpg_235[k]
                   + f_1 * fsg_85[k]
                   + f_4 * pc_x[k] * fpg_235[k];

        t_326[k] = f_1 * dpg_236[k]
                   + f_1 * fsg_86[k]
                   + f_4 * pc_x[k] * fpg_236[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, dpg_237, dpg_239, fsg_87, \
                         fsg_89, fpf0_156, fpf1_156, fpg_234, fpg_235, fpg_237, \
                         fpg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_1 * dpg_237[k]
                   + f_1 * fsg_87[k]
                   + f_4 * pc_x[k] * fpg_237[k];

        t_328[k] = f_4 * pc_y[k] * fpg_234[k];

        t_329[k] = f_1 * dpg_239[k]
                   + f_1 * fsg_89[k]
                   + f_4 * pc_x[k] * fpg_239[k];

        t_330[k] = f_2 * fpf0_156[k]
                   - f_3 * fpf1_156[k]
                   + f_4 * pc_y[k] * fpg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pc_y, fpf0_157, fpf0_158, fpf0_159, \
                         fpf1_157, fpf1_158, fpf1_159, fpg_236, fpg_237, fpg_238, \
                         fpg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_14 * fpf0_157[k]
                   - f_15 * fpf1_157[k]
                   + f_4 * pc_y[k] * fpg_236[k];

        t_332[k] = f_7 * fpf0_158[k]
                   - f_8 * fpf1_158[k]
                   + f_4 * pc_y[k] * fpg_237[k];

        t_333[k] = f_5 * fpf0_159[k]
                   - f_6 * fpf1_159[k]
                   + f_4 * pc_y[k] * fpg_238[k];

        t_334[k] = f_4 * pc_y[k] * fpg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pb_y, pc_y, pc_z, dpg_104, dpg_105, \
                         fsh0_105, fsg_75, fsh1_105, fpf0_159, fpf1_159, fpg_239, \
                         fpg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_10 * dpg_104[k]
                   + f_2 * fpf0_159[k]
                   - f_3 * fpf1_159[k]
                   + f_4 * pc_z[k] * fpg_239[k];

        t_336[k] = pb_y[k] * fsh0_105[k]
                   - f_9 * pc_y[k] * fsh1_105[k];

        t_337[k] = f_1 * fsg_75[k]
                   + f_4 * pc_y[k] * fpg_240[k];

        t_338[k] = f_10 * dpg_105[k]
                   + f_4 * pc_z[k] * fpg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_x, pb_y, pc_x, pc_y, dph0_339, dpg_243, \
                         dph1_339, fsh0_110, fsg_77, fsh1_110, \
                         fpg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_x[k] * dph0_339[k]
                   + f_0 * dpg_243[k]
                   - f_9 * pc_x[k] * dph1_339[k];

        t_340[k] = f_1 * fsg_77[k]
                   + f_4 * pc_y[k] * fpg_242[k];

        t_341[k] = pb_y[k] * fsh0_110[k]
                   - f_9 * pc_y[k] * fsh1_110[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_x, pc_x, pc_y, dph0_342, dph0_343, dpg_246, \
                         dpg_247, dph1_342, dph1_343, fsg_80, fpg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_x[k] * dph0_342[k]
                   + f_10 * dpg_246[k]
                   - f_9 * pc_x[k] * dph1_342[k];

        t_343[k] = pa_x[k] * dph0_343[k]
                   + f_10 * dpg_247[k]
                   - f_9 * pc_x[k] * dph1_343[k];

        t_344[k] = f_1 * fsg_80[k]
                   + f_4 * pc_y[k] * fpg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pb_y, pc_x, pc_y, dpg_250, dpg_251, \
                         dpg_252, fsh0_114, fsh1_114, fpg_250, fpg_251, \
                         fpg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = pb_y[k] * fsh0_114[k]
                   - f_9 * pc_y[k] * fsh1_114[k];

        t_346[k] = f_1 * dpg_250[k]
                   + f_4 * pc_x[k] * fpg_250[k];

        t_347[k] = f_1 * dpg_251[k]
                   + f_4 * pc_x[k] * fpg_251[k];

        t_348[k] = f_1 * dpg_252[k]
                   + f_4 * pc_x[k] * fpg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_x, pc_x, pc_y, dph0_351, dph0_352, \
                         dpg_254, dph1_351, dph1_352, fsg_84, fpg_249, \
                         fpg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_1 * fsg_84[k]
                   + f_4 * pc_y[k] * fpg_249[k];

        t_350[k] = f_1 * dpg_254[k]
                   + f_4 * pc_x[k] * fpg_254[k];

        t_351[k] = pa_x[k] * dph0_351[k]
                   - f_9 * pc_x[k] * dph1_351[k];

        t_352[k] = pa_x[k] * dph0_352[k]
                   - f_9 * pc_x[k] * dph1_352[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, pa_x, pc_x, pc_y, dph0_353, dph0_354, \
                         dph0_356, dph1_353, dph1_354, dph1_356, fsg_89, \
                         fpg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = pa_x[k] * dph0_353[k]
                   - f_9 * pc_x[k] * dph1_353[k];

        t_354[k] = pa_x[k] * dph0_354[k]
                   - f_9 * pc_x[k] * dph1_354[k];

        t_355[k] = f_1 * fsg_89[k]
                   + f_4 * pc_y[k] * fpg_254[k];

        t_356[k] = pa_x[k] * dph0_356[k]
                   - f_9 * pc_x[k] * dph1_356[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pa_x, pc_x, pc_y, pc_z, dph0_357, dpg_120, \
                         dpg_255, dph1_357, fsg_75, fpg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pa_x[k] * dph0_357[k]
                   + f_11 * dpg_255[k]
                   - f_9 * pc_x[k] * dph1_357[k];

        t_358[k] = f_4 * pc_y[k] * fpg_255[k];

        t_359[k] = f_10 * dpg_120[k]
                   + f_1 * fsg_75[k]
                   + f_4 * pc_z[k] * fpg_255[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_x, pc_x, pc_y, dph0_362, dpg_260, dph1_362, \
                         fpf0_170, fpf1_170, fpg_256, fpg_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_5 * fpf0_170[k]
                   - f_6 * fpf1_170[k]
                   + f_4 * pc_y[k] * fpg_256[k];

        t_361[k] = f_4 * pc_y[k] * fpg_257[k];

        t_362[k] = pa_x[k] * dph0_362[k]
                   + f_0 * dpg_260[k]
                   - f_9 * pc_x[k] * dph1_362[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pc_y, fpf0_171, fpf0_172, fpf1_171, fpf1_172, \
                         fpg_258, fpg_259, fpg_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_7 * fpf0_171[k]
                   - f_8 * fpf1_171[k]
                   + f_4 * pc_y[k] * fpg_258[k];

        t_364[k] = f_5 * fpf0_172[k]
                   - f_6 * fpf1_172[k]
                   + f_4 * pc_y[k] * fpg_259[k];

        t_365[k] = f_4 * pc_y[k] * fpg_260[k];
    }
}

static auto
compute_prim_fph_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dph0, const size_t dpg,
                                                          const size_t dph1, const size_t fsh0,
                                                          const size_t fsg, const size_t fsh1,
                                                          const size_t fpf0, const size_t fpf1,
                                                          const size_t fpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_14 = 1.5 / gamma;
    const auto f_15 = 1.5 * p / (gamma * q);
    const auto f_16 = 2.0 / q;

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
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dph0_189 = buffer.data(dph0 + 189);
    const auto *dph0_190 = buffer.data(dph0 + 190);
    const auto *dph0_192 = buffer.data(dph0 + 192);
    const auto *dph0_193 = buffer.data(dph0 + 193);
    const auto *dph0_195 = buffer.data(dph0 + 195);
    const auto *dph0_196 = buffer.data(dph0 + 196);
    const auto *dph0_197 = buffer.data(dph0 + 197);
    const auto *dph0_204 = buffer.data(dph0 + 204);
    const auto *dph0_210 = buffer.data(dph0 + 210);
    const auto *dph0_211 = buffer.data(dph0 + 211);
    const auto *dph0_213 = buffer.data(dph0 + 213);
    const auto *dph0_216 = buffer.data(dph0 + 216);
    const auto *dph0_225 = buffer.data(dph0 + 225);
    const auto *dph0_227 = buffer.data(dph0 + 227);
    const auto *dph0_228 = buffer.data(dph0 + 228);
    const auto *dph0_366 = buffer.data(dph0 + 366);
    const auto *dph0_372 = buffer.data(dph0 + 372);
    const auto *dph0_373 = buffer.data(dph0 + 373);
    const auto *dph0_374 = buffer.data(dph0 + 374);
    const auto *dph0_375 = buffer.data(dph0 + 375);
    const auto *dph0_377 = buffer.data(dph0 + 377);

    const auto *dpg_136 = buffer.data(dpg + 136);
    const auto *dpg_138 = buffer.data(dpg + 138);
    const auto *dpg_139 = buffer.data(dpg + 139);
    const auto *dpg_145 = buffer.data(dpg + 145);
    const auto *dpg_149 = buffer.data(dpg + 149);
    const auto *dpg_160 = buffer.data(dpg + 160);
    const auto *dpg_161 = buffer.data(dpg + 161);
    const auto *dpg_162 = buffer.data(dpg + 162);
    const auto *dpg_164 = buffer.data(dpg + 164);
    const auto *dpg_179 = buffer.data(dpg + 179);
    const auto *dpg_194 = buffer.data(dpg + 194);
    const auto *dpg_209 = buffer.data(dpg + 209);
    const auto *dpg_264 = buffer.data(dpg + 264);
    const auto *dpg_265 = buffer.data(dpg + 265);
    const auto *dpg_266 = buffer.data(dpg + 266);
    const auto *dpg_267 = buffer.data(dpg + 267);
    const auto *dpg_269 = buffer.data(dpg + 269);

    const auto *dph1_189 = buffer.data(dph1 + 189);
    const auto *dph1_190 = buffer.data(dph1 + 190);
    const auto *dph1_192 = buffer.data(dph1 + 192);
    const auto *dph1_193 = buffer.data(dph1 + 193);
    const auto *dph1_195 = buffer.data(dph1 + 195);
    const auto *dph1_196 = buffer.data(dph1 + 196);
    const auto *dph1_197 = buffer.data(dph1 + 197);
    const auto *dph1_204 = buffer.data(dph1 + 204);
    const auto *dph1_210 = buffer.data(dph1 + 210);
    const auto *dph1_211 = buffer.data(dph1 + 211);
    const auto *dph1_213 = buffer.data(dph1 + 213);
    const auto *dph1_216 = buffer.data(dph1 + 216);
    const auto *dph1_225 = buffer.data(dph1 + 225);
    const auto *dph1_227 = buffer.data(dph1 + 227);
    const auto *dph1_228 = buffer.data(dph1 + 228);
    const auto *dph1_366 = buffer.data(dph1 + 366);
    const auto *dph1_372 = buffer.data(dph1 + 372);
    const auto *dph1_373 = buffer.data(dph1 + 373);
    const auto *dph1_374 = buffer.data(dph1 + 374);
    const auto *dph1_375 = buffer.data(dph1 + 375);
    const auto *dph1_377 = buffer.data(dph1 + 377);

    const auto *fsh0_126 = buffer.data(fsh0 + 126);
    const auto *fsh0_127 = buffer.data(fsh0 + 127);
    const auto *fsh0_129 = buffer.data(fsh0 + 129);
    const auto *fsh0_131 = buffer.data(fsh0 + 131);
    const auto *fsh0_132 = buffer.data(fsh0 + 132);
    const auto *fsh0_134 = buffer.data(fsh0 + 134);
    const auto *fsh0_135 = buffer.data(fsh0 + 135);
    const auto *fsh0_141 = buffer.data(fsh0 + 141);
    const auto *fsh0_143 = buffer.data(fsh0 + 143);
    const auto *fsh0_144 = buffer.data(fsh0 + 144);
    const auto *fsh0_146 = buffer.data(fsh0 + 146);
    const auto *fsh0_149 = buffer.data(fsh0 + 149);
    const auto *fsh0_152 = buffer.data(fsh0 + 152);
    const auto *fsh0_156 = buffer.data(fsh0 + 156);
    const auto *fsh0_164 = buffer.data(fsh0 + 164);
    const auto *fsh0_165 = buffer.data(fsh0 + 165);
    const auto *fsh0_167 = buffer.data(fsh0 + 167);

    const auto *fsg_90 = buffer.data(fsg + 90);
    const auto *fsg_91 = buffer.data(fsg + 91);
    const auto *fsg_93 = buffer.data(fsg + 93);
    const auto *fsg_95 = buffer.data(fsg + 95);
    const auto *fsg_96 = buffer.data(fsg + 96);
    const auto *fsg_98 = buffer.data(fsg + 98);
    const auto *fsg_99 = buffer.data(fsg + 99);
    const auto *fsg_100 = buffer.data(fsg + 100);
    const auto *fsg_101 = buffer.data(fsg + 101);
    const auto *fsg_102 = buffer.data(fsg + 102);
    const auto *fsg_103 = buffer.data(fsg + 103);
    const auto *fsg_104 = buffer.data(fsg + 104);
    const auto *fsg_107 = buffer.data(fsg + 107);
    const auto *fsg_110 = buffer.data(fsg + 110);
    const auto *fsg_114 = buffer.data(fsg + 114);
    const auto *fsg_115 = buffer.data(fsg + 115);
    const auto *fsg_116 = buffer.data(fsg + 116);
    const auto *fsg_117 = buffer.data(fsg + 117);
    const auto *fsg_118 = buffer.data(fsg + 118);
    const auto *fsg_119 = buffer.data(fsg + 119);

    const auto *fsh1_126 = buffer.data(fsh1 + 126);
    const auto *fsh1_127 = buffer.data(fsh1 + 127);
    const auto *fsh1_129 = buffer.data(fsh1 + 129);
    const auto *fsh1_131 = buffer.data(fsh1 + 131);
    const auto *fsh1_132 = buffer.data(fsh1 + 132);
    const auto *fsh1_134 = buffer.data(fsh1 + 134);
    const auto *fsh1_135 = buffer.data(fsh1 + 135);
    const auto *fsh1_141 = buffer.data(fsh1 + 141);
    const auto *fsh1_143 = buffer.data(fsh1 + 143);
    const auto *fsh1_144 = buffer.data(fsh1 + 144);
    const auto *fsh1_146 = buffer.data(fsh1 + 146);
    const auto *fsh1_149 = buffer.data(fsh1 + 149);
    const auto *fsh1_152 = buffer.data(fsh1 + 152);
    const auto *fsh1_156 = buffer.data(fsh1 + 156);
    const auto *fsh1_164 = buffer.data(fsh1 + 164);
    const auto *fsh1_165 = buffer.data(fsh1 + 165);
    const auto *fsh1_167 = buffer.data(fsh1 + 167);

    const auto *fpf0_190 = buffer.data(fpf0 + 190);
    const auto *fpf0_191 = buffer.data(fpf0 + 191);
    const auto *fpf0_193 = buffer.data(fpf0 + 193);
    const auto *fpf0_195 = buffer.data(fpf0 + 195);
    const auto *fpf0_196 = buffer.data(fpf0 + 196);
    const auto *fpf0_197 = buffer.data(fpf0 + 197);
    const auto *fpf0_198 = buffer.data(fpf0 + 198);
    const auto *fpf0_199 = buffer.data(fpf0 + 199);
    const auto *fpf0_205 = buffer.data(fpf0 + 205);
    const auto *fpf0_208 = buffer.data(fpf0 + 208);
    const auto *fpf0_209 = buffer.data(fpf0 + 209);
    const auto *fpf0_222 = buffer.data(fpf0 + 222);
    const auto *fpf0_224 = buffer.data(fpf0 + 224);
    const auto *fpf0_225 = buffer.data(fpf0 + 225);
    const auto *fpf0_227 = buffer.data(fpf0 + 227);
    const auto *fpf0_228 = buffer.data(fpf0 + 228);
    const auto *fpf0_229 = buffer.data(fpf0 + 229);
    const auto *fpf0_230 = buffer.data(fpf0 + 230);
    const auto *fpf0_231 = buffer.data(fpf0 + 231);
    const auto *fpf0_232 = buffer.data(fpf0 + 232);
    const auto *fpf0_233 = buffer.data(fpf0 + 233);
    const auto *fpf0_234 = buffer.data(fpf0 + 234);
    const auto *fpf0_235 = buffer.data(fpf0 + 235);
    const auto *fpf0_236 = buffer.data(fpf0 + 236);

    const auto *fpf1_190 = buffer.data(fpf1 + 190);
    const auto *fpf1_191 = buffer.data(fpf1 + 191);
    const auto *fpf1_193 = buffer.data(fpf1 + 193);
    const auto *fpf1_195 = buffer.data(fpf1 + 195);
    const auto *fpf1_196 = buffer.data(fpf1 + 196);
    const auto *fpf1_197 = buffer.data(fpf1 + 197);
    const auto *fpf1_198 = buffer.data(fpf1 + 198);
    const auto *fpf1_199 = buffer.data(fpf1 + 199);
    const auto *fpf1_205 = buffer.data(fpf1 + 205);
    const auto *fpf1_208 = buffer.data(fpf1 + 208);
    const auto *fpf1_209 = buffer.data(fpf1 + 209);
    const auto *fpf1_222 = buffer.data(fpf1 + 222);
    const auto *fpf1_224 = buffer.data(fpf1 + 224);
    const auto *fpf1_225 = buffer.data(fpf1 + 225);
    const auto *fpf1_227 = buffer.data(fpf1 + 227);
    const auto *fpf1_228 = buffer.data(fpf1 + 228);
    const auto *fpf1_229 = buffer.data(fpf1 + 229);
    const auto *fpf1_230 = buffer.data(fpf1 + 230);
    const auto *fpf1_231 = buffer.data(fpf1 + 231);
    const auto *fpf1_232 = buffer.data(fpf1 + 232);
    const auto *fpf1_233 = buffer.data(fpf1 + 233);
    const auto *fpf1_234 = buffer.data(fpf1 + 234);
    const auto *fpf1_235 = buffer.data(fpf1 + 235);
    const auto *fpf1_236 = buffer.data(fpf1 + 236);

    const auto *fpg_264 = buffer.data(fpg + 264);
    const auto *fpg_265 = buffer.data(fpg + 265);
    const auto *fpg_266 = buffer.data(fpg + 266);
    const auto *fpg_267 = buffer.data(fpg + 267);
    const auto *fpg_269 = buffer.data(fpg + 269);
    const auto *fpg_270 = buffer.data(fpg + 270);
    const auto *fpg_271 = buffer.data(fpg + 271);
    const auto *fpg_273 = buffer.data(fpg + 273);
    const auto *fpg_280 = buffer.data(fpg + 280);
    const auto *fpg_281 = buffer.data(fpg + 281);
    const auto *fpg_282 = buffer.data(fpg + 282);
    const auto *fpg_283 = buffer.data(fpg + 283);
    const auto *fpg_284 = buffer.data(fpg + 284);
    const auto *fpg_285 = buffer.data(fpg + 285);
    const auto *fpg_286 = buffer.data(fpg + 286);
    const auto *fpg_288 = buffer.data(fpg + 288);
    const auto *fpg_290 = buffer.data(fpg + 290);
    const auto *fpg_291 = buffer.data(fpg + 291);
    const auto *fpg_293 = buffer.data(fpg + 293);
    const auto *fpg_294 = buffer.data(fpg + 294);
    const auto *fpg_295 = buffer.data(fpg + 295);
    const auto *fpg_296 = buffer.data(fpg + 296);
    const auto *fpg_297 = buffer.data(fpg + 297);
    const auto *fpg_298 = buffer.data(fpg + 298);
    const auto *fpg_299 = buffer.data(fpg + 299);
    const auto *fpg_300 = buffer.data(fpg + 300);
    const auto *fpg_301 = buffer.data(fpg + 301);
    const auto *fpg_303 = buffer.data(fpg + 303);
    const auto *fpg_305 = buffer.data(fpg + 305);
    const auto *fpg_308 = buffer.data(fpg + 308);
    const auto *fpg_309 = buffer.data(fpg + 309);
    const auto *fpg_310 = buffer.data(fpg + 310);
    const auto *fpg_311 = buffer.data(fpg + 311);
    const auto *fpg_312 = buffer.data(fpg + 312);
    const auto *fpg_313 = buffer.data(fpg + 313);
    const auto *fpg_314 = buffer.data(fpg + 314);
    const auto *fpg_325 = buffer.data(fpg + 325);
    const auto *fpg_326 = buffer.data(fpg + 326);
    const auto *fpg_327 = buffer.data(fpg + 327);
    const auto *fpg_328 = buffer.data(fpg + 328);
    const auto *fpg_329 = buffer.data(fpg + 329);
    const auto *fpg_332 = buffer.data(fpg + 332);
    const auto *fpg_334 = buffer.data(fpg + 334);
    const auto *fpg_335 = buffer.data(fpg + 335);
    const auto *fpg_337 = buffer.data(fpg + 337);
    const auto *fpg_338 = buffer.data(fpg + 338);
    const auto *fpg_339 = buffer.data(fpg + 339);
    const auto *fpg_340 = buffer.data(fpg + 340);
    const auto *fpg_341 = buffer.data(fpg + 341);
    const auto *fpg_342 = buffer.data(fpg + 342);
    const auto *fpg_343 = buffer.data(fpg + 343);
    const auto *fpg_344 = buffer.data(fpg + 344);
    const auto *fpg_345 = buffer.data(fpg + 345);
    const auto *fpg_346 = buffer.data(fpg + 346);
    const auto *fpg_347 = buffer.data(fpg + 347);
    const auto *fpg_348 = buffer.data(fpg + 348);
    const auto *fpg_349 = buffer.data(fpg + 349);
    const auto *fpg_350 = buffer.data(fpg + 350);
    const auto *fpg_351 = buffer.data(fpg + 351);

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_x, pc_x, dph0_366, dpg_264, dpg_265, \
                         dpg_266, dpg_267, dph1_366, fpg_265, fpg_266, \
                         fpg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = pa_x[k] * dph0_366[k]
                   + f_10 * dpg_264[k]
                   - f_9 * pc_x[k] * dph1_366[k];

        t_367[k] = f_1 * dpg_265[k]
                   + f_4 * pc_x[k] * fpg_265[k];

        t_368[k] = f_1 * dpg_266[k]
                   + f_4 * pc_x[k] * fpg_266[k];

        t_369[k] = f_1 * dpg_267[k]
                   + f_4 * pc_x[k] * fpg_267[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_x, pc_x, pc_y, dph0_372, dph0_373, \
                         dpg_269, dph1_372, dph1_373, fpg_264, \
                         fpg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_4 * pc_y[k] * fpg_264[k];

        t_371[k] = f_1 * dpg_269[k]
                   + f_4 * pc_x[k] * fpg_269[k];

        t_372[k] = pa_x[k] * dph0_372[k]
                   - f_9 * pc_x[k] * dph1_372[k];

        t_373[k] = pa_x[k] * dph0_373[k]
                   - f_9 * pc_x[k] * dph1_373[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, pa_x, pc_x, pc_y, dph0_374, dph0_375, \
                         dph0_377, dph1_374, dph1_375, dph1_377, \
                         fpg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_x[k] * dph0_374[k]
                   - f_9 * pc_x[k] * dph1_374[k];

        t_375[k] = pa_x[k] * dph0_375[k]
                   - f_9 * pc_x[k] * dph1_375[k];

        t_376[k] = f_4 * pc_y[k] * fpg_269[k];

        t_377[k] = pa_x[k] * dph0_377[k]
                   - f_9 * pc_x[k] * dph1_377[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pb_x, pc_x, pc_z, fsh0_126, fsh0_127, fsg_90, \
                         fsg_91, fsh1_126, fsh1_127, fpg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_x[k] * fsh0_126[k]
                   + f_11 * fsg_90[k]
                   - f_9 * pc_x[k] * fsh1_126[k];

        t_379[k] = pb_x[k] * fsh0_127[k]
                   + f_16 * fsg_91[k]
                   - f_9 * pc_x[k] * fsh1_127[k];

        t_380[k] = f_4 * pc_z[k] * fpg_270[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, pb_x, pc_x, pc_z, fsh0_129, fsh0_131, fsg_93, \
                         fsg_95, fsh1_129, fsh1_131, fpg_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = pb_x[k] * fsh0_129[k]
                   + f_0 * fsg_93[k]
                   - f_9 * pc_x[k] * fsh1_129[k];

        t_382[k] = f_4 * pc_z[k] * fpg_271[k];

        t_383[k] = pb_x[k] * fsh0_131[k]
                   + f_0 * fsg_95[k]
                   - f_9 * pc_x[k] * fsh1_131[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pb_x, pc_x, pc_z, fsh0_132, fsh0_134, fsg_96, \
                         fsg_98, fsh1_132, fsh1_134, fpg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pb_x[k] * fsh0_132[k]
                   + f_10 * fsg_96[k]
                   - f_9 * pc_x[k] * fsh1_132[k];

        t_385[k] = f_4 * pc_z[k] * fpg_273[k];

        t_386[k] = pb_x[k] * fsh0_134[k]
                   + f_10 * fsg_98[k]
                   - f_9 * pc_x[k] * fsh1_134[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pb_x, pc_x, fsh0_135, fsg_99, fsg_100, \
                         fsg_101, fsg_102, fsh1_135, fpg_280, fpg_281, \
                         fpg_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = pb_x[k] * fsh0_135[k]
                   + f_10 * fsg_99[k]
                   - f_9 * pc_x[k] * fsh1_135[k];

        t_388[k] = f_1 * fsg_100[k]
                   + f_4 * pc_x[k] * fpg_280[k];

        t_389[k] = f_1 * fsg_101[k]
                   + f_4 * pc_x[k] * fpg_281[k];

        t_390[k] = f_1 * fsg_102[k]
                   + f_4 * pc_x[k] * fpg_282[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_x, pc_x, pc_z, fsh0_141, fsg_103, \
                         fsg_104, fsh1_141, fpg_280, fpg_283, fpg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_1 * fsg_103[k]
                   + f_4 * pc_x[k] * fpg_283[k];

        t_392[k] = f_1 * fsg_104[k]
                   + f_4 * pc_x[k] * fpg_284[k];

        t_393[k] = pb_x[k] * fsh0_141[k]
                   - f_9 * pc_x[k] * fsh1_141[k];

        t_394[k] = f_4 * pc_z[k] * fpg_280[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pb_x, pc_x, pc_y, dpg_149, fsh0_143, \
                         fsh0_144, fsh0_146, fsh1_143, fsh1_144, fsh1_146, \
                         fpg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pb_x[k] * fsh0_143[k]
                   - f_9 * pc_x[k] * fsh1_143[k];

        t_396[k] = pb_x[k] * fsh0_144[k]
                   - f_9 * pc_x[k] * fsh1_144[k];

        t_397[k] = f_0 * dpg_149[k]
                   + f_4 * pc_y[k] * fpg_284[k];

        t_398[k] = pb_x[k] * fsh0_146[k]
                   - f_9 * pc_x[k] * fsh1_146[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, pc_x, pc_z, fpf0_190, fpf0_191, \
                         fpf0_193, fpf1_190, fpf1_191, fpf1_193, fpg_285, fpg_286, \
                         fpg_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_2 * fpf0_190[k]
                   - f_3 * fpf1_190[k]
                   + f_4 * pc_x[k] * fpg_285[k];

        t_400[k] = f_14 * fpf0_191[k]
                   - f_15 * fpf1_191[k]
                   + f_4 * pc_x[k] * fpg_286[k];

        t_401[k] = f_4 * pc_z[k] * fpg_285[k];

        t_402[k] = f_7 * fpf0_193[k]
                   - f_8 * fpf1_193[k]
                   + f_4 * pc_x[k] * fpg_288[k];

        t_403[k] = f_4 * pc_z[k] * fpg_286[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pc_x, pc_z, fpf0_195, fpf0_196, fpf0_198, \
                         fpf1_195, fpf1_196, fpf1_198, fpg_288, fpg_290, fpg_291, \
                         fpg_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_7 * fpf0_195[k]
                   - f_8 * fpf1_195[k]
                   + f_4 * pc_x[k] * fpg_290[k];

        t_405[k] = f_5 * fpf0_196[k]
                   - f_6 * fpf1_196[k]
                   + f_4 * pc_x[k] * fpg_291[k];

        t_406[k] = f_4 * pc_z[k] * fpg_288[k];

        t_407[k] = f_5 * fpf0_198[k]
                   - f_6 * fpf1_198[k]
                   + f_4 * pc_x[k] * fpg_293[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, t_413, pc_x, fpf0_199, fpf1_199, \
                         fpg_294, fpg_295, fpg_296, fpg_297, fpg_298, \
                         fpg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_5 * fpf0_199[k]
                   - f_6 * fpf1_199[k]
                   + f_4 * pc_x[k] * fpg_294[k];

        t_409[k] = f_4 * pc_x[k] * fpg_295[k];

        t_410[k] = f_4 * pc_x[k] * fpg_296[k];

        t_411[k] = f_4 * pc_x[k] * fpg_297[k];

        t_412[k] = f_4 * pc_x[k] * fpg_298[k];

        t_413[k] = f_4 * pc_x[k] * fpg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_y, pc_z, dpg_160, fsg_100, fpf0_196, \
                         fpf0_197, fpf1_196, fpf1_197, fpg_295, fpg_296, \
                         fpg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_0 * dpg_160[k]
                   + f_1 * fsg_100[k]
                   + f_2 * fpf0_196[k]
                   - f_3 * fpf1_196[k]
                   + f_4 * pc_y[k] * fpg_295[k];

        t_415[k] = f_4 * pc_z[k] * fpg_295[k];

        t_416[k] = f_5 * fpf0_196[k]
                   - f_6 * fpf1_196[k]
                   + f_4 * pc_z[k] * fpg_296[k];

        t_417[k] = f_7 * fpf0_197[k]
                   - f_8 * fpf1_197[k]
                   + f_4 * pc_z[k] * fpg_297[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pb_z, pc_y, pc_z, dpg_164, fsh0_126, \
                         fsh0_127, fsg_104, fsh1_126, fsh1_127, fpf0_199, fpf1_199, \
                         fpg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_0 * dpg_164[k]
                   + f_1 * fsg_104[k]
                   + f_4 * pc_y[k] * fpg_299[k];

        t_419[k] = f_2 * fpf0_199[k]
                   - f_3 * fpf1_199[k]
                   + f_4 * pc_z[k] * fpg_299[k];

        t_420[k] = pb_z[k] * fsh0_126[k]
                   - f_9 * pc_z[k] * fsh1_126[k];

        t_421[k] = pb_z[k] * fsh0_127[k]
                   - f_9 * pc_z[k] * fsh1_127[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, pb_z, pc_x, pc_z, fsh0_129, fsg_90, \
                         fsg_91, fsh1_129, fpf0_205, fpf1_205, fpg_300, fpg_301, \
                         fpg_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_1 * fsg_90[k]
                   + f_4 * pc_z[k] * fpg_300[k];

        t_423[k] = pb_z[k] * fsh0_129[k]
                   - f_9 * pc_z[k] * fsh1_129[k];

        t_424[k] = f_1 * fsg_91[k]
                   + f_4 * pc_z[k] * fpg_301[k];

        t_425[k] = f_7 * fpf0_205[k]
                   - f_8 * fpf1_205[k]
                   + f_4 * pc_x[k] * fpg_305[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pb_z, pc_x, pc_z, fsh0_132, fsg_93, fsh1_132, \
                         fpf0_208, fpf1_208, fpg_303, fpg_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = pb_z[k] * fsh0_132[k]
                   - f_9 * pc_z[k] * fsh1_132[k];

        t_427[k] = f_1 * fsg_93[k]
                   + f_4 * pc_z[k] * fpg_303[k];

        t_428[k] = f_5 * fpf0_208[k]
                   - f_6 * fpf1_208[k]
                   + f_4 * pc_x[k] * fpg_308[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, t_434, pc_x, fpf0_209, fpf1_209, \
                         fpg_309, fpg_310, fpg_311, fpg_312, fpg_313, \
                         fpg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_5 * fpf0_209[k]
                   - f_6 * fpf1_209[k]
                   + f_4 * pc_x[k] * fpg_309[k];

        t_430[k] = f_4 * pc_x[k] * fpg_310[k];

        t_431[k] = f_4 * pc_x[k] * fpg_311[k];

        t_432[k] = f_4 * pc_x[k] * fpg_312[k];

        t_433[k] = f_4 * pc_x[k] * fpg_313[k];

        t_434[k] = f_4 * pc_x[k] * fpg_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pb_z, pc_z, fsh0_141, fsh0_143, fsh0_144, \
                         fsg_100, fsg_101, fsg_102, fsh1_141, fsh1_143, fsh1_144, \
                         fpg_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = pb_z[k] * fsh0_141[k]
                   - f_9 * pc_z[k] * fsh1_141[k];

        t_436[k] = f_1 * fsg_100[k]
                   + f_4 * pc_z[k] * fpg_310[k];

        t_437[k] = pb_z[k] * fsh0_143[k]
                   + f_10 * fsg_101[k]
                   - f_9 * pc_z[k] * fsh1_143[k];

        t_438[k] = pb_z[k] * fsh0_144[k]
                   + f_0 * fsg_102[k]
                   - f_9 * pc_z[k] * fsh1_144[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, pa_z, pb_z, pc_y, pc_z, dph0_189, dpg_179, \
                         dph1_189, fsh0_146, fsg_104, fsh1_146, \
                         fpg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_0 * dpg_179[k]
                   + f_4 * pc_y[k] * fpg_314[k];

        t_440[k] = pb_z[k] * fsh0_146[k]
                   + f_11 * fsg_104[k]
                   - f_9 * pc_z[k] * fsh1_146[k];

        t_441[k] = pa_z[k] * dph0_189[k]
                   - f_9 * pc_z[k] * dph1_189[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, pa_z, pb_x, pc_x, pc_z, dph0_190, dph0_192, \
                         dph1_190, dph1_192, fsh0_149, fsg_107, \
                         fsh1_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = pa_z[k] * dph0_190[k]
                   - f_9 * pc_z[k] * dph1_190[k];

        t_443[k] = pb_x[k] * fsh0_149[k]
                   + f_16 * fsg_107[k]
                   - f_9 * pc_x[k] * fsh1_149[k];

        t_444[k] = pa_z[k] * dph0_192[k]
                   - f_9 * pc_z[k] * dph1_192[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, pa_z, pb_x, pc_x, pc_z, dph0_193, dph0_195, \
                         dpg_136, dph1_193, dph1_195, fsh0_152, fsg_110, \
                         fsh1_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pa_z[k] * dph0_193[k]
                   + f_1 * dpg_136[k]
                   - f_9 * pc_z[k] * dph1_193[k];

        t_446[k] = pb_x[k] * fsh0_152[k]
                   + f_0 * fsg_110[k]
                   - f_9 * pc_x[k] * fsh1_152[k];

        t_447[k] = pa_z[k] * dph0_195[k]
                   - f_9 * pc_z[k] * dph1_195[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, pa_z, pb_x, pc_x, pc_z, dph0_196, dph0_197, \
                         dpg_138, dpg_139, dph1_196, dph1_197, fsh0_156, fsg_114, \
                         fsh1_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pa_z[k] * dph0_196[k]
                   + f_1 * dpg_138[k]
                   - f_9 * pc_z[k] * dph1_196[k];

        t_449[k] = pa_z[k] * dph0_197[k]
                   + f_10 * dpg_139[k]
                   - f_9 * pc_z[k] * dph1_197[k];

        t_450[k] = pb_x[k] * fsh0_156[k]
                   + f_10 * fsg_114[k]
                   - f_9 * pc_x[k] * fsh1_156[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, pc_x, fsg_115, fsg_116, fsg_117, \
                         fsg_118, fsg_119, fpg_325, fpg_326, fpg_327, fpg_328, \
                         fpg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_1 * fsg_115[k]
                   + f_4 * pc_x[k] * fpg_325[k];

        t_452[k] = f_1 * fsg_116[k]
                   + f_4 * pc_x[k] * fpg_326[k];

        t_453[k] = f_1 * fsg_117[k]
                   + f_4 * pc_x[k] * fpg_327[k];

        t_454[k] = f_1 * fsg_118[k]
                   + f_4 * pc_x[k] * fpg_328[k];

        t_455[k] = f_1 * fsg_119[k]
                   + f_4 * pc_x[k] * fpg_329[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pa_z, pb_x, pc_x, pc_z, dph0_204, \
                         dpg_145, dph1_204, fsh0_164, fsh0_165, fsh1_164, fsh1_165, \
                         fpg_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = pa_z[k] * dph0_204[k]
                   - f_9 * pc_z[k] * dph1_204[k];

        t_457[k] = f_1 * dpg_145[k]
                   + f_4 * pc_z[k] * fpg_325[k];

        t_458[k] = pb_x[k] * fsh0_164[k]
                   - f_9 * pc_x[k] * fsh1_164[k];

        t_459[k] = pb_x[k] * fsh0_165[k]
                   - f_9 * pc_x[k] * fsh1_165[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pa_z, pb_x, pc_x, pc_y, pc_z, dph0_210, dpg_194, \
                         dph1_210, fsh0_167, fsh1_167, fpg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_10 * dpg_194[k]
                   + f_4 * pc_y[k] * fpg_329[k];

        t_461[k] = pb_x[k] * fsh0_167[k]
                   - f_9 * pc_x[k] * fsh1_167[k];

        t_462[k] = pa_z[k] * dph0_210[k]
                   - f_9 * pc_z[k] * dph1_210[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pa_z, pc_x, pc_z, dph0_211, dph0_213, dph1_211, \
                         dph1_213, fpf0_222, fpf1_222, fpg_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = pa_z[k] * dph0_211[k]
                   - f_9 * pc_z[k] * dph1_211[k];

        t_464[k] = f_14 * fpf0_222[k]
                   - f_15 * fpf1_222[k]
                   + f_4 * pc_x[k] * fpg_332[k];

        t_465[k] = pa_z[k] * dph0_213[k]
                   - f_9 * pc_z[k] * dph1_213[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, pa_z, pc_x, pc_z, dph0_216, dph1_216, fpf0_224, \
                         fpf0_225, fpf1_224, fpf1_225, fpg_334, \
                         fpg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_7 * fpf0_224[k]
                   - f_8 * fpf1_224[k]
                   + f_4 * pc_x[k] * fpg_334[k];

        t_467[k] = f_7 * fpf0_225[k]
                   - f_8 * fpf1_225[k]
                   + f_4 * pc_x[k] * fpg_335[k];

        t_468[k] = pa_z[k] * dph0_216[k]
                   - f_9 * pc_z[k] * dph1_216[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, pc_x, fpf0_227, fpf0_228, fpf0_229, \
                         fpf1_227, fpf1_228, fpf1_229, fpg_337, fpg_338, fpg_339, \
                         fpg_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_5 * fpf0_227[k]
                   - f_6 * fpf1_227[k]
                   + f_4 * pc_x[k] * fpg_337[k];

        t_470[k] = f_5 * fpf0_228[k]
                   - f_6 * fpf1_228[k]
                   + f_4 * pc_x[k] * fpg_338[k];

        t_471[k] = f_5 * fpf0_229[k]
                   - f_6 * fpf1_229[k]
                   + f_4 * pc_x[k] * fpg_339[k];

        t_472[k] = f_4 * pc_x[k] * fpg_340[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, t_477, pa_z, pc_x, pc_z, dph0_225, \
                         dph1_225, fpg_341, fpg_342, fpg_343, fpg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_4 * pc_x[k] * fpg_341[k];

        t_474[k] = f_4 * pc_x[k] * fpg_342[k];

        t_475[k] = f_4 * pc_x[k] * fpg_343[k];

        t_476[k] = f_4 * pc_x[k] * fpg_344[k];

        t_477[k] = pa_z[k] * dph0_225[k]
                   - f_9 * pc_z[k] * dph1_225[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pa_z, pc_z, dph0_227, dph0_228, dpg_160, \
                         dpg_161, dpg_162, dph1_227, dph1_228, \
                         fpg_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_1 * dpg_160[k]
                   + f_4 * pc_z[k] * fpg_340[k];

        t_479[k] = pa_z[k] * dph0_227[k]
                   + f_10 * dpg_161[k]
                   - f_9 * pc_z[k] * dph1_227[k];

        t_480[k] = pa_z[k] * dph0_228[k]
                   + f_0 * dpg_162[k]
                   - f_9 * pc_z[k] * dph1_228[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, pc_x, pc_y, pc_z, dpg_164, dpg_209, fsg_119, \
                         fpf0_229, fpf0_230, fpf1_229, fpf1_230, fpg_344, \
                         fpg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_10 * dpg_209[k]
                   + f_1 * fsg_119[k]
                   + f_4 * pc_y[k] * fpg_344[k];

        t_482[k] = f_1 * dpg_164[k]
                   + f_2 * fpf0_229[k]
                   - f_3 * fpf1_229[k]
                   + f_4 * pc_z[k] * fpg_344[k];

        t_483[k] = f_2 * fpf0_230[k]
                   - f_3 * fpf1_230[k]
                   + f_4 * pc_x[k] * fpg_345[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, pc_x, fpf0_231, fpf0_232, fpf0_233, fpf1_231, \
                         fpf1_232, fpf1_233, fpg_346, fpg_347, \
                         fpg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_14 * fpf0_231[k]
                   - f_15 * fpf1_231[k]
                   + f_4 * pc_x[k] * fpg_346[k];

        t_485[k] = f_14 * fpf0_232[k]
                   - f_15 * fpf1_232[k]
                   + f_4 * pc_x[k] * fpg_347[k];

        t_486[k] = f_7 * fpf0_233[k]
                   - f_8 * fpf1_233[k]
                   + f_4 * pc_x[k] * fpg_348[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, pc_x, fpf0_234, fpf0_235, fpf0_236, fpf1_234, \
                         fpf1_235, fpf1_236, fpg_349, fpg_350, \
                         fpg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_7 * fpf0_234[k]
                   - f_8 * fpf1_234[k]
                   + f_4 * pc_x[k] * fpg_349[k];

        t_488[k] = f_7 * fpf0_235[k]
                   - f_8 * fpf1_235[k]
                   + f_4 * pc_x[k] * fpg_350[k];

        t_489[k] = f_5 * fpf0_236[k]
                   - f_6 * fpf1_236[k]
                   + f_4 * pc_x[k] * fpg_351[k];
    }
}

static auto
compute_prim_fph_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t pph0, const size_t pph1,
                                                          const size_t dph0, const size_t dpg,
                                                          const size_t dph1, const size_t fsh0,
                                                          const size_t fsg, const size_t fsh1,
                                                          const size_t fpf0, const size_t fpf1,
                                                          const size_t fpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);
    const auto f_14 = 1.5 / gamma;
    const auto f_15 = 1.5 * p / (gamma * q);
    const auto f_16 = 2.0 / q;

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
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *pph0_188 = buffer.data(pph0 + 188);

    const auto *pph1_188 = buffer.data(pph1 + 188);

    const auto *dph0_314 = buffer.data(dph0 + 314);
    const auto *dph0_315 = buffer.data(dph0 + 315);
    const auto *dph0_316 = buffer.data(dph0 + 316);
    const auto *dph0_317 = buffer.data(dph0 + 317);
    const auto *dph0_318 = buffer.data(dph0 + 318);
    const auto *dph0_319 = buffer.data(dph0 + 319);
    const auto *dph0_320 = buffer.data(dph0 + 320);
    const auto *dph0_321 = buffer.data(dph0 + 321);
    const auto *dph0_322 = buffer.data(dph0 + 322);
    const auto *dph0_323 = buffer.data(dph0 + 323);
    const auto *dph0_324 = buffer.data(dph0 + 324);
    const auto *dph0_335 = buffer.data(dph0 + 335);
    const auto *dph0_357 = buffer.data(dph0 + 357);
    const auto *dph0_359 = buffer.data(dph0 + 359);
    const auto *dph0_362 = buffer.data(dph0 + 362);
    const auto *dph0_366 = buffer.data(dph0 + 366);
    const auto *dph0_372 = buffer.data(dph0 + 372);
    const auto *dph0_374 = buffer.data(dph0 + 374);
    const auto *dph0_375 = buffer.data(dph0 + 375);
    const auto *dph0_377 = buffer.data(dph0 + 377);

    const auto *dpg_175 = buffer.data(dpg + 175);
    const auto *dpg_190 = buffer.data(dpg + 190);
    const auto *dpg_205 = buffer.data(dpg + 205);
    const auto *dpg_209 = buffer.data(dpg + 209);
    const auto *dpg_220 = buffer.data(dpg + 220);
    const auto *dpg_222 = buffer.data(dpg + 222);
    const auto *dpg_223 = buffer.data(dpg + 223);
    const auto *dpg_224 = buffer.data(dpg + 224);
    const auto *dpg_225 = buffer.data(dpg + 225);
    const auto *dpg_226 = buffer.data(dpg + 226);
    const auto *dpg_227 = buffer.data(dpg + 227);
    const auto *dpg_228 = buffer.data(dpg + 228);
    const auto *dpg_229 = buffer.data(dpg + 229);
    const auto *dpg_230 = buffer.data(dpg + 230);
    const auto *dpg_239 = buffer.data(dpg + 239);
    const auto *dpg_250 = buffer.data(dpg + 250);
    const auto *dpg_252 = buffer.data(dpg + 252);
    const auto *dpg_253 = buffer.data(dpg + 253);
    const auto *dpg_254 = buffer.data(dpg + 254);
    const auto *dpg_265 = buffer.data(dpg + 265);
    const auto *dpg_267 = buffer.data(dpg + 267);
    const auto *dpg_268 = buffer.data(dpg + 268);
    const auto *dpg_269 = buffer.data(dpg + 269);

    const auto *dph1_314 = buffer.data(dph1 + 314);
    const auto *dph1_315 = buffer.data(dph1 + 315);
    const auto *dph1_316 = buffer.data(dph1 + 316);
    const auto *dph1_317 = buffer.data(dph1 + 317);
    const auto *dph1_318 = buffer.data(dph1 + 318);
    const auto *dph1_319 = buffer.data(dph1 + 319);
    const auto *dph1_320 = buffer.data(dph1 + 320);
    const auto *dph1_321 = buffer.data(dph1 + 321);
    const auto *dph1_322 = buffer.data(dph1 + 322);
    const auto *dph1_323 = buffer.data(dph1 + 323);
    const auto *dph1_324 = buffer.data(dph1 + 324);
    const auto *dph1_335 = buffer.data(dph1 + 335);
    const auto *dph1_357 = buffer.data(dph1 + 357);
    const auto *dph1_359 = buffer.data(dph1 + 359);
    const auto *dph1_362 = buffer.data(dph1 + 362);
    const auto *dph1_366 = buffer.data(dph1 + 366);
    const auto *dph1_372 = buffer.data(dph1 + 372);
    const auto *dph1_374 = buffer.data(dph1 + 374);
    const auto *dph1_375 = buffer.data(dph1 + 375);
    const auto *dph1_377 = buffer.data(dph1 + 377);

    const auto *fsh0_183 = buffer.data(fsh0 + 183);
    const auto *fsh0_185 = buffer.data(fsh0 + 185);
    const auto *fsh0_186 = buffer.data(fsh0 + 186);
    const auto *fsh0_189 = buffer.data(fsh0 + 189);
    const auto *fsh0_191 = buffer.data(fsh0 + 191);
    const auto *fsh0_192 = buffer.data(fsh0 + 192);
    const auto *fsh0_194 = buffer.data(fsh0 + 194);
    const auto *fsh0_195 = buffer.data(fsh0 + 195);
    const auto *fsh0_196 = buffer.data(fsh0 + 196);
    const auto *fsh0_198 = buffer.data(fsh0 + 198);
    const auto *fsh0_204 = buffer.data(fsh0 + 204);
    const auto *fsh0_205 = buffer.data(fsh0 + 205);
    const auto *fsh0_206 = buffer.data(fsh0 + 206);
    const auto *fsh0_207 = buffer.data(fsh0 + 207);
    const auto *fsh0_209 = buffer.data(fsh0 + 209);

    const auto *fsg_115 = buffer.data(fsg + 115);
    const auto *fsg_130 = buffer.data(fsg + 130);
    const auto *fsg_131 = buffer.data(fsg + 131);
    const auto *fsg_132 = buffer.data(fsg + 132);
    const auto *fsg_133 = buffer.data(fsg + 133);
    const auto *fsg_134 = buffer.data(fsg + 134);
    const auto *fsg_135 = buffer.data(fsg + 135);
    const auto *fsg_137 = buffer.data(fsg + 137);
    const auto *fsg_138 = buffer.data(fsg + 138);
    const auto *fsg_140 = buffer.data(fsg + 140);
    const auto *fsg_141 = buffer.data(fsg + 141);
    const auto *fsg_142 = buffer.data(fsg + 142);
    const auto *fsg_144 = buffer.data(fsg + 144);
    const auto *fsg_145 = buffer.data(fsg + 145);
    const auto *fsg_146 = buffer.data(fsg + 146);
    const auto *fsg_147 = buffer.data(fsg + 147);
    const auto *fsg_148 = buffer.data(fsg + 148);
    const auto *fsg_149 = buffer.data(fsg + 149);

    const auto *fsh1_183 = buffer.data(fsh1 + 183);
    const auto *fsh1_185 = buffer.data(fsh1 + 185);
    const auto *fsh1_186 = buffer.data(fsh1 + 186);
    const auto *fsh1_189 = buffer.data(fsh1 + 189);
    const auto *fsh1_191 = buffer.data(fsh1 + 191);
    const auto *fsh1_192 = buffer.data(fsh1 + 192);
    const auto *fsh1_194 = buffer.data(fsh1 + 194);
    const auto *fsh1_195 = buffer.data(fsh1 + 195);
    const auto *fsh1_196 = buffer.data(fsh1 + 196);
    const auto *fsh1_198 = buffer.data(fsh1 + 198);
    const auto *fsh1_204 = buffer.data(fsh1 + 204);
    const auto *fsh1_205 = buffer.data(fsh1 + 205);
    const auto *fsh1_206 = buffer.data(fsh1 + 206);
    const auto *fsh1_207 = buffer.data(fsh1 + 207);
    const auto *fsh1_209 = buffer.data(fsh1 + 209);

    const auto *fpf0_236 = buffer.data(fpf0 + 236);
    const auto *fpf0_237 = buffer.data(fpf0 + 237);
    const auto *fpf0_238 = buffer.data(fpf0 + 238);
    const auto *fpf0_239 = buffer.data(fpf0 + 239);
    const auto *fpf0_250 = buffer.data(fpf0 + 250);
    const auto *fpf0_251 = buffer.data(fpf0 + 251);
    const auto *fpf0_252 = buffer.data(fpf0 + 252);
    const auto *fpf0_253 = buffer.data(fpf0 + 253);
    const auto *fpf0_254 = buffer.data(fpf0 + 254);
    const auto *fpf0_255 = buffer.data(fpf0 + 255);
    const auto *fpf0_256 = buffer.data(fpf0 + 256);
    const auto *fpf0_257 = buffer.data(fpf0 + 257);
    const auto *fpf0_258 = buffer.data(fpf0 + 258);
    const auto *fpf0_259 = buffer.data(fpf0 + 259);
    const auto *fpf0_261 = buffer.data(fpf0 + 261);
    const auto *fpf0_263 = buffer.data(fpf0 + 263);
    const auto *fpf0_264 = buffer.data(fpf0 + 264);
    const auto *fpf0_266 = buffer.data(fpf0 + 266);
    const auto *fpf0_267 = buffer.data(fpf0 + 267);
    const auto *fpf0_268 = buffer.data(fpf0 + 268);
    const auto *fpf0_283 = buffer.data(fpf0 + 283);
    const auto *fpf0_286 = buffer.data(fpf0 + 286);
    const auto *fpf0_287 = buffer.data(fpf0 + 287);
    const auto *fpf0_290 = buffer.data(fpf0 + 290);

    const auto *fpf1_236 = buffer.data(fpf1 + 236);
    const auto *fpf1_237 = buffer.data(fpf1 + 237);
    const auto *fpf1_238 = buffer.data(fpf1 + 238);
    const auto *fpf1_239 = buffer.data(fpf1 + 239);
    const auto *fpf1_250 = buffer.data(fpf1 + 250);
    const auto *fpf1_251 = buffer.data(fpf1 + 251);
    const auto *fpf1_252 = buffer.data(fpf1 + 252);
    const auto *fpf1_253 = buffer.data(fpf1 + 253);
    const auto *fpf1_254 = buffer.data(fpf1 + 254);
    const auto *fpf1_255 = buffer.data(fpf1 + 255);
    const auto *fpf1_256 = buffer.data(fpf1 + 256);
    const auto *fpf1_257 = buffer.data(fpf1 + 257);
    const auto *fpf1_258 = buffer.data(fpf1 + 258);
    const auto *fpf1_259 = buffer.data(fpf1 + 259);
    const auto *fpf1_261 = buffer.data(fpf1 + 261);
    const auto *fpf1_263 = buffer.data(fpf1 + 263);
    const auto *fpf1_264 = buffer.data(fpf1 + 264);
    const auto *fpf1_266 = buffer.data(fpf1 + 266);
    const auto *fpf1_267 = buffer.data(fpf1 + 267);
    const auto *fpf1_268 = buffer.data(fpf1 + 268);
    const auto *fpf1_283 = buffer.data(fpf1 + 283);
    const auto *fpf1_286 = buffer.data(fpf1 + 286);
    const auto *fpf1_287 = buffer.data(fpf1 + 287);
    const auto *fpf1_290 = buffer.data(fpf1 + 290);

    const auto *fpg_352 = buffer.data(fpg + 352);
    const auto *fpg_353 = buffer.data(fpg + 353);
    const auto *fpg_354 = buffer.data(fpg + 354);
    const auto *fpg_355 = buffer.data(fpg + 355);
    const auto *fpg_356 = buffer.data(fpg + 356);
    const auto *fpg_357 = buffer.data(fpg + 357);
    const auto *fpg_358 = buffer.data(fpg + 358);
    const auto *fpg_359 = buffer.data(fpg + 359);
    const auto *fpg_370 = buffer.data(fpg + 370);
    const auto *fpg_371 = buffer.data(fpg + 371);
    const auto *fpg_372 = buffer.data(fpg + 372);
    const auto *fpg_373 = buffer.data(fpg + 373);
    const auto *fpg_374 = buffer.data(fpg + 374);
    const auto *fpg_375 = buffer.data(fpg + 375);
    const auto *fpg_376 = buffer.data(fpg + 376);
    const auto *fpg_377 = buffer.data(fpg + 377);
    const auto *fpg_378 = buffer.data(fpg + 378);
    const auto *fpg_379 = buffer.data(fpg + 379);
    const auto *fpg_380 = buffer.data(fpg + 380);
    const auto *fpg_381 = buffer.data(fpg + 381);
    const auto *fpg_382 = buffer.data(fpg + 382);
    const auto *fpg_383 = buffer.data(fpg + 383);
    const auto *fpg_384 = buffer.data(fpg + 384);
    const auto *fpg_385 = buffer.data(fpg + 385);
    const auto *fpg_386 = buffer.data(fpg + 386);
    const auto *fpg_387 = buffer.data(fpg + 387);
    const auto *fpg_388 = buffer.data(fpg + 388);
    const auto *fpg_389 = buffer.data(fpg + 389);
    const auto *fpg_391 = buffer.data(fpg + 391);
    const auto *fpg_393 = buffer.data(fpg + 393);
    const auto *fpg_394 = buffer.data(fpg + 394);
    const auto *fpg_396 = buffer.data(fpg + 396);
    const auto *fpg_397 = buffer.data(fpg + 397);
    const auto *fpg_398 = buffer.data(fpg + 398);
    const auto *fpg_400 = buffer.data(fpg + 400);
    const auto *fpg_401 = buffer.data(fpg + 401);
    const auto *fpg_402 = buffer.data(fpg + 402);
    const auto *fpg_403 = buffer.data(fpg + 403);
    const auto *fpg_404 = buffer.data(fpg + 404);
    const auto *fpg_405 = buffer.data(fpg + 405);
    const auto *fpg_407 = buffer.data(fpg + 407);
    const auto *fpg_410 = buffer.data(fpg + 410);
    const auto *fpg_415 = buffer.data(fpg + 415);
    const auto *fpg_416 = buffer.data(fpg + 416);
    const auto *fpg_417 = buffer.data(fpg + 417);
    const auto *fpg_418 = buffer.data(fpg + 418);
    const auto *fpg_419 = buffer.data(fpg + 419);
    const auto *fpg_420 = buffer.data(fpg + 420);
    const auto *fpg_422 = buffer.data(fpg + 422);
    const auto *fpg_423 = buffer.data(fpg + 423);
    const auto *fpg_425 = buffer.data(fpg + 425);
    const auto *fpg_426 = buffer.data(fpg + 426);
    const auto *fpg_427 = buffer.data(fpg + 427);
    const auto *fpg_430 = buffer.data(fpg + 430);
    const auto *fpg_431 = buffer.data(fpg + 431);
    const auto *fpg_432 = buffer.data(fpg + 432);
    const auto *fpg_433 = buffer.data(fpg + 433);
    const auto *fpg_434 = buffer.data(fpg + 434);
    const auto *fpg_435 = buffer.data(fpg + 435);

#pragma omp simd aligned(t_490, t_491, t_492, t_493, pc_x, fpf0_237, fpf0_238, fpf0_239, \
                         fpf1_237, fpf1_238, fpf1_239, fpg_352, fpg_353, fpg_354, \
                         fpg_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = f_5 * fpf0_237[k]
                   - f_6 * fpf1_237[k]
                   + f_4 * pc_x[k] * fpg_352[k];

        t_491[k] = f_5 * fpf0_238[k]
                   - f_6 * fpf1_238[k]
                   + f_4 * pc_x[k] * fpg_353[k];

        t_492[k] = f_5 * fpf0_239[k]
                   - f_6 * fpf1_239[k]
                   + f_4 * pc_x[k] * fpg_354[k];

        t_493[k] = f_4 * pc_x[k] * fpg_355[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, t_498, pc_x, pc_y, dpg_220, fpf0_236, \
                         fpf1_236, fpg_355, fpg_356, fpg_357, fpg_358, \
                         fpg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_4 * pc_x[k] * fpg_356[k];

        t_495[k] = f_4 * pc_x[k] * fpg_357[k];

        t_496[k] = f_4 * pc_x[k] * fpg_358[k];

        t_497[k] = f_4 * pc_x[k] * fpg_359[k];

        t_498[k] = f_10 * dpg_220[k]
                   + f_2 * fpf0_236[k]
                   - f_3 * fpf1_236[k]
                   + f_4 * pc_y[k] * fpg_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, dpg_175, dpg_222, dpg_223, fsg_115, \
                         fpf0_238, fpf0_239, fpf1_238, fpf1_239, fpg_355, fpg_357, \
                         fpg_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_1 * dpg_175[k]
                   + f_1 * fsg_115[k]
                   + f_4 * pc_z[k] * fpg_355[k];

        t_500[k] = f_10 * dpg_222[k]
                   + f_7 * fpf0_238[k]
                   - f_8 * fpf1_238[k]
                   + f_4 * pc_y[k] * fpg_357[k];

        t_501[k] = f_10 * dpg_223[k]
                   + f_5 * fpf0_239[k]
                   - f_6 * fpf1_239[k]
                   + f_4 * pc_y[k] * fpg_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pa_y, pc_y, pph0_188, pph1_188, dph0_314, \
                         dph0_315, dpg_224, dph1_314, dph1_315, \
                         fpg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_10 * dpg_224[k]
                   + f_4 * pc_y[k] * fpg_359[k];

        t_503[k] = f_12 * pph0_188[k]
                   - f_13 * pph1_188[k]
                   + pa_y[k] * dph0_314[k]
                   - f_9 * pc_y[k] * dph1_314[k];

        t_504[k] = pa_y[k] * dph0_315[k]
                   - f_9 * pc_y[k] * dph1_315[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pa_y, pc_y, dph0_316, dph0_317, dph0_318, \
                         dpg_225, dpg_226, dph1_316, dph1_317, \
                         dph1_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = pa_y[k] * dph0_316[k]
                   + f_1 * dpg_225[k]
                   - f_9 * pc_y[k] * dph1_316[k];

        t_506[k] = pa_y[k] * dph0_317[k]
                   - f_9 * pc_y[k] * dph1_317[k];

        t_507[k] = pa_y[k] * dph0_318[k]
                   + f_10 * dpg_226[k]
                   - f_9 * pc_y[k] * dph1_318[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pa_y, pc_y, dph0_319, dph0_320, dph0_321, \
                         dpg_227, dpg_228, dph1_319, dph1_320, \
                         dph1_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = pa_y[k] * dph0_319[k]
                   + f_1 * dpg_227[k]
                   - f_9 * pc_y[k] * dph1_319[k];

        t_509[k] = pa_y[k] * dph0_320[k]
                   - f_9 * pc_y[k] * dph1_320[k];

        t_510[k] = pa_y[k] * dph0_321[k]
                   + f_0 * dpg_228[k]
                   - f_9 * pc_y[k] * dph1_321[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pa_y, pc_y, dph0_322, dph0_323, dph0_324, \
                         dpg_229, dpg_230, dph1_322, dph1_323, \
                         dph1_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = pa_y[k] * dph0_322[k]
                   + f_10 * dpg_229[k]
                   - f_9 * pc_y[k] * dph1_322[k];

        t_512[k] = pa_y[k] * dph0_323[k]
                   + f_1 * dpg_230[k]
                   - f_9 * pc_y[k] * dph1_323[k];

        t_513[k] = pa_y[k] * dph0_324[k]
                   - f_9 * pc_y[k] * dph1_324[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, pc_x, fsg_130, fsg_131, fsg_132, \
                         fsg_133, fsg_134, fpg_370, fpg_371, fpg_372, fpg_373, \
                         fpg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_1 * fsg_130[k]
                   + f_4 * pc_x[k] * fpg_370[k];

        t_515[k] = f_1 * fsg_131[k]
                   + f_4 * pc_x[k] * fpg_371[k];

        t_516[k] = f_1 * fsg_132[k]
                   + f_4 * pc_x[k] * fpg_372[k];

        t_517[k] = f_1 * fsg_133[k]
                   + f_4 * pc_x[k] * fpg_373[k];

        t_518[k] = f_1 * fsg_134[k]
                   + f_4 * pc_x[k] * fpg_374[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, pb_x, pc_x, pc_z, dpg_190, fsh0_183, \
                         fsh0_185, fsh0_186, fsh1_183, fsh1_185, fsh1_186, \
                         fpg_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = pb_x[k] * fsh0_183[k]
                   - f_9 * pc_x[k] * fsh1_183[k];

        t_520[k] = f_10 * dpg_190[k]
                   + f_4 * pc_z[k] * fpg_370[k];

        t_521[k] = pb_x[k] * fsh0_185[k]
                   - f_9 * pc_x[k] * fsh1_185[k];

        t_522[k] = pb_x[k] * fsh0_186[k]
                   - f_9 * pc_x[k] * fsh1_186[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pa_y, pc_x, pc_y, dph0_335, dpg_239, dph1_335, \
                         fpf0_250, fpf1_250, fpg_374, fpg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_1 * dpg_239[k]
                   + f_4 * pc_y[k] * fpg_374[k];

        t_524[k] = pa_y[k] * dph0_335[k]
                   - f_9 * pc_y[k] * dph1_335[k];

        t_525[k] = f_2 * fpf0_250[k]
                   - f_3 * fpf1_250[k]
                   + f_4 * pc_x[k] * fpg_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_x, fpf0_251, fpf0_252, fpf0_253, fpf1_251, \
                         fpf1_252, fpf1_253, fpg_376, fpg_377, \
                         fpg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_14 * fpf0_251[k]
                   - f_15 * fpf1_251[k]
                   + f_4 * pc_x[k] * fpg_376[k];

        t_527[k] = f_14 * fpf0_252[k]
                   - f_15 * fpf1_252[k]
                   + f_4 * pc_x[k] * fpg_377[k];

        t_528[k] = f_7 * fpf0_253[k]
                   - f_8 * fpf1_253[k]
                   + f_4 * pc_x[k] * fpg_378[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, pc_x, fpf0_254, fpf0_255, fpf0_256, fpf1_254, \
                         fpf1_255, fpf1_256, fpg_379, fpg_380, \
                         fpg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_7 * fpf0_254[k]
                   - f_8 * fpf1_254[k]
                   + f_4 * pc_x[k] * fpg_379[k];

        t_530[k] = f_7 * fpf0_255[k]
                   - f_8 * fpf1_255[k]
                   + f_4 * pc_x[k] * fpg_380[k];

        t_531[k] = f_5 * fpf0_256[k]
                   - f_6 * fpf1_256[k]
                   + f_4 * pc_x[k] * fpg_381[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pc_x, fpf0_257, fpf0_258, fpf0_259, \
                         fpf1_257, fpf1_258, fpf1_259, fpg_382, fpg_383, fpg_384, \
                         fpg_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_5 * fpf0_257[k]
                   - f_6 * fpf1_257[k]
                   + f_4 * pc_x[k] * fpg_382[k];

        t_533[k] = f_5 * fpf0_258[k]
                   - f_6 * fpf1_258[k]
                   + f_4 * pc_x[k] * fpg_383[k];

        t_534[k] = f_5 * fpf0_259[k]
                   - f_6 * fpf1_259[k]
                   + f_4 * pc_x[k] * fpg_384[k];

        t_535[k] = f_4 * pc_x[k] * fpg_385[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, pc_x, pc_y, dpg_250, fsg_130, \
                         fpf0_256, fpf1_256, fpg_385, fpg_386, fpg_387, fpg_388, \
                         fpg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_4 * pc_x[k] * fpg_386[k];

        t_537[k] = f_4 * pc_x[k] * fpg_387[k];

        t_538[k] = f_4 * pc_x[k] * fpg_388[k];

        t_539[k] = f_4 * pc_x[k] * fpg_389[k];

        t_540[k] = f_1 * dpg_250[k]
                   + f_1 * fsg_130[k]
                   + f_2 * fpf0_256[k]
                   - f_3 * fpf1_256[k]
                   + f_4 * pc_y[k] * fpg_385[k];
    }

#pragma omp simd aligned(t_541, t_542, pc_y, pc_z, dpg_205, dpg_252, fsg_132, fpf0_258, \
                         fpf1_258, fpg_385, fpg_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_10 * dpg_205[k]
                   + f_4 * pc_z[k] * fpg_385[k];

        t_542[k] = f_1 * dpg_252[k]
                   + f_1 * fsg_132[k]
                   + f_7 * fpf0_258[k]
                   - f_8 * fpf1_258[k]
                   + f_4 * pc_y[k] * fpg_387[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pc_y, pc_z, dpg_209, dpg_253, dpg_254, fsg_133, \
                         fsg_134, fpf0_259, fpf1_259, fpg_388, \
                         fpg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_1 * dpg_253[k]
                   + f_1 * fsg_133[k]
                   + f_5 * fpf0_259[k]
                   - f_6 * fpf1_259[k]
                   + f_4 * pc_y[k] * fpg_388[k];

        t_544[k] = f_1 * dpg_254[k]
                   + f_1 * fsg_134[k]
                   + f_4 * pc_y[k] * fpg_389[k];

        t_545[k] = f_10 * dpg_209[k]
                   + f_2 * fpf0_259[k]
                   - f_3 * fpf1_259[k]
                   + f_4 * pc_z[k] * fpg_389[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, pa_y, pc_x, pc_y, dph0_357, dph0_359, dph1_357, \
                         dph1_359, fpf0_261, fpf1_261, fpg_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pa_y[k] * dph0_357[k]
                   - f_9 * pc_y[k] * dph1_357[k];

        t_547[k] = f_14 * fpf0_261[k]
                   - f_15 * fpf1_261[k]
                   + f_4 * pc_x[k] * fpg_391[k];

        t_548[k] = pa_y[k] * dph0_359[k]
                   - f_9 * pc_y[k] * dph1_359[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pa_y, pc_x, pc_y, dph0_362, dph1_362, fpf0_263, \
                         fpf0_264, fpf1_263, fpf1_264, fpg_393, \
                         fpg_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_7 * fpf0_263[k]
                   - f_8 * fpf1_263[k]
                   + f_4 * pc_x[k] * fpg_393[k];

        t_550[k] = f_7 * fpf0_264[k]
                   - f_8 * fpf1_264[k]
                   + f_4 * pc_x[k] * fpg_394[k];

        t_551[k] = pa_y[k] * dph0_362[k]
                   - f_9 * pc_y[k] * dph1_362[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pc_x, fpf0_266, fpf0_267, fpf0_268, fpf1_266, \
                         fpf1_267, fpf1_268, fpg_396, fpg_397, \
                         fpg_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_5 * fpf0_266[k]
                   - f_6 * fpf1_266[k]
                   + f_4 * pc_x[k] * fpg_396[k];

        t_553[k] = f_5 * fpf0_267[k]
                   - f_6 * fpf1_267[k]
                   + f_4 * pc_x[k] * fpg_397[k];

        t_554[k] = f_5 * fpf0_268[k]
                   - f_6 * fpf1_268[k]
                   + f_4 * pc_x[k] * fpg_398[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, t_560, pa_y, pc_x, pc_y, dph0_366, \
                         dph1_366, fpg_400, fpg_401, fpg_402, fpg_403, \
                         fpg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = pa_y[k] * dph0_366[k]
                   - f_9 * pc_y[k] * dph1_366[k];

        t_556[k] = f_4 * pc_x[k] * fpg_400[k];

        t_557[k] = f_4 * pc_x[k] * fpg_401[k];

        t_558[k] = f_4 * pc_x[k] * fpg_402[k];

        t_559[k] = f_4 * pc_x[k] * fpg_403[k];

        t_560[k] = f_4 * pc_x[k] * fpg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pa_y, pc_y, pc_z, dph0_372, dph0_374, dpg_220, \
                         dpg_265, dpg_267, dph1_372, dph1_374, fsg_130, \
                         fpg_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = pa_y[k] * dph0_372[k]
                   + f_11 * dpg_265[k]
                   - f_9 * pc_y[k] * dph1_372[k];

        t_562[k] = f_10 * dpg_220[k]
                   + f_1 * fsg_130[k]
                   + f_4 * pc_z[k] * fpg_400[k];

        t_563[k] = pa_y[k] * dph0_374[k]
                   + f_0 * dpg_267[k]
                   - f_9 * pc_y[k] * dph1_374[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_y, pc_y, dph0_375, dph0_377, dpg_268, \
                         dpg_269, dph1_375, dph1_377, fpg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = pa_y[k] * dph0_375[k]
                   + f_10 * dpg_268[k]
                   - f_9 * pc_y[k] * dph1_375[k];

        t_565[k] = f_1 * dpg_269[k]
                   + f_4 * pc_y[k] * fpg_404[k];

        t_566[k] = pa_y[k] * dph0_377[k]
                   - f_9 * pc_y[k] * dph1_377[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pb_x, pc_x, pc_y, fsh0_189, fsh0_191, fsg_135, \
                         fsg_137, fsh1_189, fsh1_191, fpg_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = pb_x[k] * fsh0_189[k]
                   + f_11 * fsg_135[k]
                   - f_9 * pc_x[k] * fsh1_189[k];

        t_568[k] = f_4 * pc_y[k] * fpg_405[k];

        t_569[k] = pb_x[k] * fsh0_191[k]
                   + f_16 * fsg_137[k]
                   - f_9 * pc_x[k] * fsh1_191[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pb_x, pc_x, pc_y, fsh0_192, fsh0_194, fsg_138, \
                         fsg_140, fsh1_192, fsh1_194, fpg_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = pb_x[k] * fsh0_192[k]
                   + f_0 * fsg_138[k]
                   - f_9 * pc_x[k] * fsh1_192[k];

        t_571[k] = f_4 * pc_y[k] * fpg_407[k];

        t_572[k] = pb_x[k] * fsh0_194[k]
                   + f_0 * fsg_140[k]
                   - f_9 * pc_x[k] * fsh1_194[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, pb_x, pc_x, pc_y, fsh0_195, fsh0_196, fsg_141, \
                         fsg_142, fsh1_195, fsh1_196, fpg_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = pb_x[k] * fsh0_195[k]
                   + f_10 * fsg_141[k]
                   - f_9 * pc_x[k] * fsh1_195[k];

        t_574[k] = pb_x[k] * fsh0_196[k]
                   + f_10 * fsg_142[k]
                   - f_9 * pc_x[k] * fsh1_196[k];

        t_575[k] = f_4 * pc_y[k] * fpg_410[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pb_x, pc_x, fsh0_198, fsg_144, fsg_145, \
                         fsg_146, fsg_147, fsh1_198, fpg_415, fpg_416, \
                         fpg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = pb_x[k] * fsh0_198[k]
                   + f_10 * fsg_144[k]
                   - f_9 * pc_x[k] * fsh1_198[k];

        t_577[k] = f_1 * fsg_145[k]
                   + f_4 * pc_x[k] * fpg_415[k];

        t_578[k] = f_1 * fsg_146[k]
                   + f_4 * pc_x[k] * fpg_416[k];

        t_579[k] = f_1 * fsg_147[k]
                   + f_4 * pc_x[k] * fpg_417[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pb_x, pc_x, fsh0_204, fsh0_205, fsg_148, \
                         fsg_149, fsh1_204, fsh1_205, fpg_418, \
                         fpg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_1 * fsg_148[k]
                   + f_4 * pc_x[k] * fpg_418[k];

        t_581[k] = f_1 * fsg_149[k]
                   + f_4 * pc_x[k] * fpg_419[k];

        t_582[k] = pb_x[k] * fsh0_204[k]
                   - f_9 * pc_x[k] * fsh1_204[k];

        t_583[k] = pb_x[k] * fsh0_205[k]
                   - f_9 * pc_x[k] * fsh1_205[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pb_x, pc_x, pc_y, fsh0_206, fsh0_207, \
                         fsh0_209, fsh1_206, fsh1_207, fsh1_209, \
                         fpg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = pb_x[k] * fsh0_206[k]
                   - f_9 * pc_x[k] * fsh1_206[k];

        t_585[k] = pb_x[k] * fsh0_207[k]
                   - f_9 * pc_x[k] * fsh1_207[k];

        t_586[k] = f_4 * pc_y[k] * fpg_419[k];

        t_587[k] = pb_x[k] * fsh0_209[k]
                   - f_9 * pc_x[k] * fsh1_209[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pb_y, pc_x, pc_y, fsh0_189, fsh0_191, \
                         fsg_135, fsh1_189, fsh1_191, fpf0_283, fpf1_283, fpg_420, \
                         fpg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pb_y[k] * fsh0_189[k]
                   - f_9 * pc_y[k] * fsh1_189[k];

        t_589[k] = f_1 * fsg_135[k]
                   + f_4 * pc_y[k] * fpg_420[k];

        t_590[k] = pb_y[k] * fsh0_191[k]
                   - f_9 * pc_y[k] * fsh1_191[k];

        t_591[k] = f_7 * fpf0_283[k]
                   - f_8 * fpf1_283[k]
                   + f_4 * pc_x[k] * fpg_423[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, pb_y, pc_x, pc_y, fsh0_194, fsg_137, fsh1_194, \
                         fpf0_286, fpf1_286, fpg_422, fpg_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_1 * fsg_137[k]
                   + f_4 * pc_y[k] * fpg_422[k];

        t_593[k] = pb_y[k] * fsh0_194[k]
                   - f_9 * pc_y[k] * fsh1_194[k];

        t_594[k] = f_5 * fpf0_286[k]
                   - f_6 * fpf1_286[k]
                   + f_4 * pc_x[k] * fpg_426[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, pb_y, pc_x, pc_y, fsh0_198, fsg_140, \
                         fsh1_198, fpf0_287, fpf1_287, fpg_425, fpg_427, \
                         fpg_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_5 * fpf0_287[k]
                   - f_6 * fpf1_287[k]
                   + f_4 * pc_x[k] * fpg_427[k];

        t_596[k] = f_1 * fsg_140[k]
                   + f_4 * pc_y[k] * fpg_425[k];

        t_597[k] = pb_y[k] * fsh0_198[k]
                   - f_9 * pc_y[k] * fsh1_198[k];

        t_598[k] = f_4 * pc_x[k] * fpg_430[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, t_602, t_603, pb_y, pc_x, pc_y, fsh0_204, \
                         fsg_145, fsh1_204, fpg_431, fpg_432, fpg_433, \
                         fpg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = f_4 * pc_x[k] * fpg_431[k];

        t_600[k] = f_4 * pc_x[k] * fpg_432[k];

        t_601[k] = f_4 * pc_x[k] * fpg_433[k];

        t_602[k] = f_4 * pc_x[k] * fpg_434[k];

        t_603[k] = pb_y[k] * fsh0_204[k]
                   + f_11 * fsg_145[k]
                   - f_9 * pc_y[k] * fsh1_204[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pb_y, pc_y, fsh0_205, fsh0_206, fsh0_207, \
                         fsg_146, fsg_147, fsg_148, fsh1_205, fsh1_206, \
                         fsh1_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pb_y[k] * fsh0_205[k]
                   + f_16 * fsg_146[k]
                   - f_9 * pc_y[k] * fsh1_205[k];

        t_605[k] = pb_y[k] * fsh0_206[k]
                   + f_0 * fsg_147[k]
                   - f_9 * pc_y[k] * fsh1_206[k];

        t_606[k] = pb_y[k] * fsh0_207[k]
                   + f_10 * fsg_148[k]
                   - f_9 * pc_y[k] * fsh1_207[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, pb_y, pc_x, pc_y, fsh0_209, fsg_149, \
                         fsh1_209, fpf0_290, fpf1_290, fpg_434, \
                         fpg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_1 * fsg_149[k]
                   + f_4 * pc_y[k] * fpg_434[k];

        t_608[k] = pb_y[k] * fsh0_209[k]
                   - f_9 * pc_y[k] * fsh1_209[k];

        t_609[k] = f_2 * fpf0_290[k]
                   - f_3 * fpf1_290[k]
                   + f_4 * pc_x[k] * fpg_435[k];

        t_610[k] = f_4 * pc_y[k] * fpg_435[k];
    }
}

static auto
compute_prim_fph_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t dpg, const size_t fsg,
                                                          const size_t fpf0, const size_t fpf1,
                                                          const size_t fpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_14 = 1.5 / gamma;
    const auto f_15 = 1.5 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpg_269 = buffer.data(dpg + 269);

    const auto *fsg_149 = buffer.data(fsg + 149);

    const auto *fpf0_292 = buffer.data(fpf0 + 292);
    const auto *fpf0_293 = buffer.data(fpf0 + 293);
    const auto *fpf0_295 = buffer.data(fpf0 + 295);
    const auto *fpf0_296 = buffer.data(fpf0 + 296);
    const auto *fpf0_297 = buffer.data(fpf0 + 297);
    const auto *fpf0_298 = buffer.data(fpf0 + 298);
    const auto *fpf0_299 = buffer.data(fpf0 + 299);

    const auto *fpf1_292 = buffer.data(fpf1 + 292);
    const auto *fpf1_293 = buffer.data(fpf1 + 293);
    const auto *fpf1_295 = buffer.data(fpf1 + 295);
    const auto *fpf1_296 = buffer.data(fpf1 + 296);
    const auto *fpf1_297 = buffer.data(fpf1 + 297);
    const auto *fpf1_298 = buffer.data(fpf1 + 298);
    const auto *fpf1_299 = buffer.data(fpf1 + 299);

    const auto *fpg_437 = buffer.data(fpg + 437);
    const auto *fpg_438 = buffer.data(fpg + 438);
    const auto *fpg_440 = buffer.data(fpg + 440);
    const auto *fpg_441 = buffer.data(fpg + 441);
    const auto *fpg_442 = buffer.data(fpg + 442);
    const auto *fpg_444 = buffer.data(fpg + 444);
    const auto *fpg_445 = buffer.data(fpg + 445);
    const auto *fpg_446 = buffer.data(fpg + 446);
    const auto *fpg_447 = buffer.data(fpg + 447);
    const auto *fpg_448 = buffer.data(fpg + 448);
    const auto *fpg_449 = buffer.data(fpg + 449);

#pragma omp simd aligned(t_611, t_612, t_613, t_614, pc_x, pc_y, fpf0_292, fpf0_293, fpf0_295, \
                         fpf1_292, fpf1_293, fpf1_295, fpg_437, fpg_438, \
                         fpg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_14 * fpf0_292[k]
                   - f_15 * fpf1_292[k]
                   + f_4 * pc_x[k] * fpg_437[k];

        t_612[k] = f_7 * fpf0_293[k]
                   - f_8 * fpf1_293[k]
                   + f_4 * pc_x[k] * fpg_438[k];

        t_613[k] = f_4 * pc_y[k] * fpg_437[k];

        t_614[k] = f_7 * fpf0_295[k]
                   - f_8 * fpf1_295[k]
                   + f_4 * pc_x[k] * fpg_440[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, pc_x, pc_y, fpf0_296, fpf0_297, fpf0_299, \
                         fpf1_296, fpf1_297, fpf1_299, fpg_440, fpg_441, fpg_442, \
                         fpg_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_5 * fpf0_296[k]
                   - f_6 * fpf1_296[k]
                   + f_4 * pc_x[k] * fpg_441[k];

        t_616[k] = f_5 * fpf0_297[k]
                   - f_6 * fpf1_297[k]
                   + f_4 * pc_x[k] * fpg_442[k];

        t_617[k] = f_4 * pc_y[k] * fpg_440[k];

        t_618[k] = f_5 * fpf0_299[k]
                   - f_6 * fpf1_299[k]
                   + f_4 * pc_x[k] * fpg_444[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, t_624, pc_x, pc_y, fpf0_296, \
                         fpf1_296, fpg_445, fpg_446, fpg_447, fpg_448, \
                         fpg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_4 * pc_x[k] * fpg_445[k];

        t_620[k] = f_4 * pc_x[k] * fpg_446[k];

        t_621[k] = f_4 * pc_x[k] * fpg_447[k];

        t_622[k] = f_4 * pc_x[k] * fpg_448[k];

        t_623[k] = f_4 * pc_x[k] * fpg_449[k];

        t_624[k] = f_2 * fpf0_296[k]
                   - f_3 * fpf1_296[k]
                   + f_4 * pc_y[k] * fpg_445[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, t_628, pc_y, fpf0_297, fpf0_298, fpf0_299, \
                         fpf1_297, fpf1_298, fpf1_299, fpg_446, fpg_447, fpg_448, \
                         fpg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_14 * fpf0_297[k]
                   - f_15 * fpf1_297[k]
                   + f_4 * pc_y[k] * fpg_446[k];

        t_626[k] = f_7 * fpf0_298[k]
                   - f_8 * fpf1_298[k]
                   + f_4 * pc_y[k] * fpg_447[k];

        t_627[k] = f_5 * fpf0_299[k]
                   - f_6 * fpf1_299[k]
                   + f_4 * pc_y[k] * fpg_448[k];

        t_628[k] = f_4 * pc_y[k] * fpg_449[k];
    }

#pragma omp simd aligned(t_629, pc_z, dpg_269, fsg_149, fpf0_299, fpf1_299, \
                         fpg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_0 * dpg_269[k]
                   + f_1 * fsg_149[k]
                   + f_2 * fpf0_299[k]
                   - f_3 * fpf1_299[k]
                   + f_4 * pc_z[k] * fpg_449[k];
    }
}

auto
compute_prim_fph_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t pph0,
                                                   const size_t pph1, const size_t dph0,
                                                   const size_t dpg, const size_t dph1,
                                                   const size_t fsh0, const size_t fsg,
                                                   const size_t fsh1, const size_t fpf0,
                                                   const size_t fpf1, const size_t fpg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fph_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, pph0,
                                                              pph1, dph0, dpg, dph1, fsh0, fsg,
                                                              fsh1, fpf0, fpf1, fpg, ncols,
                                                              gamma, p, q);

    compute_prim_fph_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, pph0,
                                                              pph1, dph0, dpg, dph1, fsh0, fsg,
                                                              fsh1, fpf0, fpf1, fpg, ncols,
                                                              gamma, p, q);

    compute_prim_fph_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, dph0,
                                                              dpg, dph1, fsh0, fsg, fsh1, fpf0,
                                                              fpf1, fpg, ncols, gamma, p, q);

    compute_prim_fph_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, dph0,
                                                              dpg, dph1, fsh0, fsg, fsh1, fpf0,
                                                              fpf1, fpg, ncols, gamma, p, q);

    compute_prim_fph_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, pph0,
                                                              pph1, dph0, dpg, dph1, fsh0, fsg,
                                                              fsh1, fpf0, fpf1, fpg, ncols,
                                                              gamma, p, q);

    compute_prim_fph_three_center_electron_repulsion_0_piece5(buffer, target, pc, dpg, fsg,
                                                              fpf0, fpf1, fpg, ncols, gamma, p,
                                                              q);
}

}  // namespace simdt3ceri
