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


#include "SimdThreeCenterElectronRepulsionVrrRecFPG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fpg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppg0, const size_t ppg1,
                                                          const size_t dpg0, const size_t dpf,
                                                          const size_t dpg1, const size_t fsg0,
                                                          const size_t fsf, const size_t fsg1,
                                                          const size_t fpd0, const size_t fpd1,
                                                          const size_t fpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppg0_70 = buffer.data(ppg0 + 70);

    const auto *ppg1_70 = buffer.data(ppg1 + 70);

    const auto *dpg0_0 = buffer.data(dpg0 + 0);
    const auto *dpg0_3 = buffer.data(dpg0 + 3);
    const auto *dpg0_5 = buffer.data(dpg0 + 5);
    const auto *dpg0_6 = buffer.data(dpg0 + 6);
    const auto *dpg0_9 = buffer.data(dpg0 + 9);
    const auto *dpg0_10 = buffer.data(dpg0 + 10);
    const auto *dpg0_14 = buffer.data(dpg0 + 14);
    const auto *dpg0_15 = buffer.data(dpg0 + 15);
    const auto *dpg0_18 = buffer.data(dpg0 + 18);
    const auto *dpg0_25 = buffer.data(dpg0 + 25);
    const auto *dpg0_30 = buffer.data(dpg0 + 30);
    const auto *dpg0_35 = buffer.data(dpg0 + 35);
    const auto *dpg0_44 = buffer.data(dpg0 + 44);
    const auto *dpg0_70 = buffer.data(dpg0 + 70);

    const auto *dpf_0 = buffer.data(dpf + 0);
    const auto *dpf_1 = buffer.data(dpf + 1);
    const auto *dpf_2 = buffer.data(dpf + 2);
    const auto *dpf_6 = buffer.data(dpf + 6);
    const auto *dpf_9 = buffer.data(dpf + 9);
    const auto *dpf_10 = buffer.data(dpf + 10);
    const auto *dpf_16 = buffer.data(dpf + 16);
    const auto *dpf_19 = buffer.data(dpf + 19);
    const auto *dpf_20 = buffer.data(dpf + 20);
    const auto *dpf_26 = buffer.data(dpf + 26);
    const auto *dpf_29 = buffer.data(dpf + 29);
    const auto *dpf_36 = buffer.data(dpf + 36);
    const auto *dpf_38 = buffer.data(dpf + 38);
    const auto *dpf_40 = buffer.data(dpf + 40);
    const auto *dpf_43 = buffer.data(dpf + 43);
    const auto *dpf_46 = buffer.data(dpf + 46);
    const auto *dpf_48 = buffer.data(dpf + 48);
    const auto *dpf_49 = buffer.data(dpf + 49);
    const auto *dpf_56 = buffer.data(dpf + 56);
    const auto *dpf_58 = buffer.data(dpf + 58);
    const auto *dpf_59 = buffer.data(dpf + 59);
    const auto *dpf_67 = buffer.data(dpf + 67);
    const auto *dpf_69 = buffer.data(dpf + 69);
    const auto *dpf_76 = buffer.data(dpf + 76);
    const auto *dpf_77 = buffer.data(dpf + 77);
    const auto *dpf_79 = buffer.data(dpf + 79);
    const auto *dpf_80 = buffer.data(dpf + 80);
    const auto *dpf_85 = buffer.data(dpf + 85);
    const auto *dpf_86 = buffer.data(dpf + 86);
    const auto *dpf_87 = buffer.data(dpf + 87);

    const auto *dpg1_0 = buffer.data(dpg1 + 0);
    const auto *dpg1_3 = buffer.data(dpg1 + 3);
    const auto *dpg1_5 = buffer.data(dpg1 + 5);
    const auto *dpg1_6 = buffer.data(dpg1 + 6);
    const auto *dpg1_9 = buffer.data(dpg1 + 9);
    const auto *dpg1_10 = buffer.data(dpg1 + 10);
    const auto *dpg1_14 = buffer.data(dpg1 + 14);
    const auto *dpg1_15 = buffer.data(dpg1 + 15);
    const auto *dpg1_18 = buffer.data(dpg1 + 18);
    const auto *dpg1_25 = buffer.data(dpg1 + 25);
    const auto *dpg1_30 = buffer.data(dpg1 + 30);
    const auto *dpg1_35 = buffer.data(dpg1 + 35);
    const auto *dpg1_44 = buffer.data(dpg1 + 44);
    const auto *dpg1_70 = buffer.data(dpg1 + 70);

    const auto *fsg0_0 = buffer.data(fsg0 + 0);
    const auto *fsg0_3 = buffer.data(fsg0 + 3);
    const auto *fsg0_5 = buffer.data(fsg0 + 5);
    const auto *fsg0_10 = buffer.data(fsg0 + 10);
    const auto *fsg0_12 = buffer.data(fsg0 + 12);
    const auto *fsg0_14 = buffer.data(fsg0 + 14);
    const auto *fsg0_18 = buffer.data(fsg0 + 18);
    const auto *fsg0_25 = buffer.data(fsg0 + 25);
    const auto *fsg0_27 = buffer.data(fsg0 + 27);
    const auto *fsg0_35 = buffer.data(fsg0 + 35);
    const auto *fsg0_41 = buffer.data(fsg0 + 41);
    const auto *fsg0_42 = buffer.data(fsg0 + 42);
    const auto *fsg0_44 = buffer.data(fsg0 + 44);

    const auto *fsf_0 = buffer.data(fsf + 0);
    const auto *fsf_1 = buffer.data(fsf + 1);
    const auto *fsf_2 = buffer.data(fsf + 2);
    const auto *fsf_3 = buffer.data(fsf + 3);
    const auto *fsf_5 = buffer.data(fsf + 5);
    const auto *fsf_6 = buffer.data(fsf + 6);
    const auto *fsf_8 = buffer.data(fsf + 8);
    const auto *fsf_9 = buffer.data(fsf + 9);
    const auto *fsf_10 = buffer.data(fsf + 10);
    const auto *fsf_11 = buffer.data(fsf + 11);
    const auto *fsf_13 = buffer.data(fsf + 13);
    const auto *fsf_16 = buffer.data(fsf + 16);
    const auto *fsf_17 = buffer.data(fsf + 17);
    const auto *fsf_18 = buffer.data(fsf + 18);
    const auto *fsf_19 = buffer.data(fsf + 19);
    const auto *fsf_20 = buffer.data(fsf + 20);
    const auto *fsf_22 = buffer.data(fsf + 22);
    const auto *fsf_25 = buffer.data(fsf + 25);
    const auto *fsf_27 = buffer.data(fsf + 27);
    const auto *fsf_28 = buffer.data(fsf + 28);
    const auto *fsf_29 = buffer.data(fsf + 29);

    const auto *fsg1_0 = buffer.data(fsg1 + 0);
    const auto *fsg1_3 = buffer.data(fsg1 + 3);
    const auto *fsg1_5 = buffer.data(fsg1 + 5);
    const auto *fsg1_10 = buffer.data(fsg1 + 10);
    const auto *fsg1_12 = buffer.data(fsg1 + 12);
    const auto *fsg1_14 = buffer.data(fsg1 + 14);
    const auto *fsg1_18 = buffer.data(fsg1 + 18);
    const auto *fsg1_25 = buffer.data(fsg1 + 25);
    const auto *fsg1_27 = buffer.data(fsg1 + 27);
    const auto *fsg1_35 = buffer.data(fsg1 + 35);
    const auto *fsg1_41 = buffer.data(fsg1 + 41);
    const auto *fsg1_42 = buffer.data(fsg1 + 42);
    const auto *fsg1_44 = buffer.data(fsg1 + 44);

    const auto *fpd0_0 = buffer.data(fpd0 + 0);
    const auto *fpd0_3 = buffer.data(fpd0 + 3);
    const auto *fpd0_5 = buffer.data(fpd0 + 5);
    const auto *fpd0_17 = buffer.data(fpd0 + 17);
    const auto *fpd0_21 = buffer.data(fpd0 + 21);
    const auto *fpd0_24 = buffer.data(fpd0 + 24);
    const auto *fpd0_27 = buffer.data(fpd0 + 27);
    const auto *fpd0_29 = buffer.data(fpd0 + 29);
    const auto *fpd0_40 = buffer.data(fpd0 + 40);
    const auto *fpd0_41 = buffer.data(fpd0 + 41);
    const auto *fpd0_48 = buffer.data(fpd0 + 48);
    const auto *fpd0_53 = buffer.data(fpd0 + 53);

    const auto *fpd1_0 = buffer.data(fpd1 + 0);
    const auto *fpd1_3 = buffer.data(fpd1 + 3);
    const auto *fpd1_5 = buffer.data(fpd1 + 5);
    const auto *fpd1_17 = buffer.data(fpd1 + 17);
    const auto *fpd1_21 = buffer.data(fpd1 + 21);
    const auto *fpd1_24 = buffer.data(fpd1 + 24);
    const auto *fpd1_27 = buffer.data(fpd1 + 27);
    const auto *fpd1_29 = buffer.data(fpd1 + 29);
    const auto *fpd1_40 = buffer.data(fpd1 + 40);
    const auto *fpd1_41 = buffer.data(fpd1 + 41);
    const auto *fpd1_48 = buffer.data(fpd1 + 48);
    const auto *fpd1_53 = buffer.data(fpd1 + 53);

    const auto *fpf_0 = buffer.data(fpf + 0);
    const auto *fpf_1 = buffer.data(fpf + 1);
    const auto *fpf_2 = buffer.data(fpf + 2);
    const auto *fpf_3 = buffer.data(fpf + 3);
    const auto *fpf_5 = buffer.data(fpf + 5);
    const auto *fpf_6 = buffer.data(fpf + 6);
    const auto *fpf_8 = buffer.data(fpf + 8);
    const auto *fpf_9 = buffer.data(fpf + 9);
    const auto *fpf_10 = buffer.data(fpf + 10);
    const auto *fpf_12 = buffer.data(fpf + 12);
    const auto *fpf_13 = buffer.data(fpf + 13);
    const auto *fpf_15 = buffer.data(fpf + 15);
    const auto *fpf_16 = buffer.data(fpf + 16);
    const auto *fpf_19 = buffer.data(fpf + 19);
    const auto *fpf_20 = buffer.data(fpf + 20);
    const auto *fpf_22 = buffer.data(fpf + 22);
    const auto *fpf_23 = buffer.data(fpf + 23);
    const auto *fpf_25 = buffer.data(fpf + 25);
    const auto *fpf_26 = buffer.data(fpf + 26);
    const auto *fpf_28 = buffer.data(fpf + 28);
    const auto *fpf_29 = buffer.data(fpf + 29);
    const auto *fpf_30 = buffer.data(fpf + 30);
    const auto *fpf_31 = buffer.data(fpf + 31);
    const auto *fpf_33 = buffer.data(fpf + 33);
    const auto *fpf_36 = buffer.data(fpf + 36);
    const auto *fpf_37 = buffer.data(fpf + 37);
    const auto *fpf_38 = buffer.data(fpf + 38);
    const auto *fpf_39 = buffer.data(fpf + 39);
    const auto *fpf_40 = buffer.data(fpf + 40);
    const auto *fpf_41 = buffer.data(fpf + 41);
    const auto *fpf_42 = buffer.data(fpf + 42);
    const auto *fpf_43 = buffer.data(fpf + 43);
    const auto *fpf_46 = buffer.data(fpf + 46);
    const auto *fpf_47 = buffer.data(fpf + 47);
    const auto *fpf_48 = buffer.data(fpf + 48);
    const auto *fpf_49 = buffer.data(fpf + 49);
    const auto *fpf_50 = buffer.data(fpf + 50);
    const auto *fpf_51 = buffer.data(fpf + 51);
    const auto *fpf_53 = buffer.data(fpf + 53);
    const auto *fpf_56 = buffer.data(fpf + 56);
    const auto *fpf_58 = buffer.data(fpf + 58);
    const auto *fpf_59 = buffer.data(fpf + 59);
    const auto *fpf_60 = buffer.data(fpf + 60);
    const auto *fpf_62 = buffer.data(fpf + 62);
    const auto *fpf_65 = buffer.data(fpf + 65);
    const auto *fpf_67 = buffer.data(fpf + 67);
    const auto *fpf_68 = buffer.data(fpf + 68);
    const auto *fpf_69 = buffer.data(fpf + 69);
    const auto *fpf_70 = buffer.data(fpf + 70);
    const auto *fpf_72 = buffer.data(fpf + 72);
    const auto *fpf_75 = buffer.data(fpf + 75);
    const auto *fpf_76 = buffer.data(fpf + 76);
    const auto *fpf_77 = buffer.data(fpf + 77);
    const auto *fpf_79 = buffer.data(fpf + 79);
    const auto *fpf_80 = buffer.data(fpf + 80);
    const auto *fpf_81 = buffer.data(fpf + 81);
    const auto *fpf_82 = buffer.data(fpf + 82);
    const auto *fpf_85 = buffer.data(fpf + 85);
    const auto *fpf_86 = buffer.data(fpf + 86);
    const auto *fpf_87 = buffer.data(fpf + 87);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dpf_0, fsf_0, fpd0_0, \
                         fpd1_0, fpf_0, fpf_1, fpf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dpf_0[k]
                 + f_1 * fsf_0[k]
                 + f_2 * fpd0_0[k]
                 - f_3 * fpd1_0[k]
                 + f_4 * pc_x[k] * fpf_0[k];

        t_1[k] = f_4 * pc_y[k] * fpf_0[k];

        t_2[k] = f_4 * pc_z[k] * fpf_0[k];

        t_3[k] = f_5 * fpd0_0[k]
                 - f_6 * fpd1_0[k]
                 + f_4 * pc_y[k] * fpf_1[k];

        t_4[k] = f_4 * pc_y[k] * fpf_2[k];

        t_5[k] = f_5 * fpd0_0[k]
                 - f_6 * fpd1_0[k]
                 + f_4 * pc_z[k] * fpf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, dpf_6, dpf_9, fsf_6, fsf_9, \
                         fpf_3, fpf_5, fpf_6, fpf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * dpf_6[k]
                 + f_1 * fsf_6[k]
                 + f_4 * pc_x[k] * fpf_6[k];

        t_7[k] = f_4 * pc_z[k] * fpf_3[k];

        t_8[k] = f_4 * pc_y[k] * fpf_5[k];

        t_9[k] = f_0 * dpf_9[k]
                 + f_1 * fsf_9[k]
                 + f_4 * pc_x[k] * fpf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pc_y, pc_z, fpd0_3, fpd0_5, fpd1_3, \
                         fpd1_5, fpf_6, fpf_8, fpf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * fpd0_3[k]
                  - f_3 * fpd1_3[k]
                  + f_4 * pc_y[k] * fpf_6[k];

        t_11[k] = f_4 * pc_z[k] * fpf_6[k];

        t_12[k] = f_5 * fpd0_5[k]
                  - f_6 * fpd1_5[k]
                  + f_4 * pc_y[k] * fpf_8[k];

        t_13[k] = f_4 * pc_y[k] * fpf_9[k];

        t_14[k] = f_2 * fpd0_5[k]
                  - f_3 * fpd1_5[k]
                  + f_4 * pc_z[k] * fpf_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pc_y, pc_z, fsg0_0, fsg0_3, fsf_0, \
                         fsf_1, fsg1_0, fsg1_3, fpf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * fsg0_0[k]
                  - f_7 * pc_y[k] * fsg1_0[k];

        t_16[k] = f_1 * fsf_0[k]
                  + f_4 * pc_y[k] * fpf_10[k];

        t_17[k] = f_4 * pc_z[k] * fpf_10[k];

        t_18[k] = pb_y[k] * fsg0_3[k]
                  + f_8 * fsf_1[k]
                  - f_7 * pc_y[k] * fsg1_3[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_y, pc_x, pc_y, pc_z, dpf_16, fsg0_5, \
                         fsf_2, fsg1_5, fpf_12, fpf_13, fpf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * fsf_2[k]
                  + f_4 * pc_y[k] * fpf_12[k];

        t_20[k] = pb_y[k] * fsg0_5[k]
                  - f_7 * pc_y[k] * fsg1_5[k];

        t_21[k] = f_0 * dpf_16[k]
                  + f_4 * pc_x[k] * fpf_16[k];

        t_22[k] = f_4 * pc_z[k] * fpf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_x, pc_y, pc_z, dpf_19, fsg0_10, \
                         fsf_5, fsf_6, fsg1_10, fpf_15, fpf_16, \
                         fpf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * fsf_5[k]
                  + f_4 * pc_y[k] * fpf_15[k];

        t_24[k] = f_0 * dpf_19[k]
                  + f_4 * pc_x[k] * fpf_19[k];

        t_25[k] = pb_y[k] * fsg0_10[k]
                  + f_9 * fsf_6[k]
                  - f_7 * pc_y[k] * fsg1_10[k];

        t_26[k] = f_4 * pc_z[k] * fpf_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pc_y, fsg0_12, fsg0_14, fsf_8, fsf_9, \
                         fsg1_12, fsg1_14, fpf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * fsg0_12[k]
                  + f_8 * fsf_8[k]
                  - f_7 * pc_y[k] * fsg1_12[k];

        t_28[k] = f_1 * fsf_9[k]
                  + f_4 * pc_y[k] * fpf_19[k];

        t_29[k] = pb_y[k] * fsg0_14[k]
                  - f_7 * pc_y[k] * fsg1_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_z, pc_y, pc_z, fsg0_0, fsg0_3, \
                         fsf_0, fsg1_0, fsg1_3, fpf_20, fpf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * fsg0_0[k]
                  - f_7 * pc_z[k] * fsg1_0[k];

        t_31[k] = f_4 * pc_y[k] * fpf_20[k];

        t_32[k] = f_1 * fsf_0[k]
                  + f_4 * pc_z[k] * fpf_20[k];

        t_33[k] = pb_z[k] * fsg0_3[k]
                  - f_7 * pc_z[k] * fsg1_3[k];

        t_34[k] = f_4 * pc_y[k] * fpf_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_z, pc_x, pc_y, pc_z, dpf_26, fsg0_5, \
                         fsf_2, fsf_3, fsg1_5, fpf_23, fpf_25, fpf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_z[k] * fsg0_5[k]
                  + f_8 * fsf_2[k]
                  - f_7 * pc_z[k] * fsg1_5[k];

        t_36[k] = f_0 * dpf_26[k]
                  + f_4 * pc_x[k] * fpf_26[k];

        t_37[k] = f_1 * fsf_3[k]
                  + f_4 * pc_z[k] * fpf_23[k];

        t_38[k] = f_4 * pc_y[k] * fpf_25[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_z, pc_x, pc_z, dpf_29, fsg0_10, fsf_6, fsg1_10, \
                         fpf_26, fpf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * dpf_29[k]
                  + f_4 * pc_x[k] * fpf_29[k];

        t_40[k] = pb_z[k] * fsg0_10[k]
                  - f_7 * pc_z[k] * fsg1_10[k];

        t_41[k] = f_1 * fsf_6[k]
                  + f_4 * pc_z[k] * fpf_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_z, pc_y, pc_z, fsg0_14, fsf_9, fsg1_14, fpd0_17, \
                         fpd1_17, fpf_28, fpf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * fpd0_17[k]
                  - f_6 * fpd1_17[k]
                  + f_4 * pc_y[k] * fpf_28[k];

        t_43[k] = f_4 * pc_y[k] * fpf_29[k];

        t_44[k] = pb_z[k] * fsg0_14[k]
                  + f_9 * fsf_9[k]
                  - f_7 * pc_z[k] * fsg1_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, dpg0_0, dpg0_3, \
                         dpf_0, dpf_1, dpg1_0, dpg1_3, fpf_30, fpf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * dpg0_0[k]
                  - f_7 * pc_y[k] * dpg1_0[k];

        t_46[k] = f_1 * dpf_0[k]
                  + f_4 * pc_y[k] * fpf_30[k];

        t_47[k] = f_4 * pc_z[k] * fpf_30[k];

        t_48[k] = pa_y[k] * dpg0_3[k]
                  + f_8 * dpf_1[k]
                  - f_7 * pc_y[k] * dpg1_3[k];

        t_49[k] = f_4 * pc_z[k] * fpf_31[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_y, pc_x, pc_y, pc_z, dpg0_5, dpf_36, dpg1_5, \
                         fsf_16, fpf_33, fpf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * dpg0_5[k]
                  - f_7 * pc_y[k] * dpg1_5[k];

        t_51[k] = f_8 * dpf_36[k]
                  + f_1 * fsf_16[k]
                  + f_4 * pc_x[k] * fpf_36[k];

        t_52[k] = f_4 * pc_z[k] * fpf_33[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pc_x, pc_y, dpg0_9, dpf_6, dpf_38, dpg1_9, \
                         fsf_18, fpd0_21, fpd1_21, fpf_36, fpf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_8 * dpf_38[k]
                  + f_1 * fsf_18[k]
                  + f_4 * pc_x[k] * fpf_38[k];

        t_54[k] = pa_y[k] * dpg0_9[k]
                  - f_7 * pc_y[k] * dpg1_9[k];

        t_55[k] = f_1 * dpf_6[k]
                  + f_2 * fpd0_21[k]
                  - f_3 * fpd1_21[k]
                  + f_4 * pc_y[k] * fpf_36[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pc_y, pc_z, dpg0_14, dpf_9, dpg1_14, \
                         fpd0_21, fpd1_21, fpf_36, fpf_37, fpf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_4 * pc_z[k] * fpf_36[k];

        t_57[k] = f_5 * fpd0_21[k]
                  - f_6 * fpd1_21[k]
                  + f_4 * pc_z[k] * fpf_37[k];

        t_58[k] = f_1 * dpf_9[k]
                  + f_4 * pc_y[k] * fpf_39[k];

        t_59[k] = pa_y[k] * dpg0_14[k]
                  - f_7 * pc_y[k] * dpg1_14[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pc_x, pc_y, pc_z, dpf_10, dpf_40, fsf_10, fpd0_24, \
                         fpd1_24, fpf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_8 * dpf_40[k]
                  + f_2 * fpd0_24[k]
                  - f_3 * fpd1_24[k]
                  + f_4 * pc_x[k] * fpf_40[k];

        t_61[k] = f_1 * dpf_10[k]
                  + f_1 * fsf_10[k]
                  + f_4 * pc_y[k] * fpf_40[k];

        t_62[k] = f_4 * pc_z[k] * fpf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_z, dpf_43, dpf_46, fpd0_24, fpd0_27, \
                         fpd1_24, fpd1_27, fpf_41, fpf_42, fpf_43, \
                         fpf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_8 * dpf_43[k]
                  + f_5 * fpd0_27[k]
                  - f_6 * fpd1_27[k]
                  + f_4 * pc_x[k] * fpf_43[k];

        t_64[k] = f_4 * pc_z[k] * fpf_41[k];

        t_65[k] = f_5 * fpd0_24[k]
                  - f_6 * fpd1_24[k]
                  + f_4 * pc_z[k] * fpf_42[k];

        t_66[k] = f_8 * dpf_46[k]
                  + f_4 * pc_x[k] * fpf_46[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_x, pc_x, pc_z, ppg0_70, ppg1_70, dpg0_70, \
                         dpf_48, dpf_49, dpg1_70, fpf_43, fpf_48, \
                         fpf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_4 * pc_z[k] * fpf_43[k];

        t_68[k] = f_8 * dpf_48[k]
                  + f_4 * pc_x[k] * fpf_48[k];

        t_69[k] = f_8 * dpf_49[k]
                  + f_4 * pc_x[k] * fpf_49[k];

        t_70[k] = f_10 * ppg0_70[k]
                  - f_11 * ppg1_70[k]
                  + pa_x[k] * dpg0_70[k]
                  - f_7 * pc_x[k] * dpg1_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_y, pc_z, dpf_19, fsf_19, fpd0_27, fpd0_29, \
                         fpd1_27, fpd1_29, fpf_46, fpf_47, fpf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_4 * pc_z[k] * fpf_46[k];

        t_72[k] = f_5 * fpd0_27[k]
                  - f_6 * fpd1_27[k]
                  + f_4 * pc_z[k] * fpf_47[k];

        t_73[k] = f_1 * dpf_19[k]
                  + f_1 * fsf_19[k]
                  + f_4 * pc_y[k] * fpf_49[k];

        t_74[k] = f_2 * fpd0_29[k]
                  - f_3 * fpd1_29[k]
                  + f_4 * pc_z[k] * fpf_49[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pb_z, pc_y, pc_z, dpg0_30, dpf_20, \
                         dpg1_30, fsg0_18, fsf_10, fsg1_18, fpf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pa_y[k] * dpg0_30[k]
                  - f_7 * pc_y[k] * dpg1_30[k];

        t_76[k] = f_1 * dpf_20[k]
                  + f_4 * pc_y[k] * fpf_50[k];

        t_77[k] = f_1 * fsf_10[k]
                  + f_4 * pc_z[k] * fpf_50[k];

        t_78[k] = pb_z[k] * fsg0_18[k]
                  - f_7 * pc_z[k] * fsg1_18[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_y, pc_x, pc_y, pc_z, dpg0_35, dpf_56, \
                         dpg1_35, fsf_11, fsf_13, fpf_51, fpf_53, \
                         fpf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * fsf_11[k]
                  + f_4 * pc_z[k] * fpf_51[k];

        t_80[k] = pa_y[k] * dpg0_35[k]
                  - f_7 * pc_y[k] * dpg1_35[k];

        t_81[k] = f_8 * dpf_56[k]
                  + f_4 * pc_x[k] * fpf_56[k];

        t_82[k] = f_1 * fsf_13[k]
                  + f_4 * pc_z[k] * fpf_53[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pb_z, pc_x, pc_z, dpf_58, dpf_59, fsg0_25, \
                         fsf_16, fsg1_25, fpf_56, fpf_58, fpf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_8 * dpf_58[k]
                  + f_4 * pc_x[k] * fpf_58[k];

        t_84[k] = f_8 * dpf_59[k]
                  + f_4 * pc_x[k] * fpf_59[k];

        t_85[k] = pb_z[k] * fsg0_25[k]
                  - f_7 * pc_z[k] * fsg1_25[k];

        t_86[k] = f_1 * fsf_16[k]
                  + f_4 * pc_z[k] * fpf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_y, pb_z, pc_y, pc_z, dpg0_44, dpf_29, dpg1_44, \
                         fsg0_27, fsf_17, fsg1_27, fpf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_z[k] * fsg0_27[k]
                  + f_8 * fsf_17[k]
                  - f_7 * pc_z[k] * fsg1_27[k];

        t_88[k] = f_1 * dpf_29[k]
                  + f_4 * pc_y[k] * fpf_59[k];

        t_89[k] = pa_y[k] * dpg0_44[k]
                  - f_7 * pc_y[k] * dpg1_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_z, pc_y, pc_z, dpg0_0, dpg0_3, \
                         dpf_0, dpg1_0, dpg1_3, fpf_60, fpf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_z[k] * dpg0_0[k]
                  - f_7 * pc_z[k] * dpg1_0[k];

        t_91[k] = f_4 * pc_y[k] * fpf_60[k];

        t_92[k] = f_1 * dpf_0[k]
                  + f_4 * pc_z[k] * fpf_60[k];

        t_93[k] = pa_z[k] * dpg0_3[k]
                  - f_7 * pc_z[k] * dpg1_3[k];

        t_94[k] = f_4 * pc_y[k] * fpf_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pa_z, pc_x, pc_z, dpg0_5, dpg0_6, dpf_2, dpf_67, \
                         dpg1_5, dpg1_6, fsf_27, fpf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pa_z[k] * dpg0_5[k]
                  + f_8 * dpf_2[k]
                  - f_7 * pc_z[k] * dpg1_5[k];

        t_96[k] = pa_z[k] * dpg0_6[k]
                  - f_7 * pc_z[k] * dpg1_6[k];

        t_97[k] = f_8 * dpf_67[k]
                  + f_1 * fsf_27[k]
                  + f_4 * pc_x[k] * fpf_67[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_z, pc_x, pc_y, pc_z, dpg0_10, dpf_69, dpg1_10, \
                         fsf_29, fpf_65, fpf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_4 * pc_y[k] * fpf_65[k];

        t_99[k] = f_8 * dpf_69[k]
                  + f_1 * fsf_29[k]
                  + f_4 * pc_x[k] * fpf_69[k];

        t_100[k] = pa_z[k] * dpg0_10[k]
                   - f_7 * pc_z[k] * dpg1_10[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_y, pc_z, dpf_9, fpd0_40, fpd0_41, \
                         fpd1_40, fpd1_41, fpf_67, fpf_68, fpf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_12 * fpd0_40[k]
                   - f_13 * fpd1_40[k]
                   + f_4 * pc_y[k] * fpf_67[k];

        t_102[k] = f_5 * fpd0_41[k]
                   - f_6 * fpd1_41[k]
                   + f_4 * pc_y[k] * fpf_68[k];

        t_103[k] = f_4 * pc_y[k] * fpf_69[k];

        t_104[k] = f_1 * dpf_9[k]
                   + f_2 * fpd0_41[k]
                   - f_3 * fpd1_41[k]
                   + f_4 * pc_z[k] * fpf_69[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_z, pc_y, pc_z, dpg0_15, dpg0_18, \
                         dpf_10, dpg1_15, dpg1_18, fsf_20, fpf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_z[k] * dpg0_15[k]
                   - f_7 * pc_z[k] * dpg1_15[k];

        t_106[k] = f_1 * fsf_20[k]
                   + f_4 * pc_y[k] * fpf_70[k];

        t_107[k] = f_1 * dpf_10[k]
                   + f_4 * pc_z[k] * fpf_70[k];

        t_108[k] = pa_z[k] * dpg0_18[k]
                   - f_7 * pc_z[k] * dpg1_18[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pb_y, pc_x, pc_y, dpf_76, dpf_77, \
                         fsg0_35, fsf_22, fsg1_35, fpf_72, fpf_76, \
                         fpf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_1 * fsf_22[k]
                   + f_4 * pc_y[k] * fpf_72[k];

        t_110[k] = pb_y[k] * fsg0_35[k]
                   - f_7 * pc_y[k] * fsg1_35[k];

        t_111[k] = f_8 * dpf_76[k]
                   + f_4 * pc_x[k] * fpf_76[k];

        t_112[k] = f_8 * dpf_77[k]
                   + f_4 * pc_x[k] * fpf_77[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_z, pc_x, pc_y, pc_z, dpg0_25, dpf_79, \
                         dpg1_25, fsf_25, fpf_75, fpf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_1 * fsf_25[k]
                   + f_4 * pc_y[k] * fpf_75[k];

        t_114[k] = f_8 * dpf_79[k]
                   + f_4 * pc_x[k] * fpf_79[k];

        t_115[k] = pa_z[k] * dpg0_25[k]
                   - f_7 * pc_z[k] * dpg1_25[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pc_y, fsg0_41, fsg0_42, fsg0_44, \
                         fsf_27, fsf_28, fsf_29, fsg1_41, fsg1_42, fsg1_44, \
                         fpf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pb_y[k] * fsg0_41[k]
                   + f_0 * fsf_27[k]
                   - f_7 * pc_y[k] * fsg1_41[k];

        t_117[k] = pb_y[k] * fsg0_42[k]
                   + f_8 * fsf_28[k]
                   - f_7 * pc_y[k] * fsg1_42[k];

        t_118[k] = f_1 * fsf_29[k]
                   + f_4 * pc_y[k] * fpf_79[k];

        t_119[k] = pb_y[k] * fsg0_44[k]
                   - f_7 * pc_y[k] * fsg1_44[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, dpf_20, dpf_80, \
                         fsf_20, fpd0_48, fpd1_48, fpf_80, fpf_81, \
                         fpf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_8 * dpf_80[k]
                   + f_2 * fpd0_48[k]
                   - f_3 * fpd1_48[k]
                   + f_4 * pc_x[k] * fpf_80[k];

        t_121[k] = f_4 * pc_y[k] * fpf_80[k];

        t_122[k] = f_1 * dpf_20[k]
                   + f_1 * fsf_20[k]
                   + f_4 * pc_z[k] * fpf_80[k];

        t_123[k] = f_5 * fpd0_48[k]
                   - f_6 * fpd1_48[k]
                   + f_4 * pc_y[k] * fpf_81[k];

        t_124[k] = f_4 * pc_y[k] * fpf_82[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_y, dpf_85, dpf_86, dpf_87, \
                         fpd0_53, fpd1_53, fpf_85, fpf_86, fpf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_8 * dpf_85[k]
                   + f_5 * fpd0_53[k]
                   - f_6 * fpd1_53[k]
                   + f_4 * pc_x[k] * fpf_85[k];

        t_126[k] = f_8 * dpf_86[k]
                   + f_4 * pc_x[k] * fpf_86[k];

        t_127[k] = f_8 * dpf_87[k]
                   + f_4 * pc_x[k] * fpf_87[k];

        t_128[k] = f_4 * pc_y[k] * fpf_85[k];
    }
}

static auto
compute_prim_fpg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppg0, const size_t ppg1,
                                                          const size_t dpg0, const size_t dpf,
                                                          const size_t dpg1, const size_t fsg0,
                                                          const size_t fsf, const size_t fsg1,
                                                          const size_t fpd0, const size_t fpd1,
                                                          const size_t fpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppg0_134 = buffer.data(ppg0 + 134);

    const auto *ppg1_134 = buffer.data(ppg1 + 134);

    const auto *dpg0_48 = buffer.data(dpg0 + 48);
    const auto *dpg0_51 = buffer.data(dpg0 + 51);
    const auto *dpg0_61 = buffer.data(dpg0 + 61);
    const auto *dpg0_63 = buffer.data(dpg0 + 63);
    const auto *dpg0_90 = buffer.data(dpg0 + 90);
    const auto *dpg0_95 = buffer.data(dpg0 + 95);
    const auto *dpg0_99 = buffer.data(dpg0 + 99);
    const auto *dpg0_120 = buffer.data(dpg0 + 120);
    const auto *dpg0_122 = buffer.data(dpg0 + 122);
    const auto *dpg0_125 = buffer.data(dpg0 + 125);
    const auto *dpg0_134 = buffer.data(dpg0 + 134);
    const auto *dpg0_150 = buffer.data(dpg0 + 150);
    const auto *dpg0_153 = buffer.data(dpg0 + 153);
    const auto *dpg0_160 = buffer.data(dpg0 + 160);
    const auto *dpg0_162 = buffer.data(dpg0 + 162);
    const auto *dpg0_163 = buffer.data(dpg0 + 163);
    const auto *dpg0_164 = buffer.data(dpg0 + 164);
    const auto *dpg0_170 = buffer.data(dpg0 + 170);
    const auto *dpg0_175 = buffer.data(dpg0 + 175);
    const auto *dpg0_177 = buffer.data(dpg0 + 177);
    const auto *dpg0_179 = buffer.data(dpg0 + 179);
    const auto *dpg0_205 = buffer.data(dpg0 + 205);
    const auto *dpg0_207 = buffer.data(dpg0 + 207);
    const auto *dpg0_208 = buffer.data(dpg0 + 208);
    const auto *dpg0_209 = buffer.data(dpg0 + 209);
    const auto *dpg0_220 = buffer.data(dpg0 + 220);
    const auto *dpg0_221 = buffer.data(dpg0 + 221);
    const auto *dpg0_222 = buffer.data(dpg0 + 222);
    const auto *dpg0_224 = buffer.data(dpg0 + 224);
    const auto *dpg0_243 = buffer.data(dpg0 + 243);

    const auto *dpf_30 = buffer.data(dpf + 30);
    const auto *dpf_36 = buffer.data(dpf + 36);
    const auto *dpf_39 = buffer.data(dpf + 39);
    const auto *dpf_40 = buffer.data(dpf + 40);
    const auto *dpf_46 = buffer.data(dpf + 46);
    const auto *dpf_50 = buffer.data(dpf + 50);
    const auto *dpf_59 = buffer.data(dpf + 59);
    const auto *dpf_60 = buffer.data(dpf + 60);
    const auto *dpf_62 = buffer.data(dpf + 62);
    const auto *dpf_66 = buffer.data(dpf + 66);
    const auto *dpf_68 = buffer.data(dpf + 68);
    const auto *dpf_69 = buffer.data(dpf + 69);
    const auto *dpf_70 = buffer.data(dpf + 70);
    const auto *dpf_72 = buffer.data(dpf + 72);
    const auto *dpf_80 = buffer.data(dpf + 80);
    const auto *dpf_82 = buffer.data(dpf + 82);
    const auto *dpf_89 = buffer.data(dpf + 89);
    const auto *dpf_90 = buffer.data(dpf + 90);
    const auto *dpf_93 = buffer.data(dpf + 93);
    const auto *dpf_96 = buffer.data(dpf + 96);
    const auto *dpf_98 = buffer.data(dpf + 98);
    const auto *dpf_99 = buffer.data(dpf + 99);
    const auto *dpf_100 = buffer.data(dpf + 100);
    const auto *dpf_103 = buffer.data(dpf + 103);
    const auto *dpf_106 = buffer.data(dpf + 106);
    const auto *dpf_108 = buffer.data(dpf + 108);
    const auto *dpf_109 = buffer.data(dpf + 109);
    const auto *dpf_115 = buffer.data(dpf + 115);
    const auto *dpf_116 = buffer.data(dpf + 116);
    const auto *dpf_118 = buffer.data(dpf + 118);
    const auto *dpf_119 = buffer.data(dpf + 119);
    const auto *dpf_127 = buffer.data(dpf + 127);
    const auto *dpf_128 = buffer.data(dpf + 128);
    const auto *dpf_130 = buffer.data(dpf + 130);
    const auto *dpf_135 = buffer.data(dpf + 135);
    const auto *dpf_136 = buffer.data(dpf + 136);
    const auto *dpf_137 = buffer.data(dpf + 137);
    const auto *dpf_138 = buffer.data(dpf + 138);
    const auto *dpf_139 = buffer.data(dpf + 139);
    const auto *dpf_143 = buffer.data(dpf + 143);
    const auto *dpf_146 = buffer.data(dpf + 146);
    const auto *dpf_147 = buffer.data(dpf + 147);
    const auto *dpf_148 = buffer.data(dpf + 148);
    const auto *dpf_149 = buffer.data(dpf + 149);
    const auto *dpf_150 = buffer.data(dpf + 150);
    const auto *dpf_155 = buffer.data(dpf + 155);
    const auto *dpf_156 = buffer.data(dpf + 156);
    const auto *dpf_157 = buffer.data(dpf + 157);
    const auto *dpf_159 = buffer.data(dpf + 159);
    const auto *dpf_163 = buffer.data(dpf + 163);
    const auto *dpf_166 = buffer.data(dpf + 166);
    const auto *dpf_167 = buffer.data(dpf + 167);
    const auto *dpf_169 = buffer.data(dpf + 169);

    const auto *dpg1_48 = buffer.data(dpg1 + 48);
    const auto *dpg1_51 = buffer.data(dpg1 + 51);
    const auto *dpg1_61 = buffer.data(dpg1 + 61);
    const auto *dpg1_63 = buffer.data(dpg1 + 63);
    const auto *dpg1_90 = buffer.data(dpg1 + 90);
    const auto *dpg1_95 = buffer.data(dpg1 + 95);
    const auto *dpg1_99 = buffer.data(dpg1 + 99);
    const auto *dpg1_120 = buffer.data(dpg1 + 120);
    const auto *dpg1_122 = buffer.data(dpg1 + 122);
    const auto *dpg1_125 = buffer.data(dpg1 + 125);
    const auto *dpg1_134 = buffer.data(dpg1 + 134);
    const auto *dpg1_150 = buffer.data(dpg1 + 150);
    const auto *dpg1_153 = buffer.data(dpg1 + 153);
    const auto *dpg1_160 = buffer.data(dpg1 + 160);
    const auto *dpg1_162 = buffer.data(dpg1 + 162);
    const auto *dpg1_163 = buffer.data(dpg1 + 163);
    const auto *dpg1_164 = buffer.data(dpg1 + 164);
    const auto *dpg1_170 = buffer.data(dpg1 + 170);
    const auto *dpg1_175 = buffer.data(dpg1 + 175);
    const auto *dpg1_177 = buffer.data(dpg1 + 177);
    const auto *dpg1_179 = buffer.data(dpg1 + 179);
    const auto *dpg1_205 = buffer.data(dpg1 + 205);
    const auto *dpg1_207 = buffer.data(dpg1 + 207);
    const auto *dpg1_208 = buffer.data(dpg1 + 208);
    const auto *dpg1_209 = buffer.data(dpg1 + 209);
    const auto *dpg1_220 = buffer.data(dpg1 + 220);
    const auto *dpg1_221 = buffer.data(dpg1 + 221);
    const auto *dpg1_222 = buffer.data(dpg1 + 222);
    const auto *dpg1_224 = buffer.data(dpg1 + 224);
    const auto *dpg1_243 = buffer.data(dpg1 + 243);

    const auto *fsg0_45 = buffer.data(fsg0 + 45);
    const auto *fsg0_48 = buffer.data(fsg0 + 48);
    const auto *fsg0_75 = buffer.data(fsg0 + 75);
    const auto *fsg0_80 = buffer.data(fsg0 + 80);

    const auto *fsf_30 = buffer.data(fsf + 30);
    const auto *fsf_31 = buffer.data(fsf + 31);
    const auto *fsf_33 = buffer.data(fsf + 33);
    const auto *fsf_36 = buffer.data(fsf + 36);
    const auto *fsf_38 = buffer.data(fsf + 38);
    const auto *fsf_39 = buffer.data(fsf + 39);
    const auto *fsf_42 = buffer.data(fsf + 42);
    const auto *fsf_47 = buffer.data(fsf + 47);
    const auto *fsf_48 = buffer.data(fsf + 48);
    const auto *fsf_50 = buffer.data(fsf + 50);
    const auto *fsf_52 = buffer.data(fsf + 52);
    const auto *fsf_55 = buffer.data(fsf + 55);
    const auto *fsf_56 = buffer.data(fsf + 56);
    const auto *fsf_57 = buffer.data(fsf + 57);
    const auto *fsf_59 = buffer.data(fsf + 59);

    const auto *fsg1_45 = buffer.data(fsg1 + 45);
    const auto *fsg1_48 = buffer.data(fsg1 + 48);
    const auto *fsg1_75 = buffer.data(fsg1 + 75);
    const auto *fsg1_80 = buffer.data(fsg1 + 80);

    const auto *fpd0_51 = buffer.data(fpd0 + 51);
    const auto *fpd0_52 = buffer.data(fpd0 + 52);
    const auto *fpd0_53 = buffer.data(fpd0 + 53);
    const auto *fpd0_54 = buffer.data(fpd0 + 54);
    const auto *fpd0_57 = buffer.data(fpd0 + 57);
    const auto *fpd0_59 = buffer.data(fpd0 + 59);
    const auto *fpd0_60 = buffer.data(fpd0 + 60);
    const auto *fpd0_75 = buffer.data(fpd0 + 75);
    const auto *fpd0_77 = buffer.data(fpd0 + 77);
    const auto *fpd0_78 = buffer.data(fpd0 + 78);
    const auto *fpd0_83 = buffer.data(fpd0 + 83);
    const auto *fpd0_87 = buffer.data(fpd0 + 87);
    const auto *fpd0_90 = buffer.data(fpd0 + 90);
    const auto *fpd0_93 = buffer.data(fpd0 + 93);
    const auto *fpd0_94 = buffer.data(fpd0 + 94);
    const auto *fpd0_95 = buffer.data(fpd0 + 95);

    const auto *fpd1_51 = buffer.data(fpd1 + 51);
    const auto *fpd1_52 = buffer.data(fpd1 + 52);
    const auto *fpd1_53 = buffer.data(fpd1 + 53);
    const auto *fpd1_54 = buffer.data(fpd1 + 54);
    const auto *fpd1_57 = buffer.data(fpd1 + 57);
    const auto *fpd1_59 = buffer.data(fpd1 + 59);
    const auto *fpd1_60 = buffer.data(fpd1 + 60);
    const auto *fpd1_75 = buffer.data(fpd1 + 75);
    const auto *fpd1_77 = buffer.data(fpd1 + 77);
    const auto *fpd1_78 = buffer.data(fpd1 + 78);
    const auto *fpd1_83 = buffer.data(fpd1 + 83);
    const auto *fpd1_87 = buffer.data(fpd1 + 87);
    const auto *fpd1_90 = buffer.data(fpd1 + 90);
    const auto *fpd1_93 = buffer.data(fpd1 + 93);
    const auto *fpd1_94 = buffer.data(fpd1 + 94);
    const auto *fpd1_95 = buffer.data(fpd1 + 95);

    const auto *fpf_86 = buffer.data(fpf + 86);
    const auto *fpf_87 = buffer.data(fpf + 87);
    const auto *fpf_88 = buffer.data(fpf + 88);
    const auto *fpf_89 = buffer.data(fpf + 89);
    const auto *fpf_90 = buffer.data(fpf + 90);
    const auto *fpf_91 = buffer.data(fpf + 91);
    const auto *fpf_92 = buffer.data(fpf + 92);
    const auto *fpf_93 = buffer.data(fpf + 93);
    const auto *fpf_96 = buffer.data(fpf + 96);
    const auto *fpf_97 = buffer.data(fpf + 97);
    const auto *fpf_98 = buffer.data(fpf + 98);
    const auto *fpf_99 = buffer.data(fpf + 99);
    const auto *fpf_100 = buffer.data(fpf + 100);
    const auto *fpf_101 = buffer.data(fpf + 101);
    const auto *fpf_102 = buffer.data(fpf + 102);
    const auto *fpf_103 = buffer.data(fpf + 103);
    const auto *fpf_106 = buffer.data(fpf + 106);
    const auto *fpf_108 = buffer.data(fpf + 108);
    const auto *fpf_109 = buffer.data(fpf + 109);
    const auto *fpf_110 = buffer.data(fpf + 110);
    const auto *fpf_111 = buffer.data(fpf + 111);
    const auto *fpf_113 = buffer.data(fpf + 113);
    const auto *fpf_116 = buffer.data(fpf + 116);
    const auto *fpf_118 = buffer.data(fpf + 118);
    const auto *fpf_119 = buffer.data(fpf + 119);
    const auto *fpf_120 = buffer.data(fpf + 120);
    const auto *fpf_122 = buffer.data(fpf + 122);
    const auto *fpf_126 = buffer.data(fpf + 126);
    const auto *fpf_127 = buffer.data(fpf + 127);
    const auto *fpf_128 = buffer.data(fpf + 128);
    const auto *fpf_129 = buffer.data(fpf + 129);
    const auto *fpf_130 = buffer.data(fpf + 130);
    const auto *fpf_132 = buffer.data(fpf + 132);
    const auto *fpf_135 = buffer.data(fpf + 135);
    const auto *fpf_136 = buffer.data(fpf + 136);
    const auto *fpf_137 = buffer.data(fpf + 137);
    const auto *fpf_138 = buffer.data(fpf + 138);
    const auto *fpf_139 = buffer.data(fpf + 139);
    const auto *fpf_140 = buffer.data(fpf + 140);
    const auto *fpf_142 = buffer.data(fpf + 142);
    const auto *fpf_143 = buffer.data(fpf + 143);
    const auto *fpf_146 = buffer.data(fpf + 146);
    const auto *fpf_147 = buffer.data(fpf + 147);
    const auto *fpf_148 = buffer.data(fpf + 148);
    const auto *fpf_149 = buffer.data(fpf + 149);
    const auto *fpf_150 = buffer.data(fpf + 150);
    const auto *fpf_151 = buffer.data(fpf + 151);
    const auto *fpf_152 = buffer.data(fpf + 152);
    const auto *fpf_155 = buffer.data(fpf + 155);
    const auto *fpf_156 = buffer.data(fpf + 156);
    const auto *fpf_157 = buffer.data(fpf + 157);
    const auto *fpf_158 = buffer.data(fpf + 158);
    const auto *fpf_159 = buffer.data(fpf + 159);
    const auto *fpf_160 = buffer.data(fpf + 160);
    const auto *fpf_162 = buffer.data(fpf + 162);
    const auto *fpf_165 = buffer.data(fpf + 165);
    const auto *fpf_166 = buffer.data(fpf + 166);
    const auto *fpf_167 = buffer.data(fpf + 167);
    const auto *fpf_169 = buffer.data(fpf + 169);

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, dpf_89, fpd0_51, fpd0_52, fpd1_51, \
                         fpd1_52, fpf_86, fpf_87, fpf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_8 * dpf_89[k]
                   + f_4 * pc_x[k] * fpf_89[k];

        t_130[k] = f_2 * fpd0_51[k]
                   - f_3 * fpd1_51[k]
                   + f_4 * pc_y[k] * fpf_86[k];

        t_131[k] = f_12 * fpd0_52[k]
                   - f_13 * fpd1_52[k]
                   + f_4 * pc_y[k] * fpf_87[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pa_x, pc_x, pc_y, ppg0_134, ppg1_134, dpg0_134, \
                         dpg1_134, fpd0_53, fpd1_53, fpf_88, fpf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_5 * fpd0_53[k]
                   - f_6 * fpd1_53[k]
                   + f_4 * pc_y[k] * fpf_88[k];

        t_133[k] = f_4 * pc_y[k] * fpf_89[k];

        t_134[k] = f_10 * ppg0_134[k]
                   - f_11 * ppg1_134[k]
                   + pa_x[k] * dpg0_134[k]
                   - f_7 * pc_x[k] * dpg1_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pc_x, pc_y, pc_z, dpf_30, dpf_90, fsf_30, \
                         fpd0_54, fpd1_54, fpf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_1 * dpf_90[k]
                   + f_1 * fsf_30[k]
                   + f_2 * fpd0_54[k]
                   - f_3 * fpd1_54[k]
                   + f_4 * pc_x[k] * fpf_90[k];

        t_136[k] = f_8 * dpf_30[k]
                   + f_4 * pc_y[k] * fpf_90[k];

        t_137[k] = f_4 * pc_z[k] * fpf_90[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_z, dpf_93, fsf_33, fpd0_54, fpd0_57, \
                         fpd1_54, fpd1_57, fpf_91, fpf_92, fpf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_1 * dpf_93[k]
                   + f_1 * fsf_33[k]
                   + f_5 * fpd0_57[k]
                   - f_6 * fpd1_57[k]
                   + f_4 * pc_x[k] * fpf_93[k];

        t_139[k] = f_4 * pc_z[k] * fpf_91[k];

        t_140[k] = f_5 * fpd0_54[k]
                   - f_6 * fpd1_54[k]
                   + f_4 * pc_z[k] * fpf_92[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pc_x, pc_z, dpf_96, dpf_98, dpf_99, \
                         fsf_36, fsf_38, fsf_39, fpf_93, fpf_96, fpf_98, \
                         fpf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_1 * dpf_96[k]
                   + f_1 * fsf_36[k]
                   + f_4 * pc_x[k] * fpf_96[k];

        t_142[k] = f_4 * pc_z[k] * fpf_93[k];

        t_143[k] = f_1 * dpf_98[k]
                   + f_1 * fsf_38[k]
                   + f_4 * pc_x[k] * fpf_98[k];

        t_144[k] = f_1 * dpf_99[k]
                   + f_1 * fsf_39[k]
                   + f_4 * pc_x[k] * fpf_99[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pc_y, pc_z, dpf_36, dpf_39, \
                         fpd0_57, fpd0_59, fpd1_57, fpd1_59, fpf_96, fpf_97, \
                         fpf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_8 * dpf_36[k]
                   + f_2 * fpd0_57[k]
                   - f_3 * fpd1_57[k]
                   + f_4 * pc_y[k] * fpf_96[k];

        t_146[k] = f_4 * pc_z[k] * fpf_96[k];

        t_147[k] = f_5 * fpd0_57[k]
                   - f_6 * fpd1_57[k]
                   + f_4 * pc_z[k] * fpf_97[k];

        t_148[k] = f_8 * dpf_39[k]
                   + f_4 * pc_y[k] * fpf_99[k];

        t_149[k] = f_2 * fpd0_59[k]
                   - f_3 * fpd1_59[k]
                   + f_4 * pc_z[k] * fpf_99[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_x, pc_x, pc_y, pc_z, dpg0_150, dpf_40, \
                         dpf_100, dpg1_150, fsf_30, fpf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_x[k] * dpg0_150[k]
                   + f_9 * dpf_100[k]
                   - f_7 * pc_x[k] * dpg1_150[k];

        t_151[k] = f_8 * dpf_40[k]
                   + f_1 * fsf_30[k]
                   + f_4 * pc_y[k] * fpf_100[k];

        t_152[k] = f_4 * pc_z[k] * fpf_100[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_x, pc_x, pc_z, dpg0_153, dpf_103, \
                         dpf_106, dpg1_153, fpd0_60, fpd1_60, fpf_101, fpf_102, \
                         fpf_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_x[k] * dpg0_153[k]
                   + f_8 * dpf_103[k]
                   - f_7 * pc_x[k] * dpg1_153[k];

        t_154[k] = f_4 * pc_z[k] * fpf_101[k];

        t_155[k] = f_5 * fpd0_60[k]
                   - f_6 * fpd1_60[k]
                   + f_4 * pc_z[k] * fpf_102[k];

        t_156[k] = f_1 * dpf_106[k]
                   + f_4 * pc_x[k] * fpf_106[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pa_x, pc_x, pc_z, dpg0_160, \
                         dpf_108, dpf_109, dpg1_160, fpf_103, fpf_106, fpf_108, \
                         fpf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_4 * pc_z[k] * fpf_103[k];

        t_158[k] = f_1 * dpf_108[k]
                   + f_4 * pc_x[k] * fpf_108[k];

        t_159[k] = f_1 * dpf_109[k]
                   + f_4 * pc_x[k] * fpf_109[k];

        t_160[k] = pa_x[k] * dpg0_160[k]
                   - f_7 * pc_x[k] * dpg1_160[k];

        t_161[k] = f_4 * pc_z[k] * fpf_106[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_x, pb_z, pc_x, pc_z, dpg0_162, \
                         dpg0_163, dpg0_164, dpg1_162, dpg1_163, dpg1_164, fsg0_45, \
                         fsg1_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_x[k] * dpg0_162[k]
                   - f_7 * pc_x[k] * dpg1_162[k];

        t_163[k] = pa_x[k] * dpg0_163[k]
                   - f_7 * pc_x[k] * dpg1_163[k];

        t_164[k] = pa_x[k] * dpg0_164[k]
                   - f_7 * pc_x[k] * dpg1_164[k];

        t_165[k] = pb_z[k] * fsg0_45[k]
                   - f_7 * pc_z[k] * fsg1_45[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, dpf_50, fsg0_48, \
                         fsf_30, fsf_31, fsg1_48, fpf_110, fpf_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_8 * dpf_50[k]
                   + f_4 * pc_y[k] * fpf_110[k];

        t_167[k] = f_1 * fsf_30[k]
                   + f_4 * pc_z[k] * fpf_110[k];

        t_168[k] = pb_z[k] * fsg0_48[k]
                   - f_7 * pc_z[k] * fsg1_48[k];

        t_169[k] = f_1 * fsf_31[k]
                   + f_4 * pc_z[k] * fpf_111[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_x, pc_x, pc_z, dpg0_170, dpf_115, \
                         dpf_116, dpf_118, dpg1_170, fsf_33, fpf_113, fpf_116, \
                         fpf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pa_x[k] * dpg0_170[k]
                   + f_8 * dpf_115[k]
                   - f_7 * pc_x[k] * dpg1_170[k];

        t_171[k] = f_1 * dpf_116[k]
                   + f_4 * pc_x[k] * fpf_116[k];

        t_172[k] = f_1 * fsf_33[k]
                   + f_4 * pc_z[k] * fpf_113[k];

        t_173[k] = f_1 * dpf_118[k]
                   + f_4 * pc_x[k] * fpf_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_x, pc_x, pc_z, dpg0_175, dpg0_177, \
                         dpf_119, dpg1_175, dpg1_177, fsf_36, fpf_116, \
                         fpf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_1 * dpf_119[k]
                   + f_4 * pc_x[k] * fpf_119[k];

        t_175[k] = pa_x[k] * dpg0_175[k]
                   - f_7 * pc_x[k] * dpg1_175[k];

        t_176[k] = f_1 * fsf_36[k]
                   + f_4 * pc_z[k] * fpf_116[k];

        t_177[k] = pa_x[k] * dpg0_177[k]
                   - f_7 * pc_x[k] * dpg1_177[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pa_x, pa_y, pc_x, pc_y, dpg0_90, \
                         dpg0_179, dpf_59, dpf_60, dpg1_90, dpg1_179, fpf_119, \
                         fpf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_8 * dpf_59[k]
                   + f_4 * pc_y[k] * fpf_119[k];

        t_179[k] = pa_x[k] * dpg0_179[k]
                   - f_7 * pc_x[k] * dpg1_179[k];

        t_180[k] = pa_y[k] * dpg0_90[k]
                   - f_7 * pc_y[k] * dpg1_90[k];

        t_181[k] = f_1 * dpf_60[k]
                   + f_4 * pc_y[k] * fpf_120[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pa_y, pa_z, pc_y, pc_z, dpg0_48, dpg0_95, \
                         dpf_30, dpf_62, dpg1_48, dpg1_95, fpf_120, \
                         fpf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_1 * dpf_30[k]
                   + f_4 * pc_z[k] * fpf_120[k];

        t_183[k] = pa_z[k] * dpg0_48[k]
                   - f_7 * pc_z[k] * dpg1_48[k];

        t_184[k] = f_1 * dpf_62[k]
                   + f_4 * pc_y[k] * fpf_122[k];

        t_185[k] = pa_y[k] * dpg0_95[k]
                   - f_7 * pc_y[k] * dpg1_95[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_z, pc_x, pc_z, dpg0_51, dpf_127, dpf_128, \
                         dpg1_51, fsf_47, fsf_48, fpf_127, fpf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_z[k] * dpg0_51[k]
                   - f_7 * pc_z[k] * dpg1_51[k];

        t_187[k] = f_1 * dpf_127[k]
                   + f_1 * fsf_47[k]
                   + f_4 * pc_x[k] * fpf_127[k];

        t_188[k] = f_1 * dpf_128[k]
                   + f_1 * fsf_48[k]
                   + f_4 * pc_x[k] * fpf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_y, pc_y, pc_z, dpg0_99, dpf_36, dpf_66, \
                         dpg1_99, fpd0_75, fpd1_75, fpf_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pa_y[k] * dpg0_99[k]
                   - f_7 * pc_y[k] * dpg1_99[k];

        t_190[k] = f_1 * dpf_66[k]
                   + f_2 * fpd0_75[k]
                   - f_3 * fpd1_75[k]
                   + f_4 * pc_y[k] * fpf_126[k];

        t_191[k] = f_1 * dpf_36[k]
                   + f_4 * pc_z[k] * fpf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pc_y, pc_z, dpf_39, dpf_68, dpf_69, fpd0_77, \
                         fpd1_77, fpf_128, fpf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_1 * dpf_68[k]
                   + f_5 * fpd0_77[k]
                   - f_6 * fpd1_77[k]
                   + f_4 * pc_y[k] * fpf_128[k];

        t_193[k] = f_1 * dpf_69[k]
                   + f_4 * pc_y[k] * fpf_129[k];

        t_194[k] = f_1 * dpf_39[k]
                   + f_2 * fpd0_77[k]
                   - f_3 * fpd1_77[k]
                   + f_4 * pc_z[k] * fpf_129[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_z, pc_x, pc_z, dpg0_61, dpg0_63, \
                         dpf_40, dpf_130, dpg1_61, dpg1_63, fpd0_78, fpd1_78, \
                         fpf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_1 * dpf_130[k]
                   + f_2 * fpd0_78[k]
                   - f_3 * fpd1_78[k]
                   + f_4 * pc_x[k] * fpf_130[k];

        t_196[k] = pa_z[k] * dpg0_61[k]
                   - f_7 * pc_z[k] * dpg1_61[k];

        t_197[k] = f_1 * dpf_40[k]
                   + f_4 * pc_z[k] * fpf_130[k];

        t_198[k] = pa_z[k] * dpg0_63[k]
                   - f_7 * pc_z[k] * dpg1_63[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, pc_x, pc_y, dpf_72, dpf_135, dpf_136, fsf_42, \
                         fpd0_83, fpd1_83, fpf_132, fpf_135, fpf_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_1 * dpf_72[k]
                   + f_1 * fsf_42[k]
                   + f_4 * pc_y[k] * fpf_132[k];

        t_200[k] = f_1 * dpf_135[k]
                   + f_5 * fpd0_83[k]
                   - f_6 * fpd1_83[k]
                   + f_4 * pc_x[k] * fpf_135[k];

        t_201[k] = f_1 * dpf_136[k]
                   + f_4 * pc_x[k] * fpf_136[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pa_x, pc_x, dpg0_205, dpf_137, dpf_138, \
                         dpf_139, dpg1_205, fpf_137, fpf_138, fpf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_1 * dpf_137[k]
                   + f_4 * pc_x[k] * fpf_137[k];

        t_203[k] = f_1 * dpf_138[k]
                   + f_4 * pc_x[k] * fpf_138[k];

        t_204[k] = f_1 * dpf_139[k]
                   + f_4 * pc_x[k] * fpf_139[k];

        t_205[k] = pa_x[k] * dpg0_205[k]
                   - f_7 * pc_x[k] * dpg1_205[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_x, pc_x, pc_z, dpg0_207, dpg0_208, \
                         dpg0_209, dpf_46, dpg1_207, dpg1_208, dpg1_209, \
                         fpf_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_1 * dpf_46[k]
                   + f_4 * pc_z[k] * fpf_136[k];

        t_207[k] = pa_x[k] * dpg0_207[k]
                   - f_7 * pc_x[k] * dpg1_207[k];

        t_208[k] = pa_x[k] * dpg0_208[k]
                   - f_7 * pc_x[k] * dpg1_208[k];

        t_209[k] = pa_x[k] * dpg0_209[k]
                   - f_7 * pc_x[k] * dpg1_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pa_y, pc_y, dpg0_120, dpg0_122, dpf_80, \
                         dpg1_120, dpg1_122, fpf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = pa_y[k] * dpg0_120[k]
                   - f_7 * pc_y[k] * dpg1_120[k];

        t_211[k] = f_1 * dpf_80[k]
                   + f_4 * pc_y[k] * fpf_140[k];

        t_212[k] = pa_y[k] * dpg0_122[k]
                   - f_7 * pc_y[k] * dpg1_122[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, pa_y, pc_x, pc_y, dpg0_125, dpf_82, dpf_143, \
                         dpg1_125, fpd0_87, fpd1_87, fpf_142, fpf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_1 * dpf_143[k]
                   + f_5 * fpd0_87[k]
                   - f_6 * fpd1_87[k]
                   + f_4 * pc_x[k] * fpf_143[k];

        t_214[k] = f_1 * dpf_82[k]
                   + f_4 * pc_y[k] * fpf_142[k];

        t_215[k] = pa_y[k] * dpg0_125[k]
                   - f_7 * pc_y[k] * dpg1_125[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, dpf_146, dpf_147, dpf_148, dpf_149, \
                         fpf_146, fpf_147, fpf_148, fpf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_1 * dpf_146[k]
                   + f_4 * pc_x[k] * fpf_146[k];

        t_217[k] = f_1 * dpf_147[k]
                   + f_4 * pc_x[k] * fpf_147[k];

        t_218[k] = f_1 * dpf_148[k]
                   + f_4 * pc_x[k] * fpf_148[k];

        t_219[k] = f_1 * dpf_149[k]
                   + f_4 * pc_x[k] * fpf_149[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pa_x, pc_x, pc_y, dpg0_220, dpg0_221, \
                         dpg0_222, dpf_89, dpg1_220, dpg1_221, dpg1_222, \
                         fpf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pa_x[k] * dpg0_220[k]
                   - f_7 * pc_x[k] * dpg1_220[k];

        t_221[k] = pa_x[k] * dpg0_221[k]
                   - f_7 * pc_x[k] * dpg1_221[k];

        t_222[k] = pa_x[k] * dpg0_222[k]
                   - f_7 * pc_x[k] * dpg1_222[k];

        t_223[k] = f_1 * dpf_89[k]
                   + f_4 * pc_y[k] * fpf_149[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_x, pc_x, pc_y, pc_z, dpg0_224, dpf_60, \
                         dpf_150, dpg1_224, fsf_50, fpd0_90, fpd1_90, \
                         fpf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_x[k] * dpg0_224[k]
                   - f_7 * pc_x[k] * dpg1_224[k];

        t_225[k] = f_1 * dpf_150[k]
                   + f_1 * fsf_50[k]
                   + f_2 * fpd0_90[k]
                   - f_3 * fpd1_90[k]
                   + f_4 * pc_x[k] * fpf_150[k];

        t_226[k] = f_4 * pc_y[k] * fpf_150[k];

        t_227[k] = f_8 * dpf_60[k]
                   + f_4 * pc_z[k] * fpf_150[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pc_x, pc_y, dpf_155, fsf_55, fpd0_90, fpd0_95, \
                         fpd1_90, fpd1_95, fpf_151, fpf_152, fpf_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_5 * fpd0_90[k]
                   - f_6 * fpd1_90[k]
                   + f_4 * pc_y[k] * fpf_151[k];

        t_229[k] = f_4 * pc_y[k] * fpf_152[k];

        t_230[k] = f_1 * dpf_155[k]
                   + f_1 * fsf_55[k]
                   + f_5 * fpd0_95[k]
                   - f_6 * fpd1_95[k]
                   + f_4 * pc_x[k] * fpf_155[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pc_x, pc_y, dpf_156, dpf_157, dpf_159, \
                         fsf_56, fsf_57, fsf_59, fpf_155, fpf_156, fpf_157, \
                         fpf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_1 * dpf_156[k]
                   + f_1 * fsf_56[k]
                   + f_4 * pc_x[k] * fpf_156[k];

        t_232[k] = f_1 * dpf_157[k]
                   + f_1 * fsf_57[k]
                   + f_4 * pc_x[k] * fpf_157[k];

        t_233[k] = f_4 * pc_y[k] * fpf_155[k];

        t_234[k] = f_1 * dpf_159[k]
                   + f_1 * fsf_59[k]
                   + f_4 * pc_x[k] * fpf_159[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_y, fpd0_93, fpd0_94, fpd0_95, fpd1_93, \
                         fpd1_94, fpd1_95, fpf_156, fpf_157, fpf_158, \
                         fpf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_2 * fpd0_93[k]
                   - f_3 * fpd1_93[k]
                   + f_4 * pc_y[k] * fpf_156[k];

        t_236[k] = f_12 * fpd0_94[k]
                   - f_13 * fpd1_94[k]
                   + f_4 * pc_y[k] * fpf_157[k];

        t_237[k] = f_5 * fpd0_95[k]
                   - f_6 * fpd1_95[k]
                   + f_4 * pc_y[k] * fpf_158[k];

        t_238[k] = f_4 * pc_y[k] * fpf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pb_y, pc_y, pc_z, dpf_69, dpf_70, \
                         fsg0_75, fsf_50, fsg1_75, fpd0_95, fpd1_95, fpf_159, \
                         fpf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_8 * dpf_69[k]
                   + f_2 * fpd0_95[k]
                   - f_3 * fpd1_95[k]
                   + f_4 * pc_z[k] * fpf_159[k];

        t_240[k] = pb_y[k] * fsg0_75[k]
                   - f_7 * pc_y[k] * fsg1_75[k];

        t_241[k] = f_1 * fsf_50[k]
                   + f_4 * pc_y[k] * fpf_160[k];

        t_242[k] = f_8 * dpf_70[k]
                   + f_4 * pc_z[k] * fpf_160[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_x, pb_y, pc_x, pc_y, dpg0_243, dpf_163, \
                         dpg1_243, fsg0_80, fsf_52, fsg1_80, fpf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pa_x[k] * dpg0_243[k]
                   + f_8 * dpf_163[k]
                   - f_7 * pc_x[k] * dpg1_243[k];

        t_244[k] = f_1 * fsf_52[k]
                   + f_4 * pc_y[k] * fpf_162[k];

        t_245[k] = pb_y[k] * fsg0_80[k]
                   - f_7 * pc_y[k] * fsg1_80[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, pc_y, dpf_166, dpf_167, dpf_169, \
                         fsf_55, fpf_165, fpf_166, fpf_167, fpf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_1 * dpf_166[k]
                   + f_4 * pc_x[k] * fpf_166[k];

        t_247[k] = f_1 * dpf_167[k]
                   + f_4 * pc_x[k] * fpf_167[k];

        t_248[k] = f_1 * fsf_55[k]
                   + f_4 * pc_y[k] * fpf_165[k];

        t_249[k] = f_1 * dpf_169[k]
                   + f_4 * pc_x[k] * fpf_169[k];
    }
}

static auto
compute_prim_fpg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppg0, const size_t ppg1,
                                                          const size_t dpg0, const size_t dpf,
                                                          const size_t dpg1, const size_t fsg0,
                                                          const size_t fsf, const size_t fsg1,
                                                          const size_t fpd0, const size_t fpd1,
                                                          const size_t fpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);

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
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppg0_134 = buffer.data(ppg0 + 134);

    const auto *ppg1_134 = buffer.data(ppg1 + 134);

    const auto *dpg0_135 = buffer.data(dpg0 + 135);
    const auto *dpg0_136 = buffer.data(dpg0 + 136);
    const auto *dpg0_138 = buffer.data(dpg0 + 138);
    const auto *dpg0_139 = buffer.data(dpg0 + 139);
    const auto *dpg0_145 = buffer.data(dpg0 + 145);
    const auto *dpg0_150 = buffer.data(dpg0 + 150);
    const auto *dpg0_151 = buffer.data(dpg0 + 151);
    const auto *dpg0_153 = buffer.data(dpg0 + 153);
    const auto *dpg0_160 = buffer.data(dpg0 + 160);
    const auto *dpg0_162 = buffer.data(dpg0 + 162);
    const auto *dpg0_224 = buffer.data(dpg0 + 224);
    const auto *dpg0_225 = buffer.data(dpg0 + 225);
    const auto *dpg0_226 = buffer.data(dpg0 + 226);
    const auto *dpg0_227 = buffer.data(dpg0 + 227);
    const auto *dpg0_228 = buffer.data(dpg0 + 228);
    const auto *dpg0_229 = buffer.data(dpg0 + 229);
    const auto *dpg0_230 = buffer.data(dpg0 + 230);
    const auto *dpg0_239 = buffer.data(dpg0 + 239);
    const auto *dpg0_250 = buffer.data(dpg0 + 250);
    const auto *dpg0_251 = buffer.data(dpg0 + 251);
    const auto *dpg0_252 = buffer.data(dpg0 + 252);
    const auto *dpg0_254 = buffer.data(dpg0 + 254);
    const auto *dpg0_255 = buffer.data(dpg0 + 255);
    const auto *dpg0_260 = buffer.data(dpg0 + 260);
    const auto *dpg0_265 = buffer.data(dpg0 + 265);
    const auto *dpg0_266 = buffer.data(dpg0 + 266);
    const auto *dpg0_267 = buffer.data(dpg0 + 267);
    const auto *dpg0_269 = buffer.data(dpg0 + 269);

    const auto *dpf_80 = buffer.data(dpf + 80);
    const auto *dpf_91 = buffer.data(dpf + 91);
    const auto *dpf_96 = buffer.data(dpf + 96);
    const auto *dpf_99 = buffer.data(dpf + 99);
    const auto *dpf_106 = buffer.data(dpf + 106);
    const auto *dpf_107 = buffer.data(dpf + 107);
    const auto *dpf_109 = buffer.data(dpf + 109);
    const auto *dpf_116 = buffer.data(dpf + 116);
    const auto *dpf_119 = buffer.data(dpf + 119);
    const auto *dpf_126 = buffer.data(dpf + 126);
    const auto *dpf_129 = buffer.data(dpf + 129);
    const auto *dpf_139 = buffer.data(dpf + 139);
    const auto *dpf_146 = buffer.data(dpf + 146);
    const auto *dpf_148 = buffer.data(dpf + 148);
    const auto *dpf_149 = buffer.data(dpf + 149);
    const auto *dpf_150 = buffer.data(dpf + 150);
    const auto *dpf_151 = buffer.data(dpf + 151);
    const auto *dpf_152 = buffer.data(dpf + 152);
    const auto *dpf_159 = buffer.data(dpf + 159);
    const auto *dpf_170 = buffer.data(dpf + 170);
    const auto *dpf_175 = buffer.data(dpf + 175);
    const auto *dpf_176 = buffer.data(dpf + 176);
    const auto *dpf_177 = buffer.data(dpf + 177);
    const auto *dpf_179 = buffer.data(dpf + 179);

    const auto *dpg1_135 = buffer.data(dpg1 + 135);
    const auto *dpg1_136 = buffer.data(dpg1 + 136);
    const auto *dpg1_138 = buffer.data(dpg1 + 138);
    const auto *dpg1_139 = buffer.data(dpg1 + 139);
    const auto *dpg1_145 = buffer.data(dpg1 + 145);
    const auto *dpg1_150 = buffer.data(dpg1 + 150);
    const auto *dpg1_151 = buffer.data(dpg1 + 151);
    const auto *dpg1_153 = buffer.data(dpg1 + 153);
    const auto *dpg1_160 = buffer.data(dpg1 + 160);
    const auto *dpg1_162 = buffer.data(dpg1 + 162);
    const auto *dpg1_224 = buffer.data(dpg1 + 224);
    const auto *dpg1_225 = buffer.data(dpg1 + 225);
    const auto *dpg1_226 = buffer.data(dpg1 + 226);
    const auto *dpg1_227 = buffer.data(dpg1 + 227);
    const auto *dpg1_228 = buffer.data(dpg1 + 228);
    const auto *dpg1_229 = buffer.data(dpg1 + 229);
    const auto *dpg1_230 = buffer.data(dpg1 + 230);
    const auto *dpg1_239 = buffer.data(dpg1 + 239);
    const auto *dpg1_250 = buffer.data(dpg1 + 250);
    const auto *dpg1_251 = buffer.data(dpg1 + 251);
    const auto *dpg1_252 = buffer.data(dpg1 + 252);
    const auto *dpg1_254 = buffer.data(dpg1 + 254);
    const auto *dpg1_255 = buffer.data(dpg1 + 255);
    const auto *dpg1_260 = buffer.data(dpg1 + 260);
    const auto *dpg1_265 = buffer.data(dpg1 + 265);
    const auto *dpg1_266 = buffer.data(dpg1 + 266);
    const auto *dpg1_267 = buffer.data(dpg1 + 267);
    const auto *dpg1_269 = buffer.data(dpg1 + 269);

    const auto *fsg0_90 = buffer.data(fsg0 + 90);
    const auto *fsg0_91 = buffer.data(fsg0 + 91);
    const auto *fsg0_93 = buffer.data(fsg0 + 93);
    const auto *fsg0_95 = buffer.data(fsg0 + 95);
    const auto *fsg0_100 = buffer.data(fsg0 + 100);
    const auto *fsg0_102 = buffer.data(fsg0 + 102);
    const auto *fsg0_104 = buffer.data(fsg0 + 104);
    const auto *fsg0_107 = buffer.data(fsg0 + 107);
    const auto *fsg0_110 = buffer.data(fsg0 + 110);
    const auto *fsg0_117 = buffer.data(fsg0 + 117);
    const auto *fsg0_119 = buffer.data(fsg0 + 119);
    const auto *fsg0_130 = buffer.data(fsg0 + 130);
    const auto *fsg0_132 = buffer.data(fsg0 + 132);

    const auto *fsf_50 = buffer.data(fsf + 50);
    const auto *fsf_59 = buffer.data(fsf + 59);
    const auto *fsf_60 = buffer.data(fsf + 60);
    const auto *fsf_61 = buffer.data(fsf + 61);
    const auto *fsf_63 = buffer.data(fsf + 63);
    const auto *fsf_65 = buffer.data(fsf + 65);
    const auto *fsf_66 = buffer.data(fsf + 66);
    const auto *fsf_67 = buffer.data(fsf + 67);
    const auto *fsf_68 = buffer.data(fsf + 68);
    const auto *fsf_69 = buffer.data(fsf + 69);
    const auto *fsf_72 = buffer.data(fsf + 72);
    const auto *fsf_75 = buffer.data(fsf + 75);
    const auto *fsf_76 = buffer.data(fsf + 76);
    const auto *fsf_77 = buffer.data(fsf + 77);
    const auto *fsf_78 = buffer.data(fsf + 78);
    const auto *fsf_79 = buffer.data(fsf + 79);
    const auto *fsf_86 = buffer.data(fsf + 86);
    const auto *fsf_87 = buffer.data(fsf + 87);
    const auto *fsf_88 = buffer.data(fsf + 88);
    const auto *fsf_89 = buffer.data(fsf + 89);

    const auto *fsg1_90 = buffer.data(fsg1 + 90);
    const auto *fsg1_91 = buffer.data(fsg1 + 91);
    const auto *fsg1_93 = buffer.data(fsg1 + 93);
    const auto *fsg1_95 = buffer.data(fsg1 + 95);
    const auto *fsg1_100 = buffer.data(fsg1 + 100);
    const auto *fsg1_102 = buffer.data(fsg1 + 102);
    const auto *fsg1_104 = buffer.data(fsg1 + 104);
    const auto *fsg1_107 = buffer.data(fsg1 + 107);
    const auto *fsg1_110 = buffer.data(fsg1 + 110);
    const auto *fsg1_117 = buffer.data(fsg1 + 117);
    const auto *fsg1_119 = buffer.data(fsg1 + 119);
    const auto *fsg1_130 = buffer.data(fsg1 + 130);
    const auto *fsg1_132 = buffer.data(fsg1 + 132);

    const auto *fpd0_102 = buffer.data(fpd0 + 102);
    const auto *fpd0_114 = buffer.data(fpd0 + 114);
    const auto *fpd0_115 = buffer.data(fpd0 + 115);
    const auto *fpd0_117 = buffer.data(fpd0 + 117);
    const auto *fpd0_119 = buffer.data(fpd0 + 119);
    const auto *fpd0_125 = buffer.data(fpd0 + 125);
    const auto *fpd0_134 = buffer.data(fpd0 + 134);
    const auto *fpd0_136 = buffer.data(fpd0 + 136);
    const auto *fpd0_137 = buffer.data(fpd0 + 137);
    const auto *fpd0_138 = buffer.data(fpd0 + 138);
    const auto *fpd0_139 = buffer.data(fpd0 + 139);
    const auto *fpd0_140 = buffer.data(fpd0 + 140);
    const auto *fpd0_141 = buffer.data(fpd0 + 141);
    const auto *fpd0_142 = buffer.data(fpd0 + 142);
    const auto *fpd0_143 = buffer.data(fpd0 + 143);

    const auto *fpd1_102 = buffer.data(fpd1 + 102);
    const auto *fpd1_114 = buffer.data(fpd1 + 114);
    const auto *fpd1_115 = buffer.data(fpd1 + 115);
    const auto *fpd1_117 = buffer.data(fpd1 + 117);
    const auto *fpd1_119 = buffer.data(fpd1 + 119);
    const auto *fpd1_125 = buffer.data(fpd1 + 125);
    const auto *fpd1_134 = buffer.data(fpd1 + 134);
    const auto *fpd1_136 = buffer.data(fpd1 + 136);
    const auto *fpd1_137 = buffer.data(fpd1 + 137);
    const auto *fpd1_138 = buffer.data(fpd1 + 138);
    const auto *fpd1_139 = buffer.data(fpd1 + 139);
    const auto *fpd1_140 = buffer.data(fpd1 + 140);
    const auto *fpd1_141 = buffer.data(fpd1 + 141);
    const auto *fpd1_142 = buffer.data(fpd1 + 142);
    const auto *fpd1_143 = buffer.data(fpd1 + 143);

    const auto *fpf_169 = buffer.data(fpf + 169);
    const auto *fpf_170 = buffer.data(fpf + 170);
    const auto *fpf_171 = buffer.data(fpf + 171);
    const auto *fpf_172 = buffer.data(fpf + 172);
    const auto *fpf_175 = buffer.data(fpf + 175);
    const auto *fpf_176 = buffer.data(fpf + 176);
    const auto *fpf_177 = buffer.data(fpf + 177);
    const auto *fpf_179 = buffer.data(fpf + 179);
    const auto *fpf_180 = buffer.data(fpf + 180);
    const auto *fpf_181 = buffer.data(fpf + 181);
    const auto *fpf_186 = buffer.data(fpf + 186);
    const auto *fpf_187 = buffer.data(fpf + 187);
    const auto *fpf_188 = buffer.data(fpf + 188);
    const auto *fpf_189 = buffer.data(fpf + 189);
    const auto *fpf_190 = buffer.data(fpf + 190);
    const auto *fpf_191 = buffer.data(fpf + 191);
    const auto *fpf_193 = buffer.data(fpf + 193);
    const auto *fpf_195 = buffer.data(fpf + 195);
    const auto *fpf_196 = buffer.data(fpf + 196);
    const auto *fpf_197 = buffer.data(fpf + 197);
    const auto *fpf_198 = buffer.data(fpf + 198);
    const auto *fpf_199 = buffer.data(fpf + 199);
    const auto *fpf_200 = buffer.data(fpf + 200);
    const auto *fpf_201 = buffer.data(fpf + 201);
    const auto *fpf_205 = buffer.data(fpf + 205);
    const auto *fpf_206 = buffer.data(fpf + 206);
    const auto *fpf_207 = buffer.data(fpf + 207);
    const auto *fpf_208 = buffer.data(fpf + 208);
    const auto *fpf_209 = buffer.data(fpf + 209);
    const auto *fpf_216 = buffer.data(fpf + 216);
    const auto *fpf_217 = buffer.data(fpf + 217);
    const auto *fpf_218 = buffer.data(fpf + 218);
    const auto *fpf_219 = buffer.data(fpf + 219);
    const auto *fpf_222 = buffer.data(fpf + 222);
    const auto *fpf_224 = buffer.data(fpf + 224);
    const auto *fpf_225 = buffer.data(fpf + 225);
    const auto *fpf_226 = buffer.data(fpf + 226);
    const auto *fpf_227 = buffer.data(fpf + 227);
    const auto *fpf_228 = buffer.data(fpf + 228);
    const auto *fpf_229 = buffer.data(fpf + 229);
    const auto *fpf_230 = buffer.data(fpf + 230);
    const auto *fpf_231 = buffer.data(fpf + 231);
    const auto *fpf_232 = buffer.data(fpf + 232);
    const auto *fpf_233 = buffer.data(fpf + 233);
    const auto *fpf_234 = buffer.data(fpf + 234);
    const auto *fpf_235 = buffer.data(fpf + 235);
    const auto *fpf_236 = buffer.data(fpf + 236);
    const auto *fpf_237 = buffer.data(fpf + 237);
    const auto *fpf_238 = buffer.data(fpf + 238);
    const auto *fpf_239 = buffer.data(fpf + 239);
    const auto *fpf_246 = buffer.data(fpf + 246);
    const auto *fpf_247 = buffer.data(fpf + 247);
    const auto *fpf_248 = buffer.data(fpf + 248);
    const auto *fpf_249 = buffer.data(fpf + 249);

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pa_x, pc_x, pc_y, dpg0_250, dpg0_251, \
                         dpg0_252, dpg1_250, dpg1_251, dpg1_252, fsf_59, \
                         fpf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = pa_x[k] * dpg0_250[k]
                   - f_7 * pc_x[k] * dpg1_250[k];

        t_251[k] = pa_x[k] * dpg0_251[k]
                   - f_7 * pc_x[k] * dpg1_251[k];

        t_252[k] = pa_x[k] * dpg0_252[k]
                   - f_7 * pc_x[k] * dpg1_252[k];

        t_253[k] = f_1 * fsf_59[k]
                   + f_4 * pc_y[k] * fpf_169[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pa_x, pc_x, pc_y, pc_z, dpg0_254, \
                         dpg0_255, dpf_80, dpf_170, dpg1_254, dpg1_255, fsf_50, \
                         fpf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = pa_x[k] * dpg0_254[k]
                   - f_7 * pc_x[k] * dpg1_254[k];

        t_255[k] = pa_x[k] * dpg0_255[k]
                   + f_9 * dpf_170[k]
                   - f_7 * pc_x[k] * dpg1_255[k];

        t_256[k] = f_4 * pc_y[k] * fpf_170[k];

        t_257[k] = f_8 * dpf_80[k]
                   + f_1 * fsf_50[k]
                   + f_4 * pc_z[k] * fpf_170[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, pa_x, pc_x, pc_y, dpg0_260, dpf_175, \
                         dpf_176, dpg1_260, fpd0_102, fpd1_102, fpf_171, fpf_172, \
                         fpf_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_5 * fpd0_102[k]
                   - f_6 * fpd1_102[k]
                   + f_4 * pc_y[k] * fpf_171[k];

        t_259[k] = f_4 * pc_y[k] * fpf_172[k];

        t_260[k] = pa_x[k] * dpg0_260[k]
                   + f_8 * dpf_175[k]
                   - f_7 * pc_x[k] * dpg1_260[k];

        t_261[k] = f_1 * dpf_176[k]
                   + f_4 * pc_x[k] * fpf_176[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, pa_x, pc_x, pc_y, dpg0_265, dpf_177, \
                         dpf_179, dpg1_265, fpf_175, fpf_177, fpf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_1 * dpf_177[k]
                   + f_4 * pc_x[k] * fpf_177[k];

        t_263[k] = f_4 * pc_y[k] * fpf_175[k];

        t_264[k] = f_1 * dpf_179[k]
                   + f_4 * pc_x[k] * fpf_179[k];

        t_265[k] = pa_x[k] * dpg0_265[k]
                   - f_7 * pc_x[k] * dpg1_265[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pa_x, pc_x, pc_y, dpg0_266, dpg0_267, \
                         dpg0_269, dpg1_266, dpg1_267, dpg1_269, \
                         fpf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pa_x[k] * dpg0_266[k]
                   - f_7 * pc_x[k] * dpg1_266[k];

        t_267[k] = pa_x[k] * dpg0_267[k]
                   - f_7 * pc_x[k] * dpg1_267[k];

        t_268[k] = f_4 * pc_y[k] * fpf_179[k];

        t_269[k] = pa_x[k] * dpg0_269[k]
                   - f_7 * pc_x[k] * dpg1_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pb_x, pc_x, pc_z, fsg0_90, fsg0_91, fsf_60, \
                         fsf_61, fsg1_90, fsg1_91, fpf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = pb_x[k] * fsg0_90[k]
                   + f_9 * fsf_60[k]
                   - f_7 * pc_x[k] * fsg1_90[k];

        t_271[k] = pb_x[k] * fsg0_91[k]
                   + f_0 * fsf_61[k]
                   - f_7 * pc_x[k] * fsg1_91[k];

        t_272[k] = f_4 * pc_z[k] * fpf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, pb_x, pc_x, pc_z, fsg0_93, fsg0_95, \
                         fsf_63, fsf_65, fsf_66, fsg1_93, fsg1_95, fpf_181, \
                         fpf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = pb_x[k] * fsg0_93[k]
                   + f_8 * fsf_63[k]
                   - f_7 * pc_x[k] * fsg1_93[k];

        t_274[k] = f_4 * pc_z[k] * fpf_181[k];

        t_275[k] = pb_x[k] * fsg0_95[k]
                   + f_8 * fsf_65[k]
                   - f_7 * pc_x[k] * fsg1_95[k];

        t_276[k] = f_1 * fsf_66[k]
                   + f_4 * pc_x[k] * fpf_186[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, pb_x, pc_x, fsg0_100, fsf_67, fsf_68, \
                         fsf_69, fsg1_100, fpf_187, fpf_188, fpf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_1 * fsf_67[k]
                   + f_4 * pc_x[k] * fpf_187[k];

        t_278[k] = f_1 * fsf_68[k]
                   + f_4 * pc_x[k] * fpf_188[k];

        t_279[k] = f_1 * fsf_69[k]
                   + f_4 * pc_x[k] * fpf_189[k];

        t_280[k] = pb_x[k] * fsg0_100[k]
                   - f_7 * pc_x[k] * fsg1_100[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pb_x, pc_x, pc_y, pc_z, dpf_99, fsg0_102, \
                         fsg0_104, fsg1_102, fsg1_104, fpf_186, \
                         fpf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_4 * pc_z[k] * fpf_186[k];

        t_282[k] = pb_x[k] * fsg0_102[k]
                   - f_7 * pc_x[k] * fsg1_102[k];

        t_283[k] = f_0 * dpf_99[k]
                   + f_4 * pc_y[k] * fpf_189[k];

        t_284[k] = pb_x[k] * fsg0_104[k]
                   - f_7 * pc_x[k] * fsg1_104[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, pc_x, pc_z, fpd0_114, fpd0_115, \
                         fpd0_117, fpd1_114, fpd1_115, fpd1_117, fpf_190, fpf_191, \
                         fpf_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_2 * fpd0_114[k]
                   - f_3 * fpd1_114[k]
                   + f_4 * pc_x[k] * fpf_190[k];

        t_286[k] = f_12 * fpd0_115[k]
                   - f_13 * fpd1_115[k]
                   + f_4 * pc_x[k] * fpf_191[k];

        t_287[k] = f_4 * pc_z[k] * fpf_190[k];

        t_288[k] = f_5 * fpd0_117[k]
                   - f_6 * fpd1_117[k]
                   + f_4 * pc_x[k] * fpf_193[k];

        t_289[k] = f_4 * pc_z[k] * fpf_191[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pc_x, fpd0_119, fpd1_119, fpf_195, \
                         fpf_196, fpf_197, fpf_198, fpf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_5 * fpd0_119[k]
                   - f_6 * fpd1_119[k]
                   + f_4 * pc_x[k] * fpf_195[k];

        t_291[k] = f_4 * pc_x[k] * fpf_196[k];

        t_292[k] = f_4 * pc_x[k] * fpf_197[k];

        t_293[k] = f_4 * pc_x[k] * fpf_198[k];

        t_294[k] = f_4 * pc_x[k] * fpf_199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pc_y, pc_z, dpf_106, dpf_109, fsf_66, \
                         fsf_69, fpd0_117, fpd1_117, fpf_196, fpf_197, \
                         fpf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_0 * dpf_106[k]
                   + f_1 * fsf_66[k]
                   + f_2 * fpd0_117[k]
                   - f_3 * fpd1_117[k]
                   + f_4 * pc_y[k] * fpf_196[k];

        t_296[k] = f_4 * pc_z[k] * fpf_196[k];

        t_297[k] = f_5 * fpd0_117[k]
                   - f_6 * fpd1_117[k]
                   + f_4 * pc_z[k] * fpf_197[k];

        t_298[k] = f_0 * dpf_109[k]
                   + f_1 * fsf_69[k]
                   + f_4 * pc_y[k] * fpf_199[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pb_z, pc_z, fsg0_90, fsg0_91, fsf_60, \
                         fsg1_90, fsg1_91, fpd0_119, fpd1_119, fpf_199, \
                         fpf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_2 * fpd0_119[k]
                   - f_3 * fpd1_119[k]
                   + f_4 * pc_z[k] * fpf_199[k];

        t_300[k] = pb_z[k] * fsg0_90[k]
                   - f_7 * pc_z[k] * fsg1_90[k];

        t_301[k] = pb_z[k] * fsg0_91[k]
                   - f_7 * pc_z[k] * fsg1_91[k];

        t_302[k] = f_1 * fsf_60[k]
                   + f_4 * pc_z[k] * fpf_200[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pb_z, pc_x, pc_z, fsg0_93, fsf_61, \
                         fsg1_93, fpd0_125, fpd1_125, fpf_201, fpf_205, \
                         fpf_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pb_z[k] * fsg0_93[k]
                   - f_7 * pc_z[k] * fsg1_93[k];

        t_304[k] = f_1 * fsf_61[k]
                   + f_4 * pc_z[k] * fpf_201[k];

        t_305[k] = f_5 * fpd0_125[k]
                   - f_6 * fpd1_125[k]
                   + f_4 * pc_x[k] * fpf_205[k];

        t_306[k] = f_4 * pc_x[k] * fpf_206[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, pb_z, pc_x, pc_z, fsg0_100, \
                         fsf_66, fsg1_100, fpf_206, fpf_207, fpf_208, \
                         fpf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_4 * pc_x[k] * fpf_207[k];

        t_308[k] = f_4 * pc_x[k] * fpf_208[k];

        t_309[k] = f_4 * pc_x[k] * fpf_209[k];

        t_310[k] = pb_z[k] * fsg0_100[k]
                   - f_7 * pc_z[k] * fsg1_100[k];

        t_311[k] = f_1 * fsf_66[k]
                   + f_4 * pc_z[k] * fpf_206[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_z, pc_y, pc_z, dpf_119, fsg0_102, fsg0_104, \
                         fsf_67, fsf_69, fsg1_102, fsg1_104, fpf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = pb_z[k] * fsg0_102[k]
                   + f_8 * fsf_67[k]
                   - f_7 * pc_z[k] * fsg1_102[k];

        t_313[k] = f_0 * dpf_119[k]
                   + f_4 * pc_y[k] * fpf_209[k];

        t_314[k] = pb_z[k] * fsg0_104[k]
                   + f_9 * fsf_69[k]
                   - f_7 * pc_z[k] * fsg1_104[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pa_z, pb_x, pc_x, pc_z, dpg0_135, dpg0_136, \
                         dpg1_135, dpg1_136, fsg0_107, fsf_72, \
                         fsg1_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * dpg0_135[k]
                   - f_7 * pc_z[k] * dpg1_135[k];

        t_316[k] = pa_z[k] * dpg0_136[k]
                   - f_7 * pc_z[k] * dpg1_136[k];

        t_317[k] = pb_x[k] * fsg0_107[k]
                   + f_0 * fsf_72[k]
                   - f_7 * pc_x[k] * fsg1_107[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pa_z, pb_x, pc_x, pc_z, dpg0_138, dpg0_139, \
                         dpf_91, dpg1_138, dpg1_139, fsg0_110, fsf_75, \
                         fsg1_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pa_z[k] * dpg0_138[k]
                   - f_7 * pc_z[k] * dpg1_138[k];

        t_319[k] = pa_z[k] * dpg0_139[k]
                   + f_1 * dpf_91[k]
                   - f_7 * pc_z[k] * dpg1_139[k];

        t_320[k] = pb_x[k] * fsg0_110[k]
                   + f_8 * fsf_75[k]
                   - f_7 * pc_x[k] * fsg1_110[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pc_x, fsf_76, fsf_77, fsf_78, fsf_79, \
                         fpf_216, fpf_217, fpf_218, fpf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_1 * fsf_76[k]
                   + f_4 * pc_x[k] * fpf_216[k];

        t_322[k] = f_1 * fsf_77[k]
                   + f_4 * pc_x[k] * fpf_217[k];

        t_323[k] = f_1 * fsf_78[k]
                   + f_4 * pc_x[k] * fpf_218[k];

        t_324[k] = f_1 * fsf_79[k]
                   + f_4 * pc_x[k] * fpf_219[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, pa_z, pb_x, pc_x, pc_z, dpg0_145, dpf_96, \
                         dpg1_145, fsg0_117, fsg1_117, fpf_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pa_z[k] * dpg0_145[k]
                   - f_7 * pc_z[k] * dpg1_145[k];

        t_326[k] = f_1 * dpf_96[k]
                   + f_4 * pc_z[k] * fpf_216[k];

        t_327[k] = pb_x[k] * fsg0_117[k]
                   - f_7 * pc_x[k] * fsg1_117[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pa_z, pb_x, pc_x, pc_y, pc_z, dpg0_150, dpf_129, \
                         dpg1_150, fsg0_119, fsg1_119, fpf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_8 * dpf_129[k]
                   + f_4 * pc_y[k] * fpf_219[k];

        t_329[k] = pb_x[k] * fsg0_119[k]
                   - f_7 * pc_x[k] * fsg1_119[k];

        t_330[k] = pa_z[k] * dpg0_150[k]
                   - f_7 * pc_z[k] * dpg1_150[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pa_z, pc_x, pc_z, dpg0_151, dpg0_153, dpg1_151, \
                         dpg1_153, fpd0_134, fpd1_134, fpf_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_z[k] * dpg0_151[k]
                   - f_7 * pc_z[k] * dpg1_151[k];

        t_332[k] = f_12 * fpd0_134[k]
                   - f_13 * fpd1_134[k]
                   + f_4 * pc_x[k] * fpf_222[k];

        t_333[k] = pa_z[k] * dpg0_153[k]
                   - f_7 * pc_z[k] * dpg1_153[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, pc_x, fpd0_136, fpd0_137, \
                         fpd1_136, fpd1_137, fpf_224, fpf_225, fpf_226, fpf_227, \
                         fpf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_5 * fpd0_136[k]
                   - f_6 * fpd1_136[k]
                   + f_4 * pc_x[k] * fpf_224[k];

        t_335[k] = f_5 * fpd0_137[k]
                   - f_6 * fpd1_137[k]
                   + f_4 * pc_x[k] * fpf_225[k];

        t_336[k] = f_4 * pc_x[k] * fpf_226[k];

        t_337[k] = f_4 * pc_x[k] * fpf_227[k];

        t_338[k] = f_4 * pc_x[k] * fpf_228[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, pa_z, pc_x, pc_z, dpg0_160, dpg0_162, \
                         dpf_106, dpf_107, dpg1_160, dpg1_162, fpf_226, \
                         fpf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_4 * pc_x[k] * fpf_229[k];

        t_340[k] = pa_z[k] * dpg0_160[k]
                   - f_7 * pc_z[k] * dpg1_160[k];

        t_341[k] = f_1 * dpf_106[k]
                   + f_4 * pc_z[k] * fpf_226[k];

        t_342[k] = pa_z[k] * dpg0_162[k]
                   + f_8 * dpf_107[k]
                   - f_7 * pc_z[k] * dpg1_162[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pc_x, pc_y, pc_z, dpf_109, dpf_139, fsf_79, \
                         fpd0_137, fpd0_138, fpd1_137, fpd1_138, fpf_229, \
                         fpf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_8 * dpf_139[k]
                   + f_1 * fsf_79[k]
                   + f_4 * pc_y[k] * fpf_229[k];

        t_344[k] = f_1 * dpf_109[k]
                   + f_2 * fpd0_137[k]
                   - f_3 * fpd1_137[k]
                   + f_4 * pc_z[k] * fpf_229[k];

        t_345[k] = f_2 * fpd0_138[k]
                   - f_3 * fpd1_138[k]
                   + f_4 * pc_x[k] * fpf_230[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, pc_x, fpd0_139, fpd0_140, fpd0_141, fpd1_139, \
                         fpd1_140, fpd1_141, fpf_231, fpf_232, \
                         fpf_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_12 * fpd0_139[k]
                   - f_13 * fpd1_139[k]
                   + f_4 * pc_x[k] * fpf_231[k];

        t_347[k] = f_12 * fpd0_140[k]
                   - f_13 * fpd1_140[k]
                   + f_4 * pc_x[k] * fpf_232[k];

        t_348[k] = f_5 * fpd0_141[k]
                   - f_6 * fpd1_141[k]
                   + f_4 * pc_x[k] * fpf_233[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, pc_x, fpd0_142, fpd0_143, \
                         fpd1_142, fpd1_143, fpf_234, fpf_235, fpf_236, fpf_237, \
                         fpf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_5 * fpd0_142[k]
                   - f_6 * fpd1_142[k]
                   + f_4 * pc_x[k] * fpf_234[k];

        t_350[k] = f_5 * fpd0_143[k]
                   - f_6 * fpd1_143[k]
                   + f_4 * pc_x[k] * fpf_235[k];

        t_351[k] = f_4 * pc_x[k] * fpf_236[k];

        t_352[k] = f_4 * pc_x[k] * fpf_237[k];

        t_353[k] = f_4 * pc_x[k] * fpf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_x, pc_y, pc_z, dpf_116, dpf_146, fsf_76, \
                         fpd0_141, fpd1_141, fpf_236, fpf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_4 * pc_x[k] * fpf_239[k];

        t_355[k] = f_8 * dpf_146[k]
                   + f_2 * fpd0_141[k]
                   - f_3 * fpd1_141[k]
                   + f_4 * pc_y[k] * fpf_236[k];

        t_356[k] = f_1 * dpf_116[k]
                   + f_1 * fsf_76[k]
                   + f_4 * pc_z[k] * fpf_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pa_y, pc_y, ppg0_134, ppg1_134, dpg0_224, \
                         dpf_148, dpf_149, dpg1_224, fpd0_143, fpd1_143, fpf_238, \
                         fpf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_8 * dpf_148[k]
                   + f_5 * fpd0_143[k]
                   - f_6 * fpd1_143[k]
                   + f_4 * pc_y[k] * fpf_238[k];

        t_358[k] = f_8 * dpf_149[k]
                   + f_4 * pc_y[k] * fpf_239[k];

        t_359[k] = f_10 * ppg0_134[k]
                   - f_11 * ppg1_134[k]
                   + pa_y[k] * dpg0_224[k]
                   - f_7 * pc_y[k] * dpg1_224[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pa_y, pc_y, dpg0_225, dpg0_226, dpg0_227, \
                         dpg0_228, dpf_150, dpf_151, dpg1_225, dpg1_226, dpg1_227, \
                         dpg1_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pa_y[k] * dpg0_225[k]
                   - f_7 * pc_y[k] * dpg1_225[k];

        t_361[k] = pa_y[k] * dpg0_226[k]
                   + f_1 * dpf_150[k]
                   - f_7 * pc_y[k] * dpg1_226[k];

        t_362[k] = pa_y[k] * dpg0_227[k]
                   - f_7 * pc_y[k] * dpg1_227[k];

        t_363[k] = pa_y[k] * dpg0_228[k]
                   + f_8 * dpf_151[k]
                   - f_7 * pc_y[k] * dpg1_228[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, pa_y, pc_x, pc_y, dpg0_229, dpg0_230, \
                         dpf_152, dpg1_229, dpg1_230, fsf_86, fsf_87, fpf_246, \
                         fpf_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = pa_y[k] * dpg0_229[k]
                   + f_1 * dpf_152[k]
                   - f_7 * pc_y[k] * dpg1_229[k];

        t_365[k] = pa_y[k] * dpg0_230[k]
                   - f_7 * pc_y[k] * dpg1_230[k];

        t_366[k] = f_1 * fsf_86[k]
                   + f_4 * pc_x[k] * fpf_246[k];

        t_367[k] = f_1 * fsf_87[k]
                   + f_4 * pc_x[k] * fpf_247[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pb_x, pc_x, pc_z, dpf_126, fsg0_130, \
                         fsf_88, fsf_89, fsg1_130, fpf_246, fpf_248, \
                         fpf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_1 * fsf_88[k]
                   + f_4 * pc_x[k] * fpf_248[k];

        t_369[k] = f_1 * fsf_89[k]
                   + f_4 * pc_x[k] * fpf_249[k];

        t_370[k] = pb_x[k] * fsg0_130[k]
                   - f_7 * pc_x[k] * fsg1_130[k];

        t_371[k] = f_8 * dpf_126[k]
                   + f_4 * pc_z[k] * fpf_246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pa_y, pb_x, pc_x, pc_y, dpg0_239, dpf_159, \
                         dpg1_239, fsg0_132, fsg1_132, fpf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = pb_x[k] * fsg0_132[k]
                   - f_7 * pc_x[k] * fsg1_132[k];

        t_373[k] = f_1 * dpf_159[k]
                   + f_4 * pc_y[k] * fpf_249[k];

        t_374[k] = pa_y[k] * dpg0_239[k]
                   - f_7 * pc_y[k] * dpg1_239[k];
    }
}

static auto
compute_prim_fpg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpg0, const size_t dpf,
                                                          const size_t dpg1, const size_t fsg0,
                                                          const size_t fsf, const size_t fsg1,
                                                          const size_t fpd0, const size_t fpd1,
                                                          const size_t fpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpg0_255 = buffer.data(dpg0 + 255);
    const auto *dpg0_257 = buffer.data(dpg0 + 257);
    const auto *dpg0_260 = buffer.data(dpg0 + 260);
    const auto *dpg0_265 = buffer.data(dpg0 + 265);
    const auto *dpg0_267 = buffer.data(dpg0 + 267);
    const auto *dpg0_269 = buffer.data(dpg0 + 269);

    const auto *dpf_136 = buffer.data(dpf + 136);
    const auto *dpf_139 = buffer.data(dpf + 139);
    const auto *dpf_146 = buffer.data(dpf + 146);
    const auto *dpf_166 = buffer.data(dpf + 166);
    const auto *dpf_168 = buffer.data(dpf + 168);
    const auto *dpf_169 = buffer.data(dpf + 169);
    const auto *dpf_176 = buffer.data(dpf + 176);
    const auto *dpf_178 = buffer.data(dpf + 178);
    const auto *dpf_179 = buffer.data(dpf + 179);

    const auto *dpg1_255 = buffer.data(dpg1 + 255);
    const auto *dpg1_257 = buffer.data(dpg1 + 257);
    const auto *dpg1_260 = buffer.data(dpg1 + 260);
    const auto *dpg1_265 = buffer.data(dpg1 + 265);
    const auto *dpg1_267 = buffer.data(dpg1 + 267);
    const auto *dpg1_269 = buffer.data(dpg1 + 269);

    const auto *fsg0_135 = buffer.data(fsg0 + 135);
    const auto *fsg0_137 = buffer.data(fsg0 + 137);
    const auto *fsg0_138 = buffer.data(fsg0 + 138);
    const auto *fsg0_140 = buffer.data(fsg0 + 140);
    const auto *fsg0_145 = buffer.data(fsg0 + 145);
    const auto *fsg0_146 = buffer.data(fsg0 + 146);
    const auto *fsg0_147 = buffer.data(fsg0 + 147);
    const auto *fsg0_149 = buffer.data(fsg0 + 149);

    const auto *fsf_86 = buffer.data(fsf + 86);
    const auto *fsf_88 = buffer.data(fsf + 88);
    const auto *fsf_89 = buffer.data(fsf + 89);
    const auto *fsf_90 = buffer.data(fsf + 90);
    const auto *fsf_92 = buffer.data(fsf + 92);
    const auto *fsf_93 = buffer.data(fsf + 93);
    const auto *fsf_95 = buffer.data(fsf + 95);
    const auto *fsf_96 = buffer.data(fsf + 96);
    const auto *fsf_97 = buffer.data(fsf + 97);
    const auto *fsf_98 = buffer.data(fsf + 98);
    const auto *fsf_99 = buffer.data(fsf + 99);

    const auto *fsg1_135 = buffer.data(fsg1 + 135);
    const auto *fsg1_137 = buffer.data(fsg1 + 137);
    const auto *fsg1_138 = buffer.data(fsg1 + 138);
    const auto *fsg1_140 = buffer.data(fsg1 + 140);
    const auto *fsg1_145 = buffer.data(fsg1 + 145);
    const auto *fsg1_146 = buffer.data(fsg1 + 146);
    const auto *fsg1_147 = buffer.data(fsg1 + 147);
    const auto *fsg1_149 = buffer.data(fsg1 + 149);

    const auto *fpd0_150 = buffer.data(fpd0 + 150);
    const auto *fpd0_151 = buffer.data(fpd0 + 151);
    const auto *fpd0_152 = buffer.data(fpd0 + 152);
    const auto *fpd0_153 = buffer.data(fpd0 + 153);
    const auto *fpd0_154 = buffer.data(fpd0 + 154);
    const auto *fpd0_155 = buffer.data(fpd0 + 155);
    const auto *fpd0_157 = buffer.data(fpd0 + 157);
    const auto *fpd0_159 = buffer.data(fpd0 + 159);
    const auto *fpd0_160 = buffer.data(fpd0 + 160);
    const auto *fpd0_171 = buffer.data(fpd0 + 171);
    const auto *fpd0_174 = buffer.data(fpd0 + 174);
    const auto *fpd0_176 = buffer.data(fpd0 + 176);
    const auto *fpd0_177 = buffer.data(fpd0 + 177);
    const auto *fpd0_178 = buffer.data(fpd0 + 178);
    const auto *fpd0_179 = buffer.data(fpd0 + 179);

    const auto *fpd1_150 = buffer.data(fpd1 + 150);
    const auto *fpd1_151 = buffer.data(fpd1 + 151);
    const auto *fpd1_152 = buffer.data(fpd1 + 152);
    const auto *fpd1_153 = buffer.data(fpd1 + 153);
    const auto *fpd1_154 = buffer.data(fpd1 + 154);
    const auto *fpd1_155 = buffer.data(fpd1 + 155);
    const auto *fpd1_157 = buffer.data(fpd1 + 157);
    const auto *fpd1_159 = buffer.data(fpd1 + 159);
    const auto *fpd1_160 = buffer.data(fpd1 + 160);
    const auto *fpd1_171 = buffer.data(fpd1 + 171);
    const auto *fpd1_174 = buffer.data(fpd1 + 174);
    const auto *fpd1_176 = buffer.data(fpd1 + 176);
    const auto *fpd1_177 = buffer.data(fpd1 + 177);
    const auto *fpd1_178 = buffer.data(fpd1 + 178);
    const auto *fpd1_179 = buffer.data(fpd1 + 179);

    const auto *fpf_250 = buffer.data(fpf + 250);
    const auto *fpf_251 = buffer.data(fpf + 251);
    const auto *fpf_252 = buffer.data(fpf + 252);
    const auto *fpf_253 = buffer.data(fpf + 253);
    const auto *fpf_254 = buffer.data(fpf + 254);
    const auto *fpf_255 = buffer.data(fpf + 255);
    const auto *fpf_256 = buffer.data(fpf + 256);
    const auto *fpf_257 = buffer.data(fpf + 257);
    const auto *fpf_258 = buffer.data(fpf + 258);
    const auto *fpf_259 = buffer.data(fpf + 259);
    const auto *fpf_261 = buffer.data(fpf + 261);
    const auto *fpf_263 = buffer.data(fpf + 263);
    const auto *fpf_264 = buffer.data(fpf + 264);
    const auto *fpf_266 = buffer.data(fpf + 266);
    const auto *fpf_267 = buffer.data(fpf + 267);
    const auto *fpf_268 = buffer.data(fpf + 268);
    const auto *fpf_269 = buffer.data(fpf + 269);
    const auto *fpf_270 = buffer.data(fpf + 270);
    const auto *fpf_272 = buffer.data(fpf + 272);
    const auto *fpf_276 = buffer.data(fpf + 276);
    const auto *fpf_277 = buffer.data(fpf + 277);
    const auto *fpf_278 = buffer.data(fpf + 278);
    const auto *fpf_279 = buffer.data(fpf + 279);
    const auto *fpf_280 = buffer.data(fpf + 280);
    const auto *fpf_282 = buffer.data(fpf + 282);
    const auto *fpf_283 = buffer.data(fpf + 283);
    const auto *fpf_286 = buffer.data(fpf + 286);
    const auto *fpf_287 = buffer.data(fpf + 287);
    const auto *fpf_288 = buffer.data(fpf + 288);
    const auto *fpf_289 = buffer.data(fpf + 289);
    const auto *fpf_290 = buffer.data(fpf + 290);
    const auto *fpf_292 = buffer.data(fpf + 292);
    const auto *fpf_293 = buffer.data(fpf + 293);
    const auto *fpf_295 = buffer.data(fpf + 295);
    const auto *fpf_296 = buffer.data(fpf + 296);
    const auto *fpf_297 = buffer.data(fpf + 297);
    const auto *fpf_298 = buffer.data(fpf + 298);
    const auto *fpf_299 = buffer.data(fpf + 299);

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, fpd0_150, fpd0_151, fpd0_152, fpd1_150, \
                         fpd1_151, fpd1_152, fpf_250, fpf_251, \
                         fpf_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_2 * fpd0_150[k]
                   - f_3 * fpd1_150[k]
                   + f_4 * pc_x[k] * fpf_250[k];

        t_376[k] = f_12 * fpd0_151[k]
                   - f_13 * fpd1_151[k]
                   + f_4 * pc_x[k] * fpf_251[k];

        t_377[k] = f_12 * fpd0_152[k]
                   - f_13 * fpd1_152[k]
                   + f_4 * pc_x[k] * fpf_252[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, fpd0_153, fpd0_154, fpd0_155, \
                         fpd1_153, fpd1_154, fpd1_155, fpf_253, fpf_254, fpf_255, \
                         fpf_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_5 * fpd0_153[k]
                   - f_6 * fpd1_153[k]
                   + f_4 * pc_x[k] * fpf_253[k];

        t_379[k] = f_5 * fpd0_154[k]
                   - f_6 * fpd1_154[k]
                   + f_4 * pc_x[k] * fpf_254[k];

        t_380[k] = f_5 * fpd0_155[k]
                   - f_6 * fpd1_155[k]
                   + f_4 * pc_x[k] * fpf_255[k];

        t_381[k] = f_4 * pc_x[k] * fpf_256[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pc_x, pc_y, dpf_166, fsf_86, fpd0_153, \
                         fpd1_153, fpf_256, fpf_257, fpf_258, fpf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_4 * pc_x[k] * fpf_257[k];

        t_383[k] = f_4 * pc_x[k] * fpf_258[k];

        t_384[k] = f_4 * pc_x[k] * fpf_259[k];

        t_385[k] = f_1 * dpf_166[k]
                   + f_1 * fsf_86[k]
                   + f_2 * fpd0_153[k]
                   - f_3 * fpd1_153[k]
                   + f_4 * pc_y[k] * fpf_256[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, pc_y, pc_z, dpf_136, dpf_168, dpf_169, fsf_88, \
                         fsf_89, fpd0_155, fpd1_155, fpf_256, fpf_258, \
                         fpf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_8 * dpf_136[k]
                   + f_4 * pc_z[k] * fpf_256[k];

        t_387[k] = f_1 * dpf_168[k]
                   + f_1 * fsf_88[k]
                   + f_5 * fpd0_155[k]
                   - f_6 * fpd1_155[k]
                   + f_4 * pc_y[k] * fpf_258[k];

        t_388[k] = f_1 * dpf_169[k]
                   + f_1 * fsf_89[k]
                   + f_4 * pc_y[k] * fpf_259[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, pa_y, pc_x, pc_y, pc_z, dpg0_255, dpf_139, \
                         dpg1_255, fpd0_155, fpd0_157, fpd1_155, fpd1_157, fpf_259, \
                         fpf_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_8 * dpf_139[k]
                   + f_2 * fpd0_155[k]
                   - f_3 * fpd1_155[k]
                   + f_4 * pc_z[k] * fpf_259[k];

        t_390[k] = pa_y[k] * dpg0_255[k]
                   - f_7 * pc_y[k] * dpg1_255[k];

        t_391[k] = f_12 * fpd0_157[k]
                   - f_13 * fpd1_157[k]
                   + f_4 * pc_x[k] * fpf_261[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, pa_y, pc_x, pc_y, dpg0_257, dpg1_257, fpd0_159, \
                         fpd0_160, fpd1_159, fpd1_160, fpf_263, \
                         fpf_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = pa_y[k] * dpg0_257[k]
                   - f_7 * pc_y[k] * dpg1_257[k];

        t_393[k] = f_5 * fpd0_159[k]
                   - f_6 * fpd1_159[k]
                   + f_4 * pc_x[k] * fpf_263[k];

        t_394[k] = f_5 * fpd0_160[k]
                   - f_6 * fpd1_160[k]
                   + f_4 * pc_x[k] * fpf_264[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, pa_y, pc_x, pc_y, dpg0_260, \
                         dpg1_260, fpf_266, fpf_267, fpf_268, fpf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_y[k] * dpg0_260[k]
                   - f_7 * pc_y[k] * dpg1_260[k];

        t_396[k] = f_4 * pc_x[k] * fpf_266[k];

        t_397[k] = f_4 * pc_x[k] * fpf_267[k];

        t_398[k] = f_4 * pc_x[k] * fpf_268[k];

        t_399[k] = f_4 * pc_x[k] * fpf_269[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pa_y, pc_y, pc_z, dpg0_265, dpg0_267, dpf_146, \
                         dpf_176, dpf_178, dpg1_265, dpg1_267, fsf_86, \
                         fpf_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = pa_y[k] * dpg0_265[k]
                   + f_9 * dpf_176[k]
                   - f_7 * pc_y[k] * dpg1_265[k];

        t_401[k] = f_8 * dpf_146[k]
                   + f_1 * fsf_86[k]
                   + f_4 * pc_z[k] * fpf_266[k];

        t_402[k] = pa_y[k] * dpg0_267[k]
                   + f_8 * dpf_178[k]
                   - f_7 * pc_y[k] * dpg1_267[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pa_y, pb_x, pc_x, pc_y, dpg0_269, \
                         dpf_179, dpg1_269, fsg0_135, fsf_90, fsg1_135, fpf_269, \
                         fpf_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_1 * dpf_179[k]
                   + f_4 * pc_y[k] * fpf_269[k];

        t_404[k] = pa_y[k] * dpg0_269[k]
                   - f_7 * pc_y[k] * dpg1_269[k];

        t_405[k] = pb_x[k] * fsg0_135[k]
                   + f_9 * fsf_90[k]
                   - f_7 * pc_x[k] * fsg1_135[k];

        t_406[k] = f_4 * pc_y[k] * fpf_270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_x, pc_x, pc_y, fsg0_137, fsg0_138, fsf_92, \
                         fsf_93, fsg1_137, fsg1_138, fpf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = pb_x[k] * fsg0_137[k]
                   + f_0 * fsf_92[k]
                   - f_7 * pc_x[k] * fsg1_137[k];

        t_408[k] = pb_x[k] * fsg0_138[k]
                   + f_8 * fsf_93[k]
                   - f_7 * pc_x[k] * fsg1_138[k];

        t_409[k] = f_4 * pc_y[k] * fpf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pb_x, pc_x, fsg0_140, fsf_95, fsf_96, \
                         fsf_97, fsf_98, fsg1_140, fpf_276, fpf_277, \
                         fpf_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pb_x[k] * fsg0_140[k]
                   + f_8 * fsf_95[k]
                   - f_7 * pc_x[k] * fsg1_140[k];

        t_411[k] = f_1 * fsf_96[k]
                   + f_4 * pc_x[k] * fpf_276[k];

        t_412[k] = f_1 * fsf_97[k]
                   + f_4 * pc_x[k] * fpf_277[k];

        t_413[k] = f_1 * fsf_98[k]
                   + f_4 * pc_x[k] * fpf_278[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, pb_x, pc_x, pc_y, fsg0_145, \
                         fsg0_146, fsg0_147, fsf_99, fsg1_145, fsg1_146, fsg1_147, \
                         fpf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_1 * fsf_99[k]
                   + f_4 * pc_x[k] * fpf_279[k];

        t_415[k] = pb_x[k] * fsg0_145[k]
                   - f_7 * pc_x[k] * fsg1_145[k];

        t_416[k] = pb_x[k] * fsg0_146[k]
                   - f_7 * pc_x[k] * fsg1_146[k];

        t_417[k] = pb_x[k] * fsg0_147[k]
                   - f_7 * pc_x[k] * fsg1_147[k];

        t_418[k] = f_4 * pc_y[k] * fpf_279[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, pb_x, pb_y, pc_x, pc_y, fsg0_135, \
                         fsg0_137, fsg0_149, fsf_90, fsg1_135, fsg1_137, fsg1_149, \
                         fpf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = pb_x[k] * fsg0_149[k]
                   - f_7 * pc_x[k] * fsg1_149[k];

        t_420[k] = pb_y[k] * fsg0_135[k]
                   - f_7 * pc_y[k] * fsg1_135[k];

        t_421[k] = f_1 * fsf_90[k]
                   + f_4 * pc_y[k] * fpf_280[k];

        t_422[k] = pb_y[k] * fsg0_137[k]
                   - f_7 * pc_y[k] * fsg1_137[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, pb_y, pc_x, pc_y, fsg0_140, fsf_92, \
                         fsg1_140, fpd0_171, fpd1_171, fpf_282, fpf_283, \
                         fpf_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_5 * fpd0_171[k]
                   - f_6 * fpd1_171[k]
                   + f_4 * pc_x[k] * fpf_283[k];

        t_424[k] = f_1 * fsf_92[k]
                   + f_4 * pc_y[k] * fpf_282[k];

        t_425[k] = pb_y[k] * fsg0_140[k]
                   - f_7 * pc_y[k] * fsg1_140[k];

        t_426[k] = f_4 * pc_x[k] * fpf_286[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pb_y, pc_x, pc_y, fsg0_145, fsf_96, \
                         fsg1_145, fpf_287, fpf_288, fpf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_4 * pc_x[k] * fpf_287[k];

        t_428[k] = f_4 * pc_x[k] * fpf_288[k];

        t_429[k] = f_4 * pc_x[k] * fpf_289[k];

        t_430[k] = pb_y[k] * fsg0_145[k]
                   + f_9 * fsf_96[k]
                   - f_7 * pc_y[k] * fsg1_145[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pb_y, pc_y, fsg0_146, fsg0_147, fsg0_149, \
                         fsf_97, fsf_98, fsf_99, fsg1_146, fsg1_147, fsg1_149, \
                         fpf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = pb_y[k] * fsg0_146[k]
                   + f_0 * fsf_97[k]
                   - f_7 * pc_y[k] * fsg1_146[k];

        t_432[k] = pb_y[k] * fsg0_147[k]
                   + f_8 * fsf_98[k]
                   - f_7 * pc_y[k] * fsg1_147[k];

        t_433[k] = f_1 * fsf_99[k]
                   + f_4 * pc_y[k] * fpf_289[k];

        t_434[k] = pb_y[k] * fsg0_149[k]
                   - f_7 * pc_y[k] * fsg1_149[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, pc_x, pc_y, fpd0_174, fpd0_176, \
                         fpd0_177, fpd1_174, fpd1_176, fpd1_177, fpf_290, fpf_292, \
                         fpf_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_2 * fpd0_174[k]
                   - f_3 * fpd1_174[k]
                   + f_4 * pc_x[k] * fpf_290[k];

        t_436[k] = f_4 * pc_y[k] * fpf_290[k];

        t_437[k] = f_12 * fpd0_176[k]
                   - f_13 * fpd1_176[k]
                   + f_4 * pc_x[k] * fpf_292[k];

        t_438[k] = f_5 * fpd0_177[k]
                   - f_6 * fpd1_177[k]
                   + f_4 * pc_x[k] * fpf_293[k];

        t_439[k] = f_4 * pc_y[k] * fpf_292[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pc_x, fpd0_179, fpd1_179, fpf_295, \
                         fpf_296, fpf_297, fpf_298, fpf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_5 * fpd0_179[k]
                   - f_6 * fpd1_179[k]
                   + f_4 * pc_x[k] * fpf_295[k];

        t_441[k] = f_4 * pc_x[k] * fpf_296[k];

        t_442[k] = f_4 * pc_x[k] * fpf_297[k];

        t_443[k] = f_4 * pc_x[k] * fpf_298[k];

        t_444[k] = f_4 * pc_x[k] * fpf_299[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_y, fpd0_177, fpd0_178, fpd0_179, \
                         fpd1_177, fpd1_178, fpd1_179, fpf_296, fpf_297, fpf_298, \
                         fpf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_2 * fpd0_177[k]
                   - f_3 * fpd1_177[k]
                   + f_4 * pc_y[k] * fpf_296[k];

        t_446[k] = f_12 * fpd0_178[k]
                   - f_13 * fpd1_178[k]
                   + f_4 * pc_y[k] * fpf_297[k];

        t_447[k] = f_5 * fpd0_179[k]
                   - f_6 * fpd1_179[k]
                   + f_4 * pc_y[k] * fpf_298[k];

        t_448[k] = f_4 * pc_y[k] * fpf_299[k];
    }

#pragma omp simd aligned(t_449, pc_z, dpf_179, fsf_99, fpd0_179, fpd1_179, \
                         fpf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_0 * dpf_179[k]
                   + f_1 * fsf_99[k]
                   + f_2 * fpd0_179[k]
                   - f_3 * fpd1_179[k]
                   + f_4 * pc_z[k] * fpf_299[k];
    }
}

auto
compute_prim_fpg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t ppg0,
                                                   const size_t ppg1, const size_t dpg0,
                                                   const size_t dpf, const size_t dpg1,
                                                   const size_t fsg0, const size_t fsf,
                                                   const size_t fsg1, const size_t fpd0,
                                                   const size_t fpd1, const size_t fpf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fpg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, ppg0,
                                                              ppg1, dpg0, dpf, dpg1, fsg0, fsf,
                                                              fsg1, fpd0, fpd1, fpf, ncols,
                                                              gamma, p, q);

    compute_prim_fpg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, ppg0,
                                                              ppg1, dpg0, dpf, dpg1, fsg0, fsf,
                                                              fsg1, fpd0, fpd1, fpf, ncols,
                                                              gamma, p, q);

    compute_prim_fpg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, ppg0,
                                                              ppg1, dpg0, dpf, dpg1, fsg0, fsf,
                                                              fsg1, fpd0, fpd1, fpf, ncols,
                                                              gamma, p, q);

    compute_prim_fpg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, dpg0,
                                                              dpf, dpg1, fsg0, fsf, fsg1, fpd0,
                                                              fpd1, fpf, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
