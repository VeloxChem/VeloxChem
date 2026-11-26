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


#include "SimdThreeCenterElectronRepulsionVrrRecDPG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_dpg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppg0, const size_t ppf,
                                                          const size_t ppg1, const size_t dsg0,
                                                          const size_t dsf, const size_t dsg1,
                                                          const size_t dpd0, const size_t dpd1,
                                                          const size_t dpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 2.0 / q;
    const auto f_9 = 1.0 / gamma;
    const auto f_10 = p / (gamma * q);

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

    const auto *ppg0_0 = buffer.data(ppg0 + 0);
    const auto *ppg0_3 = buffer.data(ppg0 + 3);
    const auto *ppg0_5 = buffer.data(ppg0 + 5);
    const auto *ppg0_6 = buffer.data(ppg0 + 6);
    const auto *ppg0_9 = buffer.data(ppg0 + 9);
    const auto *ppg0_15 = buffer.data(ppg0 + 15);
    const auto *ppg0_18 = buffer.data(ppg0 + 18);
    const auto *ppg0_30 = buffer.data(ppg0 + 30);
    const auto *ppg0_35 = buffer.data(ppg0 + 35);
    const auto *ppg0_63 = buffer.data(ppg0 + 63);
    const auto *ppg0_70 = buffer.data(ppg0 + 70);
    const auto *ppg0_72 = buffer.data(ppg0 + 72);
    const auto *ppg0_73 = buffer.data(ppg0 + 73);
    const auto *ppg0_74 = buffer.data(ppg0 + 74);
    const auto *ppg0_85 = buffer.data(ppg0 + 85);
    const auto *ppg0_87 = buffer.data(ppg0 + 87);
    const auto *ppg0_89 = buffer.data(ppg0 + 89);
    const auto *ppg0_115 = buffer.data(ppg0 + 115);
    const auto *ppg0_116 = buffer.data(ppg0 + 116);
    const auto *ppg0_117 = buffer.data(ppg0 + 117);
    const auto *ppg0_119 = buffer.data(ppg0 + 119);
    const auto *ppg0_125 = buffer.data(ppg0 + 125);

    const auto *ppf_0 = buffer.data(ppf + 0);
    const auto *ppf_6 = buffer.data(ppf + 6);
    const auto *ppf_9 = buffer.data(ppf + 9);
    const auto *ppf_10 = buffer.data(ppf + 10);
    const auto *ppf_16 = buffer.data(ppf + 16);
    const auto *ppf_19 = buffer.data(ppf + 19);
    const auto *ppf_20 = buffer.data(ppf + 20);
    const auto *ppf_26 = buffer.data(ppf + 26);
    const auto *ppf_29 = buffer.data(ppf + 29);
    const auto *ppf_33 = buffer.data(ppf + 33);
    const auto *ppf_36 = buffer.data(ppf + 36);
    const auto *ppf_38 = buffer.data(ppf + 38);
    const auto *ppf_40 = buffer.data(ppf + 40);
    const auto *ppf_43 = buffer.data(ppf + 43);
    const auto *ppf_46 = buffer.data(ppf + 46);
    const auto *ppf_48 = buffer.data(ppf + 48);
    const auto *ppf_49 = buffer.data(ppf + 49);
    const auto *ppf_56 = buffer.data(ppf + 56);
    const auto *ppf_58 = buffer.data(ppf + 58);
    const auto *ppf_59 = buffer.data(ppf + 59);
    const auto *ppf_65 = buffer.data(ppf + 65);
    const auto *ppf_67 = buffer.data(ppf + 67);
    const auto *ppf_69 = buffer.data(ppf + 69);
    const auto *ppf_76 = buffer.data(ppf + 76);
    const auto *ppf_77 = buffer.data(ppf + 77);
    const auto *ppf_79 = buffer.data(ppf + 79);
    const auto *ppf_80 = buffer.data(ppf + 80);
    const auto *ppf_85 = buffer.data(ppf + 85);
    const auto *ppf_86 = buffer.data(ppf + 86);
    const auto *ppf_87 = buffer.data(ppf + 87);

    const auto *ppg1_0 = buffer.data(ppg1 + 0);
    const auto *ppg1_3 = buffer.data(ppg1 + 3);
    const auto *ppg1_5 = buffer.data(ppg1 + 5);
    const auto *ppg1_6 = buffer.data(ppg1 + 6);
    const auto *ppg1_9 = buffer.data(ppg1 + 9);
    const auto *ppg1_15 = buffer.data(ppg1 + 15);
    const auto *ppg1_18 = buffer.data(ppg1 + 18);
    const auto *ppg1_30 = buffer.data(ppg1 + 30);
    const auto *ppg1_35 = buffer.data(ppg1 + 35);
    const auto *ppg1_63 = buffer.data(ppg1 + 63);
    const auto *ppg1_70 = buffer.data(ppg1 + 70);
    const auto *ppg1_72 = buffer.data(ppg1 + 72);
    const auto *ppg1_73 = buffer.data(ppg1 + 73);
    const auto *ppg1_74 = buffer.data(ppg1 + 74);
    const auto *ppg1_85 = buffer.data(ppg1 + 85);
    const auto *ppg1_87 = buffer.data(ppg1 + 87);
    const auto *ppg1_89 = buffer.data(ppg1 + 89);
    const auto *ppg1_115 = buffer.data(ppg1 + 115);
    const auto *ppg1_116 = buffer.data(ppg1 + 116);
    const auto *ppg1_117 = buffer.data(ppg1 + 117);
    const auto *ppg1_119 = buffer.data(ppg1 + 119);
    const auto *ppg1_125 = buffer.data(ppg1 + 125);

    const auto *dsg0_0 = buffer.data(dsg0 + 0);
    const auto *dsg0_3 = buffer.data(dsg0 + 3);
    const auto *dsg0_5 = buffer.data(dsg0 + 5);
    const auto *dsg0_10 = buffer.data(dsg0 + 10);
    const auto *dsg0_12 = buffer.data(dsg0 + 12);
    const auto *dsg0_14 = buffer.data(dsg0 + 14);
    const auto *dsg0_18 = buffer.data(dsg0 + 18);
    const auto *dsg0_35 = buffer.data(dsg0 + 35);

    const auto *dsf_0 = buffer.data(dsf + 0);
    const auto *dsf_1 = buffer.data(dsf + 1);
    const auto *dsf_2 = buffer.data(dsf + 2);
    const auto *dsf_3 = buffer.data(dsf + 3);
    const auto *dsf_5 = buffer.data(dsf + 5);
    const auto *dsf_6 = buffer.data(dsf + 6);
    const auto *dsf_8 = buffer.data(dsf + 8);
    const auto *dsf_9 = buffer.data(dsf + 9);
    const auto *dsf_10 = buffer.data(dsf + 10);
    const auto *dsf_11 = buffer.data(dsf + 11);
    const auto *dsf_13 = buffer.data(dsf + 13);
    const auto *dsf_16 = buffer.data(dsf + 16);
    const auto *dsf_18 = buffer.data(dsf + 18);
    const auto *dsf_20 = buffer.data(dsf + 20);
    const auto *dsf_22 = buffer.data(dsf + 22);
    const auto *dsf_25 = buffer.data(dsf + 25);
    const auto *dsf_27 = buffer.data(dsf + 27);
    const auto *dsf_29 = buffer.data(dsf + 29);

    const auto *dsg1_0 = buffer.data(dsg1 + 0);
    const auto *dsg1_3 = buffer.data(dsg1 + 3);
    const auto *dsg1_5 = buffer.data(dsg1 + 5);
    const auto *dsg1_10 = buffer.data(dsg1 + 10);
    const auto *dsg1_12 = buffer.data(dsg1 + 12);
    const auto *dsg1_14 = buffer.data(dsg1 + 14);
    const auto *dsg1_18 = buffer.data(dsg1 + 18);
    const auto *dsg1_35 = buffer.data(dsg1 + 35);

    const auto *dpd0_0 = buffer.data(dpd0 + 0);
    const auto *dpd0_3 = buffer.data(dpd0 + 3);
    const auto *dpd0_5 = buffer.data(dpd0 + 5);
    const auto *dpd0_17 = buffer.data(dpd0 + 17);
    const auto *dpd0_21 = buffer.data(dpd0 + 21);
    const auto *dpd0_23 = buffer.data(dpd0 + 23);
    const auto *dpd0_24 = buffer.data(dpd0 + 24);
    const auto *dpd0_39 = buffer.data(dpd0 + 39);
    const auto *dpd0_40 = buffer.data(dpd0 + 40);
    const auto *dpd0_41 = buffer.data(dpd0 + 41);
    const auto *dpd0_48 = buffer.data(dpd0 + 48);

    const auto *dpd1_0 = buffer.data(dpd1 + 0);
    const auto *dpd1_3 = buffer.data(dpd1 + 3);
    const auto *dpd1_5 = buffer.data(dpd1 + 5);
    const auto *dpd1_17 = buffer.data(dpd1 + 17);
    const auto *dpd1_21 = buffer.data(dpd1 + 21);
    const auto *dpd1_23 = buffer.data(dpd1 + 23);
    const auto *dpd1_24 = buffer.data(dpd1 + 24);
    const auto *dpd1_39 = buffer.data(dpd1 + 39);
    const auto *dpd1_40 = buffer.data(dpd1 + 40);
    const auto *dpd1_41 = buffer.data(dpd1 + 41);
    const auto *dpd1_48 = buffer.data(dpd1 + 48);

    const auto *dpf_0 = buffer.data(dpf + 0);
    const auto *dpf_1 = buffer.data(dpf + 1);
    const auto *dpf_2 = buffer.data(dpf + 2);
    const auto *dpf_3 = buffer.data(dpf + 3);
    const auto *dpf_5 = buffer.data(dpf + 5);
    const auto *dpf_6 = buffer.data(dpf + 6);
    const auto *dpf_8 = buffer.data(dpf + 8);
    const auto *dpf_9 = buffer.data(dpf + 9);
    const auto *dpf_10 = buffer.data(dpf + 10);
    const auto *dpf_12 = buffer.data(dpf + 12);
    const auto *dpf_13 = buffer.data(dpf + 13);
    const auto *dpf_15 = buffer.data(dpf + 15);
    const auto *dpf_16 = buffer.data(dpf + 16);
    const auto *dpf_19 = buffer.data(dpf + 19);
    const auto *dpf_20 = buffer.data(dpf + 20);
    const auto *dpf_22 = buffer.data(dpf + 22);
    const auto *dpf_23 = buffer.data(dpf + 23);
    const auto *dpf_25 = buffer.data(dpf + 25);
    const auto *dpf_26 = buffer.data(dpf + 26);
    const auto *dpf_28 = buffer.data(dpf + 28);
    const auto *dpf_29 = buffer.data(dpf + 29);
    const auto *dpf_30 = buffer.data(dpf + 30);
    const auto *dpf_31 = buffer.data(dpf + 31);
    const auto *dpf_33 = buffer.data(dpf + 33);
    const auto *dpf_36 = buffer.data(dpf + 36);
    const auto *dpf_37 = buffer.data(dpf + 37);
    const auto *dpf_38 = buffer.data(dpf + 38);
    const auto *dpf_39 = buffer.data(dpf + 39);
    const auto *dpf_40 = buffer.data(dpf + 40);
    const auto *dpf_41 = buffer.data(dpf + 41);
    const auto *dpf_42 = buffer.data(dpf + 42);
    const auto *dpf_43 = buffer.data(dpf + 43);
    const auto *dpf_46 = buffer.data(dpf + 46);
    const auto *dpf_48 = buffer.data(dpf + 48);
    const auto *dpf_49 = buffer.data(dpf + 49);
    const auto *dpf_50 = buffer.data(dpf + 50);
    const auto *dpf_51 = buffer.data(dpf + 51);
    const auto *dpf_53 = buffer.data(dpf + 53);
    const auto *dpf_56 = buffer.data(dpf + 56);
    const auto *dpf_58 = buffer.data(dpf + 58);
    const auto *dpf_59 = buffer.data(dpf + 59);
    const auto *dpf_60 = buffer.data(dpf + 60);
    const auto *dpf_62 = buffer.data(dpf + 62);
    const auto *dpf_65 = buffer.data(dpf + 65);
    const auto *dpf_66 = buffer.data(dpf + 66);
    const auto *dpf_67 = buffer.data(dpf + 67);
    const auto *dpf_68 = buffer.data(dpf + 68);
    const auto *dpf_69 = buffer.data(dpf + 69);
    const auto *dpf_70 = buffer.data(dpf + 70);
    const auto *dpf_72 = buffer.data(dpf + 72);
    const auto *dpf_75 = buffer.data(dpf + 75);
    const auto *dpf_76 = buffer.data(dpf + 76);
    const auto *dpf_77 = buffer.data(dpf + 77);
    const auto *dpf_79 = buffer.data(dpf + 79);
    const auto *dpf_80 = buffer.data(dpf + 80);
    const auto *dpf_81 = buffer.data(dpf + 81);
    const auto *dpf_82 = buffer.data(dpf + 82);
    const auto *dpf_85 = buffer.data(dpf + 85);
    const auto *dpf_86 = buffer.data(dpf + 86);
    const auto *dpf_87 = buffer.data(dpf + 87);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ppf_0, dsf_0, dpd0_0, \
                         dpd1_0, dpf_0, dpf_1, dpf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ppf_0[k]
                 + f_1 * dsf_0[k]
                 + f_2 * dpd0_0[k]
                 - f_3 * dpd1_0[k]
                 + f_4 * pc_x[k] * dpf_0[k];

        t_1[k] = f_4 * pc_y[k] * dpf_0[k];

        t_2[k] = f_4 * pc_z[k] * dpf_0[k];

        t_3[k] = f_5 * dpd0_0[k]
                 - f_6 * dpd1_0[k]
                 + f_4 * pc_y[k] * dpf_1[k];

        t_4[k] = f_4 * pc_y[k] * dpf_2[k];

        t_5[k] = f_5 * dpd0_0[k]
                 - f_6 * dpd1_0[k]
                 + f_4 * pc_z[k] * dpf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, ppf_6, ppf_9, dsf_6, dsf_9, \
                         dpf_3, dpf_5, dpf_6, dpf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ppf_6[k]
                 + f_1 * dsf_6[k]
                 + f_4 * pc_x[k] * dpf_6[k];

        t_7[k] = f_4 * pc_z[k] * dpf_3[k];

        t_8[k] = f_4 * pc_y[k] * dpf_5[k];

        t_9[k] = f_0 * ppf_9[k]
                 + f_1 * dsf_9[k]
                 + f_4 * pc_x[k] * dpf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pc_y, pc_z, dpd0_3, dpd0_5, dpd1_3, \
                         dpd1_5, dpf_6, dpf_8, dpf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * dpd0_3[k]
                  - f_3 * dpd1_3[k]
                  + f_4 * pc_y[k] * dpf_6[k];

        t_11[k] = f_4 * pc_z[k] * dpf_6[k];

        t_12[k] = f_5 * dpd0_5[k]
                  - f_6 * dpd1_5[k]
                  + f_4 * pc_y[k] * dpf_8[k];

        t_13[k] = f_4 * pc_y[k] * dpf_9[k];

        t_14[k] = f_2 * dpd0_5[k]
                  - f_3 * dpd1_5[k]
                  + f_4 * pc_z[k] * dpf_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pc_y, pc_z, dsg0_0, dsg0_3, dsf_0, \
                         dsf_1, dsg1_0, dsg1_3, dpf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * dsg0_0[k]
                  - f_7 * pc_y[k] * dsg1_0[k];

        t_16[k] = f_1 * dsf_0[k]
                  + f_4 * pc_y[k] * dpf_10[k];

        t_17[k] = f_4 * pc_z[k] * dpf_10[k];

        t_18[k] = pb_y[k] * dsg0_3[k]
                  + f_0 * dsf_1[k]
                  - f_7 * pc_y[k] * dsg1_3[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_y, pc_x, pc_y, pc_z, ppf_16, dsg0_5, \
                         dsf_2, dsg1_5, dpf_12, dpf_13, dpf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * dsf_2[k]
                  + f_4 * pc_y[k] * dpf_12[k];

        t_20[k] = pb_y[k] * dsg0_5[k]
                  - f_7 * pc_y[k] * dsg1_5[k];

        t_21[k] = f_0 * ppf_16[k]
                  + f_4 * pc_x[k] * dpf_16[k];

        t_22[k] = f_4 * pc_z[k] * dpf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_x, pc_y, pc_z, ppf_19, dsg0_10, \
                         dsf_5, dsf_6, dsg1_10, dpf_15, dpf_16, \
                         dpf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * dsf_5[k]
                  + f_4 * pc_y[k] * dpf_15[k];

        t_24[k] = f_0 * ppf_19[k]
                  + f_4 * pc_x[k] * dpf_19[k];

        t_25[k] = pb_y[k] * dsg0_10[k]
                  + f_8 * dsf_6[k]
                  - f_7 * pc_y[k] * dsg1_10[k];

        t_26[k] = f_4 * pc_z[k] * dpf_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pc_y, dsg0_12, dsg0_14, dsf_8, dsf_9, \
                         dsg1_12, dsg1_14, dpf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * dsg0_12[k]
                  + f_0 * dsf_8[k]
                  - f_7 * pc_y[k] * dsg1_12[k];

        t_28[k] = f_1 * dsf_9[k]
                  + f_4 * pc_y[k] * dpf_19[k];

        t_29[k] = pb_y[k] * dsg0_14[k]
                  - f_7 * pc_y[k] * dsg1_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_z, pc_y, pc_z, dsg0_0, dsg0_3, \
                         dsf_0, dsg1_0, dsg1_3, dpf_20, dpf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * dsg0_0[k]
                  - f_7 * pc_z[k] * dsg1_0[k];

        t_31[k] = f_4 * pc_y[k] * dpf_20[k];

        t_32[k] = f_1 * dsf_0[k]
                  + f_4 * pc_z[k] * dpf_20[k];

        t_33[k] = pb_z[k] * dsg0_3[k]
                  - f_7 * pc_z[k] * dsg1_3[k];

        t_34[k] = f_4 * pc_y[k] * dpf_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_z, pc_x, pc_y, pc_z, ppf_26, dsg0_5, \
                         dsf_2, dsf_3, dsg1_5, dpf_23, dpf_25, dpf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_z[k] * dsg0_5[k]
                  + f_0 * dsf_2[k]
                  - f_7 * pc_z[k] * dsg1_5[k];

        t_36[k] = f_0 * ppf_26[k]
                  + f_4 * pc_x[k] * dpf_26[k];

        t_37[k] = f_1 * dsf_3[k]
                  + f_4 * pc_z[k] * dpf_23[k];

        t_38[k] = f_4 * pc_y[k] * dpf_25[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_z, pc_x, pc_z, ppf_29, dsg0_10, dsf_6, dsg1_10, \
                         dpf_26, dpf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * ppf_29[k]
                  + f_4 * pc_x[k] * dpf_29[k];

        t_40[k] = pb_z[k] * dsg0_10[k]
                  - f_7 * pc_z[k] * dsg1_10[k];

        t_41[k] = f_1 * dsf_6[k]
                  + f_4 * pc_z[k] * dpf_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_z, pc_y, pc_z, dsg0_14, dsf_9, dsg1_14, dpd0_17, \
                         dpd1_17, dpf_28, dpf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * dpd0_17[k]
                  - f_6 * dpd1_17[k]
                  + f_4 * pc_y[k] * dpf_28[k];

        t_43[k] = f_4 * pc_y[k] * dpf_29[k];

        t_44[k] = pb_z[k] * dsg0_14[k]
                  + f_8 * dsf_9[k]
                  - f_7 * pc_z[k] * dsg1_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pc_y, pc_z, ppg0_0, ppf_0, ppg1_0, \
                         dpf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * ppg0_0[k]
                  - f_7 * pc_y[k] * ppg1_0[k];

        t_46[k] = f_1 * ppf_0[k]
                  + f_4 * pc_y[k] * dpf_30[k];

        t_47[k] = f_4 * pc_z[k] * dpf_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pa_y, pc_x, pc_y, pc_z, ppg0_5, ppf_33, ppg1_5, \
                         dsf_13, dpd0_21, dpd1_21, dpf_31, dpf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_1 * ppf_33[k]
                  + f_1 * dsf_13[k]
                  + f_5 * dpd0_21[k]
                  - f_6 * dpd1_21[k]
                  + f_4 * pc_x[k] * dpf_33[k];

        t_49[k] = f_4 * pc_z[k] * dpf_31[k];

        t_50[k] = pa_y[k] * ppg0_5[k]
                  - f_7 * pc_y[k] * ppg1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pc_x, pc_z, ppf_36, ppf_38, dsf_16, dsf_18, dpf_33, \
                         dpf_36, dpf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_1 * ppf_36[k]
                  + f_1 * dsf_16[k]
                  + f_4 * pc_x[k] * dpf_36[k];

        t_52[k] = f_4 * pc_z[k] * dpf_33[k];

        t_53[k] = f_1 * ppf_38[k]
                  + f_1 * dsf_18[k]
                  + f_4 * pc_x[k] * dpf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pc_y, pc_z, ppg0_9, ppf_6, ppg1_9, \
                         dpd0_21, dpd1_21, dpf_36, dpf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_y[k] * ppg0_9[k]
                  - f_7 * pc_y[k] * ppg1_9[k];

        t_55[k] = f_1 * ppf_6[k]
                  + f_2 * dpd0_21[k]
                  - f_3 * dpd1_21[k]
                  + f_4 * pc_y[k] * dpf_36[k];

        t_56[k] = f_4 * pc_z[k] * dpf_36[k];

        t_57[k] = f_5 * dpd0_21[k]
                  - f_6 * dpd1_21[k]
                  + f_4 * pc_z[k] * dpf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pc_x, pc_y, pc_z, ppf_9, ppf_40, dpd0_23, dpd0_24, \
                         dpd1_23, dpd1_24, dpf_39, dpf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * ppf_9[k]
                  + f_4 * pc_y[k] * dpf_39[k];

        t_59[k] = f_2 * dpd0_23[k]
                  - f_3 * dpd1_23[k]
                  + f_4 * pc_z[k] * dpf_39[k];

        t_60[k] = f_1 * ppf_40[k]
                  + f_2 * dpd0_24[k]
                  - f_3 * dpd1_24[k]
                  + f_4 * pc_x[k] * dpf_40[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pc_x, pc_y, pc_z, ppg0_63, ppf_10, \
                         ppf_43, ppg1_63, dsf_10, dpf_40, dpf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_1 * ppf_10[k]
                  + f_1 * dsf_10[k]
                  + f_4 * pc_y[k] * dpf_40[k];

        t_62[k] = f_4 * pc_z[k] * dpf_40[k];

        t_63[k] = pa_x[k] * ppg0_63[k]
                  + f_0 * ppf_43[k]
                  - f_7 * pc_x[k] * ppg1_63[k];

        t_64[k] = f_4 * pc_z[k] * dpf_41[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_x, pc_z, ppf_46, ppf_48, dpd0_24, dpd1_24, \
                         dpf_42, dpf_43, dpf_46, dpf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_5 * dpd0_24[k]
                  - f_6 * dpd1_24[k]
                  + f_4 * pc_z[k] * dpf_42[k];

        t_66[k] = f_1 * ppf_46[k]
                  + f_4 * pc_x[k] * dpf_46[k];

        t_67[k] = f_4 * pc_z[k] * dpf_43[k];

        t_68[k] = f_1 * ppf_48[k]
                  + f_4 * pc_x[k] * dpf_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pc_x, pc_z, ppg0_70, ppg0_72, ppf_49, \
                         ppg1_70, ppg1_72, dpf_46, dpf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * ppf_49[k]
                  + f_4 * pc_x[k] * dpf_49[k];

        t_70[k] = pa_x[k] * ppg0_70[k]
                  - f_7 * pc_x[k] * ppg1_70[k];

        t_71[k] = f_4 * pc_z[k] * dpf_46[k];

        t_72[k] = pa_x[k] * ppg0_72[k]
                  - f_7 * pc_x[k] * ppg1_72[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pa_y, pc_x, pc_y, ppg0_30, ppg0_73, \
                         ppg0_74, ppf_20, ppg1_30, ppg1_73, ppg1_74, \
                         dpf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_x[k] * ppg0_73[k]
                  - f_7 * pc_x[k] * ppg1_73[k];

        t_74[k] = pa_x[k] * ppg0_74[k]
                  - f_7 * pc_x[k] * ppg1_74[k];

        t_75[k] = pa_y[k] * ppg0_30[k]
                  - f_7 * pc_y[k] * ppg1_30[k];

        t_76[k] = f_1 * ppf_20[k]
                  + f_4 * pc_y[k] * dpf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pb_z, pc_y, pc_z, ppg0_35, ppg1_35, \
                         dsg0_18, dsf_10, dsf_11, dsg1_18, dpf_50, \
                         dpf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_1 * dsf_10[k]
                  + f_4 * pc_z[k] * dpf_50[k];

        t_78[k] = pb_z[k] * dsg0_18[k]
                  - f_7 * pc_z[k] * dsg1_18[k];

        t_79[k] = f_1 * dsf_11[k]
                  + f_4 * pc_z[k] * dpf_51[k];

        t_80[k] = pa_y[k] * ppg0_35[k]
                  - f_7 * pc_y[k] * ppg1_35[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_z, ppf_56, ppf_58, ppf_59, dsf_13, \
                         dpf_53, dpf_56, dpf_58, dpf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_1 * ppf_56[k]
                  + f_4 * pc_x[k] * dpf_56[k];

        t_82[k] = f_1 * dsf_13[k]
                  + f_4 * pc_z[k] * dpf_53[k];

        t_83[k] = f_1 * ppf_58[k]
                  + f_4 * pc_x[k] * dpf_58[k];

        t_84[k] = f_1 * ppf_59[k]
                  + f_4 * pc_x[k] * dpf_59[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_x, pc_x, pc_y, pc_z, ppg0_85, ppg0_87, \
                         ppf_29, ppg1_85, ppg1_87, dsf_16, dpf_56, \
                         dpf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pa_x[k] * ppg0_85[k]
                  - f_7 * pc_x[k] * ppg1_85[k];

        t_86[k] = f_1 * dsf_16[k]
                  + f_4 * pc_z[k] * dpf_56[k];

        t_87[k] = pa_x[k] * ppg0_87[k]
                  - f_7 * pc_x[k] * ppg1_87[k];

        t_88[k] = f_1 * ppf_29[k]
                  + f_4 * pc_y[k] * dpf_59[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_x, pa_z, pc_x, pc_y, pc_z, ppg0_0, \
                         ppg0_89, ppf_0, ppg1_0, ppg1_89, dpf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_x[k] * ppg0_89[k]
                  - f_7 * pc_x[k] * ppg1_89[k];

        t_90[k] = pa_z[k] * ppg0_0[k]
                  - f_7 * pc_z[k] * ppg1_0[k];

        t_91[k] = f_4 * pc_y[k] * dpf_60[k];

        t_92[k] = f_1 * ppf_0[k]
                  + f_4 * pc_z[k] * dpf_60[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pa_z, pc_x, pc_y, pc_z, ppg0_3, ppf_65, ppg1_3, \
                         dsf_25, dpd0_41, dpd1_41, dpf_62, dpf_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_z[k] * ppg0_3[k]
                  - f_7 * pc_z[k] * ppg1_3[k];

        t_94[k] = f_4 * pc_y[k] * dpf_62[k];

        t_95[k] = f_1 * ppf_65[k]
                  + f_1 * dsf_25[k]
                  + f_5 * dpd0_41[k]
                  - f_6 * dpd1_41[k]
                  + f_4 * pc_x[k] * dpf_65[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_z, pc_x, pc_y, pc_z, ppg0_6, ppf_67, ppg1_6, \
                         dsf_27, dpf_65, dpf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * ppg0_6[k]
                  - f_7 * pc_z[k] * ppg1_6[k];

        t_97[k] = f_1 * ppf_67[k]
                  + f_1 * dsf_27[k]
                  + f_4 * pc_x[k] * dpf_67[k];

        t_98[k] = f_4 * pc_y[k] * dpf_65[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pc_x, pc_y, ppf_69, dsf_29, dpd0_39, dpd0_40, \
                         dpd1_39, dpd1_40, dpf_66, dpf_67, dpf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_1 * ppf_69[k]
                  + f_1 * dsf_29[k]
                  + f_4 * pc_x[k] * dpf_69[k];

        t_100[k] = f_2 * dpd0_39[k]
                   - f_3 * dpd1_39[k]
                   + f_4 * pc_y[k] * dpf_66[k];

        t_101[k] = f_9 * dpd0_40[k]
                   - f_10 * dpd1_40[k]
                   + f_4 * pc_y[k] * dpf_67[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pa_z, pc_y, pc_z, ppg0_15, ppf_9, \
                         ppg1_15, dpd0_41, dpd1_41, dpf_68, dpf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_5 * dpd0_41[k]
                   - f_6 * dpd1_41[k]
                   + f_4 * pc_y[k] * dpf_68[k];

        t_103[k] = f_4 * pc_y[k] * dpf_69[k];

        t_104[k] = f_1 * ppf_9[k]
                   + f_2 * dpd0_41[k]
                   - f_3 * dpd1_41[k]
                   + f_4 * pc_z[k] * dpf_69[k];

        t_105[k] = pa_z[k] * ppg0_15[k]
                   - f_7 * pc_z[k] * ppg1_15[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_z, pc_y, pc_z, ppg0_18, ppf_10, \
                         ppg1_18, dsf_20, dsf_22, dpf_70, dpf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_1 * dsf_20[k]
                   + f_4 * pc_y[k] * dpf_70[k];

        t_107[k] = f_1 * ppf_10[k]
                   + f_4 * pc_z[k] * dpf_70[k];

        t_108[k] = pa_z[k] * ppg0_18[k]
                   - f_7 * pc_z[k] * ppg1_18[k];

        t_109[k] = f_1 * dsf_22[k]
                   + f_4 * pc_y[k] * dpf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_y, pc_x, pc_y, ppf_76, ppf_77, \
                         dsg0_35, dsf_25, dsg1_35, dpf_75, dpf_76, \
                         dpf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_y[k] * dsg0_35[k]
                   - f_7 * pc_y[k] * dsg1_35[k];

        t_111[k] = f_1 * ppf_76[k]
                   + f_4 * pc_x[k] * dpf_76[k];

        t_112[k] = f_1 * ppf_77[k]
                   + f_4 * pc_x[k] * dpf_77[k];

        t_113[k] = f_1 * dsf_25[k]
                   + f_4 * pc_y[k] * dpf_75[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_x, pc_x, ppg0_115, ppg0_116, ppg0_117, \
                         ppf_79, ppg1_115, ppg1_116, ppg1_117, dpf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_1 * ppf_79[k]
                   + f_4 * pc_x[k] * dpf_79[k];

        t_115[k] = pa_x[k] * ppg0_115[k]
                   - f_7 * pc_x[k] * ppg1_115[k];

        t_116[k] = pa_x[k] * ppg0_116[k]
                   - f_7 * pc_x[k] * ppg1_116[k];

        t_117[k] = pa_x[k] * ppg0_117[k]
                   - f_7 * pc_x[k] * ppg1_117[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_x, pc_x, pc_y, ppg0_119, ppf_80, \
                         ppg1_119, dsf_29, dpd0_48, dpd1_48, dpf_79, \
                         dpf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_1 * dsf_29[k]
                   + f_4 * pc_y[k] * dpf_79[k];

        t_119[k] = pa_x[k] * ppg0_119[k]
                   - f_7 * pc_x[k] * ppg1_119[k];

        t_120[k] = f_1 * ppf_80[k]
                   + f_2 * dpd0_48[k]
                   - f_3 * dpd1_48[k]
                   + f_4 * pc_x[k] * dpf_80[k];

        t_121[k] = f_4 * pc_y[k] * dpf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pc_y, pc_z, ppf_20, dsf_20, dpd0_48, dpd1_48, \
                         dpf_80, dpf_81, dpf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_1 * ppf_20[k]
                   + f_1 * dsf_20[k]
                   + f_4 * pc_z[k] * dpf_80[k];

        t_123[k] = f_5 * dpd0_48[k]
                   - f_6 * dpd1_48[k]
                   + f_4 * pc_y[k] * dpf_81[k];

        t_124[k] = f_4 * pc_y[k] * dpf_82[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_x, pc_x, pc_y, ppg0_125, ppf_85, \
                         ppf_86, ppf_87, ppg1_125, dpf_85, dpf_86, \
                         dpf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_x[k] * ppg0_125[k]
                   + f_0 * ppf_85[k]
                   - f_7 * pc_x[k] * ppg1_125[k];

        t_126[k] = f_1 * ppf_86[k]
                   + f_4 * pc_x[k] * dpf_86[k];

        t_127[k] = f_1 * ppf_87[k]
                   + f_4 * pc_x[k] * dpf_87[k];

        t_128[k] = f_4 * pc_y[k] * dpf_85[k];
    }
}

static auto
compute_prim_dpg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppg0, const size_t ppf,
                                                          const size_t ppg1, const size_t dsg0,
                                                          const size_t dsf, const size_t dsg1,
                                                          const size_t dpd0, const size_t dpd1,
                                                          const size_t dpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 2.0 / q;
    const auto f_9 = 1.0 / gamma;
    const auto f_10 = p / (gamma * q);
    const auto f_11 = 1.5 / q;

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
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppg0_46 = buffer.data(ppg0 + 46);
    const auto *ppg0_48 = buffer.data(ppg0 + 48);
    const auto *ppg0_55 = buffer.data(ppg0 + 55);
    const auto *ppg0_61 = buffer.data(ppg0 + 61);
    const auto *ppg0_63 = buffer.data(ppg0 + 63);
    const auto *ppg0_70 = buffer.data(ppg0 + 70);
    const auto *ppg0_90 = buffer.data(ppg0 + 90);
    const auto *ppg0_92 = buffer.data(ppg0 + 92);
    const auto *ppg0_94 = buffer.data(ppg0 + 94);
    const auto *ppg0_95 = buffer.data(ppg0 + 95);
    const auto *ppg0_104 = buffer.data(ppg0 + 104);
    const auto *ppg0_120 = buffer.data(ppg0 + 120);
    const auto *ppg0_122 = buffer.data(ppg0 + 122);
    const auto *ppg0_125 = buffer.data(ppg0 + 125);
    const auto *ppg0_130 = buffer.data(ppg0 + 130);
    const auto *ppg0_131 = buffer.data(ppg0 + 131);
    const auto *ppg0_132 = buffer.data(ppg0 + 132);
    const auto *ppg0_134 = buffer.data(ppg0 + 134);

    const auto *ppf_36 = buffer.data(ppf + 36);
    const auto *ppf_39 = buffer.data(ppf + 39);
    const auto *ppf_46 = buffer.data(ppf + 46);
    const auto *ppf_49 = buffer.data(ppf + 49);
    const auto *ppf_56 = buffer.data(ppf + 56);
    const auto *ppf_59 = buffer.data(ppf + 59);
    const auto *ppf_62 = buffer.data(ppf + 62);
    const auto *ppf_69 = buffer.data(ppf + 69);
    const auto *ppf_78 = buffer.data(ppf + 78);
    const auto *ppf_79 = buffer.data(ppf + 79);
    const auto *ppf_86 = buffer.data(ppf + 86);
    const auto *ppf_88 = buffer.data(ppf + 88);
    const auto *ppf_89 = buffer.data(ppf + 89);

    const auto *ppg1_46 = buffer.data(ppg1 + 46);
    const auto *ppg1_48 = buffer.data(ppg1 + 48);
    const auto *ppg1_55 = buffer.data(ppg1 + 55);
    const auto *ppg1_61 = buffer.data(ppg1 + 61);
    const auto *ppg1_63 = buffer.data(ppg1 + 63);
    const auto *ppg1_70 = buffer.data(ppg1 + 70);
    const auto *ppg1_90 = buffer.data(ppg1 + 90);
    const auto *ppg1_92 = buffer.data(ppg1 + 92);
    const auto *ppg1_94 = buffer.data(ppg1 + 94);
    const auto *ppg1_95 = buffer.data(ppg1 + 95);
    const auto *ppg1_104 = buffer.data(ppg1 + 104);
    const auto *ppg1_120 = buffer.data(ppg1 + 120);
    const auto *ppg1_122 = buffer.data(ppg1 + 122);
    const auto *ppg1_125 = buffer.data(ppg1 + 125);
    const auto *ppg1_130 = buffer.data(ppg1 + 130);
    const auto *ppg1_131 = buffer.data(ppg1 + 131);
    const auto *ppg1_132 = buffer.data(ppg1 + 132);
    const auto *ppg1_134 = buffer.data(ppg1 + 134);

    const auto *dsg0_45 = buffer.data(dsg0 + 45);
    const auto *dsg0_46 = buffer.data(dsg0 + 46);
    const auto *dsg0_48 = buffer.data(dsg0 + 48);
    const auto *dsg0_50 = buffer.data(dsg0 + 50);
    const auto *dsg0_55 = buffer.data(dsg0 + 55);
    const auto *dsg0_57 = buffer.data(dsg0 + 57);
    const auto *dsg0_59 = buffer.data(dsg0 + 59);
    const auto *dsg0_72 = buffer.data(dsg0 + 72);
    const auto *dsg0_75 = buffer.data(dsg0 + 75);
    const auto *dsg0_77 = buffer.data(dsg0 + 77);
    const auto *dsg0_78 = buffer.data(dsg0 + 78);
    const auto *dsg0_80 = buffer.data(dsg0 + 80);
    const auto *dsg0_85 = buffer.data(dsg0 + 85);
    const auto *dsg0_86 = buffer.data(dsg0 + 86);
    const auto *dsg0_87 = buffer.data(dsg0 + 87);
    const auto *dsg0_89 = buffer.data(dsg0 + 89);

    const auto *dsf_30 = buffer.data(dsf + 30);
    const auto *dsf_31 = buffer.data(dsf + 31);
    const auto *dsf_33 = buffer.data(dsf + 33);
    const auto *dsf_35 = buffer.data(dsf + 35);
    const auto *dsf_36 = buffer.data(dsf + 36);
    const auto *dsf_37 = buffer.data(dsf + 37);
    const auto *dsf_38 = buffer.data(dsf + 38);
    const auto *dsf_39 = buffer.data(dsf + 39);
    const auto *dsf_46 = buffer.data(dsf + 46);
    const auto *dsf_47 = buffer.data(dsf + 47);
    const auto *dsf_48 = buffer.data(dsf + 48);
    const auto *dsf_49 = buffer.data(dsf + 49);
    const auto *dsf_50 = buffer.data(dsf + 50);
    const auto *dsf_52 = buffer.data(dsf + 52);
    const auto *dsf_53 = buffer.data(dsf + 53);
    const auto *dsf_55 = buffer.data(dsf + 55);
    const auto *dsf_56 = buffer.data(dsf + 56);
    const auto *dsf_57 = buffer.data(dsf + 57);
    const auto *dsf_58 = buffer.data(dsf + 58);
    const auto *dsf_59 = buffer.data(dsf + 59);

    const auto *dsg1_45 = buffer.data(dsg1 + 45);
    const auto *dsg1_46 = buffer.data(dsg1 + 46);
    const auto *dsg1_48 = buffer.data(dsg1 + 48);
    const auto *dsg1_50 = buffer.data(dsg1 + 50);
    const auto *dsg1_55 = buffer.data(dsg1 + 55);
    const auto *dsg1_57 = buffer.data(dsg1 + 57);
    const auto *dsg1_59 = buffer.data(dsg1 + 59);
    const auto *dsg1_72 = buffer.data(dsg1 + 72);
    const auto *dsg1_75 = buffer.data(dsg1 + 75);
    const auto *dsg1_77 = buffer.data(dsg1 + 77);
    const auto *dsg1_78 = buffer.data(dsg1 + 78);
    const auto *dsg1_80 = buffer.data(dsg1 + 80);
    const auto *dsg1_85 = buffer.data(dsg1 + 85);
    const auto *dsg1_86 = buffer.data(dsg1 + 86);
    const auto *dsg1_87 = buffer.data(dsg1 + 87);
    const auto *dsg1_89 = buffer.data(dsg1 + 89);

    const auto *dpd0_60 = buffer.data(dpd0 + 60);
    const auto *dpd0_61 = buffer.data(dpd0 + 61);
    const auto *dpd0_63 = buffer.data(dpd0 + 63);
    const auto *dpd0_65 = buffer.data(dpd0 + 65);
    const auto *dpd0_71 = buffer.data(dpd0 + 71);
    const auto *dpd0_78 = buffer.data(dpd0 + 78);
    const auto *dpd0_80 = buffer.data(dpd0 + 80);
    const auto *dpd0_82 = buffer.data(dpd0 + 82);
    const auto *dpd0_83 = buffer.data(dpd0 + 83);
    const auto *dpd0_85 = buffer.data(dpd0 + 85);
    const auto *dpd0_87 = buffer.data(dpd0 + 87);
    const auto *dpd0_88 = buffer.data(dpd0 + 88);
    const auto *dpd0_99 = buffer.data(dpd0 + 99);
    const auto *dpd0_102 = buffer.data(dpd0 + 102);
    const auto *dpd0_104 = buffer.data(dpd0 + 104);
    const auto *dpd0_105 = buffer.data(dpd0 + 105);

    const auto *dpd1_60 = buffer.data(dpd1 + 60);
    const auto *dpd1_61 = buffer.data(dpd1 + 61);
    const auto *dpd1_63 = buffer.data(dpd1 + 63);
    const auto *dpd1_65 = buffer.data(dpd1 + 65);
    const auto *dpd1_71 = buffer.data(dpd1 + 71);
    const auto *dpd1_78 = buffer.data(dpd1 + 78);
    const auto *dpd1_80 = buffer.data(dpd1 + 80);
    const auto *dpd1_82 = buffer.data(dpd1 + 82);
    const auto *dpd1_83 = buffer.data(dpd1 + 83);
    const auto *dpd1_85 = buffer.data(dpd1 + 85);
    const auto *dpd1_87 = buffer.data(dpd1 + 87);
    const auto *dpd1_88 = buffer.data(dpd1 + 88);
    const auto *dpd1_99 = buffer.data(dpd1 + 99);
    const auto *dpd1_102 = buffer.data(dpd1 + 102);
    const auto *dpd1_104 = buffer.data(dpd1 + 104);
    const auto *dpd1_105 = buffer.data(dpd1 + 105);

    const auto *dpf_89 = buffer.data(dpf + 89);
    const auto *dpf_90 = buffer.data(dpf + 90);
    const auto *dpf_91 = buffer.data(dpf + 91);
    const auto *dpf_96 = buffer.data(dpf + 96);
    const auto *dpf_97 = buffer.data(dpf + 97);
    const auto *dpf_98 = buffer.data(dpf + 98);
    const auto *dpf_99 = buffer.data(dpf + 99);
    const auto *dpf_100 = buffer.data(dpf + 100);
    const auto *dpf_101 = buffer.data(dpf + 101);
    const auto *dpf_103 = buffer.data(dpf + 103);
    const auto *dpf_105 = buffer.data(dpf + 105);
    const auto *dpf_106 = buffer.data(dpf + 106);
    const auto *dpf_107 = buffer.data(dpf + 107);
    const auto *dpf_108 = buffer.data(dpf + 108);
    const auto *dpf_109 = buffer.data(dpf + 109);
    const auto *dpf_110 = buffer.data(dpf + 110);
    const auto *dpf_111 = buffer.data(dpf + 111);
    const auto *dpf_115 = buffer.data(dpf + 115);
    const auto *dpf_116 = buffer.data(dpf + 116);
    const auto *dpf_117 = buffer.data(dpf + 117);
    const auto *dpf_118 = buffer.data(dpf + 118);
    const auto *dpf_119 = buffer.data(dpf + 119);
    const auto *dpf_126 = buffer.data(dpf + 126);
    const auto *dpf_127 = buffer.data(dpf + 127);
    const auto *dpf_128 = buffer.data(dpf + 128);
    const auto *dpf_129 = buffer.data(dpf + 129);
    const auto *dpf_130 = buffer.data(dpf + 130);
    const auto *dpf_132 = buffer.data(dpf + 132);
    const auto *dpf_134 = buffer.data(dpf + 134);
    const auto *dpf_135 = buffer.data(dpf + 135);
    const auto *dpf_136 = buffer.data(dpf + 136);
    const auto *dpf_137 = buffer.data(dpf + 137);
    const auto *dpf_138 = buffer.data(dpf + 138);
    const auto *dpf_139 = buffer.data(dpf + 139);
    const auto *dpf_141 = buffer.data(dpf + 141);
    const auto *dpf_143 = buffer.data(dpf + 143);
    const auto *dpf_144 = buffer.data(dpf + 144);
    const auto *dpf_146 = buffer.data(dpf + 146);
    const auto *dpf_147 = buffer.data(dpf + 147);
    const auto *dpf_148 = buffer.data(dpf + 148);
    const auto *dpf_149 = buffer.data(dpf + 149);
    const auto *dpf_150 = buffer.data(dpf + 150);
    const auto *dpf_152 = buffer.data(dpf + 152);
    const auto *dpf_156 = buffer.data(dpf + 156);
    const auto *dpf_157 = buffer.data(dpf + 157);
    const auto *dpf_158 = buffer.data(dpf + 158);
    const auto *dpf_159 = buffer.data(dpf + 159);
    const auto *dpf_160 = buffer.data(dpf + 160);
    const auto *dpf_162 = buffer.data(dpf + 162);
    const auto *dpf_163 = buffer.data(dpf + 163);
    const auto *dpf_166 = buffer.data(dpf + 166);
    const auto *dpf_167 = buffer.data(dpf + 167);
    const auto *dpf_168 = buffer.data(dpf + 168);
    const auto *dpf_169 = buffer.data(dpf + 169);
    const auto *dpf_170 = buffer.data(dpf + 170);
    const auto *dpf_172 = buffer.data(dpf + 172);
    const auto *dpf_173 = buffer.data(dpf + 173);

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pa_x, pc_x, pc_y, ppg0_130, \
                         ppg0_131, ppg0_132, ppf_89, ppg1_130, ppg1_131, ppg1_132, \
                         dpf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_1 * ppf_89[k]
                   + f_4 * pc_x[k] * dpf_89[k];

        t_130[k] = pa_x[k] * ppg0_130[k]
                   - f_7 * pc_x[k] * ppg1_130[k];

        t_131[k] = pa_x[k] * ppg0_131[k]
                   - f_7 * pc_x[k] * ppg1_131[k];

        t_132[k] = pa_x[k] * ppg0_132[k]
                   - f_7 * pc_x[k] * ppg1_132[k];

        t_133[k] = f_4 * pc_y[k] * dpf_89[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pa_x, pb_x, pc_x, ppg0_134, ppg1_134, dsg0_45, \
                         dsg0_46, dsf_30, dsf_31, dsg1_45, dsg1_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_x[k] * ppg0_134[k]
                   - f_7 * pc_x[k] * ppg1_134[k];

        t_135[k] = pb_x[k] * dsg0_45[k]
                   + f_8 * dsf_30[k]
                   - f_7 * pc_x[k] * dsg1_45[k];

        t_136[k] = pb_x[k] * dsg0_46[k]
                   + f_11 * dsf_31[k]
                   - f_7 * pc_x[k] * dsg1_46[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pb_x, pc_x, pc_z, dsg0_48, dsg0_50, \
                         dsf_33, dsf_35, dsg1_48, dsg1_50, dpf_90, \
                         dpf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_4 * pc_z[k] * dpf_90[k];

        t_138[k] = pb_x[k] * dsg0_48[k]
                   + f_0 * dsf_33[k]
                   - f_7 * pc_x[k] * dsg1_48[k];

        t_139[k] = f_4 * pc_z[k] * dpf_91[k];

        t_140[k] = pb_x[k] * dsg0_50[k]
                   + f_0 * dsf_35[k]
                   - f_7 * pc_x[k] * dsg1_50[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pc_x, dsf_36, dsf_37, dsf_38, dsf_39, \
                         dpf_96, dpf_97, dpf_98, dpf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_1 * dsf_36[k]
                   + f_4 * pc_x[k] * dpf_96[k];

        t_142[k] = f_1 * dsf_37[k]
                   + f_4 * pc_x[k] * dpf_97[k];

        t_143[k] = f_1 * dsf_38[k]
                   + f_4 * pc_x[k] * dpf_98[k];

        t_144[k] = f_1 * dsf_39[k]
                   + f_4 * pc_x[k] * dpf_99[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pb_x, pc_x, pc_y, pc_z, ppf_39, dsg0_55, \
                         dsg0_57, dsg1_55, dsg1_57, dpf_96, dpf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pb_x[k] * dsg0_55[k]
                   - f_7 * pc_x[k] * dsg1_55[k];

        t_146[k] = f_4 * pc_z[k] * dpf_96[k];

        t_147[k] = pb_x[k] * dsg0_57[k]
                   - f_7 * pc_x[k] * dsg1_57[k];

        t_148[k] = f_0 * ppf_39[k]
                   + f_4 * pc_y[k] * dpf_99[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pb_x, pc_x, pc_z, dsg0_59, dsg1_59, \
                         dpd0_60, dpd0_61, dpd1_60, dpd1_61, dpf_100, \
                         dpf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pb_x[k] * dsg0_59[k]
                   - f_7 * pc_x[k] * dsg1_59[k];

        t_150[k] = f_2 * dpd0_60[k]
                   - f_3 * dpd1_60[k]
                   + f_4 * pc_x[k] * dpf_100[k];

        t_151[k] = f_9 * dpd0_61[k]
                   - f_10 * dpd1_61[k]
                   + f_4 * pc_x[k] * dpf_101[k];

        t_152[k] = f_4 * pc_z[k] * dpf_100[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pc_x, pc_z, dpd0_63, dpd0_65, \
                         dpd1_63, dpd1_65, dpf_101, dpf_103, dpf_105, dpf_106, \
                         dpf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_5 * dpd0_63[k]
                   - f_6 * dpd1_63[k]
                   + f_4 * pc_x[k] * dpf_103[k];

        t_154[k] = f_4 * pc_z[k] * dpf_101[k];

        t_155[k] = f_5 * dpd0_65[k]
                   - f_6 * dpd1_65[k]
                   + f_4 * pc_x[k] * dpf_105[k];

        t_156[k] = f_4 * pc_x[k] * dpf_106[k];

        t_157[k] = f_4 * pc_x[k] * dpf_107[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pc_x, pc_y, pc_z, ppf_46, dsf_36, \
                         dpd0_63, dpd1_63, dpf_106, dpf_107, dpf_108, \
                         dpf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_4 * pc_x[k] * dpf_108[k];

        t_159[k] = f_4 * pc_x[k] * dpf_109[k];

        t_160[k] = f_0 * ppf_46[k]
                   + f_1 * dsf_36[k]
                   + f_2 * dpd0_63[k]
                   - f_3 * dpd1_63[k]
                   + f_4 * pc_y[k] * dpf_106[k];

        t_161[k] = f_4 * pc_z[k] * dpf_106[k];

        t_162[k] = f_5 * dpd0_63[k]
                   - f_6 * dpd1_63[k]
                   + f_4 * pc_z[k] * dpf_107[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pb_z, pc_y, pc_z, ppf_49, dsg0_45, \
                         dsg0_46, dsf_39, dsg1_45, dsg1_46, dpd0_65, dpd1_65, \
                         dpf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_0 * ppf_49[k]
                   + f_1 * dsf_39[k]
                   + f_4 * pc_y[k] * dpf_109[k];

        t_164[k] = f_2 * dpd0_65[k]
                   - f_3 * dpd1_65[k]
                   + f_4 * pc_z[k] * dpf_109[k];

        t_165[k] = pb_z[k] * dsg0_45[k]
                   - f_7 * pc_z[k] * dsg1_45[k];

        t_166[k] = pb_z[k] * dsg0_46[k]
                   - f_7 * pc_z[k] * dsg1_46[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pb_z, pc_x, pc_z, dsg0_48, dsf_30, \
                         dsf_31, dsg1_48, dpd0_71, dpd1_71, dpf_110, dpf_111, \
                         dpf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_1 * dsf_30[k]
                   + f_4 * pc_z[k] * dpf_110[k];

        t_168[k] = pb_z[k] * dsg0_48[k]
                   - f_7 * pc_z[k] * dsg1_48[k];

        t_169[k] = f_1 * dsf_31[k]
                   + f_4 * pc_z[k] * dpf_111[k];

        t_170[k] = f_5 * dpd0_71[k]
                   - f_6 * dpd1_71[k]
                   + f_4 * pc_x[k] * dpf_115[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, t_175, t_176, pb_z, pc_x, pc_z, dsg0_55, \
                         dsf_36, dsg1_55, dpf_116, dpf_117, dpf_118, \
                         dpf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_4 * pc_x[k] * dpf_116[k];

        t_172[k] = f_4 * pc_x[k] * dpf_117[k];

        t_173[k] = f_4 * pc_x[k] * dpf_118[k];

        t_174[k] = f_4 * pc_x[k] * dpf_119[k];

        t_175[k] = pb_z[k] * dsg0_55[k]
                   - f_7 * pc_z[k] * dsg1_55[k];

        t_176[k] = f_1 * dsf_36[k]
                   + f_4 * pc_z[k] * dpf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_z, pc_y, pc_z, ppf_59, dsg0_57, dsg0_59, \
                         dsf_37, dsf_39, dsg1_57, dsg1_59, dpf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pb_z[k] * dsg0_57[k]
                   + f_0 * dsf_37[k]
                   - f_7 * pc_z[k] * dsg1_57[k];

        t_178[k] = f_0 * ppf_59[k]
                   + f_4 * pc_y[k] * dpf_119[k];

        t_179[k] = pb_z[k] * dsg0_59[k]
                   + f_8 * dsf_39[k]
                   - f_7 * pc_z[k] * dsg1_59[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pc_y, pc_z, ppg0_46, ppg0_48, \
                         ppg0_90, ppg0_92, ppg1_46, ppg1_48, ppg1_90, \
                         ppg1_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * ppg0_90[k]
                   - f_7 * pc_y[k] * ppg1_90[k];

        t_181[k] = pa_z[k] * ppg0_46[k]
                   - f_7 * pc_z[k] * ppg1_46[k];

        t_182[k] = pa_y[k] * ppg0_92[k]
                   - f_7 * pc_y[k] * ppg1_92[k];

        t_183[k] = pa_z[k] * ppg0_48[k]
                   - f_7 * pc_z[k] * ppg1_48[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pc_x, pc_y, ppg0_94, ppg0_95, \
                         ppf_62, ppg1_94, ppg1_95, dsf_46, dsf_47, dpf_126, \
                         dpf_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pa_y[k] * ppg0_94[k]
                   + f_1 * ppf_62[k]
                   - f_7 * pc_y[k] * ppg1_94[k];

        t_185[k] = pa_y[k] * ppg0_95[k]
                   - f_7 * pc_y[k] * ppg1_95[k];

        t_186[k] = f_1 * dsf_46[k]
                   + f_4 * pc_x[k] * dpf_126[k];

        t_187[k] = f_1 * dsf_47[k]
                   + f_4 * pc_x[k] * dpf_127[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_z, pc_x, pc_z, ppg0_55, ppf_36, \
                         ppg1_55, dsf_48, dsf_49, dpf_126, dpf_128, \
                         dpf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_1 * dsf_48[k]
                   + f_4 * pc_x[k] * dpf_128[k];

        t_189[k] = f_1 * dsf_49[k]
                   + f_4 * pc_x[k] * dpf_129[k];

        t_190[k] = pa_z[k] * ppg0_55[k]
                   - f_7 * pc_z[k] * ppg1_55[k];

        t_191[k] = f_1 * ppf_36[k]
                   + f_4 * pc_z[k] * dpf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_y, pb_x, pc_x, pc_y, ppg0_104, ppf_69, \
                         ppg1_104, dsg0_72, dsg1_72, dpf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pb_x[k] * dsg0_72[k]
                   - f_7 * pc_x[k] * dsg1_72[k];

        t_193[k] = f_1 * ppf_69[k]
                   + f_4 * pc_y[k] * dpf_129[k];

        t_194[k] = pa_y[k] * ppg0_104[k]
                   - f_7 * pc_y[k] * ppg1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_z, pc_x, pc_z, ppg0_61, ppg1_61, dpd0_78, \
                         dpd0_80, dpd1_78, dpd1_80, dpf_130, dpf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_2 * dpd0_78[k]
                   - f_3 * dpd1_78[k]
                   + f_4 * pc_x[k] * dpf_130[k];

        t_196[k] = pa_z[k] * ppg0_61[k]
                   - f_7 * pc_z[k] * ppg1_61[k];

        t_197[k] = f_9 * dpd0_80[k]
                   - f_10 * dpd1_80[k]
                   + f_4 * pc_x[k] * dpf_132[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pa_z, pc_x, pc_z, ppg0_63, ppg1_63, \
                         dpd0_82, dpd0_83, dpd1_82, dpd1_83, dpf_134, dpf_135, \
                         dpf_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_z[k] * ppg0_63[k]
                   - f_7 * pc_z[k] * ppg1_63[k];

        t_199[k] = f_5 * dpd0_82[k]
                   - f_6 * dpd1_82[k]
                   + f_4 * pc_x[k] * dpf_134[k];

        t_200[k] = f_5 * dpd0_83[k]
                   - f_6 * dpd1_83[k]
                   + f_4 * pc_x[k] * dpf_135[k];

        t_201[k] = f_4 * pc_x[k] * dpf_136[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pa_z, pc_x, pc_z, ppg0_70, ppf_46, \
                         ppg1_70, dpf_136, dpf_137, dpf_138, dpf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_4 * pc_x[k] * dpf_137[k];

        t_203[k] = f_4 * pc_x[k] * dpf_138[k];

        t_204[k] = f_4 * pc_x[k] * dpf_139[k];

        t_205[k] = pa_z[k] * ppg0_70[k]
                   - f_7 * pc_z[k] * ppg1_70[k];

        t_206[k] = f_1 * ppf_46[k]
                   + f_4 * pc_z[k] * dpf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pc_y, pc_z, ppf_49, ppf_78, ppf_79, dsf_48, \
                         dsf_49, dpd0_83, dpd1_83, dpf_138, dpf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_1 * ppf_78[k]
                   + f_1 * dsf_48[k]
                   + f_5 * dpd0_83[k]
                   - f_6 * dpd1_83[k]
                   + f_4 * pc_y[k] * dpf_138[k];

        t_208[k] = f_1 * ppf_79[k]
                   + f_1 * dsf_49[k]
                   + f_4 * pc_y[k] * dpf_139[k];

        t_209[k] = f_1 * ppf_49[k]
                   + f_2 * dpd0_83[k]
                   - f_3 * dpd1_83[k]
                   + f_4 * pc_z[k] * dpf_139[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pa_y, pc_x, pc_y, ppg0_120, ppg0_122, ppg1_120, \
                         ppg1_122, dpd0_85, dpd1_85, dpf_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = pa_y[k] * ppg0_120[k]
                   - f_7 * pc_y[k] * ppg1_120[k];

        t_211[k] = f_9 * dpd0_85[k]
                   - f_10 * dpd1_85[k]
                   + f_4 * pc_x[k] * dpf_141[k];

        t_212[k] = pa_y[k] * ppg0_122[k]
                   - f_7 * pc_y[k] * ppg1_122[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_y, pc_x, pc_y, ppg0_125, ppg1_125, \
                         dpd0_87, dpd0_88, dpd1_87, dpd1_88, dpf_143, dpf_144, \
                         dpf_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_5 * dpd0_87[k]
                   - f_6 * dpd1_87[k]
                   + f_4 * pc_x[k] * dpf_143[k];

        t_214[k] = f_5 * dpd0_88[k]
                   - f_6 * dpd1_88[k]
                   + f_4 * pc_x[k] * dpf_144[k];

        t_215[k] = pa_y[k] * ppg0_125[k]
                   - f_7 * pc_y[k] * ppg1_125[k];

        t_216[k] = f_4 * pc_x[k] * dpf_146[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pc_x, pc_y, ppf_86, dpd0_87, dpd1_87, \
                         dpf_146, dpf_147, dpf_148, dpf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_4 * pc_x[k] * dpf_147[k];

        t_218[k] = f_4 * pc_x[k] * dpf_148[k];

        t_219[k] = f_4 * pc_x[k] * dpf_149[k];

        t_220[k] = f_1 * ppf_86[k]
                   + f_2 * dpd0_87[k]
                   - f_3 * dpd1_87[k]
                   + f_4 * pc_y[k] * dpf_146[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pa_y, pc_y, pc_z, ppg0_132, ppf_56, ppf_88, \
                         ppf_89, ppg1_132, dsf_46, dpf_146, dpf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_1 * ppf_56[k]
                   + f_1 * dsf_46[k]
                   + f_4 * pc_z[k] * dpf_146[k];

        t_222[k] = pa_y[k] * ppg0_132[k]
                   + f_0 * ppf_88[k]
                   - f_7 * pc_y[k] * ppg1_132[k];

        t_223[k] = f_1 * ppf_89[k]
                   + f_4 * pc_y[k] * dpf_149[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, pa_y, pb_x, pc_x, pc_y, ppg0_134, ppg1_134, \
                         dsg0_75, dsf_50, dsg1_75, dpf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_y[k] * ppg0_134[k]
                   - f_7 * pc_y[k] * ppg1_134[k];

        t_225[k] = pb_x[k] * dsg0_75[k]
                   + f_8 * dsf_50[k]
                   - f_7 * pc_x[k] * dsg1_75[k];

        t_226[k] = f_4 * pc_y[k] * dpf_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pb_x, pc_x, pc_y, dsg0_77, dsg0_78, dsf_52, \
                         dsf_53, dsg1_77, dsg1_78, dpf_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pb_x[k] * dsg0_77[k]
                   + f_11 * dsf_52[k]
                   - f_7 * pc_x[k] * dsg1_77[k];

        t_228[k] = pb_x[k] * dsg0_78[k]
                   + f_0 * dsf_53[k]
                   - f_7 * pc_x[k] * dsg1_78[k];

        t_229[k] = f_4 * pc_y[k] * dpf_152[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pb_x, pc_x, dsg0_80, dsf_55, dsf_56, \
                         dsf_57, dsf_58, dsg1_80, dpf_156, dpf_157, \
                         dpf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = pb_x[k] * dsg0_80[k]
                   + f_0 * dsf_55[k]
                   - f_7 * pc_x[k] * dsg1_80[k];

        t_231[k] = f_1 * dsf_56[k]
                   + f_4 * pc_x[k] * dpf_156[k];

        t_232[k] = f_1 * dsf_57[k]
                   + f_4 * pc_x[k] * dpf_157[k];

        t_233[k] = f_1 * dsf_58[k]
                   + f_4 * pc_x[k] * dpf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pb_x, pc_x, pc_y, dsg0_85, \
                         dsg0_86, dsg0_87, dsf_59, dsg1_85, dsg1_86, dsg1_87, \
                         dpf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_1 * dsf_59[k]
                   + f_4 * pc_x[k] * dpf_159[k];

        t_235[k] = pb_x[k] * dsg0_85[k]
                   - f_7 * pc_x[k] * dsg1_85[k];

        t_236[k] = pb_x[k] * dsg0_86[k]
                   - f_7 * pc_x[k] * dsg1_86[k];

        t_237[k] = pb_x[k] * dsg0_87[k]
                   - f_7 * pc_x[k] * dsg1_87[k];

        t_238[k] = f_4 * pc_y[k] * dpf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pb_x, pb_y, pc_x, pc_y, dsg0_75, dsg0_77, \
                         dsg0_89, dsf_50, dsg1_75, dsg1_77, dsg1_89, \
                         dpf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = pb_x[k] * dsg0_89[k]
                   - f_7 * pc_x[k] * dsg1_89[k];

        t_240[k] = pb_y[k] * dsg0_75[k]
                   - f_7 * pc_y[k] * dsg1_75[k];

        t_241[k] = f_1 * dsf_50[k]
                   + f_4 * pc_y[k] * dpf_160[k];

        t_242[k] = pb_y[k] * dsg0_77[k]
                   - f_7 * pc_y[k] * dsg1_77[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_y, pc_x, pc_y, dsg0_80, dsf_52, \
                         dsg1_80, dpd0_99, dpd1_99, dpf_162, dpf_163, \
                         dpf_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_5 * dpd0_99[k]
                   - f_6 * dpd1_99[k]
                   + f_4 * pc_x[k] * dpf_163[k];

        t_244[k] = f_1 * dsf_52[k]
                   + f_4 * pc_y[k] * dpf_162[k];

        t_245[k] = pb_y[k] * dsg0_80[k]
                   - f_7 * pc_y[k] * dsg1_80[k];

        t_246[k] = f_4 * pc_x[k] * dpf_166[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pb_y, pc_x, pc_y, dsg0_85, dsf_56, \
                         dsg1_85, dpf_167, dpf_168, dpf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_4 * pc_x[k] * dpf_167[k];

        t_248[k] = f_4 * pc_x[k] * dpf_168[k];

        t_249[k] = f_4 * pc_x[k] * dpf_169[k];

        t_250[k] = pb_y[k] * dsg0_85[k]
                   + f_8 * dsf_56[k]
                   - f_7 * pc_y[k] * dsg1_85[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_y, dsg0_86, dsg0_87, dsg0_89, \
                         dsf_57, dsf_58, dsf_59, dsg1_86, dsg1_87, dsg1_89, \
                         dpf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * dsg0_86[k]
                   + f_11 * dsf_57[k]
                   - f_7 * pc_y[k] * dsg1_86[k];

        t_252[k] = pb_y[k] * dsg0_87[k]
                   + f_0 * dsf_58[k]
                   - f_7 * pc_y[k] * dsg1_87[k];

        t_253[k] = f_1 * dsf_59[k]
                   + f_4 * pc_y[k] * dpf_169[k];

        t_254[k] = pb_y[k] * dsg0_89[k]
                   - f_7 * pc_y[k] * dsg1_89[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pc_x, pc_y, dpd0_102, dpd0_104, \
                         dpd0_105, dpd1_102, dpd1_104, dpd1_105, dpf_170, dpf_172, \
                         dpf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_2 * dpd0_102[k]
                   - f_3 * dpd1_102[k]
                   + f_4 * pc_x[k] * dpf_170[k];

        t_256[k] = f_4 * pc_y[k] * dpf_170[k];

        t_257[k] = f_9 * dpd0_104[k]
                   - f_10 * dpd1_104[k]
                   + f_4 * pc_x[k] * dpf_172[k];

        t_258[k] = f_5 * dpd0_105[k]
                   - f_6 * dpd1_105[k]
                   + f_4 * pc_x[k] * dpf_173[k];

        t_259[k] = f_4 * pc_y[k] * dpf_172[k];
    }
}

static auto
compute_prim_dpg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t ppf, const size_t dsf,
                                                          const size_t dpd0, const size_t dpd1,
                                                          const size_t dpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_9 = 1.0 / gamma;
    const auto f_10 = p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppf_89 = buffer.data(ppf + 89);

    const auto *dsf_59 = buffer.data(dsf + 59);

    const auto *dpd0_105 = buffer.data(dpd0 + 105);
    const auto *dpd0_106 = buffer.data(dpd0 + 106);
    const auto *dpd0_107 = buffer.data(dpd0 + 107);

    const auto *dpd1_105 = buffer.data(dpd1 + 105);
    const auto *dpd1_106 = buffer.data(dpd1 + 106);
    const auto *dpd1_107 = buffer.data(dpd1 + 107);

    const auto *dpf_175 = buffer.data(dpf + 175);
    const auto *dpf_176 = buffer.data(dpf + 176);
    const auto *dpf_177 = buffer.data(dpf + 177);
    const auto *dpf_178 = buffer.data(dpf + 178);
    const auto *dpf_179 = buffer.data(dpf + 179);

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, pc_x, dpd0_107, dpd1_107, dpf_175, \
                         dpf_176, dpf_177, dpf_178, dpf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_5 * dpd0_107[k]
                   - f_6 * dpd1_107[k]
                   + f_4 * pc_x[k] * dpf_175[k];

        t_261[k] = f_4 * pc_x[k] * dpf_176[k];

        t_262[k] = f_4 * pc_x[k] * dpf_177[k];

        t_263[k] = f_4 * pc_x[k] * dpf_178[k];

        t_264[k] = f_4 * pc_x[k] * dpf_179[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pc_y, dpd0_105, dpd0_106, dpd0_107, \
                         dpd1_105, dpd1_106, dpd1_107, dpf_176, dpf_177, dpf_178, \
                         dpf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_2 * dpd0_105[k]
                   - f_3 * dpd1_105[k]
                   + f_4 * pc_y[k] * dpf_176[k];

        t_266[k] = f_9 * dpd0_106[k]
                   - f_10 * dpd1_106[k]
                   + f_4 * pc_y[k] * dpf_177[k];

        t_267[k] = f_5 * dpd0_107[k]
                   - f_6 * dpd1_107[k]
                   + f_4 * pc_y[k] * dpf_178[k];

        t_268[k] = f_4 * pc_y[k] * dpf_179[k];
    }

#pragma omp simd aligned(t_269, pc_z, ppf_89, dsf_59, dpd0_107, dpd1_107, \
                         dpf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_0 * ppf_89[k]
                   + f_1 * dsf_59[k]
                   + f_2 * dpd0_107[k]
                   - f_3 * dpd1_107[k]
                   + f_4 * pc_z[k] * dpf_179[k];
    }
}

auto
compute_prim_dpg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t ppg0,
                                                   const size_t ppf, const size_t ppg1,
                                                   const size_t dsg0, const size_t dsf,
                                                   const size_t dsg1, const size_t dpd0,
                                                   const size_t dpd1, const size_t dpf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_dpg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, ppg0,
                                                              ppf, ppg1, dsg0, dsf, dsg1, dpd0,
                                                              dpd1, dpf, ncols, gamma, p, q);

    compute_prim_dpg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, ppg0,
                                                              ppf, ppg1, dsg0, dsf, dsg1, dpd0,
                                                              dpd1, dpf, ncols, gamma, p, q);

    compute_prim_dpg_three_center_electron_repulsion_0_piece2(buffer, target, pc, ppf, dsf,
                                                              dpd0, dpd1, dpf, ncols, gamma, p,
                                                              q);
}

}  // namespace simdt3ceri
