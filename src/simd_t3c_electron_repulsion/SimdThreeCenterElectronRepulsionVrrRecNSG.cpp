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


#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_nsg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msg0,
                                                          const size_t msf, const size_t msg1,
                                                          const size_t nsd0, const size_t nsd1,
                                                          const size_t nsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 1.5 / q;

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

    const auto *msg0_0 = buffer.data(msg0 + 0);
    const auto *msg0_3 = buffer.data(msg0 + 3);
    const auto *msg0_5 = buffer.data(msg0 + 5);
    const auto *msg0_10 = buffer.data(msg0 + 10);
    const auto *msg0_14 = buffer.data(msg0 + 14);
    const auto *msg0_18 = buffer.data(msg0 + 18);
    const auto *msg0_25 = buffer.data(msg0 + 25);
    const auto *msg0_30 = buffer.data(msg0 + 30);
    const auto *msg0_35 = buffer.data(msg0 + 35);
    const auto *msg0_44 = buffer.data(msg0 + 44);
    const auto *msg0_45 = buffer.data(msg0 + 45);
    const auto *msg0_48 = buffer.data(msg0 + 48);
    const auto *msg0_55 = buffer.data(msg0 + 55);
    const auto *msg0_75 = buffer.data(msg0 + 75);
    const auto *msg0_78 = buffer.data(msg0 + 78);
    const auto *msg0_80 = buffer.data(msg0 + 80);

    const auto *msf_0 = buffer.data(msf + 0);
    const auto *msf_1 = buffer.data(msf + 1);
    const auto *msf_2 = buffer.data(msf + 2);
    const auto *msf_6 = buffer.data(msf + 6);
    const auto *msf_9 = buffer.data(msf + 9);
    const auto *msf_10 = buffer.data(msf + 10);
    const auto *msf_16 = buffer.data(msf + 16);
    const auto *msf_18 = buffer.data(msf + 18);
    const auto *msf_19 = buffer.data(msf + 19);
    const auto *msf_20 = buffer.data(msf + 20);
    const auto *msf_22 = buffer.data(msf + 22);
    const auto *msf_26 = buffer.data(msf + 26);
    const auto *msf_27 = buffer.data(msf + 27);
    const auto *msf_28 = buffer.data(msf + 28);
    const auto *msf_29 = buffer.data(msf + 29);
    const auto *msf_30 = buffer.data(msf + 30);
    const auto *msf_33 = buffer.data(msf + 33);
    const auto *msf_36 = buffer.data(msf + 36);
    const auto *msf_38 = buffer.data(msf + 38);
    const auto *msf_39 = buffer.data(msf + 39);
    const auto *msf_40 = buffer.data(msf + 40);
    const auto *msf_42 = buffer.data(msf + 42);
    const auto *msf_46 = buffer.data(msf + 46);
    const auto *msf_47 = buffer.data(msf + 47);
    const auto *msf_48 = buffer.data(msf + 48);
    const auto *msf_49 = buffer.data(msf + 49);
    const auto *msf_50 = buffer.data(msf + 50);
    const auto *msf_51 = buffer.data(msf + 51);
    const auto *msf_52 = buffer.data(msf + 52);
    const auto *msf_55 = buffer.data(msf + 55);
    const auto *msf_56 = buffer.data(msf + 56);
    const auto *msf_57 = buffer.data(msf + 57);
    const auto *msf_59 = buffer.data(msf + 59);
    const auto *msf_60 = buffer.data(msf + 60);
    const auto *msf_63 = buffer.data(msf + 63);
    const auto *msf_66 = buffer.data(msf + 66);
    const auto *msf_68 = buffer.data(msf + 68);
    const auto *msf_69 = buffer.data(msf + 69);
    const auto *msf_75 = buffer.data(msf + 75);
    const auto *msf_76 = buffer.data(msf + 76);
    const auto *msf_77 = buffer.data(msf + 77);
    const auto *msf_78 = buffer.data(msf + 78);
    const auto *msf_79 = buffer.data(msf + 79);
    const auto *msf_86 = buffer.data(msf + 86);
    const auto *msf_87 = buffer.data(msf + 87);
    const auto *msf_88 = buffer.data(msf + 88);
    const auto *msf_89 = buffer.data(msf + 89);

    const auto *msg1_0 = buffer.data(msg1 + 0);
    const auto *msg1_3 = buffer.data(msg1 + 3);
    const auto *msg1_5 = buffer.data(msg1 + 5);
    const auto *msg1_10 = buffer.data(msg1 + 10);
    const auto *msg1_14 = buffer.data(msg1 + 14);
    const auto *msg1_18 = buffer.data(msg1 + 18);
    const auto *msg1_25 = buffer.data(msg1 + 25);
    const auto *msg1_30 = buffer.data(msg1 + 30);
    const auto *msg1_35 = buffer.data(msg1 + 35);
    const auto *msg1_44 = buffer.data(msg1 + 44);
    const auto *msg1_45 = buffer.data(msg1 + 45);
    const auto *msg1_48 = buffer.data(msg1 + 48);
    const auto *msg1_55 = buffer.data(msg1 + 55);
    const auto *msg1_75 = buffer.data(msg1 + 75);
    const auto *msg1_78 = buffer.data(msg1 + 78);
    const auto *msg1_80 = buffer.data(msg1 + 80);

    const auto *nsd0_0 = buffer.data(nsd0 + 0);
    const auto *nsd0_3 = buffer.data(nsd0 + 3);
    const auto *nsd0_5 = buffer.data(nsd0 + 5);
    const auto *nsd0_9 = buffer.data(nsd0 + 9);
    const auto *nsd0_16 = buffer.data(nsd0 + 16);
    const auto *nsd0_17 = buffer.data(nsd0 + 17);
    const auto *nsd0_18 = buffer.data(nsd0 + 18);
    const auto *nsd0_21 = buffer.data(nsd0 + 21);
    const auto *nsd0_23 = buffer.data(nsd0 + 23);
    const auto *nsd0_29 = buffer.data(nsd0 + 29);
    const auto *nsd0_30 = buffer.data(nsd0 + 30);
    const auto *nsd0_33 = buffer.data(nsd0 + 33);
    const auto *nsd0_34 = buffer.data(nsd0 + 34);
    const auto *nsd0_35 = buffer.data(nsd0 + 35);
    const auto *nsd0_36 = buffer.data(nsd0 + 36);
    const auto *nsd0_39 = buffer.data(nsd0 + 39);
    const auto *nsd0_41 = buffer.data(nsd0 + 41);
    const auto *nsd0_47 = buffer.data(nsd0 + 47);

    const auto *nsd1_0 = buffer.data(nsd1 + 0);
    const auto *nsd1_3 = buffer.data(nsd1 + 3);
    const auto *nsd1_5 = buffer.data(nsd1 + 5);
    const auto *nsd1_9 = buffer.data(nsd1 + 9);
    const auto *nsd1_16 = buffer.data(nsd1 + 16);
    const auto *nsd1_17 = buffer.data(nsd1 + 17);
    const auto *nsd1_18 = buffer.data(nsd1 + 18);
    const auto *nsd1_21 = buffer.data(nsd1 + 21);
    const auto *nsd1_23 = buffer.data(nsd1 + 23);
    const auto *nsd1_29 = buffer.data(nsd1 + 29);
    const auto *nsd1_30 = buffer.data(nsd1 + 30);
    const auto *nsd1_33 = buffer.data(nsd1 + 33);
    const auto *nsd1_34 = buffer.data(nsd1 + 34);
    const auto *nsd1_35 = buffer.data(nsd1 + 35);
    const auto *nsd1_36 = buffer.data(nsd1 + 36);
    const auto *nsd1_39 = buffer.data(nsd1 + 39);
    const auto *nsd1_41 = buffer.data(nsd1 + 41);
    const auto *nsd1_47 = buffer.data(nsd1 + 47);

    const auto *nsf_0 = buffer.data(nsf + 0);
    const auto *nsf_1 = buffer.data(nsf + 1);
    const auto *nsf_2 = buffer.data(nsf + 2);
    const auto *nsf_3 = buffer.data(nsf + 3);
    const auto *nsf_5 = buffer.data(nsf + 5);
    const auto *nsf_6 = buffer.data(nsf + 6);
    const auto *nsf_8 = buffer.data(nsf + 8);
    const auto *nsf_9 = buffer.data(nsf + 9);
    const auto *nsf_10 = buffer.data(nsf + 10);
    const auto *nsf_11 = buffer.data(nsf + 11);
    const auto *nsf_13 = buffer.data(nsf + 13);
    const auto *nsf_16 = buffer.data(nsf + 16);
    const auto *nsf_17 = buffer.data(nsf + 17);
    const auto *nsf_18 = buffer.data(nsf + 18);
    const auto *nsf_19 = buffer.data(nsf + 19);
    const auto *nsf_20 = buffer.data(nsf + 20);
    const auto *nsf_22 = buffer.data(nsf + 22);
    const auto *nsf_25 = buffer.data(nsf + 25);
    const auto *nsf_26 = buffer.data(nsf + 26);
    const auto *nsf_27 = buffer.data(nsf + 27);
    const auto *nsf_28 = buffer.data(nsf + 28);
    const auto *nsf_29 = buffer.data(nsf + 29);
    const auto *nsf_30 = buffer.data(nsf + 30);
    const auto *nsf_31 = buffer.data(nsf + 31);
    const auto *nsf_32 = buffer.data(nsf + 32);
    const auto *nsf_33 = buffer.data(nsf + 33);
    const auto *nsf_36 = buffer.data(nsf + 36);
    const auto *nsf_37 = buffer.data(nsf + 37);
    const auto *nsf_38 = buffer.data(nsf + 38);
    const auto *nsf_39 = buffer.data(nsf + 39);
    const auto *nsf_40 = buffer.data(nsf + 40);
    const auto *nsf_42 = buffer.data(nsf + 42);
    const auto *nsf_46 = buffer.data(nsf + 46);
    const auto *nsf_47 = buffer.data(nsf + 47);
    const auto *nsf_48 = buffer.data(nsf + 48);
    const auto *nsf_49 = buffer.data(nsf + 49);
    const auto *nsf_50 = buffer.data(nsf + 50);
    const auto *nsf_51 = buffer.data(nsf + 51);
    const auto *nsf_52 = buffer.data(nsf + 52);
    const auto *nsf_55 = buffer.data(nsf + 55);
    const auto *nsf_56 = buffer.data(nsf + 56);
    const auto *nsf_57 = buffer.data(nsf + 57);
    const auto *nsf_58 = buffer.data(nsf + 58);
    const auto *nsf_59 = buffer.data(nsf + 59);
    const auto *nsf_60 = buffer.data(nsf + 60);
    const auto *nsf_61 = buffer.data(nsf + 61);
    const auto *nsf_62 = buffer.data(nsf + 62);
    const auto *nsf_63 = buffer.data(nsf + 63);
    const auto *nsf_66 = buffer.data(nsf + 66);
    const auto *nsf_67 = buffer.data(nsf + 67);
    const auto *nsf_68 = buffer.data(nsf + 68);
    const auto *nsf_69 = buffer.data(nsf + 69);
    const auto *nsf_70 = buffer.data(nsf + 70);
    const auto *nsf_72 = buffer.data(nsf + 72);
    const auto *nsf_75 = buffer.data(nsf + 75);
    const auto *nsf_76 = buffer.data(nsf + 76);
    const auto *nsf_77 = buffer.data(nsf + 77);
    const auto *nsf_78 = buffer.data(nsf + 78);
    const auto *nsf_79 = buffer.data(nsf + 79);
    const auto *nsf_80 = buffer.data(nsf + 80);
    const auto *nsf_82 = buffer.data(nsf + 82);
    const auto *nsf_86 = buffer.data(nsf + 86);
    const auto *nsf_87 = buffer.data(nsf + 87);
    const auto *nsf_88 = buffer.data(nsf + 88);
    const auto *nsf_89 = buffer.data(nsf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, msf_0, nsd0_0, \
                         nsd1_0, nsf_0, nsf_1, nsf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * msf_0[k]
                 + f_1 * nsd0_0[k]
                 - f_2 * nsd1_0[k]
                 + f_3 * pc_x[k] * nsf_0[k];

        t_1[k] = f_3 * pc_y[k] * nsf_0[k];

        t_2[k] = f_3 * pc_z[k] * nsf_0[k];

        t_3[k] = f_4 * nsd0_0[k]
                 - f_5 * nsd1_0[k]
                 + f_3 * pc_y[k] * nsf_1[k];

        t_4[k] = f_3 * pc_y[k] * nsf_2[k];

        t_5[k] = f_4 * nsd0_0[k]
                 - f_5 * nsd1_0[k]
                 + f_3 * pc_z[k] * nsf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, msf_6, msf_9, nsd0_3, \
                         nsd1_3, nsf_3, nsf_5, nsf_6, nsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * msf_6[k]
                 + f_3 * pc_x[k] * nsf_6[k];

        t_7[k] = f_3 * pc_z[k] * nsf_3[k];

        t_8[k] = f_3 * pc_y[k] * nsf_5[k];

        t_9[k] = f_0 * msf_9[k]
                 + f_3 * pc_x[k] * nsf_9[k];

        t_10[k] = f_1 * nsd0_3[k]
                  - f_2 * nsd1_3[k]
                  + f_3 * pc_y[k] * nsf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, msg0_0, msg1_0, \
                         nsd0_5, nsd1_5, nsf_6, nsf_8, nsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * nsf_6[k];

        t_12[k] = f_4 * nsd0_5[k]
                  - f_5 * nsd1_5[k]
                  + f_3 * pc_y[k] * nsf_8[k];

        t_13[k] = f_3 * pc_y[k] * nsf_9[k];

        t_14[k] = f_1 * nsd0_5[k]
                  - f_2 * nsd1_5[k]
                  + f_3 * pc_z[k] * nsf_9[k];

        t_15[k] = pa_y[k] * msg0_0[k]
                  - f_6 * pc_y[k] * msg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, msg0_3, msg0_5, \
                         msf_0, msf_1, msg1_3, msg1_5, nsf_10, nsf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * msf_0[k]
                  + f_3 * pc_y[k] * nsf_10[k];

        t_17[k] = f_3 * pc_z[k] * nsf_10[k];

        t_18[k] = pa_y[k] * msg0_3[k]
                  + f_8 * msf_1[k]
                  - f_6 * pc_y[k] * msg1_3[k];

        t_19[k] = f_3 * pc_z[k] * nsf_11[k];

        t_20[k] = pa_y[k] * msg0_5[k]
                  - f_6 * pc_y[k] * msg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, msf_16, msf_18, msf_19, nsf_13, \
                         nsf_16, nsf_18, nsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * msf_16[k]
                  + f_3 * pc_x[k] * nsf_16[k];

        t_22[k] = f_3 * pc_z[k] * nsf_13[k];

        t_23[k] = f_9 * msf_18[k]
                  + f_3 * pc_x[k] * nsf_18[k];

        t_24[k] = f_9 * msf_19[k]
                  + f_3 * pc_x[k] * nsf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, msf_6, msf_9, nsd0_9, nsd1_9, \
                         nsf_16, nsf_17, nsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * msf_6[k]
                  + f_1 * nsd0_9[k]
                  - f_2 * nsd1_9[k]
                  + f_3 * pc_y[k] * nsf_16[k];

        t_26[k] = f_3 * pc_z[k] * nsf_16[k];

        t_27[k] = f_4 * nsd0_9[k]
                  - f_5 * nsd1_9[k]
                  + f_3 * pc_z[k] * nsf_17[k];

        t_28[k] = f_7 * msf_9[k]
                  + f_3 * pc_y[k] * nsf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, msg0_0, msg0_14, \
                         msf_0, msg1_0, msg1_14, nsf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * msg0_14[k]
                  - f_6 * pc_y[k] * msg1_14[k];

        t_30[k] = pa_z[k] * msg0_0[k]
                  - f_6 * pc_z[k] * msg1_0[k];

        t_31[k] = f_3 * pc_y[k] * nsf_20[k];

        t_32[k] = f_7 * msf_0[k]
                  + f_3 * pc_z[k] * nsf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, msg0_3, msg0_5, \
                         msf_2, msf_26, msg1_3, msg1_5, nsf_22, \
                         nsf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * msg0_3[k]
                  - f_6 * pc_z[k] * msg1_3[k];

        t_34[k] = f_3 * pc_y[k] * nsf_22[k];

        t_35[k] = pa_z[k] * msg0_5[k]
                  + f_8 * msf_2[k]
                  - f_6 * pc_z[k] * msg1_5[k];

        t_36[k] = f_9 * msf_26[k]
                  + f_3 * pc_x[k] * nsf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, msg0_10, msf_27, \
                         msf_29, msg1_10, nsf_25, nsf_27, nsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * msf_27[k]
                  + f_3 * pc_x[k] * nsf_27[k];

        t_38[k] = f_3 * pc_y[k] * nsf_25[k];

        t_39[k] = f_9 * msf_29[k]
                  + f_3 * pc_x[k] * nsf_29[k];

        t_40[k] = pa_z[k] * msg0_10[k]
                  - f_6 * pc_z[k] * msg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, msf_9, nsd0_16, nsd0_17, nsd1_16, \
                         nsd1_17, nsf_27, nsf_28, nsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * nsd0_16[k]
                  - f_11 * nsd1_16[k]
                  + f_3 * pc_y[k] * nsf_27[k];

        t_42[k] = f_4 * nsd0_17[k]
                  - f_5 * nsd1_17[k]
                  + f_3 * pc_y[k] * nsf_28[k];

        t_43[k] = f_3 * pc_y[k] * nsf_29[k];

        t_44[k] = f_7 * msf_9[k]
                  + f_1 * nsd0_17[k]
                  - f_2 * nsd1_17[k]
                  + f_3 * pc_z[k] * nsf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, msf_10, msf_30, msf_33, \
                         nsd0_18, nsd0_21, nsd1_18, nsd1_21, nsf_30, \
                         nsf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * msf_30[k]
                  + f_1 * nsd0_18[k]
                  - f_2 * nsd1_18[k]
                  + f_3 * pc_x[k] * nsf_30[k];

        t_46[k] = f_8 * msf_10[k]
                  + f_3 * pc_y[k] * nsf_30[k];

        t_47[k] = f_3 * pc_z[k] * nsf_30[k];

        t_48[k] = f_12 * msf_33[k]
                  + f_4 * nsd0_21[k]
                  - f_5 * nsd1_21[k]
                  + f_3 * pc_x[k] * nsf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, msf_36, msf_38, nsd0_18, \
                         nsd1_18, nsf_31, nsf_32, nsf_33, nsf_36, \
                         nsf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * nsf_31[k];

        t_50[k] = f_4 * nsd0_18[k]
                  - f_5 * nsd1_18[k]
                  + f_3 * pc_z[k] * nsf_32[k];

        t_51[k] = f_12 * msf_36[k]
                  + f_3 * pc_x[k] * nsf_36[k];

        t_52[k] = f_3 * pc_z[k] * nsf_33[k];

        t_53[k] = f_12 * msf_38[k]
                  + f_3 * pc_x[k] * nsf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, msf_16, msf_19, \
                         msf_39, nsd0_21, nsd1_21, nsf_36, nsf_37, \
                         nsf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * msf_39[k]
                  + f_3 * pc_x[k] * nsf_39[k];

        t_55[k] = f_8 * msf_16[k]
                  + f_1 * nsd0_21[k]
                  - f_2 * nsd1_21[k]
                  + f_3 * pc_y[k] * nsf_36[k];

        t_56[k] = f_3 * pc_z[k] * nsf_36[k];

        t_57[k] = f_4 * nsd0_21[k]
                  - f_5 * nsd1_21[k]
                  + f_3 * pc_z[k] * nsf_37[k];

        t_58[k] = f_8 * msf_19[k]
                  + f_3 * pc_y[k] * nsf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, msg0_30, msf_10, msf_20, \
                         msg1_30, nsd0_23, nsd1_23, nsf_39, nsf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * nsd0_23[k]
                  - f_2 * nsd1_23[k]
                  + f_3 * pc_z[k] * nsf_39[k];

        t_60[k] = pa_y[k] * msg0_30[k]
                  - f_6 * pc_y[k] * msg1_30[k];

        t_61[k] = f_7 * msf_20[k]
                  + f_3 * pc_y[k] * nsf_40[k];

        t_62[k] = f_7 * msf_10[k]
                  + f_3 * pc_z[k] * nsf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, msg0_18, msg0_35, msf_22, \
                         msg1_18, msg1_35, nsf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * msg0_18[k]
                  - f_6 * pc_z[k] * msg1_18[k];

        t_64[k] = f_7 * msf_22[k]
                  + f_3 * pc_y[k] * nsf_42[k];

        t_65[k] = pa_y[k] * msg0_35[k]
                  - f_6 * pc_y[k] * msg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, msf_46, msf_47, msf_48, msf_49, nsf_46, \
                         nsf_47, nsf_48, nsf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * msf_46[k]
                  + f_3 * pc_x[k] * nsf_46[k];

        t_67[k] = f_12 * msf_47[k]
                  + f_3 * pc_x[k] * nsf_47[k];

        t_68[k] = f_12 * msf_48[k]
                  + f_3 * pc_x[k] * nsf_48[k];

        t_69[k] = f_12 * msf_49[k]
                  + f_3 * pc_x[k] * nsf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pc_y, pc_z, msg0_25, msf_16, msf_28, msg1_25, \
                         nsd0_29, nsd1_29, nsf_46, nsf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * msg0_25[k]
                  - f_6 * pc_z[k] * msg1_25[k];

        t_71[k] = f_7 * msf_16[k]
                  + f_3 * pc_z[k] * nsf_46[k];

        t_72[k] = f_7 * msf_28[k]
                  + f_4 * nsd0_29[k]
                  - f_5 * nsd1_29[k]
                  + f_3 * pc_y[k] * nsf_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, msg0_44, msf_29, msf_50, \
                         msg1_44, nsd0_30, nsd1_30, nsf_49, nsf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * msf_29[k]
                  + f_3 * pc_y[k] * nsf_49[k];

        t_74[k] = pa_y[k] * msg0_44[k]
                  - f_6 * pc_y[k] * msg1_44[k];

        t_75[k] = f_12 * msf_50[k]
                  + f_1 * nsd0_30[k]
                  - f_2 * nsd1_30[k]
                  + f_3 * pc_x[k] * nsf_50[k];

        t_76[k] = f_3 * pc_y[k] * nsf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, msf_20, nsd0_30, nsd1_30, nsf_50, \
                         nsf_51, nsf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * msf_20[k]
                  + f_3 * pc_z[k] * nsf_50[k];

        t_78[k] = f_4 * nsd0_30[k]
                  - f_5 * nsd1_30[k]
                  + f_3 * pc_y[k] * nsf_51[k];

        t_79[k] = f_3 * pc_y[k] * nsf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, msf_55, msf_56, msf_57, nsd0_35, \
                         nsd1_35, nsf_55, nsf_56, nsf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_12 * msf_55[k]
                  + f_4 * nsd0_35[k]
                  - f_5 * nsd1_35[k]
                  + f_3 * pc_x[k] * nsf_55[k];

        t_81[k] = f_12 * msf_56[k]
                  + f_3 * pc_x[k] * nsf_56[k];

        t_82[k] = f_12 * msf_57[k]
                  + f_3 * pc_x[k] * nsf_57[k];

        t_83[k] = f_3 * pc_y[k] * nsf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, msf_59, nsd0_33, nsd0_34, nsd1_33, \
                         nsd1_34, nsf_56, nsf_57, nsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * msf_59[k]
                  + f_3 * pc_x[k] * nsf_59[k];

        t_85[k] = f_1 * nsd0_33[k]
                  - f_2 * nsd1_33[k]
                  + f_3 * pc_y[k] * nsf_56[k];

        t_86[k] = f_10 * nsd0_34[k]
                  - f_11 * nsd1_34[k]
                  + f_3 * pc_y[k] * nsf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, msf_29, msf_60, nsd0_35, \
                         nsd0_36, nsd1_35, nsd1_36, nsf_58, nsf_59, \
                         nsf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * nsd0_35[k]
                  - f_5 * nsd1_35[k]
                  + f_3 * pc_y[k] * nsf_58[k];

        t_88[k] = f_3 * pc_y[k] * nsf_59[k];

        t_89[k] = f_8 * msf_29[k]
                  + f_1 * nsd0_35[k]
                  - f_2 * nsd1_35[k]
                  + f_3 * pc_z[k] * nsf_59[k];

        t_90[k] = f_13 * msf_60[k]
                  + f_1 * nsd0_36[k]
                  - f_2 * nsd1_36[k]
                  + f_3 * pc_x[k] * nsf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, msf_30, msf_63, nsd0_39, \
                         nsd1_39, nsf_60, nsf_61, nsf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_14 * msf_30[k]
                  + f_3 * pc_y[k] * nsf_60[k];

        t_92[k] = f_3 * pc_z[k] * nsf_60[k];

        t_93[k] = f_13 * msf_63[k]
                  + f_4 * nsd0_39[k]
                  - f_5 * nsd1_39[k]
                  + f_3 * pc_x[k] * nsf_63[k];

        t_94[k] = f_3 * pc_z[k] * nsf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_z, msf_66, msf_68, nsd0_36, nsd1_36, \
                         nsf_62, nsf_63, nsf_66, nsf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * nsd0_36[k]
                  - f_5 * nsd1_36[k]
                  + f_3 * pc_z[k] * nsf_62[k];

        t_96[k] = f_13 * msf_66[k]
                  + f_3 * pc_x[k] * nsf_66[k];

        t_97[k] = f_3 * pc_z[k] * nsf_63[k];

        t_98[k] = f_13 * msf_68[k]
                  + f_3 * pc_x[k] * nsf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, msf_36, msf_39, \
                         msf_69, nsd0_39, nsd1_39, nsf_66, nsf_67, \
                         nsf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * msf_69[k]
                  + f_3 * pc_x[k] * nsf_69[k];

        t_100[k] = f_14 * msf_36[k]
                   + f_1 * nsd0_39[k]
                   - f_2 * nsd1_39[k]
                   + f_3 * pc_y[k] * nsf_66[k];

        t_101[k] = f_3 * pc_z[k] * nsf_66[k];

        t_102[k] = f_4 * nsd0_39[k]
                   - f_5 * nsd1_39[k]
                   + f_3 * pc_z[k] * nsf_67[k];

        t_103[k] = f_14 * msf_39[k]
                   + f_3 * pc_y[k] * nsf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_z, pc_y, pc_z, msg0_45, msf_30, \
                         msf_40, msg1_45, nsd0_41, nsd1_41, nsf_69, \
                         nsf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * nsd0_41[k]
                   - f_2 * nsd1_41[k]
                   + f_3 * pc_z[k] * nsf_69[k];

        t_105[k] = pa_z[k] * msg0_45[k]
                   - f_6 * pc_z[k] * msg1_45[k];

        t_106[k] = f_8 * msf_40[k]
                   + f_3 * pc_y[k] * nsf_70[k];

        t_107[k] = f_7 * msf_30[k]
                   + f_3 * pc_z[k] * nsf_70[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_z, pc_x, pc_y, pc_z, msg0_48, msf_42, msf_75, \
                         msg1_48, nsd0_47, nsd1_47, nsf_72, nsf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * msg0_48[k]
                   - f_6 * pc_z[k] * msg1_48[k];

        t_109[k] = f_8 * msf_42[k]
                   + f_3 * pc_y[k] * nsf_72[k];

        t_110[k] = f_13 * msf_75[k]
                   + f_4 * nsd0_47[k]
                   - f_5 * nsd1_47[k]
                   + f_3 * pc_x[k] * nsf_75[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, msf_76, msf_77, msf_78, msf_79, \
                         nsf_76, nsf_77, nsf_78, nsf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * msf_76[k]
                   + f_3 * pc_x[k] * nsf_76[k];

        t_112[k] = f_13 * msf_77[k]
                   + f_3 * pc_x[k] * nsf_77[k];

        t_113[k] = f_13 * msf_78[k]
                   + f_3 * pc_x[k] * nsf_78[k];

        t_114[k] = f_13 * msf_79[k]
                   + f_3 * pc_x[k] * nsf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_z, pc_y, pc_z, msg0_55, msf_36, msf_48, \
                         msg1_55, nsd0_47, nsd1_47, nsf_76, nsf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * msg0_55[k]
                   - f_6 * pc_z[k] * msg1_55[k];

        t_116[k] = f_7 * msf_36[k]
                   + f_3 * pc_z[k] * nsf_76[k];

        t_117[k] = f_8 * msf_48[k]
                   + f_4 * nsd0_47[k]
                   - f_5 * nsd1_47[k]
                   + f_3 * pc_y[k] * nsf_78[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pc_y, pc_z, msg0_75, msf_39, \
                         msf_49, msf_50, msg1_75, nsd0_47, nsd1_47, nsf_79, \
                         nsf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * msf_49[k]
                   + f_3 * pc_y[k] * nsf_79[k];

        t_119[k] = f_7 * msf_39[k]
                   + f_1 * nsd0_47[k]
                   - f_2 * nsd1_47[k]
                   + f_3 * pc_z[k] * nsf_79[k];

        t_120[k] = pa_y[k] * msg0_75[k]
                   - f_6 * pc_y[k] * msg1_75[k];

        t_121[k] = f_7 * msf_50[k]
                   + f_3 * pc_y[k] * nsf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pc_y, pc_z, msg0_78, msg0_80, \
                         msf_40, msf_51, msf_52, msg1_78, msg1_80, nsf_80, \
                         nsf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * msf_40[k]
                   + f_3 * pc_z[k] * nsf_80[k];

        t_123[k] = pa_y[k] * msg0_78[k]
                   + f_8 * msf_51[k]
                   - f_6 * pc_y[k] * msg1_78[k];

        t_124[k] = f_7 * msf_52[k]
                   + f_3 * pc_y[k] * nsf_82[k];

        t_125[k] = pa_y[k] * msg0_80[k]
                   - f_6 * pc_y[k] * msg1_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, msf_86, msf_87, msf_88, msf_89, \
                         nsf_86, nsf_87, nsf_88, nsf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_13 * msf_86[k]
                   + f_3 * pc_x[k] * nsf_86[k];

        t_127[k] = f_13 * msf_87[k]
                   + f_3 * pc_x[k] * nsf_87[k];

        t_128[k] = f_13 * msf_88[k]
                   + f_3 * pc_x[k] * nsf_88[k];

        t_129[k] = f_13 * msf_89[k]
                   + f_3 * pc_x[k] * nsf_89[k];
    }
}

static auto
compute_prim_nsg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msg0,
                                                          const size_t msf, const size_t msg1,
                                                          const size_t nsd0, const size_t nsd1,
                                                          const size_t nsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_13 = 3.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msg0_89 = buffer.data(msg0 + 89);
    const auto *msg0_90 = buffer.data(msg0 + 90);
    const auto *msg0_93 = buffer.data(msg0 + 93);
    const auto *msg0_100 = buffer.data(msg0 + 100);
    const auto *msg0_135 = buffer.data(msg0 + 135);
    const auto *msg0_138 = buffer.data(msg0 + 138);
    const auto *msg0_140 = buffer.data(msg0 + 140);
    const auto *msg0_149 = buffer.data(msg0 + 149);
    const auto *msg0_150 = buffer.data(msg0 + 150);
    const auto *msg0_153 = buffer.data(msg0 + 153);
    const auto *msg0_160 = buffer.data(msg0 + 160);

    const auto *msf_46 = buffer.data(msf + 46);
    const auto *msf_50 = buffer.data(msf + 50);
    const auto *msf_56 = buffer.data(msf + 56);
    const auto *msf_58 = buffer.data(msf + 58);
    const auto *msf_59 = buffer.data(msf + 59);
    const auto *msf_60 = buffer.data(msf + 60);
    const auto *msf_66 = buffer.data(msf + 66);
    const auto *msf_69 = buffer.data(msf + 69);
    const auto *msf_70 = buffer.data(msf + 70);
    const auto *msf_72 = buffer.data(msf + 72);
    const auto *msf_76 = buffer.data(msf + 76);
    const auto *msf_78 = buffer.data(msf + 78);
    const auto *msf_79 = buffer.data(msf + 79);
    const auto *msf_80 = buffer.data(msf + 80);
    const auto *msf_82 = buffer.data(msf + 82);
    const auto *msf_86 = buffer.data(msf + 86);
    const auto *msf_88 = buffer.data(msf + 88);
    const auto *msf_89 = buffer.data(msf + 89);
    const auto *msf_90 = buffer.data(msf + 90);
    const auto *msf_91 = buffer.data(msf + 91);
    const auto *msf_92 = buffer.data(msf + 92);
    const auto *msf_95 = buffer.data(msf + 95);
    const auto *msf_96 = buffer.data(msf + 96);
    const auto *msf_97 = buffer.data(msf + 97);
    const auto *msf_98 = buffer.data(msf + 98);
    const auto *msf_99 = buffer.data(msf + 99);
    const auto *msf_100 = buffer.data(msf + 100);
    const auto *msf_103 = buffer.data(msf + 103);
    const auto *msf_106 = buffer.data(msf + 106);
    const auto *msf_108 = buffer.data(msf + 108);
    const auto *msf_109 = buffer.data(msf + 109);
    const auto *msf_110 = buffer.data(msf + 110);
    const auto *msf_112 = buffer.data(msf + 112);
    const auto *msf_115 = buffer.data(msf + 115);
    const auto *msf_116 = buffer.data(msf + 116);
    const auto *msf_117 = buffer.data(msf + 117);
    const auto *msf_118 = buffer.data(msf + 118);
    const auto *msf_119 = buffer.data(msf + 119);
    const auto *msf_120 = buffer.data(msf + 120);
    const auto *msf_123 = buffer.data(msf + 123);
    const auto *msf_125 = buffer.data(msf + 125);
    const auto *msf_126 = buffer.data(msf + 126);
    const auto *msf_127 = buffer.data(msf + 127);
    const auto *msf_128 = buffer.data(msf + 128);
    const auto *msf_129 = buffer.data(msf + 129);
    const auto *msf_136 = buffer.data(msf + 136);
    const auto *msf_137 = buffer.data(msf + 137);
    const auto *msf_138 = buffer.data(msf + 138);
    const auto *msf_139 = buffer.data(msf + 139);
    const auto *msf_140 = buffer.data(msf + 140);
    const auto *msf_145 = buffer.data(msf + 145);
    const auto *msf_146 = buffer.data(msf + 146);
    const auto *msf_147 = buffer.data(msf + 147);
    const auto *msf_149 = buffer.data(msf + 149);
    const auto *msf_150 = buffer.data(msf + 150);
    const auto *msf_153 = buffer.data(msf + 153);
    const auto *msf_156 = buffer.data(msf + 156);
    const auto *msf_158 = buffer.data(msf + 158);
    const auto *msf_159 = buffer.data(msf + 159);
    const auto *msf_165 = buffer.data(msf + 165);
    const auto *msf_166 = buffer.data(msf + 166);
    const auto *msf_167 = buffer.data(msf + 167);
    const auto *msf_168 = buffer.data(msf + 168);
    const auto *msf_169 = buffer.data(msf + 169);

    const auto *msg1_89 = buffer.data(msg1 + 89);
    const auto *msg1_90 = buffer.data(msg1 + 90);
    const auto *msg1_93 = buffer.data(msg1 + 93);
    const auto *msg1_100 = buffer.data(msg1 + 100);
    const auto *msg1_135 = buffer.data(msg1 + 135);
    const auto *msg1_138 = buffer.data(msg1 + 138);
    const auto *msg1_140 = buffer.data(msg1 + 140);
    const auto *msg1_149 = buffer.data(msg1 + 149);
    const auto *msg1_150 = buffer.data(msg1 + 150);
    const auto *msg1_153 = buffer.data(msg1 + 153);
    const auto *msg1_160 = buffer.data(msg1 + 160);

    const auto *nsd0_51 = buffer.data(nsd0 + 51);
    const auto *nsd0_53 = buffer.data(nsd0 + 53);
    const auto *nsd0_54 = buffer.data(nsd0 + 54);
    const auto *nsd0_57 = buffer.data(nsd0 + 57);
    const auto *nsd0_58 = buffer.data(nsd0 + 58);
    const auto *nsd0_59 = buffer.data(nsd0 + 59);
    const auto *nsd0_60 = buffer.data(nsd0 + 60);
    const auto *nsd0_63 = buffer.data(nsd0 + 63);
    const auto *nsd0_65 = buffer.data(nsd0 + 65);
    const auto *nsd0_71 = buffer.data(nsd0 + 71);
    const auto *nsd0_72 = buffer.data(nsd0 + 72);
    const auto *nsd0_75 = buffer.data(nsd0 + 75);
    const auto *nsd0_77 = buffer.data(nsd0 + 77);
    const auto *nsd0_81 = buffer.data(nsd0 + 81);
    const auto *nsd0_83 = buffer.data(nsd0 + 83);
    const auto *nsd0_84 = buffer.data(nsd0 + 84);
    const auto *nsd0_87 = buffer.data(nsd0 + 87);
    const auto *nsd0_88 = buffer.data(nsd0 + 88);
    const auto *nsd0_89 = buffer.data(nsd0 + 89);
    const auto *nsd0_90 = buffer.data(nsd0 + 90);
    const auto *nsd0_93 = buffer.data(nsd0 + 93);
    const auto *nsd0_95 = buffer.data(nsd0 + 95);
    const auto *nsd0_101 = buffer.data(nsd0 + 101);

    const auto *nsd1_51 = buffer.data(nsd1 + 51);
    const auto *nsd1_53 = buffer.data(nsd1 + 53);
    const auto *nsd1_54 = buffer.data(nsd1 + 54);
    const auto *nsd1_57 = buffer.data(nsd1 + 57);
    const auto *nsd1_58 = buffer.data(nsd1 + 58);
    const auto *nsd1_59 = buffer.data(nsd1 + 59);
    const auto *nsd1_60 = buffer.data(nsd1 + 60);
    const auto *nsd1_63 = buffer.data(nsd1 + 63);
    const auto *nsd1_65 = buffer.data(nsd1 + 65);
    const auto *nsd1_71 = buffer.data(nsd1 + 71);
    const auto *nsd1_72 = buffer.data(nsd1 + 72);
    const auto *nsd1_75 = buffer.data(nsd1 + 75);
    const auto *nsd1_77 = buffer.data(nsd1 + 77);
    const auto *nsd1_81 = buffer.data(nsd1 + 81);
    const auto *nsd1_83 = buffer.data(nsd1 + 83);
    const auto *nsd1_84 = buffer.data(nsd1 + 84);
    const auto *nsd1_87 = buffer.data(nsd1 + 87);
    const auto *nsd1_88 = buffer.data(nsd1 + 88);
    const auto *nsd1_89 = buffer.data(nsd1 + 89);
    const auto *nsd1_90 = buffer.data(nsd1 + 90);
    const auto *nsd1_93 = buffer.data(nsd1 + 93);
    const auto *nsd1_95 = buffer.data(nsd1 + 95);
    const auto *nsd1_101 = buffer.data(nsd1 + 101);

    const auto *nsf_86 = buffer.data(nsf + 86);
    const auto *nsf_88 = buffer.data(nsf + 88);
    const auto *nsf_89 = buffer.data(nsf + 89);
    const auto *nsf_90 = buffer.data(nsf + 90);
    const auto *nsf_91 = buffer.data(nsf + 91);
    const auto *nsf_92 = buffer.data(nsf + 92);
    const auto *nsf_95 = buffer.data(nsf + 95);
    const auto *nsf_96 = buffer.data(nsf + 96);
    const auto *nsf_97 = buffer.data(nsf + 97);
    const auto *nsf_98 = buffer.data(nsf + 98);
    const auto *nsf_99 = buffer.data(nsf + 99);
    const auto *nsf_100 = buffer.data(nsf + 100);
    const auto *nsf_101 = buffer.data(nsf + 101);
    const auto *nsf_102 = buffer.data(nsf + 102);
    const auto *nsf_103 = buffer.data(nsf + 103);
    const auto *nsf_106 = buffer.data(nsf + 106);
    const auto *nsf_107 = buffer.data(nsf + 107);
    const auto *nsf_108 = buffer.data(nsf + 108);
    const auto *nsf_109 = buffer.data(nsf + 109);
    const auto *nsf_110 = buffer.data(nsf + 110);
    const auto *nsf_112 = buffer.data(nsf + 112);
    const auto *nsf_115 = buffer.data(nsf + 115);
    const auto *nsf_116 = buffer.data(nsf + 116);
    const auto *nsf_117 = buffer.data(nsf + 117);
    const auto *nsf_118 = buffer.data(nsf + 118);
    const auto *nsf_119 = buffer.data(nsf + 119);
    const auto *nsf_120 = buffer.data(nsf + 120);
    const auto *nsf_122 = buffer.data(nsf + 122);
    const auto *nsf_123 = buffer.data(nsf + 123);
    const auto *nsf_125 = buffer.data(nsf + 125);
    const auto *nsf_126 = buffer.data(nsf + 126);
    const auto *nsf_127 = buffer.data(nsf + 127);
    const auto *nsf_128 = buffer.data(nsf + 128);
    const auto *nsf_129 = buffer.data(nsf + 129);
    const auto *nsf_130 = buffer.data(nsf + 130);
    const auto *nsf_132 = buffer.data(nsf + 132);
    const auto *nsf_136 = buffer.data(nsf + 136);
    const auto *nsf_137 = buffer.data(nsf + 137);
    const auto *nsf_138 = buffer.data(nsf + 138);
    const auto *nsf_139 = buffer.data(nsf + 139);
    const auto *nsf_140 = buffer.data(nsf + 140);
    const auto *nsf_141 = buffer.data(nsf + 141);
    const auto *nsf_142 = buffer.data(nsf + 142);
    const auto *nsf_145 = buffer.data(nsf + 145);
    const auto *nsf_146 = buffer.data(nsf + 146);
    const auto *nsf_147 = buffer.data(nsf + 147);
    const auto *nsf_148 = buffer.data(nsf + 148);
    const auto *nsf_149 = buffer.data(nsf + 149);
    const auto *nsf_150 = buffer.data(nsf + 150);
    const auto *nsf_151 = buffer.data(nsf + 151);
    const auto *nsf_152 = buffer.data(nsf + 152);
    const auto *nsf_153 = buffer.data(nsf + 153);
    const auto *nsf_156 = buffer.data(nsf + 156);
    const auto *nsf_157 = buffer.data(nsf + 157);
    const auto *nsf_158 = buffer.data(nsf + 158);
    const auto *nsf_159 = buffer.data(nsf + 159);
    const auto *nsf_160 = buffer.data(nsf + 160);
    const auto *nsf_162 = buffer.data(nsf + 162);
    const auto *nsf_165 = buffer.data(nsf + 165);
    const auto *nsf_166 = buffer.data(nsf + 166);
    const auto *nsf_167 = buffer.data(nsf + 167);
    const auto *nsf_168 = buffer.data(nsf + 168);
    const auto *nsf_169 = buffer.data(nsf + 169);

#pragma omp simd aligned(t_130, t_131, t_132, pc_y, pc_z, msf_46, msf_56, msf_58, nsd0_51, \
                         nsd0_53, nsd1_51, nsd1_53, nsf_86, nsf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * msf_56[k]
                   + f_1 * nsd0_51[k]
                   - f_2 * nsd1_51[k]
                   + f_3 * pc_y[k] * nsf_86[k];

        t_131[k] = f_8 * msf_46[k]
                   + f_3 * pc_z[k] * nsf_86[k];

        t_132[k] = f_7 * msf_58[k]
                   + f_4 * nsd0_53[k]
                   - f_5 * nsd1_53[k]
                   + f_3 * pc_y[k] * nsf_88[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pc_x, pc_y, msg0_89, msf_59, \
                         msf_90, msg1_89, nsd0_54, nsd1_54, nsf_89, \
                         nsf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * msf_59[k]
                   + f_3 * pc_y[k] * nsf_89[k];

        t_134[k] = pa_y[k] * msg0_89[k]
                   - f_6 * pc_y[k] * msg1_89[k];

        t_135[k] = f_13 * msf_90[k]
                   + f_1 * nsd0_54[k]
                   - f_2 * nsd1_54[k]
                   + f_3 * pc_x[k] * nsf_90[k];

        t_136[k] = f_3 * pc_y[k] * nsf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, msf_50, nsd0_54, nsd1_54, nsf_90, \
                         nsf_91, nsf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * msf_50[k]
                   + f_3 * pc_z[k] * nsf_90[k];

        t_138[k] = f_4 * nsd0_54[k]
                   - f_5 * nsd1_54[k]
                   + f_3 * pc_y[k] * nsf_91[k];

        t_139[k] = f_3 * pc_y[k] * nsf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, msf_95, msf_96, msf_97, \
                         nsd0_59, nsd1_59, nsf_95, nsf_96, nsf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * msf_95[k]
                   + f_4 * nsd0_59[k]
                   - f_5 * nsd1_59[k]
                   + f_3 * pc_x[k] * nsf_95[k];

        t_141[k] = f_13 * msf_96[k]
                   + f_3 * pc_x[k] * nsf_96[k];

        t_142[k] = f_13 * msf_97[k]
                   + f_3 * pc_x[k] * nsf_97[k];

        t_143[k] = f_3 * pc_y[k] * nsf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_y, msf_99, nsd0_57, nsd0_58, nsd1_57, \
                         nsd1_58, nsf_96, nsf_97, nsf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_13 * msf_99[k]
                   + f_3 * pc_x[k] * nsf_99[k];

        t_145[k] = f_1 * nsd0_57[k]
                   - f_2 * nsd1_57[k]
                   + f_3 * pc_y[k] * nsf_96[k];

        t_146[k] = f_10 * nsd0_58[k]
                   - f_11 * nsd1_58[k]
                   + f_3 * pc_y[k] * nsf_97[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, msf_59, msf_100, \
                         nsd0_59, nsd0_60, nsd1_59, nsd1_60, nsf_98, nsf_99, \
                         nsf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * nsd0_59[k]
                   - f_5 * nsd1_59[k]
                   + f_3 * pc_y[k] * nsf_98[k];

        t_148[k] = f_3 * pc_y[k] * nsf_99[k];

        t_149[k] = f_14 * msf_59[k]
                   + f_1 * nsd0_59[k]
                   - f_2 * nsd1_59[k]
                   + f_3 * pc_z[k] * nsf_99[k];

        t_150[k] = f_15 * msf_100[k]
                   + f_1 * nsd0_60[k]
                   - f_2 * nsd1_60[k]
                   + f_3 * pc_x[k] * nsf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, msf_60, msf_103, \
                         nsd0_63, nsd1_63, nsf_100, nsf_101, nsf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_16 * msf_60[k]
                   + f_3 * pc_y[k] * nsf_100[k];

        t_152[k] = f_3 * pc_z[k] * nsf_100[k];

        t_153[k] = f_15 * msf_103[k]
                   + f_4 * nsd0_63[k]
                   - f_5 * nsd1_63[k]
                   + f_3 * pc_x[k] * nsf_103[k];

        t_154[k] = f_3 * pc_z[k] * nsf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_z, msf_106, msf_108, nsd0_60, \
                         nsd1_60, nsf_102, nsf_103, nsf_106, nsf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_4 * nsd0_60[k]
                   - f_5 * nsd1_60[k]
                   + f_3 * pc_z[k] * nsf_102[k];

        t_156[k] = f_15 * msf_106[k]
                   + f_3 * pc_x[k] * nsf_106[k];

        t_157[k] = f_3 * pc_z[k] * nsf_103[k];

        t_158[k] = f_15 * msf_108[k]
                   + f_3 * pc_x[k] * nsf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, msf_66, msf_69, \
                         msf_109, nsd0_63, nsd1_63, nsf_106, nsf_107, \
                         nsf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_15 * msf_109[k]
                   + f_3 * pc_x[k] * nsf_109[k];

        t_160[k] = f_16 * msf_66[k]
                   + f_1 * nsd0_63[k]
                   - f_2 * nsd1_63[k]
                   + f_3 * pc_y[k] * nsf_106[k];

        t_161[k] = f_3 * pc_z[k] * nsf_106[k];

        t_162[k] = f_4 * nsd0_63[k]
                   - f_5 * nsd1_63[k]
                   + f_3 * pc_z[k] * nsf_107[k];

        t_163[k] = f_16 * msf_69[k]
                   + f_3 * pc_y[k] * nsf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_z, pc_y, pc_z, msg0_90, msf_60, \
                         msf_70, msg1_90, nsd0_65, nsd1_65, nsf_109, \
                         nsf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * nsd0_65[k]
                   - f_2 * nsd1_65[k]
                   + f_3 * pc_z[k] * nsf_109[k];

        t_165[k] = pa_z[k] * msg0_90[k]
                   - f_6 * pc_z[k] * msg1_90[k];

        t_166[k] = f_14 * msf_70[k]
                   + f_3 * pc_y[k] * nsf_110[k];

        t_167[k] = f_7 * msf_60[k]
                   + f_3 * pc_z[k] * nsf_110[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_z, pc_x, pc_y, pc_z, msg0_93, msf_72, \
                         msf_115, msg1_93, nsd0_71, nsd1_71, nsf_112, \
                         nsf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * msg0_93[k]
                   - f_6 * pc_z[k] * msg1_93[k];

        t_169[k] = f_14 * msf_72[k]
                   + f_3 * pc_y[k] * nsf_112[k];

        t_170[k] = f_15 * msf_115[k]
                   + f_4 * nsd0_71[k]
                   - f_5 * nsd1_71[k]
                   + f_3 * pc_x[k] * nsf_115[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pc_x, msf_116, msf_117, msf_118, msf_119, \
                         nsf_116, nsf_117, nsf_118, nsf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_15 * msf_116[k]
                   + f_3 * pc_x[k] * nsf_116[k];

        t_172[k] = f_15 * msf_117[k]
                   + f_3 * pc_x[k] * nsf_117[k];

        t_173[k] = f_15 * msf_118[k]
                   + f_3 * pc_x[k] * nsf_118[k];

        t_174[k] = f_15 * msf_119[k]
                   + f_3 * pc_x[k] * nsf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_z, pc_y, pc_z, msg0_100, msf_66, msf_78, \
                         msg1_100, nsd0_71, nsd1_71, nsf_116, nsf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * msg0_100[k]
                   - f_6 * pc_z[k] * msg1_100[k];

        t_176[k] = f_7 * msf_66[k]
                   + f_3 * pc_z[k] * nsf_116[k];

        t_177[k] = f_14 * msf_78[k]
                   + f_4 * nsd0_71[k]
                   - f_5 * nsd1_71[k]
                   + f_3 * pc_y[k] * nsf_118[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, msf_69, msf_79, msf_120, \
                         nsd0_71, nsd0_72, nsd1_71, nsd1_72, nsf_119, \
                         nsf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * msf_79[k]
                   + f_3 * pc_y[k] * nsf_119[k];

        t_179[k] = f_7 * msf_69[k]
                   + f_1 * nsd0_71[k]
                   - f_2 * nsd1_71[k]
                   + f_3 * pc_z[k] * nsf_119[k];

        t_180[k] = f_15 * msf_120[k]
                   + f_1 * nsd0_72[k]
                   - f_2 * nsd1_72[k]
                   + f_3 * pc_x[k] * nsf_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, msf_70, msf_80, msf_82, \
                         msf_123, nsd0_75, nsd1_75, nsf_120, nsf_122, \
                         nsf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * msf_80[k]
                   + f_3 * pc_y[k] * nsf_120[k];

        t_182[k] = f_8 * msf_70[k]
                   + f_3 * pc_z[k] * nsf_120[k];

        t_183[k] = f_15 * msf_123[k]
                   + f_4 * nsd0_75[k]
                   - f_5 * nsd1_75[k]
                   + f_3 * pc_x[k] * nsf_123[k];

        t_184[k] = f_8 * msf_82[k]
                   + f_3 * pc_y[k] * nsf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, msf_125, msf_126, msf_127, msf_128, \
                         nsd0_77, nsd1_77, nsf_125, nsf_126, nsf_127, \
                         nsf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_15 * msf_125[k]
                   + f_4 * nsd0_77[k]
                   - f_5 * nsd1_77[k]
                   + f_3 * pc_x[k] * nsf_125[k];

        t_186[k] = f_15 * msf_126[k]
                   + f_3 * pc_x[k] * nsf_126[k];

        t_187[k] = f_15 * msf_127[k]
                   + f_3 * pc_x[k] * nsf_127[k];

        t_188[k] = f_15 * msf_128[k]
                   + f_3 * pc_x[k] * nsf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, pc_z, msf_76, msf_86, msf_129, \
                         nsd0_75, nsd1_75, nsf_126, nsf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_15 * msf_129[k]
                   + f_3 * pc_x[k] * nsf_129[k];

        t_190[k] = f_8 * msf_86[k]
                   + f_1 * nsd0_75[k]
                   - f_2 * nsd1_75[k]
                   + f_3 * pc_y[k] * nsf_126[k];

        t_191[k] = f_8 * msf_76[k]
                   + f_3 * pc_z[k] * nsf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_y, pc_y, pc_z, msg0_135, msf_79, \
                         msf_88, msf_89, msg1_135, nsd0_77, nsd1_77, nsf_128, \
                         nsf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_8 * msf_88[k]
                   + f_4 * nsd0_77[k]
                   - f_5 * nsd1_77[k]
                   + f_3 * pc_y[k] * nsf_128[k];

        t_193[k] = f_8 * msf_89[k]
                   + f_3 * pc_y[k] * nsf_129[k];

        t_194[k] = f_8 * msf_79[k]
                   + f_1 * nsd0_77[k]
                   - f_2 * nsd1_77[k]
                   + f_3 * pc_z[k] * nsf_129[k];

        t_195[k] = pa_y[k] * msg0_135[k]
                   - f_6 * pc_y[k] * msg1_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pc_y, pc_z, msg0_138, msf_80, \
                         msf_90, msf_91, msf_92, msg1_138, nsf_130, \
                         nsf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_7 * msf_90[k]
                   + f_3 * pc_y[k] * nsf_130[k];

        t_197[k] = f_14 * msf_80[k]
                   + f_3 * pc_z[k] * nsf_130[k];

        t_198[k] = pa_y[k] * msg0_138[k]
                   + f_8 * msf_91[k]
                   - f_6 * pc_y[k] * msg1_138[k];

        t_199[k] = f_7 * msf_92[k]
                   + f_3 * pc_y[k] * nsf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pc_x, pc_y, msg0_140, msf_136, \
                         msf_137, msf_138, msg1_140, nsf_136, nsf_137, \
                         nsf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_y[k] * msg0_140[k]
                   - f_6 * pc_y[k] * msg1_140[k];

        t_201[k] = f_15 * msf_136[k]
                   + f_3 * pc_x[k] * nsf_136[k];

        t_202[k] = f_15 * msf_137[k]
                   + f_3 * pc_x[k] * nsf_137[k];

        t_203[k] = f_15 * msf_138[k]
                   + f_3 * pc_x[k] * nsf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pc_x, pc_y, pc_z, msf_86, msf_96, msf_139, \
                         nsd0_81, nsd1_81, nsf_136, nsf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_15 * msf_139[k]
                   + f_3 * pc_x[k] * nsf_139[k];

        t_205[k] = f_7 * msf_96[k]
                   + f_1 * nsd0_81[k]
                   - f_2 * nsd1_81[k]
                   + f_3 * pc_y[k] * nsf_136[k];

        t_206[k] = f_14 * msf_86[k]
                   + f_3 * pc_z[k] * nsf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_y, pc_y, msg0_149, msf_98, msf_99, msg1_149, \
                         nsd0_83, nsd1_83, nsf_138, nsf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_7 * msf_98[k]
                   + f_4 * nsd0_83[k]
                   - f_5 * nsd1_83[k]
                   + f_3 * pc_y[k] * nsf_138[k];

        t_208[k] = f_7 * msf_99[k]
                   + f_3 * pc_y[k] * nsf_139[k];

        t_209[k] = pa_y[k] * msg0_149[k]
                   - f_6 * pc_y[k] * msg1_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, pc_y, pc_z, msf_90, msf_140, \
                         nsd0_84, nsd1_84, nsf_140, nsf_141, nsf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_15 * msf_140[k]
                   + f_1 * nsd0_84[k]
                   - f_2 * nsd1_84[k]
                   + f_3 * pc_x[k] * nsf_140[k];

        t_211[k] = f_3 * pc_y[k] * nsf_140[k];

        t_212[k] = f_16 * msf_90[k]
                   + f_3 * pc_z[k] * nsf_140[k];

        t_213[k] = f_4 * nsd0_84[k]
                   - f_5 * nsd1_84[k]
                   + f_3 * pc_y[k] * nsf_141[k];

        t_214[k] = f_3 * pc_y[k] * nsf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, pc_y, msf_145, msf_146, msf_147, \
                         nsd0_89, nsd1_89, nsf_145, nsf_146, nsf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * msf_145[k]
                   + f_4 * nsd0_89[k]
                   - f_5 * nsd1_89[k]
                   + f_3 * pc_x[k] * nsf_145[k];

        t_216[k] = f_15 * msf_146[k]
                   + f_3 * pc_x[k] * nsf_146[k];

        t_217[k] = f_15 * msf_147[k]
                   + f_3 * pc_x[k] * nsf_147[k];

        t_218[k] = f_3 * pc_y[k] * nsf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pc_x, pc_y, msf_149, nsd0_87, nsd0_88, nsd1_87, \
                         nsd1_88, nsf_146, nsf_147, nsf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_15 * msf_149[k]
                   + f_3 * pc_x[k] * nsf_149[k];

        t_220[k] = f_1 * nsd0_87[k]
                   - f_2 * nsd1_87[k]
                   + f_3 * pc_y[k] * nsf_146[k];

        t_221[k] = f_10 * nsd0_88[k]
                   - f_11 * nsd1_88[k]
                   + f_3 * pc_y[k] * nsf_147[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, pc_z, msf_99, msf_150, \
                         nsd0_89, nsd0_90, nsd1_89, nsd1_90, nsf_148, nsf_149, \
                         nsf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_4 * nsd0_89[k]
                   - f_5 * nsd1_89[k]
                   + f_3 * pc_y[k] * nsf_148[k];

        t_223[k] = f_3 * pc_y[k] * nsf_149[k];

        t_224[k] = f_16 * msf_99[k]
                   + f_1 * nsd0_89[k]
                   - f_2 * nsd1_89[k]
                   + f_3 * pc_z[k] * nsf_149[k];

        t_225[k] = f_17 * msf_150[k]
                   + f_1 * nsd0_90[k]
                   - f_2 * nsd1_90[k]
                   + f_3 * pc_x[k] * nsf_150[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_x, pc_y, pc_z, msf_100, msf_153, \
                         nsd0_93, nsd1_93, nsf_150, nsf_151, nsf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_17 * msf_100[k]
                   + f_3 * pc_y[k] * nsf_150[k];

        t_227[k] = f_3 * pc_z[k] * nsf_150[k];

        t_228[k] = f_17 * msf_153[k]
                   + f_4 * nsd0_93[k]
                   - f_5 * nsd1_93[k]
                   + f_3 * pc_x[k] * nsf_153[k];

        t_229[k] = f_3 * pc_z[k] * nsf_151[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_z, msf_156, msf_158, nsd0_90, \
                         nsd1_90, nsf_152, nsf_153, nsf_156, nsf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_4 * nsd0_90[k]
                   - f_5 * nsd1_90[k]
                   + f_3 * pc_z[k] * nsf_152[k];

        t_231[k] = f_17 * msf_156[k]
                   + f_3 * pc_x[k] * nsf_156[k];

        t_232[k] = f_3 * pc_z[k] * nsf_153[k];

        t_233[k] = f_17 * msf_158[k]
                   + f_3 * pc_x[k] * nsf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pc_x, pc_y, pc_z, msf_106, \
                         msf_109, msf_159, nsd0_93, nsd1_93, nsf_156, nsf_157, \
                         nsf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_17 * msf_159[k]
                   + f_3 * pc_x[k] * nsf_159[k];

        t_235[k] = f_17 * msf_106[k]
                   + f_1 * nsd0_93[k]
                   - f_2 * nsd1_93[k]
                   + f_3 * pc_y[k] * nsf_156[k];

        t_236[k] = f_3 * pc_z[k] * nsf_156[k];

        t_237[k] = f_4 * nsd0_93[k]
                   - f_5 * nsd1_93[k]
                   + f_3 * pc_z[k] * nsf_157[k];

        t_238[k] = f_17 * msf_109[k]
                   + f_3 * pc_y[k] * nsf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_z, pc_y, pc_z, msg0_150, msf_100, \
                         msf_110, msg1_150, nsd0_95, nsd1_95, nsf_159, \
                         nsf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * nsd0_95[k]
                   - f_2 * nsd1_95[k]
                   + f_3 * pc_z[k] * nsf_159[k];

        t_240[k] = pa_z[k] * msg0_150[k]
                   - f_6 * pc_z[k] * msg1_150[k];

        t_241[k] = f_16 * msf_110[k]
                   + f_3 * pc_y[k] * nsf_160[k];

        t_242[k] = f_7 * msf_100[k]
                   + f_3 * pc_z[k] * nsf_160[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_z, pc_x, pc_y, pc_z, msg0_153, msf_112, \
                         msf_165, msg1_153, nsd0_101, nsd1_101, nsf_162, \
                         nsf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pa_z[k] * msg0_153[k]
                   - f_6 * pc_z[k] * msg1_153[k];

        t_244[k] = f_16 * msf_112[k]
                   + f_3 * pc_y[k] * nsf_162[k];

        t_245[k] = f_17 * msf_165[k]
                   + f_4 * nsd0_101[k]
                   - f_5 * nsd1_101[k]
                   + f_3 * pc_x[k] * nsf_165[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, msf_166, msf_167, msf_168, msf_169, \
                         nsf_166, nsf_167, nsf_168, nsf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_17 * msf_166[k]
                   + f_3 * pc_x[k] * nsf_166[k];

        t_247[k] = f_17 * msf_167[k]
                   + f_3 * pc_x[k] * nsf_167[k];

        t_248[k] = f_17 * msf_168[k]
                   + f_3 * pc_x[k] * nsf_168[k];

        t_249[k] = f_17 * msf_169[k]
                   + f_3 * pc_x[k] * nsf_169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pa_z, pc_y, pc_z, msg0_160, msf_106, msf_118, \
                         msg1_160, nsd0_101, nsd1_101, nsf_166, \
                         nsf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = pa_z[k] * msg0_160[k]
                   - f_6 * pc_z[k] * msg1_160[k];

        t_251[k] = f_7 * msf_106[k]
                   + f_3 * pc_z[k] * nsf_166[k];

        t_252[k] = f_16 * msf_118[k]
                   + f_4 * nsd0_101[k]
                   - f_5 * nsd1_101[k]
                   + f_3 * pc_y[k] * nsf_168[k];
    }
}

static auto
compute_prim_nsg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msg0,
                                                          const size_t msf, const size_t msg1,
                                                          const size_t nsd0, const size_t nsd1,
                                                          const size_t nsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msg0_210 = buffer.data(msg0 + 210);
    const auto *msg0_213 = buffer.data(msg0 + 213);
    const auto *msg0_215 = buffer.data(msg0 + 215);
    const auto *msg0_224 = buffer.data(msg0 + 224);
    const auto *msg0_225 = buffer.data(msg0 + 225);
    const auto *msg0_228 = buffer.data(msg0 + 228);
    const auto *msg0_235 = buffer.data(msg0 + 235);

    const auto *msf_109 = buffer.data(msf + 109);
    const auto *msf_110 = buffer.data(msf + 110);
    const auto *msf_116 = buffer.data(msf + 116);
    const auto *msf_119 = buffer.data(msf + 119);
    const auto *msf_120 = buffer.data(msf + 120);
    const auto *msf_122 = buffer.data(msf + 122);
    const auto *msf_126 = buffer.data(msf + 126);
    const auto *msf_128 = buffer.data(msf + 128);
    const auto *msf_129 = buffer.data(msf + 129);
    const auto *msf_130 = buffer.data(msf + 130);
    const auto *msf_132 = buffer.data(msf + 132);
    const auto *msf_136 = buffer.data(msf + 136);
    const auto *msf_138 = buffer.data(msf + 138);
    const auto *msf_139 = buffer.data(msf + 139);
    const auto *msf_140 = buffer.data(msf + 140);
    const auto *msf_141 = buffer.data(msf + 141);
    const auto *msf_142 = buffer.data(msf + 142);
    const auto *msf_146 = buffer.data(msf + 146);
    const auto *msf_148 = buffer.data(msf + 148);
    const auto *msf_149 = buffer.data(msf + 149);
    const auto *msf_150 = buffer.data(msf + 150);
    const auto *msf_156 = buffer.data(msf + 156);
    const auto *msf_159 = buffer.data(msf + 159);
    const auto *msf_160 = buffer.data(msf + 160);
    const auto *msf_162 = buffer.data(msf + 162);
    const auto *msf_166 = buffer.data(msf + 166);
    const auto *msf_168 = buffer.data(msf + 168);
    const auto *msf_169 = buffer.data(msf + 169);
    const auto *msf_170 = buffer.data(msf + 170);
    const auto *msf_172 = buffer.data(msf + 172);
    const auto *msf_173 = buffer.data(msf + 173);
    const auto *msf_175 = buffer.data(msf + 175);
    const auto *msf_176 = buffer.data(msf + 176);
    const auto *msf_177 = buffer.data(msf + 177);
    const auto *msf_178 = buffer.data(msf + 178);
    const auto *msf_179 = buffer.data(msf + 179);
    const auto *msf_180 = buffer.data(msf + 180);
    const auto *msf_182 = buffer.data(msf + 182);
    const auto *msf_183 = buffer.data(msf + 183);
    const auto *msf_185 = buffer.data(msf + 185);
    const auto *msf_186 = buffer.data(msf + 186);
    const auto *msf_187 = buffer.data(msf + 187);
    const auto *msf_188 = buffer.data(msf + 188);
    const auto *msf_189 = buffer.data(msf + 189);
    const auto *msf_196 = buffer.data(msf + 196);
    const auto *msf_197 = buffer.data(msf + 197);
    const auto *msf_198 = buffer.data(msf + 198);
    const auto *msf_199 = buffer.data(msf + 199);
    const auto *msf_200 = buffer.data(msf + 200);
    const auto *msf_205 = buffer.data(msf + 205);
    const auto *msf_206 = buffer.data(msf + 206);
    const auto *msf_207 = buffer.data(msf + 207);
    const auto *msf_209 = buffer.data(msf + 209);
    const auto *msf_210 = buffer.data(msf + 210);
    const auto *msf_213 = buffer.data(msf + 213);
    const auto *msf_216 = buffer.data(msf + 216);
    const auto *msf_218 = buffer.data(msf + 218);
    const auto *msf_219 = buffer.data(msf + 219);
    const auto *msf_225 = buffer.data(msf + 225);
    const auto *msf_226 = buffer.data(msf + 226);
    const auto *msf_227 = buffer.data(msf + 227);
    const auto *msf_228 = buffer.data(msf + 228);
    const auto *msf_229 = buffer.data(msf + 229);
    const auto *msf_230 = buffer.data(msf + 230);
    const auto *msf_233 = buffer.data(msf + 233);
    const auto *msf_235 = buffer.data(msf + 235);
    const auto *msf_236 = buffer.data(msf + 236);
    const auto *msf_237 = buffer.data(msf + 237);
    const auto *msf_238 = buffer.data(msf + 238);
    const auto *msf_239 = buffer.data(msf + 239);
    const auto *msf_240 = buffer.data(msf + 240);
    const auto *msf_243 = buffer.data(msf + 243);
    const auto *msf_245 = buffer.data(msf + 245);
    const auto *msf_246 = buffer.data(msf + 246);
    const auto *msf_247 = buffer.data(msf + 247);
    const auto *msf_248 = buffer.data(msf + 248);
    const auto *msf_249 = buffer.data(msf + 249);

    const auto *msg1_210 = buffer.data(msg1 + 210);
    const auto *msg1_213 = buffer.data(msg1 + 213);
    const auto *msg1_215 = buffer.data(msg1 + 215);
    const auto *msg1_224 = buffer.data(msg1 + 224);
    const auto *msg1_225 = buffer.data(msg1 + 225);
    const auto *msg1_228 = buffer.data(msg1 + 228);
    const auto *msg1_235 = buffer.data(msg1 + 235);

    const auto *nsd0_101 = buffer.data(nsd0 + 101);
    const auto *nsd0_102 = buffer.data(nsd0 + 102);
    const auto *nsd0_105 = buffer.data(nsd0 + 105);
    const auto *nsd0_107 = buffer.data(nsd0 + 107);
    const auto *nsd0_108 = buffer.data(nsd0 + 108);
    const auto *nsd0_111 = buffer.data(nsd0 + 111);
    const auto *nsd0_113 = buffer.data(nsd0 + 113);
    const auto *nsd0_117 = buffer.data(nsd0 + 117);
    const auto *nsd0_119 = buffer.data(nsd0 + 119);
    const auto *nsd0_120 = buffer.data(nsd0 + 120);
    const auto *nsd0_123 = buffer.data(nsd0 + 123);
    const auto *nsd0_124 = buffer.data(nsd0 + 124);
    const auto *nsd0_125 = buffer.data(nsd0 + 125);
    const auto *nsd0_126 = buffer.data(nsd0 + 126);
    const auto *nsd0_129 = buffer.data(nsd0 + 129);
    const auto *nsd0_131 = buffer.data(nsd0 + 131);
    const auto *nsd0_137 = buffer.data(nsd0 + 137);
    const auto *nsd0_138 = buffer.data(nsd0 + 138);
    const auto *nsd0_141 = buffer.data(nsd0 + 141);
    const auto *nsd0_143 = buffer.data(nsd0 + 143);
    const auto *nsd0_144 = buffer.data(nsd0 + 144);
    const auto *nsd0_147 = buffer.data(nsd0 + 147);
    const auto *nsd0_149 = buffer.data(nsd0 + 149);

    const auto *nsd1_101 = buffer.data(nsd1 + 101);
    const auto *nsd1_102 = buffer.data(nsd1 + 102);
    const auto *nsd1_105 = buffer.data(nsd1 + 105);
    const auto *nsd1_107 = buffer.data(nsd1 + 107);
    const auto *nsd1_108 = buffer.data(nsd1 + 108);
    const auto *nsd1_111 = buffer.data(nsd1 + 111);
    const auto *nsd1_113 = buffer.data(nsd1 + 113);
    const auto *nsd1_117 = buffer.data(nsd1 + 117);
    const auto *nsd1_119 = buffer.data(nsd1 + 119);
    const auto *nsd1_120 = buffer.data(nsd1 + 120);
    const auto *nsd1_123 = buffer.data(nsd1 + 123);
    const auto *nsd1_124 = buffer.data(nsd1 + 124);
    const auto *nsd1_125 = buffer.data(nsd1 + 125);
    const auto *nsd1_126 = buffer.data(nsd1 + 126);
    const auto *nsd1_129 = buffer.data(nsd1 + 129);
    const auto *nsd1_131 = buffer.data(nsd1 + 131);
    const auto *nsd1_137 = buffer.data(nsd1 + 137);
    const auto *nsd1_138 = buffer.data(nsd1 + 138);
    const auto *nsd1_141 = buffer.data(nsd1 + 141);
    const auto *nsd1_143 = buffer.data(nsd1 + 143);
    const auto *nsd1_144 = buffer.data(nsd1 + 144);
    const auto *nsd1_147 = buffer.data(nsd1 + 147);
    const auto *nsd1_149 = buffer.data(nsd1 + 149);

    const auto *nsf_169 = buffer.data(nsf + 169);
    const auto *nsf_170 = buffer.data(nsf + 170);
    const auto *nsf_172 = buffer.data(nsf + 172);
    const auto *nsf_173 = buffer.data(nsf + 173);
    const auto *nsf_175 = buffer.data(nsf + 175);
    const auto *nsf_176 = buffer.data(nsf + 176);
    const auto *nsf_177 = buffer.data(nsf + 177);
    const auto *nsf_178 = buffer.data(nsf + 178);
    const auto *nsf_179 = buffer.data(nsf + 179);
    const auto *nsf_180 = buffer.data(nsf + 180);
    const auto *nsf_182 = buffer.data(nsf + 182);
    const auto *nsf_183 = buffer.data(nsf + 183);
    const auto *nsf_185 = buffer.data(nsf + 185);
    const auto *nsf_186 = buffer.data(nsf + 186);
    const auto *nsf_187 = buffer.data(nsf + 187);
    const auto *nsf_188 = buffer.data(nsf + 188);
    const auto *nsf_189 = buffer.data(nsf + 189);
    const auto *nsf_190 = buffer.data(nsf + 190);
    const auto *nsf_192 = buffer.data(nsf + 192);
    const auto *nsf_196 = buffer.data(nsf + 196);
    const auto *nsf_197 = buffer.data(nsf + 197);
    const auto *nsf_198 = buffer.data(nsf + 198);
    const auto *nsf_199 = buffer.data(nsf + 199);
    const auto *nsf_200 = buffer.data(nsf + 200);
    const auto *nsf_201 = buffer.data(nsf + 201);
    const auto *nsf_202 = buffer.data(nsf + 202);
    const auto *nsf_205 = buffer.data(nsf + 205);
    const auto *nsf_206 = buffer.data(nsf + 206);
    const auto *nsf_207 = buffer.data(nsf + 207);
    const auto *nsf_208 = buffer.data(nsf + 208);
    const auto *nsf_209 = buffer.data(nsf + 209);
    const auto *nsf_210 = buffer.data(nsf + 210);
    const auto *nsf_211 = buffer.data(nsf + 211);
    const auto *nsf_212 = buffer.data(nsf + 212);
    const auto *nsf_213 = buffer.data(nsf + 213);
    const auto *nsf_216 = buffer.data(nsf + 216);
    const auto *nsf_217 = buffer.data(nsf + 217);
    const auto *nsf_218 = buffer.data(nsf + 218);
    const auto *nsf_219 = buffer.data(nsf + 219);
    const auto *nsf_220 = buffer.data(nsf + 220);
    const auto *nsf_222 = buffer.data(nsf + 222);
    const auto *nsf_225 = buffer.data(nsf + 225);
    const auto *nsf_226 = buffer.data(nsf + 226);
    const auto *nsf_227 = buffer.data(nsf + 227);
    const auto *nsf_228 = buffer.data(nsf + 228);
    const auto *nsf_229 = buffer.data(nsf + 229);
    const auto *nsf_230 = buffer.data(nsf + 230);
    const auto *nsf_232 = buffer.data(nsf + 232);
    const auto *nsf_233 = buffer.data(nsf + 233);
    const auto *nsf_235 = buffer.data(nsf + 235);
    const auto *nsf_236 = buffer.data(nsf + 236);
    const auto *nsf_237 = buffer.data(nsf + 237);
    const auto *nsf_238 = buffer.data(nsf + 238);
    const auto *nsf_239 = buffer.data(nsf + 239);
    const auto *nsf_240 = buffer.data(nsf + 240);
    const auto *nsf_242 = buffer.data(nsf + 242);
    const auto *nsf_243 = buffer.data(nsf + 243);
    const auto *nsf_245 = buffer.data(nsf + 245);
    const auto *nsf_246 = buffer.data(nsf + 246);
    const auto *nsf_247 = buffer.data(nsf + 247);
    const auto *nsf_248 = buffer.data(nsf + 248);
    const auto *nsf_249 = buffer.data(nsf + 249);

#pragma omp simd aligned(t_253, t_254, t_255, pc_x, pc_y, pc_z, msf_109, msf_119, msf_170, \
                         nsd0_101, nsd0_102, nsd1_101, nsd1_102, nsf_169, \
                         nsf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_16 * msf_119[k]
                   + f_3 * pc_y[k] * nsf_169[k];

        t_254[k] = f_7 * msf_109[k]
                   + f_1 * nsd0_101[k]
                   - f_2 * nsd1_101[k]
                   + f_3 * pc_z[k] * nsf_169[k];

        t_255[k] = f_17 * msf_170[k]
                   + f_1 * nsd0_102[k]
                   - f_2 * nsd1_102[k]
                   + f_3 * pc_x[k] * nsf_170[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_x, pc_y, pc_z, msf_110, msf_120, \
                         msf_122, msf_173, nsd0_105, nsd1_105, nsf_170, nsf_172, \
                         nsf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * msf_120[k]
                   + f_3 * pc_y[k] * nsf_170[k];

        t_257[k] = f_8 * msf_110[k]
                   + f_3 * pc_z[k] * nsf_170[k];

        t_258[k] = f_17 * msf_173[k]
                   + f_4 * nsd0_105[k]
                   - f_5 * nsd1_105[k]
                   + f_3 * pc_x[k] * nsf_173[k];

        t_259[k] = f_14 * msf_122[k]
                   + f_3 * pc_y[k] * nsf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, msf_175, msf_176, msf_177, msf_178, \
                         nsd0_107, nsd1_107, nsf_175, nsf_176, nsf_177, \
                         nsf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_17 * msf_175[k]
                   + f_4 * nsd0_107[k]
                   - f_5 * nsd1_107[k]
                   + f_3 * pc_x[k] * nsf_175[k];

        t_261[k] = f_17 * msf_176[k]
                   + f_3 * pc_x[k] * nsf_176[k];

        t_262[k] = f_17 * msf_177[k]
                   + f_3 * pc_x[k] * nsf_177[k];

        t_263[k] = f_17 * msf_178[k]
                   + f_3 * pc_x[k] * nsf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, pc_z, msf_116, msf_126, msf_179, \
                         nsd0_105, nsd1_105, nsf_176, nsf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * msf_179[k]
                   + f_3 * pc_x[k] * nsf_179[k];

        t_265[k] = f_14 * msf_126[k]
                   + f_1 * nsd0_105[k]
                   - f_2 * nsd1_105[k]
                   + f_3 * pc_y[k] * nsf_176[k];

        t_266[k] = f_8 * msf_116[k]
                   + f_3 * pc_z[k] * nsf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, msf_119, msf_128, msf_129, nsd0_107, \
                         nsd1_107, nsf_178, nsf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * msf_128[k]
                   + f_4 * nsd0_107[k]
                   - f_5 * nsd1_107[k]
                   + f_3 * pc_y[k] * nsf_178[k];

        t_268[k] = f_14 * msf_129[k]
                   + f_3 * pc_y[k] * nsf_179[k];

        t_269[k] = f_8 * msf_119[k]
                   + f_1 * nsd0_107[k]
                   - f_2 * nsd1_107[k]
                   + f_3 * pc_z[k] * nsf_179[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pc_x, pc_y, pc_z, msf_120, msf_130, msf_180, \
                         nsd0_108, nsd1_108, nsf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_17 * msf_180[k]
                   + f_1 * nsd0_108[k]
                   - f_2 * nsd1_108[k]
                   + f_3 * pc_x[k] * nsf_180[k];

        t_271[k] = f_8 * msf_130[k]
                   + f_3 * pc_y[k] * nsf_180[k];

        t_272[k] = f_14 * msf_120[k]
                   + f_3 * pc_z[k] * nsf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_x, pc_y, msf_132, msf_183, msf_185, nsd0_111, \
                         nsd0_113, nsd1_111, nsd1_113, nsf_182, nsf_183, \
                         nsf_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_17 * msf_183[k]
                   + f_4 * nsd0_111[k]
                   - f_5 * nsd1_111[k]
                   + f_3 * pc_x[k] * nsf_183[k];

        t_274[k] = f_8 * msf_132[k]
                   + f_3 * pc_y[k] * nsf_182[k];

        t_275[k] = f_17 * msf_185[k]
                   + f_4 * nsd0_113[k]
                   - f_5 * nsd1_113[k]
                   + f_3 * pc_x[k] * nsf_185[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, msf_186, msf_187, msf_188, msf_189, \
                         nsf_186, nsf_187, nsf_188, nsf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_17 * msf_186[k]
                   + f_3 * pc_x[k] * nsf_186[k];

        t_277[k] = f_17 * msf_187[k]
                   + f_3 * pc_x[k] * nsf_187[k];

        t_278[k] = f_17 * msf_188[k]
                   + f_3 * pc_x[k] * nsf_188[k];

        t_279[k] = f_17 * msf_189[k]
                   + f_3 * pc_x[k] * nsf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pc_y, pc_z, msf_126, msf_136, msf_138, nsd0_111, \
                         nsd0_113, nsd1_111, nsd1_113, nsf_186, \
                         nsf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_8 * msf_136[k]
                   + f_1 * nsd0_111[k]
                   - f_2 * nsd1_111[k]
                   + f_3 * pc_y[k] * nsf_186[k];

        t_281[k] = f_14 * msf_126[k]
                   + f_3 * pc_z[k] * nsf_186[k];

        t_282[k] = f_8 * msf_138[k]
                   + f_4 * nsd0_113[k]
                   - f_5 * nsd1_113[k]
                   + f_3 * pc_y[k] * nsf_188[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_y, pc_y, pc_z, msg0_210, msf_129, \
                         msf_139, msf_140, msg1_210, nsd0_113, nsd1_113, nsf_189, \
                         nsf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * msf_139[k]
                   + f_3 * pc_y[k] * nsf_189[k];

        t_284[k] = f_14 * msf_129[k]
                   + f_1 * nsd0_113[k]
                   - f_2 * nsd1_113[k]
                   + f_3 * pc_z[k] * nsf_189[k];

        t_285[k] = pa_y[k] * msg0_210[k]
                   - f_6 * pc_y[k] * msg1_210[k];

        t_286[k] = f_7 * msf_140[k]
                   + f_3 * pc_y[k] * nsf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pa_y, pc_y, pc_z, msg0_213, msg0_215, \
                         msf_130, msf_141, msf_142, msg1_213, msg1_215, nsf_190, \
                         nsf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_16 * msf_130[k]
                   + f_3 * pc_z[k] * nsf_190[k];

        t_288[k] = pa_y[k] * msg0_213[k]
                   + f_8 * msf_141[k]
                   - f_6 * pc_y[k] * msg1_213[k];

        t_289[k] = f_7 * msf_142[k]
                   + f_3 * pc_y[k] * nsf_192[k];

        t_290[k] = pa_y[k] * msg0_215[k]
                   - f_6 * pc_y[k] * msg1_215[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, msf_196, msf_197, msf_198, msf_199, \
                         nsf_196, nsf_197, nsf_198, nsf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_17 * msf_196[k]
                   + f_3 * pc_x[k] * nsf_196[k];

        t_292[k] = f_17 * msf_197[k]
                   + f_3 * pc_x[k] * nsf_197[k];

        t_293[k] = f_17 * msf_198[k]
                   + f_3 * pc_x[k] * nsf_198[k];

        t_294[k] = f_17 * msf_199[k]
                   + f_3 * pc_x[k] * nsf_199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, pc_y, pc_z, msf_136, msf_146, msf_148, nsd0_117, \
                         nsd0_119, nsd1_117, nsd1_119, nsf_196, \
                         nsf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_7 * msf_146[k]
                   + f_1 * nsd0_117[k]
                   - f_2 * nsd1_117[k]
                   + f_3 * pc_y[k] * nsf_196[k];

        t_296[k] = f_16 * msf_136[k]
                   + f_3 * pc_z[k] * nsf_196[k];

        t_297[k] = f_7 * msf_148[k]
                   + f_4 * nsd0_119[k]
                   - f_5 * nsd1_119[k]
                   + f_3 * pc_y[k] * nsf_198[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pc_x, pc_y, msg0_224, msf_149, \
                         msf_200, msg1_224, nsd0_120, nsd1_120, nsf_199, \
                         nsf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * msf_149[k]
                   + f_3 * pc_y[k] * nsf_199[k];

        t_299[k] = pa_y[k] * msg0_224[k]
                   - f_6 * pc_y[k] * msg1_224[k];

        t_300[k] = f_17 * msf_200[k]
                   + f_1 * nsd0_120[k]
                   - f_2 * nsd1_120[k]
                   + f_3 * pc_x[k] * nsf_200[k];

        t_301[k] = f_3 * pc_y[k] * nsf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pc_y, pc_z, msf_140, nsd0_120, nsd1_120, \
                         nsf_200, nsf_201, nsf_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_17 * msf_140[k]
                   + f_3 * pc_z[k] * nsf_200[k];

        t_303[k] = f_4 * nsd0_120[k]
                   - f_5 * nsd1_120[k]
                   + f_3 * pc_y[k] * nsf_201[k];

        t_304[k] = f_3 * pc_y[k] * nsf_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, pc_y, msf_205, msf_206, msf_207, \
                         nsd0_125, nsd1_125, nsf_205, nsf_206, \
                         nsf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_17 * msf_205[k]
                   + f_4 * nsd0_125[k]
                   - f_5 * nsd1_125[k]
                   + f_3 * pc_x[k] * nsf_205[k];

        t_306[k] = f_17 * msf_206[k]
                   + f_3 * pc_x[k] * nsf_206[k];

        t_307[k] = f_17 * msf_207[k]
                   + f_3 * pc_x[k] * nsf_207[k];

        t_308[k] = f_3 * pc_y[k] * nsf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_x, pc_y, msf_209, nsd0_123, nsd0_124, \
                         nsd1_123, nsd1_124, nsf_206, nsf_207, \
                         nsf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_17 * msf_209[k]
                   + f_3 * pc_x[k] * nsf_209[k];

        t_310[k] = f_1 * nsd0_123[k]
                   - f_2 * nsd1_123[k]
                   + f_3 * pc_y[k] * nsf_206[k];

        t_311[k] = f_10 * nsd0_124[k]
                   - f_11 * nsd1_124[k]
                   + f_3 * pc_y[k] * nsf_207[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, pc_z, msf_149, msf_210, \
                         nsd0_125, nsd0_126, nsd1_125, nsd1_126, nsf_208, nsf_209, \
                         nsf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * nsd0_125[k]
                   - f_5 * nsd1_125[k]
                   + f_3 * pc_y[k] * nsf_208[k];

        t_313[k] = f_3 * pc_y[k] * nsf_209[k];

        t_314[k] = f_17 * msf_149[k]
                   + f_1 * nsd0_125[k]
                   - f_2 * nsd1_125[k]
                   + f_3 * pc_z[k] * nsf_209[k];

        t_315[k] = f_16 * msf_210[k]
                   + f_1 * nsd0_126[k]
                   - f_2 * nsd1_126[k]
                   + f_3 * pc_x[k] * nsf_210[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_x, pc_y, pc_z, msf_150, msf_213, \
                         nsd0_129, nsd1_129, nsf_210, nsf_211, \
                         nsf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_15 * msf_150[k]
                   + f_3 * pc_y[k] * nsf_210[k];

        t_317[k] = f_3 * pc_z[k] * nsf_210[k];

        t_318[k] = f_16 * msf_213[k]
                   + f_4 * nsd0_129[k]
                   - f_5 * nsd1_129[k]
                   + f_3 * pc_x[k] * nsf_213[k];

        t_319[k] = f_3 * pc_z[k] * nsf_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, pc_z, msf_216, msf_218, nsd0_126, \
                         nsd1_126, nsf_212, nsf_213, nsf_216, nsf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_4 * nsd0_126[k]
                   - f_5 * nsd1_126[k]
                   + f_3 * pc_z[k] * nsf_212[k];

        t_321[k] = f_16 * msf_216[k]
                   + f_3 * pc_x[k] * nsf_216[k];

        t_322[k] = f_3 * pc_z[k] * nsf_213[k];

        t_323[k] = f_16 * msf_218[k]
                   + f_3 * pc_x[k] * nsf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, msf_156, \
                         msf_159, msf_219, nsd0_129, nsd1_129, nsf_216, nsf_217, \
                         nsf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_16 * msf_219[k]
                   + f_3 * pc_x[k] * nsf_219[k];

        t_325[k] = f_15 * msf_156[k]
                   + f_1 * nsd0_129[k]
                   - f_2 * nsd1_129[k]
                   + f_3 * pc_y[k] * nsf_216[k];

        t_326[k] = f_3 * pc_z[k] * nsf_216[k];

        t_327[k] = f_4 * nsd0_129[k]
                   - f_5 * nsd1_129[k]
                   + f_3 * pc_z[k] * nsf_217[k];

        t_328[k] = f_15 * msf_159[k]
                   + f_3 * pc_y[k] * nsf_219[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pa_z, pc_y, pc_z, msg0_225, msf_150, \
                         msf_160, msg1_225, nsd0_131, nsd1_131, nsf_219, \
                         nsf_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_1 * nsd0_131[k]
                   - f_2 * nsd1_131[k]
                   + f_3 * pc_z[k] * nsf_219[k];

        t_330[k] = pa_z[k] * msg0_225[k]
                   - f_6 * pc_z[k] * msg1_225[k];

        t_331[k] = f_17 * msf_160[k]
                   + f_3 * pc_y[k] * nsf_220[k];

        t_332[k] = f_7 * msf_150[k]
                   + f_3 * pc_z[k] * nsf_220[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pa_z, pc_x, pc_y, pc_z, msg0_228, msf_162, \
                         msf_225, msg1_228, nsd0_137, nsd1_137, nsf_222, \
                         nsf_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pa_z[k] * msg0_228[k]
                   - f_6 * pc_z[k] * msg1_228[k];

        t_334[k] = f_17 * msf_162[k]
                   + f_3 * pc_y[k] * nsf_222[k];

        t_335[k] = f_16 * msf_225[k]
                   + f_4 * nsd0_137[k]
                   - f_5 * nsd1_137[k]
                   + f_3 * pc_x[k] * nsf_225[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pc_x, msf_226, msf_227, msf_228, msf_229, \
                         nsf_226, nsf_227, nsf_228, nsf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_16 * msf_226[k]
                   + f_3 * pc_x[k] * nsf_226[k];

        t_337[k] = f_16 * msf_227[k]
                   + f_3 * pc_x[k] * nsf_227[k];

        t_338[k] = f_16 * msf_228[k]
                   + f_3 * pc_x[k] * nsf_228[k];

        t_339[k] = f_16 * msf_229[k]
                   + f_3 * pc_x[k] * nsf_229[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pa_z, pc_y, pc_z, msg0_235, msf_156, msf_168, \
                         msg1_235, nsd0_137, nsd1_137, nsf_226, \
                         nsf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pa_z[k] * msg0_235[k]
                   - f_6 * pc_z[k] * msg1_235[k];

        t_341[k] = f_7 * msf_156[k]
                   + f_3 * pc_z[k] * nsf_226[k];

        t_342[k] = f_17 * msf_168[k]
                   + f_4 * nsd0_137[k]
                   - f_5 * nsd1_137[k]
                   + f_3 * pc_y[k] * nsf_228[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pc_x, pc_y, pc_z, msf_159, msf_169, msf_230, \
                         nsd0_137, nsd0_138, nsd1_137, nsd1_138, nsf_229, \
                         nsf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_17 * msf_169[k]
                   + f_3 * pc_y[k] * nsf_229[k];

        t_344[k] = f_7 * msf_159[k]
                   + f_1 * nsd0_137[k]
                   - f_2 * nsd1_137[k]
                   + f_3 * pc_z[k] * nsf_229[k];

        t_345[k] = f_16 * msf_230[k]
                   + f_1 * nsd0_138[k]
                   - f_2 * nsd1_138[k]
                   + f_3 * pc_x[k] * nsf_230[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pc_x, pc_y, pc_z, msf_160, msf_170, \
                         msf_172, msf_233, nsd0_141, nsd1_141, nsf_230, nsf_232, \
                         nsf_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_16 * msf_170[k]
                   + f_3 * pc_y[k] * nsf_230[k];

        t_347[k] = f_8 * msf_160[k]
                   + f_3 * pc_z[k] * nsf_230[k];

        t_348[k] = f_16 * msf_233[k]
                   + f_4 * nsd0_141[k]
                   - f_5 * nsd1_141[k]
                   + f_3 * pc_x[k] * nsf_233[k];

        t_349[k] = f_16 * msf_172[k]
                   + f_3 * pc_y[k] * nsf_232[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pc_x, msf_235, msf_236, msf_237, msf_238, \
                         nsd0_143, nsd1_143, nsf_235, nsf_236, nsf_237, \
                         nsf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_16 * msf_235[k]
                   + f_4 * nsd0_143[k]
                   - f_5 * nsd1_143[k]
                   + f_3 * pc_x[k] * nsf_235[k];

        t_351[k] = f_16 * msf_236[k]
                   + f_3 * pc_x[k] * nsf_236[k];

        t_352[k] = f_16 * msf_237[k]
                   + f_3 * pc_x[k] * nsf_237[k];

        t_353[k] = f_16 * msf_238[k]
                   + f_3 * pc_x[k] * nsf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_x, pc_y, pc_z, msf_166, msf_176, msf_239, \
                         nsd0_141, nsd1_141, nsf_236, nsf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_16 * msf_239[k]
                   + f_3 * pc_x[k] * nsf_239[k];

        t_355[k] = f_16 * msf_176[k]
                   + f_1 * nsd0_141[k]
                   - f_2 * nsd1_141[k]
                   + f_3 * pc_y[k] * nsf_236[k];

        t_356[k] = f_8 * msf_166[k]
                   + f_3 * pc_z[k] * nsf_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pc_y, pc_z, msf_169, msf_178, msf_179, nsd0_143, \
                         nsd1_143, nsf_238, nsf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_16 * msf_178[k]
                   + f_4 * nsd0_143[k]
                   - f_5 * nsd1_143[k]
                   + f_3 * pc_y[k] * nsf_238[k];

        t_358[k] = f_16 * msf_179[k]
                   + f_3 * pc_y[k] * nsf_239[k];

        t_359[k] = f_8 * msf_169[k]
                   + f_1 * nsd0_143[k]
                   - f_2 * nsd1_143[k]
                   + f_3 * pc_z[k] * nsf_239[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pc_x, pc_y, pc_z, msf_170, msf_180, msf_240, \
                         nsd0_144, nsd1_144, nsf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_16 * msf_240[k]
                   + f_1 * nsd0_144[k]
                   - f_2 * nsd1_144[k]
                   + f_3 * pc_x[k] * nsf_240[k];

        t_361[k] = f_14 * msf_180[k]
                   + f_3 * pc_y[k] * nsf_240[k];

        t_362[k] = f_14 * msf_170[k]
                   + f_3 * pc_z[k] * nsf_240[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pc_x, pc_y, msf_182, msf_243, msf_245, nsd0_147, \
                         nsd0_149, nsd1_147, nsd1_149, nsf_242, nsf_243, \
                         nsf_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_16 * msf_243[k]
                   + f_4 * nsd0_147[k]
                   - f_5 * nsd1_147[k]
                   + f_3 * pc_x[k] * nsf_243[k];

        t_364[k] = f_14 * msf_182[k]
                   + f_3 * pc_y[k] * nsf_242[k];

        t_365[k] = f_16 * msf_245[k]
                   + f_4 * nsd0_149[k]
                   - f_5 * nsd1_149[k]
                   + f_3 * pc_x[k] * nsf_245[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pc_x, msf_246, msf_247, msf_248, msf_249, \
                         nsf_246, nsf_247, nsf_248, nsf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_16 * msf_246[k]
                   + f_3 * pc_x[k] * nsf_246[k];

        t_367[k] = f_16 * msf_247[k]
                   + f_3 * pc_x[k] * nsf_247[k];

        t_368[k] = f_16 * msf_248[k]
                   + f_3 * pc_x[k] * nsf_248[k];

        t_369[k] = f_16 * msf_249[k]
                   + f_3 * pc_x[k] * nsf_249[k];
    }
}

static auto
compute_prim_nsg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msg0,
                                                          const size_t msf, const size_t msg1,
                                                          const size_t nsd0, const size_t nsd1,
                                                          const size_t nsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_13 = 3.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msg0_300 = buffer.data(msg0 + 300);
    const auto *msg0_303 = buffer.data(msg0 + 303);
    const auto *msg0_305 = buffer.data(msg0 + 305);
    const auto *msg0_314 = buffer.data(msg0 + 314);
    const auto *msg0_315 = buffer.data(msg0 + 315);
    const auto *msg0_318 = buffer.data(msg0 + 318);
    const auto *msg0_325 = buffer.data(msg0 + 325);

    const auto *msf_176 = buffer.data(msf + 176);
    const auto *msf_179 = buffer.data(msf + 179);
    const auto *msf_180 = buffer.data(msf + 180);
    const auto *msf_186 = buffer.data(msf + 186);
    const auto *msf_188 = buffer.data(msf + 188);
    const auto *msf_189 = buffer.data(msf + 189);
    const auto *msf_190 = buffer.data(msf + 190);
    const auto *msf_192 = buffer.data(msf + 192);
    const auto *msf_196 = buffer.data(msf + 196);
    const auto *msf_198 = buffer.data(msf + 198);
    const auto *msf_199 = buffer.data(msf + 199);
    const auto *msf_200 = buffer.data(msf + 200);
    const auto *msf_201 = buffer.data(msf + 201);
    const auto *msf_202 = buffer.data(msf + 202);
    const auto *msf_206 = buffer.data(msf + 206);
    const auto *msf_208 = buffer.data(msf + 208);
    const auto *msf_209 = buffer.data(msf + 209);
    const auto *msf_210 = buffer.data(msf + 210);
    const auto *msf_216 = buffer.data(msf + 216);
    const auto *msf_219 = buffer.data(msf + 219);
    const auto *msf_220 = buffer.data(msf + 220);
    const auto *msf_222 = buffer.data(msf + 222);
    const auto *msf_226 = buffer.data(msf + 226);
    const auto *msf_228 = buffer.data(msf + 228);
    const auto *msf_229 = buffer.data(msf + 229);
    const auto *msf_230 = buffer.data(msf + 230);
    const auto *msf_232 = buffer.data(msf + 232);
    const auto *msf_236 = buffer.data(msf + 236);
    const auto *msf_238 = buffer.data(msf + 238);
    const auto *msf_239 = buffer.data(msf + 239);
    const auto *msf_240 = buffer.data(msf + 240);
    const auto *msf_242 = buffer.data(msf + 242);
    const auto *msf_246 = buffer.data(msf + 246);
    const auto *msf_248 = buffer.data(msf + 248);
    const auto *msf_249 = buffer.data(msf + 249);
    const auto *msf_250 = buffer.data(msf + 250);
    const auto *msf_252 = buffer.data(msf + 252);
    const auto *msf_253 = buffer.data(msf + 253);
    const auto *msf_255 = buffer.data(msf + 255);
    const auto *msf_256 = buffer.data(msf + 256);
    const auto *msf_257 = buffer.data(msf + 257);
    const auto *msf_258 = buffer.data(msf + 258);
    const auto *msf_259 = buffer.data(msf + 259);
    const auto *msf_266 = buffer.data(msf + 266);
    const auto *msf_267 = buffer.data(msf + 267);
    const auto *msf_268 = buffer.data(msf + 268);
    const auto *msf_269 = buffer.data(msf + 269);
    const auto *msf_270 = buffer.data(msf + 270);
    const auto *msf_275 = buffer.data(msf + 275);
    const auto *msf_276 = buffer.data(msf + 276);
    const auto *msf_277 = buffer.data(msf + 277);
    const auto *msf_279 = buffer.data(msf + 279);
    const auto *msf_280 = buffer.data(msf + 280);
    const auto *msf_283 = buffer.data(msf + 283);
    const auto *msf_286 = buffer.data(msf + 286);
    const auto *msf_288 = buffer.data(msf + 288);
    const auto *msf_289 = buffer.data(msf + 289);
    const auto *msf_295 = buffer.data(msf + 295);
    const auto *msf_296 = buffer.data(msf + 296);
    const auto *msf_297 = buffer.data(msf + 297);
    const auto *msf_298 = buffer.data(msf + 298);
    const auto *msf_299 = buffer.data(msf + 299);
    const auto *msf_300 = buffer.data(msf + 300);
    const auto *msf_303 = buffer.data(msf + 303);
    const auto *msf_305 = buffer.data(msf + 305);
    const auto *msf_306 = buffer.data(msf + 306);
    const auto *msf_307 = buffer.data(msf + 307);
    const auto *msf_308 = buffer.data(msf + 308);
    const auto *msf_309 = buffer.data(msf + 309);
    const auto *msf_310 = buffer.data(msf + 310);
    const auto *msf_313 = buffer.data(msf + 313);
    const auto *msf_315 = buffer.data(msf + 315);
    const auto *msf_316 = buffer.data(msf + 316);
    const auto *msf_317 = buffer.data(msf + 317);
    const auto *msf_318 = buffer.data(msf + 318);
    const auto *msf_319 = buffer.data(msf + 319);
    const auto *msf_320 = buffer.data(msf + 320);
    const auto *msf_323 = buffer.data(msf + 323);

    const auto *msg1_300 = buffer.data(msg1 + 300);
    const auto *msg1_303 = buffer.data(msg1 + 303);
    const auto *msg1_305 = buffer.data(msg1 + 305);
    const auto *msg1_314 = buffer.data(msg1 + 314);
    const auto *msg1_315 = buffer.data(msg1 + 315);
    const auto *msg1_318 = buffer.data(msg1 + 318);
    const auto *msg1_325 = buffer.data(msg1 + 325);

    const auto *nsd0_147 = buffer.data(nsd0 + 147);
    const auto *nsd0_149 = buffer.data(nsd0 + 149);
    const auto *nsd0_150 = buffer.data(nsd0 + 150);
    const auto *nsd0_153 = buffer.data(nsd0 + 153);
    const auto *nsd0_155 = buffer.data(nsd0 + 155);
    const auto *nsd0_159 = buffer.data(nsd0 + 159);
    const auto *nsd0_161 = buffer.data(nsd0 + 161);
    const auto *nsd0_162 = buffer.data(nsd0 + 162);
    const auto *nsd0_165 = buffer.data(nsd0 + 165);
    const auto *nsd0_166 = buffer.data(nsd0 + 166);
    const auto *nsd0_167 = buffer.data(nsd0 + 167);
    const auto *nsd0_168 = buffer.data(nsd0 + 168);
    const auto *nsd0_171 = buffer.data(nsd0 + 171);
    const auto *nsd0_173 = buffer.data(nsd0 + 173);
    const auto *nsd0_179 = buffer.data(nsd0 + 179);
    const auto *nsd0_180 = buffer.data(nsd0 + 180);
    const auto *nsd0_183 = buffer.data(nsd0 + 183);
    const auto *nsd0_185 = buffer.data(nsd0 + 185);
    const auto *nsd0_186 = buffer.data(nsd0 + 186);
    const auto *nsd0_189 = buffer.data(nsd0 + 189);
    const auto *nsd0_191 = buffer.data(nsd0 + 191);
    const auto *nsd0_192 = buffer.data(nsd0 + 192);
    const auto *nsd0_195 = buffer.data(nsd0 + 195);

    const auto *nsd1_147 = buffer.data(nsd1 + 147);
    const auto *nsd1_149 = buffer.data(nsd1 + 149);
    const auto *nsd1_150 = buffer.data(nsd1 + 150);
    const auto *nsd1_153 = buffer.data(nsd1 + 153);
    const auto *nsd1_155 = buffer.data(nsd1 + 155);
    const auto *nsd1_159 = buffer.data(nsd1 + 159);
    const auto *nsd1_161 = buffer.data(nsd1 + 161);
    const auto *nsd1_162 = buffer.data(nsd1 + 162);
    const auto *nsd1_165 = buffer.data(nsd1 + 165);
    const auto *nsd1_166 = buffer.data(nsd1 + 166);
    const auto *nsd1_167 = buffer.data(nsd1 + 167);
    const auto *nsd1_168 = buffer.data(nsd1 + 168);
    const auto *nsd1_171 = buffer.data(nsd1 + 171);
    const auto *nsd1_173 = buffer.data(nsd1 + 173);
    const auto *nsd1_179 = buffer.data(nsd1 + 179);
    const auto *nsd1_180 = buffer.data(nsd1 + 180);
    const auto *nsd1_183 = buffer.data(nsd1 + 183);
    const auto *nsd1_185 = buffer.data(nsd1 + 185);
    const auto *nsd1_186 = buffer.data(nsd1 + 186);
    const auto *nsd1_189 = buffer.data(nsd1 + 189);
    const auto *nsd1_191 = buffer.data(nsd1 + 191);
    const auto *nsd1_192 = buffer.data(nsd1 + 192);
    const auto *nsd1_195 = buffer.data(nsd1 + 195);

    const auto *nsf_246 = buffer.data(nsf + 246);
    const auto *nsf_248 = buffer.data(nsf + 248);
    const auto *nsf_249 = buffer.data(nsf + 249);
    const auto *nsf_250 = buffer.data(nsf + 250);
    const auto *nsf_252 = buffer.data(nsf + 252);
    const auto *nsf_253 = buffer.data(nsf + 253);
    const auto *nsf_255 = buffer.data(nsf + 255);
    const auto *nsf_256 = buffer.data(nsf + 256);
    const auto *nsf_257 = buffer.data(nsf + 257);
    const auto *nsf_258 = buffer.data(nsf + 258);
    const auto *nsf_259 = buffer.data(nsf + 259);
    const auto *nsf_260 = buffer.data(nsf + 260);
    const auto *nsf_262 = buffer.data(nsf + 262);
    const auto *nsf_266 = buffer.data(nsf + 266);
    const auto *nsf_267 = buffer.data(nsf + 267);
    const auto *nsf_268 = buffer.data(nsf + 268);
    const auto *nsf_269 = buffer.data(nsf + 269);
    const auto *nsf_270 = buffer.data(nsf + 270);
    const auto *nsf_271 = buffer.data(nsf + 271);
    const auto *nsf_272 = buffer.data(nsf + 272);
    const auto *nsf_275 = buffer.data(nsf + 275);
    const auto *nsf_276 = buffer.data(nsf + 276);
    const auto *nsf_277 = buffer.data(nsf + 277);
    const auto *nsf_278 = buffer.data(nsf + 278);
    const auto *nsf_279 = buffer.data(nsf + 279);
    const auto *nsf_280 = buffer.data(nsf + 280);
    const auto *nsf_281 = buffer.data(nsf + 281);
    const auto *nsf_282 = buffer.data(nsf + 282);
    const auto *nsf_283 = buffer.data(nsf + 283);
    const auto *nsf_286 = buffer.data(nsf + 286);
    const auto *nsf_287 = buffer.data(nsf + 287);
    const auto *nsf_288 = buffer.data(nsf + 288);
    const auto *nsf_289 = buffer.data(nsf + 289);
    const auto *nsf_290 = buffer.data(nsf + 290);
    const auto *nsf_292 = buffer.data(nsf + 292);
    const auto *nsf_295 = buffer.data(nsf + 295);
    const auto *nsf_296 = buffer.data(nsf + 296);
    const auto *nsf_297 = buffer.data(nsf + 297);
    const auto *nsf_298 = buffer.data(nsf + 298);
    const auto *nsf_299 = buffer.data(nsf + 299);
    const auto *nsf_300 = buffer.data(nsf + 300);
    const auto *nsf_302 = buffer.data(nsf + 302);
    const auto *nsf_303 = buffer.data(nsf + 303);
    const auto *nsf_305 = buffer.data(nsf + 305);
    const auto *nsf_306 = buffer.data(nsf + 306);
    const auto *nsf_307 = buffer.data(nsf + 307);
    const auto *nsf_308 = buffer.data(nsf + 308);
    const auto *nsf_309 = buffer.data(nsf + 309);
    const auto *nsf_310 = buffer.data(nsf + 310);
    const auto *nsf_312 = buffer.data(nsf + 312);
    const auto *nsf_313 = buffer.data(nsf + 313);
    const auto *nsf_315 = buffer.data(nsf + 315);
    const auto *nsf_316 = buffer.data(nsf + 316);
    const auto *nsf_317 = buffer.data(nsf + 317);
    const auto *nsf_318 = buffer.data(nsf + 318);
    const auto *nsf_319 = buffer.data(nsf + 319);
    const auto *nsf_320 = buffer.data(nsf + 320);
    const auto *nsf_322 = buffer.data(nsf + 322);
    const auto *nsf_323 = buffer.data(nsf + 323);

#pragma omp simd aligned(t_370, t_371, t_372, pc_y, pc_z, msf_176, msf_186, msf_188, nsd0_147, \
                         nsd0_149, nsd1_147, nsd1_149, nsf_246, \
                         nsf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_14 * msf_186[k]
                   + f_1 * nsd0_147[k]
                   - f_2 * nsd1_147[k]
                   + f_3 * pc_y[k] * nsf_246[k];

        t_371[k] = f_14 * msf_176[k]
                   + f_3 * pc_z[k] * nsf_246[k];

        t_372[k] = f_14 * msf_188[k]
                   + f_4 * nsd0_149[k]
                   - f_5 * nsd1_149[k]
                   + f_3 * pc_y[k] * nsf_248[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, pc_z, msf_179, msf_189, msf_250, \
                         nsd0_149, nsd0_150, nsd1_149, nsd1_150, nsf_249, \
                         nsf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_14 * msf_189[k]
                   + f_3 * pc_y[k] * nsf_249[k];

        t_374[k] = f_14 * msf_179[k]
                   + f_1 * nsd0_149[k]
                   - f_2 * nsd1_149[k]
                   + f_3 * pc_z[k] * nsf_249[k];

        t_375[k] = f_16 * msf_250[k]
                   + f_1 * nsd0_150[k]
                   - f_2 * nsd1_150[k]
                   + f_3 * pc_x[k] * nsf_250[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pc_x, pc_y, pc_z, msf_180, msf_190, \
                         msf_192, msf_253, nsd0_153, nsd1_153, nsf_250, nsf_252, \
                         nsf_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_8 * msf_190[k]
                   + f_3 * pc_y[k] * nsf_250[k];

        t_377[k] = f_16 * msf_180[k]
                   + f_3 * pc_z[k] * nsf_250[k];

        t_378[k] = f_16 * msf_253[k]
                   + f_4 * nsd0_153[k]
                   - f_5 * nsd1_153[k]
                   + f_3 * pc_x[k] * nsf_253[k];

        t_379[k] = f_8 * msf_192[k]
                   + f_3 * pc_y[k] * nsf_252[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, msf_255, msf_256, msf_257, msf_258, \
                         nsd0_155, nsd1_155, nsf_255, nsf_256, nsf_257, \
                         nsf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_16 * msf_255[k]
                   + f_4 * nsd0_155[k]
                   - f_5 * nsd1_155[k]
                   + f_3 * pc_x[k] * nsf_255[k];

        t_381[k] = f_16 * msf_256[k]
                   + f_3 * pc_x[k] * nsf_256[k];

        t_382[k] = f_16 * msf_257[k]
                   + f_3 * pc_x[k] * nsf_257[k];

        t_383[k] = f_16 * msf_258[k]
                   + f_3 * pc_x[k] * nsf_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pc_x, pc_y, pc_z, msf_186, msf_196, msf_259, \
                         nsd0_153, nsd1_153, nsf_256, nsf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_16 * msf_259[k]
                   + f_3 * pc_x[k] * nsf_259[k];

        t_385[k] = f_8 * msf_196[k]
                   + f_1 * nsd0_153[k]
                   - f_2 * nsd1_153[k]
                   + f_3 * pc_y[k] * nsf_256[k];

        t_386[k] = f_16 * msf_186[k]
                   + f_3 * pc_z[k] * nsf_256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_y, pc_y, pc_z, msg0_300, msf_189, \
                         msf_198, msf_199, msg1_300, nsd0_155, nsd1_155, nsf_258, \
                         nsf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_8 * msf_198[k]
                   + f_4 * nsd0_155[k]
                   - f_5 * nsd1_155[k]
                   + f_3 * pc_y[k] * nsf_258[k];

        t_388[k] = f_8 * msf_199[k]
                   + f_3 * pc_y[k] * nsf_259[k];

        t_389[k] = f_16 * msf_189[k]
                   + f_1 * nsd0_155[k]
                   - f_2 * nsd1_155[k]
                   + f_3 * pc_z[k] * nsf_259[k];

        t_390[k] = pa_y[k] * msg0_300[k]
                   - f_6 * pc_y[k] * msg1_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_y, pc_z, msg0_303, msf_190, \
                         msf_200, msf_201, msf_202, msg1_303, nsf_260, \
                         nsf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_7 * msf_200[k]
                   + f_3 * pc_y[k] * nsf_260[k];

        t_392[k] = f_17 * msf_190[k]
                   + f_3 * pc_z[k] * nsf_260[k];

        t_393[k] = pa_y[k] * msg0_303[k]
                   + f_8 * msf_201[k]
                   - f_6 * pc_y[k] * msg1_303[k];

        t_394[k] = f_7 * msf_202[k]
                   + f_3 * pc_y[k] * nsf_262[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_y, pc_x, pc_y, msg0_305, msf_266, \
                         msf_267, msf_268, msg1_305, nsf_266, nsf_267, \
                         nsf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_y[k] * msg0_305[k]
                   - f_6 * pc_y[k] * msg1_305[k];

        t_396[k] = f_16 * msf_266[k]
                   + f_3 * pc_x[k] * nsf_266[k];

        t_397[k] = f_16 * msf_267[k]
                   + f_3 * pc_x[k] * nsf_267[k];

        t_398[k] = f_16 * msf_268[k]
                   + f_3 * pc_x[k] * nsf_268[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, msf_196, msf_206, msf_269, \
                         nsd0_159, nsd1_159, nsf_266, nsf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_16 * msf_269[k]
                   + f_3 * pc_x[k] * nsf_269[k];

        t_400[k] = f_7 * msf_206[k]
                   + f_1 * nsd0_159[k]
                   - f_2 * nsd1_159[k]
                   + f_3 * pc_y[k] * nsf_266[k];

        t_401[k] = f_17 * msf_196[k]
                   + f_3 * pc_z[k] * nsf_266[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_y, pc_y, msg0_314, msf_208, msf_209, \
                         msg1_314, nsd0_161, nsd1_161, nsf_268, \
                         nsf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_7 * msf_208[k]
                   + f_4 * nsd0_161[k]
                   - f_5 * nsd1_161[k]
                   + f_3 * pc_y[k] * nsf_268[k];

        t_403[k] = f_7 * msf_209[k]
                   + f_3 * pc_y[k] * nsf_269[k];

        t_404[k] = pa_y[k] * msg0_314[k]
                   - f_6 * pc_y[k] * msg1_314[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pc_x, pc_y, pc_z, msf_200, \
                         msf_270, nsd0_162, nsd1_162, nsf_270, nsf_271, \
                         nsf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_16 * msf_270[k]
                   + f_1 * nsd0_162[k]
                   - f_2 * nsd1_162[k]
                   + f_3 * pc_x[k] * nsf_270[k];

        t_406[k] = f_3 * pc_y[k] * nsf_270[k];

        t_407[k] = f_15 * msf_200[k]
                   + f_3 * pc_z[k] * nsf_270[k];

        t_408[k] = f_4 * nsd0_162[k]
                   - f_5 * nsd1_162[k]
                   + f_3 * pc_y[k] * nsf_271[k];

        t_409[k] = f_3 * pc_y[k] * nsf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, msf_275, msf_276, msf_277, \
                         nsd0_167, nsd1_167, nsf_275, nsf_276, \
                         nsf_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_16 * msf_275[k]
                   + f_4 * nsd0_167[k]
                   - f_5 * nsd1_167[k]
                   + f_3 * pc_x[k] * nsf_275[k];

        t_411[k] = f_16 * msf_276[k]
                   + f_3 * pc_x[k] * nsf_276[k];

        t_412[k] = f_16 * msf_277[k]
                   + f_3 * pc_x[k] * nsf_277[k];

        t_413[k] = f_3 * pc_y[k] * nsf_275[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, msf_279, nsd0_165, nsd0_166, \
                         nsd1_165, nsd1_166, nsf_276, nsf_277, \
                         nsf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_16 * msf_279[k]
                   + f_3 * pc_x[k] * nsf_279[k];

        t_415[k] = f_1 * nsd0_165[k]
                   - f_2 * nsd1_165[k]
                   + f_3 * pc_y[k] * nsf_276[k];

        t_416[k] = f_10 * nsd0_166[k]
                   - f_11 * nsd1_166[k]
                   + f_3 * pc_y[k] * nsf_277[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, pc_x, pc_y, pc_z, msf_209, msf_280, \
                         nsd0_167, nsd0_168, nsd1_167, nsd1_168, nsf_278, nsf_279, \
                         nsf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * nsd0_167[k]
                   - f_5 * nsd1_167[k]
                   + f_3 * pc_y[k] * nsf_278[k];

        t_418[k] = f_3 * pc_y[k] * nsf_279[k];

        t_419[k] = f_15 * msf_209[k]
                   + f_1 * nsd0_167[k]
                   - f_2 * nsd1_167[k]
                   + f_3 * pc_z[k] * nsf_279[k];

        t_420[k] = f_14 * msf_280[k]
                   + f_1 * nsd0_168[k]
                   - f_2 * nsd1_168[k]
                   + f_3 * pc_x[k] * nsf_280[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, msf_210, msf_283, \
                         nsd0_171, nsd1_171, nsf_280, nsf_281, \
                         nsf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_13 * msf_210[k]
                   + f_3 * pc_y[k] * nsf_280[k];

        t_422[k] = f_3 * pc_z[k] * nsf_280[k];

        t_423[k] = f_14 * msf_283[k]
                   + f_4 * nsd0_171[k]
                   - f_5 * nsd1_171[k]
                   + f_3 * pc_x[k] * nsf_283[k];

        t_424[k] = f_3 * pc_z[k] * nsf_281[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_z, msf_286, msf_288, nsd0_168, \
                         nsd1_168, nsf_282, nsf_283, nsf_286, nsf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_4 * nsd0_168[k]
                   - f_5 * nsd1_168[k]
                   + f_3 * pc_z[k] * nsf_282[k];

        t_426[k] = f_14 * msf_286[k]
                   + f_3 * pc_x[k] * nsf_286[k];

        t_427[k] = f_3 * pc_z[k] * nsf_283[k];

        t_428[k] = f_14 * msf_288[k]
                   + f_3 * pc_x[k] * nsf_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, pc_x, pc_y, pc_z, msf_216, \
                         msf_219, msf_289, nsd0_171, nsd1_171, nsf_286, nsf_287, \
                         nsf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_14 * msf_289[k]
                   + f_3 * pc_x[k] * nsf_289[k];

        t_430[k] = f_13 * msf_216[k]
                   + f_1 * nsd0_171[k]
                   - f_2 * nsd1_171[k]
                   + f_3 * pc_y[k] * nsf_286[k];

        t_431[k] = f_3 * pc_z[k] * nsf_286[k];

        t_432[k] = f_4 * nsd0_171[k]
                   - f_5 * nsd1_171[k]
                   + f_3 * pc_z[k] * nsf_287[k];

        t_433[k] = f_13 * msf_219[k]
                   + f_3 * pc_y[k] * nsf_289[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pa_z, pc_y, pc_z, msg0_315, msf_210, \
                         msf_220, msg1_315, nsd0_173, nsd1_173, nsf_289, \
                         nsf_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_1 * nsd0_173[k]
                   - f_2 * nsd1_173[k]
                   + f_3 * pc_z[k] * nsf_289[k];

        t_435[k] = pa_z[k] * msg0_315[k]
                   - f_6 * pc_z[k] * msg1_315[k];

        t_436[k] = f_15 * msf_220[k]
                   + f_3 * pc_y[k] * nsf_290[k];

        t_437[k] = f_7 * msf_210[k]
                   + f_3 * pc_z[k] * nsf_290[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pa_z, pc_x, pc_y, pc_z, msg0_318, msf_222, \
                         msf_295, msg1_318, nsd0_179, nsd1_179, nsf_292, \
                         nsf_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = pa_z[k] * msg0_318[k]
                   - f_6 * pc_z[k] * msg1_318[k];

        t_439[k] = f_15 * msf_222[k]
                   + f_3 * pc_y[k] * nsf_292[k];

        t_440[k] = f_14 * msf_295[k]
                   + f_4 * nsd0_179[k]
                   - f_5 * nsd1_179[k]
                   + f_3 * pc_x[k] * nsf_295[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_x, msf_296, msf_297, msf_298, msf_299, \
                         nsf_296, nsf_297, nsf_298, nsf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_14 * msf_296[k]
                   + f_3 * pc_x[k] * nsf_296[k];

        t_442[k] = f_14 * msf_297[k]
                   + f_3 * pc_x[k] * nsf_297[k];

        t_443[k] = f_14 * msf_298[k]
                   + f_3 * pc_x[k] * nsf_298[k];

        t_444[k] = f_14 * msf_299[k]
                   + f_3 * pc_x[k] * nsf_299[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, pa_z, pc_y, pc_z, msg0_325, msf_216, msf_228, \
                         msg1_325, nsd0_179, nsd1_179, nsf_296, \
                         nsf_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pa_z[k] * msg0_325[k]
                   - f_6 * pc_z[k] * msg1_325[k];

        t_446[k] = f_7 * msf_216[k]
                   + f_3 * pc_z[k] * nsf_296[k];

        t_447[k] = f_15 * msf_228[k]
                   + f_4 * nsd0_179[k]
                   - f_5 * nsd1_179[k]
                   + f_3 * pc_y[k] * nsf_298[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, pc_x, pc_y, pc_z, msf_219, msf_229, msf_300, \
                         nsd0_179, nsd0_180, nsd1_179, nsd1_180, nsf_299, \
                         nsf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_15 * msf_229[k]
                   + f_3 * pc_y[k] * nsf_299[k];

        t_449[k] = f_7 * msf_219[k]
                   + f_1 * nsd0_179[k]
                   - f_2 * nsd1_179[k]
                   + f_3 * pc_z[k] * nsf_299[k];

        t_450[k] = f_14 * msf_300[k]
                   + f_1 * nsd0_180[k]
                   - f_2 * nsd1_180[k]
                   + f_3 * pc_x[k] * nsf_300[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pc_x, pc_y, pc_z, msf_220, msf_230, \
                         msf_232, msf_303, nsd0_183, nsd1_183, nsf_300, nsf_302, \
                         nsf_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_17 * msf_230[k]
                   + f_3 * pc_y[k] * nsf_300[k];

        t_452[k] = f_8 * msf_220[k]
                   + f_3 * pc_z[k] * nsf_300[k];

        t_453[k] = f_14 * msf_303[k]
                   + f_4 * nsd0_183[k]
                   - f_5 * nsd1_183[k]
                   + f_3 * pc_x[k] * nsf_303[k];

        t_454[k] = f_17 * msf_232[k]
                   + f_3 * pc_y[k] * nsf_302[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pc_x, msf_305, msf_306, msf_307, msf_308, \
                         nsd0_185, nsd1_185, nsf_305, nsf_306, nsf_307, \
                         nsf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_14 * msf_305[k]
                   + f_4 * nsd0_185[k]
                   - f_5 * nsd1_185[k]
                   + f_3 * pc_x[k] * nsf_305[k];

        t_456[k] = f_14 * msf_306[k]
                   + f_3 * pc_x[k] * nsf_306[k];

        t_457[k] = f_14 * msf_307[k]
                   + f_3 * pc_x[k] * nsf_307[k];

        t_458[k] = f_14 * msf_308[k]
                   + f_3 * pc_x[k] * nsf_308[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pc_x, pc_y, pc_z, msf_226, msf_236, msf_309, \
                         nsd0_183, nsd1_183, nsf_306, nsf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_14 * msf_309[k]
                   + f_3 * pc_x[k] * nsf_309[k];

        t_460[k] = f_17 * msf_236[k]
                   + f_1 * nsd0_183[k]
                   - f_2 * nsd1_183[k]
                   + f_3 * pc_y[k] * nsf_306[k];

        t_461[k] = f_8 * msf_226[k]
                   + f_3 * pc_z[k] * nsf_306[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, msf_229, msf_238, msf_239, nsd0_185, \
                         nsd1_185, nsf_308, nsf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_17 * msf_238[k]
                   + f_4 * nsd0_185[k]
                   - f_5 * nsd1_185[k]
                   + f_3 * pc_y[k] * nsf_308[k];

        t_463[k] = f_17 * msf_239[k]
                   + f_3 * pc_y[k] * nsf_309[k];

        t_464[k] = f_8 * msf_229[k]
                   + f_1 * nsd0_185[k]
                   - f_2 * nsd1_185[k]
                   + f_3 * pc_z[k] * nsf_309[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_y, pc_z, msf_230, msf_240, msf_310, \
                         nsd0_186, nsd1_186, nsf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_14 * msf_310[k]
                   + f_1 * nsd0_186[k]
                   - f_2 * nsd1_186[k]
                   + f_3 * pc_x[k] * nsf_310[k];

        t_466[k] = f_16 * msf_240[k]
                   + f_3 * pc_y[k] * nsf_310[k];

        t_467[k] = f_14 * msf_230[k]
                   + f_3 * pc_z[k] * nsf_310[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_x, pc_y, msf_242, msf_313, msf_315, nsd0_189, \
                         nsd0_191, nsd1_189, nsd1_191, nsf_312, nsf_313, \
                         nsf_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_14 * msf_313[k]
                   + f_4 * nsd0_189[k]
                   - f_5 * nsd1_189[k]
                   + f_3 * pc_x[k] * nsf_313[k];

        t_469[k] = f_16 * msf_242[k]
                   + f_3 * pc_y[k] * nsf_312[k];

        t_470[k] = f_14 * msf_315[k]
                   + f_4 * nsd0_191[k]
                   - f_5 * nsd1_191[k]
                   + f_3 * pc_x[k] * nsf_315[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, msf_316, msf_317, msf_318, msf_319, \
                         nsf_316, nsf_317, nsf_318, nsf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_14 * msf_316[k]
                   + f_3 * pc_x[k] * nsf_316[k];

        t_472[k] = f_14 * msf_317[k]
                   + f_3 * pc_x[k] * nsf_317[k];

        t_473[k] = f_14 * msf_318[k]
                   + f_3 * pc_x[k] * nsf_318[k];

        t_474[k] = f_14 * msf_319[k]
                   + f_3 * pc_x[k] * nsf_319[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, pc_y, pc_z, msf_236, msf_246, msf_248, nsd0_189, \
                         nsd0_191, nsd1_189, nsd1_191, nsf_316, \
                         nsf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_16 * msf_246[k]
                   + f_1 * nsd0_189[k]
                   - f_2 * nsd1_189[k]
                   + f_3 * pc_y[k] * nsf_316[k];

        t_476[k] = f_14 * msf_236[k]
                   + f_3 * pc_z[k] * nsf_316[k];

        t_477[k] = f_16 * msf_248[k]
                   + f_4 * nsd0_191[k]
                   - f_5 * nsd1_191[k]
                   + f_3 * pc_y[k] * nsf_318[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_x, pc_y, pc_z, msf_239, msf_249, msf_320, \
                         nsd0_191, nsd0_192, nsd1_191, nsd1_192, nsf_319, \
                         nsf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_16 * msf_249[k]
                   + f_3 * pc_y[k] * nsf_319[k];

        t_479[k] = f_14 * msf_239[k]
                   + f_1 * nsd0_191[k]
                   - f_2 * nsd1_191[k]
                   + f_3 * pc_z[k] * nsf_319[k];

        t_480[k] = f_14 * msf_320[k]
                   + f_1 * nsd0_192[k]
                   - f_2 * nsd1_192[k]
                   + f_3 * pc_x[k] * nsf_320[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, pc_z, msf_240, msf_250, \
                         msf_252, msf_323, nsd0_195, nsd1_195, nsf_320, nsf_322, \
                         nsf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_14 * msf_250[k]
                   + f_3 * pc_y[k] * nsf_320[k];

        t_482[k] = f_16 * msf_240[k]
                   + f_3 * pc_z[k] * nsf_320[k];

        t_483[k] = f_14 * msf_323[k]
                   + f_4 * nsd0_195[k]
                   - f_5 * nsd1_195[k]
                   + f_3 * pc_x[k] * nsf_323[k];

        t_484[k] = f_14 * msf_252[k]
                   + f_3 * pc_y[k] * nsf_322[k];
    }
}

static auto
compute_prim_nsg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msg0,
                                                          const size_t msf, const size_t msg1,
                                                          const size_t nsd0, const size_t nsd1,
                                                          const size_t nsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msg0_405 = buffer.data(msg0 + 405);
    const auto *msg0_408 = buffer.data(msg0 + 408);
    const auto *msg0_410 = buffer.data(msg0 + 410);
    const auto *msg0_419 = buffer.data(msg0 + 419);
    const auto *msg0_420 = buffer.data(msg0 + 420);
    const auto *msg0_423 = buffer.data(msg0 + 423);
    const auto *msg0_430 = buffer.data(msg0 + 430);

    const auto *msf_246 = buffer.data(msf + 246);
    const auto *msf_249 = buffer.data(msf + 249);
    const auto *msf_250 = buffer.data(msf + 250);
    const auto *msf_256 = buffer.data(msf + 256);
    const auto *msf_258 = buffer.data(msf + 258);
    const auto *msf_259 = buffer.data(msf + 259);
    const auto *msf_260 = buffer.data(msf + 260);
    const auto *msf_262 = buffer.data(msf + 262);
    const auto *msf_266 = buffer.data(msf + 266);
    const auto *msf_268 = buffer.data(msf + 268);
    const auto *msf_269 = buffer.data(msf + 269);
    const auto *msf_270 = buffer.data(msf + 270);
    const auto *msf_271 = buffer.data(msf + 271);
    const auto *msf_272 = buffer.data(msf + 272);
    const auto *msf_276 = buffer.data(msf + 276);
    const auto *msf_278 = buffer.data(msf + 278);
    const auto *msf_279 = buffer.data(msf + 279);
    const auto *msf_280 = buffer.data(msf + 280);
    const auto *msf_286 = buffer.data(msf + 286);
    const auto *msf_289 = buffer.data(msf + 289);
    const auto *msf_290 = buffer.data(msf + 290);
    const auto *msf_292 = buffer.data(msf + 292);
    const auto *msf_296 = buffer.data(msf + 296);
    const auto *msf_298 = buffer.data(msf + 298);
    const auto *msf_299 = buffer.data(msf + 299);
    const auto *msf_300 = buffer.data(msf + 300);
    const auto *msf_302 = buffer.data(msf + 302);
    const auto *msf_306 = buffer.data(msf + 306);
    const auto *msf_308 = buffer.data(msf + 308);
    const auto *msf_309 = buffer.data(msf + 309);
    const auto *msf_310 = buffer.data(msf + 310);
    const auto *msf_312 = buffer.data(msf + 312);
    const auto *msf_316 = buffer.data(msf + 316);
    const auto *msf_318 = buffer.data(msf + 318);
    const auto *msf_319 = buffer.data(msf + 319);
    const auto *msf_325 = buffer.data(msf + 325);
    const auto *msf_326 = buffer.data(msf + 326);
    const auto *msf_327 = buffer.data(msf + 327);
    const auto *msf_328 = buffer.data(msf + 328);
    const auto *msf_329 = buffer.data(msf + 329);
    const auto *msf_330 = buffer.data(msf + 330);
    const auto *msf_333 = buffer.data(msf + 333);
    const auto *msf_335 = buffer.data(msf + 335);
    const auto *msf_336 = buffer.data(msf + 336);
    const auto *msf_337 = buffer.data(msf + 337);
    const auto *msf_338 = buffer.data(msf + 338);
    const auto *msf_339 = buffer.data(msf + 339);
    const auto *msf_346 = buffer.data(msf + 346);
    const auto *msf_347 = buffer.data(msf + 347);
    const auto *msf_348 = buffer.data(msf + 348);
    const auto *msf_349 = buffer.data(msf + 349);
    const auto *msf_350 = buffer.data(msf + 350);
    const auto *msf_355 = buffer.data(msf + 355);
    const auto *msf_356 = buffer.data(msf + 356);
    const auto *msf_357 = buffer.data(msf + 357);
    const auto *msf_359 = buffer.data(msf + 359);
    const auto *msf_360 = buffer.data(msf + 360);
    const auto *msf_363 = buffer.data(msf + 363);
    const auto *msf_366 = buffer.data(msf + 366);
    const auto *msf_368 = buffer.data(msf + 368);
    const auto *msf_369 = buffer.data(msf + 369);
    const auto *msf_375 = buffer.data(msf + 375);
    const auto *msf_376 = buffer.data(msf + 376);
    const auto *msf_377 = buffer.data(msf + 377);
    const auto *msf_378 = buffer.data(msf + 378);
    const auto *msf_379 = buffer.data(msf + 379);
    const auto *msf_380 = buffer.data(msf + 380);
    const auto *msf_383 = buffer.data(msf + 383);
    const auto *msf_385 = buffer.data(msf + 385);
    const auto *msf_386 = buffer.data(msf + 386);
    const auto *msf_387 = buffer.data(msf + 387);
    const auto *msf_388 = buffer.data(msf + 388);
    const auto *msf_389 = buffer.data(msf + 389);
    const auto *msf_390 = buffer.data(msf + 390);
    const auto *msf_393 = buffer.data(msf + 393);
    const auto *msf_395 = buffer.data(msf + 395);
    const auto *msf_396 = buffer.data(msf + 396);
    const auto *msf_397 = buffer.data(msf + 397);
    const auto *msf_398 = buffer.data(msf + 398);
    const auto *msf_399 = buffer.data(msf + 399);
    const auto *msf_400 = buffer.data(msf + 400);

    const auto *msg1_405 = buffer.data(msg1 + 405);
    const auto *msg1_408 = buffer.data(msg1 + 408);
    const auto *msg1_410 = buffer.data(msg1 + 410);
    const auto *msg1_419 = buffer.data(msg1 + 419);
    const auto *msg1_420 = buffer.data(msg1 + 420);
    const auto *msg1_423 = buffer.data(msg1 + 423);
    const auto *msg1_430 = buffer.data(msg1 + 430);

    const auto *nsd0_195 = buffer.data(nsd0 + 195);
    const auto *nsd0_197 = buffer.data(nsd0 + 197);
    const auto *nsd0_198 = buffer.data(nsd0 + 198);
    const auto *nsd0_201 = buffer.data(nsd0 + 201);
    const auto *nsd0_203 = buffer.data(nsd0 + 203);
    const auto *nsd0_207 = buffer.data(nsd0 + 207);
    const auto *nsd0_209 = buffer.data(nsd0 + 209);
    const auto *nsd0_210 = buffer.data(nsd0 + 210);
    const auto *nsd0_213 = buffer.data(nsd0 + 213);
    const auto *nsd0_214 = buffer.data(nsd0 + 214);
    const auto *nsd0_215 = buffer.data(nsd0 + 215);
    const auto *nsd0_216 = buffer.data(nsd0 + 216);
    const auto *nsd0_219 = buffer.data(nsd0 + 219);
    const auto *nsd0_221 = buffer.data(nsd0 + 221);
    const auto *nsd0_227 = buffer.data(nsd0 + 227);
    const auto *nsd0_228 = buffer.data(nsd0 + 228);
    const auto *nsd0_231 = buffer.data(nsd0 + 231);
    const auto *nsd0_233 = buffer.data(nsd0 + 233);
    const auto *nsd0_234 = buffer.data(nsd0 + 234);
    const auto *nsd0_237 = buffer.data(nsd0 + 237);
    const auto *nsd0_239 = buffer.data(nsd0 + 239);
    const auto *nsd0_240 = buffer.data(nsd0 + 240);

    const auto *nsd1_195 = buffer.data(nsd1 + 195);
    const auto *nsd1_197 = buffer.data(nsd1 + 197);
    const auto *nsd1_198 = buffer.data(nsd1 + 198);
    const auto *nsd1_201 = buffer.data(nsd1 + 201);
    const auto *nsd1_203 = buffer.data(nsd1 + 203);
    const auto *nsd1_207 = buffer.data(nsd1 + 207);
    const auto *nsd1_209 = buffer.data(nsd1 + 209);
    const auto *nsd1_210 = buffer.data(nsd1 + 210);
    const auto *nsd1_213 = buffer.data(nsd1 + 213);
    const auto *nsd1_214 = buffer.data(nsd1 + 214);
    const auto *nsd1_215 = buffer.data(nsd1 + 215);
    const auto *nsd1_216 = buffer.data(nsd1 + 216);
    const auto *nsd1_219 = buffer.data(nsd1 + 219);
    const auto *nsd1_221 = buffer.data(nsd1 + 221);
    const auto *nsd1_227 = buffer.data(nsd1 + 227);
    const auto *nsd1_228 = buffer.data(nsd1 + 228);
    const auto *nsd1_231 = buffer.data(nsd1 + 231);
    const auto *nsd1_233 = buffer.data(nsd1 + 233);
    const auto *nsd1_234 = buffer.data(nsd1 + 234);
    const auto *nsd1_237 = buffer.data(nsd1 + 237);
    const auto *nsd1_239 = buffer.data(nsd1 + 239);
    const auto *nsd1_240 = buffer.data(nsd1 + 240);

    const auto *nsf_325 = buffer.data(nsf + 325);
    const auto *nsf_326 = buffer.data(nsf + 326);
    const auto *nsf_327 = buffer.data(nsf + 327);
    const auto *nsf_328 = buffer.data(nsf + 328);
    const auto *nsf_329 = buffer.data(nsf + 329);
    const auto *nsf_330 = buffer.data(nsf + 330);
    const auto *nsf_332 = buffer.data(nsf + 332);
    const auto *nsf_333 = buffer.data(nsf + 333);
    const auto *nsf_335 = buffer.data(nsf + 335);
    const auto *nsf_336 = buffer.data(nsf + 336);
    const auto *nsf_337 = buffer.data(nsf + 337);
    const auto *nsf_338 = buffer.data(nsf + 338);
    const auto *nsf_339 = buffer.data(nsf + 339);
    const auto *nsf_340 = buffer.data(nsf + 340);
    const auto *nsf_342 = buffer.data(nsf + 342);
    const auto *nsf_346 = buffer.data(nsf + 346);
    const auto *nsf_347 = buffer.data(nsf + 347);
    const auto *nsf_348 = buffer.data(nsf + 348);
    const auto *nsf_349 = buffer.data(nsf + 349);
    const auto *nsf_350 = buffer.data(nsf + 350);
    const auto *nsf_351 = buffer.data(nsf + 351);
    const auto *nsf_352 = buffer.data(nsf + 352);
    const auto *nsf_355 = buffer.data(nsf + 355);
    const auto *nsf_356 = buffer.data(nsf + 356);
    const auto *nsf_357 = buffer.data(nsf + 357);
    const auto *nsf_358 = buffer.data(nsf + 358);
    const auto *nsf_359 = buffer.data(nsf + 359);
    const auto *nsf_360 = buffer.data(nsf + 360);
    const auto *nsf_361 = buffer.data(nsf + 361);
    const auto *nsf_362 = buffer.data(nsf + 362);
    const auto *nsf_363 = buffer.data(nsf + 363);
    const auto *nsf_366 = buffer.data(nsf + 366);
    const auto *nsf_367 = buffer.data(nsf + 367);
    const auto *nsf_368 = buffer.data(nsf + 368);
    const auto *nsf_369 = buffer.data(nsf + 369);
    const auto *nsf_370 = buffer.data(nsf + 370);
    const auto *nsf_372 = buffer.data(nsf + 372);
    const auto *nsf_375 = buffer.data(nsf + 375);
    const auto *nsf_376 = buffer.data(nsf + 376);
    const auto *nsf_377 = buffer.data(nsf + 377);
    const auto *nsf_378 = buffer.data(nsf + 378);
    const auto *nsf_379 = buffer.data(nsf + 379);
    const auto *nsf_380 = buffer.data(nsf + 380);
    const auto *nsf_382 = buffer.data(nsf + 382);
    const auto *nsf_383 = buffer.data(nsf + 383);
    const auto *nsf_385 = buffer.data(nsf + 385);
    const auto *nsf_386 = buffer.data(nsf + 386);
    const auto *nsf_387 = buffer.data(nsf + 387);
    const auto *nsf_388 = buffer.data(nsf + 388);
    const auto *nsf_389 = buffer.data(nsf + 389);
    const auto *nsf_390 = buffer.data(nsf + 390);
    const auto *nsf_392 = buffer.data(nsf + 392);
    const auto *nsf_393 = buffer.data(nsf + 393);
    const auto *nsf_395 = buffer.data(nsf + 395);
    const auto *nsf_396 = buffer.data(nsf + 396);
    const auto *nsf_397 = buffer.data(nsf + 397);
    const auto *nsf_398 = buffer.data(nsf + 398);
    const auto *nsf_399 = buffer.data(nsf + 399);
    const auto *nsf_400 = buffer.data(nsf + 400);

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pc_x, msf_325, msf_326, msf_327, msf_328, \
                         nsd0_197, nsd1_197, nsf_325, nsf_326, nsf_327, \
                         nsf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_14 * msf_325[k]
                   + f_4 * nsd0_197[k]
                   - f_5 * nsd1_197[k]
                   + f_3 * pc_x[k] * nsf_325[k];

        t_486[k] = f_14 * msf_326[k]
                   + f_3 * pc_x[k] * nsf_326[k];

        t_487[k] = f_14 * msf_327[k]
                   + f_3 * pc_x[k] * nsf_327[k];

        t_488[k] = f_14 * msf_328[k]
                   + f_3 * pc_x[k] * nsf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_x, pc_y, pc_z, msf_246, msf_256, msf_329, \
                         nsd0_195, nsd1_195, nsf_326, nsf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_14 * msf_329[k]
                   + f_3 * pc_x[k] * nsf_329[k];

        t_490[k] = f_14 * msf_256[k]
                   + f_1 * nsd0_195[k]
                   - f_2 * nsd1_195[k]
                   + f_3 * pc_y[k] * nsf_326[k];

        t_491[k] = f_16 * msf_246[k]
                   + f_3 * pc_z[k] * nsf_326[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, pc_y, pc_z, msf_249, msf_258, msf_259, nsd0_197, \
                         nsd1_197, nsf_328, nsf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_14 * msf_258[k]
                   + f_4 * nsd0_197[k]
                   - f_5 * nsd1_197[k]
                   + f_3 * pc_y[k] * nsf_328[k];

        t_493[k] = f_14 * msf_259[k]
                   + f_3 * pc_y[k] * nsf_329[k];

        t_494[k] = f_16 * msf_249[k]
                   + f_1 * nsd0_197[k]
                   - f_2 * nsd1_197[k]
                   + f_3 * pc_z[k] * nsf_329[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, pc_x, pc_y, pc_z, msf_250, msf_260, msf_330, \
                         nsd0_198, nsd1_198, nsf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_14 * msf_330[k]
                   + f_1 * nsd0_198[k]
                   - f_2 * nsd1_198[k]
                   + f_3 * pc_x[k] * nsf_330[k];

        t_496[k] = f_8 * msf_260[k]
                   + f_3 * pc_y[k] * nsf_330[k];

        t_497[k] = f_17 * msf_250[k]
                   + f_3 * pc_z[k] * nsf_330[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_x, pc_y, msf_262, msf_333, msf_335, nsd0_201, \
                         nsd0_203, nsd1_201, nsd1_203, nsf_332, nsf_333, \
                         nsf_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_14 * msf_333[k]
                   + f_4 * nsd0_201[k]
                   - f_5 * nsd1_201[k]
                   + f_3 * pc_x[k] * nsf_333[k];

        t_499[k] = f_8 * msf_262[k]
                   + f_3 * pc_y[k] * nsf_332[k];

        t_500[k] = f_14 * msf_335[k]
                   + f_4 * nsd0_203[k]
                   - f_5 * nsd1_203[k]
                   + f_3 * pc_x[k] * nsf_335[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, pc_x, msf_336, msf_337, msf_338, msf_339, \
                         nsf_336, nsf_337, nsf_338, nsf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_14 * msf_336[k]
                   + f_3 * pc_x[k] * nsf_336[k];

        t_502[k] = f_14 * msf_337[k]
                   + f_3 * pc_x[k] * nsf_337[k];

        t_503[k] = f_14 * msf_338[k]
                   + f_3 * pc_x[k] * nsf_338[k];

        t_504[k] = f_14 * msf_339[k]
                   + f_3 * pc_x[k] * nsf_339[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pc_y, pc_z, msf_256, msf_266, msf_268, nsd0_201, \
                         nsd0_203, nsd1_201, nsd1_203, nsf_336, \
                         nsf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_8 * msf_266[k]
                   + f_1 * nsd0_201[k]
                   - f_2 * nsd1_201[k]
                   + f_3 * pc_y[k] * nsf_336[k];

        t_506[k] = f_17 * msf_256[k]
                   + f_3 * pc_z[k] * nsf_336[k];

        t_507[k] = f_8 * msf_268[k]
                   + f_4 * nsd0_203[k]
                   - f_5 * nsd1_203[k]
                   + f_3 * pc_y[k] * nsf_338[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_y, pc_y, pc_z, msg0_405, msf_259, \
                         msf_269, msf_270, msg1_405, nsd0_203, nsd1_203, nsf_339, \
                         nsf_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_8 * msf_269[k]
                   + f_3 * pc_y[k] * nsf_339[k];

        t_509[k] = f_17 * msf_259[k]
                   + f_1 * nsd0_203[k]
                   - f_2 * nsd1_203[k]
                   + f_3 * pc_z[k] * nsf_339[k];

        t_510[k] = pa_y[k] * msg0_405[k]
                   - f_6 * pc_y[k] * msg1_405[k];

        t_511[k] = f_7 * msf_270[k]
                   + f_3 * pc_y[k] * nsf_340[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pa_y, pc_y, pc_z, msg0_408, msg0_410, \
                         msf_260, msf_271, msf_272, msg1_408, msg1_410, nsf_340, \
                         nsf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_15 * msf_260[k]
                   + f_3 * pc_z[k] * nsf_340[k];

        t_513[k] = pa_y[k] * msg0_408[k]
                   + f_8 * msf_271[k]
                   - f_6 * pc_y[k] * msg1_408[k];

        t_514[k] = f_7 * msf_272[k]
                   + f_3 * pc_y[k] * nsf_342[k];

        t_515[k] = pa_y[k] * msg0_410[k]
                   - f_6 * pc_y[k] * msg1_410[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, msf_346, msf_347, msf_348, msf_349, \
                         nsf_346, nsf_347, nsf_348, nsf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_14 * msf_346[k]
                   + f_3 * pc_x[k] * nsf_346[k];

        t_517[k] = f_14 * msf_347[k]
                   + f_3 * pc_x[k] * nsf_347[k];

        t_518[k] = f_14 * msf_348[k]
                   + f_3 * pc_x[k] * nsf_348[k];

        t_519[k] = f_14 * msf_349[k]
                   + f_3 * pc_x[k] * nsf_349[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, msf_266, msf_276, msf_278, nsd0_207, \
                         nsd0_209, nsd1_207, nsd1_209, nsf_346, \
                         nsf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_7 * msf_276[k]
                   + f_1 * nsd0_207[k]
                   - f_2 * nsd1_207[k]
                   + f_3 * pc_y[k] * nsf_346[k];

        t_521[k] = f_15 * msf_266[k]
                   + f_3 * pc_z[k] * nsf_346[k];

        t_522[k] = f_7 * msf_278[k]
                   + f_4 * nsd0_209[k]
                   - f_5 * nsd1_209[k]
                   + f_3 * pc_y[k] * nsf_348[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pa_y, pc_x, pc_y, msg0_419, msf_279, \
                         msf_350, msg1_419, nsd0_210, nsd1_210, nsf_349, \
                         nsf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_7 * msf_279[k]
                   + f_3 * pc_y[k] * nsf_349[k];

        t_524[k] = pa_y[k] * msg0_419[k]
                   - f_6 * pc_y[k] * msg1_419[k];

        t_525[k] = f_14 * msf_350[k]
                   + f_1 * nsd0_210[k]
                   - f_2 * nsd1_210[k]
                   + f_3 * pc_x[k] * nsf_350[k];

        t_526[k] = f_3 * pc_y[k] * nsf_350[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pc_y, pc_z, msf_270, nsd0_210, nsd1_210, \
                         nsf_350, nsf_351, nsf_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_13 * msf_270[k]
                   + f_3 * pc_z[k] * nsf_350[k];

        t_528[k] = f_4 * nsd0_210[k]
                   - f_5 * nsd1_210[k]
                   + f_3 * pc_y[k] * nsf_351[k];

        t_529[k] = f_3 * pc_y[k] * nsf_352[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pc_x, pc_y, msf_355, msf_356, msf_357, \
                         nsd0_215, nsd1_215, nsf_355, nsf_356, \
                         nsf_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_14 * msf_355[k]
                   + f_4 * nsd0_215[k]
                   - f_5 * nsd1_215[k]
                   + f_3 * pc_x[k] * nsf_355[k];

        t_531[k] = f_14 * msf_356[k]
                   + f_3 * pc_x[k] * nsf_356[k];

        t_532[k] = f_14 * msf_357[k]
                   + f_3 * pc_x[k] * nsf_357[k];

        t_533[k] = f_3 * pc_y[k] * nsf_355[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_x, pc_y, msf_359, nsd0_213, nsd0_214, \
                         nsd1_213, nsd1_214, nsf_356, nsf_357, \
                         nsf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_14 * msf_359[k]
                   + f_3 * pc_x[k] * nsf_359[k];

        t_535[k] = f_1 * nsd0_213[k]
                   - f_2 * nsd1_213[k]
                   + f_3 * pc_y[k] * nsf_356[k];

        t_536[k] = f_10 * nsd0_214[k]
                   - f_11 * nsd1_214[k]
                   + f_3 * pc_y[k] * nsf_357[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, pc_z, msf_279, msf_360, \
                         nsd0_215, nsd0_216, nsd1_215, nsd1_216, nsf_358, nsf_359, \
                         nsf_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * nsd0_215[k]
                   - f_5 * nsd1_215[k]
                   + f_3 * pc_y[k] * nsf_358[k];

        t_538[k] = f_3 * pc_y[k] * nsf_359[k];

        t_539[k] = f_13 * msf_279[k]
                   + f_1 * nsd0_215[k]
                   - f_2 * nsd1_215[k]
                   + f_3 * pc_z[k] * nsf_359[k];

        t_540[k] = f_8 * msf_360[k]
                   + f_1 * nsd0_216[k]
                   - f_2 * nsd1_216[k]
                   + f_3 * pc_x[k] * nsf_360[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_y, pc_z, msf_280, msf_363, \
                         nsd0_219, nsd1_219, nsf_360, nsf_361, \
                         nsf_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_12 * msf_280[k]
                   + f_3 * pc_y[k] * nsf_360[k];

        t_542[k] = f_3 * pc_z[k] * nsf_360[k];

        t_543[k] = f_8 * msf_363[k]
                   + f_4 * nsd0_219[k]
                   - f_5 * nsd1_219[k]
                   + f_3 * pc_x[k] * nsf_363[k];

        t_544[k] = f_3 * pc_z[k] * nsf_361[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pc_x, pc_z, msf_366, msf_368, nsd0_216, \
                         nsd1_216, nsf_362, nsf_363, nsf_366, nsf_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * nsd0_216[k]
                   - f_5 * nsd1_216[k]
                   + f_3 * pc_z[k] * nsf_362[k];

        t_546[k] = f_8 * msf_366[k]
                   + f_3 * pc_x[k] * nsf_366[k];

        t_547[k] = f_3 * pc_z[k] * nsf_363[k];

        t_548[k] = f_8 * msf_368[k]
                   + f_3 * pc_x[k] * nsf_368[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, t_553, pc_x, pc_y, pc_z, msf_286, \
                         msf_289, msf_369, nsd0_219, nsd1_219, nsf_366, nsf_367, \
                         nsf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_8 * msf_369[k]
                   + f_3 * pc_x[k] * nsf_369[k];

        t_550[k] = f_12 * msf_286[k]
                   + f_1 * nsd0_219[k]
                   - f_2 * nsd1_219[k]
                   + f_3 * pc_y[k] * nsf_366[k];

        t_551[k] = f_3 * pc_z[k] * nsf_366[k];

        t_552[k] = f_4 * nsd0_219[k]
                   - f_5 * nsd1_219[k]
                   + f_3 * pc_z[k] * nsf_367[k];

        t_553[k] = f_12 * msf_289[k]
                   + f_3 * pc_y[k] * nsf_369[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, pa_z, pc_y, pc_z, msg0_420, msf_280, \
                         msf_290, msg1_420, nsd0_221, nsd1_221, nsf_369, \
                         nsf_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_1 * nsd0_221[k]
                   - f_2 * nsd1_221[k]
                   + f_3 * pc_z[k] * nsf_369[k];

        t_555[k] = pa_z[k] * msg0_420[k]
                   - f_6 * pc_z[k] * msg1_420[k];

        t_556[k] = f_13 * msf_290[k]
                   + f_3 * pc_y[k] * nsf_370[k];

        t_557[k] = f_7 * msf_280[k]
                   + f_3 * pc_z[k] * nsf_370[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, pa_z, pc_x, pc_y, pc_z, msg0_423, msf_292, \
                         msf_375, msg1_423, nsd0_227, nsd1_227, nsf_372, \
                         nsf_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = pa_z[k] * msg0_423[k]
                   - f_6 * pc_z[k] * msg1_423[k];

        t_559[k] = f_13 * msf_292[k]
                   + f_3 * pc_y[k] * nsf_372[k];

        t_560[k] = f_8 * msf_375[k]
                   + f_4 * nsd0_227[k]
                   - f_5 * nsd1_227[k]
                   + f_3 * pc_x[k] * nsf_375[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, pc_x, msf_376, msf_377, msf_378, msf_379, \
                         nsf_376, nsf_377, nsf_378, nsf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_8 * msf_376[k]
                   + f_3 * pc_x[k] * nsf_376[k];

        t_562[k] = f_8 * msf_377[k]
                   + f_3 * pc_x[k] * nsf_377[k];

        t_563[k] = f_8 * msf_378[k]
                   + f_3 * pc_x[k] * nsf_378[k];

        t_564[k] = f_8 * msf_379[k]
                   + f_3 * pc_x[k] * nsf_379[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, pa_z, pc_y, pc_z, msg0_430, msf_286, msf_298, \
                         msg1_430, nsd0_227, nsd1_227, nsf_376, \
                         nsf_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = pa_z[k] * msg0_430[k]
                   - f_6 * pc_z[k] * msg1_430[k];

        t_566[k] = f_7 * msf_286[k]
                   + f_3 * pc_z[k] * nsf_376[k];

        t_567[k] = f_13 * msf_298[k]
                   + f_4 * nsd0_227[k]
                   - f_5 * nsd1_227[k]
                   + f_3 * pc_y[k] * nsf_378[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, pc_x, pc_y, pc_z, msf_289, msf_299, msf_380, \
                         nsd0_227, nsd0_228, nsd1_227, nsd1_228, nsf_379, \
                         nsf_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_13 * msf_299[k]
                   + f_3 * pc_y[k] * nsf_379[k];

        t_569[k] = f_7 * msf_289[k]
                   + f_1 * nsd0_227[k]
                   - f_2 * nsd1_227[k]
                   + f_3 * pc_z[k] * nsf_379[k];

        t_570[k] = f_8 * msf_380[k]
                   + f_1 * nsd0_228[k]
                   - f_2 * nsd1_228[k]
                   + f_3 * pc_x[k] * nsf_380[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, t_574, pc_x, pc_y, pc_z, msf_290, msf_300, \
                         msf_302, msf_383, nsd0_231, nsd1_231, nsf_380, nsf_382, \
                         nsf_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_15 * msf_300[k]
                   + f_3 * pc_y[k] * nsf_380[k];

        t_572[k] = f_8 * msf_290[k]
                   + f_3 * pc_z[k] * nsf_380[k];

        t_573[k] = f_8 * msf_383[k]
                   + f_4 * nsd0_231[k]
                   - f_5 * nsd1_231[k]
                   + f_3 * pc_x[k] * nsf_383[k];

        t_574[k] = f_15 * msf_302[k]
                   + f_3 * pc_y[k] * nsf_382[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, pc_x, msf_385, msf_386, msf_387, msf_388, \
                         nsd0_233, nsd1_233, nsf_385, nsf_386, nsf_387, \
                         nsf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_8 * msf_385[k]
                   + f_4 * nsd0_233[k]
                   - f_5 * nsd1_233[k]
                   + f_3 * pc_x[k] * nsf_385[k];

        t_576[k] = f_8 * msf_386[k]
                   + f_3 * pc_x[k] * nsf_386[k];

        t_577[k] = f_8 * msf_387[k]
                   + f_3 * pc_x[k] * nsf_387[k];

        t_578[k] = f_8 * msf_388[k]
                   + f_3 * pc_x[k] * nsf_388[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pc_x, pc_y, pc_z, msf_296, msf_306, msf_389, \
                         nsd0_231, nsd1_231, nsf_386, nsf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_8 * msf_389[k]
                   + f_3 * pc_x[k] * nsf_389[k];

        t_580[k] = f_15 * msf_306[k]
                   + f_1 * nsd0_231[k]
                   - f_2 * nsd1_231[k]
                   + f_3 * pc_y[k] * nsf_386[k];

        t_581[k] = f_8 * msf_296[k]
                   + f_3 * pc_z[k] * nsf_386[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pc_y, pc_z, msf_299, msf_308, msf_309, nsd0_233, \
                         nsd1_233, nsf_388, nsf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_15 * msf_308[k]
                   + f_4 * nsd0_233[k]
                   - f_5 * nsd1_233[k]
                   + f_3 * pc_y[k] * nsf_388[k];

        t_583[k] = f_15 * msf_309[k]
                   + f_3 * pc_y[k] * nsf_389[k];

        t_584[k] = f_8 * msf_299[k]
                   + f_1 * nsd0_233[k]
                   - f_2 * nsd1_233[k]
                   + f_3 * pc_z[k] * nsf_389[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, pc_x, pc_y, pc_z, msf_300, msf_310, msf_390, \
                         nsd0_234, nsd1_234, nsf_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_8 * msf_390[k]
                   + f_1 * nsd0_234[k]
                   - f_2 * nsd1_234[k]
                   + f_3 * pc_x[k] * nsf_390[k];

        t_586[k] = f_17 * msf_310[k]
                   + f_3 * pc_y[k] * nsf_390[k];

        t_587[k] = f_14 * msf_300[k]
                   + f_3 * pc_z[k] * nsf_390[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, pc_x, pc_y, msf_312, msf_393, msf_395, nsd0_237, \
                         nsd0_239, nsd1_237, nsd1_239, nsf_392, nsf_393, \
                         nsf_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_8 * msf_393[k]
                   + f_4 * nsd0_237[k]
                   - f_5 * nsd1_237[k]
                   + f_3 * pc_x[k] * nsf_393[k];

        t_589[k] = f_17 * msf_312[k]
                   + f_3 * pc_y[k] * nsf_392[k];

        t_590[k] = f_8 * msf_395[k]
                   + f_4 * nsd0_239[k]
                   - f_5 * nsd1_239[k]
                   + f_3 * pc_x[k] * nsf_395[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, pc_x, msf_396, msf_397, msf_398, msf_399, \
                         nsf_396, nsf_397, nsf_398, nsf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_8 * msf_396[k]
                   + f_3 * pc_x[k] * nsf_396[k];

        t_592[k] = f_8 * msf_397[k]
                   + f_3 * pc_x[k] * nsf_397[k];

        t_593[k] = f_8 * msf_398[k]
                   + f_3 * pc_x[k] * nsf_398[k];

        t_594[k] = f_8 * msf_399[k]
                   + f_3 * pc_x[k] * nsf_399[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pc_y, pc_z, msf_306, msf_316, msf_318, nsd0_237, \
                         nsd0_239, nsd1_237, nsd1_239, nsf_396, \
                         nsf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_17 * msf_316[k]
                   + f_1 * nsd0_237[k]
                   - f_2 * nsd1_237[k]
                   + f_3 * pc_y[k] * nsf_396[k];

        t_596[k] = f_14 * msf_306[k]
                   + f_3 * pc_z[k] * nsf_396[k];

        t_597[k] = f_17 * msf_318[k]
                   + f_4 * nsd0_239[k]
                   - f_5 * nsd1_239[k]
                   + f_3 * pc_y[k] * nsf_398[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pc_x, pc_y, pc_z, msf_309, msf_319, msf_400, \
                         nsd0_239, nsd0_240, nsd1_239, nsd1_240, nsf_399, \
                         nsf_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_17 * msf_319[k]
                   + f_3 * pc_y[k] * nsf_399[k];

        t_599[k] = f_14 * msf_309[k]
                   + f_1 * nsd0_239[k]
                   - f_2 * nsd1_239[k]
                   + f_3 * pc_z[k] * nsf_399[k];

        t_600[k] = f_8 * msf_400[k]
                   + f_1 * nsd0_240[k]
                   - f_2 * nsd1_240[k]
                   + f_3 * pc_x[k] * nsf_400[k];
    }
}

static auto
compute_prim_nsg_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msg0,
                                                          const size_t msf, const size_t msg1,
                                                          const size_t nsd0, const size_t nsd1,
                                                          const size_t nsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msg0_525 = buffer.data(msg0 + 525);
    const auto *msg0_528 = buffer.data(msg0 + 528);
    const auto *msg0_530 = buffer.data(msg0 + 530);
    const auto *msg0_539 = buffer.data(msg0 + 539);
    const auto *msg0_540 = buffer.data(msg0 + 540);
    const auto *msg0_543 = buffer.data(msg0 + 543);
    const auto *msg0_675 = buffer.data(msg0 + 675);
    const auto *msg0_678 = buffer.data(msg0 + 678);
    const auto *msg0_685 = buffer.data(msg0 + 685);
    const auto *msg0_687 = buffer.data(msg0 + 687);
    const auto *msg0_689 = buffer.data(msg0 + 689);
    const auto *msg0_695 = buffer.data(msg0 + 695);
    const auto *msg0_700 = buffer.data(msg0 + 700);
    const auto *msg0_702 = buffer.data(msg0 + 702);
    const auto *msg0_704 = buffer.data(msg0 + 704);
    const auto *msg0_705 = buffer.data(msg0 + 705);
    const auto *msg0_708 = buffer.data(msg0 + 708);
    const auto *msg0_710 = buffer.data(msg0 + 710);
    const auto *msg0_715 = buffer.data(msg0 + 715);
    const auto *msg0_717 = buffer.data(msg0 + 717);
    const auto *msg0_719 = buffer.data(msg0 + 719);
    const auto *msg0_720 = buffer.data(msg0 + 720);

    const auto *msf_310 = buffer.data(msf + 310);
    const auto *msf_316 = buffer.data(msf + 316);
    const auto *msf_319 = buffer.data(msf + 319);
    const auto *msf_320 = buffer.data(msf + 320);
    const auto *msf_322 = buffer.data(msf + 322);
    const auto *msf_326 = buffer.data(msf + 326);
    const auto *msf_328 = buffer.data(msf + 328);
    const auto *msf_329 = buffer.data(msf + 329);
    const auto *msf_330 = buffer.data(msf + 330);
    const auto *msf_332 = buffer.data(msf + 332);
    const auto *msf_336 = buffer.data(msf + 336);
    const auto *msf_338 = buffer.data(msf + 338);
    const auto *msf_339 = buffer.data(msf + 339);
    const auto *msf_340 = buffer.data(msf + 340);
    const auto *msf_342 = buffer.data(msf + 342);
    const auto *msf_346 = buffer.data(msf + 346);
    const auto *msf_348 = buffer.data(msf + 348);
    const auto *msf_349 = buffer.data(msf + 349);
    const auto *msf_350 = buffer.data(msf + 350);
    const auto *msf_351 = buffer.data(msf + 351);
    const auto *msf_352 = buffer.data(msf + 352);
    const auto *msf_356 = buffer.data(msf + 356);
    const auto *msf_358 = buffer.data(msf + 358);
    const auto *msf_359 = buffer.data(msf + 359);
    const auto *msf_360 = buffer.data(msf + 360);
    const auto *msf_366 = buffer.data(msf + 366);
    const auto *msf_369 = buffer.data(msf + 369);
    const auto *msf_370 = buffer.data(msf + 370);
    const auto *msf_372 = buffer.data(msf + 372);
    const auto *msf_376 = buffer.data(msf + 376);
    const auto *msf_379 = buffer.data(msf + 379);
    const auto *msf_380 = buffer.data(msf + 380);
    const auto *msf_382 = buffer.data(msf + 382);
    const auto *msf_389 = buffer.data(msf + 389);
    const auto *msf_390 = buffer.data(msf + 390);
    const auto *msf_403 = buffer.data(msf + 403);
    const auto *msf_405 = buffer.data(msf + 405);
    const auto *msf_406 = buffer.data(msf + 406);
    const auto *msf_407 = buffer.data(msf + 407);
    const auto *msf_408 = buffer.data(msf + 408);
    const auto *msf_409 = buffer.data(msf + 409);
    const auto *msf_410 = buffer.data(msf + 410);
    const auto *msf_413 = buffer.data(msf + 413);
    const auto *msf_415 = buffer.data(msf + 415);
    const auto *msf_416 = buffer.data(msf + 416);
    const auto *msf_417 = buffer.data(msf + 417);
    const auto *msf_418 = buffer.data(msf + 418);
    const auto *msf_419 = buffer.data(msf + 419);
    const auto *msf_420 = buffer.data(msf + 420);
    const auto *msf_423 = buffer.data(msf + 423);
    const auto *msf_425 = buffer.data(msf + 425);
    const auto *msf_426 = buffer.data(msf + 426);
    const auto *msf_427 = buffer.data(msf + 427);
    const auto *msf_428 = buffer.data(msf + 428);
    const auto *msf_429 = buffer.data(msf + 429);
    const auto *msf_436 = buffer.data(msf + 436);
    const auto *msf_437 = buffer.data(msf + 437);
    const auto *msf_438 = buffer.data(msf + 438);
    const auto *msf_439 = buffer.data(msf + 439);
    const auto *msf_440 = buffer.data(msf + 440);
    const auto *msf_445 = buffer.data(msf + 445);
    const auto *msf_446 = buffer.data(msf + 446);
    const auto *msf_447 = buffer.data(msf + 447);
    const auto *msf_449 = buffer.data(msf + 449);
    const auto *msf_450 = buffer.data(msf + 450);
    const auto *msf_453 = buffer.data(msf + 453);
    const auto *msf_456 = buffer.data(msf + 456);
    const auto *msf_458 = buffer.data(msf + 458);
    const auto *msf_459 = buffer.data(msf + 459);
    const auto *msf_465 = buffer.data(msf + 465);
    const auto *msf_466 = buffer.data(msf + 466);
    const auto *msf_467 = buffer.data(msf + 467);
    const auto *msf_468 = buffer.data(msf + 468);
    const auto *msf_469 = buffer.data(msf + 469);
    const auto *msf_470 = buffer.data(msf + 470);
    const auto *msf_473 = buffer.data(msf + 473);
    const auto *msf_475 = buffer.data(msf + 475);
    const auto *msf_476 = buffer.data(msf + 476);
    const auto *msf_477 = buffer.data(msf + 477);
    const auto *msf_478 = buffer.data(msf + 478);
    const auto *msf_479 = buffer.data(msf + 479);
    const auto *msf_480 = buffer.data(msf + 480);

    const auto *msg1_525 = buffer.data(msg1 + 525);
    const auto *msg1_528 = buffer.data(msg1 + 528);
    const auto *msg1_530 = buffer.data(msg1 + 530);
    const auto *msg1_539 = buffer.data(msg1 + 539);
    const auto *msg1_540 = buffer.data(msg1 + 540);
    const auto *msg1_543 = buffer.data(msg1 + 543);
    const auto *msg1_675 = buffer.data(msg1 + 675);
    const auto *msg1_678 = buffer.data(msg1 + 678);
    const auto *msg1_685 = buffer.data(msg1 + 685);
    const auto *msg1_687 = buffer.data(msg1 + 687);
    const auto *msg1_689 = buffer.data(msg1 + 689);
    const auto *msg1_695 = buffer.data(msg1 + 695);
    const auto *msg1_700 = buffer.data(msg1 + 700);
    const auto *msg1_702 = buffer.data(msg1 + 702);
    const auto *msg1_704 = buffer.data(msg1 + 704);
    const auto *msg1_705 = buffer.data(msg1 + 705);
    const auto *msg1_708 = buffer.data(msg1 + 708);
    const auto *msg1_710 = buffer.data(msg1 + 710);
    const auto *msg1_715 = buffer.data(msg1 + 715);
    const auto *msg1_717 = buffer.data(msg1 + 717);
    const auto *msg1_719 = buffer.data(msg1 + 719);
    const auto *msg1_720 = buffer.data(msg1 + 720);

    const auto *nsd0_243 = buffer.data(nsd0 + 243);
    const auto *nsd0_245 = buffer.data(nsd0 + 245);
    const auto *nsd0_246 = buffer.data(nsd0 + 246);
    const auto *nsd0_249 = buffer.data(nsd0 + 249);
    const auto *nsd0_251 = buffer.data(nsd0 + 251);
    const auto *nsd0_252 = buffer.data(nsd0 + 252);
    const auto *nsd0_255 = buffer.data(nsd0 + 255);
    const auto *nsd0_257 = buffer.data(nsd0 + 257);
    const auto *nsd0_261 = buffer.data(nsd0 + 261);
    const auto *nsd0_263 = buffer.data(nsd0 + 263);
    const auto *nsd0_264 = buffer.data(nsd0 + 264);
    const auto *nsd0_267 = buffer.data(nsd0 + 267);
    const auto *nsd0_268 = buffer.data(nsd0 + 268);
    const auto *nsd0_269 = buffer.data(nsd0 + 269);
    const auto *nsd0_270 = buffer.data(nsd0 + 270);

    const auto *nsd1_243 = buffer.data(nsd1 + 243);
    const auto *nsd1_245 = buffer.data(nsd1 + 245);
    const auto *nsd1_246 = buffer.data(nsd1 + 246);
    const auto *nsd1_249 = buffer.data(nsd1 + 249);
    const auto *nsd1_251 = buffer.data(nsd1 + 251);
    const auto *nsd1_252 = buffer.data(nsd1 + 252);
    const auto *nsd1_255 = buffer.data(nsd1 + 255);
    const auto *nsd1_257 = buffer.data(nsd1 + 257);
    const auto *nsd1_261 = buffer.data(nsd1 + 261);
    const auto *nsd1_263 = buffer.data(nsd1 + 263);
    const auto *nsd1_264 = buffer.data(nsd1 + 264);
    const auto *nsd1_267 = buffer.data(nsd1 + 267);
    const auto *nsd1_268 = buffer.data(nsd1 + 268);
    const auto *nsd1_269 = buffer.data(nsd1 + 269);
    const auto *nsd1_270 = buffer.data(nsd1 + 270);

    const auto *nsf_400 = buffer.data(nsf + 400);
    const auto *nsf_402 = buffer.data(nsf + 402);
    const auto *nsf_403 = buffer.data(nsf + 403);
    const auto *nsf_405 = buffer.data(nsf + 405);
    const auto *nsf_406 = buffer.data(nsf + 406);
    const auto *nsf_407 = buffer.data(nsf + 407);
    const auto *nsf_408 = buffer.data(nsf + 408);
    const auto *nsf_409 = buffer.data(nsf + 409);
    const auto *nsf_410 = buffer.data(nsf + 410);
    const auto *nsf_412 = buffer.data(nsf + 412);
    const auto *nsf_413 = buffer.data(nsf + 413);
    const auto *nsf_415 = buffer.data(nsf + 415);
    const auto *nsf_416 = buffer.data(nsf + 416);
    const auto *nsf_417 = buffer.data(nsf + 417);
    const auto *nsf_418 = buffer.data(nsf + 418);
    const auto *nsf_419 = buffer.data(nsf + 419);
    const auto *nsf_420 = buffer.data(nsf + 420);
    const auto *nsf_422 = buffer.data(nsf + 422);
    const auto *nsf_423 = buffer.data(nsf + 423);
    const auto *nsf_425 = buffer.data(nsf + 425);
    const auto *nsf_426 = buffer.data(nsf + 426);
    const auto *nsf_427 = buffer.data(nsf + 427);
    const auto *nsf_428 = buffer.data(nsf + 428);
    const auto *nsf_429 = buffer.data(nsf + 429);
    const auto *nsf_430 = buffer.data(nsf + 430);
    const auto *nsf_432 = buffer.data(nsf + 432);
    const auto *nsf_436 = buffer.data(nsf + 436);
    const auto *nsf_437 = buffer.data(nsf + 437);
    const auto *nsf_438 = buffer.data(nsf + 438);
    const auto *nsf_439 = buffer.data(nsf + 439);
    const auto *nsf_440 = buffer.data(nsf + 440);
    const auto *nsf_441 = buffer.data(nsf + 441);
    const auto *nsf_442 = buffer.data(nsf + 442);
    const auto *nsf_445 = buffer.data(nsf + 445);
    const auto *nsf_446 = buffer.data(nsf + 446);
    const auto *nsf_447 = buffer.data(nsf + 447);
    const auto *nsf_448 = buffer.data(nsf + 448);
    const auto *nsf_449 = buffer.data(nsf + 449);
    const auto *nsf_450 = buffer.data(nsf + 450);
    const auto *nsf_451 = buffer.data(nsf + 451);
    const auto *nsf_452 = buffer.data(nsf + 452);
    const auto *nsf_453 = buffer.data(nsf + 453);
    const auto *nsf_456 = buffer.data(nsf + 456);
    const auto *nsf_458 = buffer.data(nsf + 458);
    const auto *nsf_459 = buffer.data(nsf + 459);
    const auto *nsf_460 = buffer.data(nsf + 460);
    const auto *nsf_462 = buffer.data(nsf + 462);
    const auto *nsf_466 = buffer.data(nsf + 466);
    const auto *nsf_467 = buffer.data(nsf + 467);
    const auto *nsf_468 = buffer.data(nsf + 468);
    const auto *nsf_469 = buffer.data(nsf + 469);
    const auto *nsf_470 = buffer.data(nsf + 470);
    const auto *nsf_472 = buffer.data(nsf + 472);
    const auto *nsf_476 = buffer.data(nsf + 476);
    const auto *nsf_477 = buffer.data(nsf + 477);
    const auto *nsf_478 = buffer.data(nsf + 478);
    const auto *nsf_479 = buffer.data(nsf + 479);
    const auto *nsf_480 = buffer.data(nsf + 480);

#pragma omp simd aligned(t_601, t_602, t_603, t_604, pc_x, pc_y, pc_z, msf_310, msf_320, \
                         msf_322, msf_403, nsd0_243, nsd1_243, nsf_400, nsf_402, \
                         nsf_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_16 * msf_320[k]
                   + f_3 * pc_y[k] * nsf_400[k];

        t_602[k] = f_16 * msf_310[k]
                   + f_3 * pc_z[k] * nsf_400[k];

        t_603[k] = f_8 * msf_403[k]
                   + f_4 * nsd0_243[k]
                   - f_5 * nsd1_243[k]
                   + f_3 * pc_x[k] * nsf_403[k];

        t_604[k] = f_16 * msf_322[k]
                   + f_3 * pc_y[k] * nsf_402[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pc_x, msf_405, msf_406, msf_407, msf_408, \
                         nsd0_245, nsd1_245, nsf_405, nsf_406, nsf_407, \
                         nsf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_8 * msf_405[k]
                   + f_4 * nsd0_245[k]
                   - f_5 * nsd1_245[k]
                   + f_3 * pc_x[k] * nsf_405[k];

        t_606[k] = f_8 * msf_406[k]
                   + f_3 * pc_x[k] * nsf_406[k];

        t_607[k] = f_8 * msf_407[k]
                   + f_3 * pc_x[k] * nsf_407[k];

        t_608[k] = f_8 * msf_408[k]
                   + f_3 * pc_x[k] * nsf_408[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pc_x, pc_y, pc_z, msf_316, msf_326, msf_409, \
                         nsd0_243, nsd1_243, nsf_406, nsf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_8 * msf_409[k]
                   + f_3 * pc_x[k] * nsf_409[k];

        t_610[k] = f_16 * msf_326[k]
                   + f_1 * nsd0_243[k]
                   - f_2 * nsd1_243[k]
                   + f_3 * pc_y[k] * nsf_406[k];

        t_611[k] = f_16 * msf_316[k]
                   + f_3 * pc_z[k] * nsf_406[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pc_y, pc_z, msf_319, msf_328, msf_329, nsd0_245, \
                         nsd1_245, nsf_408, nsf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_16 * msf_328[k]
                   + f_4 * nsd0_245[k]
                   - f_5 * nsd1_245[k]
                   + f_3 * pc_y[k] * nsf_408[k];

        t_613[k] = f_16 * msf_329[k]
                   + f_3 * pc_y[k] * nsf_409[k];

        t_614[k] = f_16 * msf_319[k]
                   + f_1 * nsd0_245[k]
                   - f_2 * nsd1_245[k]
                   + f_3 * pc_z[k] * nsf_409[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pc_x, pc_y, pc_z, msf_320, msf_330, msf_410, \
                         nsd0_246, nsd1_246, nsf_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_8 * msf_410[k]
                   + f_1 * nsd0_246[k]
                   - f_2 * nsd1_246[k]
                   + f_3 * pc_x[k] * nsf_410[k];

        t_616[k] = f_14 * msf_330[k]
                   + f_3 * pc_y[k] * nsf_410[k];

        t_617[k] = f_17 * msf_320[k]
                   + f_3 * pc_z[k] * nsf_410[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, pc_x, pc_y, msf_332, msf_413, msf_415, nsd0_249, \
                         nsd0_251, nsd1_249, nsd1_251, nsf_412, nsf_413, \
                         nsf_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_8 * msf_413[k]
                   + f_4 * nsd0_249[k]
                   - f_5 * nsd1_249[k]
                   + f_3 * pc_x[k] * nsf_413[k];

        t_619[k] = f_14 * msf_332[k]
                   + f_3 * pc_y[k] * nsf_412[k];

        t_620[k] = f_8 * msf_415[k]
                   + f_4 * nsd0_251[k]
                   - f_5 * nsd1_251[k]
                   + f_3 * pc_x[k] * nsf_415[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, pc_x, msf_416, msf_417, msf_418, msf_419, \
                         nsf_416, nsf_417, nsf_418, nsf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_8 * msf_416[k]
                   + f_3 * pc_x[k] * nsf_416[k];

        t_622[k] = f_8 * msf_417[k]
                   + f_3 * pc_x[k] * nsf_417[k];

        t_623[k] = f_8 * msf_418[k]
                   + f_3 * pc_x[k] * nsf_418[k];

        t_624[k] = f_8 * msf_419[k]
                   + f_3 * pc_x[k] * nsf_419[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, pc_y, pc_z, msf_326, msf_336, msf_338, nsd0_249, \
                         nsd0_251, nsd1_249, nsd1_251, nsf_416, \
                         nsf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_14 * msf_336[k]
                   + f_1 * nsd0_249[k]
                   - f_2 * nsd1_249[k]
                   + f_3 * pc_y[k] * nsf_416[k];

        t_626[k] = f_17 * msf_326[k]
                   + f_3 * pc_z[k] * nsf_416[k];

        t_627[k] = f_14 * msf_338[k]
                   + f_4 * nsd0_251[k]
                   - f_5 * nsd1_251[k]
                   + f_3 * pc_y[k] * nsf_418[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, pc_x, pc_y, pc_z, msf_329, msf_339, msf_420, \
                         nsd0_251, nsd0_252, nsd1_251, nsd1_252, nsf_419, \
                         nsf_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_14 * msf_339[k]
                   + f_3 * pc_y[k] * nsf_419[k];

        t_629[k] = f_17 * msf_329[k]
                   + f_1 * nsd0_251[k]
                   - f_2 * nsd1_251[k]
                   + f_3 * pc_z[k] * nsf_419[k];

        t_630[k] = f_8 * msf_420[k]
                   + f_1 * nsd0_252[k]
                   - f_2 * nsd1_252[k]
                   + f_3 * pc_x[k] * nsf_420[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, pc_x, pc_y, pc_z, msf_330, msf_340, \
                         msf_342, msf_423, nsd0_255, nsd1_255, nsf_420, nsf_422, \
                         nsf_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_8 * msf_340[k]
                   + f_3 * pc_y[k] * nsf_420[k];

        t_632[k] = f_15 * msf_330[k]
                   + f_3 * pc_z[k] * nsf_420[k];

        t_633[k] = f_8 * msf_423[k]
                   + f_4 * nsd0_255[k]
                   - f_5 * nsd1_255[k]
                   + f_3 * pc_x[k] * nsf_423[k];

        t_634[k] = f_8 * msf_342[k]
                   + f_3 * pc_y[k] * nsf_422[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, pc_x, msf_425, msf_426, msf_427, msf_428, \
                         nsd0_257, nsd1_257, nsf_425, nsf_426, nsf_427, \
                         nsf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_8 * msf_425[k]
                   + f_4 * nsd0_257[k]
                   - f_5 * nsd1_257[k]
                   + f_3 * pc_x[k] * nsf_425[k];

        t_636[k] = f_8 * msf_426[k]
                   + f_3 * pc_x[k] * nsf_426[k];

        t_637[k] = f_8 * msf_427[k]
                   + f_3 * pc_x[k] * nsf_427[k];

        t_638[k] = f_8 * msf_428[k]
                   + f_3 * pc_x[k] * nsf_428[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_y, pc_z, msf_336, msf_346, msf_429, \
                         nsd0_255, nsd1_255, nsf_426, nsf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_8 * msf_429[k]
                   + f_3 * pc_x[k] * nsf_429[k];

        t_640[k] = f_8 * msf_346[k]
                   + f_1 * nsd0_255[k]
                   - f_2 * nsd1_255[k]
                   + f_3 * pc_y[k] * nsf_426[k];

        t_641[k] = f_15 * msf_336[k]
                   + f_3 * pc_z[k] * nsf_426[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pa_y, pc_y, pc_z, msg0_525, msf_339, \
                         msf_348, msf_349, msg1_525, nsd0_257, nsd1_257, nsf_428, \
                         nsf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_8 * msf_348[k]
                   + f_4 * nsd0_257[k]
                   - f_5 * nsd1_257[k]
                   + f_3 * pc_y[k] * nsf_428[k];

        t_643[k] = f_8 * msf_349[k]
                   + f_3 * pc_y[k] * nsf_429[k];

        t_644[k] = f_15 * msf_339[k]
                   + f_1 * nsd0_257[k]
                   - f_2 * nsd1_257[k]
                   + f_3 * pc_z[k] * nsf_429[k];

        t_645[k] = pa_y[k] * msg0_525[k]
                   - f_6 * pc_y[k] * msg1_525[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pa_y, pc_y, pc_z, msg0_528, msf_340, \
                         msf_350, msf_351, msf_352, msg1_528, nsf_430, \
                         nsf_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_7 * msf_350[k]
                   + f_3 * pc_y[k] * nsf_430[k];

        t_647[k] = f_13 * msf_340[k]
                   + f_3 * pc_z[k] * nsf_430[k];

        t_648[k] = pa_y[k] * msg0_528[k]
                   + f_8 * msf_351[k]
                   - f_6 * pc_y[k] * msg1_528[k];

        t_649[k] = f_7 * msf_352[k]
                   + f_3 * pc_y[k] * nsf_432[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, pa_y, pc_x, pc_y, msg0_530, msf_436, \
                         msf_437, msf_438, msg1_530, nsf_436, nsf_437, \
                         nsf_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = pa_y[k] * msg0_530[k]
                   - f_6 * pc_y[k] * msg1_530[k];

        t_651[k] = f_8 * msf_436[k]
                   + f_3 * pc_x[k] * nsf_436[k];

        t_652[k] = f_8 * msf_437[k]
                   + f_3 * pc_x[k] * nsf_437[k];

        t_653[k] = f_8 * msf_438[k]
                   + f_3 * pc_x[k] * nsf_438[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_x, pc_y, pc_z, msf_346, msf_356, msf_439, \
                         nsd0_261, nsd1_261, nsf_436, nsf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_8 * msf_439[k]
                   + f_3 * pc_x[k] * nsf_439[k];

        t_655[k] = f_7 * msf_356[k]
                   + f_1 * nsd0_261[k]
                   - f_2 * nsd1_261[k]
                   + f_3 * pc_y[k] * nsf_436[k];

        t_656[k] = f_13 * msf_346[k]
                   + f_3 * pc_z[k] * nsf_436[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pa_y, pc_y, msg0_539, msf_358, msf_359, \
                         msg1_539, nsd0_263, nsd1_263, nsf_438, \
                         nsf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_7 * msf_358[k]
                   + f_4 * nsd0_263[k]
                   - f_5 * nsd1_263[k]
                   + f_3 * pc_y[k] * nsf_438[k];

        t_658[k] = f_7 * msf_359[k]
                   + f_3 * pc_y[k] * nsf_439[k];

        t_659[k] = pa_y[k] * msg0_539[k]
                   - f_6 * pc_y[k] * msg1_539[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, pc_x, pc_y, pc_z, msf_350, \
                         msf_440, nsd0_264, nsd1_264, nsf_440, nsf_441, \
                         nsf_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_8 * msf_440[k]
                   + f_1 * nsd0_264[k]
                   - f_2 * nsd1_264[k]
                   + f_3 * pc_x[k] * nsf_440[k];

        t_661[k] = f_3 * pc_y[k] * nsf_440[k];

        t_662[k] = f_12 * msf_350[k]
                   + f_3 * pc_z[k] * nsf_440[k];

        t_663[k] = f_4 * nsd0_264[k]
                   - f_5 * nsd1_264[k]
                   + f_3 * pc_y[k] * nsf_441[k];

        t_664[k] = f_3 * pc_y[k] * nsf_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pc_x, pc_y, msf_445, msf_446, msf_447, \
                         nsd0_269, nsd1_269, nsf_445, nsf_446, \
                         nsf_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_8 * msf_445[k]
                   + f_4 * nsd0_269[k]
                   - f_5 * nsd1_269[k]
                   + f_3 * pc_x[k] * nsf_445[k];

        t_666[k] = f_8 * msf_446[k]
                   + f_3 * pc_x[k] * nsf_446[k];

        t_667[k] = f_8 * msf_447[k]
                   + f_3 * pc_x[k] * nsf_447[k];

        t_668[k] = f_3 * pc_y[k] * nsf_445[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pc_x, pc_y, msf_449, nsd0_267, nsd0_268, \
                         nsd1_267, nsd1_268, nsf_446, nsf_447, \
                         nsf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_8 * msf_449[k]
                   + f_3 * pc_x[k] * nsf_449[k];

        t_670[k] = f_1 * nsd0_267[k]
                   - f_2 * nsd1_267[k]
                   + f_3 * pc_y[k] * nsf_446[k];

        t_671[k] = f_10 * nsd0_268[k]
                   - f_11 * nsd1_268[k]
                   + f_3 * pc_y[k] * nsf_447[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pa_x, pc_x, pc_y, pc_z, msg0_675, \
                         msf_359, msf_450, msg1_675, nsd0_269, nsd1_269, nsf_448, \
                         nsf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_4 * nsd0_269[k]
                   - f_5 * nsd1_269[k]
                   + f_3 * pc_y[k] * nsf_448[k];

        t_673[k] = f_3 * pc_y[k] * nsf_449[k];

        t_674[k] = f_12 * msf_359[k]
                   + f_1 * nsd0_269[k]
                   - f_2 * nsd1_269[k]
                   + f_3 * pc_z[k] * nsf_449[k];

        t_675[k] = pa_x[k] * msg0_675[k]
                   + f_16 * msf_450[k]
                   - f_6 * pc_x[k] * msg1_675[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pa_x, pc_x, pc_y, pc_z, msg0_678, \
                         msf_360, msf_453, msg1_678, nsf_450, nsf_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_9 * msf_360[k]
                   + f_3 * pc_y[k] * nsf_450[k];

        t_677[k] = f_3 * pc_z[k] * nsf_450[k];

        t_678[k] = pa_x[k] * msg0_678[k]
                   + f_8 * msf_453[k]
                   - f_6 * pc_x[k] * msg1_678[k];

        t_679[k] = f_3 * pc_z[k] * nsf_451[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_z, msf_456, msf_458, nsd0_270, \
                         nsd1_270, nsf_452, nsf_453, nsf_456, nsf_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_4 * nsd0_270[k]
                   - f_5 * nsd1_270[k]
                   + f_3 * pc_z[k] * nsf_452[k];

        t_681[k] = f_7 * msf_456[k]
                   + f_3 * pc_x[k] * nsf_456[k];

        t_682[k] = f_3 * pc_z[k] * nsf_453[k];

        t_683[k] = f_7 * msf_458[k]
                   + f_3 * pc_x[k] * nsf_458[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pa_x, pc_x, pc_z, msg0_685, msg0_687, \
                         msf_459, msg1_685, msg1_687, nsf_456, \
                         nsf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_7 * msf_459[k]
                   + f_3 * pc_x[k] * nsf_459[k];

        t_685[k] = pa_x[k] * msg0_685[k]
                   - f_6 * pc_x[k] * msg1_685[k];

        t_686[k] = f_3 * pc_z[k] * nsf_456[k];

        t_687[k] = pa_x[k] * msg0_687[k]
                   - f_6 * pc_x[k] * msg1_687[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, pa_x, pa_z, pc_x, pc_y, pc_z, msg0_540, \
                         msg0_689, msf_369, msg1_540, msg1_689, \
                         nsf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_9 * msf_369[k]
                   + f_3 * pc_y[k] * nsf_459[k];

        t_689[k] = pa_x[k] * msg0_689[k]
                   - f_6 * pc_x[k] * msg1_689[k];

        t_690[k] = pa_z[k] * msg0_540[k]
                   - f_6 * pc_z[k] * msg1_540[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, pa_z, pc_y, pc_z, msg0_543, msf_360, \
                         msf_370, msf_372, msg1_543, nsf_460, nsf_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_12 * msf_370[k]
                   + f_3 * pc_y[k] * nsf_460[k];

        t_692[k] = f_7 * msf_360[k]
                   + f_3 * pc_z[k] * nsf_460[k];

        t_693[k] = pa_z[k] * msg0_543[k]
                   - f_6 * pc_z[k] * msg1_543[k];

        t_694[k] = f_12 * msf_372[k]
                   + f_3 * pc_y[k] * nsf_462[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, pa_x, pc_x, msg0_695, msf_465, msf_466, \
                         msf_467, msf_468, msg1_695, nsf_466, nsf_467, \
                         nsf_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = pa_x[k] * msg0_695[k]
                   + f_8 * msf_465[k]
                   - f_6 * pc_x[k] * msg1_695[k];

        t_696[k] = f_7 * msf_466[k]
                   + f_3 * pc_x[k] * nsf_466[k];

        t_697[k] = f_7 * msf_467[k]
                   + f_3 * pc_x[k] * nsf_467[k];

        t_698[k] = f_7 * msf_468[k]
                   + f_3 * pc_x[k] * nsf_468[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pa_x, pc_x, pc_z, msg0_700, msg0_702, \
                         msf_366, msf_469, msg1_700, msg1_702, nsf_466, \
                         nsf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_7 * msf_469[k]
                   + f_3 * pc_x[k] * nsf_469[k];

        t_700[k] = pa_x[k] * msg0_700[k]
                   - f_6 * pc_x[k] * msg1_700[k];

        t_701[k] = f_7 * msf_366[k]
                   + f_3 * pc_z[k] * nsf_466[k];

        t_702[k] = pa_x[k] * msg0_702[k]
                   - f_6 * pc_x[k] * msg1_702[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, t_706, pa_x, pc_x, pc_y, msg0_704, msg0_705, \
                         msf_379, msf_380, msf_470, msg1_704, msg1_705, nsf_469, \
                         nsf_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_12 * msf_379[k]
                   + f_3 * pc_y[k] * nsf_469[k];

        t_704[k] = pa_x[k] * msg0_704[k]
                   - f_6 * pc_x[k] * msg1_704[k];

        t_705[k] = pa_x[k] * msg0_705[k]
                   + f_16 * msf_470[k]
                   - f_6 * pc_x[k] * msg1_705[k];

        t_706[k] = f_13 * msf_380[k]
                   + f_3 * pc_y[k] * nsf_470[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, pa_x, pc_x, pc_y, pc_z, msg0_708, msf_370, \
                         msf_382, msf_473, msg1_708, nsf_470, nsf_472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = f_8 * msf_370[k]
                   + f_3 * pc_z[k] * nsf_470[k];

        t_708[k] = pa_x[k] * msg0_708[k]
                   + f_8 * msf_473[k]
                   - f_6 * pc_x[k] * msg1_708[k];

        t_709[k] = f_13 * msf_382[k]
                   + f_3 * pc_y[k] * nsf_472[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pa_x, pc_x, msg0_710, msf_475, msf_476, \
                         msf_477, msf_478, msg1_710, nsf_476, nsf_477, \
                         nsf_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = pa_x[k] * msg0_710[k]
                   + f_8 * msf_475[k]
                   - f_6 * pc_x[k] * msg1_710[k];

        t_711[k] = f_7 * msf_476[k]
                   + f_3 * pc_x[k] * nsf_476[k];

        t_712[k] = f_7 * msf_477[k]
                   + f_3 * pc_x[k] * nsf_477[k];

        t_713[k] = f_7 * msf_478[k]
                   + f_3 * pc_x[k] * nsf_478[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, pa_x, pc_x, pc_z, msg0_715, msg0_717, \
                         msf_376, msf_479, msg1_715, msg1_717, nsf_476, \
                         nsf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_7 * msf_479[k]
                   + f_3 * pc_x[k] * nsf_479[k];

        t_715[k] = pa_x[k] * msg0_715[k]
                   - f_6 * pc_x[k] * msg1_715[k];

        t_716[k] = f_8 * msf_376[k]
                   + f_3 * pc_z[k] * nsf_476[k];

        t_717[k] = pa_x[k] * msg0_717[k]
                   - f_6 * pc_x[k] * msg1_717[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, pa_x, pc_x, pc_y, msg0_719, msg0_720, \
                         msf_389, msf_390, msf_480, msg1_719, msg1_720, nsf_479, \
                         nsf_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_13 * msf_389[k]
                   + f_3 * pc_y[k] * nsf_479[k];

        t_719[k] = pa_x[k] * msg0_719[k]
                   - f_6 * pc_x[k] * msg1_719[k];

        t_720[k] = pa_x[k] * msg0_720[k]
                   + f_16 * msf_480[k]
                   - f_6 * pc_x[k] * msg1_720[k];

        t_721[k] = f_15 * msf_390[k]
                   + f_3 * pc_y[k] * nsf_480[k];
    }
}

static auto
compute_prim_nsg_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msg0,
                                                          const size_t msf, const size_t msg1,
                                                          const size_t nsd0, const size_t nsd1,
                                                          const size_t nsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msg0_660 = buffer.data(msg0 + 660);
    const auto *msg0_665 = buffer.data(msg0 + 665);
    const auto *msg0_675 = buffer.data(msg0 + 675);
    const auto *msg0_676 = buffer.data(msg0 + 676);
    const auto *msg0_678 = buffer.data(msg0 + 678);
    const auto *msg0_723 = buffer.data(msg0 + 723);
    const auto *msg0_725 = buffer.data(msg0 + 725);
    const auto *msg0_730 = buffer.data(msg0 + 730);
    const auto *msg0_732 = buffer.data(msg0 + 732);
    const auto *msg0_734 = buffer.data(msg0 + 734);
    const auto *msg0_735 = buffer.data(msg0 + 735);
    const auto *msg0_738 = buffer.data(msg0 + 738);
    const auto *msg0_740 = buffer.data(msg0 + 740);
    const auto *msg0_745 = buffer.data(msg0 + 745);
    const auto *msg0_747 = buffer.data(msg0 + 747);
    const auto *msg0_749 = buffer.data(msg0 + 749);
    const auto *msg0_750 = buffer.data(msg0 + 750);
    const auto *msg0_753 = buffer.data(msg0 + 753);
    const auto *msg0_755 = buffer.data(msg0 + 755);
    const auto *msg0_760 = buffer.data(msg0 + 760);
    const auto *msg0_762 = buffer.data(msg0 + 762);
    const auto *msg0_764 = buffer.data(msg0 + 764);
    const auto *msg0_765 = buffer.data(msg0 + 765);
    const auto *msg0_768 = buffer.data(msg0 + 768);
    const auto *msg0_770 = buffer.data(msg0 + 770);
    const auto *msg0_775 = buffer.data(msg0 + 775);
    const auto *msg0_777 = buffer.data(msg0 + 777);
    const auto *msg0_779 = buffer.data(msg0 + 779);
    const auto *msg0_780 = buffer.data(msg0 + 780);
    const auto *msg0_783 = buffer.data(msg0 + 783);
    const auto *msg0_785 = buffer.data(msg0 + 785);
    const auto *msg0_790 = buffer.data(msg0 + 790);
    const auto *msg0_792 = buffer.data(msg0 + 792);
    const auto *msg0_794 = buffer.data(msg0 + 794);
    const auto *msg0_798 = buffer.data(msg0 + 798);
    const auto *msg0_805 = buffer.data(msg0 + 805);
    const auto *msg0_807 = buffer.data(msg0 + 807);
    const auto *msg0_809 = buffer.data(msg0 + 809);
    const auto *msg0_810 = buffer.data(msg0 + 810);
    const auto *msg0_815 = buffer.data(msg0 + 815);
    const auto *msg0_820 = buffer.data(msg0 + 820);
    const auto *msg0_821 = buffer.data(msg0 + 821);
    const auto *msg0_822 = buffer.data(msg0 + 822);
    const auto *msg0_824 = buffer.data(msg0 + 824);

    const auto *msf_380 = buffer.data(msf + 380);
    const auto *msf_386 = buffer.data(msf + 386);
    const auto *msf_390 = buffer.data(msf + 390);
    const auto *msf_392 = buffer.data(msf + 392);
    const auto *msf_396 = buffer.data(msf + 396);
    const auto *msf_399 = buffer.data(msf + 399);
    const auto *msf_400 = buffer.data(msf + 400);
    const auto *msf_402 = buffer.data(msf + 402);
    const auto *msf_406 = buffer.data(msf + 406);
    const auto *msf_409 = buffer.data(msf + 409);
    const auto *msf_410 = buffer.data(msf + 410);
    const auto *msf_412 = buffer.data(msf + 412);
    const auto *msf_416 = buffer.data(msf + 416);
    const auto *msf_419 = buffer.data(msf + 419);
    const auto *msf_420 = buffer.data(msf + 420);
    const auto *msf_422 = buffer.data(msf + 422);
    const auto *msf_426 = buffer.data(msf + 426);
    const auto *msf_429 = buffer.data(msf + 429);
    const auto *msf_430 = buffer.data(msf + 430);
    const auto *msf_432 = buffer.data(msf + 432);
    const auto *msf_436 = buffer.data(msf + 436);
    const auto *msf_439 = buffer.data(msf + 439);
    const auto *msf_440 = buffer.data(msf + 440);
    const auto *msf_442 = buffer.data(msf + 442);
    const auto *msf_449 = buffer.data(msf + 449);
    const auto *msf_456 = buffer.data(msf + 456);
    const auto *msf_459 = buffer.data(msf + 459);
    const auto *msf_483 = buffer.data(msf + 483);
    const auto *msf_485 = buffer.data(msf + 485);
    const auto *msf_486 = buffer.data(msf + 486);
    const auto *msf_487 = buffer.data(msf + 487);
    const auto *msf_488 = buffer.data(msf + 488);
    const auto *msf_489 = buffer.data(msf + 489);
    const auto *msf_490 = buffer.data(msf + 490);
    const auto *msf_493 = buffer.data(msf + 493);
    const auto *msf_495 = buffer.data(msf + 495);
    const auto *msf_496 = buffer.data(msf + 496);
    const auto *msf_497 = buffer.data(msf + 497);
    const auto *msf_498 = buffer.data(msf + 498);
    const auto *msf_499 = buffer.data(msf + 499);
    const auto *msf_500 = buffer.data(msf + 500);
    const auto *msf_503 = buffer.data(msf + 503);
    const auto *msf_505 = buffer.data(msf + 505);
    const auto *msf_506 = buffer.data(msf + 506);
    const auto *msf_507 = buffer.data(msf + 507);
    const auto *msf_508 = buffer.data(msf + 508);
    const auto *msf_509 = buffer.data(msf + 509);
    const auto *msf_510 = buffer.data(msf + 510);
    const auto *msf_513 = buffer.data(msf + 513);
    const auto *msf_515 = buffer.data(msf + 515);
    const auto *msf_516 = buffer.data(msf + 516);
    const auto *msf_517 = buffer.data(msf + 517);
    const auto *msf_518 = buffer.data(msf + 518);
    const auto *msf_519 = buffer.data(msf + 519);
    const auto *msf_520 = buffer.data(msf + 520);
    const auto *msf_523 = buffer.data(msf + 523);
    const auto *msf_525 = buffer.data(msf + 525);
    const auto *msf_526 = buffer.data(msf + 526);
    const auto *msf_527 = buffer.data(msf + 527);
    const auto *msf_528 = buffer.data(msf + 528);
    const auto *msf_529 = buffer.data(msf + 529);
    const auto *msf_533 = buffer.data(msf + 533);
    const auto *msf_536 = buffer.data(msf + 536);
    const auto *msf_537 = buffer.data(msf + 537);
    const auto *msf_538 = buffer.data(msf + 538);
    const auto *msf_539 = buffer.data(msf + 539);
    const auto *msf_540 = buffer.data(msf + 540);
    const auto *msf_545 = buffer.data(msf + 545);
    const auto *msf_546 = buffer.data(msf + 546);
    const auto *msf_547 = buffer.data(msf + 547);
    const auto *msf_549 = buffer.data(msf + 549);

    const auto *msg1_660 = buffer.data(msg1 + 660);
    const auto *msg1_665 = buffer.data(msg1 + 665);
    const auto *msg1_675 = buffer.data(msg1 + 675);
    const auto *msg1_676 = buffer.data(msg1 + 676);
    const auto *msg1_678 = buffer.data(msg1 + 678);
    const auto *msg1_723 = buffer.data(msg1 + 723);
    const auto *msg1_725 = buffer.data(msg1 + 725);
    const auto *msg1_730 = buffer.data(msg1 + 730);
    const auto *msg1_732 = buffer.data(msg1 + 732);
    const auto *msg1_734 = buffer.data(msg1 + 734);
    const auto *msg1_735 = buffer.data(msg1 + 735);
    const auto *msg1_738 = buffer.data(msg1 + 738);
    const auto *msg1_740 = buffer.data(msg1 + 740);
    const auto *msg1_745 = buffer.data(msg1 + 745);
    const auto *msg1_747 = buffer.data(msg1 + 747);
    const auto *msg1_749 = buffer.data(msg1 + 749);
    const auto *msg1_750 = buffer.data(msg1 + 750);
    const auto *msg1_753 = buffer.data(msg1 + 753);
    const auto *msg1_755 = buffer.data(msg1 + 755);
    const auto *msg1_760 = buffer.data(msg1 + 760);
    const auto *msg1_762 = buffer.data(msg1 + 762);
    const auto *msg1_764 = buffer.data(msg1 + 764);
    const auto *msg1_765 = buffer.data(msg1 + 765);
    const auto *msg1_768 = buffer.data(msg1 + 768);
    const auto *msg1_770 = buffer.data(msg1 + 770);
    const auto *msg1_775 = buffer.data(msg1 + 775);
    const auto *msg1_777 = buffer.data(msg1 + 777);
    const auto *msg1_779 = buffer.data(msg1 + 779);
    const auto *msg1_780 = buffer.data(msg1 + 780);
    const auto *msg1_783 = buffer.data(msg1 + 783);
    const auto *msg1_785 = buffer.data(msg1 + 785);
    const auto *msg1_790 = buffer.data(msg1 + 790);
    const auto *msg1_792 = buffer.data(msg1 + 792);
    const auto *msg1_794 = buffer.data(msg1 + 794);
    const auto *msg1_798 = buffer.data(msg1 + 798);
    const auto *msg1_805 = buffer.data(msg1 + 805);
    const auto *msg1_807 = buffer.data(msg1 + 807);
    const auto *msg1_809 = buffer.data(msg1 + 809);
    const auto *msg1_810 = buffer.data(msg1 + 810);
    const auto *msg1_815 = buffer.data(msg1 + 815);
    const auto *msg1_820 = buffer.data(msg1 + 820);
    const auto *msg1_821 = buffer.data(msg1 + 821);
    const auto *msg1_822 = buffer.data(msg1 + 822);
    const auto *msg1_824 = buffer.data(msg1 + 824);

    const auto *nsd0_324 = buffer.data(nsd0 + 324);
    const auto *nsd0_330 = buffer.data(nsd0 + 330);
    const auto *nsd0_331 = buffer.data(nsd0 + 331);
    const auto *nsd0_333 = buffer.data(nsd0 + 333);
    const auto *nsd0_335 = buffer.data(nsd0 + 335);
    const auto *nsd0_338 = buffer.data(nsd0 + 338);
    const auto *nsd0_340 = buffer.data(nsd0 + 340);
    const auto *nsd0_341 = buffer.data(nsd0 + 341);

    const auto *nsd1_324 = buffer.data(nsd1 + 324);
    const auto *nsd1_330 = buffer.data(nsd1 + 330);
    const auto *nsd1_331 = buffer.data(nsd1 + 331);
    const auto *nsd1_333 = buffer.data(nsd1 + 333);
    const auto *nsd1_335 = buffer.data(nsd1 + 335);
    const auto *nsd1_338 = buffer.data(nsd1 + 338);
    const auto *nsd1_340 = buffer.data(nsd1 + 340);
    const auto *nsd1_341 = buffer.data(nsd1 + 341);

    const auto *nsf_480 = buffer.data(nsf + 480);
    const auto *nsf_482 = buffer.data(nsf + 482);
    const auto *nsf_486 = buffer.data(nsf + 486);
    const auto *nsf_487 = buffer.data(nsf + 487);
    const auto *nsf_488 = buffer.data(nsf + 488);
    const auto *nsf_489 = buffer.data(nsf + 489);
    const auto *nsf_490 = buffer.data(nsf + 490);
    const auto *nsf_492 = buffer.data(nsf + 492);
    const auto *nsf_496 = buffer.data(nsf + 496);
    const auto *nsf_497 = buffer.data(nsf + 497);
    const auto *nsf_498 = buffer.data(nsf + 498);
    const auto *nsf_499 = buffer.data(nsf + 499);
    const auto *nsf_500 = buffer.data(nsf + 500);
    const auto *nsf_502 = buffer.data(nsf + 502);
    const auto *nsf_506 = buffer.data(nsf + 506);
    const auto *nsf_507 = buffer.data(nsf + 507);
    const auto *nsf_508 = buffer.data(nsf + 508);
    const auto *nsf_509 = buffer.data(nsf + 509);
    const auto *nsf_510 = buffer.data(nsf + 510);
    const auto *nsf_512 = buffer.data(nsf + 512);
    const auto *nsf_516 = buffer.data(nsf + 516);
    const auto *nsf_517 = buffer.data(nsf + 517);
    const auto *nsf_518 = buffer.data(nsf + 518);
    const auto *nsf_519 = buffer.data(nsf + 519);
    const auto *nsf_520 = buffer.data(nsf + 520);
    const auto *nsf_522 = buffer.data(nsf + 522);
    const auto *nsf_526 = buffer.data(nsf + 526);
    const auto *nsf_527 = buffer.data(nsf + 527);
    const auto *nsf_528 = buffer.data(nsf + 528);
    const auto *nsf_529 = buffer.data(nsf + 529);
    const auto *nsf_530 = buffer.data(nsf + 530);
    const auto *nsf_532 = buffer.data(nsf + 532);
    const auto *nsf_536 = buffer.data(nsf + 536);
    const auto *nsf_537 = buffer.data(nsf + 537);
    const auto *nsf_538 = buffer.data(nsf + 538);
    const auto *nsf_539 = buffer.data(nsf + 539);
    const auto *nsf_540 = buffer.data(nsf + 540);
    const auto *nsf_541 = buffer.data(nsf + 541);
    const auto *nsf_542 = buffer.data(nsf + 542);
    const auto *nsf_545 = buffer.data(nsf + 545);
    const auto *nsf_546 = buffer.data(nsf + 546);
    const auto *nsf_547 = buffer.data(nsf + 547);
    const auto *nsf_549 = buffer.data(nsf + 549);
    const auto *nsf_550 = buffer.data(nsf + 550);
    const auto *nsf_551 = buffer.data(nsf + 551);
    const auto *nsf_553 = buffer.data(nsf + 553);
    const auto *nsf_555 = buffer.data(nsf + 555);
    const auto *nsf_556 = buffer.data(nsf + 556);
    const auto *nsf_557 = buffer.data(nsf + 557);
    const auto *nsf_558 = buffer.data(nsf + 558);
    const auto *nsf_559 = buffer.data(nsf + 559);
    const auto *nsf_562 = buffer.data(nsf + 562);
    const auto *nsf_564 = buffer.data(nsf + 564);
    const auto *nsf_565 = buffer.data(nsf + 565);
    const auto *nsf_566 = buffer.data(nsf + 566);
    const auto *nsf_567 = buffer.data(nsf + 567);
    const auto *nsf_568 = buffer.data(nsf + 568);
    const auto *nsf_569 = buffer.data(nsf + 569);

#pragma omp simd aligned(t_722, t_723, t_724, pa_x, pc_x, pc_y, pc_z, msg0_723, msf_380, \
                         msf_392, msf_483, msg1_723, nsf_480, nsf_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_14 * msf_380[k]
                   + f_3 * pc_z[k] * nsf_480[k];

        t_723[k] = pa_x[k] * msg0_723[k]
                   + f_8 * msf_483[k]
                   - f_6 * pc_x[k] * msg1_723[k];

        t_724[k] = f_15 * msf_392[k]
                   + f_3 * pc_y[k] * nsf_482[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pa_x, pc_x, msg0_725, msf_485, msf_486, \
                         msf_487, msf_488, msg1_725, nsf_486, nsf_487, \
                         nsf_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = pa_x[k] * msg0_725[k]
                   + f_8 * msf_485[k]
                   - f_6 * pc_x[k] * msg1_725[k];

        t_726[k] = f_7 * msf_486[k]
                   + f_3 * pc_x[k] * nsf_486[k];

        t_727[k] = f_7 * msf_487[k]
                   + f_3 * pc_x[k] * nsf_487[k];

        t_728[k] = f_7 * msf_488[k]
                   + f_3 * pc_x[k] * nsf_488[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, pa_x, pc_x, pc_z, msg0_730, msg0_732, \
                         msf_386, msf_489, msg1_730, msg1_732, nsf_486, \
                         nsf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_7 * msf_489[k]
                   + f_3 * pc_x[k] * nsf_489[k];

        t_730[k] = pa_x[k] * msg0_730[k]
                   - f_6 * pc_x[k] * msg1_730[k];

        t_731[k] = f_14 * msf_386[k]
                   + f_3 * pc_z[k] * nsf_486[k];

        t_732[k] = pa_x[k] * msg0_732[k]
                   - f_6 * pc_x[k] * msg1_732[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, t_736, pa_x, pc_x, pc_y, msg0_734, msg0_735, \
                         msf_399, msf_400, msf_490, msg1_734, msg1_735, nsf_489, \
                         nsf_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_15 * msf_399[k]
                   + f_3 * pc_y[k] * nsf_489[k];

        t_734[k] = pa_x[k] * msg0_734[k]
                   - f_6 * pc_x[k] * msg1_734[k];

        t_735[k] = pa_x[k] * msg0_735[k]
                   + f_16 * msf_490[k]
                   - f_6 * pc_x[k] * msg1_735[k];

        t_736[k] = f_17 * msf_400[k]
                   + f_3 * pc_y[k] * nsf_490[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pa_x, pc_x, pc_y, pc_z, msg0_738, msf_390, \
                         msf_402, msf_493, msg1_738, nsf_490, nsf_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_16 * msf_390[k]
                   + f_3 * pc_z[k] * nsf_490[k];

        t_738[k] = pa_x[k] * msg0_738[k]
                   + f_8 * msf_493[k]
                   - f_6 * pc_x[k] * msg1_738[k];

        t_739[k] = f_17 * msf_402[k]
                   + f_3 * pc_y[k] * nsf_492[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pa_x, pc_x, msg0_740, msf_495, msf_496, \
                         msf_497, msf_498, msg1_740, nsf_496, nsf_497, \
                         nsf_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = pa_x[k] * msg0_740[k]
                   + f_8 * msf_495[k]
                   - f_6 * pc_x[k] * msg1_740[k];

        t_741[k] = f_7 * msf_496[k]
                   + f_3 * pc_x[k] * nsf_496[k];

        t_742[k] = f_7 * msf_497[k]
                   + f_3 * pc_x[k] * nsf_497[k];

        t_743[k] = f_7 * msf_498[k]
                   + f_3 * pc_x[k] * nsf_498[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, pa_x, pc_x, pc_z, msg0_745, msg0_747, \
                         msf_396, msf_499, msg1_745, msg1_747, nsf_496, \
                         nsf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_7 * msf_499[k]
                   + f_3 * pc_x[k] * nsf_499[k];

        t_745[k] = pa_x[k] * msg0_745[k]
                   - f_6 * pc_x[k] * msg1_745[k];

        t_746[k] = f_16 * msf_396[k]
                   + f_3 * pc_z[k] * nsf_496[k];

        t_747[k] = pa_x[k] * msg0_747[k]
                   - f_6 * pc_x[k] * msg1_747[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pa_x, pc_x, pc_y, msg0_749, msg0_750, \
                         msf_409, msf_410, msf_500, msg1_749, msg1_750, nsf_499, \
                         nsf_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_17 * msf_409[k]
                   + f_3 * pc_y[k] * nsf_499[k];

        t_749[k] = pa_x[k] * msg0_749[k]
                   - f_6 * pc_x[k] * msg1_749[k];

        t_750[k] = pa_x[k] * msg0_750[k]
                   + f_16 * msf_500[k]
                   - f_6 * pc_x[k] * msg1_750[k];

        t_751[k] = f_16 * msf_410[k]
                   + f_3 * pc_y[k] * nsf_500[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pa_x, pc_x, pc_y, pc_z, msg0_753, msf_400, \
                         msf_412, msf_503, msg1_753, nsf_500, nsf_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_17 * msf_400[k]
                   + f_3 * pc_z[k] * nsf_500[k];

        t_753[k] = pa_x[k] * msg0_753[k]
                   + f_8 * msf_503[k]
                   - f_6 * pc_x[k] * msg1_753[k];

        t_754[k] = f_16 * msf_412[k]
                   + f_3 * pc_y[k] * nsf_502[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pa_x, pc_x, msg0_755, msf_505, msf_506, \
                         msf_507, msf_508, msg1_755, nsf_506, nsf_507, \
                         nsf_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pa_x[k] * msg0_755[k]
                   + f_8 * msf_505[k]
                   - f_6 * pc_x[k] * msg1_755[k];

        t_756[k] = f_7 * msf_506[k]
                   + f_3 * pc_x[k] * nsf_506[k];

        t_757[k] = f_7 * msf_507[k]
                   + f_3 * pc_x[k] * nsf_507[k];

        t_758[k] = f_7 * msf_508[k]
                   + f_3 * pc_x[k] * nsf_508[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, t_762, pa_x, pc_x, pc_z, msg0_760, msg0_762, \
                         msf_406, msf_509, msg1_760, msg1_762, nsf_506, \
                         nsf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_7 * msf_509[k]
                   + f_3 * pc_x[k] * nsf_509[k];

        t_760[k] = pa_x[k] * msg0_760[k]
                   - f_6 * pc_x[k] * msg1_760[k];

        t_761[k] = f_17 * msf_406[k]
                   + f_3 * pc_z[k] * nsf_506[k];

        t_762[k] = pa_x[k] * msg0_762[k]
                   - f_6 * pc_x[k] * msg1_762[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, t_766, pa_x, pc_x, pc_y, msg0_764, msg0_765, \
                         msf_419, msf_420, msf_510, msg1_764, msg1_765, nsf_509, \
                         nsf_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_16 * msf_419[k]
                   + f_3 * pc_y[k] * nsf_509[k];

        t_764[k] = pa_x[k] * msg0_764[k]
                   - f_6 * pc_x[k] * msg1_764[k];

        t_765[k] = pa_x[k] * msg0_765[k]
                   + f_16 * msf_510[k]
                   - f_6 * pc_x[k] * msg1_765[k];

        t_766[k] = f_14 * msf_420[k]
                   + f_3 * pc_y[k] * nsf_510[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pa_x, pc_x, pc_y, pc_z, msg0_768, msf_410, \
                         msf_422, msf_513, msg1_768, nsf_510, nsf_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_15 * msf_410[k]
                   + f_3 * pc_z[k] * nsf_510[k];

        t_768[k] = pa_x[k] * msg0_768[k]
                   + f_8 * msf_513[k]
                   - f_6 * pc_x[k] * msg1_768[k];

        t_769[k] = f_14 * msf_422[k]
                   + f_3 * pc_y[k] * nsf_512[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, pa_x, pc_x, msg0_770, msf_515, msf_516, \
                         msf_517, msf_518, msg1_770, nsf_516, nsf_517, \
                         nsf_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = pa_x[k] * msg0_770[k]
                   + f_8 * msf_515[k]
                   - f_6 * pc_x[k] * msg1_770[k];

        t_771[k] = f_7 * msf_516[k]
                   + f_3 * pc_x[k] * nsf_516[k];

        t_772[k] = f_7 * msf_517[k]
                   + f_3 * pc_x[k] * nsf_517[k];

        t_773[k] = f_7 * msf_518[k]
                   + f_3 * pc_x[k] * nsf_518[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, pa_x, pc_x, pc_z, msg0_775, msg0_777, \
                         msf_416, msf_519, msg1_775, msg1_777, nsf_516, \
                         nsf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_7 * msf_519[k]
                   + f_3 * pc_x[k] * nsf_519[k];

        t_775[k] = pa_x[k] * msg0_775[k]
                   - f_6 * pc_x[k] * msg1_775[k];

        t_776[k] = f_15 * msf_416[k]
                   + f_3 * pc_z[k] * nsf_516[k];

        t_777[k] = pa_x[k] * msg0_777[k]
                   - f_6 * pc_x[k] * msg1_777[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, pa_x, pc_x, pc_y, msg0_779, msg0_780, \
                         msf_429, msf_430, msf_520, msg1_779, msg1_780, nsf_519, \
                         nsf_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_14 * msf_429[k]
                   + f_3 * pc_y[k] * nsf_519[k];

        t_779[k] = pa_x[k] * msg0_779[k]
                   - f_6 * pc_x[k] * msg1_779[k];

        t_780[k] = pa_x[k] * msg0_780[k]
                   + f_16 * msf_520[k]
                   - f_6 * pc_x[k] * msg1_780[k];

        t_781[k] = f_8 * msf_430[k]
                   + f_3 * pc_y[k] * nsf_520[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pa_x, pc_x, pc_y, pc_z, msg0_783, msf_420, \
                         msf_432, msf_523, msg1_783, nsf_520, nsf_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_13 * msf_420[k]
                   + f_3 * pc_z[k] * nsf_520[k];

        t_783[k] = pa_x[k] * msg0_783[k]
                   + f_8 * msf_523[k]
                   - f_6 * pc_x[k] * msg1_783[k];

        t_784[k] = f_8 * msf_432[k]
                   + f_3 * pc_y[k] * nsf_522[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, pa_x, pc_x, msg0_785, msf_525, msf_526, \
                         msf_527, msf_528, msg1_785, nsf_526, nsf_527, \
                         nsf_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = pa_x[k] * msg0_785[k]
                   + f_8 * msf_525[k]
                   - f_6 * pc_x[k] * msg1_785[k];

        t_786[k] = f_7 * msf_526[k]
                   + f_3 * pc_x[k] * nsf_526[k];

        t_787[k] = f_7 * msf_527[k]
                   + f_3 * pc_x[k] * nsf_527[k];

        t_788[k] = f_7 * msf_528[k]
                   + f_3 * pc_x[k] * nsf_528[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, pa_x, pc_x, pc_z, msg0_790, msg0_792, \
                         msf_426, msf_529, msg1_790, msg1_792, nsf_526, \
                         nsf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_7 * msf_529[k]
                   + f_3 * pc_x[k] * nsf_529[k];

        t_790[k] = pa_x[k] * msg0_790[k]
                   - f_6 * pc_x[k] * msg1_790[k];

        t_791[k] = f_13 * msf_426[k]
                   + f_3 * pc_z[k] * nsf_526[k];

        t_792[k] = pa_x[k] * msg0_792[k]
                   - f_6 * pc_x[k] * msg1_792[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, t_796, pa_x, pa_y, pc_x, pc_y, msg0_660, \
                         msg0_794, msf_439, msf_440, msg1_660, msg1_794, nsf_529, \
                         nsf_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = f_8 * msf_439[k]
                   + f_3 * pc_y[k] * nsf_529[k];

        t_794[k] = pa_x[k] * msg0_794[k]
                   - f_6 * pc_x[k] * msg1_794[k];

        t_795[k] = pa_y[k] * msg0_660[k]
                   - f_6 * pc_y[k] * msg1_660[k];

        t_796[k] = f_7 * msf_440[k]
                   + f_3 * pc_y[k] * nsf_530[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pa_x, pc_x, pc_y, pc_z, msg0_798, msf_430, \
                         msf_442, msf_533, msg1_798, nsf_530, nsf_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_12 * msf_430[k]
                   + f_3 * pc_z[k] * nsf_530[k];

        t_798[k] = pa_x[k] * msg0_798[k]
                   + f_8 * msf_533[k]
                   - f_6 * pc_x[k] * msg1_798[k];

        t_799[k] = f_7 * msf_442[k]
                   + f_3 * pc_y[k] * nsf_532[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, pa_y, pc_x, pc_y, msg0_665, msf_536, \
                         msf_537, msf_538, msg1_665, nsf_536, nsf_537, \
                         nsf_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = pa_y[k] * msg0_665[k]
                   - f_6 * pc_y[k] * msg1_665[k];

        t_801[k] = f_7 * msf_536[k]
                   + f_3 * pc_x[k] * nsf_536[k];

        t_802[k] = f_7 * msf_537[k]
                   + f_3 * pc_x[k] * nsf_537[k];

        t_803[k] = f_7 * msf_538[k]
                   + f_3 * pc_x[k] * nsf_538[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, t_807, pa_x, pc_x, pc_z, msg0_805, msg0_807, \
                         msf_436, msf_539, msg1_805, msg1_807, nsf_536, \
                         nsf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_7 * msf_539[k]
                   + f_3 * pc_x[k] * nsf_539[k];

        t_805[k] = pa_x[k] * msg0_805[k]
                   - f_6 * pc_x[k] * msg1_805[k];

        t_806[k] = f_12 * msf_436[k]
                   + f_3 * pc_z[k] * nsf_536[k];

        t_807[k] = pa_x[k] * msg0_807[k]
                   - f_6 * pc_x[k] * msg1_807[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pa_x, pc_x, pc_y, msg0_809, msg0_810, \
                         msf_449, msf_540, msg1_809, msg1_810, nsf_539, \
                         nsf_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_7 * msf_449[k]
                   + f_3 * pc_y[k] * nsf_539[k];

        t_809[k] = pa_x[k] * msg0_809[k]
                   - f_6 * pc_x[k] * msg1_809[k];

        t_810[k] = pa_x[k] * msg0_810[k]
                   + f_16 * msf_540[k]
                   - f_6 * pc_x[k] * msg1_810[k];

        t_811[k] = f_3 * pc_y[k] * nsf_540[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, pc_y, pc_z, msf_440, nsd0_324, nsd1_324, \
                         nsf_540, nsf_541, nsf_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_9 * msf_440[k]
                   + f_3 * pc_z[k] * nsf_540[k];

        t_813[k] = f_4 * nsd0_324[k]
                   - f_5 * nsd1_324[k]
                   + f_3 * pc_y[k] * nsf_541[k];

        t_814[k] = f_3 * pc_y[k] * nsf_542[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, pa_x, pc_x, pc_y, msg0_815, msf_545, \
                         msf_546, msf_547, msg1_815, nsf_545, nsf_546, \
                         nsf_547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = pa_x[k] * msg0_815[k]
                   + f_8 * msf_545[k]
                   - f_6 * pc_x[k] * msg1_815[k];

        t_816[k] = f_7 * msf_546[k]
                   + f_3 * pc_x[k] * nsf_546[k];

        t_817[k] = f_7 * msf_547[k]
                   + f_3 * pc_x[k] * nsf_547[k];

        t_818[k] = f_3 * pc_y[k] * nsf_545[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, t_823, pa_x, pc_x, pc_y, msg0_820, \
                         msg0_821, msg0_822, msf_549, msg1_820, msg1_821, msg1_822, \
                         nsf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_7 * msf_549[k]
                   + f_3 * pc_x[k] * nsf_549[k];

        t_820[k] = pa_x[k] * msg0_820[k]
                   - f_6 * pc_x[k] * msg1_820[k];

        t_821[k] = pa_x[k] * msg0_821[k]
                   - f_6 * pc_x[k] * msg1_821[k];

        t_822[k] = pa_x[k] * msg0_822[k]
                   - f_6 * pc_x[k] * msg1_822[k];

        t_823[k] = f_3 * pc_y[k] * nsf_549[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, t_827, pa_x, pc_x, pc_z, msg0_824, msg1_824, \
                         nsd0_330, nsd0_331, nsd1_330, nsd1_331, nsf_550, \
                         nsf_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = pa_x[k] * msg0_824[k]
                   - f_6 * pc_x[k] * msg1_824[k];

        t_825[k] = f_1 * nsd0_330[k]
                   - f_2 * nsd1_330[k]
                   + f_3 * pc_x[k] * nsf_550[k];

        t_826[k] = f_10 * nsd0_331[k]
                   - f_11 * nsd1_331[k]
                   + f_3 * pc_x[k] * nsf_551[k];

        t_827[k] = f_3 * pc_z[k] * nsf_550[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, t_831, t_832, pc_x, pc_z, nsd0_333, nsd0_335, \
                         nsd1_333, nsd1_335, nsf_551, nsf_553, nsf_555, nsf_556, \
                         nsf_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = f_4 * nsd0_333[k]
                   - f_5 * nsd1_333[k]
                   + f_3 * pc_x[k] * nsf_553[k];

        t_829[k] = f_3 * pc_z[k] * nsf_551[k];

        t_830[k] = f_4 * nsd0_335[k]
                   - f_5 * nsd1_335[k]
                   + f_3 * pc_x[k] * nsf_555[k];

        t_831[k] = f_3 * pc_x[k] * nsf_556[k];

        t_832[k] = f_3 * pc_x[k] * nsf_557[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, t_836, t_837, pc_x, pc_y, pc_z, msf_456, \
                         nsd0_333, nsd1_333, nsf_556, nsf_557, nsf_558, \
                         nsf_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_3 * pc_x[k] * nsf_558[k];

        t_834[k] = f_3 * pc_x[k] * nsf_559[k];

        t_835[k] = f_0 * msf_456[k]
                   + f_1 * nsd0_333[k]
                   - f_2 * nsd1_333[k]
                   + f_3 * pc_y[k] * nsf_556[k];

        t_836[k] = f_3 * pc_z[k] * nsf_556[k];

        t_837[k] = f_4 * nsd0_333[k]
                   - f_5 * nsd1_333[k]
                   + f_3 * pc_z[k] * nsf_557[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, t_841, pa_z, pc_y, pc_z, msg0_675, msg0_676, \
                         msf_459, msg1_675, msg1_676, nsd0_335, nsd1_335, \
                         nsf_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_0 * msf_459[k]
                   + f_3 * pc_y[k] * nsf_559[k];

        t_839[k] = f_1 * nsd0_335[k]
                   - f_2 * nsd1_335[k]
                   + f_3 * pc_z[k] * nsf_559[k];

        t_840[k] = pa_z[k] * msg0_675[k]
                   - f_6 * pc_z[k] * msg1_675[k];

        t_841[k] = pa_z[k] * msg0_676[k]
                   - f_6 * pc_z[k] * msg1_676[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pa_z, pc_x, pc_z, msg0_678, msg1_678, nsd0_338, \
                         nsd0_340, nsd1_338, nsd1_340, nsf_562, \
                         nsf_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_10 * nsd0_338[k]
                   - f_11 * nsd1_338[k]
                   + f_3 * pc_x[k] * nsf_562[k];

        t_843[k] = pa_z[k] * msg0_678[k]
                   - f_6 * pc_z[k] * msg1_678[k];

        t_844[k] = f_4 * nsd0_340[k]
                   - f_5 * nsd1_340[k]
                   + f_3 * pc_x[k] * nsf_564[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, t_849, pc_x, nsd0_341, nsd1_341, nsf_565, \
                         nsf_566, nsf_567, nsf_568, nsf_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_4 * nsd0_341[k]
                   - f_5 * nsd1_341[k]
                   + f_3 * pc_x[k] * nsf_565[k];

        t_846[k] = f_3 * pc_x[k] * nsf_566[k];

        t_847[k] = f_3 * pc_x[k] * nsf_567[k];

        t_848[k] = f_3 * pc_x[k] * nsf_568[k];

        t_849[k] = f_3 * pc_x[k] * nsf_569[k];
    }
}

static auto
compute_prim_nsg_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msg0,
                                                          const size_t msf, const size_t msg1,
                                                          const size_t nsd0, const size_t nsd1,
                                                          const size_t nsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.0 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msg0_685 = buffer.data(msg0 + 685);
    const auto *msg0_687 = buffer.data(msg0 + 687);
    const auto *msg0_810 = buffer.data(msg0 + 810);
    const auto *msg0_812 = buffer.data(msg0 + 812);
    const auto *msg0_815 = buffer.data(msg0 + 815);
    const auto *msg0_820 = buffer.data(msg0 + 820);

    const auto *msf_456 = buffer.data(msf + 456);
    const auto *msf_457 = buffer.data(msf + 457);
    const auto *msf_459 = buffer.data(msf + 459);
    const auto *msf_466 = buffer.data(msf + 466);
    const auto *msf_469 = buffer.data(msf + 469);
    const auto *msf_476 = buffer.data(msf + 476);
    const auto *msf_478 = buffer.data(msf + 478);
    const auto *msf_479 = buffer.data(msf + 479);
    const auto *msf_486 = buffer.data(msf + 486);
    const auto *msf_488 = buffer.data(msf + 488);
    const auto *msf_489 = buffer.data(msf + 489);
    const auto *msf_496 = buffer.data(msf + 496);
    const auto *msf_498 = buffer.data(msf + 498);
    const auto *msf_499 = buffer.data(msf + 499);
    const auto *msf_506 = buffer.data(msf + 506);
    const auto *msf_508 = buffer.data(msf + 508);
    const auto *msf_509 = buffer.data(msf + 509);
    const auto *msf_516 = buffer.data(msf + 516);
    const auto *msf_518 = buffer.data(msf + 518);
    const auto *msf_519 = buffer.data(msf + 519);
    const auto *msf_526 = buffer.data(msf + 526);
    const auto *msf_528 = buffer.data(msf + 528);
    const auto *msf_529 = buffer.data(msf + 529);
    const auto *msf_536 = buffer.data(msf + 536);
    const auto *msf_538 = buffer.data(msf + 538);
    const auto *msf_539 = buffer.data(msf + 539);
    const auto *msf_546 = buffer.data(msf + 546);

    const auto *msg1_685 = buffer.data(msg1 + 685);
    const auto *msg1_687 = buffer.data(msg1 + 687);
    const auto *msg1_810 = buffer.data(msg1 + 810);
    const auto *msg1_812 = buffer.data(msg1 + 812);
    const auto *msg1_815 = buffer.data(msg1 + 815);
    const auto *msg1_820 = buffer.data(msg1 + 820);

    const auto *nsd0_341 = buffer.data(nsd0 + 341);
    const auto *nsd0_342 = buffer.data(nsd0 + 342);
    const auto *nsd0_343 = buffer.data(nsd0 + 343);
    const auto *nsd0_344 = buffer.data(nsd0 + 344);
    const auto *nsd0_345 = buffer.data(nsd0 + 345);
    const auto *nsd0_346 = buffer.data(nsd0 + 346);
    const auto *nsd0_347 = buffer.data(nsd0 + 347);
    const auto *nsd0_348 = buffer.data(nsd0 + 348);
    const auto *nsd0_349 = buffer.data(nsd0 + 349);
    const auto *nsd0_350 = buffer.data(nsd0 + 350);
    const auto *nsd0_351 = buffer.data(nsd0 + 351);
    const auto *nsd0_352 = buffer.data(nsd0 + 352);
    const auto *nsd0_353 = buffer.data(nsd0 + 353);
    const auto *nsd0_354 = buffer.data(nsd0 + 354);
    const auto *nsd0_355 = buffer.data(nsd0 + 355);
    const auto *nsd0_356 = buffer.data(nsd0 + 356);
    const auto *nsd0_357 = buffer.data(nsd0 + 357);
    const auto *nsd0_358 = buffer.data(nsd0 + 358);
    const auto *nsd0_359 = buffer.data(nsd0 + 359);
    const auto *nsd0_360 = buffer.data(nsd0 + 360);
    const auto *nsd0_361 = buffer.data(nsd0 + 361);
    const auto *nsd0_362 = buffer.data(nsd0 + 362);
    const auto *nsd0_363 = buffer.data(nsd0 + 363);
    const auto *nsd0_364 = buffer.data(nsd0 + 364);
    const auto *nsd0_365 = buffer.data(nsd0 + 365);
    const auto *nsd0_366 = buffer.data(nsd0 + 366);
    const auto *nsd0_367 = buffer.data(nsd0 + 367);
    const auto *nsd0_368 = buffer.data(nsd0 + 368);
    const auto *nsd0_369 = buffer.data(nsd0 + 369);
    const auto *nsd0_370 = buffer.data(nsd0 + 370);
    const auto *nsd0_371 = buffer.data(nsd0 + 371);
    const auto *nsd0_372 = buffer.data(nsd0 + 372);
    const auto *nsd0_373 = buffer.data(nsd0 + 373);
    const auto *nsd0_374 = buffer.data(nsd0 + 374);
    const auto *nsd0_375 = buffer.data(nsd0 + 375);
    const auto *nsd0_376 = buffer.data(nsd0 + 376);
    const auto *nsd0_377 = buffer.data(nsd0 + 377);
    const auto *nsd0_378 = buffer.data(nsd0 + 378);
    const auto *nsd0_379 = buffer.data(nsd0 + 379);
    const auto *nsd0_380 = buffer.data(nsd0 + 380);
    const auto *nsd0_381 = buffer.data(nsd0 + 381);
    const auto *nsd0_382 = buffer.data(nsd0 + 382);
    const auto *nsd0_383 = buffer.data(nsd0 + 383);
    const auto *nsd0_385 = buffer.data(nsd0 + 385);
    const auto *nsd0_387 = buffer.data(nsd0 + 387);
    const auto *nsd0_388 = buffer.data(nsd0 + 388);

    const auto *nsd1_341 = buffer.data(nsd1 + 341);
    const auto *nsd1_342 = buffer.data(nsd1 + 342);
    const auto *nsd1_343 = buffer.data(nsd1 + 343);
    const auto *nsd1_344 = buffer.data(nsd1 + 344);
    const auto *nsd1_345 = buffer.data(nsd1 + 345);
    const auto *nsd1_346 = buffer.data(nsd1 + 346);
    const auto *nsd1_347 = buffer.data(nsd1 + 347);
    const auto *nsd1_348 = buffer.data(nsd1 + 348);
    const auto *nsd1_349 = buffer.data(nsd1 + 349);
    const auto *nsd1_350 = buffer.data(nsd1 + 350);
    const auto *nsd1_351 = buffer.data(nsd1 + 351);
    const auto *nsd1_352 = buffer.data(nsd1 + 352);
    const auto *nsd1_353 = buffer.data(nsd1 + 353);
    const auto *nsd1_354 = buffer.data(nsd1 + 354);
    const auto *nsd1_355 = buffer.data(nsd1 + 355);
    const auto *nsd1_356 = buffer.data(nsd1 + 356);
    const auto *nsd1_357 = buffer.data(nsd1 + 357);
    const auto *nsd1_358 = buffer.data(nsd1 + 358);
    const auto *nsd1_359 = buffer.data(nsd1 + 359);
    const auto *nsd1_360 = buffer.data(nsd1 + 360);
    const auto *nsd1_361 = buffer.data(nsd1 + 361);
    const auto *nsd1_362 = buffer.data(nsd1 + 362);
    const auto *nsd1_363 = buffer.data(nsd1 + 363);
    const auto *nsd1_364 = buffer.data(nsd1 + 364);
    const auto *nsd1_365 = buffer.data(nsd1 + 365);
    const auto *nsd1_366 = buffer.data(nsd1 + 366);
    const auto *nsd1_367 = buffer.data(nsd1 + 367);
    const auto *nsd1_368 = buffer.data(nsd1 + 368);
    const auto *nsd1_369 = buffer.data(nsd1 + 369);
    const auto *nsd1_370 = buffer.data(nsd1 + 370);
    const auto *nsd1_371 = buffer.data(nsd1 + 371);
    const auto *nsd1_372 = buffer.data(nsd1 + 372);
    const auto *nsd1_373 = buffer.data(nsd1 + 373);
    const auto *nsd1_374 = buffer.data(nsd1 + 374);
    const auto *nsd1_375 = buffer.data(nsd1 + 375);
    const auto *nsd1_376 = buffer.data(nsd1 + 376);
    const auto *nsd1_377 = buffer.data(nsd1 + 377);
    const auto *nsd1_378 = buffer.data(nsd1 + 378);
    const auto *nsd1_379 = buffer.data(nsd1 + 379);
    const auto *nsd1_380 = buffer.data(nsd1 + 380);
    const auto *nsd1_381 = buffer.data(nsd1 + 381);
    const auto *nsd1_382 = buffer.data(nsd1 + 382);
    const auto *nsd1_383 = buffer.data(nsd1 + 383);
    const auto *nsd1_385 = buffer.data(nsd1 + 385);
    const auto *nsd1_387 = buffer.data(nsd1 + 387);
    const auto *nsd1_388 = buffer.data(nsd1 + 388);

    const auto *nsf_566 = buffer.data(nsf + 566);
    const auto *nsf_569 = buffer.data(nsf + 569);
    const auto *nsf_570 = buffer.data(nsf + 570);
    const auto *nsf_571 = buffer.data(nsf + 571);
    const auto *nsf_572 = buffer.data(nsf + 572);
    const auto *nsf_573 = buffer.data(nsf + 573);
    const auto *nsf_574 = buffer.data(nsf + 574);
    const auto *nsf_575 = buffer.data(nsf + 575);
    const auto *nsf_576 = buffer.data(nsf + 576);
    const auto *nsf_577 = buffer.data(nsf + 577);
    const auto *nsf_578 = buffer.data(nsf + 578);
    const auto *nsf_579 = buffer.data(nsf + 579);
    const auto *nsf_580 = buffer.data(nsf + 580);
    const auto *nsf_581 = buffer.data(nsf + 581);
    const auto *nsf_582 = buffer.data(nsf + 582);
    const auto *nsf_583 = buffer.data(nsf + 583);
    const auto *nsf_584 = buffer.data(nsf + 584);
    const auto *nsf_585 = buffer.data(nsf + 585);
    const auto *nsf_586 = buffer.data(nsf + 586);
    const auto *nsf_587 = buffer.data(nsf + 587);
    const auto *nsf_588 = buffer.data(nsf + 588);
    const auto *nsf_589 = buffer.data(nsf + 589);
    const auto *nsf_590 = buffer.data(nsf + 590);
    const auto *nsf_591 = buffer.data(nsf + 591);
    const auto *nsf_592 = buffer.data(nsf + 592);
    const auto *nsf_593 = buffer.data(nsf + 593);
    const auto *nsf_594 = buffer.data(nsf + 594);
    const auto *nsf_595 = buffer.data(nsf + 595);
    const auto *nsf_596 = buffer.data(nsf + 596);
    const auto *nsf_597 = buffer.data(nsf + 597);
    const auto *nsf_598 = buffer.data(nsf + 598);
    const auto *nsf_599 = buffer.data(nsf + 599);
    const auto *nsf_600 = buffer.data(nsf + 600);
    const auto *nsf_601 = buffer.data(nsf + 601);
    const auto *nsf_602 = buffer.data(nsf + 602);
    const auto *nsf_603 = buffer.data(nsf + 603);
    const auto *nsf_604 = buffer.data(nsf + 604);
    const auto *nsf_605 = buffer.data(nsf + 605);
    const auto *nsf_606 = buffer.data(nsf + 606);
    const auto *nsf_607 = buffer.data(nsf + 607);
    const auto *nsf_608 = buffer.data(nsf + 608);
    const auto *nsf_609 = buffer.data(nsf + 609);
    const auto *nsf_610 = buffer.data(nsf + 610);
    const auto *nsf_611 = buffer.data(nsf + 611);
    const auto *nsf_612 = buffer.data(nsf + 612);
    const auto *nsf_613 = buffer.data(nsf + 613);
    const auto *nsf_614 = buffer.data(nsf + 614);
    const auto *nsf_615 = buffer.data(nsf + 615);
    const auto *nsf_616 = buffer.data(nsf + 616);
    const auto *nsf_617 = buffer.data(nsf + 617);
    const auto *nsf_618 = buffer.data(nsf + 618);
    const auto *nsf_619 = buffer.data(nsf + 619);
    const auto *nsf_620 = buffer.data(nsf + 620);
    const auto *nsf_621 = buffer.data(nsf + 621);
    const auto *nsf_622 = buffer.data(nsf + 622);
    const auto *nsf_623 = buffer.data(nsf + 623);
    const auto *nsf_624 = buffer.data(nsf + 624);
    const auto *nsf_625 = buffer.data(nsf + 625);
    const auto *nsf_626 = buffer.data(nsf + 626);
    const auto *nsf_627 = buffer.data(nsf + 627);
    const auto *nsf_628 = buffer.data(nsf + 628);
    const auto *nsf_629 = buffer.data(nsf + 629);
    const auto *nsf_630 = buffer.data(nsf + 630);
    const auto *nsf_631 = buffer.data(nsf + 631);
    const auto *nsf_632 = buffer.data(nsf + 632);
    const auto *nsf_633 = buffer.data(nsf + 633);
    const auto *nsf_634 = buffer.data(nsf + 634);
    const auto *nsf_635 = buffer.data(nsf + 635);
    const auto *nsf_636 = buffer.data(nsf + 636);
    const auto *nsf_637 = buffer.data(nsf + 637);
    const auto *nsf_638 = buffer.data(nsf + 638);
    const auto *nsf_639 = buffer.data(nsf + 639);
    const auto *nsf_641 = buffer.data(nsf + 641);
    const auto *nsf_643 = buffer.data(nsf + 643);
    const auto *nsf_644 = buffer.data(nsf + 644);
    const auto *nsf_646 = buffer.data(nsf + 646);
    const auto *nsf_647 = buffer.data(nsf + 647);
    const auto *nsf_648 = buffer.data(nsf + 648);
    const auto *nsf_649 = buffer.data(nsf + 649);

#pragma omp simd aligned(t_850, t_851, t_852, t_853, pa_z, pc_y, pc_z, msg0_685, msg0_687, \
                         msf_456, msf_457, msf_469, msg1_685, msg1_687, nsf_566, \
                         nsf_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = pa_z[k] * msg0_685[k]
                   - f_6 * pc_z[k] * msg1_685[k];

        t_851[k] = f_7 * msf_456[k]
                   + f_3 * pc_z[k] * nsf_566[k];

        t_852[k] = pa_z[k] * msg0_687[k]
                   + f_8 * msf_457[k]
                   - f_6 * pc_z[k] * msg1_687[k];

        t_853[k] = f_9 * msf_469[k]
                   + f_3 * pc_y[k] * nsf_569[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, pc_x, pc_z, msf_459, nsd0_341, nsd0_342, \
                         nsd0_343, nsd1_341, nsd1_342, nsd1_343, nsf_569, nsf_570, \
                         nsf_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_7 * msf_459[k]
                   + f_1 * nsd0_341[k]
                   - f_2 * nsd1_341[k]
                   + f_3 * pc_z[k] * nsf_569[k];

        t_855[k] = f_1 * nsd0_342[k]
                   - f_2 * nsd1_342[k]
                   + f_3 * pc_x[k] * nsf_570[k];

        t_856[k] = f_10 * nsd0_343[k]
                   - f_11 * nsd1_343[k]
                   + f_3 * pc_x[k] * nsf_571[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, pc_x, nsd0_344, nsd0_345, nsd0_346, nsd1_344, \
                         nsd1_345, nsd1_346, nsf_572, nsf_573, \
                         nsf_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = f_10 * nsd0_344[k]
                   - f_11 * nsd1_344[k]
                   + f_3 * pc_x[k] * nsf_572[k];

        t_858[k] = f_4 * nsd0_345[k]
                   - f_5 * nsd1_345[k]
                   + f_3 * pc_x[k] * nsf_573[k];

        t_859[k] = f_4 * nsd0_346[k]
                   - f_5 * nsd1_346[k]
                   + f_3 * pc_x[k] * nsf_574[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, t_864, pc_x, nsd0_347, nsd1_347, nsf_575, \
                         nsf_576, nsf_577, nsf_578, nsf_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_4 * nsd0_347[k]
                   - f_5 * nsd1_347[k]
                   + f_3 * pc_x[k] * nsf_575[k];

        t_861[k] = f_3 * pc_x[k] * nsf_576[k];

        t_862[k] = f_3 * pc_x[k] * nsf_577[k];

        t_863[k] = f_3 * pc_x[k] * nsf_578[k];

        t_864[k] = f_3 * pc_x[k] * nsf_579[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, pc_y, pc_z, msf_466, msf_476, msf_478, nsd0_345, \
                         nsd0_347, nsd1_345, nsd1_347, nsf_576, \
                         nsf_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_12 * msf_476[k]
                   + f_1 * nsd0_345[k]
                   - f_2 * nsd1_345[k]
                   + f_3 * pc_y[k] * nsf_576[k];

        t_866[k] = f_8 * msf_466[k]
                   + f_3 * pc_z[k] * nsf_576[k];

        t_867[k] = f_12 * msf_478[k]
                   + f_4 * nsd0_347[k]
                   - f_5 * nsd1_347[k]
                   + f_3 * pc_y[k] * nsf_578[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, pc_x, pc_y, pc_z, msf_469, msf_479, nsd0_347, \
                         nsd0_348, nsd1_347, nsd1_348, nsf_579, \
                         nsf_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = f_12 * msf_479[k]
                   + f_3 * pc_y[k] * nsf_579[k];

        t_869[k] = f_8 * msf_469[k]
                   + f_1 * nsd0_347[k]
                   - f_2 * nsd1_347[k]
                   + f_3 * pc_z[k] * nsf_579[k];

        t_870[k] = f_1 * nsd0_348[k]
                   - f_2 * nsd1_348[k]
                   + f_3 * pc_x[k] * nsf_580[k];
    }

#pragma omp simd aligned(t_871, t_872, t_873, pc_x, nsd0_349, nsd0_350, nsd0_351, nsd1_349, \
                         nsd1_350, nsd1_351, nsf_581, nsf_582, \
                         nsf_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_871[k] = f_10 * nsd0_349[k]
                   - f_11 * nsd1_349[k]
                   + f_3 * pc_x[k] * nsf_581[k];

        t_872[k] = f_10 * nsd0_350[k]
                   - f_11 * nsd1_350[k]
                   + f_3 * pc_x[k] * nsf_582[k];

        t_873[k] = f_4 * nsd0_351[k]
                   - f_5 * nsd1_351[k]
                   + f_3 * pc_x[k] * nsf_583[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, t_878, pc_x, nsd0_352, nsd0_353, \
                         nsd1_352, nsd1_353, nsf_584, nsf_585, nsf_586, nsf_587, \
                         nsf_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_4 * nsd0_352[k]
                   - f_5 * nsd1_352[k]
                   + f_3 * pc_x[k] * nsf_584[k];

        t_875[k] = f_4 * nsd0_353[k]
                   - f_5 * nsd1_353[k]
                   + f_3 * pc_x[k] * nsf_585[k];

        t_876[k] = f_3 * pc_x[k] * nsf_586[k];

        t_877[k] = f_3 * pc_x[k] * nsf_587[k];

        t_878[k] = f_3 * pc_x[k] * nsf_588[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, pc_z, msf_476, msf_486, nsd0_351, \
                         nsd1_351, nsf_586, nsf_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_3 * pc_x[k] * nsf_589[k];

        t_880[k] = f_13 * msf_486[k]
                   + f_1 * nsd0_351[k]
                   - f_2 * nsd1_351[k]
                   + f_3 * pc_y[k] * nsf_586[k];

        t_881[k] = f_14 * msf_476[k]
                   + f_3 * pc_z[k] * nsf_586[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, pc_y, pc_z, msf_479, msf_488, msf_489, nsd0_353, \
                         nsd1_353, nsf_588, nsf_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_13 * msf_488[k]
                   + f_4 * nsd0_353[k]
                   - f_5 * nsd1_353[k]
                   + f_3 * pc_y[k] * nsf_588[k];

        t_883[k] = f_13 * msf_489[k]
                   + f_3 * pc_y[k] * nsf_589[k];

        t_884[k] = f_14 * msf_479[k]
                   + f_1 * nsd0_353[k]
                   - f_2 * nsd1_353[k]
                   + f_3 * pc_z[k] * nsf_589[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, pc_x, nsd0_354, nsd0_355, nsd0_356, nsd1_354, \
                         nsd1_355, nsd1_356, nsf_590, nsf_591, \
                         nsf_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = f_1 * nsd0_354[k]
                   - f_2 * nsd1_354[k]
                   + f_3 * pc_x[k] * nsf_590[k];

        t_886[k] = f_10 * nsd0_355[k]
                   - f_11 * nsd1_355[k]
                   + f_3 * pc_x[k] * nsf_591[k];

        t_887[k] = f_10 * nsd0_356[k]
                   - f_11 * nsd1_356[k]
                   + f_3 * pc_x[k] * nsf_592[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pc_x, nsd0_357, nsd0_358, nsd0_359, \
                         nsd1_357, nsd1_358, nsd1_359, nsf_593, nsf_594, nsf_595, \
                         nsf_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_4 * nsd0_357[k]
                   - f_5 * nsd1_357[k]
                   + f_3 * pc_x[k] * nsf_593[k];

        t_889[k] = f_4 * nsd0_358[k]
                   - f_5 * nsd1_358[k]
                   + f_3 * pc_x[k] * nsf_594[k];

        t_890[k] = f_4 * nsd0_359[k]
                   - f_5 * nsd1_359[k]
                   + f_3 * pc_x[k] * nsf_595[k];

        t_891[k] = f_3 * pc_x[k] * nsf_596[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, t_895, t_896, pc_x, pc_y, pc_z, msf_486, \
                         msf_496, nsd0_357, nsd1_357, nsf_596, nsf_597, nsf_598, \
                         nsf_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_3 * pc_x[k] * nsf_597[k];

        t_893[k] = f_3 * pc_x[k] * nsf_598[k];

        t_894[k] = f_3 * pc_x[k] * nsf_599[k];

        t_895[k] = f_15 * msf_496[k]
                   + f_1 * nsd0_357[k]
                   - f_2 * nsd1_357[k]
                   + f_3 * pc_y[k] * nsf_596[k];

        t_896[k] = f_16 * msf_486[k]
                   + f_3 * pc_z[k] * nsf_596[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, pc_y, pc_z, msf_489, msf_498, msf_499, nsd0_359, \
                         nsd1_359, nsf_598, nsf_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = f_15 * msf_498[k]
                   + f_4 * nsd0_359[k]
                   - f_5 * nsd1_359[k]
                   + f_3 * pc_y[k] * nsf_598[k];

        t_898[k] = f_15 * msf_499[k]
                   + f_3 * pc_y[k] * nsf_599[k];

        t_899[k] = f_16 * msf_489[k]
                   + f_1 * nsd0_359[k]
                   - f_2 * nsd1_359[k]
                   + f_3 * pc_z[k] * nsf_599[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, pc_x, nsd0_360, nsd0_361, nsd0_362, nsd1_360, \
                         nsd1_361, nsd1_362, nsf_600, nsf_601, \
                         nsf_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_1 * nsd0_360[k]
                   - f_2 * nsd1_360[k]
                   + f_3 * pc_x[k] * nsf_600[k];

        t_901[k] = f_10 * nsd0_361[k]
                   - f_11 * nsd1_361[k]
                   + f_3 * pc_x[k] * nsf_601[k];

        t_902[k] = f_10 * nsd0_362[k]
                   - f_11 * nsd1_362[k]
                   + f_3 * pc_x[k] * nsf_602[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, t_906, pc_x, nsd0_363, nsd0_364, nsd0_365, \
                         nsd1_363, nsd1_364, nsd1_365, nsf_603, nsf_604, nsf_605, \
                         nsf_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_4 * nsd0_363[k]
                   - f_5 * nsd1_363[k]
                   + f_3 * pc_x[k] * nsf_603[k];

        t_904[k] = f_4 * nsd0_364[k]
                   - f_5 * nsd1_364[k]
                   + f_3 * pc_x[k] * nsf_604[k];

        t_905[k] = f_4 * nsd0_365[k]
                   - f_5 * nsd1_365[k]
                   + f_3 * pc_x[k] * nsf_605[k];

        t_906[k] = f_3 * pc_x[k] * nsf_606[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, t_910, t_911, pc_x, pc_y, pc_z, msf_496, \
                         msf_506, nsd0_363, nsd1_363, nsf_606, nsf_607, nsf_608, \
                         nsf_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = f_3 * pc_x[k] * nsf_607[k];

        t_908[k] = f_3 * pc_x[k] * nsf_608[k];

        t_909[k] = f_3 * pc_x[k] * nsf_609[k];

        t_910[k] = f_17 * msf_506[k]
                   + f_1 * nsd0_363[k]
                   - f_2 * nsd1_363[k]
                   + f_3 * pc_y[k] * nsf_606[k];

        t_911[k] = f_17 * msf_496[k]
                   + f_3 * pc_z[k] * nsf_606[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pc_y, pc_z, msf_499, msf_508, msf_509, nsd0_365, \
                         nsd1_365, nsf_608, nsf_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_17 * msf_508[k]
                   + f_4 * nsd0_365[k]
                   - f_5 * nsd1_365[k]
                   + f_3 * pc_y[k] * nsf_608[k];

        t_913[k] = f_17 * msf_509[k]
                   + f_3 * pc_y[k] * nsf_609[k];

        t_914[k] = f_17 * msf_499[k]
                   + f_1 * nsd0_365[k]
                   - f_2 * nsd1_365[k]
                   + f_3 * pc_z[k] * nsf_609[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pc_x, nsd0_366, nsd0_367, nsd0_368, nsd1_366, \
                         nsd1_367, nsd1_368, nsf_610, nsf_611, \
                         nsf_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_1 * nsd0_366[k]
                   - f_2 * nsd1_366[k]
                   + f_3 * pc_x[k] * nsf_610[k];

        t_916[k] = f_10 * nsd0_367[k]
                   - f_11 * nsd1_367[k]
                   + f_3 * pc_x[k] * nsf_611[k];

        t_917[k] = f_10 * nsd0_368[k]
                   - f_11 * nsd1_368[k]
                   + f_3 * pc_x[k] * nsf_612[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, t_921, pc_x, nsd0_369, nsd0_370, nsd0_371, \
                         nsd1_369, nsd1_370, nsd1_371, nsf_613, nsf_614, nsf_615, \
                         nsf_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_4 * nsd0_369[k]
                   - f_5 * nsd1_369[k]
                   + f_3 * pc_x[k] * nsf_613[k];

        t_919[k] = f_4 * nsd0_370[k]
                   - f_5 * nsd1_370[k]
                   + f_3 * pc_x[k] * nsf_614[k];

        t_920[k] = f_4 * nsd0_371[k]
                   - f_5 * nsd1_371[k]
                   + f_3 * pc_x[k] * nsf_615[k];

        t_921[k] = f_3 * pc_x[k] * nsf_616[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, t_926, pc_x, pc_y, pc_z, msf_506, \
                         msf_516, nsd0_369, nsd1_369, nsf_616, nsf_617, nsf_618, \
                         nsf_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_3 * pc_x[k] * nsf_617[k];

        t_923[k] = f_3 * pc_x[k] * nsf_618[k];

        t_924[k] = f_3 * pc_x[k] * nsf_619[k];

        t_925[k] = f_16 * msf_516[k]
                   + f_1 * nsd0_369[k]
                   - f_2 * nsd1_369[k]
                   + f_3 * pc_y[k] * nsf_616[k];

        t_926[k] = f_15 * msf_506[k]
                   + f_3 * pc_z[k] * nsf_616[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, pc_y, pc_z, msf_509, msf_518, msf_519, nsd0_371, \
                         nsd1_371, nsf_618, nsf_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_16 * msf_518[k]
                   + f_4 * nsd0_371[k]
                   - f_5 * nsd1_371[k]
                   + f_3 * pc_y[k] * nsf_618[k];

        t_928[k] = f_16 * msf_519[k]
                   + f_3 * pc_y[k] * nsf_619[k];

        t_929[k] = f_15 * msf_509[k]
                   + f_1 * nsd0_371[k]
                   - f_2 * nsd1_371[k]
                   + f_3 * pc_z[k] * nsf_619[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, pc_x, nsd0_372, nsd0_373, nsd0_374, nsd1_372, \
                         nsd1_373, nsd1_374, nsf_620, nsf_621, \
                         nsf_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_1 * nsd0_372[k]
                   - f_2 * nsd1_372[k]
                   + f_3 * pc_x[k] * nsf_620[k];

        t_931[k] = f_10 * nsd0_373[k]
                   - f_11 * nsd1_373[k]
                   + f_3 * pc_x[k] * nsf_621[k];

        t_932[k] = f_10 * nsd0_374[k]
                   - f_11 * nsd1_374[k]
                   + f_3 * pc_x[k] * nsf_622[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, pc_x, nsd0_375, nsd0_376, nsd0_377, \
                         nsd1_375, nsd1_376, nsd1_377, nsf_623, nsf_624, nsf_625, \
                         nsf_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_4 * nsd0_375[k]
                   - f_5 * nsd1_375[k]
                   + f_3 * pc_x[k] * nsf_623[k];

        t_934[k] = f_4 * nsd0_376[k]
                   - f_5 * nsd1_376[k]
                   + f_3 * pc_x[k] * nsf_624[k];

        t_935[k] = f_4 * nsd0_377[k]
                   - f_5 * nsd1_377[k]
                   + f_3 * pc_x[k] * nsf_625[k];

        t_936[k] = f_3 * pc_x[k] * nsf_626[k];
    }

#pragma omp simd aligned(t_937, t_938, t_939, t_940, t_941, pc_x, pc_y, pc_z, msf_516, \
                         msf_526, nsd0_375, nsd1_375, nsf_626, nsf_627, nsf_628, \
                         nsf_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_937[k] = f_3 * pc_x[k] * nsf_627[k];

        t_938[k] = f_3 * pc_x[k] * nsf_628[k];

        t_939[k] = f_3 * pc_x[k] * nsf_629[k];

        t_940[k] = f_14 * msf_526[k]
                   + f_1 * nsd0_375[k]
                   - f_2 * nsd1_375[k]
                   + f_3 * pc_y[k] * nsf_626[k];

        t_941[k] = f_13 * msf_516[k]
                   + f_3 * pc_z[k] * nsf_626[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, pc_y, pc_z, msf_519, msf_528, msf_529, nsd0_377, \
                         nsd1_377, nsf_628, nsf_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_14 * msf_528[k]
                   + f_4 * nsd0_377[k]
                   - f_5 * nsd1_377[k]
                   + f_3 * pc_y[k] * nsf_628[k];

        t_943[k] = f_14 * msf_529[k]
                   + f_3 * pc_y[k] * nsf_629[k];

        t_944[k] = f_13 * msf_519[k]
                   + f_1 * nsd0_377[k]
                   - f_2 * nsd1_377[k]
                   + f_3 * pc_z[k] * nsf_629[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, pc_x, nsd0_378, nsd0_379, nsd0_380, nsd1_378, \
                         nsd1_379, nsd1_380, nsf_630, nsf_631, \
                         nsf_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = f_1 * nsd0_378[k]
                   - f_2 * nsd1_378[k]
                   + f_3 * pc_x[k] * nsf_630[k];

        t_946[k] = f_10 * nsd0_379[k]
                   - f_11 * nsd1_379[k]
                   + f_3 * pc_x[k] * nsf_631[k];

        t_947[k] = f_10 * nsd0_380[k]
                   - f_11 * nsd1_380[k]
                   + f_3 * pc_x[k] * nsf_632[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, t_951, pc_x, nsd0_381, nsd0_382, nsd0_383, \
                         nsd1_381, nsd1_382, nsd1_383, nsf_633, nsf_634, nsf_635, \
                         nsf_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = f_4 * nsd0_381[k]
                   - f_5 * nsd1_381[k]
                   + f_3 * pc_x[k] * nsf_633[k];

        t_949[k] = f_4 * nsd0_382[k]
                   - f_5 * nsd1_382[k]
                   + f_3 * pc_x[k] * nsf_634[k];

        t_950[k] = f_4 * nsd0_383[k]
                   - f_5 * nsd1_383[k]
                   + f_3 * pc_x[k] * nsf_635[k];

        t_951[k] = f_3 * pc_x[k] * nsf_636[k];
    }

#pragma omp simd aligned(t_952, t_953, t_954, t_955, t_956, pc_x, pc_y, pc_z, msf_526, \
                         msf_536, nsd0_381, nsd1_381, nsf_636, nsf_637, nsf_638, \
                         nsf_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_952[k] = f_3 * pc_x[k] * nsf_637[k];

        t_953[k] = f_3 * pc_x[k] * nsf_638[k];

        t_954[k] = f_3 * pc_x[k] * nsf_639[k];

        t_955[k] = f_8 * msf_536[k]
                   + f_1 * nsd0_381[k]
                   - f_2 * nsd1_381[k]
                   + f_3 * pc_y[k] * nsf_636[k];

        t_956[k] = f_12 * msf_526[k]
                   + f_3 * pc_z[k] * nsf_636[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pa_y, pc_y, pc_z, msg0_810, msf_529, \
                         msf_538, msf_539, msg1_810, nsd0_383, nsd1_383, nsf_638, \
                         nsf_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_8 * msf_538[k]
                   + f_4 * nsd0_383[k]
                   - f_5 * nsd1_383[k]
                   + f_3 * pc_y[k] * nsf_638[k];

        t_958[k] = f_8 * msf_539[k]
                   + f_3 * pc_y[k] * nsf_639[k];

        t_959[k] = f_12 * msf_529[k]
                   + f_1 * nsd0_383[k]
                   - f_2 * nsd1_383[k]
                   + f_3 * pc_z[k] * nsf_639[k];

        t_960[k] = pa_y[k] * msg0_810[k]
                   - f_6 * pc_y[k] * msg1_810[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pa_y, pc_x, pc_y, msg0_812, msg1_812, nsd0_385, \
                         nsd0_387, nsd1_385, nsd1_387, nsf_641, \
                         nsf_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = f_10 * nsd0_385[k]
                   - f_11 * nsd1_385[k]
                   + f_3 * pc_x[k] * nsf_641[k];

        t_962[k] = pa_y[k] * msg0_812[k]
                   - f_6 * pc_y[k] * msg1_812[k];

        t_963[k] = f_4 * nsd0_387[k]
                   - f_5 * nsd1_387[k]
                   + f_3 * pc_x[k] * nsf_643[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, t_968, pa_y, pc_x, pc_y, msg0_815, \
                         msg1_815, nsd0_388, nsd1_388, nsf_644, nsf_646, nsf_647, \
                         nsf_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = f_4 * nsd0_388[k]
                   - f_5 * nsd1_388[k]
                   + f_3 * pc_x[k] * nsf_644[k];

        t_965[k] = pa_y[k] * msg0_815[k]
                   - f_6 * pc_y[k] * msg1_815[k];

        t_966[k] = f_3 * pc_x[k] * nsf_646[k];

        t_967[k] = f_3 * pc_x[k] * nsf_647[k];

        t_968[k] = f_3 * pc_x[k] * nsf_648[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pa_y, pc_x, pc_y, pc_z, msg0_820, msf_536, \
                         msf_546, msg1_820, nsf_646, nsf_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_3 * pc_x[k] * nsf_649[k];

        t_970[k] = pa_y[k] * msg0_820[k]
                   + f_16 * msf_546[k]
                   - f_6 * pc_y[k] * msg1_820[k];

        t_971[k] = f_9 * msf_536[k]
                   + f_3 * pc_z[k] * nsf_646[k];
    }
}

static auto
compute_prim_nsg_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msg0,
                                                          const size_t msf, const size_t msg1,
                                                          const size_t nsd0, const size_t nsd1,
                                                          const size_t nsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msg0_822 = buffer.data(msg0 + 822);
    const auto *msg0_824 = buffer.data(msg0 + 824);

    const auto *msf_548 = buffer.data(msf + 548);
    const auto *msf_549 = buffer.data(msf + 549);

    const auto *msg1_822 = buffer.data(msg1 + 822);
    const auto *msg1_824 = buffer.data(msg1 + 824);

    const auto *nsd0_390 = buffer.data(nsd0 + 390);
    const auto *nsd0_392 = buffer.data(nsd0 + 392);
    const auto *nsd0_393 = buffer.data(nsd0 + 393);
    const auto *nsd0_394 = buffer.data(nsd0 + 394);
    const auto *nsd0_395 = buffer.data(nsd0 + 395);

    const auto *nsd1_390 = buffer.data(nsd1 + 390);
    const auto *nsd1_392 = buffer.data(nsd1 + 392);
    const auto *nsd1_393 = buffer.data(nsd1 + 393);
    const auto *nsd1_394 = buffer.data(nsd1 + 394);
    const auto *nsd1_395 = buffer.data(nsd1 + 395);

    const auto *nsf_649 = buffer.data(nsf + 649);
    const auto *nsf_650 = buffer.data(nsf + 650);
    const auto *nsf_652 = buffer.data(nsf + 652);
    const auto *nsf_653 = buffer.data(nsf + 653);
    const auto *nsf_655 = buffer.data(nsf + 655);
    const auto *nsf_656 = buffer.data(nsf + 656);
    const auto *nsf_657 = buffer.data(nsf + 657);
    const auto *nsf_658 = buffer.data(nsf + 658);
    const auto *nsf_659 = buffer.data(nsf + 659);

#pragma omp simd aligned(t_972, t_973, t_974, pa_y, pc_y, msg0_822, msg0_824, msf_548, \
                         msf_549, msg1_822, msg1_824, nsf_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = pa_y[k] * msg0_822[k]
                   + f_8 * msf_548[k]
                   - f_6 * pc_y[k] * msg1_822[k];

        t_973[k] = f_7 * msf_549[k]
                   + f_3 * pc_y[k] * nsf_649[k];

        t_974[k] = pa_y[k] * msg0_824[k]
                   - f_6 * pc_y[k] * msg1_824[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, pc_x, pc_y, nsd0_390, nsd0_392, \
                         nsd0_393, nsd1_390, nsd1_392, nsd1_393, nsf_650, nsf_652, \
                         nsf_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = f_1 * nsd0_390[k]
                   - f_2 * nsd1_390[k]
                   + f_3 * pc_x[k] * nsf_650[k];

        t_976[k] = f_3 * pc_y[k] * nsf_650[k];

        t_977[k] = f_10 * nsd0_392[k]
                   - f_11 * nsd1_392[k]
                   + f_3 * pc_x[k] * nsf_652[k];

        t_978[k] = f_4 * nsd0_393[k]
                   - f_5 * nsd1_393[k]
                   + f_3 * pc_x[k] * nsf_653[k];

        t_979[k] = f_3 * pc_y[k] * nsf_652[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, t_984, pc_x, nsd0_395, nsd1_395, nsf_655, \
                         nsf_656, nsf_657, nsf_658, nsf_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = f_4 * nsd0_395[k]
                   - f_5 * nsd1_395[k]
                   + f_3 * pc_x[k] * nsf_655[k];

        t_981[k] = f_3 * pc_x[k] * nsf_656[k];

        t_982[k] = f_3 * pc_x[k] * nsf_657[k];

        t_983[k] = f_3 * pc_x[k] * nsf_658[k];

        t_984[k] = f_3 * pc_x[k] * nsf_659[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, pc_y, nsd0_393, nsd0_394, nsd0_395, \
                         nsd1_393, nsd1_394, nsd1_395, nsf_656, nsf_657, nsf_658, \
                         nsf_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_1 * nsd0_393[k]
                   - f_2 * nsd1_393[k]
                   + f_3 * pc_y[k] * nsf_656[k];

        t_986[k] = f_10 * nsd0_394[k]
                   - f_11 * nsd1_394[k]
                   + f_3 * pc_y[k] * nsf_657[k];

        t_987[k] = f_4 * nsd0_395[k]
                   - f_5 * nsd1_395[k]
                   + f_3 * pc_y[k] * nsf_658[k];

        t_988[k] = f_3 * pc_y[k] * nsf_659[k];
    }

#pragma omp simd aligned(t_989, pc_z, msf_549, nsd0_395, nsd1_395, \
                         nsf_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_0 * msf_549[k]
                   + f_1 * nsd0_395[k]
                   - f_2 * nsd1_395[k]
                   + f_3 * pc_z[k] * nsf_659[k];
    }
}

auto
compute_prim_nsg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t msg0, const size_t msf,
                                                   const size_t msg1, const size_t nsd0,
                                                   const size_t nsd1, const size_t nsf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_nsg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, msg0, msf,
                                                              msg1, nsd0, nsd1, nsf, ncols,
                                                              gamma, p, q);

    compute_prim_nsg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, msg0, msf,
                                                              msg1, nsd0, nsd1, nsf, ncols,
                                                              gamma, p, q);

    compute_prim_nsg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, msg0, msf,
                                                              msg1, nsd0, nsd1, nsf, ncols,
                                                              gamma, p, q);

    compute_prim_nsg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, msg0, msf,
                                                              msg1, nsd0, nsd1, nsf, ncols,
                                                              gamma, p, q);

    compute_prim_nsg_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, msg0, msf,
                                                              msg1, nsd0, nsd1, nsf, ncols,
                                                              gamma, p, q);

    compute_prim_nsg_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, msg0, msf,
                                                              msg1, nsd0, nsd1, nsf, ncols,
                                                              gamma, p, q);

    compute_prim_nsg_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, msg0, msf,
                                                              msg1, nsd0, nsd1, nsf, ncols,
                                                              gamma, p, q);

    compute_prim_nsg_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, msg0, msf,
                                                              msg1, nsd0, nsd1, nsf, ncols,
                                                              gamma, p, q);

    compute_prim_nsg_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, msg0, msf,
                                                              msg1, nsd0, nsd1, nsf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
