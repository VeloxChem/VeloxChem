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


#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_nsf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msf0,
                                                          const size_t msd, const size_t msf1,
                                                          const size_t nsp0, const size_t nsp1,
                                                          const size_t nsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msf0_0 = buffer.data(msf0 + 0);
    const auto *msf0_6 = buffer.data(msf0 + 6);
    const auto *msf0_9 = buffer.data(msf0 + 9);
    const auto *msf0_16 = buffer.data(msf0 + 16);
    const auto *msf0_20 = buffer.data(msf0 + 20);
    const auto *msf0_29 = buffer.data(msf0 + 29);
    const auto *msf0_30 = buffer.data(msf0 + 30);
    const auto *msf0_36 = buffer.data(msf0 + 36);
    const auto *msf0_50 = buffer.data(msf0 + 50);
    const auto *msf0_59 = buffer.data(msf0 + 59);
    const auto *msf0_60 = buffer.data(msf0 + 60);
    const auto *msf0_66 = buffer.data(msf0 + 66);
    const auto *msf0_90 = buffer.data(msf0 + 90);

    const auto *msd_0 = buffer.data(msd + 0);
    const auto *msd_3 = buffer.data(msd + 3);
    const auto *msd_5 = buffer.data(msd + 5);
    const auto *msd_6 = buffer.data(msd + 6);
    const auto *msd_9 = buffer.data(msd + 9);
    const auto *msd_11 = buffer.data(msd + 11);
    const auto *msd_12 = buffer.data(msd + 12);
    const auto *msd_15 = buffer.data(msd + 15);
    const auto *msd_17 = buffer.data(msd + 17);
    const auto *msd_18 = buffer.data(msd + 18);
    const auto *msd_21 = buffer.data(msd + 21);
    const auto *msd_23 = buffer.data(msd + 23);
    const auto *msd_24 = buffer.data(msd + 24);
    const auto *msd_27 = buffer.data(msd + 27);
    const auto *msd_28 = buffer.data(msd + 28);
    const auto *msd_29 = buffer.data(msd + 29);
    const auto *msd_30 = buffer.data(msd + 30);
    const auto *msd_33 = buffer.data(msd + 33);
    const auto *msd_35 = buffer.data(msd + 35);
    const auto *msd_36 = buffer.data(msd + 36);
    const auto *msd_39 = buffer.data(msd + 39);
    const auto *msd_41 = buffer.data(msd + 41);
    const auto *msd_42 = buffer.data(msd + 42);
    const auto *msd_45 = buffer.data(msd + 45);
    const auto *msd_46 = buffer.data(msd + 46);
    const auto *msd_47 = buffer.data(msd + 47);
    const auto *msd_48 = buffer.data(msd + 48);
    const auto *msd_51 = buffer.data(msd + 51);
    const auto *msd_52 = buffer.data(msd + 52);
    const auto *msd_53 = buffer.data(msd + 53);
    const auto *msd_54 = buffer.data(msd + 54);
    const auto *msd_57 = buffer.data(msd + 57);
    const auto *msd_59 = buffer.data(msd + 59);
    const auto *msd_60 = buffer.data(msd + 60);
    const auto *msd_63 = buffer.data(msd + 63);
    const auto *msd_65 = buffer.data(msd + 65);
    const auto *msd_69 = buffer.data(msd + 69);
    const auto *msd_70 = buffer.data(msd + 70);
    const auto *msd_71 = buffer.data(msd + 71);
    const auto *msd_72 = buffer.data(msd + 72);
    const auto *msd_75 = buffer.data(msd + 75);
    const auto *msd_76 = buffer.data(msd + 76);
    const auto *msd_77 = buffer.data(msd + 77);

    const auto *msf1_0 = buffer.data(msf1 + 0);
    const auto *msf1_6 = buffer.data(msf1 + 6);
    const auto *msf1_9 = buffer.data(msf1 + 9);
    const auto *msf1_16 = buffer.data(msf1 + 16);
    const auto *msf1_20 = buffer.data(msf1 + 20);
    const auto *msf1_29 = buffer.data(msf1 + 29);
    const auto *msf1_30 = buffer.data(msf1 + 30);
    const auto *msf1_36 = buffer.data(msf1 + 36);
    const auto *msf1_50 = buffer.data(msf1 + 50);
    const auto *msf1_59 = buffer.data(msf1 + 59);
    const auto *msf1_60 = buffer.data(msf1 + 60);
    const auto *msf1_66 = buffer.data(msf1 + 66);
    const auto *msf1_90 = buffer.data(msf1 + 90);

    const auto *nsp0_0 = buffer.data(nsp0 + 0);
    const auto *nsp0_1 = buffer.data(nsp0 + 1);
    const auto *nsp0_2 = buffer.data(nsp0 + 2);
    const auto *nsp0_4 = buffer.data(nsp0 + 4);
    const auto *nsp0_8 = buffer.data(nsp0 + 8);
    const auto *nsp0_9 = buffer.data(nsp0 + 9);
    const auto *nsp0_10 = buffer.data(nsp0 + 10);
    const auto *nsp0_11 = buffer.data(nsp0 + 11);
    const auto *nsp0_15 = buffer.data(nsp0 + 15);
    const auto *nsp0_16 = buffer.data(nsp0 + 16);
    const auto *nsp0_17 = buffer.data(nsp0 + 17);
    const auto *nsp0_18 = buffer.data(nsp0 + 18);
    const auto *nsp0_19 = buffer.data(nsp0 + 19);
    const auto *nsp0_20 = buffer.data(nsp0 + 20);
    const auto *nsp0_23 = buffer.data(nsp0 + 23);
    const auto *nsp0_25 = buffer.data(nsp0 + 25);
    const auto *nsp0_27 = buffer.data(nsp0 + 27);
    const auto *nsp0_28 = buffer.data(nsp0 + 28);
    const auto *nsp0_29 = buffer.data(nsp0 + 29);
    const auto *nsp0_30 = buffer.data(nsp0 + 30);
    const auto *nsp0_31 = buffer.data(nsp0 + 31);
    const auto *nsp0_32 = buffer.data(nsp0 + 32);
    const auto *nsp0_35 = buffer.data(nsp0 + 35);
    const auto *nsp0_36 = buffer.data(nsp0 + 36);
    const auto *nsp0_37 = buffer.data(nsp0 + 37);
    const auto *nsp0_38 = buffer.data(nsp0 + 38);

    const auto *nsp1_0 = buffer.data(nsp1 + 0);
    const auto *nsp1_1 = buffer.data(nsp1 + 1);
    const auto *nsp1_2 = buffer.data(nsp1 + 2);
    const auto *nsp1_4 = buffer.data(nsp1 + 4);
    const auto *nsp1_8 = buffer.data(nsp1 + 8);
    const auto *nsp1_9 = buffer.data(nsp1 + 9);
    const auto *nsp1_10 = buffer.data(nsp1 + 10);
    const auto *nsp1_11 = buffer.data(nsp1 + 11);
    const auto *nsp1_15 = buffer.data(nsp1 + 15);
    const auto *nsp1_16 = buffer.data(nsp1 + 16);
    const auto *nsp1_17 = buffer.data(nsp1 + 17);
    const auto *nsp1_18 = buffer.data(nsp1 + 18);
    const auto *nsp1_19 = buffer.data(nsp1 + 19);
    const auto *nsp1_20 = buffer.data(nsp1 + 20);
    const auto *nsp1_23 = buffer.data(nsp1 + 23);
    const auto *nsp1_25 = buffer.data(nsp1 + 25);
    const auto *nsp1_27 = buffer.data(nsp1 + 27);
    const auto *nsp1_28 = buffer.data(nsp1 + 28);
    const auto *nsp1_29 = buffer.data(nsp1 + 29);
    const auto *nsp1_30 = buffer.data(nsp1 + 30);
    const auto *nsp1_31 = buffer.data(nsp1 + 31);
    const auto *nsp1_32 = buffer.data(nsp1 + 32);
    const auto *nsp1_35 = buffer.data(nsp1 + 35);
    const auto *nsp1_36 = buffer.data(nsp1 + 36);
    const auto *nsp1_37 = buffer.data(nsp1 + 37);
    const auto *nsp1_38 = buffer.data(nsp1 + 38);

    const auto *nsd_0 = buffer.data(nsd + 0);
    const auto *nsd_2 = buffer.data(nsd + 2);
    const auto *nsd_3 = buffer.data(nsd + 3);
    const auto *nsd_5 = buffer.data(nsd + 5);
    const auto *nsd_6 = buffer.data(nsd + 6);
    const auto *nsd_7 = buffer.data(nsd + 7);
    const auto *nsd_9 = buffer.data(nsd + 9);
    const auto *nsd_11 = buffer.data(nsd + 11);
    const auto *nsd_12 = buffer.data(nsd + 12);
    const auto *nsd_14 = buffer.data(nsd + 14);
    const auto *nsd_15 = buffer.data(nsd + 15);
    const auto *nsd_16 = buffer.data(nsd + 16);
    const auto *nsd_17 = buffer.data(nsd + 17);
    const auto *nsd_18 = buffer.data(nsd + 18);
    const auto *nsd_19 = buffer.data(nsd + 19);
    const auto *nsd_21 = buffer.data(nsd + 21);
    const auto *nsd_23 = buffer.data(nsd + 23);
    const auto *nsd_24 = buffer.data(nsd + 24);
    const auto *nsd_27 = buffer.data(nsd + 27);
    const auto *nsd_28 = buffer.data(nsd + 28);
    const auto *nsd_29 = buffer.data(nsd + 29);
    const auto *nsd_30 = buffer.data(nsd + 30);
    const auto *nsd_32 = buffer.data(nsd + 32);
    const auto *nsd_33 = buffer.data(nsd + 33);
    const auto *nsd_34 = buffer.data(nsd + 34);
    const auto *nsd_35 = buffer.data(nsd + 35);
    const auto *nsd_36 = buffer.data(nsd + 36);
    const auto *nsd_37 = buffer.data(nsd + 37);
    const auto *nsd_39 = buffer.data(nsd + 39);
    const auto *nsd_41 = buffer.data(nsd + 41);
    const auto *nsd_42 = buffer.data(nsd + 42);
    const auto *nsd_45 = buffer.data(nsd + 45);
    const auto *nsd_46 = buffer.data(nsd + 46);
    const auto *nsd_47 = buffer.data(nsd + 47);
    const auto *nsd_48 = buffer.data(nsd + 48);
    const auto *nsd_51 = buffer.data(nsd + 51);
    const auto *nsd_52 = buffer.data(nsd + 52);
    const auto *nsd_53 = buffer.data(nsd + 53);
    const auto *nsd_54 = buffer.data(nsd + 54);
    const auto *nsd_56 = buffer.data(nsd + 56);
    const auto *nsd_57 = buffer.data(nsd + 57);
    const auto *nsd_58 = buffer.data(nsd + 58);
    const auto *nsd_59 = buffer.data(nsd + 59);
    const auto *nsd_60 = buffer.data(nsd + 60);
    const auto *nsd_61 = buffer.data(nsd + 61);
    const auto *nsd_63 = buffer.data(nsd + 63);
    const auto *nsd_65 = buffer.data(nsd + 65);
    const auto *nsd_66 = buffer.data(nsd + 66);
    const auto *nsd_69 = buffer.data(nsd + 69);
    const auto *nsd_70 = buffer.data(nsd + 70);
    const auto *nsd_71 = buffer.data(nsd + 71);
    const auto *nsd_72 = buffer.data(nsd + 72);
    const auto *nsd_75 = buffer.data(nsd + 75);
    const auto *nsd_76 = buffer.data(nsd + 76);
    const auto *nsd_77 = buffer.data(nsd + 77);
    const auto *nsd_78 = buffer.data(nsd + 78);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, msd_0, msd_3, nsp0_0, \
                         nsp1_0, nsd_0, nsd_2, nsd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * msd_0[k]
                 + f_1 * nsp0_0[k]
                 - f_2 * nsp1_0[k]
                 + f_3 * pc_x[k] * nsd_0[k];

        t_1[k] = f_3 * pc_y[k] * nsd_0[k];

        t_2[k] = f_3 * pc_z[k] * nsd_0[k];

        t_3[k] = f_0 * msd_3[k]
                 + f_3 * pc_x[k] * nsd_3[k];

        t_4[k] = f_3 * pc_y[k] * nsd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, msd_5, nsp0_1, nsp0_2, \
                         nsp1_1, nsp1_2, nsd_3, nsd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * msd_5[k]
                 + f_3 * pc_x[k] * nsd_5[k];

        t_6[k] = f_1 * nsp0_1[k]
                 - f_2 * nsp1_1[k]
                 + f_3 * pc_y[k] * nsd_3[k];

        t_7[k] = f_3 * pc_z[k] * nsd_3[k];

        t_8[k] = f_3 * pc_y[k] * nsd_5[k];

        t_9[k] = f_1 * nsp0_2[k]
                 - f_2 * nsp1_2[k]
                 + f_3 * pc_z[k] * nsd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, msf0_0, msd_0, \
                         msd_9, msf1_0, nsd_6, nsd_7, nsd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * msf0_0[k]
                  - f_4 * pc_y[k] * msf1_0[k];

        t_11[k] = f_5 * msd_0[k]
                  + f_3 * pc_y[k] * nsd_6[k];

        t_12[k] = f_3 * pc_z[k] * nsd_6[k];

        t_13[k] = f_6 * msd_9[k]
                  + f_3 * pc_x[k] * nsd_9[k];

        t_14[k] = f_3 * pc_z[k] * nsd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, msd_3, msd_5, msd_11, \
                         nsp0_4, nsp1_4, nsd_9, nsd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * msd_11[k]
                  + f_3 * pc_x[k] * nsd_11[k];

        t_16[k] = f_5 * msd_3[k]
                  + f_1 * nsp0_4[k]
                  - f_2 * nsp1_4[k]
                  + f_3 * pc_y[k] * nsd_9[k];

        t_17[k] = f_3 * pc_z[k] * nsd_9[k];

        t_18[k] = f_5 * msd_5[k]
                  + f_3 * pc_y[k] * nsd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, msf0_0, msf0_9, \
                         msd_0, msf1_0, msf1_9, nsd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * msf0_9[k]
                  - f_4 * pc_y[k] * msf1_9[k];

        t_20[k] = pa_z[k] * msf0_0[k]
                  - f_4 * pc_z[k] * msf1_0[k];

        t_21[k] = f_3 * pc_y[k] * nsd_12[k];

        t_22[k] = f_5 * msd_0[k]
                  + f_3 * pc_z[k] * nsd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, msf0_6, msd_15, \
                         msd_17, msf1_6, nsd_14, nsd_15, nsd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * msd_15[k]
                  + f_3 * pc_x[k] * nsd_15[k];

        t_24[k] = f_3 * pc_y[k] * nsd_14[k];

        t_25[k] = f_6 * msd_17[k]
                  + f_3 * pc_x[k] * nsd_17[k];

        t_26[k] = pa_z[k] * msf0_6[k]
                  - f_4 * pc_z[k] * msf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, msd_5, msd_18, nsp0_8, \
                         nsp0_9, nsp1_8, nsp1_9, nsd_16, nsd_17, \
                         nsd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * nsp0_8[k]
                  - f_8 * nsp1_8[k]
                  + f_3 * pc_y[k] * nsd_16[k];

        t_28[k] = f_3 * pc_y[k] * nsd_17[k];

        t_29[k] = f_5 * msd_5[k]
                  + f_1 * nsp0_8[k]
                  - f_2 * nsp1_8[k]
                  + f_3 * pc_z[k] * nsd_17[k];

        t_30[k] = f_9 * msd_18[k]
                  + f_1 * nsp0_9[k]
                  - f_2 * nsp1_9[k]
                  + f_3 * pc_x[k] * nsd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, msd_6, msd_21, \
                         msd_23, nsd_18, nsd_19, nsd_21, nsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * msd_6[k]
                  + f_3 * pc_y[k] * nsd_18[k];

        t_32[k] = f_3 * pc_z[k] * nsd_18[k];

        t_33[k] = f_9 * msd_21[k]
                  + f_3 * pc_x[k] * nsd_21[k];

        t_34[k] = f_3 * pc_z[k] * nsd_19[k];

        t_35[k] = f_9 * msd_23[k]
                  + f_3 * pc_x[k] * nsd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, msd_9, msd_11, nsp0_10, nsp0_11, \
                         nsp1_10, nsp1_11, nsd_21, nsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_10 * msd_9[k]
                  + f_1 * nsp0_10[k]
                  - f_2 * nsp1_10[k]
                  + f_3 * pc_y[k] * nsd_21[k];

        t_37[k] = f_3 * pc_z[k] * nsd_21[k];

        t_38[k] = f_10 * msd_11[k]
                  + f_3 * pc_y[k] * nsd_23[k];

        t_39[k] = f_1 * nsp0_11[k]
                  - f_2 * nsp1_11[k]
                  + f_3 * pc_z[k] * nsd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, msf0_20, msd_6, \
                         msd_12, msd_27, msf1_20, nsd_24, nsd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * msf0_20[k]
                  - f_4 * pc_y[k] * msf1_20[k];

        t_41[k] = f_5 * msd_12[k]
                  + f_3 * pc_y[k] * nsd_24[k];

        t_42[k] = f_5 * msd_6[k]
                  + f_3 * pc_z[k] * nsd_24[k];

        t_43[k] = f_9 * msd_27[k]
                  + f_3 * pc_x[k] * nsd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_x, pc_z, msf0_16, msd_9, msd_28, \
                         msd_29, msf1_16, nsd_27, nsd_28, nsd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * msd_28[k]
                  + f_3 * pc_x[k] * nsd_28[k];

        t_45[k] = f_9 * msd_29[k]
                  + f_3 * pc_x[k] * nsd_29[k];

        t_46[k] = pa_z[k] * msf0_16[k]
                  - f_4 * pc_z[k] * msf1_16[k];

        t_47[k] = f_5 * msd_9[k]
                  + f_3 * pc_z[k] * nsd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pc_x, pc_y, msf0_29, msd_17, msd_30, \
                         msf1_29, nsp0_15, nsp1_15, nsd_29, nsd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * msd_17[k]
                  + f_3 * pc_y[k] * nsd_29[k];

        t_49[k] = pa_y[k] * msf0_29[k]
                  - f_4 * pc_y[k] * msf1_29[k];

        t_50[k] = f_9 * msd_30[k]
                  + f_1 * nsp0_15[k]
                  - f_2 * nsp1_15[k]
                  + f_3 * pc_x[k] * nsd_30[k];

        t_51[k] = f_3 * pc_y[k] * nsd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, msd_12, msd_33, msd_35, \
                         nsd_30, nsd_32, nsd_33, nsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * msd_12[k]
                  + f_3 * pc_z[k] * nsd_30[k];

        t_53[k] = f_9 * msd_33[k]
                  + f_3 * pc_x[k] * nsd_33[k];

        t_54[k] = f_3 * pc_y[k] * nsd_32[k];

        t_55[k] = f_9 * msd_35[k]
                  + f_3 * pc_x[k] * nsd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_y, pc_z, msd_17, nsp0_16, nsp0_17, \
                         nsp1_16, nsp1_17, nsd_33, nsd_34, nsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * nsp0_16[k]
                  - f_2 * nsp1_16[k]
                  + f_3 * pc_y[k] * nsd_33[k];

        t_57[k] = f_7 * nsp0_17[k]
                  - f_8 * nsp1_17[k]
                  + f_3 * pc_y[k] * nsd_34[k];

        t_58[k] = f_3 * pc_y[k] * nsd_35[k];

        t_59[k] = f_10 * msd_17[k]
                  + f_1 * nsp0_17[k]
                  - f_2 * nsp1_17[k]
                  + f_3 * pc_z[k] * nsd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pc_x, pc_y, pc_z, msd_18, msd_36, \
                         msd_39, nsp0_18, nsp1_18, nsd_36, nsd_37, \
                         nsd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * msd_36[k]
                  + f_1 * nsp0_18[k]
                  - f_2 * nsp1_18[k]
                  + f_3 * pc_x[k] * nsd_36[k];

        t_61[k] = f_12 * msd_18[k]
                  + f_3 * pc_y[k] * nsd_36[k];

        t_62[k] = f_3 * pc_z[k] * nsd_36[k];

        t_63[k] = f_11 * msd_39[k]
                  + f_3 * pc_x[k] * nsd_39[k];

        t_64[k] = f_3 * pc_z[k] * nsd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_x, pc_y, pc_z, msd_21, msd_23, msd_41, \
                         nsp0_19, nsp1_19, nsd_39, nsd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * msd_41[k]
                  + f_3 * pc_x[k] * nsd_41[k];

        t_66[k] = f_12 * msd_21[k]
                  + f_1 * nsp0_19[k]
                  - f_2 * nsp1_19[k]
                  + f_3 * pc_y[k] * nsd_39[k];

        t_67[k] = f_3 * pc_z[k] * nsd_39[k];

        t_68[k] = f_12 * msd_23[k]
                  + f_3 * pc_y[k] * nsd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, msf0_30, msd_18, msd_24, \
                         msf1_30, nsp0_20, nsp1_20, nsd_41, nsd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * nsp0_20[k]
                  - f_2 * nsp1_20[k]
                  + f_3 * pc_z[k] * nsd_41[k];

        t_70[k] = pa_z[k] * msf0_30[k]
                  - f_4 * pc_z[k] * msf1_30[k];

        t_71[k] = f_10 * msd_24[k]
                  + f_3 * pc_y[k] * nsd_42[k];

        t_72[k] = f_5 * msd_18[k]
                  + f_3 * pc_z[k] * nsd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_z, msf0_36, msd_45, msd_46, \
                         msd_47, msf1_36, nsd_45, nsd_46, nsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_11 * msd_45[k]
                  + f_3 * pc_x[k] * nsd_45[k];

        t_74[k] = f_11 * msd_46[k]
                  + f_3 * pc_x[k] * nsd_46[k];

        t_75[k] = f_11 * msd_47[k]
                  + f_3 * pc_x[k] * nsd_47[k];

        t_76[k] = pa_z[k] * msf0_36[k]
                  - f_4 * pc_z[k] * msf1_36[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, msf0_50, msd_21, msd_23, \
                         msd_29, msf1_50, nsp0_23, nsp1_23, nsd_45, \
                         nsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * msd_21[k]
                  + f_3 * pc_z[k] * nsd_45[k];

        t_78[k] = f_10 * msd_29[k]
                  + f_3 * pc_y[k] * nsd_47[k];

        t_79[k] = f_5 * msd_23[k]
                  + f_1 * nsp0_23[k]
                  - f_2 * nsp1_23[k]
                  + f_3 * pc_z[k] * nsd_47[k];

        t_80[k] = pa_y[k] * msf0_50[k]
                  - f_4 * pc_y[k] * msf1_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, msd_24, msd_30, msd_51, \
                         msd_52, nsd_48, nsd_51, nsd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * msd_30[k]
                  + f_3 * pc_y[k] * nsd_48[k];

        t_82[k] = f_10 * msd_24[k]
                  + f_3 * pc_z[k] * nsd_48[k];

        t_83[k] = f_11 * msd_51[k]
                  + f_3 * pc_x[k] * nsd_51[k];

        t_84[k] = f_11 * msd_52[k]
                  + f_3 * pc_x[k] * nsd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, msd_27, msd_33, msd_35, \
                         msd_53, nsp0_25, nsp1_25, nsd_51, nsd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_11 * msd_53[k]
                  + f_3 * pc_x[k] * nsd_53[k];

        t_86[k] = f_5 * msd_33[k]
                  + f_1 * nsp0_25[k]
                  - f_2 * nsp1_25[k]
                  + f_3 * pc_y[k] * nsd_51[k];

        t_87[k] = f_10 * msd_27[k]
                  + f_3 * pc_z[k] * nsd_51[k];

        t_88[k] = f_5 * msd_35[k]
                  + f_3 * pc_y[k] * nsd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_x, pc_y, pc_z, msf0_59, msd_30, \
                         msd_54, msf1_59, nsp0_27, nsp1_27, nsd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * msf0_59[k]
                  - f_4 * pc_y[k] * msf1_59[k];

        t_90[k] = f_11 * msd_54[k]
                  + f_1 * nsp0_27[k]
                  - f_2 * nsp1_27[k]
                  + f_3 * pc_x[k] * nsd_54[k];

        t_91[k] = f_3 * pc_y[k] * nsd_54[k];

        t_92[k] = f_12 * msd_30[k]
                  + f_3 * pc_z[k] * nsd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_y, msd_57, msd_59, nsp0_28, nsp1_28, \
                         nsd_56, nsd_57, nsd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_11 * msd_57[k]
                  + f_3 * pc_x[k] * nsd_57[k];

        t_94[k] = f_3 * pc_y[k] * nsd_56[k];

        t_95[k] = f_11 * msd_59[k]
                  + f_3 * pc_x[k] * nsd_59[k];

        t_96[k] = f_1 * nsp0_28[k]
                  - f_2 * nsp1_28[k]
                  + f_3 * pc_y[k] * nsd_57[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, pc_z, msd_35, msd_60, nsp0_29, \
                         nsp0_30, nsp1_29, nsp1_30, nsd_58, nsd_59, \
                         nsd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_7 * nsp0_29[k]
                  - f_8 * nsp1_29[k]
                  + f_3 * pc_y[k] * nsd_58[k];

        t_98[k] = f_3 * pc_y[k] * nsd_59[k];

        t_99[k] = f_12 * msd_35[k]
                  + f_1 * nsp0_29[k]
                  - f_2 * nsp1_29[k]
                  + f_3 * pc_z[k] * nsd_59[k];

        t_100[k] = f_13 * msd_60[k]
                   + f_1 * nsp0_30[k]
                   - f_2 * nsp1_30[k]
                   + f_3 * pc_x[k] * nsd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pc_x, pc_y, pc_z, msd_36, msd_63, \
                         msd_65, nsd_60, nsd_61, nsd_63, nsd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_14 * msd_36[k]
                   + f_3 * pc_y[k] * nsd_60[k];

        t_102[k] = f_3 * pc_z[k] * nsd_60[k];

        t_103[k] = f_13 * msd_63[k]
                   + f_3 * pc_x[k] * nsd_63[k];

        t_104[k] = f_3 * pc_z[k] * nsd_61[k];

        t_105[k] = f_13 * msd_65[k]
                   + f_3 * pc_x[k] * nsd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_y, pc_z, msd_39, msd_41, nsp0_31, \
                         nsp0_32, nsp1_31, nsp1_32, nsd_63, nsd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_14 * msd_39[k]
                   + f_1 * nsp0_31[k]
                   - f_2 * nsp1_31[k]
                   + f_3 * pc_y[k] * nsd_63[k];

        t_107[k] = f_3 * pc_z[k] * nsd_63[k];

        t_108[k] = f_14 * msd_41[k]
                   + f_3 * pc_y[k] * nsd_65[k];

        t_109[k] = f_1 * nsp0_32[k]
                   - f_2 * nsp1_32[k]
                   + f_3 * pc_z[k] * nsd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pc_x, pc_y, pc_z, msf0_60, msd_36, \
                         msd_42, msd_69, msf1_60, nsd_66, nsd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * msf0_60[k]
                   - f_4 * pc_z[k] * msf1_60[k];

        t_111[k] = f_12 * msd_42[k]
                   + f_3 * pc_y[k] * nsd_66[k];

        t_112[k] = f_5 * msd_36[k]
                   + f_3 * pc_z[k] * nsd_66[k];

        t_113[k] = f_13 * msd_69[k]
                   + f_3 * pc_x[k] * nsd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pc_x, pc_z, msf0_66, msd_39, \
                         msd_70, msd_71, msf1_66, nsd_69, nsd_70, \
                         nsd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * msd_70[k]
                   + f_3 * pc_x[k] * nsd_70[k];

        t_115[k] = f_13 * msd_71[k]
                   + f_3 * pc_x[k] * nsd_71[k];

        t_116[k] = pa_z[k] * msf0_66[k]
                   - f_4 * pc_z[k] * msf1_66[k];

        t_117[k] = f_5 * msd_39[k]
                   + f_3 * pc_z[k] * nsd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_y, pc_z, msd_41, msd_47, msd_72, \
                         nsp0_35, nsp0_36, nsp1_35, nsp1_36, nsd_71, \
                         nsd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_12 * msd_47[k]
                   + f_3 * pc_y[k] * nsd_71[k];

        t_119[k] = f_5 * msd_41[k]
                   + f_1 * nsp0_35[k]
                   - f_2 * nsp1_35[k]
                   + f_3 * pc_z[k] * nsd_71[k];

        t_120[k] = f_13 * msd_72[k]
                   + f_1 * nsp0_36[k]
                   - f_2 * nsp1_36[k]
                   + f_3 * pc_x[k] * nsd_72[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, msd_42, msd_48, msd_75, \
                         msd_76, nsd_72, nsd_75, nsd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * msd_48[k]
                   + f_3 * pc_y[k] * nsd_72[k];

        t_122[k] = f_10 * msd_42[k]
                   + f_3 * pc_z[k] * nsd_72[k];

        t_123[k] = f_13 * msd_75[k]
                   + f_3 * pc_x[k] * nsd_75[k];

        t_124[k] = f_13 * msd_76[k]
                   + f_3 * pc_x[k] * nsd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_y, pc_z, msd_45, msd_51, msd_53, \
                         msd_77, nsp0_37, nsp1_37, nsd_75, nsd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_13 * msd_77[k]
                   + f_3 * pc_x[k] * nsd_77[k];

        t_126[k] = f_10 * msd_51[k]
                   + f_1 * nsp0_37[k]
                   - f_2 * nsp1_37[k]
                   + f_3 * pc_y[k] * nsd_75[k];

        t_127[k] = f_10 * msd_45[k]
                   + f_3 * pc_z[k] * nsd_75[k];

        t_128[k] = f_10 * msd_53[k]
                   + f_3 * pc_y[k] * nsd_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pc_y, pc_z, msf0_90, msd_47, \
                         msd_48, msd_54, msf1_90, nsp0_38, nsp1_38, nsd_77, \
                         nsd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * msd_47[k]
                   + f_1 * nsp0_38[k]
                   - f_2 * nsp1_38[k]
                   + f_3 * pc_z[k] * nsd_77[k];

        t_130[k] = pa_y[k] * msf0_90[k]
                   - f_4 * pc_y[k] * msf1_90[k];

        t_131[k] = f_5 * msd_54[k]
                   + f_3 * pc_y[k] * nsd_78[k];

        t_132[k] = f_12 * msd_48[k]
                   + f_3 * pc_z[k] * nsd_78[k];
    }
}

static auto
compute_prim_nsf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msf0,
                                                          const size_t msd, const size_t msf1,
                                                          const size_t nsp0, const size_t nsp1,
                                                          const size_t nsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_10 = 1.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msf0_99 = buffer.data(msf0 + 99);
    const auto *msf0_100 = buffer.data(msf0 + 100);
    const auto *msf0_106 = buffer.data(msf0 + 106);
    const auto *msf0_140 = buffer.data(msf0 + 140);
    const auto *msf0_149 = buffer.data(msf0 + 149);
    const auto *msf0_150 = buffer.data(msf0 + 150);
    const auto *msf0_156 = buffer.data(msf0 + 156);

    const auto *msd_51 = buffer.data(msd + 51);
    const auto *msd_54 = buffer.data(msd + 54);
    const auto *msd_57 = buffer.data(msd + 57);
    const auto *msd_59 = buffer.data(msd + 59);
    const auto *msd_60 = buffer.data(msd + 60);
    const auto *msd_63 = buffer.data(msd + 63);
    const auto *msd_65 = buffer.data(msd + 65);
    const auto *msd_66 = buffer.data(msd + 66);
    const auto *msd_69 = buffer.data(msd + 69);
    const auto *msd_71 = buffer.data(msd + 71);
    const auto *msd_72 = buffer.data(msd + 72);
    const auto *msd_75 = buffer.data(msd + 75);
    const auto *msd_77 = buffer.data(msd + 77);
    const auto *msd_78 = buffer.data(msd + 78);
    const auto *msd_81 = buffer.data(msd + 81);
    const auto *msd_82 = buffer.data(msd + 82);
    const auto *msd_83 = buffer.data(msd + 83);
    const auto *msd_84 = buffer.data(msd + 84);
    const auto *msd_87 = buffer.data(msd + 87);
    const auto *msd_89 = buffer.data(msd + 89);
    const auto *msd_90 = buffer.data(msd + 90);
    const auto *msd_93 = buffer.data(msd + 93);
    const auto *msd_95 = buffer.data(msd + 95);
    const auto *msd_96 = buffer.data(msd + 96);
    const auto *msd_99 = buffer.data(msd + 99);
    const auto *msd_100 = buffer.data(msd + 100);
    const auto *msd_101 = buffer.data(msd + 101);
    const auto *msd_102 = buffer.data(msd + 102);
    const auto *msd_105 = buffer.data(msd + 105);
    const auto *msd_106 = buffer.data(msd + 106);
    const auto *msd_107 = buffer.data(msd + 107);
    const auto *msd_108 = buffer.data(msd + 108);
    const auto *msd_111 = buffer.data(msd + 111);
    const auto *msd_112 = buffer.data(msd + 112);
    const auto *msd_113 = buffer.data(msd + 113);
    const auto *msd_114 = buffer.data(msd + 114);
    const auto *msd_117 = buffer.data(msd + 117);
    const auto *msd_118 = buffer.data(msd + 118);
    const auto *msd_119 = buffer.data(msd + 119);
    const auto *msd_120 = buffer.data(msd + 120);
    const auto *msd_123 = buffer.data(msd + 123);
    const auto *msd_125 = buffer.data(msd + 125);
    const auto *msd_126 = buffer.data(msd + 126);
    const auto *msd_129 = buffer.data(msd + 129);
    const auto *msd_131 = buffer.data(msd + 131);
    const auto *msd_135 = buffer.data(msd + 135);
    const auto *msd_136 = buffer.data(msd + 136);
    const auto *msd_137 = buffer.data(msd + 137);
    const auto *msd_138 = buffer.data(msd + 138);
    const auto *msd_141 = buffer.data(msd + 141);
    const auto *msd_142 = buffer.data(msd + 142);
    const auto *msd_143 = buffer.data(msd + 143);
    const auto *msd_144 = buffer.data(msd + 144);
    const auto *msd_147 = buffer.data(msd + 147);
    const auto *msd_148 = buffer.data(msd + 148);
    const auto *msd_149 = buffer.data(msd + 149);
    const auto *msd_150 = buffer.data(msd + 150);
    const auto *msd_153 = buffer.data(msd + 153);
    const auto *msd_154 = buffer.data(msd + 154);
    const auto *msd_155 = buffer.data(msd + 155);

    const auto *msf1_99 = buffer.data(msf1 + 99);
    const auto *msf1_100 = buffer.data(msf1 + 100);
    const auto *msf1_106 = buffer.data(msf1 + 106);
    const auto *msf1_140 = buffer.data(msf1 + 140);
    const auto *msf1_149 = buffer.data(msf1 + 149);
    const auto *msf1_150 = buffer.data(msf1 + 150);
    const auto *msf1_156 = buffer.data(msf1 + 156);

    const auto *nsp0_40 = buffer.data(nsp0 + 40);
    const auto *nsp0_42 = buffer.data(nsp0 + 42);
    const auto *nsp0_43 = buffer.data(nsp0 + 43);
    const auto *nsp0_44 = buffer.data(nsp0 + 44);
    const auto *nsp0_45 = buffer.data(nsp0 + 45);
    const auto *nsp0_46 = buffer.data(nsp0 + 46);
    const auto *nsp0_47 = buffer.data(nsp0 + 47);
    const auto *nsp0_50 = buffer.data(nsp0 + 50);
    const auto *nsp0_51 = buffer.data(nsp0 + 51);
    const auto *nsp0_52 = buffer.data(nsp0 + 52);
    const auto *nsp0_53 = buffer.data(nsp0 + 53);
    const auto *nsp0_54 = buffer.data(nsp0 + 54);
    const auto *nsp0_55 = buffer.data(nsp0 + 55);
    const auto *nsp0_56 = buffer.data(nsp0 + 56);
    const auto *nsp0_58 = buffer.data(nsp0 + 58);
    const auto *nsp0_60 = buffer.data(nsp0 + 60);
    const auto *nsp0_61 = buffer.data(nsp0 + 61);
    const auto *nsp0_62 = buffer.data(nsp0 + 62);
    const auto *nsp0_63 = buffer.data(nsp0 + 63);
    const auto *nsp0_64 = buffer.data(nsp0 + 64);
    const auto *nsp0_65 = buffer.data(nsp0 + 65);
    const auto *nsp0_68 = buffer.data(nsp0 + 68);
    const auto *nsp0_69 = buffer.data(nsp0 + 69);
    const auto *nsp0_70 = buffer.data(nsp0 + 70);
    const auto *nsp0_71 = buffer.data(nsp0 + 71);
    const auto *nsp0_72 = buffer.data(nsp0 + 72);
    const auto *nsp0_73 = buffer.data(nsp0 + 73);
    const auto *nsp0_74 = buffer.data(nsp0 + 74);
    const auto *nsp0_75 = buffer.data(nsp0 + 75);

    const auto *nsp1_40 = buffer.data(nsp1 + 40);
    const auto *nsp1_42 = buffer.data(nsp1 + 42);
    const auto *nsp1_43 = buffer.data(nsp1 + 43);
    const auto *nsp1_44 = buffer.data(nsp1 + 44);
    const auto *nsp1_45 = buffer.data(nsp1 + 45);
    const auto *nsp1_46 = buffer.data(nsp1 + 46);
    const auto *nsp1_47 = buffer.data(nsp1 + 47);
    const auto *nsp1_50 = buffer.data(nsp1 + 50);
    const auto *nsp1_51 = buffer.data(nsp1 + 51);
    const auto *nsp1_52 = buffer.data(nsp1 + 52);
    const auto *nsp1_53 = buffer.data(nsp1 + 53);
    const auto *nsp1_54 = buffer.data(nsp1 + 54);
    const auto *nsp1_55 = buffer.data(nsp1 + 55);
    const auto *nsp1_56 = buffer.data(nsp1 + 56);
    const auto *nsp1_58 = buffer.data(nsp1 + 58);
    const auto *nsp1_60 = buffer.data(nsp1 + 60);
    const auto *nsp1_61 = buffer.data(nsp1 + 61);
    const auto *nsp1_62 = buffer.data(nsp1 + 62);
    const auto *nsp1_63 = buffer.data(nsp1 + 63);
    const auto *nsp1_64 = buffer.data(nsp1 + 64);
    const auto *nsp1_65 = buffer.data(nsp1 + 65);
    const auto *nsp1_68 = buffer.data(nsp1 + 68);
    const auto *nsp1_69 = buffer.data(nsp1 + 69);
    const auto *nsp1_70 = buffer.data(nsp1 + 70);
    const auto *nsp1_71 = buffer.data(nsp1 + 71);
    const auto *nsp1_72 = buffer.data(nsp1 + 72);
    const auto *nsp1_73 = buffer.data(nsp1 + 73);
    const auto *nsp1_74 = buffer.data(nsp1 + 74);
    const auto *nsp1_75 = buffer.data(nsp1 + 75);

    const auto *nsd_81 = buffer.data(nsd + 81);
    const auto *nsd_82 = buffer.data(nsd + 82);
    const auto *nsd_83 = buffer.data(nsd + 83);
    const auto *nsd_84 = buffer.data(nsd + 84);
    const auto *nsd_86 = buffer.data(nsd + 86);
    const auto *nsd_87 = buffer.data(nsd + 87);
    const auto *nsd_88 = buffer.data(nsd + 88);
    const auto *nsd_89 = buffer.data(nsd + 89);
    const auto *nsd_90 = buffer.data(nsd + 90);
    const auto *nsd_91 = buffer.data(nsd + 91);
    const auto *nsd_93 = buffer.data(nsd + 93);
    const auto *nsd_95 = buffer.data(nsd + 95);
    const auto *nsd_96 = buffer.data(nsd + 96);
    const auto *nsd_99 = buffer.data(nsd + 99);
    const auto *nsd_100 = buffer.data(nsd + 100);
    const auto *nsd_101 = buffer.data(nsd + 101);
    const auto *nsd_102 = buffer.data(nsd + 102);
    const auto *nsd_105 = buffer.data(nsd + 105);
    const auto *nsd_106 = buffer.data(nsd + 106);
    const auto *nsd_107 = buffer.data(nsd + 107);
    const auto *nsd_108 = buffer.data(nsd + 108);
    const auto *nsd_111 = buffer.data(nsd + 111);
    const auto *nsd_112 = buffer.data(nsd + 112);
    const auto *nsd_113 = buffer.data(nsd + 113);
    const auto *nsd_114 = buffer.data(nsd + 114);
    const auto *nsd_117 = buffer.data(nsd + 117);
    const auto *nsd_118 = buffer.data(nsd + 118);
    const auto *nsd_119 = buffer.data(nsd + 119);
    const auto *nsd_120 = buffer.data(nsd + 120);
    const auto *nsd_122 = buffer.data(nsd + 122);
    const auto *nsd_123 = buffer.data(nsd + 123);
    const auto *nsd_124 = buffer.data(nsd + 124);
    const auto *nsd_125 = buffer.data(nsd + 125);
    const auto *nsd_126 = buffer.data(nsd + 126);
    const auto *nsd_127 = buffer.data(nsd + 127);
    const auto *nsd_129 = buffer.data(nsd + 129);
    const auto *nsd_131 = buffer.data(nsd + 131);
    const auto *nsd_132 = buffer.data(nsd + 132);
    const auto *nsd_135 = buffer.data(nsd + 135);
    const auto *nsd_136 = buffer.data(nsd + 136);
    const auto *nsd_137 = buffer.data(nsd + 137);
    const auto *nsd_138 = buffer.data(nsd + 138);
    const auto *nsd_141 = buffer.data(nsd + 141);
    const auto *nsd_142 = buffer.data(nsd + 142);
    const auto *nsd_143 = buffer.data(nsd + 143);
    const auto *nsd_144 = buffer.data(nsd + 144);
    const auto *nsd_147 = buffer.data(nsd + 147);
    const auto *nsd_148 = buffer.data(nsd + 148);
    const auto *nsd_149 = buffer.data(nsd + 149);
    const auto *nsd_150 = buffer.data(nsd + 150);
    const auto *nsd_153 = buffer.data(nsd + 153);
    const auto *nsd_154 = buffer.data(nsd + 154);
    const auto *nsd_155 = buffer.data(nsd + 155);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, msd_57, msd_81, msd_82, \
                         msd_83, nsp0_40, nsp1_40, nsd_81, nsd_82, \
                         nsd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_13 * msd_81[k]
                   + f_3 * pc_x[k] * nsd_81[k];

        t_134[k] = f_13 * msd_82[k]
                   + f_3 * pc_x[k] * nsd_82[k];

        t_135[k] = f_13 * msd_83[k]
                   + f_3 * pc_x[k] * nsd_83[k];

        t_136[k] = f_5 * msd_57[k]
                   + f_1 * nsp0_40[k]
                   - f_2 * nsp1_40[k]
                   + f_3 * pc_y[k] * nsd_81[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pc_y, pc_z, msf0_99, msd_51, msd_59, \
                         msf1_99, nsd_81, nsd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_12 * msd_51[k]
                   + f_3 * pc_z[k] * nsd_81[k];

        t_138[k] = f_5 * msd_59[k]
                   + f_3 * pc_y[k] * nsd_83[k];

        t_139[k] = pa_y[k] * msf0_99[k]
                   - f_4 * pc_y[k] * msf1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, msd_54, msd_84, \
                         msd_87, nsp0_42, nsp1_42, nsd_84, nsd_86, \
                         nsd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * msd_84[k]
                   + f_1 * nsp0_42[k]
                   - f_2 * nsp1_42[k]
                   + f_3 * pc_x[k] * nsd_84[k];

        t_141[k] = f_3 * pc_y[k] * nsd_84[k];

        t_142[k] = f_14 * msd_54[k]
                   + f_3 * pc_z[k] * nsd_84[k];

        t_143[k] = f_13 * msd_87[k]
                   + f_3 * pc_x[k] * nsd_87[k];

        t_144[k] = f_3 * pc_y[k] * nsd_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, msd_89, nsp0_43, nsp0_44, \
                         nsp1_43, nsp1_44, nsd_87, nsd_88, nsd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * msd_89[k]
                   + f_3 * pc_x[k] * nsd_89[k];

        t_146[k] = f_1 * nsp0_43[k]
                   - f_2 * nsp1_43[k]
                   + f_3 * pc_y[k] * nsd_87[k];

        t_147[k] = f_7 * nsp0_44[k]
                   - f_8 * nsp1_44[k]
                   + f_3 * pc_y[k] * nsd_88[k];

        t_148[k] = f_3 * pc_y[k] * nsd_89[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pc_x, pc_y, pc_z, msd_59, msd_60, msd_90, \
                         nsp0_44, nsp0_45, nsp1_44, nsp1_45, nsd_89, \
                         nsd_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_14 * msd_59[k]
                   + f_1 * nsp0_44[k]
                   - f_2 * nsp1_44[k]
                   + f_3 * pc_z[k] * nsd_89[k];

        t_150[k] = f_15 * msd_90[k]
                   + f_1 * nsp0_45[k]
                   - f_2 * nsp1_45[k]
                   + f_3 * pc_x[k] * nsd_90[k];

        t_151[k] = f_15 * msd_60[k]
                   + f_3 * pc_y[k] * nsd_90[k];

        t_152[k] = f_3 * pc_z[k] * nsd_90[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, msd_63, msd_93, \
                         msd_95, nsp0_46, nsp1_46, nsd_91, nsd_93, \
                         nsd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_15 * msd_93[k]
                   + f_3 * pc_x[k] * nsd_93[k];

        t_154[k] = f_3 * pc_z[k] * nsd_91[k];

        t_155[k] = f_15 * msd_95[k]
                   + f_3 * pc_x[k] * nsd_95[k];

        t_156[k] = f_15 * msd_63[k]
                   + f_1 * nsp0_46[k]
                   - f_2 * nsp1_46[k]
                   + f_3 * pc_y[k] * nsd_93[k];

        t_157[k] = f_3 * pc_z[k] * nsd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_z, pc_y, pc_z, msf0_100, msd_65, \
                         msd_66, msf1_100, nsp0_47, nsp1_47, nsd_95, \
                         nsd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_15 * msd_65[k]
                   + f_3 * pc_y[k] * nsd_95[k];

        t_159[k] = f_1 * nsp0_47[k]
                   - f_2 * nsp1_47[k]
                   + f_3 * pc_z[k] * nsd_95[k];

        t_160[k] = pa_z[k] * msf0_100[k]
                   - f_4 * pc_z[k] * msf1_100[k];

        t_161[k] = f_14 * msd_66[k]
                   + f_3 * pc_y[k] * nsd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, msd_60, msd_99, msd_100, \
                         msd_101, nsd_96, nsd_99, nsd_100, nsd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * msd_60[k]
                   + f_3 * pc_z[k] * nsd_96[k];

        t_163[k] = f_15 * msd_99[k]
                   + f_3 * pc_x[k] * nsd_99[k];

        t_164[k] = f_15 * msd_100[k]
                   + f_3 * pc_x[k] * nsd_100[k];

        t_165[k] = f_15 * msd_101[k]
                   + f_3 * pc_x[k] * nsd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, pc_y, pc_z, msf0_106, msd_63, \
                         msd_65, msd_71, msf1_106, nsp0_50, nsp1_50, nsd_99, \
                         nsd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_z[k] * msf0_106[k]
                   - f_4 * pc_z[k] * msf1_106[k];

        t_167[k] = f_5 * msd_63[k]
                   + f_3 * pc_z[k] * nsd_99[k];

        t_168[k] = f_14 * msd_71[k]
                   + f_3 * pc_y[k] * nsd_101[k];

        t_169[k] = f_5 * msd_65[k]
                   + f_1 * nsp0_50[k]
                   - f_2 * nsp1_50[k]
                   + f_3 * pc_z[k] * nsd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, msd_66, msd_72, \
                         msd_102, msd_105, nsp0_51, nsp1_51, nsd_102, \
                         nsd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_15 * msd_102[k]
                   + f_1 * nsp0_51[k]
                   - f_2 * nsp1_51[k]
                   + f_3 * pc_x[k] * nsd_102[k];

        t_171[k] = f_12 * msd_72[k]
                   + f_3 * pc_y[k] * nsd_102[k];

        t_172[k] = f_10 * msd_66[k]
                   + f_3 * pc_z[k] * nsd_102[k];

        t_173[k] = f_15 * msd_105[k]
                   + f_3 * pc_x[k] * nsd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, msd_69, msd_75, \
                         msd_106, msd_107, nsp0_52, nsp1_52, nsd_105, nsd_106, \
                         nsd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_15 * msd_106[k]
                   + f_3 * pc_x[k] * nsd_106[k];

        t_175[k] = f_15 * msd_107[k]
                   + f_3 * pc_x[k] * nsd_107[k];

        t_176[k] = f_12 * msd_75[k]
                   + f_1 * nsp0_52[k]
                   - f_2 * nsp1_52[k]
                   + f_3 * pc_y[k] * nsd_105[k];

        t_177[k] = f_10 * msd_69[k]
                   + f_3 * pc_z[k] * nsd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, msd_71, msd_77, msd_108, \
                         nsp0_53, nsp0_54, nsp1_53, nsp1_54, nsd_107, \
                         nsd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_12 * msd_77[k]
                   + f_3 * pc_y[k] * nsd_107[k];

        t_179[k] = f_10 * msd_71[k]
                   + f_1 * nsp0_53[k]
                   - f_2 * nsp1_53[k]
                   + f_3 * pc_z[k] * nsd_107[k];

        t_180[k] = f_15 * msd_108[k]
                   + f_1 * nsp0_54[k]
                   - f_2 * nsp1_54[k]
                   + f_3 * pc_x[k] * nsd_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, msd_72, msd_78, \
                         msd_111, msd_112, nsd_108, nsd_111, nsd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_10 * msd_78[k]
                   + f_3 * pc_y[k] * nsd_108[k];

        t_182[k] = f_12 * msd_72[k]
                   + f_3 * pc_z[k] * nsd_108[k];

        t_183[k] = f_15 * msd_111[k]
                   + f_3 * pc_x[k] * nsd_111[k];

        t_184[k] = f_15 * msd_112[k]
                   + f_3 * pc_x[k] * nsd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, msd_75, msd_81, msd_83, \
                         msd_113, nsp0_55, nsp1_55, nsd_111, nsd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_15 * msd_113[k]
                   + f_3 * pc_x[k] * nsd_113[k];

        t_186[k] = f_10 * msd_81[k]
                   + f_1 * nsp0_55[k]
                   - f_2 * nsp1_55[k]
                   + f_3 * pc_y[k] * nsd_111[k];

        t_187[k] = f_12 * msd_75[k]
                   + f_3 * pc_z[k] * nsd_111[k];

        t_188[k] = f_10 * msd_83[k]
                   + f_3 * pc_y[k] * nsd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_y, pc_y, pc_z, msf0_140, msd_77, \
                         msd_78, msd_84, msf1_140, nsp0_56, nsp1_56, nsd_113, \
                         nsd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_12 * msd_77[k]
                   + f_1 * nsp0_56[k]
                   - f_2 * nsp1_56[k]
                   + f_3 * pc_z[k] * nsd_113[k];

        t_190[k] = pa_y[k] * msf0_140[k]
                   - f_4 * pc_y[k] * msf1_140[k];

        t_191[k] = f_5 * msd_84[k]
                   + f_3 * pc_y[k] * nsd_114[k];

        t_192[k] = f_14 * msd_78[k]
                   + f_3 * pc_z[k] * nsd_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, msd_87, msd_117, msd_118, \
                         msd_119, nsp0_58, nsp1_58, nsd_117, nsd_118, \
                         nsd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_15 * msd_117[k]
                   + f_3 * pc_x[k] * nsd_117[k];

        t_194[k] = f_15 * msd_118[k]
                   + f_3 * pc_x[k] * nsd_118[k];

        t_195[k] = f_15 * msd_119[k]
                   + f_3 * pc_x[k] * nsd_119[k];

        t_196[k] = f_5 * msd_87[k]
                   + f_1 * nsp0_58[k]
                   - f_2 * nsp1_58[k]
                   + f_3 * pc_y[k] * nsd_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pc_y, pc_z, msf0_149, msd_81, msd_89, \
                         msf1_149, nsd_117, nsd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_14 * msd_81[k]
                   + f_3 * pc_z[k] * nsd_117[k];

        t_198[k] = f_5 * msd_89[k]
                   + f_3 * pc_y[k] * nsd_119[k];

        t_199[k] = pa_y[k] * msf0_149[k]
                   - f_4 * pc_y[k] * msf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pc_x, pc_y, pc_z, msd_84, msd_120, \
                         msd_123, nsp0_60, nsp1_60, nsd_120, nsd_122, \
                         nsd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_15 * msd_120[k]
                   + f_1 * nsp0_60[k]
                   - f_2 * nsp1_60[k]
                   + f_3 * pc_x[k] * nsd_120[k];

        t_201[k] = f_3 * pc_y[k] * nsd_120[k];

        t_202[k] = f_15 * msd_84[k]
                   + f_3 * pc_z[k] * nsd_120[k];

        t_203[k] = f_15 * msd_123[k]
                   + f_3 * pc_x[k] * nsd_123[k];

        t_204[k] = f_3 * pc_y[k] * nsd_122[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_x, pc_y, msd_125, nsp0_61, nsp0_62, \
                         nsp1_61, nsp1_62, nsd_123, nsd_124, nsd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_15 * msd_125[k]
                   + f_3 * pc_x[k] * nsd_125[k];

        t_206[k] = f_1 * nsp0_61[k]
                   - f_2 * nsp1_61[k]
                   + f_3 * pc_y[k] * nsd_123[k];

        t_207[k] = f_7 * nsp0_62[k]
                   - f_8 * nsp1_62[k]
                   + f_3 * pc_y[k] * nsd_124[k];

        t_208[k] = f_3 * pc_y[k] * nsd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pc_x, pc_y, pc_z, msd_89, msd_90, \
                         msd_126, nsp0_62, nsp0_63, nsp1_62, nsp1_63, nsd_125, \
                         nsd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_15 * msd_89[k]
                   + f_1 * nsp0_62[k]
                   - f_2 * nsp1_62[k]
                   + f_3 * pc_z[k] * nsd_125[k];

        t_210[k] = f_14 * msd_126[k]
                   + f_1 * nsp0_63[k]
                   - f_2 * nsp1_63[k]
                   + f_3 * pc_x[k] * nsd_126[k];

        t_211[k] = f_13 * msd_90[k]
                   + f_3 * pc_y[k] * nsd_126[k];

        t_212[k] = f_3 * pc_z[k] * nsd_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, pc_x, pc_y, pc_z, msd_93, msd_129, \
                         msd_131, nsp0_64, nsp1_64, nsd_127, nsd_129, \
                         nsd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_14 * msd_129[k]
                   + f_3 * pc_x[k] * nsd_129[k];

        t_214[k] = f_3 * pc_z[k] * nsd_127[k];

        t_215[k] = f_14 * msd_131[k]
                   + f_3 * pc_x[k] * nsd_131[k];

        t_216[k] = f_13 * msd_93[k]
                   + f_1 * nsp0_64[k]
                   - f_2 * nsp1_64[k]
                   + f_3 * pc_y[k] * nsd_129[k];

        t_217[k] = f_3 * pc_z[k] * nsd_129[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pa_z, pc_y, pc_z, msf0_150, msd_95, \
                         msd_96, msf1_150, nsp0_65, nsp1_65, nsd_131, \
                         nsd_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_13 * msd_95[k]
                   + f_3 * pc_y[k] * nsd_131[k];

        t_219[k] = f_1 * nsp0_65[k]
                   - f_2 * nsp1_65[k]
                   + f_3 * pc_z[k] * nsd_131[k];

        t_220[k] = pa_z[k] * msf0_150[k]
                   - f_4 * pc_z[k] * msf1_150[k];

        t_221[k] = f_15 * msd_96[k]
                   + f_3 * pc_y[k] * nsd_132[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_z, msd_90, msd_135, msd_136, \
                         msd_137, nsd_132, nsd_135, nsd_136, nsd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_5 * msd_90[k]
                   + f_3 * pc_z[k] * nsd_132[k];

        t_223[k] = f_14 * msd_135[k]
                   + f_3 * pc_x[k] * nsd_135[k];

        t_224[k] = f_14 * msd_136[k]
                   + f_3 * pc_x[k] * nsd_136[k];

        t_225[k] = f_14 * msd_137[k]
                   + f_3 * pc_x[k] * nsd_137[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pa_z, pc_y, pc_z, msf0_156, msd_93, \
                         msd_95, msd_101, msf1_156, nsp0_68, nsp1_68, nsd_135, \
                         nsd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pa_z[k] * msf0_156[k]
                   - f_4 * pc_z[k] * msf1_156[k];

        t_227[k] = f_5 * msd_93[k]
                   + f_3 * pc_z[k] * nsd_135[k];

        t_228[k] = f_15 * msd_101[k]
                   + f_3 * pc_y[k] * nsd_137[k];

        t_229[k] = f_5 * msd_95[k]
                   + f_1 * nsp0_68[k]
                   - f_2 * nsp1_68[k]
                   + f_3 * pc_z[k] * nsd_137[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, pc_z, msd_96, msd_102, \
                         msd_138, msd_141, nsp0_69, nsp1_69, nsd_138, \
                         nsd_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_14 * msd_138[k]
                   + f_1 * nsp0_69[k]
                   - f_2 * nsp1_69[k]
                   + f_3 * pc_x[k] * nsd_138[k];

        t_231[k] = f_14 * msd_102[k]
                   + f_3 * pc_y[k] * nsd_138[k];

        t_232[k] = f_10 * msd_96[k]
                   + f_3 * pc_z[k] * nsd_138[k];

        t_233[k] = f_14 * msd_141[k]
                   + f_3 * pc_x[k] * nsd_141[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, msd_99, msd_105, \
                         msd_142, msd_143, nsp0_70, nsp1_70, nsd_141, nsd_142, \
                         nsd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_14 * msd_142[k]
                   + f_3 * pc_x[k] * nsd_142[k];

        t_235[k] = f_14 * msd_143[k]
                   + f_3 * pc_x[k] * nsd_143[k];

        t_236[k] = f_14 * msd_105[k]
                   + f_1 * nsp0_70[k]
                   - f_2 * nsp1_70[k]
                   + f_3 * pc_y[k] * nsd_141[k];

        t_237[k] = f_10 * msd_99[k]
                   + f_3 * pc_z[k] * nsd_141[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pc_x, pc_y, pc_z, msd_101, msd_107, msd_144, \
                         nsp0_71, nsp0_72, nsp1_71, nsp1_72, nsd_143, \
                         nsd_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_14 * msd_107[k]
                   + f_3 * pc_y[k] * nsd_143[k];

        t_239[k] = f_10 * msd_101[k]
                   + f_1 * nsp0_71[k]
                   - f_2 * nsp1_71[k]
                   + f_3 * pc_z[k] * nsd_143[k];

        t_240[k] = f_14 * msd_144[k]
                   + f_1 * nsp0_72[k]
                   - f_2 * nsp1_72[k]
                   + f_3 * pc_x[k] * nsd_144[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pc_x, pc_y, pc_z, msd_102, msd_108, \
                         msd_147, msd_148, nsd_144, nsd_147, nsd_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_12 * msd_108[k]
                   + f_3 * pc_y[k] * nsd_144[k];

        t_242[k] = f_12 * msd_102[k]
                   + f_3 * pc_z[k] * nsd_144[k];

        t_243[k] = f_14 * msd_147[k]
                   + f_3 * pc_x[k] * nsd_147[k];

        t_244[k] = f_14 * msd_148[k]
                   + f_3 * pc_x[k] * nsd_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, pc_y, pc_z, msd_105, msd_111, \
                         msd_113, msd_149, nsp0_73, nsp1_73, nsd_147, \
                         nsd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_14 * msd_149[k]
                   + f_3 * pc_x[k] * nsd_149[k];

        t_246[k] = f_12 * msd_111[k]
                   + f_1 * nsp0_73[k]
                   - f_2 * nsp1_73[k]
                   + f_3 * pc_y[k] * nsd_147[k];

        t_247[k] = f_12 * msd_105[k]
                   + f_3 * pc_z[k] * nsd_147[k];

        t_248[k] = f_12 * msd_113[k]
                   + f_3 * pc_y[k] * nsd_149[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_x, pc_y, pc_z, msd_107, msd_114, msd_150, \
                         nsp0_74, nsp0_75, nsp1_74, nsp1_75, nsd_149, \
                         nsd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_12 * msd_107[k]
                   + f_1 * nsp0_74[k]
                   - f_2 * nsp1_74[k]
                   + f_3 * pc_z[k] * nsd_149[k];

        t_250[k] = f_14 * msd_150[k]
                   + f_1 * nsp0_75[k]
                   - f_2 * nsp1_75[k]
                   + f_3 * pc_x[k] * nsd_150[k];

        t_251[k] = f_10 * msd_114[k]
                   + f_3 * pc_y[k] * nsd_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, pc_z, msd_108, msd_153, msd_154, \
                         msd_155, nsd_150, nsd_153, nsd_154, nsd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_14 * msd_108[k]
                   + f_3 * pc_z[k] * nsd_150[k];

        t_253[k] = f_14 * msd_153[k]
                   + f_3 * pc_x[k] * nsd_153[k];

        t_254[k] = f_14 * msd_154[k]
                   + f_3 * pc_x[k] * nsd_154[k];

        t_255[k] = f_14 * msd_155[k]
                   + f_3 * pc_x[k] * nsd_155[k];
    }
}

static auto
compute_prim_nsf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msf0,
                                                          const size_t msd, const size_t msf1,
                                                          const size_t nsp0, const size_t nsp1,
                                                          const size_t nsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msf0_200 = buffer.data(msf0 + 200);
    const auto *msf0_209 = buffer.data(msf0 + 209);
    const auto *msf0_210 = buffer.data(msf0 + 210);
    const auto *msf0_216 = buffer.data(msf0 + 216);
    const auto *msf0_270 = buffer.data(msf0 + 270);
    const auto *msf0_279 = buffer.data(msf0 + 279);
    const auto *msf0_280 = buffer.data(msf0 + 280);
    const auto *msf0_286 = buffer.data(msf0 + 286);

    const auto *msd_111 = buffer.data(msd + 111);
    const auto *msd_113 = buffer.data(msd + 113);
    const auto *msd_114 = buffer.data(msd + 114);
    const auto *msd_117 = buffer.data(msd + 117);
    const auto *msd_119 = buffer.data(msd + 119);
    const auto *msd_120 = buffer.data(msd + 120);
    const auto *msd_123 = buffer.data(msd + 123);
    const auto *msd_125 = buffer.data(msd + 125);
    const auto *msd_126 = buffer.data(msd + 126);
    const auto *msd_129 = buffer.data(msd + 129);
    const auto *msd_131 = buffer.data(msd + 131);
    const auto *msd_132 = buffer.data(msd + 132);
    const auto *msd_135 = buffer.data(msd + 135);
    const auto *msd_137 = buffer.data(msd + 137);
    const auto *msd_138 = buffer.data(msd + 138);
    const auto *msd_141 = buffer.data(msd + 141);
    const auto *msd_143 = buffer.data(msd + 143);
    const auto *msd_144 = buffer.data(msd + 144);
    const auto *msd_147 = buffer.data(msd + 147);
    const auto *msd_149 = buffer.data(msd + 149);
    const auto *msd_150 = buffer.data(msd + 150);
    const auto *msd_153 = buffer.data(msd + 153);
    const auto *msd_155 = buffer.data(msd + 155);
    const auto *msd_156 = buffer.data(msd + 156);
    const auto *msd_159 = buffer.data(msd + 159);
    const auto *msd_160 = buffer.data(msd + 160);
    const auto *msd_161 = buffer.data(msd + 161);
    const auto *msd_162 = buffer.data(msd + 162);
    const auto *msd_165 = buffer.data(msd + 165);
    const auto *msd_167 = buffer.data(msd + 167);
    const auto *msd_168 = buffer.data(msd + 168);
    const auto *msd_171 = buffer.data(msd + 171);
    const auto *msd_173 = buffer.data(msd + 173);
    const auto *msd_174 = buffer.data(msd + 174);
    const auto *msd_177 = buffer.data(msd + 177);
    const auto *msd_178 = buffer.data(msd + 178);
    const auto *msd_179 = buffer.data(msd + 179);
    const auto *msd_180 = buffer.data(msd + 180);
    const auto *msd_183 = buffer.data(msd + 183);
    const auto *msd_184 = buffer.data(msd + 184);
    const auto *msd_185 = buffer.data(msd + 185);
    const auto *msd_186 = buffer.data(msd + 186);
    const auto *msd_189 = buffer.data(msd + 189);
    const auto *msd_190 = buffer.data(msd + 190);
    const auto *msd_191 = buffer.data(msd + 191);
    const auto *msd_192 = buffer.data(msd + 192);
    const auto *msd_195 = buffer.data(msd + 195);
    const auto *msd_196 = buffer.data(msd + 196);
    const auto *msd_197 = buffer.data(msd + 197);
    const auto *msd_198 = buffer.data(msd + 198);
    const auto *msd_201 = buffer.data(msd + 201);
    const auto *msd_202 = buffer.data(msd + 202);
    const auto *msd_203 = buffer.data(msd + 203);
    const auto *msd_207 = buffer.data(msd + 207);
    const auto *msd_208 = buffer.data(msd + 208);
    const auto *msd_209 = buffer.data(msd + 209);
    const auto *msd_210 = buffer.data(msd + 210);
    const auto *msd_213 = buffer.data(msd + 213);
    const auto *msd_215 = buffer.data(msd + 215);
    const auto *msd_216 = buffer.data(msd + 216);
    const auto *msd_219 = buffer.data(msd + 219);
    const auto *msd_221 = buffer.data(msd + 221);
    const auto *msd_225 = buffer.data(msd + 225);
    const auto *msd_226 = buffer.data(msd + 226);
    const auto *msd_227 = buffer.data(msd + 227);

    const auto *msf1_200 = buffer.data(msf1 + 200);
    const auto *msf1_209 = buffer.data(msf1 + 209);
    const auto *msf1_210 = buffer.data(msf1 + 210);
    const auto *msf1_216 = buffer.data(msf1 + 216);
    const auto *msf1_270 = buffer.data(msf1 + 270);
    const auto *msf1_279 = buffer.data(msf1 + 279);
    const auto *msf1_280 = buffer.data(msf1 + 280);
    const auto *msf1_286 = buffer.data(msf1 + 286);

    const auto *nsp0_76 = buffer.data(nsp0 + 76);
    const auto *nsp0_77 = buffer.data(nsp0 + 77);
    const auto *nsp0_79 = buffer.data(nsp0 + 79);
    const auto *nsp0_81 = buffer.data(nsp0 + 81);
    const auto *nsp0_82 = buffer.data(nsp0 + 82);
    const auto *nsp0_83 = buffer.data(nsp0 + 83);
    const auto *nsp0_84 = buffer.data(nsp0 + 84);
    const auto *nsp0_85 = buffer.data(nsp0 + 85);
    const auto *nsp0_86 = buffer.data(nsp0 + 86);
    const auto *nsp0_89 = buffer.data(nsp0 + 89);
    const auto *nsp0_90 = buffer.data(nsp0 + 90);
    const auto *nsp0_91 = buffer.data(nsp0 + 91);
    const auto *nsp0_92 = buffer.data(nsp0 + 92);
    const auto *nsp0_93 = buffer.data(nsp0 + 93);
    const auto *nsp0_94 = buffer.data(nsp0 + 94);
    const auto *nsp0_95 = buffer.data(nsp0 + 95);
    const auto *nsp0_96 = buffer.data(nsp0 + 96);
    const auto *nsp0_97 = buffer.data(nsp0 + 97);
    const auto *nsp0_98 = buffer.data(nsp0 + 98);
    const auto *nsp0_99 = buffer.data(nsp0 + 99);
    const auto *nsp0_100 = buffer.data(nsp0 + 100);
    const auto *nsp0_101 = buffer.data(nsp0 + 101);
    const auto *nsp0_103 = buffer.data(nsp0 + 103);
    const auto *nsp0_105 = buffer.data(nsp0 + 105);
    const auto *nsp0_106 = buffer.data(nsp0 + 106);
    const auto *nsp0_107 = buffer.data(nsp0 + 107);
    const auto *nsp0_108 = buffer.data(nsp0 + 108);
    const auto *nsp0_109 = buffer.data(nsp0 + 109);
    const auto *nsp0_110 = buffer.data(nsp0 + 110);
    const auto *nsp0_113 = buffer.data(nsp0 + 113);

    const auto *nsp1_76 = buffer.data(nsp1 + 76);
    const auto *nsp1_77 = buffer.data(nsp1 + 77);
    const auto *nsp1_79 = buffer.data(nsp1 + 79);
    const auto *nsp1_81 = buffer.data(nsp1 + 81);
    const auto *nsp1_82 = buffer.data(nsp1 + 82);
    const auto *nsp1_83 = buffer.data(nsp1 + 83);
    const auto *nsp1_84 = buffer.data(nsp1 + 84);
    const auto *nsp1_85 = buffer.data(nsp1 + 85);
    const auto *nsp1_86 = buffer.data(nsp1 + 86);
    const auto *nsp1_89 = buffer.data(nsp1 + 89);
    const auto *nsp1_90 = buffer.data(nsp1 + 90);
    const auto *nsp1_91 = buffer.data(nsp1 + 91);
    const auto *nsp1_92 = buffer.data(nsp1 + 92);
    const auto *nsp1_93 = buffer.data(nsp1 + 93);
    const auto *nsp1_94 = buffer.data(nsp1 + 94);
    const auto *nsp1_95 = buffer.data(nsp1 + 95);
    const auto *nsp1_96 = buffer.data(nsp1 + 96);
    const auto *nsp1_97 = buffer.data(nsp1 + 97);
    const auto *nsp1_98 = buffer.data(nsp1 + 98);
    const auto *nsp1_99 = buffer.data(nsp1 + 99);
    const auto *nsp1_100 = buffer.data(nsp1 + 100);
    const auto *nsp1_101 = buffer.data(nsp1 + 101);
    const auto *nsp1_103 = buffer.data(nsp1 + 103);
    const auto *nsp1_105 = buffer.data(nsp1 + 105);
    const auto *nsp1_106 = buffer.data(nsp1 + 106);
    const auto *nsp1_107 = buffer.data(nsp1 + 107);
    const auto *nsp1_108 = buffer.data(nsp1 + 108);
    const auto *nsp1_109 = buffer.data(nsp1 + 109);
    const auto *nsp1_110 = buffer.data(nsp1 + 110);
    const auto *nsp1_113 = buffer.data(nsp1 + 113);

    const auto *nsd_153 = buffer.data(nsd + 153);
    const auto *nsd_155 = buffer.data(nsd + 155);
    const auto *nsd_156 = buffer.data(nsd + 156);
    const auto *nsd_159 = buffer.data(nsd + 159);
    const auto *nsd_160 = buffer.data(nsd + 160);
    const auto *nsd_161 = buffer.data(nsd + 161);
    const auto *nsd_162 = buffer.data(nsd + 162);
    const auto *nsd_164 = buffer.data(nsd + 164);
    const auto *nsd_165 = buffer.data(nsd + 165);
    const auto *nsd_166 = buffer.data(nsd + 166);
    const auto *nsd_167 = buffer.data(nsd + 167);
    const auto *nsd_168 = buffer.data(nsd + 168);
    const auto *nsd_169 = buffer.data(nsd + 169);
    const auto *nsd_171 = buffer.data(nsd + 171);
    const auto *nsd_173 = buffer.data(nsd + 173);
    const auto *nsd_174 = buffer.data(nsd + 174);
    const auto *nsd_177 = buffer.data(nsd + 177);
    const auto *nsd_178 = buffer.data(nsd + 178);
    const auto *nsd_179 = buffer.data(nsd + 179);
    const auto *nsd_180 = buffer.data(nsd + 180);
    const auto *nsd_183 = buffer.data(nsd + 183);
    const auto *nsd_184 = buffer.data(nsd + 184);
    const auto *nsd_185 = buffer.data(nsd + 185);
    const auto *nsd_186 = buffer.data(nsd + 186);
    const auto *nsd_189 = buffer.data(nsd + 189);
    const auto *nsd_190 = buffer.data(nsd + 190);
    const auto *nsd_191 = buffer.data(nsd + 191);
    const auto *nsd_192 = buffer.data(nsd + 192);
    const auto *nsd_195 = buffer.data(nsd + 195);
    const auto *nsd_196 = buffer.data(nsd + 196);
    const auto *nsd_197 = buffer.data(nsd + 197);
    const auto *nsd_198 = buffer.data(nsd + 198);
    const auto *nsd_201 = buffer.data(nsd + 201);
    const auto *nsd_202 = buffer.data(nsd + 202);
    const auto *nsd_203 = buffer.data(nsd + 203);
    const auto *nsd_204 = buffer.data(nsd + 204);
    const auto *nsd_207 = buffer.data(nsd + 207);
    const auto *nsd_208 = buffer.data(nsd + 208);
    const auto *nsd_209 = buffer.data(nsd + 209);
    const auto *nsd_210 = buffer.data(nsd + 210);
    const auto *nsd_212 = buffer.data(nsd + 212);
    const auto *nsd_213 = buffer.data(nsd + 213);
    const auto *nsd_214 = buffer.data(nsd + 214);
    const auto *nsd_215 = buffer.data(nsd + 215);
    const auto *nsd_216 = buffer.data(nsd + 216);
    const auto *nsd_217 = buffer.data(nsd + 217);
    const auto *nsd_219 = buffer.data(nsd + 219);
    const auto *nsd_221 = buffer.data(nsd + 221);
    const auto *nsd_222 = buffer.data(nsd + 222);
    const auto *nsd_225 = buffer.data(nsd + 225);
    const auto *nsd_226 = buffer.data(nsd + 226);
    const auto *nsd_227 = buffer.data(nsd + 227);

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_y, pc_z, msd_111, msd_113, msd_117, \
                         msd_119, nsp0_76, nsp0_77, nsp1_76, nsp1_77, nsd_153, \
                         nsd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_10 * msd_117[k]
                   + f_1 * nsp0_76[k]
                   - f_2 * nsp1_76[k]
                   + f_3 * pc_y[k] * nsd_153[k];

        t_257[k] = f_14 * msd_111[k]
                   + f_3 * pc_z[k] * nsd_153[k];

        t_258[k] = f_10 * msd_119[k]
                   + f_3 * pc_y[k] * nsd_155[k];

        t_259[k] = f_14 * msd_113[k]
                   + f_1 * nsp0_77[k]
                   - f_2 * nsp1_77[k]
                   + f_3 * pc_z[k] * nsd_155[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, pc_z, msf0_200, \
                         msd_114, msd_120, msd_159, msf1_200, nsd_156, \
                         nsd_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * msf0_200[k]
                   - f_4 * pc_y[k] * msf1_200[k];

        t_261[k] = f_5 * msd_120[k]
                   + f_3 * pc_y[k] * nsd_156[k];

        t_262[k] = f_15 * msd_114[k]
                   + f_3 * pc_z[k] * nsd_156[k];

        t_263[k] = f_14 * msd_159[k]
                   + f_3 * pc_x[k] * nsd_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, pc_z, msd_117, msd_123, \
                         msd_160, msd_161, nsp0_79, nsp1_79, nsd_159, nsd_160, \
                         nsd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_14 * msd_160[k]
                   + f_3 * pc_x[k] * nsd_160[k];

        t_265[k] = f_14 * msd_161[k]
                   + f_3 * pc_x[k] * nsd_161[k];

        t_266[k] = f_5 * msd_123[k]
                   + f_1 * nsp0_79[k]
                   - f_2 * nsp1_79[k]
                   + f_3 * pc_y[k] * nsd_159[k];

        t_267[k] = f_15 * msd_117[k]
                   + f_3 * pc_z[k] * nsd_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_y, pc_x, pc_y, msf0_209, msd_125, \
                         msd_162, msf1_209, nsp0_81, nsp1_81, nsd_161, \
                         nsd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_5 * msd_125[k]
                   + f_3 * pc_y[k] * nsd_161[k];

        t_269[k] = pa_y[k] * msf0_209[k]
                   - f_4 * pc_y[k] * msf1_209[k];

        t_270[k] = f_14 * msd_162[k]
                   + f_1 * nsp0_81[k]
                   - f_2 * nsp1_81[k]
                   + f_3 * pc_x[k] * nsd_162[k];

        t_271[k] = f_3 * pc_y[k] * nsd_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, msd_120, msd_165, \
                         msd_167, nsd_162, nsd_164, nsd_165, nsd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_13 * msd_120[k]
                   + f_3 * pc_z[k] * nsd_162[k];

        t_273[k] = f_14 * msd_165[k]
                   + f_3 * pc_x[k] * nsd_165[k];

        t_274[k] = f_3 * pc_y[k] * nsd_164[k];

        t_275[k] = f_14 * msd_167[k]
                   + f_3 * pc_x[k] * nsd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, msd_125, nsp0_82, nsp0_83, \
                         nsp1_82, nsp1_83, nsd_165, nsd_166, nsd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * nsp0_82[k]
                   - f_2 * nsp1_82[k]
                   + f_3 * pc_y[k] * nsd_165[k];

        t_277[k] = f_7 * nsp0_83[k]
                   - f_8 * nsp1_83[k]
                   + f_3 * pc_y[k] * nsd_166[k];

        t_278[k] = f_3 * pc_y[k] * nsd_167[k];

        t_279[k] = f_13 * msd_125[k]
                   + f_1 * nsp0_83[k]
                   - f_2 * nsp1_83[k]
                   + f_3 * pc_z[k] * nsd_167[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pc_x, pc_y, pc_z, msd_126, \
                         msd_168, msd_171, nsp0_84, nsp1_84, nsd_168, nsd_169, \
                         nsd_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_12 * msd_168[k]
                   + f_1 * nsp0_84[k]
                   - f_2 * nsp1_84[k]
                   + f_3 * pc_x[k] * nsd_168[k];

        t_281[k] = f_11 * msd_126[k]
                   + f_3 * pc_y[k] * nsd_168[k];

        t_282[k] = f_3 * pc_z[k] * nsd_168[k];

        t_283[k] = f_12 * msd_171[k]
                   + f_3 * pc_x[k] * nsd_171[k];

        t_284[k] = f_3 * pc_z[k] * nsd_169[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pc_x, pc_y, pc_z, msd_129, msd_131, \
                         msd_173, nsp0_85, nsp1_85, nsd_171, nsd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_12 * msd_173[k]
                   + f_3 * pc_x[k] * nsd_173[k];

        t_286[k] = f_11 * msd_129[k]
                   + f_1 * nsp0_85[k]
                   - f_2 * nsp1_85[k]
                   + f_3 * pc_y[k] * nsd_171[k];

        t_287[k] = f_3 * pc_z[k] * nsd_171[k];

        t_288[k] = f_11 * msd_131[k]
                   + f_3 * pc_y[k] * nsd_173[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_z, pc_y, pc_z, msf0_210, msd_126, \
                         msd_132, msf1_210, nsp0_86, nsp1_86, nsd_173, \
                         nsd_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_1 * nsp0_86[k]
                   - f_2 * nsp1_86[k]
                   + f_3 * pc_z[k] * nsd_173[k];

        t_290[k] = pa_z[k] * msf0_210[k]
                   - f_4 * pc_z[k] * msf1_210[k];

        t_291[k] = f_13 * msd_132[k]
                   + f_3 * pc_y[k] * nsd_174[k];

        t_292[k] = f_5 * msd_126[k]
                   + f_3 * pc_z[k] * nsd_174[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_z, pc_x, pc_z, msf0_216, msd_177, \
                         msd_178, msd_179, msf1_216, nsd_177, nsd_178, \
                         nsd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_12 * msd_177[k]
                   + f_3 * pc_x[k] * nsd_177[k];

        t_294[k] = f_12 * msd_178[k]
                   + f_3 * pc_x[k] * nsd_178[k];

        t_295[k] = f_12 * msd_179[k]
                   + f_3 * pc_x[k] * nsd_179[k];

        t_296[k] = pa_z[k] * msf0_216[k]
                   - f_4 * pc_z[k] * msf1_216[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pc_y, pc_z, msd_129, msd_131, msd_137, nsp0_89, \
                         nsp1_89, nsd_177, nsd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_5 * msd_129[k]
                   + f_3 * pc_z[k] * nsd_177[k];

        t_298[k] = f_13 * msd_137[k]
                   + f_3 * pc_y[k] * nsd_179[k];

        t_299[k] = f_5 * msd_131[k]
                   + f_1 * nsp0_89[k]
                   - f_2 * nsp1_89[k]
                   + f_3 * pc_z[k] * nsd_179[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, msd_132, msd_138, \
                         msd_180, msd_183, nsp0_90, nsp1_90, nsd_180, \
                         nsd_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_12 * msd_180[k]
                   + f_1 * nsp0_90[k]
                   - f_2 * nsp1_90[k]
                   + f_3 * pc_x[k] * nsd_180[k];

        t_301[k] = f_15 * msd_138[k]
                   + f_3 * pc_y[k] * nsd_180[k];

        t_302[k] = f_10 * msd_132[k]
                   + f_3 * pc_z[k] * nsd_180[k];

        t_303[k] = f_12 * msd_183[k]
                   + f_3 * pc_x[k] * nsd_183[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, pc_z, msd_135, msd_141, \
                         msd_184, msd_185, nsp0_91, nsp1_91, nsd_183, nsd_184, \
                         nsd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_12 * msd_184[k]
                   + f_3 * pc_x[k] * nsd_184[k];

        t_305[k] = f_12 * msd_185[k]
                   + f_3 * pc_x[k] * nsd_185[k];

        t_306[k] = f_15 * msd_141[k]
                   + f_1 * nsp0_91[k]
                   - f_2 * nsp1_91[k]
                   + f_3 * pc_y[k] * nsd_183[k];

        t_307[k] = f_10 * msd_135[k]
                   + f_3 * pc_z[k] * nsd_183[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pc_x, pc_y, pc_z, msd_137, msd_143, msd_186, \
                         nsp0_92, nsp0_93, nsp1_92, nsp1_93, nsd_185, \
                         nsd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_15 * msd_143[k]
                   + f_3 * pc_y[k] * nsd_185[k];

        t_309[k] = f_10 * msd_137[k]
                   + f_1 * nsp0_92[k]
                   - f_2 * nsp1_92[k]
                   + f_3 * pc_z[k] * nsd_185[k];

        t_310[k] = f_12 * msd_186[k]
                   + f_1 * nsp0_93[k]
                   - f_2 * nsp1_93[k]
                   + f_3 * pc_x[k] * nsd_186[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_x, pc_y, pc_z, msd_138, msd_144, \
                         msd_189, msd_190, nsd_186, nsd_189, nsd_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_14 * msd_144[k]
                   + f_3 * pc_y[k] * nsd_186[k];

        t_312[k] = f_12 * msd_138[k]
                   + f_3 * pc_z[k] * nsd_186[k];

        t_313[k] = f_12 * msd_189[k]
                   + f_3 * pc_x[k] * nsd_189[k];

        t_314[k] = f_12 * msd_190[k]
                   + f_3 * pc_x[k] * nsd_190[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, msd_141, msd_147, \
                         msd_149, msd_191, nsp0_94, nsp1_94, nsd_189, \
                         nsd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_12 * msd_191[k]
                   + f_3 * pc_x[k] * nsd_191[k];

        t_316[k] = f_14 * msd_147[k]
                   + f_1 * nsp0_94[k]
                   - f_2 * nsp1_94[k]
                   + f_3 * pc_y[k] * nsd_189[k];

        t_317[k] = f_12 * msd_141[k]
                   + f_3 * pc_z[k] * nsd_189[k];

        t_318[k] = f_14 * msd_149[k]
                   + f_3 * pc_y[k] * nsd_191[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pc_x, pc_y, pc_z, msd_143, msd_150, msd_192, \
                         nsp0_95, nsp0_96, nsp1_95, nsp1_96, nsd_191, \
                         nsd_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_12 * msd_143[k]
                   + f_1 * nsp0_95[k]
                   - f_2 * nsp1_95[k]
                   + f_3 * pc_z[k] * nsd_191[k];

        t_320[k] = f_12 * msd_192[k]
                   + f_1 * nsp0_96[k]
                   - f_2 * nsp1_96[k]
                   + f_3 * pc_x[k] * nsd_192[k];

        t_321[k] = f_12 * msd_150[k]
                   + f_3 * pc_y[k] * nsd_192[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, pc_z, msd_144, msd_195, msd_196, \
                         msd_197, nsd_192, nsd_195, nsd_196, nsd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_14 * msd_144[k]
                   + f_3 * pc_z[k] * nsd_192[k];

        t_323[k] = f_12 * msd_195[k]
                   + f_3 * pc_x[k] * nsd_195[k];

        t_324[k] = f_12 * msd_196[k]
                   + f_3 * pc_x[k] * nsd_196[k];

        t_325[k] = f_12 * msd_197[k]
                   + f_3 * pc_x[k] * nsd_197[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pc_y, pc_z, msd_147, msd_149, msd_153, \
                         msd_155, nsp0_97, nsp0_98, nsp1_97, nsp1_98, nsd_195, \
                         nsd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_12 * msd_153[k]
                   + f_1 * nsp0_97[k]
                   - f_2 * nsp1_97[k]
                   + f_3 * pc_y[k] * nsd_195[k];

        t_327[k] = f_14 * msd_147[k]
                   + f_3 * pc_z[k] * nsd_195[k];

        t_328[k] = f_12 * msd_155[k]
                   + f_3 * pc_y[k] * nsd_197[k];

        t_329[k] = f_14 * msd_149[k]
                   + f_1 * nsp0_98[k]
                   - f_2 * nsp1_98[k]
                   + f_3 * pc_z[k] * nsd_197[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pc_x, pc_y, pc_z, msd_150, msd_156, \
                         msd_198, msd_201, nsp0_99, nsp1_99, nsd_198, \
                         nsd_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_12 * msd_198[k]
                   + f_1 * nsp0_99[k]
                   - f_2 * nsp1_99[k]
                   + f_3 * pc_x[k] * nsd_198[k];

        t_331[k] = f_10 * msd_156[k]
                   + f_3 * pc_y[k] * nsd_198[k];

        t_332[k] = f_15 * msd_150[k]
                   + f_3 * pc_z[k] * nsd_198[k];

        t_333[k] = f_12 * msd_201[k]
                   + f_3 * pc_x[k] * nsd_201[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pc_x, pc_y, pc_z, msd_153, msd_159, \
                         msd_202, msd_203, nsp0_100, nsp1_100, nsd_201, nsd_202, \
                         nsd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_12 * msd_202[k]
                   + f_3 * pc_x[k] * nsd_202[k];

        t_335[k] = f_12 * msd_203[k]
                   + f_3 * pc_x[k] * nsd_203[k];

        t_336[k] = f_10 * msd_159[k]
                   + f_1 * nsp0_100[k]
                   - f_2 * nsp1_100[k]
                   + f_3 * pc_y[k] * nsd_201[k];

        t_337[k] = f_15 * msd_153[k]
                   + f_3 * pc_z[k] * nsd_201[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pa_y, pc_y, pc_z, msf0_270, msd_155, \
                         msd_161, msd_162, msf1_270, nsp0_101, nsp1_101, nsd_203, \
                         nsd_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_10 * msd_161[k]
                   + f_3 * pc_y[k] * nsd_203[k];

        t_339[k] = f_15 * msd_155[k]
                   + f_1 * nsp0_101[k]
                   - f_2 * nsp1_101[k]
                   + f_3 * pc_z[k] * nsd_203[k];

        t_340[k] = pa_y[k] * msf0_270[k]
                   - f_4 * pc_y[k] * msf1_270[k];

        t_341[k] = f_5 * msd_162[k]
                   + f_3 * pc_y[k] * nsd_204[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_z, msd_156, msd_207, msd_208, \
                         msd_209, nsd_204, nsd_207, nsd_208, nsd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_13 * msd_156[k]
                   + f_3 * pc_z[k] * nsd_204[k];

        t_343[k] = f_12 * msd_207[k]
                   + f_3 * pc_x[k] * nsd_207[k];

        t_344[k] = f_12 * msd_208[k]
                   + f_3 * pc_x[k] * nsd_208[k];

        t_345[k] = f_12 * msd_209[k]
                   + f_3 * pc_x[k] * nsd_209[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pa_y, pc_y, pc_z, msf0_279, msd_159, \
                         msd_165, msd_167, msf1_279, nsp0_103, nsp1_103, nsd_207, \
                         nsd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_5 * msd_165[k]
                   + f_1 * nsp0_103[k]
                   - f_2 * nsp1_103[k]
                   + f_3 * pc_y[k] * nsd_207[k];

        t_347[k] = f_13 * msd_159[k]
                   + f_3 * pc_z[k] * nsd_207[k];

        t_348[k] = f_5 * msd_167[k]
                   + f_3 * pc_y[k] * nsd_209[k];

        t_349[k] = pa_y[k] * msf0_279[k]
                   - f_4 * pc_y[k] * msf1_279[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, pc_x, pc_y, pc_z, msd_162, \
                         msd_210, msd_213, nsp0_105, nsp1_105, nsd_210, nsd_212, \
                         nsd_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_12 * msd_210[k]
                   + f_1 * nsp0_105[k]
                   - f_2 * nsp1_105[k]
                   + f_3 * pc_x[k] * nsd_210[k];

        t_351[k] = f_3 * pc_y[k] * nsd_210[k];

        t_352[k] = f_11 * msd_162[k]
                   + f_3 * pc_z[k] * nsd_210[k];

        t_353[k] = f_12 * msd_213[k]
                   + f_3 * pc_x[k] * nsd_213[k];

        t_354[k] = f_3 * pc_y[k] * nsd_212[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pc_x, pc_y, msd_215, nsp0_106, nsp0_107, \
                         nsp1_106, nsp1_107, nsd_213, nsd_214, \
                         nsd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_12 * msd_215[k]
                   + f_3 * pc_x[k] * nsd_215[k];

        t_356[k] = f_1 * nsp0_106[k]
                   - f_2 * nsp1_106[k]
                   + f_3 * pc_y[k] * nsd_213[k];

        t_357[k] = f_7 * nsp0_107[k]
                   - f_8 * nsp1_107[k]
                   + f_3 * pc_y[k] * nsd_214[k];

        t_358[k] = f_3 * pc_y[k] * nsd_215[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pc_x, pc_y, pc_z, msd_167, msd_168, \
                         msd_216, nsp0_107, nsp0_108, nsp1_107, nsp1_108, nsd_215, \
                         nsd_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_11 * msd_167[k]
                   + f_1 * nsp0_107[k]
                   - f_2 * nsp1_107[k]
                   + f_3 * pc_z[k] * nsd_215[k];

        t_360[k] = f_10 * msd_216[k]
                   + f_1 * nsp0_108[k]
                   - f_2 * nsp1_108[k]
                   + f_3 * pc_x[k] * nsd_216[k];

        t_361[k] = f_9 * msd_168[k]
                   + f_3 * pc_y[k] * nsd_216[k];

        t_362[k] = f_3 * pc_z[k] * nsd_216[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, pc_x, pc_y, pc_z, msd_171, \
                         msd_219, msd_221, nsp0_109, nsp1_109, nsd_217, nsd_219, \
                         nsd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_10 * msd_219[k]
                   + f_3 * pc_x[k] * nsd_219[k];

        t_364[k] = f_3 * pc_z[k] * nsd_217[k];

        t_365[k] = f_10 * msd_221[k]
                   + f_3 * pc_x[k] * nsd_221[k];

        t_366[k] = f_9 * msd_171[k]
                   + f_1 * nsp0_109[k]
                   - f_2 * nsp1_109[k]
                   + f_3 * pc_y[k] * nsd_219[k];

        t_367[k] = f_3 * pc_z[k] * nsd_219[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pa_z, pc_y, pc_z, msf0_280, msd_173, \
                         msd_174, msf1_280, nsp0_110, nsp1_110, nsd_221, \
                         nsd_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_9 * msd_173[k]
                   + f_3 * pc_y[k] * nsd_221[k];

        t_369[k] = f_1 * nsp0_110[k]
                   - f_2 * nsp1_110[k]
                   + f_3 * pc_z[k] * nsd_221[k];

        t_370[k] = pa_z[k] * msf0_280[k]
                   - f_4 * pc_z[k] * msf1_280[k];

        t_371[k] = f_11 * msd_174[k]
                   + f_3 * pc_y[k] * nsd_222[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, pc_x, pc_z, msd_168, msd_225, msd_226, \
                         msd_227, nsd_222, nsd_225, nsd_226, nsd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_5 * msd_168[k]
                   + f_3 * pc_z[k] * nsd_222[k];

        t_373[k] = f_10 * msd_225[k]
                   + f_3 * pc_x[k] * nsd_225[k];

        t_374[k] = f_10 * msd_226[k]
                   + f_3 * pc_x[k] * nsd_226[k];

        t_375[k] = f_10 * msd_227[k]
                   + f_3 * pc_x[k] * nsd_227[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pa_z, pc_y, pc_z, msf0_286, msd_171, \
                         msd_173, msd_179, msf1_286, nsp0_113, nsp1_113, nsd_225, \
                         nsd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pa_z[k] * msf0_286[k]
                   - f_4 * pc_z[k] * msf1_286[k];

        t_377[k] = f_5 * msd_171[k]
                   + f_3 * pc_z[k] * nsd_225[k];

        t_378[k] = f_11 * msd_179[k]
                   + f_3 * pc_y[k] * nsd_227[k];

        t_379[k] = f_5 * msd_173[k]
                   + f_1 * nsp0_113[k]
                   - f_2 * nsp1_113[k]
                   + f_3 * pc_z[k] * nsd_227[k];
    }
}

static auto
compute_prim_nsf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msf0,
                                                          const size_t msd, const size_t msf1,
                                                          const size_t nsp0, const size_t nsp1,
                                                          const size_t nsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msf0_350 = buffer.data(msf0 + 350);
    const auto *msf0_359 = buffer.data(msf0 + 359);
    const auto *msf0_360 = buffer.data(msf0 + 360);
    const auto *msf0_450 = buffer.data(msf0 + 450);
    const auto *msf0_456 = buffer.data(msf0 + 456);
    const auto *msf0_459 = buffer.data(msf0 + 459);
    const auto *msf0_466 = buffer.data(msf0 + 466);
    const auto *msf0_469 = buffer.data(msf0 + 469);
    const auto *msf0_470 = buffer.data(msf0 + 470);
    const auto *msf0_476 = buffer.data(msf0 + 476);
    const auto *msf0_479 = buffer.data(msf0 + 479);
    const auto *msf0_480 = buffer.data(msf0 + 480);
    const auto *msf0_486 = buffer.data(msf0 + 486);
    const auto *msf0_489 = buffer.data(msf0 + 489);
    const auto *msf0_490 = buffer.data(msf0 + 490);
    const auto *msf0_496 = buffer.data(msf0 + 496);
    const auto *msf0_499 = buffer.data(msf0 + 499);
    const auto *msf0_500 = buffer.data(msf0 + 500);

    const auto *msd_174 = buffer.data(msd + 174);
    const auto *msd_177 = buffer.data(msd + 177);
    const auto *msd_179 = buffer.data(msd + 179);
    const auto *msd_180 = buffer.data(msd + 180);
    const auto *msd_183 = buffer.data(msd + 183);
    const auto *msd_185 = buffer.data(msd + 185);
    const auto *msd_186 = buffer.data(msd + 186);
    const auto *msd_189 = buffer.data(msd + 189);
    const auto *msd_191 = buffer.data(msd + 191);
    const auto *msd_192 = buffer.data(msd + 192);
    const auto *msd_195 = buffer.data(msd + 195);
    const auto *msd_197 = buffer.data(msd + 197);
    const auto *msd_198 = buffer.data(msd + 198);
    const auto *msd_201 = buffer.data(msd + 201);
    const auto *msd_203 = buffer.data(msd + 203);
    const auto *msd_204 = buffer.data(msd + 204);
    const auto *msd_207 = buffer.data(msd + 207);
    const auto *msd_209 = buffer.data(msd + 209);
    const auto *msd_210 = buffer.data(msd + 210);
    const auto *msd_213 = buffer.data(msd + 213);
    const auto *msd_215 = buffer.data(msd + 215);
    const auto *msd_216 = buffer.data(msd + 216);
    const auto *msd_219 = buffer.data(msd + 219);
    const auto *msd_221 = buffer.data(msd + 221);
    const auto *msd_222 = buffer.data(msd + 222);
    const auto *msd_225 = buffer.data(msd + 225);
    const auto *msd_227 = buffer.data(msd + 227);
    const auto *msd_228 = buffer.data(msd + 228);
    const auto *msd_231 = buffer.data(msd + 231);
    const auto *msd_232 = buffer.data(msd + 232);
    const auto *msd_233 = buffer.data(msd + 233);
    const auto *msd_234 = buffer.data(msd + 234);
    const auto *msd_237 = buffer.data(msd + 237);
    const auto *msd_238 = buffer.data(msd + 238);
    const auto *msd_239 = buffer.data(msd + 239);
    const auto *msd_240 = buffer.data(msd + 240);
    const auto *msd_243 = buffer.data(msd + 243);
    const auto *msd_244 = buffer.data(msd + 244);
    const auto *msd_245 = buffer.data(msd + 245);
    const auto *msd_246 = buffer.data(msd + 246);
    const auto *msd_249 = buffer.data(msd + 249);
    const auto *msd_250 = buffer.data(msd + 250);
    const auto *msd_251 = buffer.data(msd + 251);
    const auto *msd_252 = buffer.data(msd + 252);
    const auto *msd_255 = buffer.data(msd + 255);
    const auto *msd_256 = buffer.data(msd + 256);
    const auto *msd_257 = buffer.data(msd + 257);
    const auto *msd_261 = buffer.data(msd + 261);
    const auto *msd_262 = buffer.data(msd + 262);
    const auto *msd_263 = buffer.data(msd + 263);
    const auto *msd_264 = buffer.data(msd + 264);
    const auto *msd_267 = buffer.data(msd + 267);
    const auto *msd_269 = buffer.data(msd + 269);
    const auto *msd_270 = buffer.data(msd + 270);
    const auto *msd_273 = buffer.data(msd + 273);
    const auto *msd_275 = buffer.data(msd + 275);
    const auto *msd_279 = buffer.data(msd + 279);
    const auto *msd_280 = buffer.data(msd + 280);
    const auto *msd_281 = buffer.data(msd + 281);
    const auto *msd_282 = buffer.data(msd + 282);
    const auto *msd_285 = buffer.data(msd + 285);
    const auto *msd_286 = buffer.data(msd + 286);
    const auto *msd_287 = buffer.data(msd + 287);
    const auto *msd_288 = buffer.data(msd + 288);
    const auto *msd_291 = buffer.data(msd + 291);
    const auto *msd_292 = buffer.data(msd + 292);
    const auto *msd_293 = buffer.data(msd + 293);
    const auto *msd_294 = buffer.data(msd + 294);
    const auto *msd_297 = buffer.data(msd + 297);
    const auto *msd_298 = buffer.data(msd + 298);
    const auto *msd_299 = buffer.data(msd + 299);
    const auto *msd_300 = buffer.data(msd + 300);
    const auto *msd_303 = buffer.data(msd + 303);

    const auto *msf1_350 = buffer.data(msf1 + 350);
    const auto *msf1_359 = buffer.data(msf1 + 359);
    const auto *msf1_360 = buffer.data(msf1 + 360);
    const auto *msf1_450 = buffer.data(msf1 + 450);
    const auto *msf1_456 = buffer.data(msf1 + 456);
    const auto *msf1_459 = buffer.data(msf1 + 459);
    const auto *msf1_466 = buffer.data(msf1 + 466);
    const auto *msf1_469 = buffer.data(msf1 + 469);
    const auto *msf1_470 = buffer.data(msf1 + 470);
    const auto *msf1_476 = buffer.data(msf1 + 476);
    const auto *msf1_479 = buffer.data(msf1 + 479);
    const auto *msf1_480 = buffer.data(msf1 + 480);
    const auto *msf1_486 = buffer.data(msf1 + 486);
    const auto *msf1_489 = buffer.data(msf1 + 489);
    const auto *msf1_490 = buffer.data(msf1 + 490);
    const auto *msf1_496 = buffer.data(msf1 + 496);
    const auto *msf1_499 = buffer.data(msf1 + 499);
    const auto *msf1_500 = buffer.data(msf1 + 500);

    const auto *nsp0_114 = buffer.data(nsp0 + 114);
    const auto *nsp0_115 = buffer.data(nsp0 + 115);
    const auto *nsp0_116 = buffer.data(nsp0 + 116);
    const auto *nsp0_117 = buffer.data(nsp0 + 117);
    const auto *nsp0_118 = buffer.data(nsp0 + 118);
    const auto *nsp0_119 = buffer.data(nsp0 + 119);
    const auto *nsp0_120 = buffer.data(nsp0 + 120);
    const auto *nsp0_121 = buffer.data(nsp0 + 121);
    const auto *nsp0_122 = buffer.data(nsp0 + 122);
    const auto *nsp0_123 = buffer.data(nsp0 + 123);
    const auto *nsp0_124 = buffer.data(nsp0 + 124);
    const auto *nsp0_125 = buffer.data(nsp0 + 125);
    const auto *nsp0_126 = buffer.data(nsp0 + 126);
    const auto *nsp0_127 = buffer.data(nsp0 + 127);
    const auto *nsp0_128 = buffer.data(nsp0 + 128);
    const auto *nsp0_130 = buffer.data(nsp0 + 130);
    const auto *nsp0_132 = buffer.data(nsp0 + 132);
    const auto *nsp0_133 = buffer.data(nsp0 + 133);
    const auto *nsp0_134 = buffer.data(nsp0 + 134);

    const auto *nsp1_114 = buffer.data(nsp1 + 114);
    const auto *nsp1_115 = buffer.data(nsp1 + 115);
    const auto *nsp1_116 = buffer.data(nsp1 + 116);
    const auto *nsp1_117 = buffer.data(nsp1 + 117);
    const auto *nsp1_118 = buffer.data(nsp1 + 118);
    const auto *nsp1_119 = buffer.data(nsp1 + 119);
    const auto *nsp1_120 = buffer.data(nsp1 + 120);
    const auto *nsp1_121 = buffer.data(nsp1 + 121);
    const auto *nsp1_122 = buffer.data(nsp1 + 122);
    const auto *nsp1_123 = buffer.data(nsp1 + 123);
    const auto *nsp1_124 = buffer.data(nsp1 + 124);
    const auto *nsp1_125 = buffer.data(nsp1 + 125);
    const auto *nsp1_126 = buffer.data(nsp1 + 126);
    const auto *nsp1_127 = buffer.data(nsp1 + 127);
    const auto *nsp1_128 = buffer.data(nsp1 + 128);
    const auto *nsp1_130 = buffer.data(nsp1 + 130);
    const auto *nsp1_132 = buffer.data(nsp1 + 132);
    const auto *nsp1_133 = buffer.data(nsp1 + 133);
    const auto *nsp1_134 = buffer.data(nsp1 + 134);

    const auto *nsd_228 = buffer.data(nsd + 228);
    const auto *nsd_231 = buffer.data(nsd + 231);
    const auto *nsd_232 = buffer.data(nsd + 232);
    const auto *nsd_233 = buffer.data(nsd + 233);
    const auto *nsd_234 = buffer.data(nsd + 234);
    const auto *nsd_237 = buffer.data(nsd + 237);
    const auto *nsd_238 = buffer.data(nsd + 238);
    const auto *nsd_239 = buffer.data(nsd + 239);
    const auto *nsd_240 = buffer.data(nsd + 240);
    const auto *nsd_243 = buffer.data(nsd + 243);
    const auto *nsd_244 = buffer.data(nsd + 244);
    const auto *nsd_245 = buffer.data(nsd + 245);
    const auto *nsd_246 = buffer.data(nsd + 246);
    const auto *nsd_249 = buffer.data(nsd + 249);
    const auto *nsd_250 = buffer.data(nsd + 250);
    const auto *nsd_251 = buffer.data(nsd + 251);
    const auto *nsd_252 = buffer.data(nsd + 252);
    const auto *nsd_255 = buffer.data(nsd + 255);
    const auto *nsd_256 = buffer.data(nsd + 256);
    const auto *nsd_257 = buffer.data(nsd + 257);
    const auto *nsd_258 = buffer.data(nsd + 258);
    const auto *nsd_261 = buffer.data(nsd + 261);
    const auto *nsd_262 = buffer.data(nsd + 262);
    const auto *nsd_263 = buffer.data(nsd + 263);
    const auto *nsd_264 = buffer.data(nsd + 264);
    const auto *nsd_266 = buffer.data(nsd + 266);
    const auto *nsd_267 = buffer.data(nsd + 267);
    const auto *nsd_268 = buffer.data(nsd + 268);
    const auto *nsd_269 = buffer.data(nsd + 269);
    const auto *nsd_270 = buffer.data(nsd + 270);
    const auto *nsd_271 = buffer.data(nsd + 271);
    const auto *nsd_273 = buffer.data(nsd + 273);
    const auto *nsd_275 = buffer.data(nsd + 275);
    const auto *nsd_276 = buffer.data(nsd + 276);
    const auto *nsd_279 = buffer.data(nsd + 279);
    const auto *nsd_280 = buffer.data(nsd + 280);
    const auto *nsd_281 = buffer.data(nsd + 281);
    const auto *nsd_282 = buffer.data(nsd + 282);
    const auto *nsd_285 = buffer.data(nsd + 285);
    const auto *nsd_286 = buffer.data(nsd + 286);
    const auto *nsd_287 = buffer.data(nsd + 287);
    const auto *nsd_288 = buffer.data(nsd + 288);
    const auto *nsd_291 = buffer.data(nsd + 291);
    const auto *nsd_292 = buffer.data(nsd + 292);
    const auto *nsd_293 = buffer.data(nsd + 293);
    const auto *nsd_294 = buffer.data(nsd + 294);
    const auto *nsd_297 = buffer.data(nsd + 297);
    const auto *nsd_298 = buffer.data(nsd + 298);
    const auto *nsd_299 = buffer.data(nsd + 299);
    const auto *nsd_300 = buffer.data(nsd + 300);
    const auto *nsd_303 = buffer.data(nsd + 303);

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, pc_y, pc_z, msd_174, msd_180, \
                         msd_228, msd_231, nsp0_114, nsp1_114, nsd_228, \
                         nsd_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_10 * msd_228[k]
                   + f_1 * nsp0_114[k]
                   - f_2 * nsp1_114[k]
                   + f_3 * pc_x[k] * nsd_228[k];

        t_381[k] = f_13 * msd_180[k]
                   + f_3 * pc_y[k] * nsd_228[k];

        t_382[k] = f_10 * msd_174[k]
                   + f_3 * pc_z[k] * nsd_228[k];

        t_383[k] = f_10 * msd_231[k]
                   + f_3 * pc_x[k] * nsd_231[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pc_x, pc_y, pc_z, msd_177, msd_183, \
                         msd_232, msd_233, nsp0_115, nsp1_115, nsd_231, nsd_232, \
                         nsd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_10 * msd_232[k]
                   + f_3 * pc_x[k] * nsd_232[k];

        t_385[k] = f_10 * msd_233[k]
                   + f_3 * pc_x[k] * nsd_233[k];

        t_386[k] = f_13 * msd_183[k]
                   + f_1 * nsp0_115[k]
                   - f_2 * nsp1_115[k]
                   + f_3 * pc_y[k] * nsd_231[k];

        t_387[k] = f_10 * msd_177[k]
                   + f_3 * pc_z[k] * nsd_231[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_x, pc_y, pc_z, msd_179, msd_185, msd_234, \
                         nsp0_116, nsp0_117, nsp1_116, nsp1_117, nsd_233, \
                         nsd_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_13 * msd_185[k]
                   + f_3 * pc_y[k] * nsd_233[k];

        t_389[k] = f_10 * msd_179[k]
                   + f_1 * nsp0_116[k]
                   - f_2 * nsp1_116[k]
                   + f_3 * pc_z[k] * nsd_233[k];

        t_390[k] = f_10 * msd_234[k]
                   + f_1 * nsp0_117[k]
                   - f_2 * nsp1_117[k]
                   + f_3 * pc_x[k] * nsd_234[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pc_x, pc_y, pc_z, msd_180, msd_186, \
                         msd_237, msd_238, nsd_234, nsd_237, nsd_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_15 * msd_186[k]
                   + f_3 * pc_y[k] * nsd_234[k];

        t_392[k] = f_12 * msd_180[k]
                   + f_3 * pc_z[k] * nsd_234[k];

        t_393[k] = f_10 * msd_237[k]
                   + f_3 * pc_x[k] * nsd_237[k];

        t_394[k] = f_10 * msd_238[k]
                   + f_3 * pc_x[k] * nsd_238[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pc_x, pc_y, pc_z, msd_183, msd_189, \
                         msd_191, msd_239, nsp0_118, nsp1_118, nsd_237, \
                         nsd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_10 * msd_239[k]
                   + f_3 * pc_x[k] * nsd_239[k];

        t_396[k] = f_15 * msd_189[k]
                   + f_1 * nsp0_118[k]
                   - f_2 * nsp1_118[k]
                   + f_3 * pc_y[k] * nsd_237[k];

        t_397[k] = f_12 * msd_183[k]
                   + f_3 * pc_z[k] * nsd_237[k];

        t_398[k] = f_15 * msd_191[k]
                   + f_3 * pc_y[k] * nsd_239[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, msd_185, msd_192, msd_240, \
                         nsp0_119, nsp0_120, nsp1_119, nsp1_120, nsd_239, \
                         nsd_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_12 * msd_185[k]
                   + f_1 * nsp0_119[k]
                   - f_2 * nsp1_119[k]
                   + f_3 * pc_z[k] * nsd_239[k];

        t_400[k] = f_10 * msd_240[k]
                   + f_1 * nsp0_120[k]
                   - f_2 * nsp1_120[k]
                   + f_3 * pc_x[k] * nsd_240[k];

        t_401[k] = f_14 * msd_192[k]
                   + f_3 * pc_y[k] * nsd_240[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, pc_x, pc_z, msd_186, msd_243, msd_244, \
                         msd_245, nsd_240, nsd_243, nsd_244, nsd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_14 * msd_186[k]
                   + f_3 * pc_z[k] * nsd_240[k];

        t_403[k] = f_10 * msd_243[k]
                   + f_3 * pc_x[k] * nsd_243[k];

        t_404[k] = f_10 * msd_244[k]
                   + f_3 * pc_x[k] * nsd_244[k];

        t_405[k] = f_10 * msd_245[k]
                   + f_3 * pc_x[k] * nsd_245[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pc_y, pc_z, msd_189, msd_191, msd_195, \
                         msd_197, nsp0_121, nsp0_122, nsp1_121, nsp1_122, nsd_243, \
                         nsd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = f_14 * msd_195[k]
                   + f_1 * nsp0_121[k]
                   - f_2 * nsp1_121[k]
                   + f_3 * pc_y[k] * nsd_243[k];

        t_407[k] = f_14 * msd_189[k]
                   + f_3 * pc_z[k] * nsd_243[k];

        t_408[k] = f_14 * msd_197[k]
                   + f_3 * pc_y[k] * nsd_245[k];

        t_409[k] = f_14 * msd_191[k]
                   + f_1 * nsp0_122[k]
                   - f_2 * nsp1_122[k]
                   + f_3 * pc_z[k] * nsd_245[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, pc_z, msd_192, msd_198, \
                         msd_246, msd_249, nsp0_123, nsp1_123, nsd_246, \
                         nsd_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_10 * msd_246[k]
                   + f_1 * nsp0_123[k]
                   - f_2 * nsp1_123[k]
                   + f_3 * pc_x[k] * nsd_246[k];

        t_411[k] = f_12 * msd_198[k]
                   + f_3 * pc_y[k] * nsd_246[k];

        t_412[k] = f_15 * msd_192[k]
                   + f_3 * pc_z[k] * nsd_246[k];

        t_413[k] = f_10 * msd_249[k]
                   + f_3 * pc_x[k] * nsd_249[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_x, pc_y, pc_z, msd_195, msd_201, \
                         msd_250, msd_251, nsp0_124, nsp1_124, nsd_249, nsd_250, \
                         nsd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_10 * msd_250[k]
                   + f_3 * pc_x[k] * nsd_250[k];

        t_415[k] = f_10 * msd_251[k]
                   + f_3 * pc_x[k] * nsd_251[k];

        t_416[k] = f_12 * msd_201[k]
                   + f_1 * nsp0_124[k]
                   - f_2 * nsp1_124[k]
                   + f_3 * pc_y[k] * nsd_249[k];

        t_417[k] = f_15 * msd_195[k]
                   + f_3 * pc_z[k] * nsd_249[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, pc_x, pc_y, pc_z, msd_197, msd_203, msd_252, \
                         nsp0_125, nsp0_126, nsp1_125, nsp1_126, nsd_251, \
                         nsd_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_12 * msd_203[k]
                   + f_3 * pc_y[k] * nsd_251[k];

        t_419[k] = f_15 * msd_197[k]
                   + f_1 * nsp0_125[k]
                   - f_2 * nsp1_125[k]
                   + f_3 * pc_z[k] * nsd_251[k];

        t_420[k] = f_10 * msd_252[k]
                   + f_1 * nsp0_126[k]
                   - f_2 * nsp1_126[k]
                   + f_3 * pc_x[k] * nsd_252[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, msd_198, msd_204, \
                         msd_255, msd_256, nsd_252, nsd_255, nsd_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_10 * msd_204[k]
                   + f_3 * pc_y[k] * nsd_252[k];

        t_422[k] = f_13 * msd_198[k]
                   + f_3 * pc_z[k] * nsd_252[k];

        t_423[k] = f_10 * msd_255[k]
                   + f_3 * pc_x[k] * nsd_255[k];

        t_424[k] = f_10 * msd_256[k]
                   + f_3 * pc_x[k] * nsd_256[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, pc_z, msd_201, msd_207, \
                         msd_209, msd_257, nsp0_127, nsp1_127, nsd_255, \
                         nsd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_10 * msd_257[k]
                   + f_3 * pc_x[k] * nsd_257[k];

        t_426[k] = f_10 * msd_207[k]
                   + f_1 * nsp0_127[k]
                   - f_2 * nsp1_127[k]
                   + f_3 * pc_y[k] * nsd_255[k];

        t_427[k] = f_13 * msd_201[k]
                   + f_3 * pc_z[k] * nsd_255[k];

        t_428[k] = f_10 * msd_209[k]
                   + f_3 * pc_y[k] * nsd_257[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pa_y, pc_y, pc_z, msf0_350, msd_203, \
                         msd_204, msd_210, msf1_350, nsp0_128, nsp1_128, nsd_257, \
                         nsd_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_13 * msd_203[k]
                   + f_1 * nsp0_128[k]
                   - f_2 * nsp1_128[k]
                   + f_3 * pc_z[k] * nsd_257[k];

        t_430[k] = pa_y[k] * msf0_350[k]
                   - f_4 * pc_y[k] * msf1_350[k];

        t_431[k] = f_5 * msd_210[k]
                   + f_3 * pc_y[k] * nsd_258[k];

        t_432[k] = f_11 * msd_204[k]
                   + f_3 * pc_z[k] * nsd_258[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, msd_213, msd_261, msd_262, \
                         msd_263, nsp0_130, nsp1_130, nsd_261, nsd_262, \
                         nsd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_10 * msd_261[k]
                   + f_3 * pc_x[k] * nsd_261[k];

        t_434[k] = f_10 * msd_262[k]
                   + f_3 * pc_x[k] * nsd_262[k];

        t_435[k] = f_10 * msd_263[k]
                   + f_3 * pc_x[k] * nsd_263[k];

        t_436[k] = f_5 * msd_213[k]
                   + f_1 * nsp0_130[k]
                   - f_2 * nsp1_130[k]
                   + f_3 * pc_y[k] * nsd_261[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pa_y, pc_y, pc_z, msf0_359, msd_207, msd_215, \
                         msf1_359, nsd_261, nsd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_11 * msd_207[k]
                   + f_3 * pc_z[k] * nsd_261[k];

        t_438[k] = f_5 * msd_215[k]
                   + f_3 * pc_y[k] * nsd_263[k];

        t_439[k] = pa_y[k] * msf0_359[k]
                   - f_4 * pc_y[k] * msf1_359[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pc_x, pc_y, pc_z, msd_210, \
                         msd_264, msd_267, nsp0_132, nsp1_132, nsd_264, nsd_266, \
                         nsd_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_10 * msd_264[k]
                   + f_1 * nsp0_132[k]
                   - f_2 * nsp1_132[k]
                   + f_3 * pc_x[k] * nsd_264[k];

        t_441[k] = f_3 * pc_y[k] * nsd_264[k];

        t_442[k] = f_9 * msd_210[k]
                   + f_3 * pc_z[k] * nsd_264[k];

        t_443[k] = f_10 * msd_267[k]
                   + f_3 * pc_x[k] * nsd_267[k];

        t_444[k] = f_3 * pc_y[k] * nsd_266[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_x, pc_y, msd_269, nsp0_133, nsp0_134, \
                         nsp1_133, nsp1_134, nsd_267, nsd_268, \
                         nsd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_10 * msd_269[k]
                   + f_3 * pc_x[k] * nsd_269[k];

        t_446[k] = f_1 * nsp0_133[k]
                   - f_2 * nsp1_133[k]
                   + f_3 * pc_y[k] * nsd_267[k];

        t_447[k] = f_7 * nsp0_134[k]
                   - f_8 * nsp1_134[k]
                   + f_3 * pc_y[k] * nsd_268[k];

        t_448[k] = f_3 * pc_y[k] * nsd_269[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pa_x, pc_x, pc_y, pc_z, msf0_450, msd_215, \
                         msd_216, msd_270, msf1_450, nsp0_134, nsp1_134, nsd_269, \
                         nsd_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_9 * msd_215[k]
                   + f_1 * nsp0_134[k]
                   - f_2 * nsp1_134[k]
                   + f_3 * pc_z[k] * nsd_269[k];

        t_450[k] = pa_x[k] * msf0_450[k]
                   + f_12 * msd_270[k]
                   - f_4 * pc_x[k] * msf1_450[k];

        t_451[k] = f_6 * msd_216[k]
                   + f_3 * pc_y[k] * nsd_270[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, t_456, pa_x, pc_x, pc_z, msf0_456, \
                         msd_273, msd_275, msf1_456, nsd_270, nsd_271, nsd_273, \
                         nsd_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_3 * pc_z[k] * nsd_270[k];

        t_453[k] = f_5 * msd_273[k]
                   + f_3 * pc_x[k] * nsd_273[k];

        t_454[k] = f_3 * pc_z[k] * nsd_271[k];

        t_455[k] = f_5 * msd_275[k]
                   + f_3 * pc_x[k] * nsd_275[k];

        t_456[k] = pa_x[k] * msf0_456[k]
                   - f_4 * pc_x[k] * msf1_456[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pa_x, pa_z, pc_x, pc_y, pc_z, msf0_360, \
                         msf0_459, msd_221, msf1_360, msf1_459, nsd_273, \
                         nsd_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * nsd_273[k];

        t_458[k] = f_6 * msd_221[k]
                   + f_3 * pc_y[k] * nsd_275[k];

        t_459[k] = pa_x[k] * msf0_459[k]
                   - f_4 * pc_x[k] * msf1_459[k];

        t_460[k] = pa_z[k] * msf0_360[k]
                   - f_4 * pc_z[k] * msf1_360[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pc_x, pc_y, pc_z, msd_216, msd_222, \
                         msd_279, msd_280, nsd_276, nsd_279, nsd_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_9 * msd_222[k]
                   + f_3 * pc_y[k] * nsd_276[k];

        t_462[k] = f_5 * msd_216[k]
                   + f_3 * pc_z[k] * nsd_276[k];

        t_463[k] = f_5 * msd_279[k]
                   + f_3 * pc_x[k] * nsd_279[k];

        t_464[k] = f_5 * msd_280[k]
                   + f_3 * pc_x[k] * nsd_280[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, pa_x, pc_x, pc_y, pc_z, msf0_466, \
                         msd_219, msd_227, msd_281, msf1_466, nsd_279, \
                         nsd_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_5 * msd_281[k]
                   + f_3 * pc_x[k] * nsd_281[k];

        t_466[k] = pa_x[k] * msf0_466[k]
                   - f_4 * pc_x[k] * msf1_466[k];

        t_467[k] = f_5 * msd_219[k]
                   + f_3 * pc_z[k] * nsd_279[k];

        t_468[k] = f_9 * msd_227[k]
                   + f_3 * pc_y[k] * nsd_281[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, pa_x, pc_x, pc_y, pc_z, msf0_469, \
                         msf0_470, msd_222, msd_228, msd_282, msf1_469, msf1_470, \
                         nsd_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = pa_x[k] * msf0_469[k]
                   - f_4 * pc_x[k] * msf1_469[k];

        t_470[k] = pa_x[k] * msf0_470[k]
                   + f_12 * msd_282[k]
                   - f_4 * pc_x[k] * msf1_470[k];

        t_471[k] = f_11 * msd_228[k]
                   + f_3 * pc_y[k] * nsd_282[k];

        t_472[k] = f_10 * msd_222[k]
                   + f_3 * pc_z[k] * nsd_282[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, pa_x, pc_x, msf0_476, msd_285, msd_286, \
                         msd_287, msf1_476, nsd_285, nsd_286, nsd_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_5 * msd_285[k]
                   + f_3 * pc_x[k] * nsd_285[k];

        t_474[k] = f_5 * msd_286[k]
                   + f_3 * pc_x[k] * nsd_286[k];

        t_475[k] = f_5 * msd_287[k]
                   + f_3 * pc_x[k] * nsd_287[k];

        t_476[k] = pa_x[k] * msf0_476[k]
                   - f_4 * pc_x[k] * msf1_476[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, pa_x, pc_x, pc_y, pc_z, msf0_479, msd_225, \
                         msd_233, msf1_479, nsd_285, nsd_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_10 * msd_225[k]
                   + f_3 * pc_z[k] * nsd_285[k];

        t_478[k] = f_11 * msd_233[k]
                   + f_3 * pc_y[k] * nsd_287[k];

        t_479[k] = pa_x[k] * msf0_479[k]
                   - f_4 * pc_x[k] * msf1_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, pa_x, pc_x, pc_y, pc_z, msf0_480, \
                         msd_228, msd_234, msd_288, msd_291, msf1_480, nsd_288, \
                         nsd_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = pa_x[k] * msf0_480[k]
                   + f_12 * msd_288[k]
                   - f_4 * pc_x[k] * msf1_480[k];

        t_481[k] = f_13 * msd_234[k]
                   + f_3 * pc_y[k] * nsd_288[k];

        t_482[k] = f_12 * msd_228[k]
                   + f_3 * pc_z[k] * nsd_288[k];

        t_483[k] = f_5 * msd_291[k]
                   + f_3 * pc_x[k] * nsd_291[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pa_x, pc_x, pc_z, msf0_486, msd_231, \
                         msd_292, msd_293, msf1_486, nsd_291, nsd_292, \
                         nsd_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_5 * msd_292[k]
                   + f_3 * pc_x[k] * nsd_292[k];

        t_485[k] = f_5 * msd_293[k]
                   + f_3 * pc_x[k] * nsd_293[k];

        t_486[k] = pa_x[k] * msf0_486[k]
                   - f_4 * pc_x[k] * msf1_486[k];

        t_487[k] = f_12 * msd_231[k]
                   + f_3 * pc_z[k] * nsd_291[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pa_x, pc_x, pc_y, msf0_489, msf0_490, \
                         msd_239, msd_240, msd_294, msf1_489, msf1_490, nsd_293, \
                         nsd_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_13 * msd_239[k]
                   + f_3 * pc_y[k] * nsd_293[k];

        t_489[k] = pa_x[k] * msf0_489[k]
                   - f_4 * pc_x[k] * msf1_489[k];

        t_490[k] = pa_x[k] * msf0_490[k]
                   + f_12 * msd_294[k]
                   - f_4 * pc_x[k] * msf1_490[k];

        t_491[k] = f_15 * msd_240[k]
                   + f_3 * pc_y[k] * nsd_294[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pc_x, pc_z, msd_234, msd_297, msd_298, \
                         msd_299, nsd_294, nsd_297, nsd_298, nsd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_14 * msd_234[k]
                   + f_3 * pc_z[k] * nsd_294[k];

        t_493[k] = f_5 * msd_297[k]
                   + f_3 * pc_x[k] * nsd_297[k];

        t_494[k] = f_5 * msd_298[k]
                   + f_3 * pc_x[k] * nsd_298[k];

        t_495[k] = f_5 * msd_299[k]
                   + f_3 * pc_x[k] * nsd_299[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pa_x, pc_x, pc_y, pc_z, msf0_496, \
                         msf0_499, msd_237, msd_245, msf1_496, msf1_499, nsd_297, \
                         nsd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = pa_x[k] * msf0_496[k]
                   - f_4 * pc_x[k] * msf1_496[k];

        t_497[k] = f_14 * msd_237[k]
                   + f_3 * pc_z[k] * nsd_297[k];

        t_498[k] = f_15 * msd_245[k]
                   + f_3 * pc_y[k] * nsd_299[k];

        t_499[k] = pa_x[k] * msf0_499[k]
                   - f_4 * pc_x[k] * msf1_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_x, pc_x, pc_y, pc_z, msf0_500, \
                         msd_240, msd_246, msd_300, msd_303, msf1_500, nsd_300, \
                         nsd_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pa_x[k] * msf0_500[k]
                   + f_12 * msd_300[k]
                   - f_4 * pc_x[k] * msf1_500[k];

        t_501[k] = f_14 * msd_246[k]
                   + f_3 * pc_y[k] * nsd_300[k];

        t_502[k] = f_15 * msd_240[k]
                   + f_3 * pc_z[k] * nsd_300[k];

        t_503[k] = f_5 * msd_303[k]
                   + f_3 * pc_x[k] * nsd_303[k];
    }
}

static auto
compute_prim_nsf_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msf0,
                                                          const size_t msd, const size_t msf1,
                                                          const size_t nsp0, const size_t nsp1,
                                                          const size_t nsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msf0_440 = buffer.data(msf0 + 440);
    const auto *msf0_450 = buffer.data(msf0 + 450);
    const auto *msf0_451 = buffer.data(msf0 + 451);
    const auto *msf0_456 = buffer.data(msf0 + 456);
    const auto *msf0_506 = buffer.data(msf0 + 506);
    const auto *msf0_509 = buffer.data(msf0 + 509);
    const auto *msf0_510 = buffer.data(msf0 + 510);
    const auto *msf0_516 = buffer.data(msf0 + 516);
    const auto *msf0_519 = buffer.data(msf0 + 519);
    const auto *msf0_520 = buffer.data(msf0 + 520);
    const auto *msf0_526 = buffer.data(msf0 + 526);
    const auto *msf0_529 = buffer.data(msf0 + 529);
    const auto *msf0_536 = buffer.data(msf0 + 536);
    const auto *msf0_539 = buffer.data(msf0 + 539);
    const auto *msf0_540 = buffer.data(msf0 + 540);
    const auto *msf0_546 = buffer.data(msf0 + 546);
    const auto *msf0_547 = buffer.data(msf0 + 547);
    const auto *msf0_549 = buffer.data(msf0 + 549);

    const auto *msd_243 = buffer.data(msd + 243);
    const auto *msd_246 = buffer.data(msd + 246);
    const auto *msd_249 = buffer.data(msd + 249);
    const auto *msd_251 = buffer.data(msd + 251);
    const auto *msd_252 = buffer.data(msd + 252);
    const auto *msd_255 = buffer.data(msd + 255);
    const auto *msd_257 = buffer.data(msd + 257);
    const auto *msd_258 = buffer.data(msd + 258);
    const auto *msd_261 = buffer.data(msd + 261);
    const auto *msd_263 = buffer.data(msd + 263);
    const auto *msd_264 = buffer.data(msd + 264);
    const auto *msd_269 = buffer.data(msd + 269);
    const auto *msd_273 = buffer.data(msd + 273);
    const auto *msd_275 = buffer.data(msd + 275);
    const auto *msd_279 = buffer.data(msd + 279);
    const auto *msd_281 = buffer.data(msd + 281);
    const auto *msd_285 = buffer.data(msd + 285);
    const auto *msd_287 = buffer.data(msd + 287);
    const auto *msd_291 = buffer.data(msd + 291);
    const auto *msd_293 = buffer.data(msd + 293);
    const auto *msd_297 = buffer.data(msd + 297);
    const auto *msd_299 = buffer.data(msd + 299);
    const auto *msd_303 = buffer.data(msd + 303);
    const auto *msd_304 = buffer.data(msd + 304);
    const auto *msd_305 = buffer.data(msd + 305);
    const auto *msd_306 = buffer.data(msd + 306);
    const auto *msd_309 = buffer.data(msd + 309);
    const auto *msd_310 = buffer.data(msd + 310);
    const auto *msd_311 = buffer.data(msd + 311);
    const auto *msd_312 = buffer.data(msd + 312);
    const auto *msd_315 = buffer.data(msd + 315);
    const auto *msd_316 = buffer.data(msd + 316);
    const auto *msd_317 = buffer.data(msd + 317);
    const auto *msd_321 = buffer.data(msd + 321);
    const auto *msd_322 = buffer.data(msd + 322);
    const auto *msd_323 = buffer.data(msd + 323);
    const auto *msd_324 = buffer.data(msd + 324);
    const auto *msd_327 = buffer.data(msd + 327);
    const auto *msd_329 = buffer.data(msd + 329);

    const auto *msf1_440 = buffer.data(msf1 + 440);
    const auto *msf1_450 = buffer.data(msf1 + 450);
    const auto *msf1_451 = buffer.data(msf1 + 451);
    const auto *msf1_456 = buffer.data(msf1 + 456);
    const auto *msf1_506 = buffer.data(msf1 + 506);
    const auto *msf1_509 = buffer.data(msf1 + 509);
    const auto *msf1_510 = buffer.data(msf1 + 510);
    const auto *msf1_516 = buffer.data(msf1 + 516);
    const auto *msf1_519 = buffer.data(msf1 + 519);
    const auto *msf1_520 = buffer.data(msf1 + 520);
    const auto *msf1_526 = buffer.data(msf1 + 526);
    const auto *msf1_529 = buffer.data(msf1 + 529);
    const auto *msf1_536 = buffer.data(msf1 + 536);
    const auto *msf1_539 = buffer.data(msf1 + 539);
    const auto *msf1_540 = buffer.data(msf1 + 540);
    const auto *msf1_546 = buffer.data(msf1 + 546);
    const auto *msf1_547 = buffer.data(msf1 + 547);
    const auto *msf1_549 = buffer.data(msf1 + 549);

    const auto *nsp0_165 = buffer.data(nsp0 + 165);
    const auto *nsp0_166 = buffer.data(nsp0 + 166);
    const auto *nsp0_167 = buffer.data(nsp0 + 167);
    const auto *nsp0_170 = buffer.data(nsp0 + 170);
    const auto *nsp0_171 = buffer.data(nsp0 + 171);
    const auto *nsp0_172 = buffer.data(nsp0 + 172);
    const auto *nsp0_173 = buffer.data(nsp0 + 173);
    const auto *nsp0_174 = buffer.data(nsp0 + 174);
    const auto *nsp0_175 = buffer.data(nsp0 + 175);
    const auto *nsp0_176 = buffer.data(nsp0 + 176);
    const auto *nsp0_177 = buffer.data(nsp0 + 177);
    const auto *nsp0_178 = buffer.data(nsp0 + 178);
    const auto *nsp0_179 = buffer.data(nsp0 + 179);
    const auto *nsp0_180 = buffer.data(nsp0 + 180);
    const auto *nsp0_181 = buffer.data(nsp0 + 181);
    const auto *nsp0_182 = buffer.data(nsp0 + 182);
    const auto *nsp0_183 = buffer.data(nsp0 + 183);
    const auto *nsp0_184 = buffer.data(nsp0 + 184);
    const auto *nsp0_185 = buffer.data(nsp0 + 185);
    const auto *nsp0_186 = buffer.data(nsp0 + 186);
    const auto *nsp0_187 = buffer.data(nsp0 + 187);
    const auto *nsp0_188 = buffer.data(nsp0 + 188);
    const auto *nsp0_189 = buffer.data(nsp0 + 189);
    const auto *nsp0_190 = buffer.data(nsp0 + 190);
    const auto *nsp0_191 = buffer.data(nsp0 + 191);

    const auto *nsp1_165 = buffer.data(nsp1 + 165);
    const auto *nsp1_166 = buffer.data(nsp1 + 166);
    const auto *nsp1_167 = buffer.data(nsp1 + 167);
    const auto *nsp1_170 = buffer.data(nsp1 + 170);
    const auto *nsp1_171 = buffer.data(nsp1 + 171);
    const auto *nsp1_172 = buffer.data(nsp1 + 172);
    const auto *nsp1_173 = buffer.data(nsp1 + 173);
    const auto *nsp1_174 = buffer.data(nsp1 + 174);
    const auto *nsp1_175 = buffer.data(nsp1 + 175);
    const auto *nsp1_176 = buffer.data(nsp1 + 176);
    const auto *nsp1_177 = buffer.data(nsp1 + 177);
    const auto *nsp1_178 = buffer.data(nsp1 + 178);
    const auto *nsp1_179 = buffer.data(nsp1 + 179);
    const auto *nsp1_180 = buffer.data(nsp1 + 180);
    const auto *nsp1_181 = buffer.data(nsp1 + 181);
    const auto *nsp1_182 = buffer.data(nsp1 + 182);
    const auto *nsp1_183 = buffer.data(nsp1 + 183);
    const auto *nsp1_184 = buffer.data(nsp1 + 184);
    const auto *nsp1_185 = buffer.data(nsp1 + 185);
    const auto *nsp1_186 = buffer.data(nsp1 + 186);
    const auto *nsp1_187 = buffer.data(nsp1 + 187);
    const auto *nsp1_188 = buffer.data(nsp1 + 188);
    const auto *nsp1_189 = buffer.data(nsp1 + 189);
    const auto *nsp1_190 = buffer.data(nsp1 + 190);
    const auto *nsp1_191 = buffer.data(nsp1 + 191);

    const auto *nsd_303 = buffer.data(nsd + 303);
    const auto *nsd_304 = buffer.data(nsd + 304);
    const auto *nsd_305 = buffer.data(nsd + 305);
    const auto *nsd_306 = buffer.data(nsd + 306);
    const auto *nsd_309 = buffer.data(nsd + 309);
    const auto *nsd_310 = buffer.data(nsd + 310);
    const auto *nsd_311 = buffer.data(nsd + 311);
    const auto *nsd_312 = buffer.data(nsd + 312);
    const auto *nsd_315 = buffer.data(nsd + 315);
    const auto *nsd_316 = buffer.data(nsd + 316);
    const auto *nsd_317 = buffer.data(nsd + 317);
    const auto *nsd_318 = buffer.data(nsd + 318);
    const auto *nsd_321 = buffer.data(nsd + 321);
    const auto *nsd_322 = buffer.data(nsd + 322);
    const auto *nsd_323 = buffer.data(nsd + 323);
    const auto *nsd_324 = buffer.data(nsd + 324);
    const auto *nsd_326 = buffer.data(nsd + 326);
    const auto *nsd_327 = buffer.data(nsd + 327);
    const auto *nsd_329 = buffer.data(nsd + 329);
    const auto *nsd_330 = buffer.data(nsd + 330);
    const auto *nsd_331 = buffer.data(nsd + 331);
    const auto *nsd_333 = buffer.data(nsd + 333);
    const auto *nsd_334 = buffer.data(nsd + 334);
    const auto *nsd_335 = buffer.data(nsd + 335);
    const auto *nsd_338 = buffer.data(nsd + 338);
    const auto *nsd_339 = buffer.data(nsd + 339);
    const auto *nsd_340 = buffer.data(nsd + 340);
    const auto *nsd_341 = buffer.data(nsd + 341);
    const auto *nsd_342 = buffer.data(nsd + 342);
    const auto *nsd_343 = buffer.data(nsd + 343);
    const auto *nsd_344 = buffer.data(nsd + 344);
    const auto *nsd_345 = buffer.data(nsd + 345);
    const auto *nsd_346 = buffer.data(nsd + 346);
    const auto *nsd_347 = buffer.data(nsd + 347);
    const auto *nsd_348 = buffer.data(nsd + 348);
    const auto *nsd_349 = buffer.data(nsd + 349);
    const auto *nsd_350 = buffer.data(nsd + 350);
    const auto *nsd_351 = buffer.data(nsd + 351);
    const auto *nsd_352 = buffer.data(nsd + 352);
    const auto *nsd_353 = buffer.data(nsd + 353);
    const auto *nsd_354 = buffer.data(nsd + 354);
    const auto *nsd_355 = buffer.data(nsd + 355);
    const auto *nsd_356 = buffer.data(nsd + 356);
    const auto *nsd_357 = buffer.data(nsd + 357);
    const auto *nsd_358 = buffer.data(nsd + 358);
    const auto *nsd_359 = buffer.data(nsd + 359);
    const auto *nsd_360 = buffer.data(nsd + 360);
    const auto *nsd_361 = buffer.data(nsd + 361);
    const auto *nsd_362 = buffer.data(nsd + 362);
    const auto *nsd_363 = buffer.data(nsd + 363);
    const auto *nsd_364 = buffer.data(nsd + 364);
    const auto *nsd_365 = buffer.data(nsd + 365);
    const auto *nsd_366 = buffer.data(nsd + 366);
    const auto *nsd_367 = buffer.data(nsd + 367);
    const auto *nsd_368 = buffer.data(nsd + 368);
    const auto *nsd_369 = buffer.data(nsd + 369);
    const auto *nsd_370 = buffer.data(nsd + 370);
    const auto *nsd_371 = buffer.data(nsd + 371);
    const auto *nsd_372 = buffer.data(nsd + 372);
    const auto *nsd_373 = buffer.data(nsd + 373);
    const auto *nsd_374 = buffer.data(nsd + 374);
    const auto *nsd_375 = buffer.data(nsd + 375);
    const auto *nsd_376 = buffer.data(nsd + 376);
    const auto *nsd_377 = buffer.data(nsd + 377);
    const auto *nsd_378 = buffer.data(nsd + 378);
    const auto *nsd_379 = buffer.data(nsd + 379);
    const auto *nsd_380 = buffer.data(nsd + 380);
    const auto *nsd_381 = buffer.data(nsd + 381);

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pa_x, pc_x, pc_z, msf0_506, msd_243, \
                         msd_304, msd_305, msf1_506, nsd_303, nsd_304, \
                         nsd_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_5 * msd_304[k]
                   + f_3 * pc_x[k] * nsd_304[k];

        t_505[k] = f_5 * msd_305[k]
                   + f_3 * pc_x[k] * nsd_305[k];

        t_506[k] = pa_x[k] * msf0_506[k]
                   - f_4 * pc_x[k] * msf1_506[k];

        t_507[k] = f_15 * msd_243[k]
                   + f_3 * pc_z[k] * nsd_303[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_x, pc_x, pc_y, msf0_509, msf0_510, \
                         msd_251, msd_252, msd_306, msf1_509, msf1_510, nsd_305, \
                         nsd_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_14 * msd_251[k]
                   + f_3 * pc_y[k] * nsd_305[k];

        t_509[k] = pa_x[k] * msf0_509[k]
                   - f_4 * pc_x[k] * msf1_509[k];

        t_510[k] = pa_x[k] * msf0_510[k]
                   + f_12 * msd_306[k]
                   - f_4 * pc_x[k] * msf1_510[k];

        t_511[k] = f_12 * msd_252[k]
                   + f_3 * pc_y[k] * nsd_306[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_z, msd_246, msd_309, msd_310, \
                         msd_311, nsd_306, nsd_309, nsd_310, nsd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_13 * msd_246[k]
                   + f_3 * pc_z[k] * nsd_306[k];

        t_513[k] = f_5 * msd_309[k]
                   + f_3 * pc_x[k] * nsd_309[k];

        t_514[k] = f_5 * msd_310[k]
                   + f_3 * pc_x[k] * nsd_310[k];

        t_515[k] = f_5 * msd_311[k]
                   + f_3 * pc_x[k] * nsd_311[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pa_x, pc_x, pc_y, pc_z, msf0_516, \
                         msf0_519, msd_249, msd_257, msf1_516, msf1_519, nsd_309, \
                         nsd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = pa_x[k] * msf0_516[k]
                   - f_4 * pc_x[k] * msf1_516[k];

        t_517[k] = f_13 * msd_249[k]
                   + f_3 * pc_z[k] * nsd_309[k];

        t_518[k] = f_12 * msd_257[k]
                   + f_3 * pc_y[k] * nsd_311[k];

        t_519[k] = pa_x[k] * msf0_519[k]
                   - f_4 * pc_x[k] * msf1_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, pa_x, pc_x, pc_y, pc_z, msf0_520, \
                         msd_252, msd_258, msd_312, msd_315, msf1_520, nsd_312, \
                         nsd_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = pa_x[k] * msf0_520[k]
                   + f_12 * msd_312[k]
                   - f_4 * pc_x[k] * msf1_520[k];

        t_521[k] = f_10 * msd_258[k]
                   + f_3 * pc_y[k] * nsd_312[k];

        t_522[k] = f_11 * msd_252[k]
                   + f_3 * pc_z[k] * nsd_312[k];

        t_523[k] = f_5 * msd_315[k]
                   + f_3 * pc_x[k] * nsd_315[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, pa_x, pc_x, pc_z, msf0_526, msd_255, \
                         msd_316, msd_317, msf1_526, nsd_315, nsd_316, \
                         nsd_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_5 * msd_316[k]
                   + f_3 * pc_x[k] * nsd_316[k];

        t_525[k] = f_5 * msd_317[k]
                   + f_3 * pc_x[k] * nsd_317[k];

        t_526[k] = pa_x[k] * msf0_526[k]
                   - f_4 * pc_x[k] * msf1_526[k];

        t_527[k] = f_11 * msd_255[k]
                   + f_3 * pc_z[k] * nsd_315[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, pa_x, pa_y, pc_x, pc_y, msf0_440, \
                         msf0_529, msd_263, msd_264, msf1_440, msf1_529, nsd_317, \
                         nsd_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_10 * msd_263[k]
                   + f_3 * pc_y[k] * nsd_317[k];

        t_529[k] = pa_x[k] * msf0_529[k]
                   - f_4 * pc_x[k] * msf1_529[k];

        t_530[k] = pa_y[k] * msf0_440[k]
                   - f_4 * pc_y[k] * msf1_440[k];

        t_531[k] = f_5 * msd_264[k]
                   + f_3 * pc_y[k] * nsd_318[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pc_x, pc_z, msd_258, msd_321, msd_322, \
                         msd_323, nsd_318, nsd_321, nsd_322, nsd_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_9 * msd_258[k]
                   + f_3 * pc_z[k] * nsd_318[k];

        t_533[k] = f_5 * msd_321[k]
                   + f_3 * pc_x[k] * nsd_321[k];

        t_534[k] = f_5 * msd_322[k]
                   + f_3 * pc_x[k] * nsd_322[k];

        t_535[k] = f_5 * msd_323[k]
                   + f_3 * pc_x[k] * nsd_323[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pa_x, pc_x, pc_y, pc_z, msf0_536, \
                         msf0_539, msd_261, msd_269, msf1_536, msf1_539, nsd_321, \
                         nsd_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = pa_x[k] * msf0_536[k]
                   - f_4 * pc_x[k] * msf1_536[k];

        t_537[k] = f_9 * msd_261[k]
                   + f_3 * pc_z[k] * nsd_321[k];

        t_538[k] = f_5 * msd_269[k]
                   + f_3 * pc_y[k] * nsd_323[k];

        t_539[k] = pa_x[k] * msf0_539[k]
                   - f_4 * pc_x[k] * msf1_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pa_x, pc_x, pc_y, pc_z, msf0_540, \
                         msd_264, msd_324, msd_327, msf1_540, nsd_324, \
                         nsd_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = pa_x[k] * msf0_540[k]
                   + f_12 * msd_324[k]
                   - f_4 * pc_x[k] * msf1_540[k];

        t_541[k] = f_3 * pc_y[k] * nsd_324[k];

        t_542[k] = f_6 * msd_264[k]
                   + f_3 * pc_z[k] * nsd_324[k];

        t_543[k] = f_5 * msd_327[k]
                   + f_3 * pc_x[k] * nsd_327[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, t_548, pa_x, pc_x, pc_y, msf0_546, \
                         msf0_547, msd_329, msf1_546, msf1_547, nsd_326, \
                         nsd_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_3 * pc_y[k] * nsd_326[k];

        t_545[k] = f_5 * msd_329[k]
                   + f_3 * pc_x[k] * nsd_329[k];

        t_546[k] = pa_x[k] * msf0_546[k]
                   - f_4 * pc_x[k] * msf1_546[k];

        t_547[k] = pa_x[k] * msf0_547[k]
                   - f_4 * pc_x[k] * msf1_547[k];

        t_548[k] = f_3 * pc_y[k] * nsd_329[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, pa_x, pc_x, pc_z, msf0_549, msf1_549, \
                         nsp0_165, nsp0_166, nsp1_165, nsp1_166, nsd_330, \
                         nsd_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = pa_x[k] * msf0_549[k]
                   - f_4 * pc_x[k] * msf1_549[k];

        t_550[k] = f_1 * nsp0_165[k]
                   - f_2 * nsp1_165[k]
                   + f_3 * pc_x[k] * nsd_330[k];

        t_551[k] = f_7 * nsp0_166[k]
                   - f_8 * nsp1_166[k]
                   + f_3 * pc_x[k] * nsd_331[k];

        t_552[k] = f_3 * pc_z[k] * nsd_330[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, t_557, t_558, pc_x, pc_y, pc_z, msd_273, \
                         msd_275, nsp0_166, nsp1_166, nsd_333, nsd_334, \
                         nsd_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_3 * pc_x[k] * nsd_333[k];

        t_554[k] = f_3 * pc_x[k] * nsd_334[k];

        t_555[k] = f_3 * pc_x[k] * nsd_335[k];

        t_556[k] = f_0 * msd_273[k]
                   + f_1 * nsp0_166[k]
                   - f_2 * nsp1_166[k]
                   + f_3 * pc_y[k] * nsd_333[k];

        t_557[k] = f_3 * pc_z[k] * nsd_333[k];

        t_558[k] = f_0 * msd_275[k]
                   + f_3 * pc_y[k] * nsd_335[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, pa_z, pc_z, msf0_450, msf0_451, msf1_450, \
                         msf1_451, nsp0_167, nsp1_167, nsd_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_1 * nsp0_167[k]
                   - f_2 * nsp1_167[k]
                   + f_3 * pc_z[k] * nsd_335[k];

        t_560[k] = pa_z[k] * msf0_450[k]
                   - f_4 * pc_z[k] * msf1_450[k];

        t_561[k] = pa_z[k] * msf0_451[k]
                   - f_4 * pc_z[k] * msf1_451[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, pa_z, pc_x, pc_z, msf0_456, \
                         msf1_456, nsp0_170, nsp1_170, nsd_338, nsd_339, nsd_340, \
                         nsd_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_7 * nsp0_170[k]
                   - f_8 * nsp1_170[k]
                   + f_3 * pc_x[k] * nsd_338[k];

        t_563[k] = f_3 * pc_x[k] * nsd_339[k];

        t_564[k] = f_3 * pc_x[k] * nsd_340[k];

        t_565[k] = f_3 * pc_x[k] * nsd_341[k];

        t_566[k] = pa_z[k] * msf0_456[k]
                   - f_4 * pc_z[k] * msf1_456[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pc_y, pc_z, msd_273, msd_275, msd_281, nsp0_170, \
                         nsp1_170, nsd_339, nsd_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_5 * msd_273[k]
                   + f_3 * pc_z[k] * nsd_339[k];

        t_568[k] = f_6 * msd_281[k]
                   + f_3 * pc_y[k] * nsd_341[k];

        t_569[k] = f_5 * msd_275[k]
                   + f_1 * nsp0_170[k]
                   - f_2 * nsp1_170[k]
                   + f_3 * pc_z[k] * nsd_341[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, pc_x, nsp0_171, nsp0_172, nsp0_173, \
                         nsp1_171, nsp1_172, nsp1_173, nsd_342, nsd_343, nsd_344, \
                         nsd_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_1 * nsp0_171[k]
                   - f_2 * nsp1_171[k]
                   + f_3 * pc_x[k] * nsd_342[k];

        t_571[k] = f_7 * nsp0_172[k]
                   - f_8 * nsp1_172[k]
                   + f_3 * pc_x[k] * nsd_343[k];

        t_572[k] = f_7 * nsp0_173[k]
                   - f_8 * nsp1_173[k]
                   + f_3 * pc_x[k] * nsd_344[k];

        t_573[k] = f_3 * pc_x[k] * nsd_345[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, pc_x, pc_y, pc_z, msd_279, \
                         msd_285, msd_287, nsp0_172, nsp1_172, nsd_345, nsd_346, \
                         nsd_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_3 * pc_x[k] * nsd_346[k];

        t_575[k] = f_3 * pc_x[k] * nsd_347[k];

        t_576[k] = f_9 * msd_285[k]
                   + f_1 * nsp0_172[k]
                   - f_2 * nsp1_172[k]
                   + f_3 * pc_y[k] * nsd_345[k];

        t_577[k] = f_10 * msd_279[k]
                   + f_3 * pc_z[k] * nsd_345[k];

        t_578[k] = f_9 * msd_287[k]
                   + f_3 * pc_y[k] * nsd_347[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pc_x, pc_z, msd_281, nsp0_173, nsp0_174, \
                         nsp0_175, nsp1_173, nsp1_174, nsp1_175, nsd_347, nsd_348, \
                         nsd_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_10 * msd_281[k]
                   + f_1 * nsp0_173[k]
                   - f_2 * nsp1_173[k]
                   + f_3 * pc_z[k] * nsd_347[k];

        t_580[k] = f_1 * nsp0_174[k]
                   - f_2 * nsp1_174[k]
                   + f_3 * pc_x[k] * nsd_348[k];

        t_581[k] = f_7 * nsp0_175[k]
                   - f_8 * nsp1_175[k]
                   + f_3 * pc_x[k] * nsd_349[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, t_586, pc_x, pc_y, msd_291, nsp0_175, \
                         nsp0_176, nsp1_175, nsp1_176, nsd_350, nsd_351, nsd_352, \
                         nsd_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_7 * nsp0_176[k]
                   - f_8 * nsp1_176[k]
                   + f_3 * pc_x[k] * nsd_350[k];

        t_583[k] = f_3 * pc_x[k] * nsd_351[k];

        t_584[k] = f_3 * pc_x[k] * nsd_352[k];

        t_585[k] = f_3 * pc_x[k] * nsd_353[k];

        t_586[k] = f_11 * msd_291[k]
                   + f_1 * nsp0_175[k]
                   - f_2 * nsp1_175[k]
                   + f_3 * pc_y[k] * nsd_351[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pc_y, pc_z, msd_285, msd_287, msd_293, nsp0_176, \
                         nsp1_176, nsd_351, nsd_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_12 * msd_285[k]
                   + f_3 * pc_z[k] * nsd_351[k];

        t_588[k] = f_11 * msd_293[k]
                   + f_3 * pc_y[k] * nsd_353[k];

        t_589[k] = f_12 * msd_287[k]
                   + f_1 * nsp0_176[k]
                   - f_2 * nsp1_176[k]
                   + f_3 * pc_z[k] * nsd_353[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pc_x, nsp0_177, nsp0_178, nsp0_179, \
                         nsp1_177, nsp1_178, nsp1_179, nsd_354, nsd_355, nsd_356, \
                         nsd_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_1 * nsp0_177[k]
                   - f_2 * nsp1_177[k]
                   + f_3 * pc_x[k] * nsd_354[k];

        t_591[k] = f_7 * nsp0_178[k]
                   - f_8 * nsp1_178[k]
                   + f_3 * pc_x[k] * nsd_355[k];

        t_592[k] = f_7 * nsp0_179[k]
                   - f_8 * nsp1_179[k]
                   + f_3 * pc_x[k] * nsd_356[k];

        t_593[k] = f_3 * pc_x[k] * nsd_357[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, pc_x, pc_y, pc_z, msd_291, \
                         msd_297, msd_299, nsp0_178, nsp1_178, nsd_357, nsd_358, \
                         nsd_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_3 * pc_x[k] * nsd_358[k];

        t_595[k] = f_3 * pc_x[k] * nsd_359[k];

        t_596[k] = f_13 * msd_297[k]
                   + f_1 * nsp0_178[k]
                   - f_2 * nsp1_178[k]
                   + f_3 * pc_y[k] * nsd_357[k];

        t_597[k] = f_14 * msd_291[k]
                   + f_3 * pc_z[k] * nsd_357[k];

        t_598[k] = f_13 * msd_299[k]
                   + f_3 * pc_y[k] * nsd_359[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, pc_x, pc_z, msd_293, nsp0_179, nsp0_180, \
                         nsp0_181, nsp1_179, nsp1_180, nsp1_181, nsd_359, nsd_360, \
                         nsd_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = f_14 * msd_293[k]
                   + f_1 * nsp0_179[k]
                   - f_2 * nsp1_179[k]
                   + f_3 * pc_z[k] * nsd_359[k];

        t_600[k] = f_1 * nsp0_180[k]
                   - f_2 * nsp1_180[k]
                   + f_3 * pc_x[k] * nsd_360[k];

        t_601[k] = f_7 * nsp0_181[k]
                   - f_8 * nsp1_181[k]
                   + f_3 * pc_x[k] * nsd_361[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, pc_x, pc_y, msd_303, nsp0_181, \
                         nsp0_182, nsp1_181, nsp1_182, nsd_362, nsd_363, nsd_364, \
                         nsd_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_7 * nsp0_182[k]
                   - f_8 * nsp1_182[k]
                   + f_3 * pc_x[k] * nsd_362[k];

        t_603[k] = f_3 * pc_x[k] * nsd_363[k];

        t_604[k] = f_3 * pc_x[k] * nsd_364[k];

        t_605[k] = f_3 * pc_x[k] * nsd_365[k];

        t_606[k] = f_15 * msd_303[k]
                   + f_1 * nsp0_181[k]
                   - f_2 * nsp1_181[k]
                   + f_3 * pc_y[k] * nsd_363[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, pc_z, msd_297, msd_299, msd_305, nsp0_182, \
                         nsp1_182, nsd_363, nsd_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_15 * msd_297[k]
                   + f_3 * pc_z[k] * nsd_363[k];

        t_608[k] = f_15 * msd_305[k]
                   + f_3 * pc_y[k] * nsd_365[k];

        t_609[k] = f_15 * msd_299[k]
                   + f_1 * nsp0_182[k]
                   - f_2 * nsp1_182[k]
                   + f_3 * pc_z[k] * nsd_365[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pc_x, nsp0_183, nsp0_184, nsp0_185, \
                         nsp1_183, nsp1_184, nsp1_185, nsd_366, nsd_367, nsd_368, \
                         nsd_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_1 * nsp0_183[k]
                   - f_2 * nsp1_183[k]
                   + f_3 * pc_x[k] * nsd_366[k];

        t_611[k] = f_7 * nsp0_184[k]
                   - f_8 * nsp1_184[k]
                   + f_3 * pc_x[k] * nsd_367[k];

        t_612[k] = f_7 * nsp0_185[k]
                   - f_8 * nsp1_185[k]
                   + f_3 * pc_x[k] * nsd_368[k];

        t_613[k] = f_3 * pc_x[k] * nsd_369[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pc_x, pc_y, pc_z, msd_303, \
                         msd_309, msd_311, nsp0_184, nsp1_184, nsd_369, nsd_370, \
                         nsd_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_3 * pc_x[k] * nsd_370[k];

        t_615[k] = f_3 * pc_x[k] * nsd_371[k];

        t_616[k] = f_14 * msd_309[k]
                   + f_1 * nsp0_184[k]
                   - f_2 * nsp1_184[k]
                   + f_3 * pc_y[k] * nsd_369[k];

        t_617[k] = f_13 * msd_303[k]
                   + f_3 * pc_z[k] * nsd_369[k];

        t_618[k] = f_14 * msd_311[k]
                   + f_3 * pc_y[k] * nsd_371[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, pc_x, pc_z, msd_305, nsp0_185, nsp0_186, \
                         nsp0_187, nsp1_185, nsp1_186, nsp1_187, nsd_371, nsd_372, \
                         nsd_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_13 * msd_305[k]
                   + f_1 * nsp0_185[k]
                   - f_2 * nsp1_185[k]
                   + f_3 * pc_z[k] * nsd_371[k];

        t_620[k] = f_1 * nsp0_186[k]
                   - f_2 * nsp1_186[k]
                   + f_3 * pc_x[k] * nsd_372[k];

        t_621[k] = f_7 * nsp0_187[k]
                   - f_8 * nsp1_187[k]
                   + f_3 * pc_x[k] * nsd_373[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, pc_x, pc_y, msd_315, nsp0_187, \
                         nsp0_188, nsp1_187, nsp1_188, nsd_374, nsd_375, nsd_376, \
                         nsd_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_7 * nsp0_188[k]
                   - f_8 * nsp1_188[k]
                   + f_3 * pc_x[k] * nsd_374[k];

        t_623[k] = f_3 * pc_x[k] * nsd_375[k];

        t_624[k] = f_3 * pc_x[k] * nsd_376[k];

        t_625[k] = f_3 * pc_x[k] * nsd_377[k];

        t_626[k] = f_12 * msd_315[k]
                   + f_1 * nsp0_187[k]
                   - f_2 * nsp1_187[k]
                   + f_3 * pc_y[k] * nsd_375[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pc_y, pc_z, msd_309, msd_311, msd_317, nsp0_188, \
                         nsp1_188, nsd_375, nsd_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_11 * msd_309[k]
                   + f_3 * pc_z[k] * nsd_375[k];

        t_628[k] = f_12 * msd_317[k]
                   + f_3 * pc_y[k] * nsd_377[k];

        t_629[k] = f_11 * msd_311[k]
                   + f_1 * nsp0_188[k]
                   - f_2 * nsp1_188[k]
                   + f_3 * pc_z[k] * nsd_377[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, nsp0_189, nsp0_190, nsp0_191, \
                         nsp1_189, nsp1_190, nsp1_191, nsd_378, nsd_379, nsd_380, \
                         nsd_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_1 * nsp0_189[k]
                   - f_2 * nsp1_189[k]
                   + f_3 * pc_x[k] * nsd_378[k];

        t_631[k] = f_7 * nsp0_190[k]
                   - f_8 * nsp1_190[k]
                   + f_3 * pc_x[k] * nsd_379[k];

        t_632[k] = f_7 * nsp0_191[k]
                   - f_8 * nsp1_191[k]
                   + f_3 * pc_x[k] * nsd_380[k];

        t_633[k] = f_3 * pc_x[k] * nsd_381[k];
    }
}

static auto
compute_prim_nsf_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msf0,
                                                          const size_t msd, const size_t msf1,
                                                          const size_t nsp0, const size_t nsp1,
                                                          const size_t nsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msf0_540 = buffer.data(msf0 + 540);
    const auto *msf0_542 = buffer.data(msf0 + 542);
    const auto *msf0_546 = buffer.data(msf0 + 546);
    const auto *msf0_549 = buffer.data(msf0 + 549);

    const auto *msd_315 = buffer.data(msd + 315);
    const auto *msd_317 = buffer.data(msd + 317);
    const auto *msd_321 = buffer.data(msd + 321);
    const auto *msd_323 = buffer.data(msd + 323);
    const auto *msd_327 = buffer.data(msd + 327);
    const auto *msd_329 = buffer.data(msd + 329);

    const auto *msf1_540 = buffer.data(msf1 + 540);
    const auto *msf1_542 = buffer.data(msf1 + 542);
    const auto *msf1_546 = buffer.data(msf1 + 546);
    const auto *msf1_549 = buffer.data(msf1 + 549);

    const auto *nsp0_190 = buffer.data(nsp0 + 190);
    const auto *nsp0_191 = buffer.data(nsp0 + 191);
    const auto *nsp0_193 = buffer.data(nsp0 + 193);
    const auto *nsp0_195 = buffer.data(nsp0 + 195);
    const auto *nsp0_196 = buffer.data(nsp0 + 196);
    const auto *nsp0_197 = buffer.data(nsp0 + 197);

    const auto *nsp1_190 = buffer.data(nsp1 + 190);
    const auto *nsp1_191 = buffer.data(nsp1 + 191);
    const auto *nsp1_193 = buffer.data(nsp1 + 193);
    const auto *nsp1_195 = buffer.data(nsp1 + 195);
    const auto *nsp1_196 = buffer.data(nsp1 + 196);
    const auto *nsp1_197 = buffer.data(nsp1 + 197);

    const auto *nsd_381 = buffer.data(nsd + 381);
    const auto *nsd_382 = buffer.data(nsd + 382);
    const auto *nsd_383 = buffer.data(nsd + 383);
    const auto *nsd_385 = buffer.data(nsd + 385);
    const auto *nsd_387 = buffer.data(nsd + 387);
    const auto *nsd_388 = buffer.data(nsd + 388);
    const auto *nsd_389 = buffer.data(nsd + 389);
    const auto *nsd_390 = buffer.data(nsd + 390);
    const auto *nsd_392 = buffer.data(nsd + 392);
    const auto *nsd_393 = buffer.data(nsd + 393);
    const auto *nsd_394 = buffer.data(nsd + 394);
    const auto *nsd_395 = buffer.data(nsd + 395);

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, pc_x, pc_y, pc_z, msd_315, \
                         msd_321, msd_323, nsp0_190, nsp1_190, nsd_381, nsd_382, \
                         nsd_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_3 * pc_x[k] * nsd_382[k];

        t_635[k] = f_3 * pc_x[k] * nsd_383[k];

        t_636[k] = f_10 * msd_321[k]
                   + f_1 * nsp0_190[k]
                   - f_2 * nsp1_190[k]
                   + f_3 * pc_y[k] * nsd_381[k];

        t_637[k] = f_9 * msd_315[k]
                   + f_3 * pc_z[k] * nsd_381[k];

        t_638[k] = f_10 * msd_323[k]
                   + f_3 * pc_y[k] * nsd_383[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pa_y, pc_x, pc_y, pc_z, msf0_540, msd_317, \
                         msf1_540, nsp0_191, nsp0_193, nsp1_191, nsp1_193, nsd_383, \
                         nsd_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_9 * msd_317[k]
                   + f_1 * nsp0_191[k]
                   - f_2 * nsp1_191[k]
                   + f_3 * pc_z[k] * nsd_383[k];

        t_640[k] = pa_y[k] * msf0_540[k]
                   - f_4 * pc_y[k] * msf1_540[k];

        t_641[k] = f_7 * nsp0_193[k]
                   - f_8 * nsp1_193[k]
                   + f_3 * pc_x[k] * nsd_385[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, t_646, pa_y, pc_x, pc_y, msf0_542, \
                         msf0_546, msd_327, msf1_542, msf1_546, nsd_387, nsd_388, \
                         nsd_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = pa_y[k] * msf0_542[k]
                   - f_4 * pc_y[k] * msf1_542[k];

        t_643[k] = f_3 * pc_x[k] * nsd_387[k];

        t_644[k] = f_3 * pc_x[k] * nsd_388[k];

        t_645[k] = f_3 * pc_x[k] * nsd_389[k];

        t_646[k] = pa_y[k] * msf0_546[k]
                   + f_12 * msd_327[k]
                   - f_4 * pc_y[k] * msf1_546[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pa_y, pc_y, pc_z, msf0_549, msd_321, msd_329, \
                         msf1_549, nsd_387, nsd_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_6 * msd_321[k]
                   + f_3 * pc_z[k] * nsd_387[k];

        t_648[k] = f_5 * msd_329[k]
                   + f_3 * pc_y[k] * nsd_389[k];

        t_649[k] = pa_y[k] * msf0_549[k]
                   - f_4 * pc_y[k] * msf1_549[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, pc_x, pc_y, nsp0_195, nsp0_197, \
                         nsp1_195, nsp1_197, nsd_390, nsd_392, nsd_393, \
                         nsd_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_1 * nsp0_195[k]
                   - f_2 * nsp1_195[k]
                   + f_3 * pc_x[k] * nsd_390[k];

        t_651[k] = f_3 * pc_y[k] * nsd_390[k];

        t_652[k] = f_7 * nsp0_197[k]
                   - f_8 * nsp1_197[k]
                   + f_3 * pc_x[k] * nsd_392[k];

        t_653[k] = f_3 * pc_x[k] * nsd_393[k];

        t_654[k] = f_3 * pc_x[k] * nsd_394[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, t_658, t_659, pc_x, pc_y, pc_z, msd_329, \
                         nsp0_196, nsp0_197, nsp1_196, nsp1_197, nsd_393, nsd_394, \
                         nsd_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_3 * pc_x[k] * nsd_395[k];

        t_656[k] = f_1 * nsp0_196[k]
                   - f_2 * nsp1_196[k]
                   + f_3 * pc_y[k] * nsd_393[k];

        t_657[k] = f_7 * nsp0_197[k]
                   - f_8 * nsp1_197[k]
                   + f_3 * pc_y[k] * nsd_394[k];

        t_658[k] = f_3 * pc_y[k] * nsd_395[k];

        t_659[k] = f_0 * msd_329[k]
                   + f_1 * nsp0_197[k]
                   - f_2 * nsp1_197[k]
                   + f_3 * pc_z[k] * nsd_395[k];
    }
}

auto
compute_prim_nsf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t msf0, const size_t msd,
                                                   const size_t msf1, const size_t nsp0,
                                                   const size_t nsp1, const size_t nsd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_nsf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, msf0, msd,
                                                              msf1, nsp0, nsp1, nsd, ncols,
                                                              gamma, p, q);

    compute_prim_nsf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, msf0, msd,
                                                              msf1, nsp0, nsp1, nsd, ncols,
                                                              gamma, p, q);

    compute_prim_nsf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, msf0, msd,
                                                              msf1, nsp0, nsp1, nsd, ncols,
                                                              gamma, p, q);

    compute_prim_nsf_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, msf0, msd,
                                                              msf1, nsp0, nsp1, nsd, ncols,
                                                              gamma, p, q);

    compute_prim_nsf_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, msf0, msd,
                                                              msf1, nsp0, nsp1, nsd, ncols,
                                                              gamma, p, q);

    compute_prim_nsf_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, msf0, msd,
                                                              msf1, nsp0, nsp1, nsd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
