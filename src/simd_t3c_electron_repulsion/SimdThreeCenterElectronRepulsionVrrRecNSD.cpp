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


#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_nsd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msd0,
                                                          const size_t msp, const size_t msd1,
                                                          const size_t nss0, const size_t nss1,
                                                          const size_t nsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 4.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msd0_0 = buffer.data(msd0 + 0);
    const auto *msd0_3 = buffer.data(msd0 + 3);
    const auto *msd0_5 = buffer.data(msd0 + 5);
    const auto *msd0_9 = buffer.data(msd0 + 9);
    const auto *msd0_12 = buffer.data(msd0 + 12);
    const auto *msd0_17 = buffer.data(msd0 + 17);
    const auto *msd0_18 = buffer.data(msd0 + 18);
    const auto *msd0_21 = buffer.data(msd0 + 21);
    const auto *msd0_30 = buffer.data(msd0 + 30);
    const auto *msd0_35 = buffer.data(msd0 + 35);
    const auto *msd0_36 = buffer.data(msd0 + 36);
    const auto *msd0_39 = buffer.data(msd0 + 39);
    const auto *msd0_54 = buffer.data(msd0 + 54);
    const auto *msd0_59 = buffer.data(msd0 + 59);
    const auto *msd0_60 = buffer.data(msd0 + 60);
    const auto *msd0_63 = buffer.data(msd0 + 63);
    const auto *msd0_84 = buffer.data(msd0 + 84);
    const auto *msd0_89 = buffer.data(msd0 + 89);

    const auto *msp_0 = buffer.data(msp + 0);
    const auto *msp_1 = buffer.data(msp + 1);
    const auto *msp_2 = buffer.data(msp + 2);
    const auto *msp_4 = buffer.data(msp + 4);
    const auto *msp_8 = buffer.data(msp + 8);
    const auto *msp_9 = buffer.data(msp + 9);
    const auto *msp_10 = buffer.data(msp + 10);
    const auto *msp_11 = buffer.data(msp + 11);
    const auto *msp_13 = buffer.data(msp + 13);
    const auto *msp_14 = buffer.data(msp + 14);
    const auto *msp_15 = buffer.data(msp + 15);
    const auto *msp_16 = buffer.data(msp + 16);
    const auto *msp_17 = buffer.data(msp + 17);
    const auto *msp_18 = buffer.data(msp + 18);
    const auto *msp_19 = buffer.data(msp + 19);
    const auto *msp_20 = buffer.data(msp + 20);
    const auto *msp_22 = buffer.data(msp + 22);
    const auto *msp_23 = buffer.data(msp + 23);
    const auto *msp_25 = buffer.data(msp + 25);
    const auto *msp_26 = buffer.data(msp + 26);
    const auto *msp_27 = buffer.data(msp + 27);
    const auto *msp_28 = buffer.data(msp + 28);
    const auto *msp_29 = buffer.data(msp + 29);
    const auto *msp_30 = buffer.data(msp + 30);
    const auto *msp_31 = buffer.data(msp + 31);
    const auto *msp_32 = buffer.data(msp + 32);
    const auto *msp_34 = buffer.data(msp + 34);
    const auto *msp_35 = buffer.data(msp + 35);
    const auto *msp_36 = buffer.data(msp + 36);
    const auto *msp_37 = buffer.data(msp + 37);
    const auto *msp_38 = buffer.data(msp + 38);
    const auto *msp_40 = buffer.data(msp + 40);
    const auto *msp_41 = buffer.data(msp + 41);
    const auto *msp_42 = buffer.data(msp + 42);
    const auto *msp_43 = buffer.data(msp + 43);
    const auto *msp_44 = buffer.data(msp + 44);
    const auto *msp_45 = buffer.data(msp + 45);
    const auto *msp_46 = buffer.data(msp + 46);
    const auto *msp_49 = buffer.data(msp + 49);
    const auto *msp_50 = buffer.data(msp + 50);
    const auto *msp_51 = buffer.data(msp + 51);
    const auto *msp_52 = buffer.data(msp + 52);
    const auto *msp_53 = buffer.data(msp + 53);
    const auto *msp_54 = buffer.data(msp + 54);
    const auto *msp_55 = buffer.data(msp + 55);
    const auto *msp_56 = buffer.data(msp + 56);
    const auto *msp_58 = buffer.data(msp + 58);
    const auto *msp_59 = buffer.data(msp + 59);
    const auto *msp_60 = buffer.data(msp + 60);
    const auto *msp_62 = buffer.data(msp + 62);

    const auto *msd1_0 = buffer.data(msd1 + 0);
    const auto *msd1_3 = buffer.data(msd1 + 3);
    const auto *msd1_5 = buffer.data(msd1 + 5);
    const auto *msd1_9 = buffer.data(msd1 + 9);
    const auto *msd1_12 = buffer.data(msd1 + 12);
    const auto *msd1_17 = buffer.data(msd1 + 17);
    const auto *msd1_18 = buffer.data(msd1 + 18);
    const auto *msd1_21 = buffer.data(msd1 + 21);
    const auto *msd1_30 = buffer.data(msd1 + 30);
    const auto *msd1_35 = buffer.data(msd1 + 35);
    const auto *msd1_36 = buffer.data(msd1 + 36);
    const auto *msd1_39 = buffer.data(msd1 + 39);
    const auto *msd1_54 = buffer.data(msd1 + 54);
    const auto *msd1_59 = buffer.data(msd1 + 59);
    const auto *msd1_60 = buffer.data(msd1 + 60);
    const auto *msd1_63 = buffer.data(msd1 + 63);
    const auto *msd1_84 = buffer.data(msd1 + 84);
    const auto *msd1_89 = buffer.data(msd1 + 89);

    const auto *nss0_0 = buffer.data(nss0 + 0);
    const auto *nss0_1 = buffer.data(nss0 + 1);
    const auto *nss0_2 = buffer.data(nss0 + 2);
    const auto *nss0_3 = buffer.data(nss0 + 3);
    const auto *nss0_5 = buffer.data(nss0 + 5);
    const auto *nss0_6 = buffer.data(nss0 + 6);
    const auto *nss0_7 = buffer.data(nss0 + 7);
    const auto *nss0_8 = buffer.data(nss0 + 8);
    const auto *nss0_9 = buffer.data(nss0 + 9);
    const auto *nss0_10 = buffer.data(nss0 + 10);
    const auto *nss0_11 = buffer.data(nss0 + 11);
    const auto *nss0_12 = buffer.data(nss0 + 12);
    const auto *nss0_13 = buffer.data(nss0 + 13);
    const auto *nss0_14 = buffer.data(nss0 + 14);
    const auto *nss0_15 = buffer.data(nss0 + 15);
    const auto *nss0_16 = buffer.data(nss0 + 16);
    const auto *nss0_17 = buffer.data(nss0 + 17);
    const auto *nss0_18 = buffer.data(nss0 + 18);
    const auto *nss0_19 = buffer.data(nss0 + 19);
    const auto *nss0_20 = buffer.data(nss0 + 20);

    const auto *nss1_0 = buffer.data(nss1 + 0);
    const auto *nss1_1 = buffer.data(nss1 + 1);
    const auto *nss1_2 = buffer.data(nss1 + 2);
    const auto *nss1_3 = buffer.data(nss1 + 3);
    const auto *nss1_5 = buffer.data(nss1 + 5);
    const auto *nss1_6 = buffer.data(nss1 + 6);
    const auto *nss1_7 = buffer.data(nss1 + 7);
    const auto *nss1_8 = buffer.data(nss1 + 8);
    const auto *nss1_9 = buffer.data(nss1 + 9);
    const auto *nss1_10 = buffer.data(nss1 + 10);
    const auto *nss1_11 = buffer.data(nss1 + 11);
    const auto *nss1_12 = buffer.data(nss1 + 12);
    const auto *nss1_13 = buffer.data(nss1 + 13);
    const auto *nss1_14 = buffer.data(nss1 + 14);
    const auto *nss1_15 = buffer.data(nss1 + 15);
    const auto *nss1_16 = buffer.data(nss1 + 16);
    const auto *nss1_17 = buffer.data(nss1 + 17);
    const auto *nss1_18 = buffer.data(nss1 + 18);
    const auto *nss1_19 = buffer.data(nss1 + 19);
    const auto *nss1_20 = buffer.data(nss1 + 20);

    const auto *nsp_0 = buffer.data(nsp + 0);
    const auto *nsp_1 = buffer.data(nsp + 1);
    const auto *nsp_2 = buffer.data(nsp + 2);
    const auto *nsp_3 = buffer.data(nsp + 3);
    const auto *nsp_4 = buffer.data(nsp + 4);
    const auto *nsp_6 = buffer.data(nsp + 6);
    const auto *nsp_8 = buffer.data(nsp + 8);
    const auto *nsp_9 = buffer.data(nsp + 9);
    const auto *nsp_10 = buffer.data(nsp + 10);
    const auto *nsp_11 = buffer.data(nsp + 11);
    const auto *nsp_13 = buffer.data(nsp + 13);
    const auto *nsp_14 = buffer.data(nsp + 14);
    const auto *nsp_15 = buffer.data(nsp + 15);
    const auto *nsp_16 = buffer.data(nsp + 16);
    const auto *nsp_17 = buffer.data(nsp + 17);
    const auto *nsp_18 = buffer.data(nsp + 18);
    const auto *nsp_19 = buffer.data(nsp + 19);
    const auto *nsp_20 = buffer.data(nsp + 20);
    const auto *nsp_22 = buffer.data(nsp + 22);
    const auto *nsp_23 = buffer.data(nsp + 23);
    const auto *nsp_25 = buffer.data(nsp + 25);
    const auto *nsp_26 = buffer.data(nsp + 26);
    const auto *nsp_27 = buffer.data(nsp + 27);
    const auto *nsp_28 = buffer.data(nsp + 28);
    const auto *nsp_29 = buffer.data(nsp + 29);
    const auto *nsp_30 = buffer.data(nsp + 30);
    const auto *nsp_31 = buffer.data(nsp + 31);
    const auto *nsp_32 = buffer.data(nsp + 32);
    const auto *nsp_34 = buffer.data(nsp + 34);
    const auto *nsp_35 = buffer.data(nsp + 35);
    const auto *nsp_36 = buffer.data(nsp + 36);
    const auto *nsp_37 = buffer.data(nsp + 37);
    const auto *nsp_38 = buffer.data(nsp + 38);
    const auto *nsp_40 = buffer.data(nsp + 40);
    const auto *nsp_41 = buffer.data(nsp + 41);
    const auto *nsp_42 = buffer.data(nsp + 42);
    const auto *nsp_43 = buffer.data(nsp + 43);
    const auto *nsp_44 = buffer.data(nsp + 44);
    const auto *nsp_45 = buffer.data(nsp + 45);
    const auto *nsp_46 = buffer.data(nsp + 46);
    const auto *nsp_47 = buffer.data(nsp + 47);
    const auto *nsp_49 = buffer.data(nsp + 49);
    const auto *nsp_50 = buffer.data(nsp + 50);
    const auto *nsp_51 = buffer.data(nsp + 51);
    const auto *nsp_52 = buffer.data(nsp + 52);
    const auto *nsp_53 = buffer.data(nsp + 53);
    const auto *nsp_54 = buffer.data(nsp + 54);
    const auto *nsp_55 = buffer.data(nsp + 55);
    const auto *nsp_56 = buffer.data(nsp + 56);
    const auto *nsp_58 = buffer.data(nsp + 58);
    const auto *nsp_59 = buffer.data(nsp + 59);
    const auto *nsp_60 = buffer.data(nsp + 60);
    const auto *nsp_61 = buffer.data(nsp + 61);
    const auto *nsp_62 = buffer.data(nsp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, msp_0, nss0_0, \
                         nss1_0, nsp_0, nsp_1, nsp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * msp_0[k]
                 + f_1 * nss0_0[k]
                 - f_2 * nss1_0[k]
                 + f_3 * pc_x[k] * nsp_0[k];

        t_1[k] = f_3 * pc_y[k] * nsp_0[k];

        t_2[k] = f_3 * pc_z[k] * nsp_0[k];

        t_3[k] = f_1 * nss0_0[k]
                 - f_2 * nss1_0[k]
                 + f_3 * pc_y[k] * nsp_1[k];

        t_4[k] = f_3 * pc_y[k] * nsp_2[k];

        t_5[k] = f_1 * nss0_0[k]
                 - f_2 * nss1_0[k]
                 + f_3 * pc_z[k] * nsp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, msd0_0, msp_1, msp_4, \
                         msd1_0, nss0_1, nss1_1, nsp_3, nsp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * msd0_0[k]
                 - f_4 * pc_y[k] * msd1_0[k];

        t_7[k] = f_5 * msp_4[k]
                 + f_3 * pc_x[k] * nsp_4[k];

        t_8[k] = f_3 * pc_z[k] * nsp_3[k];

        t_9[k] = f_6 * msp_1[k]
                 + f_1 * nss0_1[k]
                 - f_2 * nss1_1[k]
                 + f_3 * pc_y[k] * nsp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, msd0_0, msd0_5, \
                         msd1_0, msd1_5, nsp_4, nsp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * nsp_4[k];

        t_11[k] = pa_y[k] * msd0_5[k]
                  - f_4 * pc_y[k] * msd1_5[k];

        t_12[k] = pa_z[k] * msd0_0[k]
                  - f_4 * pc_z[k] * msd1_0[k];

        t_13[k] = f_3 * pc_y[k] * nsp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, msd0_3, msp_2, msp_8, \
                         msd1_3, nss0_2, nss1_2, nsp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * msp_8[k]
                  + f_3 * pc_x[k] * nsp_8[k];

        t_15[k] = pa_z[k] * msd0_3[k]
                  - f_4 * pc_z[k] * msd1_3[k];

        t_16[k] = f_3 * pc_y[k] * nsp_8[k];

        t_17[k] = f_6 * msp_2[k]
                  + f_1 * nss0_2[k]
                  - f_2 * nss1_2[k]
                  + f_3 * pc_z[k] * nsp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, msp_4, msp_9, msp_10, \
                         nss0_3, nss1_3, nsp_9, nsp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * msp_9[k]
                  + f_1 * nss0_3[k]
                  - f_2 * nss1_3[k]
                  + f_3 * pc_x[k] * nsp_9[k];

        t_19[k] = f_7 * msp_10[k]
                  + f_3 * pc_x[k] * nsp_10[k];

        t_20[k] = f_3 * pc_z[k] * nsp_9[k];

        t_21[k] = f_8 * msp_4[k]
                  + f_1 * nss0_3[k]
                  - f_2 * nss1_3[k]
                  + f_3 * pc_y[k] * nsp_10[k];

        t_22[k] = f_3 * pc_z[k] * nsp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, msd0_12, msp_13, msd1_12, \
                         nss0_3, nss1_3, nsp_11, nsp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * nss0_3[k]
                  - f_2 * nss1_3[k]
                  + f_3 * pc_z[k] * nsp_11[k];

        t_24[k] = pa_y[k] * msd0_12[k]
                  - f_4 * pc_y[k] * msd1_12[k];

        t_25[k] = f_7 * msp_13[k]
                  + f_3 * pc_x[k] * nsp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pc_x, pc_y, pc_z, msd0_9, \
                         msd0_17, msp_8, msp_14, msd1_9, msd1_17, \
                         nsp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * msp_14[k]
                  + f_3 * pc_x[k] * nsp_14[k];

        t_27[k] = pa_z[k] * msd0_9[k]
                  - f_4 * pc_z[k] * msd1_9[k];

        t_28[k] = f_6 * msp_8[k]
                  + f_3 * pc_y[k] * nsp_14[k];

        t_29[k] = pa_y[k] * msd0_17[k]
                  - f_4 * pc_y[k] * msd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, msp_15, msp_17, nss0_5, \
                         nss1_5, nsp_15, nsp_16, nsp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * msp_15[k]
                  + f_1 * nss0_5[k]
                  - f_2 * nss1_5[k]
                  + f_3 * pc_x[k] * nsp_15[k];

        t_31[k] = f_3 * pc_y[k] * nsp_15[k];

        t_32[k] = f_7 * msp_17[k]
                  + f_3 * pc_x[k] * nsp_17[k];

        t_33[k] = f_1 * nss0_5[k]
                  - f_2 * nss1_5[k]
                  + f_3 * pc_y[k] * nsp_16[k];

        t_34[k] = f_3 * pc_y[k] * nsp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_z, msp_8, msp_18, msp_19, nss0_5, \
                         nss0_6, nss1_5, nss1_6, nsp_17, nsp_18, \
                         nsp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * msp_8[k]
                  + f_1 * nss0_5[k]
                  - f_2 * nss1_5[k]
                  + f_3 * pc_z[k] * nsp_17[k];

        t_36[k] = f_9 * msp_18[k]
                  + f_1 * nss0_6[k]
                  - f_2 * nss1_6[k]
                  + f_3 * pc_x[k] * nsp_18[k];

        t_37[k] = f_9 * msp_19[k]
                  + f_3 * pc_x[k] * nsp_19[k];

        t_38[k] = f_3 * pc_z[k] * nsp_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pc_y, pc_z, msd0_18, msp_10, msd1_18, \
                         nss0_6, nss1_6, nsp_19, nsp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * msp_10[k]
                  + f_1 * nss0_6[k]
                  - f_2 * nss1_6[k]
                  + f_3 * pc_y[k] * nsp_19[k];

        t_40[k] = f_3 * pc_z[k] * nsp_19[k];

        t_41[k] = f_1 * nss0_6[k]
                  - f_2 * nss1_6[k]
                  + f_3 * pc_z[k] * nsp_20[k];

        t_42[k] = pa_z[k] * msd0_18[k]
                  - f_4 * pc_z[k] * msd1_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pc_x, pc_y, pc_z, msd0_21, msp_14, \
                         msp_22, msp_23, msd1_21, nsp_22, nsp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_9 * msp_22[k]
                  + f_3 * pc_x[k] * nsp_22[k];

        t_44[k] = f_9 * msp_23[k]
                  + f_3 * pc_x[k] * nsp_23[k];

        t_45[k] = pa_z[k] * msd0_21[k]
                  - f_4 * pc_z[k] * msd1_21[k];

        t_46[k] = f_8 * msp_14[k]
                  + f_3 * pc_y[k] * nsp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pc_x, pc_y, pc_z, msd0_30, msp_11, msp_25, \
                         msd1_30, nss0_7, nss1_7, nsp_23, nsp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * msp_11[k]
                  + f_1 * nss0_7[k]
                  - f_2 * nss1_7[k]
                  + f_3 * pc_z[k] * nsp_23[k];

        t_48[k] = pa_y[k] * msd0_30[k]
                  - f_4 * pc_y[k] * msd1_30[k];

        t_49[k] = f_9 * msp_25[k]
                  + f_3 * pc_x[k] * nsp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_x, pc_y, msd0_35, msp_16, msp_17, \
                         msp_26, msd1_35, nss0_8, nss1_8, nsp_25, \
                         nsp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * msp_26[k]
                  + f_3 * pc_x[k] * nsp_26[k];

        t_51[k] = f_6 * msp_16[k]
                  + f_1 * nss0_8[k]
                  - f_2 * nss1_8[k]
                  + f_3 * pc_y[k] * nsp_25[k];

        t_52[k] = f_6 * msp_17[k]
                  + f_3 * pc_y[k] * nsp_26[k];

        t_53[k] = pa_y[k] * msd0_35[k]
                  - f_4 * pc_y[k] * msd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, msp_27, msp_29, nss0_9, \
                         nss1_9, nsp_27, nsp_28, nsp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * msp_27[k]
                  + f_1 * nss0_9[k]
                  - f_2 * nss1_9[k]
                  + f_3 * pc_x[k] * nsp_27[k];

        t_55[k] = f_3 * pc_y[k] * nsp_27[k];

        t_56[k] = f_9 * msp_29[k]
                  + f_3 * pc_x[k] * nsp_29[k];

        t_57[k] = f_1 * nss0_9[k]
                  - f_2 * nss1_9[k]
                  + f_3 * pc_y[k] * nsp_28[k];

        t_58[k] = f_3 * pc_y[k] * nsp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_z, msp_17, msp_30, msp_31, nss0_9, \
                         nss0_10, nss1_9, nss1_10, nsp_29, nsp_30, \
                         nsp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_10 * msp_17[k]
                  + f_1 * nss0_9[k]
                  - f_2 * nss1_9[k]
                  + f_3 * pc_z[k] * nsp_29[k];

        t_60[k] = f_11 * msp_30[k]
                  + f_1 * nss0_10[k]
                  - f_2 * nss1_10[k]
                  + f_3 * pc_x[k] * nsp_30[k];

        t_61[k] = f_11 * msp_31[k]
                  + f_3 * pc_x[k] * nsp_31[k];

        t_62[k] = f_3 * pc_z[k] * nsp_30[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pc_y, pc_z, msd0_36, msp_19, msd1_36, \
                         nss0_10, nss1_10, nsp_31, nsp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_12 * msp_19[k]
                  + f_1 * nss0_10[k]
                  - f_2 * nss1_10[k]
                  + f_3 * pc_y[k] * nsp_31[k];

        t_64[k] = f_3 * pc_z[k] * nsp_31[k];

        t_65[k] = f_1 * nss0_10[k]
                  - f_2 * nss1_10[k]
                  + f_3 * pc_z[k] * nsp_32[k];

        t_66[k] = pa_z[k] * msd0_36[k]
                  - f_4 * pc_z[k] * msd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pc_x, pc_y, pc_z, msd0_39, msp_23, \
                         msp_34, msp_35, msd1_39, nsp_34, nsp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_11 * msp_34[k]
                  + f_3 * pc_x[k] * nsp_34[k];

        t_68[k] = f_11 * msp_35[k]
                  + f_3 * pc_x[k] * nsp_35[k];

        t_69[k] = pa_z[k] * msd0_39[k]
                  - f_4 * pc_z[k] * msd1_39[k];

        t_70[k] = f_10 * msp_23[k]
                  + f_3 * pc_y[k] * nsp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pc_x, pc_z, msp_20, msp_36, msp_37, nss0_11, \
                         nss0_12, nss1_11, nss1_12, nsp_35, nsp_36, \
                         nsp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_6 * msp_20[k]
                  + f_1 * nss0_11[k]
                  - f_2 * nss1_11[k]
                  + f_3 * pc_z[k] * nsp_35[k];

        t_72[k] = f_11 * msp_36[k]
                  + f_1 * nss0_12[k]
                  - f_2 * nss1_12[k]
                  + f_3 * pc_x[k] * nsp_36[k];

        t_73[k] = f_11 * msp_37[k]
                  + f_3 * pc_x[k] * nsp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, msp_23, msp_25, msp_26, \
                         msp_38, nss0_12, nss1_12, nsp_37, nsp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_11 * msp_38[k]
                  + f_3 * pc_x[k] * nsp_38[k];

        t_75[k] = f_8 * msp_25[k]
                  + f_1 * nss0_12[k]
                  - f_2 * nss1_12[k]
                  + f_3 * pc_y[k] * nsp_37[k];

        t_76[k] = f_8 * msp_26[k]
                  + f_3 * pc_y[k] * nsp_38[k];

        t_77[k] = f_8 * msp_23[k]
                  + f_1 * nss0_12[k]
                  - f_2 * nss1_12[k]
                  + f_3 * pc_z[k] * nsp_38[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_y, pc_x, pc_y, msd0_54, msp_28, msp_40, \
                         msp_41, msd1_54, nss0_13, nss1_13, nsp_40, \
                         nsp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * msd0_54[k]
                  - f_4 * pc_y[k] * msd1_54[k];

        t_79[k] = f_11 * msp_40[k]
                  + f_3 * pc_x[k] * nsp_40[k];

        t_80[k] = f_11 * msp_41[k]
                  + f_3 * pc_x[k] * nsp_41[k];

        t_81[k] = f_6 * msp_28[k]
                  + f_1 * nss0_13[k]
                  - f_2 * nss1_13[k]
                  + f_3 * pc_y[k] * nsp_40[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_y, pc_x, pc_y, msd0_59, msp_29, msp_42, \
                         msd1_59, nss0_14, nss1_14, nsp_41, nsp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * msp_29[k]
                  + f_3 * pc_y[k] * nsp_41[k];

        t_83[k] = pa_y[k] * msd0_59[k]
                  - f_4 * pc_y[k] * msd1_59[k];

        t_84[k] = f_11 * msp_42[k]
                  + f_1 * nss0_14[k]
                  - f_2 * nss1_14[k]
                  + f_3 * pc_x[k] * nsp_42[k];

        t_85[k] = f_3 * pc_y[k] * nsp_42[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, msp_29, msp_44, nss0_14, \
                         nss1_14, nsp_43, nsp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * msp_44[k]
                  + f_3 * pc_x[k] * nsp_44[k];

        t_87[k] = f_1 * nss0_14[k]
                  - f_2 * nss1_14[k]
                  + f_3 * pc_y[k] * nsp_43[k];

        t_88[k] = f_3 * pc_y[k] * nsp_44[k];

        t_89[k] = f_12 * msp_29[k]
                  + f_1 * nss0_14[k]
                  - f_2 * nss1_14[k]
                  + f_3 * pc_z[k] * nsp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, msp_31, msp_45, \
                         msp_46, nss0_15, nss1_15, nsp_45, nsp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_13 * msp_45[k]
                  + f_1 * nss0_15[k]
                  - f_2 * nss1_15[k]
                  + f_3 * pc_x[k] * nsp_45[k];

        t_91[k] = f_13 * msp_46[k]
                  + f_3 * pc_x[k] * nsp_46[k];

        t_92[k] = f_3 * pc_z[k] * nsp_45[k];

        t_93[k] = f_13 * msp_31[k]
                  + f_1 * nss0_15[k]
                  - f_2 * nss1_15[k]
                  + f_3 * pc_y[k] * nsp_46[k];

        t_94[k] = f_3 * pc_z[k] * nsp_46[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_z, pc_x, pc_z, msd0_60, msp_49, msp_50, \
                         msd1_60, nss0_15, nss1_15, nsp_47, nsp_49, \
                         nsp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_1 * nss0_15[k]
                  - f_2 * nss1_15[k]
                  + f_3 * pc_z[k] * nsp_47[k];

        t_96[k] = pa_z[k] * msd0_60[k]
                  - f_4 * pc_z[k] * msd1_60[k];

        t_97[k] = f_13 * msp_49[k]
                  + f_3 * pc_x[k] * nsp_49[k];

        t_98[k] = f_13 * msp_50[k]
                  + f_3 * pc_x[k] * nsp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, msd0_63, msp_32, msp_35, \
                         msd1_63, nss0_16, nss1_16, nsp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * msd0_63[k]
                  - f_4 * pc_z[k] * msd1_63[k];

        t_100[k] = f_12 * msp_35[k]
                   + f_3 * pc_y[k] * nsp_50[k];

        t_101[k] = f_6 * msp_32[k]
                   + f_1 * nss0_16[k]
                   - f_2 * nss1_16[k]
                   + f_3 * pc_z[k] * nsp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, msp_37, msp_51, msp_52, \
                         msp_53, nss0_17, nss1_17, nsp_51, nsp_52, \
                         nsp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_13 * msp_51[k]
                   + f_1 * nss0_17[k]
                   - f_2 * nss1_17[k]
                   + f_3 * pc_x[k] * nsp_51[k];

        t_103[k] = f_13 * msp_52[k]
                   + f_3 * pc_x[k] * nsp_52[k];

        t_104[k] = f_13 * msp_53[k]
                   + f_3 * pc_x[k] * nsp_53[k];

        t_105[k] = f_10 * msp_37[k]
                   + f_1 * nss0_17[k]
                   - f_2 * nss1_17[k]
                   + f_3 * pc_y[k] * nsp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, msp_35, msp_38, msp_54, \
                         nss0_17, nss0_18, nss1_17, nss1_18, nsp_53, \
                         nsp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * msp_38[k]
                   + f_3 * pc_y[k] * nsp_53[k];

        t_107[k] = f_8 * msp_35[k]
                   + f_1 * nss0_17[k]
                   - f_2 * nss1_17[k]
                   + f_3 * pc_z[k] * nsp_53[k];

        t_108[k] = f_13 * msp_54[k]
                   + f_1 * nss0_18[k]
                   - f_2 * nss1_18[k]
                   + f_3 * pc_x[k] * nsp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, msp_40, msp_41, msp_55, \
                         msp_56, nss0_18, nss1_18, nsp_55, nsp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_13 * msp_55[k]
                   + f_3 * pc_x[k] * nsp_55[k];

        t_110[k] = f_13 * msp_56[k]
                   + f_3 * pc_x[k] * nsp_56[k];

        t_111[k] = f_8 * msp_40[k]
                   + f_1 * nss0_18[k]
                   - f_2 * nss1_18[k]
                   + f_3 * pc_y[k] * nsp_55[k];

        t_112[k] = f_8 * msp_41[k]
                   + f_3 * pc_y[k] * nsp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_y, pc_x, pc_y, pc_z, msd0_84, msp_38, msp_58, \
                         msd1_84, nss0_18, nss1_18, nsp_56, nsp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * msp_38[k]
                   + f_1 * nss0_18[k]
                   - f_2 * nss1_18[k]
                   + f_3 * pc_z[k] * nsp_56[k];

        t_114[k] = pa_y[k] * msd0_84[k]
                   - f_4 * pc_y[k] * msd1_84[k];

        t_115[k] = f_13 * msp_58[k]
                   + f_3 * pc_x[k] * nsp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_y, pc_x, pc_y, msd0_89, msp_43, \
                         msp_44, msp_59, msd1_89, nss0_19, nss1_19, nsp_58, \
                         nsp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_13 * msp_59[k]
                   + f_3 * pc_x[k] * nsp_59[k];

        t_117[k] = f_6 * msp_43[k]
                   + f_1 * nss0_19[k]
                   - f_2 * nss1_19[k]
                   + f_3 * pc_y[k] * nsp_58[k];

        t_118[k] = f_6 * msp_44[k]
                   + f_3 * pc_y[k] * nsp_59[k];

        t_119[k] = pa_y[k] * msd0_89[k]
                   - f_4 * pc_y[k] * msd1_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, msp_60, msp_62, \
                         nss0_20, nss1_20, nsp_60, nsp_61, nsp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_13 * msp_60[k]
                   + f_1 * nss0_20[k]
                   - f_2 * nss1_20[k]
                   + f_3 * pc_x[k] * nsp_60[k];

        t_121[k] = f_3 * pc_y[k] * nsp_60[k];

        t_122[k] = f_13 * msp_62[k]
                   + f_3 * pc_x[k] * nsp_62[k];

        t_123[k] = f_1 * nss0_20[k]
                   - f_2 * nss1_20[k]
                   + f_3 * pc_y[k] * nsp_61[k];

        t_124[k] = f_3 * pc_y[k] * nsp_62[k];
    }
}

static auto
compute_prim_nsd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msd0,
                                                          const size_t msp, const size_t msd1,
                                                          const size_t nss0, const size_t nss1,
                                                          const size_t nsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msd0_90 = buffer.data(msd0 + 90);
    const auto *msd0_93 = buffer.data(msd0 + 93);
    const auto *msd0_120 = buffer.data(msd0 + 120);
    const auto *msd0_125 = buffer.data(msd0 + 125);
    const auto *msd0_126 = buffer.data(msd0 + 126);
    const auto *msd0_129 = buffer.data(msd0 + 129);
    const auto *msd0_162 = buffer.data(msd0 + 162);
    const auto *msd0_167 = buffer.data(msd0 + 167);
    const auto *msd0_168 = buffer.data(msd0 + 168);
    const auto *msd0_171 = buffer.data(msd0 + 171);

    const auto *msp_44 = buffer.data(msp + 44);
    const auto *msp_46 = buffer.data(msp + 46);
    const auto *msp_47 = buffer.data(msp + 47);
    const auto *msp_50 = buffer.data(msp + 50);
    const auto *msp_52 = buffer.data(msp + 52);
    const auto *msp_53 = buffer.data(msp + 53);
    const auto *msp_55 = buffer.data(msp + 55);
    const auto *msp_56 = buffer.data(msp + 56);
    const auto *msp_58 = buffer.data(msp + 58);
    const auto *msp_59 = buffer.data(msp + 59);
    const auto *msp_61 = buffer.data(msp + 61);
    const auto *msp_62 = buffer.data(msp + 62);
    const auto *msp_63 = buffer.data(msp + 63);
    const auto *msp_64 = buffer.data(msp + 64);
    const auto *msp_65 = buffer.data(msp + 65);
    const auto *msp_67 = buffer.data(msp + 67);
    const auto *msp_68 = buffer.data(msp + 68);
    const auto *msp_69 = buffer.data(msp + 69);
    const auto *msp_70 = buffer.data(msp + 70);
    const auto *msp_71 = buffer.data(msp + 71);
    const auto *msp_72 = buffer.data(msp + 72);
    const auto *msp_73 = buffer.data(msp + 73);
    const auto *msp_74 = buffer.data(msp + 74);
    const auto *msp_75 = buffer.data(msp + 75);
    const auto *msp_76 = buffer.data(msp + 76);
    const auto *msp_77 = buffer.data(msp + 77);
    const auto *msp_79 = buffer.data(msp + 79);
    const auto *msp_80 = buffer.data(msp + 80);
    const auto *msp_81 = buffer.data(msp + 81);
    const auto *msp_82 = buffer.data(msp + 82);
    const auto *msp_83 = buffer.data(msp + 83);
    const auto *msp_84 = buffer.data(msp + 84);
    const auto *msp_85 = buffer.data(msp + 85);
    const auto *msp_86 = buffer.data(msp + 86);
    const auto *msp_88 = buffer.data(msp + 88);
    const auto *msp_89 = buffer.data(msp + 89);
    const auto *msp_90 = buffer.data(msp + 90);
    const auto *msp_91 = buffer.data(msp + 91);
    const auto *msp_92 = buffer.data(msp + 92);
    const auto *msp_93 = buffer.data(msp + 93);
    const auto *msp_94 = buffer.data(msp + 94);
    const auto *msp_95 = buffer.data(msp + 95);
    const auto *msp_96 = buffer.data(msp + 96);
    const auto *msp_97 = buffer.data(msp + 97);
    const auto *msp_98 = buffer.data(msp + 98);
    const auto *msp_99 = buffer.data(msp + 99);
    const auto *msp_100 = buffer.data(msp + 100);
    const auto *msp_101 = buffer.data(msp + 101);
    const auto *msp_103 = buffer.data(msp + 103);
    const auto *msp_104 = buffer.data(msp + 104);
    const auto *msp_105 = buffer.data(msp + 105);
    const auto *msp_107 = buffer.data(msp + 107);
    const auto *msp_108 = buffer.data(msp + 108);
    const auto *msp_109 = buffer.data(msp + 109);
    const auto *msp_112 = buffer.data(msp + 112);
    const auto *msp_113 = buffer.data(msp + 113);
    const auto *msp_114 = buffer.data(msp + 114);
    const auto *msp_115 = buffer.data(msp + 115);
    const auto *msp_116 = buffer.data(msp + 116);
    const auto *msp_117 = buffer.data(msp + 117);
    const auto *msp_118 = buffer.data(msp + 118);
    const auto *msp_119 = buffer.data(msp + 119);
    const auto *msp_120 = buffer.data(msp + 120);
    const auto *msp_121 = buffer.data(msp + 121);

    const auto *msd1_90 = buffer.data(msd1 + 90);
    const auto *msd1_93 = buffer.data(msd1 + 93);
    const auto *msd1_120 = buffer.data(msd1 + 120);
    const auto *msd1_125 = buffer.data(msd1 + 125);
    const auto *msd1_126 = buffer.data(msd1 + 126);
    const auto *msd1_129 = buffer.data(msd1 + 129);
    const auto *msd1_162 = buffer.data(msd1 + 162);
    const auto *msd1_167 = buffer.data(msd1 + 167);
    const auto *msd1_168 = buffer.data(msd1 + 168);
    const auto *msd1_171 = buffer.data(msd1 + 171);

    const auto *nss0_20 = buffer.data(nss0 + 20);
    const auto *nss0_21 = buffer.data(nss0 + 21);
    const auto *nss0_22 = buffer.data(nss0 + 22);
    const auto *nss0_23 = buffer.data(nss0 + 23);
    const auto *nss0_24 = buffer.data(nss0 + 24);
    const auto *nss0_25 = buffer.data(nss0 + 25);
    const auto *nss0_26 = buffer.data(nss0 + 26);
    const auto *nss0_27 = buffer.data(nss0 + 27);
    const auto *nss0_28 = buffer.data(nss0 + 28);
    const auto *nss0_29 = buffer.data(nss0 + 29);
    const auto *nss0_30 = buffer.data(nss0 + 30);
    const auto *nss0_31 = buffer.data(nss0 + 31);
    const auto *nss0_32 = buffer.data(nss0 + 32);
    const auto *nss0_33 = buffer.data(nss0 + 33);
    const auto *nss0_34 = buffer.data(nss0 + 34);
    const auto *nss0_35 = buffer.data(nss0 + 35);
    const auto *nss0_36 = buffer.data(nss0 + 36);
    const auto *nss0_37 = buffer.data(nss0 + 37);
    const auto *nss0_38 = buffer.data(nss0 + 38);
    const auto *nss0_39 = buffer.data(nss0 + 39);
    const auto *nss0_40 = buffer.data(nss0 + 40);

    const auto *nss1_20 = buffer.data(nss1 + 20);
    const auto *nss1_21 = buffer.data(nss1 + 21);
    const auto *nss1_22 = buffer.data(nss1 + 22);
    const auto *nss1_23 = buffer.data(nss1 + 23);
    const auto *nss1_24 = buffer.data(nss1 + 24);
    const auto *nss1_25 = buffer.data(nss1 + 25);
    const auto *nss1_26 = buffer.data(nss1 + 26);
    const auto *nss1_27 = buffer.data(nss1 + 27);
    const auto *nss1_28 = buffer.data(nss1 + 28);
    const auto *nss1_29 = buffer.data(nss1 + 29);
    const auto *nss1_30 = buffer.data(nss1 + 30);
    const auto *nss1_31 = buffer.data(nss1 + 31);
    const auto *nss1_32 = buffer.data(nss1 + 32);
    const auto *nss1_33 = buffer.data(nss1 + 33);
    const auto *nss1_34 = buffer.data(nss1 + 34);
    const auto *nss1_35 = buffer.data(nss1 + 35);
    const auto *nss1_36 = buffer.data(nss1 + 36);
    const auto *nss1_37 = buffer.data(nss1 + 37);
    const auto *nss1_38 = buffer.data(nss1 + 38);
    const auto *nss1_39 = buffer.data(nss1 + 39);
    const auto *nss1_40 = buffer.data(nss1 + 40);

    const auto *nsp_62 = buffer.data(nsp + 62);
    const auto *nsp_63 = buffer.data(nsp + 63);
    const auto *nsp_64 = buffer.data(nsp + 64);
    const auto *nsp_65 = buffer.data(nsp + 65);
    const auto *nsp_67 = buffer.data(nsp + 67);
    const auto *nsp_68 = buffer.data(nsp + 68);
    const auto *nsp_69 = buffer.data(nsp + 69);
    const auto *nsp_70 = buffer.data(nsp + 70);
    const auto *nsp_71 = buffer.data(nsp + 71);
    const auto *nsp_72 = buffer.data(nsp + 72);
    const auto *nsp_73 = buffer.data(nsp + 73);
    const auto *nsp_74 = buffer.data(nsp + 74);
    const auto *nsp_75 = buffer.data(nsp + 75);
    const auto *nsp_76 = buffer.data(nsp + 76);
    const auto *nsp_77 = buffer.data(nsp + 77);
    const auto *nsp_79 = buffer.data(nsp + 79);
    const auto *nsp_80 = buffer.data(nsp + 80);
    const auto *nsp_81 = buffer.data(nsp + 81);
    const auto *nsp_82 = buffer.data(nsp + 82);
    const auto *nsp_83 = buffer.data(nsp + 83);
    const auto *nsp_84 = buffer.data(nsp + 84);
    const auto *nsp_85 = buffer.data(nsp + 85);
    const auto *nsp_86 = buffer.data(nsp + 86);
    const auto *nsp_88 = buffer.data(nsp + 88);
    const auto *nsp_89 = buffer.data(nsp + 89);
    const auto *nsp_90 = buffer.data(nsp + 90);
    const auto *nsp_91 = buffer.data(nsp + 91);
    const auto *nsp_92 = buffer.data(nsp + 92);
    const auto *nsp_93 = buffer.data(nsp + 93);
    const auto *nsp_94 = buffer.data(nsp + 94);
    const auto *nsp_95 = buffer.data(nsp + 95);
    const auto *nsp_96 = buffer.data(nsp + 96);
    const auto *nsp_97 = buffer.data(nsp + 97);
    const auto *nsp_98 = buffer.data(nsp + 98);
    const auto *nsp_99 = buffer.data(nsp + 99);
    const auto *nsp_100 = buffer.data(nsp + 100);
    const auto *nsp_101 = buffer.data(nsp + 101);
    const auto *nsp_103 = buffer.data(nsp + 103);
    const auto *nsp_104 = buffer.data(nsp + 104);
    const auto *nsp_105 = buffer.data(nsp + 105);
    const auto *nsp_106 = buffer.data(nsp + 106);
    const auto *nsp_107 = buffer.data(nsp + 107);
    const auto *nsp_108 = buffer.data(nsp + 108);
    const auto *nsp_109 = buffer.data(nsp + 109);
    const auto *nsp_110 = buffer.data(nsp + 110);
    const auto *nsp_112 = buffer.data(nsp + 112);
    const auto *nsp_113 = buffer.data(nsp + 113);
    const auto *nsp_114 = buffer.data(nsp + 114);
    const auto *nsp_115 = buffer.data(nsp + 115);
    const auto *nsp_116 = buffer.data(nsp + 116);
    const auto *nsp_117 = buffer.data(nsp + 117);
    const auto *nsp_118 = buffer.data(nsp + 118);
    const auto *nsp_119 = buffer.data(nsp + 119);
    const auto *nsp_120 = buffer.data(nsp + 120);
    const auto *nsp_121 = buffer.data(nsp + 121);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_z, msp_44, msp_63, msp_64, \
                         nss0_20, nss0_21, nss1_20, nss1_21, nsp_62, nsp_63, \
                         nsp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_13 * msp_44[k]
                   + f_1 * nss0_20[k]
                   - f_2 * nss1_20[k]
                   + f_3 * pc_z[k] * nsp_62[k];

        t_126[k] = f_12 * msp_63[k]
                   + f_1 * nss0_21[k]
                   - f_2 * nss1_21[k]
                   + f_3 * pc_x[k] * nsp_63[k];

        t_127[k] = f_12 * msp_64[k]
                   + f_3 * pc_x[k] * nsp_64[k];

        t_128[k] = f_3 * pc_z[k] * nsp_63[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_z, pc_y, pc_z, msd0_90, msp_46, \
                         msd1_90, nss0_21, nss1_21, nsp_64, nsp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * msp_46[k]
                   + f_1 * nss0_21[k]
                   - f_2 * nss1_21[k]
                   + f_3 * pc_y[k] * nsp_64[k];

        t_130[k] = f_3 * pc_z[k] * nsp_64[k];

        t_131[k] = f_1 * nss0_21[k]
                   - f_2 * nss1_21[k]
                   + f_3 * pc_z[k] * nsp_65[k];

        t_132[k] = pa_z[k] * msd0_90[k]
                   - f_4 * pc_z[k] * msd1_90[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_z, pc_x, pc_y, pc_z, msd0_93, msp_50, \
                         msp_67, msp_68, msd1_93, nsp_67, nsp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_12 * msp_67[k]
                   + f_3 * pc_x[k] * nsp_67[k];

        t_134[k] = f_12 * msp_68[k]
                   + f_3 * pc_x[k] * nsp_68[k];

        t_135[k] = pa_z[k] * msd0_93[k]
                   - f_4 * pc_z[k] * msd1_93[k];

        t_136[k] = f_13 * msp_50[k]
                   + f_3 * pc_y[k] * nsp_68[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_x, pc_z, msp_47, msp_69, msp_70, nss0_22, \
                         nss0_23, nss1_22, nss1_23, nsp_68, nsp_69, \
                         nsp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_6 * msp_47[k]
                   + f_1 * nss0_22[k]
                   - f_2 * nss1_22[k]
                   + f_3 * pc_z[k] * nsp_68[k];

        t_138[k] = f_12 * msp_69[k]
                   + f_1 * nss0_23[k]
                   - f_2 * nss1_23[k]
                   + f_3 * pc_x[k] * nsp_69[k];

        t_139[k] = f_12 * msp_70[k]
                   + f_3 * pc_x[k] * nsp_70[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, pc_z, msp_50, msp_52, msp_53, \
                         msp_71, nss0_23, nss1_23, nsp_70, nsp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_12 * msp_71[k]
                   + f_3 * pc_x[k] * nsp_71[k];

        t_141[k] = f_12 * msp_52[k]
                   + f_1 * nss0_23[k]
                   - f_2 * nss1_23[k]
                   + f_3 * pc_y[k] * nsp_70[k];

        t_142[k] = f_12 * msp_53[k]
                   + f_3 * pc_y[k] * nsp_71[k];

        t_143[k] = f_8 * msp_50[k]
                   + f_1 * nss0_23[k]
                   - f_2 * nss1_23[k]
                   + f_3 * pc_z[k] * nsp_71[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pc_x, pc_y, msp_55, msp_72, msp_73, \
                         msp_74, nss0_24, nss1_24, nsp_72, nsp_73, \
                         nsp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_12 * msp_72[k]
                   + f_1 * nss0_24[k]
                   - f_2 * nss1_24[k]
                   + f_3 * pc_x[k] * nsp_72[k];

        t_145[k] = f_12 * msp_73[k]
                   + f_3 * pc_x[k] * nsp_73[k];

        t_146[k] = f_12 * msp_74[k]
                   + f_3 * pc_x[k] * nsp_74[k];

        t_147[k] = f_10 * msp_55[k]
                   + f_1 * nss0_24[k]
                   - f_2 * nss1_24[k]
                   + f_3 * pc_y[k] * nsp_73[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pc_x, pc_y, pc_z, msp_53, msp_56, msp_75, \
                         nss0_24, nss0_25, nss1_24, nss1_25, nsp_74, \
                         nsp_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * msp_56[k]
                   + f_3 * pc_y[k] * nsp_74[k];

        t_149[k] = f_10 * msp_53[k]
                   + f_1 * nss0_24[k]
                   - f_2 * nss1_24[k]
                   + f_3 * pc_z[k] * nsp_74[k];

        t_150[k] = f_12 * msp_75[k]
                   + f_1 * nss0_25[k]
                   - f_2 * nss1_25[k]
                   + f_3 * pc_x[k] * nsp_75[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, msp_58, msp_59, msp_76, \
                         msp_77, nss0_25, nss1_25, nsp_76, nsp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_12 * msp_76[k]
                   + f_3 * pc_x[k] * nsp_76[k];

        t_152[k] = f_12 * msp_77[k]
                   + f_3 * pc_x[k] * nsp_77[k];

        t_153[k] = f_8 * msp_58[k]
                   + f_1 * nss0_25[k]
                   - f_2 * nss1_25[k]
                   + f_3 * pc_y[k] * nsp_76[k];

        t_154[k] = f_8 * msp_59[k]
                   + f_3 * pc_y[k] * nsp_77[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pa_y, pc_x, pc_y, pc_z, msd0_120, msp_56, \
                         msp_79, msd1_120, nss0_25, nss1_25, nsp_77, \
                         nsp_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_12 * msp_56[k]
                   + f_1 * nss0_25[k]
                   - f_2 * nss1_25[k]
                   + f_3 * pc_z[k] * nsp_77[k];

        t_156[k] = pa_y[k] * msd0_120[k]
                   - f_4 * pc_y[k] * msd1_120[k];

        t_157[k] = f_12 * msp_79[k]
                   + f_3 * pc_x[k] * nsp_79[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_y, pc_x, pc_y, msd0_125, msp_61, \
                         msp_62, msp_80, msd1_125, nss0_26, nss1_26, nsp_79, \
                         nsp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_12 * msp_80[k]
                   + f_3 * pc_x[k] * nsp_80[k];

        t_159[k] = f_6 * msp_61[k]
                   + f_1 * nss0_26[k]
                   - f_2 * nss1_26[k]
                   + f_3 * pc_y[k] * nsp_79[k];

        t_160[k] = f_6 * msp_62[k]
                   + f_3 * pc_y[k] * nsp_80[k];

        t_161[k] = pa_y[k] * msd0_125[k]
                   - f_4 * pc_y[k] * msd1_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pc_x, pc_y, msp_81, msp_83, \
                         nss0_27, nss1_27, nsp_81, nsp_82, nsp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_12 * msp_81[k]
                   + f_1 * nss0_27[k]
                   - f_2 * nss1_27[k]
                   + f_3 * pc_x[k] * nsp_81[k];

        t_163[k] = f_3 * pc_y[k] * nsp_81[k];

        t_164[k] = f_12 * msp_83[k]
                   + f_3 * pc_x[k] * nsp_83[k];

        t_165[k] = f_1 * nss0_27[k]
                   - f_2 * nss1_27[k]
                   + f_3 * pc_y[k] * nsp_82[k];

        t_166[k] = f_3 * pc_y[k] * nsp_83[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pc_x, pc_z, msp_62, msp_84, msp_85, \
                         nss0_27, nss0_28, nss1_27, nss1_28, nsp_83, nsp_84, \
                         nsp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_11 * msp_62[k]
                   + f_1 * nss0_27[k]
                   - f_2 * nss1_27[k]
                   + f_3 * pc_z[k] * nsp_83[k];

        t_168[k] = f_10 * msp_84[k]
                   + f_1 * nss0_28[k]
                   - f_2 * nss1_28[k]
                   + f_3 * pc_x[k] * nsp_84[k];

        t_169[k] = f_10 * msp_85[k]
                   + f_3 * pc_x[k] * nsp_85[k];

        t_170[k] = f_3 * pc_z[k] * nsp_84[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_z, pc_y, pc_z, msd0_126, msp_64, \
                         msd1_126, nss0_28, nss1_28, nsp_85, nsp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_9 * msp_64[k]
                   + f_1 * nss0_28[k]
                   - f_2 * nss1_28[k]
                   + f_3 * pc_y[k] * nsp_85[k];

        t_172[k] = f_3 * pc_z[k] * nsp_85[k];

        t_173[k] = f_1 * nss0_28[k]
                   - f_2 * nss1_28[k]
                   + f_3 * pc_z[k] * nsp_86[k];

        t_174[k] = pa_z[k] * msd0_126[k]
                   - f_4 * pc_z[k] * msd1_126[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_z, pc_x, pc_y, pc_z, msd0_129, msp_68, \
                         msp_88, msp_89, msd1_129, nsp_88, nsp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * msp_88[k]
                   + f_3 * pc_x[k] * nsp_88[k];

        t_176[k] = f_10 * msp_89[k]
                   + f_3 * pc_x[k] * nsp_89[k];

        t_177[k] = pa_z[k] * msd0_129[k]
                   - f_4 * pc_z[k] * msd1_129[k];

        t_178[k] = f_11 * msp_68[k]
                   + f_3 * pc_y[k] * nsp_89[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pc_x, pc_z, msp_65, msp_90, msp_91, nss0_29, \
                         nss0_30, nss1_29, nss1_30, nsp_89, nsp_90, \
                         nsp_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_6 * msp_65[k]
                   + f_1 * nss0_29[k]
                   - f_2 * nss1_29[k]
                   + f_3 * pc_z[k] * nsp_89[k];

        t_180[k] = f_10 * msp_90[k]
                   + f_1 * nss0_30[k]
                   - f_2 * nss1_30[k]
                   + f_3 * pc_x[k] * nsp_90[k];

        t_181[k] = f_10 * msp_91[k]
                   + f_3 * pc_x[k] * nsp_91[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pc_x, pc_y, pc_z, msp_68, msp_70, msp_71, \
                         msp_92, nss0_30, nss1_30, nsp_91, nsp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_10 * msp_92[k]
                   + f_3 * pc_x[k] * nsp_92[k];

        t_183[k] = f_13 * msp_70[k]
                   + f_1 * nss0_30[k]
                   - f_2 * nss1_30[k]
                   + f_3 * pc_y[k] * nsp_91[k];

        t_184[k] = f_13 * msp_71[k]
                   + f_3 * pc_y[k] * nsp_92[k];

        t_185[k] = f_8 * msp_68[k]
                   + f_1 * nss0_30[k]
                   - f_2 * nss1_30[k]
                   + f_3 * pc_z[k] * nsp_92[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, pc_y, msp_73, msp_93, msp_94, \
                         msp_95, nss0_31, nss1_31, nsp_93, nsp_94, \
                         nsp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_10 * msp_93[k]
                   + f_1 * nss0_31[k]
                   - f_2 * nss1_31[k]
                   + f_3 * pc_x[k] * nsp_93[k];

        t_187[k] = f_10 * msp_94[k]
                   + f_3 * pc_x[k] * nsp_94[k];

        t_188[k] = f_10 * msp_95[k]
                   + f_3 * pc_x[k] * nsp_95[k];

        t_189[k] = f_12 * msp_73[k]
                   + f_1 * nss0_31[k]
                   - f_2 * nss1_31[k]
                   + f_3 * pc_y[k] * nsp_94[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_y, pc_z, msp_71, msp_74, msp_96, \
                         nss0_31, nss0_32, nss1_31, nss1_32, nsp_95, \
                         nsp_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_12 * msp_74[k]
                   + f_3 * pc_y[k] * nsp_95[k];

        t_191[k] = f_10 * msp_71[k]
                   + f_1 * nss0_31[k]
                   - f_2 * nss1_31[k]
                   + f_3 * pc_z[k] * nsp_95[k];

        t_192[k] = f_10 * msp_96[k]
                   + f_1 * nss0_32[k]
                   - f_2 * nss1_32[k]
                   + f_3 * pc_x[k] * nsp_96[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, msp_76, msp_77, msp_97, \
                         msp_98, nss0_32, nss1_32, nsp_97, nsp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_10 * msp_97[k]
                   + f_3 * pc_x[k] * nsp_97[k];

        t_194[k] = f_10 * msp_98[k]
                   + f_3 * pc_x[k] * nsp_98[k];

        t_195[k] = f_10 * msp_76[k]
                   + f_1 * nss0_32[k]
                   - f_2 * nss1_32[k]
                   + f_3 * pc_y[k] * nsp_97[k];

        t_196[k] = f_10 * msp_77[k]
                   + f_3 * pc_y[k] * nsp_98[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pc_x, pc_z, msp_74, msp_99, msp_100, nss0_32, \
                         nss0_33, nss1_32, nss1_33, nsp_98, nsp_99, \
                         nsp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * msp_74[k]
                   + f_1 * nss0_32[k]
                   - f_2 * nss1_32[k]
                   + f_3 * pc_z[k] * nsp_98[k];

        t_198[k] = f_10 * msp_99[k]
                   + f_1 * nss0_33[k]
                   - f_2 * nss1_33[k]
                   + f_3 * pc_x[k] * nsp_99[k];

        t_199[k] = f_10 * msp_100[k]
                   + f_3 * pc_x[k] * nsp_100[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, msp_77, msp_79, msp_80, \
                         msp_101, nss0_33, nss1_33, nsp_100, nsp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_10 * msp_101[k]
                   + f_3 * pc_x[k] * nsp_101[k];

        t_201[k] = f_8 * msp_79[k]
                   + f_1 * nss0_33[k]
                   - f_2 * nss1_33[k]
                   + f_3 * pc_y[k] * nsp_100[k];

        t_202[k] = f_8 * msp_80[k]
                   + f_3 * pc_y[k] * nsp_101[k];

        t_203[k] = f_13 * msp_77[k]
                   + f_1 * nss0_33[k]
                   - f_2 * nss1_33[k]
                   + f_3 * pc_z[k] * nsp_101[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_x, pc_y, msd0_162, msp_82, \
                         msp_103, msp_104, msd1_162, nss0_34, nss1_34, nsp_103, \
                         nsp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * msd0_162[k]
                   - f_4 * pc_y[k] * msd1_162[k];

        t_205[k] = f_10 * msp_103[k]
                   + f_3 * pc_x[k] * nsp_103[k];

        t_206[k] = f_10 * msp_104[k]
                   + f_3 * pc_x[k] * nsp_104[k];

        t_207[k] = f_6 * msp_82[k]
                   + f_1 * nss0_34[k]
                   - f_2 * nss1_34[k]
                   + f_3 * pc_y[k] * nsp_103[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_y, pc_x, pc_y, msd0_167, msp_83, \
                         msp_105, msd1_167, nss0_35, nss1_35, nsp_104, \
                         nsp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_6 * msp_83[k]
                   + f_3 * pc_y[k] * nsp_104[k];

        t_209[k] = pa_y[k] * msd0_167[k]
                   - f_4 * pc_y[k] * msd1_167[k];

        t_210[k] = f_10 * msp_105[k]
                   + f_1 * nss0_35[k]
                   - f_2 * nss1_35[k]
                   + f_3 * pc_x[k] * nsp_105[k];

        t_211[k] = f_3 * pc_y[k] * nsp_105[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, msp_83, msp_107, \
                         nss0_35, nss1_35, nsp_106, nsp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_10 * msp_107[k]
                   + f_3 * pc_x[k] * nsp_107[k];

        t_213[k] = f_1 * nss0_35[k]
                   - f_2 * nss1_35[k]
                   + f_3 * pc_y[k] * nsp_106[k];

        t_214[k] = f_3 * pc_y[k] * nsp_107[k];

        t_215[k] = f_9 * msp_83[k]
                   + f_1 * nss0_35[k]
                   - f_2 * nss1_35[k]
                   + f_3 * pc_z[k] * nsp_107[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, pc_x, pc_y, pc_z, msp_85, msp_108, \
                         msp_109, nss0_36, nss1_36, nsp_108, nsp_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_8 * msp_108[k]
                   + f_1 * nss0_36[k]
                   - f_2 * nss1_36[k]
                   + f_3 * pc_x[k] * nsp_108[k];

        t_217[k] = f_8 * msp_109[k]
                   + f_3 * pc_x[k] * nsp_109[k];

        t_218[k] = f_3 * pc_z[k] * nsp_108[k];

        t_219[k] = f_7 * msp_85[k]
                   + f_1 * nss0_36[k]
                   - f_2 * nss1_36[k]
                   + f_3 * pc_y[k] * nsp_109[k];

        t_220[k] = f_3 * pc_z[k] * nsp_109[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_z, pc_x, pc_z, msd0_168, msp_112, \
                         msp_113, msd1_168, nss0_36, nss1_36, nsp_110, nsp_112, \
                         nsp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_1 * nss0_36[k]
                   - f_2 * nss1_36[k]
                   + f_3 * pc_z[k] * nsp_110[k];

        t_222[k] = pa_z[k] * msd0_168[k]
                   - f_4 * pc_z[k] * msd1_168[k];

        t_223[k] = f_8 * msp_112[k]
                   + f_3 * pc_x[k] * nsp_112[k];

        t_224[k] = f_8 * msp_113[k]
                   + f_3 * pc_x[k] * nsp_113[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_z, pc_y, pc_z, msd0_171, msp_86, msp_89, \
                         msd1_171, nss0_37, nss1_37, nsp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = pa_z[k] * msd0_171[k]
                   - f_4 * pc_z[k] * msd1_171[k];

        t_226[k] = f_9 * msp_89[k]
                   + f_3 * pc_y[k] * nsp_113[k];

        t_227[k] = f_6 * msp_86[k]
                   + f_1 * nss0_37[k]
                   - f_2 * nss1_37[k]
                   + f_3 * pc_z[k] * nsp_113[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pc_x, pc_y, msp_91, msp_114, msp_115, \
                         msp_116, nss0_38, nss1_38, nsp_114, nsp_115, \
                         nsp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_8 * msp_114[k]
                   + f_1 * nss0_38[k]
                   - f_2 * nss1_38[k]
                   + f_3 * pc_x[k] * nsp_114[k];

        t_229[k] = f_8 * msp_115[k]
                   + f_3 * pc_x[k] * nsp_115[k];

        t_230[k] = f_8 * msp_116[k]
                   + f_3 * pc_x[k] * nsp_116[k];

        t_231[k] = f_11 * msp_91[k]
                   + f_1 * nss0_38[k]
                   - f_2 * nss1_38[k]
                   + f_3 * pc_y[k] * nsp_115[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_y, pc_z, msp_89, msp_92, msp_117, \
                         nss0_38, nss0_39, nss1_38, nss1_39, nsp_116, \
                         nsp_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_11 * msp_92[k]
                   + f_3 * pc_y[k] * nsp_116[k];

        t_233[k] = f_8 * msp_89[k]
                   + f_1 * nss0_38[k]
                   - f_2 * nss1_38[k]
                   + f_3 * pc_z[k] * nsp_116[k];

        t_234[k] = f_8 * msp_117[k]
                   + f_1 * nss0_39[k]
                   - f_2 * nss1_39[k]
                   + f_3 * pc_x[k] * nsp_117[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, msp_94, msp_95, msp_118, \
                         msp_119, nss0_39, nss1_39, nsp_118, nsp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_8 * msp_118[k]
                   + f_3 * pc_x[k] * nsp_118[k];

        t_236[k] = f_8 * msp_119[k]
                   + f_3 * pc_x[k] * nsp_119[k];

        t_237[k] = f_13 * msp_94[k]
                   + f_1 * nss0_39[k]
                   - f_2 * nss1_39[k]
                   + f_3 * pc_y[k] * nsp_118[k];

        t_238[k] = f_13 * msp_95[k]
                   + f_3 * pc_y[k] * nsp_119[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, msp_92, msp_120, msp_121, nss0_39, \
                         nss0_40, nss1_39, nss1_40, nsp_119, nsp_120, \
                         nsp_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_10 * msp_92[k]
                   + f_1 * nss0_39[k]
                   - f_2 * nss1_39[k]
                   + f_3 * pc_z[k] * nsp_119[k];

        t_240[k] = f_8 * msp_120[k]
                   + f_1 * nss0_40[k]
                   - f_2 * nss1_40[k]
                   + f_3 * pc_x[k] * nsp_120[k];

        t_241[k] = f_8 * msp_121[k]
                   + f_3 * pc_x[k] * nsp_121[k];
    }
}

static auto
compute_prim_nsd_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msd0,
                                                          const size_t msp, const size_t msd1,
                                                          const size_t nss0, const size_t nss1,
                                                          const size_t nsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 4.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msd0_210 = buffer.data(msd0 + 210);
    const auto *msd0_215 = buffer.data(msd0 + 215);
    const auto *msd0_216 = buffer.data(msd0 + 216);
    const auto *msd0_264 = buffer.data(msd0 + 264);
    const auto *msd0_270 = buffer.data(msd0 + 270);
    const auto *msd0_273 = buffer.data(msd0 + 273);
    const auto *msd0_275 = buffer.data(msd0 + 275);
    const auto *msd0_279 = buffer.data(msd0 + 279);
    const auto *msd0_281 = buffer.data(msd0 + 281);
    const auto *msd0_282 = buffer.data(msd0 + 282);
    const auto *msd0_285 = buffer.data(msd0 + 285);
    const auto *msd0_287 = buffer.data(msd0 + 287);
    const auto *msd0_288 = buffer.data(msd0 + 288);
    const auto *msd0_291 = buffer.data(msd0 + 291);
    const auto *msd0_293 = buffer.data(msd0 + 293);
    const auto *msd0_294 = buffer.data(msd0 + 294);
    const auto *msd0_297 = buffer.data(msd0 + 297);
    const auto *msd0_299 = buffer.data(msd0 + 299);
    const auto *msd0_300 = buffer.data(msd0 + 300);
    const auto *msd0_303 = buffer.data(msd0 + 303);
    const auto *msd0_305 = buffer.data(msd0 + 305);
    const auto *msd0_306 = buffer.data(msd0 + 306);
    const auto *msd0_309 = buffer.data(msd0 + 309);
    const auto *msd0_311 = buffer.data(msd0 + 311);
    const auto *msd0_312 = buffer.data(msd0 + 312);
    const auto *msd0_315 = buffer.data(msd0 + 315);
    const auto *msd0_317 = buffer.data(msd0 + 317);
    const auto *msd0_321 = buffer.data(msd0 + 321);
    const auto *msd0_323 = buffer.data(msd0 + 323);
    const auto *msd0_324 = buffer.data(msd0 + 324);
    const auto *msd0_327 = buffer.data(msd0 + 327);
    const auto *msd0_329 = buffer.data(msd0 + 329);

    const auto *msp_95 = buffer.data(msp + 95);
    const auto *msp_97 = buffer.data(msp + 97);
    const auto *msp_98 = buffer.data(msp + 98);
    const auto *msp_100 = buffer.data(msp + 100);
    const auto *msp_101 = buffer.data(msp + 101);
    const auto *msp_103 = buffer.data(msp + 103);
    const auto *msp_104 = buffer.data(msp + 104);
    const auto *msp_106 = buffer.data(msp + 106);
    const auto *msp_107 = buffer.data(msp + 107);
    const auto *msp_113 = buffer.data(msp + 113);
    const auto *msp_116 = buffer.data(msp + 116);
    const auto *msp_119 = buffer.data(msp + 119);
    const auto *msp_122 = buffer.data(msp + 122);
    const auto *msp_123 = buffer.data(msp + 123);
    const auto *msp_124 = buffer.data(msp + 124);
    const auto *msp_125 = buffer.data(msp + 125);
    const auto *msp_126 = buffer.data(msp + 126);
    const auto *msp_127 = buffer.data(msp + 127);
    const auto *msp_128 = buffer.data(msp + 128);
    const auto *msp_130 = buffer.data(msp + 130);
    const auto *msp_131 = buffer.data(msp + 131);
    const auto *msp_132 = buffer.data(msp + 132);
    const auto *msp_134 = buffer.data(msp + 134);
    const auto *msp_135 = buffer.data(msp + 135);
    const auto *msp_136 = buffer.data(msp + 136);
    const auto *msp_137 = buffer.data(msp + 137);
    const auto *msp_139 = buffer.data(msp + 139);
    const auto *msp_140 = buffer.data(msp + 140);
    const auto *msp_141 = buffer.data(msp + 141);
    const auto *msp_142 = buffer.data(msp + 142);
    const auto *msp_143 = buffer.data(msp + 143);
    const auto *msp_144 = buffer.data(msp + 144);
    const auto *msp_145 = buffer.data(msp + 145);
    const auto *msp_146 = buffer.data(msp + 146);
    const auto *msp_147 = buffer.data(msp + 147);
    const auto *msp_148 = buffer.data(msp + 148);
    const auto *msp_149 = buffer.data(msp + 149);
    const auto *msp_150 = buffer.data(msp + 150);
    const auto *msp_151 = buffer.data(msp + 151);
    const auto *msp_152 = buffer.data(msp + 152);
    const auto *msp_153 = buffer.data(msp + 153);
    const auto *msp_154 = buffer.data(msp + 154);
    const auto *msp_155 = buffer.data(msp + 155);
    const auto *msp_156 = buffer.data(msp + 156);
    const auto *msp_157 = buffer.data(msp + 157);
    const auto *msp_158 = buffer.data(msp + 158);
    const auto *msp_160 = buffer.data(msp + 160);
    const auto *msp_161 = buffer.data(msp + 161);
    const auto *msp_162 = buffer.data(msp + 162);
    const auto *msp_164 = buffer.data(msp + 164);

    const auto *msd1_210 = buffer.data(msd1 + 210);
    const auto *msd1_215 = buffer.data(msd1 + 215);
    const auto *msd1_216 = buffer.data(msd1 + 216);
    const auto *msd1_264 = buffer.data(msd1 + 264);
    const auto *msd1_270 = buffer.data(msd1 + 270);
    const auto *msd1_273 = buffer.data(msd1 + 273);
    const auto *msd1_275 = buffer.data(msd1 + 275);
    const auto *msd1_279 = buffer.data(msd1 + 279);
    const auto *msd1_281 = buffer.data(msd1 + 281);
    const auto *msd1_282 = buffer.data(msd1 + 282);
    const auto *msd1_285 = buffer.data(msd1 + 285);
    const auto *msd1_287 = buffer.data(msd1 + 287);
    const auto *msd1_288 = buffer.data(msd1 + 288);
    const auto *msd1_291 = buffer.data(msd1 + 291);
    const auto *msd1_293 = buffer.data(msd1 + 293);
    const auto *msd1_294 = buffer.data(msd1 + 294);
    const auto *msd1_297 = buffer.data(msd1 + 297);
    const auto *msd1_299 = buffer.data(msd1 + 299);
    const auto *msd1_300 = buffer.data(msd1 + 300);
    const auto *msd1_303 = buffer.data(msd1 + 303);
    const auto *msd1_305 = buffer.data(msd1 + 305);
    const auto *msd1_306 = buffer.data(msd1 + 306);
    const auto *msd1_309 = buffer.data(msd1 + 309);
    const auto *msd1_311 = buffer.data(msd1 + 311);
    const auto *msd1_312 = buffer.data(msd1 + 312);
    const auto *msd1_315 = buffer.data(msd1 + 315);
    const auto *msd1_317 = buffer.data(msd1 + 317);
    const auto *msd1_321 = buffer.data(msd1 + 321);
    const auto *msd1_323 = buffer.data(msd1 + 323);
    const auto *msd1_324 = buffer.data(msd1 + 324);
    const auto *msd1_327 = buffer.data(msd1 + 327);
    const auto *msd1_329 = buffer.data(msd1 + 329);

    const auto *nss0_40 = buffer.data(nss0 + 40);
    const auto *nss0_41 = buffer.data(nss0 + 41);
    const auto *nss0_42 = buffer.data(nss0 + 42);
    const auto *nss0_43 = buffer.data(nss0 + 43);
    const auto *nss0_44 = buffer.data(nss0 + 44);
    const auto *nss0_55 = buffer.data(nss0 + 55);
    const auto *nss0_56 = buffer.data(nss0 + 56);
    const auto *nss0_57 = buffer.data(nss0 + 57);
    const auto *nss0_58 = buffer.data(nss0 + 58);
    const auto *nss0_59 = buffer.data(nss0 + 59);
    const auto *nss0_60 = buffer.data(nss0 + 60);
    const auto *nss0_61 = buffer.data(nss0 + 61);

    const auto *nss1_40 = buffer.data(nss1 + 40);
    const auto *nss1_41 = buffer.data(nss1 + 41);
    const auto *nss1_42 = buffer.data(nss1 + 42);
    const auto *nss1_43 = buffer.data(nss1 + 43);
    const auto *nss1_44 = buffer.data(nss1 + 44);
    const auto *nss1_55 = buffer.data(nss1 + 55);
    const auto *nss1_56 = buffer.data(nss1 + 56);
    const auto *nss1_57 = buffer.data(nss1 + 57);
    const auto *nss1_58 = buffer.data(nss1 + 58);
    const auto *nss1_59 = buffer.data(nss1 + 59);
    const auto *nss1_60 = buffer.data(nss1 + 60);
    const auto *nss1_61 = buffer.data(nss1 + 61);

    const auto *nsp_121 = buffer.data(nsp + 121);
    const auto *nsp_122 = buffer.data(nsp + 122);
    const auto *nsp_123 = buffer.data(nsp + 123);
    const auto *nsp_124 = buffer.data(nsp + 124);
    const auto *nsp_125 = buffer.data(nsp + 125);
    const auto *nsp_126 = buffer.data(nsp + 126);
    const auto *nsp_127 = buffer.data(nsp + 127);
    const auto *nsp_128 = buffer.data(nsp + 128);
    const auto *nsp_130 = buffer.data(nsp + 130);
    const auto *nsp_131 = buffer.data(nsp + 131);
    const auto *nsp_132 = buffer.data(nsp + 132);
    const auto *nsp_133 = buffer.data(nsp + 133);
    const auto *nsp_134 = buffer.data(nsp + 134);
    const auto *nsp_135 = buffer.data(nsp + 135);
    const auto *nsp_136 = buffer.data(nsp + 136);
    const auto *nsp_139 = buffer.data(nsp + 139);
    const auto *nsp_140 = buffer.data(nsp + 140);
    const auto *nsp_142 = buffer.data(nsp + 142);
    const auto *nsp_143 = buffer.data(nsp + 143);
    const auto *nsp_145 = buffer.data(nsp + 145);
    const auto *nsp_146 = buffer.data(nsp + 146);
    const auto *nsp_148 = buffer.data(nsp + 148);
    const auto *nsp_149 = buffer.data(nsp + 149);
    const auto *nsp_151 = buffer.data(nsp + 151);
    const auto *nsp_152 = buffer.data(nsp + 152);
    const auto *nsp_154 = buffer.data(nsp + 154);
    const auto *nsp_155 = buffer.data(nsp + 155);
    const auto *nsp_157 = buffer.data(nsp + 157);
    const auto *nsp_158 = buffer.data(nsp + 158);
    const auto *nsp_160 = buffer.data(nsp + 160);
    const auto *nsp_161 = buffer.data(nsp + 161);
    const auto *nsp_162 = buffer.data(nsp + 162);
    const auto *nsp_164 = buffer.data(nsp + 164);
    const auto *nsp_165 = buffer.data(nsp + 165);
    const auto *nsp_166 = buffer.data(nsp + 166);
    const auto *nsp_167 = buffer.data(nsp + 167);
    const auto *nsp_169 = buffer.data(nsp + 169);
    const auto *nsp_170 = buffer.data(nsp + 170);
    const auto *nsp_171 = buffer.data(nsp + 171);
    const auto *nsp_172 = buffer.data(nsp + 172);
    const auto *nsp_173 = buffer.data(nsp + 173);
    const auto *nsp_174 = buffer.data(nsp + 174);
    const auto *nsp_175 = buffer.data(nsp + 175);
    const auto *nsp_176 = buffer.data(nsp + 176);
    const auto *nsp_177 = buffer.data(nsp + 177);
    const auto *nsp_178 = buffer.data(nsp + 178);
    const auto *nsp_179 = buffer.data(nsp + 179);
    const auto *nsp_180 = buffer.data(nsp + 180);
    const auto *nsp_181 = buffer.data(nsp + 181);
    const auto *nsp_182 = buffer.data(nsp + 182);
    const auto *nsp_183 = buffer.data(nsp + 183);
    const auto *nsp_184 = buffer.data(nsp + 184);
    const auto *nsp_185 = buffer.data(nsp + 185);

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pc_x, pc_y, pc_z, msp_95, msp_97, msp_98, \
                         msp_122, nss0_40, nss1_40, nsp_121, nsp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_8 * msp_122[k]
                   + f_3 * pc_x[k] * nsp_122[k];

        t_243[k] = f_12 * msp_97[k]
                   + f_1 * nss0_40[k]
                   - f_2 * nss1_40[k]
                   + f_3 * pc_y[k] * nsp_121[k];

        t_244[k] = f_12 * msp_98[k]
                   + f_3 * pc_y[k] * nsp_122[k];

        t_245[k] = f_12 * msp_95[k]
                   + f_1 * nss0_40[k]
                   - f_2 * nss1_40[k]
                   + f_3 * pc_z[k] * nsp_122[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, pc_y, msp_100, msp_123, msp_124, \
                         msp_125, nss0_41, nss1_41, nsp_123, nsp_124, \
                         nsp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_8 * msp_123[k]
                   + f_1 * nss0_41[k]
                   - f_2 * nss1_41[k]
                   + f_3 * pc_x[k] * nsp_123[k];

        t_247[k] = f_8 * msp_124[k]
                   + f_3 * pc_x[k] * nsp_124[k];

        t_248[k] = f_8 * msp_125[k]
                   + f_3 * pc_x[k] * nsp_125[k];

        t_249[k] = f_10 * msp_100[k]
                   + f_1 * nss0_41[k]
                   - f_2 * nss1_41[k]
                   + f_3 * pc_y[k] * nsp_124[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, msp_98, msp_101, msp_126, \
                         nss0_41, nss0_42, nss1_41, nss1_42, nsp_125, \
                         nsp_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_10 * msp_101[k]
                   + f_3 * pc_y[k] * nsp_125[k];

        t_251[k] = f_13 * msp_98[k]
                   + f_1 * nss0_41[k]
                   - f_2 * nss1_41[k]
                   + f_3 * pc_z[k] * nsp_125[k];

        t_252[k] = f_8 * msp_126[k]
                   + f_1 * nss0_42[k]
                   - f_2 * nss1_42[k]
                   + f_3 * pc_x[k] * nsp_126[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, msp_103, msp_104, msp_127, \
                         msp_128, nss0_42, nss1_42, nsp_127, nsp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_8 * msp_127[k]
                   + f_3 * pc_x[k] * nsp_127[k];

        t_254[k] = f_8 * msp_128[k]
                   + f_3 * pc_x[k] * nsp_128[k];

        t_255[k] = f_8 * msp_103[k]
                   + f_1 * nss0_42[k]
                   - f_2 * nss1_42[k]
                   + f_3 * pc_y[k] * nsp_127[k];

        t_256[k] = f_8 * msp_104[k]
                   + f_3 * pc_y[k] * nsp_128[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_y, pc_x, pc_y, pc_z, msd0_210, msp_101, \
                         msp_130, msd1_210, nss0_42, nss1_42, nsp_128, \
                         nsp_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_11 * msp_101[k]
                   + f_1 * nss0_42[k]
                   - f_2 * nss1_42[k]
                   + f_3 * pc_z[k] * nsp_128[k];

        t_258[k] = pa_y[k] * msd0_210[k]
                   - f_4 * pc_y[k] * msd1_210[k];

        t_259[k] = f_8 * msp_130[k]
                   + f_3 * pc_x[k] * nsp_130[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, msd0_215, msp_106, \
                         msp_107, msp_131, msd1_215, nss0_43, nss1_43, nsp_130, \
                         nsp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_8 * msp_131[k]
                   + f_3 * pc_x[k] * nsp_131[k];

        t_261[k] = f_6 * msp_106[k]
                   + f_1 * nss0_43[k]
                   - f_2 * nss1_43[k]
                   + f_3 * pc_y[k] * nsp_130[k];

        t_262[k] = f_6 * msp_107[k]
                   + f_3 * pc_y[k] * nsp_131[k];

        t_263[k] = pa_y[k] * msd0_215[k]
                   - f_4 * pc_y[k] * msd1_215[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, pc_x, pc_y, msp_132, msp_134, \
                         nss0_44, nss1_44, nsp_132, nsp_133, nsp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_8 * msp_132[k]
                   + f_1 * nss0_44[k]
                   - f_2 * nss1_44[k]
                   + f_3 * pc_x[k] * nsp_132[k];

        t_265[k] = f_3 * pc_y[k] * nsp_132[k];

        t_266[k] = f_8 * msp_134[k]
                   + f_3 * pc_x[k] * nsp_134[k];

        t_267[k] = f_1 * nss0_44[k]
                   - f_2 * nss1_44[k]
                   + f_3 * pc_y[k] * nsp_133[k];

        t_268[k] = f_3 * pc_y[k] * nsp_134[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_x, pc_x, pc_z, msd0_270, msp_107, msp_135, \
                         msp_136, msd1_270, nss0_44, nss1_44, nsp_134, \
                         nsp_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_7 * msp_107[k]
                   + f_1 * nss0_44[k]
                   - f_2 * nss1_44[k]
                   + f_3 * pc_z[k] * nsp_134[k];

        t_270[k] = pa_x[k] * msd0_270[k]
                   + f_8 * msp_135[k]
                   - f_4 * pc_x[k] * msd1_270[k];

        t_271[k] = f_6 * msp_136[k]
                   + f_3 * pc_x[k] * nsp_136[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_x, pc_x, pc_z, msd0_273, msd0_275, \
                         msd1_273, msd1_275, nsp_135, nsp_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_3 * pc_z[k] * nsp_135[k];

        t_273[k] = pa_x[k] * msd0_273[k]
                   - f_4 * pc_x[k] * msd1_273[k];

        t_274[k] = f_3 * pc_z[k] * nsp_136[k];

        t_275[k] = pa_x[k] * msd0_275[k]
                   - f_4 * pc_x[k] * msd1_275[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_x, pa_z, pc_x, pc_z, msd0_216, \
                         msd0_279, msp_139, msp_140, msd1_216, msd1_279, nsp_139, \
                         nsp_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pa_z[k] * msd0_216[k]
                   - f_4 * pc_z[k] * msd1_216[k];

        t_277[k] = f_6 * msp_139[k]
                   + f_3 * pc_x[k] * nsp_139[k];

        t_278[k] = f_6 * msp_140[k]
                   + f_3 * pc_x[k] * nsp_140[k];

        t_279[k] = pa_x[k] * msd0_279[k]
                   - f_4 * pc_x[k] * msd1_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_x, pc_x, pc_y, msd0_281, msd0_282, \
                         msp_113, msp_141, msp_142, msd1_281, msd1_282, nsp_140, \
                         nsp_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_7 * msp_113[k]
                   + f_3 * pc_y[k] * nsp_140[k];

        t_281[k] = pa_x[k] * msd0_281[k]
                   - f_4 * pc_x[k] * msd1_281[k];

        t_282[k] = pa_x[k] * msd0_282[k]
                   + f_8 * msp_141[k]
                   - f_4 * pc_x[k] * msd1_282[k];

        t_283[k] = f_6 * msp_142[k]
                   + f_3 * pc_x[k] * nsp_142[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pa_x, pc_x, pc_y, msd0_285, msd0_287, \
                         msp_116, msp_143, msd1_285, msd1_287, \
                         nsp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_6 * msp_143[k]
                   + f_3 * pc_x[k] * nsp_143[k];

        t_285[k] = pa_x[k] * msd0_285[k]
                   - f_4 * pc_x[k] * msd1_285[k];

        t_286[k] = f_9 * msp_116[k]
                   + f_3 * pc_y[k] * nsp_143[k];

        t_287[k] = pa_x[k] * msd0_287[k]
                   - f_4 * pc_x[k] * msd1_287[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pa_x, pc_x, msd0_288, msd0_291, msp_144, \
                         msp_145, msp_146, msd1_288, msd1_291, nsp_145, \
                         nsp_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = pa_x[k] * msd0_288[k]
                   + f_8 * msp_144[k]
                   - f_4 * pc_x[k] * msd1_288[k];

        t_289[k] = f_6 * msp_145[k]
                   + f_3 * pc_x[k] * nsp_145[k];

        t_290[k] = f_6 * msp_146[k]
                   + f_3 * pc_x[k] * nsp_146[k];

        t_291[k] = pa_x[k] * msd0_291[k]
                   - f_4 * pc_x[k] * msd1_291[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_x, pc_x, pc_y, msd0_293, msd0_294, \
                         msp_119, msp_147, msp_148, msd1_293, msd1_294, nsp_146, \
                         nsp_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_11 * msp_119[k]
                   + f_3 * pc_y[k] * nsp_146[k];

        t_293[k] = pa_x[k] * msd0_293[k]
                   - f_4 * pc_x[k] * msd1_293[k];

        t_294[k] = pa_x[k] * msd0_294[k]
                   + f_8 * msp_147[k]
                   - f_4 * pc_x[k] * msd1_294[k];

        t_295[k] = f_6 * msp_148[k]
                   + f_3 * pc_x[k] * nsp_148[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_x, pc_x, pc_y, msd0_297, msd0_299, \
                         msp_122, msp_149, msd1_297, msd1_299, \
                         nsp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_6 * msp_149[k]
                   + f_3 * pc_x[k] * nsp_149[k];

        t_297[k] = pa_x[k] * msd0_297[k]
                   - f_4 * pc_x[k] * msd1_297[k];

        t_298[k] = f_13 * msp_122[k]
                   + f_3 * pc_y[k] * nsp_149[k];

        t_299[k] = pa_x[k] * msd0_299[k]
                   - f_4 * pc_x[k] * msd1_299[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_x, pc_x, msd0_300, msd0_303, msp_150, \
                         msp_151, msp_152, msd1_300, msd1_303, nsp_151, \
                         nsp_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = pa_x[k] * msd0_300[k]
                   + f_8 * msp_150[k]
                   - f_4 * pc_x[k] * msd1_300[k];

        t_301[k] = f_6 * msp_151[k]
                   + f_3 * pc_x[k] * nsp_151[k];

        t_302[k] = f_6 * msp_152[k]
                   + f_3 * pc_x[k] * nsp_152[k];

        t_303[k] = pa_x[k] * msd0_303[k]
                   - f_4 * pc_x[k] * msd1_303[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pa_x, pc_x, pc_y, msd0_305, msd0_306, \
                         msp_125, msp_153, msp_154, msd1_305, msd1_306, nsp_152, \
                         nsp_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_12 * msp_125[k]
                   + f_3 * pc_y[k] * nsp_152[k];

        t_305[k] = pa_x[k] * msd0_305[k]
                   - f_4 * pc_x[k] * msd1_305[k];

        t_306[k] = pa_x[k] * msd0_306[k]
                   + f_8 * msp_153[k]
                   - f_4 * pc_x[k] * msd1_306[k];

        t_307[k] = f_6 * msp_154[k]
                   + f_3 * pc_x[k] * nsp_154[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_x, pc_x, pc_y, msd0_309, msd0_311, \
                         msp_128, msp_155, msd1_309, msd1_311, \
                         nsp_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_6 * msp_155[k]
                   + f_3 * pc_x[k] * nsp_155[k];

        t_309[k] = pa_x[k] * msd0_309[k]
                   - f_4 * pc_x[k] * msd1_309[k];

        t_310[k] = f_10 * msp_128[k]
                   + f_3 * pc_y[k] * nsp_155[k];

        t_311[k] = pa_x[k] * msd0_311[k]
                   - f_4 * pc_x[k] * msd1_311[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pa_x, pc_x, msd0_312, msd0_315, msp_156, \
                         msp_157, msp_158, msd1_312, msd1_315, nsp_157, \
                         nsp_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = pa_x[k] * msd0_312[k]
                   + f_8 * msp_156[k]
                   - f_4 * pc_x[k] * msd1_312[k];

        t_313[k] = f_6 * msp_157[k]
                   + f_3 * pc_x[k] * nsp_157[k];

        t_314[k] = f_6 * msp_158[k]
                   + f_3 * pc_x[k] * nsp_158[k];

        t_315[k] = pa_x[k] * msd0_315[k]
                   - f_4 * pc_x[k] * msd1_315[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_x, pa_y, pc_x, pc_y, msd0_264, \
                         msd0_317, msp_131, msp_160, msd1_264, msd1_317, nsp_158, \
                         nsp_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_8 * msp_131[k]
                   + f_3 * pc_y[k] * nsp_158[k];

        t_317[k] = pa_x[k] * msd0_317[k]
                   - f_4 * pc_x[k] * msd1_317[k];

        t_318[k] = pa_y[k] * msd0_264[k]
                   - f_4 * pc_y[k] * msd1_264[k];

        t_319[k] = f_6 * msp_160[k]
                   + f_3 * pc_x[k] * nsp_160[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_x, pc_x, pc_y, msd0_321, msd0_323, \
                         msp_134, msp_161, msd1_321, msd1_323, \
                         nsp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_6 * msp_161[k]
                   + f_3 * pc_x[k] * nsp_161[k];

        t_321[k] = pa_x[k] * msd0_321[k]
                   - f_4 * pc_x[k] * msd1_321[k];

        t_322[k] = f_6 * msp_134[k]
                   + f_3 * pc_y[k] * nsp_161[k];

        t_323[k] = pa_x[k] * msd0_323[k]
                   - f_4 * pc_x[k] * msd1_323[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pa_x, pc_x, pc_y, msd0_324, \
                         msd0_327, msp_162, msp_164, msd1_324, msd1_327, nsp_162, \
                         nsp_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = pa_x[k] * msd0_324[k]
                   + f_8 * msp_162[k]
                   - f_4 * pc_x[k] * msd1_324[k];

        t_325[k] = f_3 * pc_y[k] * nsp_162[k];

        t_326[k] = f_6 * msp_164[k]
                   + f_3 * pc_x[k] * nsp_164[k];

        t_327[k] = pa_x[k] * msd0_327[k]
                   - f_4 * pc_x[k] * msd1_327[k];

        t_328[k] = f_3 * pc_y[k] * nsp_164[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, t_333, pa_x, pc_x, pc_y, msd0_329, \
                         msp_136, msd1_329, nss0_55, nss1_55, nsp_165, nsp_166, \
                         nsp_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = pa_x[k] * msd0_329[k]
                   - f_4 * pc_x[k] * msd1_329[k];

        t_330[k] = f_1 * nss0_55[k]
                   - f_2 * nss1_55[k]
                   + f_3 * pc_x[k] * nsp_165[k];

        t_331[k] = f_3 * pc_x[k] * nsp_166[k];

        t_332[k] = f_3 * pc_x[k] * nsp_167[k];

        t_333[k] = f_0 * msp_136[k]
                   + f_1 * nss0_55[k]
                   - f_2 * nss1_55[k]
                   + f_3 * pc_y[k] * nsp_166[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, pa_z, pc_x, pc_z, msd0_270, \
                         msd1_270, nss0_55, nss1_55, nsp_166, nsp_167, nsp_169, \
                         nsp_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_3 * pc_z[k] * nsp_166[k];

        t_335[k] = f_1 * nss0_55[k]
                   - f_2 * nss1_55[k]
                   + f_3 * pc_z[k] * nsp_167[k];

        t_336[k] = pa_z[k] * msd0_270[k]
                   - f_4 * pc_z[k] * msd1_270[k];

        t_337[k] = f_3 * pc_x[k] * nsp_169[k];

        t_338[k] = f_3 * pc_x[k] * nsp_170[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_z, pc_y, pc_z, msd0_273, msp_137, msp_140, \
                         msd1_273, nss0_56, nss1_56, nsp_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * msd0_273[k]
                   - f_4 * pc_z[k] * msd1_273[k];

        t_340[k] = f_5 * msp_140[k]
                   + f_3 * pc_y[k] * nsp_170[k];

        t_341[k] = f_6 * msp_137[k]
                   + f_1 * nss0_56[k]
                   - f_2 * nss1_56[k]
                   + f_3 * pc_z[k] * nsp_170[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, pc_x, pc_y, msp_142, msp_143, \
                         nss0_57, nss1_57, nsp_171, nsp_172, nsp_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_1 * nss0_57[k]
                   - f_2 * nss1_57[k]
                   + f_3 * pc_x[k] * nsp_171[k];

        t_343[k] = f_3 * pc_x[k] * nsp_172[k];

        t_344[k] = f_3 * pc_x[k] * nsp_173[k];

        t_345[k] = f_7 * msp_142[k]
                   + f_1 * nss0_57[k]
                   - f_2 * nss1_57[k]
                   + f_3 * pc_y[k] * nsp_172[k];

        t_346[k] = f_7 * msp_143[k]
                   + f_3 * pc_y[k] * nsp_173[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, pc_x, pc_z, msp_140, nss0_57, nss0_58, \
                         nss1_57, nss1_58, nsp_173, nsp_174, nsp_175, \
                         nsp_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_8 * msp_140[k]
                   + f_1 * nss0_57[k]
                   - f_2 * nss1_57[k]
                   + f_3 * pc_z[k] * nsp_173[k];

        t_348[k] = f_1 * nss0_58[k]
                   - f_2 * nss1_58[k]
                   + f_3 * pc_x[k] * nsp_174[k];

        t_349[k] = f_3 * pc_x[k] * nsp_175[k];

        t_350[k] = f_3 * pc_x[k] * nsp_176[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pc_y, pc_z, msp_143, msp_145, msp_146, nss0_58, \
                         nss1_58, nsp_175, nsp_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_9 * msp_145[k]
                   + f_1 * nss0_58[k]
                   - f_2 * nss1_58[k]
                   + f_3 * pc_y[k] * nsp_175[k];

        t_352[k] = f_9 * msp_146[k]
                   + f_3 * pc_y[k] * nsp_176[k];

        t_353[k] = f_10 * msp_143[k]
                   + f_1 * nss0_58[k]
                   - f_2 * nss1_58[k]
                   + f_3 * pc_z[k] * nsp_176[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, t_358, pc_x, pc_y, msp_148, msp_149, \
                         nss0_59, nss1_59, nsp_177, nsp_178, nsp_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_1 * nss0_59[k]
                   - f_2 * nss1_59[k]
                   + f_3 * pc_x[k] * nsp_177[k];

        t_355[k] = f_3 * pc_x[k] * nsp_178[k];

        t_356[k] = f_3 * pc_x[k] * nsp_179[k];

        t_357[k] = f_11 * msp_148[k]
                   + f_1 * nss0_59[k]
                   - f_2 * nss1_59[k]
                   + f_3 * pc_y[k] * nsp_178[k];

        t_358[k] = f_11 * msp_149[k]
                   + f_3 * pc_y[k] * nsp_179[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pc_x, pc_z, msp_146, nss0_59, nss0_60, \
                         nss1_59, nss1_60, nsp_179, nsp_180, nsp_181, \
                         nsp_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * msp_146[k]
                   + f_1 * nss0_59[k]
                   - f_2 * nss1_59[k]
                   + f_3 * pc_z[k] * nsp_179[k];

        t_360[k] = f_1 * nss0_60[k]
                   - f_2 * nss1_60[k]
                   + f_3 * pc_x[k] * nsp_180[k];

        t_361[k] = f_3 * pc_x[k] * nsp_181[k];

        t_362[k] = f_3 * pc_x[k] * nsp_182[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pc_y, pc_z, msp_149, msp_151, msp_152, nss0_60, \
                         nss1_60, nsp_181, nsp_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_13 * msp_151[k]
                   + f_1 * nss0_60[k]
                   - f_2 * nss1_60[k]
                   + f_3 * pc_y[k] * nsp_181[k];

        t_364[k] = f_13 * msp_152[k]
                   + f_3 * pc_y[k] * nsp_182[k];

        t_365[k] = f_13 * msp_149[k]
                   + f_1 * nss0_60[k]
                   - f_2 * nss1_60[k]
                   + f_3 * pc_z[k] * nsp_182[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, t_370, pc_x, pc_y, msp_154, msp_155, \
                         nss0_61, nss1_61, nsp_183, nsp_184, nsp_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_1 * nss0_61[k]
                   - f_2 * nss1_61[k]
                   + f_3 * pc_x[k] * nsp_183[k];

        t_367[k] = f_3 * pc_x[k] * nsp_184[k];

        t_368[k] = f_3 * pc_x[k] * nsp_185[k];

        t_369[k] = f_12 * msp_154[k]
                   + f_1 * nss0_61[k]
                   - f_2 * nss1_61[k]
                   + f_3 * pc_y[k] * nsp_184[k];

        t_370[k] = f_12 * msp_155[k]
                   + f_3 * pc_y[k] * nsp_185[k];
    }
}

static auto
compute_prim_nsd_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msd0,
                                                          const size_t msp, const size_t msd1,
                                                          const size_t nss0, const size_t nss1,
                                                          const size_t nsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msd0_324 = buffer.data(msd0 + 324);
    const auto *msd0_327 = buffer.data(msd0 + 327);
    const auto *msd0_329 = buffer.data(msd0 + 329);

    const auto *msp_152 = buffer.data(msp + 152);
    const auto *msp_155 = buffer.data(msp + 155);
    const auto *msp_157 = buffer.data(msp + 157);
    const auto *msp_158 = buffer.data(msp + 158);
    const auto *msp_160 = buffer.data(msp + 160);
    const auto *msp_161 = buffer.data(msp + 161);
    const auto *msp_163 = buffer.data(msp + 163);
    const auto *msp_164 = buffer.data(msp + 164);

    const auto *msd1_324 = buffer.data(msd1 + 324);
    const auto *msd1_327 = buffer.data(msd1 + 327);
    const auto *msd1_329 = buffer.data(msd1 + 329);

    const auto *nss0_61 = buffer.data(nss0 + 61);
    const auto *nss0_62 = buffer.data(nss0 + 62);
    const auto *nss0_63 = buffer.data(nss0 + 63);
    const auto *nss0_65 = buffer.data(nss0 + 65);

    const auto *nss1_61 = buffer.data(nss1 + 61);
    const auto *nss1_62 = buffer.data(nss1 + 62);
    const auto *nss1_63 = buffer.data(nss1 + 63);
    const auto *nss1_65 = buffer.data(nss1 + 65);

    const auto *nsp_185 = buffer.data(nsp + 185);
    const auto *nsp_186 = buffer.data(nsp + 186);
    const auto *nsp_187 = buffer.data(nsp + 187);
    const auto *nsp_188 = buffer.data(nsp + 188);
    const auto *nsp_189 = buffer.data(nsp + 189);
    const auto *nsp_190 = buffer.data(nsp + 190);
    const auto *nsp_191 = buffer.data(nsp + 191);
    const auto *nsp_193 = buffer.data(nsp + 193);
    const auto *nsp_194 = buffer.data(nsp + 194);
    const auto *nsp_195 = buffer.data(nsp + 195);
    const auto *nsp_196 = buffer.data(nsp + 196);
    const auto *nsp_197 = buffer.data(nsp + 197);

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pc_x, pc_z, msp_152, nss0_61, nss0_62, \
                         nss1_61, nss1_62, nsp_185, nsp_186, nsp_187, \
                         nsp_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_11 * msp_152[k]
                   + f_1 * nss0_61[k]
                   - f_2 * nss1_61[k]
                   + f_3 * pc_z[k] * nsp_185[k];

        t_372[k] = f_1 * nss0_62[k]
                   - f_2 * nss1_62[k]
                   + f_3 * pc_x[k] * nsp_186[k];

        t_373[k] = f_3 * pc_x[k] * nsp_187[k];

        t_374[k] = f_3 * pc_x[k] * nsp_188[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_y, pc_z, msp_155, msp_157, msp_158, nss0_62, \
                         nss1_62, nsp_187, nsp_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_10 * msp_157[k]
                   + f_1 * nss0_62[k]
                   - f_2 * nss1_62[k]
                   + f_3 * pc_y[k] * nsp_187[k];

        t_376[k] = f_10 * msp_158[k]
                   + f_3 * pc_y[k] * nsp_188[k];

        t_377[k] = f_9 * msp_155[k]
                   + f_1 * nss0_62[k]
                   - f_2 * nss1_62[k]
                   + f_3 * pc_z[k] * nsp_188[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, pc_x, pc_y, msp_160, msp_161, \
                         nss0_63, nss1_63, nsp_189, nsp_190, nsp_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_1 * nss0_63[k]
                   - f_2 * nss1_63[k]
                   + f_3 * pc_x[k] * nsp_189[k];

        t_379[k] = f_3 * pc_x[k] * nsp_190[k];

        t_380[k] = f_3 * pc_x[k] * nsp_191[k];

        t_381[k] = f_8 * msp_160[k]
                   + f_1 * nss0_63[k]
                   - f_2 * nss1_63[k]
                   + f_3 * pc_y[k] * nsp_190[k];

        t_382[k] = f_8 * msp_161[k]
                   + f_3 * pc_y[k] * nsp_191[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, pa_y, pc_x, pc_y, pc_z, msd0_324, \
                         msp_158, msd1_324, nss0_63, nss1_63, nsp_191, nsp_193, \
                         nsp_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_7 * msp_158[k]
                   + f_1 * nss0_63[k]
                   - f_2 * nss1_63[k]
                   + f_3 * pc_z[k] * nsp_191[k];

        t_384[k] = pa_y[k] * msd0_324[k]
                   - f_4 * pc_y[k] * msd1_324[k];

        t_385[k] = f_3 * pc_x[k] * nsp_193[k];

        t_386[k] = f_3 * pc_x[k] * nsp_194[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, pa_y, pc_y, msd0_327, msd0_329, msp_163, \
                         msp_164, msd1_327, msd1_329, nsp_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = pa_y[k] * msd0_327[k]
                   + f_8 * msp_163[k]
                   - f_4 * pc_y[k] * msd1_327[k];

        t_388[k] = f_6 * msp_164[k]
                   + f_3 * pc_y[k] * nsp_194[k];

        t_389[k] = pa_y[k] * msd0_329[k]
                   - f_4 * pc_y[k] * msd1_329[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, t_395, pc_x, pc_y, pc_z, msp_164, \
                         nss0_65, nss1_65, nsp_195, nsp_196, nsp_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_1 * nss0_65[k]
                   - f_2 * nss1_65[k]
                   + f_3 * pc_x[k] * nsp_195[k];

        t_391[k] = f_3 * pc_x[k] * nsp_196[k];

        t_392[k] = f_3 * pc_x[k] * nsp_197[k];

        t_393[k] = f_1 * nss0_65[k]
                   - f_2 * nss1_65[k]
                   + f_3 * pc_y[k] * nsp_196[k];

        t_394[k] = f_3 * pc_y[k] * nsp_197[k];

        t_395[k] = f_0 * msp_164[k]
                   + f_1 * nss0_65[k]
                   - f_2 * nss1_65[k]
                   + f_3 * pc_z[k] * nsp_197[k];
    }
}

auto
compute_prim_nsd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t msd0, const size_t msp,
                                                   const size_t msd1, const size_t nss0,
                                                   const size_t nss1, const size_t nsp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_nsd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, msd0, msp,
                                                              msd1, nss0, nss1, nsp, ncols,
                                                              gamma, p, q);

    compute_prim_nsd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, msd0, msp,
                                                              msd1, nss0, nss1, nsp, ncols,
                                                              gamma, p, q);

    compute_prim_nsd_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, msd0, msp,
                                                              msd1, nss0, nss1, nsp, ncols,
                                                              gamma, p, q);

    compute_prim_nsd_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, msd0, msp,
                                                              msd1, nss0, nss1, nsp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
