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


#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_psk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ssk0, const size_t ssi,
                                                   const size_t ssk1, const size_t psh0,
                                                   const size_t psh1, const size_t psi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 0.5 / gamma;
    const auto f_4 = 0.5 * p / (gamma * q);
    const auto f_5 = 1.0 / gamma;
    const auto f_6 = p / (gamma * q);
    const auto f_7 = 1.5 / gamma;
    const auto f_8 = 1.5 * p / (gamma * q);
    const auto f_9 = 2.0 / gamma;
    const auto f_10 = 2.0 * p / (gamma * q);
    const auto f_11 = 0.5 / q;
    const auto f_12 = 2.5 / gamma;
    const auto f_13 = 2.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssk0_0 = buffer.data(ssk0 + 0);
    const auto *ssk0_3 = buffer.data(ssk0 + 3);
    const auto *ssk0_5 = buffer.data(ssk0 + 5);
    const auto *ssk0_6 = buffer.data(ssk0 + 6);
    const auto *ssk0_9 = buffer.data(ssk0 + 9);
    const auto *ssk0_10 = buffer.data(ssk0 + 10);
    const auto *ssk0_14 = buffer.data(ssk0 + 14);
    const auto *ssk0_15 = buffer.data(ssk0 + 15);
    const auto *ssk0_20 = buffer.data(ssk0 + 20);
    const auto *ssk0_28 = buffer.data(ssk0 + 28);
    const auto *ssk0_30 = buffer.data(ssk0 + 30);
    const auto *ssk0_31 = buffer.data(ssk0 + 31);
    const auto *ssk0_32 = buffer.data(ssk0 + 32);
    const auto *ssk0_33 = buffer.data(ssk0 + 33);
    const auto *ssk0_35 = buffer.data(ssk0 + 35);

    const auto *ssi_0 = buffer.data(ssi + 0);
    const auto *ssi_21 = buffer.data(ssi + 21);
    const auto *ssi_23 = buffer.data(ssi + 23);
    const auto *ssi_24 = buffer.data(ssi + 24);
    const auto *ssi_25 = buffer.data(ssi + 25);
    const auto *ssi_27 = buffer.data(ssi + 27);

    const auto *ssk1_0 = buffer.data(ssk1 + 0);
    const auto *ssk1_3 = buffer.data(ssk1 + 3);
    const auto *ssk1_5 = buffer.data(ssk1 + 5);
    const auto *ssk1_6 = buffer.data(ssk1 + 6);
    const auto *ssk1_9 = buffer.data(ssk1 + 9);
    const auto *ssk1_10 = buffer.data(ssk1 + 10);
    const auto *ssk1_14 = buffer.data(ssk1 + 14);
    const auto *ssk1_15 = buffer.data(ssk1 + 15);
    const auto *ssk1_20 = buffer.data(ssk1 + 20);
    const auto *ssk1_28 = buffer.data(ssk1 + 28);
    const auto *ssk1_30 = buffer.data(ssk1 + 30);
    const auto *ssk1_31 = buffer.data(ssk1 + 31);
    const auto *ssk1_32 = buffer.data(ssk1 + 32);
    const auto *ssk1_33 = buffer.data(ssk1 + 33);
    const auto *ssk1_35 = buffer.data(ssk1 + 35);

    const auto *psh0_0 = buffer.data(psh0 + 0);
    const auto *psh0_1 = buffer.data(psh0 + 1);
    const auto *psh0_2 = buffer.data(psh0 + 2);
    const auto *psh0_3 = buffer.data(psh0 + 3);
    const auto *psh0_5 = buffer.data(psh0 + 5);
    const auto *psh0_6 = buffer.data(psh0 + 6);
    const auto *psh0_8 = buffer.data(psh0 + 8);
    const auto *psh0_9 = buffer.data(psh0 + 9);
    const auto *psh0_22 = buffer.data(psh0 + 22);
    const auto *psh0_24 = buffer.data(psh0 + 24);
    const auto *psh0_27 = buffer.data(psh0 + 27);
    const auto *psh0_29 = buffer.data(psh0 + 29);
    const auto *psh0_31 = buffer.data(psh0 + 31);
    const auto *psh0_33 = buffer.data(psh0 + 33);
    const auto *psh0_34 = buffer.data(psh0 + 34);
    const auto *psh0_36 = buffer.data(psh0 + 36);
    const auto *psh0_37 = buffer.data(psh0 + 37);
    const auto *psh0_38 = buffer.data(psh0 + 38);
    const auto *psh0_39 = buffer.data(psh0 + 39);
    const auto *psh0_40 = buffer.data(psh0 + 40);
    const auto *psh0_44 = buffer.data(psh0 + 44);
    const auto *psh0_47 = buffer.data(psh0 + 47);
    const auto *psh0_49 = buffer.data(psh0 + 49);
    const auto *psh0_51 = buffer.data(psh0 + 51);
    const auto *psh0_53 = buffer.data(psh0 + 53);
    const auto *psh0_54 = buffer.data(psh0 + 54);
    const auto *psh0_56 = buffer.data(psh0 + 56);
    const auto *psh0_58 = buffer.data(psh0 + 58);
    const auto *psh0_59 = buffer.data(psh0 + 59);
    const auto *psh0_60 = buffer.data(psh0 + 60);
    const auto *psh0_61 = buffer.data(psh0 + 61);
    const auto *psh0_62 = buffer.data(psh0 + 62);

    const auto *psh1_0 = buffer.data(psh1 + 0);
    const auto *psh1_1 = buffer.data(psh1 + 1);
    const auto *psh1_2 = buffer.data(psh1 + 2);
    const auto *psh1_3 = buffer.data(psh1 + 3);
    const auto *psh1_5 = buffer.data(psh1 + 5);
    const auto *psh1_6 = buffer.data(psh1 + 6);
    const auto *psh1_8 = buffer.data(psh1 + 8);
    const auto *psh1_9 = buffer.data(psh1 + 9);
    const auto *psh1_22 = buffer.data(psh1 + 22);
    const auto *psh1_24 = buffer.data(psh1 + 24);
    const auto *psh1_27 = buffer.data(psh1 + 27);
    const auto *psh1_29 = buffer.data(psh1 + 29);
    const auto *psh1_31 = buffer.data(psh1 + 31);
    const auto *psh1_33 = buffer.data(psh1 + 33);
    const auto *psh1_34 = buffer.data(psh1 + 34);
    const auto *psh1_36 = buffer.data(psh1 + 36);
    const auto *psh1_37 = buffer.data(psh1 + 37);
    const auto *psh1_38 = buffer.data(psh1 + 38);
    const auto *psh1_39 = buffer.data(psh1 + 39);
    const auto *psh1_40 = buffer.data(psh1 + 40);
    const auto *psh1_44 = buffer.data(psh1 + 44);
    const auto *psh1_47 = buffer.data(psh1 + 47);
    const auto *psh1_49 = buffer.data(psh1 + 49);
    const auto *psh1_51 = buffer.data(psh1 + 51);
    const auto *psh1_53 = buffer.data(psh1 + 53);
    const auto *psh1_54 = buffer.data(psh1 + 54);
    const auto *psh1_56 = buffer.data(psh1 + 56);
    const auto *psh1_58 = buffer.data(psh1 + 58);
    const auto *psh1_59 = buffer.data(psh1 + 59);
    const auto *psh1_60 = buffer.data(psh1 + 60);
    const auto *psh1_61 = buffer.data(psh1 + 61);
    const auto *psh1_62 = buffer.data(psh1 + 62);

    const auto *psi_0 = buffer.data(psi + 0);
    const auto *psi_1 = buffer.data(psi + 1);
    const auto *psi_2 = buffer.data(psi + 2);
    const auto *psi_3 = buffer.data(psi + 3);
    const auto *psi_5 = buffer.data(psi + 5);
    const auto *psi_6 = buffer.data(psi + 6);
    const auto *psi_8 = buffer.data(psi + 8);
    const auto *psi_9 = buffer.data(psi + 9);
    const auto *psi_10 = buffer.data(psi + 10);
    const auto *psi_12 = buffer.data(psi + 12);
    const auto *psi_13 = buffer.data(psi + 13);
    const auto *psi_14 = buffer.data(psi + 14);
    const auto *psi_15 = buffer.data(psi + 15);
    const auto *psi_20 = buffer.data(psi + 20);
    const auto *psi_21 = buffer.data(psi + 21);
    const auto *psi_23 = buffer.data(psi + 23);
    const auto *psi_24 = buffer.data(psi + 24);
    const auto *psi_25 = buffer.data(psi + 25);
    const auto *psi_27 = buffer.data(psi + 27);
    const auto *psi_28 = buffer.data(psi + 28);
    const auto *psi_29 = buffer.data(psi + 29);
    const auto *psi_31 = buffer.data(psi + 31);
    const auto *psi_34 = buffer.data(psi + 34);
    const auto *psi_36 = buffer.data(psi + 36);
    const auto *psi_38 = buffer.data(psi + 38);
    const auto *psi_40 = buffer.data(psi + 40);
    const auto *psi_41 = buffer.data(psi + 41);
    const auto *psi_43 = buffer.data(psi + 43);
    const auto *psi_45 = buffer.data(psi + 45);
    const auto *psi_46 = buffer.data(psi + 46);
    const auto *psi_47 = buffer.data(psi + 47);
    const auto *psi_49 = buffer.data(psi + 49);
    const auto *psi_50 = buffer.data(psi + 50);
    const auto *psi_51 = buffer.data(psi + 51);
    const auto *psi_52 = buffer.data(psi + 52);
    const auto *psi_53 = buffer.data(psi + 53);
    const auto *psi_54 = buffer.data(psi + 54);
    const auto *psi_55 = buffer.data(psi + 55);
    const auto *psi_56 = buffer.data(psi + 56);
    const auto *psi_58 = buffer.data(psi + 58);
    const auto *psi_61 = buffer.data(psi + 61);
    const auto *psi_63 = buffer.data(psi + 63);
    const auto *psi_65 = buffer.data(psi + 65);
    const auto *psi_67 = buffer.data(psi + 67);
    const auto *psi_68 = buffer.data(psi + 68);
    const auto *psi_70 = buffer.data(psi + 70);
    const auto *psi_72 = buffer.data(psi + 72);
    const auto *psi_73 = buffer.data(psi + 73);
    const auto *psi_74 = buffer.data(psi + 74);
    const auto *psi_76 = buffer.data(psi + 76);
    const auto *psi_77 = buffer.data(psi + 77);
    const auto *psi_78 = buffer.data(psi + 78);
    const auto *psi_79 = buffer.data(psi + 79);
    const auto *psi_80 = buffer.data(psi + 80);
    const auto *psi_81 = buffer.data(psi + 81);
    const auto *psi_82 = buffer.data(psi + 82);
    const auto *psi_83 = buffer.data(psi + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pc_x, pc_y, pc_z, ssk0_0, ssi_0, ssk1_0, \
                         psh0_0, psh1_0, psi_0, psi_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pa_x[k] * ssk0_0[k]
                 + f_0 * ssi_0[k]
                 - f_1 * pc_x[k] * ssk1_0[k];

        t_1[k] = f_2 * pc_y[k] * psi_0[k];

        t_2[k] = f_2 * pc_z[k] * psi_0[k];

        t_3[k] = f_3 * psh0_0[k]
                 - f_4 * psh1_0[k]
                 + f_2 * pc_y[k] * psi_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_y, pc_z, psh0_0, psh0_1, psh1_0, psh1_1, \
                         psi_2, psi_3, psi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * psi_2[k];

        t_5[k] = f_3 * psh0_0[k]
                 - f_4 * psh1_0[k]
                 + f_2 * pc_z[k] * psi_2[k];

        t_6[k] = f_5 * psh0_1[k]
                 - f_6 * psh1_1[k]
                 + f_2 * pc_y[k] * psi_3[k];

        t_7[k] = f_2 * pc_z[k] * psi_3[k];

        t_8[k] = f_2 * pc_y[k] * psi_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pc_y, pc_z, psh0_2, psh0_3, psh0_5, psh1_2, \
                         psh1_3, psh1_5, psi_5, psi_6, psi_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * psh0_2[k]
                 - f_6 * psh1_2[k]
                 + f_2 * pc_z[k] * psi_5[k];

        t_10[k] = f_7 * psh0_3[k]
                  - f_8 * psh1_3[k]
                  + f_2 * pc_y[k] * psi_6[k];

        t_11[k] = f_2 * pc_z[k] * psi_6[k];

        t_12[k] = f_3 * psh0_5[k]
                  - f_4 * psh1_5[k]
                  + f_2 * pc_y[k] * psi_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pc_y, pc_z, psh0_5, psh0_6, psh0_8, \
                         psh1_5, psh1_6, psh1_8, psi_9, psi_10, \
                         psi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * pc_y[k] * psi_9[k];

        t_14[k] = f_7 * psh0_5[k]
                  - f_8 * psh1_5[k]
                  + f_2 * pc_z[k] * psi_9[k];

        t_15[k] = f_9 * psh0_6[k]
                  - f_10 * psh1_6[k]
                  + f_2 * pc_y[k] * psi_10[k];

        t_16[k] = f_2 * pc_z[k] * psi_10[k];

        t_17[k] = f_5 * psh0_8[k]
                  - f_6 * psh1_8[k]
                  + f_2 * pc_y[k] * psi_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, ssi_21, psh0_9, \
                         psh1_9, psi_13, psi_14, psi_15, psi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * psh0_9[k]
                  - f_4 * psh1_9[k]
                  + f_2 * pc_y[k] * psi_13[k];

        t_19[k] = f_2 * pc_y[k] * psi_14[k];

        t_20[k] = f_9 * psh0_9[k]
                  - f_10 * psh1_9[k]
                  + f_2 * pc_z[k] * psi_14[k];

        t_21[k] = f_11 * ssi_21[k]
                  + f_2 * pc_x[k] * psi_21[k];

        t_22[k] = f_2 * pc_z[k] * psi_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, pc_y, ssi_23, ssi_24, ssi_25, \
                         ssi_27, psi_20, psi_23, psi_24, psi_25, \
                         psi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_11 * ssi_23[k]
                  + f_2 * pc_x[k] * psi_23[k];

        t_24[k] = f_11 * ssi_24[k]
                  + f_2 * pc_x[k] * psi_24[k];

        t_25[k] = f_11 * ssi_25[k]
                  + f_2 * pc_x[k] * psi_25[k];

        t_26[k] = f_2 * pc_y[k] * psi_20[k];

        t_27[k] = f_11 * ssi_27[k]
                  + f_2 * pc_x[k] * psi_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pc_x, pc_z, ssk0_28, ssk0_30, ssk0_31, \
                         ssk1_28, ssk1_30, ssk1_31, psi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_x[k] * ssk0_28[k]
                  - f_1 * pc_x[k] * ssk1_28[k];

        t_29[k] = f_2 * pc_z[k] * psi_21[k];

        t_30[k] = pa_x[k] * ssk0_30[k]
                  - f_1 * pc_x[k] * ssk1_30[k];

        t_31[k] = pa_x[k] * ssk0_31[k]
                  - f_1 * pc_x[k] * ssk1_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pc_x, pc_y, ssk0_32, ssk0_33, ssk0_35, \
                         ssk1_32, ssk1_33, ssk1_35, psi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pa_x[k] * ssk0_32[k]
                  - f_1 * pc_x[k] * ssk1_32[k];

        t_33[k] = pa_x[k] * ssk0_33[k]
                  - f_1 * pc_x[k] * ssk1_33[k];

        t_34[k] = f_2 * pc_y[k] * psi_27[k];

        t_35[k] = pa_x[k] * ssk0_35[k]
                  - f_1 * pc_x[k] * ssk1_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_y, pc_x, pc_y, pc_z, ssk0_0, ssk1_0, psh0_22, \
                         psh1_22, psi_28, psi_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_y[k] * ssk0_0[k]
                  - f_1 * pc_y[k] * ssk1_0[k];

        t_37[k] = f_12 * psh0_22[k]
                  - f_13 * psh1_22[k]
                  + f_2 * pc_x[k] * psi_29[k];

        t_38[k] = f_2 * pc_z[k] * psi_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pc_x, pc_y, pc_z, ssk0_5, ssk1_5, psh0_24, \
                         psh1_24, psi_29, psi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * psh0_24[k]
                  - f_10 * psh1_24[k]
                  + f_2 * pc_x[k] * psi_31[k];

        t_40[k] = f_2 * pc_z[k] * psi_29[k];

        t_41[k] = pa_y[k] * ssk0_5[k]
                  - f_1 * pc_y[k] * ssk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pc_x, pc_z, psh0_27, psh0_29, psh1_27, psh1_29, \
                         psi_31, psi_34, psi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_7 * psh0_27[k]
                  - f_8 * psh1_27[k]
                  + f_2 * pc_x[k] * psi_34[k];

        t_43[k] = f_2 * pc_z[k] * psi_31[k];

        t_44[k] = f_7 * psh0_29[k]
                  - f_8 * psh1_29[k]
                  + f_2 * pc_x[k] * psi_36[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pc_x, pc_y, pc_z, ssk0_9, ssk1_9, psh0_31, \
                         psh1_31, psi_34, psi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * ssk0_9[k]
                  - f_1 * pc_y[k] * ssk1_9[k];

        t_46[k] = f_5 * psh0_31[k]
                  - f_6 * psh1_31[k]
                  + f_2 * pc_x[k] * psi_38[k];

        t_47[k] = f_2 * pc_z[k] * psi_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pa_y, pc_x, pc_y, ssk0_14, ssk1_14, psh0_33, \
                         psh0_34, psh1_33, psh1_34, psi_40, psi_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * psh0_33[k]
                  - f_6 * psh1_33[k]
                  + f_2 * pc_x[k] * psi_40[k];

        t_49[k] = f_5 * psh0_34[k]
                  - f_6 * psh1_34[k]
                  + f_2 * pc_x[k] * psi_41[k];

        t_50[k] = pa_y[k] * ssk0_14[k]
                  - f_1 * pc_y[k] * ssk1_14[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pc_x, pc_z, psh0_36, psh0_38, psh0_39, \
                         psh1_36, psh1_38, psh1_39, psi_38, psi_43, psi_45, \
                         psi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * psh0_36[k]
                  - f_4 * psh1_36[k]
                  + f_2 * pc_x[k] * psi_43[k];

        t_52[k] = f_2 * pc_z[k] * psi_38[k];

        t_53[k] = f_3 * psh0_38[k]
                  - f_4 * psh1_38[k]
                  + f_2 * pc_x[k] * psi_45[k];

        t_54[k] = f_3 * psh0_39[k]
                  - f_4 * psh1_39[k]
                  + f_2 * pc_x[k] * psi_46[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_y, pc_x, pc_y, ssk0_20, ssk1_20, \
                         psh0_40, psh1_40, psi_47, psi_49, psi_50, \
                         psi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_3 * psh0_40[k]
                  - f_4 * psh1_40[k]
                  + f_2 * pc_x[k] * psi_47[k];

        t_56[k] = pa_y[k] * ssk0_20[k]
                  - f_1 * pc_y[k] * ssk1_20[k];

        t_57[k] = f_2 * pc_x[k] * psi_49[k];

        t_58[k] = f_2 * pc_x[k] * psi_50[k];

        t_59[k] = f_2 * pc_x[k] * psi_51[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_y, pc_x, pc_y, ssk0_28, ssi_21, \
                         ssk1_28, psi_52, psi_53, psi_54, psi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * pc_x[k] * psi_52[k];

        t_61[k] = f_2 * pc_x[k] * psi_53[k];

        t_62[k] = f_2 * pc_x[k] * psi_54[k];

        t_63[k] = f_2 * pc_x[k] * psi_55[k];

        t_64[k] = pa_y[k] * ssk0_28[k]
                  + f_0 * ssi_21[k]
                  - f_1 * pc_y[k] * ssk1_28[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_z, psh0_36, psh0_37, psh0_38, psh1_36, \
                         psh1_37, psh1_38, psi_49, psi_50, psi_51, \
                         psi_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_2 * pc_z[k] * psi_49[k];

        t_66[k] = f_3 * psh0_36[k]
                  - f_4 * psh1_36[k]
                  + f_2 * pc_z[k] * psi_50[k];

        t_67[k] = f_5 * psh0_37[k]
                  - f_6 * psh1_37[k]
                  + f_2 * pc_z[k] * psi_51[k];

        t_68[k] = f_7 * psh0_38[k]
                  - f_8 * psh1_38[k]
                  + f_2 * pc_z[k] * psi_52[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, ssk0_35, ssi_27, ssk1_35, \
                         psh0_39, psh1_39, psi_53, psi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_9 * psh0_39[k]
                  - f_10 * psh1_39[k]
                  + f_2 * pc_z[k] * psi_53[k];

        t_70[k] = f_11 * ssi_27[k]
                  + f_2 * pc_y[k] * psi_55[k];

        t_71[k] = pa_y[k] * ssk0_35[k]
                  - f_1 * pc_y[k] * ssk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pc_x, pc_y, pc_z, ssk0_0, ssk0_3, \
                         ssk1_0, ssk1_3, psh0_44, psh1_44, psi_56, \
                         psi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * ssk0_0[k]
                  - f_1 * pc_z[k] * ssk1_0[k];

        t_73[k] = f_2 * pc_y[k] * psi_56[k];

        t_74[k] = f_12 * psh0_44[k]
                  - f_13 * psh1_44[k]
                  + f_2 * pc_x[k] * psi_58[k];

        t_75[k] = pa_z[k] * ssk0_3[k]
                  - f_1 * pc_z[k] * ssk1_3[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_z, pc_x, pc_y, pc_z, ssk0_6, ssk1_6, psh0_47, \
                         psh1_47, psi_58, psi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * pc_y[k] * psi_58[k];

        t_77[k] = f_9 * psh0_47[k]
                  - f_10 * psh1_47[k]
                  + f_2 * pc_x[k] * psi_61[k];

        t_78[k] = pa_z[k] * ssk0_6[k]
                  - f_1 * pc_z[k] * ssk1_6[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pc_x, pc_y, psh0_49, psh0_51, psh1_49, psh1_51, \
                         psi_61, psi_63, psi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_7 * psh0_49[k]
                  - f_8 * psh1_49[k]
                  + f_2 * pc_x[k] * psi_63[k];

        t_80[k] = f_2 * pc_y[k] * psi_61[k];

        t_81[k] = f_7 * psh0_51[k]
                  - f_8 * psh1_51[k]
                  + f_2 * pc_x[k] * psi_65[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_z, pc_x, pc_z, ssk0_10, ssk1_10, psh0_53, \
                         psh0_54, psh1_53, psh1_54, psi_67, psi_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_z[k] * ssk0_10[k]
                  - f_1 * pc_z[k] * ssk1_10[k];

        t_83[k] = f_5 * psh0_53[k]
                  - f_6 * psh1_53[k]
                  + f_2 * pc_x[k] * psi_67[k];

        t_84[k] = f_5 * psh0_54[k]
                  - f_6 * psh1_54[k]
                  + f_2 * pc_x[k] * psi_68[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pa_z, pc_x, pc_y, pc_z, ssk0_15, ssk1_15, psh0_56, \
                         psh1_56, psi_65, psi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_2 * pc_y[k] * psi_65[k];

        t_86[k] = f_5 * psh0_56[k]
                  - f_6 * psh1_56[k]
                  + f_2 * pc_x[k] * psi_70[k];

        t_87[k] = pa_z[k] * ssk0_15[k]
                  - f_1 * pc_z[k] * ssk1_15[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_x, pc_y, psh0_58, psh0_59, psh0_60, \
                         psh1_58, psh1_59, psh1_60, psi_70, psi_72, psi_73, \
                         psi_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * psh0_58[k]
                  - f_4 * psh1_58[k]
                  + f_2 * pc_x[k] * psi_72[k];

        t_89[k] = f_3 * psh0_59[k]
                  - f_4 * psh1_59[k]
                  + f_2 * pc_x[k] * psi_73[k];

        t_90[k] = f_3 * psh0_60[k]
                  - f_4 * psh1_60[k]
                  + f_2 * pc_x[k] * psi_74[k];

        t_91[k] = f_2 * pc_y[k] * psi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, t_97, pc_x, psh0_62, psh1_62, psi_76, \
                         psi_77, psi_78, psi_79, psi_80, psi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_3 * psh0_62[k]
                  - f_4 * psh1_62[k]
                  + f_2 * pc_x[k] * psi_76[k];

        t_93[k] = f_2 * pc_x[k] * psi_77[k];

        t_94[k] = f_2 * pc_x[k] * psi_78[k];

        t_95[k] = f_2 * pc_x[k] * psi_79[k];

        t_96[k] = f_2 * pc_x[k] * psi_80[k];

        t_97[k] = f_2 * pc_x[k] * psi_81[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pa_z, pc_x, pc_y, pc_z, ssk0_28, ssk1_28, \
                         psh0_58, psh1_58, psi_78, psi_82, psi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_2 * pc_x[k] * psi_82[k];

        t_99[k] = f_2 * pc_x[k] * psi_83[k];

        t_100[k] = pa_z[k] * ssk0_28[k]
                   - f_1 * pc_z[k] * ssk1_28[k];

        t_101[k] = f_12 * psh0_58[k]
                   - f_13 * psh1_58[k]
                   + f_2 * pc_y[k] * psi_78[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, psh0_59, psh0_60, psh0_61, psh1_59, \
                         psh1_60, psh1_61, psi_79, psi_80, psi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * psh0_59[k]
                   - f_10 * psh1_59[k]
                   + f_2 * pc_y[k] * psi_79[k];

        t_103[k] = f_7 * psh0_60[k]
                   - f_8 * psh1_60[k]
                   + f_2 * pc_y[k] * psi_80[k];

        t_104[k] = f_5 * psh0_61[k]
                   - f_6 * psh1_61[k]
                   + f_2 * pc_y[k] * psi_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_z, pc_y, pc_z, ssk0_35, ssi_27, ssk1_35, \
                         psh0_62, psh1_62, psi_82, psi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_3 * psh0_62[k]
                   - f_4 * psh1_62[k]
                   + f_2 * pc_y[k] * psi_82[k];

        t_106[k] = f_2 * pc_y[k] * psi_83[k];

        t_107[k] = pa_z[k] * ssk0_35[k]
                   + f_0 * ssi_27[k]
                   - f_1 * pc_z[k] * ssk1_35[k];
    }
}

}  // namespace simdt3ceri
