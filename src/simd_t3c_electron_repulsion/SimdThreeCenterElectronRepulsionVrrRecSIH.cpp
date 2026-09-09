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


#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sih_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shh0,
                                                          const size_t shg, const size_t shh1,
                                                          const size_t sif0, const size_t sif1,
                                                          const size_t sig, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_13 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shh0_0 = buffer.data(shh0 + 0);
    const auto *shh0_3 = buffer.data(shh0 + 3);
    const auto *shh0_5 = buffer.data(shh0 + 5);
    const auto *shh0_6 = buffer.data(shh0 + 6);
    const auto *shh0_9 = buffer.data(shh0 + 9);
    const auto *shh0_15 = buffer.data(shh0 + 15);
    const auto *shh0_20 = buffer.data(shh0 + 20);
    const auto *shh0_24 = buffer.data(shh0 + 24);
    const auto *shh0_27 = buffer.data(shh0 + 27);
    const auto *shh0_36 = buffer.data(shh0 + 36);
    const auto *shh0_42 = buffer.data(shh0 + 42);
    const auto *shh0_47 = buffer.data(shh0 + 47);
    const auto *shh0_51 = buffer.data(shh0 + 51);
    const auto *shh0_62 = buffer.data(shh0 + 62);

    const auto *shg_0 = buffer.data(shg + 0);
    const auto *shg_1 = buffer.data(shg + 1);
    const auto *shg_2 = buffer.data(shg + 2);
    const auto *shg_3 = buffer.data(shg + 3);
    const auto *shg_5 = buffer.data(shg + 5);
    const auto *shg_6 = buffer.data(shg + 6);
    const auto *shg_9 = buffer.data(shg + 9);
    const auto *shg_10 = buffer.data(shg + 10);
    const auto *shg_11 = buffer.data(shg + 11);
    const auto *shg_12 = buffer.data(shg + 12);
    const auto *shg_13 = buffer.data(shg + 13);
    const auto *shg_14 = buffer.data(shg + 14);
    const auto *shg_15 = buffer.data(shg + 15);
    const auto *shg_17 = buffer.data(shg + 17);
    const auto *shg_18 = buffer.data(shg + 18);
    const auto *shg_20 = buffer.data(shg + 20);
    const auto *shg_25 = buffer.data(shg + 25);
    const auto *shg_26 = buffer.data(shg + 26);
    const auto *shg_27 = buffer.data(shg + 27);
    const auto *shg_28 = buffer.data(shg + 28);
    const auto *shg_29 = buffer.data(shg + 29);
    const auto *shg_30 = buffer.data(shg + 30);
    const auto *shg_32 = buffer.data(shg + 32);
    const auto *shg_33 = buffer.data(shg + 33);
    const auto *shg_35 = buffer.data(shg + 35);
    const auto *shg_40 = buffer.data(shg + 40);
    const auto *shg_41 = buffer.data(shg + 41);
    const auto *shg_42 = buffer.data(shg + 42);
    const auto *shg_43 = buffer.data(shg + 43);
    const auto *shg_44 = buffer.data(shg + 44);
    const auto *shg_45 = buffer.data(shg + 45);
    const auto *shg_48 = buffer.data(shg + 48);
    const auto *shg_50 = buffer.data(shg + 50);
    const auto *shg_51 = buffer.data(shg + 51);
    const auto *shg_54 = buffer.data(shg + 54);
    const auto *shg_55 = buffer.data(shg + 55);
    const auto *shg_56 = buffer.data(shg + 56);
    const auto *shg_57 = buffer.data(shg + 57);
    const auto *shg_58 = buffer.data(shg + 58);
    const auto *shg_59 = buffer.data(shg + 59);
    const auto *shg_70 = buffer.data(shg + 70);
    const auto *shg_71 = buffer.data(shg + 71);
    const auto *shg_72 = buffer.data(shg + 72);
    const auto *shg_73 = buffer.data(shg + 73);
    const auto *shg_74 = buffer.data(shg + 74);
    const auto *shg_75 = buffer.data(shg + 75);
    const auto *shg_78 = buffer.data(shg + 78);
    const auto *shg_80 = buffer.data(shg + 80);
    const auto *shg_81 = buffer.data(shg + 81);
    const auto *shg_84 = buffer.data(shg + 84);
    const auto *shg_85 = buffer.data(shg + 85);
    const auto *shg_86 = buffer.data(shg + 86);
    const auto *shg_87 = buffer.data(shg + 87);
    const auto *shg_88 = buffer.data(shg + 88);
    const auto *shg_89 = buffer.data(shg + 89);

    const auto *shh1_0 = buffer.data(shh1 + 0);
    const auto *shh1_3 = buffer.data(shh1 + 3);
    const auto *shh1_5 = buffer.data(shh1 + 5);
    const auto *shh1_6 = buffer.data(shh1 + 6);
    const auto *shh1_9 = buffer.data(shh1 + 9);
    const auto *shh1_15 = buffer.data(shh1 + 15);
    const auto *shh1_20 = buffer.data(shh1 + 20);
    const auto *shh1_24 = buffer.data(shh1 + 24);
    const auto *shh1_27 = buffer.data(shh1 + 27);
    const auto *shh1_36 = buffer.data(shh1 + 36);
    const auto *shh1_42 = buffer.data(shh1 + 42);
    const auto *shh1_47 = buffer.data(shh1 + 47);
    const auto *shh1_51 = buffer.data(shh1 + 51);
    const auto *shh1_62 = buffer.data(shh1 + 62);

    const auto *sif0_0 = buffer.data(sif0 + 0);
    const auto *sif0_3 = buffer.data(sif0 + 3);
    const auto *sif0_5 = buffer.data(sif0 + 5);
    const auto *sif0_6 = buffer.data(sif0 + 6);
    const auto *sif0_8 = buffer.data(sif0 + 8);
    const auto *sif0_9 = buffer.data(sif0 + 9);
    const auto *sif0_16 = buffer.data(sif0 + 16);
    const auto *sif0_18 = buffer.data(sif0 + 18);
    const auto *sif0_19 = buffer.data(sif0 + 19);
    const auto *sif0_28 = buffer.data(sif0 + 28);
    const auto *sif0_29 = buffer.data(sif0 + 29);
    const auto *sif0_30 = buffer.data(sif0 + 30);
    const auto *sif0_33 = buffer.data(sif0 + 33);
    const auto *sif0_35 = buffer.data(sif0 + 35);
    const auto *sif0_36 = buffer.data(sif0 + 36);
    const auto *sif0_38 = buffer.data(sif0 + 38);
    const auto *sif0_39 = buffer.data(sif0 + 39);
    const auto *sif0_48 = buffer.data(sif0 + 48);
    const auto *sif0_49 = buffer.data(sif0 + 49);
    const auto *sif0_50 = buffer.data(sif0 + 50);
    const auto *sif0_53 = buffer.data(sif0 + 53);
    const auto *sif0_55 = buffer.data(sif0 + 55);
    const auto *sif0_56 = buffer.data(sif0 + 56);
    const auto *sif0_58 = buffer.data(sif0 + 58);
    const auto *sif0_59 = buffer.data(sif0 + 59);

    const auto *sif1_0 = buffer.data(sif1 + 0);
    const auto *sif1_3 = buffer.data(sif1 + 3);
    const auto *sif1_5 = buffer.data(sif1 + 5);
    const auto *sif1_6 = buffer.data(sif1 + 6);
    const auto *sif1_8 = buffer.data(sif1 + 8);
    const auto *sif1_9 = buffer.data(sif1 + 9);
    const auto *sif1_16 = buffer.data(sif1 + 16);
    const auto *sif1_18 = buffer.data(sif1 + 18);
    const auto *sif1_19 = buffer.data(sif1 + 19);
    const auto *sif1_28 = buffer.data(sif1 + 28);
    const auto *sif1_29 = buffer.data(sif1 + 29);
    const auto *sif1_30 = buffer.data(sif1 + 30);
    const auto *sif1_33 = buffer.data(sif1 + 33);
    const auto *sif1_35 = buffer.data(sif1 + 35);
    const auto *sif1_36 = buffer.data(sif1 + 36);
    const auto *sif1_38 = buffer.data(sif1 + 38);
    const auto *sif1_39 = buffer.data(sif1 + 39);
    const auto *sif1_48 = buffer.data(sif1 + 48);
    const auto *sif1_49 = buffer.data(sif1 + 49);
    const auto *sif1_50 = buffer.data(sif1 + 50);
    const auto *sif1_53 = buffer.data(sif1 + 53);
    const auto *sif1_55 = buffer.data(sif1 + 55);
    const auto *sif1_56 = buffer.data(sif1 + 56);
    const auto *sif1_58 = buffer.data(sif1 + 58);
    const auto *sif1_59 = buffer.data(sif1 + 59);

    const auto *sig_0 = buffer.data(sig + 0);
    const auto *sig_2 = buffer.data(sig + 2);
    const auto *sig_3 = buffer.data(sig + 3);
    const auto *sig_5 = buffer.data(sig + 5);
    const auto *sig_6 = buffer.data(sig + 6);
    const auto *sig_9 = buffer.data(sig + 9);
    const auto *sig_10 = buffer.data(sig + 10);
    const auto *sig_11 = buffer.data(sig + 11);
    const auto *sig_12 = buffer.data(sig + 12);
    const auto *sig_13 = buffer.data(sig + 13);
    const auto *sig_14 = buffer.data(sig + 14);
    const auto *sig_15 = buffer.data(sig + 15);
    const auto *sig_17 = buffer.data(sig + 17);
    const auto *sig_18 = buffer.data(sig + 18);
    const auto *sig_20 = buffer.data(sig + 20);
    const auto *sig_25 = buffer.data(sig + 25);
    const auto *sig_26 = buffer.data(sig + 26);
    const auto *sig_27 = buffer.data(sig + 27);
    const auto *sig_28 = buffer.data(sig + 28);
    const auto *sig_29 = buffer.data(sig + 29);
    const auto *sig_30 = buffer.data(sig + 30);
    const auto *sig_32 = buffer.data(sig + 32);
    const auto *sig_33 = buffer.data(sig + 33);
    const auto *sig_35 = buffer.data(sig + 35);
    const auto *sig_40 = buffer.data(sig + 40);
    const auto *sig_41 = buffer.data(sig + 41);
    const auto *sig_42 = buffer.data(sig + 42);
    const auto *sig_43 = buffer.data(sig + 43);
    const auto *sig_44 = buffer.data(sig + 44);
    const auto *sig_45 = buffer.data(sig + 45);
    const auto *sig_47 = buffer.data(sig + 47);
    const auto *sig_48 = buffer.data(sig + 48);
    const auto *sig_50 = buffer.data(sig + 50);
    const auto *sig_51 = buffer.data(sig + 51);
    const auto *sig_54 = buffer.data(sig + 54);
    const auto *sig_55 = buffer.data(sig + 55);
    const auto *sig_56 = buffer.data(sig + 56);
    const auto *sig_57 = buffer.data(sig + 57);
    const auto *sig_58 = buffer.data(sig + 58);
    const auto *sig_59 = buffer.data(sig + 59);
    const auto *sig_60 = buffer.data(sig + 60);
    const auto *sig_62 = buffer.data(sig + 62);
    const auto *sig_63 = buffer.data(sig + 63);
    const auto *sig_65 = buffer.data(sig + 65);
    const auto *sig_70 = buffer.data(sig + 70);
    const auto *sig_71 = buffer.data(sig + 71);
    const auto *sig_72 = buffer.data(sig + 72);
    const auto *sig_73 = buffer.data(sig + 73);
    const auto *sig_74 = buffer.data(sig + 74);
    const auto *sig_75 = buffer.data(sig + 75);
    const auto *sig_77 = buffer.data(sig + 77);
    const auto *sig_78 = buffer.data(sig + 78);
    const auto *sig_80 = buffer.data(sig + 80);
    const auto *sig_81 = buffer.data(sig + 81);
    const auto *sig_84 = buffer.data(sig + 84);
    const auto *sig_85 = buffer.data(sig + 85);
    const auto *sig_86 = buffer.data(sig + 86);
    const auto *sig_87 = buffer.data(sig + 87);
    const auto *sig_88 = buffer.data(sig + 88);
    const auto *sig_89 = buffer.data(sig + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, shg_0, shg_3, sif0_0, sif0_3, \
                         sif1_0, sif1_3, sig_0, sig_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * shg_0[k]
                 + f_1 * sif0_0[k]
                 - f_2 * sif1_0[k]
                 + f_3 * pc_x[k] * sig_0[k];

        t_1[k] = f_3 * pc_y[k] * sig_0[k];

        t_2[k] = f_3 * pc_z[k] * sig_0[k];

        t_3[k] = f_0 * shg_3[k]
                 + f_4 * sif0_3[k]
                 - f_5 * sif1_3[k]
                 + f_3 * pc_x[k] * sig_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, shg_5, shg_6, sif0_5, sif0_6, sif1_5, \
                         sif1_6, sig_2, sig_5, sig_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sig_2[k];

        t_5[k] = f_0 * shg_5[k]
                 + f_4 * sif0_5[k]
                 - f_5 * sif1_5[k]
                 + f_3 * pc_x[k] * sig_5[k];

        t_6[k] = f_0 * shg_6[k]
                 + f_6 * sif0_6[k]
                 - f_7 * sif1_6[k]
                 + f_3 * pc_x[k] * sig_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, shg_9, shg_10, sif0_9, sif1_9, \
                         sig_3, sig_5, sig_9, sig_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sig_3[k];

        t_8[k] = f_3 * pc_y[k] * sig_5[k];

        t_9[k] = f_0 * shg_9[k]
                 + f_6 * sif0_9[k]
                 - f_7 * sif1_9[k]
                 + f_3 * pc_x[k] * sig_9[k];

        t_10[k] = f_0 * shg_10[k]
                  + f_3 * pc_x[k] * sig_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, shg_11, shg_12, shg_13, shg_14, sig_11, \
                         sig_12, sig_13, sig_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * shg_11[k]
                  + f_3 * pc_x[k] * sig_11[k];

        t_12[k] = f_0 * shg_12[k]
                  + f_3 * pc_x[k] * sig_12[k];

        t_13[k] = f_0 * shg_13[k]
                  + f_3 * pc_x[k] * sig_13[k];

        t_14[k] = f_0 * shg_14[k]
                  + f_3 * pc_x[k] * sig_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, sif0_6, sif0_8, sif0_9, sif1_6, \
                         sif1_8, sif1_9, sig_10, sig_12, sig_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * sif0_6[k]
                  - f_2 * sif1_6[k]
                  + f_3 * pc_y[k] * sig_10[k];

        t_16[k] = f_3 * pc_z[k] * sig_10[k];

        t_17[k] = f_4 * sif0_8[k]
                  - f_5 * sif1_8[k]
                  + f_3 * pc_y[k] * sig_12[k];

        t_18[k] = f_6 * sif0_9[k]
                  - f_7 * sif1_9[k]
                  + f_3 * pc_y[k] * sig_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, shh0_0, shg_0, \
                         shh1_0, sif0_9, sif1_9, sig_14, sig_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sig_14[k];

        t_20[k] = f_1 * sif0_9[k]
                  - f_2 * sif1_9[k]
                  + f_3 * pc_z[k] * sig_14[k];

        t_21[k] = pb_y[k] * shh0_0[k]
                  - f_8 * pc_y[k] * shh1_0[k];

        t_22[k] = f_9 * shg_0[k]
                  + f_3 * pc_y[k] * sig_15[k];

        t_23[k] = f_3 * pc_z[k] * sig_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_y, shh0_3, shh0_5, shh0_6, shg_1, \
                         shg_2, shg_3, shh1_3, shh1_5, shh1_6, sig_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * shh0_3[k]
                  + f_10 * shg_1[k]
                  - f_8 * pc_y[k] * shh1_3[k];

        t_25[k] = f_9 * shg_2[k]
                  + f_3 * pc_y[k] * sig_17[k];

        t_26[k] = pb_y[k] * shh0_5[k]
                  - f_8 * pc_y[k] * shh1_5[k];

        t_27[k] = pb_y[k] * shh0_6[k]
                  + f_11 * shg_3[k]
                  - f_8 * pc_y[k] * shh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_x, pc_y, pc_z, shh0_9, shg_5, \
                         shg_25, shh1_9, sig_18, sig_20, sig_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * sig_18[k];

        t_29[k] = f_9 * shg_5[k]
                  + f_3 * pc_y[k] * sig_20[k];

        t_30[k] = pb_y[k] * shh0_9[k]
                  - f_8 * pc_y[k] * shh1_9[k];

        t_31[k] = f_12 * shg_25[k]
                  + f_3 * pc_x[k] * sig_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, shg_26, shg_27, shg_28, shg_29, sig_26, \
                         sig_27, sig_28, sig_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * shg_26[k]
                  + f_3 * pc_x[k] * sig_26[k];

        t_33[k] = f_12 * shg_27[k]
                  + f_3 * pc_x[k] * sig_27[k];

        t_34[k] = f_12 * shg_28[k]
                  + f_3 * pc_x[k] * sig_28[k];

        t_35[k] = f_12 * shg_29[k]
                  + f_3 * pc_x[k] * sig_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pc_y, pc_z, shg_10, shg_12, sif0_16, sif0_18, \
                         sif1_16, sif1_18, sig_25, sig_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * shg_10[k]
                  + f_1 * sif0_16[k]
                  - f_2 * sif1_16[k]
                  + f_3 * pc_y[k] * sig_25[k];

        t_37[k] = f_3 * pc_z[k] * sig_25[k];

        t_38[k] = f_9 * shg_12[k]
                  + f_4 * sif0_18[k]
                  - f_5 * sif1_18[k]
                  + f_3 * pc_y[k] * sig_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, shh0_20, shg_13, shg_14, shh1_20, \
                         sif0_19, sif1_19, sig_28, sig_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * shg_13[k]
                  + f_6 * sif0_19[k]
                  - f_7 * sif1_19[k]
                  + f_3 * pc_y[k] * sig_28[k];

        t_40[k] = f_9 * shg_14[k]
                  + f_3 * pc_y[k] * sig_29[k];

        t_41[k] = pb_y[k] * shh0_20[k]
                  - f_8 * pc_y[k] * shh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, shh0_0, shh0_3, \
                         shg_0, shh1_0, shh1_3, sig_30, sig_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * shh0_0[k]
                  - f_8 * pc_z[k] * shh1_0[k];

        t_43[k] = f_3 * pc_y[k] * sig_30[k];

        t_44[k] = f_9 * shg_0[k]
                  + f_3 * pc_z[k] * sig_30[k];

        t_45[k] = pb_z[k] * shh0_3[k]
                  - f_8 * pc_z[k] * shh1_3[k];

        t_46[k] = f_3 * pc_y[k] * sig_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, shh0_5, shh0_6, shg_2, \
                         shg_3, shh1_5, shh1_6, sig_33, sig_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * shh0_5[k]
                  + f_10 * shg_2[k]
                  - f_8 * pc_z[k] * shh1_5[k];

        t_48[k] = pb_z[k] * shh0_6[k]
                  - f_8 * pc_z[k] * shh1_6[k];

        t_49[k] = f_9 * shg_3[k]
                  + f_3 * pc_z[k] * sig_33[k];

        t_50[k] = f_3 * pc_y[k] * sig_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, shh0_9, shg_5, shg_40, \
                         shg_41, shg_42, shh1_9, sig_40, sig_41, \
                         sig_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * shh0_9[k]
                  + f_11 * shg_5[k]
                  - f_8 * pc_z[k] * shh1_9[k];

        t_52[k] = f_12 * shg_40[k]
                  + f_3 * pc_x[k] * sig_40[k];

        t_53[k] = f_12 * shg_41[k]
                  + f_3 * pc_x[k] * sig_41[k];

        t_54[k] = f_12 * shg_42[k]
                  + f_3 * pc_x[k] * sig_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_z, shh0_15, shg_10, shg_43, \
                         shg_44, shh1_15, sig_40, sig_43, sig_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_12 * shg_43[k]
                  + f_3 * pc_x[k] * sig_43[k];

        t_56[k] = f_12 * shg_44[k]
                  + f_3 * pc_x[k] * sig_44[k];

        t_57[k] = pb_z[k] * shh0_15[k]
                  - f_8 * pc_z[k] * shh1_15[k];

        t_58[k] = f_9 * shg_10[k]
                  + f_3 * pc_z[k] * sig_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_y, pc_z, shg_14, sif0_28, sif0_29, \
                         sif1_28, sif1_29, sig_42, sig_43, sig_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_4 * sif0_28[k]
                  - f_5 * sif1_28[k]
                  + f_3 * pc_y[k] * sig_42[k];

        t_60[k] = f_6 * sif0_29[k]
                  - f_7 * sif1_29[k]
                  + f_3 * pc_y[k] * sig_43[k];

        t_61[k] = f_3 * pc_y[k] * sig_44[k];

        t_62[k] = f_9 * shg_14[k]
                  + f_1 * sif0_29[k]
                  - f_2 * sif1_29[k]
                  + f_3 * pc_z[k] * sig_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, shg_15, shg_45, shg_48, \
                         sif0_30, sif0_33, sif1_30, sif1_33, sig_45, \
                         sig_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_13 * shg_45[k]
                  + f_1 * sif0_30[k]
                  - f_2 * sif1_30[k]
                  + f_3 * pc_x[k] * sig_45[k];

        t_64[k] = f_10 * shg_15[k]
                  + f_3 * pc_y[k] * sig_45[k];

        t_65[k] = f_3 * pc_z[k] * sig_45[k];

        t_66[k] = f_13 * shg_48[k]
                  + f_4 * sif0_33[k]
                  - f_5 * sif1_33[k]
                  + f_3 * pc_x[k] * sig_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_x, pc_y, shg_17, shg_50, shg_51, sif0_35, \
                         sif0_36, sif1_35, sif1_36, sig_47, sig_50, \
                         sig_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * shg_17[k]
                  + f_3 * pc_y[k] * sig_47[k];

        t_68[k] = f_13 * shg_50[k]
                  + f_4 * sif0_35[k]
                  - f_5 * sif1_35[k]
                  + f_3 * pc_x[k] * sig_50[k];

        t_69[k] = f_13 * shg_51[k]
                  + f_6 * sif0_36[k]
                  - f_7 * sif1_36[k]
                  + f_3 * pc_x[k] * sig_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pc_x, pc_y, pc_z, shg_20, shg_54, shg_55, \
                         sif0_39, sif1_39, sig_48, sig_50, sig_54, \
                         sig_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * sig_48[k];

        t_71[k] = f_10 * shg_20[k]
                  + f_3 * pc_y[k] * sig_50[k];

        t_72[k] = f_13 * shg_54[k]
                  + f_6 * sif0_39[k]
                  - f_7 * sif1_39[k]
                  + f_3 * pc_x[k] * sig_54[k];

        t_73[k] = f_13 * shg_55[k]
                  + f_3 * pc_x[k] * sig_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, shg_56, shg_57, shg_58, shg_59, sig_56, \
                         sig_57, sig_58, sig_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * shg_56[k]
                  + f_3 * pc_x[k] * sig_56[k];

        t_75[k] = f_13 * shg_57[k]
                  + f_3 * pc_x[k] * sig_57[k];

        t_76[k] = f_13 * shg_58[k]
                  + f_3 * pc_x[k] * sig_58[k];

        t_77[k] = f_13 * shg_59[k]
                  + f_3 * pc_x[k] * sig_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, shg_25, shg_27, sif0_36, sif0_38, \
                         sif1_36, sif1_38, sig_55, sig_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * shg_25[k]
                  + f_1 * sif0_36[k]
                  - f_2 * sif1_36[k]
                  + f_3 * pc_y[k] * sig_55[k];

        t_79[k] = f_3 * pc_z[k] * sig_55[k];

        t_80[k] = f_10 * shg_27[k]
                  + f_4 * sif0_38[k]
                  - f_5 * sif1_38[k]
                  + f_3 * pc_y[k] * sig_57[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pc_y, pc_z, shh0_42, shg_28, shg_29, \
                         shh1_42, sif0_39, sif1_39, sig_58, sig_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * shg_28[k]
                  + f_6 * sif0_39[k]
                  - f_7 * sif1_39[k]
                  + f_3 * pc_y[k] * sig_58[k];

        t_82[k] = f_10 * shg_29[k]
                  + f_3 * pc_y[k] * sig_59[k];

        t_83[k] = f_1 * sif0_39[k]
                  - f_2 * sif1_39[k]
                  + f_3 * pc_z[k] * sig_59[k];

        t_84[k] = pb_y[k] * shh0_42[k]
                  - f_8 * pc_y[k] * shh1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, shh0_24, shg_15, shg_30, \
                         shg_32, shh1_24, sig_60, sig_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_9 * shg_30[k]
                  + f_3 * pc_y[k] * sig_60[k];

        t_86[k] = f_9 * shg_15[k]
                  + f_3 * pc_z[k] * sig_60[k];

        t_87[k] = pb_z[k] * shh0_24[k]
                  - f_8 * pc_z[k] * shh1_24[k];

        t_88[k] = f_9 * shg_32[k]
                  + f_3 * pc_y[k] * sig_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, shh0_27, shh0_47, \
                         shg_18, shg_35, shh1_27, shh1_47, sig_63, \
                         sig_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * shh0_47[k]
                  - f_8 * pc_y[k] * shh1_47[k];

        t_90[k] = pb_z[k] * shh0_27[k]
                  - f_8 * pc_z[k] * shh1_27[k];

        t_91[k] = f_9 * shg_18[k]
                  + f_3 * pc_z[k] * sig_63[k];

        t_92[k] = f_9 * shg_35[k]
                  + f_3 * pc_y[k] * sig_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_x, pc_y, shh0_51, shg_70, shg_71, \
                         shg_72, shh1_51, sig_70, sig_71, sig_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * shh0_51[k]
                  - f_8 * pc_y[k] * shh1_51[k];

        t_94[k] = f_13 * shg_70[k]
                  + f_3 * pc_x[k] * sig_70[k];

        t_95[k] = f_13 * shg_71[k]
                  + f_3 * pc_x[k] * sig_71[k];

        t_96[k] = f_13 * shg_72[k]
                  + f_3 * pc_x[k] * sig_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_z, pc_x, pc_z, shh0_36, shg_25, shg_73, \
                         shg_74, shh1_36, sig_70, sig_73, sig_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_13 * shg_73[k]
                  + f_3 * pc_x[k] * sig_73[k];

        t_98[k] = f_13 * shg_74[k]
                  + f_3 * pc_x[k] * sig_74[k];

        t_99[k] = pb_z[k] * shh0_36[k]
                  - f_8 * pc_z[k] * shh1_36[k];

        t_100[k] = f_9 * shg_25[k]
                   + f_3 * pc_z[k] * sig_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, shg_42, shg_43, shg_44, sif0_48, sif0_49, \
                         sif1_48, sif1_49, sig_72, sig_73, sig_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * shg_42[k]
                   + f_4 * sif0_48[k]
                   - f_5 * sif1_48[k]
                   + f_3 * pc_y[k] * sig_72[k];

        t_102[k] = f_9 * shg_43[k]
                   + f_6 * sif0_49[k]
                   - f_7 * sif1_49[k]
                   + f_3 * pc_y[k] * sig_73[k];

        t_103[k] = f_9 * shg_44[k]
                   + f_3 * pc_y[k] * sig_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, pc_x, pc_y, pc_z, shh0_62, shg_30, \
                         shg_75, shh1_62, sif0_50, sif1_50, sig_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_y[k] * shh0_62[k]
                   - f_8 * pc_y[k] * shh1_62[k];

        t_105[k] = f_13 * shg_75[k]
                   + f_1 * sif0_50[k]
                   - f_2 * sif1_50[k]
                   + f_3 * pc_x[k] * sig_75[k];

        t_106[k] = f_3 * pc_y[k] * sig_75[k];

        t_107[k] = f_10 * shg_30[k]
                   + f_3 * pc_z[k] * sig_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pc_y, shg_78, shg_80, sif0_53, sif0_55, \
                         sif1_53, sif1_55, sig_77, sig_78, sig_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_13 * shg_78[k]
                   + f_4 * sif0_53[k]
                   - f_5 * sif1_53[k]
                   + f_3 * pc_x[k] * sig_78[k];

        t_109[k] = f_3 * pc_y[k] * sig_77[k];

        t_110[k] = f_13 * shg_80[k]
                   + f_4 * sif0_55[k]
                   - f_5 * sif1_55[k]
                   + f_3 * pc_x[k] * sig_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_y, pc_z, shg_33, shg_81, sif0_56, \
                         sif1_56, sig_78, sig_80, sig_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * shg_81[k]
                   + f_6 * sif0_56[k]
                   - f_7 * sif1_56[k]
                   + f_3 * pc_x[k] * sig_81[k];

        t_112[k] = f_10 * shg_33[k]
                   + f_3 * pc_z[k] * sig_78[k];

        t_113[k] = f_3 * pc_y[k] * sig_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, shg_84, shg_85, shg_86, shg_87, \
                         sif0_59, sif1_59, sig_84, sig_85, sig_86, \
                         sig_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * shg_84[k]
                   + f_6 * sif0_59[k]
                   - f_7 * sif1_59[k]
                   + f_3 * pc_x[k] * sig_84[k];

        t_115[k] = f_13 * shg_85[k]
                   + f_3 * pc_x[k] * sig_85[k];

        t_116[k] = f_13 * shg_86[k]
                   + f_3 * pc_x[k] * sig_86[k];

        t_117[k] = f_13 * shg_87[k]
                   + f_3 * pc_x[k] * sig_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, pc_z, shg_40, shg_88, shg_89, \
                         sif0_56, sif1_56, sig_85, sig_88, sig_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_13 * shg_88[k]
                   + f_3 * pc_x[k] * sig_88[k];

        t_119[k] = f_13 * shg_89[k]
                   + f_3 * pc_x[k] * sig_89[k];

        t_120[k] = f_1 * sif0_56[k]
                   - f_2 * sif1_56[k]
                   + f_3 * pc_y[k] * sig_85[k];

        t_121[k] = f_10 * shg_40[k]
                   + f_3 * pc_z[k] * sig_85[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, shg_44, sif0_58, sif0_59, \
                         sif1_58, sif1_59, sig_87, sig_88, sig_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_4 * sif0_58[k]
                   - f_5 * sif1_58[k]
                   + f_3 * pc_y[k] * sig_87[k];

        t_123[k] = f_6 * sif0_59[k]
                   - f_7 * sif1_59[k]
                   + f_3 * pc_y[k] * sig_88[k];

        t_124[k] = f_3 * pc_y[k] * sig_89[k];

        t_125[k] = f_10 * shg_44[k]
                   + f_1 * sif0_59[k]
                   - f_2 * sif1_59[k]
                   + f_3 * pc_z[k] * sig_89[k];
    }
}

static auto
compute_prim_sih_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shh0,
                                                          const size_t shg, const size_t shh1,
                                                          const size_t sif0, const size_t sif1,
                                                          const size_t sig, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 2.0 / q;

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
    auto *t_242 = buffer.data(target + 242);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shh0_63 = buffer.data(shh0 + 63);
    const auto *shh0_66 = buffer.data(shh0 + 66);
    const auto *shh0_69 = buffer.data(shh0 + 69);
    const auto *shh0_78 = buffer.data(shh0 + 78);
    const auto *shh0_105 = buffer.data(shh0 + 105);
    const auto *shh0_108 = buffer.data(shh0 + 108);
    const auto *shh0_110 = buffer.data(shh0 + 110);
    const auto *shh0_111 = buffer.data(shh0 + 111);
    const auto *shh0_114 = buffer.data(shh0 + 114);
    const auto *shh0_125 = buffer.data(shh0 + 125);
    const auto *shh0_126 = buffer.data(shh0 + 126);
    const auto *shh0_129 = buffer.data(shh0 + 129);
    const auto *shh0_132 = buffer.data(shh0 + 132);

    const auto *shg_45 = buffer.data(shg + 45);
    const auto *shg_47 = buffer.data(shg + 47);
    const auto *shg_48 = buffer.data(shg + 48);
    const auto *shg_50 = buffer.data(shg + 50);
    const auto *shg_55 = buffer.data(shg + 55);
    const auto *shg_57 = buffer.data(shg + 57);
    const auto *shg_58 = buffer.data(shg + 58);
    const auto *shg_59 = buffer.data(shg + 59);
    const auto *shg_60 = buffer.data(shg + 60);
    const auto *shg_62 = buffer.data(shg + 62);
    const auto *shg_63 = buffer.data(shg + 63);
    const auto *shg_65 = buffer.data(shg + 65);
    const auto *shg_70 = buffer.data(shg + 70);
    const auto *shg_72 = buffer.data(shg + 72);
    const auto *shg_73 = buffer.data(shg + 73);
    const auto *shg_74 = buffer.data(shg + 74);
    const auto *shg_75 = buffer.data(shg + 75);
    const auto *shg_76 = buffer.data(shg + 76);
    const auto *shg_77 = buffer.data(shg + 77);
    const auto *shg_78 = buffer.data(shg + 78);
    const auto *shg_80 = buffer.data(shg + 80);
    const auto *shg_85 = buffer.data(shg + 85);
    const auto *shg_87 = buffer.data(shg + 87);
    const auto *shg_88 = buffer.data(shg + 88);
    const auto *shg_89 = buffer.data(shg + 89);
    const auto *shg_90 = buffer.data(shg + 90);
    const auto *shg_92 = buffer.data(shg + 92);
    const auto *shg_93 = buffer.data(shg + 93);
    const auto *shg_95 = buffer.data(shg + 95);
    const auto *shg_96 = buffer.data(shg + 96);
    const auto *shg_99 = buffer.data(shg + 99);
    const auto *shg_100 = buffer.data(shg + 100);
    const auto *shg_101 = buffer.data(shg + 101);
    const auto *shg_102 = buffer.data(shg + 102);
    const auto *shg_103 = buffer.data(shg + 103);
    const auto *shg_104 = buffer.data(shg + 104);
    const auto *shg_105 = buffer.data(shg + 105);
    const auto *shg_107 = buffer.data(shg + 107);
    const auto *shg_110 = buffer.data(shg + 110);
    const auto *shg_114 = buffer.data(shg + 114);
    const auto *shg_115 = buffer.data(shg + 115);
    const auto *shg_116 = buffer.data(shg + 116);
    const auto *shg_117 = buffer.data(shg + 117);
    const auto *shg_118 = buffer.data(shg + 118);
    const auto *shg_119 = buffer.data(shg + 119);
    const auto *shg_130 = buffer.data(shg + 130);
    const auto *shg_131 = buffer.data(shg + 131);
    const auto *shg_132 = buffer.data(shg + 132);
    const auto *shg_133 = buffer.data(shg + 133);
    const auto *shg_134 = buffer.data(shg + 134);
    const auto *shg_135 = buffer.data(shg + 135);
    const auto *shg_138 = buffer.data(shg + 138);
    const auto *shg_140 = buffer.data(shg + 140);
    const auto *shg_141 = buffer.data(shg + 141);
    const auto *shg_144 = buffer.data(shg + 144);
    const auto *shg_145 = buffer.data(shg + 145);
    const auto *shg_146 = buffer.data(shg + 146);
    const auto *shg_147 = buffer.data(shg + 147);
    const auto *shg_148 = buffer.data(shg + 148);
    const auto *shg_149 = buffer.data(shg + 149);
    const auto *shg_150 = buffer.data(shg + 150);
    const auto *shg_153 = buffer.data(shg + 153);
    const auto *shg_155 = buffer.data(shg + 155);
    const auto *shg_156 = buffer.data(shg + 156);
    const auto *shg_159 = buffer.data(shg + 159);
    const auto *shg_160 = buffer.data(shg + 160);
    const auto *shg_161 = buffer.data(shg + 161);
    const auto *shg_162 = buffer.data(shg + 162);
    const auto *shg_163 = buffer.data(shg + 163);
    const auto *shg_164 = buffer.data(shg + 164);
    const auto *shg_170 = buffer.data(shg + 170);
    const auto *shg_174 = buffer.data(shg + 174);
    const auto *shg_175 = buffer.data(shg + 175);
    const auto *shg_176 = buffer.data(shg + 176);

    const auto *shh1_63 = buffer.data(shh1 + 63);
    const auto *shh1_66 = buffer.data(shh1 + 66);
    const auto *shh1_69 = buffer.data(shh1 + 69);
    const auto *shh1_78 = buffer.data(shh1 + 78);
    const auto *shh1_105 = buffer.data(shh1 + 105);
    const auto *shh1_108 = buffer.data(shh1 + 108);
    const auto *shh1_110 = buffer.data(shh1 + 110);
    const auto *shh1_111 = buffer.data(shh1 + 111);
    const auto *shh1_114 = buffer.data(shh1 + 114);
    const auto *shh1_125 = buffer.data(shh1 + 125);
    const auto *shh1_126 = buffer.data(shh1 + 126);
    const auto *shh1_129 = buffer.data(shh1 + 129);
    const auto *shh1_132 = buffer.data(shh1 + 132);

    const auto *sif0_60 = buffer.data(sif0 + 60);
    const auto *sif0_63 = buffer.data(sif0 + 63);
    const auto *sif0_65 = buffer.data(sif0 + 65);
    const auto *sif0_66 = buffer.data(sif0 + 66);
    const auto *sif0_68 = buffer.data(sif0 + 68);
    const auto *sif0_69 = buffer.data(sif0 + 69);
    const auto *sif0_75 = buffer.data(sif0 + 75);
    const auto *sif0_78 = buffer.data(sif0 + 78);
    const auto *sif0_79 = buffer.data(sif0 + 79);
    const auto *sif0_86 = buffer.data(sif0 + 86);
    const auto *sif0_88 = buffer.data(sif0 + 88);
    const auto *sif0_89 = buffer.data(sif0 + 89);
    const auto *sif0_90 = buffer.data(sif0 + 90);
    const auto *sif0_93 = buffer.data(sif0 + 93);
    const auto *sif0_95 = buffer.data(sif0 + 95);
    const auto *sif0_96 = buffer.data(sif0 + 96);
    const auto *sif0_98 = buffer.data(sif0 + 98);
    const auto *sif0_99 = buffer.data(sif0 + 99);
    const auto *sif0_100 = buffer.data(sif0 + 100);
    const auto *sif0_103 = buffer.data(sif0 + 103);
    const auto *sif0_105 = buffer.data(sif0 + 105);
    const auto *sif0_106 = buffer.data(sif0 + 106);
    const auto *sif0_108 = buffer.data(sif0 + 108);
    const auto *sif0_109 = buffer.data(sif0 + 109);
    const auto *sif0_115 = buffer.data(sif0 + 115);
    const auto *sif0_119 = buffer.data(sif0 + 119);

    const auto *sif1_60 = buffer.data(sif1 + 60);
    const auto *sif1_63 = buffer.data(sif1 + 63);
    const auto *sif1_65 = buffer.data(sif1 + 65);
    const auto *sif1_66 = buffer.data(sif1 + 66);
    const auto *sif1_68 = buffer.data(sif1 + 68);
    const auto *sif1_69 = buffer.data(sif1 + 69);
    const auto *sif1_75 = buffer.data(sif1 + 75);
    const auto *sif1_78 = buffer.data(sif1 + 78);
    const auto *sif1_79 = buffer.data(sif1 + 79);
    const auto *sif1_86 = buffer.data(sif1 + 86);
    const auto *sif1_88 = buffer.data(sif1 + 88);
    const auto *sif1_89 = buffer.data(sif1 + 89);
    const auto *sif1_90 = buffer.data(sif1 + 90);
    const auto *sif1_93 = buffer.data(sif1 + 93);
    const auto *sif1_95 = buffer.data(sif1 + 95);
    const auto *sif1_96 = buffer.data(sif1 + 96);
    const auto *sif1_98 = buffer.data(sif1 + 98);
    const auto *sif1_99 = buffer.data(sif1 + 99);
    const auto *sif1_100 = buffer.data(sif1 + 100);
    const auto *sif1_103 = buffer.data(sif1 + 103);
    const auto *sif1_105 = buffer.data(sif1 + 105);
    const auto *sif1_106 = buffer.data(sif1 + 106);
    const auto *sif1_108 = buffer.data(sif1 + 108);
    const auto *sif1_109 = buffer.data(sif1 + 109);
    const auto *sif1_115 = buffer.data(sif1 + 115);
    const auto *sif1_119 = buffer.data(sif1 + 119);

    const auto *sig_90 = buffer.data(sig + 90);
    const auto *sig_92 = buffer.data(sig + 92);
    const auto *sig_93 = buffer.data(sig + 93);
    const auto *sig_95 = buffer.data(sig + 95);
    const auto *sig_96 = buffer.data(sig + 96);
    const auto *sig_99 = buffer.data(sig + 99);
    const auto *sig_100 = buffer.data(sig + 100);
    const auto *sig_101 = buffer.data(sig + 101);
    const auto *sig_102 = buffer.data(sig + 102);
    const auto *sig_103 = buffer.data(sig + 103);
    const auto *sig_104 = buffer.data(sig + 104);
    const auto *sig_105 = buffer.data(sig + 105);
    const auto *sig_107 = buffer.data(sig + 107);
    const auto *sig_108 = buffer.data(sig + 108);
    const auto *sig_110 = buffer.data(sig + 110);
    const auto *sig_114 = buffer.data(sig + 114);
    const auto *sig_115 = buffer.data(sig + 115);
    const auto *sig_116 = buffer.data(sig + 116);
    const auto *sig_117 = buffer.data(sig + 117);
    const auto *sig_118 = buffer.data(sig + 118);
    const auto *sig_119 = buffer.data(sig + 119);
    const auto *sig_120 = buffer.data(sig + 120);
    const auto *sig_122 = buffer.data(sig + 122);
    const auto *sig_123 = buffer.data(sig + 123);
    const auto *sig_125 = buffer.data(sig + 125);
    const auto *sig_130 = buffer.data(sig + 130);
    const auto *sig_131 = buffer.data(sig + 131);
    const auto *sig_132 = buffer.data(sig + 132);
    const auto *sig_133 = buffer.data(sig + 133);
    const auto *sig_134 = buffer.data(sig + 134);
    const auto *sig_135 = buffer.data(sig + 135);
    const auto *sig_137 = buffer.data(sig + 137);
    const auto *sig_138 = buffer.data(sig + 138);
    const auto *sig_140 = buffer.data(sig + 140);
    const auto *sig_141 = buffer.data(sig + 141);
    const auto *sig_144 = buffer.data(sig + 144);
    const auto *sig_145 = buffer.data(sig + 145);
    const auto *sig_146 = buffer.data(sig + 146);
    const auto *sig_147 = buffer.data(sig + 147);
    const auto *sig_148 = buffer.data(sig + 148);
    const auto *sig_149 = buffer.data(sig + 149);
    const auto *sig_150 = buffer.data(sig + 150);
    const auto *sig_152 = buffer.data(sig + 152);
    const auto *sig_153 = buffer.data(sig + 153);
    const auto *sig_155 = buffer.data(sig + 155);
    const auto *sig_156 = buffer.data(sig + 156);
    const auto *sig_159 = buffer.data(sig + 159);
    const auto *sig_160 = buffer.data(sig + 160);
    const auto *sig_161 = buffer.data(sig + 161);
    const auto *sig_162 = buffer.data(sig + 162);
    const auto *sig_163 = buffer.data(sig + 163);
    const auto *sig_164 = buffer.data(sig + 164);
    const auto *sig_165 = buffer.data(sig + 165);
    const auto *sig_167 = buffer.data(sig + 167);
    const auto *sig_168 = buffer.data(sig + 168);
    const auto *sig_170 = buffer.data(sig + 170);
    const auto *sig_174 = buffer.data(sig + 174);
    const auto *sig_175 = buffer.data(sig + 175);
    const auto *sig_176 = buffer.data(sig + 176);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, shg_45, shg_90, shg_93, \
                         sif0_60, sif0_63, sif1_60, sif1_63, sig_90, \
                         sig_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_11 * shg_90[k]
                   + f_1 * sif0_60[k]
                   - f_2 * sif1_60[k]
                   + f_3 * pc_x[k] * sig_90[k];

        t_127[k] = f_11 * shg_45[k]
                   + f_3 * pc_y[k] * sig_90[k];

        t_128[k] = f_3 * pc_z[k] * sig_90[k];

        t_129[k] = f_11 * shg_93[k]
                   + f_4 * sif0_63[k]
                   - f_5 * sif1_63[k]
                   + f_3 * pc_x[k] * sig_93[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pc_x, pc_y, shg_47, shg_95, shg_96, sif0_65, \
                         sif0_66, sif1_65, sif1_66, sig_92, sig_95, \
                         sig_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_11 * shg_47[k]
                   + f_3 * pc_y[k] * sig_92[k];

        t_131[k] = f_11 * shg_95[k]
                   + f_4 * sif0_65[k]
                   - f_5 * sif1_65[k]
                   + f_3 * pc_x[k] * sig_95[k];

        t_132[k] = f_11 * shg_96[k]
                   + f_6 * sif0_66[k]
                   - f_7 * sif1_66[k]
                   + f_3 * pc_x[k] * sig_96[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, pc_z, shg_50, shg_99, \
                         shg_100, sif0_69, sif1_69, sig_93, sig_95, sig_99, \
                         sig_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_z[k] * sig_93[k];

        t_134[k] = f_11 * shg_50[k]
                   + f_3 * pc_y[k] * sig_95[k];

        t_135[k] = f_11 * shg_99[k]
                   + f_6 * sif0_69[k]
                   - f_7 * sif1_69[k]
                   + f_3 * pc_x[k] * sig_99[k];

        t_136[k] = f_11 * shg_100[k]
                   + f_3 * pc_x[k] * sig_100[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, shg_101, shg_102, shg_103, shg_104, \
                         sig_101, sig_102, sig_103, sig_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_11 * shg_101[k]
                   + f_3 * pc_x[k] * sig_101[k];

        t_138[k] = f_11 * shg_102[k]
                   + f_3 * pc_x[k] * sig_102[k];

        t_139[k] = f_11 * shg_103[k]
                   + f_3 * pc_x[k] * sig_103[k];

        t_140[k] = f_11 * shg_104[k]
                   + f_3 * pc_x[k] * sig_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_y, pc_z, shg_55, shg_57, sif0_66, sif0_68, \
                         sif1_66, sif1_68, sig_100, sig_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_11 * shg_55[k]
                   + f_1 * sif0_66[k]
                   - f_2 * sif1_66[k]
                   + f_3 * pc_y[k] * sig_100[k];

        t_142[k] = f_3 * pc_z[k] * sig_100[k];

        t_143[k] = f_11 * shg_57[k]
                   + f_4 * sif0_68[k]
                   - f_5 * sif1_68[k]
                   + f_3 * pc_y[k] * sig_102[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_z, pc_y, pc_z, shh0_63, shg_58, \
                         shg_59, shh1_63, sif0_69, sif1_69, sig_103, \
                         sig_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_11 * shg_58[k]
                   + f_6 * sif0_69[k]
                   - f_7 * sif1_69[k]
                   + f_3 * pc_y[k] * sig_103[k];

        t_145[k] = f_11 * shg_59[k]
                   + f_3 * pc_y[k] * sig_104[k];

        t_146[k] = f_1 * sif0_69[k]
                   - f_2 * sif1_69[k]
                   + f_3 * pc_z[k] * sig_104[k];

        t_147[k] = pb_z[k] * shh0_63[k]
                   - f_8 * pc_z[k] * shh1_63[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_z, pc_y, pc_z, shh0_66, shg_45, \
                         shg_60, shg_62, shh1_66, sig_105, sig_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * shg_60[k]
                   + f_3 * pc_y[k] * sig_105[k];

        t_149[k] = f_9 * shg_45[k]
                   + f_3 * pc_z[k] * sig_105[k];

        t_150[k] = pb_z[k] * shh0_66[k]
                   - f_8 * pc_z[k] * shh1_66[k];

        t_151[k] = f_10 * shg_62[k]
                   + f_3 * pc_y[k] * sig_107[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_z, pc_x, pc_z, shh0_69, shg_48, shg_110, \
                         shh1_69, sif0_75, sif1_75, sig_108, sig_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_11 * shg_110[k]
                   + f_4 * sif0_75[k]
                   - f_5 * sif1_75[k]
                   + f_3 * pc_x[k] * sig_110[k];

        t_153[k] = pb_z[k] * shh0_69[k]
                   - f_8 * pc_z[k] * shh1_69[k];

        t_154[k] = f_9 * shg_48[k]
                   + f_3 * pc_z[k] * sig_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_y, shg_65, shg_114, shg_115, \
                         shg_116, sif0_79, sif1_79, sig_110, sig_114, sig_115, \
                         sig_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_10 * shg_65[k]
                   + f_3 * pc_y[k] * sig_110[k];

        t_156[k] = f_11 * shg_114[k]
                   + f_6 * sif0_79[k]
                   - f_7 * sif1_79[k]
                   + f_3 * pc_x[k] * sig_114[k];

        t_157[k] = f_11 * shg_115[k]
                   + f_3 * pc_x[k] * sig_115[k];

        t_158[k] = f_11 * shg_116[k]
                   + f_3 * pc_x[k] * sig_116[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_z, pc_x, pc_z, shh0_78, shg_117, \
                         shg_118, shg_119, shh1_78, sig_117, sig_118, \
                         sig_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_11 * shg_117[k]
                   + f_3 * pc_x[k] * sig_117[k];

        t_160[k] = f_11 * shg_118[k]
                   + f_3 * pc_x[k] * sig_118[k];

        t_161[k] = f_11 * shg_119[k]
                   + f_3 * pc_x[k] * sig_119[k];

        t_162[k] = pb_z[k] * shh0_78[k]
                   - f_8 * pc_z[k] * shh1_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pc_y, pc_z, shg_55, shg_72, shg_73, sif0_78, \
                         sif0_79, sif1_78, sif1_79, sig_115, sig_117, \
                         sig_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * shg_55[k]
                   + f_3 * pc_z[k] * sig_115[k];

        t_164[k] = f_10 * shg_72[k]
                   + f_4 * sif0_78[k]
                   - f_5 * sif1_78[k]
                   + f_3 * pc_y[k] * sig_117[k];

        t_165[k] = f_10 * shg_73[k]
                   + f_6 * sif0_79[k]
                   - f_7 * sif1_79[k]
                   + f_3 * pc_y[k] * sig_118[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_y, pc_y, pc_z, shh0_105, shg_59, \
                         shg_74, shg_75, shh1_105, sif0_79, sif1_79, sig_119, \
                         sig_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * shg_74[k]
                   + f_3 * pc_y[k] * sig_119[k];

        t_167[k] = f_9 * shg_59[k]
                   + f_1 * sif0_79[k]
                   - f_2 * sif1_79[k]
                   + f_3 * pc_z[k] * sig_119[k];

        t_168[k] = pb_y[k] * shh0_105[k]
                   - f_8 * pc_y[k] * shh1_105[k];

        t_169[k] = f_9 * shg_75[k]
                   + f_3 * pc_y[k] * sig_120[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_y, pc_y, pc_z, shh0_108, shh0_110, \
                         shg_60, shg_76, shg_77, shh1_108, shh1_110, sig_120, \
                         sig_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * shg_60[k]
                   + f_3 * pc_z[k] * sig_120[k];

        t_171[k] = pb_y[k] * shh0_108[k]
                   + f_10 * shg_76[k]
                   - f_8 * pc_y[k] * shh1_108[k];

        t_172[k] = f_9 * shg_77[k]
                   + f_3 * pc_y[k] * sig_122[k];

        t_173[k] = pb_y[k] * shh0_110[k]
                   - f_8 * pc_y[k] * shh1_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_y, pc_y, pc_z, shh0_111, shh0_114, \
                         shg_63, shg_78, shg_80, shh1_111, shh1_114, sig_123, \
                         sig_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_y[k] * shh0_111[k]
                   + f_11 * shg_78[k]
                   - f_8 * pc_y[k] * shh1_111[k];

        t_175[k] = f_10 * shg_63[k]
                   + f_3 * pc_z[k] * sig_123[k];

        t_176[k] = f_9 * shg_80[k]
                   + f_3 * pc_y[k] * sig_125[k];

        t_177[k] = pb_y[k] * shh0_114[k]
                   - f_8 * pc_y[k] * shh1_114[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pc_x, shg_130, shg_131, shg_132, \
                         shg_133, shg_134, sig_130, sig_131, sig_132, sig_133, \
                         sig_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_11 * shg_130[k]
                   + f_3 * pc_x[k] * sig_130[k];

        t_179[k] = f_11 * shg_131[k]
                   + f_3 * pc_x[k] * sig_131[k];

        t_180[k] = f_11 * shg_132[k]
                   + f_3 * pc_x[k] * sig_132[k];

        t_181[k] = f_11 * shg_133[k]
                   + f_3 * pc_x[k] * sig_133[k];

        t_182[k] = f_11 * shg_134[k]
                   + f_3 * pc_x[k] * sig_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, shg_70, shg_85, shg_87, sif0_86, \
                         sif0_88, sif1_86, sif1_88, sig_130, sig_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * shg_85[k]
                   + f_1 * sif0_86[k]
                   - f_2 * sif1_86[k]
                   + f_3 * pc_y[k] * sig_130[k];

        t_184[k] = f_10 * shg_70[k]
                   + f_3 * pc_z[k] * sig_130[k];

        t_185[k] = f_9 * shg_87[k]
                   + f_4 * sif0_88[k]
                   - f_5 * sif1_88[k]
                   + f_3 * pc_y[k] * sig_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_y, pc_y, shh0_125, shg_88, shg_89, shh1_125, \
                         sif0_89, sif1_89, sig_133, sig_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * shg_88[k]
                   + f_6 * sif0_89[k]
                   - f_7 * sif1_89[k]
                   + f_3 * pc_y[k] * sig_133[k];

        t_187[k] = f_9 * shg_89[k]
                   + f_3 * pc_y[k] * sig_134[k];

        t_188[k] = pb_y[k] * shh0_125[k]
                   - f_8 * pc_y[k] * shh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_x, pc_y, pc_z, shg_75, shg_135, \
                         shg_138, sif0_90, sif0_93, sif1_90, sif1_93, sig_135, \
                         sig_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_11 * shg_135[k]
                   + f_1 * sif0_90[k]
                   - f_2 * sif1_90[k]
                   + f_3 * pc_x[k] * sig_135[k];

        t_190[k] = f_3 * pc_y[k] * sig_135[k];

        t_191[k] = f_11 * shg_75[k]
                   + f_3 * pc_z[k] * sig_135[k];

        t_192[k] = f_11 * shg_138[k]
                   + f_4 * sif0_93[k]
                   - f_5 * sif1_93[k]
                   + f_3 * pc_x[k] * sig_138[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, shg_140, shg_141, sif0_95, sif0_96, \
                         sif1_95, sif1_96, sig_137, sig_140, sig_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * sig_137[k];

        t_194[k] = f_11 * shg_140[k]
                   + f_4 * sif0_95[k]
                   - f_5 * sif1_95[k]
                   + f_3 * pc_x[k] * sig_140[k];

        t_195[k] = f_11 * shg_141[k]
                   + f_6 * sif0_96[k]
                   - f_7 * sif1_96[k]
                   + f_3 * pc_x[k] * sig_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pc_x, pc_y, pc_z, shg_78, shg_144, \
                         shg_145, sif0_99, sif1_99, sig_138, sig_140, sig_144, \
                         sig_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_11 * shg_78[k]
                   + f_3 * pc_z[k] * sig_138[k];

        t_197[k] = f_3 * pc_y[k] * sig_140[k];

        t_198[k] = f_11 * shg_144[k]
                   + f_6 * sif0_99[k]
                   - f_7 * sif1_99[k]
                   + f_3 * pc_x[k] * sig_144[k];

        t_199[k] = f_11 * shg_145[k]
                   + f_3 * pc_x[k] * sig_145[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, shg_146, shg_147, shg_148, shg_149, \
                         sig_146, sig_147, sig_148, sig_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_11 * shg_146[k]
                   + f_3 * pc_x[k] * sig_146[k];

        t_201[k] = f_11 * shg_147[k]
                   + f_3 * pc_x[k] * sig_147[k];

        t_202[k] = f_11 * shg_148[k]
                   + f_3 * pc_x[k] * sig_148[k];

        t_203[k] = f_11 * shg_149[k]
                   + f_3 * pc_x[k] * sig_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_y, pc_z, shg_85, sif0_96, sif0_98, \
                         sif0_99, sif1_96, sif1_98, sif1_99, sig_145, sig_147, \
                         sig_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * sif0_96[k]
                   - f_2 * sif1_96[k]
                   + f_3 * pc_y[k] * sig_145[k];

        t_205[k] = f_11 * shg_85[k]
                   + f_3 * pc_z[k] * sig_145[k];

        t_206[k] = f_4 * sif0_98[k]
                   - f_5 * sif1_98[k]
                   + f_3 * pc_y[k] * sig_147[k];

        t_207[k] = f_6 * sif0_99[k]
                   - f_7 * sif1_99[k]
                   + f_3 * pc_y[k] * sig_148[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_x, pc_y, pc_z, shg_89, shg_90, \
                         shg_150, sif0_99, sif0_100, sif1_99, sif1_100, sig_149, \
                         sig_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * sig_149[k];

        t_209[k] = f_11 * shg_89[k]
                   + f_1 * sif0_99[k]
                   - f_2 * sif1_99[k]
                   + f_3 * pc_z[k] * sig_149[k];

        t_210[k] = f_10 * shg_150[k]
                   + f_1 * sif0_100[k]
                   - f_2 * sif1_100[k]
                   + f_3 * pc_x[k] * sig_150[k];

        t_211[k] = f_13 * shg_90[k]
                   + f_3 * pc_y[k] * sig_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, shg_92, shg_153, sif0_103, \
                         sif1_103, sig_150, sig_152, sig_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_z[k] * sig_150[k];

        t_213[k] = f_10 * shg_153[k]
                   + f_4 * sif0_103[k]
                   - f_5 * sif1_103[k]
                   + f_3 * pc_x[k] * sig_153[k];

        t_214[k] = f_13 * shg_92[k]
                   + f_3 * pc_y[k] * sig_152[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pc_x, pc_z, shg_155, shg_156, sif0_105, \
                         sif0_106, sif1_105, sif1_106, sig_153, sig_155, \
                         sig_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_10 * shg_155[k]
                   + f_4 * sif0_105[k]
                   - f_5 * sif1_105[k]
                   + f_3 * pc_x[k] * sig_155[k];

        t_216[k] = f_10 * shg_156[k]
                   + f_6 * sif0_106[k]
                   - f_7 * sif1_106[k]
                   + f_3 * pc_x[k] * sig_156[k];

        t_217[k] = f_3 * pc_z[k] * sig_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, shg_95, shg_159, shg_160, \
                         shg_161, sif0_109, sif1_109, sig_155, sig_159, sig_160, \
                         sig_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_13 * shg_95[k]
                   + f_3 * pc_y[k] * sig_155[k];

        t_219[k] = f_10 * shg_159[k]
                   + f_6 * sif0_109[k]
                   - f_7 * sif1_109[k]
                   + f_3 * pc_x[k] * sig_159[k];

        t_220[k] = f_10 * shg_160[k]
                   + f_3 * pc_x[k] * sig_160[k];

        t_221[k] = f_10 * shg_161[k]
                   + f_3 * pc_x[k] * sig_161[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, shg_100, shg_162, shg_163, \
                         shg_164, sif0_106, sif1_106, sig_160, sig_162, sig_163, \
                         sig_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_10 * shg_162[k]
                   + f_3 * pc_x[k] * sig_162[k];

        t_223[k] = f_10 * shg_163[k]
                   + f_3 * pc_x[k] * sig_163[k];

        t_224[k] = f_10 * shg_164[k]
                   + f_3 * pc_x[k] * sig_164[k];

        t_225[k] = f_13 * shg_100[k]
                   + f_1 * sif0_106[k]
                   - f_2 * sif1_106[k]
                   + f_3 * pc_y[k] * sig_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_y, pc_z, shg_102, shg_103, sif0_108, \
                         sif0_109, sif1_108, sif1_109, sig_160, sig_162, \
                         sig_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * sig_160[k];

        t_227[k] = f_13 * shg_102[k]
                   + f_4 * sif0_108[k]
                   - f_5 * sif1_108[k]
                   + f_3 * pc_y[k] * sig_162[k];

        t_228[k] = f_13 * shg_103[k]
                   + f_6 * sif0_109[k]
                   - f_7 * sif1_109[k]
                   + f_3 * pc_y[k] * sig_163[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_z, pc_y, pc_z, shh0_126, shg_104, \
                         shg_105, shh1_126, sif0_109, sif1_109, sig_164, \
                         sig_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_13 * shg_104[k]
                   + f_3 * pc_y[k] * sig_164[k];

        t_230[k] = f_1 * sif0_109[k]
                   - f_2 * sif1_109[k]
                   + f_3 * pc_z[k] * sig_164[k];

        t_231[k] = pb_z[k] * shh0_126[k]
                   - f_8 * pc_z[k] * shh1_126[k];

        t_232[k] = f_11 * shg_105[k]
                   + f_3 * pc_y[k] * sig_165[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_z, pc_y, pc_z, shh0_129, shg_90, shg_107, \
                         shh1_129, sig_165, sig_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_9 * shg_90[k]
                   + f_3 * pc_z[k] * sig_165[k];

        t_234[k] = pb_z[k] * shh0_129[k]
                   - f_8 * pc_z[k] * shh1_129[k];

        t_235[k] = f_11 * shg_107[k]
                   + f_3 * pc_y[k] * sig_167[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pb_z, pc_x, pc_z, shh0_132, shg_93, shg_170, \
                         shh1_132, sif0_115, sif1_115, sig_168, \
                         sig_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_10 * shg_170[k]
                   + f_4 * sif0_115[k]
                   - f_5 * sif1_115[k]
                   + f_3 * pc_x[k] * sig_170[k];

        t_237[k] = pb_z[k] * shh0_132[k]
                   - f_8 * pc_z[k] * shh1_132[k];

        t_238[k] = f_9 * shg_93[k]
                   + f_3 * pc_z[k] * sig_168[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pc_x, pc_y, shg_110, shg_174, shg_175, \
                         shg_176, sif0_119, sif1_119, sig_170, sig_174, sig_175, \
                         sig_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_11 * shg_110[k]
                   + f_3 * pc_y[k] * sig_170[k];

        t_240[k] = f_10 * shg_174[k]
                   + f_6 * sif0_119[k]
                   - f_7 * sif1_119[k]
                   + f_3 * pc_x[k] * sig_174[k];

        t_241[k] = f_10 * shg_175[k]
                   + f_3 * pc_x[k] * sig_175[k];

        t_242[k] = f_10 * shg_176[k]
                   + f_3 * pc_x[k] * sig_176[k];
    }
}

static auto
compute_prim_sih_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shh0,
                                                          const size_t shg, const size_t shh1,
                                                          const size_t sif0, const size_t sif1,
                                                          const size_t sig, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_13 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shh0_141 = buffer.data(shh0 + 141);
    const auto *shh0_189 = buffer.data(shh0 + 189);
    const auto *shh0_192 = buffer.data(shh0 + 192);
    const auto *shh0_194 = buffer.data(shh0 + 194);
    const auto *shh0_195 = buffer.data(shh0 + 195);
    const auto *shh0_198 = buffer.data(shh0 + 198);
    const auto *shh0_209 = buffer.data(shh0 + 209);
    const auto *shh0_210 = buffer.data(shh0 + 210);
    const auto *shh0_213 = buffer.data(shh0 + 213);
    const auto *shh0_216 = buffer.data(shh0 + 216);
    const auto *shh0_315 = buffer.data(shh0 + 315);
    const auto *shh0_318 = buffer.data(shh0 + 318);
    const auto *shh0_320 = buffer.data(shh0 + 320);
    const auto *shh0_321 = buffer.data(shh0 + 321);
    const auto *shh0_324 = buffer.data(shh0 + 324);
    const auto *shh0_330 = buffer.data(shh0 + 330);
    const auto *shh0_332 = buffer.data(shh0 + 332);
    const auto *shh0_333 = buffer.data(shh0 + 333);
    const auto *shh0_335 = buffer.data(shh0 + 335);
    const auto *shh0_341 = buffer.data(shh0 + 341);
    const auto *shh0_345 = buffer.data(shh0 + 345);
    const auto *shh0_351 = buffer.data(shh0 + 351);
    const auto *shh0_353 = buffer.data(shh0 + 353);
    const auto *shh0_354 = buffer.data(shh0 + 354);
    const auto *shh0_356 = buffer.data(shh0 + 356);
    const auto *shh0_357 = buffer.data(shh0 + 357);
    const auto *shh0_360 = buffer.data(shh0 + 360);
    const auto *shh0_362 = buffer.data(shh0 + 362);

    const auto *shg_100 = buffer.data(shg + 100);
    const auto *shg_104 = buffer.data(shg + 104);
    const auto *shg_105 = buffer.data(shg + 105);
    const auto *shg_108 = buffer.data(shg + 108);
    const auto *shg_115 = buffer.data(shg + 115);
    const auto *shg_117 = buffer.data(shg + 117);
    const auto *shg_118 = buffer.data(shg + 118);
    const auto *shg_119 = buffer.data(shg + 119);
    const auto *shg_120 = buffer.data(shg + 120);
    const auto *shg_122 = buffer.data(shg + 122);
    const auto *shg_123 = buffer.data(shg + 123);
    const auto *shg_125 = buffer.data(shg + 125);
    const auto *shg_130 = buffer.data(shg + 130);
    const auto *shg_132 = buffer.data(shg + 132);
    const auto *shg_133 = buffer.data(shg + 133);
    const auto *shg_134 = buffer.data(shg + 134);
    const auto *shg_135 = buffer.data(shg + 135);
    const auto *shg_136 = buffer.data(shg + 136);
    const auto *shg_137 = buffer.data(shg + 137);
    const auto *shg_138 = buffer.data(shg + 138);
    const auto *shg_140 = buffer.data(shg + 140);
    const auto *shg_145 = buffer.data(shg + 145);
    const auto *shg_147 = buffer.data(shg + 147);
    const auto *shg_148 = buffer.data(shg + 148);
    const auto *shg_149 = buffer.data(shg + 149);
    const auto *shg_150 = buffer.data(shg + 150);
    const auto *shg_152 = buffer.data(shg + 152);
    const auto *shg_153 = buffer.data(shg + 153);
    const auto *shg_155 = buffer.data(shg + 155);
    const auto *shg_160 = buffer.data(shg + 160);
    const auto *shg_164 = buffer.data(shg + 164);
    const auto *shg_165 = buffer.data(shg + 165);
    const auto *shg_167 = buffer.data(shg + 167);
    const auto *shg_170 = buffer.data(shg + 170);
    const auto *shg_177 = buffer.data(shg + 177);
    const auto *shg_178 = buffer.data(shg + 178);
    const auto *shg_179 = buffer.data(shg + 179);
    const auto *shg_180 = buffer.data(shg + 180);
    const auto *shg_182 = buffer.data(shg + 182);
    const auto *shg_183 = buffer.data(shg + 183);
    const auto *shg_185 = buffer.data(shg + 185);
    const auto *shg_186 = buffer.data(shg + 186);
    const auto *shg_189 = buffer.data(shg + 189);
    const auto *shg_190 = buffer.data(shg + 190);
    const auto *shg_191 = buffer.data(shg + 191);
    const auto *shg_192 = buffer.data(shg + 192);
    const auto *shg_193 = buffer.data(shg + 193);
    const auto *shg_194 = buffer.data(shg + 194);
    const auto *shg_205 = buffer.data(shg + 205);
    const auto *shg_206 = buffer.data(shg + 206);
    const auto *shg_207 = buffer.data(shg + 207);
    const auto *shg_208 = buffer.data(shg + 208);
    const auto *shg_209 = buffer.data(shg + 209);
    const auto *shg_210 = buffer.data(shg + 210);
    const auto *shg_213 = buffer.data(shg + 213);
    const auto *shg_215 = buffer.data(shg + 215);
    const auto *shg_216 = buffer.data(shg + 216);
    const auto *shg_219 = buffer.data(shg + 219);
    const auto *shg_220 = buffer.data(shg + 220);
    const auto *shg_221 = buffer.data(shg + 221);
    const auto *shg_222 = buffer.data(shg + 222);
    const auto *shg_223 = buffer.data(shg + 223);
    const auto *shg_224 = buffer.data(shg + 224);
    const auto *shg_225 = buffer.data(shg + 225);
    const auto *shg_228 = buffer.data(shg + 228);
    const auto *shg_230 = buffer.data(shg + 230);
    const auto *shg_231 = buffer.data(shg + 231);
    const auto *shg_234 = buffer.data(shg + 234);
    const auto *shg_235 = buffer.data(shg + 235);
    const auto *shg_236 = buffer.data(shg + 236);
    const auto *shg_237 = buffer.data(shg + 237);
    const auto *shg_238 = buffer.data(shg + 238);
    const auto *shg_239 = buffer.data(shg + 239);
    const auto *shg_245 = buffer.data(shg + 245);
    const auto *shg_249 = buffer.data(shg + 249);
    const auto *shg_250 = buffer.data(shg + 250);
    const auto *shg_251 = buffer.data(shg + 251);
    const auto *shg_252 = buffer.data(shg + 252);
    const auto *shg_253 = buffer.data(shg + 253);
    const auto *shg_254 = buffer.data(shg + 254);
    const auto *shg_255 = buffer.data(shg + 255);
    const auto *shg_258 = buffer.data(shg + 258);
    const auto *shg_260 = buffer.data(shg + 260);

    const auto *shh1_141 = buffer.data(shh1 + 141);
    const auto *shh1_189 = buffer.data(shh1 + 189);
    const auto *shh1_192 = buffer.data(shh1 + 192);
    const auto *shh1_194 = buffer.data(shh1 + 194);
    const auto *shh1_195 = buffer.data(shh1 + 195);
    const auto *shh1_198 = buffer.data(shh1 + 198);
    const auto *shh1_209 = buffer.data(shh1 + 209);
    const auto *shh1_210 = buffer.data(shh1 + 210);
    const auto *shh1_213 = buffer.data(shh1 + 213);
    const auto *shh1_216 = buffer.data(shh1 + 216);
    const auto *shh1_315 = buffer.data(shh1 + 315);
    const auto *shh1_318 = buffer.data(shh1 + 318);
    const auto *shh1_320 = buffer.data(shh1 + 320);
    const auto *shh1_321 = buffer.data(shh1 + 321);
    const auto *shh1_324 = buffer.data(shh1 + 324);
    const auto *shh1_330 = buffer.data(shh1 + 330);
    const auto *shh1_332 = buffer.data(shh1 + 332);
    const auto *shh1_333 = buffer.data(shh1 + 333);
    const auto *shh1_335 = buffer.data(shh1 + 335);
    const auto *shh1_341 = buffer.data(shh1 + 341);
    const auto *shh1_345 = buffer.data(shh1 + 345);
    const auto *shh1_351 = buffer.data(shh1 + 351);
    const auto *shh1_353 = buffer.data(shh1 + 353);
    const auto *shh1_354 = buffer.data(shh1 + 354);
    const auto *shh1_356 = buffer.data(shh1 + 356);
    const auto *shh1_357 = buffer.data(shh1 + 357);
    const auto *shh1_360 = buffer.data(shh1 + 360);
    const auto *shh1_362 = buffer.data(shh1 + 362);

    const auto *sif0_118 = buffer.data(sif0 + 118);
    const auto *sif0_119 = buffer.data(sif0 + 119);
    const auto *sif0_120 = buffer.data(sif0 + 120);
    const auto *sif0_123 = buffer.data(sif0 + 123);
    const auto *sif0_125 = buffer.data(sif0 + 125);
    const auto *sif0_126 = buffer.data(sif0 + 126);
    const auto *sif0_128 = buffer.data(sif0 + 128);
    const auto *sif0_129 = buffer.data(sif0 + 129);
    const auto *sif0_136 = buffer.data(sif0 + 136);
    const auto *sif0_138 = buffer.data(sif0 + 138);
    const auto *sif0_139 = buffer.data(sif0 + 139);
    const auto *sif0_140 = buffer.data(sif0 + 140);
    const auto *sif0_143 = buffer.data(sif0 + 143);
    const auto *sif0_145 = buffer.data(sif0 + 145);
    const auto *sif0_146 = buffer.data(sif0 + 146);
    const auto *sif0_148 = buffer.data(sif0 + 148);
    const auto *sif0_149 = buffer.data(sif0 + 149);

    const auto *sif1_118 = buffer.data(sif1 + 118);
    const auto *sif1_119 = buffer.data(sif1 + 119);
    const auto *sif1_120 = buffer.data(sif1 + 120);
    const auto *sif1_123 = buffer.data(sif1 + 123);
    const auto *sif1_125 = buffer.data(sif1 + 125);
    const auto *sif1_126 = buffer.data(sif1 + 126);
    const auto *sif1_128 = buffer.data(sif1 + 128);
    const auto *sif1_129 = buffer.data(sif1 + 129);
    const auto *sif1_136 = buffer.data(sif1 + 136);
    const auto *sif1_138 = buffer.data(sif1 + 138);
    const auto *sif1_139 = buffer.data(sif1 + 139);
    const auto *sif1_140 = buffer.data(sif1 + 140);
    const auto *sif1_143 = buffer.data(sif1 + 143);
    const auto *sif1_145 = buffer.data(sif1 + 145);
    const auto *sif1_146 = buffer.data(sif1 + 146);
    const auto *sif1_148 = buffer.data(sif1 + 148);
    const auto *sif1_149 = buffer.data(sif1 + 149);

    const auto *sig_175 = buffer.data(sig + 175);
    const auto *sig_177 = buffer.data(sig + 177);
    const auto *sig_178 = buffer.data(sig + 178);
    const auto *sig_179 = buffer.data(sig + 179);
    const auto *sig_180 = buffer.data(sig + 180);
    const auto *sig_182 = buffer.data(sig + 182);
    const auto *sig_183 = buffer.data(sig + 183);
    const auto *sig_185 = buffer.data(sig + 185);
    const auto *sig_186 = buffer.data(sig + 186);
    const auto *sig_189 = buffer.data(sig + 189);
    const auto *sig_190 = buffer.data(sig + 190);
    const auto *sig_191 = buffer.data(sig + 191);
    const auto *sig_192 = buffer.data(sig + 192);
    const auto *sig_193 = buffer.data(sig + 193);
    const auto *sig_194 = buffer.data(sig + 194);
    const auto *sig_195 = buffer.data(sig + 195);
    const auto *sig_197 = buffer.data(sig + 197);
    const auto *sig_198 = buffer.data(sig + 198);
    const auto *sig_200 = buffer.data(sig + 200);
    const auto *sig_205 = buffer.data(sig + 205);
    const auto *sig_206 = buffer.data(sig + 206);
    const auto *sig_207 = buffer.data(sig + 207);
    const auto *sig_208 = buffer.data(sig + 208);
    const auto *sig_209 = buffer.data(sig + 209);
    const auto *sig_210 = buffer.data(sig + 210);
    const auto *sig_212 = buffer.data(sig + 212);
    const auto *sig_213 = buffer.data(sig + 213);
    const auto *sig_215 = buffer.data(sig + 215);
    const auto *sig_216 = buffer.data(sig + 216);
    const auto *sig_219 = buffer.data(sig + 219);
    const auto *sig_220 = buffer.data(sig + 220);
    const auto *sig_221 = buffer.data(sig + 221);
    const auto *sig_222 = buffer.data(sig + 222);
    const auto *sig_223 = buffer.data(sig + 223);
    const auto *sig_224 = buffer.data(sig + 224);
    const auto *sig_225 = buffer.data(sig + 225);
    const auto *sig_227 = buffer.data(sig + 227);
    const auto *sig_228 = buffer.data(sig + 228);
    const auto *sig_230 = buffer.data(sig + 230);
    const auto *sig_235 = buffer.data(sig + 235);
    const auto *sig_236 = buffer.data(sig + 236);
    const auto *sig_237 = buffer.data(sig + 237);
    const auto *sig_238 = buffer.data(sig + 238);
    const auto *sig_239 = buffer.data(sig + 239);
    const auto *sig_240 = buffer.data(sig + 240);
    const auto *sig_242 = buffer.data(sig + 242);
    const auto *sig_243 = buffer.data(sig + 243);
    const auto *sig_245 = buffer.data(sig + 245);
    const auto *sig_250 = buffer.data(sig + 250);
    const auto *sig_251 = buffer.data(sig + 251);
    const auto *sig_252 = buffer.data(sig + 252);
    const auto *sig_253 = buffer.data(sig + 253);
    const auto *sig_254 = buffer.data(sig + 254);
    const auto *sig_255 = buffer.data(sig + 255);
    const auto *sig_257 = buffer.data(sig + 257);

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_z, pc_x, pc_z, shh0_141, shg_177, \
                         shg_178, shg_179, shh1_141, sig_177, sig_178, \
                         sig_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_10 * shg_177[k]
                   + f_3 * pc_x[k] * sig_177[k];

        t_244[k] = f_10 * shg_178[k]
                   + f_3 * pc_x[k] * sig_178[k];

        t_245[k] = f_10 * shg_179[k]
                   + f_3 * pc_x[k] * sig_179[k];

        t_246[k] = pb_z[k] * shh0_141[k]
                   - f_8 * pc_z[k] * shh1_141[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, pc_z, shg_100, shg_117, shg_118, sif0_118, \
                         sif0_119, sif1_118, sif1_119, sig_175, sig_177, \
                         sig_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_9 * shg_100[k]
                   + f_3 * pc_z[k] * sig_175[k];

        t_248[k] = f_11 * shg_117[k]
                   + f_4 * sif0_118[k]
                   - f_5 * sif1_118[k]
                   + f_3 * pc_y[k] * sig_177[k];

        t_249[k] = f_11 * shg_118[k]
                   + f_6 * sif0_119[k]
                   - f_7 * sif1_119[k]
                   + f_3 * pc_y[k] * sig_178[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, shg_104, shg_119, shg_180, \
                         sif0_119, sif0_120, sif1_119, sif1_120, sig_179, \
                         sig_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_11 * shg_119[k]
                   + f_3 * pc_y[k] * sig_179[k];

        t_251[k] = f_9 * shg_104[k]
                   + f_1 * sif0_119[k]
                   - f_2 * sif1_119[k]
                   + f_3 * pc_z[k] * sig_179[k];

        t_252[k] = f_10 * shg_180[k]
                   + f_1 * sif0_120[k]
                   - f_2 * sif1_120[k]
                   + f_3 * pc_x[k] * sig_180[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, pc_z, shg_105, shg_120, \
                         shg_122, shg_183, sif0_123, sif1_123, sig_180, sig_182, \
                         sig_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * shg_120[k]
                   + f_3 * pc_y[k] * sig_180[k];

        t_254[k] = f_10 * shg_105[k]
                   + f_3 * pc_z[k] * sig_180[k];

        t_255[k] = f_10 * shg_183[k]
                   + f_4 * sif0_123[k]
                   - f_5 * sif1_123[k]
                   + f_3 * pc_x[k] * sig_183[k];

        t_256[k] = f_10 * shg_122[k]
                   + f_3 * pc_y[k] * sig_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, shg_108, shg_185, shg_186, sif0_125, \
                         sif0_126, sif1_125, sif1_126, sig_183, sig_185, \
                         sig_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_10 * shg_185[k]
                   + f_4 * sif0_125[k]
                   - f_5 * sif1_125[k]
                   + f_3 * pc_x[k] * sig_185[k];

        t_258[k] = f_10 * shg_186[k]
                   + f_6 * sif0_126[k]
                   - f_7 * sif1_126[k]
                   + f_3 * pc_x[k] * sig_186[k];

        t_259[k] = f_10 * shg_108[k]
                   + f_3 * pc_z[k] * sig_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, shg_125, shg_189, shg_190, \
                         shg_191, sif0_129, sif1_129, sig_185, sig_189, sig_190, \
                         sig_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * shg_125[k]
                   + f_3 * pc_y[k] * sig_185[k];

        t_261[k] = f_10 * shg_189[k]
                   + f_6 * sif0_129[k]
                   - f_7 * sif1_129[k]
                   + f_3 * pc_x[k] * sig_189[k];

        t_262[k] = f_10 * shg_190[k]
                   + f_3 * pc_x[k] * sig_190[k];

        t_263[k] = f_10 * shg_191[k]
                   + f_3 * pc_x[k] * sig_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, shg_130, shg_192, shg_193, \
                         shg_194, sif0_126, sif1_126, sig_190, sig_192, sig_193, \
                         sig_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * shg_192[k]
                   + f_3 * pc_x[k] * sig_192[k];

        t_265[k] = f_10 * shg_193[k]
                   + f_3 * pc_x[k] * sig_193[k];

        t_266[k] = f_10 * shg_194[k]
                   + f_3 * pc_x[k] * sig_194[k];

        t_267[k] = f_10 * shg_130[k]
                   + f_1 * sif0_126[k]
                   - f_2 * sif1_126[k]
                   + f_3 * pc_y[k] * sig_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, shg_115, shg_132, shg_133, sif0_128, \
                         sif0_129, sif1_128, sif1_129, sig_190, sig_192, \
                         sig_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * shg_115[k]
                   + f_3 * pc_z[k] * sig_190[k];

        t_269[k] = f_10 * shg_132[k]
                   + f_4 * sif0_128[k]
                   - f_5 * sif1_128[k]
                   + f_3 * pc_y[k] * sig_192[k];

        t_270[k] = f_10 * shg_133[k]
                   + f_6 * sif0_129[k]
                   - f_7 * sif1_129[k]
                   + f_3 * pc_y[k] * sig_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pb_y, pc_y, pc_z, shh0_189, shg_119, \
                         shg_134, shg_135, shh1_189, sif0_129, sif1_129, sig_194, \
                         sig_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * shg_134[k]
                   + f_3 * pc_y[k] * sig_194[k];

        t_272[k] = f_10 * shg_119[k]
                   + f_1 * sif0_129[k]
                   - f_2 * sif1_129[k]
                   + f_3 * pc_z[k] * sig_194[k];

        t_273[k] = pb_y[k] * shh0_189[k]
                   - f_8 * pc_y[k] * shh1_189[k];

        t_274[k] = f_9 * shg_135[k]
                   + f_3 * pc_y[k] * sig_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pb_y, pc_y, pc_z, shh0_192, shh0_194, \
                         shg_120, shg_136, shg_137, shh1_192, shh1_194, sig_195, \
                         sig_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * shg_120[k]
                   + f_3 * pc_z[k] * sig_195[k];

        t_276[k] = pb_y[k] * shh0_192[k]
                   + f_10 * shg_136[k]
                   - f_8 * pc_y[k] * shh1_192[k];

        t_277[k] = f_9 * shg_137[k]
                   + f_3 * pc_y[k] * sig_197[k];

        t_278[k] = pb_y[k] * shh0_194[k]
                   - f_8 * pc_y[k] * shh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_y, pc_y, pc_z, shh0_195, shh0_198, \
                         shg_123, shg_138, shg_140, shh1_195, shh1_198, sig_198, \
                         sig_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_y[k] * shh0_195[k]
                   + f_11 * shg_138[k]
                   - f_8 * pc_y[k] * shh1_195[k];

        t_280[k] = f_11 * shg_123[k]
                   + f_3 * pc_z[k] * sig_198[k];

        t_281[k] = f_9 * shg_140[k]
                   + f_3 * pc_y[k] * sig_200[k];

        t_282[k] = pb_y[k] * shh0_198[k]
                   - f_8 * pc_y[k] * shh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, shg_205, shg_206, shg_207, \
                         shg_208, shg_209, sig_205, sig_206, sig_207, sig_208, \
                         sig_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_10 * shg_205[k]
                   + f_3 * pc_x[k] * sig_205[k];

        t_284[k] = f_10 * shg_206[k]
                   + f_3 * pc_x[k] * sig_206[k];

        t_285[k] = f_10 * shg_207[k]
                   + f_3 * pc_x[k] * sig_207[k];

        t_286[k] = f_10 * shg_208[k]
                   + f_3 * pc_x[k] * sig_208[k];

        t_287[k] = f_10 * shg_209[k]
                   + f_3 * pc_x[k] * sig_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, shg_130, shg_145, shg_147, sif0_136, \
                         sif0_138, sif1_136, sif1_138, sig_205, \
                         sig_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * shg_145[k]
                   + f_1 * sif0_136[k]
                   - f_2 * sif1_136[k]
                   + f_3 * pc_y[k] * sig_205[k];

        t_289[k] = f_11 * shg_130[k]
                   + f_3 * pc_z[k] * sig_205[k];

        t_290[k] = f_9 * shg_147[k]
                   + f_4 * sif0_138[k]
                   - f_5 * sif1_138[k]
                   + f_3 * pc_y[k] * sig_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_y, pc_y, shh0_209, shg_148, shg_149, \
                         shh1_209, sif0_139, sif1_139, sig_208, \
                         sig_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * shg_148[k]
                   + f_6 * sif0_139[k]
                   - f_7 * sif1_139[k]
                   + f_3 * pc_y[k] * sig_208[k];

        t_292[k] = f_9 * shg_149[k]
                   + f_3 * pc_y[k] * sig_209[k];

        t_293[k] = pb_y[k] * shh0_209[k]
                   - f_8 * pc_y[k] * shh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pc_x, pc_y, pc_z, shg_135, shg_210, \
                         shg_213, sif0_140, sif0_143, sif1_140, sif1_143, sig_210, \
                         sig_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_10 * shg_210[k]
                   + f_1 * sif0_140[k]
                   - f_2 * sif1_140[k]
                   + f_3 * pc_x[k] * sig_210[k];

        t_295[k] = f_3 * pc_y[k] * sig_210[k];

        t_296[k] = f_13 * shg_135[k]
                   + f_3 * pc_z[k] * sig_210[k];

        t_297[k] = f_10 * shg_213[k]
                   + f_4 * sif0_143[k]
                   - f_5 * sif1_143[k]
                   + f_3 * pc_x[k] * sig_213[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pc_x, pc_y, shg_215, shg_216, sif0_145, \
                         sif0_146, sif1_145, sif1_146, sig_212, sig_215, \
                         sig_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_3 * pc_y[k] * sig_212[k];

        t_299[k] = f_10 * shg_215[k]
                   + f_4 * sif0_145[k]
                   - f_5 * sif1_145[k]
                   + f_3 * pc_x[k] * sig_215[k];

        t_300[k] = f_10 * shg_216[k]
                   + f_6 * sif0_146[k]
                   - f_7 * sif1_146[k]
                   + f_3 * pc_x[k] * sig_216[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_x, pc_y, pc_z, shg_138, shg_219, \
                         shg_220, sif0_149, sif1_149, sig_213, sig_215, sig_219, \
                         sig_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * shg_138[k]
                   + f_3 * pc_z[k] * sig_213[k];

        t_302[k] = f_3 * pc_y[k] * sig_215[k];

        t_303[k] = f_10 * shg_219[k]
                   + f_6 * sif0_149[k]
                   - f_7 * sif1_149[k]
                   + f_3 * pc_x[k] * sig_219[k];

        t_304[k] = f_10 * shg_220[k]
                   + f_3 * pc_x[k] * sig_220[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, shg_221, shg_222, shg_223, shg_224, \
                         sig_221, sig_222, sig_223, sig_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_10 * shg_221[k]
                   + f_3 * pc_x[k] * sig_221[k];

        t_306[k] = f_10 * shg_222[k]
                   + f_3 * pc_x[k] * sig_222[k];

        t_307[k] = f_10 * shg_223[k]
                   + f_3 * pc_x[k] * sig_223[k];

        t_308[k] = f_10 * shg_224[k]
                   + f_3 * pc_x[k] * sig_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pc_y, pc_z, shg_145, sif0_146, sif0_148, \
                         sif0_149, sif1_146, sif1_148, sif1_149, sig_220, sig_222, \
                         sig_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * sif0_146[k]
                   - f_2 * sif1_146[k]
                   + f_3 * pc_y[k] * sig_220[k];

        t_310[k] = f_13 * shg_145[k]
                   + f_3 * pc_z[k] * sig_220[k];

        t_311[k] = f_4 * sif0_148[k]
                   - f_5 * sif1_148[k]
                   + f_3 * pc_y[k] * sig_222[k];

        t_312[k] = f_6 * sif0_149[k]
                   - f_7 * sif1_149[k]
                   + f_3 * pc_y[k] * sig_223[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pb_x, pc_x, pc_y, pc_z, shh0_315, shg_149, \
                         shg_225, shh1_315, sif0_149, sif1_149, \
                         sig_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * sig_224[k];

        t_314[k] = f_13 * shg_149[k]
                   + f_1 * sif0_149[k]
                   - f_2 * sif1_149[k]
                   + f_3 * pc_z[k] * sig_224[k];

        t_315[k] = pb_x[k] * shh0_315[k]
                   + f_12 * shg_225[k]
                   - f_8 * pc_x[k] * shh1_315[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pb_x, pc_x, pc_y, pc_z, shh0_318, \
                         shg_150, shg_152, shg_228, shh1_318, sig_225, \
                         sig_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_12 * shg_150[k]
                   + f_3 * pc_y[k] * sig_225[k];

        t_317[k] = f_3 * pc_z[k] * sig_225[k];

        t_318[k] = pb_x[k] * shh0_318[k]
                   + f_11 * shg_228[k]
                   - f_8 * pc_x[k] * shh1_318[k];

        t_319[k] = f_12 * shg_152[k]
                   + f_3 * pc_y[k] * sig_227[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pb_x, pc_x, pc_z, shh0_320, shh0_321, shg_230, \
                         shg_231, shh1_320, shh1_321, sig_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = pb_x[k] * shh0_320[k]
                   + f_11 * shg_230[k]
                   - f_8 * pc_x[k] * shh1_320[k];

        t_321[k] = pb_x[k] * shh0_321[k]
                   + f_10 * shg_231[k]
                   - f_8 * pc_x[k] * shh1_321[k];

        t_322[k] = f_3 * pc_z[k] * sig_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pb_x, pc_x, pc_y, shh0_324, shg_155, \
                         shg_234, shg_235, shg_236, shh1_324, sig_230, sig_235, \
                         sig_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_12 * shg_155[k]
                   + f_3 * pc_y[k] * sig_230[k];

        t_324[k] = pb_x[k] * shh0_324[k]
                   + f_10 * shg_234[k]
                   - f_8 * pc_x[k] * shh1_324[k];

        t_325[k] = f_9 * shg_235[k]
                   + f_3 * pc_x[k] * sig_235[k];

        t_326[k] = f_9 * shg_236[k]
                   + f_3 * pc_x[k] * sig_236[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pb_x, pc_x, shh0_330, shg_237, shg_238, \
                         shg_239, shh1_330, sig_237, sig_238, sig_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_9 * shg_237[k]
                   + f_3 * pc_x[k] * sig_237[k];

        t_328[k] = f_9 * shg_238[k]
                   + f_3 * pc_x[k] * sig_238[k];

        t_329[k] = f_9 * shg_239[k]
                   + f_3 * pc_x[k] * sig_239[k];

        t_330[k] = pb_x[k] * shh0_330[k]
                   - f_8 * pc_x[k] * shh1_330[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_x, pc_x, pc_y, pc_z, shh0_332, \
                         shh0_333, shg_164, shh1_332, shh1_333, sig_235, \
                         sig_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * sig_235[k];

        t_332[k] = pb_x[k] * shh0_332[k]
                   - f_8 * pc_x[k] * shh1_332[k];

        t_333[k] = pb_x[k] * shh0_333[k]
                   - f_8 * pc_x[k] * shh1_333[k];

        t_334[k] = f_12 * shg_164[k]
                   + f_3 * pc_y[k] * sig_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pb_x, pb_z, pc_x, pc_y, pc_z, shh0_210, \
                         shh0_335, shg_150, shg_165, shh1_210, shh1_335, \
                         sig_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pb_x[k] * shh0_335[k]
                   - f_8 * pc_x[k] * shh1_335[k];

        t_336[k] = pb_z[k] * shh0_210[k]
                   - f_8 * pc_z[k] * shh1_210[k];

        t_337[k] = f_13 * shg_165[k]
                   + f_3 * pc_y[k] * sig_240[k];

        t_338[k] = f_9 * shg_150[k]
                   + f_3 * pc_z[k] * sig_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pb_x, pb_z, pc_x, pc_y, pc_z, shh0_213, \
                         shh0_341, shg_167, shg_245, shh1_213, shh1_341, \
                         sig_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pb_z[k] * shh0_213[k]
                   - f_8 * pc_z[k] * shh1_213[k];

        t_340[k] = f_13 * shg_167[k]
                   + f_3 * pc_y[k] * sig_242[k];

        t_341[k] = pb_x[k] * shh0_341[k]
                   + f_11 * shg_245[k]
                   - f_8 * pc_x[k] * shh1_341[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pb_z, pc_y, pc_z, shh0_216, shg_153, shg_170, \
                         shh1_216, sig_243, sig_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pb_z[k] * shh0_216[k]
                   - f_8 * pc_z[k] * shh1_216[k];

        t_343[k] = f_9 * shg_153[k]
                   + f_3 * pc_z[k] * sig_243[k];

        t_344[k] = f_13 * shg_170[k]
                   + f_3 * pc_y[k] * sig_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pb_x, pc_x, shh0_345, shg_249, shg_250, \
                         shg_251, shg_252, shh1_345, sig_250, sig_251, \
                         sig_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = pb_x[k] * shh0_345[k]
                   + f_10 * shg_249[k]
                   - f_8 * pc_x[k] * shh1_345[k];

        t_346[k] = f_9 * shg_250[k]
                   + f_3 * pc_x[k] * sig_250[k];

        t_347[k] = f_9 * shg_251[k]
                   + f_3 * pc_x[k] * sig_251[k];

        t_348[k] = f_9 * shg_252[k]
                   + f_3 * pc_x[k] * sig_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pb_x, pc_x, pc_z, shh0_351, shg_160, \
                         shg_253, shg_254, shh1_351, sig_250, sig_253, \
                         sig_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_9 * shg_253[k]
                   + f_3 * pc_x[k] * sig_253[k];

        t_350[k] = f_9 * shg_254[k]
                   + f_3 * pc_x[k] * sig_254[k];

        t_351[k] = pb_x[k] * shh0_351[k]
                   - f_8 * pc_x[k] * shh1_351[k];

        t_352[k] = f_9 * shg_160[k]
                   + f_3 * pc_z[k] * sig_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, pb_x, pc_x, pc_y, shh0_353, shh0_354, \
                         shh0_356, shg_179, shh1_353, shh1_354, shh1_356, \
                         sig_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = pb_x[k] * shh0_353[k]
                   - f_8 * pc_x[k] * shh1_353[k];

        t_354[k] = pb_x[k] * shh0_354[k]
                   - f_8 * pc_x[k] * shh1_354[k];

        t_355[k] = f_13 * shg_179[k]
                   + f_3 * pc_y[k] * sig_254[k];

        t_356[k] = pb_x[k] * shh0_356[k]
                   - f_8 * pc_x[k] * shh1_356[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pb_x, pc_x, pc_y, pc_z, shh0_357, shg_165, \
                         shg_180, shg_255, shh1_357, sig_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pb_x[k] * shh0_357[k]
                   + f_12 * shg_255[k]
                   - f_8 * pc_x[k] * shh1_357[k];

        t_358[k] = f_11 * shg_180[k]
                   + f_3 * pc_y[k] * sig_255[k];

        t_359[k] = f_10 * shg_165[k]
                   + f_3 * pc_z[k] * sig_255[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pb_x, pc_x, pc_y, shh0_360, shh0_362, shg_182, \
                         shg_258, shg_260, shh1_360, shh1_362, \
                         sig_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pb_x[k] * shh0_360[k]
                   + f_11 * shg_258[k]
                   - f_8 * pc_x[k] * shh1_360[k];

        t_361[k] = f_11 * shg_182[k]
                   + f_3 * pc_y[k] * sig_257[k];

        t_362[k] = pb_x[k] * shh0_362[k]
                   + f_11 * shg_260[k]
                   - f_8 * pc_x[k] * shh1_362[k];
    }
}

static auto
compute_prim_sih_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shh0,
                                                          const size_t shg, const size_t shh1,
                                                          const size_t sif0, const size_t sif1,
                                                          const size_t sig, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_13 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shh0_294 = buffer.data(shh0 + 294);
    const auto *shh0_299 = buffer.data(shh0 + 299);
    const auto *shh0_303 = buffer.data(shh0 + 303);
    const auto *shh0_315 = buffer.data(shh0 + 315);
    const auto *shh0_318 = buffer.data(shh0 + 318);
    const auto *shh0_321 = buffer.data(shh0 + 321);
    const auto *shh0_330 = buffer.data(shh0 + 330);
    const auto *shh0_332 = buffer.data(shh0 + 332);
    const auto *shh0_333 = buffer.data(shh0 + 333);
    const auto *shh0_363 = buffer.data(shh0 + 363);
    const auto *shh0_366 = buffer.data(shh0 + 366);
    const auto *shh0_372 = buffer.data(shh0 + 372);
    const auto *shh0_374 = buffer.data(shh0 + 374);
    const auto *shh0_375 = buffer.data(shh0 + 375);
    const auto *shh0_377 = buffer.data(shh0 + 377);
    const auto *shh0_378 = buffer.data(shh0 + 378);
    const auto *shh0_381 = buffer.data(shh0 + 381);
    const auto *shh0_383 = buffer.data(shh0 + 383);
    const auto *shh0_384 = buffer.data(shh0 + 384);
    const auto *shh0_387 = buffer.data(shh0 + 387);
    const auto *shh0_393 = buffer.data(shh0 + 393);
    const auto *shh0_395 = buffer.data(shh0 + 395);
    const auto *shh0_396 = buffer.data(shh0 + 396);
    const auto *shh0_398 = buffer.data(shh0 + 398);
    const auto *shh0_402 = buffer.data(shh0 + 402);
    const auto *shh0_405 = buffer.data(shh0 + 405);
    const auto *shh0_414 = buffer.data(shh0 + 414);
    const auto *shh0_416 = buffer.data(shh0 + 416);
    const auto *shh0_417 = buffer.data(shh0 + 417);
    const auto *shh0_419 = buffer.data(shh0 + 419);
    const auto *shh0_420 = buffer.data(shh0 + 420);
    const auto *shh0_423 = buffer.data(shh0 + 423);
    const auto *shh0_425 = buffer.data(shh0 + 425);
    const auto *shh0_426 = buffer.data(shh0 + 426);
    const auto *shh0_429 = buffer.data(shh0 + 429);
    const auto *shh0_435 = buffer.data(shh0 + 435);
    const auto *shh0_437 = buffer.data(shh0 + 437);
    const auto *shh0_438 = buffer.data(shh0 + 438);
    const auto *shh0_440 = buffer.data(shh0 + 440);

    const auto *shg_168 = buffer.data(shg + 168);
    const auto *shg_175 = buffer.data(shg + 175);
    const auto *shg_180 = buffer.data(shg + 180);
    const auto *shg_183 = buffer.data(shg + 183);
    const auto *shg_185 = buffer.data(shg + 185);
    const auto *shg_190 = buffer.data(shg + 190);
    const auto *shg_194 = buffer.data(shg + 194);
    const auto *shg_195 = buffer.data(shg + 195);
    const auto *shg_197 = buffer.data(shg + 197);
    const auto *shg_198 = buffer.data(shg + 198);
    const auto *shg_200 = buffer.data(shg + 200);
    const auto *shg_205 = buffer.data(shg + 205);
    const auto *shg_209 = buffer.data(shg + 209);
    const auto *shg_210 = buffer.data(shg + 210);
    const auto *shg_212 = buffer.data(shg + 212);
    const auto *shg_213 = buffer.data(shg + 213);
    const auto *shg_215 = buffer.data(shg + 215);
    const auto *shg_220 = buffer.data(shg + 220);
    const auto *shg_224 = buffer.data(shg + 224);
    const auto *shg_225 = buffer.data(shg + 225);
    const auto *shg_227 = buffer.data(shg + 227);
    const auto *shg_228 = buffer.data(shg + 228);
    const auto *shg_230 = buffer.data(shg + 230);
    const auto *shg_235 = buffer.data(shg + 235);
    const auto *shg_236 = buffer.data(shg + 236);
    const auto *shg_237 = buffer.data(shg + 237);
    const auto *shg_238 = buffer.data(shg + 238);
    const auto *shg_239 = buffer.data(shg + 239);
    const auto *shg_240 = buffer.data(shg + 240);
    const auto *shg_242 = buffer.data(shg + 242);
    const auto *shg_245 = buffer.data(shg + 245);
    const auto *shg_254 = buffer.data(shg + 254);
    const auto *shg_255 = buffer.data(shg + 255);
    const auto *shg_257 = buffer.data(shg + 257);
    const auto *shg_261 = buffer.data(shg + 261);
    const auto *shg_264 = buffer.data(shg + 264);
    const auto *shg_265 = buffer.data(shg + 265);
    const auto *shg_266 = buffer.data(shg + 266);
    const auto *shg_267 = buffer.data(shg + 267);
    const auto *shg_268 = buffer.data(shg + 268);
    const auto *shg_269 = buffer.data(shg + 269);
    const auto *shg_270 = buffer.data(shg + 270);
    const auto *shg_273 = buffer.data(shg + 273);
    const auto *shg_275 = buffer.data(shg + 275);
    const auto *shg_276 = buffer.data(shg + 276);
    const auto *shg_279 = buffer.data(shg + 279);
    const auto *shg_280 = buffer.data(shg + 280);
    const auto *shg_281 = buffer.data(shg + 281);
    const auto *shg_282 = buffer.data(shg + 282);
    const auto *shg_283 = buffer.data(shg + 283);
    const auto *shg_284 = buffer.data(shg + 284);
    const auto *shg_288 = buffer.data(shg + 288);
    const auto *shg_291 = buffer.data(shg + 291);
    const auto *shg_295 = buffer.data(shg + 295);
    const auto *shg_296 = buffer.data(shg + 296);
    const auto *shg_297 = buffer.data(shg + 297);
    const auto *shg_298 = buffer.data(shg + 298);
    const auto *shg_299 = buffer.data(shg + 299);
    const auto *shg_300 = buffer.data(shg + 300);
    const auto *shg_303 = buffer.data(shg + 303);
    const auto *shg_305 = buffer.data(shg + 305);
    const auto *shg_306 = buffer.data(shg + 306);
    const auto *shg_309 = buffer.data(shg + 309);
    const auto *shg_310 = buffer.data(shg + 310);
    const auto *shg_311 = buffer.data(shg + 311);
    const auto *shg_312 = buffer.data(shg + 312);
    const auto *shg_313 = buffer.data(shg + 313);
    const auto *shg_314 = buffer.data(shg + 314);

    const auto *shh1_294 = buffer.data(shh1 + 294);
    const auto *shh1_299 = buffer.data(shh1 + 299);
    const auto *shh1_303 = buffer.data(shh1 + 303);
    const auto *shh1_315 = buffer.data(shh1 + 315);
    const auto *shh1_318 = buffer.data(shh1 + 318);
    const auto *shh1_321 = buffer.data(shh1 + 321);
    const auto *shh1_330 = buffer.data(shh1 + 330);
    const auto *shh1_332 = buffer.data(shh1 + 332);
    const auto *shh1_333 = buffer.data(shh1 + 333);
    const auto *shh1_363 = buffer.data(shh1 + 363);
    const auto *shh1_366 = buffer.data(shh1 + 366);
    const auto *shh1_372 = buffer.data(shh1 + 372);
    const auto *shh1_374 = buffer.data(shh1 + 374);
    const auto *shh1_375 = buffer.data(shh1 + 375);
    const auto *shh1_377 = buffer.data(shh1 + 377);
    const auto *shh1_378 = buffer.data(shh1 + 378);
    const auto *shh1_381 = buffer.data(shh1 + 381);
    const auto *shh1_383 = buffer.data(shh1 + 383);
    const auto *shh1_384 = buffer.data(shh1 + 384);
    const auto *shh1_387 = buffer.data(shh1 + 387);
    const auto *shh1_393 = buffer.data(shh1 + 393);
    const auto *shh1_395 = buffer.data(shh1 + 395);
    const auto *shh1_396 = buffer.data(shh1 + 396);
    const auto *shh1_398 = buffer.data(shh1 + 398);
    const auto *shh1_402 = buffer.data(shh1 + 402);
    const auto *shh1_405 = buffer.data(shh1 + 405);
    const auto *shh1_414 = buffer.data(shh1 + 414);
    const auto *shh1_416 = buffer.data(shh1 + 416);
    const auto *shh1_417 = buffer.data(shh1 + 417);
    const auto *shh1_419 = buffer.data(shh1 + 419);
    const auto *shh1_420 = buffer.data(shh1 + 420);
    const auto *shh1_423 = buffer.data(shh1 + 423);
    const auto *shh1_425 = buffer.data(shh1 + 425);
    const auto *shh1_426 = buffer.data(shh1 + 426);
    const auto *shh1_429 = buffer.data(shh1 + 429);
    const auto *shh1_435 = buffer.data(shh1 + 435);
    const auto *shh1_437 = buffer.data(shh1 + 437);
    const auto *shh1_438 = buffer.data(shh1 + 438);
    const auto *shh1_440 = buffer.data(shh1 + 440);

    const auto *sif0_210 = buffer.data(sif0 + 210);
    const auto *sif0_213 = buffer.data(sif0 + 213);
    const auto *sif0_215 = buffer.data(sif0 + 215);
    const auto *sif0_216 = buffer.data(sif0 + 216);
    const auto *sif0_218 = buffer.data(sif0 + 218);
    const auto *sif0_219 = buffer.data(sif0 + 219);
    const auto *sif0_225 = buffer.data(sif0 + 225);
    const auto *sif0_229 = buffer.data(sif0 + 229);
    const auto *sif0_230 = buffer.data(sif0 + 230);
    const auto *sif0_233 = buffer.data(sif0 + 233);

    const auto *sif1_210 = buffer.data(sif1 + 210);
    const auto *sif1_213 = buffer.data(sif1 + 213);
    const auto *sif1_215 = buffer.data(sif1 + 215);
    const auto *sif1_216 = buffer.data(sif1 + 216);
    const auto *sif1_218 = buffer.data(sif1 + 218);
    const auto *sif1_219 = buffer.data(sif1 + 219);
    const auto *sif1_225 = buffer.data(sif1 + 225);
    const auto *sif1_229 = buffer.data(sif1 + 229);
    const auto *sif1_230 = buffer.data(sif1 + 230);
    const auto *sif1_233 = buffer.data(sif1 + 233);

    const auto *sig_258 = buffer.data(sig + 258);
    const auto *sig_260 = buffer.data(sig + 260);
    const auto *sig_265 = buffer.data(sig + 265);
    const auto *sig_266 = buffer.data(sig + 266);
    const auto *sig_267 = buffer.data(sig + 267);
    const auto *sig_268 = buffer.data(sig + 268);
    const auto *sig_269 = buffer.data(sig + 269);
    const auto *sig_270 = buffer.data(sig + 270);
    const auto *sig_272 = buffer.data(sig + 272);
    const auto *sig_273 = buffer.data(sig + 273);
    const auto *sig_275 = buffer.data(sig + 275);
    const auto *sig_280 = buffer.data(sig + 280);
    const auto *sig_281 = buffer.data(sig + 281);
    const auto *sig_282 = buffer.data(sig + 282);
    const auto *sig_283 = buffer.data(sig + 283);
    const auto *sig_284 = buffer.data(sig + 284);
    const auto *sig_285 = buffer.data(sig + 285);
    const auto *sig_287 = buffer.data(sig + 287);
    const auto *sig_288 = buffer.data(sig + 288);
    const auto *sig_290 = buffer.data(sig + 290);
    const auto *sig_295 = buffer.data(sig + 295);
    const auto *sig_296 = buffer.data(sig + 296);
    const auto *sig_297 = buffer.data(sig + 297);
    const auto *sig_298 = buffer.data(sig + 298);
    const auto *sig_299 = buffer.data(sig + 299);
    const auto *sig_300 = buffer.data(sig + 300);
    const auto *sig_302 = buffer.data(sig + 302);
    const auto *sig_303 = buffer.data(sig + 303);
    const auto *sig_305 = buffer.data(sig + 305);
    const auto *sig_310 = buffer.data(sig + 310);
    const auto *sig_311 = buffer.data(sig + 311);
    const auto *sig_312 = buffer.data(sig + 312);
    const auto *sig_313 = buffer.data(sig + 313);
    const auto *sig_314 = buffer.data(sig + 314);
    const auto *sig_315 = buffer.data(sig + 315);
    const auto *sig_317 = buffer.data(sig + 317);
    const auto *sig_318 = buffer.data(sig + 318);
    const auto *sig_320 = buffer.data(sig + 320);
    const auto *sig_321 = buffer.data(sig + 321);
    const auto *sig_324 = buffer.data(sig + 324);
    const auto *sig_325 = buffer.data(sig + 325);
    const auto *sig_326 = buffer.data(sig + 326);
    const auto *sig_327 = buffer.data(sig + 327);
    const auto *sig_328 = buffer.data(sig + 328);
    const auto *sig_329 = buffer.data(sig + 329);
    const auto *sig_330 = buffer.data(sig + 330);
    const auto *sig_332 = buffer.data(sig + 332);
    const auto *sig_333 = buffer.data(sig + 333);
    const auto *sig_335 = buffer.data(sig + 335);
    const auto *sig_339 = buffer.data(sig + 339);
    const auto *sig_340 = buffer.data(sig + 340);
    const auto *sig_341 = buffer.data(sig + 341);
    const auto *sig_342 = buffer.data(sig + 342);
    const auto *sig_343 = buffer.data(sig + 343);
    const auto *sig_344 = buffer.data(sig + 344);
    const auto *sig_345 = buffer.data(sig + 345);
    const auto *sig_347 = buffer.data(sig + 347);
    const auto *sig_348 = buffer.data(sig + 348);

#pragma omp simd aligned(t_363, t_364, t_365, pb_x, pc_x, pc_y, pc_z, shh0_363, shg_168, \
                         shg_185, shg_261, shh1_363, sig_258, sig_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pb_x[k] * shh0_363[k]
                   + f_10 * shg_261[k]
                   - f_8 * pc_x[k] * shh1_363[k];

        t_364[k] = f_10 * shg_168[k]
                   + f_3 * pc_z[k] * sig_258[k];

        t_365[k] = f_11 * shg_185[k]
                   + f_3 * pc_y[k] * sig_260[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_x, pc_x, shh0_366, shg_264, shg_265, \
                         shg_266, shg_267, shh1_366, sig_265, sig_266, \
                         sig_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = pb_x[k] * shh0_366[k]
                   + f_10 * shg_264[k]
                   - f_8 * pc_x[k] * shh1_366[k];

        t_367[k] = f_9 * shg_265[k]
                   + f_3 * pc_x[k] * sig_265[k];

        t_368[k] = f_9 * shg_266[k]
                   + f_3 * pc_x[k] * sig_266[k];

        t_369[k] = f_9 * shg_267[k]
                   + f_3 * pc_x[k] * sig_267[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pb_x, pc_x, pc_z, shh0_372, shg_175, \
                         shg_268, shg_269, shh1_372, sig_265, sig_268, \
                         sig_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_9 * shg_268[k]
                   + f_3 * pc_x[k] * sig_268[k];

        t_371[k] = f_9 * shg_269[k]
                   + f_3 * pc_x[k] * sig_269[k];

        t_372[k] = pb_x[k] * shh0_372[k]
                   - f_8 * pc_x[k] * shh1_372[k];

        t_373[k] = f_10 * shg_175[k]
                   + f_3 * pc_z[k] * sig_265[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, pb_x, pc_x, pc_y, shh0_374, shh0_375, \
                         shh0_377, shg_194, shh1_374, shh1_375, shh1_377, \
                         sig_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pb_x[k] * shh0_374[k]
                   - f_8 * pc_x[k] * shh1_374[k];

        t_375[k] = pb_x[k] * shh0_375[k]
                   - f_8 * pc_x[k] * shh1_375[k];

        t_376[k] = f_11 * shg_194[k]
                   + f_3 * pc_y[k] * sig_269[k];

        t_377[k] = pb_x[k] * shh0_377[k]
                   - f_8 * pc_x[k] * shh1_377[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pb_x, pc_x, pc_y, pc_z, shh0_378, shg_180, \
                         shg_195, shg_270, shh1_378, sig_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_x[k] * shh0_378[k]
                   + f_12 * shg_270[k]
                   - f_8 * pc_x[k] * shh1_378[k];

        t_379[k] = f_10 * shg_195[k]
                   + f_3 * pc_y[k] * sig_270[k];

        t_380[k] = f_11 * shg_180[k]
                   + f_3 * pc_z[k] * sig_270[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, pb_x, pc_x, pc_y, shh0_381, shh0_383, shg_197, \
                         shg_273, shg_275, shh1_381, shh1_383, \
                         sig_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = pb_x[k] * shh0_381[k]
                   + f_11 * shg_273[k]
                   - f_8 * pc_x[k] * shh1_381[k];

        t_382[k] = f_10 * shg_197[k]
                   + f_3 * pc_y[k] * sig_272[k];

        t_383[k] = pb_x[k] * shh0_383[k]
                   + f_11 * shg_275[k]
                   - f_8 * pc_x[k] * shh1_383[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pb_x, pc_x, pc_y, pc_z, shh0_384, shg_183, \
                         shg_200, shg_276, shh1_384, sig_273, sig_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pb_x[k] * shh0_384[k]
                   + f_10 * shg_276[k]
                   - f_8 * pc_x[k] * shh1_384[k];

        t_385[k] = f_11 * shg_183[k]
                   + f_3 * pc_z[k] * sig_273[k];

        t_386[k] = f_10 * shg_200[k]
                   + f_3 * pc_y[k] * sig_275[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pb_x, pc_x, shh0_387, shg_279, shg_280, \
                         shg_281, shg_282, shh1_387, sig_280, sig_281, \
                         sig_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = pb_x[k] * shh0_387[k]
                   + f_10 * shg_279[k]
                   - f_8 * pc_x[k] * shh1_387[k];

        t_388[k] = f_9 * shg_280[k]
                   + f_3 * pc_x[k] * sig_280[k];

        t_389[k] = f_9 * shg_281[k]
                   + f_3 * pc_x[k] * sig_281[k];

        t_390[k] = f_9 * shg_282[k]
                   + f_3 * pc_x[k] * sig_282[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_x, pc_x, pc_z, shh0_393, shg_190, \
                         shg_283, shg_284, shh1_393, sig_280, sig_283, \
                         sig_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_9 * shg_283[k]
                   + f_3 * pc_x[k] * sig_283[k];

        t_392[k] = f_9 * shg_284[k]
                   + f_3 * pc_x[k] * sig_284[k];

        t_393[k] = pb_x[k] * shh0_393[k]
                   - f_8 * pc_x[k] * shh1_393[k];

        t_394[k] = f_11 * shg_190[k]
                   + f_3 * pc_z[k] * sig_280[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pb_x, pc_x, pc_y, shh0_395, shh0_396, \
                         shh0_398, shg_209, shh1_395, shh1_396, shh1_398, \
                         sig_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pb_x[k] * shh0_395[k]
                   - f_8 * pc_x[k] * shh1_395[k];

        t_396[k] = pb_x[k] * shh0_396[k]
                   - f_8 * pc_x[k] * shh1_396[k];

        t_397[k] = f_10 * shg_209[k]
                   + f_3 * pc_y[k] * sig_284[k];

        t_398[k] = pb_x[k] * shh0_398[k]
                   - f_8 * pc_x[k] * shh1_398[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pb_y, pc_y, pc_z, shh0_294, shg_195, shg_210, \
                         shh1_294, sig_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = pb_y[k] * shh0_294[k]
                   - f_8 * pc_y[k] * shh1_294[k];

        t_400[k] = f_9 * shg_210[k]
                   + f_3 * pc_y[k] * sig_285[k];

        t_401[k] = f_13 * shg_195[k]
                   + f_3 * pc_z[k] * sig_285[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pb_x, pb_y, pc_x, pc_y, shh0_299, shh0_402, \
                         shg_212, shg_288, shh1_299, shh1_402, \
                         sig_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pb_x[k] * shh0_402[k]
                   + f_11 * shg_288[k]
                   - f_8 * pc_x[k] * shh1_402[k];

        t_403[k] = f_9 * shg_212[k]
                   + f_3 * pc_y[k] * sig_287[k];

        t_404[k] = pb_y[k] * shh0_299[k]
                   - f_8 * pc_y[k] * shh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, pb_x, pc_x, pc_y, pc_z, shh0_405, shg_198, \
                         shg_215, shg_291, shh1_405, sig_288, sig_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pb_x[k] * shh0_405[k]
                   + f_10 * shg_291[k]
                   - f_8 * pc_x[k] * shh1_405[k];

        t_406[k] = f_13 * shg_198[k]
                   + f_3 * pc_z[k] * sig_288[k];

        t_407[k] = f_9 * shg_215[k]
                   + f_3 * pc_y[k] * sig_290[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pb_y, pc_x, pc_y, shh0_303, shg_295, \
                         shg_296, shg_297, shh1_303, sig_295, sig_296, \
                         sig_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = pb_y[k] * shh0_303[k]
                   - f_8 * pc_y[k] * shh1_303[k];

        t_409[k] = f_9 * shg_295[k]
                   + f_3 * pc_x[k] * sig_295[k];

        t_410[k] = f_9 * shg_296[k]
                   + f_3 * pc_x[k] * sig_296[k];

        t_411[k] = f_9 * shg_297[k]
                   + f_3 * pc_x[k] * sig_297[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pb_x, pc_x, pc_z, shh0_414, shg_205, \
                         shg_298, shg_299, shh1_414, sig_295, sig_298, \
                         sig_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_9 * shg_298[k]
                   + f_3 * pc_x[k] * sig_298[k];

        t_413[k] = f_9 * shg_299[k]
                   + f_3 * pc_x[k] * sig_299[k];

        t_414[k] = pb_x[k] * shh0_414[k]
                   - f_8 * pc_x[k] * shh1_414[k];

        t_415[k] = f_13 * shg_205[k]
                   + f_3 * pc_z[k] * sig_295[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pb_x, pc_x, pc_y, shh0_416, shh0_417, \
                         shh0_419, shg_224, shh1_416, shh1_417, shh1_419, \
                         sig_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = pb_x[k] * shh0_416[k]
                   - f_8 * pc_x[k] * shh1_416[k];

        t_417[k] = pb_x[k] * shh0_417[k]
                   - f_8 * pc_x[k] * shh1_417[k];

        t_418[k] = f_9 * shg_224[k]
                   + f_3 * pc_y[k] * sig_299[k];

        t_419[k] = pb_x[k] * shh0_419[k]
                   - f_8 * pc_x[k] * shh1_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pb_x, pc_x, pc_y, pc_z, shh0_420, \
                         shh0_423, shg_210, shg_300, shg_303, shh1_420, shh1_423, \
                         sig_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pb_x[k] * shh0_420[k]
                   + f_12 * shg_300[k]
                   - f_8 * pc_x[k] * shh1_420[k];

        t_421[k] = f_3 * pc_y[k] * sig_300[k];

        t_422[k] = f_12 * shg_210[k]
                   + f_3 * pc_z[k] * sig_300[k];

        t_423[k] = pb_x[k] * shh0_423[k]
                   + f_11 * shg_303[k]
                   - f_8 * pc_x[k] * shh1_423[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_x, pc_x, pc_y, shh0_425, shh0_426, shg_305, \
                         shg_306, shh1_425, shh1_426, sig_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_y[k] * sig_302[k];

        t_425[k] = pb_x[k] * shh0_425[k]
                   + f_11 * shg_305[k]
                   - f_8 * pc_x[k] * shh1_425[k];

        t_426[k] = pb_x[k] * shh0_426[k]
                   + f_10 * shg_306[k]
                   - f_8 * pc_x[k] * shh1_426[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pb_x, pc_x, pc_y, pc_z, shh0_429, \
                         shg_213, shg_309, shg_310, shh1_429, sig_303, sig_305, \
                         sig_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_12 * shg_213[k]
                   + f_3 * pc_z[k] * sig_303[k];

        t_428[k] = f_3 * pc_y[k] * sig_305[k];

        t_429[k] = pb_x[k] * shh0_429[k]
                   + f_10 * shg_309[k]
                   - f_8 * pc_x[k] * shh1_429[k];

        t_430[k] = f_9 * shg_310[k]
                   + f_3 * pc_x[k] * sig_310[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, shg_311, shg_312, shg_313, shg_314, \
                         sig_311, sig_312, sig_313, sig_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_9 * shg_311[k]
                   + f_3 * pc_x[k] * sig_311[k];

        t_432[k] = f_9 * shg_312[k]
                   + f_3 * pc_x[k] * sig_312[k];

        t_433[k] = f_9 * shg_313[k]
                   + f_3 * pc_x[k] * sig_313[k];

        t_434[k] = f_9 * shg_314[k]
                   + f_3 * pc_x[k] * sig_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pb_x, pc_x, pc_z, shh0_435, shh0_437, \
                         shh0_438, shg_220, shh1_435, shh1_437, shh1_438, \
                         sig_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = pb_x[k] * shh0_435[k]
                   - f_8 * pc_x[k] * shh1_435[k];

        t_436[k] = f_12 * shg_220[k]
                   + f_3 * pc_z[k] * sig_310[k];

        t_437[k] = pb_x[k] * shh0_437[k]
                   - f_8 * pc_x[k] * shh1_437[k];

        t_438[k] = pb_x[k] * shh0_438[k]
                   - f_8 * pc_x[k] * shh1_438[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, pb_x, pc_x, pc_y, pc_z, shh0_440, \
                         shg_225, shh1_440, sif0_210, sif1_210, sig_314, \
                         sig_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_3 * pc_y[k] * sig_314[k];

        t_440[k] = pb_x[k] * shh0_440[k]
                   - f_8 * pc_x[k] * shh1_440[k];

        t_441[k] = f_1 * sif0_210[k]
                   - f_2 * sif1_210[k]
                   + f_3 * pc_x[k] * sig_315[k];

        t_442[k] = f_0 * shg_225[k]
                   + f_3 * pc_y[k] * sig_315[k];

        t_443[k] = f_3 * pc_z[k] * sig_315[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_x, pc_y, shg_227, sif0_213, sif0_215, \
                         sif1_213, sif1_215, sig_317, sig_318, \
                         sig_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_4 * sif0_213[k]
                   - f_5 * sif1_213[k]
                   + f_3 * pc_x[k] * sig_318[k];

        t_445[k] = f_0 * shg_227[k]
                   + f_3 * pc_y[k] * sig_317[k];

        t_446[k] = f_4 * sif0_215[k]
                   - f_5 * sif1_215[k]
                   + f_3 * pc_x[k] * sig_320[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pc_x, pc_y, pc_z, shg_230, sif0_216, \
                         sif0_219, sif1_216, sif1_219, sig_318, sig_320, sig_321, \
                         sig_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_6 * sif0_216[k]
                   - f_7 * sif1_216[k]
                   + f_3 * pc_x[k] * sig_321[k];

        t_448[k] = f_3 * pc_z[k] * sig_318[k];

        t_449[k] = f_0 * shg_230[k]
                   + f_3 * pc_y[k] * sig_320[k];

        t_450[k] = f_6 * sif0_219[k]
                   - f_7 * sif1_219[k]
                   + f_3 * pc_x[k] * sig_324[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, t_456, pc_x, pc_y, shg_235, \
                         sif0_216, sif1_216, sig_325, sig_326, sig_327, sig_328, \
                         sig_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_3 * pc_x[k] * sig_325[k];

        t_452[k] = f_3 * pc_x[k] * sig_326[k];

        t_453[k] = f_3 * pc_x[k] * sig_327[k];

        t_454[k] = f_3 * pc_x[k] * sig_328[k];

        t_455[k] = f_3 * pc_x[k] * sig_329[k];

        t_456[k] = f_0 * shg_235[k]
                   + f_1 * sif0_216[k]
                   - f_2 * sif1_216[k]
                   + f_3 * pc_y[k] * sig_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, pc_y, pc_z, shg_237, shg_238, sif0_218, \
                         sif0_219, sif1_218, sif1_219, sig_325, sig_327, \
                         sig_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * sig_325[k];

        t_458[k] = f_0 * shg_237[k]
                   + f_4 * sif0_218[k]
                   - f_5 * sif1_218[k]
                   + f_3 * pc_y[k] * sig_327[k];

        t_459[k] = f_0 * shg_238[k]
                   + f_6 * sif0_219[k]
                   - f_7 * sif1_219[k]
                   + f_3 * pc_y[k] * sig_328[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pb_z, pc_y, pc_z, shh0_315, shg_239, \
                         shg_240, shh1_315, sif0_219, sif1_219, sig_329, \
                         sig_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_0 * shg_239[k]
                   + f_3 * pc_y[k] * sig_329[k];

        t_461[k] = f_1 * sif0_219[k]
                   - f_2 * sif1_219[k]
                   + f_3 * pc_z[k] * sig_329[k];

        t_462[k] = pb_z[k] * shh0_315[k]
                   - f_8 * pc_z[k] * shh1_315[k];

        t_463[k] = f_12 * shg_240[k]
                   + f_3 * pc_y[k] * sig_330[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pb_z, pc_y, pc_z, shh0_318, shg_225, shg_242, \
                         shh1_318, sig_330, sig_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_9 * shg_225[k]
                   + f_3 * pc_z[k] * sig_330[k];

        t_465[k] = pb_z[k] * shh0_318[k]
                   - f_8 * pc_z[k] * shh1_318[k];

        t_466[k] = f_12 * shg_242[k]
                   + f_3 * pc_y[k] * sig_332[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, pb_z, pc_x, pc_y, pc_z, shh0_321, \
                         shg_228, shg_245, shh1_321, sif0_225, sif1_225, sig_333, \
                         sig_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_4 * sif0_225[k]
                   - f_5 * sif1_225[k]
                   + f_3 * pc_x[k] * sig_335[k];

        t_468[k] = pb_z[k] * shh0_321[k]
                   - f_8 * pc_z[k] * shh1_321[k];

        t_469[k] = f_9 * shg_228[k]
                   + f_3 * pc_z[k] * sig_333[k];

        t_470[k] = f_12 * shg_245[k]
                   + f_3 * pc_y[k] * sig_335[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, t_475, t_476, pc_x, sif0_229, sif1_229, \
                         sig_339, sig_340, sig_341, sig_342, sig_343, \
                         sig_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_6 * sif0_229[k]
                   - f_7 * sif1_229[k]
                   + f_3 * pc_x[k] * sig_339[k];

        t_472[k] = f_3 * pc_x[k] * sig_340[k];

        t_473[k] = f_3 * pc_x[k] * sig_341[k];

        t_474[k] = f_3 * pc_x[k] * sig_342[k];

        t_475[k] = f_3 * pc_x[k] * sig_343[k];

        t_476[k] = f_3 * pc_x[k] * sig_344[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pb_z, pc_z, shh0_330, shh0_332, shh0_333, \
                         shg_235, shg_236, shg_237, shh1_330, shh1_332, shh1_333, \
                         sig_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = pb_z[k] * shh0_330[k]
                   - f_8 * pc_z[k] * shh1_330[k];

        t_478[k] = f_9 * shg_235[k]
                   + f_3 * pc_z[k] * sig_340[k];

        t_479[k] = pb_z[k] * shh0_332[k]
                   + f_10 * shg_236[k]
                   - f_8 * pc_z[k] * shh1_332[k];

        t_480[k] = pb_z[k] * shh0_333[k]
                   + f_11 * shg_237[k]
                   - f_8 * pc_z[k] * shh1_333[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, pc_z, shg_239, shg_254, \
                         shg_255, sif0_229, sif0_230, sif1_229, sif1_230, sig_344, \
                         sig_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_12 * shg_254[k]
                   + f_3 * pc_y[k] * sig_344[k];

        t_482[k] = f_9 * shg_239[k]
                   + f_1 * sif0_229[k]
                   - f_2 * sif1_229[k]
                   + f_3 * pc_z[k] * sig_344[k];

        t_483[k] = f_1 * sif0_230[k]
                   - f_2 * sif1_230[k]
                   + f_3 * pc_x[k] * sig_345[k];

        t_484[k] = f_13 * shg_255[k]
                   + f_3 * pc_y[k] * sig_345[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_y, pc_z, shg_240, shg_257, sif0_233, \
                         sif1_233, sig_345, sig_347, sig_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_10 * shg_240[k]
                   + f_3 * pc_z[k] * sig_345[k];

        t_486[k] = f_4 * sif0_233[k]
                   - f_5 * sif1_233[k]
                   + f_3 * pc_x[k] * sig_348[k];

        t_487[k] = f_13 * shg_257[k]
                   + f_3 * pc_y[k] * sig_347[k];
    }
}

static auto
compute_prim_sih_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shh0,
                                                          const size_t shg, const size_t shh1,
                                                          const size_t sif0, const size_t sif1,
                                                          const size_t sig, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_13 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shh0_420 = buffer.data(shh0 + 420);
    const auto *shh0_425 = buffer.data(shh0 + 425);
    const auto *shh0_429 = buffer.data(shh0 + 429);
    const auto *shh0_435 = buffer.data(shh0 + 435);
    const auto *shh0_437 = buffer.data(shh0 + 437);
    const auto *shh0_438 = buffer.data(shh0 + 438);
    const auto *shh0_440 = buffer.data(shh0 + 440);

    const auto *shg_243 = buffer.data(shg + 243);
    const auto *shg_250 = buffer.data(shg + 250);
    const auto *shg_254 = buffer.data(shg + 254);
    const auto *shg_255 = buffer.data(shg + 255);
    const auto *shg_258 = buffer.data(shg + 258);
    const auto *shg_260 = buffer.data(shg + 260);
    const auto *shg_265 = buffer.data(shg + 265);
    const auto *shg_267 = buffer.data(shg + 267);
    const auto *shg_268 = buffer.data(shg + 268);
    const auto *shg_269 = buffer.data(shg + 269);
    const auto *shg_270 = buffer.data(shg + 270);
    const auto *shg_272 = buffer.data(shg + 272);
    const auto *shg_273 = buffer.data(shg + 273);
    const auto *shg_275 = buffer.data(shg + 275);
    const auto *shg_280 = buffer.data(shg + 280);
    const auto *shg_282 = buffer.data(shg + 282);
    const auto *shg_283 = buffer.data(shg + 283);
    const auto *shg_284 = buffer.data(shg + 284);
    const auto *shg_285 = buffer.data(shg + 285);
    const auto *shg_287 = buffer.data(shg + 287);
    const auto *shg_288 = buffer.data(shg + 288);
    const auto *shg_290 = buffer.data(shg + 290);
    const auto *shg_295 = buffer.data(shg + 295);
    const auto *shg_297 = buffer.data(shg + 297);
    const auto *shg_298 = buffer.data(shg + 298);
    const auto *shg_299 = buffer.data(shg + 299);
    const auto *shg_300 = buffer.data(shg + 300);
    const auto *shg_302 = buffer.data(shg + 302);
    const auto *shg_303 = buffer.data(shg + 303);
    const auto *shg_305 = buffer.data(shg + 305);
    const auto *shg_310 = buffer.data(shg + 310);
    const auto *shg_312 = buffer.data(shg + 312);
    const auto *shg_313 = buffer.data(shg + 313);
    const auto *shg_314 = buffer.data(shg + 314);

    const auto *shh1_420 = buffer.data(shh1 + 420);
    const auto *shh1_425 = buffer.data(shh1 + 425);
    const auto *shh1_429 = buffer.data(shh1 + 429);
    const auto *shh1_435 = buffer.data(shh1 + 435);
    const auto *shh1_437 = buffer.data(shh1 + 437);
    const auto *shh1_438 = buffer.data(shh1 + 438);
    const auto *shh1_440 = buffer.data(shh1 + 440);

    const auto *sif0_235 = buffer.data(sif0 + 235);
    const auto *sif0_236 = buffer.data(sif0 + 236);
    const auto *sif0_238 = buffer.data(sif0 + 238);
    const auto *sif0_239 = buffer.data(sif0 + 239);
    const auto *sif0_240 = buffer.data(sif0 + 240);
    const auto *sif0_243 = buffer.data(sif0 + 243);
    const auto *sif0_245 = buffer.data(sif0 + 245);
    const auto *sif0_246 = buffer.data(sif0 + 246);
    const auto *sif0_248 = buffer.data(sif0 + 248);
    const auto *sif0_249 = buffer.data(sif0 + 249);
    const auto *sif0_250 = buffer.data(sif0 + 250);
    const auto *sif0_253 = buffer.data(sif0 + 253);
    const auto *sif0_255 = buffer.data(sif0 + 255);
    const auto *sif0_256 = buffer.data(sif0 + 256);
    const auto *sif0_258 = buffer.data(sif0 + 258);
    const auto *sif0_259 = buffer.data(sif0 + 259);
    const auto *sif0_263 = buffer.data(sif0 + 263);
    const auto *sif0_266 = buffer.data(sif0 + 266);
    const auto *sif0_270 = buffer.data(sif0 + 270);
    const auto *sif0_273 = buffer.data(sif0 + 273);
    const auto *sif0_275 = buffer.data(sif0 + 275);
    const auto *sif0_276 = buffer.data(sif0 + 276);
    const auto *sif0_278 = buffer.data(sif0 + 278);
    const auto *sif0_279 = buffer.data(sif0 + 279);

    const auto *sif1_235 = buffer.data(sif1 + 235);
    const auto *sif1_236 = buffer.data(sif1 + 236);
    const auto *sif1_238 = buffer.data(sif1 + 238);
    const auto *sif1_239 = buffer.data(sif1 + 239);
    const auto *sif1_240 = buffer.data(sif1 + 240);
    const auto *sif1_243 = buffer.data(sif1 + 243);
    const auto *sif1_245 = buffer.data(sif1 + 245);
    const auto *sif1_246 = buffer.data(sif1 + 246);
    const auto *sif1_248 = buffer.data(sif1 + 248);
    const auto *sif1_249 = buffer.data(sif1 + 249);
    const auto *sif1_250 = buffer.data(sif1 + 250);
    const auto *sif1_253 = buffer.data(sif1 + 253);
    const auto *sif1_255 = buffer.data(sif1 + 255);
    const auto *sif1_256 = buffer.data(sif1 + 256);
    const auto *sif1_258 = buffer.data(sif1 + 258);
    const auto *sif1_259 = buffer.data(sif1 + 259);
    const auto *sif1_263 = buffer.data(sif1 + 263);
    const auto *sif1_266 = buffer.data(sif1 + 266);
    const auto *sif1_270 = buffer.data(sif1 + 270);
    const auto *sif1_273 = buffer.data(sif1 + 273);
    const auto *sif1_275 = buffer.data(sif1 + 275);
    const auto *sif1_276 = buffer.data(sif1 + 276);
    const auto *sif1_278 = buffer.data(sif1 + 278);
    const auto *sif1_279 = buffer.data(sif1 + 279);

    const auto *sig_348 = buffer.data(sig + 348);
    const auto *sig_350 = buffer.data(sig + 350);
    const auto *sig_351 = buffer.data(sig + 351);
    const auto *sig_354 = buffer.data(sig + 354);
    const auto *sig_355 = buffer.data(sig + 355);
    const auto *sig_356 = buffer.data(sig + 356);
    const auto *sig_357 = buffer.data(sig + 357);
    const auto *sig_358 = buffer.data(sig + 358);
    const auto *sig_359 = buffer.data(sig + 359);
    const auto *sig_360 = buffer.data(sig + 360);
    const auto *sig_362 = buffer.data(sig + 362);
    const auto *sig_363 = buffer.data(sig + 363);
    const auto *sig_365 = buffer.data(sig + 365);
    const auto *sig_366 = buffer.data(sig + 366);
    const auto *sig_369 = buffer.data(sig + 369);
    const auto *sig_370 = buffer.data(sig + 370);
    const auto *sig_371 = buffer.data(sig + 371);
    const auto *sig_372 = buffer.data(sig + 372);
    const auto *sig_373 = buffer.data(sig + 373);
    const auto *sig_374 = buffer.data(sig + 374);
    const auto *sig_375 = buffer.data(sig + 375);
    const auto *sig_377 = buffer.data(sig + 377);
    const auto *sig_378 = buffer.data(sig + 378);
    const auto *sig_380 = buffer.data(sig + 380);
    const auto *sig_381 = buffer.data(sig + 381);
    const auto *sig_384 = buffer.data(sig + 384);
    const auto *sig_385 = buffer.data(sig + 385);
    const auto *sig_386 = buffer.data(sig + 386);
    const auto *sig_387 = buffer.data(sig + 387);
    const auto *sig_388 = buffer.data(sig + 388);
    const auto *sig_389 = buffer.data(sig + 389);
    const auto *sig_390 = buffer.data(sig + 390);
    const auto *sig_392 = buffer.data(sig + 392);
    const auto *sig_393 = buffer.data(sig + 393);
    const auto *sig_395 = buffer.data(sig + 395);
    const auto *sig_396 = buffer.data(sig + 396);
    const auto *sig_400 = buffer.data(sig + 400);
    const auto *sig_401 = buffer.data(sig + 401);
    const auto *sig_402 = buffer.data(sig + 402);
    const auto *sig_403 = buffer.data(sig + 403);
    const auto *sig_404 = buffer.data(sig + 404);
    const auto *sig_405 = buffer.data(sig + 405);
    const auto *sig_407 = buffer.data(sig + 407);
    const auto *sig_408 = buffer.data(sig + 408);
    const auto *sig_410 = buffer.data(sig + 410);
    const auto *sig_411 = buffer.data(sig + 411);
    const auto *sig_414 = buffer.data(sig + 414);
    const auto *sig_415 = buffer.data(sig + 415);
    const auto *sig_416 = buffer.data(sig + 416);
    const auto *sig_417 = buffer.data(sig + 417);
    const auto *sig_418 = buffer.data(sig + 418);
    const auto *sig_419 = buffer.data(sig + 419);

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pc_x, pc_y, pc_z, shg_243, shg_260, \
                         sif0_235, sif0_236, sif1_235, sif1_236, sig_348, sig_350, \
                         sig_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_4 * sif0_235[k]
                   - f_5 * sif1_235[k]
                   + f_3 * pc_x[k] * sig_350[k];

        t_489[k] = f_6 * sif0_236[k]
                   - f_7 * sif1_236[k]
                   + f_3 * pc_x[k] * sig_351[k];

        t_490[k] = f_10 * shg_243[k]
                   + f_3 * pc_z[k] * sig_348[k];

        t_491[k] = f_13 * shg_260[k]
                   + f_3 * pc_y[k] * sig_350[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, t_497, pc_x, sif0_239, sif1_239, \
                         sig_354, sig_355, sig_356, sig_357, sig_358, \
                         sig_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_6 * sif0_239[k]
                   - f_7 * sif1_239[k]
                   + f_3 * pc_x[k] * sig_354[k];

        t_493[k] = f_3 * pc_x[k] * sig_355[k];

        t_494[k] = f_3 * pc_x[k] * sig_356[k];

        t_495[k] = f_3 * pc_x[k] * sig_357[k];

        t_496[k] = f_3 * pc_x[k] * sig_358[k];

        t_497[k] = f_3 * pc_x[k] * sig_359[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, pc_z, shg_250, shg_265, shg_267, sif0_236, \
                         sif0_238, sif1_236, sif1_238, sig_355, \
                         sig_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * shg_265[k]
                   + f_1 * sif0_236[k]
                   - f_2 * sif1_236[k]
                   + f_3 * pc_y[k] * sig_355[k];

        t_499[k] = f_10 * shg_250[k]
                   + f_3 * pc_z[k] * sig_355[k];

        t_500[k] = f_13 * shg_267[k]
                   + f_4 * sif0_238[k]
                   - f_5 * sif1_238[k]
                   + f_3 * pc_y[k] * sig_357[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pc_y, pc_z, shg_254, shg_268, shg_269, sif0_239, \
                         sif1_239, sig_358, sig_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * shg_268[k]
                   + f_6 * sif0_239[k]
                   - f_7 * sif1_239[k]
                   + f_3 * pc_y[k] * sig_358[k];

        t_502[k] = f_13 * shg_269[k]
                   + f_3 * pc_y[k] * sig_359[k];

        t_503[k] = f_10 * shg_254[k]
                   + f_1 * sif0_239[k]
                   - f_2 * sif1_239[k]
                   + f_3 * pc_z[k] * sig_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pc_x, pc_y, pc_z, shg_255, shg_270, \
                         sif0_240, sif0_243, sif1_240, sif1_243, sig_360, \
                         sig_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_1 * sif0_240[k]
                   - f_2 * sif1_240[k]
                   + f_3 * pc_x[k] * sig_360[k];

        t_505[k] = f_11 * shg_270[k]
                   + f_3 * pc_y[k] * sig_360[k];

        t_506[k] = f_11 * shg_255[k]
                   + f_3 * pc_z[k] * sig_360[k];

        t_507[k] = f_4 * sif0_243[k]
                   - f_5 * sif1_243[k]
                   + f_3 * pc_x[k] * sig_363[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pc_x, pc_y, shg_272, sif0_245, sif0_246, \
                         sif1_245, sif1_246, sig_362, sig_365, \
                         sig_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_11 * shg_272[k]
                   + f_3 * pc_y[k] * sig_362[k];

        t_509[k] = f_4 * sif0_245[k]
                   - f_5 * sif1_245[k]
                   + f_3 * pc_x[k] * sig_365[k];

        t_510[k] = f_6 * sif0_246[k]
                   - f_7 * sif1_246[k]
                   + f_3 * pc_x[k] * sig_366[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pc_x, pc_y, pc_z, shg_258, shg_275, \
                         sif0_249, sif1_249, sig_363, sig_365, sig_369, \
                         sig_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_11 * shg_258[k]
                   + f_3 * pc_z[k] * sig_363[k];

        t_512[k] = f_11 * shg_275[k]
                   + f_3 * pc_y[k] * sig_365[k];

        t_513[k] = f_6 * sif0_249[k]
                   - f_7 * sif1_249[k]
                   + f_3 * pc_x[k] * sig_369[k];

        t_514[k] = f_3 * pc_x[k] * sig_370[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, pc_x, pc_y, shg_280, sif0_246, \
                         sif1_246, sig_370, sig_371, sig_372, sig_373, \
                         sig_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_3 * pc_x[k] * sig_371[k];

        t_516[k] = f_3 * pc_x[k] * sig_372[k];

        t_517[k] = f_3 * pc_x[k] * sig_373[k];

        t_518[k] = f_3 * pc_x[k] * sig_374[k];

        t_519[k] = f_11 * shg_280[k]
                   + f_1 * sif0_246[k]
                   - f_2 * sif1_246[k]
                   + f_3 * pc_y[k] * sig_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, shg_265, shg_282, shg_283, sif0_248, \
                         sif0_249, sif1_248, sif1_249, sig_370, sig_372, \
                         sig_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * shg_265[k]
                   + f_3 * pc_z[k] * sig_370[k];

        t_521[k] = f_11 * shg_282[k]
                   + f_4 * sif0_248[k]
                   - f_5 * sif1_248[k]
                   + f_3 * pc_y[k] * sig_372[k];

        t_522[k] = f_11 * shg_283[k]
                   + f_6 * sif0_249[k]
                   - f_7 * sif1_249[k]
                   + f_3 * pc_y[k] * sig_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pc_x, pc_y, pc_z, shg_269, shg_284, \
                         shg_285, sif0_249, sif0_250, sif1_249, sif1_250, sig_374, \
                         sig_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * shg_284[k]
                   + f_3 * pc_y[k] * sig_374[k];

        t_524[k] = f_11 * shg_269[k]
                   + f_1 * sif0_249[k]
                   - f_2 * sif1_249[k]
                   + f_3 * pc_z[k] * sig_374[k];

        t_525[k] = f_1 * sif0_250[k]
                   - f_2 * sif1_250[k]
                   + f_3 * pc_x[k] * sig_375[k];

        t_526[k] = f_10 * shg_285[k]
                   + f_3 * pc_y[k] * sig_375[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pc_x, pc_y, pc_z, shg_270, shg_287, sif0_253, \
                         sif1_253, sig_375, sig_377, sig_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_13 * shg_270[k]
                   + f_3 * pc_z[k] * sig_375[k];

        t_528[k] = f_4 * sif0_253[k]
                   - f_5 * sif1_253[k]
                   + f_3 * pc_x[k] * sig_378[k];

        t_529[k] = f_10 * shg_287[k]
                   + f_3 * pc_y[k] * sig_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pc_x, pc_y, pc_z, shg_273, shg_290, \
                         sif0_255, sif0_256, sif1_255, sif1_256, sig_378, sig_380, \
                         sig_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_4 * sif0_255[k]
                   - f_5 * sif1_255[k]
                   + f_3 * pc_x[k] * sig_380[k];

        t_531[k] = f_6 * sif0_256[k]
                   - f_7 * sif1_256[k]
                   + f_3 * pc_x[k] * sig_381[k];

        t_532[k] = f_13 * shg_273[k]
                   + f_3 * pc_z[k] * sig_378[k];

        t_533[k] = f_10 * shg_290[k]
                   + f_3 * pc_y[k] * sig_380[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, t_539, pc_x, sif0_259, sif1_259, \
                         sig_384, sig_385, sig_386, sig_387, sig_388, \
                         sig_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_6 * sif0_259[k]
                   - f_7 * sif1_259[k]
                   + f_3 * pc_x[k] * sig_384[k];

        t_535[k] = f_3 * pc_x[k] * sig_385[k];

        t_536[k] = f_3 * pc_x[k] * sig_386[k];

        t_537[k] = f_3 * pc_x[k] * sig_387[k];

        t_538[k] = f_3 * pc_x[k] * sig_388[k];

        t_539[k] = f_3 * pc_x[k] * sig_389[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pc_y, pc_z, shg_280, shg_295, shg_297, sif0_256, \
                         sif0_258, sif1_256, sif1_258, sig_385, \
                         sig_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_10 * shg_295[k]
                   + f_1 * sif0_256[k]
                   - f_2 * sif1_256[k]
                   + f_3 * pc_y[k] * sig_385[k];

        t_541[k] = f_13 * shg_280[k]
                   + f_3 * pc_z[k] * sig_385[k];

        t_542[k] = f_10 * shg_297[k]
                   + f_4 * sif0_258[k]
                   - f_5 * sif1_258[k]
                   + f_3 * pc_y[k] * sig_387[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, pb_y, pc_y, pc_z, shh0_420, shg_284, \
                         shg_298, shg_299, shh1_420, sif0_259, sif1_259, sig_388, \
                         sig_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_10 * shg_298[k]
                   + f_6 * sif0_259[k]
                   - f_7 * sif1_259[k]
                   + f_3 * pc_y[k] * sig_388[k];

        t_544[k] = f_10 * shg_299[k]
                   + f_3 * pc_y[k] * sig_389[k];

        t_545[k] = f_13 * shg_284[k]
                   + f_1 * sif0_259[k]
                   - f_2 * sif1_259[k]
                   + f_3 * pc_z[k] * sig_389[k];

        t_546[k] = pb_y[k] * shh0_420[k]
                   - f_8 * pc_y[k] * shh1_420[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, pc_x, pc_y, pc_z, shg_285, shg_300, \
                         shg_302, sif0_263, sif1_263, sig_390, sig_392, \
                         sig_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_9 * shg_300[k]
                   + f_3 * pc_y[k] * sig_390[k];

        t_548[k] = f_12 * shg_285[k]
                   + f_3 * pc_z[k] * sig_390[k];

        t_549[k] = f_4 * sif0_263[k]
                   - f_5 * sif1_263[k]
                   + f_3 * pc_x[k] * sig_393[k];

        t_550[k] = f_9 * shg_302[k]
                   + f_3 * pc_y[k] * sig_392[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pb_y, pc_x, pc_y, pc_z, shh0_425, shg_288, \
                         shh1_425, sif0_266, sif1_266, sig_393, \
                         sig_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = pb_y[k] * shh0_425[k]
                   - f_8 * pc_y[k] * shh1_425[k];

        t_552[k] = f_6 * sif0_266[k]
                   - f_7 * sif1_266[k]
                   + f_3 * pc_x[k] * sig_396[k];

        t_553[k] = f_12 * shg_288[k]
                   + f_3 * pc_z[k] * sig_393[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, t_558, pb_y, pc_x, pc_y, shh0_429, \
                         shg_305, shh1_429, sig_395, sig_400, sig_401, \
                         sig_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_9 * shg_305[k]
                   + f_3 * pc_y[k] * sig_395[k];

        t_555[k] = pb_y[k] * shh0_429[k]
                   - f_8 * pc_y[k] * shh1_429[k];

        t_556[k] = f_3 * pc_x[k] * sig_400[k];

        t_557[k] = f_3 * pc_x[k] * sig_401[k];

        t_558[k] = f_3 * pc_x[k] * sig_402[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pb_y, pc_x, pc_y, pc_z, shh0_435, \
                         shg_295, shg_310, shh1_435, sig_400, sig_403, \
                         sig_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_3 * pc_x[k] * sig_403[k];

        t_560[k] = f_3 * pc_x[k] * sig_404[k];

        t_561[k] = pb_y[k] * shh0_435[k]
                   + f_12 * shg_310[k]
                   - f_8 * pc_y[k] * shh1_435[k];

        t_562[k] = f_12 * shg_295[k]
                   + f_3 * pc_z[k] * sig_400[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, pb_y, pc_y, shh0_437, shh0_438, shh0_440, \
                         shg_312, shg_313, shg_314, shh1_437, shh1_438, shh1_440, \
                         sig_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = pb_y[k] * shh0_437[k]
                   + f_11 * shg_312[k]
                   - f_8 * pc_y[k] * shh1_437[k];

        t_564[k] = pb_y[k] * shh0_438[k]
                   + f_10 * shg_313[k]
                   - f_8 * pc_y[k] * shh1_438[k];

        t_565[k] = f_9 * shg_314[k]
                   + f_3 * pc_y[k] * sig_404[k];

        t_566[k] = pb_y[k] * shh0_440[k]
                   - f_8 * pc_y[k] * shh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, pc_x, pc_y, pc_z, shg_300, \
                         sif0_270, sif0_273, sif1_270, sif1_273, sig_405, sig_407, \
                         sig_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_1 * sif0_270[k]
                   - f_2 * sif1_270[k]
                   + f_3 * pc_x[k] * sig_405[k];

        t_568[k] = f_3 * pc_y[k] * sig_405[k];

        t_569[k] = f_0 * shg_300[k]
                   + f_3 * pc_z[k] * sig_405[k];

        t_570[k] = f_4 * sif0_273[k]
                   - f_5 * sif1_273[k]
                   + f_3 * pc_x[k] * sig_408[k];

        t_571[k] = f_3 * pc_y[k] * sig_407[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, pc_z, shg_303, sif0_275, \
                         sif0_276, sif1_275, sif1_276, sig_408, sig_410, \
                         sig_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_4 * sif0_275[k]
                   - f_5 * sif1_275[k]
                   + f_3 * pc_x[k] * sig_410[k];

        t_573[k] = f_6 * sif0_276[k]
                   - f_7 * sif1_276[k]
                   + f_3 * pc_x[k] * sig_411[k];

        t_574[k] = f_0 * shg_303[k]
                   + f_3 * pc_z[k] * sig_408[k];

        t_575[k] = f_3 * pc_y[k] * sig_410[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, t_581, pc_x, sif0_279, sif1_279, \
                         sig_414, sig_415, sig_416, sig_417, sig_418, \
                         sig_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_6 * sif0_279[k]
                   - f_7 * sif1_279[k]
                   + f_3 * pc_x[k] * sig_414[k];

        t_577[k] = f_3 * pc_x[k] * sig_415[k];

        t_578[k] = f_3 * pc_x[k] * sig_416[k];

        t_579[k] = f_3 * pc_x[k] * sig_417[k];

        t_580[k] = f_3 * pc_x[k] * sig_418[k];

        t_581[k] = f_3 * pc_x[k] * sig_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pc_y, pc_z, shg_310, sif0_276, sif0_278, \
                         sif0_279, sif1_276, sif1_278, sif1_279, sig_415, sig_417, \
                         sig_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_1 * sif0_276[k]
                   - f_2 * sif1_276[k]
                   + f_3 * pc_y[k] * sig_415[k];

        t_583[k] = f_0 * shg_310[k]
                   + f_3 * pc_z[k] * sig_415[k];

        t_584[k] = f_4 * sif0_278[k]
                   - f_5 * sif1_278[k]
                   + f_3 * pc_y[k] * sig_417[k];

        t_585[k] = f_6 * sif0_279[k]
                   - f_7 * sif1_279[k]
                   + f_3 * pc_y[k] * sig_418[k];
    }

#pragma omp simd aligned(t_586, t_587, pc_y, pc_z, shg_314, sif0_279, sif1_279, \
                         sig_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_3 * pc_y[k] * sig_419[k];

        t_587[k] = f_0 * shg_314[k]
                   + f_1 * sif0_279[k]
                   - f_2 * sif1_279[k]
                   + f_3 * pc_z[k] * sig_419[k];
    }
}

auto
compute_prim_sih_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t shh0, const size_t shg,
                                                   const size_t shh1, const size_t sif0,
                                                   const size_t sif1, const size_t sig,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sih_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, shh0, shg,
                                                              shh1, sif0, sif1, sig, ncols,
                                                              gamma, p, q);

    compute_prim_sih_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, shh0, shg,
                                                              shh1, sif0, sif1, sig, ncols,
                                                              gamma, p, q);

    compute_prim_sih_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, shh0, shg,
                                                              shh1, sif0, sif1, sig, ncols,
                                                              gamma, p, q);

    compute_prim_sih_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, shh0, shg,
                                                              shh1, sif0, sif1, sig, ncols,
                                                              gamma, p, q);

    compute_prim_sih_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, shh0, shg,
                                                              shh1, sif0, sif1, sig, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
