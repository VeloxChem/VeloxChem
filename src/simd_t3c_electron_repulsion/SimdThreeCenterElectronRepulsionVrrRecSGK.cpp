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


#include "SimdThreeCenterElectronRepulsionVrrRecSGK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sgk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfk0,
                                                          const size_t sfi, const size_t sfk1,
                                                          const size_t sgh0, const size_t sgh1,
                                                          const size_t sgi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfk0_0 = buffer.data(sfk0 + 0);
    const auto *sfk0_3 = buffer.data(sfk0 + 3);
    const auto *sfk0_5 = buffer.data(sfk0 + 5);
    const auto *sfk0_6 = buffer.data(sfk0 + 6);
    const auto *sfk0_9 = buffer.data(sfk0 + 9);
    const auto *sfk0_10 = buffer.data(sfk0 + 10);
    const auto *sfk0_12 = buffer.data(sfk0 + 12);
    const auto *sfk0_14 = buffer.data(sfk0 + 14);
    const auto *sfk0_15 = buffer.data(sfk0 + 15);
    const auto *sfk0_17 = buffer.data(sfk0 + 17);
    const auto *sfk0_18 = buffer.data(sfk0 + 18);
    const auto *sfk0_20 = buffer.data(sfk0 + 20);
    const auto *sfk0_28 = buffer.data(sfk0 + 28);
    const auto *sfk0_35 = buffer.data(sfk0 + 35);

    const auto *sfi_0 = buffer.data(sfi + 0);
    const auto *sfi_1 = buffer.data(sfi + 1);
    const auto *sfi_2 = buffer.data(sfi + 2);
    const auto *sfi_3 = buffer.data(sfi + 3);
    const auto *sfi_5 = buffer.data(sfi + 5);
    const auto *sfi_6 = buffer.data(sfi + 6);
    const auto *sfi_7 = buffer.data(sfi + 7);
    const auto *sfi_8 = buffer.data(sfi + 8);
    const auto *sfi_9 = buffer.data(sfi + 9);
    const auto *sfi_10 = buffer.data(sfi + 10);
    const auto *sfi_11 = buffer.data(sfi + 11);
    const auto *sfi_12 = buffer.data(sfi + 12);
    const auto *sfi_13 = buffer.data(sfi + 13);
    const auto *sfi_14 = buffer.data(sfi + 14);
    const auto *sfi_15 = buffer.data(sfi + 15);
    const auto *sfi_17 = buffer.data(sfi + 17);
    const auto *sfi_18 = buffer.data(sfi + 18);
    const auto *sfi_20 = buffer.data(sfi + 20);
    const auto *sfi_21 = buffer.data(sfi + 21);
    const auto *sfi_22 = buffer.data(sfi + 22);
    const auto *sfi_23 = buffer.data(sfi + 23);
    const auto *sfi_24 = buffer.data(sfi + 24);
    const auto *sfi_25 = buffer.data(sfi + 25);
    const auto *sfi_26 = buffer.data(sfi + 26);
    const auto *sfi_27 = buffer.data(sfi + 27);
    const auto *sfi_28 = buffer.data(sfi + 28);
    const auto *sfi_30 = buffer.data(sfi + 30);
    const auto *sfi_33 = buffer.data(sfi + 33);
    const auto *sfi_37 = buffer.data(sfi + 37);
    const auto *sfi_49 = buffer.data(sfi + 49);
    const auto *sfi_50 = buffer.data(sfi + 50);
    const auto *sfi_51 = buffer.data(sfi + 51);
    const auto *sfi_52 = buffer.data(sfi + 52);
    const auto *sfi_53 = buffer.data(sfi + 53);
    const auto *sfi_54 = buffer.data(sfi + 54);
    const auto *sfi_55 = buffer.data(sfi + 55);
    const auto *sfi_77 = buffer.data(sfi + 77);
    const auto *sfi_78 = buffer.data(sfi + 78);
    const auto *sfi_79 = buffer.data(sfi + 79);
    const auto *sfi_80 = buffer.data(sfi + 80);
    const auto *sfi_81 = buffer.data(sfi + 81);
    const auto *sfi_82 = buffer.data(sfi + 82);
    const auto *sfi_83 = buffer.data(sfi + 83);
    const auto *sfi_84 = buffer.data(sfi + 84);
    const auto *sfi_87 = buffer.data(sfi + 87);
    const auto *sfi_89 = buffer.data(sfi + 89);
    const auto *sfi_90 = buffer.data(sfi + 90);
    const auto *sfi_93 = buffer.data(sfi + 93);
    const auto *sfi_94 = buffer.data(sfi + 94);
    const auto *sfi_96 = buffer.data(sfi + 96);

    const auto *sfk1_0 = buffer.data(sfk1 + 0);
    const auto *sfk1_3 = buffer.data(sfk1 + 3);
    const auto *sfk1_5 = buffer.data(sfk1 + 5);
    const auto *sfk1_6 = buffer.data(sfk1 + 6);
    const auto *sfk1_9 = buffer.data(sfk1 + 9);
    const auto *sfk1_10 = buffer.data(sfk1 + 10);
    const auto *sfk1_12 = buffer.data(sfk1 + 12);
    const auto *sfk1_14 = buffer.data(sfk1 + 14);
    const auto *sfk1_15 = buffer.data(sfk1 + 15);
    const auto *sfk1_17 = buffer.data(sfk1 + 17);
    const auto *sfk1_18 = buffer.data(sfk1 + 18);
    const auto *sfk1_20 = buffer.data(sfk1 + 20);
    const auto *sfk1_28 = buffer.data(sfk1 + 28);
    const auto *sfk1_35 = buffer.data(sfk1 + 35);

    const auto *sgh0_0 = buffer.data(sgh0 + 0);
    const auto *sgh0_3 = buffer.data(sgh0 + 3);
    const auto *sgh0_5 = buffer.data(sgh0 + 5);
    const auto *sgh0_6 = buffer.data(sgh0 + 6);
    const auto *sgh0_9 = buffer.data(sgh0 + 9);
    const auto *sgh0_10 = buffer.data(sgh0 + 10);
    const auto *sgh0_12 = buffer.data(sgh0 + 12);
    const auto *sgh0_14 = buffer.data(sgh0 + 14);
    const auto *sgh0_15 = buffer.data(sgh0 + 15);
    const auto *sgh0_17 = buffer.data(sgh0 + 17);
    const auto *sgh0_18 = buffer.data(sgh0 + 18);
    const auto *sgh0_19 = buffer.data(sgh0 + 19);
    const auto *sgh0_20 = buffer.data(sgh0 + 20);
    const auto *sgh0_36 = buffer.data(sgh0 + 36);
    const auto *sgh0_38 = buffer.data(sgh0 + 38);
    const auto *sgh0_39 = buffer.data(sgh0 + 39);
    const auto *sgh0_40 = buffer.data(sgh0 + 40);
    const auto *sgh0_41 = buffer.data(sgh0 + 41);
    const auto *sgh0_59 = buffer.data(sgh0 + 59);
    const auto *sgh0_60 = buffer.data(sgh0 + 60);
    const auto *sgh0_61 = buffer.data(sgh0 + 61);
    const auto *sgh0_62 = buffer.data(sgh0 + 62);
    const auto *sgh0_63 = buffer.data(sgh0 + 63);
    const auto *sgh0_66 = buffer.data(sgh0 + 66);
    const auto *sgh0_68 = buffer.data(sgh0 + 68);
    const auto *sgh0_69 = buffer.data(sgh0 + 69);
    const auto *sgh0_72 = buffer.data(sgh0 + 72);
    const auto *sgh0_73 = buffer.data(sgh0 + 73);
    const auto *sgh0_75 = buffer.data(sgh0 + 75);

    const auto *sgh1_0 = buffer.data(sgh1 + 0);
    const auto *sgh1_3 = buffer.data(sgh1 + 3);
    const auto *sgh1_5 = buffer.data(sgh1 + 5);
    const auto *sgh1_6 = buffer.data(sgh1 + 6);
    const auto *sgh1_9 = buffer.data(sgh1 + 9);
    const auto *sgh1_10 = buffer.data(sgh1 + 10);
    const auto *sgh1_12 = buffer.data(sgh1 + 12);
    const auto *sgh1_14 = buffer.data(sgh1 + 14);
    const auto *sgh1_15 = buffer.data(sgh1 + 15);
    const auto *sgh1_17 = buffer.data(sgh1 + 17);
    const auto *sgh1_18 = buffer.data(sgh1 + 18);
    const auto *sgh1_19 = buffer.data(sgh1 + 19);
    const auto *sgh1_20 = buffer.data(sgh1 + 20);
    const auto *sgh1_36 = buffer.data(sgh1 + 36);
    const auto *sgh1_38 = buffer.data(sgh1 + 38);
    const auto *sgh1_39 = buffer.data(sgh1 + 39);
    const auto *sgh1_40 = buffer.data(sgh1 + 40);
    const auto *sgh1_41 = buffer.data(sgh1 + 41);
    const auto *sgh1_59 = buffer.data(sgh1 + 59);
    const auto *sgh1_60 = buffer.data(sgh1 + 60);
    const auto *sgh1_61 = buffer.data(sgh1 + 61);
    const auto *sgh1_62 = buffer.data(sgh1 + 62);
    const auto *sgh1_63 = buffer.data(sgh1 + 63);
    const auto *sgh1_66 = buffer.data(sgh1 + 66);
    const auto *sgh1_68 = buffer.data(sgh1 + 68);
    const auto *sgh1_69 = buffer.data(sgh1 + 69);
    const auto *sgh1_72 = buffer.data(sgh1 + 72);
    const auto *sgh1_73 = buffer.data(sgh1 + 73);
    const auto *sgh1_75 = buffer.data(sgh1 + 75);

    const auto *sgi_0 = buffer.data(sgi + 0);
    const auto *sgi_2 = buffer.data(sgi + 2);
    const auto *sgi_3 = buffer.data(sgi + 3);
    const auto *sgi_5 = buffer.data(sgi + 5);
    const auto *sgi_6 = buffer.data(sgi + 6);
    const auto *sgi_9 = buffer.data(sgi + 9);
    const auto *sgi_10 = buffer.data(sgi + 10);
    const auto *sgi_12 = buffer.data(sgi + 12);
    const auto *sgi_14 = buffer.data(sgi + 14);
    const auto *sgi_15 = buffer.data(sgi + 15);
    const auto *sgi_17 = buffer.data(sgi + 17);
    const auto *sgi_18 = buffer.data(sgi + 18);
    const auto *sgi_20 = buffer.data(sgi + 20);
    const auto *sgi_21 = buffer.data(sgi + 21);
    const auto *sgi_22 = buffer.data(sgi + 22);
    const auto *sgi_23 = buffer.data(sgi + 23);
    const auto *sgi_24 = buffer.data(sgi + 24);
    const auto *sgi_25 = buffer.data(sgi + 25);
    const auto *sgi_26 = buffer.data(sgi + 26);
    const auto *sgi_27 = buffer.data(sgi + 27);
    const auto *sgi_28 = buffer.data(sgi + 28);
    const auto *sgi_30 = buffer.data(sgi + 30);
    const auto *sgi_31 = buffer.data(sgi + 31);
    const auto *sgi_33 = buffer.data(sgi + 33);
    const auto *sgi_34 = buffer.data(sgi + 34);
    const auto *sgi_37 = buffer.data(sgi + 37);
    const auto *sgi_38 = buffer.data(sgi + 38);
    const auto *sgi_42 = buffer.data(sgi + 42);
    const auto *sgi_49 = buffer.data(sgi + 49);
    const auto *sgi_50 = buffer.data(sgi + 50);
    const auto *sgi_51 = buffer.data(sgi + 51);
    const auto *sgi_52 = buffer.data(sgi + 52);
    const auto *sgi_53 = buffer.data(sgi + 53);
    const auto *sgi_54 = buffer.data(sgi + 54);
    const auto *sgi_55 = buffer.data(sgi + 55);
    const auto *sgi_56 = buffer.data(sgi + 56);
    const auto *sgi_58 = buffer.data(sgi + 58);
    const auto *sgi_59 = buffer.data(sgi + 59);
    const auto *sgi_61 = buffer.data(sgi + 61);
    const auto *sgi_62 = buffer.data(sgi + 62);
    const auto *sgi_65 = buffer.data(sgi + 65);
    const auto *sgi_66 = buffer.data(sgi + 66);
    const auto *sgi_70 = buffer.data(sgi + 70);
    const auto *sgi_77 = buffer.data(sgi + 77);
    const auto *sgi_78 = buffer.data(sgi + 78);
    const auto *sgi_79 = buffer.data(sgi + 79);
    const auto *sgi_80 = buffer.data(sgi + 80);
    const auto *sgi_81 = buffer.data(sgi + 81);
    const auto *sgi_82 = buffer.data(sgi + 82);
    const auto *sgi_83 = buffer.data(sgi + 83);
    const auto *sgi_84 = buffer.data(sgi + 84);
    const auto *sgi_86 = buffer.data(sgi + 86);
    const auto *sgi_87 = buffer.data(sgi + 87);
    const auto *sgi_89 = buffer.data(sgi + 89);
    const auto *sgi_90 = buffer.data(sgi + 90);
    const auto *sgi_93 = buffer.data(sgi + 93);
    const auto *sgi_94 = buffer.data(sgi + 94);
    const auto *sgi_96 = buffer.data(sgi + 96);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sfi_0, sfi_3, sgh0_0, sgh0_3, \
                         sgh1_0, sgh1_3, sgi_0, sgi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfi_0[k]
                 + f_1 * sgh0_0[k]
                 - f_2 * sgh1_0[k]
                 + f_3 * pc_x[k] * sgi_0[k];

        t_1[k] = f_3 * pc_y[k] * sgi_0[k];

        t_2[k] = f_3 * pc_z[k] * sgi_0[k];

        t_3[k] = f_0 * sfi_3[k]
                 + f_4 * sgh0_3[k]
                 - f_5 * sgh1_3[k]
                 + f_3 * pc_x[k] * sgi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sfi_5, sfi_6, sgh0_5, sgh0_6, sgh1_5, \
                         sgh1_6, sgi_2, sgi_5, sgi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sgi_2[k];

        t_5[k] = f_0 * sfi_5[k]
                 + f_4 * sgh0_5[k]
                 - f_5 * sgh1_5[k]
                 + f_3 * pc_x[k] * sgi_5[k];

        t_6[k] = f_0 * sfi_6[k]
                 + f_6 * sgh0_6[k]
                 - f_7 * sgh1_6[k]
                 + f_3 * pc_x[k] * sgi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sfi_9, sgh0_9, sgh1_9, sgi_3, sgi_5, \
                         sgi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sgi_3[k];

        t_8[k] = f_3 * pc_y[k] * sgi_5[k];

        t_9[k] = f_0 * sfi_9[k]
                 + f_6 * sgh0_9[k]
                 - f_7 * sgh1_9[k]
                 + f_3 * pc_x[k] * sgi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sfi_10, sfi_12, sgh0_10, sgh0_12, \
                         sgh1_10, sgh1_12, sgi_6, sgi_10, sgi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sfi_10[k]
                  + f_8 * sgh0_10[k]
                  - f_9 * sgh1_10[k]
                  + f_3 * pc_x[k] * sgi_10[k];

        t_11[k] = f_3 * pc_z[k] * sgi_6[k];

        t_12[k] = f_0 * sfi_12[k]
                  + f_8 * sgh0_12[k]
                  - f_9 * sgh1_12[k]
                  + f_3 * pc_x[k] * sgi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, sfi_14, sfi_15, sgh0_14, sgh0_15, \
                         sgh1_14, sgh1_15, sgi_9, sgi_14, sgi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sgi_9[k];

        t_14[k] = f_0 * sfi_14[k]
                  + f_8 * sgh0_14[k]
                  - f_9 * sgh1_14[k]
                  + f_3 * pc_x[k] * sgi_14[k];

        t_15[k] = f_0 * sfi_15[k]
                  + f_10 * sgh0_15[k]
                  - f_11 * sgh1_15[k]
                  + f_3 * pc_x[k] * sgi_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, sfi_17, sfi_18, sgh0_17, sgh0_18, \
                         sgh1_17, sgh1_18, sgi_10, sgi_17, sgi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sgi_10[k];

        t_17[k] = f_0 * sfi_17[k]
                  + f_10 * sgh0_17[k]
                  - f_11 * sgh1_17[k]
                  + f_3 * pc_x[k] * sgi_17[k];

        t_18[k] = f_0 * sfi_18[k]
                  + f_10 * sgh0_18[k]
                  - f_11 * sgh1_18[k]
                  + f_3 * pc_x[k] * sgi_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pc_x, pc_y, sfi_20, sfi_21, sfi_22, sgh0_20, \
                         sgh1_20, sgi_14, sgi_20, sgi_21, sgi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sgi_14[k];

        t_20[k] = f_0 * sfi_20[k]
                  + f_10 * sgh0_20[k]
                  - f_11 * sgh1_20[k]
                  + f_3 * pc_x[k] * sgi_20[k];

        t_21[k] = f_0 * sfi_21[k]
                  + f_3 * pc_x[k] * sgi_21[k];

        t_22[k] = f_0 * sfi_22[k]
                  + f_3 * pc_x[k] * sgi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, sfi_23, sfi_24, sfi_25, sfi_26, \
                         sfi_27, sgi_23, sgi_24, sgi_25, sgi_26, \
                         sgi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * sfi_23[k]
                  + f_3 * pc_x[k] * sgi_23[k];

        t_24[k] = f_0 * sfi_24[k]
                  + f_3 * pc_x[k] * sgi_24[k];

        t_25[k] = f_0 * sfi_25[k]
                  + f_3 * pc_x[k] * sgi_25[k];

        t_26[k] = f_0 * sfi_26[k]
                  + f_3 * pc_x[k] * sgi_26[k];

        t_27[k] = f_0 * sfi_27[k]
                  + f_3 * pc_x[k] * sgi_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, sgh0_15, sgh0_17, sgh0_18, \
                         sgh1_15, sgh1_17, sgh1_18, sgi_21, sgi_23, \
                         sgi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * sgh0_15[k]
                  - f_2 * sgh1_15[k]
                  + f_3 * pc_y[k] * sgi_21[k];

        t_29[k] = f_3 * pc_z[k] * sgi_21[k];

        t_30[k] = f_4 * sgh0_17[k]
                  - f_5 * sgh1_17[k]
                  + f_3 * pc_y[k] * sgi_23[k];

        t_31[k] = f_6 * sgh0_18[k]
                  - f_7 * sgh1_18[k]
                  + f_3 * pc_y[k] * sgi_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, sgh0_19, sgh0_20, sgh1_19, \
                         sgh1_20, sgi_25, sgi_26, sgi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * sgh0_19[k]
                  - f_9 * sgh1_19[k]
                  + f_3 * pc_y[k] * sgi_25[k];

        t_33[k] = f_10 * sgh0_20[k]
                  - f_11 * sgh1_20[k]
                  + f_3 * pc_y[k] * sgi_26[k];

        t_34[k] = f_3 * pc_y[k] * sgi_27[k];

        t_35[k] = f_1 * sgh0_20[k]
                  - f_2 * sgh1_20[k]
                  + f_3 * pc_z[k] * sgi_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, sfk0_0, sfk0_3, sfi_0, \
                         sfi_1, sfk1_0, sfk1_3, sgi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * sfk0_0[k]
                  - f_12 * pc_y[k] * sfk1_0[k];

        t_37[k] = f_13 * sfi_0[k]
                  + f_3 * pc_y[k] * sgi_28[k];

        t_38[k] = f_3 * pc_z[k] * sgi_28[k];

        t_39[k] = pb_y[k] * sfk0_3[k]
                  + f_14 * sfi_1[k]
                  - f_12 * pc_y[k] * sfk1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_y, pc_z, sfk0_5, sfk0_6, sfi_2, \
                         sfi_3, sfk1_5, sfk1_6, sgi_30, sgi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * sfi_2[k]
                  + f_3 * pc_y[k] * sgi_30[k];

        t_41[k] = pb_y[k] * sfk0_5[k]
                  - f_12 * pc_y[k] * sfk1_5[k];

        t_42[k] = pb_y[k] * sfk0_6[k]
                  + f_15 * sfi_3[k]
                  - f_12 * pc_y[k] * sfk1_6[k];

        t_43[k] = f_3 * pc_z[k] * sgi_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pc_y, pc_z, sfk0_9, sfk0_10, sfi_5, \
                         sfi_6, sfk1_9, sfk1_10, sgi_33, sgi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_13 * sfi_5[k]
                  + f_3 * pc_y[k] * sgi_33[k];

        t_45[k] = pb_y[k] * sfk0_9[k]
                  - f_12 * pc_y[k] * sfk1_9[k];

        t_46[k] = pb_y[k] * sfk0_10[k]
                  + f_0 * sfi_6[k]
                  - f_12 * pc_y[k] * sfk1_10[k];

        t_47[k] = f_3 * pc_z[k] * sgi_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_y, pc_y, sfk0_12, sfk0_14, sfk0_15, sfi_8, \
                         sfi_9, sfi_10, sfk1_12, sfk1_14, sfk1_15, \
                         sgi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * sfk0_12[k]
                  + f_14 * sfi_8[k]
                  - f_12 * pc_y[k] * sfk1_12[k];

        t_49[k] = f_13 * sfi_9[k]
                  + f_3 * pc_y[k] * sgi_37[k];

        t_50[k] = pb_y[k] * sfk0_14[k]
                  - f_12 * pc_y[k] * sfk1_14[k];

        t_51[k] = pb_y[k] * sfk0_15[k]
                  + f_16 * sfi_10[k]
                  - f_12 * pc_y[k] * sfk1_15[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, pc_y, pc_z, sfk0_17, sfk0_18, sfi_12, \
                         sfi_13, sfi_14, sfk1_17, sfk1_18, sgi_38, \
                         sgi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * sgi_38[k];

        t_53[k] = pb_y[k] * sfk0_17[k]
                  + f_15 * sfi_12[k]
                  - f_12 * pc_y[k] * sfk1_17[k];

        t_54[k] = pb_y[k] * sfk0_18[k]
                  + f_14 * sfi_13[k]
                  - f_12 * pc_y[k] * sfk1_18[k];

        t_55[k] = f_13 * sfi_14[k]
                  + f_3 * pc_y[k] * sgi_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pc_x, pc_y, sfk0_20, sfi_49, sfi_50, \
                         sfi_51, sfk1_20, sgi_49, sgi_50, sgi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_y[k] * sfk0_20[k]
                  - f_12 * pc_y[k] * sfk1_20[k];

        t_57[k] = f_15 * sfi_49[k]
                  + f_3 * pc_x[k] * sgi_49[k];

        t_58[k] = f_15 * sfi_50[k]
                  + f_3 * pc_x[k] * sgi_50[k];

        t_59[k] = f_15 * sfi_51[k]
                  + f_3 * pc_x[k] * sgi_51[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, sfi_52, sfi_53, sfi_54, sfi_55, sgi_52, \
                         sgi_53, sgi_54, sgi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_15 * sfi_52[k]
                  + f_3 * pc_x[k] * sgi_52[k];

        t_61[k] = f_15 * sfi_53[k]
                  + f_3 * pc_x[k] * sgi_53[k];

        t_62[k] = f_15 * sfi_54[k]
                  + f_3 * pc_x[k] * sgi_54[k];

        t_63[k] = f_15 * sfi_55[k]
                  + f_3 * pc_x[k] * sgi_55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pc_y, pc_z, sfi_21, sfi_23, sgh0_36, sgh0_38, \
                         sgh1_36, sgh1_38, sgi_49, sgi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_13 * sfi_21[k]
                  + f_1 * sgh0_36[k]
                  - f_2 * sgh1_36[k]
                  + f_3 * pc_y[k] * sgi_49[k];

        t_65[k] = f_3 * pc_z[k] * sgi_49[k];

        t_66[k] = f_13 * sfi_23[k]
                  + f_4 * sgh0_38[k]
                  - f_5 * sgh1_38[k]
                  + f_3 * pc_y[k] * sgi_51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, sfi_24, sfi_25, sfi_26, sgh0_39, sgh0_40, \
                         sgh0_41, sgh1_39, sgh1_40, sgh1_41, sgi_52, sgi_53, \
                         sgi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_13 * sfi_24[k]
                  + f_6 * sgh0_39[k]
                  - f_7 * sgh1_39[k]
                  + f_3 * pc_y[k] * sgi_52[k];

        t_68[k] = f_13 * sfi_25[k]
                  + f_8 * sgh0_40[k]
                  - f_9 * sgh1_40[k]
                  + f_3 * pc_y[k] * sgi_53[k];

        t_69[k] = f_13 * sfi_26[k]
                  + f_10 * sgh0_41[k]
                  - f_11 * sgh1_41[k]
                  + f_3 * pc_y[k] * sgi_54[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_y, pb_z, pc_y, pc_z, sfk0_0, sfk0_35, \
                         sfi_27, sfk1_0, sfk1_35, sgi_55, sgi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * sfi_27[k]
                  + f_3 * pc_y[k] * sgi_55[k];

        t_71[k] = pb_y[k] * sfk0_35[k]
                  - f_12 * pc_y[k] * sfk1_35[k];

        t_72[k] = pb_z[k] * sfk0_0[k]
                  - f_12 * pc_z[k] * sfk1_0[k];

        t_73[k] = f_3 * pc_y[k] * sgi_56[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_y, pc_z, sfk0_3, sfk0_5, sfi_0, \
                         sfi_2, sfk1_3, sfk1_5, sgi_56, sgi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * sfi_0[k]
                  + f_3 * pc_z[k] * sgi_56[k];

        t_75[k] = pb_z[k] * sfk0_3[k]
                  - f_12 * pc_z[k] * sfk1_3[k];

        t_76[k] = f_3 * pc_y[k] * sgi_58[k];

        t_77[k] = pb_z[k] * sfk0_5[k]
                  + f_14 * sfi_2[k]
                  - f_12 * pc_z[k] * sfk1_5[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_z, pc_y, pc_z, sfk0_6, sfk0_9, sfi_3, \
                         sfi_5, sfk1_6, sfk1_9, sgi_59, sgi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_z[k] * sfk0_6[k]
                  - f_12 * pc_z[k] * sfk1_6[k];

        t_79[k] = f_13 * sfi_3[k]
                  + f_3 * pc_z[k] * sgi_59[k];

        t_80[k] = f_3 * pc_y[k] * sgi_61[k];

        t_81[k] = pb_z[k] * sfk0_9[k]
                  + f_15 * sfi_5[k]
                  - f_12 * pc_z[k] * sfk1_9[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_z, pc_y, pc_z, sfk0_10, sfk0_12, sfi_6, \
                         sfi_7, sfk1_10, sfk1_12, sgi_62, sgi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * sfk0_10[k]
                  - f_12 * pc_z[k] * sfk1_10[k];

        t_83[k] = f_13 * sfi_6[k]
                  + f_3 * pc_z[k] * sgi_62[k];

        t_84[k] = pb_z[k] * sfk0_12[k]
                  + f_14 * sfi_7[k]
                  - f_12 * pc_z[k] * sfk1_12[k];

        t_85[k] = f_3 * pc_y[k] * sgi_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_z, pc_z, sfk0_14, sfk0_15, sfk0_17, sfi_9, \
                         sfi_10, sfi_11, sfk1_14, sfk1_15, sfk1_17, \
                         sgi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * sfk0_14[k]
                  + f_0 * sfi_9[k]
                  - f_12 * pc_z[k] * sfk1_14[k];

        t_87[k] = pb_z[k] * sfk0_15[k]
                  - f_12 * pc_z[k] * sfk1_15[k];

        t_88[k] = f_13 * sfi_10[k]
                  + f_3 * pc_z[k] * sgi_66[k];

        t_89[k] = pb_z[k] * sfk0_17[k]
                  + f_14 * sfi_11[k]
                  - f_12 * pc_z[k] * sfk1_17[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_z, pc_y, pc_z, sfk0_18, sfk0_20, sfi_12, sfi_14, \
                         sfk1_18, sfk1_20, sgi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * sfk0_18[k]
                  + f_15 * sfi_12[k]
                  - f_12 * pc_z[k] * sfk1_18[k];

        t_91[k] = f_3 * pc_y[k] * sgi_70[k];

        t_92[k] = pb_z[k] * sfk0_20[k]
                  + f_16 * sfi_14[k]
                  - f_12 * pc_z[k] * sfk1_20[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pc_x, sfi_77, sfi_78, sfi_79, sfi_80, \
                         sfi_81, sgi_77, sgi_78, sgi_79, sgi_80, \
                         sgi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_15 * sfi_77[k]
                  + f_3 * pc_x[k] * sgi_77[k];

        t_94[k] = f_15 * sfi_78[k]
                  + f_3 * pc_x[k] * sgi_78[k];

        t_95[k] = f_15 * sfi_79[k]
                  + f_3 * pc_x[k] * sgi_79[k];

        t_96[k] = f_15 * sfi_80[k]
                  + f_3 * pc_x[k] * sgi_80[k];

        t_97[k] = f_15 * sfi_81[k]
                  + f_3 * pc_x[k] * sgi_81[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_z, sfk0_28, sfi_21, sfi_82, \
                         sfi_83, sfk1_28, sgi_77, sgi_82, sgi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_15 * sfi_82[k]
                  + f_3 * pc_x[k] * sgi_82[k];

        t_99[k] = f_15 * sfi_83[k]
                  + f_3 * pc_x[k] * sgi_83[k];

        t_100[k] = pb_z[k] * sfk0_28[k]
                   - f_12 * pc_z[k] * sfk1_28[k];

        t_101[k] = f_13 * sfi_21[k]
                   + f_3 * pc_z[k] * sgi_77[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, sgh0_59, sgh0_60, sgh0_61, sgh1_59, \
                         sgh1_60, sgh1_61, sgi_79, sgi_80, sgi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * sgh0_59[k]
                   - f_5 * sgh1_59[k]
                   + f_3 * pc_y[k] * sgi_79[k];

        t_103[k] = f_6 * sgh0_60[k]
                   - f_7 * sgh1_60[k]
                   + f_3 * pc_y[k] * sgi_80[k];

        t_104[k] = f_8 * sgh0_61[k]
                   - f_9 * sgh1_61[k]
                   + f_3 * pc_y[k] * sgi_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_x, pc_y, pc_z, sfi_27, sfi_84, \
                         sgh0_62, sgh0_63, sgh1_62, sgh1_63, sgi_82, sgi_83, \
                         sgi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_10 * sgh0_62[k]
                   - f_11 * sgh1_62[k]
                   + f_3 * pc_y[k] * sgi_82[k];

        t_106[k] = f_3 * pc_y[k] * sgi_83[k];

        t_107[k] = f_13 * sfi_27[k]
                   + f_1 * sgh0_62[k]
                   - f_2 * sgh1_62[k]
                   + f_3 * pc_z[k] * sgi_83[k];

        t_108[k] = f_14 * sfi_84[k]
                   + f_1 * sgh0_63[k]
                   - f_2 * sgh1_63[k]
                   + f_3 * pc_x[k] * sgi_84[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, pc_z, sfi_28, sfi_30, sfi_87, \
                         sgh0_66, sgh1_66, sgi_84, sgi_86, sgi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_14 * sfi_28[k]
                   + f_3 * pc_y[k] * sgi_84[k];

        t_110[k] = f_3 * pc_z[k] * sgi_84[k];

        t_111[k] = f_14 * sfi_87[k]
                   + f_4 * sgh0_66[k]
                   - f_5 * sgh1_66[k]
                   + f_3 * pc_x[k] * sgi_87[k];

        t_112[k] = f_14 * sfi_30[k]
                   + f_3 * pc_y[k] * sgi_86[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pc_x, pc_z, sfi_89, sfi_90, sgh0_68, sgh0_69, \
                         sgh1_68, sgh1_69, sgi_87, sgi_89, sgi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_14 * sfi_89[k]
                   + f_4 * sgh0_68[k]
                   - f_5 * sgh1_68[k]
                   + f_3 * pc_x[k] * sgi_89[k];

        t_114[k] = f_14 * sfi_90[k]
                   + f_6 * sgh0_69[k]
                   - f_7 * sgh1_69[k]
                   + f_3 * pc_x[k] * sgi_90[k];

        t_115[k] = f_3 * pc_z[k] * sgi_87[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pc_x, pc_y, sfi_33, sfi_93, sfi_94, sgh0_72, \
                         sgh0_73, sgh1_72, sgh1_73, sgi_89, sgi_93, \
                         sgi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_14 * sfi_33[k]
                   + f_3 * pc_y[k] * sgi_89[k];

        t_117[k] = f_14 * sfi_93[k]
                   + f_6 * sgh0_72[k]
                   - f_7 * sgh1_72[k]
                   + f_3 * pc_x[k] * sgi_93[k];

        t_118[k] = f_14 * sfi_94[k]
                   + f_8 * sgh0_73[k]
                   - f_9 * sgh1_73[k]
                   + f_3 * pc_x[k] * sgi_94[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, sfi_37, sfi_96, sgh0_75, \
                         sgh1_75, sgi_90, sgi_93, sgi_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * pc_z[k] * sgi_90[k];

        t_120[k] = f_14 * sfi_96[k]
                   + f_8 * sgh0_75[k]
                   - f_9 * sgh1_75[k]
                   + f_3 * pc_x[k] * sgi_96[k];

        t_121[k] = f_14 * sfi_37[k]
                   + f_3 * pc_y[k] * sgi_93[k];
    }
}

static auto
compute_prim_sgk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfk0,
                                                          const size_t sfi, const size_t sfk1,
                                                          const size_t sgh0, const size_t sgh1,
                                                          const size_t sgi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfk0_39 = buffer.data(sfk0 + 39);
    const auto *sfk0_42 = buffer.data(sfk0 + 42);
    const auto *sfk0_46 = buffer.data(sfk0 + 46);
    const auto *sfk0_51 = buffer.data(sfk0 + 51);
    const auto *sfk0_64 = buffer.data(sfk0 + 64);
    const auto *sfk0_72 = buffer.data(sfk0 + 72);
    const auto *sfk0_77 = buffer.data(sfk0 + 77);
    const auto *sfk0_81 = buffer.data(sfk0 + 81);
    const auto *sfk0_84 = buffer.data(sfk0 + 84);
    const auto *sfk0_86 = buffer.data(sfk0 + 86);
    const auto *sfk0_89 = buffer.data(sfk0 + 89);
    const auto *sfk0_90 = buffer.data(sfk0 + 90);
    const auto *sfk0_92 = buffer.data(sfk0 + 92);
    const auto *sfk0_107 = buffer.data(sfk0 + 107);
    const auto *sfk0_216 = buffer.data(sfk0 + 216);
    const auto *sfk0_219 = buffer.data(sfk0 + 219);
    const auto *sfk0_221 = buffer.data(sfk0 + 221);
    const auto *sfk0_222 = buffer.data(sfk0 + 222);
    const auto *sfk0_225 = buffer.data(sfk0 + 225);
    const auto *sfk0_226 = buffer.data(sfk0 + 226);
    const auto *sfk0_228 = buffer.data(sfk0 + 228);
    const auto *sfk0_230 = buffer.data(sfk0 + 230);
    const auto *sfk0_231 = buffer.data(sfk0 + 231);
    const auto *sfk0_233 = buffer.data(sfk0 + 233);
    const auto *sfk0_234 = buffer.data(sfk0 + 234);
    const auto *sfk0_236 = buffer.data(sfk0 + 236);

    const auto *sfi_28 = buffer.data(sfi + 28);
    const auto *sfi_31 = buffer.data(sfi + 31);
    const auto *sfi_34 = buffer.data(sfi + 34);
    const auto *sfi_38 = buffer.data(sfi + 38);
    const auto *sfi_42 = buffer.data(sfi + 42);
    const auto *sfi_49 = buffer.data(sfi + 49);
    const auto *sfi_51 = buffer.data(sfi + 51);
    const auto *sfi_52 = buffer.data(sfi + 52);
    const auto *sfi_53 = buffer.data(sfi + 53);
    const auto *sfi_54 = buffer.data(sfi + 54);
    const auto *sfi_55 = buffer.data(sfi + 55);
    const auto *sfi_56 = buffer.data(sfi + 56);
    const auto *sfi_58 = buffer.data(sfi + 58);
    const auto *sfi_59 = buffer.data(sfi + 59);
    const auto *sfi_61 = buffer.data(sfi + 61);
    const auto *sfi_62 = buffer.data(sfi + 62);
    const auto *sfi_64 = buffer.data(sfi + 64);
    const auto *sfi_65 = buffer.data(sfi + 65);
    const auto *sfi_66 = buffer.data(sfi + 66);
    const auto *sfi_68 = buffer.data(sfi + 68);
    const auto *sfi_69 = buffer.data(sfi + 69);
    const auto *sfi_70 = buffer.data(sfi + 70);
    const auto *sfi_77 = buffer.data(sfi + 77);
    const auto *sfi_79 = buffer.data(sfi + 79);
    const auto *sfi_80 = buffer.data(sfi + 80);
    const auto *sfi_81 = buffer.data(sfi + 81);
    const auto *sfi_82 = buffer.data(sfi + 82);
    const auto *sfi_83 = buffer.data(sfi + 83);
    const auto *sfi_84 = buffer.data(sfi + 84);
    const auto *sfi_86 = buffer.data(sfi + 86);
    const auto *sfi_89 = buffer.data(sfi + 89);
    const auto *sfi_93 = buffer.data(sfi + 93);
    const auto *sfi_98 = buffer.data(sfi + 98);
    const auto *sfi_99 = buffer.data(sfi + 99);
    const auto *sfi_101 = buffer.data(sfi + 101);
    const auto *sfi_102 = buffer.data(sfi + 102);
    const auto *sfi_104 = buffer.data(sfi + 104);
    const auto *sfi_105 = buffer.data(sfi + 105);
    const auto *sfi_106 = buffer.data(sfi + 106);
    const auto *sfi_107 = buffer.data(sfi + 107);
    const auto *sfi_108 = buffer.data(sfi + 108);
    const auto *sfi_109 = buffer.data(sfi + 109);
    const auto *sfi_110 = buffer.data(sfi + 110);
    const auto *sfi_111 = buffer.data(sfi + 111);
    const auto *sfi_133 = buffer.data(sfi + 133);
    const auto *sfi_134 = buffer.data(sfi + 134);
    const auto *sfi_135 = buffer.data(sfi + 135);
    const auto *sfi_136 = buffer.data(sfi + 136);
    const auto *sfi_137 = buffer.data(sfi + 137);
    const auto *sfi_138 = buffer.data(sfi + 138);
    const auto *sfi_139 = buffer.data(sfi + 139);
    const auto *sfi_140 = buffer.data(sfi + 140);
    const auto *sfi_143 = buffer.data(sfi + 143);
    const auto *sfi_145 = buffer.data(sfi + 145);
    const auto *sfi_146 = buffer.data(sfi + 146);
    const auto *sfi_149 = buffer.data(sfi + 149);
    const auto *sfi_150 = buffer.data(sfi + 150);
    const auto *sfi_152 = buffer.data(sfi + 152);
    const auto *sfi_154 = buffer.data(sfi + 154);
    const auto *sfi_155 = buffer.data(sfi + 155);
    const auto *sfi_157 = buffer.data(sfi + 157);
    const auto *sfi_158 = buffer.data(sfi + 158);
    const auto *sfi_160 = buffer.data(sfi + 160);
    const auto *sfi_161 = buffer.data(sfi + 161);
    const auto *sfi_162 = buffer.data(sfi + 162);
    const auto *sfi_163 = buffer.data(sfi + 163);
    const auto *sfi_164 = buffer.data(sfi + 164);
    const auto *sfi_165 = buffer.data(sfi + 165);
    const auto *sfi_166 = buffer.data(sfi + 166);
    const auto *sfi_167 = buffer.data(sfi + 167);
    const auto *sfi_168 = buffer.data(sfi + 168);
    const auto *sfi_171 = buffer.data(sfi + 171);
    const auto *sfi_173 = buffer.data(sfi + 173);
    const auto *sfi_174 = buffer.data(sfi + 174);
    const auto *sfi_177 = buffer.data(sfi + 177);
    const auto *sfi_178 = buffer.data(sfi + 178);
    const auto *sfi_180 = buffer.data(sfi + 180);
    const auto *sfi_182 = buffer.data(sfi + 182);
    const auto *sfi_183 = buffer.data(sfi + 183);
    const auto *sfi_185 = buffer.data(sfi + 185);
    const auto *sfi_186 = buffer.data(sfi + 186);
    const auto *sfi_188 = buffer.data(sfi + 188);
    const auto *sfi_189 = buffer.data(sfi + 189);
    const auto *sfi_190 = buffer.data(sfi + 190);

    const auto *sfk1_39 = buffer.data(sfk1 + 39);
    const auto *sfk1_42 = buffer.data(sfk1 + 42);
    const auto *sfk1_46 = buffer.data(sfk1 + 46);
    const auto *sfk1_51 = buffer.data(sfk1 + 51);
    const auto *sfk1_64 = buffer.data(sfk1 + 64);
    const auto *sfk1_72 = buffer.data(sfk1 + 72);
    const auto *sfk1_77 = buffer.data(sfk1 + 77);
    const auto *sfk1_81 = buffer.data(sfk1 + 81);
    const auto *sfk1_84 = buffer.data(sfk1 + 84);
    const auto *sfk1_86 = buffer.data(sfk1 + 86);
    const auto *sfk1_89 = buffer.data(sfk1 + 89);
    const auto *sfk1_90 = buffer.data(sfk1 + 90);
    const auto *sfk1_92 = buffer.data(sfk1 + 92);
    const auto *sfk1_107 = buffer.data(sfk1 + 107);
    const auto *sfk1_216 = buffer.data(sfk1 + 216);
    const auto *sfk1_219 = buffer.data(sfk1 + 219);
    const auto *sfk1_221 = buffer.data(sfk1 + 221);
    const auto *sfk1_222 = buffer.data(sfk1 + 222);
    const auto *sfk1_225 = buffer.data(sfk1 + 225);
    const auto *sfk1_226 = buffer.data(sfk1 + 226);
    const auto *sfk1_228 = buffer.data(sfk1 + 228);
    const auto *sfk1_230 = buffer.data(sfk1 + 230);
    const auto *sfk1_231 = buffer.data(sfk1 + 231);
    const auto *sfk1_233 = buffer.data(sfk1 + 233);
    const auto *sfk1_234 = buffer.data(sfk1 + 234);
    const auto *sfk1_236 = buffer.data(sfk1 + 236);

    const auto *sgh0_77 = buffer.data(sgh0 + 77);
    const auto *sgh0_78 = buffer.data(sgh0 + 78);
    const auto *sgh0_80 = buffer.data(sgh0 + 80);
    const auto *sgh0_81 = buffer.data(sgh0 + 81);
    const auto *sgh0_82 = buffer.data(sgh0 + 82);
    const auto *sgh0_83 = buffer.data(sgh0 + 83);
    const auto *sgh0_101 = buffer.data(sgh0 + 101);
    const auto *sgh0_102 = buffer.data(sgh0 + 102);
    const auto *sgh0_103 = buffer.data(sgh0 + 103);
    const auto *sgh0_104 = buffer.data(sgh0 + 104);
    const auto *sgh0_105 = buffer.data(sgh0 + 105);
    const auto *sgh0_108 = buffer.data(sgh0 + 108);
    const auto *sgh0_110 = buffer.data(sgh0 + 110);
    const auto *sgh0_111 = buffer.data(sgh0 + 111);
    const auto *sgh0_114 = buffer.data(sgh0 + 114);
    const auto *sgh0_115 = buffer.data(sgh0 + 115);
    const auto *sgh0_117 = buffer.data(sgh0 + 117);
    const auto *sgh0_119 = buffer.data(sgh0 + 119);
    const auto *sgh0_120 = buffer.data(sgh0 + 120);
    const auto *sgh0_122 = buffer.data(sgh0 + 122);
    const auto *sgh0_123 = buffer.data(sgh0 + 123);
    const auto *sgh0_124 = buffer.data(sgh0 + 124);
    const auto *sgh0_125 = buffer.data(sgh0 + 125);

    const auto *sgh1_77 = buffer.data(sgh1 + 77);
    const auto *sgh1_78 = buffer.data(sgh1 + 78);
    const auto *sgh1_80 = buffer.data(sgh1 + 80);
    const auto *sgh1_81 = buffer.data(sgh1 + 81);
    const auto *sgh1_82 = buffer.data(sgh1 + 82);
    const auto *sgh1_83 = buffer.data(sgh1 + 83);
    const auto *sgh1_101 = buffer.data(sgh1 + 101);
    const auto *sgh1_102 = buffer.data(sgh1 + 102);
    const auto *sgh1_103 = buffer.data(sgh1 + 103);
    const auto *sgh1_104 = buffer.data(sgh1 + 104);
    const auto *sgh1_105 = buffer.data(sgh1 + 105);
    const auto *sgh1_108 = buffer.data(sgh1 + 108);
    const auto *sgh1_110 = buffer.data(sgh1 + 110);
    const auto *sgh1_111 = buffer.data(sgh1 + 111);
    const auto *sgh1_114 = buffer.data(sgh1 + 114);
    const auto *sgh1_115 = buffer.data(sgh1 + 115);
    const auto *sgh1_117 = buffer.data(sgh1 + 117);
    const auto *sgh1_119 = buffer.data(sgh1 + 119);
    const auto *sgh1_120 = buffer.data(sgh1 + 120);
    const auto *sgh1_122 = buffer.data(sgh1 + 122);
    const auto *sgh1_123 = buffer.data(sgh1 + 123);
    const auto *sgh1_124 = buffer.data(sgh1 + 124);
    const auto *sgh1_125 = buffer.data(sgh1 + 125);

    const auto *sgi_94 = buffer.data(sgi + 94);
    const auto *sgi_98 = buffer.data(sgi + 98);
    const auto *sgi_99 = buffer.data(sgi + 99);
    const auto *sgi_101 = buffer.data(sgi + 101);
    const auto *sgi_102 = buffer.data(sgi + 102);
    const auto *sgi_104 = buffer.data(sgi + 104);
    const auto *sgi_105 = buffer.data(sgi + 105);
    const auto *sgi_106 = buffer.data(sgi + 106);
    const auto *sgi_107 = buffer.data(sgi + 107);
    const auto *sgi_108 = buffer.data(sgi + 108);
    const auto *sgi_109 = buffer.data(sgi + 109);
    const auto *sgi_110 = buffer.data(sgi + 110);
    const auto *sgi_111 = buffer.data(sgi + 111);
    const auto *sgi_112 = buffer.data(sgi + 112);
    const auto *sgi_114 = buffer.data(sgi + 114);
    const auto *sgi_115 = buffer.data(sgi + 115);
    const auto *sgi_117 = buffer.data(sgi + 117);
    const auto *sgi_118 = buffer.data(sgi + 118);
    const auto *sgi_121 = buffer.data(sgi + 121);
    const auto *sgi_122 = buffer.data(sgi + 122);
    const auto *sgi_126 = buffer.data(sgi + 126);
    const auto *sgi_133 = buffer.data(sgi + 133);
    const auto *sgi_134 = buffer.data(sgi + 134);
    const auto *sgi_135 = buffer.data(sgi + 135);
    const auto *sgi_136 = buffer.data(sgi + 136);
    const auto *sgi_137 = buffer.data(sgi + 137);
    const auto *sgi_138 = buffer.data(sgi + 138);
    const auto *sgi_139 = buffer.data(sgi + 139);
    const auto *sgi_140 = buffer.data(sgi + 140);
    const auto *sgi_142 = buffer.data(sgi + 142);
    const auto *sgi_143 = buffer.data(sgi + 143);
    const auto *sgi_145 = buffer.data(sgi + 145);
    const auto *sgi_146 = buffer.data(sgi + 146);
    const auto *sgi_149 = buffer.data(sgi + 149);
    const auto *sgi_150 = buffer.data(sgi + 150);
    const auto *sgi_152 = buffer.data(sgi + 152);
    const auto *sgi_154 = buffer.data(sgi + 154);
    const auto *sgi_155 = buffer.data(sgi + 155);
    const auto *sgi_157 = buffer.data(sgi + 157);
    const auto *sgi_158 = buffer.data(sgi + 158);
    const auto *sgi_160 = buffer.data(sgi + 160);
    const auto *sgi_161 = buffer.data(sgi + 161);
    const auto *sgi_162 = buffer.data(sgi + 162);
    const auto *sgi_163 = buffer.data(sgi + 163);
    const auto *sgi_164 = buffer.data(sgi + 164);
    const auto *sgi_165 = buffer.data(sgi + 165);
    const auto *sgi_166 = buffer.data(sgi + 166);
    const auto *sgi_167 = buffer.data(sgi + 167);
    const auto *sgi_168 = buffer.data(sgi + 168);
    const auto *sgi_170 = buffer.data(sgi + 170);
    const auto *sgi_171 = buffer.data(sgi + 171);
    const auto *sgi_173 = buffer.data(sgi + 173);
    const auto *sgi_174 = buffer.data(sgi + 174);
    const auto *sgi_177 = buffer.data(sgi + 177);
    const auto *sgi_178 = buffer.data(sgi + 178);
    const auto *sgi_182 = buffer.data(sgi + 182);
    const auto *sgi_189 = buffer.data(sgi + 189);
    const auto *sgi_190 = buffer.data(sgi + 190);

#pragma omp simd aligned(t_122, t_123, t_124, pc_x, pc_z, sfi_98, sfi_99, sgh0_77, sgh0_78, \
                         sgh1_77, sgh1_78, sgi_94, sgi_98, sgi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_14 * sfi_98[k]
                   + f_8 * sgh0_77[k]
                   - f_9 * sgh1_77[k]
                   + f_3 * pc_x[k] * sgi_98[k];

        t_123[k] = f_14 * sfi_99[k]
                   + f_10 * sgh0_78[k]
                   - f_11 * sgh1_78[k]
                   + f_3 * pc_x[k] * sgi_99[k];

        t_124[k] = f_3 * pc_z[k] * sgi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_y, sfi_42, sfi_101, sfi_102, sgh0_80, \
                         sgh0_81, sgh1_80, sgh1_81, sgi_98, sgi_101, \
                         sgi_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_14 * sfi_101[k]
                   + f_10 * sgh0_80[k]
                   - f_11 * sgh1_80[k]
                   + f_3 * pc_x[k] * sgi_101[k];

        t_126[k] = f_14 * sfi_102[k]
                   + f_10 * sgh0_81[k]
                   - f_11 * sgh1_81[k]
                   + f_3 * pc_x[k] * sgi_102[k];

        t_127[k] = f_14 * sfi_42[k]
                   + f_3 * pc_y[k] * sgi_98[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, sfi_104, sfi_105, sfi_106, sfi_107, \
                         sgh0_83, sgh1_83, sgi_104, sgi_105, sgi_106, \
                         sgi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_14 * sfi_104[k]
                   + f_10 * sgh0_83[k]
                   - f_11 * sgh1_83[k]
                   + f_3 * pc_x[k] * sgi_104[k];

        t_129[k] = f_14 * sfi_105[k]
                   + f_3 * pc_x[k] * sgi_105[k];

        t_130[k] = f_14 * sfi_106[k]
                   + f_3 * pc_x[k] * sgi_106[k];

        t_131[k] = f_14 * sfi_107[k]
                   + f_3 * pc_x[k] * sgi_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, sfi_108, sfi_109, sfi_110, sfi_111, \
                         sgi_108, sgi_109, sgi_110, sgi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_14 * sfi_108[k]
                   + f_3 * pc_x[k] * sgi_108[k];

        t_133[k] = f_14 * sfi_109[k]
                   + f_3 * pc_x[k] * sgi_109[k];

        t_134[k] = f_14 * sfi_110[k]
                   + f_3 * pc_x[k] * sgi_110[k];

        t_135[k] = f_14 * sfi_111[k]
                   + f_3 * pc_x[k] * sgi_111[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, pc_z, sfi_49, sfi_51, sgh0_78, sgh0_80, \
                         sgh1_78, sgh1_80, sgi_105, sgi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_14 * sfi_49[k]
                   + f_1 * sgh0_78[k]
                   - f_2 * sgh1_78[k]
                   + f_3 * pc_y[k] * sgi_105[k];

        t_137[k] = f_3 * pc_z[k] * sgi_105[k];

        t_138[k] = f_14 * sfi_51[k]
                   + f_4 * sgh0_80[k]
                   - f_5 * sgh1_80[k]
                   + f_3 * pc_y[k] * sgi_107[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pc_y, sfi_52, sfi_53, sfi_54, sgh0_81, sgh0_82, \
                         sgh0_83, sgh1_81, sgh1_82, sgh1_83, sgi_108, sgi_109, \
                         sgi_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_14 * sfi_52[k]
                   + f_6 * sgh0_81[k]
                   - f_7 * sgh1_81[k]
                   + f_3 * pc_y[k] * sgi_108[k];

        t_140[k] = f_14 * sfi_53[k]
                   + f_8 * sgh0_82[k]
                   - f_9 * sgh1_82[k]
                   + f_3 * pc_y[k] * sgi_109[k];

        t_141[k] = f_14 * sfi_54[k]
                   + f_10 * sgh0_83[k]
                   - f_11 * sgh1_83[k]
                   + f_3 * pc_y[k] * sgi_110[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_y, pc_y, pc_z, sfk0_72, sfi_55, \
                         sfi_56, sfk1_72, sgh0_83, sgh1_83, sgi_111, \
                         sgi_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * sfi_55[k]
                   + f_3 * pc_y[k] * sgi_111[k];

        t_143[k] = f_1 * sgh0_83[k]
                   - f_2 * sgh1_83[k]
                   + f_3 * pc_z[k] * sgi_111[k];

        t_144[k] = pb_y[k] * sfk0_72[k]
                   - f_12 * pc_y[k] * sfk1_72[k];

        t_145[k] = f_13 * sfi_56[k]
                   + f_3 * pc_y[k] * sgi_112[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_y, pb_z, pc_y, pc_z, sfk0_39, sfk0_77, \
                         sfi_28, sfi_58, sfk1_39, sfk1_77, sgi_112, \
                         sgi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_13 * sfi_28[k]
                   + f_3 * pc_z[k] * sgi_112[k];

        t_147[k] = pb_z[k] * sfk0_39[k]
                   - f_12 * pc_z[k] * sfk1_39[k];

        t_148[k] = f_13 * sfi_58[k]
                   + f_3 * pc_y[k] * sgi_114[k];

        t_149[k] = pb_y[k] * sfk0_77[k]
                   - f_12 * pc_y[k] * sfk1_77[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, pb_z, pc_y, pc_z, sfk0_42, sfk0_81, \
                         sfi_31, sfi_61, sfk1_42, sfk1_81, sgi_115, \
                         sgi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_z[k] * sfk0_42[k]
                   - f_12 * pc_z[k] * sfk1_42[k];

        t_151[k] = f_13 * sfi_31[k]
                   + f_3 * pc_z[k] * sgi_115[k];

        t_152[k] = f_13 * sfi_61[k]
                   + f_3 * pc_y[k] * sgi_117[k];

        t_153[k] = pb_y[k] * sfk0_81[k]
                   - f_12 * pc_y[k] * sfk1_81[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pb_z, pc_y, pc_z, sfk0_46, sfk0_84, \
                         sfi_34, sfi_64, sfk1_46, sfk1_84, sgi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * sfk0_46[k]
                   - f_12 * pc_z[k] * sfk1_46[k];

        t_155[k] = f_13 * sfi_34[k]
                   + f_3 * pc_z[k] * sgi_118[k];

        t_156[k] = pb_y[k] * sfk0_84[k]
                   + f_14 * sfi_64[k]
                   - f_12 * pc_y[k] * sfk1_84[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pb_y, pb_z, pc_y, pc_z, sfk0_51, sfk0_86, \
                         sfi_38, sfi_65, sfk1_51, sfk1_86, sgi_121, \
                         sgi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_13 * sfi_65[k]
                   + f_3 * pc_y[k] * sgi_121[k];

        t_158[k] = pb_y[k] * sfk0_86[k]
                   - f_12 * pc_y[k] * sfk1_86[k];

        t_159[k] = pb_z[k] * sfk0_51[k]
                   - f_12 * pc_z[k] * sfk1_51[k];

        t_160[k] = f_13 * sfi_38[k]
                   + f_3 * pc_z[k] * sgi_122[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_y, pc_y, sfk0_89, sfk0_90, sfk0_92, \
                         sfi_68, sfi_69, sfi_70, sfk1_89, sfk1_90, sfk1_92, \
                         sgi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * sfk0_89[k]
                   + f_15 * sfi_68[k]
                   - f_12 * pc_y[k] * sfk1_89[k];

        t_162[k] = pb_y[k] * sfk0_90[k]
                   + f_14 * sfi_69[k]
                   - f_12 * pc_y[k] * sfk1_90[k];

        t_163[k] = f_13 * sfi_70[k]
                   + f_3 * pc_y[k] * sgi_126[k];

        t_164[k] = pb_y[k] * sfk0_92[k]
                   - f_12 * pc_y[k] * sfk1_92[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pc_x, sfi_133, sfi_134, sfi_135, \
                         sfi_136, sfi_137, sgi_133, sgi_134, sgi_135, sgi_136, \
                         sgi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_14 * sfi_133[k]
                   + f_3 * pc_x[k] * sgi_133[k];

        t_166[k] = f_14 * sfi_134[k]
                   + f_3 * pc_x[k] * sgi_134[k];

        t_167[k] = f_14 * sfi_135[k]
                   + f_3 * pc_x[k] * sgi_135[k];

        t_168[k] = f_14 * sfi_136[k]
                   + f_3 * pc_x[k] * sgi_136[k];

        t_169[k] = f_14 * sfi_137[k]
                   + f_3 * pc_x[k] * sgi_137[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_z, pc_x, pc_z, sfk0_64, sfi_49, \
                         sfi_138, sfi_139, sfk1_64, sgi_133, sgi_138, \
                         sgi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_14 * sfi_138[k]
                   + f_3 * pc_x[k] * sgi_138[k];

        t_171[k] = f_14 * sfi_139[k]
                   + f_3 * pc_x[k] * sgi_139[k];

        t_172[k] = pb_z[k] * sfk0_64[k]
                   - f_12 * pc_z[k] * sfk1_64[k];

        t_173[k] = f_13 * sfi_49[k]
                   + f_3 * pc_z[k] * sgi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, sfi_79, sfi_80, sfi_81, sgh0_101, \
                         sgh0_102, sgh0_103, sgh1_101, sgh1_102, sgh1_103, sgi_135, sgi_136, \
                         sgi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * sfi_79[k]
                   + f_4 * sgh0_101[k]
                   - f_5 * sgh1_101[k]
                   + f_3 * pc_y[k] * sgi_135[k];

        t_175[k] = f_13 * sfi_80[k]
                   + f_6 * sgh0_102[k]
                   - f_7 * sgh1_102[k]
                   + f_3 * pc_y[k] * sgi_136[k];

        t_176[k] = f_13 * sfi_81[k]
                   + f_8 * sgh0_103[k]
                   - f_9 * sgh1_103[k]
                   + f_3 * pc_y[k] * sgi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_y, pc_y, sfk0_107, sfi_82, sfi_83, sfk1_107, \
                         sgh0_104, sgh1_104, sgi_138, sgi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * sfi_82[k]
                   + f_10 * sgh0_104[k]
                   - f_11 * sgh1_104[k]
                   + f_3 * pc_y[k] * sgi_138[k];

        t_178[k] = f_13 * sfi_83[k]
                   + f_3 * pc_y[k] * sgi_139[k];

        t_179[k] = pb_y[k] * sfk0_107[k]
                   - f_12 * pc_y[k] * sfk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, sfi_56, sfi_140, \
                         sfi_143, sgh0_105, sgh0_108, sgh1_105, sgh1_108, sgi_140, \
                         sgi_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_14 * sfi_140[k]
                   + f_1 * sgh0_105[k]
                   - f_2 * sgh1_105[k]
                   + f_3 * pc_x[k] * sgi_140[k];

        t_181[k] = f_3 * pc_y[k] * sgi_140[k];

        t_182[k] = f_14 * sfi_56[k]
                   + f_3 * pc_z[k] * sgi_140[k];

        t_183[k] = f_14 * sfi_143[k]
                   + f_4 * sgh0_108[k]
                   - f_5 * sgh1_108[k]
                   + f_3 * pc_x[k] * sgi_143[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, sfi_145, sfi_146, sgh0_110, \
                         sgh0_111, sgh1_110, sgh1_111, sgi_142, sgi_145, \
                         sgi_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_y[k] * sgi_142[k];

        t_185[k] = f_14 * sfi_145[k]
                   + f_4 * sgh0_110[k]
                   - f_5 * sgh1_110[k]
                   + f_3 * pc_x[k] * sgi_145[k];

        t_186[k] = f_14 * sfi_146[k]
                   + f_6 * sgh0_111[k]
                   - f_7 * sgh1_111[k]
                   + f_3 * pc_x[k] * sgi_146[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pc_x, pc_y, pc_z, sfi_59, sfi_149, sgh0_114, \
                         sgh1_114, sgi_143, sgi_145, sgi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_14 * sfi_59[k]
                   + f_3 * pc_z[k] * sgi_143[k];

        t_188[k] = f_3 * pc_y[k] * sgi_145[k];

        t_189[k] = f_14 * sfi_149[k]
                   + f_6 * sgh0_114[k]
                   - f_7 * sgh1_114[k]
                   + f_3 * pc_x[k] * sgi_149[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_z, sfi_62, sfi_150, sfi_152, sgh0_115, \
                         sgh0_117, sgh1_115, sgh1_117, sgi_146, sgi_150, \
                         sgi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_14 * sfi_150[k]
                   + f_8 * sgh0_115[k]
                   - f_9 * sgh1_115[k]
                   + f_3 * pc_x[k] * sgi_150[k];

        t_191[k] = f_14 * sfi_62[k]
                   + f_3 * pc_z[k] * sgi_146[k];

        t_192[k] = f_14 * sfi_152[k]
                   + f_8 * sgh0_117[k]
                   - f_9 * sgh1_117[k]
                   + f_3 * pc_x[k] * sgi_152[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, sfi_154, sfi_155, sgh0_119, \
                         sgh0_120, sgh1_119, sgh1_120, sgi_149, sgi_154, \
                         sgi_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * sgi_149[k];

        t_194[k] = f_14 * sfi_154[k]
                   + f_8 * sgh0_119[k]
                   - f_9 * sgh1_119[k]
                   + f_3 * pc_x[k] * sgi_154[k];

        t_195[k] = f_14 * sfi_155[k]
                   + f_10 * sgh0_120[k]
                   - f_11 * sgh1_120[k]
                   + f_3 * pc_x[k] * sgi_155[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pc_x, pc_z, sfi_66, sfi_157, sfi_158, sgh0_122, \
                         sgh0_123, sgh1_122, sgh1_123, sgi_150, sgi_157, \
                         sgi_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_14 * sfi_66[k]
                   + f_3 * pc_z[k] * sgi_150[k];

        t_197[k] = f_14 * sfi_157[k]
                   + f_10 * sgh0_122[k]
                   - f_11 * sgh1_122[k]
                   + f_3 * pc_x[k] * sgi_157[k];

        t_198[k] = f_14 * sfi_158[k]
                   + f_10 * sgh0_123[k]
                   - f_11 * sgh1_123[k]
                   + f_3 * pc_x[k] * sgi_158[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pc_x, pc_y, sfi_160, sfi_161, sfi_162, \
                         sgh0_125, sgh1_125, sgi_154, sgi_160, sgi_161, \
                         sgi_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_3 * pc_y[k] * sgi_154[k];

        t_200[k] = f_14 * sfi_160[k]
                   + f_10 * sgh0_125[k]
                   - f_11 * sgh1_125[k]
                   + f_3 * pc_x[k] * sgi_160[k];

        t_201[k] = f_14 * sfi_161[k]
                   + f_3 * pc_x[k] * sgi_161[k];

        t_202[k] = f_14 * sfi_162[k]
                   + f_3 * pc_x[k] * sgi_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pc_x, sfi_163, sfi_164, sfi_165, \
                         sfi_166, sfi_167, sgi_163, sgi_164, sgi_165, sgi_166, \
                         sgi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_14 * sfi_163[k]
                   + f_3 * pc_x[k] * sgi_163[k];

        t_204[k] = f_14 * sfi_164[k]
                   + f_3 * pc_x[k] * sgi_164[k];

        t_205[k] = f_14 * sfi_165[k]
                   + f_3 * pc_x[k] * sgi_165[k];

        t_206[k] = f_14 * sfi_166[k]
                   + f_3 * pc_x[k] * sgi_166[k];

        t_207[k] = f_14 * sfi_167[k]
                   + f_3 * pc_x[k] * sgi_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_y, pc_z, sfi_77, sgh0_120, sgh0_122, \
                         sgh0_123, sgh1_120, sgh1_122, sgh1_123, sgi_161, sgi_163, \
                         sgi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * sgh0_120[k]
                   - f_2 * sgh1_120[k]
                   + f_3 * pc_y[k] * sgi_161[k];

        t_209[k] = f_14 * sfi_77[k]
                   + f_3 * pc_z[k] * sgi_161[k];

        t_210[k] = f_4 * sgh0_122[k]
                   - f_5 * sgh1_122[k]
                   + f_3 * pc_y[k] * sgi_163[k];

        t_211[k] = f_6 * sgh0_123[k]
                   - f_7 * sgh1_123[k]
                   + f_3 * pc_y[k] * sgi_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_y, pc_z, sfi_83, sgh0_124, sgh0_125, \
                         sgh1_124, sgh1_125, sgi_165, sgi_166, \
                         sgi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * sgh0_124[k]
                   - f_9 * sgh1_124[k]
                   + f_3 * pc_y[k] * sgi_165[k];

        t_213[k] = f_10 * sgh0_125[k]
                   - f_11 * sgh1_125[k]
                   + f_3 * pc_y[k] * sgi_166[k];

        t_214[k] = f_3 * pc_y[k] * sgi_167[k];

        t_215[k] = f_14 * sfi_83[k]
                   + f_1 * sgh0_125[k]
                   - f_2 * sgh1_125[k]
                   + f_3 * pc_z[k] * sgi_167[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pb_x, pc_x, pc_y, pc_z, sfk0_216, \
                         sfk0_219, sfi_84, sfi_168, sfi_171, sfk1_216, sfk1_219, \
                         sgi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pb_x[k] * sfk0_216[k]
                   + f_17 * sfi_168[k]
                   - f_12 * pc_x[k] * sfk1_216[k];

        t_217[k] = f_15 * sfi_84[k]
                   + f_3 * pc_y[k] * sgi_168[k];

        t_218[k] = f_3 * pc_z[k] * sgi_168[k];

        t_219[k] = pb_x[k] * sfk0_219[k]
                   + f_16 * sfi_171[k]
                   - f_12 * pc_x[k] * sfk1_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pb_x, pc_x, pc_y, sfk0_221, sfk0_222, sfi_86, \
                         sfi_173, sfi_174, sfk1_221, sfk1_222, \
                         sgi_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * sfi_86[k]
                   + f_3 * pc_y[k] * sgi_170[k];

        t_221[k] = pb_x[k] * sfk0_221[k]
                   + f_16 * sfi_173[k]
                   - f_12 * pc_x[k] * sfk1_221[k];

        t_222[k] = pb_x[k] * sfk0_222[k]
                   + f_0 * sfi_174[k]
                   - f_12 * pc_x[k] * sfk1_222[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, pb_x, pc_x, pc_y, pc_z, sfk0_225, sfi_89, \
                         sfi_177, sfk1_225, sgi_171, sgi_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_z[k] * sgi_171[k];

        t_224[k] = f_15 * sfi_89[k]
                   + f_3 * pc_y[k] * sgi_173[k];

        t_225[k] = pb_x[k] * sfk0_225[k]
                   + f_0 * sfi_177[k]
                   - f_12 * pc_x[k] * sfk1_225[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pb_x, pc_x, pc_z, sfk0_226, sfk0_228, sfi_178, \
                         sfi_180, sfk1_226, sfk1_228, sgi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pb_x[k] * sfk0_226[k]
                   + f_15 * sfi_178[k]
                   - f_12 * pc_x[k] * sfk1_226[k];

        t_227[k] = f_3 * pc_z[k] * sgi_174[k];

        t_228[k] = pb_x[k] * sfk0_228[k]
                   + f_15 * sfi_180[k]
                   - f_12 * pc_x[k] * sfk1_228[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pb_x, pc_x, pc_y, sfk0_230, sfk0_231, sfi_93, \
                         sfi_182, sfi_183, sfk1_230, sfk1_231, \
                         sgi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_15 * sfi_93[k]
                   + f_3 * pc_y[k] * sgi_177[k];

        t_230[k] = pb_x[k] * sfk0_230[k]
                   + f_15 * sfi_182[k]
                   - f_12 * pc_x[k] * sfk1_230[k];

        t_231[k] = pb_x[k] * sfk0_231[k]
                   + f_14 * sfi_183[k]
                   - f_12 * pc_x[k] * sfk1_231[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pb_x, pc_x, pc_z, sfk0_233, sfk0_234, sfi_185, \
                         sfi_186, sfk1_233, sfk1_234, sgi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * pc_z[k] * sgi_178[k];

        t_233[k] = pb_x[k] * sfk0_233[k]
                   + f_14 * sfi_185[k]
                   - f_12 * pc_x[k] * sfk1_233[k];

        t_234[k] = pb_x[k] * sfk0_234[k]
                   + f_14 * sfi_186[k]
                   - f_12 * pc_x[k] * sfk1_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pb_x, pc_x, pc_y, sfk0_236, sfi_98, \
                         sfi_188, sfi_189, sfi_190, sfk1_236, sgi_182, sgi_189, \
                         sgi_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_15 * sfi_98[k]
                   + f_3 * pc_y[k] * sgi_182[k];

        t_236[k] = pb_x[k] * sfk0_236[k]
                   + f_14 * sfi_188[k]
                   - f_12 * pc_x[k] * sfk1_236[k];

        t_237[k] = f_13 * sfi_189[k]
                   + f_3 * pc_x[k] * sgi_189[k];

        t_238[k] = f_13 * sfi_190[k]
                   + f_3 * pc_x[k] * sgi_190[k];
    }
}

static auto
compute_prim_sgk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfk0,
                                                          const size_t sfi, const size_t sfk1,
                                                          const size_t sgh0, const size_t sgh1,
                                                          const size_t sgi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfk0_108 = buffer.data(sfk0 + 108);
    const auto *sfk0_111 = buffer.data(sfk0 + 111);
    const auto *sfk0_114 = buffer.data(sfk0 + 114);
    const auto *sfk0_118 = buffer.data(sfk0 + 118);
    const auto *sfk0_123 = buffer.data(sfk0 + 123);
    const auto *sfk0_180 = buffer.data(sfk0 + 180);
    const auto *sfk0_185 = buffer.data(sfk0 + 185);
    const auto *sfk0_189 = buffer.data(sfk0 + 189);
    const auto *sfk0_194 = buffer.data(sfk0 + 194);
    const auto *sfk0_200 = buffer.data(sfk0 + 200);
    const auto *sfk0_244 = buffer.data(sfk0 + 244);
    const auto *sfk0_246 = buffer.data(sfk0 + 246);
    const auto *sfk0_247 = buffer.data(sfk0 + 247);
    const auto *sfk0_248 = buffer.data(sfk0 + 248);
    const auto *sfk0_249 = buffer.data(sfk0 + 249);
    const auto *sfk0_251 = buffer.data(sfk0 + 251);
    const auto *sfk0_257 = buffer.data(sfk0 + 257);
    const auto *sfk0_261 = buffer.data(sfk0 + 261);
    const auto *sfk0_264 = buffer.data(sfk0 + 264);
    const auto *sfk0_266 = buffer.data(sfk0 + 266);
    const auto *sfk0_269 = buffer.data(sfk0 + 269);
    const auto *sfk0_270 = buffer.data(sfk0 + 270);
    const auto *sfk0_272 = buffer.data(sfk0 + 272);
    const auto *sfk0_280 = buffer.data(sfk0 + 280);
    const auto *sfk0_282 = buffer.data(sfk0 + 282);
    const auto *sfk0_283 = buffer.data(sfk0 + 283);
    const auto *sfk0_284 = buffer.data(sfk0 + 284);
    const auto *sfk0_285 = buffer.data(sfk0 + 285);
    const auto *sfk0_287 = buffer.data(sfk0 + 287);
    const auto *sfk0_291 = buffer.data(sfk0 + 291);
    const auto *sfk0_294 = buffer.data(sfk0 + 294);
    const auto *sfk0_298 = buffer.data(sfk0 + 298);
    const auto *sfk0_300 = buffer.data(sfk0 + 300);
    const auto *sfk0_303 = buffer.data(sfk0 + 303);
    const auto *sfk0_305 = buffer.data(sfk0 + 305);
    const auto *sfk0_306 = buffer.data(sfk0 + 306);
    const auto *sfk0_316 = buffer.data(sfk0 + 316);
    const auto *sfk0_318 = buffer.data(sfk0 + 318);
    const auto *sfk0_319 = buffer.data(sfk0 + 319);
    const auto *sfk0_320 = buffer.data(sfk0 + 320);
    const auto *sfk0_321 = buffer.data(sfk0 + 321);
    const auto *sfk0_323 = buffer.data(sfk0 + 323);
    const auto *sfk0_324 = buffer.data(sfk0 + 324);
    const auto *sfk0_327 = buffer.data(sfk0 + 327);
    const auto *sfk0_329 = buffer.data(sfk0 + 329);
    const auto *sfk0_330 = buffer.data(sfk0 + 330);
    const auto *sfk0_333 = buffer.data(sfk0 + 333);
    const auto *sfk0_334 = buffer.data(sfk0 + 334);
    const auto *sfk0_336 = buffer.data(sfk0 + 336);
    const auto *sfk0_338 = buffer.data(sfk0 + 338);
    const auto *sfk0_339 = buffer.data(sfk0 + 339);
    const auto *sfk0_341 = buffer.data(sfk0 + 341);
    const auto *sfk0_342 = buffer.data(sfk0 + 342);
    const auto *sfk0_344 = buffer.data(sfk0 + 344);
    const auto *sfk0_352 = buffer.data(sfk0 + 352);
    const auto *sfk0_354 = buffer.data(sfk0 + 354);
    const auto *sfk0_355 = buffer.data(sfk0 + 355);
    const auto *sfk0_356 = buffer.data(sfk0 + 356);
    const auto *sfk0_357 = buffer.data(sfk0 + 357);
    const auto *sfk0_359 = buffer.data(sfk0 + 359);

    const auto *sfi_84 = buffer.data(sfi + 84);
    const auto *sfi_87 = buffer.data(sfi + 87);
    const auto *sfi_90 = buffer.data(sfi + 90);
    const auto *sfi_94 = buffer.data(sfi + 94);
    const auto *sfi_105 = buffer.data(sfi + 105);
    const auto *sfi_111 = buffer.data(sfi + 111);
    const auto *sfi_112 = buffer.data(sfi + 112);
    const auto *sfi_114 = buffer.data(sfi + 114);
    const auto *sfi_115 = buffer.data(sfi + 115);
    const auto *sfi_117 = buffer.data(sfi + 117);
    const auto *sfi_118 = buffer.data(sfi + 118);
    const auto *sfi_121 = buffer.data(sfi + 121);
    const auto *sfi_122 = buffer.data(sfi + 122);
    const auto *sfi_126 = buffer.data(sfi + 126);
    const auto *sfi_133 = buffer.data(sfi + 133);
    const auto *sfi_139 = buffer.data(sfi + 139);
    const auto *sfi_140 = buffer.data(sfi + 140);
    const auto *sfi_142 = buffer.data(sfi + 142);
    const auto *sfi_143 = buffer.data(sfi + 143);
    const auto *sfi_145 = buffer.data(sfi + 145);
    const auto *sfi_146 = buffer.data(sfi + 146);
    const auto *sfi_149 = buffer.data(sfi + 149);
    const auto *sfi_150 = buffer.data(sfi + 150);
    const auto *sfi_154 = buffer.data(sfi + 154);
    const auto *sfi_161 = buffer.data(sfi + 161);
    const auto *sfi_167 = buffer.data(sfi + 167);
    const auto *sfi_168 = buffer.data(sfi + 168);
    const auto *sfi_191 = buffer.data(sfi + 191);
    const auto *sfi_192 = buffer.data(sfi + 192);
    const auto *sfi_193 = buffer.data(sfi + 193);
    const auto *sfi_194 = buffer.data(sfi + 194);
    const auto *sfi_195 = buffer.data(sfi + 195);
    const auto *sfi_201 = buffer.data(sfi + 201);
    const auto *sfi_205 = buffer.data(sfi + 205);
    const auto *sfi_208 = buffer.data(sfi + 208);
    const auto *sfi_210 = buffer.data(sfi + 210);
    const auto *sfi_213 = buffer.data(sfi + 213);
    const auto *sfi_214 = buffer.data(sfi + 214);
    const auto *sfi_216 = buffer.data(sfi + 216);
    const auto *sfi_217 = buffer.data(sfi + 217);
    const auto *sfi_218 = buffer.data(sfi + 218);
    const auto *sfi_219 = buffer.data(sfi + 219);
    const auto *sfi_220 = buffer.data(sfi + 220);
    const auto *sfi_221 = buffer.data(sfi + 221);
    const auto *sfi_222 = buffer.data(sfi + 222);
    const auto *sfi_223 = buffer.data(sfi + 223);
    const auto *sfi_227 = buffer.data(sfi + 227);
    const auto *sfi_230 = buffer.data(sfi + 230);
    const auto *sfi_234 = buffer.data(sfi + 234);
    const auto *sfi_236 = buffer.data(sfi + 236);
    const auto *sfi_239 = buffer.data(sfi + 239);
    const auto *sfi_241 = buffer.data(sfi + 241);
    const auto *sfi_242 = buffer.data(sfi + 242);
    const auto *sfi_245 = buffer.data(sfi + 245);
    const auto *sfi_246 = buffer.data(sfi + 246);
    const auto *sfi_247 = buffer.data(sfi + 247);
    const auto *sfi_248 = buffer.data(sfi + 248);
    const auto *sfi_249 = buffer.data(sfi + 249);
    const auto *sfi_250 = buffer.data(sfi + 250);
    const auto *sfi_251 = buffer.data(sfi + 251);
    const auto *sfi_252 = buffer.data(sfi + 252);
    const auto *sfi_255 = buffer.data(sfi + 255);
    const auto *sfi_257 = buffer.data(sfi + 257);
    const auto *sfi_258 = buffer.data(sfi + 258);
    const auto *sfi_261 = buffer.data(sfi + 261);
    const auto *sfi_262 = buffer.data(sfi + 262);
    const auto *sfi_264 = buffer.data(sfi + 264);
    const auto *sfi_266 = buffer.data(sfi + 266);
    const auto *sfi_267 = buffer.data(sfi + 267);
    const auto *sfi_269 = buffer.data(sfi + 269);
    const auto *sfi_270 = buffer.data(sfi + 270);
    const auto *sfi_272 = buffer.data(sfi + 272);
    const auto *sfi_273 = buffer.data(sfi + 273);
    const auto *sfi_274 = buffer.data(sfi + 274);
    const auto *sfi_275 = buffer.data(sfi + 275);
    const auto *sfi_276 = buffer.data(sfi + 276);
    const auto *sfi_277 = buffer.data(sfi + 277);
    const auto *sfi_278 = buffer.data(sfi + 278);
    const auto *sfi_279 = buffer.data(sfi + 279);

    const auto *sfk1_108 = buffer.data(sfk1 + 108);
    const auto *sfk1_111 = buffer.data(sfk1 + 111);
    const auto *sfk1_114 = buffer.data(sfk1 + 114);
    const auto *sfk1_118 = buffer.data(sfk1 + 118);
    const auto *sfk1_123 = buffer.data(sfk1 + 123);
    const auto *sfk1_180 = buffer.data(sfk1 + 180);
    const auto *sfk1_185 = buffer.data(sfk1 + 185);
    const auto *sfk1_189 = buffer.data(sfk1 + 189);
    const auto *sfk1_194 = buffer.data(sfk1 + 194);
    const auto *sfk1_200 = buffer.data(sfk1 + 200);
    const auto *sfk1_244 = buffer.data(sfk1 + 244);
    const auto *sfk1_246 = buffer.data(sfk1 + 246);
    const auto *sfk1_247 = buffer.data(sfk1 + 247);
    const auto *sfk1_248 = buffer.data(sfk1 + 248);
    const auto *sfk1_249 = buffer.data(sfk1 + 249);
    const auto *sfk1_251 = buffer.data(sfk1 + 251);
    const auto *sfk1_257 = buffer.data(sfk1 + 257);
    const auto *sfk1_261 = buffer.data(sfk1 + 261);
    const auto *sfk1_264 = buffer.data(sfk1 + 264);
    const auto *sfk1_266 = buffer.data(sfk1 + 266);
    const auto *sfk1_269 = buffer.data(sfk1 + 269);
    const auto *sfk1_270 = buffer.data(sfk1 + 270);
    const auto *sfk1_272 = buffer.data(sfk1 + 272);
    const auto *sfk1_280 = buffer.data(sfk1 + 280);
    const auto *sfk1_282 = buffer.data(sfk1 + 282);
    const auto *sfk1_283 = buffer.data(sfk1 + 283);
    const auto *sfk1_284 = buffer.data(sfk1 + 284);
    const auto *sfk1_285 = buffer.data(sfk1 + 285);
    const auto *sfk1_287 = buffer.data(sfk1 + 287);
    const auto *sfk1_291 = buffer.data(sfk1 + 291);
    const auto *sfk1_294 = buffer.data(sfk1 + 294);
    const auto *sfk1_298 = buffer.data(sfk1 + 298);
    const auto *sfk1_300 = buffer.data(sfk1 + 300);
    const auto *sfk1_303 = buffer.data(sfk1 + 303);
    const auto *sfk1_305 = buffer.data(sfk1 + 305);
    const auto *sfk1_306 = buffer.data(sfk1 + 306);
    const auto *sfk1_316 = buffer.data(sfk1 + 316);
    const auto *sfk1_318 = buffer.data(sfk1 + 318);
    const auto *sfk1_319 = buffer.data(sfk1 + 319);
    const auto *sfk1_320 = buffer.data(sfk1 + 320);
    const auto *sfk1_321 = buffer.data(sfk1 + 321);
    const auto *sfk1_323 = buffer.data(sfk1 + 323);
    const auto *sfk1_324 = buffer.data(sfk1 + 324);
    const auto *sfk1_327 = buffer.data(sfk1 + 327);
    const auto *sfk1_329 = buffer.data(sfk1 + 329);
    const auto *sfk1_330 = buffer.data(sfk1 + 330);
    const auto *sfk1_333 = buffer.data(sfk1 + 333);
    const auto *sfk1_334 = buffer.data(sfk1 + 334);
    const auto *sfk1_336 = buffer.data(sfk1 + 336);
    const auto *sfk1_338 = buffer.data(sfk1 + 338);
    const auto *sfk1_339 = buffer.data(sfk1 + 339);
    const auto *sfk1_341 = buffer.data(sfk1 + 341);
    const auto *sfk1_342 = buffer.data(sfk1 + 342);
    const auto *sfk1_344 = buffer.data(sfk1 + 344);
    const auto *sfk1_352 = buffer.data(sfk1 + 352);
    const auto *sfk1_354 = buffer.data(sfk1 + 354);
    const auto *sfk1_355 = buffer.data(sfk1 + 355);
    const auto *sfk1_356 = buffer.data(sfk1 + 356);
    const auto *sfk1_357 = buffer.data(sfk1 + 357);
    const auto *sfk1_359 = buffer.data(sfk1 + 359);

    const auto *sgh0_210 = buffer.data(sgh0 + 210);
    const auto *sgh0_213 = buffer.data(sgh0 + 213);

    const auto *sgh1_210 = buffer.data(sgh1 + 210);
    const auto *sgh1_213 = buffer.data(sgh1 + 213);

    const auto *sgi_189 = buffer.data(sgi + 189);
    const auto *sgi_191 = buffer.data(sgi + 191);
    const auto *sgi_192 = buffer.data(sgi + 192);
    const auto *sgi_193 = buffer.data(sgi + 193);
    const auto *sgi_194 = buffer.data(sgi + 194);
    const auto *sgi_195 = buffer.data(sgi + 195);
    const auto *sgi_196 = buffer.data(sgi + 196);
    const auto *sgi_198 = buffer.data(sgi + 198);
    const auto *sgi_199 = buffer.data(sgi + 199);
    const auto *sgi_201 = buffer.data(sgi + 201);
    const auto *sgi_202 = buffer.data(sgi + 202);
    const auto *sgi_205 = buffer.data(sgi + 205);
    const auto *sgi_206 = buffer.data(sgi + 206);
    const auto *sgi_210 = buffer.data(sgi + 210);
    const auto *sgi_217 = buffer.data(sgi + 217);
    const auto *sgi_218 = buffer.data(sgi + 218);
    const auto *sgi_219 = buffer.data(sgi + 219);
    const auto *sgi_220 = buffer.data(sgi + 220);
    const auto *sgi_221 = buffer.data(sgi + 221);
    const auto *sgi_222 = buffer.data(sgi + 222);
    const auto *sgi_223 = buffer.data(sgi + 223);
    const auto *sgi_224 = buffer.data(sgi + 224);
    const auto *sgi_226 = buffer.data(sgi + 226);
    const auto *sgi_227 = buffer.data(sgi + 227);
    const auto *sgi_229 = buffer.data(sgi + 229);
    const auto *sgi_230 = buffer.data(sgi + 230);
    const auto *sgi_233 = buffer.data(sgi + 233);
    const auto *sgi_234 = buffer.data(sgi + 234);
    const auto *sgi_238 = buffer.data(sgi + 238);
    const auto *sgi_245 = buffer.data(sgi + 245);
    const auto *sgi_246 = buffer.data(sgi + 246);
    const auto *sgi_247 = buffer.data(sgi + 247);
    const auto *sgi_248 = buffer.data(sgi + 248);
    const auto *sgi_249 = buffer.data(sgi + 249);
    const auto *sgi_250 = buffer.data(sgi + 250);
    const auto *sgi_251 = buffer.data(sgi + 251);
    const auto *sgi_252 = buffer.data(sgi + 252);
    const auto *sgi_254 = buffer.data(sgi + 254);
    const auto *sgi_255 = buffer.data(sgi + 255);
    const auto *sgi_257 = buffer.data(sgi + 257);
    const auto *sgi_258 = buffer.data(sgi + 258);
    const auto *sgi_261 = buffer.data(sgi + 261);
    const auto *sgi_262 = buffer.data(sgi + 262);
    const auto *sgi_266 = buffer.data(sgi + 266);
    const auto *sgi_273 = buffer.data(sgi + 273);
    const auto *sgi_274 = buffer.data(sgi + 274);
    const auto *sgi_275 = buffer.data(sgi + 275);
    const auto *sgi_276 = buffer.data(sgi + 276);
    const auto *sgi_277 = buffer.data(sgi + 277);
    const auto *sgi_278 = buffer.data(sgi + 278);
    const auto *sgi_279 = buffer.data(sgi + 279);
    const auto *sgi_280 = buffer.data(sgi + 280);
    const auto *sgi_283 = buffer.data(sgi + 283);

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pc_x, sfi_191, sfi_192, sfi_193, \
                         sfi_194, sfi_195, sgi_191, sgi_192, sgi_193, sgi_194, \
                         sgi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_13 * sfi_191[k]
                   + f_3 * pc_x[k] * sgi_191[k];

        t_240[k] = f_13 * sfi_192[k]
                   + f_3 * pc_x[k] * sgi_192[k];

        t_241[k] = f_13 * sfi_193[k]
                   + f_3 * pc_x[k] * sgi_193[k];

        t_242[k] = f_13 * sfi_194[k]
                   + f_3 * pc_x[k] * sgi_194[k];

        t_243[k] = f_13 * sfi_195[k]
                   + f_3 * pc_x[k] * sgi_195[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pc_x, pc_z, sfk0_244, sfk0_246, \
                         sfk0_247, sfk1_244, sfk1_246, sfk1_247, \
                         sgi_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = pb_x[k] * sfk0_244[k]
                   - f_12 * pc_x[k] * sfk1_244[k];

        t_245[k] = f_3 * pc_z[k] * sgi_189[k];

        t_246[k] = pb_x[k] * sfk0_246[k]
                   - f_12 * pc_x[k] * sfk1_246[k];

        t_247[k] = pb_x[k] * sfk0_247[k]
                   - f_12 * pc_x[k] * sfk1_247[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pb_x, pc_x, pc_y, sfk0_248, sfk0_249, \
                         sfk0_251, sfi_111, sfk1_248, sfk1_249, sfk1_251, \
                         sgi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pb_x[k] * sfk0_248[k]
                   - f_12 * pc_x[k] * sfk1_248[k];

        t_249[k] = pb_x[k] * sfk0_249[k]
                   - f_12 * pc_x[k] * sfk1_249[k];

        t_250[k] = f_15 * sfi_111[k]
                   + f_3 * pc_y[k] * sgi_195[k];

        t_251[k] = pb_x[k] * sfk0_251[k]
                   - f_12 * pc_x[k] * sfk1_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_z, pc_y, pc_z, sfk0_108, sfk0_111, \
                         sfi_84, sfi_112, sfk1_108, sfk1_111, sgi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pb_z[k] * sfk0_108[k]
                   - f_12 * pc_z[k] * sfk1_108[k];

        t_253[k] = f_14 * sfi_112[k]
                   + f_3 * pc_y[k] * sgi_196[k];

        t_254[k] = f_13 * sfi_84[k]
                   + f_3 * pc_z[k] * sgi_196[k];

        t_255[k] = pb_z[k] * sfk0_111[k]
                   - f_12 * pc_z[k] * sfk1_111[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pb_x, pb_z, pc_x, pc_y, pc_z, sfk0_114, \
                         sfk0_257, sfi_114, sfi_201, sfk1_114, sfk1_257, \
                         sgi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * sfi_114[k]
                   + f_3 * pc_y[k] * sgi_198[k];

        t_257[k] = pb_x[k] * sfk0_257[k]
                   + f_16 * sfi_201[k]
                   - f_12 * pc_x[k] * sfk1_257[k];

        t_258[k] = pb_z[k] * sfk0_114[k]
                   - f_12 * pc_z[k] * sfk1_114[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, pb_x, pc_x, pc_y, pc_z, sfk0_261, sfi_87, \
                         sfi_117, sfi_205, sfk1_261, sgi_199, sgi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_13 * sfi_87[k]
                   + f_3 * pc_z[k] * sgi_199[k];

        t_260[k] = f_14 * sfi_117[k]
                   + f_3 * pc_y[k] * sgi_201[k];

        t_261[k] = pb_x[k] * sfk0_261[k]
                   + f_0 * sfi_205[k]
                   - f_12 * pc_x[k] * sfk1_261[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, pb_x, pb_z, pc_x, pc_z, sfk0_118, sfk0_264, \
                         sfi_90, sfi_208, sfk1_118, sfk1_264, sgi_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = pb_z[k] * sfk0_118[k]
                   - f_12 * pc_z[k] * sfk1_118[k];

        t_263[k] = f_13 * sfi_90[k]
                   + f_3 * pc_z[k] * sgi_202[k];

        t_264[k] = pb_x[k] * sfk0_264[k]
                   + f_15 * sfi_208[k]
                   - f_12 * pc_x[k] * sfk1_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pb_x, pb_z, pc_x, pc_y, pc_z, sfk0_123, \
                         sfk0_266, sfi_121, sfi_210, sfk1_123, sfk1_266, \
                         sgi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_14 * sfi_121[k]
                   + f_3 * pc_y[k] * sgi_205[k];

        t_266[k] = pb_x[k] * sfk0_266[k]
                   + f_15 * sfi_210[k]
                   - f_12 * pc_x[k] * sfk1_266[k];

        t_267[k] = pb_z[k] * sfk0_123[k]
                   - f_12 * pc_z[k] * sfk1_123[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pb_x, pc_x, pc_z, sfk0_269, sfk0_270, sfi_94, \
                         sfi_213, sfi_214, sfk1_269, sfk1_270, \
                         sgi_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_13 * sfi_94[k]
                   + f_3 * pc_z[k] * sgi_206[k];

        t_269[k] = pb_x[k] * sfk0_269[k]
                   + f_14 * sfi_213[k]
                   - f_12 * pc_x[k] * sfk1_269[k];

        t_270[k] = pb_x[k] * sfk0_270[k]
                   + f_14 * sfi_214[k]
                   - f_12 * pc_x[k] * sfk1_270[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pb_x, pc_x, pc_y, sfk0_272, sfi_126, \
                         sfi_216, sfi_217, sfi_218, sfk1_272, sgi_210, sgi_217, \
                         sgi_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_14 * sfi_126[k]
                   + f_3 * pc_y[k] * sgi_210[k];

        t_272[k] = pb_x[k] * sfk0_272[k]
                   + f_14 * sfi_216[k]
                   - f_12 * pc_x[k] * sfk1_272[k];

        t_273[k] = f_13 * sfi_217[k]
                   + f_3 * pc_x[k] * sgi_217[k];

        t_274[k] = f_13 * sfi_218[k]
                   + f_3 * pc_x[k] * sgi_218[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, pc_x, sfi_219, sfi_220, sfi_221, \
                         sfi_222, sfi_223, sgi_219, sgi_220, sgi_221, sgi_222, \
                         sgi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_13 * sfi_219[k]
                   + f_3 * pc_x[k] * sgi_219[k];

        t_276[k] = f_13 * sfi_220[k]
                   + f_3 * pc_x[k] * sgi_220[k];

        t_277[k] = f_13 * sfi_221[k]
                   + f_3 * pc_x[k] * sgi_221[k];

        t_278[k] = f_13 * sfi_222[k]
                   + f_3 * pc_x[k] * sgi_222[k];

        t_279[k] = f_13 * sfi_223[k]
                   + f_3 * pc_x[k] * sgi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pb_x, pc_x, pc_z, sfk0_280, sfk0_282, \
                         sfk0_283, sfi_105, sfk1_280, sfk1_282, sfk1_283, \
                         sgi_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pb_x[k] * sfk0_280[k]
                   - f_12 * pc_x[k] * sfk1_280[k];

        t_281[k] = f_13 * sfi_105[k]
                   + f_3 * pc_z[k] * sgi_217[k];

        t_282[k] = pb_x[k] * sfk0_282[k]
                   - f_12 * pc_x[k] * sfk1_282[k];

        t_283[k] = pb_x[k] * sfk0_283[k]
                   - f_12 * pc_x[k] * sfk1_283[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pb_x, pc_x, pc_y, sfk0_284, sfk0_285, \
                         sfk0_287, sfi_139, sfk1_284, sfk1_285, sfk1_287, \
                         sgi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = pb_x[k] * sfk0_284[k]
                   - f_12 * pc_x[k] * sfk1_284[k];

        t_285[k] = pb_x[k] * sfk0_285[k]
                   - f_12 * pc_x[k] * sfk1_285[k];

        t_286[k] = f_14 * sfi_139[k]
                   + f_3 * pc_y[k] * sgi_223[k];

        t_287[k] = pb_x[k] * sfk0_287[k]
                   - f_12 * pc_x[k] * sfk1_287[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pb_y, pc_y, pc_z, sfk0_180, sfi_112, sfi_140, \
                         sfk1_180, sgi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = pb_y[k] * sfk0_180[k]
                   - f_12 * pc_y[k] * sfk1_180[k];

        t_289[k] = f_13 * sfi_140[k]
                   + f_3 * pc_y[k] * sgi_224[k];

        t_290[k] = f_14 * sfi_112[k]
                   + f_3 * pc_z[k] * sgi_224[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_x, pb_y, pc_x, pc_y, sfk0_185, sfk0_291, \
                         sfi_142, sfi_227, sfk1_185, sfk1_291, \
                         sgi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = pb_x[k] * sfk0_291[k]
                   + f_16 * sfi_227[k]
                   - f_12 * pc_x[k] * sfk1_291[k];

        t_292[k] = f_13 * sfi_142[k]
                   + f_3 * pc_y[k] * sgi_226[k];

        t_293[k] = pb_y[k] * sfk0_185[k]
                   - f_12 * pc_y[k] * sfk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pb_x, pc_x, pc_y, pc_z, sfk0_294, sfi_115, \
                         sfi_145, sfi_230, sfk1_294, sgi_227, sgi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pb_x[k] * sfk0_294[k]
                   + f_0 * sfi_230[k]
                   - f_12 * pc_x[k] * sfk1_294[k];

        t_295[k] = f_14 * sfi_115[k]
                   + f_3 * pc_z[k] * sgi_227[k];

        t_296[k] = f_13 * sfi_145[k]
                   + f_3 * pc_y[k] * sgi_229[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_x, pb_y, pc_x, pc_y, pc_z, sfk0_189, \
                         sfk0_298, sfi_118, sfi_234, sfk1_189, sfk1_298, \
                         sgi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = pb_y[k] * sfk0_189[k]
                   - f_12 * pc_y[k] * sfk1_189[k];

        t_298[k] = pb_x[k] * sfk0_298[k]
                   + f_15 * sfi_234[k]
                   - f_12 * pc_x[k] * sfk1_298[k];

        t_299[k] = f_14 * sfi_118[k]
                   + f_3 * pc_z[k] * sgi_230[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pb_x, pb_y, pc_x, pc_y, sfk0_194, sfk0_300, \
                         sfi_149, sfi_236, sfk1_194, sfk1_300, \
                         sgi_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = pb_x[k] * sfk0_300[k]
                   + f_15 * sfi_236[k]
                   - f_12 * pc_x[k] * sfk1_300[k];

        t_301[k] = f_13 * sfi_149[k]
                   + f_3 * pc_y[k] * sgi_233[k];

        t_302[k] = pb_y[k] * sfk0_194[k]
                   - f_12 * pc_y[k] * sfk1_194[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, pb_x, pc_x, pc_z, sfk0_303, sfk0_305, sfi_122, \
                         sfi_239, sfi_241, sfk1_303, sfk1_305, \
                         sgi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pb_x[k] * sfk0_303[k]
                   + f_14 * sfi_239[k]
                   - f_12 * pc_x[k] * sfk1_303[k];

        t_304[k] = f_14 * sfi_122[k]
                   + f_3 * pc_z[k] * sgi_234[k];

        t_305[k] = pb_x[k] * sfk0_305[k]
                   + f_14 * sfi_241[k]
                   - f_12 * pc_x[k] * sfk1_305[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, pb_x, pb_y, pc_x, pc_y, sfk0_200, sfk0_306, \
                         sfi_154, sfi_242, sfk1_200, sfk1_306, \
                         sgi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = pb_x[k] * sfk0_306[k]
                   + f_14 * sfi_242[k]
                   - f_12 * pc_x[k] * sfk1_306[k];

        t_307[k] = f_13 * sfi_154[k]
                   + f_3 * pc_y[k] * sgi_238[k];

        t_308[k] = pb_y[k] * sfk0_200[k]
                   - f_12 * pc_y[k] * sfk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, sfi_245, sfi_246, sfi_247, \
                         sfi_248, sfi_249, sgi_245, sgi_246, sgi_247, sgi_248, \
                         sgi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_13 * sfi_245[k]
                   + f_3 * pc_x[k] * sgi_245[k];

        t_310[k] = f_13 * sfi_246[k]
                   + f_3 * pc_x[k] * sgi_246[k];

        t_311[k] = f_13 * sfi_247[k]
                   + f_3 * pc_x[k] * sgi_247[k];

        t_312[k] = f_13 * sfi_248[k]
                   + f_3 * pc_x[k] * sgi_248[k];

        t_313[k] = f_13 * sfi_249[k]
                   + f_3 * pc_x[k] * sgi_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pb_x, pc_x, pc_z, sfk0_316, sfi_133, \
                         sfi_250, sfi_251, sfk1_316, sgi_245, sgi_250, \
                         sgi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_13 * sfi_250[k]
                   + f_3 * pc_x[k] * sgi_250[k];

        t_315[k] = f_13 * sfi_251[k]
                   + f_3 * pc_x[k] * sgi_251[k];

        t_316[k] = pb_x[k] * sfk0_316[k]
                   - f_12 * pc_x[k] * sfk1_316[k];

        t_317[k] = f_14 * sfi_133[k]
                   + f_3 * pc_z[k] * sgi_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pb_x, pc_x, sfk0_318, sfk0_319, sfk0_320, \
                         sfk0_321, sfk1_318, sfk1_319, sfk1_320, \
                         sfk1_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_x[k] * sfk0_318[k]
                   - f_12 * pc_x[k] * sfk1_318[k];

        t_319[k] = pb_x[k] * sfk0_319[k]
                   - f_12 * pc_x[k] * sfk1_319[k];

        t_320[k] = pb_x[k] * sfk0_320[k]
                   - f_12 * pc_x[k] * sfk1_320[k];

        t_321[k] = pb_x[k] * sfk0_321[k]
                   - f_12 * pc_x[k] * sfk1_321[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pb_x, pc_x, pc_y, sfk0_323, sfk0_324, \
                         sfi_167, sfi_252, sfk1_323, sfk1_324, sgi_251, \
                         sgi_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_13 * sfi_167[k]
                   + f_3 * pc_y[k] * sgi_251[k];

        t_323[k] = pb_x[k] * sfk0_323[k]
                   - f_12 * pc_x[k] * sfk1_323[k];

        t_324[k] = pb_x[k] * sfk0_324[k]
                   + f_17 * sfi_252[k]
                   - f_12 * pc_x[k] * sfk1_324[k];

        t_325[k] = f_3 * pc_y[k] * sgi_252[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, pb_x, pc_x, pc_y, pc_z, sfk0_327, sfi_140, \
                         sfi_255, sfk1_327, sgi_252, sgi_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_15 * sfi_140[k]
                   + f_3 * pc_z[k] * sgi_252[k];

        t_327[k] = pb_x[k] * sfk0_327[k]
                   + f_16 * sfi_255[k]
                   - f_12 * pc_x[k] * sfk1_327[k];

        t_328[k] = f_3 * pc_y[k] * sgi_254[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pb_x, pc_x, pc_z, sfk0_329, sfk0_330, sfi_143, \
                         sfi_257, sfi_258, sfk1_329, sfk1_330, \
                         sgi_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = pb_x[k] * sfk0_329[k]
                   + f_16 * sfi_257[k]
                   - f_12 * pc_x[k] * sfk1_329[k];

        t_330[k] = pb_x[k] * sfk0_330[k]
                   + f_0 * sfi_258[k]
                   - f_12 * pc_x[k] * sfk1_330[k];

        t_331[k] = f_15 * sfi_143[k]
                   + f_3 * pc_z[k] * sgi_255[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pb_x, pc_x, pc_y, sfk0_333, sfk0_334, sfi_261, \
                         sfi_262, sfk1_333, sfk1_334, sgi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_3 * pc_y[k] * sgi_257[k];

        t_333[k] = pb_x[k] * sfk0_333[k]
                   + f_0 * sfi_261[k]
                   - f_12 * pc_x[k] * sfk1_333[k];

        t_334[k] = pb_x[k] * sfk0_334[k]
                   + f_15 * sfi_262[k]
                   - f_12 * pc_x[k] * sfk1_334[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pb_x, pc_x, pc_y, pc_z, sfk0_336, sfi_146, \
                         sfi_264, sfk1_336, sgi_258, sgi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_15 * sfi_146[k]
                   + f_3 * pc_z[k] * sgi_258[k];

        t_336[k] = pb_x[k] * sfk0_336[k]
                   + f_15 * sfi_264[k]
                   - f_12 * pc_x[k] * sfk1_336[k];

        t_337[k] = f_3 * pc_y[k] * sgi_261[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_x, pc_x, pc_z, sfk0_338, sfk0_339, sfi_150, \
                         sfi_266, sfi_267, sfk1_338, sfk1_339, \
                         sgi_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pb_x[k] * sfk0_338[k]
                   + f_15 * sfi_266[k]
                   - f_12 * pc_x[k] * sfk1_338[k];

        t_339[k] = pb_x[k] * sfk0_339[k]
                   + f_14 * sfi_267[k]
                   - f_12 * pc_x[k] * sfk1_339[k];

        t_340[k] = f_15 * sfi_150[k]
                   + f_3 * pc_z[k] * sgi_262[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pb_x, pc_x, pc_y, sfk0_341, sfk0_342, sfi_269, \
                         sfi_270, sfk1_341, sfk1_342, sgi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = pb_x[k] * sfk0_341[k]
                   + f_14 * sfi_269[k]
                   - f_12 * pc_x[k] * sfk1_341[k];

        t_342[k] = pb_x[k] * sfk0_342[k]
                   + f_14 * sfi_270[k]
                   - f_12 * pc_x[k] * sfk1_342[k];

        t_343[k] = f_3 * pc_y[k] * sgi_266[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pb_x, pc_x, sfk0_344, sfi_272, sfi_273, \
                         sfi_274, sfi_275, sfk1_344, sgi_273, sgi_274, \
                         sgi_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = pb_x[k] * sfk0_344[k]
                   + f_14 * sfi_272[k]
                   - f_12 * pc_x[k] * sfk1_344[k];

        t_345[k] = f_13 * sfi_273[k]
                   + f_3 * pc_x[k] * sgi_273[k];

        t_346[k] = f_13 * sfi_274[k]
                   + f_3 * pc_x[k] * sgi_274[k];

        t_347[k] = f_13 * sfi_275[k]
                   + f_3 * pc_x[k] * sgi_275[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pc_x, sfi_276, sfi_277, sfi_278, sfi_279, \
                         sgi_276, sgi_277, sgi_278, sgi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_13 * sfi_276[k]
                   + f_3 * pc_x[k] * sgi_276[k];

        t_349[k] = f_13 * sfi_277[k]
                   + f_3 * pc_x[k] * sgi_277[k];

        t_350[k] = f_13 * sfi_278[k]
                   + f_3 * pc_x[k] * sgi_278[k];

        t_351[k] = f_13 * sfi_279[k]
                   + f_3 * pc_x[k] * sgi_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pb_x, pc_x, pc_z, sfk0_352, sfk0_354, \
                         sfk0_355, sfi_161, sfk1_352, sfk1_354, sfk1_355, \
                         sgi_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = pb_x[k] * sfk0_352[k]
                   - f_12 * pc_x[k] * sfk1_352[k];

        t_353[k] = f_15 * sfi_161[k]
                   + f_3 * pc_z[k] * sgi_273[k];

        t_354[k] = pb_x[k] * sfk0_354[k]
                   - f_12 * pc_x[k] * sfk1_354[k];

        t_355[k] = pb_x[k] * sfk0_355[k]
                   - f_12 * pc_x[k] * sfk1_355[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pb_x, pc_x, pc_y, sfk0_356, sfk0_357, \
                         sfk0_359, sfk1_356, sfk1_357, sfk1_359, \
                         sgi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pb_x[k] * sfk0_356[k]
                   - f_12 * pc_x[k] * sfk1_356[k];

        t_357[k] = pb_x[k] * sfk0_357[k]
                   - f_12 * pc_x[k] * sfk1_357[k];

        t_358[k] = f_3 * pc_y[k] * sgi_279[k];

        t_359[k] = pb_x[k] * sfk0_359[k]
                   - f_12 * pc_x[k] * sfk1_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pc_x, pc_y, pc_z, sfi_168, sgh0_210, \
                         sgh0_213, sgh1_210, sgh1_213, sgi_280, \
                         sgi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_1 * sgh0_210[k]
                   - f_2 * sgh1_210[k]
                   + f_3 * pc_x[k] * sgi_280[k];

        t_361[k] = f_0 * sfi_168[k]
                   + f_3 * pc_y[k] * sgi_280[k];

        t_362[k] = f_3 * pc_z[k] * sgi_280[k];

        t_363[k] = f_4 * sgh0_213[k]
                   - f_5 * sgh1_213[k]
                   + f_3 * pc_x[k] * sgi_283[k];
    }
}

static auto
compute_prim_sgk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfk0,
                                                          const size_t sfi, const size_t sfk1,
                                                          const size_t sgh0, const size_t sgh1,
                                                          const size_t sgi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfk0_216 = buffer.data(sfk0 + 216);
    const auto *sfk0_219 = buffer.data(sfk0 + 219);
    const auto *sfk0_222 = buffer.data(sfk0 + 222);
    const auto *sfk0_226 = buffer.data(sfk0 + 226);
    const auto *sfk0_231 = buffer.data(sfk0 + 231);
    const auto *sfk0_244 = buffer.data(sfk0 + 244);
    const auto *sfk0_246 = buffer.data(sfk0 + 246);
    const auto *sfk0_247 = buffer.data(sfk0 + 247);
    const auto *sfk0_248 = buffer.data(sfk0 + 248);
    const auto *sfk0_249 = buffer.data(sfk0 + 249);
    const auto *sfk0_324 = buffer.data(sfk0 + 324);
    const auto *sfk0_329 = buffer.data(sfk0 + 329);
    const auto *sfk0_333 = buffer.data(sfk0 + 333);
    const auto *sfk0_338 = buffer.data(sfk0 + 338);

    const auto *sfi_168 = buffer.data(sfi + 168);
    const auto *sfi_170 = buffer.data(sfi + 170);
    const auto *sfi_171 = buffer.data(sfi + 171);
    const auto *sfi_173 = buffer.data(sfi + 173);
    const auto *sfi_174 = buffer.data(sfi + 174);
    const auto *sfi_177 = buffer.data(sfi + 177);
    const auto *sfi_178 = buffer.data(sfi + 178);
    const auto *sfi_182 = buffer.data(sfi + 182);
    const auto *sfi_189 = buffer.data(sfi + 189);
    const auto *sfi_190 = buffer.data(sfi + 190);
    const auto *sfi_191 = buffer.data(sfi + 191);
    const auto *sfi_192 = buffer.data(sfi + 192);
    const auto *sfi_193 = buffer.data(sfi + 193);
    const auto *sfi_194 = buffer.data(sfi + 194);
    const auto *sfi_195 = buffer.data(sfi + 195);
    const auto *sfi_196 = buffer.data(sfi + 196);
    const auto *sfi_198 = buffer.data(sfi + 198);
    const auto *sfi_199 = buffer.data(sfi + 199);
    const auto *sfi_201 = buffer.data(sfi + 201);
    const auto *sfi_202 = buffer.data(sfi + 202);
    const auto *sfi_205 = buffer.data(sfi + 205);
    const auto *sfi_206 = buffer.data(sfi + 206);
    const auto *sfi_210 = buffer.data(sfi + 210);
    const auto *sfi_217 = buffer.data(sfi + 217);
    const auto *sfi_223 = buffer.data(sfi + 223);
    const auto *sfi_224 = buffer.data(sfi + 224);
    const auto *sfi_226 = buffer.data(sfi + 226);
    const auto *sfi_227 = buffer.data(sfi + 227);
    const auto *sfi_229 = buffer.data(sfi + 229);
    const auto *sfi_230 = buffer.data(sfi + 230);
    const auto *sfi_233 = buffer.data(sfi + 233);
    const auto *sfi_234 = buffer.data(sfi + 234);
    const auto *sfi_238 = buffer.data(sfi + 238);
    const auto *sfi_245 = buffer.data(sfi + 245);
    const auto *sfi_247 = buffer.data(sfi + 247);
    const auto *sfi_248 = buffer.data(sfi + 248);
    const auto *sfi_249 = buffer.data(sfi + 249);
    const auto *sfi_250 = buffer.data(sfi + 250);
    const auto *sfi_251 = buffer.data(sfi + 251);
    const auto *sfi_252 = buffer.data(sfi + 252);
    const auto *sfi_254 = buffer.data(sfi + 254);
    const auto *sfi_257 = buffer.data(sfi + 257);
    const auto *sfi_261 = buffer.data(sfi + 261);

    const auto *sfk1_216 = buffer.data(sfk1 + 216);
    const auto *sfk1_219 = buffer.data(sfk1 + 219);
    const auto *sfk1_222 = buffer.data(sfk1 + 222);
    const auto *sfk1_226 = buffer.data(sfk1 + 226);
    const auto *sfk1_231 = buffer.data(sfk1 + 231);
    const auto *sfk1_244 = buffer.data(sfk1 + 244);
    const auto *sfk1_246 = buffer.data(sfk1 + 246);
    const auto *sfk1_247 = buffer.data(sfk1 + 247);
    const auto *sfk1_248 = buffer.data(sfk1 + 248);
    const auto *sfk1_249 = buffer.data(sfk1 + 249);
    const auto *sfk1_324 = buffer.data(sfk1 + 324);
    const auto *sfk1_329 = buffer.data(sfk1 + 329);
    const auto *sfk1_333 = buffer.data(sfk1 + 333);
    const auto *sfk1_338 = buffer.data(sfk1 + 338);

    const auto *sgh0_215 = buffer.data(sgh0 + 215);
    const auto *sgh0_216 = buffer.data(sgh0 + 216);
    const auto *sgh0_219 = buffer.data(sgh0 + 219);
    const auto *sgh0_220 = buffer.data(sgh0 + 220);
    const auto *sgh0_222 = buffer.data(sgh0 + 222);
    const auto *sgh0_224 = buffer.data(sgh0 + 224);
    const auto *sgh0_225 = buffer.data(sgh0 + 225);
    const auto *sgh0_227 = buffer.data(sgh0 + 227);
    const auto *sgh0_228 = buffer.data(sgh0 + 228);
    const auto *sgh0_229 = buffer.data(sgh0 + 229);
    const auto *sgh0_230 = buffer.data(sgh0 + 230);
    const auto *sgh0_236 = buffer.data(sgh0 + 236);
    const auto *sgh0_240 = buffer.data(sgh0 + 240);
    const auto *sgh0_243 = buffer.data(sgh0 + 243);
    const auto *sgh0_245 = buffer.data(sgh0 + 245);
    const auto *sgh0_248 = buffer.data(sgh0 + 248);
    const auto *sgh0_249 = buffer.data(sgh0 + 249);
    const auto *sgh0_251 = buffer.data(sgh0 + 251);
    const auto *sgh0_252 = buffer.data(sgh0 + 252);
    const auto *sgh0_255 = buffer.data(sgh0 + 255);
    const auto *sgh0_257 = buffer.data(sgh0 + 257);
    const auto *sgh0_258 = buffer.data(sgh0 + 258);
    const auto *sgh0_261 = buffer.data(sgh0 + 261);
    const auto *sgh0_262 = buffer.data(sgh0 + 262);
    const auto *sgh0_264 = buffer.data(sgh0 + 264);
    const auto *sgh0_266 = buffer.data(sgh0 + 266);
    const auto *sgh0_267 = buffer.data(sgh0 + 267);
    const auto *sgh0_269 = buffer.data(sgh0 + 269);
    const auto *sgh0_270 = buffer.data(sgh0 + 270);
    const auto *sgh0_271 = buffer.data(sgh0 + 271);
    const auto *sgh0_272 = buffer.data(sgh0 + 272);
    const auto *sgh0_276 = buffer.data(sgh0 + 276);
    const auto *sgh0_279 = buffer.data(sgh0 + 279);
    const auto *sgh0_283 = buffer.data(sgh0 + 283);
    const auto *sgh0_285 = buffer.data(sgh0 + 285);
    const auto *sgh0_288 = buffer.data(sgh0 + 288);

    const auto *sgh1_215 = buffer.data(sgh1 + 215);
    const auto *sgh1_216 = buffer.data(sgh1 + 216);
    const auto *sgh1_219 = buffer.data(sgh1 + 219);
    const auto *sgh1_220 = buffer.data(sgh1 + 220);
    const auto *sgh1_222 = buffer.data(sgh1 + 222);
    const auto *sgh1_224 = buffer.data(sgh1 + 224);
    const auto *sgh1_225 = buffer.data(sgh1 + 225);
    const auto *sgh1_227 = buffer.data(sgh1 + 227);
    const auto *sgh1_228 = buffer.data(sgh1 + 228);
    const auto *sgh1_229 = buffer.data(sgh1 + 229);
    const auto *sgh1_230 = buffer.data(sgh1 + 230);
    const auto *sgh1_236 = buffer.data(sgh1 + 236);
    const auto *sgh1_240 = buffer.data(sgh1 + 240);
    const auto *sgh1_243 = buffer.data(sgh1 + 243);
    const auto *sgh1_245 = buffer.data(sgh1 + 245);
    const auto *sgh1_248 = buffer.data(sgh1 + 248);
    const auto *sgh1_249 = buffer.data(sgh1 + 249);
    const auto *sgh1_251 = buffer.data(sgh1 + 251);
    const auto *sgh1_252 = buffer.data(sgh1 + 252);
    const auto *sgh1_255 = buffer.data(sgh1 + 255);
    const auto *sgh1_257 = buffer.data(sgh1 + 257);
    const auto *sgh1_258 = buffer.data(sgh1 + 258);
    const auto *sgh1_261 = buffer.data(sgh1 + 261);
    const auto *sgh1_262 = buffer.data(sgh1 + 262);
    const auto *sgh1_264 = buffer.data(sgh1 + 264);
    const auto *sgh1_266 = buffer.data(sgh1 + 266);
    const auto *sgh1_267 = buffer.data(sgh1 + 267);
    const auto *sgh1_269 = buffer.data(sgh1 + 269);
    const auto *sgh1_270 = buffer.data(sgh1 + 270);
    const auto *sgh1_271 = buffer.data(sgh1 + 271);
    const auto *sgh1_272 = buffer.data(sgh1 + 272);
    const auto *sgh1_276 = buffer.data(sgh1 + 276);
    const auto *sgh1_279 = buffer.data(sgh1 + 279);
    const auto *sgh1_283 = buffer.data(sgh1 + 283);
    const auto *sgh1_285 = buffer.data(sgh1 + 285);
    const auto *sgh1_288 = buffer.data(sgh1 + 288);

    const auto *sgi_282 = buffer.data(sgi + 282);
    const auto *sgi_283 = buffer.data(sgi + 283);
    const auto *sgi_285 = buffer.data(sgi + 285);
    const auto *sgi_286 = buffer.data(sgi + 286);
    const auto *sgi_289 = buffer.data(sgi + 289);
    const auto *sgi_290 = buffer.data(sgi + 290);
    const auto *sgi_292 = buffer.data(sgi + 292);
    const auto *sgi_294 = buffer.data(sgi + 294);
    const auto *sgi_295 = buffer.data(sgi + 295);
    const auto *sgi_297 = buffer.data(sgi + 297);
    const auto *sgi_298 = buffer.data(sgi + 298);
    const auto *sgi_300 = buffer.data(sgi + 300);
    const auto *sgi_301 = buffer.data(sgi + 301);
    const auto *sgi_302 = buffer.data(sgi + 302);
    const auto *sgi_303 = buffer.data(sgi + 303);
    const auto *sgi_304 = buffer.data(sgi + 304);
    const auto *sgi_305 = buffer.data(sgi + 305);
    const auto *sgi_306 = buffer.data(sgi + 306);
    const auto *sgi_307 = buffer.data(sgi + 307);
    const auto *sgi_308 = buffer.data(sgi + 308);
    const auto *sgi_310 = buffer.data(sgi + 310);
    const auto *sgi_311 = buffer.data(sgi + 311);
    const auto *sgi_313 = buffer.data(sgi + 313);
    const auto *sgi_314 = buffer.data(sgi + 314);
    const auto *sgi_317 = buffer.data(sgi + 317);
    const auto *sgi_318 = buffer.data(sgi + 318);
    const auto *sgi_320 = buffer.data(sgi + 320);
    const auto *sgi_322 = buffer.data(sgi + 322);
    const auto *sgi_325 = buffer.data(sgi + 325);
    const auto *sgi_326 = buffer.data(sgi + 326);
    const auto *sgi_328 = buffer.data(sgi + 328);
    const auto *sgi_329 = buffer.data(sgi + 329);
    const auto *sgi_330 = buffer.data(sgi + 330);
    const auto *sgi_331 = buffer.data(sgi + 331);
    const auto *sgi_332 = buffer.data(sgi + 332);
    const auto *sgi_333 = buffer.data(sgi + 333);
    const auto *sgi_334 = buffer.data(sgi + 334);
    const auto *sgi_335 = buffer.data(sgi + 335);
    const auto *sgi_336 = buffer.data(sgi + 336);
    const auto *sgi_338 = buffer.data(sgi + 338);
    const auto *sgi_339 = buffer.data(sgi + 339);
    const auto *sgi_341 = buffer.data(sgi + 341);
    const auto *sgi_342 = buffer.data(sgi + 342);
    const auto *sgi_345 = buffer.data(sgi + 345);
    const auto *sgi_346 = buffer.data(sgi + 346);
    const auto *sgi_348 = buffer.data(sgi + 348);
    const auto *sgi_350 = buffer.data(sgi + 350);
    const auto *sgi_351 = buffer.data(sgi + 351);
    const auto *sgi_353 = buffer.data(sgi + 353);
    const auto *sgi_354 = buffer.data(sgi + 354);
    const auto *sgi_356 = buffer.data(sgi + 356);
    const auto *sgi_357 = buffer.data(sgi + 357);
    const auto *sgi_358 = buffer.data(sgi + 358);
    const auto *sgi_359 = buffer.data(sgi + 359);
    const auto *sgi_360 = buffer.data(sgi + 360);
    const auto *sgi_361 = buffer.data(sgi + 361);
    const auto *sgi_362 = buffer.data(sgi + 362);
    const auto *sgi_363 = buffer.data(sgi + 363);
    const auto *sgi_364 = buffer.data(sgi + 364);
    const auto *sgi_366 = buffer.data(sgi + 366);
    const auto *sgi_367 = buffer.data(sgi + 367);
    const auto *sgi_369 = buffer.data(sgi + 369);
    const auto *sgi_370 = buffer.data(sgi + 370);
    const auto *sgi_373 = buffer.data(sgi + 373);
    const auto *sgi_374 = buffer.data(sgi + 374);
    const auto *sgi_376 = buffer.data(sgi + 376);
    const auto *sgi_379 = buffer.data(sgi + 379);

#pragma omp simd aligned(t_364, t_365, t_366, t_367, pc_x, pc_y, pc_z, sfi_170, sgh0_215, \
                         sgh0_216, sgh1_215, sgh1_216, sgi_282, sgi_283, sgi_285, \
                         sgi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_0 * sfi_170[k]
                   + f_3 * pc_y[k] * sgi_282[k];

        t_365[k] = f_4 * sgh0_215[k]
                   - f_5 * sgh1_215[k]
                   + f_3 * pc_x[k] * sgi_285[k];

        t_366[k] = f_6 * sgh0_216[k]
                   - f_7 * sgh1_216[k]
                   + f_3 * pc_x[k] * sgi_286[k];

        t_367[k] = f_3 * pc_z[k] * sgi_283[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pc_x, pc_y, pc_z, sfi_173, sgh0_219, \
                         sgh0_220, sgh1_219, sgh1_220, sgi_285, sgi_286, sgi_289, \
                         sgi_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_0 * sfi_173[k]
                   + f_3 * pc_y[k] * sgi_285[k];

        t_369[k] = f_6 * sgh0_219[k]
                   - f_7 * sgh1_219[k]
                   + f_3 * pc_x[k] * sgi_289[k];

        t_370[k] = f_8 * sgh0_220[k]
                   - f_9 * sgh1_220[k]
                   + f_3 * pc_x[k] * sgi_290[k];

        t_371[k] = f_3 * pc_z[k] * sgi_286[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_x, pc_y, sfi_177, sgh0_222, sgh0_224, \
                         sgh1_222, sgh1_224, sgi_289, sgi_292, \
                         sgi_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_8 * sgh0_222[k]
                   - f_9 * sgh1_222[k]
                   + f_3 * pc_x[k] * sgi_292[k];

        t_373[k] = f_0 * sfi_177[k]
                   + f_3 * pc_y[k] * sgi_289[k];

        t_374[k] = f_8 * sgh0_224[k]
                   - f_9 * sgh1_224[k]
                   + f_3 * pc_x[k] * sgi_294[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, pc_x, pc_z, sgh0_225, sgh0_227, sgh0_228, \
                         sgh1_225, sgh1_227, sgh1_228, sgi_290, sgi_295, sgi_297, \
                         sgi_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_10 * sgh0_225[k]
                   - f_11 * sgh1_225[k]
                   + f_3 * pc_x[k] * sgi_295[k];

        t_376[k] = f_3 * pc_z[k] * sgi_290[k];

        t_377[k] = f_10 * sgh0_227[k]
                   - f_11 * sgh1_227[k]
                   + f_3 * pc_x[k] * sgi_297[k];

        t_378[k] = f_10 * sgh0_228[k]
                   - f_11 * sgh1_228[k]
                   + f_3 * pc_x[k] * sgi_298[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, pc_x, pc_y, sfi_182, sgh0_230, \
                         sgh1_230, sgi_294, sgi_300, sgi_301, sgi_302, \
                         sgi_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_0 * sfi_182[k]
                   + f_3 * pc_y[k] * sgi_294[k];

        t_380[k] = f_10 * sgh0_230[k]
                   - f_11 * sgh1_230[k]
                   + f_3 * pc_x[k] * sgi_300[k];

        t_381[k] = f_3 * pc_x[k] * sgi_301[k];

        t_382[k] = f_3 * pc_x[k] * sgi_302[k];

        t_383[k] = f_3 * pc_x[k] * sgi_303[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, t_388, pc_x, pc_y, sfi_189, sgh0_225, \
                         sgh1_225, sgi_301, sgi_304, sgi_305, sgi_306, \
                         sgi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_3 * pc_x[k] * sgi_304[k];

        t_385[k] = f_3 * pc_x[k] * sgi_305[k];

        t_386[k] = f_3 * pc_x[k] * sgi_306[k];

        t_387[k] = f_3 * pc_x[k] * sgi_307[k];

        t_388[k] = f_0 * sfi_189[k]
                   + f_1 * sgh0_225[k]
                   - f_2 * sgh1_225[k]
                   + f_3 * pc_y[k] * sgi_301[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, pc_y, pc_z, sfi_191, sfi_192, sgh0_227, \
                         sgh0_228, sgh1_227, sgh1_228, sgi_301, sgi_303, \
                         sgi_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_3 * pc_z[k] * sgi_301[k];

        t_390[k] = f_0 * sfi_191[k]
                   + f_4 * sgh0_227[k]
                   - f_5 * sgh1_227[k]
                   + f_3 * pc_y[k] * sgi_303[k];

        t_391[k] = f_0 * sfi_192[k]
                   + f_6 * sgh0_228[k]
                   - f_7 * sgh1_228[k]
                   + f_3 * pc_y[k] * sgi_304[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pc_y, pc_z, sfi_193, sfi_194, sfi_195, \
                         sgh0_229, sgh0_230, sgh1_229, sgh1_230, sgi_305, sgi_306, \
                         sgi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_0 * sfi_193[k]
                   + f_8 * sgh0_229[k]
                   - f_9 * sgh1_229[k]
                   + f_3 * pc_y[k] * sgi_305[k];

        t_393[k] = f_0 * sfi_194[k]
                   + f_10 * sgh0_230[k]
                   - f_11 * sgh1_230[k]
                   + f_3 * pc_y[k] * sgi_306[k];

        t_394[k] = f_0 * sfi_195[k]
                   + f_3 * pc_y[k] * sgi_307[k];

        t_395[k] = f_1 * sgh0_230[k]
                   - f_2 * sgh1_230[k]
                   + f_3 * pc_z[k] * sgi_307[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pb_z, pc_y, pc_z, sfk0_216, sfk0_219, \
                         sfi_168, sfi_196, sfk1_216, sfk1_219, \
                         sgi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = pb_z[k] * sfk0_216[k]
                   - f_12 * pc_z[k] * sfk1_216[k];

        t_397[k] = f_15 * sfi_196[k]
                   + f_3 * pc_y[k] * sgi_308[k];

        t_398[k] = f_13 * sfi_168[k]
                   + f_3 * pc_z[k] * sgi_308[k];

        t_399[k] = pb_z[k] * sfk0_219[k]
                   - f_12 * pc_z[k] * sfk1_219[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pb_z, pc_x, pc_y, pc_z, sfk0_222, sfi_198, \
                         sfk1_222, sgh0_236, sgh1_236, sgi_310, \
                         sgi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_15 * sfi_198[k]
                   + f_3 * pc_y[k] * sgi_310[k];

        t_401[k] = f_4 * sgh0_236[k]
                   - f_5 * sgh1_236[k]
                   + f_3 * pc_x[k] * sgi_313[k];

        t_402[k] = pb_z[k] * sfk0_222[k]
                   - f_12 * pc_z[k] * sfk1_222[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, pc_x, pc_y, pc_z, sfi_171, sfi_201, sgh0_240, \
                         sgh1_240, sgi_311, sgi_313, sgi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_13 * sfi_171[k]
                   + f_3 * pc_z[k] * sgi_311[k];

        t_404[k] = f_15 * sfi_201[k]
                   + f_3 * pc_y[k] * sgi_313[k];

        t_405[k] = f_6 * sgh0_240[k]
                   - f_7 * sgh1_240[k]
                   + f_3 * pc_x[k] * sgi_317[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, pb_z, pc_x, pc_z, sfk0_226, sfi_174, sfk1_226, \
                         sgh0_243, sgh1_243, sgi_314, sgi_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pb_z[k] * sfk0_226[k]
                   - f_12 * pc_z[k] * sfk1_226[k];

        t_407[k] = f_13 * sfi_174[k]
                   + f_3 * pc_z[k] * sgi_314[k];

        t_408[k] = f_8 * sgh0_243[k]
                   - f_9 * sgh1_243[k]
                   + f_3 * pc_x[k] * sgi_320[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, pb_z, pc_x, pc_y, pc_z, sfk0_231, sfi_205, \
                         sfk1_231, sgh0_245, sgh1_245, sgi_317, \
                         sgi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_15 * sfi_205[k]
                   + f_3 * pc_y[k] * sgi_317[k];

        t_410[k] = f_8 * sgh0_245[k]
                   - f_9 * sgh1_245[k]
                   + f_3 * pc_x[k] * sgi_322[k];

        t_411[k] = pb_z[k] * sfk0_231[k]
                   - f_12 * pc_z[k] * sfk1_231[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, pc_x, pc_z, sfi_178, sgh0_248, sgh0_249, \
                         sgh1_248, sgh1_249, sgi_318, sgi_325, \
                         sgi_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_13 * sfi_178[k]
                   + f_3 * pc_z[k] * sgi_318[k];

        t_413[k] = f_10 * sgh0_248[k]
                   - f_11 * sgh1_248[k]
                   + f_3 * pc_x[k] * sgi_325[k];

        t_414[k] = f_10 * sgh0_249[k]
                   - f_11 * sgh1_249[k]
                   + f_3 * pc_x[k] * sgi_326[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, pc_x, pc_y, sfi_210, sgh0_251, \
                         sgh1_251, sgi_322, sgi_328, sgi_329, sgi_330, \
                         sgi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_15 * sfi_210[k]
                   + f_3 * pc_y[k] * sgi_322[k];

        t_416[k] = f_10 * sgh0_251[k]
                   - f_11 * sgh1_251[k]
                   + f_3 * pc_x[k] * sgi_328[k];

        t_417[k] = f_3 * pc_x[k] * sgi_329[k];

        t_418[k] = f_3 * pc_x[k] * sgi_330[k];

        t_419[k] = f_3 * pc_x[k] * sgi_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pb_z, pc_x, pc_z, sfk0_244, \
                         sfk1_244, sgi_332, sgi_333, sgi_334, sgi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_3 * pc_x[k] * sgi_332[k];

        t_421[k] = f_3 * pc_x[k] * sgi_333[k];

        t_422[k] = f_3 * pc_x[k] * sgi_334[k];

        t_423[k] = f_3 * pc_x[k] * sgi_335[k];

        t_424[k] = pb_z[k] * sfk0_244[k]
                   - f_12 * pc_z[k] * sfk1_244[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, pb_z, pc_z, sfk0_246, sfk0_247, sfi_189, \
                         sfi_190, sfi_191, sfk1_246, sfk1_247, \
                         sgi_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_13 * sfi_189[k]
                   + f_3 * pc_z[k] * sgi_329[k];

        t_426[k] = pb_z[k] * sfk0_246[k]
                   + f_14 * sfi_190[k]
                   - f_12 * pc_z[k] * sfk1_246[k];

        t_427[k] = pb_z[k] * sfk0_247[k]
                   + f_15 * sfi_191[k]
                   - f_12 * pc_z[k] * sfk1_247[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, pb_z, pc_y, pc_z, sfk0_248, sfk0_249, sfi_192, \
                         sfi_193, sfi_223, sfk1_248, sfk1_249, \
                         sgi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = pb_z[k] * sfk0_248[k]
                   + f_0 * sfi_192[k]
                   - f_12 * pc_z[k] * sfk1_248[k];

        t_429[k] = pb_z[k] * sfk0_249[k]
                   + f_16 * sfi_193[k]
                   - f_12 * pc_z[k] * sfk1_249[k];

        t_430[k] = f_15 * sfi_223[k]
                   + f_3 * pc_y[k] * sgi_335[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, pc_y, pc_z, sfi_195, sfi_196, \
                         sfi_224, sgh0_251, sgh0_252, sgh1_251, sgh1_252, sgi_335, \
                         sgi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_13 * sfi_195[k]
                   + f_1 * sgh0_251[k]
                   - f_2 * sgh1_251[k]
                   + f_3 * pc_z[k] * sgi_335[k];

        t_432[k] = f_1 * sgh0_252[k]
                   - f_2 * sgh1_252[k]
                   + f_3 * pc_x[k] * sgi_336[k];

        t_433[k] = f_14 * sfi_224[k]
                   + f_3 * pc_y[k] * sgi_336[k];

        t_434[k] = f_14 * sfi_196[k]
                   + f_3 * pc_z[k] * sgi_336[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, pc_x, pc_y, sfi_226, sgh0_255, sgh0_257, \
                         sgh1_255, sgh1_257, sgi_338, sgi_339, \
                         sgi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_4 * sgh0_255[k]
                   - f_5 * sgh1_255[k]
                   + f_3 * pc_x[k] * sgi_339[k];

        t_436[k] = f_14 * sfi_226[k]
                   + f_3 * pc_y[k] * sgi_338[k];

        t_437[k] = f_4 * sgh0_257[k]
                   - f_5 * sgh1_257[k]
                   + f_3 * pc_x[k] * sgi_341[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pc_x, pc_y, pc_z, sfi_199, sfi_229, sgh0_258, \
                         sgh1_258, sgi_339, sgi_341, sgi_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_6 * sgh0_258[k]
                   - f_7 * sgh1_258[k]
                   + f_3 * pc_x[k] * sgi_342[k];

        t_439[k] = f_14 * sfi_199[k]
                   + f_3 * pc_z[k] * sgi_339[k];

        t_440[k] = f_14 * sfi_229[k]
                   + f_3 * pc_y[k] * sgi_341[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_x, pc_z, sfi_202, sgh0_261, sgh0_262, \
                         sgh1_261, sgh1_262, sgi_342, sgi_345, \
                         sgi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_6 * sgh0_261[k]
                   - f_7 * sgh1_261[k]
                   + f_3 * pc_x[k] * sgi_345[k];

        t_442[k] = f_8 * sgh0_262[k]
                   - f_9 * sgh1_262[k]
                   + f_3 * pc_x[k] * sgi_346[k];

        t_443[k] = f_14 * sfi_202[k]
                   + f_3 * pc_z[k] * sgi_342[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_x, pc_y, sfi_233, sgh0_264, sgh0_266, \
                         sgh1_264, sgh1_266, sgi_345, sgi_348, \
                         sgi_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_8 * sgh0_264[k]
                   - f_9 * sgh1_264[k]
                   + f_3 * pc_x[k] * sgi_348[k];

        t_445[k] = f_14 * sfi_233[k]
                   + f_3 * pc_y[k] * sgi_345[k];

        t_446[k] = f_8 * sgh0_266[k]
                   - f_9 * sgh1_266[k]
                   + f_3 * pc_x[k] * sgi_350[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pc_x, pc_z, sfi_206, sgh0_267, sgh0_269, \
                         sgh1_267, sgh1_269, sgi_346, sgi_351, \
                         sgi_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_10 * sgh0_267[k]
                   - f_11 * sgh1_267[k]
                   + f_3 * pc_x[k] * sgi_351[k];

        t_448[k] = f_14 * sfi_206[k]
                   + f_3 * pc_z[k] * sgi_346[k];

        t_449[k] = f_10 * sgh0_269[k]
                   - f_11 * sgh1_269[k]
                   + f_3 * pc_x[k] * sgi_353[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, pc_x, pc_y, sfi_238, sgh0_270, sgh0_272, \
                         sgh1_270, sgh1_272, sgi_350, sgi_354, sgi_356, \
                         sgi_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_10 * sgh0_270[k]
                   - f_11 * sgh1_270[k]
                   + f_3 * pc_x[k] * sgi_354[k];

        t_451[k] = f_14 * sfi_238[k]
                   + f_3 * pc_y[k] * sgi_350[k];

        t_452[k] = f_10 * sgh0_272[k]
                   - f_11 * sgh1_272[k]
                   + f_3 * pc_x[k] * sgi_356[k];

        t_453[k] = f_3 * pc_x[k] * sgi_357[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, t_459, pc_x, sgi_358, sgi_359, \
                         sgi_360, sgi_361, sgi_362, sgi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_3 * pc_x[k] * sgi_358[k];

        t_455[k] = f_3 * pc_x[k] * sgi_359[k];

        t_456[k] = f_3 * pc_x[k] * sgi_360[k];

        t_457[k] = f_3 * pc_x[k] * sgi_361[k];

        t_458[k] = f_3 * pc_x[k] * sgi_362[k];

        t_459[k] = f_3 * pc_x[k] * sgi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, sfi_217, sfi_245, sfi_247, sgh0_267, \
                         sgh0_269, sgh1_267, sgh1_269, sgi_357, \
                         sgi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * sfi_245[k]
                   + f_1 * sgh0_267[k]
                   - f_2 * sgh1_267[k]
                   + f_3 * pc_y[k] * sgi_357[k];

        t_461[k] = f_14 * sfi_217[k]
                   + f_3 * pc_z[k] * sgi_357[k];

        t_462[k] = f_14 * sfi_247[k]
                   + f_4 * sgh0_269[k]
                   - f_5 * sgh1_269[k]
                   + f_3 * pc_y[k] * sgi_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, sfi_248, sfi_249, sfi_250, sgh0_270, \
                         sgh0_271, sgh0_272, sgh1_270, sgh1_271, sgh1_272, sgi_360, sgi_361, \
                         sgi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * sfi_248[k]
                   + f_6 * sgh0_270[k]
                   - f_7 * sgh1_270[k]
                   + f_3 * pc_y[k] * sgi_360[k];

        t_464[k] = f_14 * sfi_249[k]
                   + f_8 * sgh0_271[k]
                   - f_9 * sgh1_271[k]
                   + f_3 * pc_y[k] * sgi_361[k];

        t_465[k] = f_14 * sfi_250[k]
                   + f_10 * sgh0_272[k]
                   - f_11 * sgh1_272[k]
                   + f_3 * pc_y[k] * sgi_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_y, pc_y, pc_z, sfk0_324, sfi_223, \
                         sfi_251, sfi_252, sfk1_324, sgh0_272, sgh1_272, sgi_363, \
                         sgi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * sfi_251[k]
                   + f_3 * pc_y[k] * sgi_363[k];

        t_467[k] = f_14 * sfi_223[k]
                   + f_1 * sgh0_272[k]
                   - f_2 * sgh1_272[k]
                   + f_3 * pc_z[k] * sgi_363[k];

        t_468[k] = pb_y[k] * sfk0_324[k]
                   - f_12 * pc_y[k] * sfk1_324[k];

        t_469[k] = f_13 * sfi_252[k]
                   + f_3 * pc_y[k] * sgi_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_x, pc_y, pc_z, sfi_224, sfi_254, sgh0_276, \
                         sgh1_276, sgi_364, sgi_366, sgi_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * sfi_224[k]
                   + f_3 * pc_z[k] * sgi_364[k];

        t_471[k] = f_4 * sgh0_276[k]
                   - f_5 * sgh1_276[k]
                   + f_3 * pc_x[k] * sgi_367[k];

        t_472[k] = f_13 * sfi_254[k]
                   + f_3 * pc_y[k] * sgi_366[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pb_y, pc_x, pc_y, pc_z, sfk0_329, sfi_227, \
                         sfk1_329, sgh0_279, sgh1_279, sgi_367, \
                         sgi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = pb_y[k] * sfk0_329[k]
                   - f_12 * pc_y[k] * sfk1_329[k];

        t_474[k] = f_6 * sgh0_279[k]
                   - f_7 * sgh1_279[k]
                   + f_3 * pc_x[k] * sgi_370[k];

        t_475[k] = f_15 * sfi_227[k]
                   + f_3 * pc_z[k] * sgi_367[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pb_y, pc_x, pc_y, sfk0_333, sfi_257, sfk1_333, \
                         sgh0_283, sgh1_283, sgi_369, sgi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_13 * sfi_257[k]
                   + f_3 * pc_y[k] * sgi_369[k];

        t_477[k] = pb_y[k] * sfk0_333[k]
                   - f_12 * pc_y[k] * sfk1_333[k];

        t_478[k] = f_8 * sgh0_283[k]
                   - f_9 * sgh1_283[k]
                   + f_3 * pc_x[k] * sgi_374[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, pc_z, sfi_230, sfi_261, sgh0_285, \
                         sgh1_285, sgi_370, sgi_373, sgi_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_15 * sfi_230[k]
                   + f_3 * pc_z[k] * sgi_370[k];

        t_480[k] = f_8 * sgh0_285[k]
                   - f_9 * sgh1_285[k]
                   + f_3 * pc_x[k] * sgi_376[k];

        t_481[k] = f_13 * sfi_261[k]
                   + f_3 * pc_y[k] * sgi_373[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pb_y, pc_x, pc_y, pc_z, sfk0_338, sfi_234, \
                         sfk1_338, sgh0_288, sgh1_288, sgi_374, \
                         sgi_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = pb_y[k] * sfk0_338[k]
                   - f_12 * pc_y[k] * sfk1_338[k];

        t_483[k] = f_10 * sgh0_288[k]
                   - f_11 * sgh1_288[k]
                   + f_3 * pc_x[k] * sgi_379[k];

        t_484[k] = f_15 * sfi_234[k]
                   + f_3 * pc_z[k] * sgi_374[k];
    }
}

static auto
compute_prim_sgk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfk0,
                                                          const size_t sfi, const size_t sfk1,
                                                          const size_t sgh0, const size_t sgh1,
                                                          const size_t sgi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 3.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfk0_344 = buffer.data(sfk0 + 344);
    const auto *sfk0_352 = buffer.data(sfk0 + 352);
    const auto *sfk0_354 = buffer.data(sfk0 + 354);
    const auto *sfk0_355 = buffer.data(sfk0 + 355);
    const auto *sfk0_356 = buffer.data(sfk0 + 356);
    const auto *sfk0_357 = buffer.data(sfk0 + 357);
    const auto *sfk0_359 = buffer.data(sfk0 + 359);

    const auto *sfi_245 = buffer.data(sfi + 245);
    const auto *sfi_252 = buffer.data(sfi + 252);
    const auto *sfi_255 = buffer.data(sfi + 255);
    const auto *sfi_258 = buffer.data(sfi + 258);
    const auto *sfi_262 = buffer.data(sfi + 262);
    const auto *sfi_266 = buffer.data(sfi + 266);
    const auto *sfi_273 = buffer.data(sfi + 273);
    const auto *sfi_275 = buffer.data(sfi + 275);
    const auto *sfi_276 = buffer.data(sfi + 276);
    const auto *sfi_277 = buffer.data(sfi + 277);
    const auto *sfi_278 = buffer.data(sfi + 278);
    const auto *sfi_279 = buffer.data(sfi + 279);

    const auto *sfk1_344 = buffer.data(sfk1 + 344);
    const auto *sfk1_352 = buffer.data(sfk1 + 352);
    const auto *sfk1_354 = buffer.data(sfk1 + 354);
    const auto *sfk1_355 = buffer.data(sfk1 + 355);
    const auto *sfk1_356 = buffer.data(sfk1 + 356);
    const auto *sfk1_357 = buffer.data(sfk1 + 357);
    const auto *sfk1_359 = buffer.data(sfk1 + 359);

    const auto *sgh0_290 = buffer.data(sgh0 + 290);
    const auto *sgh0_291 = buffer.data(sgh0 + 291);
    const auto *sgh0_294 = buffer.data(sgh0 + 294);
    const auto *sgh0_297 = buffer.data(sgh0 + 297);
    const auto *sgh0_299 = buffer.data(sgh0 + 299);
    const auto *sgh0_300 = buffer.data(sgh0 + 300);
    const auto *sgh0_303 = buffer.data(sgh0 + 303);
    const auto *sgh0_304 = buffer.data(sgh0 + 304);
    const auto *sgh0_306 = buffer.data(sgh0 + 306);
    const auto *sgh0_308 = buffer.data(sgh0 + 308);
    const auto *sgh0_309 = buffer.data(sgh0 + 309);
    const auto *sgh0_311 = buffer.data(sgh0 + 311);
    const auto *sgh0_312 = buffer.data(sgh0 + 312);
    const auto *sgh0_313 = buffer.data(sgh0 + 313);
    const auto *sgh0_314 = buffer.data(sgh0 + 314);

    const auto *sgh1_290 = buffer.data(sgh1 + 290);
    const auto *sgh1_291 = buffer.data(sgh1 + 291);
    const auto *sgh1_294 = buffer.data(sgh1 + 294);
    const auto *sgh1_297 = buffer.data(sgh1 + 297);
    const auto *sgh1_299 = buffer.data(sgh1 + 299);
    const auto *sgh1_300 = buffer.data(sgh1 + 300);
    const auto *sgh1_303 = buffer.data(sgh1 + 303);
    const auto *sgh1_304 = buffer.data(sgh1 + 304);
    const auto *sgh1_306 = buffer.data(sgh1 + 306);
    const auto *sgh1_308 = buffer.data(sgh1 + 308);
    const auto *sgh1_309 = buffer.data(sgh1 + 309);
    const auto *sgh1_311 = buffer.data(sgh1 + 311);
    const auto *sgh1_312 = buffer.data(sgh1 + 312);
    const auto *sgh1_313 = buffer.data(sgh1 + 313);
    const auto *sgh1_314 = buffer.data(sgh1 + 314);

    const auto *sgi_378 = buffer.data(sgi + 378);
    const auto *sgi_381 = buffer.data(sgi + 381);
    const auto *sgi_382 = buffer.data(sgi + 382);
    const auto *sgi_385 = buffer.data(sgi + 385);
    const auto *sgi_386 = buffer.data(sgi + 386);
    const auto *sgi_387 = buffer.data(sgi + 387);
    const auto *sgi_388 = buffer.data(sgi + 388);
    const auto *sgi_389 = buffer.data(sgi + 389);
    const auto *sgi_390 = buffer.data(sgi + 390);
    const auto *sgi_391 = buffer.data(sgi + 391);
    const auto *sgi_392 = buffer.data(sgi + 392);
    const auto *sgi_394 = buffer.data(sgi + 394);
    const auto *sgi_395 = buffer.data(sgi + 395);
    const auto *sgi_397 = buffer.data(sgi + 397);
    const auto *sgi_398 = buffer.data(sgi + 398);
    const auto *sgi_401 = buffer.data(sgi + 401);
    const auto *sgi_402 = buffer.data(sgi + 402);
    const auto *sgi_404 = buffer.data(sgi + 404);
    const auto *sgi_406 = buffer.data(sgi + 406);
    const auto *sgi_407 = buffer.data(sgi + 407);
    const auto *sgi_409 = buffer.data(sgi + 409);
    const auto *sgi_410 = buffer.data(sgi + 410);
    const auto *sgi_412 = buffer.data(sgi + 412);
    const auto *sgi_413 = buffer.data(sgi + 413);
    const auto *sgi_414 = buffer.data(sgi + 414);
    const auto *sgi_415 = buffer.data(sgi + 415);
    const auto *sgi_416 = buffer.data(sgi + 416);
    const auto *sgi_417 = buffer.data(sgi + 417);
    const auto *sgi_418 = buffer.data(sgi + 418);
    const auto *sgi_419 = buffer.data(sgi + 419);

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_y, sfi_266, sgh0_290, sgh0_291, \
                         sgh1_290, sgh1_291, sgi_378, sgi_381, \
                         sgi_382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_10 * sgh0_290[k]
                   - f_11 * sgh1_290[k]
                   + f_3 * pc_x[k] * sgi_381[k];

        t_486[k] = f_10 * sgh0_291[k]
                   - f_11 * sgh1_291[k]
                   + f_3 * pc_x[k] * sgi_382[k];

        t_487[k] = f_13 * sfi_266[k]
                   + f_3 * pc_y[k] * sgi_378[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, t_492, t_493, pb_y, pc_x, pc_y, sfk0_344, \
                         sfk1_344, sgi_385, sgi_386, sgi_387, sgi_388, \
                         sgi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = pb_y[k] * sfk0_344[k]
                   - f_12 * pc_y[k] * sfk1_344[k];

        t_489[k] = f_3 * pc_x[k] * sgi_385[k];

        t_490[k] = f_3 * pc_x[k] * sgi_386[k];

        t_491[k] = f_3 * pc_x[k] * sgi_387[k];

        t_492[k] = f_3 * pc_x[k] * sgi_388[k];

        t_493[k] = f_3 * pc_x[k] * sgi_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pb_y, pc_x, pc_y, pc_z, sfk0_352, \
                         sfi_245, sfi_273, sfk1_352, sgi_385, sgi_390, \
                         sgi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_3 * pc_x[k] * sgi_390[k];

        t_495[k] = f_3 * pc_x[k] * sgi_391[k];

        t_496[k] = pb_y[k] * sfk0_352[k]
                   + f_17 * sfi_273[k]
                   - f_12 * pc_y[k] * sfk1_352[k];

        t_497[k] = f_15 * sfi_245[k]
                   + f_3 * pc_z[k] * sgi_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pb_y, pc_y, sfk0_354, sfk0_355, sfk0_356, \
                         sfi_275, sfi_276, sfi_277, sfk1_354, sfk1_355, \
                         sfk1_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = pb_y[k] * sfk0_354[k]
                   + f_16 * sfi_275[k]
                   - f_12 * pc_y[k] * sfk1_354[k];

        t_499[k] = pb_y[k] * sfk0_355[k]
                   + f_0 * sfi_276[k]
                   - f_12 * pc_y[k] * sfk1_355[k];

        t_500[k] = pb_y[k] * sfk0_356[k]
                   + f_15 * sfi_277[k]
                   - f_12 * pc_y[k] * sfk1_356[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pb_y, pc_y, sfk0_357, sfk0_359, sfi_278, \
                         sfi_279, sfk1_357, sfk1_359, sgi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = pb_y[k] * sfk0_357[k]
                   + f_14 * sfi_278[k]
                   - f_12 * pc_y[k] * sfk1_357[k];

        t_502[k] = f_13 * sfi_279[k]
                   + f_3 * pc_y[k] * sgi_391[k];

        t_503[k] = pb_y[k] * sfk0_359[k]
                   - f_12 * pc_y[k] * sfk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, sfi_252, \
                         sgh0_294, sgh0_297, sgh1_294, sgh1_297, sgi_392, sgi_394, \
                         sgi_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_1 * sgh0_294[k]
                   - f_2 * sgh1_294[k]
                   + f_3 * pc_x[k] * sgi_392[k];

        t_505[k] = f_3 * pc_y[k] * sgi_392[k];

        t_506[k] = f_0 * sfi_252[k]
                   + f_3 * pc_z[k] * sgi_392[k];

        t_507[k] = f_4 * sgh0_297[k]
                   - f_5 * sgh1_297[k]
                   + f_3 * pc_x[k] * sgi_395[k];

        t_508[k] = f_3 * pc_y[k] * sgi_394[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, t_512, pc_x, pc_y, pc_z, sfi_255, sgh0_299, \
                         sgh0_300, sgh1_299, sgh1_300, sgi_395, sgi_397, \
                         sgi_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_4 * sgh0_299[k]
                   - f_5 * sgh1_299[k]
                   + f_3 * pc_x[k] * sgi_397[k];

        t_510[k] = f_6 * sgh0_300[k]
                   - f_7 * sgh1_300[k]
                   + f_3 * pc_x[k] * sgi_398[k];

        t_511[k] = f_0 * sfi_255[k]
                   + f_3 * pc_z[k] * sgi_395[k];

        t_512[k] = f_3 * pc_y[k] * sgi_397[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_z, sfi_258, sgh0_303, sgh0_304, \
                         sgh1_303, sgh1_304, sgi_398, sgi_401, \
                         sgi_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_6 * sgh0_303[k]
                   - f_7 * sgh1_303[k]
                   + f_3 * pc_x[k] * sgi_401[k];

        t_514[k] = f_8 * sgh0_304[k]
                   - f_9 * sgh1_304[k]
                   + f_3 * pc_x[k] * sgi_402[k];

        t_515[k] = f_0 * sfi_258[k]
                   + f_3 * pc_z[k] * sgi_398[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, sgh0_306, sgh0_308, sgh0_309, \
                         sgh1_306, sgh1_308, sgh1_309, sgi_401, sgi_404, sgi_406, \
                         sgi_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_8 * sgh0_306[k]
                   - f_9 * sgh1_306[k]
                   + f_3 * pc_x[k] * sgi_404[k];

        t_517[k] = f_3 * pc_y[k] * sgi_401[k];

        t_518[k] = f_8 * sgh0_308[k]
                   - f_9 * sgh1_308[k]
                   + f_3 * pc_x[k] * sgi_406[k];

        t_519[k] = f_10 * sgh0_309[k]
                   - f_11 * sgh1_309[k]
                   + f_3 * pc_x[k] * sgi_407[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, pc_x, pc_y, pc_z, sfi_262, sgh0_311, \
                         sgh0_312, sgh1_311, sgh1_312, sgi_402, sgi_406, sgi_409, \
                         sgi_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_0 * sfi_262[k]
                   + f_3 * pc_z[k] * sgi_402[k];

        t_521[k] = f_10 * sgh0_311[k]
                   - f_11 * sgh1_311[k]
                   + f_3 * pc_x[k] * sgi_409[k];

        t_522[k] = f_10 * sgh0_312[k]
                   - f_11 * sgh1_312[k]
                   + f_3 * pc_x[k] * sgi_410[k];

        t_523[k] = f_3 * pc_y[k] * sgi_406[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, t_528, t_529, pc_x, sgh0_314, sgh1_314, \
                         sgi_412, sgi_413, sgi_414, sgi_415, sgi_416, \
                         sgi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_10 * sgh0_314[k]
                   - f_11 * sgh1_314[k]
                   + f_3 * pc_x[k] * sgi_412[k];

        t_525[k] = f_3 * pc_x[k] * sgi_413[k];

        t_526[k] = f_3 * pc_x[k] * sgi_414[k];

        t_527[k] = f_3 * pc_x[k] * sgi_415[k];

        t_528[k] = f_3 * pc_x[k] * sgi_416[k];

        t_529[k] = f_3 * pc_x[k] * sgi_417[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pc_x, pc_y, pc_z, sfi_273, sgh0_309, \
                         sgh1_309, sgi_413, sgi_418, sgi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_3 * pc_x[k] * sgi_418[k];

        t_531[k] = f_3 * pc_x[k] * sgi_419[k];

        t_532[k] = f_1 * sgh0_309[k]
                   - f_2 * sgh1_309[k]
                   + f_3 * pc_y[k] * sgi_413[k];

        t_533[k] = f_0 * sfi_273[k]
                   + f_3 * pc_z[k] * sgi_413[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_y, sgh0_311, sgh0_312, sgh0_313, sgh1_311, \
                         sgh1_312, sgh1_313, sgi_415, sgi_416, \
                         sgi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_4 * sgh0_311[k]
                   - f_5 * sgh1_311[k]
                   + f_3 * pc_y[k] * sgi_415[k];

        t_535[k] = f_6 * sgh0_312[k]
                   - f_7 * sgh1_312[k]
                   + f_3 * pc_y[k] * sgi_416[k];

        t_536[k] = f_8 * sgh0_313[k]
                   - f_9 * sgh1_313[k]
                   + f_3 * pc_y[k] * sgi_417[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pc_y, pc_z, sfi_279, sgh0_314, sgh1_314, \
                         sgi_418, sgi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_10 * sgh0_314[k]
                   - f_11 * sgh1_314[k]
                   + f_3 * pc_y[k] * sgi_418[k];

        t_538[k] = f_3 * pc_y[k] * sgi_419[k];

        t_539[k] = f_0 * sfi_279[k]
                   + f_1 * sgh0_314[k]
                   - f_2 * sgh1_314[k]
                   + f_3 * pc_z[k] * sgi_419[k];
    }
}

auto
compute_prim_sgk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sfk0, const size_t sfi,
                                                   const size_t sfk1, const size_t sgh0,
                                                   const size_t sgh1, const size_t sgi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sgk_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sfk0, sfi,
                                                              sfk1, sgh0, sgh1, sgi, ncols,
                                                              gamma, p, q);

    compute_prim_sgk_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sfk0, sfi,
                                                              sfk1, sgh0, sgh1, sgi, ncols,
                                                              gamma, p, q);

    compute_prim_sgk_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sfk0, sfi,
                                                              sfk1, sgh0, sgh1, sgi, ncols,
                                                              gamma, p, q);

    compute_prim_sgk_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sfk0, sfi,
                                                              sfk1, sgh0, sgh1, sgi, ncols,
                                                              gamma, p, q);

    compute_prim_sgk_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sfk0, sfi,
                                                              sfk1, sgh0, sgh1, sgi, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
