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


#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sgf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sff0,
                                                          const size_t sfd, const size_t sff1,
                                                          const size_t sgp0, const size_t sgp1,
                                                          const size_t sgd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 1.0 / q;

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
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);
    auto *t_135 = buffer.data(target + 135);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sff0_0 = buffer.data(sff0 + 0);
    const auto *sff0_6 = buffer.data(sff0 + 6);
    const auto *sff0_9 = buffer.data(sff0 + 9);
    const auto *sff0_16 = buffer.data(sff0 + 16);
    const auto *sff0_20 = buffer.data(sff0 + 20);
    const auto *sff0_29 = buffer.data(sff0 + 29);
    const auto *sff0_30 = buffer.data(sff0 + 30);
    const auto *sff0_50 = buffer.data(sff0 + 50);
    const auto *sff0_60 = buffer.data(sff0 + 60);
    const auto *sff0_66 = buffer.data(sff0 + 66);
    const auto *sff0_69 = buffer.data(sff0 + 69);
    const auto *sff0_76 = buffer.data(sff0 + 76);
    const auto *sff0_79 = buffer.data(sff0 + 79);
    const auto *sff0_86 = buffer.data(sff0 + 86);
    const auto *sff0_89 = buffer.data(sff0 + 89);
    const auto *sff0_90 = buffer.data(sff0 + 90);
    const auto *sff0_96 = buffer.data(sff0 + 96);
    const auto *sff0_99 = buffer.data(sff0 + 99);

    const auto *sfd_0 = buffer.data(sfd + 0);
    const auto *sfd_3 = buffer.data(sfd + 3);
    const auto *sfd_4 = buffer.data(sfd + 4);
    const auto *sfd_5 = buffer.data(sfd + 5);
    const auto *sfd_6 = buffer.data(sfd + 6);
    const auto *sfd_9 = buffer.data(sfd + 9);
    const auto *sfd_10 = buffer.data(sfd + 10);
    const auto *sfd_11 = buffer.data(sfd + 11);
    const auto *sfd_12 = buffer.data(sfd + 12);
    const auto *sfd_15 = buffer.data(sfd + 15);
    const auto *sfd_16 = buffer.data(sfd + 16);
    const auto *sfd_17 = buffer.data(sfd + 17);
    const auto *sfd_18 = buffer.data(sfd + 18);
    const auto *sfd_21 = buffer.data(sfd + 21);
    const auto *sfd_22 = buffer.data(sfd + 22);
    const auto *sfd_23 = buffer.data(sfd + 23);
    const auto *sfd_24 = buffer.data(sfd + 24);
    const auto *sfd_27 = buffer.data(sfd + 27);
    const auto *sfd_28 = buffer.data(sfd + 28);
    const auto *sfd_29 = buffer.data(sfd + 29);
    const auto *sfd_30 = buffer.data(sfd + 30);
    const auto *sfd_33 = buffer.data(sfd + 33);
    const auto *sfd_34 = buffer.data(sfd + 34);
    const auto *sfd_35 = buffer.data(sfd + 35);
    const auto *sfd_36 = buffer.data(sfd + 36);
    const auto *sfd_39 = buffer.data(sfd + 39);
    const auto *sfd_40 = buffer.data(sfd + 40);
    const auto *sfd_41 = buffer.data(sfd + 41);
    const auto *sfd_42 = buffer.data(sfd + 42);
    const auto *sfd_45 = buffer.data(sfd + 45);
    const auto *sfd_46 = buffer.data(sfd + 46);
    const auto *sfd_47 = buffer.data(sfd + 47);
    const auto *sfd_48 = buffer.data(sfd + 48);
    const auto *sfd_51 = buffer.data(sfd + 51);
    const auto *sfd_52 = buffer.data(sfd + 52);
    const auto *sfd_53 = buffer.data(sfd + 53);
    const auto *sfd_54 = buffer.data(sfd + 54);
    const auto *sfd_57 = buffer.data(sfd + 57);
    const auto *sfd_58 = buffer.data(sfd + 58);
    const auto *sfd_59 = buffer.data(sfd + 59);

    const auto *sff1_0 = buffer.data(sff1 + 0);
    const auto *sff1_6 = buffer.data(sff1 + 6);
    const auto *sff1_9 = buffer.data(sff1 + 9);
    const auto *sff1_16 = buffer.data(sff1 + 16);
    const auto *sff1_20 = buffer.data(sff1 + 20);
    const auto *sff1_29 = buffer.data(sff1 + 29);
    const auto *sff1_30 = buffer.data(sff1 + 30);
    const auto *sff1_50 = buffer.data(sff1 + 50);
    const auto *sff1_60 = buffer.data(sff1 + 60);
    const auto *sff1_66 = buffer.data(sff1 + 66);
    const auto *sff1_69 = buffer.data(sff1 + 69);
    const auto *sff1_76 = buffer.data(sff1 + 76);
    const auto *sff1_79 = buffer.data(sff1 + 79);
    const auto *sff1_86 = buffer.data(sff1 + 86);
    const auto *sff1_89 = buffer.data(sff1 + 89);
    const auto *sff1_90 = buffer.data(sff1 + 90);
    const auto *sff1_96 = buffer.data(sff1 + 96);
    const auto *sff1_99 = buffer.data(sff1 + 99);

    const auto *sgp0_0 = buffer.data(sgp0 + 0);
    const auto *sgp0_1 = buffer.data(sgp0 + 1);
    const auto *sgp0_2 = buffer.data(sgp0 + 2);
    const auto *sgp0_4 = buffer.data(sgp0 + 4);
    const auto *sgp0_8 = buffer.data(sgp0 + 8);
    const auto *sgp0_9 = buffer.data(sgp0 + 9);
    const auto *sgp0_10 = buffer.data(sgp0 + 10);
    const auto *sgp0_11 = buffer.data(sgp0 + 11);
    const auto *sgp0_15 = buffer.data(sgp0 + 15);
    const auto *sgp0_16 = buffer.data(sgp0 + 16);
    const auto *sgp0_17 = buffer.data(sgp0 + 17);
    const auto *sgp0_30 = buffer.data(sgp0 + 30);
    const auto *sgp0_31 = buffer.data(sgp0 + 31);
    const auto *sgp0_32 = buffer.data(sgp0 + 32);
    const auto *sgp0_35 = buffer.data(sgp0 + 35);
    const auto *sgp0_36 = buffer.data(sgp0 + 36);
    const auto *sgp0_37 = buffer.data(sgp0 + 37);
    const auto *sgp0_38 = buffer.data(sgp0 + 38);

    const auto *sgp1_0 = buffer.data(sgp1 + 0);
    const auto *sgp1_1 = buffer.data(sgp1 + 1);
    const auto *sgp1_2 = buffer.data(sgp1 + 2);
    const auto *sgp1_4 = buffer.data(sgp1 + 4);
    const auto *sgp1_8 = buffer.data(sgp1 + 8);
    const auto *sgp1_9 = buffer.data(sgp1 + 9);
    const auto *sgp1_10 = buffer.data(sgp1 + 10);
    const auto *sgp1_11 = buffer.data(sgp1 + 11);
    const auto *sgp1_15 = buffer.data(sgp1 + 15);
    const auto *sgp1_16 = buffer.data(sgp1 + 16);
    const auto *sgp1_17 = buffer.data(sgp1 + 17);
    const auto *sgp1_30 = buffer.data(sgp1 + 30);
    const auto *sgp1_31 = buffer.data(sgp1 + 31);
    const auto *sgp1_32 = buffer.data(sgp1 + 32);
    const auto *sgp1_35 = buffer.data(sgp1 + 35);
    const auto *sgp1_36 = buffer.data(sgp1 + 36);
    const auto *sgp1_37 = buffer.data(sgp1 + 37);
    const auto *sgp1_38 = buffer.data(sgp1 + 38);

    const auto *sgd_0 = buffer.data(sgd + 0);
    const auto *sgd_3 = buffer.data(sgd + 3);
    const auto *sgd_4 = buffer.data(sgd + 4);
    const auto *sgd_5 = buffer.data(sgd + 5);
    const auto *sgd_6 = buffer.data(sgd + 6);
    const auto *sgd_9 = buffer.data(sgd + 9);
    const auto *sgd_10 = buffer.data(sgd + 10);
    const auto *sgd_11 = buffer.data(sgd + 11);
    const auto *sgd_12 = buffer.data(sgd + 12);
    const auto *sgd_15 = buffer.data(sgd + 15);
    const auto *sgd_16 = buffer.data(sgd + 16);
    const auto *sgd_17 = buffer.data(sgd + 17);
    const auto *sgd_18 = buffer.data(sgd + 18);
    const auto *sgd_21 = buffer.data(sgd + 21);
    const auto *sgd_22 = buffer.data(sgd + 22);
    const auto *sgd_23 = buffer.data(sgd + 23);
    const auto *sgd_24 = buffer.data(sgd + 24);
    const auto *sgd_27 = buffer.data(sgd + 27);
    const auto *sgd_28 = buffer.data(sgd + 28);
    const auto *sgd_29 = buffer.data(sgd + 29);
    const auto *sgd_30 = buffer.data(sgd + 30);
    const auto *sgd_33 = buffer.data(sgd + 33);
    const auto *sgd_34 = buffer.data(sgd + 34);
    const auto *sgd_35 = buffer.data(sgd + 35);
    const auto *sgd_36 = buffer.data(sgd + 36);
    const auto *sgd_39 = buffer.data(sgd + 39);
    const auto *sgd_40 = buffer.data(sgd + 40);
    const auto *sgd_41 = buffer.data(sgd + 41);
    const auto *sgd_42 = buffer.data(sgd + 42);
    const auto *sgd_45 = buffer.data(sgd + 45);
    const auto *sgd_46 = buffer.data(sgd + 46);
    const auto *sgd_47 = buffer.data(sgd + 47);
    const auto *sgd_48 = buffer.data(sgd + 48);
    const auto *sgd_51 = buffer.data(sgd + 51);
    const auto *sgd_52 = buffer.data(sgd + 52);
    const auto *sgd_53 = buffer.data(sgd + 53);
    const auto *sgd_54 = buffer.data(sgd + 54);
    const auto *sgd_57 = buffer.data(sgd + 57);
    const auto *sgd_58 = buffer.data(sgd + 58);
    const auto *sgd_59 = buffer.data(sgd + 59);
    const auto *sgd_60 = buffer.data(sgd + 60);
    const auto *sgd_63 = buffer.data(sgd + 63);
    const auto *sgd_64 = buffer.data(sgd + 64);
    const auto *sgd_65 = buffer.data(sgd + 65);
    const auto *sgd_66 = buffer.data(sgd + 66);
    const auto *sgd_69 = buffer.data(sgd + 69);
    const auto *sgd_70 = buffer.data(sgd + 70);
    const auto *sgd_71 = buffer.data(sgd + 71);
    const auto *sgd_72 = buffer.data(sgd + 72);
    const auto *sgd_75 = buffer.data(sgd + 75);
    const auto *sgd_76 = buffer.data(sgd + 76);
    const auto *sgd_77 = buffer.data(sgd + 77);
    const auto *sgd_78 = buffer.data(sgd + 78);
    const auto *sgd_81 = buffer.data(sgd + 81);
    const auto *sgd_82 = buffer.data(sgd + 82);
    const auto *sgd_83 = buffer.data(sgd + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, sfd_0, sfd_3, sfd_4, \
                         sgp0_0, sgp1_0, sgd_0, sgd_3, sgd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfd_0[k]
                 + f_1 * sgp0_0[k]
                 - f_2 * sgp1_0[k]
                 + f_3 * pc_x[k] * sgd_0[k];

        t_1[k] = f_3 * pc_y[k] * sgd_0[k];

        t_2[k] = f_3 * pc_z[k] * sgd_0[k];

        t_3[k] = f_0 * sfd_3[k]
                 + f_3 * pc_x[k] * sgd_3[k];

        t_4[k] = f_0 * sfd_4[k]
                 + f_3 * pc_x[k] * sgd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, sfd_5, sgp0_1, sgp0_2, \
                         sgp1_1, sgp1_2, sgd_3, sgd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * sfd_5[k]
                 + f_3 * pc_x[k] * sgd_5[k];

        t_6[k] = f_1 * sgp0_1[k]
                 - f_2 * sgp1_1[k]
                 + f_3 * pc_y[k] * sgd_3[k];

        t_7[k] = f_3 * pc_z[k] * sgd_3[k];

        t_8[k] = f_3 * pc_y[k] * sgd_5[k];

        t_9[k] = f_1 * sgp0_2[k]
                 - f_2 * sgp1_2[k]
                 + f_3 * pc_z[k] * sgd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, sff0_0, sfd_0, sfd_9, \
                         sff1_0, sgd_6, sgd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * sff0_0[k]
                  - f_4 * pc_y[k] * sff1_0[k];

        t_11[k] = f_5 * sfd_0[k]
                  + f_3 * pc_y[k] * sgd_6[k];

        t_12[k] = f_3 * pc_z[k] * sgd_6[k];

        t_13[k] = f_6 * sfd_9[k]
                  + f_3 * pc_x[k] * sgd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sfd_3, sfd_10, sfd_11, \
                         sgp0_4, sgp1_4, sgd_9, sgd_10, sgd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * sfd_10[k]
                  + f_3 * pc_x[k] * sgd_10[k];

        t_15[k] = f_6 * sfd_11[k]
                  + f_3 * pc_x[k] * sgd_11[k];

        t_16[k] = f_5 * sfd_3[k]
                  + f_1 * sgp0_4[k]
                  - f_2 * sgp1_4[k]
                  + f_3 * pc_y[k] * sgd_9[k];

        t_17[k] = f_3 * pc_z[k] * sgd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, pc_y, pc_z, sff0_0, sff0_9, \
                         sfd_5, sff1_0, sff1_9, sgd_11, sgd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * sfd_5[k]
                  + f_3 * pc_y[k] * sgd_11[k];

        t_19[k] = pb_y[k] * sff0_9[k]
                  - f_4 * pc_y[k] * sff1_9[k];

        t_20[k] = pb_z[k] * sff0_0[k]
                  - f_4 * pc_z[k] * sff1_0[k];

        t_21[k] = f_3 * pc_y[k] * sgd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, sfd_0, sfd_15, sfd_16, sfd_17, \
                         sgd_12, sgd_15, sgd_16, sgd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * sfd_0[k]
                  + f_3 * pc_z[k] * sgd_12[k];

        t_23[k] = f_6 * sfd_15[k]
                  + f_3 * pc_x[k] * sgd_15[k];

        t_24[k] = f_6 * sfd_16[k]
                  + f_3 * pc_x[k] * sgd_16[k];

        t_25[k] = f_6 * sfd_17[k]
                  + f_3 * pc_x[k] * sgd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, sff0_6, sfd_3, sfd_5, \
                         sff1_6, sgp0_8, sgp1_8, sgd_15, sgd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * sff0_6[k]
                  - f_4 * pc_z[k] * sff1_6[k];

        t_27[k] = f_5 * sfd_3[k]
                  + f_3 * pc_z[k] * sgd_15[k];

        t_28[k] = f_3 * pc_y[k] * sgd_17[k];

        t_29[k] = f_5 * sfd_5[k]
                  + f_1 * sgp0_8[k]
                  - f_2 * sgp1_8[k]
                  + f_3 * pc_z[k] * sgd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pc_x, pc_y, pc_z, sfd_6, sfd_18, sfd_21, \
                         sgp0_9, sgp1_9, sgd_18, sgd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * sfd_18[k]
                  + f_1 * sgp0_9[k]
                  - f_2 * sgp1_9[k]
                  + f_3 * pc_x[k] * sgd_18[k];

        t_31[k] = f_7 * sfd_6[k]
                  + f_3 * pc_y[k] * sgd_18[k];

        t_32[k] = f_3 * pc_z[k] * sgd_18[k];

        t_33[k] = f_7 * sfd_21[k]
                  + f_3 * pc_x[k] * sgd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pc_x, pc_y, pc_z, sfd_9, sfd_22, sfd_23, \
                         sgp0_10, sgp1_10, sgd_21, sgd_22, sgd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * sfd_22[k]
                  + f_3 * pc_x[k] * sgd_22[k];

        t_35[k] = f_7 * sfd_23[k]
                  + f_3 * pc_x[k] * sgd_23[k];

        t_36[k] = f_7 * sfd_9[k]
                  + f_1 * sgp0_10[k]
                  - f_2 * sgp1_10[k]
                  + f_3 * pc_y[k] * sgd_21[k];

        t_37[k] = f_3 * pc_z[k] * sgd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, sff0_20, sfd_11, sfd_12, \
                         sff1_20, sgp0_11, sgp1_11, sgd_23, sgd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_7 * sfd_11[k]
                  + f_3 * pc_y[k] * sgd_23[k];

        t_39[k] = f_1 * sgp0_11[k]
                  - f_2 * sgp1_11[k]
                  + f_3 * pc_z[k] * sgd_23[k];

        t_40[k] = pb_y[k] * sff0_20[k]
                  - f_4 * pc_y[k] * sff1_20[k];

        t_41[k] = f_5 * sfd_12[k]
                  + f_3 * pc_y[k] * sgd_24[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pc_x, pc_z, sfd_6, sfd_27, sfd_28, sfd_29, \
                         sgd_24, sgd_27, sgd_28, sgd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * sfd_6[k]
                  + f_3 * pc_z[k] * sgd_24[k];

        t_43[k] = f_7 * sfd_27[k]
                  + f_3 * pc_x[k] * sgd_27[k];

        t_44[k] = f_7 * sfd_28[k]
                  + f_3 * pc_x[k] * sgd_28[k];

        t_45[k] = f_7 * sfd_29[k]
                  + f_3 * pc_x[k] * sgd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, pc_y, pc_z, sff0_16, sff0_29, \
                         sfd_9, sfd_17, sff1_16, sff1_29, sgd_27, \
                         sgd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * sff0_16[k]
                  - f_4 * pc_z[k] * sff1_16[k];

        t_47[k] = f_5 * sfd_9[k]
                  + f_3 * pc_z[k] * sgd_27[k];

        t_48[k] = f_5 * sfd_17[k]
                  + f_3 * pc_y[k] * sgd_29[k];

        t_49[k] = pb_y[k] * sff0_29[k]
                  - f_4 * pc_y[k] * sff1_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, sfd_12, sfd_30, sfd_33, \
                         sgp0_15, sgp1_15, sgd_30, sgd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_7 * sfd_30[k]
                  + f_1 * sgp0_15[k]
                  - f_2 * sgp1_15[k]
                  + f_3 * pc_x[k] * sgd_30[k];

        t_51[k] = f_3 * pc_y[k] * sgd_30[k];

        t_52[k] = f_7 * sfd_12[k]
                  + f_3 * pc_z[k] * sgd_30[k];

        t_53[k] = f_7 * sfd_33[k]
                  + f_3 * pc_x[k] * sgd_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, sfd_15, sfd_34, \
                         sfd_35, sgp0_16, sgp1_16, sgd_33, sgd_34, \
                         sgd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * sfd_34[k]
                  + f_3 * pc_x[k] * sgd_34[k];

        t_55[k] = f_7 * sfd_35[k]
                  + f_3 * pc_x[k] * sgd_35[k];

        t_56[k] = f_1 * sgp0_16[k]
                  - f_2 * sgp1_16[k]
                  + f_3 * pc_y[k] * sgd_33[k];

        t_57[k] = f_7 * sfd_15[k]
                  + f_3 * pc_z[k] * sgd_33[k];

        t_58[k] = f_3 * pc_y[k] * sgd_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_x, pc_x, pc_y, pc_z, sff0_60, sfd_17, sfd_18, \
                         sfd_36, sff1_60, sgp0_17, sgp1_17, sgd_35, \
                         sgd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_7 * sfd_17[k]
                  + f_1 * sgp0_17[k]
                  - f_2 * sgp1_17[k]
                  + f_3 * pc_z[k] * sgd_35[k];

        t_60[k] = pb_x[k] * sff0_60[k]
                  + f_6 * sfd_36[k]
                  - f_4 * pc_x[k] * sff1_60[k];

        t_61[k] = f_6 * sfd_18[k]
                  + f_3 * pc_y[k] * sgd_36[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_z, sfd_39, sfd_40, sfd_41, sgd_36, \
                         sgd_39, sgd_40, sgd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_3 * pc_z[k] * sgd_36[k];

        t_63[k] = f_5 * sfd_39[k]
                  + f_3 * pc_x[k] * sgd_39[k];

        t_64[k] = f_5 * sfd_40[k]
                  + f_3 * pc_x[k] * sgd_40[k];

        t_65[k] = f_5 * sfd_41[k]
                  + f_3 * pc_x[k] * sgd_41[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_x, pc_x, pc_y, pc_z, sff0_66, sff0_69, \
                         sfd_23, sff1_66, sff1_69, sgd_39, sgd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_x[k] * sff0_66[k]
                  - f_4 * pc_x[k] * sff1_66[k];

        t_67[k] = f_3 * pc_z[k] * sgd_39[k];

        t_68[k] = f_6 * sfd_23[k]
                  + f_3 * pc_y[k] * sgd_41[k];

        t_69[k] = pb_x[k] * sff0_69[k]
                  - f_4 * pc_x[k] * sff1_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_z, pc_x, pc_y, pc_z, sff0_30, sfd_18, \
                         sfd_24, sfd_45, sff1_30, sgd_42, sgd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_z[k] * sff0_30[k]
                  - f_4 * pc_z[k] * sff1_30[k];

        t_71[k] = f_7 * sfd_24[k]
                  + f_3 * pc_y[k] * sgd_42[k];

        t_72[k] = f_5 * sfd_18[k]
                  + f_3 * pc_z[k] * sgd_42[k];

        t_73[k] = f_5 * sfd_45[k]
                  + f_3 * pc_x[k] * sgd_45[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_x, pc_x, pc_z, sff0_76, sfd_21, sfd_46, \
                         sfd_47, sff1_76, sgd_45, sgd_46, sgd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_5 * sfd_46[k]
                  + f_3 * pc_x[k] * sgd_46[k];

        t_75[k] = f_5 * sfd_47[k]
                  + f_3 * pc_x[k] * sgd_47[k];

        t_76[k] = pb_x[k] * sff0_76[k]
                  - f_4 * pc_x[k] * sff1_76[k];

        t_77[k] = f_5 * sfd_21[k]
                  + f_3 * pc_z[k] * sgd_45[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_x, pb_y, pc_x, pc_y, sff0_50, sff0_79, \
                         sfd_29, sfd_30, sff1_50, sff1_79, sgd_47, \
                         sgd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_7 * sfd_29[k]
                  + f_3 * pc_y[k] * sgd_47[k];

        t_79[k] = pb_x[k] * sff0_79[k]
                  - f_4 * pc_x[k] * sff1_79[k];

        t_80[k] = pb_y[k] * sff0_50[k]
                  - f_4 * pc_y[k] * sff1_50[k];

        t_81[k] = f_5 * sfd_30[k]
                  + f_3 * pc_y[k] * sgd_48[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pc_x, pc_z, sfd_24, sfd_51, sfd_52, sfd_53, \
                         sgd_48, sgd_51, sgd_52, sgd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_7 * sfd_24[k]
                  + f_3 * pc_z[k] * sgd_48[k];

        t_83[k] = f_5 * sfd_51[k]
                  + f_3 * pc_x[k] * sgd_51[k];

        t_84[k] = f_5 * sfd_52[k]
                  + f_3 * pc_x[k] * sgd_52[k];

        t_85[k] = f_5 * sfd_53[k]
                  + f_3 * pc_x[k] * sgd_53[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_x, pc_x, pc_y, pc_z, sff0_86, sff0_89, \
                         sfd_27, sfd_35, sff1_86, sff1_89, sgd_51, \
                         sgd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_x[k] * sff0_86[k]
                  - f_4 * pc_x[k] * sff1_86[k];

        t_87[k] = f_7 * sfd_27[k]
                  + f_3 * pc_z[k] * sgd_51[k];

        t_88[k] = f_5 * sfd_35[k]
                  + f_3 * pc_y[k] * sgd_53[k];

        t_89[k] = pb_x[k] * sff0_89[k]
                  - f_4 * pc_x[k] * sff1_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_x, pc_x, pc_y, pc_z, sff0_90, sfd_30, \
                         sfd_54, sfd_57, sff1_90, sgd_54, sgd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_x[k] * sff0_90[k]
                  + f_6 * sfd_54[k]
                  - f_4 * pc_x[k] * sff1_90[k];

        t_91[k] = f_3 * pc_y[k] * sgd_54[k];

        t_92[k] = f_6 * sfd_30[k]
                  + f_3 * pc_z[k] * sgd_54[k];

        t_93[k] = f_5 * sfd_57[k]
                  + f_3 * pc_x[k] * sgd_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_x, pc_x, pc_z, sff0_96, sfd_33, sfd_58, \
                         sfd_59, sff1_96, sgd_57, sgd_58, sgd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * sfd_58[k]
                  + f_3 * pc_x[k] * sgd_58[k];

        t_95[k] = f_5 * sfd_59[k]
                  + f_3 * pc_x[k] * sgd_59[k];

        t_96[k] = pb_x[k] * sff0_96[k]
                  - f_4 * pc_x[k] * sff1_96[k];

        t_97[k] = f_6 * sfd_33[k]
                  + f_3 * pc_z[k] * sgd_57[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pb_x, pc_x, pc_y, pc_z, sff0_99, \
                         sfd_36, sff1_99, sgp0_30, sgp1_30, sgd_59, \
                         sgd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_3 * pc_y[k] * sgd_59[k];

        t_99[k] = pb_x[k] * sff0_99[k]
                  - f_4 * pc_x[k] * sff1_99[k];

        t_100[k] = f_1 * sgp0_30[k]
                   - f_2 * sgp1_30[k]
                   + f_3 * pc_x[k] * sgd_60[k];

        t_101[k] = f_0 * sfd_36[k]
                   + f_3 * pc_y[k] * sgd_60[k];

        t_102[k] = f_3 * pc_z[k] * sgd_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, t_108, pc_x, pc_y, pc_z, sfd_39, \
                         sfd_41, sgp0_31, sgp1_31, sgd_63, sgd_64, \
                         sgd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_3 * pc_x[k] * sgd_63[k];

        t_104[k] = f_3 * pc_x[k] * sgd_64[k];

        t_105[k] = f_3 * pc_x[k] * sgd_65[k];

        t_106[k] = f_0 * sfd_39[k]
                   + f_1 * sgp0_31[k]
                   - f_2 * sgp1_31[k]
                   + f_3 * pc_y[k] * sgd_63[k];

        t_107[k] = f_3 * pc_z[k] * sgd_63[k];

        t_108[k] = f_0 * sfd_41[k]
                   + f_3 * pc_y[k] * sgd_65[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pb_z, pc_y, pc_z, sff0_60, sfd_36, \
                         sfd_42, sff1_60, sgp0_32, sgp1_32, sgd_65, \
                         sgd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_1 * sgp0_32[k]
                   - f_2 * sgp1_32[k]
                   + f_3 * pc_z[k] * sgd_65[k];

        t_110[k] = pb_z[k] * sff0_60[k]
                   - f_4 * pc_z[k] * sff1_60[k];

        t_111[k] = f_6 * sfd_42[k]
                   + f_3 * pc_y[k] * sgd_66[k];

        t_112[k] = f_5 * sfd_36[k]
                   + f_3 * pc_z[k] * sgd_66[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pb_z, pc_x, pc_z, sff0_66, sfd_39, \
                         sff1_66, sgd_69, sgd_70, sgd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_3 * pc_x[k] * sgd_69[k];

        t_114[k] = f_3 * pc_x[k] * sgd_70[k];

        t_115[k] = f_3 * pc_x[k] * sgd_71[k];

        t_116[k] = pb_z[k] * sff0_66[k]
                   - f_4 * pc_z[k] * sff1_66[k];

        t_117[k] = f_5 * sfd_39[k]
                   + f_3 * pc_z[k] * sgd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, pc_z, sfd_41, sfd_47, sfd_48, \
                         sgp0_35, sgp0_36, sgp1_35, sgp1_36, sgd_71, \
                         sgd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_6 * sfd_47[k]
                   + f_3 * pc_y[k] * sgd_71[k];

        t_119[k] = f_5 * sfd_41[k]
                   + f_1 * sgp0_35[k]
                   - f_2 * sgp1_35[k]
                   + f_3 * pc_z[k] * sgd_71[k];

        t_120[k] = f_1 * sgp0_36[k]
                   - f_2 * sgp1_36[k]
                   + f_3 * pc_x[k] * sgd_72[k];

        t_121[k] = f_7 * sfd_48[k]
                   + f_3 * pc_y[k] * sgd_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, pc_x, pc_y, pc_z, sfd_42, sfd_51, \
                         sgp0_37, sgp1_37, sgd_72, sgd_75, sgd_76, \
                         sgd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_7 * sfd_42[k]
                   + f_3 * pc_z[k] * sgd_72[k];

        t_123[k] = f_3 * pc_x[k] * sgd_75[k];

        t_124[k] = f_3 * pc_x[k] * sgd_76[k];

        t_125[k] = f_3 * pc_x[k] * sgd_77[k];

        t_126[k] = f_7 * sfd_51[k]
                   + f_1 * sgp0_37[k]
                   - f_2 * sgp1_37[k]
                   + f_3 * pc_y[k] * sgd_75[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pb_y, pc_y, pc_z, sff0_90, sfd_45, \
                         sfd_47, sfd_53, sff1_90, sgp0_38, sgp1_38, sgd_75, \
                         sgd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_7 * sfd_45[k]
                   + f_3 * pc_z[k] * sgd_75[k];

        t_128[k] = f_7 * sfd_53[k]
                   + f_3 * pc_y[k] * sgd_77[k];

        t_129[k] = f_7 * sfd_47[k]
                   + f_1 * sgp0_38[k]
                   - f_2 * sgp1_38[k]
                   + f_3 * pc_z[k] * sgd_77[k];

        t_130[k] = pb_y[k] * sff0_90[k]
                   - f_4 * pc_y[k] * sff1_90[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, pc_x, pc_y, pc_z, sfd_48, sfd_54, \
                         sgd_78, sgd_81, sgd_82, sgd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_5 * sfd_54[k]
                   + f_3 * pc_y[k] * sgd_78[k];

        t_132[k] = f_6 * sfd_48[k]
                   + f_3 * pc_z[k] * sgd_78[k];

        t_133[k] = f_3 * pc_x[k] * sgd_81[k];

        t_134[k] = f_3 * pc_x[k] * sgd_82[k];

        t_135[k] = f_3 * pc_x[k] * sgd_83[k];
    }
}

static auto
compute_prim_sgf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sff0,
                                                          const size_t sfd, const size_t sff1,
                                                          const size_t sgp0, const size_t sgp1,
                                                          const size_t sgd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sff0_96 = buffer.data(sff0 + 96);
    const auto *sff0_99 = buffer.data(sff0 + 99);

    const auto *sfd_51 = buffer.data(sfd + 51);
    const auto *sfd_54 = buffer.data(sfd + 54);
    const auto *sfd_57 = buffer.data(sfd + 57);
    const auto *sfd_59 = buffer.data(sfd + 59);

    const auto *sff1_96 = buffer.data(sff1 + 96);
    const auto *sff1_99 = buffer.data(sff1 + 99);

    const auto *sgp0_42 = buffer.data(sgp0 + 42);
    const auto *sgp0_43 = buffer.data(sgp0 + 43);
    const auto *sgp0_44 = buffer.data(sgp0 + 44);

    const auto *sgp1_42 = buffer.data(sgp1 + 42);
    const auto *sgp1_43 = buffer.data(sgp1 + 43);
    const auto *sgp1_44 = buffer.data(sgp1 + 44);

    const auto *sgd_81 = buffer.data(sgd + 81);
    const auto *sgd_83 = buffer.data(sgd + 83);
    const auto *sgd_84 = buffer.data(sgd + 84);
    const auto *sgd_87 = buffer.data(sgd + 87);
    const auto *sgd_88 = buffer.data(sgd + 88);
    const auto *sgd_89 = buffer.data(sgd + 89);

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pb_y, pc_y, pc_z, sff0_96, sff0_99, \
                         sfd_51, sfd_57, sfd_59, sff1_96, sff1_99, sgd_81, \
                         sgd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_y[k] * sff0_96[k]
                   + f_6 * sfd_57[k]
                   - f_4 * pc_y[k] * sff1_96[k];

        t_137[k] = f_6 * sfd_51[k]
                   + f_3 * pc_z[k] * sgd_81[k];

        t_138[k] = f_5 * sfd_59[k]
                   + f_3 * pc_y[k] * sgd_83[k];

        t_139[k] = pb_y[k] * sff0_99[k]
                   - f_4 * pc_y[k] * sff1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, t_145, pc_x, pc_y, pc_z, sfd_54, \
                         sgp0_42, sgp1_42, sgd_84, sgd_87, sgd_88, \
                         sgd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_1 * sgp0_42[k]
                   - f_2 * sgp1_42[k]
                   + f_3 * pc_x[k] * sgd_84[k];

        t_141[k] = f_3 * pc_y[k] * sgd_84[k];

        t_142[k] = f_0 * sfd_54[k]
                   + f_3 * pc_z[k] * sgd_84[k];

        t_143[k] = f_3 * pc_x[k] * sgd_87[k];

        t_144[k] = f_3 * pc_x[k] * sgd_88[k];

        t_145[k] = f_3 * pc_x[k] * sgd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_y, pc_z, sfd_57, sfd_59, sgp0_43, \
                         sgp0_44, sgp1_43, sgp1_44, sgd_87, sgd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * sgp0_43[k]
                   - f_2 * sgp1_43[k]
                   + f_3 * pc_y[k] * sgd_87[k];

        t_147[k] = f_0 * sfd_57[k]
                   + f_3 * pc_z[k] * sgd_87[k];

        t_148[k] = f_3 * pc_y[k] * sgd_89[k];

        t_149[k] = f_0 * sfd_59[k]
                   + f_1 * sgp0_44[k]
                   - f_2 * sgp1_44[k]
                   + f_3 * pc_z[k] * sgd_89[k];
    }
}

auto
compute_prim_sgf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sff0, const size_t sfd,
                                                   const size_t sff1, const size_t sgp0,
                                                   const size_t sgp1, const size_t sgd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sgf_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sff0, sfd,
                                                              sff1, sgp0, sgp1, sgd, ncols,
                                                              gamma, p, q);

    compute_prim_sgf_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sff0, sfd,
                                                              sff1, sgp0, sgp1, sgd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
