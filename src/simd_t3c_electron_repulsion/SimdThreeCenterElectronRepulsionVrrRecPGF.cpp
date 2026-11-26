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


#include "SimdThreeCenterElectronRepulsionVrrRecPGF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pgf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgf0, const size_t sgd,
                                                          const size_t sgf1, const size_t pdf0,
                                                          const size_t pdf1, const size_t pff0,
                                                          const size_t pfd, const size_t pff1,
                                                          const size_t pgp0, const size_t pgp1,
                                                          const size_t pgd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);
    const auto f_9 = 1.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgf0_100 = buffer.data(sgf0 + 100);
    const auto *sgf0_106 = buffer.data(sgf0 + 106);
    const auto *sgf0_109 = buffer.data(sgf0 + 109);
    const auto *sgf0_116 = buffer.data(sgf0 + 116);
    const auto *sgf0_119 = buffer.data(sgf0 + 119);
    const auto *sgf0_120 = buffer.data(sgf0 + 120);

    const auto *sgd_0 = buffer.data(sgd + 0);
    const auto *sgd_3 = buffer.data(sgd + 3);
    const auto *sgd_5 = buffer.data(sgd + 5);
    const auto *sgd_9 = buffer.data(sgd + 9);
    const auto *sgd_17 = buffer.data(sgd + 17);
    const auto *sgd_21 = buffer.data(sgd + 21);
    const auto *sgd_23 = buffer.data(sgd + 23);
    const auto *sgd_33 = buffer.data(sgd + 33);
    const auto *sgd_35 = buffer.data(sgd + 35);
    const auto *sgd_36 = buffer.data(sgd + 36);
    const auto *sgd_39 = buffer.data(sgd + 39);
    const auto *sgd_41 = buffer.data(sgd + 41);
    const auto *sgd_47 = buffer.data(sgd + 47);
    const auto *sgd_51 = buffer.data(sgd + 51);
    const auto *sgd_54 = buffer.data(sgd + 54);
    const auto *sgd_57 = buffer.data(sgd + 57);
    const auto *sgd_59 = buffer.data(sgd + 59);
    const auto *sgd_60 = buffer.data(sgd + 60);
    const auto *sgd_63 = buffer.data(sgd + 63);
    const auto *sgd_65 = buffer.data(sgd + 65);
    const auto *sgd_69 = buffer.data(sgd + 69);
    const auto *sgd_71 = buffer.data(sgd + 71);
    const auto *sgd_72 = buffer.data(sgd + 72);

    const auto *sgf1_100 = buffer.data(sgf1 + 100);
    const auto *sgf1_106 = buffer.data(sgf1 + 106);
    const auto *sgf1_109 = buffer.data(sgf1 + 109);
    const auto *sgf1_116 = buffer.data(sgf1 + 116);
    const auto *sgf1_119 = buffer.data(sgf1 + 119);
    const auto *sgf1_120 = buffer.data(sgf1 + 120);

    const auto *pdf0_0 = buffer.data(pdf0 + 0);

    const auto *pdf1_0 = buffer.data(pdf1 + 0);

    const auto *pff0_0 = buffer.data(pff0 + 0);
    const auto *pff0_3 = buffer.data(pff0 + 3);
    const auto *pff0_5 = buffer.data(pff0 + 5);
    const auto *pff0_6 = buffer.data(pff0 + 6);
    const auto *pff0_9 = buffer.data(pff0 + 9);
    const auto *pff0_10 = buffer.data(pff0 + 10);
    const auto *pff0_13 = buffer.data(pff0 + 13);
    const auto *pff0_16 = buffer.data(pff0 + 16);
    const auto *pff0_20 = buffer.data(pff0 + 20);
    const auto *pff0_25 = buffer.data(pff0 + 25);
    const auto *pff0_29 = buffer.data(pff0 + 29);
    const auto *pff0_30 = buffer.data(pff0 + 30);
    const auto *pff0_33 = buffer.data(pff0 + 33);
    const auto *pff0_50 = buffer.data(pff0 + 50);
    const auto *pff0_55 = buffer.data(pff0 + 55);
    const auto *pff0_60 = buffer.data(pff0 + 60);

    const auto *pfd_0 = buffer.data(pfd + 0);
    const auto *pfd_2 = buffer.data(pfd + 2);
    const auto *pfd_3 = buffer.data(pfd + 3);
    const auto *pfd_5 = buffer.data(pfd + 5);
    const auto *pfd_6 = buffer.data(pfd + 6);
    const auto *pfd_8 = buffer.data(pfd + 8);
    const auto *pfd_9 = buffer.data(pfd + 9);
    const auto *pfd_11 = buffer.data(pfd + 11);
    const auto *pfd_12 = buffer.data(pfd + 12);
    const auto *pfd_14 = buffer.data(pfd + 14);
    const auto *pfd_15 = buffer.data(pfd + 15);
    const auto *pfd_17 = buffer.data(pfd + 17);
    const auto *pfd_18 = buffer.data(pfd + 18);
    const auto *pfd_20 = buffer.data(pfd + 20);
    const auto *pfd_21 = buffer.data(pfd + 21);
    const auto *pfd_23 = buffer.data(pfd + 23);
    const auto *pfd_24 = buffer.data(pfd + 24);
    const auto *pfd_26 = buffer.data(pfd + 26);
    const auto *pfd_27 = buffer.data(pfd + 27);
    const auto *pfd_29 = buffer.data(pfd + 29);
    const auto *pfd_30 = buffer.data(pfd + 30);
    const auto *pfd_32 = buffer.data(pfd + 32);
    const auto *pfd_33 = buffer.data(pfd + 33);
    const auto *pfd_35 = buffer.data(pfd + 35);
    const auto *pfd_36 = buffer.data(pfd + 36);
    const auto *pfd_38 = buffer.data(pfd + 38);
    const auto *pfd_39 = buffer.data(pfd + 39);
    const auto *pfd_41 = buffer.data(pfd + 41);
    const auto *pfd_42 = buffer.data(pfd + 42);
    const auto *pfd_44 = buffer.data(pfd + 44);
    const auto *pfd_47 = buffer.data(pfd + 47);
    const auto *pfd_48 = buffer.data(pfd + 48);
    const auto *pfd_51 = buffer.data(pfd + 51);
    const auto *pfd_54 = buffer.data(pfd + 54);
    const auto *pfd_57 = buffer.data(pfd + 57);
    const auto *pfd_59 = buffer.data(pfd + 59);

    const auto *pff1_0 = buffer.data(pff1 + 0);
    const auto *pff1_3 = buffer.data(pff1 + 3);
    const auto *pff1_5 = buffer.data(pff1 + 5);
    const auto *pff1_6 = buffer.data(pff1 + 6);
    const auto *pff1_9 = buffer.data(pff1 + 9);
    const auto *pff1_10 = buffer.data(pff1 + 10);
    const auto *pff1_13 = buffer.data(pff1 + 13);
    const auto *pff1_16 = buffer.data(pff1 + 16);
    const auto *pff1_20 = buffer.data(pff1 + 20);
    const auto *pff1_25 = buffer.data(pff1 + 25);
    const auto *pff1_29 = buffer.data(pff1 + 29);
    const auto *pff1_30 = buffer.data(pff1 + 30);
    const auto *pff1_33 = buffer.data(pff1 + 33);
    const auto *pff1_50 = buffer.data(pff1 + 50);
    const auto *pff1_55 = buffer.data(pff1 + 55);
    const auto *pff1_60 = buffer.data(pff1 + 60);

    const auto *pgp0_0 = buffer.data(pgp0 + 0);
    const auto *pgp0_1 = buffer.data(pgp0 + 1);
    const auto *pgp0_2 = buffer.data(pgp0 + 2);
    const auto *pgp0_4 = buffer.data(pgp0 + 4);
    const auto *pgp0_8 = buffer.data(pgp0 + 8);
    const auto *pgp0_10 = buffer.data(pgp0 + 10);
    const auto *pgp0_11 = buffer.data(pgp0 + 11);
    const auto *pgp0_16 = buffer.data(pgp0 + 16);
    const auto *pgp0_17 = buffer.data(pgp0 + 17);
    const auto *pgp0_18 = buffer.data(pgp0 + 18);
    const auto *pgp0_19 = buffer.data(pgp0 + 19);
    const auto *pgp0_20 = buffer.data(pgp0 + 20);
    const auto *pgp0_22 = buffer.data(pgp0 + 22);
    const auto *pgp0_23 = buffer.data(pgp0 + 23);
    const auto *pgp0_25 = buffer.data(pgp0 + 25);
    const auto *pgp0_26 = buffer.data(pgp0 + 26);
    const auto *pgp0_27 = buffer.data(pgp0 + 27);
    const auto *pgp0_28 = buffer.data(pgp0 + 28);
    const auto *pgp0_29 = buffer.data(pgp0 + 29);

    const auto *pgp1_0 = buffer.data(pgp1 + 0);
    const auto *pgp1_1 = buffer.data(pgp1 + 1);
    const auto *pgp1_2 = buffer.data(pgp1 + 2);
    const auto *pgp1_4 = buffer.data(pgp1 + 4);
    const auto *pgp1_8 = buffer.data(pgp1 + 8);
    const auto *pgp1_10 = buffer.data(pgp1 + 10);
    const auto *pgp1_11 = buffer.data(pgp1 + 11);
    const auto *pgp1_16 = buffer.data(pgp1 + 16);
    const auto *pgp1_17 = buffer.data(pgp1 + 17);
    const auto *pgp1_18 = buffer.data(pgp1 + 18);
    const auto *pgp1_19 = buffer.data(pgp1 + 19);
    const auto *pgp1_20 = buffer.data(pgp1 + 20);
    const auto *pgp1_22 = buffer.data(pgp1 + 22);
    const auto *pgp1_23 = buffer.data(pgp1 + 23);
    const auto *pgp1_25 = buffer.data(pgp1 + 25);
    const auto *pgp1_26 = buffer.data(pgp1 + 26);
    const auto *pgp1_27 = buffer.data(pgp1 + 27);
    const auto *pgp1_28 = buffer.data(pgp1 + 28);
    const auto *pgp1_29 = buffer.data(pgp1 + 29);

    const auto *pgd_0 = buffer.data(pgd + 0);
    const auto *pgd_2 = buffer.data(pgd + 2);
    const auto *pgd_3 = buffer.data(pgd + 3);
    const auto *pgd_5 = buffer.data(pgd + 5);
    const auto *pgd_6 = buffer.data(pgd + 6);
    const auto *pgd_8 = buffer.data(pgd + 8);
    const auto *pgd_9 = buffer.data(pgd + 9);
    const auto *pgd_11 = buffer.data(pgd + 11);
    const auto *pgd_12 = buffer.data(pgd + 12);
    const auto *pgd_14 = buffer.data(pgd + 14);
    const auto *pgd_15 = buffer.data(pgd + 15);
    const auto *pgd_17 = buffer.data(pgd + 17);
    const auto *pgd_18 = buffer.data(pgd + 18);
    const auto *pgd_20 = buffer.data(pgd + 20);
    const auto *pgd_21 = buffer.data(pgd + 21);
    const auto *pgd_23 = buffer.data(pgd + 23);
    const auto *pgd_24 = buffer.data(pgd + 24);
    const auto *pgd_26 = buffer.data(pgd + 26);
    const auto *pgd_27 = buffer.data(pgd + 27);
    const auto *pgd_29 = buffer.data(pgd + 29);
    const auto *pgd_30 = buffer.data(pgd + 30);
    const auto *pgd_32 = buffer.data(pgd + 32);
    const auto *pgd_33 = buffer.data(pgd + 33);
    const auto *pgd_35 = buffer.data(pgd + 35);
    const auto *pgd_36 = buffer.data(pgd + 36);
    const auto *pgd_38 = buffer.data(pgd + 38);
    const auto *pgd_39 = buffer.data(pgd + 39);
    const auto *pgd_41 = buffer.data(pgd + 41);
    const auto *pgd_42 = buffer.data(pgd + 42);
    const auto *pgd_44 = buffer.data(pgd + 44);
    const auto *pgd_45 = buffer.data(pgd + 45);
    const auto *pgd_47 = buffer.data(pgd + 47);
    const auto *pgd_48 = buffer.data(pgd + 48);
    const auto *pgd_50 = buffer.data(pgd + 50);
    const auto *pgd_51 = buffer.data(pgd + 51);
    const auto *pgd_53 = buffer.data(pgd + 53);
    const auto *pgd_54 = buffer.data(pgd + 54);
    const auto *pgd_56 = buffer.data(pgd + 56);
    const auto *pgd_57 = buffer.data(pgd + 57);
    const auto *pgd_59 = buffer.data(pgd + 59);
    const auto *pgd_60 = buffer.data(pgd + 60);
    const auto *pgd_62 = buffer.data(pgd + 62);
    const auto *pgd_63 = buffer.data(pgd + 63);
    const auto *pgd_65 = buffer.data(pgd + 65);
    const auto *pgd_66 = buffer.data(pgd + 66);
    const auto *pgd_68 = buffer.data(pgd + 68);
    const auto *pgd_69 = buffer.data(pgd + 69);
    const auto *pgd_71 = buffer.data(pgd + 71);
    const auto *pgd_72 = buffer.data(pgd + 72);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sgd_0, sgd_3, pfd_0, pfd_3, \
                         pgp0_0, pgp1_0, pgd_0, pgd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgd_0[k]
                 + f_1 * pfd_0[k]
                 + f_2 * pgp0_0[k]
                 - f_3 * pgp1_0[k]
                 + f_4 * pc_x[k] * pgd_0[k];

        t_1[k] = f_4 * pc_y[k] * pgd_0[k];

        t_2[k] = f_4 * pc_z[k] * pgd_0[k];

        t_3[k] = f_0 * sgd_3[k]
                 + f_1 * pfd_3[k]
                 + f_4 * pc_x[k] * pgd_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, sgd_5, pfd_5, pgp0_1, \
                         pgp1_1, pgd_2, pgd_3, pgd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * pc_y[k] * pgd_2[k];

        t_5[k] = f_0 * sgd_5[k]
                 + f_1 * pfd_5[k]
                 + f_4 * pc_x[k] * pgd_5[k];

        t_6[k] = f_2 * pgp0_1[k]
                 - f_3 * pgp1_1[k]
                 + f_4 * pc_y[k] * pgd_3[k];

        t_7[k] = f_4 * pc_z[k] * pgd_3[k];

        t_8[k] = f_4 * pc_y[k] * pgd_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pc_y, pc_z, pff0_0, pfd_0, pff1_0, \
                         pgp0_2, pgp1_2, pgd_5, pgd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * pgp0_2[k]
                 - f_3 * pgp1_2[k]
                 + f_4 * pc_z[k] * pgd_5[k];

        t_10[k] = pb_y[k] * pff0_0[k]
                  - f_5 * pc_y[k] * pff1_0[k];

        t_11[k] = f_0 * pfd_0[k]
                  + f_4 * pc_y[k] * pgd_6[k];

        t_12[k] = f_4 * pc_z[k] * pgd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_y, pc_x, pc_y, sgd_9, pff0_5, pfd_2, pfd_9, \
                         pff1_5, pgd_8, pgd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * sgd_9[k]
                  + f_6 * pfd_9[k]
                  + f_4 * pc_x[k] * pgd_9[k];

        t_14[k] = f_0 * pfd_2[k]
                  + f_4 * pc_y[k] * pgd_8[k];

        t_15[k] = pb_y[k] * pff0_5[k]
                  - f_5 * pc_y[k] * pff1_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, pc_y, pc_z, pff0_9, pfd_3, pfd_5, \
                         pff1_9, pgp0_4, pgp1_4, pgd_9, pgd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pfd_3[k]
                  + f_2 * pgp0_4[k]
                  - f_3 * pgp1_4[k]
                  + f_4 * pc_y[k] * pgd_9[k];

        t_17[k] = f_4 * pc_z[k] * pgd_9[k];

        t_18[k] = f_0 * pfd_5[k]
                  + f_4 * pc_y[k] * pgd_11[k];

        t_19[k] = pb_y[k] * pff0_9[k]
                  - f_5 * pc_y[k] * pff1_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pb_z, pc_y, pc_z, pff0_0, pff0_3, \
                         pfd_0, pff1_0, pff1_3, pgd_12, pgd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_z[k] * pff0_0[k]
                  - f_5 * pc_z[k] * pff1_0[k];

        t_21[k] = f_4 * pc_y[k] * pgd_12[k];

        t_22[k] = f_0 * pfd_0[k]
                  + f_4 * pc_z[k] * pgd_12[k];

        t_23[k] = pb_z[k] * pff0_3[k]
                  - f_5 * pc_z[k] * pff1_3[k];

        t_24[k] = f_4 * pc_y[k] * pgd_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_z, pc_x, pc_y, pc_z, sgd_17, pff0_6, \
                         pfd_3, pfd_17, pff1_6, pgd_15, pgd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * sgd_17[k]
                  + f_6 * pfd_17[k]
                  + f_4 * pc_x[k] * pgd_17[k];

        t_26[k] = pb_z[k] * pff0_6[k]
                  - f_5 * pc_z[k] * pff1_6[k];

        t_27[k] = f_0 * pfd_3[k]
                  + f_4 * pc_z[k] * pgd_15[k];

        t_28[k] = f_4 * pc_y[k] * pgd_17[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_y, pc_y, pc_z, pdf0_0, pdf1_0, pff0_10, pfd_5, \
                         pfd_6, pff1_10, pgp0_8, pgp1_8, pgd_17, \
                         pgd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * pfd_5[k]
                  + f_2 * pgp0_8[k]
                  - f_3 * pgp1_8[k]
                  + f_4 * pc_z[k] * pgd_17[k];

        t_30[k] = f_7 * pdf0_0[k]
                  - f_8 * pdf1_0[k]
                  + pb_y[k] * pff0_10[k]
                  - f_5 * pc_y[k] * pff1_10[k];

        t_31[k] = f_9 * pfd_6[k]
                  + f_4 * pc_y[k] * pgd_18[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, sgd_21, sgd_23, pfd_8, \
                         pfd_21, pfd_23, pgd_18, pgd_20, pgd_21, \
                         pgd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_4 * pc_z[k] * pgd_18[k];

        t_33[k] = f_0 * sgd_21[k]
                  + f_9 * pfd_21[k]
                  + f_4 * pc_x[k] * pgd_21[k];

        t_34[k] = f_9 * pfd_8[k]
                  + f_4 * pc_y[k] * pgd_20[k];

        t_35[k] = f_0 * sgd_23[k]
                  + f_9 * pfd_23[k]
                  + f_4 * pc_x[k] * pgd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, pfd_9, pfd_11, pgp0_10, pgp0_11, \
                         pgp1_10, pgp1_11, pgd_21, pgd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * pfd_9[k]
                  + f_2 * pgp0_10[k]
                  - f_3 * pgp1_10[k]
                  + f_4 * pc_y[k] * pgd_21[k];

        t_37[k] = f_4 * pc_z[k] * pgd_21[k];

        t_38[k] = f_9 * pfd_11[k]
                  + f_4 * pc_y[k] * pgd_23[k];

        t_39[k] = f_2 * pgp0_11[k]
                  - f_3 * pgp1_11[k]
                  + f_4 * pc_z[k] * pgd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pb_z, pc_y, pc_z, pff0_13, pff0_20, \
                         pfd_6, pfd_12, pff1_13, pff1_20, pgd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * pff0_20[k]
                  - f_5 * pc_y[k] * pff1_20[k];

        t_41[k] = f_0 * pfd_12[k]
                  + f_4 * pc_y[k] * pgd_24[k];

        t_42[k] = f_0 * pfd_6[k]
                  + f_4 * pc_z[k] * pgd_24[k];

        t_43[k] = pb_z[k] * pff0_13[k]
                  - f_5 * pc_z[k] * pff1_13[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pb_z, pc_y, pc_z, pff0_16, pff0_25, \
                         pfd_9, pfd_14, pff1_16, pff1_25, pgd_26, \
                         pgd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * pfd_14[k]
                  + f_4 * pc_y[k] * pgd_26[k];

        t_45[k] = pb_y[k] * pff0_25[k]
                  - f_5 * pc_y[k] * pff1_25[k];

        t_46[k] = pb_z[k] * pff0_16[k]
                  - f_5 * pc_z[k] * pff1_16[k];

        t_47[k] = f_0 * pfd_9[k]
                  + f_4 * pc_z[k] * pgd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_y, pb_z, pc_y, pc_z, pdf0_0, pdf1_0, pff0_20, \
                         pff0_29, pfd_17, pff1_20, pff1_29, pgd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * pfd_17[k]
                  + f_4 * pc_y[k] * pgd_29[k];

        t_49[k] = pb_y[k] * pff0_29[k]
                  - f_5 * pc_y[k] * pff1_29[k];

        t_50[k] = f_7 * pdf0_0[k]
                  - f_8 * pdf1_0[k]
                  + pb_z[k] * pff0_20[k]
                  - f_5 * pc_z[k] * pff1_20[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pc_x, pc_y, pc_z, sgd_33, pfd_12, pfd_33, \
                         pgd_30, pgd_32, pgd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_4 * pc_y[k] * pgd_30[k];

        t_52[k] = f_9 * pfd_12[k]
                  + f_4 * pc_z[k] * pgd_30[k];

        t_53[k] = f_0 * sgd_33[k]
                  + f_9 * pfd_33[k]
                  + f_4 * pc_x[k] * pgd_33[k];

        t_54[k] = f_4 * pc_y[k] * pgd_32[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, sgd_35, pfd_15, pfd_35, \
                         pgp0_16, pgp1_16, pgd_33, pgd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * sgd_35[k]
                  + f_9 * pfd_35[k]
                  + f_4 * pc_x[k] * pgd_35[k];

        t_56[k] = f_2 * pgp0_16[k]
                  - f_3 * pgp1_16[k]
                  + f_4 * pc_y[k] * pgd_33[k];

        t_57[k] = f_9 * pfd_15[k]
                  + f_4 * pc_z[k] * pgd_33[k];

        t_58[k] = f_4 * pc_y[k] * pgd_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pc_x, pc_y, pc_z, sgd_36, pfd_17, pfd_18, pfd_36, \
                         pgp0_17, pgp0_18, pgp1_17, pgp1_18, pgd_35, \
                         pgd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_9 * pfd_17[k]
                  + f_2 * pgp0_17[k]
                  - f_3 * pgp1_17[k]
                  + f_4 * pc_z[k] * pgd_35[k];

        t_60[k] = f_0 * sgd_36[k]
                  + f_0 * pfd_36[k]
                  + f_2 * pgp0_18[k]
                  - f_3 * pgp1_18[k]
                  + f_4 * pc_x[k] * pgd_36[k];

        t_61[k] = f_6 * pfd_18[k]
                  + f_4 * pc_y[k] * pgd_36[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, sgd_39, sgd_41, pfd_20, \
                         pfd_39, pfd_41, pgd_36, pgd_38, pgd_39, \
                         pgd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_4 * pc_z[k] * pgd_36[k];

        t_63[k] = f_0 * sgd_39[k]
                  + f_0 * pfd_39[k]
                  + f_4 * pc_x[k] * pgd_39[k];

        t_64[k] = f_6 * pfd_20[k]
                  + f_4 * pc_y[k] * pgd_38[k];

        t_65[k] = f_0 * sgd_41[k]
                  + f_0 * pfd_41[k]
                  + f_4 * pc_x[k] * pgd_41[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_y, pc_z, pfd_21, pfd_23, pgp0_19, pgp0_20, \
                         pgp1_19, pgp1_20, pgd_39, pgd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_6 * pfd_21[k]
                  + f_2 * pgp0_19[k]
                  - f_3 * pgp1_19[k]
                  + f_4 * pc_y[k] * pgd_39[k];

        t_67[k] = f_4 * pc_z[k] * pgd_39[k];

        t_68[k] = f_6 * pfd_23[k]
                  + f_4 * pc_y[k] * pgd_41[k];

        t_69[k] = f_2 * pgp0_20[k]
                  - f_3 * pgp1_20[k]
                  + f_4 * pc_z[k] * pgd_41[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_z, pc_y, pc_z, pff0_30, pff0_33, pfd_18, \
                         pfd_24, pff1_30, pff1_33, pgd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_z[k] * pff0_30[k]
                  - f_5 * pc_z[k] * pff1_30[k];

        t_71[k] = f_9 * pfd_24[k]
                  + f_4 * pc_y[k] * pgd_42[k];

        t_72[k] = f_0 * pfd_18[k]
                  + f_4 * pc_z[k] * pgd_42[k];

        t_73[k] = pb_z[k] * pff0_33[k]
                  - f_5 * pc_z[k] * pff1_33[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pc_x, pc_y, sgd_47, pfd_26, pfd_27, pfd_47, \
                         pgp0_22, pgp1_22, pgd_44, pgd_45, pgd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_9 * pfd_26[k]
                  + f_4 * pc_y[k] * pgd_44[k];

        t_75[k] = f_0 * sgd_47[k]
                  + f_0 * pfd_47[k]
                  + f_4 * pc_x[k] * pgd_47[k];

        t_76[k] = f_9 * pfd_27[k]
                  + f_2 * pgp0_22[k]
                  - f_3 * pgp1_22[k]
                  + f_4 * pc_y[k] * pgd_45[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pb_y, pc_y, pc_z, pff0_50, pfd_21, pfd_23, \
                         pfd_29, pff1_50, pgp0_23, pgp1_23, pgd_45, \
                         pgd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * pfd_21[k]
                  + f_4 * pc_z[k] * pgd_45[k];

        t_78[k] = f_9 * pfd_29[k]
                  + f_4 * pc_y[k] * pgd_47[k];

        t_79[k] = f_0 * pfd_23[k]
                  + f_2 * pgp0_23[k]
                  - f_3 * pgp1_23[k]
                  + f_4 * pc_z[k] * pgd_47[k];

        t_80[k] = pb_y[k] * pff0_50[k]
                  - f_5 * pc_y[k] * pff1_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, sgd_51, pfd_24, pfd_30, \
                         pfd_32, pfd_51, pgd_48, pgd_50, pgd_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_0 * pfd_30[k]
                  + f_4 * pc_y[k] * pgd_48[k];

        t_82[k] = f_9 * pfd_24[k]
                  + f_4 * pc_z[k] * pgd_48[k];

        t_83[k] = f_0 * sgd_51[k]
                  + f_0 * pfd_51[k]
                  + f_4 * pc_x[k] * pgd_51[k];

        t_84[k] = f_0 * pfd_32[k]
                  + f_4 * pc_y[k] * pgd_50[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_y, pc_y, pc_z, pff0_55, pfd_27, pfd_33, \
                         pfd_35, pff1_55, pgp0_25, pgp1_25, pgd_51, \
                         pgd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_y[k] * pff0_55[k]
                  - f_5 * pc_y[k] * pff1_55[k];

        t_86[k] = f_0 * pfd_33[k]
                  + f_2 * pgp0_25[k]
                  - f_3 * pgp1_25[k]
                  + f_4 * pc_y[k] * pgd_51[k];

        t_87[k] = f_9 * pfd_27[k]
                  + f_4 * pc_z[k] * pgd_51[k];

        t_88[k] = f_0 * pfd_35[k]
                  + f_4 * pc_y[k] * pgd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_y, pc_z, sgd_54, pfd_29, pfd_54, pgp0_26, \
                         pgp0_27, pgp1_26, pgp1_27, pgd_53, pgd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_9 * pfd_29[k]
                  + f_2 * pgp0_26[k]
                  - f_3 * pgp1_26[k]
                  + f_4 * pc_z[k] * pgd_53[k];

        t_90[k] = f_0 * sgd_54[k]
                  + f_0 * pfd_54[k]
                  + f_2 * pgp0_27[k]
                  - f_3 * pgp1_27[k]
                  + f_4 * pc_x[k] * pgd_54[k];

        t_91[k] = f_4 * pc_y[k] * pgd_54[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, sgd_57, sgd_59, pfd_30, \
                         pfd_57, pfd_59, pgd_54, pgd_56, pgd_57, \
                         pgd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_6 * pfd_30[k]
                  + f_4 * pc_z[k] * pgd_54[k];

        t_93[k] = f_0 * sgd_57[k]
                  + f_0 * pfd_57[k]
                  + f_4 * pc_x[k] * pgd_57[k];

        t_94[k] = f_4 * pc_y[k] * pgd_56[k];

        t_95[k] = f_0 * sgd_59[k]
                  + f_0 * pfd_59[k]
                  + f_4 * pc_x[k] * pgd_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_y, pc_z, pfd_33, pfd_35, pgp0_28, pgp0_29, \
                         pgp1_28, pgp1_29, pgd_57, pgd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_2 * pgp0_28[k]
                  - f_3 * pgp1_28[k]
                  + f_4 * pc_y[k] * pgd_57[k];

        t_97[k] = f_6 * pfd_33[k]
                  + f_4 * pc_z[k] * pgd_57[k];

        t_98[k] = f_4 * pc_y[k] * pgd_59[k];

        t_99[k] = f_6 * pfd_35[k]
                  + f_2 * pgp0_29[k]
                  - f_3 * pgp1_29[k]
                  + f_4 * pc_z[k] * pgd_59[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pc_x, pc_y, pc_z, sgf0_100, sgd_60, \
                         sgd_63, sgf1_100, pfd_36, pgd_60, pgd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_x[k] * sgf0_100[k]
                   + f_6 * sgd_60[k]
                   - f_5 * pc_x[k] * sgf1_100[k];

        t_101[k] = f_1 * pfd_36[k]
                   + f_4 * pc_y[k] * pgd_60[k];

        t_102[k] = f_4 * pc_z[k] * pgd_60[k];

        t_103[k] = f_0 * sgd_63[k]
                   + f_4 * pc_x[k] * pgd_63[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pc_x, pc_y, pc_z, sgf0_106, sgd_65, \
                         sgf1_106, pfd_38, pgd_62, pgd_63, pgd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * pfd_38[k]
                   + f_4 * pc_y[k] * pgd_62[k];

        t_105[k] = f_0 * sgd_65[k]
                   + f_4 * pc_x[k] * pgd_65[k];

        t_106[k] = pa_x[k] * sgf0_106[k]
                   - f_5 * pc_x[k] * sgf1_106[k];

        t_107[k] = f_4 * pc_z[k] * pgd_63[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_x, pb_z, pc_x, pc_y, pc_z, sgf0_109, \
                         sgf1_109, pff0_60, pfd_41, pff1_60, pgd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_1 * pfd_41[k]
                   + f_4 * pc_y[k] * pgd_65[k];

        t_109[k] = pa_x[k] * sgf0_109[k]
                   - f_5 * pc_x[k] * sgf1_109[k];

        t_110[k] = pb_z[k] * pff0_60[k]
                   - f_5 * pc_z[k] * pff1_60[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, pc_z, sgd_69, pfd_36, pfd_42, \
                         pfd_44, pgd_66, pgd_68, pgd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_6 * pfd_42[k]
                   + f_4 * pc_y[k] * pgd_66[k];

        t_112[k] = f_0 * pfd_36[k]
                   + f_4 * pc_z[k] * pgd_66[k];

        t_113[k] = f_0 * sgd_69[k]
                   + f_4 * pc_x[k] * pgd_69[k];

        t_114[k] = f_6 * pfd_44[k]
                   + f_4 * pc_y[k] * pgd_68[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pa_x, pc_x, pc_y, pc_z, sgf0_116, sgd_71, \
                         sgf1_116, pfd_39, pfd_47, pgd_69, pgd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_0 * sgd_71[k]
                   + f_4 * pc_x[k] * pgd_71[k];

        t_116[k] = pa_x[k] * sgf0_116[k]
                   - f_5 * pc_x[k] * sgf1_116[k];

        t_117[k] = f_0 * pfd_39[k]
                   + f_4 * pc_z[k] * pgd_69[k];

        t_118[k] = f_6 * pfd_47[k]
                   + f_4 * pc_y[k] * pgd_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_x, pc_x, pc_y, pc_z, sgf0_119, \
                         sgf0_120, sgd_72, sgf1_119, sgf1_120, pfd_42, pfd_48, \
                         pgd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pa_x[k] * sgf0_119[k]
                   - f_5 * pc_x[k] * sgf1_119[k];

        t_120[k] = pa_x[k] * sgf0_120[k]
                   + f_6 * sgd_72[k]
                   - f_5 * pc_x[k] * sgf1_120[k];

        t_121[k] = f_9 * pfd_48[k]
                   + f_4 * pc_y[k] * pgd_72[k];

        t_122[k] = f_9 * pfd_42[k]
                   + f_4 * pc_z[k] * pgd_72[k];
    }
}

static auto
compute_prim_pgf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgf0, const size_t sgd,
                                                          const size_t sgf1, const size_t pdf0,
                                                          const size_t pdf1, const size_t pff0,
                                                          const size_t pfd, const size_t pff1,
                                                          const size_t pgp0, const size_t pgp1,
                                                          const size_t pgd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);
    const auto f_9 = 1.0 / q;
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = 1.0 / p;
    const auto f_13 = gamma / (p * q);

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
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgf0_0 = buffer.data(sgf0 + 0);
    const auto *sgf0_1 = buffer.data(sgf0 + 1);
    const auto *sgf0_6 = buffer.data(sgf0 + 6);
    const auto *sgf0_9 = buffer.data(sgf0 + 9);
    const auto *sgf0_20 = buffer.data(sgf0 + 20);
    const auto *sgf0_29 = buffer.data(sgf0 + 29);
    const auto *sgf0_50 = buffer.data(sgf0 + 50);
    const auto *sgf0_51 = buffer.data(sgf0 + 51);
    const auto *sgf0_56 = buffer.data(sgf0 + 56);
    const auto *sgf0_59 = buffer.data(sgf0 + 59);
    const auto *sgf0_90 = buffer.data(sgf0 + 90);
    const auto *sgf0_91 = buffer.data(sgf0 + 91);
    const auto *sgf0_126 = buffer.data(sgf0 + 126);
    const auto *sgf0_129 = buffer.data(sgf0 + 129);
    const auto *sgf0_136 = buffer.data(sgf0 + 136);
    const auto *sgf0_139 = buffer.data(sgf0 + 139);
    const auto *sgf0_140 = buffer.data(sgf0 + 140);
    const auto *sgf0_146 = buffer.data(sgf0 + 146);
    const auto *sgf0_149 = buffer.data(sgf0 + 149);

    const auto *sgd_0 = buffer.data(sgd + 0);
    const auto *sgd_3 = buffer.data(sgd + 3);
    const auto *sgd_5 = buffer.data(sgd + 5);
    const auto *sgd_11 = buffer.data(sgd + 11);
    const auto *sgd_17 = buffer.data(sgd + 17);
    const auto *sgd_23 = buffer.data(sgd + 23);
    const auto *sgd_29 = buffer.data(sgd + 29);
    const auto *sgd_30 = buffer.data(sgd + 30);
    const auto *sgd_33 = buffer.data(sgd + 33);
    const auto *sgd_35 = buffer.data(sgd + 35);
    const auto *sgd_54 = buffer.data(sgd + 54);
    const auto *sgd_75 = buffer.data(sgd + 75);
    const auto *sgd_77 = buffer.data(sgd + 77);
    const auto *sgd_81 = buffer.data(sgd + 81);
    const auto *sgd_83 = buffer.data(sgd + 83);
    const auto *sgd_84 = buffer.data(sgd + 84);
    const auto *sgd_87 = buffer.data(sgd + 87);
    const auto *sgd_89 = buffer.data(sgd + 89);

    const auto *sgf1_0 = buffer.data(sgf1 + 0);
    const auto *sgf1_1 = buffer.data(sgf1 + 1);
    const auto *sgf1_6 = buffer.data(sgf1 + 6);
    const auto *sgf1_9 = buffer.data(sgf1 + 9);
    const auto *sgf1_20 = buffer.data(sgf1 + 20);
    const auto *sgf1_29 = buffer.data(sgf1 + 29);
    const auto *sgf1_50 = buffer.data(sgf1 + 50);
    const auto *sgf1_51 = buffer.data(sgf1 + 51);
    const auto *sgf1_56 = buffer.data(sgf1 + 56);
    const auto *sgf1_59 = buffer.data(sgf1 + 59);
    const auto *sgf1_90 = buffer.data(sgf1 + 90);
    const auto *sgf1_91 = buffer.data(sgf1 + 91);
    const auto *sgf1_126 = buffer.data(sgf1 + 126);
    const auto *sgf1_129 = buffer.data(sgf1 + 129);
    const auto *sgf1_136 = buffer.data(sgf1 + 136);
    const auto *sgf1_139 = buffer.data(sgf1 + 139);
    const auto *sgf1_140 = buffer.data(sgf1 + 140);
    const auto *sgf1_146 = buffer.data(sgf1 + 146);
    const auto *sgf1_149 = buffer.data(sgf1 + 149);

    const auto *pdf0_76 = buffer.data(pdf0 + 76);
    const auto *pdf0_96 = buffer.data(pdf0 + 96);

    const auto *pdf1_76 = buffer.data(pdf1 + 76);
    const auto *pdf1_96 = buffer.data(pdf1 + 96);

    const auto *pff0_90 = buffer.data(pff0 + 90);
    const auto *pff0_101 = buffer.data(pff0 + 101);
    const auto *pff0_106 = buffer.data(pff0 + 106);
    const auto *pff0_111 = buffer.data(pff0 + 111);
    const auto *pff0_116 = buffer.data(pff0 + 116);
    const auto *pff0_130 = buffer.data(pff0 + 130);
    const auto *pff0_131 = buffer.data(pff0 + 131);
    const auto *pff0_136 = buffer.data(pff0 + 136);
    const auto *pff0_160 = buffer.data(pff0 + 160);
    const auto *pff0_161 = buffer.data(pff0 + 161);
    const auto *pff0_166 = buffer.data(pff0 + 166);
    const auto *pff0_168 = buffer.data(pff0 + 168);
    const auto *pff0_169 = buffer.data(pff0 + 169);
    const auto *pff0_176 = buffer.data(pff0 + 176);
    const auto *pff0_178 = buffer.data(pff0 + 178);
    const auto *pff0_179 = buffer.data(pff0 + 179);
    const auto *pff0_181 = buffer.data(pff0 + 181);
    const auto *pff0_186 = buffer.data(pff0 + 186);
    const auto *pff0_188 = buffer.data(pff0 + 188);
    const auto *pff0_189 = buffer.data(pff0 + 189);
    const auto *pff0_196 = buffer.data(pff0 + 196);

    const auto *pfd_45 = buffer.data(pfd + 45);
    const auto *pfd_48 = buffer.data(pfd + 48);
    const auto *pfd_50 = buffer.data(pfd + 50);
    const auto *pfd_51 = buffer.data(pfd + 51);
    const auto *pfd_53 = buffer.data(pfd + 53);
    const auto *pfd_54 = buffer.data(pfd + 54);
    const auto *pfd_56 = buffer.data(pfd + 56);
    const auto *pfd_57 = buffer.data(pfd + 57);
    const auto *pfd_59 = buffer.data(pfd + 59);
    const auto *pfd_60 = buffer.data(pfd + 60);
    const auto *pfd_63 = buffer.data(pfd + 63);
    const auto *pfd_64 = buffer.data(pfd + 64);
    const auto *pfd_65 = buffer.data(pfd + 65);
    const auto *pfd_66 = buffer.data(pfd + 66);
    const auto *pfd_67 = buffer.data(pfd + 67);
    const auto *pfd_69 = buffer.data(pfd + 69);
    const auto *pfd_70 = buffer.data(pfd + 70);
    const auto *pfd_71 = buffer.data(pfd + 71);
    const auto *pfd_72 = buffer.data(pfd + 72);
    const auto *pfd_75 = buffer.data(pfd + 75);
    const auto *pfd_76 = buffer.data(pfd + 76);
    const auto *pfd_77 = buffer.data(pfd + 77);
    const auto *pfd_78 = buffer.data(pfd + 78);
    const auto *pfd_79 = buffer.data(pfd + 79);
    const auto *pfd_81 = buffer.data(pfd + 81);
    const auto *pfd_82 = buffer.data(pfd + 82);
    const auto *pfd_83 = buffer.data(pfd + 83);
    const auto *pfd_84 = buffer.data(pfd + 84);
    const auto *pfd_87 = buffer.data(pfd + 87);
    const auto *pfd_88 = buffer.data(pfd + 88);
    const auto *pfd_89 = buffer.data(pfd + 89);
    const auto *pfd_90 = buffer.data(pfd + 90);
    const auto *pfd_93 = buffer.data(pfd + 93);
    const auto *pfd_94 = buffer.data(pfd + 94);
    const auto *pfd_95 = buffer.data(pfd + 95);
    const auto *pfd_96 = buffer.data(pfd + 96);
    const auto *pfd_97 = buffer.data(pfd + 97);
    const auto *pfd_99 = buffer.data(pfd + 99);
    const auto *pfd_100 = buffer.data(pfd + 100);
    const auto *pfd_101 = buffer.data(pfd + 101);
    const auto *pfd_105 = buffer.data(pfd + 105);
    const auto *pfd_106 = buffer.data(pfd + 106);
    const auto *pfd_107 = buffer.data(pfd + 107);
    const auto *pfd_108 = buffer.data(pfd + 108);
    const auto *pfd_109 = buffer.data(pfd + 109);
    const auto *pfd_111 = buffer.data(pfd + 111);
    const auto *pfd_112 = buffer.data(pfd + 112);
    const auto *pfd_113 = buffer.data(pfd + 113);
    const auto *pfd_117 = buffer.data(pfd + 117);
    const auto *pfd_118 = buffer.data(pfd + 118);
    const auto *pfd_119 = buffer.data(pfd + 119);

    const auto *pff1_90 = buffer.data(pff1 + 90);
    const auto *pff1_101 = buffer.data(pff1 + 101);
    const auto *pff1_106 = buffer.data(pff1 + 106);
    const auto *pff1_111 = buffer.data(pff1 + 111);
    const auto *pff1_116 = buffer.data(pff1 + 116);
    const auto *pff1_130 = buffer.data(pff1 + 130);
    const auto *pff1_131 = buffer.data(pff1 + 131);
    const auto *pff1_136 = buffer.data(pff1 + 136);
    const auto *pff1_160 = buffer.data(pff1 + 160);
    const auto *pff1_161 = buffer.data(pff1 + 161);
    const auto *pff1_166 = buffer.data(pff1 + 166);
    const auto *pff1_168 = buffer.data(pff1 + 168);
    const auto *pff1_169 = buffer.data(pff1 + 169);
    const auto *pff1_176 = buffer.data(pff1 + 176);
    const auto *pff1_178 = buffer.data(pff1 + 178);
    const auto *pff1_179 = buffer.data(pff1 + 179);
    const auto *pff1_181 = buffer.data(pff1 + 181);
    const auto *pff1_186 = buffer.data(pff1 + 186);
    const auto *pff1_188 = buffer.data(pff1 + 188);
    const auto *pff1_189 = buffer.data(pff1 + 189);
    const auto *pff1_196 = buffer.data(pff1 + 196);

    const auto *pgp0_48 = buffer.data(pgp0 + 48);
    const auto *pgp0_49 = buffer.data(pgp0 + 49);
    const auto *pgp0_50 = buffer.data(pgp0 + 50);
    const auto *pgp0_54 = buffer.data(pgp0 + 54);
    const auto *pgp0_55 = buffer.data(pgp0 + 55);
    const auto *pgp0_56 = buffer.data(pgp0 + 56);
    const auto *pgp0_57 = buffer.data(pgp0 + 57);
    const auto *pgp0_59 = buffer.data(pgp0 + 59);
    const auto *pgp0_69 = buffer.data(pgp0 + 69);

    const auto *pgp1_48 = buffer.data(pgp1 + 48);
    const auto *pgp1_49 = buffer.data(pgp1 + 49);
    const auto *pgp1_50 = buffer.data(pgp1 + 50);
    const auto *pgp1_54 = buffer.data(pgp1 + 54);
    const auto *pgp1_55 = buffer.data(pgp1 + 55);
    const auto *pgp1_56 = buffer.data(pgp1 + 56);
    const auto *pgp1_57 = buffer.data(pgp1 + 57);
    const auto *pgp1_59 = buffer.data(pgp1 + 59);
    const auto *pgp1_69 = buffer.data(pgp1 + 69);

    const auto *pgd_74 = buffer.data(pgd + 74);
    const auto *pgd_75 = buffer.data(pgd + 75);
    const auto *pgd_77 = buffer.data(pgd + 77);
    const auto *pgd_78 = buffer.data(pgd + 78);
    const auto *pgd_80 = buffer.data(pgd + 80);
    const auto *pgd_81 = buffer.data(pgd + 81);
    const auto *pgd_83 = buffer.data(pgd + 83);
    const auto *pgd_84 = buffer.data(pgd + 84);
    const auto *pgd_86 = buffer.data(pgd + 86);
    const auto *pgd_87 = buffer.data(pgd + 87);
    const auto *pgd_89 = buffer.data(pgd + 89);
    const auto *pgd_90 = buffer.data(pgd + 90);
    const auto *pgd_93 = buffer.data(pgd + 93);
    const auto *pgd_94 = buffer.data(pgd + 94);
    const auto *pgd_95 = buffer.data(pgd + 95);
    const auto *pgd_96 = buffer.data(pgd + 96);
    const auto *pgd_97 = buffer.data(pgd + 97);
    const auto *pgd_99 = buffer.data(pgd + 99);
    const auto *pgd_100 = buffer.data(pgd + 100);
    const auto *pgd_101 = buffer.data(pgd + 101);
    const auto *pgd_102 = buffer.data(pgd + 102);
    const auto *pgd_105 = buffer.data(pgd + 105);
    const auto *pgd_106 = buffer.data(pgd + 106);
    const auto *pgd_107 = buffer.data(pgd + 107);
    const auto *pgd_108 = buffer.data(pgd + 108);
    const auto *pgd_109 = buffer.data(pgd + 109);
    const auto *pgd_111 = buffer.data(pgd + 111);
    const auto *pgd_112 = buffer.data(pgd + 112);
    const auto *pgd_113 = buffer.data(pgd + 113);
    const auto *pgd_114 = buffer.data(pgd + 114);
    const auto *pgd_117 = buffer.data(pgd + 117);
    const auto *pgd_118 = buffer.data(pgd + 118);
    const auto *pgd_119 = buffer.data(pgd + 119);
    const auto *pgd_120 = buffer.data(pgd + 120);
    const auto *pgd_123 = buffer.data(pgd + 123);
    const auto *pgd_124 = buffer.data(pgd + 124);
    const auto *pgd_125 = buffer.data(pgd + 125);
    const auto *pgd_126 = buffer.data(pgd + 126);
    const auto *pgd_129 = buffer.data(pgd + 129);
    const auto *pgd_130 = buffer.data(pgd + 130);
    const auto *pgd_131 = buffer.data(pgd + 131);
    const auto *pgd_132 = buffer.data(pgd + 132);
    const auto *pgd_135 = buffer.data(pgd + 135);
    const auto *pgd_136 = buffer.data(pgd + 136);
    const auto *pgd_137 = buffer.data(pgd + 137);
    const auto *pgd_138 = buffer.data(pgd + 138);
    const auto *pgd_141 = buffer.data(pgd + 141);
    const auto *pgd_142 = buffer.data(pgd + 142);
    const auto *pgd_143 = buffer.data(pgd + 143);
    const auto *pgd_144 = buffer.data(pgd + 144);
    const auto *pgd_147 = buffer.data(pgd + 147);
    const auto *pgd_148 = buffer.data(pgd + 148);
    const auto *pgd_149 = buffer.data(pgd + 149);

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_x, pc_x, pc_y, sgf0_126, sgd_75, \
                         sgd_77, sgf1_126, pfd_50, pgd_74, pgd_75, \
                         pgd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * sgd_75[k]
                   + f_4 * pc_x[k] * pgd_75[k];

        t_124[k] = f_9 * pfd_50[k]
                   + f_4 * pc_y[k] * pgd_74[k];

        t_125[k] = f_0 * sgd_77[k]
                   + f_4 * pc_x[k] * pgd_77[k];

        t_126[k] = pa_x[k] * sgf0_126[k]
                   - f_5 * pc_x[k] * sgf1_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pa_x, pc_x, pc_y, pc_z, sgf0_129, sgf1_129, \
                         pfd_45, pfd_53, pgd_75, pgd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_9 * pfd_45[k]
                   + f_4 * pc_z[k] * pgd_75[k];

        t_128[k] = f_9 * pfd_53[k]
                   + f_4 * pc_y[k] * pgd_77[k];

        t_129[k] = pa_x[k] * sgf0_129[k]
                   - f_5 * pc_x[k] * sgf1_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pb_y, pc_x, pc_y, pc_z, sgd_81, pff0_90, \
                         pfd_48, pfd_54, pff1_90, pgd_78, pgd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_y[k] * pff0_90[k]
                   - f_5 * pc_y[k] * pff1_90[k];

        t_131[k] = f_0 * pfd_54[k]
                   + f_4 * pc_y[k] * pgd_78[k];

        t_132[k] = f_6 * pfd_48[k]
                   + f_4 * pc_z[k] * pgd_78[k];

        t_133[k] = f_0 * sgd_81[k]
                   + f_4 * pc_x[k] * pgd_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pc_x, pc_y, pc_z, sgf0_136, sgd_83, \
                         sgf1_136, pfd_51, pfd_56, pgd_80, pgd_81, \
                         pgd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_0 * pfd_56[k]
                   + f_4 * pc_y[k] * pgd_80[k];

        t_135[k] = f_0 * sgd_83[k]
                   + f_4 * pc_x[k] * pgd_83[k];

        t_136[k] = pa_x[k] * sgf0_136[k]
                   - f_5 * pc_x[k] * sgf1_136[k];

        t_137[k] = f_6 * pfd_51[k]
                   + f_4 * pc_z[k] * pgd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_x, pc_x, pc_y, sgf0_139, sgf0_140, \
                         sgd_84, sgf1_139, sgf1_140, pfd_59, pgd_83, \
                         pgd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * pfd_59[k]
                   + f_4 * pc_y[k] * pgd_83[k];

        t_139[k] = pa_x[k] * sgf0_139[k]
                   - f_5 * pc_x[k] * sgf1_139[k];

        t_140[k] = pa_x[k] * sgf0_140[k]
                   + f_6 * sgd_84[k]
                   - f_5 * pc_x[k] * sgf1_140[k];

        t_141[k] = f_4 * pc_y[k] * pgd_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_y, pc_z, sgd_87, sgd_89, pfd_54, \
                         pgd_84, pgd_86, pgd_87, pgd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_1 * pfd_54[k]
                   + f_4 * pc_z[k] * pgd_84[k];

        t_143[k] = f_0 * sgd_87[k]
                   + f_4 * pc_x[k] * pgd_87[k];

        t_144[k] = f_4 * pc_y[k] * pgd_86[k];

        t_145[k] = f_0 * sgd_89[k]
                   + f_4 * pc_x[k] * pgd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pc_x, pc_y, pc_z, sgf0_146, \
                         sgf0_149, sgf1_146, sgf1_149, pfd_57, pgd_87, \
                         pgd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_x[k] * sgf0_146[k]
                   - f_5 * pc_x[k] * sgf1_146[k];

        t_147[k] = f_1 * pfd_57[k]
                   + f_4 * pc_z[k] * pgd_87[k];

        t_148[k] = f_4 * pc_y[k] * pgd_89[k];

        t_149[k] = pa_x[k] * sgf0_149[k]
                   - f_5 * pc_x[k] * sgf1_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_y, pc_x, pc_y, pc_z, sgf0_0, sgf0_1, \
                         sgd_0, sgf1_0, sgf1_1, pfd_63, pgd_90, \
                         pgd_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_y[k] * sgf0_0[k]
                   - f_5 * pc_y[k] * sgf1_0[k];

        t_151[k] = pa_y[k] * sgf0_1[k]
                   + f_0 * sgd_0[k]
                   - f_5 * pc_y[k] * sgf1_1[k];

        t_152[k] = f_4 * pc_z[k] * pgd_90[k];

        t_153[k] = f_1 * pfd_63[k]
                   + f_4 * pc_x[k] * pgd_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_y, pc_x, pc_y, pc_z, sgf0_6, sgd_3, \
                         sgf1_6, pfd_64, pfd_65, pgd_93, pgd_94, \
                         pgd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_1 * pfd_64[k]
                   + f_4 * pc_x[k] * pgd_94[k];

        t_155[k] = f_1 * pfd_65[k]
                   + f_4 * pc_x[k] * pgd_95[k];

        t_156[k] = pa_y[k] * sgf0_6[k]
                   + f_6 * sgd_3[k]
                   - f_5 * pc_y[k] * sgf1_6[k];

        t_157[k] = f_4 * pc_z[k] * pgd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pc_x, pc_y, sgf0_9, sgd_5, sgf1_9, pfd_66, \
                         pgp0_48, pgp1_48, pgd_95, pgd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_0 * sgd_5[k]
                   + f_4 * pc_y[k] * pgd_95[k];

        t_159[k] = pa_y[k] * sgf0_9[k]
                   - f_5 * pc_y[k] * sgf1_9[k];

        t_160[k] = f_6 * pfd_66[k]
                   + f_2 * pgp0_48[k]
                   - f_3 * pgp1_48[k]
                   + f_4 * pc_x[k] * pgd_96[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pc_x, pc_z, pfd_67, pfd_69, pfd_70, \
                         pgp0_49, pgp1_49, pgd_96, pgd_97, pgd_99, \
                         pgd_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_6 * pfd_67[k]
                   + f_10 * pgp0_49[k]
                   - f_11 * pgp1_49[k]
                   + f_4 * pc_x[k] * pgd_97[k];

        t_162[k] = f_4 * pc_z[k] * pgd_96[k];

        t_163[k] = f_6 * pfd_69[k]
                   + f_4 * pc_x[k] * pgd_99[k];

        t_164[k] = f_6 * pfd_70[k]
                   + f_4 * pc_x[k] * pgd_100[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pb_x, pc_x, pc_z, pdf0_76, pdf1_76, pff0_116, \
                         pfd_71, pff1_116, pgd_99, pgd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_6 * pfd_71[k]
                   + f_4 * pc_x[k] * pgd_101[k];

        t_166[k] = f_12 * pdf0_76[k]
                   - f_13 * pdf1_76[k]
                   + pb_x[k] * pff0_116[k]
                   - f_5 * pc_x[k] * pff1_116[k];

        t_167[k] = f_4 * pc_z[k] * pgd_99[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_y, pc_y, pc_z, sgf0_20, sgd_11, sgf1_20, \
                         pfd_65, pgp0_50, pgp1_50, pgd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_0 * sgd_11[k]
                   + f_0 * pfd_65[k]
                   + f_4 * pc_y[k] * pgd_101[k];

        t_169[k] = f_2 * pgp0_50[k]
                   - f_3 * pgp1_50[k]
                   + f_4 * pc_z[k] * pgd_101[k];

        t_170[k] = pa_y[k] * sgf0_20[k]
                   - f_5 * pc_y[k] * sgf1_20[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pb_z, pc_x, pc_z, pff0_101, pfd_60, \
                         pfd_75, pfd_76, pff1_101, pgd_102, pgd_105, \
                         pgd_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pb_z[k] * pff0_101[k]
                   - f_5 * pc_z[k] * pff1_101[k];

        t_172[k] = f_0 * pfd_60[k]
                   + f_4 * pc_z[k] * pgd_102[k];

        t_173[k] = f_6 * pfd_75[k]
                   + f_4 * pc_x[k] * pgd_105[k];

        t_174[k] = f_6 * pfd_76[k]
                   + f_4 * pc_x[k] * pgd_106[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pb_z, pc_x, pc_y, pc_z, sgd_17, pff0_106, \
                         pfd_63, pfd_77, pff1_106, pgd_105, pgd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_6 * pfd_77[k]
                   + f_4 * pc_x[k] * pgd_107[k];

        t_176[k] = pb_z[k] * pff0_106[k]
                   - f_5 * pc_z[k] * pff1_106[k];

        t_177[k] = f_0 * pfd_63[k]
                   + f_4 * pc_z[k] * pgd_105[k];

        t_178[k] = f_0 * sgd_17[k]
                   + f_4 * pc_y[k] * pgd_107[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pa_y, pc_x, pc_y, sgf0_29, sgf1_29, pfd_78, \
                         pfd_79, pgp0_54, pgp0_55, pgp1_54, pgp1_55, pgd_108, \
                         pgd_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pa_y[k] * sgf0_29[k]
                   - f_5 * pc_y[k] * sgf1_29[k];

        t_180[k] = f_9 * pfd_78[k]
                   + f_2 * pgp0_54[k]
                   - f_3 * pgp1_54[k]
                   + f_4 * pc_x[k] * pgd_108[k];

        t_181[k] = f_9 * pfd_79[k]
                   + f_10 * pgp0_55[k]
                   - f_11 * pgp1_55[k]
                   + f_4 * pc_x[k] * pgd_109[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pc_x, pc_z, pfd_81, pfd_82, pfd_83, \
                         pgd_108, pgd_111, pgd_112, pgd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_4 * pc_z[k] * pgd_108[k];

        t_183[k] = f_9 * pfd_81[k]
                   + f_4 * pc_x[k] * pgd_111[k];

        t_184[k] = f_9 * pfd_82[k]
                   + f_4 * pc_x[k] * pgd_112[k];

        t_185[k] = f_9 * pfd_83[k]
                   + f_4 * pc_x[k] * pgd_113[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_x, pc_x, pc_y, pc_z, sgd_23, pdf0_96, \
                         pdf1_96, pff0_136, pfd_71, pff1_136, pgd_111, \
                         pgd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_7 * pdf0_96[k]
                   - f_8 * pdf1_96[k]
                   + pb_x[k] * pff0_136[k]
                   - f_5 * pc_x[k] * pff1_136[k];

        t_187[k] = f_4 * pc_z[k] * pgd_111[k];

        t_188[k] = f_0 * sgd_23[k]
                   + f_9 * pfd_71[k]
                   + f_4 * pc_y[k] * pgd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_z, pc_x, pc_z, pff0_111, pfd_84, pff1_111, \
                         pgp0_56, pgp0_57, pgp1_56, pgp1_57, pgd_113, \
                         pgd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_2 * pgp0_56[k]
                   - f_3 * pgp1_56[k]
                   + f_4 * pc_z[k] * pgd_113[k];

        t_190[k] = f_9 * pfd_84[k]
                   + f_2 * pgp0_57[k]
                   - f_3 * pgp1_57[k]
                   + f_4 * pc_x[k] * pgd_114[k];

        t_191[k] = pb_z[k] * pff0_111[k]
                   - f_5 * pc_z[k] * pff1_111[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_x, pc_z, pfd_66, pfd_87, pfd_88, \
                         pfd_89, pgd_114, pgd_117, pgd_118, pgd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_0 * pfd_66[k]
                   + f_4 * pc_z[k] * pgd_114[k];

        t_193[k] = f_9 * pfd_87[k]
                   + f_4 * pc_x[k] * pgd_117[k];

        t_194[k] = f_9 * pfd_88[k]
                   + f_4 * pc_x[k] * pgd_118[k];

        t_195[k] = f_9 * pfd_89[k]
                   + f_4 * pc_x[k] * pgd_119[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pb_z, pc_y, pc_z, sgd_29, pff0_116, pfd_69, \
                         pfd_77, pff1_116, pgd_117, pgd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_z[k] * pff0_116[k]
                   - f_5 * pc_z[k] * pff1_116[k];

        t_197[k] = f_0 * pfd_69[k]
                   + f_4 * pc_z[k] * pgd_117[k];

        t_198[k] = f_0 * sgd_29[k]
                   + f_0 * pfd_77[k]
                   + f_4 * pc_y[k] * pgd_119[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, pa_y, pc_y, pc_z, sgf0_50, sgf0_51, sgd_30, \
                         sgf1_50, sgf1_51, pfd_71, pgp0_59, pgp1_59, \
                         pgd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * pfd_71[k]
                   + f_2 * pgp0_59[k]
                   - f_3 * pgp1_59[k]
                   + f_4 * pc_z[k] * pgd_119[k];

        t_200[k] = pa_y[k] * sgf0_50[k]
                   - f_5 * pc_y[k] * sgf1_50[k];

        t_201[k] = pa_y[k] * sgf0_51[k]
                   + f_0 * sgd_30[k]
                   - f_5 * pc_y[k] * sgf1_51[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pc_x, pc_z, pfd_72, pfd_93, pfd_94, \
                         pfd_95, pgd_120, pgd_123, pgd_124, pgd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_9 * pfd_72[k]
                   + f_4 * pc_z[k] * pgd_120[k];

        t_203[k] = f_9 * pfd_93[k]
                   + f_4 * pc_x[k] * pgd_123[k];

        t_204[k] = f_9 * pfd_94[k]
                   + f_4 * pc_x[k] * pgd_124[k];

        t_205[k] = f_9 * pfd_95[k]
                   + f_4 * pc_x[k] * pgd_125[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_y, pc_y, pc_z, sgf0_56, sgf0_59, \
                         sgd_33, sgd_35, sgf1_56, sgf1_59, pfd_75, pgd_123, \
                         pgd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pa_y[k] * sgf0_56[k]
                   + f_6 * sgd_33[k]
                   - f_5 * pc_y[k] * sgf1_56[k];

        t_207[k] = f_9 * pfd_75[k]
                   + f_4 * pc_z[k] * pgd_123[k];

        t_208[k] = f_0 * sgd_35[k]
                   + f_4 * pc_y[k] * pgd_125[k];

        t_209[k] = pa_y[k] * sgf0_59[k]
                   - f_5 * pc_y[k] * sgf1_59[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pb_x, pc_x, pc_z, pff0_160, pff0_161, \
                         pfd_96, pfd_97, pfd_99, pff1_160, pff1_161, pgd_126, \
                         pgd_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = pb_x[k] * pff0_160[k]
                   + f_6 * pfd_96[k]
                   - f_5 * pc_x[k] * pff1_160[k];

        t_211[k] = pb_x[k] * pff0_161[k]
                   + f_9 * pfd_97[k]
                   - f_5 * pc_x[k] * pff1_161[k];

        t_212[k] = f_4 * pc_z[k] * pgd_126[k];

        t_213[k] = f_0 * pfd_99[k]
                   + f_4 * pc_x[k] * pgd_129[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pb_x, pc_x, pc_z, pff0_166, pfd_100, \
                         pfd_101, pff1_166, pgd_129, pgd_130, pgd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_0 * pfd_100[k]
                   + f_4 * pc_x[k] * pgd_130[k];

        t_215[k] = f_0 * pfd_101[k]
                   + f_4 * pc_x[k] * pgd_131[k];

        t_216[k] = pb_x[k] * pff0_166[k]
                   - f_5 * pc_x[k] * pff1_166[k];

        t_217[k] = f_4 * pc_z[k] * pgd_129[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pb_x, pb_z, pc_x, pc_z, pff0_130, \
                         pff0_131, pff0_168, pff0_169, pff1_130, pff1_131, pff1_168, \
                         pff1_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = pb_x[k] * pff0_168[k]
                   - f_5 * pc_x[k] * pff1_168[k];

        t_219[k] = pb_x[k] * pff0_169[k]
                   - f_5 * pc_x[k] * pff1_169[k];

        t_220[k] = pb_z[k] * pff0_130[k]
                   - f_5 * pc_z[k] * pff1_130[k];

        t_221[k] = pb_z[k] * pff0_131[k]
                   - f_5 * pc_z[k] * pff1_131[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_z, pfd_78, pfd_105, pfd_106, \
                         pfd_107, pgd_132, pgd_135, pgd_136, pgd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * pfd_78[k]
                   + f_4 * pc_z[k] * pgd_132[k];

        t_223[k] = f_0 * pfd_105[k]
                   + f_4 * pc_x[k] * pgd_135[k];

        t_224[k] = f_0 * pfd_106[k]
                   + f_4 * pc_x[k] * pgd_136[k];

        t_225[k] = f_0 * pfd_107[k]
                   + f_4 * pc_x[k] * pgd_137[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pb_x, pc_x, pc_z, pff0_176, pff0_178, \
                         pff0_179, pfd_81, pff1_176, pff1_178, pff1_179, \
                         pgd_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pb_x[k] * pff0_176[k]
                   - f_5 * pc_x[k] * pff1_176[k];

        t_227[k] = f_0 * pfd_81[k]
                   + f_4 * pc_z[k] * pgd_135[k];

        t_228[k] = pb_x[k] * pff0_178[k]
                   - f_5 * pc_x[k] * pff1_178[k];

        t_229[k] = pb_x[k] * pff0_179[k]
                   - f_5 * pc_x[k] * pff1_179[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pb_x, pc_x, pc_z, pff0_181, pfd_84, pfd_108, \
                         pfd_109, pff1_181, pgp0_69, pgp1_69, pgd_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_0 * pfd_108[k]
                   + f_2 * pgp0_69[k]
                   - f_3 * pgp1_69[k]
                   + f_4 * pc_x[k] * pgd_138[k];

        t_231[k] = pb_x[k] * pff0_181[k]
                   + f_9 * pfd_109[k]
                   - f_5 * pc_x[k] * pff1_181[k];

        t_232[k] = f_9 * pfd_84[k]
                   + f_4 * pc_z[k] * pgd_138[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_x, pc_x, pff0_186, pfd_111, pfd_112, \
                         pfd_113, pff1_186, pgd_141, pgd_142, pgd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_0 * pfd_111[k]
                   + f_4 * pc_x[k] * pgd_141[k];

        t_234[k] = f_0 * pfd_112[k]
                   + f_4 * pc_x[k] * pgd_142[k];

        t_235[k] = f_0 * pfd_113[k]
                   + f_4 * pc_x[k] * pgd_143[k];

        t_236[k] = pb_x[k] * pff0_186[k]
                   - f_5 * pc_x[k] * pff1_186[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pb_x, pc_x, pc_z, pff0_188, pff0_189, pfd_87, \
                         pff1_188, pff1_189, pgd_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_9 * pfd_87[k]
                   + f_4 * pc_z[k] * pgd_141[k];

        t_238[k] = pb_x[k] * pff0_188[k]
                   - f_5 * pc_x[k] * pff1_188[k];

        t_239[k] = pb_x[k] * pff0_189[k]
                   - f_5 * pc_x[k] * pff1_189[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pa_y, pc_y, pc_z, sgf0_90, sgf0_91, sgd_54, \
                         sgf1_90, sgf1_91, pfd_90, pgd_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pa_y[k] * sgf0_90[k]
                   - f_5 * pc_y[k] * sgf1_90[k];

        t_241[k] = pa_y[k] * sgf0_91[k]
                   + f_0 * sgd_54[k]
                   - f_5 * pc_y[k] * sgf1_91[k];

        t_242[k] = f_6 * pfd_90[k]
                   + f_4 * pc_z[k] * pgd_144[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_x, pc_x, pff0_196, pfd_117, pfd_118, \
                         pfd_119, pff1_196, pgd_147, pgd_148, pgd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_0 * pfd_117[k]
                   + f_4 * pc_x[k] * pgd_147[k];

        t_244[k] = f_0 * pfd_118[k]
                   + f_4 * pc_x[k] * pgd_148[k];

        t_245[k] = f_0 * pfd_119[k]
                   + f_4 * pc_x[k] * pgd_149[k];

        t_246[k] = pb_x[k] * pff0_196[k]
                   - f_5 * pc_x[k] * pff1_196[k];
    }
}

static auto
compute_prim_pgf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgf0, const size_t sgd,
                                                          const size_t sgf1, const size_t pdf0,
                                                          const size_t pdf1, const size_t pff0,
                                                          const size_t pfd, const size_t pff1,
                                                          const size_t pgp0, const size_t pgp1,
                                                          const size_t pgd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);
    const auto f_9 = 1.0 / q;
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = 1.0 / p;
    const auto f_13 = gamma / (p * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgf0_0 = buffer.data(sgf0 + 0);
    const auto *sgf0_2 = buffer.data(sgf0 + 2);
    const auto *sgf0_6 = buffer.data(sgf0 + 6);
    const auto *sgf0_9 = buffer.data(sgf0 + 9);
    const auto *sgf0_10 = buffer.data(sgf0 + 10);
    const auto *sgf0_16 = buffer.data(sgf0 + 16);
    const auto *sgf0_17 = buffer.data(sgf0 + 17);
    const auto *sgf0_30 = buffer.data(sgf0 + 30);
    const auto *sgf0_32 = buffer.data(sgf0 + 32);
    const auto *sgf0_36 = buffer.data(sgf0 + 36);
    const auto *sgf0_37 = buffer.data(sgf0 + 37);
    const auto *sgf0_39 = buffer.data(sgf0 + 39);
    const auto *sgf0_60 = buffer.data(sgf0 + 60);
    const auto *sgf0_62 = buffer.data(sgf0 + 62);
    const auto *sgf0_66 = buffer.data(sgf0 + 66);
    const auto *sgf0_99 = buffer.data(sgf0 + 99);
    const auto *sgf0_140 = buffer.data(sgf0 + 140);
    const auto *sgf0_146 = buffer.data(sgf0 + 146);
    const auto *sgf0_149 = buffer.data(sgf0 + 149);

    const auto *sgd_0 = buffer.data(sgd + 0);
    const auto *sgd_5 = buffer.data(sgd + 5);
    const auto *sgd_9 = buffer.data(sgd + 9);
    const auto *sgd_18 = buffer.data(sgd + 18);
    const auto *sgd_21 = buffer.data(sgd + 21);
    const auto *sgd_23 = buffer.data(sgd + 23);
    const auto *sgd_36 = buffer.data(sgd + 36);
    const auto *sgd_59 = buffer.data(sgd + 59);
    const auto *sgd_63 = buffer.data(sgd + 63);
    const auto *sgd_65 = buffer.data(sgd + 65);
    const auto *sgd_71 = buffer.data(sgd + 71);
    const auto *sgd_77 = buffer.data(sgd + 77);
    const auto *sgd_81 = buffer.data(sgd + 81);
    const auto *sgd_83 = buffer.data(sgd + 83);
    const auto *sgd_87 = buffer.data(sgd + 87);
    const auto *sgd_89 = buffer.data(sgd + 89);

    const auto *sgf1_0 = buffer.data(sgf1 + 0);
    const auto *sgf1_2 = buffer.data(sgf1 + 2);
    const auto *sgf1_6 = buffer.data(sgf1 + 6);
    const auto *sgf1_9 = buffer.data(sgf1 + 9);
    const auto *sgf1_10 = buffer.data(sgf1 + 10);
    const auto *sgf1_16 = buffer.data(sgf1 + 16);
    const auto *sgf1_17 = buffer.data(sgf1 + 17);
    const auto *sgf1_30 = buffer.data(sgf1 + 30);
    const auto *sgf1_32 = buffer.data(sgf1 + 32);
    const auto *sgf1_36 = buffer.data(sgf1 + 36);
    const auto *sgf1_37 = buffer.data(sgf1 + 37);
    const auto *sgf1_39 = buffer.data(sgf1 + 39);
    const auto *sgf1_60 = buffer.data(sgf1 + 60);
    const auto *sgf1_62 = buffer.data(sgf1 + 62);
    const auto *sgf1_66 = buffer.data(sgf1 + 66);
    const auto *sgf1_99 = buffer.data(sgf1 + 99);
    const auto *sgf1_140 = buffer.data(sgf1 + 140);
    const auto *sgf1_146 = buffer.data(sgf1 + 146);
    const auto *sgf1_149 = buffer.data(sgf1 + 149);

    const auto *pdf0_96 = buffer.data(pdf0 + 96);
    const auto *pdf0_149 = buffer.data(pdf0 + 149);
    const auto *pdf0_179 = buffer.data(pdf0 + 179);

    const auto *pdf1_96 = buffer.data(pdf1 + 96);
    const auto *pdf1_149 = buffer.data(pdf1 + 149);
    const auto *pdf1_179 = buffer.data(pdf1 + 179);

    const auto *pff0_160 = buffer.data(pff0 + 160);
    const auto *pff0_161 = buffer.data(pff0 + 161);
    const auto *pff0_166 = buffer.data(pff0 + 166);
    const auto *pff0_176 = buffer.data(pff0 + 176);
    const auto *pff0_202 = buffer.data(pff0 + 202);
    const auto *pff0_209 = buffer.data(pff0 + 209);
    const auto *pff0_220 = buffer.data(pff0 + 220);
    const auto *pff0_222 = buffer.data(pff0 + 222);
    const auto *pff0_229 = buffer.data(pff0 + 229);
    const auto *pff0_259 = buffer.data(pff0 + 259);
    const auto *pff0_267 = buffer.data(pff0 + 267);

    const auto *pfd_93 = buffer.data(pfd + 93);
    const auto *pfd_96 = buffer.data(pfd + 96);
    const auto *pfd_99 = buffer.data(pfd + 99);
    const auto *pfd_101 = buffer.data(pfd + 101);
    const auto *pfd_102 = buffer.data(pfd + 102);
    const auto *pfd_105 = buffer.data(pfd + 105);
    const auto *pfd_107 = buffer.data(pfd + 107);
    const auto *pfd_108 = buffer.data(pfd + 108);
    const auto *pfd_111 = buffer.data(pfd + 111);
    const auto *pfd_113 = buffer.data(pfd + 113);
    const auto *pfd_114 = buffer.data(pfd + 114);
    const auto *pfd_117 = buffer.data(pfd + 117);
    const auto *pfd_119 = buffer.data(pfd + 119);
    const auto *pfd_120 = buffer.data(pfd + 120);
    const auto *pfd_123 = buffer.data(pfd + 123);
    const auto *pfd_124 = buffer.data(pfd + 124);
    const auto *pfd_125 = buffer.data(pfd + 125);
    const auto *pfd_126 = buffer.data(pfd + 126);
    const auto *pfd_129 = buffer.data(pfd + 129);
    const auto *pfd_130 = buffer.data(pfd + 130);
    const auto *pfd_131 = buffer.data(pfd + 131);
    const auto *pfd_132 = buffer.data(pfd + 132);
    const auto *pfd_134 = buffer.data(pfd + 134);
    const auto *pfd_135 = buffer.data(pfd + 135);
    const auto *pfd_136 = buffer.data(pfd + 136);
    const auto *pfd_137 = buffer.data(pfd + 137);
    const auto *pfd_138 = buffer.data(pfd + 138);
    const auto *pfd_141 = buffer.data(pfd + 141);
    const auto *pfd_142 = buffer.data(pfd + 142);
    const auto *pfd_143 = buffer.data(pfd + 143);
    const auto *pfd_147 = buffer.data(pfd + 147);
    const auto *pfd_148 = buffer.data(pfd + 148);
    const auto *pfd_149 = buffer.data(pfd + 149);
    const auto *pfd_150 = buffer.data(pfd + 150);
    const auto *pfd_152 = buffer.data(pfd + 152);
    const auto *pfd_153 = buffer.data(pfd + 153);
    const auto *pfd_154 = buffer.data(pfd + 154);
    const auto *pfd_155 = buffer.data(pfd + 155);
    const auto *pfd_159 = buffer.data(pfd + 159);
    const auto *pfd_160 = buffer.data(pfd + 160);
    const auto *pfd_161 = buffer.data(pfd + 161);

    const auto *pff1_160 = buffer.data(pff1 + 160);
    const auto *pff1_161 = buffer.data(pff1 + 161);
    const auto *pff1_166 = buffer.data(pff1 + 166);
    const auto *pff1_176 = buffer.data(pff1 + 176);
    const auto *pff1_202 = buffer.data(pff1 + 202);
    const auto *pff1_209 = buffer.data(pff1 + 209);
    const auto *pff1_220 = buffer.data(pff1 + 220);
    const auto *pff1_222 = buffer.data(pff1 + 222);
    const auto *pff1_229 = buffer.data(pff1 + 229);
    const auto *pff1_259 = buffer.data(pff1 + 259);
    const auto *pff1_267 = buffer.data(pff1 + 267);

    const auto *pgp0_75 = buffer.data(pgp0 + 75);
    const auto *pgp0_76 = buffer.data(pgp0 + 76);
    const auto *pgp0_77 = buffer.data(pgp0 + 77);
    const auto *pgp0_80 = buffer.data(pgp0 + 80);
    const auto *pgp0_81 = buffer.data(pgp0 + 81);
    const auto *pgp0_82 = buffer.data(pgp0 + 82);
    const auto *pgp0_83 = buffer.data(pgp0 + 83);
    const auto *pgp0_84 = buffer.data(pgp0 + 84);
    const auto *pgp0_85 = buffer.data(pgp0 + 85);
    const auto *pgp0_86 = buffer.data(pgp0 + 86);
    const auto *pgp0_88 = buffer.data(pgp0 + 88);
    const auto *pgp0_92 = buffer.data(pgp0 + 92);
    const auto *pgp0_96 = buffer.data(pgp0 + 96);
    const auto *pgp0_97 = buffer.data(pgp0 + 97);
    const auto *pgp0_98 = buffer.data(pgp0 + 98);
    const auto *pgp0_103 = buffer.data(pgp0 + 103);
    const auto *pgp0_104 = buffer.data(pgp0 + 104);
    const auto *pgp0_105 = buffer.data(pgp0 + 105);
    const auto *pgp0_106 = buffer.data(pgp0 + 106);
    const auto *pgp0_107 = buffer.data(pgp0 + 107);

    const auto *pgp1_75 = buffer.data(pgp1 + 75);
    const auto *pgp1_76 = buffer.data(pgp1 + 76);
    const auto *pgp1_77 = buffer.data(pgp1 + 77);
    const auto *pgp1_80 = buffer.data(pgp1 + 80);
    const auto *pgp1_81 = buffer.data(pgp1 + 81);
    const auto *pgp1_82 = buffer.data(pgp1 + 82);
    const auto *pgp1_83 = buffer.data(pgp1 + 83);
    const auto *pgp1_84 = buffer.data(pgp1 + 84);
    const auto *pgp1_85 = buffer.data(pgp1 + 85);
    const auto *pgp1_86 = buffer.data(pgp1 + 86);
    const auto *pgp1_88 = buffer.data(pgp1 + 88);
    const auto *pgp1_92 = buffer.data(pgp1 + 92);
    const auto *pgp1_96 = buffer.data(pgp1 + 96);
    const auto *pgp1_97 = buffer.data(pgp1 + 97);
    const auto *pgp1_98 = buffer.data(pgp1 + 98);
    const auto *pgp1_103 = buffer.data(pgp1 + 103);
    const auto *pgp1_104 = buffer.data(pgp1 + 104);
    const auto *pgp1_105 = buffer.data(pgp1 + 105);
    const auto *pgp1_106 = buffer.data(pgp1 + 106);
    const auto *pgp1_107 = buffer.data(pgp1 + 107);

    const auto *pgd_147 = buffer.data(pgd + 147);
    const auto *pgd_149 = buffer.data(pgd + 149);
    const auto *pgd_150 = buffer.data(pgd + 150);
    const auto *pgd_151 = buffer.data(pgd + 151);
    const auto *pgd_153 = buffer.data(pgd + 153);
    const auto *pgd_154 = buffer.data(pgd + 154);
    const auto *pgd_155 = buffer.data(pgd + 155);
    const auto *pgd_156 = buffer.data(pgd + 156);
    const auto *pgd_159 = buffer.data(pgd + 159);
    const auto *pgd_160 = buffer.data(pgd + 160);
    const auto *pgd_161 = buffer.data(pgd + 161);
    const auto *pgd_162 = buffer.data(pgd + 162);
    const auto *pgd_163 = buffer.data(pgd + 163);
    const auto *pgd_165 = buffer.data(pgd + 165);
    const auto *pgd_166 = buffer.data(pgd + 166);
    const auto *pgd_167 = buffer.data(pgd + 167);
    const auto *pgd_168 = buffer.data(pgd + 168);
    const auto *pgd_169 = buffer.data(pgd + 169);
    const auto *pgd_171 = buffer.data(pgd + 171);
    const auto *pgd_172 = buffer.data(pgd + 172);
    const auto *pgd_173 = buffer.data(pgd + 173);
    const auto *pgd_174 = buffer.data(pgd + 174);
    const auto *pgd_175 = buffer.data(pgd + 175);
    const auto *pgd_177 = buffer.data(pgd + 177);
    const auto *pgd_178 = buffer.data(pgd + 178);
    const auto *pgd_179 = buffer.data(pgd + 179);
    const auto *pgd_180 = buffer.data(pgd + 180);
    const auto *pgd_183 = buffer.data(pgd + 183);
    const auto *pgd_184 = buffer.data(pgd + 184);
    const auto *pgd_185 = buffer.data(pgd + 185);
    const auto *pgd_186 = buffer.data(pgd + 186);
    const auto *pgd_189 = buffer.data(pgd + 189);
    const auto *pgd_190 = buffer.data(pgd + 190);
    const auto *pgd_191 = buffer.data(pgd + 191);
    const auto *pgd_192 = buffer.data(pgd + 192);
    const auto *pgd_194 = buffer.data(pgd + 194);
    const auto *pgd_195 = buffer.data(pgd + 195);
    const auto *pgd_196 = buffer.data(pgd + 196);
    const auto *pgd_197 = buffer.data(pgd + 197);
    const auto *pgd_198 = buffer.data(pgd + 198);
    const auto *pgd_201 = buffer.data(pgd + 201);
    const auto *pgd_202 = buffer.data(pgd + 202);
    const auto *pgd_203 = buffer.data(pgd + 203);
    const auto *pgd_204 = buffer.data(pgd + 204);
    const auto *pgd_207 = buffer.data(pgd + 207);
    const auto *pgd_208 = buffer.data(pgd + 208);
    const auto *pgd_209 = buffer.data(pgd + 209);
    const auto *pgd_210 = buffer.data(pgd + 210);
    const auto *pgd_212 = buffer.data(pgd + 212);
    const auto *pgd_213 = buffer.data(pgd + 213);
    const auto *pgd_214 = buffer.data(pgd + 214);
    const auto *pgd_215 = buffer.data(pgd + 215);
    const auto *pgd_216 = buffer.data(pgd + 216);
    const auto *pgd_219 = buffer.data(pgd + 219);
    const auto *pgd_220 = buffer.data(pgd + 220);
    const auto *pgd_221 = buffer.data(pgd + 221);

#pragma omp simd aligned(t_247, t_248, t_249, pa_y, pc_y, pc_z, sgf0_99, sgd_59, sgf1_99, \
                         pfd_93, pgd_147, pgd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_6 * pfd_93[k]
                   + f_4 * pc_z[k] * pgd_147[k];

        t_248[k] = f_0 * sgd_59[k]
                   + f_4 * pc_y[k] * pgd_149[k];

        t_249[k] = pa_y[k] * sgf0_99[k]
                   - f_5 * pc_y[k] * sgf1_99[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, pc_x, pc_z, pgp0_75, pgp0_76, \
                         pgp1_75, pgp1_76, pgd_150, pgd_151, pgd_153, \
                         pgd_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_2 * pgp0_75[k]
                   - f_3 * pgp1_75[k]
                   + f_4 * pc_x[k] * pgd_150[k];

        t_251[k] = f_10 * pgp0_76[k]
                   - f_11 * pgp1_76[k]
                   + f_4 * pc_x[k] * pgd_151[k];

        t_252[k] = f_4 * pc_z[k] * pgd_150[k];

        t_253[k] = f_4 * pc_x[k] * pgd_153[k];

        t_254[k] = f_4 * pc_x[k] * pgd_154[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pc_x, pc_y, pc_z, sgd_63, sgd_65, pfd_99, \
                         pfd_101, pgp0_76, pgp1_76, pgd_153, pgd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_4 * pc_x[k] * pgd_155[k];

        t_256[k] = f_0 * sgd_63[k]
                   + f_1 * pfd_99[k]
                   + f_2 * pgp0_76[k]
                   - f_3 * pgp1_76[k]
                   + f_4 * pc_y[k] * pgd_153[k];

        t_257[k] = f_4 * pc_z[k] * pgd_153[k];

        t_258[k] = f_0 * sgd_65[k]
                   + f_1 * pfd_101[k]
                   + f_4 * pc_y[k] * pgd_155[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pb_z, pc_z, pff0_160, pff0_161, pfd_96, \
                         pff1_160, pff1_161, pgp0_77, pgp1_77, pgd_155, \
                         pgd_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_2 * pgp0_77[k]
                   - f_3 * pgp1_77[k]
                   + f_4 * pc_z[k] * pgd_155[k];

        t_260[k] = pb_z[k] * pff0_160[k]
                   - f_5 * pc_z[k] * pff1_160[k];

        t_261[k] = pb_z[k] * pff0_161[k]
                   - f_5 * pc_z[k] * pff1_161[k];

        t_262[k] = f_0 * pfd_96[k]
                   + f_4 * pc_z[k] * pgd_156[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, pb_z, pc_x, pc_z, pff0_166, \
                         pfd_99, pff1_166, pgd_159, pgd_160, pgd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_4 * pc_x[k] * pgd_159[k];

        t_264[k] = f_4 * pc_x[k] * pgd_160[k];

        t_265[k] = f_4 * pc_x[k] * pgd_161[k];

        t_266[k] = pb_z[k] * pff0_166[k]
                   - f_5 * pc_z[k] * pff1_166[k];

        t_267[k] = f_0 * pfd_99[k]
                   + f_4 * pc_z[k] * pgd_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_x, pc_y, pc_z, sgd_71, pfd_101, pfd_107, \
                         pgp0_80, pgp0_81, pgp1_80, pgp1_81, pgd_161, \
                         pgd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_0 * sgd_71[k]
                   + f_6 * pfd_107[k]
                   + f_4 * pc_y[k] * pgd_161[k];

        t_269[k] = f_0 * pfd_101[k]
                   + f_2 * pgp0_80[k]
                   - f_3 * pgp1_80[k]
                   + f_4 * pc_z[k] * pgd_161[k];

        t_270[k] = f_2 * pgp0_81[k]
                   - f_3 * pgp1_81[k]
                   + f_4 * pc_x[k] * pgd_162[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, pc_x, pc_z, pfd_102, pgp0_82, \
                         pgp1_82, pgd_162, pgd_163, pgd_165, pgd_166, \
                         pgd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * pgp0_82[k]
                   - f_11 * pgp1_82[k]
                   + f_4 * pc_x[k] * pgd_163[k];

        t_272[k] = f_9 * pfd_102[k]
                   + f_4 * pc_z[k] * pgd_162[k];

        t_273[k] = f_4 * pc_x[k] * pgd_165[k];

        t_274[k] = f_4 * pc_x[k] * pgd_166[k];

        t_275[k] = f_4 * pc_x[k] * pgd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, pb_z, pc_y, pc_z, sgd_77, pdf0_96, pdf1_96, \
                         pff0_176, pfd_105, pfd_113, pff1_176, pgd_165, \
                         pgd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_7 * pdf0_96[k]
                   - f_8 * pdf1_96[k]
                   + pb_z[k] * pff0_176[k]
                   - f_5 * pc_z[k] * pff1_176[k];

        t_277[k] = f_9 * pfd_105[k]
                   + f_4 * pc_z[k] * pgd_165[k];

        t_278[k] = f_0 * sgd_77[k]
                   + f_9 * pfd_113[k]
                   + f_4 * pc_y[k] * pgd_167[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pc_x, pc_z, pfd_107, pgp0_83, pgp0_84, pgp0_85, \
                         pgp1_83, pgp1_84, pgp1_85, pgd_167, pgd_168, \
                         pgd_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_9 * pfd_107[k]
                   + f_2 * pgp0_83[k]
                   - f_3 * pgp1_83[k]
                   + f_4 * pc_z[k] * pgd_167[k];

        t_280[k] = f_2 * pgp0_84[k]
                   - f_3 * pgp1_84[k]
                   + f_4 * pc_x[k] * pgd_168[k];

        t_281[k] = f_10 * pgp0_85[k]
                   - f_11 * pgp1_85[k]
                   + f_4 * pc_x[k] * pgd_169[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pc_x, pc_z, pfd_108, pgd_168, pgd_171, \
                         pgd_172, pgd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_6 * pfd_108[k]
                   + f_4 * pc_z[k] * pgd_168[k];

        t_283[k] = f_4 * pc_x[k] * pgd_171[k];

        t_284[k] = f_4 * pc_x[k] * pgd_172[k];

        t_285[k] = f_4 * pc_x[k] * pgd_173[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, pc_y, pc_z, sgd_81, sgd_83, pfd_111, pfd_117, \
                         pfd_119, pgp0_85, pgp1_85, pgd_171, pgd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * sgd_81[k]
                   + f_0 * pfd_117[k]
                   + f_2 * pgp0_85[k]
                   - f_3 * pgp1_85[k]
                   + f_4 * pc_y[k] * pgd_171[k];

        t_287[k] = f_6 * pfd_111[k]
                   + f_4 * pc_z[k] * pgd_171[k];

        t_288[k] = f_0 * sgd_83[k]
                   + f_0 * pfd_119[k]
                   + f_4 * pc_y[k] * pgd_173[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, pa_y, pc_x, pc_y, pc_z, sgf0_140, sgf1_140, \
                         pfd_113, pgp0_86, pgp0_88, pgp1_86, pgp1_88, pgd_173, \
                         pgd_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_6 * pfd_113[k]
                   + f_2 * pgp0_86[k]
                   - f_3 * pgp1_86[k]
                   + f_4 * pc_z[k] * pgd_173[k];

        t_290[k] = pa_y[k] * sgf0_140[k]
                   - f_5 * pc_y[k] * sgf1_140[k];

        t_291[k] = f_10 * pgp0_88[k]
                   - f_11 * pgp1_88[k]
                   + f_4 * pc_x[k] * pgd_175[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_z, pfd_114, pgd_174, pgd_177, \
                         pgd_178, pgd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_1 * pfd_114[k]
                   + f_4 * pc_z[k] * pgd_174[k];

        t_293[k] = f_4 * pc_x[k] * pgd_177[k];

        t_294[k] = f_4 * pc_x[k] * pgd_178[k];

        t_295[k] = f_4 * pc_x[k] * pgd_179[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_y, pc_y, pc_z, sgf0_146, sgf0_149, \
                         sgd_87, sgd_89, sgf1_146, sgf1_149, pfd_117, pgd_177, \
                         pgd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = pa_y[k] * sgf0_146[k]
                   + f_6 * sgd_87[k]
                   - f_5 * pc_y[k] * sgf1_146[k];

        t_297[k] = f_1 * pfd_117[k]
                   + f_4 * pc_z[k] * pgd_177[k];

        t_298[k] = f_0 * sgd_89[k]
                   + f_4 * pc_y[k] * pgd_179[k];

        t_299[k] = pa_y[k] * sgf0_149[k]
                   - f_5 * pc_y[k] * sgf1_149[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_z, pc_x, pc_y, pc_z, sgf0_0, sgf0_2, \
                         sgd_0, sgf1_0, sgf1_2, pfd_123, pgd_180, \
                         pgd_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = pa_z[k] * sgf0_0[k]
                   - f_5 * pc_z[k] * sgf1_0[k];

        t_301[k] = f_4 * pc_y[k] * pgd_180[k];

        t_302[k] = pa_z[k] * sgf0_2[k]
                   + f_0 * sgd_0[k]
                   - f_5 * pc_z[k] * sgf1_2[k];

        t_303[k] = f_1 * pfd_123[k]
                   + f_4 * pc_x[k] * pgd_183[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pa_z, pc_x, pc_y, pc_z, sgf0_6, sgf1_6, \
                         pfd_124, pfd_125, pgp0_92, pgp1_92, pgd_184, \
                         pgd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_1 * pfd_124[k]
                   + f_4 * pc_x[k] * pgd_184[k];

        t_305[k] = f_1 * pfd_125[k]
                   + f_4 * pc_x[k] * pgd_185[k];

        t_306[k] = pa_z[k] * sgf0_6[k]
                   - f_5 * pc_z[k] * sgf1_6[k];

        t_307[k] = f_10 * pgp0_92[k]
                   - f_11 * pgp1_92[k]
                   + f_4 * pc_y[k] * pgd_184[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_z, pc_y, pc_z, sgf0_9, sgf0_10, sgd_5, \
                         sgf1_9, sgf1_10, pfd_120, pgd_185, pgd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_4 * pc_y[k] * pgd_185[k];

        t_309[k] = pa_z[k] * sgf0_9[k]
                   + f_6 * sgd_5[k]
                   - f_5 * pc_z[k] * sgf1_9[k];

        t_310[k] = pa_z[k] * sgf0_10[k]
                   - f_5 * pc_z[k] * sgf1_10[k];

        t_311[k] = f_0 * pfd_120[k]
                   + f_4 * pc_y[k] * pgd_186[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pb_y, pc_x, pc_y, pff0_202, pfd_129, \
                         pfd_130, pfd_131, pff1_202, pgd_189, pgd_190, \
                         pgd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = pb_y[k] * pff0_202[k]
                   - f_5 * pc_y[k] * pff1_202[k];

        t_313[k] = f_6 * pfd_129[k]
                   + f_4 * pc_x[k] * pgd_189[k];

        t_314[k] = f_6 * pfd_130[k]
                   + f_4 * pc_x[k] * pgd_190[k];

        t_315[k] = f_6 * pfd_131[k]
                   + f_4 * pc_x[k] * pgd_191[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_z, pc_y, pc_z, sgf0_16, sgf0_17, sgd_9, \
                         sgf1_16, sgf1_17, pfd_125, pgd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_z[k] * sgf0_16[k]
                   - f_5 * pc_z[k] * sgf1_16[k];

        t_317[k] = pa_z[k] * sgf0_17[k]
                   + f_0 * sgd_9[k]
                   - f_5 * pc_z[k] * sgf1_17[k];

        t_318[k] = f_0 * pfd_125[k]
                   + f_4 * pc_y[k] * pgd_191[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pb_y, pc_x, pc_y, pff0_209, pfd_132, pff1_209, \
                         pgp0_96, pgp1_96, pgd_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = pb_y[k] * pff0_209[k]
                   - f_5 * pc_y[k] * pff1_209[k];

        t_320[k] = f_6 * pfd_132[k]
                   + f_2 * pgp0_96[k]
                   - f_3 * pgp1_96[k]
                   + f_4 * pc_x[k] * pgd_192[k];

        t_321[k] = f_4 * pc_y[k] * pgd_192[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, pfd_134, pfd_135, pfd_136, pfd_137, \
                         pgp0_98, pgp1_98, pgd_194, pgd_195, pgd_196, \
                         pgd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_6 * pfd_134[k]
                   + f_10 * pgp0_98[k]
                   - f_11 * pgp1_98[k]
                   + f_4 * pc_x[k] * pgd_194[k];

        t_323[k] = f_6 * pfd_135[k]
                   + f_4 * pc_x[k] * pgd_195[k];

        t_324[k] = f_6 * pfd_136[k]
                   + f_4 * pc_x[k] * pgd_196[k];

        t_325[k] = f_6 * pfd_137[k]
                   + f_4 * pc_x[k] * pgd_197[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, pc_y, pgp0_97, pgp0_98, pgp1_97, pgp1_98, \
                         pgd_195, pgd_196, pgd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_2 * pgp0_97[k]
                   - f_3 * pgp1_97[k]
                   + f_4 * pc_y[k] * pgd_195[k];

        t_327[k] = f_10 * pgp0_98[k]
                   - f_11 * pgp1_98[k]
                   + f_4 * pc_y[k] * pgd_196[k];

        t_328[k] = f_4 * pc_y[k] * pgd_197[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pa_z, pb_x, pc_x, pc_y, pc_z, sgf0_30, sgf1_30, \
                         pdf0_149, pdf1_149, pff0_229, pfd_126, pff1_229, \
                         pgd_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_12 * pdf0_149[k]
                   - f_13 * pdf1_149[k]
                   + pb_x[k] * pff0_229[k]
                   - f_5 * pc_x[k] * pff1_229[k];

        t_330[k] = pa_z[k] * sgf0_30[k]
                   - f_5 * pc_z[k] * sgf1_30[k];

        t_331[k] = f_9 * pfd_126[k]
                   + f_4 * pc_y[k] * pgd_198[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pa_z, pc_x, pc_z, sgf0_32, sgd_18, \
                         sgf1_32, pfd_141, pfd_142, pfd_143, pgd_201, pgd_202, \
                         pgd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pa_z[k] * sgf0_32[k]
                   + f_0 * sgd_18[k]
                   - f_5 * pc_z[k] * sgf1_32[k];

        t_333[k] = f_9 * pfd_141[k]
                   + f_4 * pc_x[k] * pgd_201[k];

        t_334[k] = f_9 * pfd_142[k]
                   + f_4 * pc_x[k] * pgd_202[k];

        t_335[k] = f_9 * pfd_143[k]
                   + f_4 * pc_x[k] * pgd_203[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pa_z, pc_y, pc_z, sgf0_36, sgf0_37, sgd_21, \
                         sgf1_36, sgf1_37, pfd_131, pgd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_z[k] * sgf0_36[k]
                   - f_5 * pc_z[k] * sgf1_36[k];

        t_337[k] = pa_z[k] * sgf0_37[k]
                   + f_0 * sgd_21[k]
                   - f_5 * pc_z[k] * sgf1_37[k];

        t_338[k] = f_9 * pfd_131[k]
                   + f_4 * pc_y[k] * pgd_203[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_z, pb_y, pc_y, pc_z, sgf0_39, sgd_23, \
                         sgf1_39, pff0_220, pfd_132, pff1_220, \
                         pgd_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * sgf0_39[k]
                   + f_6 * sgd_23[k]
                   - f_5 * pc_z[k] * sgf1_39[k];

        t_340[k] = pb_y[k] * pff0_220[k]
                   - f_5 * pc_y[k] * pff1_220[k];

        t_341[k] = f_0 * pfd_132[k]
                   + f_4 * pc_y[k] * pgd_204[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pb_y, pc_x, pc_y, pff0_222, pfd_147, \
                         pfd_148, pfd_149, pff1_222, pgd_207, pgd_208, \
                         pgd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pb_y[k] * pff0_222[k]
                   - f_5 * pc_y[k] * pff1_222[k];

        t_343[k] = f_9 * pfd_147[k]
                   + f_4 * pc_x[k] * pgd_207[k];

        t_344[k] = f_9 * pfd_148[k]
                   + f_4 * pc_x[k] * pgd_208[k];

        t_345[k] = f_9 * pfd_149[k]
                   + f_4 * pc_x[k] * pgd_209[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, pc_y, pfd_135, pfd_136, pfd_137, pgp0_103, \
                         pgp0_104, pgp1_103, pgp1_104, pgd_207, pgd_208, \
                         pgd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_0 * pfd_135[k]
                   + f_2 * pgp0_103[k]
                   - f_3 * pgp1_103[k]
                   + f_4 * pc_y[k] * pgd_207[k];

        t_347[k] = f_0 * pfd_136[k]
                   + f_10 * pgp0_104[k]
                   - f_11 * pgp1_104[k]
                   + f_4 * pc_y[k] * pgd_208[k];

        t_348[k] = f_0 * pfd_137[k]
                   + f_4 * pc_y[k] * pgd_209[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, pb_y, pc_x, pc_y, pff0_229, pfd_150, pff1_229, \
                         pgp0_105, pgp1_105, pgd_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = pb_y[k] * pff0_229[k]
                   - f_5 * pc_y[k] * pff1_229[k];

        t_350[k] = f_9 * pfd_150[k]
                   + f_2 * pgp0_105[k]
                   - f_3 * pgp1_105[k]
                   + f_4 * pc_x[k] * pgd_210[k];

        t_351[k] = f_4 * pc_y[k] * pgd_210[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_x, pfd_152, pfd_153, pfd_154, pfd_155, \
                         pgp0_107, pgp1_107, pgd_212, pgd_213, pgd_214, \
                         pgd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_9 * pfd_152[k]
                   + f_10 * pgp0_107[k]
                   - f_11 * pgp1_107[k]
                   + f_4 * pc_x[k] * pgd_212[k];

        t_353[k] = f_9 * pfd_153[k]
                   + f_4 * pc_x[k] * pgd_213[k];

        t_354[k] = f_9 * pfd_154[k]
                   + f_4 * pc_x[k] * pgd_214[k];

        t_355[k] = f_9 * pfd_155[k]
                   + f_4 * pc_x[k] * pgd_215[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, pgp0_106, pgp0_107, pgp1_106, pgp1_107, \
                         pgd_213, pgd_214, pgd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_2 * pgp0_106[k]
                   - f_3 * pgp1_106[k]
                   + f_4 * pc_y[k] * pgd_213[k];

        t_357[k] = f_10 * pgp0_107[k]
                   - f_11 * pgp1_107[k]
                   + f_4 * pc_y[k] * pgd_214[k];

        t_358[k] = f_4 * pc_y[k] * pgd_215[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pa_z, pb_x, pc_x, pc_y, pc_z, sgf0_60, sgf1_60, \
                         pdf0_179, pdf1_179, pff0_259, pfd_138, pff1_259, \
                         pgd_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_7 * pdf0_179[k]
                   - f_8 * pdf1_179[k]
                   + pb_x[k] * pff0_259[k]
                   - f_5 * pc_x[k] * pff1_259[k];

        t_360[k] = pa_z[k] * sgf0_60[k]
                   - f_5 * pc_z[k] * sgf1_60[k];

        t_361[k] = f_6 * pfd_138[k]
                   + f_4 * pc_y[k] * pgd_216[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_z, pc_x, pc_z, sgf0_62, sgd_36, \
                         sgf1_62, pfd_159, pfd_160, pfd_161, pgd_219, pgd_220, \
                         pgd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = pa_z[k] * sgf0_62[k]
                   + f_0 * sgd_36[k]
                   - f_5 * pc_z[k] * sgf1_62[k];

        t_363[k] = f_0 * pfd_159[k]
                   + f_4 * pc_x[k] * pgd_219[k];

        t_364[k] = f_0 * pfd_160[k]
                   + f_4 * pc_x[k] * pgd_220[k];

        t_365[k] = f_0 * pfd_161[k]
                   + f_4 * pc_x[k] * pgd_221[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, pa_z, pb_x, pc_x, pc_y, pc_z, sgf0_66, sgf1_66, \
                         pff0_267, pfd_143, pff1_267, pgd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = pa_z[k] * sgf0_66[k]
                   - f_5 * pc_z[k] * sgf1_66[k];

        t_367[k] = pb_x[k] * pff0_267[k]
                   - f_5 * pc_x[k] * pff1_267[k];

        t_368[k] = f_6 * pfd_143[k]
                   + f_4 * pc_y[k] * pgd_221[k];
    }
}

static auto
compute_prim_pgf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgf0, const size_t sgd,
                                                          const size_t sgf1, const size_t pdf0,
                                                          const size_t pdf1, const size_t pff0,
                                                          const size_t pfd, const size_t pff1,
                                                          const size_t pgp0, const size_t pgp1,
                                                          const size_t pgd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);
    const auto f_9 = 1.0 / q;
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = 1.0 / p;
    const auto f_13 = gamma / (p * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgf0_100 = buffer.data(sgf0 + 100);
    const auto *sgf0_106 = buffer.data(sgf0 + 106);
    const auto *sgf0_107 = buffer.data(sgf0 + 107);
    const auto *sgf0_109 = buffer.data(sgf0 + 109);

    const auto *sgd_63 = buffer.data(sgd + 63);
    const auto *sgd_65 = buffer.data(sgd + 65);
    const auto *sgd_89 = buffer.data(sgd + 89);

    const auto *sgf1_100 = buffer.data(sgf1 + 100);
    const auto *sgf1_106 = buffer.data(sgf1 + 106);
    const auto *sgf1_107 = buffer.data(sgf1 + 107);
    const auto *sgf1_109 = buffer.data(sgf1 + 109);

    const auto *pdf0_169 = buffer.data(pdf0 + 169);
    const auto *pdf0_179 = buffer.data(pdf0 + 179);

    const auto *pdf1_169 = buffer.data(pdf1 + 169);
    const auto *pdf1_179 = buffer.data(pdf1 + 179);

    const auto *pff0_250 = buffer.data(pff0 + 250);
    const auto *pff0_252 = buffer.data(pff0 + 252);
    const auto *pff0_269 = buffer.data(pff0 + 269);
    const auto *pff0_272 = buffer.data(pff0 + 272);
    const auto *pff0_276 = buffer.data(pff0 + 276);
    const auto *pff0_277 = buffer.data(pff0 + 277);
    const auto *pff0_279 = buffer.data(pff0 + 279);
    const auto *pff0_286 = buffer.data(pff0 + 286);
    const auto *pff0_287 = buffer.data(pff0 + 287);
    const auto *pff0_289 = buffer.data(pff0 + 289);
    const auto *pff0_290 = buffer.data(pff0 + 290);
    const auto *pff0_292 = buffer.data(pff0 + 292);
    const auto *pff0_296 = buffer.data(pff0 + 296);
    const auto *pff0_297 = buffer.data(pff0 + 297);
    const auto *pff0_299 = buffer.data(pff0 + 299);

    const auto *pfd_144 = buffer.data(pfd + 144);
    const auto *pfd_149 = buffer.data(pfd + 149);
    const auto *pfd_150 = buffer.data(pfd + 150);
    const auto *pfd_155 = buffer.data(pfd + 155);
    const auto *pfd_156 = buffer.data(pfd + 156);
    const auto *pfd_161 = buffer.data(pfd + 161);
    const auto *pfd_162 = buffer.data(pfd + 162);
    const auto *pfd_164 = buffer.data(pfd + 164);
    const auto *pfd_165 = buffer.data(pfd + 165);
    const auto *pfd_166 = buffer.data(pfd + 166);
    const auto *pfd_167 = buffer.data(pfd + 167);
    const auto *pfd_168 = buffer.data(pfd + 168);
    const auto *pfd_171 = buffer.data(pfd + 171);
    const auto *pfd_172 = buffer.data(pfd + 172);
    const auto *pfd_173 = buffer.data(pfd + 173);
    const auto *pfd_174 = buffer.data(pfd + 174);
    const auto *pfd_176 = buffer.data(pfd + 176);
    const auto *pfd_177 = buffer.data(pfd + 177);
    const auto *pfd_178 = buffer.data(pfd + 178);
    const auto *pfd_179 = buffer.data(pfd + 179);

    const auto *pff1_250 = buffer.data(pff1 + 250);
    const auto *pff1_252 = buffer.data(pff1 + 252);
    const auto *pff1_269 = buffer.data(pff1 + 269);
    const auto *pff1_272 = buffer.data(pff1 + 272);
    const auto *pff1_276 = buffer.data(pff1 + 276);
    const auto *pff1_277 = buffer.data(pff1 + 277);
    const auto *pff1_279 = buffer.data(pff1 + 279);
    const auto *pff1_286 = buffer.data(pff1 + 286);
    const auto *pff1_287 = buffer.data(pff1 + 287);
    const auto *pff1_289 = buffer.data(pff1 + 289);
    const auto *pff1_290 = buffer.data(pff1 + 290);
    const auto *pff1_292 = buffer.data(pff1 + 292);
    const auto *pff1_296 = buffer.data(pff1 + 296);
    const auto *pff1_297 = buffer.data(pff1 + 297);
    const auto *pff1_299 = buffer.data(pff1 + 299);

    const auto *pgp0_111 = buffer.data(pgp0 + 111);
    const auto *pgp0_122 = buffer.data(pgp0 + 122);
    const auto *pgp0_123 = buffer.data(pgp0 + 123);
    const auto *pgp0_124 = buffer.data(pgp0 + 124);
    const auto *pgp0_125 = buffer.data(pgp0 + 125);
    const auto *pgp0_126 = buffer.data(pgp0 + 126);
    const auto *pgp0_127 = buffer.data(pgp0 + 127);
    const auto *pgp0_128 = buffer.data(pgp0 + 128);
    const auto *pgp0_132 = buffer.data(pgp0 + 132);
    const auto *pgp0_133 = buffer.data(pgp0 + 133);
    const auto *pgp0_134 = buffer.data(pgp0 + 134);

    const auto *pgp1_111 = buffer.data(pgp1 + 111);
    const auto *pgp1_122 = buffer.data(pgp1 + 122);
    const auto *pgp1_123 = buffer.data(pgp1 + 123);
    const auto *pgp1_124 = buffer.data(pgp1 + 124);
    const auto *pgp1_125 = buffer.data(pgp1 + 125);
    const auto *pgp1_126 = buffer.data(pgp1 + 126);
    const auto *pgp1_127 = buffer.data(pgp1 + 127);
    const auto *pgp1_128 = buffer.data(pgp1 + 128);
    const auto *pgp1_132 = buffer.data(pgp1 + 132);
    const auto *pgp1_133 = buffer.data(pgp1 + 133);
    const auto *pgp1_134 = buffer.data(pgp1 + 134);

    const auto *pgd_222 = buffer.data(pgd + 222);
    const auto *pgd_225 = buffer.data(pgd + 225);
    const auto *pgd_226 = buffer.data(pgd + 226);
    const auto *pgd_227 = buffer.data(pgd + 227);
    const auto *pgd_228 = buffer.data(pgd + 228);
    const auto *pgd_231 = buffer.data(pgd + 231);
    const auto *pgd_232 = buffer.data(pgd + 232);
    const auto *pgd_233 = buffer.data(pgd + 233);
    const auto *pgd_234 = buffer.data(pgd + 234);
    const auto *pgd_237 = buffer.data(pgd + 237);
    const auto *pgd_238 = buffer.data(pgd + 238);
    const auto *pgd_239 = buffer.data(pgd + 239);
    const auto *pgd_240 = buffer.data(pgd + 240);
    const auto *pgd_242 = buffer.data(pgd + 242);
    const auto *pgd_243 = buffer.data(pgd + 243);
    const auto *pgd_244 = buffer.data(pgd + 244);
    const auto *pgd_245 = buffer.data(pgd + 245);
    const auto *pgd_246 = buffer.data(pgd + 246);
    const auto *pgd_248 = buffer.data(pgd + 248);
    const auto *pgd_249 = buffer.data(pgd + 249);
    const auto *pgd_250 = buffer.data(pgd + 250);
    const auto *pgd_251 = buffer.data(pgd + 251);
    const auto *pgd_252 = buffer.data(pgd + 252);
    const auto *pgd_254 = buffer.data(pgd + 254);
    const auto *pgd_255 = buffer.data(pgd + 255);
    const auto *pgd_256 = buffer.data(pgd + 256);
    const auto *pgd_257 = buffer.data(pgd + 257);
    const auto *pgd_258 = buffer.data(pgd + 258);
    const auto *pgd_261 = buffer.data(pgd + 261);
    const auto *pgd_262 = buffer.data(pgd + 262);
    const auto *pgd_263 = buffer.data(pgd + 263);
    const auto *pgd_264 = buffer.data(pgd + 264);
    const auto *pgd_266 = buffer.data(pgd + 266);
    const auto *pgd_267 = buffer.data(pgd + 267);
    const auto *pgd_268 = buffer.data(pgd + 268);
    const auto *pgd_269 = buffer.data(pgd + 269);

#pragma omp simd aligned(t_369, t_370, t_371, pb_x, pc_x, pc_y, pff0_269, pfd_144, pfd_162, \
                         pff1_269, pgp0_111, pgp1_111, pgd_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pb_x[k] * pff0_269[k]
                   - f_5 * pc_x[k] * pff1_269[k];

        t_370[k] = f_0 * pfd_162[k]
                   + f_2 * pgp0_111[k]
                   - f_3 * pgp1_111[k]
                   + f_4 * pc_x[k] * pgd_222[k];

        t_371[k] = f_9 * pfd_144[k]
                   + f_4 * pc_y[k] * pgd_222[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, pb_x, pc_x, pff0_272, pfd_164, pfd_165, \
                         pfd_166, pfd_167, pff1_272, pgd_225, pgd_226, \
                         pgd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = pb_x[k] * pff0_272[k]
                   + f_9 * pfd_164[k]
                   - f_5 * pc_x[k] * pff1_272[k];

        t_373[k] = f_0 * pfd_165[k]
                   + f_4 * pc_x[k] * pgd_225[k];

        t_374[k] = f_0 * pfd_166[k]
                   + f_4 * pc_x[k] * pgd_226[k];

        t_375[k] = f_0 * pfd_167[k]
                   + f_4 * pc_x[k] * pgd_227[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pb_x, pc_x, pc_y, pff0_276, pff0_277, \
                         pff0_279, pfd_149, pff1_276, pff1_277, pff1_279, \
                         pgd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pb_x[k] * pff0_276[k]
                   - f_5 * pc_x[k] * pff1_276[k];

        t_377[k] = pb_x[k] * pff0_277[k]
                   - f_5 * pc_x[k] * pff1_277[k];

        t_378[k] = f_9 * pfd_149[k]
                   + f_4 * pc_y[k] * pgd_227[k];

        t_379[k] = pb_x[k] * pff0_279[k]
                   - f_5 * pc_x[k] * pff1_279[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pb_y, pc_x, pc_y, pff0_250, pff0_252, \
                         pfd_150, pfd_171, pff1_250, pff1_252, pgd_228, \
                         pgd_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = pb_y[k] * pff0_250[k]
                   - f_5 * pc_y[k] * pff1_250[k];

        t_381[k] = f_0 * pfd_150[k]
                   + f_4 * pc_y[k] * pgd_228[k];

        t_382[k] = pb_y[k] * pff0_252[k]
                   - f_5 * pc_y[k] * pff1_252[k];

        t_383[k] = f_0 * pfd_171[k]
                   + f_4 * pc_x[k] * pgd_231[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pb_x, pc_x, pff0_286, pff0_287, pfd_172, \
                         pfd_173, pff1_286, pff1_287, pgd_232, \
                         pgd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_0 * pfd_172[k]
                   + f_4 * pc_x[k] * pgd_232[k];

        t_385[k] = f_0 * pfd_173[k]
                   + f_4 * pc_x[k] * pgd_233[k];

        t_386[k] = pb_x[k] * pff0_286[k]
                   - f_5 * pc_x[k] * pff1_286[k];

        t_387[k] = pb_x[k] * pff0_287[k]
                   - f_5 * pc_x[k] * pff1_287[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pb_x, pc_x, pc_y, pff0_289, pff0_290, \
                         pfd_155, pfd_174, pff1_289, pff1_290, pgd_233, \
                         pgd_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_0 * pfd_155[k]
                   + f_4 * pc_y[k] * pgd_233[k];

        t_389[k] = pb_x[k] * pff0_289[k]
                   - f_5 * pc_x[k] * pff1_289[k];

        t_390[k] = pb_x[k] * pff0_290[k]
                   + f_6 * pfd_174[k]
                   - f_5 * pc_x[k] * pff1_290[k];

        t_391[k] = f_4 * pc_y[k] * pgd_234[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pb_x, pc_x, pff0_292, pfd_176, pfd_177, \
                         pfd_178, pfd_179, pff1_292, pgd_237, pgd_238, \
                         pgd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = pb_x[k] * pff0_292[k]
                   + f_9 * pfd_176[k]
                   - f_5 * pc_x[k] * pff1_292[k];

        t_393[k] = f_0 * pfd_177[k]
                   + f_4 * pc_x[k] * pgd_237[k];

        t_394[k] = f_0 * pfd_178[k]
                   + f_4 * pc_x[k] * pgd_238[k];

        t_395[k] = f_0 * pfd_179[k]
                   + f_4 * pc_x[k] * pgd_239[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pb_x, pc_x, pc_y, pff0_296, pff0_297, \
                         pff0_299, pff1_296, pff1_297, pff1_299, \
                         pgd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = pb_x[k] * pff0_296[k]
                   - f_5 * pc_x[k] * pff1_296[k];

        t_397[k] = pb_x[k] * pff0_297[k]
                   - f_5 * pc_x[k] * pff1_297[k];

        t_398[k] = f_4 * pc_y[k] * pgd_239[k];

        t_399[k] = pb_x[k] * pff0_299[k]
                   - f_5 * pc_x[k] * pff1_299[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, pa_z, pc_x, pc_y, pc_z, sgf0_100, \
                         sgf1_100, pfd_156, pgp0_122, pgp1_122, pgd_240, pgd_242, \
                         pgd_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = pa_z[k] * sgf0_100[k]
                   - f_5 * pc_z[k] * sgf1_100[k];

        t_401[k] = f_1 * pfd_156[k]
                   + f_4 * pc_y[k] * pgd_240[k];

        t_402[k] = f_10 * pgp0_122[k]
                   - f_11 * pgp1_122[k]
                   + f_4 * pc_x[k] * pgd_242[k];

        t_403[k] = f_4 * pc_x[k] * pgd_243[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_z, pc_x, pc_z, sgf0_106, sgf0_107, \
                         sgd_63, sgf1_106, sgf1_107, pgd_244, pgd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_4 * pc_x[k] * pgd_244[k];

        t_405[k] = f_4 * pc_x[k] * pgd_245[k];

        t_406[k] = pa_z[k] * sgf0_106[k]
                   - f_5 * pc_z[k] * sgf1_106[k];

        t_407[k] = pa_z[k] * sgf0_107[k]
                   + f_0 * sgd_63[k]
                   - f_5 * pc_z[k] * sgf1_107[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pa_z, pc_x, pc_y, pc_z, sgf0_109, sgd_65, \
                         sgf1_109, pfd_161, pgp0_123, pgp1_123, pgd_245, \
                         pgd_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_1 * pfd_161[k]
                   + f_4 * pc_y[k] * pgd_245[k];

        t_409[k] = pa_z[k] * sgf0_109[k]
                   + f_6 * sgd_65[k]
                   - f_5 * pc_z[k] * sgf1_109[k];

        t_410[k] = f_2 * pgp0_123[k]
                   - f_3 * pgp1_123[k]
                   + f_4 * pc_x[k] * pgd_246[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, t_415, pc_x, pc_y, pfd_162, pgp0_125, \
                         pgp1_125, pgd_246, pgd_248, pgd_249, pgd_250, \
                         pgd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_6 * pfd_162[k]
                   + f_4 * pc_y[k] * pgd_246[k];

        t_412[k] = f_10 * pgp0_125[k]
                   - f_11 * pgp1_125[k]
                   + f_4 * pc_x[k] * pgd_248[k];

        t_413[k] = f_4 * pc_x[k] * pgd_249[k];

        t_414[k] = f_4 * pc_x[k] * pgd_250[k];

        t_415[k] = f_4 * pc_x[k] * pgd_251[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, pc_y, pfd_165, pfd_166, pfd_167, pgp0_124, \
                         pgp0_125, pgp1_124, pgp1_125, pgd_249, pgd_250, \
                         pgd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * pfd_165[k]
                   + f_2 * pgp0_124[k]
                   - f_3 * pgp1_124[k]
                   + f_4 * pc_y[k] * pgd_249[k];

        t_417[k] = f_6 * pfd_166[k]
                   + f_10 * pgp0_125[k]
                   - f_11 * pgp1_125[k]
                   + f_4 * pc_y[k] * pgd_250[k];

        t_418[k] = f_6 * pfd_167[k]
                   + f_4 * pc_y[k] * pgd_251[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, pb_y, pc_x, pc_y, pdf0_169, pdf1_169, pff0_279, \
                         pfd_168, pff1_279, pgp0_126, pgp1_126, \
                         pgd_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_12 * pdf0_169[k]
                   - f_13 * pdf1_169[k]
                   + pb_y[k] * pff0_279[k]
                   - f_5 * pc_y[k] * pff1_279[k];

        t_420[k] = f_2 * pgp0_126[k]
                   - f_3 * pgp1_126[k]
                   + f_4 * pc_x[k] * pgd_252[k];

        t_421[k] = f_9 * pfd_168[k]
                   + f_4 * pc_y[k] * pgd_252[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, pc_x, pc_y, pfd_171, pgp0_127, \
                         pgp0_128, pgp1_127, pgp1_128, pgd_254, pgd_255, pgd_256, \
                         pgd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_10 * pgp0_128[k]
                   - f_11 * pgp1_128[k]
                   + f_4 * pc_x[k] * pgd_254[k];

        t_423[k] = f_4 * pc_x[k] * pgd_255[k];

        t_424[k] = f_4 * pc_x[k] * pgd_256[k];

        t_425[k] = f_4 * pc_x[k] * pgd_257[k];

        t_426[k] = f_9 * pfd_171[k]
                   + f_2 * pgp0_127[k]
                   - f_3 * pgp1_127[k]
                   + f_4 * pc_y[k] * pgd_255[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pb_y, pc_y, pdf0_179, pdf1_179, pff0_289, \
                         pfd_172, pfd_173, pff1_289, pgp0_128, pgp1_128, pgd_256, \
                         pgd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_9 * pfd_172[k]
                   + f_10 * pgp0_128[k]
                   - f_11 * pgp1_128[k]
                   + f_4 * pc_y[k] * pgd_256[k];

        t_428[k] = f_9 * pfd_173[k]
                   + f_4 * pc_y[k] * pgd_257[k];

        t_429[k] = f_7 * pdf0_179[k]
                   - f_8 * pdf1_179[k]
                   + pb_y[k] * pff0_289[k]
                   - f_5 * pc_y[k] * pff1_289[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, pb_y, pc_x, pc_y, pff0_290, \
                         pff0_292, pfd_174, pff1_290, pff1_292, pgd_258, pgd_261, \
                         pgd_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = pb_y[k] * pff0_290[k]
                   - f_5 * pc_y[k] * pff1_290[k];

        t_431[k] = f_0 * pfd_174[k]
                   + f_4 * pc_y[k] * pgd_258[k];

        t_432[k] = pb_y[k] * pff0_292[k]
                   - f_5 * pc_y[k] * pff1_292[k];

        t_433[k] = f_4 * pc_x[k] * pgd_261[k];

        t_434[k] = f_4 * pc_x[k] * pgd_262[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pb_y, pc_x, pc_y, pff0_296, pff0_297, \
                         pfd_177, pfd_178, pfd_179, pff1_296, pff1_297, \
                         pgd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_4 * pc_x[k] * pgd_263[k];

        t_436[k] = pb_y[k] * pff0_296[k]
                   + f_6 * pfd_177[k]
                   - f_5 * pc_y[k] * pff1_296[k];

        t_437[k] = pb_y[k] * pff0_297[k]
                   + f_9 * pfd_178[k]
                   - f_5 * pc_y[k] * pff1_297[k];

        t_438[k] = f_0 * pfd_179[k]
                   + f_4 * pc_y[k] * pgd_263[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pb_y, pc_x, pc_y, pff0_299, pff1_299, \
                         pgp0_132, pgp0_134, pgp1_132, pgp1_134, pgd_264, \
                         pgd_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = pb_y[k] * pff0_299[k]
                   - f_5 * pc_y[k] * pff1_299[k];

        t_440[k] = f_2 * pgp0_132[k]
                   - f_3 * pgp1_132[k]
                   + f_4 * pc_x[k] * pgd_264[k];

        t_441[k] = f_4 * pc_y[k] * pgd_264[k];

        t_442[k] = f_10 * pgp0_134[k]
                   - f_11 * pgp1_134[k]
                   + f_4 * pc_x[k] * pgd_266[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, t_447, t_448, pc_x, pc_y, pgp0_133, \
                         pgp0_134, pgp1_133, pgp1_134, pgd_267, pgd_268, \
                         pgd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_4 * pc_x[k] * pgd_267[k];

        t_444[k] = f_4 * pc_x[k] * pgd_268[k];

        t_445[k] = f_4 * pc_x[k] * pgd_269[k];

        t_446[k] = f_2 * pgp0_133[k]
                   - f_3 * pgp1_133[k]
                   + f_4 * pc_y[k] * pgd_267[k];

        t_447[k] = f_10 * pgp0_134[k]
                   - f_11 * pgp1_134[k]
                   + f_4 * pc_y[k] * pgd_268[k];

        t_448[k] = f_4 * pc_y[k] * pgd_269[k];
    }

#pragma omp simd aligned(t_449, pc_z, sgd_89, pfd_179, pgp0_134, pgp1_134, \
                         pgd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_0 * sgd_89[k]
                   + f_1 * pfd_179[k]
                   + f_2 * pgp0_134[k]
                   - f_3 * pgp1_134[k]
                   + f_4 * pc_z[k] * pgd_269[k];
    }
}

auto
compute_prim_pgf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sgf0,
                                                   const size_t sgd, const size_t sgf1,
                                                   const size_t pdf0, const size_t pdf1,
                                                   const size_t pff0, const size_t pfd,
                                                   const size_t pff1, const size_t pgp0,
                                                   const size_t pgp1, const size_t pgd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pgf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sgf0,
                                                              sgd, sgf1, pdf0, pdf1, pff0, pfd,
                                                              pff1, pgp0, pgp1, pgd, ncols,
                                                              gamma, p, q);

    compute_prim_pgf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sgf0,
                                                              sgd, sgf1, pdf0, pdf1, pff0, pfd,
                                                              pff1, pgp0, pgp1, pgd, ncols,
                                                              gamma, p, q);

    compute_prim_pgf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sgf0,
                                                              sgd, sgf1, pdf0, pdf1, pff0, pfd,
                                                              pff1, pgp0, pgp1, pgd, ncols,
                                                              gamma, p, q);

    compute_prim_pgf_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, sgf0,
                                                              sgd, sgf1, pdf0, pdf1, pff0, pfd,
                                                              pff1, pgp0, pgp1, pgd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
