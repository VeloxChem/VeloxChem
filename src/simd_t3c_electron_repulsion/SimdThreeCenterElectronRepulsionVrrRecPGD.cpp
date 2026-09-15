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


#include "SimdThreeCenterElectronRepulsionVrrRecPGD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pgd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgd0, const size_t sgp,
                                                          const size_t sgd1, const size_t pdd0,
                                                          const size_t pdd1, const size_t pfd0,
                                                          const size_t pfp, const size_t pfd1,
                                                          const size_t pgs0, const size_t pgs1,
                                                          const size_t pgp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 0.5 / p;
    const auto f_7 = 0.5 * gamma / (p * q);
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / p;
    const auto f_11 = gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgd0_0 = buffer.data(sgd0 + 0);
    const auto *sgd0_3 = buffer.data(sgd0 + 3);
    const auto *sgd0_5 = buffer.data(sgd0 + 5);
    const auto *sgd0_12 = buffer.data(sgd0 + 12);
    const auto *sgd0_17 = buffer.data(sgd0 + 17);
    const auto *sgd0_30 = buffer.data(sgd0 + 30);
    const auto *sgd0_33 = buffer.data(sgd0 + 33);
    const auto *sgd0_60 = buffer.data(sgd0 + 60);
    const auto *sgd0_63 = buffer.data(sgd0 + 63);
    const auto *sgd0_65 = buffer.data(sgd0 + 65);
    const auto *sgd0_69 = buffer.data(sgd0 + 69);
    const auto *sgd0_71 = buffer.data(sgd0 + 71);
    const auto *sgd0_72 = buffer.data(sgd0 + 72);
    const auto *sgd0_75 = buffer.data(sgd0 + 75);
    const auto *sgd0_77 = buffer.data(sgd0 + 77);
    const auto *sgd0_81 = buffer.data(sgd0 + 81);
    const auto *sgd0_83 = buffer.data(sgd0 + 83);
    const auto *sgd0_84 = buffer.data(sgd0 + 84);
    const auto *sgd0_87 = buffer.data(sgd0 + 87);
    const auto *sgd0_89 = buffer.data(sgd0 + 89);

    const auto *sgp_0 = buffer.data(sgp + 0);
    const auto *sgp_1 = buffer.data(sgp + 1);
    const auto *sgp_16 = buffer.data(sgp + 16);
    const auto *sgp_18 = buffer.data(sgp + 18);
    const auto *sgp_27 = buffer.data(sgp + 27);
    const auto *sgp_30 = buffer.data(sgp + 30);
    const auto *sgp_36 = buffer.data(sgp + 36);
    const auto *sgp_42 = buffer.data(sgp + 42);

    const auto *sgd1_0 = buffer.data(sgd1 + 0);
    const auto *sgd1_3 = buffer.data(sgd1 + 3);
    const auto *sgd1_5 = buffer.data(sgd1 + 5);
    const auto *sgd1_12 = buffer.data(sgd1 + 12);
    const auto *sgd1_17 = buffer.data(sgd1 + 17);
    const auto *sgd1_30 = buffer.data(sgd1 + 30);
    const auto *sgd1_33 = buffer.data(sgd1 + 33);
    const auto *sgd1_60 = buffer.data(sgd1 + 60);
    const auto *sgd1_63 = buffer.data(sgd1 + 63);
    const auto *sgd1_65 = buffer.data(sgd1 + 65);
    const auto *sgd1_69 = buffer.data(sgd1 + 69);
    const auto *sgd1_71 = buffer.data(sgd1 + 71);
    const auto *sgd1_72 = buffer.data(sgd1 + 72);
    const auto *sgd1_75 = buffer.data(sgd1 + 75);
    const auto *sgd1_77 = buffer.data(sgd1 + 77);
    const auto *sgd1_81 = buffer.data(sgd1 + 81);
    const auto *sgd1_83 = buffer.data(sgd1 + 83);
    const auto *sgd1_84 = buffer.data(sgd1 + 84);
    const auto *sgd1_87 = buffer.data(sgd1 + 87);
    const auto *sgd1_89 = buffer.data(sgd1 + 89);

    const auto *pdd0_0 = buffer.data(pdd0 + 0);
    const auto *pdd0_45 = buffer.data(pdd0 + 45);
    const auto *pdd0_57 = buffer.data(pdd0 + 57);

    const auto *pdd1_0 = buffer.data(pdd1 + 0);
    const auto *pdd1_45 = buffer.data(pdd1 + 45);
    const auto *pdd1_57 = buffer.data(pdd1 + 57);

    const auto *pfd0_0 = buffer.data(pfd0 + 0);
    const auto *pfd0_3 = buffer.data(pfd0 + 3);
    const auto *pfd0_5 = buffer.data(pfd0 + 5);
    const auto *pfd0_6 = buffer.data(pfd0 + 6);
    const auto *pfd0_9 = buffer.data(pfd0 + 9);
    const auto *pfd0_12 = buffer.data(pfd0 + 12);
    const auto *pfd0_17 = buffer.data(pfd0 + 17);
    const auto *pfd0_18 = buffer.data(pfd0 + 18);
    const auto *pfd0_30 = buffer.data(pfd0 + 30);
    const auto *pfd0_36 = buffer.data(pfd0 + 36);
    const auto *pfd0_54 = buffer.data(pfd0 + 54);
    const auto *pfd0_63 = buffer.data(pfd0 + 63);
    const auto *pfd0_69 = buffer.data(pfd0 + 69);
    const auto *pfd0_81 = buffer.data(pfd0 + 81);

    const auto *pfp_0 = buffer.data(pfp + 0);
    const auto *pfp_1 = buffer.data(pfp + 1);
    const auto *pfp_2 = buffer.data(pfp + 2);
    const auto *pfp_3 = buffer.data(pfp + 3);
    const auto *pfp_4 = buffer.data(pfp + 4);
    const auto *pfp_5 = buffer.data(pfp + 5);
    const auto *pfp_6 = buffer.data(pfp + 6);
    const auto *pfp_8 = buffer.data(pfp + 8);
    const auto *pfp_9 = buffer.data(pfp + 9);
    const auto *pfp_10 = buffer.data(pfp + 10);
    const auto *pfp_11 = buffer.data(pfp + 11);
    const auto *pfp_12 = buffer.data(pfp + 12);
    const auto *pfp_13 = buffer.data(pfp + 13);
    const auto *pfp_14 = buffer.data(pfp + 14);
    const auto *pfp_15 = buffer.data(pfp + 15);
    const auto *pfp_16 = buffer.data(pfp + 16);
    const auto *pfp_17 = buffer.data(pfp + 17);
    const auto *pfp_18 = buffer.data(pfp + 18);
    const auto *pfp_20 = buffer.data(pfp + 20);
    const auto *pfp_21 = buffer.data(pfp + 21);
    const auto *pfp_23 = buffer.data(pfp + 23);
    const auto *pfp_24 = buffer.data(pfp + 24);
    const auto *pfp_26 = buffer.data(pfp + 26);
    const auto *pfp_27 = buffer.data(pfp + 27);
    const auto *pfp_29 = buffer.data(pfp + 29);
    const auto *pfp_31 = buffer.data(pfp + 31);
    const auto *pfp_32 = buffer.data(pfp + 32);
    const auto *pfp_33 = buffer.data(pfp + 33);
    const auto *pfp_34 = buffer.data(pfp + 34);
    const auto *pfp_35 = buffer.data(pfp + 35);
    const auto *pfp_37 = buffer.data(pfp + 37);
    const auto *pfp_38 = buffer.data(pfp + 38);
    const auto *pfp_39 = buffer.data(pfp + 39);
    const auto *pfp_40 = buffer.data(pfp + 40);
    const auto *pfp_41 = buffer.data(pfp + 41);
    const auto *pfp_42 = buffer.data(pfp + 42);
    const auto *pfp_43 = buffer.data(pfp + 43);
    const auto *pfp_44 = buffer.data(pfp + 44);
    const auto *pfp_46 = buffer.data(pfp + 46);
    const auto *pfp_47 = buffer.data(pfp + 47);

    const auto *pfd1_0 = buffer.data(pfd1 + 0);
    const auto *pfd1_3 = buffer.data(pfd1 + 3);
    const auto *pfd1_5 = buffer.data(pfd1 + 5);
    const auto *pfd1_6 = buffer.data(pfd1 + 6);
    const auto *pfd1_9 = buffer.data(pfd1 + 9);
    const auto *pfd1_12 = buffer.data(pfd1 + 12);
    const auto *pfd1_17 = buffer.data(pfd1 + 17);
    const auto *pfd1_18 = buffer.data(pfd1 + 18);
    const auto *pfd1_30 = buffer.data(pfd1 + 30);
    const auto *pfd1_36 = buffer.data(pfd1 + 36);
    const auto *pfd1_54 = buffer.data(pfd1 + 54);
    const auto *pfd1_63 = buffer.data(pfd1 + 63);
    const auto *pfd1_69 = buffer.data(pfd1 + 69);
    const auto *pfd1_81 = buffer.data(pfd1 + 81);

    const auto *pgs0_0 = buffer.data(pgs0 + 0);
    const auto *pgs0_1 = buffer.data(pgs0 + 1);
    const auto *pgs0_2 = buffer.data(pgs0 + 2);
    const auto *pgs0_3 = buffer.data(pgs0 + 3);
    const auto *pgs0_5 = buffer.data(pgs0 + 5);
    const auto *pgs0_6 = buffer.data(pgs0 + 6);
    const auto *pgs0_7 = buffer.data(pgs0 + 7);
    const auto *pgs0_8 = buffer.data(pgs0 + 8);
    const auto *pgs0_9 = buffer.data(pgs0 + 9);
    const auto *pgs0_16 = buffer.data(pgs0 + 16);
    const auto *pgs0_18 = buffer.data(pgs0 + 18);
    const auto *pgs0_19 = buffer.data(pgs0 + 19);

    const auto *pgs1_0 = buffer.data(pgs1 + 0);
    const auto *pgs1_1 = buffer.data(pgs1 + 1);
    const auto *pgs1_2 = buffer.data(pgs1 + 2);
    const auto *pgs1_3 = buffer.data(pgs1 + 3);
    const auto *pgs1_5 = buffer.data(pgs1 + 5);
    const auto *pgs1_6 = buffer.data(pgs1 + 6);
    const auto *pgs1_7 = buffer.data(pgs1 + 7);
    const auto *pgs1_8 = buffer.data(pgs1 + 8);
    const auto *pgs1_9 = buffer.data(pgs1 + 9);
    const auto *pgs1_16 = buffer.data(pgs1 + 16);
    const auto *pgs1_18 = buffer.data(pgs1 + 18);
    const auto *pgs1_19 = buffer.data(pgs1 + 19);

    const auto *pgp_0 = buffer.data(pgp + 0);
    const auto *pgp_1 = buffer.data(pgp + 1);
    const auto *pgp_2 = buffer.data(pgp + 2);
    const auto *pgp_3 = buffer.data(pgp + 3);
    const auto *pgp_4 = buffer.data(pgp + 4);
    const auto *pgp_5 = buffer.data(pgp + 5);
    const auto *pgp_6 = buffer.data(pgp + 6);
    const auto *pgp_8 = buffer.data(pgp + 8);
    const auto *pgp_9 = buffer.data(pgp + 9);
    const auto *pgp_10 = buffer.data(pgp + 10);
    const auto *pgp_11 = buffer.data(pgp + 11);
    const auto *pgp_12 = buffer.data(pgp + 12);
    const auto *pgp_14 = buffer.data(pgp + 14);
    const auto *pgp_15 = buffer.data(pgp + 15);
    const auto *pgp_16 = buffer.data(pgp + 16);
    const auto *pgp_17 = buffer.data(pgp + 17);
    const auto *pgp_18 = buffer.data(pgp + 18);
    const auto *pgp_19 = buffer.data(pgp + 19);
    const auto *pgp_20 = buffer.data(pgp + 20);
    const auto *pgp_21 = buffer.data(pgp + 21);
    const auto *pgp_22 = buffer.data(pgp + 22);
    const auto *pgp_23 = buffer.data(pgp + 23);
    const auto *pgp_24 = buffer.data(pgp + 24);
    const auto *pgp_25 = buffer.data(pgp + 25);
    const auto *pgp_26 = buffer.data(pgp + 26);
    const auto *pgp_27 = buffer.data(pgp + 27);
    const auto *pgp_28 = buffer.data(pgp + 28);
    const auto *pgp_29 = buffer.data(pgp + 29);
    const auto *pgp_30 = buffer.data(pgp + 30);
    const auto *pgp_32 = buffer.data(pgp + 32);
    const auto *pgp_33 = buffer.data(pgp + 33);
    const auto *pgp_35 = buffer.data(pgp + 35);
    const auto *pgp_36 = buffer.data(pgp + 36);
    const auto *pgp_38 = buffer.data(pgp + 38);
    const auto *pgp_39 = buffer.data(pgp + 39);
    const auto *pgp_41 = buffer.data(pgp + 41);
    const auto *pgp_42 = buffer.data(pgp + 42);
    const auto *pgp_44 = buffer.data(pgp + 44);
    const auto *pgp_46 = buffer.data(pgp + 46);
    const auto *pgp_47 = buffer.data(pgp + 47);
    const auto *pgp_48 = buffer.data(pgp + 48);
    const auto *pgp_49 = buffer.data(pgp + 49);
    const auto *pgp_50 = buffer.data(pgp + 50);
    const auto *pgp_52 = buffer.data(pgp + 52);
    const auto *pgp_53 = buffer.data(pgp + 53);
    const auto *pgp_54 = buffer.data(pgp + 54);
    const auto *pgp_55 = buffer.data(pgp + 55);
    const auto *pgp_56 = buffer.data(pgp + 56);
    const auto *pgp_57 = buffer.data(pgp + 57);
    const auto *pgp_58 = buffer.data(pgp + 58);
    const auto *pgp_59 = buffer.data(pgp + 59);
    const auto *pgp_61 = buffer.data(pgp + 61);
    const auto *pgp_62 = buffer.data(pgp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sgp_0, pfp_0, pgs0_0, \
                         pgs1_0, pgp_0, pgp_1, pgp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgp_0[k]
                 + f_1 * pfp_0[k]
                 + f_2 * pgs0_0[k]
                 - f_3 * pgs1_0[k]
                 + f_4 * pc_x[k] * pgp_0[k];

        t_1[k] = f_4 * pc_y[k] * pgp_0[k];

        t_2[k] = f_4 * pc_z[k] * pgp_0[k];

        t_3[k] = f_2 * pgs0_0[k]
                 - f_3 * pgs1_0[k]
                 + f_4 * pc_y[k] * pgp_1[k];

        t_4[k] = f_4 * pc_y[k] * pgp_2[k];

        t_5[k] = f_2 * pgs0_0[k]
                 - f_3 * pgs1_0[k]
                 + f_4 * pc_z[k] * pgp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pc_y, pc_z, pfd0_0, pfp_0, pfp_1, pfd1_0, \
                         pgs0_1, pgs1_1, pgp_3, pgp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * pfd0_0[k]
                 - f_5 * pc_y[k] * pfd1_0[k];

        t_7[k] = f_0 * pfp_0[k]
                 + f_4 * pc_y[k] * pgp_3[k];

        t_8[k] = f_4 * pc_z[k] * pgp_3[k];

        t_9[k] = f_0 * pfp_1[k]
                 + f_2 * pgs0_1[k]
                 - f_3 * pgs1_1[k]
                 + f_4 * pc_y[k] * pgp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pb_z, pc_y, pc_z, pfd0_0, pfd0_5, \
                         pfp_2, pfd1_0, pfd1_5, pgp_5, pgp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * pfp_2[k]
                  + f_4 * pc_y[k] * pgp_5[k];

        t_11[k] = pb_y[k] * pfd0_5[k]
                  - f_5 * pc_y[k] * pfd1_5[k];

        t_12[k] = pb_z[k] * pfd0_0[k]
                  - f_5 * pc_z[k] * pfd1_0[k];

        t_13[k] = f_4 * pc_y[k] * pgp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_z, pc_y, pc_z, pfd0_3, pfp_0, pfp_2, \
                         pfd1_3, pgs0_2, pgs1_2, pgp_6, pgp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * pfp_0[k]
                  + f_4 * pc_z[k] * pgp_6[k];

        t_15[k] = pb_z[k] * pfd0_3[k]
                  - f_5 * pc_z[k] * pfd1_3[k];

        t_16[k] = f_4 * pc_y[k] * pgp_8[k];

        t_17[k] = f_0 * pfp_2[k]
                  + f_2 * pgs0_2[k]
                  - f_3 * pgs1_2[k]
                  + f_4 * pc_z[k] * pgp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, pc_y, pc_z, pdd0_0, pdd1_0, pfd0_6, pfp_3, \
                         pfd1_6, pgp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_6 * pdd0_0[k]
                  - f_7 * pdd1_0[k]
                  + pb_y[k] * pfd0_6[k]
                  - f_5 * pc_y[k] * pfd1_6[k];

        t_19[k] = f_8 * pfp_3[k]
                  + f_4 * pc_y[k] * pgp_9[k];

        t_20[k] = f_4 * pc_z[k] * pgp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pc_y, pc_z, pfd0_12, pfp_4, pfp_5, \
                         pfd1_12, pgs0_3, pgs1_3, pgp_10, pgp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * pfp_4[k]
                  + f_2 * pgs0_3[k]
                  - f_3 * pgs1_3[k]
                  + f_4 * pc_y[k] * pgp_10[k];

        t_22[k] = f_8 * pfp_5[k]
                  + f_4 * pc_y[k] * pgp_11[k];

        t_23[k] = f_2 * pgs0_3[k]
                  - f_3 * pgs1_3[k]
                  + f_4 * pc_z[k] * pgp_11[k];

        t_24[k] = pb_y[k] * pfd0_12[k]
                  - f_5 * pc_y[k] * pfd1_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_z, pc_y, pc_z, pfd0_9, pfp_3, pfp_6, \
                         pfp_8, pfd1_9, pgp_12, pgp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * pfp_6[k]
                  + f_4 * pc_y[k] * pgp_12[k];

        t_26[k] = f_0 * pfp_3[k]
                  + f_4 * pc_z[k] * pgp_12[k];

        t_27[k] = pb_z[k] * pfd0_9[k]
                  - f_5 * pc_z[k] * pfd1_9[k];

        t_28[k] = f_0 * pfp_8[k]
                  + f_4 * pc_y[k] * pgp_14[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, pdd0_0, pdd1_0, \
                         pfd0_12, pfd0_17, pfp_6, pfd1_12, pfd1_17, \
                         pgp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * pfd0_17[k]
                  - f_5 * pc_y[k] * pfd1_17[k];

        t_30[k] = f_6 * pdd0_0[k]
                  - f_7 * pdd1_0[k]
                  + pb_z[k] * pfd0_12[k]
                  - f_5 * pc_z[k] * pfd1_12[k];

        t_31[k] = f_4 * pc_y[k] * pgp_15[k];

        t_32[k] = f_8 * pfp_6[k]
                  + f_4 * pc_z[k] * pgp_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pc_y, pc_z, pfp_8, pgs0_5, pgs1_5, pgp_16, \
                         pgp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_2 * pgs0_5[k]
                  - f_3 * pgs1_5[k]
                  + f_4 * pc_y[k] * pgp_16[k];

        t_34[k] = f_4 * pc_y[k] * pgp_17[k];

        t_35[k] = f_8 * pfp_8[k]
                  + f_2 * pgs0_5[k]
                  - f_3 * pgs1_5[k]
                  + f_4 * pc_z[k] * pgp_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_x, pc_y, pc_z, sgp_18, pfp_9, pfp_10, \
                         pfp_18, pgs0_6, pgs1_6, pgp_18, pgp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * sgp_18[k]
                  + f_0 * pfp_18[k]
                  + f_2 * pgs0_6[k]
                  - f_3 * pgs1_6[k]
                  + f_4 * pc_x[k] * pgp_18[k];

        t_37[k] = f_9 * pfp_9[k]
                  + f_4 * pc_y[k] * pgp_18[k];

        t_38[k] = f_4 * pc_z[k] * pgp_18[k];

        t_39[k] = f_9 * pfp_10[k]
                  + f_2 * pgs0_6[k]
                  - f_3 * pgs1_6[k]
                  + f_4 * pc_y[k] * pgp_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_z, pc_y, pc_z, pfd0_18, pfp_11, pfp_12, \
                         pfd1_18, pgs0_6, pgs1_6, pgp_20, pgp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * pfp_11[k]
                  + f_4 * pc_y[k] * pgp_20[k];

        t_41[k] = f_2 * pgs0_6[k]
                  - f_3 * pgs1_6[k]
                  + f_4 * pc_z[k] * pgp_20[k];

        t_42[k] = pb_z[k] * pfd0_18[k]
                  - f_5 * pc_z[k] * pfd1_18[k];

        t_43[k] = f_8 * pfp_12[k]
                  + f_4 * pc_y[k] * pgp_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pc_y, pc_z, pfp_9, pfp_11, pfp_13, pfp_14, \
                         pgs0_7, pgs1_7, pgp_21, pgp_22, pgp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * pfp_9[k]
                  + f_4 * pc_z[k] * pgp_21[k];

        t_45[k] = f_8 * pfp_13[k]
                  + f_2 * pgs0_7[k]
                  - f_3 * pgs1_7[k]
                  + f_4 * pc_y[k] * pgp_22[k];

        t_46[k] = f_8 * pfp_14[k]
                  + f_4 * pc_y[k] * pgp_23[k];

        t_47[k] = f_0 * pfp_11[k]
                  + f_2 * pgs0_7[k]
                  - f_3 * pgs1_7[k]
                  + f_4 * pc_z[k] * pgp_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_y, pc_y, pc_z, pfd0_30, pfp_12, pfp_15, \
                         pfp_16, pfd1_30, pgs0_8, pgs1_8, pgp_24, \
                         pgp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * pfd0_30[k]
                  - f_5 * pc_y[k] * pfd1_30[k];

        t_49[k] = f_0 * pfp_15[k]
                  + f_4 * pc_y[k] * pgp_24[k];

        t_50[k] = f_8 * pfp_12[k]
                  + f_4 * pc_z[k] * pgp_24[k];

        t_51[k] = f_0 * pfp_16[k]
                  + f_2 * pgs0_8[k]
                  - f_3 * pgs1_8[k]
                  + f_4 * pc_y[k] * pgp_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pc_x, pc_y, pc_z, sgp_27, pfp_14, pfp_17, pfp_27, \
                         pgs0_8, pgs0_9, pgs1_8, pgs1_9, pgp_26, \
                         pgp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * pfp_17[k]
                  + f_4 * pc_y[k] * pgp_26[k];

        t_53[k] = f_8 * pfp_14[k]
                  + f_2 * pgs0_8[k]
                  - f_3 * pgs1_8[k]
                  + f_4 * pc_z[k] * pgp_26[k];

        t_54[k] = f_0 * sgp_27[k]
                  + f_0 * pfp_27[k]
                  + f_2 * pgs0_9[k]
                  - f_3 * pgs1_9[k]
                  + f_4 * pc_x[k] * pgp_27[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pc_y, pc_z, pfp_15, pfp_17, pgs0_9, \
                         pgs1_9, pgp_27, pgp_28, pgp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_4 * pc_y[k] * pgp_27[k];

        t_56[k] = f_9 * pfp_15[k]
                  + f_4 * pc_z[k] * pgp_27[k];

        t_57[k] = f_2 * pgs0_9[k]
                  - f_3 * pgs1_9[k]
                  + f_4 * pc_y[k] * pgp_28[k];

        t_58[k] = f_4 * pc_y[k] * pgp_29[k];

        t_59[k] = f_9 * pfp_17[k]
                  + f_2 * pgs0_9[k]
                  - f_3 * pgs1_9[k]
                  + f_4 * pc_z[k] * pgp_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_x, pc_x, pc_y, pc_z, sgd0_60, sgd0_63, \
                         sgp_30, sgd1_60, sgd1_63, pfp_18, pgp_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_x[k] * sgd0_60[k]
                  + f_8 * sgp_30[k]
                  - f_5 * pc_x[k] * sgd1_60[k];

        t_61[k] = f_1 * pfp_18[k]
                  + f_4 * pc_y[k] * pgp_30[k];

        t_62[k] = f_4 * pc_z[k] * pgp_30[k];

        t_63[k] = pa_x[k] * sgd0_63[k]
                  - f_5 * pc_x[k] * sgd1_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_x, pb_z, pc_x, pc_y, pc_z, sgd0_65, sgd1_65, \
                         pfd0_36, pfp_20, pfd1_36, pgp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * pfp_20[k]
                  + f_4 * pc_y[k] * pgp_32[k];

        t_65[k] = pa_x[k] * sgd0_65[k]
                  - f_5 * pc_x[k] * sgd1_65[k];

        t_66[k] = pb_z[k] * pfd0_36[k]
                  - f_5 * pc_z[k] * pfd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_x, pc_x, pc_y, pc_z, sgd0_69, sgd1_69, \
                         pfp_18, pfp_21, pfp_23, pgp_33, pgp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_9 * pfp_21[k]
                  + f_4 * pc_y[k] * pgp_33[k];

        t_68[k] = f_0 * pfp_18[k]
                  + f_4 * pc_z[k] * pgp_33[k];

        t_69[k] = pa_x[k] * sgd0_69[k]
                  - f_5 * pc_x[k] * sgd1_69[k];

        t_70[k] = f_9 * pfp_23[k]
                  + f_4 * pc_y[k] * pgp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pc_x, pc_y, pc_z, sgd0_71, sgd0_72, \
                         sgp_36, sgd1_71, sgd1_72, pfp_21, pfp_24, \
                         pgp_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_x[k] * sgd0_71[k]
                  - f_5 * pc_x[k] * sgd1_71[k];

        t_72[k] = pa_x[k] * sgd0_72[k]
                  + f_8 * sgp_36[k]
                  - f_5 * pc_x[k] * sgd1_72[k];

        t_73[k] = f_8 * pfp_24[k]
                  + f_4 * pc_y[k] * pgp_36[k];

        t_74[k] = f_8 * pfp_21[k]
                  + f_4 * pc_z[k] * pgp_36[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_x, pb_y, pc_x, pc_y, sgd0_75, sgd0_77, \
                         sgd1_75, sgd1_77, pfd0_54, pfp_26, pfd1_54, \
                         pgp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pa_x[k] * sgd0_75[k]
                  - f_5 * pc_x[k] * sgd1_75[k];

        t_76[k] = f_8 * pfp_26[k]
                  + f_4 * pc_y[k] * pgp_38[k];

        t_77[k] = pa_x[k] * sgd0_77[k]
                  - f_5 * pc_x[k] * sgd1_77[k];

        t_78[k] = pb_y[k] * pfd0_54[k]
                  - f_5 * pc_y[k] * pfd1_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_x, pc_x, pc_y, pc_z, sgd0_81, sgd1_81, \
                         pfp_24, pfp_27, pfp_29, pgp_39, pgp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_0 * pfp_27[k]
                  + f_4 * pc_y[k] * pgp_39[k];

        t_80[k] = f_9 * pfp_24[k]
                  + f_4 * pc_z[k] * pgp_39[k];

        t_81[k] = pa_x[k] * sgd0_81[k]
                  - f_5 * pc_x[k] * sgd1_81[k];

        t_82[k] = f_0 * pfp_29[k]
                  + f_4 * pc_y[k] * pgp_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pc_x, pc_y, pc_z, sgd0_83, sgd0_84, \
                         sgp_42, sgd1_83, sgd1_84, pfp_27, pgp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pa_x[k] * sgd0_83[k]
                  - f_5 * pc_x[k] * sgd1_83[k];

        t_84[k] = pa_x[k] * sgd0_84[k]
                  + f_8 * sgp_42[k]
                  - f_5 * pc_x[k] * sgd1_84[k];

        t_85[k] = f_4 * pc_y[k] * pgp_42[k];

        t_86[k] = f_1 * pfp_27[k]
                  + f_4 * pc_z[k] * pgp_42[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pa_y, pc_x, pc_y, sgd0_0, sgd0_87, \
                         sgd0_89, sgd1_0, sgd1_87, sgd1_89, pgp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_x[k] * sgd0_87[k]
                  - f_5 * pc_x[k] * sgd1_87[k];

        t_88[k] = f_4 * pc_y[k] * pgp_44[k];

        t_89[k] = pa_x[k] * sgd0_89[k]
                  - f_5 * pc_x[k] * sgd1_89[k];

        t_90[k] = pa_y[k] * sgd0_0[k]
                  - f_5 * pc_y[k] * sgd1_0[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, sgd0_3, sgp_1, \
                         sgd1_3, pfp_31, pfp_32, pgp_46, pgp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_1 * pfp_31[k]
                  + f_4 * pc_x[k] * pgp_46[k];

        t_92[k] = f_1 * pfp_32[k]
                  + f_4 * pc_x[k] * pgp_47[k];

        t_93[k] = pa_y[k] * sgd0_3[k]
                  + f_8 * sgp_1[k]
                  - f_5 * pc_y[k] * sgd1_3[k];

        t_94[k] = f_4 * pc_z[k] * pgp_46[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pa_y, pc_x, pc_y, sgd0_5, sgd1_5, pfp_33, pfp_34, \
                         pgs0_16, pgs1_16, pgp_48, pgp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pa_y[k] * sgd0_5[k]
                  - f_5 * pc_y[k] * sgd1_5[k];

        t_96[k] = f_9 * pfp_33[k]
                  + f_2 * pgs0_16[k]
                  - f_3 * pgs1_16[k]
                  + f_4 * pc_x[k] * pgp_48[k];

        t_97[k] = f_9 * pfp_34[k]
                  + f_4 * pc_x[k] * pgp_49[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_x, pc_x, pc_z, pdd0_45, pdd1_45, \
                         pfd0_69, pfp_35, pfd1_69, pgs0_16, pgs1_16, pgp_49, \
                         pgp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_9 * pfp_35[k]
                  + f_4 * pc_x[k] * pgp_50[k];

        t_99[k] = f_10 * pdd0_45[k]
                  - f_11 * pdd1_45[k]
                  + pb_x[k] * pfd0_69[k]
                  - f_5 * pc_x[k] * pfd1_69[k];

        t_100[k] = f_4 * pc_z[k] * pgp_49[k];

        t_101[k] = f_2 * pgs0_16[k]
                   - f_3 * pgs1_16[k]
                   + f_4 * pc_z[k] * pgp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_x, pc_y, sgd0_12, sgd1_12, pfp_37, \
                         pfp_38, pgp_52, pgp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = pa_y[k] * sgd0_12[k]
                   - f_5 * pc_y[k] * sgd1_12[k];

        t_103[k] = f_9 * pfp_37[k]
                   + f_4 * pc_x[k] * pgp_52[k];

        t_104[k] = f_9 * pfp_38[k]
                   + f_4 * pc_x[k] * pgp_53[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_y, pb_z, pc_y, pc_z, sgd0_17, sgd1_17, \
                         pfd0_63, pfp_31, pfd1_63, pgp_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pb_z[k] * pfd0_63[k]
                   - f_5 * pc_z[k] * pfd1_63[k];

        t_106[k] = f_0 * pfp_31[k]
                   + f_4 * pc_z[k] * pgp_52[k];

        t_107[k] = pa_y[k] * sgd0_17[k]
                   - f_5 * pc_y[k] * sgd1_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pfp_39, pfp_40, pfp_41, pgs0_18, pgs1_18, \
                         pgp_54, pgp_55, pgp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_8 * pfp_39[k]
                   + f_2 * pgs0_18[k]
                   - f_3 * pgs1_18[k]
                   + f_4 * pc_x[k] * pgp_54[k];

        t_109[k] = f_8 * pfp_40[k]
                   + f_4 * pc_x[k] * pgp_55[k];

        t_110[k] = f_8 * pfp_41[k]
                   + f_4 * pc_x[k] * pgp_56[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_x, pc_x, pc_z, pdd0_57, pdd1_57, pfd0_81, \
                         pfd1_81, pgs0_18, pgs1_18, pgp_55, pgp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_6 * pdd0_57[k]
                   - f_7 * pdd1_57[k]
                   + pb_x[k] * pfd0_81[k]
                   - f_5 * pc_x[k] * pfd1_81[k];

        t_112[k] = f_4 * pc_z[k] * pgp_55[k];

        t_113[k] = f_2 * pgs0_18[k]
                   - f_3 * pgs1_18[k]
                   + f_4 * pc_z[k] * pgp_56[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pc_x, pfp_42, pfp_43, pfp_44, pgs0_19, pgs1_19, \
                         pgp_57, pgp_58, pgp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_8 * pfp_42[k]
                   + f_2 * pgs0_19[k]
                   - f_3 * pgs1_19[k]
                   + f_4 * pc_x[k] * pgp_57[k];

        t_115[k] = f_8 * pfp_43[k]
                   + f_4 * pc_x[k] * pgp_58[k];

        t_116[k] = f_8 * pfp_44[k]
                   + f_4 * pc_x[k] * pgp_59[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pb_z, pc_z, pfd0_69, pfp_34, pfp_35, pfd1_69, \
                         pgs0_19, pgs1_19, pgp_58, pgp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pb_z[k] * pfd0_69[k]
                   - f_5 * pc_z[k] * pfd1_69[k];

        t_118[k] = f_0 * pfp_34[k]
                   + f_4 * pc_z[k] * pgp_58[k];

        t_119[k] = f_0 * pfp_35[k]
                   + f_2 * pgs0_19[k]
                   - f_3 * pgs1_19[k]
                   + f_4 * pc_z[k] * pgp_59[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, pc_x, pc_y, sgd0_30, sgd0_33, \
                         sgp_16, sgd1_30, sgd1_33, pfp_46, pfp_47, pgp_61, \
                         pgp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_y[k] * sgd0_30[k]
                   - f_5 * pc_y[k] * sgd1_30[k];

        t_121[k] = f_8 * pfp_46[k]
                   + f_4 * pc_x[k] * pgp_61[k];

        t_122[k] = f_8 * pfp_47[k]
                   + f_4 * pc_x[k] * pgp_62[k];

        t_123[k] = pa_y[k] * sgd0_33[k]
                   + f_8 * sgp_16[k]
                   - f_5 * pc_y[k] * sgd1_33[k];
    }
}

static auto
compute_prim_pgd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgd0, const size_t sgp,
                                                          const size_t sgd1, const size_t pdd0,
                                                          const size_t pdd1, const size_t pfd0,
                                                          const size_t pfp, const size_t pfd1,
                                                          const size_t pgs0, const size_t pgs1,
                                                          const size_t pgp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 0.5 / p;
    const auto f_7 = 0.5 * gamma / (p * q);
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / p;
    const auto f_11 = gamma / (p * q);

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
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgd0_0 = buffer.data(sgd0 + 0);
    const auto *sgd0_3 = buffer.data(sgd0 + 3);
    const auto *sgd0_5 = buffer.data(sgd0 + 5);
    const auto *sgd0_6 = buffer.data(sgd0 + 6);
    const auto *sgd0_9 = buffer.data(sgd0 + 9);
    const auto *sgd0_18 = buffer.data(sgd0 + 18);
    const auto *sgd0_21 = buffer.data(sgd0 + 21);
    const auto *sgd0_23 = buffer.data(sgd0 + 23);
    const auto *sgd0_35 = buffer.data(sgd0 + 35);
    const auto *sgd0_36 = buffer.data(sgd0 + 36);
    const auto *sgd0_39 = buffer.data(sgd0 + 39);
    const auto *sgd0_54 = buffer.data(sgd0 + 54);
    const auto *sgd0_59 = buffer.data(sgd0 + 59);
    const auto *sgd0_60 = buffer.data(sgd0 + 60);
    const auto *sgd0_63 = buffer.data(sgd0 + 63);
    const auto *sgd0_65 = buffer.data(sgd0 + 65);
    const auto *sgd0_84 = buffer.data(sgd0 + 84);
    const auto *sgd0_87 = buffer.data(sgd0 + 87);
    const auto *sgd0_89 = buffer.data(sgd0 + 89);

    const auto *sgp_2 = buffer.data(sgp + 2);
    const auto *sgp_11 = buffer.data(sgp + 11);
    const auto *sgp_31 = buffer.data(sgp + 31);
    const auto *sgp_32 = buffer.data(sgp + 32);
    const auto *sgp_40 = buffer.data(sgp + 40);
    const auto *sgp_43 = buffer.data(sgp + 43);

    const auto *sgd1_0 = buffer.data(sgd1 + 0);
    const auto *sgd1_3 = buffer.data(sgd1 + 3);
    const auto *sgd1_5 = buffer.data(sgd1 + 5);
    const auto *sgd1_6 = buffer.data(sgd1 + 6);
    const auto *sgd1_9 = buffer.data(sgd1 + 9);
    const auto *sgd1_18 = buffer.data(sgd1 + 18);
    const auto *sgd1_21 = buffer.data(sgd1 + 21);
    const auto *sgd1_23 = buffer.data(sgd1 + 23);
    const auto *sgd1_35 = buffer.data(sgd1 + 35);
    const auto *sgd1_36 = buffer.data(sgd1 + 36);
    const auto *sgd1_39 = buffer.data(sgd1 + 39);
    const auto *sgd1_54 = buffer.data(sgd1 + 54);
    const auto *sgd1_59 = buffer.data(sgd1 + 59);
    const auto *sgd1_60 = buffer.data(sgd1 + 60);
    const auto *sgd1_63 = buffer.data(sgd1 + 63);
    const auto *sgd1_65 = buffer.data(sgd1 + 65);
    const auto *sgd1_84 = buffer.data(sgd1 + 84);
    const auto *sgd1_87 = buffer.data(sgd1 + 87);
    const auto *sgd1_89 = buffer.data(sgd1 + 89);

    const auto *pdd0_57 = buffer.data(pdd0 + 57);
    const auto *pdd0_89 = buffer.data(pdd0 + 89);
    const auto *pdd0_107 = buffer.data(pdd0 + 107);

    const auto *pdd1_57 = buffer.data(pdd1 + 57);
    const auto *pdd1_89 = buffer.data(pdd1 + 89);
    const auto *pdd1_107 = buffer.data(pdd1 + 107);

    const auto *pfd0_78 = buffer.data(pfd0 + 78);
    const auto *pfd0_96 = buffer.data(pfd0 + 96);
    const auto *pfd0_99 = buffer.data(pfd0 + 99);
    const auto *pfd0_101 = buffer.data(pfd0 + 101);
    const auto *pfd0_105 = buffer.data(pfd0 + 105);
    const auto *pfd0_107 = buffer.data(pfd0 + 107);
    const auto *pfd0_111 = buffer.data(pfd0 + 111);
    const auto *pfd0_113 = buffer.data(pfd0 + 113);
    const auto *pfd0_117 = buffer.data(pfd0 + 117);
    const auto *pfd0_125 = buffer.data(pfd0 + 125);
    const auto *pfd0_132 = buffer.data(pfd0 + 132);
    const auto *pfd0_137 = buffer.data(pfd0 + 137);
    const auto *pfd0_150 = buffer.data(pfd0 + 150);
    const auto *pfd0_155 = buffer.data(pfd0 + 155);
    const auto *pfd0_161 = buffer.data(pfd0 + 161);
    const auto *pfd0_165 = buffer.data(pfd0 + 165);
    const auto *pfd0_167 = buffer.data(pfd0 + 167);
    const auto *pfd0_171 = buffer.data(pfd0 + 171);
    const auto *pfd0_173 = buffer.data(pfd0 + 173);
    const auto *pfd0_174 = buffer.data(pfd0 + 174);
    const auto *pfd0_177 = buffer.data(pfd0 + 177);
    const auto *pfd0_179 = buffer.data(pfd0 + 179);

    const auto *pfp_37 = buffer.data(pfp + 37);
    const auto *pfp_40 = buffer.data(pfp + 40);
    const auto *pfp_43 = buffer.data(pfp + 43);
    const auto *pfp_46 = buffer.data(pfp + 46);
    const auto *pfp_48 = buffer.data(pfp + 48);
    const auto *pfp_49 = buffer.data(pfp + 49);
    const auto *pfp_50 = buffer.data(pfp + 50);
    const auto *pfp_52 = buffer.data(pfp + 52);
    const auto *pfp_53 = buffer.data(pfp + 53);
    const auto *pfp_54 = buffer.data(pfp + 54);
    const auto *pfp_55 = buffer.data(pfp + 55);
    const auto *pfp_56 = buffer.data(pfp + 56);
    const auto *pfp_58 = buffer.data(pfp + 58);
    const auto *pfp_59 = buffer.data(pfp + 59);
    const auto *pfp_61 = buffer.data(pfp + 61);
    const auto *pfp_62 = buffer.data(pfp + 62);
    const auto *pfp_64 = buffer.data(pfp + 64);
    const auto *pfp_65 = buffer.data(pfp + 65);
    const auto *pfp_66 = buffer.data(pfp + 66);
    const auto *pfp_67 = buffer.data(pfp + 67);
    const auto *pfp_68 = buffer.data(pfp + 68);
    const auto *pfp_70 = buffer.data(pfp + 70);
    const auto *pfp_71 = buffer.data(pfp + 71);
    const auto *pfp_73 = buffer.data(pfp + 73);
    const auto *pfp_74 = buffer.data(pfp + 74);
    const auto *pfp_75 = buffer.data(pfp + 75);
    const auto *pfp_76 = buffer.data(pfp + 76);
    const auto *pfp_77 = buffer.data(pfp + 77);
    const auto *pfp_79 = buffer.data(pfp + 79);
    const auto *pfp_80 = buffer.data(pfp + 80);
    const auto *pfp_81 = buffer.data(pfp + 81);
    const auto *pfp_82 = buffer.data(pfp + 82);
    const auto *pfp_83 = buffer.data(pfp + 83);
    const auto *pfp_85 = buffer.data(pfp + 85);
    const auto *pfp_86 = buffer.data(pfp + 86);
    const auto *pfp_87 = buffer.data(pfp + 87);
    const auto *pfp_88 = buffer.data(pfp + 88);
    const auto *pfp_89 = buffer.data(pfp + 89);

    const auto *pfd1_78 = buffer.data(pfd1 + 78);
    const auto *pfd1_96 = buffer.data(pfd1 + 96);
    const auto *pfd1_99 = buffer.data(pfd1 + 99);
    const auto *pfd1_101 = buffer.data(pfd1 + 101);
    const auto *pfd1_105 = buffer.data(pfd1 + 105);
    const auto *pfd1_107 = buffer.data(pfd1 + 107);
    const auto *pfd1_111 = buffer.data(pfd1 + 111);
    const auto *pfd1_113 = buffer.data(pfd1 + 113);
    const auto *pfd1_117 = buffer.data(pfd1 + 117);
    const auto *pfd1_125 = buffer.data(pfd1 + 125);
    const auto *pfd1_132 = buffer.data(pfd1 + 132);
    const auto *pfd1_137 = buffer.data(pfd1 + 137);
    const auto *pfd1_150 = buffer.data(pfd1 + 150);
    const auto *pfd1_155 = buffer.data(pfd1 + 155);
    const auto *pfd1_161 = buffer.data(pfd1 + 161);
    const auto *pfd1_165 = buffer.data(pfd1 + 165);
    const auto *pfd1_167 = buffer.data(pfd1 + 167);
    const auto *pfd1_171 = buffer.data(pfd1 + 171);
    const auto *pfd1_173 = buffer.data(pfd1 + 173);
    const auto *pfd1_174 = buffer.data(pfd1 + 174);
    const auto *pfd1_177 = buffer.data(pfd1 + 177);
    const auto *pfd1_179 = buffer.data(pfd1 + 179);

    const auto *pgs0_23 = buffer.data(pgs0 + 23);
    const auto *pgs0_25 = buffer.data(pgs0 + 25);
    const auto *pgs0_26 = buffer.data(pgs0 + 26);
    const auto *pgs0_27 = buffer.data(pgs0 + 27);
    const auto *pgs0_28 = buffer.data(pgs0 + 28);
    const auto *pgs0_32 = buffer.data(pgs0 + 32);
    const auto *pgs0_34 = buffer.data(pgs0 + 34);
    const auto *pgs0_35 = buffer.data(pgs0 + 35);
    const auto *pgs0_37 = buffer.data(pgs0 + 37);
    const auto *pgs0_41 = buffer.data(pgs0 + 41);

    const auto *pgs1_23 = buffer.data(pgs1 + 23);
    const auto *pgs1_25 = buffer.data(pgs1 + 25);
    const auto *pgs1_26 = buffer.data(pgs1 + 26);
    const auto *pgs1_27 = buffer.data(pgs1 + 27);
    const auto *pgs1_28 = buffer.data(pgs1 + 28);
    const auto *pgs1_32 = buffer.data(pgs1 + 32);
    const auto *pgs1_34 = buffer.data(pgs1 + 34);
    const auto *pgs1_35 = buffer.data(pgs1 + 35);
    const auto *pgs1_37 = buffer.data(pgs1 + 37);
    const auto *pgs1_41 = buffer.data(pgs1 + 41);

    const auto *pgp_61 = buffer.data(pgp + 61);
    const auto *pgp_64 = buffer.data(pgp + 64);
    const auto *pgp_65 = buffer.data(pgp + 65);
    const auto *pgp_67 = buffer.data(pgp + 67);
    const auto *pgp_68 = buffer.data(pgp + 68);
    const auto *pgp_69 = buffer.data(pgp + 69);
    const auto *pgp_70 = buffer.data(pgp + 70);
    const auto *pgp_71 = buffer.data(pgp + 71);
    const auto *pgp_73 = buffer.data(pgp + 73);
    const auto *pgp_74 = buffer.data(pgp + 74);
    const auto *pgp_75 = buffer.data(pgp + 75);
    const auto *pgp_76 = buffer.data(pgp + 76);
    const auto *pgp_77 = buffer.data(pgp + 77);
    const auto *pgp_79 = buffer.data(pgp + 79);
    const auto *pgp_80 = buffer.data(pgp + 80);
    const auto *pgp_81 = buffer.data(pgp + 81);
    const auto *pgp_82 = buffer.data(pgp + 82);
    const auto *pgp_83 = buffer.data(pgp + 83);
    const auto *pgp_84 = buffer.data(pgp + 84);
    const auto *pgp_85 = buffer.data(pgp + 85);
    const auto *pgp_86 = buffer.data(pgp + 86);
    const auto *pgp_88 = buffer.data(pgp + 88);
    const auto *pgp_89 = buffer.data(pgp + 89);
    const auto *pgp_91 = buffer.data(pgp + 91);
    const auto *pgp_92 = buffer.data(pgp + 92);
    const auto *pgp_94 = buffer.data(pgp + 94);
    const auto *pgp_95 = buffer.data(pgp + 95);
    const auto *pgp_96 = buffer.data(pgp + 96);
    const auto *pgp_97 = buffer.data(pgp + 97);
    const auto *pgp_98 = buffer.data(pgp + 98);
    const auto *pgp_100 = buffer.data(pgp + 100);
    const auto *pgp_101 = buffer.data(pgp + 101);
    const auto *pgp_103 = buffer.data(pgp + 103);
    const auto *pgp_104 = buffer.data(pgp + 104);
    const auto *pgp_105 = buffer.data(pgp + 105);
    const auto *pgp_106 = buffer.data(pgp + 106);
    const auto *pgp_107 = buffer.data(pgp + 107);
    const auto *pgp_109 = buffer.data(pgp + 109);
    const auto *pgp_110 = buffer.data(pgp + 110);
    const auto *pgp_111 = buffer.data(pgp + 111);
    const auto *pgp_112 = buffer.data(pgp + 112);
    const auto *pgp_113 = buffer.data(pgp + 113);
    const auto *pgp_115 = buffer.data(pgp + 115);
    const auto *pgp_116 = buffer.data(pgp + 116);
    const auto *pgp_118 = buffer.data(pgp + 118);
    const auto *pgp_119 = buffer.data(pgp + 119);
    const auto *pgp_121 = buffer.data(pgp + 121);
    const auto *pgp_122 = buffer.data(pgp + 122);
    const auto *pgp_123 = buffer.data(pgp + 123);
    const auto *pgp_124 = buffer.data(pgp + 124);
    const auto *pgp_125 = buffer.data(pgp + 125);

#pragma omp simd aligned(t_124, t_125, t_126, pa_y, pb_x, pc_x, pc_y, pc_z, sgd0_35, sgd1_35, \
                         pfd0_96, pfp_37, pfp_48, pfd1_96, pgp_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_8 * pfp_37[k]
                   + f_4 * pc_z[k] * pgp_61[k];

        t_125[k] = pa_y[k] * sgd0_35[k]
                   - f_5 * pc_y[k] * sgd1_35[k];

        t_126[k] = pb_x[k] * pfd0_96[k]
                   + f_8 * pfp_48[k]
                   - f_5 * pc_x[k] * pfd1_96[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pb_x, pc_x, pc_z, pfd0_99, \
                         pfd0_101, pfp_49, pfp_50, pfd1_99, pfd1_101, pgp_64, \
                         pgp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_0 * pfp_49[k]
                   + f_4 * pc_x[k] * pgp_64[k];

        t_128[k] = f_0 * pfp_50[k]
                   + f_4 * pc_x[k] * pgp_65[k];

        t_129[k] = pb_x[k] * pfd0_99[k]
                   - f_5 * pc_x[k] * pfd1_99[k];

        t_130[k] = f_4 * pc_z[k] * pgp_64[k];

        t_131[k] = pb_x[k] * pfd0_101[k]
                   - f_5 * pc_x[k] * pfd1_101[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, pb_z, pc_x, pc_z, pfd0_78, \
                         pfd0_105, pfp_52, pfp_53, pfd1_78, pfd1_105, pgp_67, \
                         pgp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_z[k] * pfd0_78[k]
                   - f_5 * pc_z[k] * pfd1_78[k];

        t_133[k] = f_0 * pfp_52[k]
                   + f_4 * pc_x[k] * pgp_67[k];

        t_134[k] = f_0 * pfp_53[k]
                   + f_4 * pc_x[k] * pgp_68[k];

        t_135[k] = pb_x[k] * pfd0_105[k]
                   - f_5 * pc_x[k] * pfd1_105[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pb_x, pc_x, pc_z, pfd0_107, pfp_40, pfp_54, \
                         pfd1_107, pgs0_23, pgs1_23, pgp_67, pgp_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_0 * pfp_40[k]
                   + f_4 * pc_z[k] * pgp_67[k];

        t_137[k] = pb_x[k] * pfd0_107[k]
                   - f_5 * pc_x[k] * pfd1_107[k];

        t_138[k] = f_0 * pfp_54[k]
                   + f_2 * pgs0_23[k]
                   - f_3 * pgs1_23[k]
                   + f_4 * pc_x[k] * pgp_69[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_x, pc_x, pc_z, pfd0_111, pfp_43, \
                         pfp_55, pfp_56, pfd1_111, pgp_70, pgp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_0 * pfp_55[k]
                   + f_4 * pc_x[k] * pgp_70[k];

        t_140[k] = f_0 * pfp_56[k]
                   + f_4 * pc_x[k] * pgp_71[k];

        t_141[k] = pb_x[k] * pfd0_111[k]
                   - f_5 * pc_x[k] * pfd1_111[k];

        t_142[k] = f_8 * pfp_43[k]
                   + f_4 * pc_z[k] * pgp_70[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_y, pb_x, pc_x, pc_y, sgd0_54, sgd1_54, \
                         pfd0_113, pfp_58, pfp_59, pfd1_113, pgp_73, \
                         pgp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = pb_x[k] * pfd0_113[k]
                   - f_5 * pc_x[k] * pfd1_113[k];

        t_144[k] = pa_y[k] * sgd0_54[k]
                   - f_5 * pc_y[k] * sgd1_54[k];

        t_145[k] = f_0 * pfp_58[k]
                   + f_4 * pc_x[k] * pgp_73[k];

        t_146[k] = f_0 * pfp_59[k]
                   + f_4 * pc_x[k] * pgp_74[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pa_y, pb_x, pc_x, pc_y, pc_z, sgd0_59, sgd1_59, \
                         pfd0_117, pfp_46, pfd1_117, pgp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pb_x[k] * pfd0_117[k]
                   - f_5 * pc_x[k] * pfd1_117[k];

        t_148[k] = f_9 * pfp_46[k]
                   + f_4 * pc_z[k] * pgp_73[k];

        t_149[k] = pa_y[k] * sgd0_59[k]
                   - f_5 * pc_y[k] * sgd1_59[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, pc_x, pc_y, pc_z, sgp_31, \
                         pfp_49, pgs0_25, pgs1_25, pgp_75, pgp_76, \
                         pgp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_2 * pgs0_25[k]
                   - f_3 * pgs1_25[k]
                   + f_4 * pc_x[k] * pgp_75[k];

        t_151[k] = f_4 * pc_x[k] * pgp_76[k];

        t_152[k] = f_4 * pc_x[k] * pgp_77[k];

        t_153[k] = f_0 * sgp_31[k]
                   + f_1 * pfp_49[k]
                   + f_2 * pgs0_25[k]
                   - f_3 * pgs1_25[k]
                   + f_4 * pc_y[k] * pgp_76[k];

        t_154[k] = f_4 * pc_z[k] * pgp_76[k];

        t_155[k] = f_2 * pgs0_25[k]
                   - f_3 * pgs1_25[k]
                   + f_4 * pc_z[k] * pgp_77[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pb_z, pc_x, pc_z, pfd0_96, \
                         pfd0_99, pfp_49, pfd1_96, pfd1_99, pgp_79, \
                         pgp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pb_z[k] * pfd0_96[k]
                   - f_5 * pc_z[k] * pfd1_96[k];

        t_157[k] = f_4 * pc_x[k] * pgp_79[k];

        t_158[k] = f_4 * pc_x[k] * pgp_80[k];

        t_159[k] = pb_z[k] * pfd0_99[k]
                   - f_5 * pc_z[k] * pfd1_99[k];

        t_160[k] = f_0 * pfp_49[k]
                   + f_4 * pc_z[k] * pgp_79[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pc_x, pc_z, pfp_50, pgs0_26, pgs0_27, \
                         pgs1_26, pgs1_27, pgp_80, pgp_81, pgp_82, \
                         pgp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_0 * pfp_50[k]
                   + f_2 * pgs0_26[k]
                   - f_3 * pgs1_26[k]
                   + f_4 * pc_z[k] * pgp_80[k];

        t_162[k] = f_2 * pgs0_27[k]
                   - f_3 * pgs1_27[k]
                   + f_4 * pc_x[k] * pgp_81[k];

        t_163[k] = f_4 * pc_x[k] * pgp_82[k];

        t_164[k] = f_4 * pc_x[k] * pgp_83[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pb_z, pc_z, pdd0_57, pdd1_57, pfd0_105, pfp_52, \
                         pfp_53, pfd1_105, pgs0_27, pgs1_27, pgp_82, \
                         pgp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_6 * pdd0_57[k]
                   - f_7 * pdd1_57[k]
                   + pb_z[k] * pfd0_105[k]
                   - f_5 * pc_z[k] * pfd1_105[k];

        t_166[k] = f_8 * pfp_52[k]
                   + f_4 * pc_z[k] * pgp_82[k];

        t_167[k] = f_8 * pfp_53[k]
                   + f_2 * pgs0_27[k]
                   - f_3 * pgs1_27[k]
                   + f_4 * pc_z[k] * pgp_83[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pc_x, pc_y, pc_z, sgp_40, pfp_55, \
                         pfp_58, pgs0_28, pgs1_28, pgp_84, pgp_85, \
                         pgp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_2 * pgs0_28[k]
                   - f_3 * pgs1_28[k]
                   + f_4 * pc_x[k] * pgp_84[k];

        t_169[k] = f_4 * pc_x[k] * pgp_85[k];

        t_170[k] = f_4 * pc_x[k] * pgp_86[k];

        t_171[k] = f_0 * sgp_40[k]
                   + f_0 * pfp_58[k]
                   + f_2 * pgs0_28[k]
                   - f_3 * pgs1_28[k]
                   + f_4 * pc_y[k] * pgp_85[k];

        t_172[k] = f_9 * pfp_55[k]
                   + f_4 * pc_z[k] * pgp_85[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_y, pc_x, pc_y, pc_z, sgd0_84, sgd1_84, \
                         pfp_56, pgs0_28, pgs1_28, pgp_86, pgp_88, \
                         pgp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_9 * pfp_56[k]
                   + f_2 * pgs0_28[k]
                   - f_3 * pgs1_28[k]
                   + f_4 * pc_z[k] * pgp_86[k];

        t_174[k] = pa_y[k] * sgd0_84[k]
                   - f_5 * pc_y[k] * sgd1_84[k];

        t_175[k] = f_4 * pc_x[k] * pgp_88[k];

        t_176[k] = f_4 * pc_x[k] * pgp_89[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pc_y, pc_z, sgd0_87, sgd0_89, sgp_43, \
                         sgd1_87, sgd1_89, pfp_58, pgp_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pa_y[k] * sgd0_87[k]
                   + f_8 * sgp_43[k]
                   - f_5 * pc_y[k] * sgd1_87[k];

        t_178[k] = f_1 * pfp_58[k]
                   + f_4 * pc_z[k] * pgp_88[k];

        t_179[k] = pa_y[k] * sgd0_89[k]
                   - f_5 * pc_y[k] * sgd1_89[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pc_x, pc_z, sgd0_0, sgd0_3, sgd1_0, \
                         sgd1_3, pfp_61, pfp_62, pgp_91, pgp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_z[k] * sgd0_0[k]
                   - f_5 * pc_z[k] * sgd1_0[k];

        t_181[k] = f_1 * pfp_61[k]
                   + f_4 * pc_x[k] * pgp_91[k];

        t_182[k] = f_1 * pfp_62[k]
                   + f_4 * pc_x[k] * pgp_92[k];

        t_183[k] = pa_z[k] * sgd0_3[k]
                   - f_5 * pc_z[k] * sgd1_3[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_z, pc_x, pc_y, pc_z, sgd0_5, sgd0_6, \
                         sgp_2, sgd1_5, sgd1_6, pfp_64, pgp_92, \
                         pgp_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_4 * pc_y[k] * pgp_92[k];

        t_185[k] = pa_z[k] * sgd0_5[k]
                   + f_8 * sgp_2[k]
                   - f_5 * pc_z[k] * sgd1_5[k];

        t_186[k] = pa_z[k] * sgd0_6[k]
                   - f_5 * pc_z[k] * sgd1_6[k];

        t_187[k] = f_9 * pfp_64[k]
                   + f_4 * pc_x[k] * pgp_94[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_z, pb_y, pc_x, pc_y, pc_z, sgd0_9, \
                         sgd1_9, pfd0_125, pfp_62, pfp_65, pfd1_125, \
                         pgp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_9 * pfp_65[k]
                   + f_4 * pc_x[k] * pgp_95[k];

        t_189[k] = pa_z[k] * sgd0_9[k]
                   - f_5 * pc_z[k] * sgd1_9[k];

        t_190[k] = f_0 * pfp_62[k]
                   + f_4 * pc_y[k] * pgp_95[k];

        t_191[k] = pb_y[k] * pfd0_125[k]
                   - f_5 * pc_y[k] * pfd1_125[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, pc_x, pc_y, pfp_66, pfp_67, \
                         pfp_68, pgs0_32, pgs1_32, pgp_96, pgp_97, \
                         pgp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_9 * pfp_66[k]
                   + f_2 * pgs0_32[k]
                   - f_3 * pgs1_32[k]
                   + f_4 * pc_x[k] * pgp_96[k];

        t_193[k] = f_9 * pfp_67[k]
                   + f_4 * pc_x[k] * pgp_97[k];

        t_194[k] = f_9 * pfp_68[k]
                   + f_4 * pc_x[k] * pgp_98[k];

        t_195[k] = f_2 * pgs0_32[k]
                   - f_3 * pgs1_32[k]
                   + f_4 * pc_y[k] * pgp_97[k];

        t_196[k] = f_4 * pc_y[k] * pgp_98[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_z, pb_x, pc_x, pc_z, sgd0_18, sgd1_18, \
                         pdd0_89, pdd1_89, pfd0_137, pfp_70, pfd1_137, \
                         pgp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_10 * pdd0_89[k]
                   - f_11 * pdd1_89[k]
                   + pb_x[k] * pfd0_137[k]
                   - f_5 * pc_x[k] * pfd1_137[k];

        t_198[k] = pa_z[k] * sgd0_18[k]
                   - f_5 * pc_z[k] * sgd1_18[k];

        t_199[k] = f_8 * pfp_70[k]
                   + f_4 * pc_x[k] * pgp_100[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_z, pc_x, pc_y, pc_z, sgd0_21, sgd0_23, \
                         sgp_11, sgd1_21, sgd1_23, pfp_65, pfp_71, \
                         pgp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_8 * pfp_71[k]
                   + f_4 * pc_x[k] * pgp_101[k];

        t_201[k] = pa_z[k] * sgd0_21[k]
                   - f_5 * pc_z[k] * sgd1_21[k];

        t_202[k] = f_8 * pfp_65[k]
                   + f_4 * pc_y[k] * pgp_101[k];

        t_203[k] = pa_z[k] * sgd0_23[k]
                   + f_8 * sgp_11[k]
                   - f_5 * pc_z[k] * sgd1_23[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pc_x, pc_y, pfd0_132, pfp_67, \
                         pfp_73, pfp_74, pfd1_132, pgs0_34, pgs1_34, pgp_103, \
                         pgp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_y[k] * pfd0_132[k]
                   - f_5 * pc_y[k] * pfd1_132[k];

        t_205[k] = f_8 * pfp_73[k]
                   + f_4 * pc_x[k] * pgp_103[k];

        t_206[k] = f_8 * pfp_74[k]
                   + f_4 * pc_x[k] * pgp_104[k];

        t_207[k] = f_0 * pfp_67[k]
                   + f_2 * pgs0_34[k]
                   - f_3 * pgs1_34[k]
                   + f_4 * pc_y[k] * pgp_103[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pb_y, pc_x, pc_y, pfd0_137, pfp_68, pfp_75, \
                         pfd1_137, pgs0_35, pgs1_35, pgp_104, pgp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_0 * pfp_68[k]
                   + f_4 * pc_y[k] * pgp_104[k];

        t_209[k] = pb_y[k] * pfd0_137[k]
                   - f_5 * pc_y[k] * pfd1_137[k];

        t_210[k] = f_8 * pfp_75[k]
                   + f_2 * pgs0_35[k]
                   - f_3 * pgs1_35[k]
                   + f_4 * pc_x[k] * pgp_105[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, pc_x, pc_y, pfp_76, pfp_77, pgs0_35, \
                         pgs1_35, pgp_106, pgp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_8 * pfp_76[k]
                   + f_4 * pc_x[k] * pgp_106[k];

        t_212[k] = f_8 * pfp_77[k]
                   + f_4 * pc_x[k] * pgp_107[k];

        t_213[k] = f_2 * pgs0_35[k]
                   - f_3 * pgs1_35[k]
                   + f_4 * pc_y[k] * pgp_106[k];

        t_214[k] = f_4 * pc_y[k] * pgp_107[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pa_z, pb_x, pc_x, pc_z, sgd0_36, sgd1_36, \
                         pdd0_107, pdd1_107, pfd0_155, pfp_79, pfd1_155, \
                         pgp_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_6 * pdd0_107[k]
                   - f_7 * pdd1_107[k]
                   + pb_x[k] * pfd0_155[k]
                   - f_5 * pc_x[k] * pfd1_155[k];

        t_216[k] = pa_z[k] * sgd0_36[k]
                   - f_5 * pc_z[k] * sgd1_36[k];

        t_217[k] = f_0 * pfp_79[k]
                   + f_4 * pc_x[k] * pgp_109[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pa_z, pb_x, pc_x, pc_y, pc_z, sgd0_39, \
                         sgd1_39, pfd0_161, pfp_71, pfp_80, pfd1_161, \
                         pgp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * pfp_80[k]
                   + f_4 * pc_x[k] * pgp_110[k];

        t_219[k] = pa_z[k] * sgd0_39[k]
                   - f_5 * pc_z[k] * sgd1_39[k];

        t_220[k] = f_9 * pfp_71[k]
                   + f_4 * pc_y[k] * pgp_110[k];

        t_221[k] = pb_x[k] * pfd0_161[k]
                   - f_5 * pc_x[k] * pfd1_161[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pb_x, pc_x, pfd0_165, pfp_81, pfp_82, \
                         pfp_83, pfd1_165, pgs0_37, pgs1_37, pgp_111, pgp_112, \
                         pgp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_0 * pfp_81[k]
                   + f_2 * pgs0_37[k]
                   - f_3 * pgs1_37[k]
                   + f_4 * pc_x[k] * pgp_111[k];

        t_223[k] = f_0 * pfp_82[k]
                   + f_4 * pc_x[k] * pgp_112[k];

        t_224[k] = f_0 * pfp_83[k]
                   + f_4 * pc_x[k] * pgp_113[k];

        t_225[k] = pb_x[k] * pfd0_165[k]
                   - f_5 * pc_x[k] * pfd1_165[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pb_x, pb_y, pc_x, pc_y, pfd0_150, \
                         pfd0_167, pfp_74, pfp_85, pfd1_150, pfd1_167, pgp_113, \
                         pgp_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_8 * pfp_74[k]
                   + f_4 * pc_y[k] * pgp_113[k];

        t_227[k] = pb_x[k] * pfd0_167[k]
                   - f_5 * pc_x[k] * pfd1_167[k];

        t_228[k] = pb_y[k] * pfd0_150[k]
                   - f_5 * pc_y[k] * pfd1_150[k];

        t_229[k] = f_0 * pfp_85[k]
                   + f_4 * pc_x[k] * pgp_115[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pb_x, pc_x, pc_y, pfd0_171, pfd0_173, \
                         pfp_77, pfp_86, pfd1_171, pfd1_173, pgp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_0 * pfp_86[k]
                   + f_4 * pc_x[k] * pgp_116[k];

        t_231[k] = pb_x[k] * pfd0_171[k]
                   - f_5 * pc_x[k] * pfd1_171[k];

        t_232[k] = f_0 * pfp_77[k]
                   + f_4 * pc_y[k] * pgp_116[k];

        t_233[k] = pb_x[k] * pfd0_173[k]
                   - f_5 * pc_x[k] * pfd1_173[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pb_x, pc_x, pfd0_174, pfd0_177, pfp_87, \
                         pfp_88, pfp_89, pfd1_174, pfd1_177, pgp_118, \
                         pgp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pb_x[k] * pfd0_174[k]
                   + f_8 * pfp_87[k]
                   - f_5 * pc_x[k] * pfd1_174[k];

        t_235[k] = f_0 * pfp_88[k]
                   + f_4 * pc_x[k] * pgp_118[k];

        t_236[k] = f_0 * pfp_89[k]
                   + f_4 * pc_x[k] * pgp_119[k];

        t_237[k] = pb_x[k] * pfd0_177[k]
                   - f_5 * pc_x[k] * pfd1_177[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_z, pb_x, pc_x, pc_y, pc_z, sgd0_60, \
                         sgd1_60, pfd0_179, pfd1_179, pgp_119, \
                         pgp_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_4 * pc_y[k] * pgp_119[k];

        t_239[k] = pb_x[k] * pfd0_179[k]
                   - f_5 * pc_x[k] * pfd1_179[k];

        t_240[k] = pa_z[k] * sgd0_60[k]
                   - f_5 * pc_z[k] * sgd1_60[k];

        t_241[k] = f_4 * pc_x[k] * pgp_121[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pa_z, pc_x, pc_y, pc_z, sgd0_63, sgd0_65, \
                         sgp_32, sgd1_63, sgd1_65, pfp_80, pgp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_4 * pc_x[k] * pgp_122[k];

        t_243[k] = pa_z[k] * sgd0_63[k]
                   - f_5 * pc_z[k] * sgd1_63[k];

        t_244[k] = f_1 * pfp_80[k]
                   + f_4 * pc_y[k] * pgp_122[k];

        t_245[k] = pa_z[k] * sgd0_65[k]
                   + f_8 * sgp_32[k]
                   - f_5 * pc_z[k] * sgd1_65[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, pc_x, pc_y, pfp_82, pfp_83, \
                         pgs0_41, pgs1_41, pgp_123, pgp_124, pgp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_2 * pgs0_41[k]
                   - f_3 * pgs1_41[k]
                   + f_4 * pc_x[k] * pgp_123[k];

        t_247[k] = f_4 * pc_x[k] * pgp_124[k];

        t_248[k] = f_4 * pc_x[k] * pgp_125[k];

        t_249[k] = f_9 * pfp_82[k]
                   + f_2 * pgs0_41[k]
                   - f_3 * pgs1_41[k]
                   + f_4 * pc_y[k] * pgp_124[k];

        t_250[k] = f_9 * pfp_83[k]
                   + f_4 * pc_y[k] * pgp_125[k];
    }
}

static auto
compute_prim_pgd_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgp,
                                                          const size_t pdd0, const size_t pdd1,
                                                          const size_t pfd0, const size_t pfp,
                                                          const size_t pfd1, const size_t pgs0,
                                                          const size_t pgs1, const size_t pgp,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 0.5 / p;
    const auto f_7 = 0.5 * gamma / (p * q);
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / p;
    const auto f_11 = gamma / (p * q);

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgp_44 = buffer.data(sgp + 44);

    const auto *pdd0_101 = buffer.data(pdd0 + 101);
    const auto *pdd0_107 = buffer.data(pdd0 + 107);

    const auto *pdd1_101 = buffer.data(pdd1 + 101);
    const auto *pdd1_107 = buffer.data(pdd1 + 107);

    const auto *pfd0_167 = buffer.data(pfd0 + 167);
    const auto *pfd0_173 = buffer.data(pfd0 + 173);
    const auto *pfd0_174 = buffer.data(pfd0 + 174);
    const auto *pfd0_177 = buffer.data(pfd0 + 177);
    const auto *pfd0_179 = buffer.data(pfd0 + 179);

    const auto *pfp_85 = buffer.data(pfp + 85);
    const auto *pfp_86 = buffer.data(pfp + 86);
    const auto *pfp_88 = buffer.data(pfp + 88);
    const auto *pfp_89 = buffer.data(pfp + 89);

    const auto *pfd1_167 = buffer.data(pfd1 + 167);
    const auto *pfd1_173 = buffer.data(pfd1 + 173);
    const auto *pfd1_174 = buffer.data(pfd1 + 174);
    const auto *pfd1_177 = buffer.data(pfd1 + 177);
    const auto *pfd1_179 = buffer.data(pfd1 + 179);

    const auto *pgs0_42 = buffer.data(pgs0 + 42);
    const auto *pgs0_44 = buffer.data(pgs0 + 44);

    const auto *pgs1_42 = buffer.data(pgs1 + 42);
    const auto *pgs1_44 = buffer.data(pgs1 + 44);

    const auto *pgp_126 = buffer.data(pgp + 126);
    const auto *pgp_127 = buffer.data(pgp + 127);
    const auto *pgp_128 = buffer.data(pgp + 128);
    const auto *pgp_130 = buffer.data(pgp + 130);
    const auto *pgp_131 = buffer.data(pgp + 131);
    const auto *pgp_132 = buffer.data(pgp + 132);
    const auto *pgp_133 = buffer.data(pgp + 133);
    const auto *pgp_134 = buffer.data(pgp + 134);

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_x, pc_y, pdd0_101, pdd1_101, \
                         pfd0_167, pfd1_167, pgs0_42, pgs1_42, pgp_126, pgp_127, \
                         pgp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_10 * pdd0_101[k]
                   - f_11 * pdd1_101[k]
                   + pb_y[k] * pfd0_167[k]
                   - f_5 * pc_y[k] * pfd1_167[k];

        t_252[k] = f_2 * pgs0_42[k]
                   - f_3 * pgs1_42[k]
                   + f_4 * pc_x[k] * pgp_126[k];

        t_253[k] = f_4 * pc_x[k] * pgp_127[k];

        t_254[k] = f_4 * pc_x[k] * pgp_128[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pb_y, pc_y, pdd0_107, pdd1_107, pfd0_173, \
                         pfp_85, pfp_86, pfd1_173, pgs0_42, pgs1_42, pgp_127, \
                         pgp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_8 * pfp_85[k]
                   + f_2 * pgs0_42[k]
                   - f_3 * pgs1_42[k]
                   + f_4 * pc_y[k] * pgp_127[k];

        t_256[k] = f_8 * pfp_86[k]
                   + f_4 * pc_y[k] * pgp_128[k];

        t_257[k] = f_6 * pdd0_107[k]
                   - f_7 * pdd1_107[k]
                   + pb_y[k] * pfd0_173[k]
                   - f_5 * pc_y[k] * pfd1_173[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, pb_y, pc_x, pc_y, pfd0_174, \
                         pfd0_177, pfp_88, pfp_89, pfd1_174, pfd1_177, pgp_130, \
                         pgp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pb_y[k] * pfd0_174[k]
                   - f_5 * pc_y[k] * pfd1_174[k];

        t_259[k] = f_4 * pc_x[k] * pgp_130[k];

        t_260[k] = f_4 * pc_x[k] * pgp_131[k];

        t_261[k] = pb_y[k] * pfd0_177[k]
                   + f_8 * pfp_88[k]
                   - f_5 * pc_y[k] * pfd1_177[k];

        t_262[k] = f_0 * pfp_89[k]
                   + f_4 * pc_y[k] * pgp_131[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, t_268, pb_y, pc_x, pc_y, pfd0_179, \
                         pfd1_179, pgs0_44, pgs1_44, pgp_132, pgp_133, \
                         pgp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = pb_y[k] * pfd0_179[k]
                   - f_5 * pc_y[k] * pfd1_179[k];

        t_264[k] = f_2 * pgs0_44[k]
                   - f_3 * pgs1_44[k]
                   + f_4 * pc_x[k] * pgp_132[k];

        t_265[k] = f_4 * pc_x[k] * pgp_133[k];

        t_266[k] = f_4 * pc_x[k] * pgp_134[k];

        t_267[k] = f_2 * pgs0_44[k]
                   - f_3 * pgs1_44[k]
                   + f_4 * pc_y[k] * pgp_133[k];

        t_268[k] = f_4 * pc_y[k] * pgp_134[k];
    }

#pragma omp simd aligned(t_269, pc_z, sgp_44, pfp_89, pgs0_44, pgs1_44, \
                         pgp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_0 * sgp_44[k]
                   + f_1 * pfp_89[k]
                   + f_2 * pgs0_44[k]
                   - f_3 * pgs1_44[k]
                   + f_4 * pc_z[k] * pgp_134[k];
    }
}

auto
compute_prim_pgd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sgd0,
                                                   const size_t sgp, const size_t sgd1,
                                                   const size_t pdd0, const size_t pdd1,
                                                   const size_t pfd0, const size_t pfp,
                                                   const size_t pfd1, const size_t pgs0,
                                                   const size_t pgs1, const size_t pgp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pgd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sgd0,
                                                              sgp, sgd1, pdd0, pdd1, pfd0, pfp,
                                                              pfd1, pgs0, pgs1, pgp, ncols,
                                                              gamma, p, q);

    compute_prim_pgd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sgd0,
                                                              sgp, sgd1, pdd0, pdd1, pfd0, pfp,
                                                              pfd1, pgs0, pgs1, pgp, ncols,
                                                              gamma, p, q);

    compute_prim_pgd_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sgp, pdd0,
                                                              pdd1, pfd0, pfp, pfd1, pgs0, pgs1,
                                                              pgp, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
