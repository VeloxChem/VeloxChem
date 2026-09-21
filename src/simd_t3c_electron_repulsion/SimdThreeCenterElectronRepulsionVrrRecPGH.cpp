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


#include "SimdThreeCenterElectronRepulsionVrrRecPGH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pgh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgg,
                                                          const size_t pdh0, const size_t pdh1,
                                                          const size_t pfh0, const size_t pfg,
                                                          const size_t pfh1, const size_t pgf0,
                                                          const size_t pgf1, const size_t pgg,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgg_0 = buffer.data(sgg + 0);
    const auto *sgg_10 = buffer.data(sgg + 10);
    const auto *sgg_12 = buffer.data(sgg + 12);
    const auto *sgg_14 = buffer.data(sgg + 14);
    const auto *sgg_25 = buffer.data(sgg + 25);
    const auto *sgg_27 = buffer.data(sgg + 27);
    const auto *sgg_42 = buffer.data(sgg + 42);
    const auto *sgg_44 = buffer.data(sgg + 44);
    const auto *sgg_55 = buffer.data(sgg + 55);
    const auto *sgg_57 = buffer.data(sgg + 57);
    const auto *sgg_59 = buffer.data(sgg + 59);
    const auto *sgg_72 = buffer.data(sgg + 72);
    const auto *sgg_85 = buffer.data(sgg + 85);
    const auto *sgg_87 = buffer.data(sgg + 87);
    const auto *sgg_89 = buffer.data(sgg + 89);

    const auto *pdh0_0 = buffer.data(pdh0 + 0);

    const auto *pdh1_0 = buffer.data(pdh1 + 0);

    const auto *pfh0_0 = buffer.data(pfh0 + 0);
    const auto *pfh0_3 = buffer.data(pfh0 + 3);
    const auto *pfh0_5 = buffer.data(pfh0 + 5);
    const auto *pfh0_6 = buffer.data(pfh0 + 6);
    const auto *pfh0_9 = buffer.data(pfh0 + 9);
    const auto *pfh0_10 = buffer.data(pfh0 + 10);
    const auto *pfh0_14 = buffer.data(pfh0 + 14);
    const auto *pfh0_15 = buffer.data(pfh0 + 15);
    const auto *pfh0_20 = buffer.data(pfh0 + 20);
    const auto *pfh0_21 = buffer.data(pfh0 + 21);
    const auto *pfh0_24 = buffer.data(pfh0 + 24);
    const auto *pfh0_27 = buffer.data(pfh0 + 27);
    const auto *pfh0_31 = buffer.data(pfh0 + 31);
    const auto *pfh0_36 = buffer.data(pfh0 + 36);
    const auto *pfh0_42 = buffer.data(pfh0 + 42);
    const auto *pfh0_47 = buffer.data(pfh0 + 47);
    const auto *pfh0_51 = buffer.data(pfh0 + 51);
    const auto *pfh0_56 = buffer.data(pfh0 + 56);
    const auto *pfh0_62 = buffer.data(pfh0 + 62);

    const auto *pfg_0 = buffer.data(pfg + 0);
    const auto *pfg_1 = buffer.data(pfg + 1);
    const auto *pfg_2 = buffer.data(pfg + 2);
    const auto *pfg_3 = buffer.data(pfg + 3);
    const auto *pfg_5 = buffer.data(pfg + 5);
    const auto *pfg_6 = buffer.data(pfg + 6);
    const auto *pfg_9 = buffer.data(pfg + 9);
    const auto *pfg_10 = buffer.data(pfg + 10);
    const auto *pfg_12 = buffer.data(pfg + 12);
    const auto *pfg_13 = buffer.data(pfg + 13);
    const auto *pfg_14 = buffer.data(pfg + 14);
    const auto *pfg_15 = buffer.data(pfg + 15);
    const auto *pfg_16 = buffer.data(pfg + 16);
    const auto *pfg_17 = buffer.data(pfg + 17);
    const auto *pfg_18 = buffer.data(pfg + 18);
    const auto *pfg_20 = buffer.data(pfg + 20);
    const auto *pfg_21 = buffer.data(pfg + 21);
    const auto *pfg_24 = buffer.data(pfg + 24);
    const auto *pfg_25 = buffer.data(pfg + 25);
    const auto *pfg_27 = buffer.data(pfg + 27);
    const auto *pfg_28 = buffer.data(pfg + 28);
    const auto *pfg_29 = buffer.data(pfg + 29);
    const auto *pfg_30 = buffer.data(pfg + 30);
    const auto *pfg_32 = buffer.data(pfg + 32);
    const auto *pfg_33 = buffer.data(pfg + 33);
    const auto *pfg_35 = buffer.data(pfg + 35);
    const auto *pfg_36 = buffer.data(pfg + 36);
    const auto *pfg_39 = buffer.data(pfg + 39);
    const auto *pfg_40 = buffer.data(pfg + 40);
    const auto *pfg_42 = buffer.data(pfg + 42);
    const auto *pfg_43 = buffer.data(pfg + 43);
    const auto *pfg_44 = buffer.data(pfg + 44);
    const auto *pfg_55 = buffer.data(pfg + 55);
    const auto *pfg_57 = buffer.data(pfg + 57);
    const auto *pfg_59 = buffer.data(pfg + 59);
    const auto *pfg_72 = buffer.data(pfg + 72);
    const auto *pfg_85 = buffer.data(pfg + 85);
    const auto *pfg_87 = buffer.data(pfg + 87);
    const auto *pfg_89 = buffer.data(pfg + 89);

    const auto *pfh1_0 = buffer.data(pfh1 + 0);
    const auto *pfh1_3 = buffer.data(pfh1 + 3);
    const auto *pfh1_5 = buffer.data(pfh1 + 5);
    const auto *pfh1_6 = buffer.data(pfh1 + 6);
    const auto *pfh1_9 = buffer.data(pfh1 + 9);
    const auto *pfh1_10 = buffer.data(pfh1 + 10);
    const auto *pfh1_14 = buffer.data(pfh1 + 14);
    const auto *pfh1_15 = buffer.data(pfh1 + 15);
    const auto *pfh1_20 = buffer.data(pfh1 + 20);
    const auto *pfh1_21 = buffer.data(pfh1 + 21);
    const auto *pfh1_24 = buffer.data(pfh1 + 24);
    const auto *pfh1_27 = buffer.data(pfh1 + 27);
    const auto *pfh1_31 = buffer.data(pfh1 + 31);
    const auto *pfh1_36 = buffer.data(pfh1 + 36);
    const auto *pfh1_42 = buffer.data(pfh1 + 42);
    const auto *pfh1_47 = buffer.data(pfh1 + 47);
    const auto *pfh1_51 = buffer.data(pfh1 + 51);
    const auto *pfh1_56 = buffer.data(pfh1 + 56);
    const auto *pfh1_62 = buffer.data(pfh1 + 62);

    const auto *pgf0_0 = buffer.data(pgf0 + 0);
    const auto *pgf0_1 = buffer.data(pgf0 + 1);
    const auto *pgf0_2 = buffer.data(pgf0 + 2);
    const auto *pgf0_6 = buffer.data(pgf0 + 6);
    const auto *pgf0_8 = buffer.data(pgf0 + 8);
    const auto *pgf0_9 = buffer.data(pgf0 + 9);
    const auto *pgf0_16 = buffer.data(pgf0 + 16);
    const auto *pgf0_18 = buffer.data(pgf0 + 18);
    const auto *pgf0_19 = buffer.data(pgf0 + 19);
    const auto *pgf0_28 = buffer.data(pgf0 + 28);
    const auto *pgf0_29 = buffer.data(pgf0 + 29);
    const auto *pgf0_30 = buffer.data(pgf0 + 30);
    const auto *pgf0_31 = buffer.data(pgf0 + 31);
    const auto *pgf0_32 = buffer.data(pgf0 + 32);
    const auto *pgf0_36 = buffer.data(pgf0 + 36);
    const auto *pgf0_38 = buffer.data(pgf0 + 38);
    const auto *pgf0_39 = buffer.data(pgf0 + 39);
    const auto *pgf0_48 = buffer.data(pgf0 + 48);
    const auto *pgf0_49 = buffer.data(pgf0 + 49);
    const auto *pgf0_50 = buffer.data(pgf0 + 50);
    const auto *pgf0_51 = buffer.data(pgf0 + 51);
    const auto *pgf0_52 = buffer.data(pgf0 + 52);
    const auto *pgf0_56 = buffer.data(pgf0 + 56);
    const auto *pgf0_58 = buffer.data(pgf0 + 58);
    const auto *pgf0_59 = buffer.data(pgf0 + 59);

    const auto *pgf1_0 = buffer.data(pgf1 + 0);
    const auto *pgf1_1 = buffer.data(pgf1 + 1);
    const auto *pgf1_2 = buffer.data(pgf1 + 2);
    const auto *pgf1_6 = buffer.data(pgf1 + 6);
    const auto *pgf1_8 = buffer.data(pgf1 + 8);
    const auto *pgf1_9 = buffer.data(pgf1 + 9);
    const auto *pgf1_16 = buffer.data(pgf1 + 16);
    const auto *pgf1_18 = buffer.data(pgf1 + 18);
    const auto *pgf1_19 = buffer.data(pgf1 + 19);
    const auto *pgf1_28 = buffer.data(pgf1 + 28);
    const auto *pgf1_29 = buffer.data(pgf1 + 29);
    const auto *pgf1_30 = buffer.data(pgf1 + 30);
    const auto *pgf1_31 = buffer.data(pgf1 + 31);
    const auto *pgf1_32 = buffer.data(pgf1 + 32);
    const auto *pgf1_36 = buffer.data(pgf1 + 36);
    const auto *pgf1_38 = buffer.data(pgf1 + 38);
    const auto *pgf1_39 = buffer.data(pgf1 + 39);
    const auto *pgf1_48 = buffer.data(pgf1 + 48);
    const auto *pgf1_49 = buffer.data(pgf1 + 49);
    const auto *pgf1_50 = buffer.data(pgf1 + 50);
    const auto *pgf1_51 = buffer.data(pgf1 + 51);
    const auto *pgf1_52 = buffer.data(pgf1 + 52);
    const auto *pgf1_56 = buffer.data(pgf1 + 56);
    const auto *pgf1_58 = buffer.data(pgf1 + 58);
    const auto *pgf1_59 = buffer.data(pgf1 + 59);

    const auto *pgg_0 = buffer.data(pgg + 0);
    const auto *pgg_1 = buffer.data(pgg + 1);
    const auto *pgg_2 = buffer.data(pgg + 2);
    const auto *pgg_3 = buffer.data(pgg + 3);
    const auto *pgg_5 = buffer.data(pgg + 5);
    const auto *pgg_6 = buffer.data(pgg + 6);
    const auto *pgg_9 = buffer.data(pgg + 9);
    const auto *pgg_10 = buffer.data(pgg + 10);
    const auto *pgg_12 = buffer.data(pgg + 12);
    const auto *pgg_13 = buffer.data(pgg + 13);
    const auto *pgg_14 = buffer.data(pgg + 14);
    const auto *pgg_15 = buffer.data(pgg + 15);
    const auto *pgg_17 = buffer.data(pgg + 17);
    const auto *pgg_18 = buffer.data(pgg + 18);
    const auto *pgg_20 = buffer.data(pgg + 20);
    const auto *pgg_21 = buffer.data(pgg + 21);
    const auto *pgg_24 = buffer.data(pgg + 24);
    const auto *pgg_25 = buffer.data(pgg + 25);
    const auto *pgg_27 = buffer.data(pgg + 27);
    const auto *pgg_28 = buffer.data(pgg + 28);
    const auto *pgg_29 = buffer.data(pgg + 29);
    const auto *pgg_30 = buffer.data(pgg + 30);
    const auto *pgg_32 = buffer.data(pgg + 32);
    const auto *pgg_33 = buffer.data(pgg + 33);
    const auto *pgg_35 = buffer.data(pgg + 35);
    const auto *pgg_36 = buffer.data(pgg + 36);
    const auto *pgg_39 = buffer.data(pgg + 39);
    const auto *pgg_40 = buffer.data(pgg + 40);
    const auto *pgg_42 = buffer.data(pgg + 42);
    const auto *pgg_43 = buffer.data(pgg + 43);
    const auto *pgg_44 = buffer.data(pgg + 44);
    const auto *pgg_45 = buffer.data(pgg + 45);
    const auto *pgg_46 = buffer.data(pgg + 46);
    const auto *pgg_47 = buffer.data(pgg + 47);
    const auto *pgg_48 = buffer.data(pgg + 48);
    const auto *pgg_50 = buffer.data(pgg + 50);
    const auto *pgg_51 = buffer.data(pgg + 51);
    const auto *pgg_54 = buffer.data(pgg + 54);
    const auto *pgg_55 = buffer.data(pgg + 55);
    const auto *pgg_57 = buffer.data(pgg + 57);
    const auto *pgg_58 = buffer.data(pgg + 58);
    const auto *pgg_59 = buffer.data(pgg + 59);
    const auto *pgg_60 = buffer.data(pgg + 60);
    const auto *pgg_62 = buffer.data(pgg + 62);
    const auto *pgg_63 = buffer.data(pgg + 63);
    const auto *pgg_65 = buffer.data(pgg + 65);
    const auto *pgg_66 = buffer.data(pgg + 66);
    const auto *pgg_69 = buffer.data(pgg + 69);
    const auto *pgg_70 = buffer.data(pgg + 70);
    const auto *pgg_72 = buffer.data(pgg + 72);
    const auto *pgg_73 = buffer.data(pgg + 73);
    const auto *pgg_74 = buffer.data(pgg + 74);
    const auto *pgg_75 = buffer.data(pgg + 75);
    const auto *pgg_76 = buffer.data(pgg + 76);
    const auto *pgg_77 = buffer.data(pgg + 77);
    const auto *pgg_78 = buffer.data(pgg + 78);
    const auto *pgg_80 = buffer.data(pgg + 80);
    const auto *pgg_81 = buffer.data(pgg + 81);
    const auto *pgg_84 = buffer.data(pgg + 84);
    const auto *pgg_85 = buffer.data(pgg + 85);
    const auto *pgg_87 = buffer.data(pgg + 87);
    const auto *pgg_88 = buffer.data(pgg + 88);
    const auto *pgg_89 = buffer.data(pgg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sgg_0, pfg_0, pgf0_0, \
                         pgf1_0, pgg_0, pgg_1, pgg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgg_0[k]
                 + f_1 * pfg_0[k]
                 + f_2 * pgf0_0[k]
                 - f_3 * pgf1_0[k]
                 + f_4 * pc_x[k] * pgg_0[k];

        t_1[k] = f_4 * pc_y[k] * pgg_0[k];

        t_2[k] = f_4 * pc_z[k] * pgg_0[k];

        t_3[k] = f_5 * pgf0_0[k]
                 - f_6 * pgf1_0[k]
                 + f_4 * pc_y[k] * pgg_1[k];

        t_4[k] = f_4 * pc_y[k] * pgg_2[k];

        t_5[k] = f_5 * pgf0_0[k]
                 - f_6 * pgf1_0[k]
                 + f_4 * pc_z[k] * pgg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_y, pc_z, pgf0_1, pgf0_2, pgf1_1, pgf1_2, \
                         pgg_3, pgg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * pgf0_1[k]
                 - f_8 * pgf1_1[k]
                 + f_4 * pc_y[k] * pgg_3[k];

        t_7[k] = f_4 * pc_z[k] * pgg_3[k];

        t_8[k] = f_4 * pc_y[k] * pgg_5[k];

        t_9[k] = f_7 * pgf0_2[k]
                 - f_8 * pgf1_2[k]
                 + f_4 * pc_z[k] * pgg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pc_x, pc_y, pc_z, sgg_10, sgg_12, pfg_10, \
                         pfg_12, pgg_6, pgg_9, pgg_10, pgg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sgg_10[k]
                  + f_1 * pfg_10[k]
                  + f_4 * pc_x[k] * pgg_10[k];

        t_11[k] = f_4 * pc_z[k] * pgg_6[k];

        t_12[k] = f_0 * sgg_12[k]
                  + f_1 * pfg_12[k]
                  + f_4 * pc_x[k] * pgg_12[k];

        t_13[k] = f_4 * pc_y[k] * pgg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sgg_14, pfg_14, pgf0_6, \
                         pgf0_8, pgf1_6, pgf1_8, pgg_10, pgg_12, \
                         pgg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * sgg_14[k]
                  + f_1 * pfg_14[k]
                  + f_4 * pc_x[k] * pgg_14[k];

        t_15[k] = f_2 * pgf0_6[k]
                  - f_3 * pgf1_6[k]
                  + f_4 * pc_y[k] * pgg_10[k];

        t_16[k] = f_4 * pc_z[k] * pgg_10[k];

        t_17[k] = f_7 * pgf0_8[k]
                  - f_8 * pgf1_8[k]
                  + f_4 * pc_y[k] * pgg_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pb_y, pc_y, pc_z, pfh0_0, pfg_0, \
                         pfh1_0, pgf0_9, pgf1_9, pgg_13, pgg_14, \
                         pgg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * pgf0_9[k]
                  - f_6 * pgf1_9[k]
                  + f_4 * pc_y[k] * pgg_13[k];

        t_19[k] = f_4 * pc_y[k] * pgg_14[k];

        t_20[k] = f_2 * pgf0_9[k]
                  - f_3 * pgf1_9[k]
                  + f_4 * pc_z[k] * pgg_14[k];

        t_21[k] = pb_y[k] * pfh0_0[k]
                  - f_9 * pc_y[k] * pfh1_0[k];

        t_22[k] = f_0 * pfg_0[k]
                  + f_4 * pc_y[k] * pgg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_y, pc_z, pfh0_3, pfh0_5, pfg_1, \
                         pfg_2, pfh1_3, pfh1_5, pgg_15, pgg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * pc_z[k] * pgg_15[k];

        t_24[k] = pb_y[k] * pfh0_3[k]
                  + f_10 * pfg_1[k]
                  - f_9 * pc_y[k] * pfh1_3[k];

        t_25[k] = f_0 * pfg_2[k]
                  + f_4 * pc_y[k] * pgg_17[k];

        t_26[k] = pb_y[k] * pfh0_5[k]
                  - f_9 * pc_y[k] * pfh1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_y, pc_y, pc_z, pfh0_6, pfh0_9, pfg_3, \
                         pfg_5, pfh1_6, pfh1_9, pgg_18, pgg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * pfh0_6[k]
                  + f_11 * pfg_3[k]
                  - f_9 * pc_y[k] * pfh1_6[k];

        t_28[k] = f_4 * pc_z[k] * pgg_18[k];

        t_29[k] = f_0 * pfg_5[k]
                  + f_4 * pc_y[k] * pgg_20[k];

        t_30[k] = pb_y[k] * pfh0_9[k]
                  - f_9 * pc_y[k] * pfh1_9[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pc_x, pc_y, pc_z, sgg_25, sgg_27, pfg_9, \
                         pfg_25, pfg_27, pgg_21, pgg_24, pgg_25, \
                         pgg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * sgg_25[k]
                  + f_11 * pfg_25[k]
                  + f_4 * pc_x[k] * pgg_25[k];

        t_32[k] = f_4 * pc_z[k] * pgg_21[k];

        t_33[k] = f_0 * sgg_27[k]
                  + f_11 * pfg_27[k]
                  + f_4 * pc_x[k] * pgg_27[k];

        t_34[k] = f_0 * pfg_9[k]
                  + f_4 * pc_y[k] * pgg_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pc_y, pc_z, pfh0_14, pfg_10, pfh1_14, \
                         pgf0_16, pgf1_16, pgg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_y[k] * pfh0_14[k]
                  - f_9 * pc_y[k] * pfh1_14[k];

        t_36[k] = f_0 * pfg_10[k]
                  + f_2 * pgf0_16[k]
                  - f_3 * pgf1_16[k]
                  + f_4 * pc_y[k] * pgg_25[k];

        t_37[k] = f_4 * pc_z[k] * pgg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pc_y, pfg_12, pfg_13, pfg_14, pgf0_18, pgf0_19, \
                         pgf1_18, pgf1_19, pgg_27, pgg_28, pgg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * pfg_12[k]
                  + f_7 * pgf0_18[k]
                  - f_8 * pgf1_18[k]
                  + f_4 * pc_y[k] * pgg_27[k];

        t_39[k] = f_0 * pfg_13[k]
                  + f_5 * pgf0_19[k]
                  - f_6 * pgf1_19[k]
                  + f_4 * pc_y[k] * pgg_28[k];

        t_40[k] = f_0 * pfg_14[k]
                  + f_4 * pc_y[k] * pgg_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pb_y, pb_z, pc_y, pc_z, pfh0_0, pfh0_20, \
                         pfg_0, pfh1_0, pfh1_20, pgg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pb_y[k] * pfh0_20[k]
                  - f_9 * pc_y[k] * pfh1_20[k];

        t_42[k] = pb_z[k] * pfh0_0[k]
                  - f_9 * pc_z[k] * pfh1_0[k];

        t_43[k] = f_4 * pc_y[k] * pgg_30[k];

        t_44[k] = f_0 * pfg_0[k]
                  + f_4 * pc_z[k] * pgg_30[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_z, pc_y, pc_z, pfh0_3, pfh0_5, pfh0_6, \
                         pfg_2, pfh1_3, pfh1_5, pfh1_6, pgg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_z[k] * pfh0_3[k]
                  - f_9 * pc_z[k] * pfh1_3[k];

        t_46[k] = f_4 * pc_y[k] * pgg_32[k];

        t_47[k] = pb_z[k] * pfh0_5[k]
                  + f_10 * pfg_2[k]
                  - f_9 * pc_z[k] * pfh1_5[k];

        t_48[k] = pb_z[k] * pfh0_6[k]
                  - f_9 * pc_z[k] * pfh1_6[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_z, pc_y, pc_z, pfh0_9, pfh0_10, pfg_3, \
                         pfg_5, pfh1_9, pfh1_10, pgg_33, pgg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_0 * pfg_3[k]
                  + f_4 * pc_z[k] * pgg_33[k];

        t_50[k] = f_4 * pc_y[k] * pgg_35[k];

        t_51[k] = pb_z[k] * pfh0_9[k]
                  + f_11 * pfg_5[k]
                  - f_9 * pc_z[k] * pfh1_9[k];

        t_52[k] = pb_z[k] * pfh0_10[k]
                  - f_9 * pc_z[k] * pfh1_10[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, sgg_42, sgg_44, pfg_6, \
                         pfg_42, pfg_44, pgg_36, pgg_39, pgg_42, \
                         pgg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * pfg_6[k]
                  + f_4 * pc_z[k] * pgg_36[k];

        t_54[k] = f_0 * sgg_42[k]
                  + f_11 * pfg_42[k]
                  + f_4 * pc_x[k] * pgg_42[k];

        t_55[k] = f_4 * pc_y[k] * pgg_39[k];

        t_56[k] = f_0 * sgg_44[k]
                  + f_11 * pfg_44[k]
                  + f_4 * pc_x[k] * pgg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_z, pc_y, pc_z, pfh0_15, pfg_10, pfh1_15, \
                         pgf0_28, pgf1_28, pgg_40, pgg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_z[k] * pfh0_15[k]
                  - f_9 * pc_z[k] * pfh1_15[k];

        t_58[k] = f_0 * pfg_10[k]
                  + f_4 * pc_z[k] * pgg_40[k];

        t_59[k] = f_7 * pgf0_28[k]
                  - f_8 * pgf1_28[k]
                  + f_4 * pc_y[k] * pgg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pb_y, pc_y, pc_z, pdh0_0, pdh1_0, pfh0_21, \
                         pfg_14, pfh1_21, pgf0_29, pgf1_29, pgg_43, \
                         pgg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_5 * pgf0_29[k]
                  - f_6 * pgf1_29[k]
                  + f_4 * pc_y[k] * pgg_43[k];

        t_61[k] = f_4 * pc_y[k] * pgg_44[k];

        t_62[k] = f_0 * pfg_14[k]
                  + f_2 * pgf0_29[k]
                  - f_3 * pgf1_29[k]
                  + f_4 * pc_z[k] * pgg_44[k];

        t_63[k] = f_12 * pdh0_0[k]
                  - f_13 * pdh1_0[k]
                  + pb_y[k] * pfh0_21[k]
                  - f_9 * pc_y[k] * pfh1_21[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pc_y, pc_z, pfg_15, pfg_16, pfg_17, \
                         pgf0_30, pgf1_30, pgg_45, pgg_46, pgg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * pfg_15[k]
                  + f_4 * pc_y[k] * pgg_45[k];

        t_65[k] = f_4 * pc_z[k] * pgg_45[k];

        t_66[k] = f_10 * pfg_16[k]
                  + f_5 * pgf0_30[k]
                  - f_6 * pgf1_30[k]
                  + f_4 * pc_y[k] * pgg_46[k];

        t_67[k] = f_10 * pfg_17[k]
                  + f_4 * pc_y[k] * pgg_47[k];

        t_68[k] = f_5 * pgf0_30[k]
                  - f_6 * pgf1_30[k]
                  + f_4 * pc_z[k] * pgg_47[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pc_y, pc_z, pfg_18, pfg_20, pgf0_31, pgf0_32, \
                         pgf1_31, pgf1_32, pgg_48, pgg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * pfg_18[k]
                  + f_7 * pgf0_31[k]
                  - f_8 * pgf1_31[k]
                  + f_4 * pc_y[k] * pgg_48[k];

        t_70[k] = f_4 * pc_z[k] * pgg_48[k];

        t_71[k] = f_10 * pfg_20[k]
                  + f_4 * pc_y[k] * pgg_50[k];

        t_72[k] = f_7 * pgf0_32[k]
                  - f_8 * pgf1_32[k]
                  + f_4 * pc_z[k] * pgg_50[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, pc_y, pc_z, sgg_55, sgg_57, pfg_24, \
                         pfg_55, pfg_57, pgg_51, pgg_54, pgg_55, \
                         pgg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_0 * sgg_55[k]
                  + f_10 * pfg_55[k]
                  + f_4 * pc_x[k] * pgg_55[k];

        t_74[k] = f_4 * pc_z[k] * pgg_51[k];

        t_75[k] = f_0 * sgg_57[k]
                  + f_10 * pfg_57[k]
                  + f_4 * pc_x[k] * pgg_57[k];

        t_76[k] = f_10 * pfg_24[k]
                  + f_4 * pc_y[k] * pgg_54[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_x, pc_y, pc_z, sgg_59, pfg_25, pfg_59, pgf0_36, \
                         pgf1_36, pgg_55, pgg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * sgg_59[k]
                  + f_10 * pfg_59[k]
                  + f_4 * pc_x[k] * pgg_59[k];

        t_78[k] = f_10 * pfg_25[k]
                  + f_2 * pgf0_36[k]
                  - f_3 * pgf1_36[k]
                  + f_4 * pc_y[k] * pgg_55[k];

        t_79[k] = f_4 * pc_z[k] * pgg_55[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, pfg_27, pfg_28, pfg_29, pgf0_38, \
                         pgf0_39, pgf1_38, pgf1_39, pgg_57, pgg_58, \
                         pgg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_10 * pfg_27[k]
                  + f_7 * pgf0_38[k]
                  - f_8 * pgf1_38[k]
                  + f_4 * pc_y[k] * pgg_57[k];

        t_81[k] = f_10 * pfg_28[k]
                  + f_5 * pgf0_39[k]
                  - f_6 * pgf1_39[k]
                  + f_4 * pc_y[k] * pgg_58[k];

        t_82[k] = f_10 * pfg_29[k]
                  + f_4 * pc_y[k] * pgg_59[k];

        t_83[k] = f_2 * pgf0_39[k]
                  - f_3 * pgf1_39[k]
                  + f_4 * pc_z[k] * pgg_59[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_y, pb_z, pc_y, pc_z, pfh0_24, pfh0_42, \
                         pfg_15, pfg_30, pfh1_24, pfh1_42, pgg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pb_y[k] * pfh0_42[k]
                  - f_9 * pc_y[k] * pfh1_42[k];

        t_85[k] = f_0 * pfg_30[k]
                  + f_4 * pc_y[k] * pgg_60[k];

        t_86[k] = f_0 * pfg_15[k]
                  + f_4 * pc_z[k] * pgg_60[k];

        t_87[k] = pb_z[k] * pfh0_24[k]
                  - f_9 * pc_z[k] * pfh1_24[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_y, pb_z, pc_y, pc_z, pfh0_27, pfh0_47, \
                         pfg_18, pfg_32, pfh1_27, pfh1_47, pgg_62, \
                         pgg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * pfg_32[k]
                  + f_4 * pc_y[k] * pgg_62[k];

        t_89[k] = pb_y[k] * pfh0_47[k]
                  - f_9 * pc_y[k] * pfh1_47[k];

        t_90[k] = pb_z[k] * pfh0_27[k]
                  - f_9 * pc_z[k] * pfh1_27[k];

        t_91[k] = f_0 * pfg_18[k]
                  + f_4 * pc_z[k] * pgg_63[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_y, pb_z, pc_y, pc_z, pfh0_31, pfh0_51, \
                         pfg_21, pfg_35, pfh1_31, pfh1_51, pgg_65, \
                         pgg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_0 * pfg_35[k]
                  + f_4 * pc_y[k] * pgg_65[k];

        t_93[k] = pb_y[k] * pfh0_51[k]
                  - f_9 * pc_y[k] * pfh1_51[k];

        t_94[k] = pb_z[k] * pfh0_31[k]
                  - f_9 * pc_z[k] * pfh1_31[k];

        t_95[k] = f_0 * pfg_21[k]
                  + f_4 * pc_z[k] * pgg_66[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, pc_x, pc_y, sgg_72, pfh0_56, pfg_39, pfg_72, \
                         pfh1_56, pgg_69, pgg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_0 * sgg_72[k]
                  + f_10 * pfg_72[k]
                  + f_4 * pc_x[k] * pgg_72[k];

        t_97[k] = f_0 * pfg_39[k]
                  + f_4 * pc_y[k] * pgg_69[k];

        t_98[k] = pb_y[k] * pfh0_56[k]
                  - f_9 * pc_y[k] * pfh1_56[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_z, pc_y, pc_z, pfh0_36, pfg_25, pfg_42, \
                         pfh1_36, pgf0_48, pgf1_48, pgg_70, pgg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pb_z[k] * pfh0_36[k]
                  - f_9 * pc_z[k] * pfh1_36[k];

        t_100[k] = f_0 * pfg_25[k]
                   + f_4 * pc_z[k] * pgg_70[k];

        t_101[k] = f_0 * pfg_42[k]
                   + f_7 * pgf0_48[k]
                   - f_8 * pgf1_48[k]
                   + f_4 * pc_y[k] * pgg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_y, pc_y, pfh0_62, pfg_43, pfg_44, pfh1_62, \
                         pgf0_49, pgf1_49, pgg_73, pgg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_0 * pfg_43[k]
                   + f_5 * pgf0_49[k]
                   - f_6 * pgf1_49[k]
                   + f_4 * pc_y[k] * pgg_73[k];

        t_103[k] = f_0 * pfg_44[k]
                   + f_4 * pc_y[k] * pgg_74[k];

        t_104[k] = pb_y[k] * pfh0_62[k]
                   - f_9 * pc_y[k] * pfh1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_z, pc_y, pc_z, pdh0_0, pdh1_0, \
                         pfh0_42, pfg_30, pfh1_42, pgf0_50, pgf1_50, pgg_75, \
                         pgg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * pdh0_0[k]
                   - f_13 * pdh1_0[k]
                   + pb_z[k] * pfh0_42[k]
                   - f_9 * pc_z[k] * pfh1_42[k];

        t_106[k] = f_4 * pc_y[k] * pgg_75[k];

        t_107[k] = f_10 * pfg_30[k]
                   + f_4 * pc_z[k] * pgg_75[k];

        t_108[k] = f_5 * pgf0_50[k]
                   - f_6 * pgf1_50[k]
                   + f_4 * pc_y[k] * pgg_76[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pc_y, pc_z, pfg_32, pfg_33, \
                         pgf0_50, pgf0_51, pgf1_50, pgf1_51, pgg_77, pgg_78, \
                         pgg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_4 * pc_y[k] * pgg_77[k];

        t_110[k] = f_10 * pfg_32[k]
                   + f_5 * pgf0_50[k]
                   - f_6 * pgf1_50[k]
                   + f_4 * pc_z[k] * pgg_77[k];

        t_111[k] = f_7 * pgf0_51[k]
                   - f_8 * pgf1_51[k]
                   + f_4 * pc_y[k] * pgg_78[k];

        t_112[k] = f_10 * pfg_33[k]
                   + f_4 * pc_z[k] * pgg_78[k];

        t_113[k] = f_4 * pc_y[k] * pgg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pc_x, pc_z, sgg_85, pfg_35, pfg_36, pfg_85, \
                         pgf0_52, pgf1_52, pgg_80, pgg_81, pgg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_10 * pfg_35[k]
                   + f_7 * pgf0_52[k]
                   - f_8 * pgf1_52[k]
                   + f_4 * pc_z[k] * pgg_80[k];

        t_115[k] = f_0 * sgg_85[k]
                   + f_10 * pfg_85[k]
                   + f_4 * pc_x[k] * pgg_85[k];

        t_116[k] = f_10 * pfg_36[k]
                   + f_4 * pc_z[k] * pgg_81[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pc_x, pc_y, sgg_87, sgg_89, pfg_87, \
                         pfg_89, pgf0_56, pgf1_56, pgg_84, pgg_85, pgg_87, \
                         pgg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * sgg_87[k]
                   + f_10 * pfg_87[k]
                   + f_4 * pc_x[k] * pgg_87[k];

        t_118[k] = f_4 * pc_y[k] * pgg_84[k];

        t_119[k] = f_0 * sgg_89[k]
                   + f_10 * pfg_89[k]
                   + f_4 * pc_x[k] * pgg_89[k];

        t_120[k] = f_2 * pgf0_56[k]
                   - f_3 * pgf1_56[k]
                   + f_4 * pc_y[k] * pgg_85[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_y, pc_z, pfg_40, pgf0_58, pgf0_59, \
                         pgf1_58, pgf1_59, pgg_85, pgg_87, pgg_88, \
                         pgg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * pfg_40[k]
                   + f_4 * pc_z[k] * pgg_85[k];

        t_122[k] = f_7 * pgf0_58[k]
                   - f_8 * pgf1_58[k]
                   + f_4 * pc_y[k] * pgg_87[k];

        t_123[k] = f_5 * pgf0_59[k]
                   - f_6 * pgf1_59[k]
                   + f_4 * pc_y[k] * pgg_88[k];

        t_124[k] = f_4 * pc_y[k] * pgg_89[k];
    }
}

static auto
compute_prim_pgh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgh0, const size_t sgg,
                                                          const size_t sgh1, const size_t pfh0,
                                                          const size_t pfg, const size_t pfh1,
                                                          const size_t pgf0, const size_t pgf1,
                                                          const size_t pgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_210 = buffer.data(sgh0 + 210);
    const auto *sgh0_213 = buffer.data(sgh0 + 213);
    const auto *sgh0_216 = buffer.data(sgh0 + 216);
    const auto *sgh0_225 = buffer.data(sgh0 + 225);
    const auto *sgh0_227 = buffer.data(sgh0 + 227);
    const auto *sgh0_228 = buffer.data(sgh0 + 228);
    const auto *sgh0_230 = buffer.data(sgh0 + 230);
    const auto *sgh0_236 = buffer.data(sgh0 + 236);
    const auto *sgh0_240 = buffer.data(sgh0 + 240);

    const auto *sgg_90 = buffer.data(sgg + 90);
    const auto *sgg_100 = buffer.data(sgg + 100);
    const auto *sgg_102 = buffer.data(sgg + 102);
    const auto *sgg_104 = buffer.data(sgg + 104);
    const auto *sgg_117 = buffer.data(sgg + 117);
    const auto *sgg_119 = buffer.data(sgg + 119);
    const auto *sgg_130 = buffer.data(sgg + 130);
    const auto *sgg_132 = buffer.data(sgg + 132);
    const auto *sgg_135 = buffer.data(sgg + 135);
    const auto *sgg_145 = buffer.data(sgg + 145);
    const auto *sgg_147 = buffer.data(sgg + 147);
    const auto *sgg_149 = buffer.data(sgg + 149);
    const auto *sgg_150 = buffer.data(sgg + 150);
    const auto *sgg_153 = buffer.data(sgg + 153);
    const auto *sgg_156 = buffer.data(sgg + 156);
    const auto *sgg_160 = buffer.data(sgg + 160);
    const auto *sgg_162 = buffer.data(sgg + 162);
    const auto *sgg_164 = buffer.data(sgg + 164);
    const auto *sgg_170 = buffer.data(sgg + 170);
    const auto *sgg_174 = buffer.data(sgg + 174);

    const auto *sgh1_210 = buffer.data(sgh1 + 210);
    const auto *sgh1_213 = buffer.data(sgh1 + 213);
    const auto *sgh1_216 = buffer.data(sgh1 + 216);
    const auto *sgh1_225 = buffer.data(sgh1 + 225);
    const auto *sgh1_227 = buffer.data(sgh1 + 227);
    const auto *sgh1_228 = buffer.data(sgh1 + 228);
    const auto *sgh1_230 = buffer.data(sgh1 + 230);
    const auto *sgh1_236 = buffer.data(sgh1 + 236);
    const auto *sgh1_240 = buffer.data(sgh1 + 240);

    const auto *pfh0_63 = buffer.data(pfh0 + 63);
    const auto *pfh0_66 = buffer.data(pfh0 + 66);
    const auto *pfh0_69 = buffer.data(pfh0 + 69);
    const auto *pfh0_73 = buffer.data(pfh0 + 73);
    const auto *pfh0_105 = buffer.data(pfh0 + 105);
    const auto *pfh0_110 = buffer.data(pfh0 + 110);
    const auto *pfh0_114 = buffer.data(pfh0 + 114);
    const auto *pfh0_119 = buffer.data(pfh0 + 119);
    const auto *pfh0_126 = buffer.data(pfh0 + 126);
    const auto *pfh0_129 = buffer.data(pfh0 + 129);
    const auto *pfh0_132 = buffer.data(pfh0 + 132);

    const auto *pfg_44 = buffer.data(pfg + 44);
    const auto *pfg_45 = buffer.data(pfg + 45);
    const auto *pfg_46 = buffer.data(pfg + 46);
    const auto *pfg_47 = buffer.data(pfg + 47);
    const auto *pfg_48 = buffer.data(pfg + 48);
    const auto *pfg_50 = buffer.data(pfg + 50);
    const auto *pfg_51 = buffer.data(pfg + 51);
    const auto *pfg_54 = buffer.data(pfg + 54);
    const auto *pfg_55 = buffer.data(pfg + 55);
    const auto *pfg_57 = buffer.data(pfg + 57);
    const auto *pfg_58 = buffer.data(pfg + 58);
    const auto *pfg_59 = buffer.data(pfg + 59);
    const auto *pfg_60 = buffer.data(pfg + 60);
    const auto *pfg_62 = buffer.data(pfg + 62);
    const auto *pfg_63 = buffer.data(pfg + 63);
    const auto *pfg_65 = buffer.data(pfg + 65);
    const auto *pfg_66 = buffer.data(pfg + 66);
    const auto *pfg_69 = buffer.data(pfg + 69);
    const auto *pfg_70 = buffer.data(pfg + 70);
    const auto *pfg_72 = buffer.data(pfg + 72);
    const auto *pfg_73 = buffer.data(pfg + 73);
    const auto *pfg_74 = buffer.data(pfg + 74);
    const auto *pfg_75 = buffer.data(pfg + 75);
    const auto *pfg_76 = buffer.data(pfg + 76);
    const auto *pfg_77 = buffer.data(pfg + 77);
    const auto *pfg_78 = buffer.data(pfg + 78);
    const auto *pfg_80 = buffer.data(pfg + 80);
    const auto *pfg_81 = buffer.data(pfg + 81);
    const auto *pfg_84 = buffer.data(pfg + 84);
    const auto *pfg_85 = buffer.data(pfg + 85);
    const auto *pfg_87 = buffer.data(pfg + 87);
    const auto *pfg_88 = buffer.data(pfg + 88);
    const auto *pfg_89 = buffer.data(pfg + 89);
    const auto *pfg_90 = buffer.data(pfg + 90);
    const auto *pfg_92 = buffer.data(pfg + 92);
    const auto *pfg_93 = buffer.data(pfg + 93);
    const auto *pfg_95 = buffer.data(pfg + 95);
    const auto *pfg_99 = buffer.data(pfg + 99);
    const auto *pfg_100 = buffer.data(pfg + 100);
    const auto *pfg_102 = buffer.data(pfg + 102);
    const auto *pfg_104 = buffer.data(pfg + 104);
    const auto *pfg_105 = buffer.data(pfg + 105);
    const auto *pfg_107 = buffer.data(pfg + 107);
    const auto *pfg_110 = buffer.data(pfg + 110);
    const auto *pfg_117 = buffer.data(pfg + 117);
    const auto *pfg_119 = buffer.data(pfg + 119);
    const auto *pfg_130 = buffer.data(pfg + 130);
    const auto *pfg_132 = buffer.data(pfg + 132);
    const auto *pfg_135 = buffer.data(pfg + 135);
    const auto *pfg_145 = buffer.data(pfg + 145);
    const auto *pfg_147 = buffer.data(pfg + 147);
    const auto *pfg_149 = buffer.data(pfg + 149);

    const auto *pfh1_63 = buffer.data(pfh1 + 63);
    const auto *pfh1_66 = buffer.data(pfh1 + 66);
    const auto *pfh1_69 = buffer.data(pfh1 + 69);
    const auto *pfh1_73 = buffer.data(pfh1 + 73);
    const auto *pfh1_105 = buffer.data(pfh1 + 105);
    const auto *pfh1_110 = buffer.data(pfh1 + 110);
    const auto *pfh1_114 = buffer.data(pfh1 + 114);
    const auto *pfh1_119 = buffer.data(pfh1 + 119);
    const auto *pfh1_126 = buffer.data(pfh1 + 126);
    const auto *pfh1_129 = buffer.data(pfh1 + 129);
    const auto *pfh1_132 = buffer.data(pfh1 + 132);

    const auto *pgf0_59 = buffer.data(pgf0 + 59);
    const auto *pgf0_60 = buffer.data(pgf0 + 60);
    const auto *pgf0_61 = buffer.data(pgf0 + 61);
    const auto *pgf0_62 = buffer.data(pgf0 + 62);
    const auto *pgf0_66 = buffer.data(pgf0 + 66);
    const auto *pgf0_68 = buffer.data(pgf0 + 68);
    const auto *pgf0_69 = buffer.data(pgf0 + 69);
    const auto *pgf0_70 = buffer.data(pgf0 + 70);
    const auto *pgf0_72 = buffer.data(pgf0 + 72);
    const auto *pgf0_76 = buffer.data(pgf0 + 76);
    const auto *pgf0_78 = buffer.data(pgf0 + 78);
    const auto *pgf0_79 = buffer.data(pgf0 + 79);
    const auto *pgf0_80 = buffer.data(pgf0 + 80);
    const auto *pgf0_81 = buffer.data(pgf0 + 81);
    const auto *pgf0_86 = buffer.data(pgf0 + 86);
    const auto *pgf0_88 = buffer.data(pgf0 + 88);
    const auto *pgf0_89 = buffer.data(pgf0 + 89);
    const auto *pgf0_90 = buffer.data(pgf0 + 90);
    const auto *pgf0_91 = buffer.data(pgf0 + 91);
    const auto *pgf0_92 = buffer.data(pgf0 + 92);
    const auto *pgf0_96 = buffer.data(pgf0 + 96);
    const auto *pgf0_98 = buffer.data(pgf0 + 98);
    const auto *pgf0_99 = buffer.data(pgf0 + 99);
    const auto *pgf0_100 = buffer.data(pgf0 + 100);
    const auto *pgf0_102 = buffer.data(pgf0 + 102);

    const auto *pgf1_59 = buffer.data(pgf1 + 59);
    const auto *pgf1_60 = buffer.data(pgf1 + 60);
    const auto *pgf1_61 = buffer.data(pgf1 + 61);
    const auto *pgf1_62 = buffer.data(pgf1 + 62);
    const auto *pgf1_66 = buffer.data(pgf1 + 66);
    const auto *pgf1_68 = buffer.data(pgf1 + 68);
    const auto *pgf1_69 = buffer.data(pgf1 + 69);
    const auto *pgf1_70 = buffer.data(pgf1 + 70);
    const auto *pgf1_72 = buffer.data(pgf1 + 72);
    const auto *pgf1_76 = buffer.data(pgf1 + 76);
    const auto *pgf1_78 = buffer.data(pgf1 + 78);
    const auto *pgf1_79 = buffer.data(pgf1 + 79);
    const auto *pgf1_80 = buffer.data(pgf1 + 80);
    const auto *pgf1_81 = buffer.data(pgf1 + 81);
    const auto *pgf1_86 = buffer.data(pgf1 + 86);
    const auto *pgf1_88 = buffer.data(pgf1 + 88);
    const auto *pgf1_89 = buffer.data(pgf1 + 89);
    const auto *pgf1_90 = buffer.data(pgf1 + 90);
    const auto *pgf1_91 = buffer.data(pgf1 + 91);
    const auto *pgf1_92 = buffer.data(pgf1 + 92);
    const auto *pgf1_96 = buffer.data(pgf1 + 96);
    const auto *pgf1_98 = buffer.data(pgf1 + 98);
    const auto *pgf1_99 = buffer.data(pgf1 + 99);
    const auto *pgf1_100 = buffer.data(pgf1 + 100);
    const auto *pgf1_102 = buffer.data(pgf1 + 102);

    const auto *pgg_89 = buffer.data(pgg + 89);
    const auto *pgg_90 = buffer.data(pgg + 90);
    const auto *pgg_91 = buffer.data(pgg + 91);
    const auto *pgg_92 = buffer.data(pgg + 92);
    const auto *pgg_93 = buffer.data(pgg + 93);
    const auto *pgg_95 = buffer.data(pgg + 95);
    const auto *pgg_96 = buffer.data(pgg + 96);
    const auto *pgg_99 = buffer.data(pgg + 99);
    const auto *pgg_100 = buffer.data(pgg + 100);
    const auto *pgg_102 = buffer.data(pgg + 102);
    const auto *pgg_103 = buffer.data(pgg + 103);
    const auto *pgg_104 = buffer.data(pgg + 104);
    const auto *pgg_105 = buffer.data(pgg + 105);
    const auto *pgg_107 = buffer.data(pgg + 107);
    const auto *pgg_108 = buffer.data(pgg + 108);
    const auto *pgg_110 = buffer.data(pgg + 110);
    const auto *pgg_111 = buffer.data(pgg + 111);
    const auto *pgg_114 = buffer.data(pgg + 114);
    const auto *pgg_115 = buffer.data(pgg + 115);
    const auto *pgg_117 = buffer.data(pgg + 117);
    const auto *pgg_118 = buffer.data(pgg + 118);
    const auto *pgg_119 = buffer.data(pgg + 119);
    const auto *pgg_120 = buffer.data(pgg + 120);
    const auto *pgg_121 = buffer.data(pgg + 121);
    const auto *pgg_122 = buffer.data(pgg + 122);
    const auto *pgg_123 = buffer.data(pgg + 123);
    const auto *pgg_125 = buffer.data(pgg + 125);
    const auto *pgg_126 = buffer.data(pgg + 126);
    const auto *pgg_129 = buffer.data(pgg + 129);
    const auto *pgg_130 = buffer.data(pgg + 130);
    const auto *pgg_132 = buffer.data(pgg + 132);
    const auto *pgg_133 = buffer.data(pgg + 133);
    const auto *pgg_134 = buffer.data(pgg + 134);
    const auto *pgg_135 = buffer.data(pgg + 135);
    const auto *pgg_136 = buffer.data(pgg + 136);
    const auto *pgg_137 = buffer.data(pgg + 137);
    const auto *pgg_138 = buffer.data(pgg + 138);
    const auto *pgg_140 = buffer.data(pgg + 140);
    const auto *pgg_141 = buffer.data(pgg + 141);
    const auto *pgg_144 = buffer.data(pgg + 144);
    const auto *pgg_145 = buffer.data(pgg + 145);
    const auto *pgg_147 = buffer.data(pgg + 147);
    const auto *pgg_148 = buffer.data(pgg + 148);
    const auto *pgg_149 = buffer.data(pgg + 149);
    const auto *pgg_150 = buffer.data(pgg + 150);
    const auto *pgg_152 = buffer.data(pgg + 152);
    const auto *pgg_153 = buffer.data(pgg + 153);
    const auto *pgg_155 = buffer.data(pgg + 155);
    const auto *pgg_156 = buffer.data(pgg + 156);
    const auto *pgg_159 = buffer.data(pgg + 159);
    const auto *pgg_160 = buffer.data(pgg + 160);
    const auto *pgg_162 = buffer.data(pgg + 162);
    const auto *pgg_164 = buffer.data(pgg + 164);
    const auto *pgg_165 = buffer.data(pgg + 165);
    const auto *pgg_167 = buffer.data(pgg + 167);
    const auto *pgg_168 = buffer.data(pgg + 168);
    const auto *pgg_170 = buffer.data(pgg + 170);

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_y, pc_z, sgg_90, pfg_44, pfg_45, \
                         pfg_90, pgf0_59, pgf0_60, pgf1_59, pgf1_60, pgg_89, \
                         pgg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_10 * pfg_44[k]
                   + f_2 * pgf0_59[k]
                   - f_3 * pgf1_59[k]
                   + f_4 * pc_z[k] * pgg_89[k];

        t_126[k] = f_0 * sgg_90[k]
                   + f_0 * pfg_90[k]
                   + f_2 * pgf0_60[k]
                   - f_3 * pgf1_60[k]
                   + f_4 * pc_x[k] * pgg_90[k];

        t_127[k] = f_11 * pfg_45[k]
                   + f_4 * pc_y[k] * pgg_90[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_y, pc_z, pfg_46, pfg_47, pgf0_60, \
                         pgf1_60, pgg_90, pgg_91, pgg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_4 * pc_z[k] * pgg_90[k];

        t_129[k] = f_11 * pfg_46[k]
                   + f_5 * pgf0_60[k]
                   - f_6 * pgf1_60[k]
                   + f_4 * pc_y[k] * pgg_91[k];

        t_130[k] = f_11 * pfg_47[k]
                   + f_4 * pc_y[k] * pgg_92[k];

        t_131[k] = f_5 * pgf0_60[k]
                   - f_6 * pgf1_60[k]
                   + f_4 * pc_z[k] * pgg_92[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_y, pc_z, pfg_48, pfg_50, pgf0_61, \
                         pgf0_62, pgf1_61, pgf1_62, pgg_93, pgg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_11 * pfg_48[k]
                   + f_7 * pgf0_61[k]
                   - f_8 * pgf1_61[k]
                   + f_4 * pc_y[k] * pgg_93[k];

        t_133[k] = f_4 * pc_z[k] * pgg_93[k];

        t_134[k] = f_11 * pfg_50[k]
                   + f_4 * pc_y[k] * pgg_95[k];

        t_135[k] = f_7 * pgf0_62[k]
                   - f_8 * pgf1_62[k]
                   + f_4 * pc_z[k] * pgg_95[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pc_x, pc_y, pc_z, sgg_100, sgg_102, \
                         pfg_54, pfg_100, pfg_102, pgg_96, pgg_99, pgg_100, \
                         pgg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_0 * sgg_100[k]
                   + f_0 * pfg_100[k]
                   + f_4 * pc_x[k] * pgg_100[k];

        t_137[k] = f_4 * pc_z[k] * pgg_96[k];

        t_138[k] = f_0 * sgg_102[k]
                   + f_0 * pfg_102[k]
                   + f_4 * pc_x[k] * pgg_102[k];

        t_139[k] = f_11 * pfg_54[k]
                   + f_4 * pc_y[k] * pgg_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pc_x, pc_y, pc_z, sgg_104, pfg_55, pfg_104, \
                         pgf0_66, pgf1_66, pgg_100, pgg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * sgg_104[k]
                   + f_0 * pfg_104[k]
                   + f_4 * pc_x[k] * pgg_104[k];

        t_141[k] = f_11 * pfg_55[k]
                   + f_2 * pgf0_66[k]
                   - f_3 * pgf1_66[k]
                   + f_4 * pc_y[k] * pgg_100[k];

        t_142[k] = f_4 * pc_z[k] * pgg_100[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_y, pc_z, pfg_57, pfg_58, pfg_59, \
                         pgf0_68, pgf0_69, pgf1_68, pgf1_69, pgg_102, pgg_103, \
                         pgg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_11 * pfg_57[k]
                   + f_7 * pgf0_68[k]
                   - f_8 * pgf1_68[k]
                   + f_4 * pc_y[k] * pgg_102[k];

        t_144[k] = f_11 * pfg_58[k]
                   + f_5 * pgf0_69[k]
                   - f_6 * pgf1_69[k]
                   + f_4 * pc_y[k] * pgg_103[k];

        t_145[k] = f_11 * pfg_59[k]
                   + f_4 * pc_y[k] * pgg_104[k];

        t_146[k] = f_2 * pgf0_69[k]
                   - f_3 * pgf1_69[k]
                   + f_4 * pc_z[k] * pgg_104[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pb_z, pc_y, pc_z, pfh0_63, pfh0_66, \
                         pfg_45, pfg_60, pfh1_63, pfh1_66, pgg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pb_z[k] * pfh0_63[k]
                   - f_9 * pc_z[k] * pfh1_63[k];

        t_148[k] = f_10 * pfg_60[k]
                   + f_4 * pc_y[k] * pgg_105[k];

        t_149[k] = f_0 * pfg_45[k]
                   + f_4 * pc_z[k] * pgg_105[k];

        t_150[k] = pb_z[k] * pfh0_66[k]
                   - f_9 * pc_z[k] * pfh1_66[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pb_z, pc_y, pc_z, pfh0_69, pfg_47, \
                         pfg_48, pfg_62, pfh1_69, pgf0_70, pgf1_70, pgg_107, \
                         pgg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_10 * pfg_62[k]
                   + f_4 * pc_y[k] * pgg_107[k];

        t_152[k] = f_0 * pfg_47[k]
                   + f_5 * pgf0_70[k]
                   - f_6 * pgf1_70[k]
                   + f_4 * pc_z[k] * pgg_107[k];

        t_153[k] = pb_z[k] * pfh0_69[k]
                   - f_9 * pc_z[k] * pfh1_69[k];

        t_154[k] = f_0 * pfg_48[k]
                   + f_4 * pc_z[k] * pgg_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pb_z, pc_y, pc_z, pfh0_73, pfg_50, \
                         pfg_51, pfg_65, pfh1_73, pgf0_72, pgf1_72, pgg_110, \
                         pgg_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_10 * pfg_65[k]
                   + f_4 * pc_y[k] * pgg_110[k];

        t_156[k] = f_0 * pfg_50[k]
                   + f_7 * pgf0_72[k]
                   - f_8 * pgf1_72[k]
                   + f_4 * pc_z[k] * pgg_110[k];

        t_157[k] = pb_z[k] * pfh0_73[k]
                   - f_9 * pc_z[k] * pfh1_73[k];

        t_158[k] = f_0 * pfg_51[k]
                   + f_4 * pc_z[k] * pgg_111[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, sgg_117, sgg_119, pfg_69, pfg_117, \
                         pfg_119, pgg_114, pgg_117, pgg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_0 * sgg_117[k]
                   + f_0 * pfg_117[k]
                   + f_4 * pc_x[k] * pgg_117[k];

        t_160[k] = f_10 * pfg_69[k]
                   + f_4 * pc_y[k] * pgg_114[k];

        t_161[k] = f_0 * sgg_119[k]
                   + f_0 * pfg_119[k]
                   + f_4 * pc_x[k] * pgg_119[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pc_y, pc_z, pfg_55, pfg_70, pfg_72, pgf0_76, \
                         pgf0_78, pgf1_76, pgf1_78, pgg_115, pgg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_10 * pfg_70[k]
                   + f_2 * pgf0_76[k]
                   - f_3 * pgf1_76[k]
                   + f_4 * pc_y[k] * pgg_115[k];

        t_163[k] = f_0 * pfg_55[k]
                   + f_4 * pc_z[k] * pgg_115[k];

        t_164[k] = f_10 * pfg_72[k]
                   + f_7 * pgf0_78[k]
                   - f_8 * pgf1_78[k]
                   + f_4 * pc_y[k] * pgg_117[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pb_y, pc_y, pc_z, pfh0_105, pfg_59, \
                         pfg_73, pfg_74, pfh1_105, pgf0_79, pgf1_79, pgg_118, \
                         pgg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_10 * pfg_73[k]
                   + f_5 * pgf0_79[k]
                   - f_6 * pgf1_79[k]
                   + f_4 * pc_y[k] * pgg_118[k];

        t_166[k] = f_10 * pfg_74[k]
                   + f_4 * pc_y[k] * pgg_119[k];

        t_167[k] = f_0 * pfg_59[k]
                   + f_2 * pgf0_79[k]
                   - f_3 * pgf1_79[k]
                   + f_4 * pc_z[k] * pgg_119[k];

        t_168[k] = pb_y[k] * pfh0_105[k]
                   - f_9 * pc_y[k] * pfh1_105[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pc_y, pc_z, pfg_60, pfg_75, pfg_76, \
                         pfg_77, pgf0_80, pgf1_80, pgg_120, pgg_121, \
                         pgg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_0 * pfg_75[k]
                   + f_4 * pc_y[k] * pgg_120[k];

        t_170[k] = f_10 * pfg_60[k]
                   + f_4 * pc_z[k] * pgg_120[k];

        t_171[k] = f_0 * pfg_76[k]
                   + f_5 * pgf0_80[k]
                   - f_6 * pgf1_80[k]
                   + f_4 * pc_y[k] * pgg_121[k];

        t_172[k] = f_0 * pfg_77[k]
                   + f_4 * pc_y[k] * pgg_122[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pb_y, pc_y, pc_z, pfh0_110, pfg_63, \
                         pfg_78, pfg_80, pfh1_110, pgf0_81, pgf1_81, pgg_123, \
                         pgg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = pb_y[k] * pfh0_110[k]
                   - f_9 * pc_y[k] * pfh1_110[k];

        t_174[k] = f_0 * pfg_78[k]
                   + f_7 * pgf0_81[k]
                   - f_8 * pgf1_81[k]
                   + f_4 * pc_y[k] * pgg_123[k];

        t_175[k] = f_10 * pfg_63[k]
                   + f_4 * pc_z[k] * pgg_123[k];

        t_176[k] = f_0 * pfg_80[k]
                   + f_4 * pc_y[k] * pgg_125[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_y, pc_x, pc_y, pc_z, sgg_130, pfh0_114, \
                         pfg_66, pfg_130, pfh1_114, pgg_126, pgg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pb_y[k] * pfh0_114[k]
                   - f_9 * pc_y[k] * pfh1_114[k];

        t_178[k] = f_0 * sgg_130[k]
                   + f_0 * pfg_130[k]
                   + f_4 * pc_x[k] * pgg_130[k];

        t_179[k] = f_10 * pfg_66[k]
                   + f_4 * pc_z[k] * pgg_126[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pb_y, pc_x, pc_y, sgg_132, pfh0_119, pfg_84, \
                         pfg_132, pfh1_119, pgg_129, pgg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_0 * sgg_132[k]
                   + f_0 * pfg_132[k]
                   + f_4 * pc_x[k] * pgg_132[k];

        t_181[k] = f_0 * pfg_84[k]
                   + f_4 * pc_y[k] * pgg_129[k];

        t_182[k] = pb_y[k] * pfh0_119[k]
                   - f_9 * pc_y[k] * pfh1_119[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, pfg_70, pfg_85, pfg_87, pgf0_86, \
                         pgf0_88, pgf1_86, pgf1_88, pgg_130, pgg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * pfg_85[k]
                   + f_2 * pgf0_86[k]
                   - f_3 * pgf1_86[k]
                   + f_4 * pc_y[k] * pgg_130[k];

        t_184[k] = f_10 * pfg_70[k]
                   + f_4 * pc_z[k] * pgg_130[k];

        t_185[k] = f_0 * pfg_87[k]
                   + f_7 * pgf0_88[k]
                   - f_8 * pgf1_88[k]
                   + f_4 * pc_y[k] * pgg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pc_y, pc_z, pfg_74, pfg_88, pfg_89, pgf0_89, \
                         pgf1_89, pgg_133, pgg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_0 * pfg_88[k]
                   + f_5 * pgf0_89[k]
                   - f_6 * pgf1_89[k]
                   + f_4 * pc_y[k] * pgg_133[k];

        t_187[k] = f_0 * pfg_89[k]
                   + f_4 * pc_y[k] * pgg_134[k];

        t_188[k] = f_10 * pfg_74[k]
                   + f_2 * pgf0_89[k]
                   - f_3 * pgf1_89[k]
                   + f_4 * pc_z[k] * pgg_134[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, sgg_135, pfg_75, \
                         pfg_135, pgf0_90, pgf1_90, pgg_135, pgg_136, \
                         pgg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_0 * sgg_135[k]
                   + f_0 * pfg_135[k]
                   + f_2 * pgf0_90[k]
                   - f_3 * pgf1_90[k]
                   + f_4 * pc_x[k] * pgg_135[k];

        t_190[k] = f_4 * pc_y[k] * pgg_135[k];

        t_191[k] = f_11 * pfg_75[k]
                   + f_4 * pc_z[k] * pgg_135[k];

        t_192[k] = f_5 * pgf0_90[k]
                   - f_6 * pgf1_90[k]
                   + f_4 * pc_y[k] * pgg_136[k];

        t_193[k] = f_4 * pc_y[k] * pgg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_y, pc_z, pfg_77, pfg_78, pgf0_90, \
                         pgf0_91, pgf1_90, pgf1_91, pgg_137, pgg_138, \
                         pgg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_11 * pfg_77[k]
                   + f_5 * pgf0_90[k]
                   - f_6 * pgf1_90[k]
                   + f_4 * pc_z[k] * pgg_137[k];

        t_195[k] = f_7 * pgf0_91[k]
                   - f_8 * pgf1_91[k]
                   + f_4 * pc_y[k] * pgg_138[k];

        t_196[k] = f_11 * pfg_78[k]
                   + f_4 * pc_z[k] * pgg_138[k];

        t_197[k] = f_4 * pc_y[k] * pgg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pc_x, pc_z, sgg_145, pfg_80, pfg_81, pfg_145, \
                         pgf0_92, pgf1_92, pgg_140, pgg_141, pgg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_11 * pfg_80[k]
                   + f_7 * pgf0_92[k]
                   - f_8 * pgf1_92[k]
                   + f_4 * pc_z[k] * pgg_140[k];

        t_199[k] = f_0 * sgg_145[k]
                   + f_0 * pfg_145[k]
                   + f_4 * pc_x[k] * pgg_145[k];

        t_200[k] = f_11 * pfg_81[k]
                   + f_4 * pc_z[k] * pgg_141[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pc_x, pc_y, sgg_147, sgg_149, pfg_147, \
                         pfg_149, pgf0_96, pgf1_96, pgg_144, pgg_145, pgg_147, \
                         pgg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_0 * sgg_147[k]
                   + f_0 * pfg_147[k]
                   + f_4 * pc_x[k] * pgg_147[k];

        t_202[k] = f_4 * pc_y[k] * pgg_144[k];

        t_203[k] = f_0 * sgg_149[k]
                   + f_0 * pfg_149[k]
                   + f_4 * pc_x[k] * pgg_149[k];

        t_204[k] = f_2 * pgf0_96[k]
                   - f_3 * pgf1_96[k]
                   + f_4 * pc_y[k] * pgg_145[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_y, pc_z, pfg_85, pgf0_98, pgf0_99, \
                         pgf1_98, pgf1_99, pgg_145, pgg_147, pgg_148, \
                         pgg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_11 * pfg_85[k]
                   + f_4 * pc_z[k] * pgg_145[k];

        t_206[k] = f_7 * pgf0_98[k]
                   - f_8 * pgf1_98[k]
                   + f_4 * pc_y[k] * pgg_147[k];

        t_207[k] = f_5 * pgf0_99[k]
                   - f_6 * pgf1_99[k]
                   + f_4 * pc_y[k] * pgg_148[k];

        t_208[k] = f_4 * pc_y[k] * pgg_149[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pa_x, pc_x, pc_y, pc_z, sgh0_210, sgg_150, \
                         sgh1_210, pfg_89, pfg_90, pgf0_99, pgf1_99, pgg_149, \
                         pgg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_11 * pfg_89[k]
                   + f_2 * pgf0_99[k]
                   - f_3 * pgf1_99[k]
                   + f_4 * pc_z[k] * pgg_149[k];

        t_210[k] = pa_x[k] * sgh0_210[k]
                   + f_14 * sgg_150[k]
                   - f_9 * pc_x[k] * sgh1_210[k];

        t_211[k] = f_1 * pfg_90[k]
                   + f_4 * pc_y[k] * pgg_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_x, pc_x, pc_y, pc_z, sgh0_213, \
                         sgg_153, sgh1_213, pfg_92, pgf0_100, pgf1_100, pgg_150, \
                         pgg_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_4 * pc_z[k] * pgg_150[k];

        t_213[k] = pa_x[k] * sgh0_213[k]
                   + f_11 * sgg_153[k]
                   - f_9 * pc_x[k] * sgh1_213[k];

        t_214[k] = f_1 * pfg_92[k]
                   + f_4 * pc_y[k] * pgg_152[k];

        t_215[k] = f_5 * pgf0_100[k]
                   - f_6 * pgf1_100[k]
                   + f_4 * pc_z[k] * pgg_152[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pa_x, pc_x, pc_y, pc_z, sgh0_216, \
                         sgg_156, sgh1_216, pfg_95, pgf0_102, pgf1_102, pgg_153, \
                         pgg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pa_x[k] * sgh0_216[k]
                   + f_10 * sgg_156[k]
                   - f_9 * pc_x[k] * sgh1_216[k];

        t_217[k] = f_4 * pc_z[k] * pgg_153[k];

        t_218[k] = f_1 * pfg_95[k]
                   + f_4 * pc_y[k] * pgg_155[k];

        t_219[k] = f_7 * pgf0_102[k]
                   - f_8 * pgf1_102[k]
                   + f_4 * pc_z[k] * pgg_155[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pc_x, pc_y, pc_z, sgg_160, sgg_162, \
                         pfg_99, pgg_156, pgg_159, pgg_160, pgg_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_0 * sgg_160[k]
                   + f_4 * pc_x[k] * pgg_160[k];

        t_221[k] = f_4 * pc_z[k] * pgg_156[k];

        t_222[k] = f_0 * sgg_162[k]
                   + f_4 * pc_x[k] * pgg_162[k];

        t_223[k] = f_1 * pfg_99[k]
                   + f_4 * pc_y[k] * pgg_159[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_x, pc_x, pc_z, sgh0_225, sgh0_227, \
                         sgg_164, sgh1_225, sgh1_227, pgg_160, \
                         pgg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_0 * sgg_164[k]
                   + f_4 * pc_x[k] * pgg_164[k];

        t_225[k] = pa_x[k] * sgh0_225[k]
                   - f_9 * pc_x[k] * sgh1_225[k];

        t_226[k] = f_4 * pc_z[k] * pgg_160[k];

        t_227[k] = pa_x[k] * sgh0_227[k]
                   - f_9 * pc_x[k] * sgh1_227[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pa_x, pc_x, pc_y, sgh0_228, sgh0_230, sgh1_228, \
                         sgh1_230, pfg_104, pgg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = pa_x[k] * sgh0_228[k]
                   - f_9 * pc_x[k] * sgh1_228[k];

        t_229[k] = f_1 * pfg_104[k]
                   + f_4 * pc_y[k] * pgg_164[k];

        t_230[k] = pa_x[k] * sgh0_230[k]
                   - f_9 * pc_x[k] * sgh1_230[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_z, pc_y, pc_z, pfh0_126, pfh0_129, \
                         pfg_90, pfg_105, pfh1_126, pfh1_129, pgg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = pb_z[k] * pfh0_126[k]
                   - f_9 * pc_z[k] * pfh1_126[k];

        t_232[k] = f_11 * pfg_105[k]
                   + f_4 * pc_y[k] * pgg_165[k];

        t_233[k] = f_0 * pfg_90[k]
                   + f_4 * pc_z[k] * pgg_165[k];

        t_234[k] = pb_z[k] * pfh0_129[k]
                   - f_9 * pc_z[k] * pfh1_129[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pa_x, pb_z, pc_x, pc_y, pc_z, sgh0_236, sgg_170, \
                         sgh1_236, pfh0_132, pfg_107, pfh1_132, \
                         pgg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_11 * pfg_107[k]
                   + f_4 * pc_y[k] * pgg_167[k];

        t_236[k] = pa_x[k] * sgh0_236[k]
                   + f_11 * sgg_170[k]
                   - f_9 * pc_x[k] * sgh1_236[k];

        t_237[k] = pb_z[k] * pfh0_132[k]
                   - f_9 * pc_z[k] * pfh1_132[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pa_x, pc_x, pc_y, pc_z, sgh0_240, sgg_174, \
                         sgh1_240, pfg_93, pfg_110, pgg_168, pgg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * pfg_93[k]
                   + f_4 * pc_z[k] * pgg_168[k];

        t_239[k] = f_11 * pfg_110[k]
                   + f_4 * pc_y[k] * pgg_170[k];

        t_240[k] = pa_x[k] * sgh0_240[k]
                   + f_10 * sgg_174[k]
                   - f_9 * pc_x[k] * sgh1_240[k];
    }
}

static auto
compute_prim_pgh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgh0, const size_t sgg,
                                                          const size_t sgh1, const size_t pdh0,
                                                          const size_t pdh1, const size_t pfh0,
                                                          const size_t pfg, const size_t pfh1,
                                                          const size_t pgf0, const size_t pgf1,
                                                          const size_t pgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);
    const auto f_17 = 1.0 / p;
    const auto f_18 = gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_0 = buffer.data(sgh0 + 0);
    const auto *sgh0_1 = buffer.data(sgh0 + 1);
    const auto *sgh0_3 = buffer.data(sgh0 + 3);
    const auto *sgh0_5 = buffer.data(sgh0 + 5);
    const auto *sgh0_6 = buffer.data(sgh0 + 6);
    const auto *sgh0_8 = buffer.data(sgh0 + 8);
    const auto *sgh0_9 = buffer.data(sgh0 + 9);
    const auto *sgh0_15 = buffer.data(sgh0 + 15);
    const auto *sgh0_20 = buffer.data(sgh0 + 20);
    const auto *sgh0_42 = buffer.data(sgh0 + 42);
    const auto *sgh0_47 = buffer.data(sgh0 + 47);
    const auto *sgh0_246 = buffer.data(sgh0 + 246);
    const auto *sgh0_248 = buffer.data(sgh0 + 248);
    const auto *sgh0_249 = buffer.data(sgh0 + 249);
    const auto *sgh0_251 = buffer.data(sgh0 + 251);
    const auto *sgh0_252 = buffer.data(sgh0 + 252);
    const auto *sgh0_255 = buffer.data(sgh0 + 255);
    const auto *sgh0_257 = buffer.data(sgh0 + 257);
    const auto *sgh0_258 = buffer.data(sgh0 + 258);
    const auto *sgh0_261 = buffer.data(sgh0 + 261);
    const auto *sgh0_267 = buffer.data(sgh0 + 267);
    const auto *sgh0_269 = buffer.data(sgh0 + 269);
    const auto *sgh0_270 = buffer.data(sgh0 + 270);
    const auto *sgh0_272 = buffer.data(sgh0 + 272);
    const auto *sgh0_276 = buffer.data(sgh0 + 276);
    const auto *sgh0_279 = buffer.data(sgh0 + 279);
    const auto *sgh0_288 = buffer.data(sgh0 + 288);
    const auto *sgh0_290 = buffer.data(sgh0 + 290);
    const auto *sgh0_291 = buffer.data(sgh0 + 291);
    const auto *sgh0_293 = buffer.data(sgh0 + 293);
    const auto *sgh0_294 = buffer.data(sgh0 + 294);
    const auto *sgh0_299 = buffer.data(sgh0 + 299);
    const auto *sgh0_303 = buffer.data(sgh0 + 303);
    const auto *sgh0_309 = buffer.data(sgh0 + 309);
    const auto *sgh0_311 = buffer.data(sgh0 + 311);
    const auto *sgh0_312 = buffer.data(sgh0 + 312);
    const auto *sgh0_314 = buffer.data(sgh0 + 314);

    const auto *sgg_0 = buffer.data(sgg + 0);
    const auto *sgg_1 = buffer.data(sgg + 1);
    const auto *sgg_3 = buffer.data(sgg + 3);
    const auto *sgg_5 = buffer.data(sgg + 5);
    const auto *sgg_10 = buffer.data(sgg + 10);
    const auto *sgg_14 = buffer.data(sgg + 14);
    const auto *sgg_29 = buffer.data(sgg + 29);
    const auto *sgg_175 = buffer.data(sgg + 175);
    const auto *sgg_177 = buffer.data(sgg + 177);
    const auto *sgg_179 = buffer.data(sgg + 179);
    const auto *sgg_180 = buffer.data(sgg + 180);
    const auto *sgg_183 = buffer.data(sgg + 183);
    const auto *sgg_185 = buffer.data(sgg + 185);
    const auto *sgg_186 = buffer.data(sgg + 186);
    const auto *sgg_189 = buffer.data(sgg + 189);
    const auto *sgg_190 = buffer.data(sgg + 190);
    const auto *sgg_192 = buffer.data(sgg + 192);
    const auto *sgg_194 = buffer.data(sgg + 194);
    const auto *sgg_198 = buffer.data(sgg + 198);
    const auto *sgg_201 = buffer.data(sgg + 201);
    const auto *sgg_205 = buffer.data(sgg + 205);
    const auto *sgg_207 = buffer.data(sgg + 207);
    const auto *sgg_209 = buffer.data(sgg + 209);
    const auto *sgg_210 = buffer.data(sgg + 210);
    const auto *sgg_215 = buffer.data(sgg + 215);
    const auto *sgg_219 = buffer.data(sgg + 219);
    const auto *sgg_220 = buffer.data(sgg + 220);
    const auto *sgg_222 = buffer.data(sgg + 222);
    const auto *sgg_224 = buffer.data(sgg + 224);

    const auto *sgh1_0 = buffer.data(sgh1 + 0);
    const auto *sgh1_1 = buffer.data(sgh1 + 1);
    const auto *sgh1_3 = buffer.data(sgh1 + 3);
    const auto *sgh1_5 = buffer.data(sgh1 + 5);
    const auto *sgh1_6 = buffer.data(sgh1 + 6);
    const auto *sgh1_8 = buffer.data(sgh1 + 8);
    const auto *sgh1_9 = buffer.data(sgh1 + 9);
    const auto *sgh1_15 = buffer.data(sgh1 + 15);
    const auto *sgh1_20 = buffer.data(sgh1 + 20);
    const auto *sgh1_42 = buffer.data(sgh1 + 42);
    const auto *sgh1_47 = buffer.data(sgh1 + 47);
    const auto *sgh1_246 = buffer.data(sgh1 + 246);
    const auto *sgh1_248 = buffer.data(sgh1 + 248);
    const auto *sgh1_249 = buffer.data(sgh1 + 249);
    const auto *sgh1_251 = buffer.data(sgh1 + 251);
    const auto *sgh1_252 = buffer.data(sgh1 + 252);
    const auto *sgh1_255 = buffer.data(sgh1 + 255);
    const auto *sgh1_257 = buffer.data(sgh1 + 257);
    const auto *sgh1_258 = buffer.data(sgh1 + 258);
    const auto *sgh1_261 = buffer.data(sgh1 + 261);
    const auto *sgh1_267 = buffer.data(sgh1 + 267);
    const auto *sgh1_269 = buffer.data(sgh1 + 269);
    const auto *sgh1_270 = buffer.data(sgh1 + 270);
    const auto *sgh1_272 = buffer.data(sgh1 + 272);
    const auto *sgh1_276 = buffer.data(sgh1 + 276);
    const auto *sgh1_279 = buffer.data(sgh1 + 279);
    const auto *sgh1_288 = buffer.data(sgh1 + 288);
    const auto *sgh1_290 = buffer.data(sgh1 + 290);
    const auto *sgh1_291 = buffer.data(sgh1 + 291);
    const auto *sgh1_293 = buffer.data(sgh1 + 293);
    const auto *sgh1_294 = buffer.data(sgh1 + 294);
    const auto *sgh1_299 = buffer.data(sgh1 + 299);
    const auto *sgh1_303 = buffer.data(sgh1 + 303);
    const auto *sgh1_309 = buffer.data(sgh1 + 309);
    const auto *sgh1_311 = buffer.data(sgh1 + 311);
    const auto *sgh1_312 = buffer.data(sgh1 + 312);
    const auto *sgh1_314 = buffer.data(sgh1 + 314);

    const auto *pdh0_162 = buffer.data(pdh0 + 162);

    const auto *pdh1_162 = buffer.data(pdh1 + 162);

    const auto *pfh0_189 = buffer.data(pfh0 + 189);
    const auto *pfh0_194 = buffer.data(pfh0 + 194);
    const auto *pfh0_198 = buffer.data(pfh0 + 198);
    const auto *pfh0_211 = buffer.data(pfh0 + 211);
    const auto *pfh0_213 = buffer.data(pfh0 + 213);
    const auto *pfh0_246 = buffer.data(pfh0 + 246);

    const auto *pfg_96 = buffer.data(pfg + 96);
    const auto *pfg_100 = buffer.data(pfg + 100);
    const auto *pfg_105 = buffer.data(pfg + 105);
    const auto *pfg_108 = buffer.data(pfg + 108);
    const auto *pfg_111 = buffer.data(pfg + 111);
    const auto *pfg_114 = buffer.data(pfg + 114);
    const auto *pfg_115 = buffer.data(pfg + 115);
    const auto *pfg_119 = buffer.data(pfg + 119);
    const auto *pfg_120 = buffer.data(pfg + 120);
    const auto *pfg_122 = buffer.data(pfg + 122);
    const auto *pfg_123 = buffer.data(pfg + 123);
    const auto *pfg_125 = buffer.data(pfg + 125);
    const auto *pfg_126 = buffer.data(pfg + 126);
    const auto *pfg_129 = buffer.data(pfg + 129);
    const auto *pfg_130 = buffer.data(pfg + 130);
    const auto *pfg_134 = buffer.data(pfg + 134);
    const auto *pfg_135 = buffer.data(pfg + 135);
    const auto *pfg_137 = buffer.data(pfg + 137);
    const auto *pfg_138 = buffer.data(pfg + 138);
    const auto *pfg_140 = buffer.data(pfg + 140);
    const auto *pfg_141 = buffer.data(pfg + 141);
    const auto *pfg_144 = buffer.data(pfg + 144);
    const auto *pfg_145 = buffer.data(pfg + 145);
    const auto *pfg_149 = buffer.data(pfg + 149);
    const auto *pfg_150 = buffer.data(pfg + 150);
    const auto *pfg_151 = buffer.data(pfg + 151);
    const auto *pfg_160 = buffer.data(pfg + 160);
    const auto *pfg_161 = buffer.data(pfg + 161);
    const auto *pfg_162 = buffer.data(pfg + 162);
    const auto *pfg_163 = buffer.data(pfg + 163);
    const auto *pfg_164 = buffer.data(pfg + 164);
    const auto *pfg_165 = buffer.data(pfg + 165);
    const auto *pfg_166 = buffer.data(pfg + 166);
    const auto *pfg_168 = buffer.data(pfg + 168);
    const auto *pfg_170 = buffer.data(pfg + 170);
    const auto *pfg_171 = buffer.data(pfg + 171);
    const auto *pfg_173 = buffer.data(pfg + 173);
    const auto *pfg_174 = buffer.data(pfg + 174);
    const auto *pfg_175 = buffer.data(pfg + 175);
    const auto *pfg_176 = buffer.data(pfg + 176);
    const auto *pfg_177 = buffer.data(pfg + 177);
    const auto *pfg_178 = buffer.data(pfg + 178);
    const auto *pfg_179 = buffer.data(pfg + 179);

    const auto *pfh1_189 = buffer.data(pfh1 + 189);
    const auto *pfh1_194 = buffer.data(pfh1 + 194);
    const auto *pfh1_198 = buffer.data(pfh1 + 198);
    const auto *pfh1_211 = buffer.data(pfh1 + 211);
    const auto *pfh1_213 = buffer.data(pfh1 + 213);
    const auto *pfh1_246 = buffer.data(pfh1 + 246);

    const auto *pgf0_140 = buffer.data(pgf0 + 140);
    const auto *pgf0_141 = buffer.data(pgf0 + 141);
    const auto *pgf0_156 = buffer.data(pgf0 + 156);
    const auto *pgf0_157 = buffer.data(pgf0 + 157);
    const auto *pgf0_160 = buffer.data(pgf0 + 160);
    const auto *pgf0_161 = buffer.data(pgf0 + 161);
    const auto *pgf0_163 = buffer.data(pgf0 + 163);
    const auto *pgf0_165 = buffer.data(pgf0 + 165);
    const auto *pgf0_166 = buffer.data(pgf0 + 166);
    const auto *pgf0_167 = buffer.data(pgf0 + 167);
    const auto *pgf0_168 = buffer.data(pgf0 + 168);
    const auto *pgf0_169 = buffer.data(pgf0 + 169);

    const auto *pgf1_140 = buffer.data(pgf1 + 140);
    const auto *pgf1_141 = buffer.data(pgf1 + 141);
    const auto *pgf1_156 = buffer.data(pgf1 + 156);
    const auto *pgf1_157 = buffer.data(pgf1 + 157);
    const auto *pgf1_160 = buffer.data(pgf1 + 160);
    const auto *pgf1_161 = buffer.data(pgf1 + 161);
    const auto *pgf1_163 = buffer.data(pgf1 + 163);
    const auto *pgf1_165 = buffer.data(pgf1 + 165);
    const auto *pgf1_166 = buffer.data(pgf1 + 166);
    const auto *pgf1_167 = buffer.data(pgf1 + 167);
    const auto *pgf1_168 = buffer.data(pgf1 + 168);
    const auto *pgf1_169 = buffer.data(pgf1 + 169);

    const auto *pgg_171 = buffer.data(pgg + 171);
    const auto *pgg_174 = buffer.data(pgg + 174);
    const auto *pgg_175 = buffer.data(pgg + 175);
    const auto *pgg_177 = buffer.data(pgg + 177);
    const auto *pgg_179 = buffer.data(pgg + 179);
    const auto *pgg_180 = buffer.data(pgg + 180);
    const auto *pgg_182 = buffer.data(pgg + 182);
    const auto *pgg_183 = buffer.data(pgg + 183);
    const auto *pgg_185 = buffer.data(pgg + 185);
    const auto *pgg_186 = buffer.data(pgg + 186);
    const auto *pgg_189 = buffer.data(pgg + 189);
    const auto *pgg_190 = buffer.data(pgg + 190);
    const auto *pgg_192 = buffer.data(pgg + 192);
    const auto *pgg_194 = buffer.data(pgg + 194);
    const auto *pgg_195 = buffer.data(pgg + 195);
    const auto *pgg_197 = buffer.data(pgg + 197);
    const auto *pgg_198 = buffer.data(pgg + 198);
    const auto *pgg_200 = buffer.data(pgg + 200);
    const auto *pgg_201 = buffer.data(pgg + 201);
    const auto *pgg_204 = buffer.data(pgg + 204);
    const auto *pgg_205 = buffer.data(pgg + 205);
    const auto *pgg_207 = buffer.data(pgg + 207);
    const auto *pgg_209 = buffer.data(pgg + 209);
    const auto *pgg_210 = buffer.data(pgg + 210);
    const auto *pgg_211 = buffer.data(pgg + 211);
    const auto *pgg_212 = buffer.data(pgg + 212);
    const auto *pgg_213 = buffer.data(pgg + 213);
    const auto *pgg_215 = buffer.data(pgg + 215);
    const auto *pgg_216 = buffer.data(pgg + 216);
    const auto *pgg_219 = buffer.data(pgg + 219);
    const auto *pgg_220 = buffer.data(pgg + 220);
    const auto *pgg_222 = buffer.data(pgg + 222);
    const auto *pgg_224 = buffer.data(pgg + 224);
    const auto *pgg_225 = buffer.data(pgg + 225);
    const auto *pgg_226 = buffer.data(pgg + 226);
    const auto *pgg_228 = buffer.data(pgg + 228);
    const auto *pgg_235 = buffer.data(pgg + 235);
    const auto *pgg_236 = buffer.data(pgg + 236);
    const auto *pgg_237 = buffer.data(pgg + 237);
    const auto *pgg_238 = buffer.data(pgg + 238);
    const auto *pgg_239 = buffer.data(pgg + 239);
    const auto *pgg_240 = buffer.data(pgg + 240);
    const auto *pgg_241 = buffer.data(pgg + 241);
    const auto *pgg_243 = buffer.data(pgg + 243);
    const auto *pgg_245 = buffer.data(pgg + 245);
    const auto *pgg_246 = buffer.data(pgg + 246);
    const auto *pgg_248 = buffer.data(pgg + 248);
    const auto *pgg_249 = buffer.data(pgg + 249);
    const auto *pgg_250 = buffer.data(pgg + 250);
    const auto *pgg_251 = buffer.data(pgg + 251);
    const auto *pgg_252 = buffer.data(pgg + 252);
    const auto *pgg_253 = buffer.data(pgg + 253);
    const auto *pgg_254 = buffer.data(pgg + 254);
    const auto *pgg_255 = buffer.data(pgg + 255);
    const auto *pgg_256 = buffer.data(pgg + 256);

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pc_x, pc_y, pc_z, sgg_175, sgg_177, \
                         pfg_96, pfg_114, pgg_171, pgg_174, pgg_175, \
                         pgg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_0 * sgg_175[k]
                   + f_4 * pc_x[k] * pgg_175[k];

        t_242[k] = f_0 * pfg_96[k]
                   + f_4 * pc_z[k] * pgg_171[k];

        t_243[k] = f_0 * sgg_177[k]
                   + f_4 * pc_x[k] * pgg_177[k];

        t_244[k] = f_11 * pfg_114[k]
                   + f_4 * pc_y[k] * pgg_174[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pa_x, pc_x, pc_z, sgh0_246, sgh0_248, \
                         sgg_179, sgh1_246, sgh1_248, pfg_100, pgg_175, \
                         pgg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_0 * sgg_179[k]
                   + f_4 * pc_x[k] * pgg_179[k];

        t_246[k] = pa_x[k] * sgh0_246[k]
                   - f_9 * pc_x[k] * sgh1_246[k];

        t_247[k] = f_0 * pfg_100[k]
                   + f_4 * pc_z[k] * pgg_175[k];

        t_248[k] = pa_x[k] * sgh0_248[k]
                   - f_9 * pc_x[k] * sgh1_248[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pa_x, pc_x, pc_y, sgh0_249, sgh0_251, \
                         sgh0_252, sgg_180, sgh1_249, sgh1_251, sgh1_252, pfg_119, \
                         pgg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = pa_x[k] * sgh0_249[k]
                   - f_9 * pc_x[k] * sgh1_249[k];

        t_250[k] = f_11 * pfg_119[k]
                   + f_4 * pc_y[k] * pgg_179[k];

        t_251[k] = pa_x[k] * sgh0_251[k]
                   - f_9 * pc_x[k] * sgh1_251[k];

        t_252[k] = pa_x[k] * sgh0_252[k]
                   + f_14 * sgg_180[k]
                   - f_9 * pc_x[k] * sgh1_252[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pa_x, pc_x, pc_y, pc_z, sgh0_255, \
                         sgg_183, sgh1_255, pfg_105, pfg_120, pfg_122, pgg_180, \
                         pgg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * pfg_120[k]
                   + f_4 * pc_y[k] * pgg_180[k];

        t_254[k] = f_10 * pfg_105[k]
                   + f_4 * pc_z[k] * pgg_180[k];

        t_255[k] = pa_x[k] * sgh0_255[k]
                   + f_11 * sgg_183[k]
                   - f_9 * pc_x[k] * sgh1_255[k];

        t_256[k] = f_10 * pfg_122[k]
                   + f_4 * pc_y[k] * pgg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_x, pc_x, pc_z, sgh0_257, sgh0_258, sgg_185, \
                         sgg_186, sgh1_257, sgh1_258, pfg_108, \
                         pgg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = pa_x[k] * sgh0_257[k]
                   + f_11 * sgg_185[k]
                   - f_9 * pc_x[k] * sgh1_257[k];

        t_258[k] = pa_x[k] * sgh0_258[k]
                   + f_10 * sgg_186[k]
                   - f_9 * pc_x[k] * sgh1_258[k];

        t_259[k] = f_10 * pfg_108[k]
                   + f_4 * pc_z[k] * pgg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_x, pc_x, pc_y, sgh0_261, sgg_189, sgg_190, \
                         sgh1_261, pfg_125, pgg_185, pgg_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * pfg_125[k]
                   + f_4 * pc_y[k] * pgg_185[k];

        t_261[k] = pa_x[k] * sgh0_261[k]
                   + f_10 * sgg_189[k]
                   - f_9 * pc_x[k] * sgh1_261[k];

        t_262[k] = f_0 * sgg_190[k]
                   + f_4 * pc_x[k] * pgg_190[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pc_x, pc_y, pc_z, sgg_192, sgg_194, \
                         pfg_111, pfg_129, pgg_186, pgg_189, pgg_192, \
                         pgg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_10 * pfg_111[k]
                   + f_4 * pc_z[k] * pgg_186[k];

        t_264[k] = f_0 * sgg_192[k]
                   + f_4 * pc_x[k] * pgg_192[k];

        t_265[k] = f_10 * pfg_129[k]
                   + f_4 * pc_y[k] * pgg_189[k];

        t_266[k] = f_0 * sgg_194[k]
                   + f_4 * pc_x[k] * pgg_194[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pa_x, pc_x, pc_z, sgh0_267, sgh0_269, \
                         sgh0_270, sgh1_267, sgh1_269, sgh1_270, pfg_115, \
                         pgg_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pa_x[k] * sgh0_267[k]
                   - f_9 * pc_x[k] * sgh1_267[k];

        t_268[k] = f_10 * pfg_115[k]
                   + f_4 * pc_z[k] * pgg_190[k];

        t_269[k] = pa_x[k] * sgh0_269[k]
                   - f_9 * pc_x[k] * sgh1_269[k];

        t_270[k] = pa_x[k] * sgh0_270[k]
                   - f_9 * pc_x[k] * sgh1_270[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_x, pb_y, pc_x, pc_y, sgh0_272, \
                         sgh1_272, pfh0_189, pfg_134, pfg_135, pfh1_189, pgg_194, \
                         pgg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * pfg_134[k]
                   + f_4 * pc_y[k] * pgg_194[k];

        t_272[k] = pa_x[k] * sgh0_272[k]
                   - f_9 * pc_x[k] * sgh1_272[k];

        t_273[k] = pb_y[k] * pfh0_189[k]
                   - f_9 * pc_y[k] * pfh1_189[k];

        t_274[k] = f_0 * pfg_135[k]
                   + f_4 * pc_y[k] * pgg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pa_x, pc_x, pc_y, pc_z, sgh0_276, sgg_198, \
                         sgh1_276, pfg_120, pfg_137, pgg_195, pgg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * pfg_120[k]
                   + f_4 * pc_z[k] * pgg_195[k];

        t_276[k] = pa_x[k] * sgh0_276[k]
                   + f_11 * sgg_198[k]
                   - f_9 * pc_x[k] * sgh1_276[k];

        t_277[k] = f_0 * pfg_137[k]
                   + f_4 * pc_y[k] * pgg_197[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pa_x, pb_y, pc_x, pc_y, pc_z, sgh0_279, sgg_201, \
                         sgh1_279, pfh0_194, pfg_123, pfh1_194, \
                         pgg_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = pb_y[k] * pfh0_194[k]
                   - f_9 * pc_y[k] * pfh1_194[k];

        t_279[k] = pa_x[k] * sgh0_279[k]
                   + f_10 * sgg_201[k]
                   - f_9 * pc_x[k] * sgh1_279[k];

        t_280[k] = f_11 * pfg_123[k]
                   + f_4 * pc_z[k] * pgg_198[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pb_y, pc_x, pc_y, pc_z, sgg_205, \
                         pfh0_198, pfg_126, pfg_140, pfh1_198, pgg_200, pgg_201, \
                         pgg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_0 * pfg_140[k]
                   + f_4 * pc_y[k] * pgg_200[k];

        t_282[k] = pb_y[k] * pfh0_198[k]
                   - f_9 * pc_y[k] * pfh1_198[k];

        t_283[k] = f_0 * sgg_205[k]
                   + f_4 * pc_x[k] * pgg_205[k];

        t_284[k] = f_11 * pfg_126[k]
                   + f_4 * pc_z[k] * pgg_201[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pa_x, pc_x, pc_y, sgh0_288, sgg_207, \
                         sgg_209, sgh1_288, pfg_144, pgg_204, pgg_207, \
                         pgg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_0 * sgg_207[k]
                   + f_4 * pc_x[k] * pgg_207[k];

        t_286[k] = f_0 * pfg_144[k]
                   + f_4 * pc_y[k] * pgg_204[k];

        t_287[k] = f_0 * sgg_209[k]
                   + f_4 * pc_x[k] * pgg_209[k];

        t_288[k] = pa_x[k] * sgh0_288[k]
                   - f_9 * pc_x[k] * sgh1_288[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_x, pc_x, pc_y, pc_z, sgh0_290, \
                         sgh0_291, sgh1_290, sgh1_291, pfg_130, pfg_149, pgg_205, \
                         pgg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_11 * pfg_130[k]
                   + f_4 * pc_z[k] * pgg_205[k];

        t_290[k] = pa_x[k] * sgh0_290[k]
                   - f_9 * pc_x[k] * sgh1_290[k];

        t_291[k] = pa_x[k] * sgh0_291[k]
                   - f_9 * pc_x[k] * sgh1_291[k];

        t_292[k] = f_0 * pfg_149[k]
                   + f_4 * pc_y[k] * pgg_209[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_x, pc_x, pc_y, pc_z, sgh0_293, \
                         sgh0_294, sgg_210, sgh1_293, sgh1_294, pfg_135, \
                         pgg_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pa_x[k] * sgh0_293[k]
                   - f_9 * pc_x[k] * sgh1_293[k];

        t_294[k] = pa_x[k] * sgh0_294[k]
                   + f_14 * sgg_210[k]
                   - f_9 * pc_x[k] * sgh1_294[k];

        t_295[k] = f_4 * pc_y[k] * pgg_210[k];

        t_296[k] = f_1 * pfg_135[k]
                   + f_4 * pc_z[k] * pgg_210[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pa_x, pc_x, pc_y, sgh0_299, sgg_215, sgh1_299, \
                         pgf0_140, pgf1_140, pgg_211, pgg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_5 * pgf0_140[k]
                   - f_6 * pgf1_140[k]
                   + f_4 * pc_y[k] * pgg_211[k];

        t_298[k] = f_4 * pc_y[k] * pgg_212[k];

        t_299[k] = pa_x[k] * sgh0_299[k]
                   + f_11 * sgg_215[k]
                   - f_9 * pc_x[k] * sgh1_299[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_x, pc_x, pc_y, pc_z, sgh0_303, \
                         sgg_219, sgh1_303, pfg_138, pgf0_141, pgf1_141, pgg_213, \
                         pgg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_7 * pgf0_141[k]
                   - f_8 * pgf1_141[k]
                   + f_4 * pc_y[k] * pgg_213[k];

        t_301[k] = f_1 * pfg_138[k]
                   + f_4 * pc_z[k] * pgg_213[k];

        t_302[k] = f_4 * pc_y[k] * pgg_215[k];

        t_303[k] = pa_x[k] * sgh0_303[k]
                   + f_10 * sgg_219[k]
                   - f_9 * pc_x[k] * sgh1_303[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, pc_z, sgg_220, sgg_222, \
                         pfg_141, pgg_216, pgg_219, pgg_220, pgg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_0 * sgg_220[k]
                   + f_4 * pc_x[k] * pgg_220[k];

        t_305[k] = f_1 * pfg_141[k]
                   + f_4 * pc_z[k] * pgg_216[k];

        t_306[k] = f_0 * sgg_222[k]
                   + f_4 * pc_x[k] * pgg_222[k];

        t_307[k] = f_4 * pc_y[k] * pgg_219[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_x, pc_x, pc_z, sgh0_309, sgh0_311, \
                         sgg_224, sgh1_309, sgh1_311, pfg_145, pgg_220, \
                         pgg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_0 * sgg_224[k]
                   + f_4 * pc_x[k] * pgg_224[k];

        t_309[k] = pa_x[k] * sgh0_309[k]
                   - f_9 * pc_x[k] * sgh1_309[k];

        t_310[k] = f_1 * pfg_145[k]
                   + f_4 * pc_z[k] * pgg_220[k];

        t_311[k] = pa_x[k] * sgh0_311[k]
                   - f_9 * pc_x[k] * sgh1_311[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pa_x, pa_y, pc_x, pc_y, sgh0_0, sgh0_312, \
                         sgh0_314, sgh1_0, sgh1_312, sgh1_314, \
                         pgg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = pa_x[k] * sgh0_312[k]
                   - f_9 * pc_x[k] * sgh1_312[k];

        t_313[k] = f_4 * pc_y[k] * pgg_224[k];

        t_314[k] = pa_x[k] * sgh0_314[k]
                   - f_9 * pc_x[k] * sgh1_314[k];

        t_315[k] = pa_y[k] * sgh0_0[k]
                   - f_9 * pc_y[k] * sgh1_0[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_y, pc_y, pc_z, sgh0_1, sgh0_3, sgg_0, \
                         sgg_1, sgh1_1, sgh1_3, pgg_225, pgg_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_y[k] * sgh0_1[k]
                   + f_0 * sgg_0[k]
                   - f_9 * pc_y[k] * sgh1_1[k];

        t_317[k] = f_4 * pc_z[k] * pgg_225[k];

        t_318[k] = pa_y[k] * sgh0_3[k]
                   + f_10 * sgg_1[k]
                   - f_9 * pc_y[k] * sgh1_3[k];

        t_319[k] = f_4 * pc_z[k] * pgg_226[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_y, pc_y, pc_z, sgh0_5, sgh0_6, sgh0_8, \
                         sgg_3, sgg_5, sgh1_5, sgh1_6, sgh1_8, \
                         pgg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = pa_y[k] * sgh0_5[k]
                   - f_9 * pc_y[k] * sgh1_5[k];

        t_321[k] = pa_y[k] * sgh0_6[k]
                   + f_11 * sgg_3[k]
                   - f_9 * pc_y[k] * sgh1_6[k];

        t_322[k] = f_4 * pc_z[k] * pgg_228[k];

        t_323[k] = pa_y[k] * sgh0_8[k]
                   + f_0 * sgg_5[k]
                   - f_9 * pc_y[k] * sgh1_8[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_y, pc_x, pc_y, sgh0_9, sgh1_9, \
                         pfg_160, pfg_161, pfg_162, pgg_235, pgg_236, \
                         pgg_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = pa_y[k] * sgh0_9[k]
                   - f_9 * pc_y[k] * sgh1_9[k];

        t_325[k] = f_1 * pfg_160[k]
                   + f_4 * pc_x[k] * pgg_235[k];

        t_326[k] = f_1 * pfg_161[k]
                   + f_4 * pc_x[k] * pgg_236[k];

        t_327[k] = f_1 * pfg_162[k]
                   + f_4 * pc_x[k] * pgg_237[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_y, pc_x, pc_y, pc_z, sgh0_15, sgg_10, \
                         sgh1_15, pfg_163, pfg_164, pgg_235, pgg_238, \
                         pgg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_1 * pfg_163[k]
                   + f_4 * pc_x[k] * pgg_238[k];

        t_329[k] = f_1 * pfg_164[k]
                   + f_4 * pc_x[k] * pgg_239[k];

        t_330[k] = pa_y[k] * sgh0_15[k]
                   + f_14 * sgg_10[k]
                   - f_9 * pc_y[k] * sgh1_15[k];

        t_331[k] = f_4 * pc_z[k] * pgg_235[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pc_y, pc_z, sgg_14, pgf0_156, pgf0_157, \
                         pgf1_156, pgf1_157, pgg_236, pgg_237, \
                         pgg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_5 * pgf0_156[k]
                   - f_6 * pgf1_156[k]
                   + f_4 * pc_z[k] * pgg_236[k];

        t_333[k] = f_7 * pgf0_157[k]
                   - f_8 * pgf1_157[k]
                   + f_4 * pc_z[k] * pgg_237[k];

        t_334[k] = f_0 * sgg_14[k]
                   + f_4 * pc_y[k] * pgg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pa_y, pc_x, pc_y, sgh0_20, sgh1_20, pfg_165, \
                         pfg_166, pgf0_160, pgf0_161, pgf1_160, pgf1_161, pgg_240, \
                         pgg_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pa_y[k] * sgh0_20[k]
                   - f_9 * pc_y[k] * sgh1_20[k];

        t_336[k] = f_11 * pfg_165[k]
                   + f_2 * pgf0_160[k]
                   - f_3 * pgf1_160[k]
                   + f_4 * pc_x[k] * pgg_240[k];

        t_337[k] = f_11 * pfg_166[k]
                   + f_15 * pgf0_161[k]
                   - f_16 * pgf1_161[k]
                   + f_4 * pc_x[k] * pgg_241[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pc_x, pc_z, pfg_168, pfg_170, pgf0_163, \
                         pgf0_165, pgf1_163, pgf1_165, pgg_240, pgg_241, pgg_243, \
                         pgg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_4 * pc_z[k] * pgg_240[k];

        t_339[k] = f_11 * pfg_168[k]
                   + f_7 * pgf0_163[k]
                   - f_8 * pgf1_163[k]
                   + f_4 * pc_x[k] * pgg_243[k];

        t_340[k] = f_4 * pc_z[k] * pgg_241[k];

        t_341[k] = f_11 * pfg_170[k]
                   + f_7 * pgf0_165[k]
                   - f_8 * pgf1_165[k]
                   + f_4 * pc_x[k] * pgg_245[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_z, pfg_171, pfg_173, pgf0_166, \
                         pgf0_168, pgf1_166, pgf1_168, pgg_243, pgg_246, \
                         pgg_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_11 * pfg_171[k]
                   + f_5 * pgf0_166[k]
                   - f_6 * pgf1_166[k]
                   + f_4 * pc_x[k] * pgg_246[k];

        t_343[k] = f_4 * pc_z[k] * pgg_243[k];

        t_344[k] = f_11 * pfg_173[k]
                   + f_5 * pgf0_168[k]
                   - f_6 * pgf1_168[k]
                   + f_4 * pc_x[k] * pgg_248[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pc_x, pfg_174, pfg_175, pfg_176, pfg_177, \
                         pgf0_169, pgf1_169, pgg_249, pgg_250, pgg_251, \
                         pgg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_11 * pfg_174[k]
                   + f_5 * pgf0_169[k]
                   - f_6 * pgf1_169[k]
                   + f_4 * pc_x[k] * pgg_249[k];

        t_346[k] = f_11 * pfg_175[k]
                   + f_4 * pc_x[k] * pgg_250[k];

        t_347[k] = f_11 * pfg_176[k]
                   + f_4 * pc_x[k] * pgg_251[k];

        t_348[k] = f_11 * pfg_177[k]
                   + f_4 * pc_x[k] * pgg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pb_x, pc_x, pc_z, pdh0_162, pdh1_162, \
                         pfh0_246, pfg_178, pfg_179, pfh1_246, pgg_250, pgg_253, \
                         pgg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_11 * pfg_178[k]
                   + f_4 * pc_x[k] * pgg_253[k];

        t_350[k] = f_11 * pfg_179[k]
                   + f_4 * pc_x[k] * pgg_254[k];

        t_351[k] = f_17 * pdh0_162[k]
                   - f_18 * pdh1_162[k]
                   + pb_x[k] * pfh0_246[k]
                   - f_9 * pc_x[k] * pfh1_246[k];

        t_352[k] = f_4 * pc_z[k] * pgg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, pc_z, sgg_29, pfg_164, pgf0_166, pgf0_167, \
                         pgf1_166, pgf1_167, pgg_251, pgg_252, \
                         pgg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_5 * pgf0_166[k]
                   - f_6 * pgf1_166[k]
                   + f_4 * pc_z[k] * pgg_251[k];

        t_354[k] = f_7 * pgf0_167[k]
                   - f_8 * pgf1_167[k]
                   + f_4 * pc_z[k] * pgg_252[k];

        t_355[k] = f_0 * sgg_29[k]
                   + f_0 * pfg_164[k]
                   + f_4 * pc_y[k] * pgg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pa_y, pb_z, pc_y, pc_z, sgh0_42, sgh1_42, \
                         pfh0_211, pfh1_211, pgf0_169, pgf1_169, \
                         pgg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_2 * pgf0_169[k]
                   - f_3 * pgf1_169[k]
                   + f_4 * pc_z[k] * pgg_254[k];

        t_357[k] = pa_y[k] * sgh0_42[k]
                   - f_9 * pc_y[k] * sgh1_42[k];

        t_358[k] = pb_z[k] * pfh0_211[k]
                   - f_9 * pc_z[k] * pfh1_211[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pa_y, pb_z, pc_y, pc_z, sgh0_47, sgh1_47, \
                         pfh0_213, pfg_150, pfg_151, pfh1_213, pgg_255, \
                         pgg_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_0 * pfg_150[k]
                   + f_4 * pc_z[k] * pgg_255[k];

        t_360[k] = pb_z[k] * pfh0_213[k]
                   - f_9 * pc_z[k] * pfh1_213[k];

        t_361[k] = f_0 * pfg_151[k]
                   + f_4 * pc_z[k] * pgg_256[k];

        t_362[k] = pa_y[k] * sgh0_47[k]
                   - f_9 * pc_y[k] * sgh1_47[k];
    }
}

static auto
compute_prim_pgh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgh0, const size_t sgg,
                                                          const size_t sgh1, const size_t pdh0,
                                                          const size_t pdh1, const size_t pfh0,
                                                          const size_t pfg, const size_t pfh1,
                                                          const size_t pgf0, const size_t pgf1,
                                                          const size_t pgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);
    const auto f_14 = 2.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_50 = buffer.data(sgh0 + 50);
    const auto *sgh0_51 = buffer.data(sgh0 + 51);
    const auto *sgh0_59 = buffer.data(sgh0 + 59);
    const auto *sgh0_60 = buffer.data(sgh0 + 60);
    const auto *sgh0_62 = buffer.data(sgh0 + 62);
    const auto *sgh0_105 = buffer.data(sgh0 + 105);
    const auto *sgh0_106 = buffer.data(sgh0 + 106);
    const auto *sgh0_108 = buffer.data(sgh0 + 108);
    const auto *sgh0_110 = buffer.data(sgh0 + 110);
    const auto *sgh0_111 = buffer.data(sgh0 + 111);
    const auto *sgh0_113 = buffer.data(sgh0 + 113);
    const auto *sgh0_114 = buffer.data(sgh0 + 114);
    const auto *sgh0_120 = buffer.data(sgh0 + 120);
    const auto *sgh0_122 = buffer.data(sgh0 + 122);
    const auto *sgh0_123 = buffer.data(sgh0 + 123);
    const auto *sgh0_125 = buffer.data(sgh0 + 125);

    const auto *sgg_35 = buffer.data(sgg + 35);
    const auto *sgg_42 = buffer.data(sgg + 42);
    const auto *sgg_43 = buffer.data(sgg + 43);
    const auto *sgg_44 = buffer.data(sgg + 44);
    const auto *sgg_59 = buffer.data(sgg + 59);
    const auto *sgg_74 = buffer.data(sgg + 74);
    const auto *sgg_75 = buffer.data(sgg + 75);
    const auto *sgg_76 = buffer.data(sgg + 76);
    const auto *sgg_78 = buffer.data(sgg + 78);
    const auto *sgg_80 = buffer.data(sgg + 80);
    const auto *sgg_85 = buffer.data(sgg + 85);
    const auto *sgg_87 = buffer.data(sgg + 87);
    const auto *sgg_88 = buffer.data(sgg + 88);
    const auto *sgg_89 = buffer.data(sgg + 89);

    const auto *sgh1_50 = buffer.data(sgh1 + 50);
    const auto *sgh1_51 = buffer.data(sgh1 + 51);
    const auto *sgh1_59 = buffer.data(sgh1 + 59);
    const auto *sgh1_60 = buffer.data(sgh1 + 60);
    const auto *sgh1_62 = buffer.data(sgh1 + 62);
    const auto *sgh1_105 = buffer.data(sgh1 + 105);
    const auto *sgh1_106 = buffer.data(sgh1 + 106);
    const auto *sgh1_108 = buffer.data(sgh1 + 108);
    const auto *sgh1_110 = buffer.data(sgh1 + 110);
    const auto *sgh1_111 = buffer.data(sgh1 + 111);
    const auto *sgh1_113 = buffer.data(sgh1 + 113);
    const auto *sgh1_114 = buffer.data(sgh1 + 114);
    const auto *sgh1_120 = buffer.data(sgh1 + 120);
    const auto *sgh1_122 = buffer.data(sgh1 + 122);
    const auto *sgh1_123 = buffer.data(sgh1 + 123);
    const auto *sgh1_125 = buffer.data(sgh1 + 125);

    const auto *pdh0_204 = buffer.data(pdh0 + 204);

    const auto *pdh1_204 = buffer.data(pdh1 + 204);

    const auto *pfh0_216 = buffer.data(pfh0 + 216);
    const auto *pfh0_225 = buffer.data(pfh0 + 225);
    const auto *pfh0_232 = buffer.data(pfh0 + 232);
    const auto *pfh0_234 = buffer.data(pfh0 + 234);
    const auto *pfh0_237 = buffer.data(pfh0 + 237);
    const auto *pfh0_246 = buffer.data(pfh0 + 246);
    const auto *pfh0_273 = buffer.data(pfh0 + 273);
    const auto *pfh0_274 = buffer.data(pfh0 + 274);
    const auto *pfh0_276 = buffer.data(pfh0 + 276);
    const auto *pfh0_279 = buffer.data(pfh0 + 279);
    const auto *pfh0_288 = buffer.data(pfh0 + 288);
    const auto *pfh0_336 = buffer.data(pfh0 + 336);
    const auto *pfh0_337 = buffer.data(pfh0 + 337);
    const auto *pfh0_339 = buffer.data(pfh0 + 339);
    const auto *pfh0_341 = buffer.data(pfh0 + 341);
    const auto *pfh0_342 = buffer.data(pfh0 + 342);
    const auto *pfh0_344 = buffer.data(pfh0 + 344);
    const auto *pfh0_345 = buffer.data(pfh0 + 345);
    const auto *pfh0_351 = buffer.data(pfh0 + 351);
    const auto *pfh0_353 = buffer.data(pfh0 + 353);
    const auto *pfh0_354 = buffer.data(pfh0 + 354);
    const auto *pfh0_355 = buffer.data(pfh0 + 355);
    const auto *pfh0_356 = buffer.data(pfh0 + 356);
    const auto *pfh0_362 = buffer.data(pfh0 + 362);
    const auto *pfh0_365 = buffer.data(pfh0 + 365);
    const auto *pfh0_366 = buffer.data(pfh0 + 366);
    const auto *pfh0_372 = buffer.data(pfh0 + 372);

    const auto *pfg_153 = buffer.data(pfg + 153);
    const auto *pfg_160 = buffer.data(pfg + 160);
    const auto *pfg_165 = buffer.data(pfg + 165);
    const auto *pfg_166 = buffer.data(pfg + 166);
    const auto *pfg_168 = buffer.data(pfg + 168);
    const auto *pfg_175 = buffer.data(pfg + 175);
    const auto *pfg_176 = buffer.data(pfg + 176);
    const auto *pfg_177 = buffer.data(pfg + 177);
    const auto *pfg_179 = buffer.data(pfg + 179);
    const auto *pfg_180 = buffer.data(pfg + 180);
    const auto *pfg_181 = buffer.data(pfg + 181);
    const auto *pfg_183 = buffer.data(pfg + 183);
    const auto *pfg_190 = buffer.data(pfg + 190);
    const auto *pfg_191 = buffer.data(pfg + 191);
    const auto *pfg_192 = buffer.data(pfg + 192);
    const auto *pfg_193 = buffer.data(pfg + 193);
    const auto *pfg_194 = buffer.data(pfg + 194);
    const auto *pfg_195 = buffer.data(pfg + 195);
    const auto *pfg_196 = buffer.data(pfg + 196);
    const auto *pfg_198 = buffer.data(pfg + 198);
    const auto *pfg_200 = buffer.data(pfg + 200);
    const auto *pfg_201 = buffer.data(pfg + 201);
    const auto *pfg_203 = buffer.data(pfg + 203);
    const auto *pfg_204 = buffer.data(pfg + 204);
    const auto *pfg_205 = buffer.data(pfg + 205);
    const auto *pfg_206 = buffer.data(pfg + 206);
    const auto *pfg_207 = buffer.data(pfg + 207);
    const auto *pfg_208 = buffer.data(pfg + 208);
    const auto *pfg_209 = buffer.data(pfg + 209);
    const auto *pfg_210 = buffer.data(pfg + 210);
    const auto *pfg_215 = buffer.data(pfg + 215);
    const auto *pfg_218 = buffer.data(pfg + 218);
    const auto *pfg_219 = buffer.data(pfg + 219);
    const auto *pfg_220 = buffer.data(pfg + 220);
    const auto *pfg_221 = buffer.data(pfg + 221);
    const auto *pfg_222 = buffer.data(pfg + 222);
    const auto *pfg_223 = buffer.data(pfg + 223);
    const auto *pfg_224 = buffer.data(pfg + 224);
    const auto *pfg_235 = buffer.data(pfg + 235);
    const auto *pfg_236 = buffer.data(pfg + 236);
    const auto *pfg_237 = buffer.data(pfg + 237);
    const auto *pfg_238 = buffer.data(pfg + 238);
    const auto *pfg_239 = buffer.data(pfg + 239);
    const auto *pfg_240 = buffer.data(pfg + 240);
    const auto *pfg_241 = buffer.data(pfg + 241);
    const auto *pfg_243 = buffer.data(pfg + 243);
    const auto *pfg_245 = buffer.data(pfg + 245);
    const auto *pfg_246 = buffer.data(pfg + 246);
    const auto *pfg_248 = buffer.data(pfg + 248);
    const auto *pfg_249 = buffer.data(pfg + 249);
    const auto *pfg_250 = buffer.data(pfg + 250);
    const auto *pfg_251 = buffer.data(pfg + 251);
    const auto *pfg_252 = buffer.data(pfg + 252);
    const auto *pfg_253 = buffer.data(pfg + 253);
    const auto *pfg_254 = buffer.data(pfg + 254);
    const auto *pfg_260 = buffer.data(pfg + 260);
    const auto *pfg_263 = buffer.data(pfg + 263);
    const auto *pfg_264 = buffer.data(pfg + 264);
    const auto *pfg_265 = buffer.data(pfg + 265);
    const auto *pfg_266 = buffer.data(pfg + 266);
    const auto *pfg_267 = buffer.data(pfg + 267);
    const auto *pfg_268 = buffer.data(pfg + 268);
    const auto *pfg_269 = buffer.data(pfg + 269);

    const auto *pfh1_216 = buffer.data(pfh1 + 216);
    const auto *pfh1_225 = buffer.data(pfh1 + 225);
    const auto *pfh1_232 = buffer.data(pfh1 + 232);
    const auto *pfh1_234 = buffer.data(pfh1 + 234);
    const auto *pfh1_237 = buffer.data(pfh1 + 237);
    const auto *pfh1_246 = buffer.data(pfh1 + 246);
    const auto *pfh1_273 = buffer.data(pfh1 + 273);
    const auto *pfh1_274 = buffer.data(pfh1 + 274);
    const auto *pfh1_276 = buffer.data(pfh1 + 276);
    const auto *pfh1_279 = buffer.data(pfh1 + 279);
    const auto *pfh1_288 = buffer.data(pfh1 + 288);
    const auto *pfh1_336 = buffer.data(pfh1 + 336);
    const auto *pfh1_337 = buffer.data(pfh1 + 337);
    const auto *pfh1_339 = buffer.data(pfh1 + 339);
    const auto *pfh1_341 = buffer.data(pfh1 + 341);
    const auto *pfh1_342 = buffer.data(pfh1 + 342);
    const auto *pfh1_344 = buffer.data(pfh1 + 344);
    const auto *pfh1_345 = buffer.data(pfh1 + 345);
    const auto *pfh1_351 = buffer.data(pfh1 + 351);
    const auto *pfh1_353 = buffer.data(pfh1 + 353);
    const auto *pfh1_354 = buffer.data(pfh1 + 354);
    const auto *pfh1_355 = buffer.data(pfh1 + 355);
    const auto *pfh1_356 = buffer.data(pfh1 + 356);
    const auto *pfh1_362 = buffer.data(pfh1 + 362);
    const auto *pfh1_365 = buffer.data(pfh1 + 365);
    const auto *pfh1_366 = buffer.data(pfh1 + 366);
    const auto *pfh1_372 = buffer.data(pfh1 + 372);

    const auto *pgf0_180 = buffer.data(pgf0 + 180);
    const auto *pgf0_181 = buffer.data(pgf0 + 181);
    const auto *pgf0_183 = buffer.data(pgf0 + 183);
    const auto *pgf0_185 = buffer.data(pgf0 + 185);
    const auto *pgf0_186 = buffer.data(pgf0 + 186);
    const auto *pgf0_187 = buffer.data(pgf0 + 187);
    const auto *pgf0_188 = buffer.data(pgf0 + 188);
    const auto *pgf0_189 = buffer.data(pgf0 + 189);
    const auto *pgf0_190 = buffer.data(pgf0 + 190);
    const auto *pgf0_195 = buffer.data(pgf0 + 195);
    const auto *pgf0_196 = buffer.data(pgf0 + 196);
    const auto *pgf0_197 = buffer.data(pgf0 + 197);
    const auto *pgf0_198 = buffer.data(pgf0 + 198);
    const auto *pgf0_199 = buffer.data(pgf0 + 199);

    const auto *pgf1_180 = buffer.data(pgf1 + 180);
    const auto *pgf1_181 = buffer.data(pgf1 + 181);
    const auto *pgf1_183 = buffer.data(pgf1 + 183);
    const auto *pgf1_185 = buffer.data(pgf1 + 185);
    const auto *pgf1_186 = buffer.data(pgf1 + 186);
    const auto *pgf1_187 = buffer.data(pgf1 + 187);
    const auto *pgf1_188 = buffer.data(pgf1 + 188);
    const auto *pgf1_189 = buffer.data(pgf1 + 189);
    const auto *pgf1_190 = buffer.data(pgf1 + 190);
    const auto *pgf1_195 = buffer.data(pgf1 + 195);
    const auto *pgf1_196 = buffer.data(pgf1 + 196);
    const auto *pgf1_197 = buffer.data(pgf1 + 197);
    const auto *pgf1_198 = buffer.data(pgf1 + 198);
    const auto *pgf1_199 = buffer.data(pgf1 + 199);

    const auto *pgg_258 = buffer.data(pgg + 258);
    const auto *pgg_265 = buffer.data(pgg + 265);
    const auto *pgg_266 = buffer.data(pgg + 266);
    const auto *pgg_267 = buffer.data(pgg + 267);
    const auto *pgg_268 = buffer.data(pgg + 268);
    const auto *pgg_269 = buffer.data(pgg + 269);
    const auto *pgg_270 = buffer.data(pgg + 270);
    const auto *pgg_271 = buffer.data(pgg + 271);
    const auto *pgg_273 = buffer.data(pgg + 273);
    const auto *pgg_275 = buffer.data(pgg + 275);
    const auto *pgg_276 = buffer.data(pgg + 276);
    const auto *pgg_278 = buffer.data(pgg + 278);
    const auto *pgg_279 = buffer.data(pgg + 279);
    const auto *pgg_280 = buffer.data(pgg + 280);
    const auto *pgg_281 = buffer.data(pgg + 281);
    const auto *pgg_282 = buffer.data(pgg + 282);
    const auto *pgg_283 = buffer.data(pgg + 283);
    const auto *pgg_284 = buffer.data(pgg + 284);
    const auto *pgg_285 = buffer.data(pgg + 285);
    const auto *pgg_286 = buffer.data(pgg + 286);
    const auto *pgg_288 = buffer.data(pgg + 288);
    const auto *pgg_290 = buffer.data(pgg + 290);
    const auto *pgg_293 = buffer.data(pgg + 293);
    const auto *pgg_294 = buffer.data(pgg + 294);
    const auto *pgg_295 = buffer.data(pgg + 295);
    const auto *pgg_296 = buffer.data(pgg + 296);
    const auto *pgg_297 = buffer.data(pgg + 297);
    const auto *pgg_298 = buffer.data(pgg + 298);
    const auto *pgg_299 = buffer.data(pgg + 299);
    const auto *pgg_300 = buffer.data(pgg + 300);
    const auto *pgg_301 = buffer.data(pgg + 301);
    const auto *pgg_303 = buffer.data(pgg + 303);
    const auto *pgg_310 = buffer.data(pgg + 310);
    const auto *pgg_311 = buffer.data(pgg + 311);
    const auto *pgg_312 = buffer.data(pgg + 312);
    const auto *pgg_313 = buffer.data(pgg + 313);
    const auto *pgg_314 = buffer.data(pgg + 314);
    const auto *pgg_315 = buffer.data(pgg + 315);
    const auto *pgg_316 = buffer.data(pgg + 316);
    const auto *pgg_318 = buffer.data(pgg + 318);
    const auto *pgg_325 = buffer.data(pgg + 325);
    const auto *pgg_326 = buffer.data(pgg + 326);
    const auto *pgg_327 = buffer.data(pgg + 327);
    const auto *pgg_328 = buffer.data(pgg + 328);
    const auto *pgg_329 = buffer.data(pgg + 329);
    const auto *pgg_330 = buffer.data(pgg + 330);
    const auto *pgg_331 = buffer.data(pgg + 331);
    const auto *pgg_333 = buffer.data(pgg + 333);
    const auto *pgg_340 = buffer.data(pgg + 340);
    const auto *pgg_341 = buffer.data(pgg + 341);
    const auto *pgg_342 = buffer.data(pgg + 342);
    const auto *pgg_343 = buffer.data(pgg + 343);
    const auto *pgg_344 = buffer.data(pgg + 344);

#pragma omp simd aligned(t_363, t_364, t_365, pa_y, pb_z, pc_y, pc_z, sgh0_50, sgg_35, \
                         sgh1_50, pfh0_216, pfg_153, pfh1_216, \
                         pgg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pb_z[k] * pfh0_216[k]
                   - f_9 * pc_z[k] * pfh1_216[k];

        t_364[k] = f_0 * pfg_153[k]
                   + f_4 * pc_z[k] * pgg_258[k];

        t_365[k] = pa_y[k] * sgh0_50[k]
                   + f_0 * sgg_35[k]
                   - f_9 * pc_y[k] * sgh1_50[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pc_x, pc_y, sgh0_51, sgh1_51, \
                         pfg_190, pfg_191, pfg_192, pgg_265, pgg_266, \
                         pgg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = pa_y[k] * sgh0_51[k]
                   - f_9 * pc_y[k] * sgh1_51[k];

        t_367[k] = f_11 * pfg_190[k]
                   + f_4 * pc_x[k] * pgg_265[k];

        t_368[k] = f_11 * pfg_191[k]
                   + f_4 * pc_x[k] * pgg_266[k];

        t_369[k] = f_11 * pfg_192[k]
                   + f_4 * pc_x[k] * pgg_267[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pb_z, pc_x, pc_z, pfh0_225, pfg_160, \
                         pfg_193, pfg_194, pfh1_225, pgg_265, pgg_268, \
                         pgg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_11 * pfg_193[k]
                   + f_4 * pc_x[k] * pgg_268[k];

        t_371[k] = f_11 * pfg_194[k]
                   + f_4 * pc_x[k] * pgg_269[k];

        t_372[k] = pb_z[k] * pfh0_225[k]
                   - f_9 * pc_z[k] * pfh1_225[k];

        t_373[k] = f_0 * pfg_160[k]
                   + f_4 * pc_z[k] * pgg_265[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, pa_y, pc_y, sgh0_59, sgh0_60, sgh0_62, \
                         sgg_42, sgg_43, sgg_44, sgh1_59, sgh1_60, sgh1_62, \
                         pgg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * sgh0_59[k]
                   + f_11 * sgg_42[k]
                   - f_9 * pc_y[k] * sgh1_59[k];

        t_375[k] = pa_y[k] * sgh0_60[k]
                   + f_10 * sgg_43[k]
                   - f_9 * pc_y[k] * sgh1_60[k];

        t_376[k] = f_0 * sgg_44[k]
                   + f_4 * pc_y[k] * pgg_269[k];

        t_377[k] = pa_y[k] * sgh0_62[k]
                   - f_9 * pc_y[k] * sgh1_62[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pc_x, pc_z, pfg_195, pfg_196, pgf0_180, \
                         pgf0_181, pgf1_180, pgf1_181, pgg_270, \
                         pgg_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_10 * pfg_195[k]
                   + f_2 * pgf0_180[k]
                   - f_3 * pgf1_180[k]
                   + f_4 * pc_x[k] * pgg_270[k];

        t_379[k] = f_10 * pfg_196[k]
                   + f_15 * pgf0_181[k]
                   - f_16 * pgf1_181[k]
                   + f_4 * pc_x[k] * pgg_271[k];

        t_380[k] = f_4 * pc_z[k] * pgg_270[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, pc_x, pc_z, pfg_198, pfg_200, pgf0_183, \
                         pgf0_185, pgf1_183, pgf1_185, pgg_271, pgg_273, \
                         pgg_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_10 * pfg_198[k]
                   + f_7 * pgf0_183[k]
                   - f_8 * pgf1_183[k]
                   + f_4 * pc_x[k] * pgg_273[k];

        t_382[k] = f_4 * pc_z[k] * pgg_271[k];

        t_383[k] = f_10 * pfg_200[k]
                   + f_7 * pgf0_185[k]
                   - f_8 * pgf1_185[k]
                   + f_4 * pc_x[k] * pgg_275[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pc_x, pc_z, pfg_201, pfg_203, pgf0_186, \
                         pgf0_188, pgf1_186, pgf1_188, pgg_273, pgg_276, \
                         pgg_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_10 * pfg_201[k]
                   + f_5 * pgf0_186[k]
                   - f_6 * pgf1_186[k]
                   + f_4 * pc_x[k] * pgg_276[k];

        t_385[k] = f_4 * pc_z[k] * pgg_273[k];

        t_386[k] = f_10 * pfg_203[k]
                   + f_5 * pgf0_188[k]
                   - f_6 * pgf1_188[k]
                   + f_4 * pc_x[k] * pgg_278[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pc_x, pfg_204, pfg_205, pfg_206, pfg_207, \
                         pgf0_189, pgf1_189, pgg_279, pgg_280, pgg_281, \
                         pgg_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_10 * pfg_204[k]
                   + f_5 * pgf0_189[k]
                   - f_6 * pgf1_189[k]
                   + f_4 * pc_x[k] * pgg_279[k];

        t_388[k] = f_10 * pfg_205[k]
                   + f_4 * pc_x[k] * pgg_280[k];

        t_389[k] = f_10 * pfg_206[k]
                   + f_4 * pc_x[k] * pgg_281[k];

        t_390[k] = f_10 * pfg_207[k]
                   + f_4 * pc_x[k] * pgg_282[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_x, pc_x, pc_z, pdh0_204, pdh1_204, \
                         pfh0_288, pfg_208, pfg_209, pfh1_288, pgg_280, pgg_283, \
                         pgg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_10 * pfg_208[k]
                   + f_4 * pc_x[k] * pgg_283[k];

        t_392[k] = f_10 * pfg_209[k]
                   + f_4 * pc_x[k] * pgg_284[k];

        t_393[k] = f_12 * pdh0_204[k]
                   - f_13 * pdh1_204[k]
                   + pb_x[k] * pfh0_288[k]
                   - f_9 * pc_x[k] * pfh1_288[k];

        t_394[k] = f_4 * pc_z[k] * pgg_280[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_y, pc_z, sgg_59, pfg_179, pgf0_186, pgf0_187, \
                         pgf1_186, pgf1_187, pgg_281, pgg_282, \
                         pgg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_5 * pgf0_186[k]
                   - f_6 * pgf1_186[k]
                   + f_4 * pc_z[k] * pgg_281[k];

        t_396[k] = f_7 * pgf0_187[k]
                   - f_8 * pgf1_187[k]
                   + f_4 * pc_z[k] * pgg_282[k];

        t_397[k] = f_0 * sgg_59[k]
                   + f_10 * pfg_179[k]
                   + f_4 * pc_y[k] * pgg_284[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pb_z, pc_x, pc_z, pfh0_232, pfg_210, pfh1_232, \
                         pgf0_189, pgf0_190, pgf1_189, pgf1_190, pgg_284, \
                         pgg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_2 * pgf0_189[k]
                   - f_3 * pgf1_189[k]
                   + f_4 * pc_z[k] * pgg_284[k];

        t_399[k] = f_10 * pfg_210[k]
                   + f_2 * pgf0_190[k]
                   - f_3 * pgf1_190[k]
                   + f_4 * pc_x[k] * pgg_285[k];

        t_400[k] = pb_z[k] * pfh0_232[k]
                   - f_9 * pc_z[k] * pfh1_232[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pb_z, pc_z, pfh0_234, pfg_165, pfg_166, \
                         pfh1_234, pgg_285, pgg_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_0 * pfg_165[k]
                   + f_4 * pc_z[k] * pgg_285[k];

        t_402[k] = pb_z[k] * pfh0_234[k]
                   - f_9 * pc_z[k] * pfh1_234[k];

        t_403[k] = f_0 * pfg_166[k]
                   + f_4 * pc_z[k] * pgg_286[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pb_z, pc_x, pc_z, pfh0_237, pfg_168, pfg_215, \
                         pfh1_237, pgf0_195, pgf1_195, pgg_288, \
                         pgg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_10 * pfg_215[k]
                   + f_7 * pgf0_195[k]
                   - f_8 * pgf1_195[k]
                   + f_4 * pc_x[k] * pgg_290[k];

        t_405[k] = pb_z[k] * pfh0_237[k]
                   - f_9 * pc_z[k] * pfh1_237[k];

        t_406[k] = f_0 * pfg_168[k]
                   + f_4 * pc_z[k] * pgg_288[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pc_x, pfg_218, pfg_219, pfg_220, pgf0_198, \
                         pgf0_199, pgf1_198, pgf1_199, pgg_293, pgg_294, \
                         pgg_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_10 * pfg_218[k]
                   + f_5 * pgf0_198[k]
                   - f_6 * pgf1_198[k]
                   + f_4 * pc_x[k] * pgg_293[k];

        t_408[k] = f_10 * pfg_219[k]
                   + f_5 * pgf0_199[k]
                   - f_6 * pgf1_199[k]
                   + f_4 * pc_x[k] * pgg_294[k];

        t_409[k] = f_10 * pfg_220[k]
                   + f_4 * pc_x[k] * pgg_295[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pfg_221, pfg_222, pfg_223, pfg_224, \
                         pgg_296, pgg_297, pgg_298, pgg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_10 * pfg_221[k]
                   + f_4 * pc_x[k] * pgg_296[k];

        t_411[k] = f_10 * pfg_222[k]
                   + f_4 * pc_x[k] * pgg_297[k];

        t_412[k] = f_10 * pfg_223[k]
                   + f_4 * pc_x[k] * pgg_298[k];

        t_413[k] = f_10 * pfg_224[k]
                   + f_4 * pc_x[k] * pgg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pb_z, pc_z, pfh0_246, pfg_175, pfg_176, \
                         pfh1_246, pgf0_196, pgf1_196, pgg_295, \
                         pgg_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = pb_z[k] * pfh0_246[k]
                   - f_9 * pc_z[k] * pfh1_246[k];

        t_415[k] = f_0 * pfg_175[k]
                   + f_4 * pc_z[k] * pgg_295[k];

        t_416[k] = f_0 * pfg_176[k]
                   + f_5 * pgf0_196[k]
                   - f_6 * pgf1_196[k]
                   + f_4 * pc_z[k] * pgg_296[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_y, pc_z, sgg_74, pfg_177, pfg_179, pfg_194, \
                         pgf0_197, pgf0_199, pgf1_197, pgf1_199, pgg_297, \
                         pgg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_0 * pfg_177[k]
                   + f_7 * pgf0_197[k]
                   - f_8 * pgf1_197[k]
                   + f_4 * pc_z[k] * pgg_297[k];

        t_418[k] = f_0 * sgg_74[k]
                   + f_0 * pfg_194[k]
                   + f_4 * pc_y[k] * pgg_299[k];

        t_419[k] = f_0 * pfg_179[k]
                   + f_2 * pgf0_199[k]
                   - f_3 * pgf1_199[k]
                   + f_4 * pc_z[k] * pgg_299[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pa_y, pc_y, pc_z, sgh0_105, sgh0_106, sgg_75, \
                         sgh1_105, sgh1_106, pfg_180, pgg_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pa_y[k] * sgh0_105[k]
                   - f_9 * pc_y[k] * sgh1_105[k];

        t_421[k] = pa_y[k] * sgh0_106[k]
                   + f_0 * sgg_75[k]
                   - f_9 * pc_y[k] * sgh1_106[k];

        t_422[k] = f_10 * pfg_180[k]
                   + f_4 * pc_z[k] * pgg_300[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pa_y, pc_y, pc_z, sgh0_108, sgh0_110, sgg_76, \
                         sgh1_108, sgh1_110, pfg_181, pgg_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = pa_y[k] * sgh0_108[k]
                   + f_10 * sgg_76[k]
                   - f_9 * pc_y[k] * sgh1_108[k];

        t_424[k] = f_10 * pfg_181[k]
                   + f_4 * pc_z[k] * pgg_301[k];

        t_425[k] = pa_y[k] * sgh0_110[k]
                   - f_9 * pc_y[k] * sgh1_110[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pa_y, pc_y, pc_z, sgh0_111, sgh0_113, sgg_78, \
                         sgg_80, sgh1_111, sgh1_113, pfg_183, pgg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = pa_y[k] * sgh0_111[k]
                   + f_11 * sgg_78[k]
                   - f_9 * pc_y[k] * sgh1_111[k];

        t_427[k] = f_10 * pfg_183[k]
                   + f_4 * pc_z[k] * pgg_303[k];

        t_428[k] = pa_y[k] * sgh0_113[k]
                   + f_0 * sgg_80[k]
                   - f_9 * pc_y[k] * sgh1_113[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pa_y, pc_x, pc_y, sgh0_114, sgh1_114, \
                         pfg_235, pfg_236, pfg_237, pgg_310, pgg_311, \
                         pgg_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = pa_y[k] * sgh0_114[k]
                   - f_9 * pc_y[k] * sgh1_114[k];

        t_430[k] = f_10 * pfg_235[k]
                   + f_4 * pc_x[k] * pgg_310[k];

        t_431[k] = f_10 * pfg_236[k]
                   + f_4 * pc_x[k] * pgg_311[k];

        t_432[k] = f_10 * pfg_237[k]
                   + f_4 * pc_x[k] * pgg_312[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, pa_y, pc_x, pc_y, sgh0_120, sgg_85, sgh1_120, \
                         pfg_238, pfg_239, pgg_313, pgg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_10 * pfg_238[k]
                   + f_4 * pc_x[k] * pgg_313[k];

        t_434[k] = f_10 * pfg_239[k]
                   + f_4 * pc_x[k] * pgg_314[k];

        t_435[k] = pa_y[k] * sgh0_120[k]
                   + f_14 * sgg_85[k]
                   - f_9 * pc_y[k] * sgh1_120[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, pa_y, pc_y, pc_z, sgh0_122, sgh0_123, sgg_87, \
                         sgg_88, sgh1_122, sgh1_123, pfg_190, pgg_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_10 * pfg_190[k]
                   + f_4 * pc_z[k] * pgg_310[k];

        t_437[k] = pa_y[k] * sgh0_122[k]
                   + f_11 * sgg_87[k]
                   - f_9 * pc_y[k] * sgh1_122[k];

        t_438[k] = pa_y[k] * sgh0_123[k]
                   + f_10 * sgg_88[k]
                   - f_9 * pc_y[k] * sgh1_123[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, pa_y, pb_x, pc_x, pc_y, sgh0_125, sgg_89, \
                         sgh1_125, pfh0_336, pfg_240, pfh1_336, \
                         pgg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_0 * sgg_89[k]
                   + f_4 * pc_y[k] * pgg_314[k];

        t_440[k] = pa_y[k] * sgh0_125[k]
                   - f_9 * pc_y[k] * sgh1_125[k];

        t_441[k] = pb_x[k] * pfh0_336[k]
                   + f_14 * pfg_240[k]
                   - f_9 * pc_x[k] * pfh1_336[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, pb_x, pc_x, pc_z, pfh0_337, pfh0_339, \
                         pfg_241, pfg_243, pfh1_337, pfh1_339, pgg_315, \
                         pgg_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = pb_x[k] * pfh0_337[k]
                   + f_1 * pfg_241[k]
                   - f_9 * pc_x[k] * pfh1_337[k];

        t_443[k] = f_4 * pc_z[k] * pgg_315[k];

        t_444[k] = pb_x[k] * pfh0_339[k]
                   + f_11 * pfg_243[k]
                   - f_9 * pc_x[k] * pfh1_339[k];

        t_445[k] = f_4 * pc_z[k] * pgg_316[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pb_x, pc_x, pc_z, pfh0_341, pfh0_342, pfg_245, \
                         pfg_246, pfh1_341, pfh1_342, pgg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = pb_x[k] * pfh0_341[k]
                   + f_11 * pfg_245[k]
                   - f_9 * pc_x[k] * pfh1_341[k];

        t_447[k] = pb_x[k] * pfh0_342[k]
                   + f_10 * pfg_246[k]
                   - f_9 * pc_x[k] * pfh1_342[k];

        t_448[k] = f_4 * pc_z[k] * pgg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pb_x, pc_x, pfh0_344, pfh0_345, pfg_248, \
                         pfg_249, pfg_250, pfg_251, pfh1_344, pfh1_345, pgg_325, \
                         pgg_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = pb_x[k] * pfh0_344[k]
                   + f_10 * pfg_248[k]
                   - f_9 * pc_x[k] * pfh1_344[k];

        t_450[k] = pb_x[k] * pfh0_345[k]
                   + f_10 * pfg_249[k]
                   - f_9 * pc_x[k] * pfh1_345[k];

        t_451[k] = f_0 * pfg_250[k]
                   + f_4 * pc_x[k] * pgg_325[k];

        t_452[k] = f_0 * pfg_251[k]
                   + f_4 * pc_x[k] * pgg_326[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pb_x, pc_x, pfh0_351, pfg_252, pfg_253, \
                         pfg_254, pfh1_351, pgg_327, pgg_328, pgg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_0 * pfg_252[k]
                   + f_4 * pc_x[k] * pgg_327[k];

        t_454[k] = f_0 * pfg_253[k]
                   + f_4 * pc_x[k] * pgg_328[k];

        t_455[k] = f_0 * pfg_254[k]
                   + f_4 * pc_x[k] * pgg_329[k];

        t_456[k] = pb_x[k] * pfh0_351[k]
                   - f_9 * pc_x[k] * pfh1_351[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pb_x, pc_x, pc_z, pfh0_353, pfh0_354, \
                         pfh0_355, pfh1_353, pfh1_354, pfh1_355, \
                         pgg_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_4 * pc_z[k] * pgg_325[k];

        t_458[k] = pb_x[k] * pfh0_353[k]
                   - f_9 * pc_x[k] * pfh1_353[k];

        t_459[k] = pb_x[k] * pfh0_354[k]
                   - f_9 * pc_x[k] * pfh1_354[k];

        t_460[k] = pb_x[k] * pfh0_355[k]
                   - f_9 * pc_x[k] * pfh1_355[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pb_x, pb_z, pc_x, pc_z, pfh0_273, \
                         pfh0_274, pfh0_356, pfg_195, pfh1_273, pfh1_274, pfh1_356, \
                         pgg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = pb_x[k] * pfh0_356[k]
                   - f_9 * pc_x[k] * pfh1_356[k];

        t_462[k] = pb_z[k] * pfh0_273[k]
                   - f_9 * pc_z[k] * pfh1_273[k];

        t_463[k] = pb_z[k] * pfh0_274[k]
                   - f_9 * pc_z[k] * pfh1_274[k];

        t_464[k] = f_0 * pfg_195[k]
                   + f_4 * pc_z[k] * pgg_330[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pb_x, pb_z, pc_x, pc_z, pfh0_276, pfh0_362, \
                         pfg_196, pfg_260, pfh1_276, pfh1_362, \
                         pgg_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pb_z[k] * pfh0_276[k]
                   - f_9 * pc_z[k] * pfh1_276[k];

        t_466[k] = f_0 * pfg_196[k]
                   + f_4 * pc_z[k] * pgg_331[k];

        t_467[k] = pb_x[k] * pfh0_362[k]
                   + f_11 * pfg_260[k]
                   - f_9 * pc_x[k] * pfh1_362[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pb_x, pb_z, pc_x, pc_z, pfh0_279, pfh0_365, \
                         pfg_198, pfg_263, pfh1_279, pfh1_365, \
                         pgg_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = pb_z[k] * pfh0_279[k]
                   - f_9 * pc_z[k] * pfh1_279[k];

        t_469[k] = f_0 * pfg_198[k]
                   + f_4 * pc_z[k] * pgg_333[k];

        t_470[k] = pb_x[k] * pfh0_365[k]
                   + f_10 * pfg_263[k]
                   - f_9 * pc_x[k] * pfh1_365[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pb_x, pc_x, pfh0_366, pfg_264, pfg_265, \
                         pfg_266, pfg_267, pfh1_366, pgg_340, pgg_341, \
                         pgg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pb_x[k] * pfh0_366[k]
                   + f_10 * pfg_264[k]
                   - f_9 * pc_x[k] * pfh1_366[k];

        t_472[k] = f_0 * pfg_265[k]
                   + f_4 * pc_x[k] * pgg_340[k];

        t_473[k] = f_0 * pfg_266[k]
                   + f_4 * pc_x[k] * pgg_341[k];

        t_474[k] = f_0 * pfg_267[k]
                   + f_4 * pc_x[k] * pgg_342[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pb_x, pc_x, pc_z, pfh0_372, pfg_205, \
                         pfg_268, pfg_269, pfh1_372, pgg_340, pgg_343, \
                         pgg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_0 * pfg_268[k]
                   + f_4 * pc_x[k] * pgg_343[k];

        t_476[k] = f_0 * pfg_269[k]
                   + f_4 * pc_x[k] * pgg_344[k];

        t_477[k] = pb_x[k] * pfh0_372[k]
                   - f_9 * pc_x[k] * pfh1_372[k];

        t_478[k] = f_0 * pfg_205[k]
                   + f_4 * pc_z[k] * pgg_340[k];
    }
}

static auto
compute_prim_pgh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgh0, const size_t sgg,
                                                          const size_t sgh1, const size_t pdh0,
                                                          const size_t pdh1, const size_t pfh0,
                                                          const size_t pfg, const size_t pfh1,
                                                          const size_t pgf0, const size_t pgf1,
                                                          const size_t pgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_189 = buffer.data(sgh0 + 189);
    const auto *sgh0_190 = buffer.data(sgh0 + 190);
    const auto *sgh0_192 = buffer.data(sgh0 + 192);
    const auto *sgh0_194 = buffer.data(sgh0 + 194);
    const auto *sgh0_195 = buffer.data(sgh0 + 195);
    const auto *sgh0_197 = buffer.data(sgh0 + 197);
    const auto *sgh0_198 = buffer.data(sgh0 + 198);
    const auto *sgh0_209 = buffer.data(sgh0 + 209);

    const auto *sgg_135 = buffer.data(sgg + 135);
    const auto *sgg_136 = buffer.data(sgg + 136);
    const auto *sgg_138 = buffer.data(sgg + 138);
    const auto *sgg_140 = buffer.data(sgg + 140);
    const auto *sgg_149 = buffer.data(sgg + 149);
    const auto *sgg_160 = buffer.data(sgg + 160);
    const auto *sgg_164 = buffer.data(sgg + 164);
    const auto *sgg_179 = buffer.data(sgg + 179);
    const auto *sgg_194 = buffer.data(sgg + 194);

    const auto *sgh1_189 = buffer.data(sgh1 + 189);
    const auto *sgh1_190 = buffer.data(sgh1 + 190);
    const auto *sgh1_192 = buffer.data(sgh1 + 192);
    const auto *sgh1_194 = buffer.data(sgh1 + 194);
    const auto *sgh1_195 = buffer.data(sgh1 + 195);
    const auto *sgh1_197 = buffer.data(sgh1 + 197);
    const auto *sgh1_198 = buffer.data(sgh1 + 198);
    const auto *sgh1_209 = buffer.data(sgh1 + 209);

    const auto *pdh0_204 = buffer.data(pdh0 + 204);

    const auto *pdh1_204 = buffer.data(pdh1 + 204);

    const auto *pfh0_336 = buffer.data(pfh0 + 336);
    const auto *pfh0_337 = buffer.data(pfh0 + 337);
    const auto *pfh0_339 = buffer.data(pfh0 + 339);
    const auto *pfh0_342 = buffer.data(pfh0 + 342);
    const auto *pfh0_351 = buffer.data(pfh0 + 351);
    const auto *pfh0_353 = buffer.data(pfh0 + 353);
    const auto *pfh0_354 = buffer.data(pfh0 + 354);
    const auto *pfh0_372 = buffer.data(pfh0 + 372);
    const auto *pfh0_374 = buffer.data(pfh0 + 374);
    const auto *pfh0_375 = buffer.data(pfh0 + 375);
    const auto *pfh0_376 = buffer.data(pfh0 + 376);
    const auto *pfh0_377 = buffer.data(pfh0 + 377);
    const auto *pfh0_379 = buffer.data(pfh0 + 379);
    const auto *pfh0_381 = buffer.data(pfh0 + 381);
    const auto *pfh0_384 = buffer.data(pfh0 + 384);
    const auto *pfh0_386 = buffer.data(pfh0 + 386);
    const auto *pfh0_393 = buffer.data(pfh0 + 393);
    const auto *pfh0_395 = buffer.data(pfh0 + 395);
    const auto *pfh0_396 = buffer.data(pfh0 + 396);
    const auto *pfh0_397 = buffer.data(pfh0 + 397);
    const auto *pfh0_398 = buffer.data(pfh0 + 398);
    const auto *pfh0_414 = buffer.data(pfh0 + 414);
    const auto *pfh0_416 = buffer.data(pfh0 + 416);
    const auto *pfh0_417 = buffer.data(pfh0 + 417);

    const auto *pfg_210 = buffer.data(pfg + 210);
    const auto *pfg_211 = buffer.data(pfg + 211);
    const auto *pfg_213 = buffer.data(pfg + 213);
    const auto *pfg_220 = buffer.data(pfg + 220);
    const auto *pfg_225 = buffer.data(pfg + 225);
    const auto *pfg_226 = buffer.data(pfg + 226);
    const auto *pfg_228 = buffer.data(pfg + 228);
    const auto *pfg_235 = buffer.data(pfg + 235);
    const auto *pfg_240 = buffer.data(pfg + 240);
    const auto *pfg_241 = buffer.data(pfg + 241);
    const auto *pfg_243 = buffer.data(pfg + 243);
    const auto *pfg_250 = buffer.data(pfg + 250);
    const auto *pfg_251 = buffer.data(pfg + 251);
    const auto *pfg_252 = buffer.data(pfg + 252);
    const auto *pfg_254 = buffer.data(pfg + 254);
    const auto *pfg_255 = buffer.data(pfg + 255);
    const auto *pfg_256 = buffer.data(pfg + 256);
    const auto *pfg_258 = buffer.data(pfg + 258);
    const auto *pfg_265 = buffer.data(pfg + 265);
    const auto *pfg_266 = buffer.data(pfg + 266);
    const auto *pfg_267 = buffer.data(pfg + 267);
    const auto *pfg_269 = buffer.data(pfg + 269);
    const auto *pfg_270 = buffer.data(pfg + 270);
    const auto *pfg_271 = buffer.data(pfg + 271);
    const auto *pfg_273 = buffer.data(pfg + 273);
    const auto *pfg_275 = buffer.data(pfg + 275);
    const auto *pfg_276 = buffer.data(pfg + 276);
    const auto *pfg_278 = buffer.data(pfg + 278);
    const auto *pfg_279 = buffer.data(pfg + 279);
    const auto *pfg_280 = buffer.data(pfg + 280);
    const auto *pfg_281 = buffer.data(pfg + 281);
    const auto *pfg_282 = buffer.data(pfg + 282);
    const auto *pfg_283 = buffer.data(pfg + 283);
    const auto *pfg_284 = buffer.data(pfg + 284);
    const auto *pfg_295 = buffer.data(pfg + 295);
    const auto *pfg_296 = buffer.data(pfg + 296);
    const auto *pfg_297 = buffer.data(pfg + 297);
    const auto *pfg_298 = buffer.data(pfg + 298);
    const auto *pfg_299 = buffer.data(pfg + 299);

    const auto *pfh1_336 = buffer.data(pfh1 + 336);
    const auto *pfh1_337 = buffer.data(pfh1 + 337);
    const auto *pfh1_339 = buffer.data(pfh1 + 339);
    const auto *pfh1_342 = buffer.data(pfh1 + 342);
    const auto *pfh1_351 = buffer.data(pfh1 + 351);
    const auto *pfh1_353 = buffer.data(pfh1 + 353);
    const auto *pfh1_354 = buffer.data(pfh1 + 354);
    const auto *pfh1_372 = buffer.data(pfh1 + 372);
    const auto *pfh1_374 = buffer.data(pfh1 + 374);
    const auto *pfh1_375 = buffer.data(pfh1 + 375);
    const auto *pfh1_376 = buffer.data(pfh1 + 376);
    const auto *pfh1_377 = buffer.data(pfh1 + 377);
    const auto *pfh1_379 = buffer.data(pfh1 + 379);
    const auto *pfh1_381 = buffer.data(pfh1 + 381);
    const auto *pfh1_384 = buffer.data(pfh1 + 384);
    const auto *pfh1_386 = buffer.data(pfh1 + 386);
    const auto *pfh1_393 = buffer.data(pfh1 + 393);
    const auto *pfh1_395 = buffer.data(pfh1 + 395);
    const auto *pfh1_396 = buffer.data(pfh1 + 396);
    const auto *pfh1_397 = buffer.data(pfh1 + 397);
    const auto *pfh1_398 = buffer.data(pfh1 + 398);
    const auto *pfh1_414 = buffer.data(pfh1 + 414);
    const auto *pfh1_416 = buffer.data(pfh1 + 416);
    const auto *pfh1_417 = buffer.data(pfh1 + 417);

    const auto *pgf0_230 = buffer.data(pgf0 + 230);
    const auto *pgf0_235 = buffer.data(pgf0 + 235);
    const auto *pgf0_239 = buffer.data(pgf0 + 239);
    const auto *pgf0_250 = buffer.data(pgf0 + 250);
    const auto *pgf0_251 = buffer.data(pgf0 + 251);
    const auto *pgf0_253 = buffer.data(pgf0 + 253);
    const auto *pgf0_255 = buffer.data(pgf0 + 255);
    const auto *pgf0_256 = buffer.data(pgf0 + 256);
    const auto *pgf0_257 = buffer.data(pgf0 + 257);
    const auto *pgf0_258 = buffer.data(pgf0 + 258);
    const auto *pgf0_259 = buffer.data(pgf0 + 259);
    const auto *pgf0_265 = buffer.data(pgf0 + 265);
    const auto *pgf0_268 = buffer.data(pgf0 + 268);
    const auto *pgf0_269 = buffer.data(pgf0 + 269);
    const auto *pgf0_270 = buffer.data(pgf0 + 270);
    const auto *pgf0_271 = buffer.data(pgf0 + 271);
    const auto *pgf0_273 = buffer.data(pgf0 + 273);
    const auto *pgf0_275 = buffer.data(pgf0 + 275);
    const auto *pgf0_276 = buffer.data(pgf0 + 276);
    const auto *pgf0_277 = buffer.data(pgf0 + 277);
    const auto *pgf0_278 = buffer.data(pgf0 + 278);
    const auto *pgf0_279 = buffer.data(pgf0 + 279);
    const auto *pgf0_280 = buffer.data(pgf0 + 280);
    const auto *pgf0_281 = buffer.data(pgf0 + 281);
    const auto *pgf0_283 = buffer.data(pgf0 + 283);
    const auto *pgf0_285 = buffer.data(pgf0 + 285);
    const auto *pgf0_286 = buffer.data(pgf0 + 286);
    const auto *pgf0_288 = buffer.data(pgf0 + 288);
    const auto *pgf0_289 = buffer.data(pgf0 + 289);

    const auto *pgf1_230 = buffer.data(pgf1 + 230);
    const auto *pgf1_235 = buffer.data(pgf1 + 235);
    const auto *pgf1_239 = buffer.data(pgf1 + 239);
    const auto *pgf1_250 = buffer.data(pgf1 + 250);
    const auto *pgf1_251 = buffer.data(pgf1 + 251);
    const auto *pgf1_253 = buffer.data(pgf1 + 253);
    const auto *pgf1_255 = buffer.data(pgf1 + 255);
    const auto *pgf1_256 = buffer.data(pgf1 + 256);
    const auto *pgf1_257 = buffer.data(pgf1 + 257);
    const auto *pgf1_258 = buffer.data(pgf1 + 258);
    const auto *pgf1_259 = buffer.data(pgf1 + 259);
    const auto *pgf1_265 = buffer.data(pgf1 + 265);
    const auto *pgf1_268 = buffer.data(pgf1 + 268);
    const auto *pgf1_269 = buffer.data(pgf1 + 269);
    const auto *pgf1_270 = buffer.data(pgf1 + 270);
    const auto *pgf1_271 = buffer.data(pgf1 + 271);
    const auto *pgf1_273 = buffer.data(pgf1 + 273);
    const auto *pgf1_275 = buffer.data(pgf1 + 275);
    const auto *pgf1_276 = buffer.data(pgf1 + 276);
    const auto *pgf1_277 = buffer.data(pgf1 + 277);
    const auto *pgf1_278 = buffer.data(pgf1 + 278);
    const auto *pgf1_279 = buffer.data(pgf1 + 279);
    const auto *pgf1_280 = buffer.data(pgf1 + 280);
    const auto *pgf1_281 = buffer.data(pgf1 + 281);
    const auto *pgf1_283 = buffer.data(pgf1 + 283);
    const auto *pgf1_285 = buffer.data(pgf1 + 285);
    const auto *pgf1_286 = buffer.data(pgf1 + 286);
    const auto *pgf1_288 = buffer.data(pgf1 + 288);
    const auto *pgf1_289 = buffer.data(pgf1 + 289);

    const auto *pgg_345 = buffer.data(pgg + 345);
    const auto *pgg_346 = buffer.data(pgg + 346);
    const auto *pgg_348 = buffer.data(pgg + 348);
    const auto *pgg_350 = buffer.data(pgg + 350);
    const auto *pgg_354 = buffer.data(pgg + 354);
    const auto *pgg_355 = buffer.data(pgg + 355);
    const auto *pgg_356 = buffer.data(pgg + 356);
    const auto *pgg_357 = buffer.data(pgg + 357);
    const auto *pgg_358 = buffer.data(pgg + 358);
    const auto *pgg_359 = buffer.data(pgg + 359);
    const auto *pgg_360 = buffer.data(pgg + 360);
    const auto *pgg_361 = buffer.data(pgg + 361);
    const auto *pgg_363 = buffer.data(pgg + 363);
    const auto *pgg_370 = buffer.data(pgg + 370);
    const auto *pgg_371 = buffer.data(pgg + 371);
    const auto *pgg_372 = buffer.data(pgg + 372);
    const auto *pgg_373 = buffer.data(pgg + 373);
    const auto *pgg_374 = buffer.data(pgg + 374);
    const auto *pgg_375 = buffer.data(pgg + 375);
    const auto *pgg_376 = buffer.data(pgg + 376);
    const auto *pgg_378 = buffer.data(pgg + 378);
    const auto *pgg_380 = buffer.data(pgg + 380);
    const auto *pgg_381 = buffer.data(pgg + 381);
    const auto *pgg_383 = buffer.data(pgg + 383);
    const auto *pgg_384 = buffer.data(pgg + 384);
    const auto *pgg_385 = buffer.data(pgg + 385);
    const auto *pgg_386 = buffer.data(pgg + 386);
    const auto *pgg_387 = buffer.data(pgg + 387);
    const auto *pgg_388 = buffer.data(pgg + 388);
    const auto *pgg_389 = buffer.data(pgg + 389);
    const auto *pgg_390 = buffer.data(pgg + 390);
    const auto *pgg_391 = buffer.data(pgg + 391);
    const auto *pgg_393 = buffer.data(pgg + 393);
    const auto *pgg_395 = buffer.data(pgg + 395);
    const auto *pgg_398 = buffer.data(pgg + 398);
    const auto *pgg_399 = buffer.data(pgg + 399);
    const auto *pgg_400 = buffer.data(pgg + 400);
    const auto *pgg_401 = buffer.data(pgg + 401);
    const auto *pgg_402 = buffer.data(pgg + 402);
    const auto *pgg_403 = buffer.data(pgg + 403);
    const auto *pgg_404 = buffer.data(pgg + 404);
    const auto *pgg_405 = buffer.data(pgg + 405);
    const auto *pgg_406 = buffer.data(pgg + 406);
    const auto *pgg_408 = buffer.data(pgg + 408);
    const auto *pgg_410 = buffer.data(pgg + 410);
    const auto *pgg_411 = buffer.data(pgg + 411);
    const auto *pgg_413 = buffer.data(pgg + 413);
    const auto *pgg_414 = buffer.data(pgg + 414);
    const auto *pgg_415 = buffer.data(pgg + 415);
    const auto *pgg_416 = buffer.data(pgg + 416);
    const auto *pgg_417 = buffer.data(pgg + 417);
    const auto *pgg_418 = buffer.data(pgg + 418);
    const auto *pgg_419 = buffer.data(pgg + 419);
    const auto *pgg_420 = buffer.data(pgg + 420);
    const auto *pgg_421 = buffer.data(pgg + 421);
    const auto *pgg_423 = buffer.data(pgg + 423);
    const auto *pgg_425 = buffer.data(pgg + 425);
    const auto *pgg_426 = buffer.data(pgg + 426);
    const auto *pgg_428 = buffer.data(pgg + 428);
    const auto *pgg_429 = buffer.data(pgg + 429);
    const auto *pgg_430 = buffer.data(pgg + 430);
    const auto *pgg_431 = buffer.data(pgg + 431);
    const auto *pgg_432 = buffer.data(pgg + 432);
    const auto *pgg_433 = buffer.data(pgg + 433);
    const auto *pgg_434 = buffer.data(pgg + 434);

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pb_x, pc_x, pfh0_374, pfh0_375, pfh0_376, \
                         pfh0_377, pfh1_374, pfh1_375, pfh1_376, \
                         pfh1_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = pb_x[k] * pfh0_374[k]
                   - f_9 * pc_x[k] * pfh1_374[k];

        t_480[k] = pb_x[k] * pfh0_375[k]
                   - f_9 * pc_x[k] * pfh1_375[k];

        t_481[k] = pb_x[k] * pfh0_376[k]
                   - f_9 * pc_x[k] * pfh1_376[k];

        t_482[k] = pb_x[k] * pfh0_377[k]
                   - f_9 * pc_x[k] * pfh1_377[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, pb_x, pc_x, pc_z, pfh0_379, pfg_210, pfg_270, \
                         pfg_271, pfh1_379, pgf0_230, pgf1_230, \
                         pgg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_0 * pfg_270[k]
                   + f_2 * pgf0_230[k]
                   - f_3 * pgf1_230[k]
                   + f_4 * pc_x[k] * pgg_345[k];

        t_484[k] = pb_x[k] * pfh0_379[k]
                   + f_1 * pfg_271[k]
                   - f_9 * pc_x[k] * pfh1_379[k];

        t_485[k] = f_10 * pfg_210[k]
                   + f_4 * pc_z[k] * pgg_345[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pb_x, pc_x, pc_z, pfh0_381, pfg_211, pfg_273, \
                         pfg_275, pfh1_381, pgf0_235, pgf1_235, pgg_346, \
                         pgg_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = pb_x[k] * pfh0_381[k]
                   + f_11 * pfg_273[k]
                   - f_9 * pc_x[k] * pfh1_381[k];

        t_487[k] = f_10 * pfg_211[k]
                   + f_4 * pc_z[k] * pgg_346[k];

        t_488[k] = f_0 * pfg_275[k]
                   + f_7 * pgf0_235[k]
                   - f_8 * pgf1_235[k]
                   + f_4 * pc_x[k] * pgg_350[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pb_x, pc_x, pc_z, pfh0_384, pfh0_386, pfg_213, \
                         pfg_276, pfg_278, pfh1_384, pfh1_386, \
                         pgg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = pb_x[k] * pfh0_384[k]
                   + f_10 * pfg_276[k]
                   - f_9 * pc_x[k] * pfh1_384[k];

        t_490[k] = f_10 * pfg_213[k]
                   + f_4 * pc_z[k] * pgg_348[k];

        t_491[k] = pb_x[k] * pfh0_386[k]
                   + f_10 * pfg_278[k]
                   - f_9 * pc_x[k] * pfh1_386[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pc_x, pfg_279, pfg_280, pfg_281, pfg_282, \
                         pgf0_239, pgf1_239, pgg_354, pgg_355, pgg_356, \
                         pgg_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_0 * pfg_279[k]
                   + f_5 * pgf0_239[k]
                   - f_6 * pgf1_239[k]
                   + f_4 * pc_x[k] * pgg_354[k];

        t_493[k] = f_0 * pfg_280[k]
                   + f_4 * pc_x[k] * pgg_355[k];

        t_494[k] = f_0 * pfg_281[k]
                   + f_4 * pc_x[k] * pgg_356[k];

        t_495[k] = f_0 * pfg_282[k]
                   + f_4 * pc_x[k] * pgg_357[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_x, pc_x, pc_z, pfh0_393, pfg_220, \
                         pfg_283, pfg_284, pfh1_393, pgg_355, pgg_358, \
                         pgg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_0 * pfg_283[k]
                   + f_4 * pc_x[k] * pgg_358[k];

        t_497[k] = f_0 * pfg_284[k]
                   + f_4 * pc_x[k] * pgg_359[k];

        t_498[k] = pb_x[k] * pfh0_393[k]
                   - f_9 * pc_x[k] * pfh1_393[k];

        t_499[k] = f_10 * pfg_220[k]
                   + f_4 * pc_z[k] * pgg_355[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pb_x, pc_x, pfh0_395, pfh0_396, pfh0_397, \
                         pfh0_398, pfh1_395, pfh1_396, pfh1_397, \
                         pfh1_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pb_x[k] * pfh0_395[k]
                   - f_9 * pc_x[k] * pfh1_395[k];

        t_501[k] = pb_x[k] * pfh0_396[k]
                   - f_9 * pc_x[k] * pfh1_396[k];

        t_502[k] = pb_x[k] * pfh0_397[k]
                   - f_9 * pc_x[k] * pfh1_397[k];

        t_503[k] = pb_x[k] * pfh0_398[k]
                   - f_9 * pc_x[k] * pfh1_398[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pa_y, pc_y, pc_z, sgh0_189, sgh0_190, sgg_135, \
                         sgh1_189, sgh1_190, pfg_225, pgg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = pa_y[k] * sgh0_189[k]
                   - f_9 * pc_y[k] * sgh1_189[k];

        t_505[k] = pa_y[k] * sgh0_190[k]
                   + f_0 * sgg_135[k]
                   - f_9 * pc_y[k] * sgh1_190[k];

        t_506[k] = f_11 * pfg_225[k]
                   + f_4 * pc_z[k] * pgg_360[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, pa_y, pc_y, pc_z, sgh0_192, sgh0_194, sgg_136, \
                         sgh1_192, sgh1_194, pfg_226, pgg_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = pa_y[k] * sgh0_192[k]
                   + f_10 * sgg_136[k]
                   - f_9 * pc_y[k] * sgh1_192[k];

        t_508[k] = f_11 * pfg_226[k]
                   + f_4 * pc_z[k] * pgg_361[k];

        t_509[k] = pa_y[k] * sgh0_194[k]
                   - f_9 * pc_y[k] * sgh1_194[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, pa_y, pc_y, pc_z, sgh0_195, sgh0_197, sgg_138, \
                         sgg_140, sgh1_195, sgh1_197, pfg_228, \
                         pgg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = pa_y[k] * sgh0_195[k]
                   + f_11 * sgg_138[k]
                   - f_9 * pc_y[k] * sgh1_195[k];

        t_511[k] = f_11 * pfg_228[k]
                   + f_4 * pc_z[k] * pgg_363[k];

        t_512[k] = pa_y[k] * sgh0_197[k]
                   + f_0 * sgg_140[k]
                   - f_9 * pc_y[k] * sgh1_197[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, t_516, pa_y, pc_x, pc_y, sgh0_198, sgh1_198, \
                         pfg_295, pfg_296, pfg_297, pgg_370, pgg_371, \
                         pgg_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = pa_y[k] * sgh0_198[k]
                   - f_9 * pc_y[k] * sgh1_198[k];

        t_514[k] = f_0 * pfg_295[k]
                   + f_4 * pc_x[k] * pgg_370[k];

        t_515[k] = f_0 * pfg_296[k]
                   + f_4 * pc_x[k] * pgg_371[k];

        t_516[k] = f_0 * pfg_297[k]
                   + f_4 * pc_x[k] * pgg_372[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pb_x, pc_x, pc_z, pfh0_414, pfg_235, \
                         pfg_298, pfg_299, pfh1_414, pgg_370, pgg_373, \
                         pgg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_0 * pfg_298[k]
                   + f_4 * pc_x[k] * pgg_373[k];

        t_518[k] = f_0 * pfg_299[k]
                   + f_4 * pc_x[k] * pgg_374[k];

        t_519[k] = pb_x[k] * pfh0_414[k]
                   - f_9 * pc_x[k] * pfh1_414[k];

        t_520[k] = f_11 * pfg_235[k]
                   + f_4 * pc_z[k] * pgg_370[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pa_y, pb_x, pc_x, pc_y, sgh0_209, \
                         sgg_149, sgh1_209, pfh0_416, pfh0_417, pfh1_416, pfh1_417, \
                         pgg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = pb_x[k] * pfh0_416[k]
                   - f_9 * pc_x[k] * pfh1_416[k];

        t_522[k] = pb_x[k] * pfh0_417[k]
                   - f_9 * pc_x[k] * pfh1_417[k];

        t_523[k] = f_0 * sgg_149[k]
                   + f_4 * pc_y[k] * pgg_374[k];

        t_524[k] = pa_y[k] * sgh0_209[k]
                   - f_9 * pc_y[k] * sgh1_209[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, pc_z, pgf0_250, pgf0_251, \
                         pgf0_253, pgf1_250, pgf1_251, pgf1_253, pgg_375, pgg_376, \
                         pgg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_2 * pgf0_250[k]
                   - f_3 * pgf1_250[k]
                   + f_4 * pc_x[k] * pgg_375[k];

        t_526[k] = f_15 * pgf0_251[k]
                   - f_16 * pgf1_251[k]
                   + f_4 * pc_x[k] * pgg_376[k];

        t_527[k] = f_4 * pc_z[k] * pgg_375[k];

        t_528[k] = f_7 * pgf0_253[k]
                   - f_8 * pgf1_253[k]
                   + f_4 * pc_x[k] * pgg_378[k];

        t_529[k] = f_4 * pc_z[k] * pgg_376[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pc_x, pc_z, pgf0_255, pgf0_256, pgf0_258, \
                         pgf1_255, pgf1_256, pgf1_258, pgg_378, pgg_380, pgg_381, \
                         pgg_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_7 * pgf0_255[k]
                   - f_8 * pgf1_255[k]
                   + f_4 * pc_x[k] * pgg_380[k];

        t_531[k] = f_5 * pgf0_256[k]
                   - f_6 * pgf1_256[k]
                   + f_4 * pc_x[k] * pgg_381[k];

        t_532[k] = f_4 * pc_z[k] * pgg_378[k];

        t_533[k] = f_5 * pgf0_258[k]
                   - f_6 * pgf1_258[k]
                   + f_4 * pc_x[k] * pgg_383[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, t_539, pc_x, pgf0_259, pgf1_259, \
                         pgg_384, pgg_385, pgg_386, pgg_387, pgg_388, \
                         pgg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_5 * pgf0_259[k]
                   - f_6 * pgf1_259[k]
                   + f_4 * pc_x[k] * pgg_384[k];

        t_535[k] = f_4 * pc_x[k] * pgg_385[k];

        t_536[k] = f_4 * pc_x[k] * pgg_386[k];

        t_537[k] = f_4 * pc_x[k] * pgg_387[k];

        t_538[k] = f_4 * pc_x[k] * pgg_388[k];

        t_539[k] = f_4 * pc_x[k] * pgg_389[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pc_y, pc_z, sgg_160, pfg_250, pgf0_256, \
                         pgf0_257, pgf1_256, pgf1_257, pgg_385, pgg_386, \
                         pgg_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_0 * sgg_160[k]
                   + f_1 * pfg_250[k]
                   + f_2 * pgf0_256[k]
                   - f_3 * pgf1_256[k]
                   + f_4 * pc_y[k] * pgg_385[k];

        t_541[k] = f_4 * pc_z[k] * pgg_385[k];

        t_542[k] = f_5 * pgf0_256[k]
                   - f_6 * pgf1_256[k]
                   + f_4 * pc_z[k] * pgg_386[k];

        t_543[k] = f_7 * pgf0_257[k]
                   - f_8 * pgf1_257[k]
                   + f_4 * pc_z[k] * pgg_387[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_z, pc_y, pc_z, sgg_164, pfh0_336, \
                         pfh0_337, pfg_254, pfh1_336, pfh1_337, pgf0_259, pgf1_259, \
                         pgg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_0 * sgg_164[k]
                   + f_1 * pfg_254[k]
                   + f_4 * pc_y[k] * pgg_389[k];

        t_545[k] = f_2 * pgf0_259[k]
                   - f_3 * pgf1_259[k]
                   + f_4 * pc_z[k] * pgg_389[k];

        t_546[k] = pb_z[k] * pfh0_336[k]
                   - f_9 * pc_z[k] * pfh1_336[k];

        t_547[k] = pb_z[k] * pfh0_337[k]
                   - f_9 * pc_z[k] * pfh1_337[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pb_z, pc_x, pc_z, pfh0_339, pfg_240, \
                         pfg_241, pfh1_339, pgf0_265, pgf1_265, pgg_390, pgg_391, \
                         pgg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_0 * pfg_240[k]
                   + f_4 * pc_z[k] * pgg_390[k];

        t_549[k] = pb_z[k] * pfh0_339[k]
                   - f_9 * pc_z[k] * pfh1_339[k];

        t_550[k] = f_0 * pfg_241[k]
                   + f_4 * pc_z[k] * pgg_391[k];

        t_551[k] = f_7 * pgf0_265[k]
                   - f_8 * pgf1_265[k]
                   + f_4 * pc_x[k] * pgg_395[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pb_z, pc_x, pc_z, pfh0_342, pfg_243, pfh1_342, \
                         pgf0_268, pgf1_268, pgg_393, pgg_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pb_z[k] * pfh0_342[k]
                   - f_9 * pc_z[k] * pfh1_342[k];

        t_553[k] = f_0 * pfg_243[k]
                   + f_4 * pc_z[k] * pgg_393[k];

        t_554[k] = f_5 * pgf0_268[k]
                   - f_6 * pgf1_268[k]
                   + f_4 * pc_x[k] * pgg_398[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, t_560, pc_x, pgf0_269, pgf1_269, \
                         pgg_399, pgg_400, pgg_401, pgg_402, pgg_403, \
                         pgg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_5 * pgf0_269[k]
                   - f_6 * pgf1_269[k]
                   + f_4 * pc_x[k] * pgg_399[k];

        t_556[k] = f_4 * pc_x[k] * pgg_400[k];

        t_557[k] = f_4 * pc_x[k] * pgg_401[k];

        t_558[k] = f_4 * pc_x[k] * pgg_402[k];

        t_559[k] = f_4 * pc_x[k] * pgg_403[k];

        t_560[k] = f_4 * pc_x[k] * pgg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, pb_z, pc_z, pfh0_351, pfh0_353, pfh0_354, \
                         pfg_250, pfg_251, pfg_252, pfh1_351, pfh1_353, pfh1_354, \
                         pgg_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = pb_z[k] * pfh0_351[k]
                   - f_9 * pc_z[k] * pfh1_351[k];

        t_562[k] = f_0 * pfg_250[k]
                   + f_4 * pc_z[k] * pgg_400[k];

        t_563[k] = pb_z[k] * pfh0_353[k]
                   + f_10 * pfg_251[k]
                   - f_9 * pc_z[k] * pfh1_353[k];

        t_564[k] = pb_z[k] * pfh0_354[k]
                   + f_11 * pfg_252[k]
                   - f_9 * pc_z[k] * pfh1_354[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, pc_x, pc_y, pc_z, sgg_179, pfg_254, pfg_269, \
                         pgf0_269, pgf0_270, pgf1_269, pgf1_270, pgg_404, \
                         pgg_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = f_0 * sgg_179[k]
                   + f_11 * pfg_269[k]
                   + f_4 * pc_y[k] * pgg_404[k];

        t_566[k] = f_0 * pfg_254[k]
                   + f_2 * pgf0_269[k]
                   - f_3 * pgf1_269[k]
                   + f_4 * pc_z[k] * pgg_404[k];

        t_567[k] = f_2 * pgf0_270[k]
                   - f_3 * pgf1_270[k]
                   + f_4 * pc_x[k] * pgg_405[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, pc_x, pc_z, pfg_255, pfg_256, pgf0_271, \
                         pgf0_273, pgf1_271, pgf1_273, pgg_405, pgg_406, \
                         pgg_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_15 * pgf0_271[k]
                   - f_16 * pgf1_271[k]
                   + f_4 * pc_x[k] * pgg_406[k];

        t_569[k] = f_10 * pfg_255[k]
                   + f_4 * pc_z[k] * pgg_405[k];

        t_570[k] = f_7 * pgf0_273[k]
                   - f_8 * pgf1_273[k]
                   + f_4 * pc_x[k] * pgg_408[k];

        t_571[k] = f_10 * pfg_256[k]
                   + f_4 * pc_z[k] * pgg_406[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pc_x, pc_z, pfg_258, pgf0_275, pgf0_276, \
                         pgf1_275, pgf1_276, pgg_408, pgg_410, \
                         pgg_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_7 * pgf0_275[k]
                   - f_8 * pgf1_275[k]
                   + f_4 * pc_x[k] * pgg_410[k];

        t_573[k] = f_5 * pgf0_276[k]
                   - f_6 * pgf1_276[k]
                   + f_4 * pc_x[k] * pgg_411[k];

        t_574[k] = f_10 * pfg_258[k]
                   + f_4 * pc_z[k] * pgg_408[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, pc_x, pgf0_278, pgf0_279, \
                         pgf1_278, pgf1_279, pgg_413, pgg_414, pgg_415, pgg_416, \
                         pgg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_5 * pgf0_278[k]
                   - f_6 * pgf1_278[k]
                   + f_4 * pc_x[k] * pgg_413[k];

        t_576[k] = f_5 * pgf0_279[k]
                   - f_6 * pgf1_279[k]
                   + f_4 * pc_x[k] * pgg_414[k];

        t_577[k] = f_4 * pc_x[k] * pgg_415[k];

        t_578[k] = f_4 * pc_x[k] * pgg_416[k];

        t_579[k] = f_4 * pc_x[k] * pgg_417[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pb_z, pc_x, pc_z, pdh0_204, pdh1_204, \
                         pfh0_372, pfg_265, pfh1_372, pgg_415, pgg_418, \
                         pgg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_4 * pc_x[k] * pgg_418[k];

        t_581[k] = f_4 * pc_x[k] * pgg_419[k];

        t_582[k] = f_12 * pdh0_204[k]
                   - f_13 * pdh1_204[k]
                   + pb_z[k] * pfh0_372[k]
                   - f_9 * pc_z[k] * pfh1_372[k];

        t_583[k] = f_10 * pfg_265[k]
                   + f_4 * pc_z[k] * pgg_415[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pc_y, pc_z, sgg_194, pfg_266, pfg_267, pfg_284, \
                         pgf0_276, pgf0_277, pgf1_276, pgf1_277, pgg_416, pgg_417, \
                         pgg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_10 * pfg_266[k]
                   + f_5 * pgf0_276[k]
                   - f_6 * pgf1_276[k]
                   + f_4 * pc_z[k] * pgg_416[k];

        t_585[k] = f_10 * pfg_267[k]
                   + f_7 * pgf0_277[k]
                   - f_8 * pgf1_277[k]
                   + f_4 * pc_z[k] * pgg_417[k];

        t_586[k] = f_0 * sgg_194[k]
                   + f_10 * pfg_284[k]
                   + f_4 * pc_y[k] * pgg_419[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pc_x, pc_z, pfg_269, pgf0_279, pgf0_280, \
                         pgf0_281, pgf1_279, pgf1_280, pgf1_281, pgg_419, pgg_420, \
                         pgg_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_10 * pfg_269[k]
                   + f_2 * pgf0_279[k]
                   - f_3 * pgf1_279[k]
                   + f_4 * pc_z[k] * pgg_419[k];

        t_588[k] = f_2 * pgf0_280[k]
                   - f_3 * pgf1_280[k]
                   + f_4 * pc_x[k] * pgg_420[k];

        t_589[k] = f_15 * pgf0_281[k]
                   - f_16 * pgf1_281[k]
                   + f_4 * pc_x[k] * pgg_421[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pc_x, pc_z, pfg_270, pfg_271, pgf0_283, \
                         pgf0_285, pgf1_283, pgf1_285, pgg_420, pgg_421, pgg_423, \
                         pgg_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_11 * pfg_270[k]
                   + f_4 * pc_z[k] * pgg_420[k];

        t_591[k] = f_7 * pgf0_283[k]
                   - f_8 * pgf1_283[k]
                   + f_4 * pc_x[k] * pgg_423[k];

        t_592[k] = f_11 * pfg_271[k]
                   + f_4 * pc_z[k] * pgg_421[k];

        t_593[k] = f_7 * pgf0_285[k]
                   - f_8 * pgf1_285[k]
                   + f_4 * pc_x[k] * pgg_425[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pc_x, pc_z, pfg_273, pgf0_286, pgf0_288, \
                         pgf1_286, pgf1_288, pgg_423, pgg_426, \
                         pgg_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_5 * pgf0_286[k]
                   - f_6 * pgf1_286[k]
                   + f_4 * pc_x[k] * pgg_426[k];

        t_595[k] = f_11 * pfg_273[k]
                   + f_4 * pc_z[k] * pgg_423[k];

        t_596[k] = f_5 * pgf0_288[k]
                   - f_6 * pgf1_288[k]
                   + f_4 * pc_x[k] * pgg_428[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, t_602, pc_x, pgf0_289, pgf1_289, \
                         pgg_429, pgg_430, pgg_431, pgg_432, pgg_433, \
                         pgg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = f_5 * pgf0_289[k]
                   - f_6 * pgf1_289[k]
                   + f_4 * pc_x[k] * pgg_429[k];

        t_598[k] = f_4 * pc_x[k] * pgg_430[k];

        t_599[k] = f_4 * pc_x[k] * pgg_431[k];

        t_600[k] = f_4 * pc_x[k] * pgg_432[k];

        t_601[k] = f_4 * pc_x[k] * pgg_433[k];

        t_602[k] = f_4 * pc_x[k] * pgg_434[k];
    }
}

static auto
compute_prim_pgh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgh0, const size_t sgg,
                                                          const size_t sgh1, const size_t pdh0,
                                                          const size_t pdh1, const size_t pfh0,
                                                          const size_t pfg, const size_t pfh1,
                                                          const size_t pgf0, const size_t pgf1,
                                                          const size_t pgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);
    const auto f_17 = 1.0 / p;
    const auto f_18 = gamma / (p * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_0 = buffer.data(sgh0 + 0);
    const auto *sgh0_2 = buffer.data(sgh0 + 2);
    const auto *sgh0_3 = buffer.data(sgh0 + 3);
    const auto *sgh0_5 = buffer.data(sgh0 + 5);
    const auto *sgh0_6 = buffer.data(sgh0 + 6);
    const auto *sgh0_7 = buffer.data(sgh0 + 7);
    const auto *sgh0_9 = buffer.data(sgh0 + 9);
    const auto *sgh0_15 = buffer.data(sgh0 + 15);
    const auto *sgh0_20 = buffer.data(sgh0 + 20);
    const auto *sgh0_21 = buffer.data(sgh0 + 21);
    const auto *sgh0_24 = buffer.data(sgh0 + 24);
    const auto *sgh0_27 = buffer.data(sgh0 + 27);
    const auto *sgh0_28 = buffer.data(sgh0 + 28);
    const auto *sgh0_36 = buffer.data(sgh0 + 36);
    const auto *sgh0_37 = buffer.data(sgh0 + 37);
    const auto *sgh0_38 = buffer.data(sgh0 + 38);
    const auto *sgh0_39 = buffer.data(sgh0 + 39);
    const auto *sgh0_63 = buffer.data(sgh0 + 63);
    const auto *sgh0_65 = buffer.data(sgh0 + 65);
    const auto *sgh0_66 = buffer.data(sgh0 + 66);
    const auto *sgh0_68 = buffer.data(sgh0 + 68);
    const auto *sgh0_69 = buffer.data(sgh0 + 69);
    const auto *sgh0_70 = buffer.data(sgh0 + 70);
    const auto *sgh0_72 = buffer.data(sgh0 + 72);
    const auto *sgh0_78 = buffer.data(sgh0 + 78);
    const auto *sgh0_79 = buffer.data(sgh0 + 79);
    const auto *sgh0_80 = buffer.data(sgh0 + 80);
    const auto *sgh0_81 = buffer.data(sgh0 + 81);
    const auto *sgh0_83 = buffer.data(sgh0 + 83);
    const auto *sgh0_294 = buffer.data(sgh0 + 294);
    const auto *sgh0_299 = buffer.data(sgh0 + 299);
    const auto *sgh0_303 = buffer.data(sgh0 + 303);
    const auto *sgh0_309 = buffer.data(sgh0 + 309);
    const auto *sgh0_311 = buffer.data(sgh0 + 311);
    const auto *sgh0_312 = buffer.data(sgh0 + 312);
    const auto *sgh0_314 = buffer.data(sgh0 + 314);

    const auto *sgg_0 = buffer.data(sgg + 0);
    const auto *sgg_2 = buffer.data(sgg + 2);
    const auto *sgg_3 = buffer.data(sgg + 3);
    const auto *sgg_5 = buffer.data(sgg + 5);
    const auto *sgg_14 = buffer.data(sgg + 14);
    const auto *sgg_18 = buffer.data(sgg + 18);
    const auto *sgg_25 = buffer.data(sgg + 25);
    const auto *sgg_26 = buffer.data(sgg + 26);
    const auto *sgg_27 = buffer.data(sgg + 27);
    const auto *sgg_45 = buffer.data(sgg + 45);
    const auto *sgg_47 = buffer.data(sgg + 47);
    const auto *sgg_48 = buffer.data(sgg + 48);
    const auto *sgg_50 = buffer.data(sgg + 50);
    const auto *sgg_55 = buffer.data(sgg + 55);
    const auto *sgg_56 = buffer.data(sgg + 56);
    const auto *sgg_57 = buffer.data(sgg + 57);
    const auto *sgg_59 = buffer.data(sgg + 59);
    const auto *sgg_205 = buffer.data(sgg + 205);
    const auto *sgg_209 = buffer.data(sgg + 209);
    const auto *sgg_220 = buffer.data(sgg + 220);
    const auto *sgg_222 = buffer.data(sgg + 222);
    const auto *sgg_223 = buffer.data(sgg + 223);
    const auto *sgg_224 = buffer.data(sgg + 224);

    const auto *sgh1_0 = buffer.data(sgh1 + 0);
    const auto *sgh1_2 = buffer.data(sgh1 + 2);
    const auto *sgh1_3 = buffer.data(sgh1 + 3);
    const auto *sgh1_5 = buffer.data(sgh1 + 5);
    const auto *sgh1_6 = buffer.data(sgh1 + 6);
    const auto *sgh1_7 = buffer.data(sgh1 + 7);
    const auto *sgh1_9 = buffer.data(sgh1 + 9);
    const auto *sgh1_15 = buffer.data(sgh1 + 15);
    const auto *sgh1_20 = buffer.data(sgh1 + 20);
    const auto *sgh1_21 = buffer.data(sgh1 + 21);
    const auto *sgh1_24 = buffer.data(sgh1 + 24);
    const auto *sgh1_27 = buffer.data(sgh1 + 27);
    const auto *sgh1_28 = buffer.data(sgh1 + 28);
    const auto *sgh1_36 = buffer.data(sgh1 + 36);
    const auto *sgh1_37 = buffer.data(sgh1 + 37);
    const auto *sgh1_38 = buffer.data(sgh1 + 38);
    const auto *sgh1_39 = buffer.data(sgh1 + 39);
    const auto *sgh1_63 = buffer.data(sgh1 + 63);
    const auto *sgh1_65 = buffer.data(sgh1 + 65);
    const auto *sgh1_66 = buffer.data(sgh1 + 66);
    const auto *sgh1_68 = buffer.data(sgh1 + 68);
    const auto *sgh1_69 = buffer.data(sgh1 + 69);
    const auto *sgh1_70 = buffer.data(sgh1 + 70);
    const auto *sgh1_72 = buffer.data(sgh1 + 72);
    const auto *sgh1_78 = buffer.data(sgh1 + 78);
    const auto *sgh1_79 = buffer.data(sgh1 + 79);
    const auto *sgh1_80 = buffer.data(sgh1 + 80);
    const auto *sgh1_81 = buffer.data(sgh1 + 81);
    const auto *sgh1_83 = buffer.data(sgh1 + 83);
    const auto *sgh1_294 = buffer.data(sgh1 + 294);
    const auto *sgh1_299 = buffer.data(sgh1 + 299);
    const auto *sgh1_303 = buffer.data(sgh1 + 303);
    const auto *sgh1_309 = buffer.data(sgh1 + 309);
    const auto *sgh1_311 = buffer.data(sgh1 + 311);
    const auto *sgh1_312 = buffer.data(sgh1 + 312);
    const auto *sgh1_314 = buffer.data(sgh1 + 314);

    const auto *pdh0_314 = buffer.data(pdh0 + 314);

    const auto *pdh1_314 = buffer.data(pdh1 + 314);

    const auto *pfh0_422 = buffer.data(pfh0 + 422);
    const auto *pfh0_425 = buffer.data(pfh0 + 425);
    const auto *pfh0_429 = buffer.data(pfh0 + 429);
    const auto *pfh0_440 = buffer.data(pfh0 + 440);
    const auto *pfh0_462 = buffer.data(pfh0 + 462);
    const auto *pfh0_464 = buffer.data(pfh0 + 464);
    const auto *pfh0_467 = buffer.data(pfh0 + 467);
    const auto *pfh0_482 = buffer.data(pfh0 + 482);

    const auto *pfg_280 = buffer.data(pfg + 280);
    const auto *pfg_281 = buffer.data(pfg + 281);
    const auto *pfg_282 = buffer.data(pfg + 282);
    const auto *pfg_284 = buffer.data(pfg + 284);
    const auto *pfg_285 = buffer.data(pfg + 285);
    const auto *pfg_286 = buffer.data(pfg + 286);
    const auto *pfg_288 = buffer.data(pfg + 288);
    const auto *pfg_295 = buffer.data(pfg + 295);
    const auto *pfg_299 = buffer.data(pfg + 299);
    const auto *pfg_300 = buffer.data(pfg + 300);
    const auto *pfg_302 = buffer.data(pfg + 302);
    const auto *pfg_305 = buffer.data(pfg + 305);
    const auto *pfg_310 = buffer.data(pfg + 310);
    const auto *pfg_311 = buffer.data(pfg + 311);
    const auto *pfg_312 = buffer.data(pfg + 312);
    const auto *pfg_313 = buffer.data(pfg + 313);
    const auto *pfg_314 = buffer.data(pfg + 314);
    const auto *pfg_315 = buffer.data(pfg + 315);
    const auto *pfg_317 = buffer.data(pfg + 317);
    const auto *pfg_320 = buffer.data(pfg + 320);
    const auto *pfg_325 = buffer.data(pfg + 325);
    const auto *pfg_326 = buffer.data(pfg + 326);
    const auto *pfg_327 = buffer.data(pfg + 327);
    const auto *pfg_328 = buffer.data(pfg + 328);
    const auto *pfg_329 = buffer.data(pfg + 329);
    const auto *pfg_330 = buffer.data(pfg + 330);
    const auto *pfg_332 = buffer.data(pfg + 332);
    const auto *pfg_333 = buffer.data(pfg + 333);
    const auto *pfg_335 = buffer.data(pfg + 335);
    const auto *pfg_336 = buffer.data(pfg + 336);
    const auto *pfg_337 = buffer.data(pfg + 337);
    const auto *pfg_339 = buffer.data(pfg + 339);
    const auto *pfg_340 = buffer.data(pfg + 340);
    const auto *pfg_341 = buffer.data(pfg + 341);
    const auto *pfg_342 = buffer.data(pfg + 342);
    const auto *pfg_343 = buffer.data(pfg + 343);
    const auto *pfg_344 = buffer.data(pfg + 344);
    const auto *pfg_355 = buffer.data(pfg + 355);
    const auto *pfg_356 = buffer.data(pfg + 356);
    const auto *pfg_357 = buffer.data(pfg + 357);
    const auto *pfg_358 = buffer.data(pfg + 358);
    const auto *pfg_359 = buffer.data(pfg + 359);
    const auto *pfg_363 = buffer.data(pfg + 363);

    const auto *pfh1_422 = buffer.data(pfh1 + 422);
    const auto *pfh1_425 = buffer.data(pfh1 + 425);
    const auto *pfh1_429 = buffer.data(pfh1 + 429);
    const auto *pfh1_440 = buffer.data(pfh1 + 440);
    const auto *pfh1_462 = buffer.data(pfh1 + 462);
    const auto *pfh1_464 = buffer.data(pfh1 + 464);
    const auto *pfh1_467 = buffer.data(pfh1 + 467);
    const auto *pfh1_482 = buffer.data(pfh1 + 482);

    const auto *pgf0_286 = buffer.data(pgf0 + 286);
    const auto *pgf0_287 = buffer.data(pgf0 + 287);
    const auto *pgf0_289 = buffer.data(pgf0 + 289);
    const auto *pgf0_291 = buffer.data(pgf0 + 291);
    const auto *pgf0_293 = buffer.data(pgf0 + 293);
    const auto *pgf0_296 = buffer.data(pgf0 + 296);
    const auto *pgf0_298 = buffer.data(pgf0 + 298);
    const auto *pgf0_307 = buffer.data(pgf0 + 307);
    const auto *pgf0_308 = buffer.data(pgf0 + 308);
    const auto *pgf0_309 = buffer.data(pgf0 + 309);
    const auto *pgf0_320 = buffer.data(pgf0 + 320);
    const auto *pgf0_322 = buffer.data(pgf0 + 322);
    const auto *pgf0_323 = buffer.data(pgf0 + 323);
    const auto *pgf0_325 = buffer.data(pgf0 + 325);
    const auto *pgf0_326 = buffer.data(pgf0 + 326);
    const auto *pgf0_327 = buffer.data(pgf0 + 327);
    const auto *pgf0_328 = buffer.data(pgf0 + 328);
    const auto *pgf0_329 = buffer.data(pgf0 + 329);
    const auto *pgf0_343 = buffer.data(pgf0 + 343);

    const auto *pgf1_286 = buffer.data(pgf1 + 286);
    const auto *pgf1_287 = buffer.data(pgf1 + 287);
    const auto *pgf1_289 = buffer.data(pgf1 + 289);
    const auto *pgf1_291 = buffer.data(pgf1 + 291);
    const auto *pgf1_293 = buffer.data(pgf1 + 293);
    const auto *pgf1_296 = buffer.data(pgf1 + 296);
    const auto *pgf1_298 = buffer.data(pgf1 + 298);
    const auto *pgf1_307 = buffer.data(pgf1 + 307);
    const auto *pgf1_308 = buffer.data(pgf1 + 308);
    const auto *pgf1_309 = buffer.data(pgf1 + 309);
    const auto *pgf1_320 = buffer.data(pgf1 + 320);
    const auto *pgf1_322 = buffer.data(pgf1 + 322);
    const auto *pgf1_323 = buffer.data(pgf1 + 323);
    const auto *pgf1_325 = buffer.data(pgf1 + 325);
    const auto *pgf1_326 = buffer.data(pgf1 + 326);
    const auto *pgf1_327 = buffer.data(pgf1 + 327);
    const auto *pgf1_328 = buffer.data(pgf1 + 328);
    const auto *pgf1_329 = buffer.data(pgf1 + 329);
    const auto *pgf1_343 = buffer.data(pgf1 + 343);

    const auto *pgg_430 = buffer.data(pgg + 430);
    const auto *pgg_431 = buffer.data(pgg + 431);
    const auto *pgg_432 = buffer.data(pgg + 432);
    const auto *pgg_434 = buffer.data(pgg + 434);
    const auto *pgg_435 = buffer.data(pgg + 435);
    const auto *pgg_436 = buffer.data(pgg + 436);
    const auto *pgg_438 = buffer.data(pgg + 438);
    const auto *pgg_441 = buffer.data(pgg + 441);
    const auto *pgg_443 = buffer.data(pgg + 443);
    const auto *pgg_445 = buffer.data(pgg + 445);
    const auto *pgg_446 = buffer.data(pgg + 446);
    const auto *pgg_447 = buffer.data(pgg + 447);
    const auto *pgg_448 = buffer.data(pgg + 448);
    const auto *pgg_449 = buffer.data(pgg + 449);
    const auto *pgg_450 = buffer.data(pgg + 450);
    const auto *pgg_452 = buffer.data(pgg + 452);
    const auto *pgg_455 = buffer.data(pgg + 455);
    const auto *pgg_460 = buffer.data(pgg + 460);
    const auto *pgg_461 = buffer.data(pgg + 461);
    const auto *pgg_462 = buffer.data(pgg + 462);
    const auto *pgg_463 = buffer.data(pgg + 463);
    const auto *pgg_464 = buffer.data(pgg + 464);
    const auto *pgg_465 = buffer.data(pgg + 465);
    const auto *pgg_467 = buffer.data(pgg + 467);
    const auto *pgg_470 = buffer.data(pgg + 470);
    const auto *pgg_475 = buffer.data(pgg + 475);
    const auto *pgg_476 = buffer.data(pgg + 476);
    const auto *pgg_477 = buffer.data(pgg + 477);
    const auto *pgg_478 = buffer.data(pgg + 478);
    const auto *pgg_479 = buffer.data(pgg + 479);
    const auto *pgg_480 = buffer.data(pgg + 480);
    const auto *pgg_482 = buffer.data(pgg + 482);
    const auto *pgg_483 = buffer.data(pgg + 483);
    const auto *pgg_485 = buffer.data(pgg + 485);
    const auto *pgg_486 = buffer.data(pgg + 486);
    const auto *pgg_487 = buffer.data(pgg + 487);
    const auto *pgg_489 = buffer.data(pgg + 489);
    const auto *pgg_490 = buffer.data(pgg + 490);
    const auto *pgg_491 = buffer.data(pgg + 491);
    const auto *pgg_492 = buffer.data(pgg + 492);
    const auto *pgg_493 = buffer.data(pgg + 493);
    const auto *pgg_494 = buffer.data(pgg + 494);
    const auto *pgg_495 = buffer.data(pgg + 495);
    const auto *pgg_497 = buffer.data(pgg + 497);
    const auto *pgg_500 = buffer.data(pgg + 500);
    const auto *pgg_505 = buffer.data(pgg + 505);
    const auto *pgg_506 = buffer.data(pgg + 506);
    const auto *pgg_507 = buffer.data(pgg + 507);
    const auto *pgg_508 = buffer.data(pgg + 508);
    const auto *pgg_509 = buffer.data(pgg + 509);
    const auto *pgg_510 = buffer.data(pgg + 510);
    const auto *pgg_512 = buffer.data(pgg + 512);
    const auto *pgg_513 = buffer.data(pgg + 513);

#pragma omp simd aligned(t_603, t_604, t_605, pc_y, pc_z, sgg_205, pfg_280, pfg_281, pfg_295, \
                         pgf0_286, pgf1_286, pgg_430, pgg_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = f_0 * sgg_205[k]
                   + f_0 * pfg_295[k]
                   + f_2 * pgf0_286[k]
                   - f_3 * pgf1_286[k]
                   + f_4 * pc_y[k] * pgg_430[k];

        t_604[k] = f_11 * pfg_280[k]
                   + f_4 * pc_z[k] * pgg_430[k];

        t_605[k] = f_11 * pfg_281[k]
                   + f_5 * pgf0_286[k]
                   - f_6 * pgf1_286[k]
                   + f_4 * pc_z[k] * pgg_431[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, pc_y, pc_z, sgg_209, pfg_282, pfg_284, pfg_299, \
                         pgf0_287, pgf0_289, pgf1_287, pgf1_289, pgg_432, \
                         pgg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = f_11 * pfg_282[k]
                   + f_7 * pgf0_287[k]
                   - f_8 * pgf1_287[k]
                   + f_4 * pc_z[k] * pgg_432[k];

        t_607[k] = f_0 * sgg_209[k]
                   + f_0 * pfg_299[k]
                   + f_4 * pc_y[k] * pgg_434[k];

        t_608[k] = f_11 * pfg_284[k]
                   + f_2 * pgf0_289[k]
                   - f_3 * pgf1_289[k]
                   + f_4 * pc_z[k] * pgg_434[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pa_y, pc_x, pc_y, pc_z, sgh0_294, sgh1_294, \
                         pfg_285, pgf0_291, pgf1_291, pgg_435, \
                         pgg_436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = pa_y[k] * sgh0_294[k]
                   - f_9 * pc_y[k] * sgh1_294[k];

        t_610[k] = f_15 * pgf0_291[k]
                   - f_16 * pgf1_291[k]
                   + f_4 * pc_x[k] * pgg_436[k];

        t_611[k] = f_1 * pfg_285[k]
                   + f_4 * pc_z[k] * pgg_435[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pa_y, pc_x, pc_y, pc_z, sgh0_299, sgh1_299, \
                         pfg_286, pgf0_293, pgf1_293, pgg_436, \
                         pgg_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_7 * pgf0_293[k]
                   - f_8 * pgf1_293[k]
                   + f_4 * pc_x[k] * pgg_438[k];

        t_613[k] = f_1 * pfg_286[k]
                   + f_4 * pc_z[k] * pgg_436[k];

        t_614[k] = pa_y[k] * sgh0_299[k]
                   - f_9 * pc_y[k] * sgh1_299[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pc_x, pc_z, pfg_288, pgf0_296, pgf0_298, \
                         pgf1_296, pgf1_298, pgg_438, pgg_441, \
                         pgg_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_5 * pgf0_296[k]
                   - f_6 * pgf1_296[k]
                   + f_4 * pc_x[k] * pgg_441[k];

        t_616[k] = f_1 * pfg_288[k]
                   + f_4 * pc_z[k] * pgg_438[k];

        t_617[k] = f_5 * pgf0_298[k]
                   - f_6 * pgf1_298[k]
                   + f_4 * pc_x[k] * pgg_443[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, t_622, t_623, pa_y, pc_x, pc_y, sgh0_303, \
                         sgh1_303, pgg_445, pgg_446, pgg_447, pgg_448, \
                         pgg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = pa_y[k] * sgh0_303[k]
                   - f_9 * pc_y[k] * sgh1_303[k];

        t_619[k] = f_4 * pc_x[k] * pgg_445[k];

        t_620[k] = f_4 * pc_x[k] * pgg_446[k];

        t_621[k] = f_4 * pc_x[k] * pgg_447[k];

        t_622[k] = f_4 * pc_x[k] * pgg_448[k];

        t_623[k] = f_4 * pc_x[k] * pgg_449[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_y, pc_y, pc_z, sgh0_309, sgh0_311, sgg_220, \
                         sgg_222, sgh1_309, sgh1_311, pfg_295, \
                         pgg_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = pa_y[k] * sgh0_309[k]
                   + f_14 * sgg_220[k]
                   - f_9 * pc_y[k] * sgh1_309[k];

        t_625[k] = f_1 * pfg_295[k]
                   + f_4 * pc_z[k] * pgg_445[k];

        t_626[k] = pa_y[k] * sgh0_311[k]
                   + f_11 * sgg_222[k]
                   - f_9 * pc_y[k] * sgh1_311[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_y, pc_y, sgh0_312, sgh0_314, sgg_223, \
                         sgg_224, sgh1_312, sgh1_314, pgg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = pa_y[k] * sgh0_312[k]
                   + f_10 * sgg_223[k]
                   - f_9 * pc_y[k] * sgh1_312[k];

        t_628[k] = f_0 * sgg_224[k]
                   + f_4 * pc_y[k] * pgg_449[k];

        t_629[k] = pa_y[k] * sgh0_314[k]
                   - f_9 * pc_y[k] * sgh1_314[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pa_z, pc_y, pc_z, sgh0_0, sgh0_2, sgh0_3, \
                         sgg_0, sgh1_0, sgh1_2, sgh1_3, pgg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = pa_z[k] * sgh0_0[k]
                   - f_9 * pc_z[k] * sgh1_0[k];

        t_631[k] = f_4 * pc_y[k] * pgg_450[k];

        t_632[k] = pa_z[k] * sgh0_2[k]
                   + f_0 * sgg_0[k]
                   - f_9 * pc_z[k] * sgh1_2[k];

        t_633[k] = pa_z[k] * sgh0_3[k]
                   - f_9 * pc_z[k] * sgh1_3[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pa_z, pc_y, pc_z, sgh0_5, sgh0_6, sgh0_7, \
                         sgg_2, sgg_3, sgh1_5, sgh1_6, sgh1_7, \
                         pgg_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_4 * pc_y[k] * pgg_452[k];

        t_635[k] = pa_z[k] * sgh0_5[k]
                   + f_10 * sgg_2[k]
                   - f_9 * pc_z[k] * sgh1_5[k];

        t_636[k] = pa_z[k] * sgh0_6[k]
                   - f_9 * pc_z[k] * sgh1_6[k];

        t_637[k] = pa_z[k] * sgh0_7[k]
                   + f_0 * sgg_3[k]
                   - f_9 * pc_z[k] * sgh1_7[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pa_z, pc_x, pc_y, pc_z, sgh0_9, sgg_5, \
                         sgh1_9, pfg_310, pfg_311, pgg_455, pgg_460, \
                         pgg_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_4 * pc_y[k] * pgg_455[k];

        t_639[k] = pa_z[k] * sgh0_9[k]
                   + f_11 * sgg_5[k]
                   - f_9 * pc_z[k] * sgh1_9[k];

        t_640[k] = f_1 * pfg_310[k]
                   + f_4 * pc_x[k] * pgg_460[k];

        t_641[k] = f_1 * pfg_311[k]
                   + f_4 * pc_x[k] * pgg_461[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pa_z, pc_x, pc_z, sgh0_15, sgh1_15, \
                         pfg_312, pfg_313, pfg_314, pgg_462, pgg_463, \
                         pgg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_1 * pfg_312[k]
                   + f_4 * pc_x[k] * pgg_462[k];

        t_643[k] = f_1 * pfg_313[k]
                   + f_4 * pc_x[k] * pgg_463[k];

        t_644[k] = f_1 * pfg_314[k]
                   + f_4 * pc_x[k] * pgg_464[k];

        t_645[k] = pa_z[k] * sgh0_15[k]
                   - f_9 * pc_z[k] * sgh1_15[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pc_y, pgf0_307, pgf0_308, pgf0_309, \
                         pgf1_307, pgf1_308, pgf1_309, pgg_461, pgg_462, pgg_463, \
                         pgg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * pgf0_307[k]
                   - f_16 * pgf1_307[k]
                   + f_4 * pc_y[k] * pgg_461[k];

        t_647[k] = f_7 * pgf0_308[k]
                   - f_8 * pgf1_308[k]
                   + f_4 * pc_y[k] * pgg_462[k];

        t_648[k] = f_5 * pgf0_309[k]
                   - f_6 * pgf1_309[k]
                   + f_4 * pc_y[k] * pgg_463[k];

        t_649[k] = f_4 * pc_y[k] * pgg_464[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pa_z, pc_y, pc_z, sgh0_20, sgh0_21, sgg_14, \
                         sgh1_20, sgh1_21, pfg_300, pgg_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = pa_z[k] * sgh0_20[k]
                   + f_14 * sgg_14[k]
                   - f_9 * pc_z[k] * sgh1_20[k];

        t_651[k] = pa_z[k] * sgh0_21[k]
                   - f_9 * pc_z[k] * sgh1_21[k];

        t_652[k] = f_0 * pfg_300[k]
                   + f_4 * pc_y[k] * pgg_465[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, t_656, pa_z, pb_y, pc_y, pc_z, sgh0_24, sgh1_24, \
                         pfh0_422, pfh0_425, pfg_302, pfh1_422, pfh1_425, \
                         pgg_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = pb_y[k] * pfh0_422[k]
                   - f_9 * pc_y[k] * pfh1_422[k];

        t_654[k] = pa_z[k] * sgh0_24[k]
                   - f_9 * pc_z[k] * sgh1_24[k];

        t_655[k] = f_0 * pfg_302[k]
                   + f_4 * pc_y[k] * pgg_467[k];

        t_656[k] = pb_y[k] * pfh0_425[k]
                   - f_9 * pc_y[k] * pfh1_425[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pa_z, pc_y, pc_z, sgh0_27, sgh0_28, sgg_18, \
                         sgh1_27, sgh1_28, pfg_305, pgg_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = pa_z[k] * sgh0_27[k]
                   - f_9 * pc_z[k] * sgh1_27[k];

        t_658[k] = pa_z[k] * sgh0_28[k]
                   + f_0 * sgg_18[k]
                   - f_9 * pc_z[k] * sgh1_28[k];

        t_659[k] = f_0 * pfg_305[k]
                   + f_4 * pc_y[k] * pgg_470[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, pb_y, pc_x, pc_y, pfh0_429, pfg_325, \
                         pfg_326, pfg_327, pfh1_429, pgg_475, pgg_476, \
                         pgg_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = pb_y[k] * pfh0_429[k]
                   - f_9 * pc_y[k] * pfh1_429[k];

        t_661[k] = f_11 * pfg_325[k]
                   + f_4 * pc_x[k] * pgg_475[k];

        t_662[k] = f_11 * pfg_326[k]
                   + f_4 * pc_x[k] * pgg_476[k];

        t_663[k] = f_11 * pfg_327[k]
                   + f_4 * pc_x[k] * pgg_477[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, pa_z, pc_x, pc_z, sgh0_36, sgh0_37, \
                         sgg_25, sgh1_36, sgh1_37, pfg_328, pfg_329, pgg_478, \
                         pgg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_11 * pfg_328[k]
                   + f_4 * pc_x[k] * pgg_478[k];

        t_665[k] = f_11 * pfg_329[k]
                   + f_4 * pc_x[k] * pgg_479[k];

        t_666[k] = pa_z[k] * sgh0_36[k]
                   - f_9 * pc_z[k] * sgh1_36[k];

        t_667[k] = pa_z[k] * sgh0_37[k]
                   + f_0 * sgg_25[k]
                   - f_9 * pc_z[k] * sgh1_37[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, pa_z, pc_y, pc_z, sgh0_38, sgh0_39, sgg_26, \
                         sgg_27, sgh1_38, sgh1_39, pfg_314, pgg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = pa_z[k] * sgh0_38[k]
                   + f_10 * sgg_26[k]
                   - f_9 * pc_z[k] * sgh1_38[k];

        t_669[k] = pa_z[k] * sgh0_39[k]
                   + f_11 * sgg_27[k]
                   - f_9 * pc_z[k] * sgh1_39[k];

        t_670[k] = f_0 * pfg_314[k]
                   + f_4 * pc_y[k] * pgg_479[k];
    }

#pragma omp simd aligned(t_671, t_672, t_673, pb_y, pc_x, pc_y, pfh0_440, pfg_330, pfh1_440, \
                         pgf0_320, pgf1_320, pgg_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_671[k] = pb_y[k] * pfh0_440[k]
                   - f_9 * pc_y[k] * pfh1_440[k];

        t_672[k] = f_11 * pfg_330[k]
                   + f_2 * pgf0_320[k]
                   - f_3 * pgf1_320[k]
                   + f_4 * pc_x[k] * pgg_480[k];

        t_673[k] = f_4 * pc_y[k] * pgg_480[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, pc_x, pc_y, pfg_332, pfg_333, pgf0_322, \
                         pgf0_323, pgf1_322, pgf1_323, pgg_482, \
                         pgg_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_11 * pfg_332[k]
                   + f_15 * pgf0_322[k]
                   - f_16 * pgf1_322[k]
                   + f_4 * pc_x[k] * pgg_482[k];

        t_675[k] = f_11 * pfg_333[k]
                   + f_7 * pgf0_323[k]
                   - f_8 * pgf1_323[k]
                   + f_4 * pc_x[k] * pgg_483[k];

        t_676[k] = f_4 * pc_y[k] * pgg_482[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pfg_335, pfg_336, pfg_337, pgf0_325, \
                         pgf0_326, pgf0_327, pgf1_325, pgf1_326, pgf1_327, pgg_485, pgg_486, \
                         pgg_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_11 * pfg_335[k]
                   + f_7 * pgf0_325[k]
                   - f_8 * pgf1_325[k]
                   + f_4 * pc_x[k] * pgg_485[k];

        t_678[k] = f_11 * pfg_336[k]
                   + f_5 * pgf0_326[k]
                   - f_6 * pgf1_326[k]
                   + f_4 * pc_x[k] * pgg_486[k];

        t_679[k] = f_11 * pfg_337[k]
                   + f_5 * pgf0_327[k]
                   - f_6 * pgf1_327[k]
                   + f_4 * pc_x[k] * pgg_487[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_y, pfg_339, pfg_340, pfg_341, \
                         pgf0_329, pgf1_329, pgg_485, pgg_489, pgg_490, \
                         pgg_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_4 * pc_y[k] * pgg_485[k];

        t_681[k] = f_11 * pfg_339[k]
                   + f_5 * pgf0_329[k]
                   - f_6 * pgf1_329[k]
                   + f_4 * pc_x[k] * pgg_489[k];

        t_682[k] = f_11 * pfg_340[k]
                   + f_4 * pc_x[k] * pgg_490[k];

        t_683[k] = f_11 * pfg_341[k]
                   + f_4 * pc_x[k] * pgg_491[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, pc_x, pc_y, pfg_342, pfg_343, pfg_344, \
                         pgf0_326, pgf1_326, pgg_490, pgg_492, pgg_493, \
                         pgg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_11 * pfg_342[k]
                   + f_4 * pc_x[k] * pgg_492[k];

        t_685[k] = f_11 * pfg_343[k]
                   + f_4 * pc_x[k] * pgg_493[k];

        t_686[k] = f_11 * pfg_344[k]
                   + f_4 * pc_x[k] * pgg_494[k];

        t_687[k] = f_2 * pgf0_326[k]
                   - f_3 * pgf1_326[k]
                   + f_4 * pc_y[k] * pgg_490[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, t_691, pc_y, pgf0_327, pgf0_328, pgf0_329, \
                         pgf1_327, pgf1_328, pgf1_329, pgg_491, pgg_492, pgg_493, \
                         pgg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_15 * pgf0_327[k]
                   - f_16 * pgf1_327[k]
                   + f_4 * pc_y[k] * pgg_491[k];

        t_689[k] = f_7 * pgf0_328[k]
                   - f_8 * pgf1_328[k]
                   + f_4 * pc_y[k] * pgg_492[k];

        t_690[k] = f_5 * pgf0_329[k]
                   - f_6 * pgf1_329[k]
                   + f_4 * pc_y[k] * pgg_493[k];

        t_691[k] = f_4 * pc_y[k] * pgg_494[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, pa_z, pb_x, pc_x, pc_y, pc_z, sgh0_63, sgh1_63, \
                         pdh0_314, pdh1_314, pfh0_482, pfg_315, pfh1_482, \
                         pgg_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_17 * pdh0_314[k]
                   - f_18 * pdh1_314[k]
                   + pb_x[k] * pfh0_482[k]
                   - f_9 * pc_x[k] * pfh1_482[k];

        t_693[k] = pa_z[k] * sgh0_63[k]
                   - f_9 * pc_z[k] * sgh1_63[k];

        t_694[k] = f_10 * pfg_315[k]
                   + f_4 * pc_y[k] * pgg_495[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, pa_z, pc_y, pc_z, sgh0_65, sgh0_66, sgg_45, \
                         sgh1_65, sgh1_66, pfg_317, pgg_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = pa_z[k] * sgh0_65[k]
                   + f_0 * sgg_45[k]
                   - f_9 * pc_z[k] * sgh1_65[k];

        t_696[k] = pa_z[k] * sgh0_66[k]
                   - f_9 * pc_z[k] * sgh1_66[k];

        t_697[k] = f_10 * pfg_317[k]
                   + f_4 * pc_y[k] * pgg_497[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pa_z, pc_z, sgh0_68, sgh0_69, sgh0_70, sgg_47, \
                         sgg_48, sgh1_68, sgh1_69, sgh1_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = pa_z[k] * sgh0_68[k]
                   + f_10 * sgg_47[k]
                   - f_9 * pc_z[k] * sgh1_68[k];

        t_699[k] = pa_z[k] * sgh0_69[k]
                   - f_9 * pc_z[k] * sgh1_69[k];

        t_700[k] = pa_z[k] * sgh0_70[k]
                   + f_0 * sgg_48[k]
                   - f_9 * pc_z[k] * sgh1_70[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, pa_z, pc_x, pc_y, pc_z, sgh0_72, sgg_50, \
                         sgh1_72, pfg_320, pfg_355, pgg_500, pgg_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_10 * pfg_320[k]
                   + f_4 * pc_y[k] * pgg_500[k];

        t_702[k] = pa_z[k] * sgh0_72[k]
                   + f_11 * sgg_50[k]
                   - f_9 * pc_z[k] * sgh1_72[k];

        t_703[k] = f_10 * pfg_355[k]
                   + f_4 * pc_x[k] * pgg_505[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, pc_x, pfg_356, pfg_357, pfg_358, pfg_359, \
                         pgg_506, pgg_507, pgg_508, pgg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = f_10 * pfg_356[k]
                   + f_4 * pc_x[k] * pgg_506[k];

        t_705[k] = f_10 * pfg_357[k]
                   + f_4 * pc_x[k] * pgg_507[k];

        t_706[k] = f_10 * pfg_358[k]
                   + f_4 * pc_x[k] * pgg_508[k];

        t_707[k] = f_10 * pfg_359[k]
                   + f_4 * pc_x[k] * pgg_509[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, pa_z, pc_z, sgh0_78, sgh0_79, sgh0_80, sgg_55, \
                         sgg_56, sgh1_78, sgh1_79, sgh1_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = pa_z[k] * sgh0_78[k]
                   - f_9 * pc_z[k] * sgh1_78[k];

        t_709[k] = pa_z[k] * sgh0_79[k]
                   + f_0 * sgg_55[k]
                   - f_9 * pc_z[k] * sgh1_79[k];

        t_710[k] = pa_z[k] * sgh0_80[k]
                   + f_10 * sgg_56[k]
                   - f_9 * pc_z[k] * sgh1_80[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, pa_z, pc_y, pc_z, sgh0_81, sgh0_83, sgg_57, \
                         sgg_59, sgh1_81, sgh1_83, pfg_329, pgg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = pa_z[k] * sgh0_81[k]
                   + f_11 * sgg_57[k]
                   - f_9 * pc_z[k] * sgh1_81[k];

        t_712[k] = f_10 * pfg_329[k]
                   + f_4 * pc_y[k] * pgg_509[k];

        t_713[k] = pa_z[k] * sgh0_83[k]
                   + f_14 * sgg_59[k]
                   - f_9 * pc_z[k] * sgh1_83[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pb_y, pc_y, pfh0_462, pfh0_464, pfg_330, \
                         pfh1_462, pfh1_464, pgg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pb_y[k] * pfh0_462[k]
                   - f_9 * pc_y[k] * pfh1_462[k];

        t_715[k] = f_0 * pfg_330[k]
                   + f_4 * pc_y[k] * pgg_510[k];

        t_716[k] = pb_y[k] * pfh0_464[k]
                   - f_9 * pc_y[k] * pfh1_464[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pb_y, pc_x, pc_y, pfh0_467, pfg_332, pfg_363, \
                         pfh1_467, pgf0_343, pgf1_343, pgg_512, \
                         pgg_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_10 * pfg_363[k]
                   + f_7 * pgf0_343[k]
                   - f_8 * pgf1_343[k]
                   + f_4 * pc_x[k] * pgg_513[k];

        t_718[k] = f_0 * pfg_332[k]
                   + f_4 * pc_y[k] * pgg_512[k];

        t_719[k] = pb_y[k] * pfh0_467[k]
                   - f_9 * pc_y[k] * pfh1_467[k];
    }
}

static auto
compute_prim_pgh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgh0, const size_t sgg,
                                                          const size_t sgh1, const size_t pdh0,
                                                          const size_t pdh1, const size_t pfh0,
                                                          const size_t pfg, const size_t pfh1,
                                                          const size_t pgf0, const size_t pgf1,
                                                          const size_t pgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);
    const auto f_14 = 2.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);

    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_126 = buffer.data(sgh0 + 126);
    const auto *sgh0_128 = buffer.data(sgh0 + 128);
    const auto *sgh0_129 = buffer.data(sgh0 + 129);
    const auto *sgh0_131 = buffer.data(sgh0 + 131);
    const auto *sgh0_132 = buffer.data(sgh0 + 132);
    const auto *sgh0_133 = buffer.data(sgh0 + 133);
    const auto *sgh0_135 = buffer.data(sgh0 + 135);
    const auto *sgh0_141 = buffer.data(sgh0 + 141);

    const auto *sgg_90 = buffer.data(sgg + 90);
    const auto *sgg_92 = buffer.data(sgg + 92);
    const auto *sgg_93 = buffer.data(sgg + 93);
    const auto *sgg_95 = buffer.data(sgg + 95);

    const auto *sgh1_126 = buffer.data(sgh1 + 126);
    const auto *sgh1_128 = buffer.data(sgh1 + 128);
    const auto *sgh1_129 = buffer.data(sgh1 + 129);
    const auto *sgh1_131 = buffer.data(sgh1 + 131);
    const auto *sgh1_132 = buffer.data(sgh1 + 132);
    const auto *sgh1_133 = buffer.data(sgh1 + 133);
    const auto *sgh1_135 = buffer.data(sgh1 + 135);
    const auto *sgh1_141 = buffer.data(sgh1 + 141);

    const auto *pdh0_377 = buffer.data(pdh0 + 377);

    const auto *pdh1_377 = buffer.data(pdh1 + 377);

    const auto *pfh0_469 = buffer.data(pfh0 + 469);
    const auto *pfh0_471 = buffer.data(pfh0 + 471);
    const auto *pfh0_482 = buffer.data(pfh0 + 482);
    const auto *pfh0_525 = buffer.data(pfh0 + 525);
    const auto *pfh0_527 = buffer.data(pfh0 + 527);
    const auto *pfh0_530 = buffer.data(pfh0 + 530);
    const auto *pfh0_534 = buffer.data(pfh0 + 534);
    const auto *pfh0_545 = buffer.data(pfh0 + 545);
    const auto *pfh0_562 = buffer.data(pfh0 + 562);
    const auto *pfh0_563 = buffer.data(pfh0 + 563);
    const auto *pfh0_564 = buffer.data(pfh0 + 564);
    const auto *pfh0_566 = buffer.data(pfh0 + 566);
    const auto *pfh0_569 = buffer.data(pfh0 + 569);
    const auto *pfh0_572 = buffer.data(pfh0 + 572);
    const auto *pfh0_574 = buffer.data(pfh0 + 574);
    const auto *pfh0_576 = buffer.data(pfh0 + 576);
    const auto *pfh0_582 = buffer.data(pfh0 + 582);
    const auto *pfh0_583 = buffer.data(pfh0 + 583);
    const auto *pfh0_584 = buffer.data(pfh0 + 584);
    const auto *pfh0_585 = buffer.data(pfh0 + 585);
    const auto *pfh0_587 = buffer.data(pfh0 + 587);
    const auto *pfh0_591 = buffer.data(pfh0 + 591);
    const auto *pfh0_594 = buffer.data(pfh0 + 594);
    const auto *pfh0_595 = buffer.data(pfh0 + 595);
    const auto *pfh0_603 = buffer.data(pfh0 + 603);
    const auto *pfh0_604 = buffer.data(pfh0 + 604);
    const auto *pfh0_605 = buffer.data(pfh0 + 605);
    const auto *pfh0_606 = buffer.data(pfh0 + 606);
    const auto *pfh0_608 = buffer.data(pfh0 + 608);
    const auto *pfh0_609 = buffer.data(pfh0 + 609);
    const auto *pfh0_611 = buffer.data(pfh0 + 611);
    const auto *pfh0_612 = buffer.data(pfh0 + 612);
    const auto *pfh0_614 = buffer.data(pfh0 + 614);
    const auto *pfh0_615 = buffer.data(pfh0 + 615);
    const auto *pfh0_616 = buffer.data(pfh0 + 616);
    const auto *pfh0_618 = buffer.data(pfh0 + 618);
    const auto *pfh0_624 = buffer.data(pfh0 + 624);
    const auto *pfh0_625 = buffer.data(pfh0 + 625);

    const auto *pfg_334 = buffer.data(pfg + 334);
    const auto *pfg_335 = buffer.data(pfg + 335);
    const auto *pfg_340 = buffer.data(pfg + 340);
    const auto *pfg_341 = buffer.data(pfg + 341);
    const auto *pfg_342 = buffer.data(pfg + 342);
    const auto *pfg_343 = buffer.data(pfg + 343);
    const auto *pfg_344 = buffer.data(pfg + 344);
    const auto *pfg_345 = buffer.data(pfg + 345);
    const auto *pfg_347 = buffer.data(pfg + 347);
    const auto *pfg_350 = buffer.data(pfg + 350);
    const auto *pfg_359 = buffer.data(pfg + 359);
    const auto *pfg_360 = buffer.data(pfg + 360);
    const auto *pfg_362 = buffer.data(pfg + 362);
    const auto *pfg_365 = buffer.data(pfg + 365);
    const auto *pfg_366 = buffer.data(pfg + 366);
    const auto *pfg_370 = buffer.data(pfg + 370);
    const auto *pfg_371 = buffer.data(pfg + 371);
    const auto *pfg_372 = buffer.data(pfg + 372);
    const auto *pfg_373 = buffer.data(pfg + 373);
    const auto *pfg_374 = buffer.data(pfg + 374);
    const auto *pfg_375 = buffer.data(pfg + 375);
    const auto *pfg_377 = buffer.data(pfg + 377);
    const auto *pfg_378 = buffer.data(pfg + 378);
    const auto *pfg_380 = buffer.data(pfg + 380);
    const auto *pfg_381 = buffer.data(pfg + 381);
    const auto *pfg_382 = buffer.data(pfg + 382);
    const auto *pfg_384 = buffer.data(pfg + 384);
    const auto *pfg_385 = buffer.data(pfg + 385);
    const auto *pfg_386 = buffer.data(pfg + 386);
    const auto *pfg_387 = buffer.data(pfg + 387);
    const auto *pfg_388 = buffer.data(pfg + 388);
    const auto *pfg_389 = buffer.data(pfg + 389);
    const auto *pfg_400 = buffer.data(pfg + 400);
    const auto *pfg_401 = buffer.data(pfg + 401);
    const auto *pfg_402 = buffer.data(pfg + 402);
    const auto *pfg_403 = buffer.data(pfg + 403);
    const auto *pfg_404 = buffer.data(pfg + 404);
    const auto *pfg_405 = buffer.data(pfg + 405);
    const auto *pfg_407 = buffer.data(pfg + 407);
    const auto *pfg_408 = buffer.data(pfg + 408);
    const auto *pfg_410 = buffer.data(pfg + 410);
    const auto *pfg_411 = buffer.data(pfg + 411);
    const auto *pfg_412 = buffer.data(pfg + 412);
    const auto *pfg_414 = buffer.data(pfg + 414);
    const auto *pfg_415 = buffer.data(pfg + 415);
    const auto *pfg_416 = buffer.data(pfg + 416);
    const auto *pfg_417 = buffer.data(pfg + 417);
    const auto *pfg_418 = buffer.data(pfg + 418);
    const auto *pfg_419 = buffer.data(pfg + 419);
    const auto *pfg_423 = buffer.data(pfg + 423);
    const auto *pfg_426 = buffer.data(pfg + 426);
    const auto *pfg_427 = buffer.data(pfg + 427);
    const auto *pfg_430 = buffer.data(pfg + 430);
    const auto *pfg_431 = buffer.data(pfg + 431);
    const auto *pfg_432 = buffer.data(pfg + 432);
    const auto *pfg_433 = buffer.data(pfg + 433);
    const auto *pfg_434 = buffer.data(pfg + 434);
    const auto *pfg_435 = buffer.data(pfg + 435);
    const auto *pfg_437 = buffer.data(pfg + 437);
    const auto *pfg_438 = buffer.data(pfg + 438);
    const auto *pfg_440 = buffer.data(pfg + 440);
    const auto *pfg_441 = buffer.data(pfg + 441);
    const auto *pfg_442 = buffer.data(pfg + 442);
    const auto *pfg_444 = buffer.data(pfg + 444);
    const auto *pfg_445 = buffer.data(pfg + 445);
    const auto *pfg_446 = buffer.data(pfg + 446);
    const auto *pfg_447 = buffer.data(pfg + 447);
    const auto *pfg_448 = buffer.data(pfg + 448);
    const auto *pfg_449 = buffer.data(pfg + 449);

    const auto *pfh1_469 = buffer.data(pfh1 + 469);
    const auto *pfh1_471 = buffer.data(pfh1 + 471);
    const auto *pfh1_482 = buffer.data(pfh1 + 482);
    const auto *pfh1_525 = buffer.data(pfh1 + 525);
    const auto *pfh1_527 = buffer.data(pfh1 + 527);
    const auto *pfh1_530 = buffer.data(pfh1 + 530);
    const auto *pfh1_534 = buffer.data(pfh1 + 534);
    const auto *pfh1_545 = buffer.data(pfh1 + 545);
    const auto *pfh1_562 = buffer.data(pfh1 + 562);
    const auto *pfh1_563 = buffer.data(pfh1 + 563);
    const auto *pfh1_564 = buffer.data(pfh1 + 564);
    const auto *pfh1_566 = buffer.data(pfh1 + 566);
    const auto *pfh1_569 = buffer.data(pfh1 + 569);
    const auto *pfh1_572 = buffer.data(pfh1 + 572);
    const auto *pfh1_574 = buffer.data(pfh1 + 574);
    const auto *pfh1_576 = buffer.data(pfh1 + 576);
    const auto *pfh1_582 = buffer.data(pfh1 + 582);
    const auto *pfh1_583 = buffer.data(pfh1 + 583);
    const auto *pfh1_584 = buffer.data(pfh1 + 584);
    const auto *pfh1_585 = buffer.data(pfh1 + 585);
    const auto *pfh1_587 = buffer.data(pfh1 + 587);
    const auto *pfh1_591 = buffer.data(pfh1 + 591);
    const auto *pfh1_594 = buffer.data(pfh1 + 594);
    const auto *pfh1_595 = buffer.data(pfh1 + 595);
    const auto *pfh1_603 = buffer.data(pfh1 + 603);
    const auto *pfh1_604 = buffer.data(pfh1 + 604);
    const auto *pfh1_605 = buffer.data(pfh1 + 605);
    const auto *pfh1_606 = buffer.data(pfh1 + 606);
    const auto *pfh1_608 = buffer.data(pfh1 + 608);
    const auto *pfh1_609 = buffer.data(pfh1 + 609);
    const auto *pfh1_611 = buffer.data(pfh1 + 611);
    const auto *pfh1_612 = buffer.data(pfh1 + 612);
    const auto *pfh1_614 = buffer.data(pfh1 + 614);
    const auto *pfh1_615 = buffer.data(pfh1 + 615);
    const auto *pfh1_616 = buffer.data(pfh1 + 616);
    const auto *pfh1_618 = buffer.data(pfh1 + 618);
    const auto *pfh1_624 = buffer.data(pfh1 + 624);
    const auto *pfh1_625 = buffer.data(pfh1 + 625);

    const auto *pgf0_346 = buffer.data(pgf0 + 346);
    const auto *pgf0_347 = buffer.data(pgf0 + 347);
    const auto *pgf0_348 = buffer.data(pgf0 + 348);
    const auto *pgf0_349 = buffer.data(pgf0 + 349);
    const auto *pgf0_350 = buffer.data(pgf0 + 350);
    const auto *pgf0_352 = buffer.data(pgf0 + 352);
    const auto *pgf0_353 = buffer.data(pgf0 + 353);
    const auto *pgf0_355 = buffer.data(pgf0 + 355);
    const auto *pgf0_356 = buffer.data(pgf0 + 356);
    const auto *pgf0_357 = buffer.data(pgf0 + 357);
    const auto *pgf0_358 = buffer.data(pgf0 + 358);
    const auto *pgf0_359 = buffer.data(pgf0 + 359);
    const auto *pgf0_370 = buffer.data(pgf0 + 370);
    const auto *pgf0_373 = buffer.data(pgf0 + 373);
    const auto *pgf0_376 = buffer.data(pgf0 + 376);

    const auto *pgf1_346 = buffer.data(pgf1 + 346);
    const auto *pgf1_347 = buffer.data(pgf1 + 347);
    const auto *pgf1_348 = buffer.data(pgf1 + 348);
    const auto *pgf1_349 = buffer.data(pgf1 + 349);
    const auto *pgf1_350 = buffer.data(pgf1 + 350);
    const auto *pgf1_352 = buffer.data(pgf1 + 352);
    const auto *pgf1_353 = buffer.data(pgf1 + 353);
    const auto *pgf1_355 = buffer.data(pgf1 + 355);
    const auto *pgf1_356 = buffer.data(pgf1 + 356);
    const auto *pgf1_357 = buffer.data(pgf1 + 357);
    const auto *pgf1_358 = buffer.data(pgf1 + 358);
    const auto *pgf1_359 = buffer.data(pgf1 + 359);
    const auto *pgf1_370 = buffer.data(pgf1 + 370);
    const auto *pgf1_373 = buffer.data(pgf1 + 373);
    const auto *pgf1_376 = buffer.data(pgf1 + 376);

    const auto *pgg_515 = buffer.data(pgg + 515);
    const auto *pgg_516 = buffer.data(pgg + 516);
    const auto *pgg_520 = buffer.data(pgg + 520);
    const auto *pgg_521 = buffer.data(pgg + 521);
    const auto *pgg_522 = buffer.data(pgg + 522);
    const auto *pgg_523 = buffer.data(pgg + 523);
    const auto *pgg_524 = buffer.data(pgg + 524);
    const auto *pgg_525 = buffer.data(pgg + 525);
    const auto *pgg_527 = buffer.data(pgg + 527);
    const auto *pgg_528 = buffer.data(pgg + 528);
    const auto *pgg_530 = buffer.data(pgg + 530);
    const auto *pgg_531 = buffer.data(pgg + 531);
    const auto *pgg_532 = buffer.data(pgg + 532);
    const auto *pgg_534 = buffer.data(pgg + 534);
    const auto *pgg_535 = buffer.data(pgg + 535);
    const auto *pgg_536 = buffer.data(pgg + 536);
    const auto *pgg_537 = buffer.data(pgg + 537);
    const auto *pgg_538 = buffer.data(pgg + 538);
    const auto *pgg_539 = buffer.data(pgg + 539);
    const auto *pgg_540 = buffer.data(pgg + 540);
    const auto *pgg_542 = buffer.data(pgg + 542);
    const auto *pgg_545 = buffer.data(pgg + 545);
    const auto *pgg_550 = buffer.data(pgg + 550);
    const auto *pgg_551 = buffer.data(pgg + 551);
    const auto *pgg_552 = buffer.data(pgg + 552);
    const auto *pgg_553 = buffer.data(pgg + 553);
    const auto *pgg_554 = buffer.data(pgg + 554);
    const auto *pgg_555 = buffer.data(pgg + 555);
    const auto *pgg_557 = buffer.data(pgg + 557);
    const auto *pgg_558 = buffer.data(pgg + 558);
    const auto *pgg_560 = buffer.data(pgg + 560);
    const auto *pgg_561 = buffer.data(pgg + 561);
    const auto *pgg_565 = buffer.data(pgg + 565);
    const auto *pgg_566 = buffer.data(pgg + 566);
    const auto *pgg_567 = buffer.data(pgg + 567);
    const auto *pgg_568 = buffer.data(pgg + 568);
    const auto *pgg_569 = buffer.data(pgg + 569);
    const auto *pgg_570 = buffer.data(pgg + 570);
    const auto *pgg_572 = buffer.data(pgg + 572);
    const auto *pgg_575 = buffer.data(pgg + 575);
    const auto *pgg_580 = buffer.data(pgg + 580);
    const auto *pgg_581 = buffer.data(pgg + 581);
    const auto *pgg_582 = buffer.data(pgg + 582);
    const auto *pgg_583 = buffer.data(pgg + 583);
    const auto *pgg_584 = buffer.data(pgg + 584);
    const auto *pgg_585 = buffer.data(pgg + 585);
    const auto *pgg_587 = buffer.data(pgg + 587);
    const auto *pgg_590 = buffer.data(pgg + 590);
    const auto *pgg_595 = buffer.data(pgg + 595);
    const auto *pgg_596 = buffer.data(pgg + 596);
    const auto *pgg_597 = buffer.data(pgg + 597);
    const auto *pgg_598 = buffer.data(pgg + 598);
    const auto *pgg_599 = buffer.data(pgg + 599);

#pragma omp simd aligned(t_720, t_721, t_722, pb_y, pc_x, pc_y, pfh0_469, pfg_334, pfg_335, \
                         pfg_366, pfh1_469, pgf0_346, pgf1_346, pgg_515, \
                         pgg_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_10 * pfg_366[k]
                   + f_5 * pgf0_346[k]
                   - f_6 * pgf1_346[k]
                   + f_4 * pc_x[k] * pgg_516[k];

        t_721[k] = pb_y[k] * pfh0_469[k]
                   + f_10 * pfg_334[k]
                   - f_9 * pc_y[k] * pfh1_469[k];

        t_722[k] = f_0 * pfg_335[k]
                   + f_4 * pc_y[k] * pgg_515[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, t_726, pb_y, pc_x, pc_y, pfh0_471, pfg_370, \
                         pfg_371, pfg_372, pfh1_471, pgg_520, pgg_521, \
                         pgg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = pb_y[k] * pfh0_471[k]
                   - f_9 * pc_y[k] * pfh1_471[k];

        t_724[k] = f_10 * pfg_370[k]
                   + f_4 * pc_x[k] * pgg_520[k];

        t_725[k] = f_10 * pfg_371[k]
                   + f_4 * pc_x[k] * pgg_521[k];

        t_726[k] = f_10 * pfg_372[k]
                   + f_4 * pc_x[k] * pgg_522[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pc_x, pc_y, pfg_340, pfg_373, pfg_374, pgf0_346, \
                         pgf1_346, pgg_520, pgg_523, pgg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_10 * pfg_373[k]
                   + f_4 * pc_x[k] * pgg_523[k];

        t_728[k] = f_10 * pfg_374[k]
                   + f_4 * pc_x[k] * pgg_524[k];

        t_729[k] = f_0 * pfg_340[k]
                   + f_2 * pgf0_346[k]
                   - f_3 * pgf1_346[k]
                   + f_4 * pc_y[k] * pgg_520[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pc_y, pfg_341, pfg_342, pfg_343, pgf0_347, \
                         pgf0_348, pgf0_349, pgf1_347, pgf1_348, pgf1_349, pgg_521, pgg_522, \
                         pgg_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_0 * pfg_341[k]
                   + f_15 * pgf0_347[k]
                   - f_16 * pgf1_347[k]
                   + f_4 * pc_y[k] * pgg_521[k];

        t_731[k] = f_0 * pfg_342[k]
                   + f_7 * pgf0_348[k]
                   - f_8 * pgf1_348[k]
                   + f_4 * pc_y[k] * pgg_522[k];

        t_732[k] = f_0 * pfg_343[k]
                   + f_5 * pgf0_349[k]
                   - f_6 * pgf1_349[k]
                   + f_4 * pc_y[k] * pgg_523[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, t_736, pb_y, pc_x, pc_y, pfh0_482, pfg_344, \
                         pfg_375, pfh1_482, pgf0_350, pgf1_350, pgg_524, \
                         pgg_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_0 * pfg_344[k]
                   + f_4 * pc_y[k] * pgg_524[k];

        t_734[k] = pb_y[k] * pfh0_482[k]
                   - f_9 * pc_y[k] * pfh1_482[k];

        t_735[k] = f_10 * pfg_375[k]
                   + f_2 * pgf0_350[k]
                   - f_3 * pgf1_350[k]
                   + f_4 * pc_x[k] * pgg_525[k];

        t_736[k] = f_4 * pc_y[k] * pgg_525[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pc_x, pc_y, pfg_377, pfg_378, pgf0_352, \
                         pgf0_353, pgf1_352, pgf1_353, pgg_527, \
                         pgg_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_10 * pfg_377[k]
                   + f_15 * pgf0_352[k]
                   - f_16 * pgf1_352[k]
                   + f_4 * pc_x[k] * pgg_527[k];

        t_738[k] = f_10 * pfg_378[k]
                   + f_7 * pgf0_353[k]
                   - f_8 * pgf1_353[k]
                   + f_4 * pc_x[k] * pgg_528[k];

        t_739[k] = f_4 * pc_y[k] * pgg_527[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, pc_x, pfg_380, pfg_381, pfg_382, pgf0_355, \
                         pgf0_356, pgf0_357, pgf1_355, pgf1_356, pgf1_357, pgg_530, pgg_531, \
                         pgg_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_10 * pfg_380[k]
                   + f_7 * pgf0_355[k]
                   - f_8 * pgf1_355[k]
                   + f_4 * pc_x[k] * pgg_530[k];

        t_741[k] = f_10 * pfg_381[k]
                   + f_5 * pgf0_356[k]
                   - f_6 * pgf1_356[k]
                   + f_4 * pc_x[k] * pgg_531[k];

        t_742[k] = f_10 * pfg_382[k]
                   + f_5 * pgf0_357[k]
                   - f_6 * pgf1_357[k]
                   + f_4 * pc_x[k] * pgg_532[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, pc_x, pc_y, pfg_384, pfg_385, pfg_386, \
                         pgf0_359, pgf1_359, pgg_530, pgg_534, pgg_535, \
                         pgg_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_4 * pc_y[k] * pgg_530[k];

        t_744[k] = f_10 * pfg_384[k]
                   + f_5 * pgf0_359[k]
                   - f_6 * pgf1_359[k]
                   + f_4 * pc_x[k] * pgg_534[k];

        t_745[k] = f_10 * pfg_385[k]
                   + f_4 * pc_x[k] * pgg_535[k];

        t_746[k] = f_10 * pfg_386[k]
                   + f_4 * pc_x[k] * pgg_536[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, t_750, pc_x, pc_y, pfg_387, pfg_388, pfg_389, \
                         pgf0_356, pgf1_356, pgg_535, pgg_537, pgg_538, \
                         pgg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_10 * pfg_387[k]
                   + f_4 * pc_x[k] * pgg_537[k];

        t_748[k] = f_10 * pfg_388[k]
                   + f_4 * pc_x[k] * pgg_538[k];

        t_749[k] = f_10 * pfg_389[k]
                   + f_4 * pc_x[k] * pgg_539[k];

        t_750[k] = f_2 * pgf0_356[k]
                   - f_3 * pgf1_356[k]
                   + f_4 * pc_y[k] * pgg_535[k];
    }

#pragma omp simd aligned(t_751, t_752, t_753, t_754, pc_y, pgf0_357, pgf0_358, pgf0_359, \
                         pgf1_357, pgf1_358, pgf1_359, pgg_536, pgg_537, pgg_538, \
                         pgg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_751[k] = f_15 * pgf0_357[k]
                   - f_16 * pgf1_357[k]
                   + f_4 * pc_y[k] * pgg_536[k];

        t_752[k] = f_7 * pgf0_358[k]
                   - f_8 * pgf1_358[k]
                   + f_4 * pc_y[k] * pgg_537[k];

        t_753[k] = f_5 * pgf0_359[k]
                   - f_6 * pgf1_359[k]
                   + f_4 * pc_y[k] * pgg_538[k];

        t_754[k] = f_4 * pc_y[k] * pgg_539[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, pa_z, pb_x, pc_x, pc_y, pc_z, sgh0_126, \
                         sgh1_126, pdh0_377, pdh1_377, pfh0_545, pfg_345, pfh1_545, \
                         pgg_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_12 * pdh0_377[k]
                   - f_13 * pdh1_377[k]
                   + pb_x[k] * pfh0_545[k]
                   - f_9 * pc_x[k] * pfh1_545[k];

        t_756[k] = pa_z[k] * sgh0_126[k]
                   - f_9 * pc_z[k] * sgh1_126[k];

        t_757[k] = f_11 * pfg_345[k]
                   + f_4 * pc_y[k] * pgg_540[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pa_z, pc_y, pc_z, sgh0_128, sgh0_129, sgg_90, \
                         sgh1_128, sgh1_129, pfg_347, pgg_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = pa_z[k] * sgh0_128[k]
                   + f_0 * sgg_90[k]
                   - f_9 * pc_z[k] * sgh1_128[k];

        t_759[k] = pa_z[k] * sgh0_129[k]
                   - f_9 * pc_z[k] * sgh1_129[k];

        t_760[k] = f_11 * pfg_347[k]
                   + f_4 * pc_y[k] * pgg_542[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pa_z, pc_z, sgh0_131, sgh0_132, sgh0_133, \
                         sgg_92, sgg_93, sgh1_131, sgh1_132, sgh1_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = pa_z[k] * sgh0_131[k]
                   + f_10 * sgg_92[k]
                   - f_9 * pc_z[k] * sgh1_131[k];

        t_762[k] = pa_z[k] * sgh0_132[k]
                   - f_9 * pc_z[k] * sgh1_132[k];

        t_763[k] = pa_z[k] * sgh0_133[k]
                   + f_0 * sgg_93[k]
                   - f_9 * pc_z[k] * sgh1_133[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pa_z, pc_x, pc_y, pc_z, sgh0_135, sgg_95, \
                         sgh1_135, pfg_350, pfg_400, pgg_545, pgg_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_11 * pfg_350[k]
                   + f_4 * pc_y[k] * pgg_545[k];

        t_765[k] = pa_z[k] * sgh0_135[k]
                   + f_11 * sgg_95[k]
                   - f_9 * pc_z[k] * sgh1_135[k];

        t_766[k] = f_0 * pfg_400[k]
                   + f_4 * pc_x[k] * pgg_550[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, t_770, pc_x, pfg_401, pfg_402, pfg_403, pfg_404, \
                         pgg_551, pgg_552, pgg_553, pgg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_0 * pfg_401[k]
                   + f_4 * pc_x[k] * pgg_551[k];

        t_768[k] = f_0 * pfg_402[k]
                   + f_4 * pc_x[k] * pgg_552[k];

        t_769[k] = f_0 * pfg_403[k]
                   + f_4 * pc_x[k] * pgg_553[k];

        t_770[k] = f_0 * pfg_404[k]
                   + f_4 * pc_x[k] * pgg_554[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, pa_z, pb_x, pc_x, pc_z, sgh0_141, \
                         sgh1_141, pfh0_562, pfh0_563, pfh0_564, pfh1_562, pfh1_563, \
                         pfh1_564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = pa_z[k] * sgh0_141[k]
                   - f_9 * pc_z[k] * sgh1_141[k];

        t_772[k] = pb_x[k] * pfh0_562[k]
                   - f_9 * pc_x[k] * pfh1_562[k];

        t_773[k] = pb_x[k] * pfh0_563[k]
                   - f_9 * pc_x[k] * pfh1_563[k];

        t_774[k] = pb_x[k] * pfh0_564[k]
                   - f_9 * pc_x[k] * pfh1_564[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, pb_x, pc_x, pc_y, pfh0_566, pfg_359, \
                         pfg_360, pfg_405, pfh1_566, pgf0_370, pgf1_370, pgg_554, \
                         pgg_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_11 * pfg_359[k]
                   + f_4 * pc_y[k] * pgg_554[k];

        t_776[k] = pb_x[k] * pfh0_566[k]
                   - f_9 * pc_x[k] * pfh1_566[k];

        t_777[k] = f_0 * pfg_405[k]
                   + f_2 * pgf0_370[k]
                   - f_3 * pgf1_370[k]
                   + f_4 * pc_x[k] * pgg_555[k];

        t_778[k] = f_10 * pfg_360[k]
                   + f_4 * pc_y[k] * pgg_555[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pb_x, pc_x, pc_y, pfh0_569, pfg_362, pfg_407, \
                         pfg_408, pfh1_569, pgf0_373, pgf1_373, pgg_557, \
                         pgg_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = pb_x[k] * pfh0_569[k]
                   + f_1 * pfg_407[k]
                   - f_9 * pc_x[k] * pfh1_569[k];

        t_780[k] = f_0 * pfg_408[k]
                   + f_7 * pgf0_373[k]
                   - f_8 * pgf1_373[k]
                   + f_4 * pc_x[k] * pgg_558[k];

        t_781[k] = f_10 * pfg_362[k]
                   + f_4 * pc_y[k] * pgg_557[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pb_x, pc_x, pfh0_572, pfh0_574, pfg_410, \
                         pfg_411, pfg_412, pfh1_572, pfh1_574, pgf0_376, pgf1_376, \
                         pgg_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = pb_x[k] * pfh0_572[k]
                   + f_11 * pfg_410[k]
                   - f_9 * pc_x[k] * pfh1_572[k];

        t_783[k] = f_0 * pfg_411[k]
                   + f_5 * pgf0_376[k]
                   - f_6 * pgf1_376[k]
                   + f_4 * pc_x[k] * pgg_561[k];

        t_784[k] = pb_x[k] * pfh0_574[k]
                   + f_10 * pfg_412[k]
                   - f_9 * pc_x[k] * pfh1_574[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, pb_x, pc_x, pc_y, pfh0_576, pfg_365, \
                         pfg_414, pfg_415, pfg_416, pfh1_576, pgg_560, pgg_565, \
                         pgg_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_10 * pfg_365[k]
                   + f_4 * pc_y[k] * pgg_560[k];

        t_786[k] = pb_x[k] * pfh0_576[k]
                   + f_10 * pfg_414[k]
                   - f_9 * pc_x[k] * pfh1_576[k];

        t_787[k] = f_0 * pfg_415[k]
                   + f_4 * pc_x[k] * pgg_565[k];

        t_788[k] = f_0 * pfg_416[k]
                   + f_4 * pc_x[k] * pgg_566[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, t_792, pb_x, pc_x, pfh0_582, pfg_417, pfg_418, \
                         pfg_419, pfh1_582, pgg_567, pgg_568, pgg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_0 * pfg_417[k]
                   + f_4 * pc_x[k] * pgg_567[k];

        t_790[k] = f_0 * pfg_418[k]
                   + f_4 * pc_x[k] * pgg_568[k];

        t_791[k] = f_0 * pfg_419[k]
                   + f_4 * pc_x[k] * pgg_569[k];

        t_792[k] = pb_x[k] * pfh0_582[k]
                   - f_9 * pc_x[k] * pfh1_582[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, t_796, pb_x, pc_x, pc_y, pfh0_583, pfh0_584, \
                         pfh0_585, pfg_374, pfh1_583, pfh1_584, pfh1_585, \
                         pgg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = pb_x[k] * pfh0_583[k]
                   - f_9 * pc_x[k] * pfh1_583[k];

        t_794[k] = pb_x[k] * pfh0_584[k]
                   - f_9 * pc_x[k] * pfh1_584[k];

        t_795[k] = pb_x[k] * pfh0_585[k]
                   - f_9 * pc_x[k] * pfh1_585[k];

        t_796[k] = f_10 * pfg_374[k]
                   + f_4 * pc_y[k] * pgg_569[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, t_800, pb_x, pb_y, pc_x, pc_y, pfh0_525, \
                         pfh0_527, pfh0_587, pfg_375, pfh1_525, pfh1_527, pfh1_587, \
                         pgg_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = pb_x[k] * pfh0_587[k]
                   - f_9 * pc_x[k] * pfh1_587[k];

        t_798[k] = pb_y[k] * pfh0_525[k]
                   - f_9 * pc_y[k] * pfh1_525[k];

        t_799[k] = f_0 * pfg_375[k]
                   + f_4 * pc_y[k] * pgg_570[k];

        t_800[k] = pb_y[k] * pfh0_527[k]
                   - f_9 * pc_y[k] * pfh1_527[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, pb_x, pb_y, pc_x, pc_y, pfh0_530, pfh0_591, \
                         pfg_377, pfg_423, pfh1_530, pfh1_591, \
                         pgg_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = pb_x[k] * pfh0_591[k]
                   + f_11 * pfg_423[k]
                   - f_9 * pc_x[k] * pfh1_591[k];

        t_802[k] = f_0 * pfg_377[k]
                   + f_4 * pc_y[k] * pgg_572[k];

        t_803[k] = pb_y[k] * pfh0_530[k]
                   - f_9 * pc_y[k] * pfh1_530[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pb_x, pc_x, pc_y, pfh0_594, pfh0_595, pfg_380, \
                         pfg_426, pfg_427, pfh1_594, pfh1_595, \
                         pgg_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = pb_x[k] * pfh0_594[k]
                   + f_10 * pfg_426[k]
                   - f_9 * pc_x[k] * pfh1_594[k];

        t_805[k] = pb_x[k] * pfh0_595[k]
                   + f_10 * pfg_427[k]
                   - f_9 * pc_x[k] * pfh1_595[k];

        t_806[k] = f_0 * pfg_380[k]
                   + f_4 * pc_y[k] * pgg_575[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, pb_y, pc_x, pc_y, pfh0_534, pfg_430, \
                         pfg_431, pfg_432, pfh1_534, pgg_580, pgg_581, \
                         pgg_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = pb_y[k] * pfh0_534[k]
                   - f_9 * pc_y[k] * pfh1_534[k];

        t_808[k] = f_0 * pfg_430[k]
                   + f_4 * pc_x[k] * pgg_580[k];

        t_809[k] = f_0 * pfg_431[k]
                   + f_4 * pc_x[k] * pgg_581[k];

        t_810[k] = f_0 * pfg_432[k]
                   + f_4 * pc_x[k] * pgg_582[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, pb_x, pc_x, pfh0_603, pfh0_604, pfg_433, \
                         pfg_434, pfh1_603, pfh1_604, pgg_583, \
                         pgg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_0 * pfg_433[k]
                   + f_4 * pc_x[k] * pgg_583[k];

        t_812[k] = f_0 * pfg_434[k]
                   + f_4 * pc_x[k] * pgg_584[k];

        t_813[k] = pb_x[k] * pfh0_603[k]
                   - f_9 * pc_x[k] * pfh1_603[k];

        t_814[k] = pb_x[k] * pfh0_604[k]
                   - f_9 * pc_x[k] * pfh1_604[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, pb_x, pc_x, pc_y, pfh0_605, pfh0_606, \
                         pfh0_608, pfg_389, pfh1_605, pfh1_606, pfh1_608, \
                         pgg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = pb_x[k] * pfh0_605[k]
                   - f_9 * pc_x[k] * pfh1_605[k];

        t_816[k] = pb_x[k] * pfh0_606[k]
                   - f_9 * pc_x[k] * pfh1_606[k];

        t_817[k] = f_0 * pfg_389[k]
                   + f_4 * pc_y[k] * pgg_584[k];

        t_818[k] = pb_x[k] * pfh0_608[k]
                   - f_9 * pc_x[k] * pfh1_608[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pb_x, pc_x, pc_y, pfh0_609, pfh0_611, pfg_435, \
                         pfg_437, pfh1_609, pfh1_611, pgg_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = pb_x[k] * pfh0_609[k]
                   + f_14 * pfg_435[k]
                   - f_9 * pc_x[k] * pfh1_609[k];

        t_820[k] = f_4 * pc_y[k] * pgg_585[k];

        t_821[k] = pb_x[k] * pfh0_611[k]
                   + f_1 * pfg_437[k]
                   - f_9 * pc_x[k] * pfh1_611[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pb_x, pc_x, pc_y, pfh0_612, pfh0_614, pfg_438, \
                         pfg_440, pfh1_612, pfh1_614, pgg_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pb_x[k] * pfh0_612[k]
                   + f_11 * pfg_438[k]
                   - f_9 * pc_x[k] * pfh1_612[k];

        t_823[k] = f_4 * pc_y[k] * pgg_587[k];

        t_824[k] = pb_x[k] * pfh0_614[k]
                   + f_11 * pfg_440[k]
                   - f_9 * pc_x[k] * pfh1_614[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, pb_x, pc_x, pc_y, pfh0_615, pfh0_616, pfg_441, \
                         pfg_442, pfh1_615, pfh1_616, pgg_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = pb_x[k] * pfh0_615[k]
                   + f_10 * pfg_441[k]
                   - f_9 * pc_x[k] * pfh1_615[k];

        t_826[k] = pb_x[k] * pfh0_616[k]
                   + f_10 * pfg_442[k]
                   - f_9 * pc_x[k] * pfh1_616[k];

        t_827[k] = f_4 * pc_y[k] * pgg_590[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, t_831, pb_x, pc_x, pfh0_618, pfg_444, pfg_445, \
                         pfg_446, pfg_447, pfh1_618, pgg_595, pgg_596, \
                         pgg_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = pb_x[k] * pfh0_618[k]
                   + f_10 * pfg_444[k]
                   - f_9 * pc_x[k] * pfh1_618[k];

        t_829[k] = f_0 * pfg_445[k]
                   + f_4 * pc_x[k] * pgg_595[k];

        t_830[k] = f_0 * pfg_446[k]
                   + f_4 * pc_x[k] * pgg_596[k];

        t_831[k] = f_0 * pfg_447[k]
                   + f_4 * pc_x[k] * pgg_597[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, pb_x, pc_x, pfh0_624, pfh0_625, pfg_448, \
                         pfg_449, pfh1_624, pfh1_625, pgg_598, \
                         pgg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_0 * pfg_448[k]
                   + f_4 * pc_x[k] * pgg_598[k];

        t_833[k] = f_0 * pfg_449[k]
                   + f_4 * pc_x[k] * pgg_599[k];

        t_834[k] = pb_x[k] * pfh0_624[k]
                   - f_9 * pc_x[k] * pfh1_624[k];

        t_835[k] = pb_x[k] * pfh0_625[k]
                   - f_9 * pc_x[k] * pfh1_625[k];
    }
}

static auto
compute_prim_pgh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgh0, const size_t sgg,
                                                          const size_t sgh1, const size_t pdh0,
                                                          const size_t pdh1, const size_t pfh0,
                                                          const size_t pfg, const size_t pfh1,
                                                          const size_t pgf0, const size_t pgf1,
                                                          const size_t pgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);
    const auto f_14 = 2.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);
    const auto f_17 = 1.0 / p;
    const auto f_18 = gamma / (p * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_210 = buffer.data(sgh0 + 210);
    const auto *sgh0_213 = buffer.data(sgh0 + 213);
    const auto *sgh0_216 = buffer.data(sgh0 + 216);
    const auto *sgh0_225 = buffer.data(sgh0 + 225);
    const auto *sgh0_226 = buffer.data(sgh0 + 226);
    const auto *sgh0_227 = buffer.data(sgh0 + 227);
    const auto *sgh0_228 = buffer.data(sgh0 + 228);
    const auto *sgh0_230 = buffer.data(sgh0 + 230);

    const auto *sgg_160 = buffer.data(sgg + 160);
    const auto *sgg_161 = buffer.data(sgg + 161);
    const auto *sgg_162 = buffer.data(sgg + 162);
    const auto *sgg_164 = buffer.data(sgg + 164);
    const auto *sgg_224 = buffer.data(sgg + 224);

    const auto *sgh1_210 = buffer.data(sgh1 + 210);
    const auto *sgh1_213 = buffer.data(sgh1 + 213);
    const auto *sgh1_216 = buffer.data(sgh1 + 216);
    const auto *sgh1_225 = buffer.data(sgh1 + 225);
    const auto *sgh1_226 = buffer.data(sgh1 + 226);
    const auto *sgh1_227 = buffer.data(sgh1 + 227);
    const auto *sgh1_228 = buffer.data(sgh1 + 228);
    const auto *sgh1_230 = buffer.data(sgh1 + 230);

    const auto *pdh0_356 = buffer.data(pdh0 + 356);
    const auto *pdh0_377 = buffer.data(pdh0 + 377);

    const auto *pdh1_356 = buffer.data(pdh1 + 356);
    const auto *pdh1_377 = buffer.data(pdh1 + 377);

    const auto *pfh0_587 = buffer.data(pfh0 + 587);
    const auto *pfh0_608 = buffer.data(pfh0 + 608);
    const auto *pfh0_609 = buffer.data(pfh0 + 609);
    const auto *pfh0_611 = buffer.data(pfh0 + 611);
    const auto *pfh0_614 = buffer.data(pfh0 + 614);
    const auto *pfh0_618 = buffer.data(pfh0 + 618);
    const auto *pfh0_624 = buffer.data(pfh0 + 624);
    const auto *pfh0_625 = buffer.data(pfh0 + 625);
    const auto *pfh0_626 = buffer.data(pfh0 + 626);
    const auto *pfh0_627 = buffer.data(pfh0 + 627);
    const auto *pfh0_629 = buffer.data(pfh0 + 629);

    const auto *pfg_390 = buffer.data(pfg + 390);
    const auto *pfg_392 = buffer.data(pfg + 392);
    const auto *pfg_395 = buffer.data(pfg + 395);
    const auto *pfg_404 = buffer.data(pfg + 404);
    const auto *pfg_405 = buffer.data(pfg + 405);
    const auto *pfg_407 = buffer.data(pfg + 407);
    const auto *pfg_410 = buffer.data(pfg + 410);
    const auto *pfg_415 = buffer.data(pfg + 415);
    const auto *pfg_416 = buffer.data(pfg + 416);
    const auto *pfg_417 = buffer.data(pfg + 417);
    const auto *pfg_418 = buffer.data(pfg + 418);
    const auto *pfg_419 = buffer.data(pfg + 419);
    const auto *pfg_420 = buffer.data(pfg + 420);
    const auto *pfg_422 = buffer.data(pfg + 422);
    const auto *pfg_425 = buffer.data(pfg + 425);
    const auto *pfg_430 = buffer.data(pfg + 430);
    const auto *pfg_431 = buffer.data(pfg + 431);
    const auto *pfg_432 = buffer.data(pfg + 432);
    const auto *pfg_433 = buffer.data(pfg + 433);
    const auto *pfg_434 = buffer.data(pfg + 434);
    const auto *pfg_435 = buffer.data(pfg + 435);
    const auto *pfg_437 = buffer.data(pfg + 437);
    const auto *pfg_440 = buffer.data(pfg + 440);
    const auto *pfg_445 = buffer.data(pfg + 445);
    const auto *pfg_446 = buffer.data(pfg + 446);
    const auto *pfg_447 = buffer.data(pfg + 447);
    const auto *pfg_448 = buffer.data(pfg + 448);
    const auto *pfg_449 = buffer.data(pfg + 449);

    const auto *pfh1_587 = buffer.data(pfh1 + 587);
    const auto *pfh1_608 = buffer.data(pfh1 + 608);
    const auto *pfh1_609 = buffer.data(pfh1 + 609);
    const auto *pfh1_611 = buffer.data(pfh1 + 611);
    const auto *pfh1_614 = buffer.data(pfh1 + 614);
    const auto *pfh1_618 = buffer.data(pfh1 + 618);
    const auto *pfh1_624 = buffer.data(pfh1 + 624);
    const auto *pfh1_625 = buffer.data(pfh1 + 625);
    const auto *pfh1_626 = buffer.data(pfh1 + 626);
    const auto *pfh1_627 = buffer.data(pfh1 + 627);
    const auto *pfh1_629 = buffer.data(pfh1 + 629);

    const auto *pgf0_402 = buffer.data(pgf0 + 402);
    const auto *pgf0_405 = buffer.data(pgf0 + 405);
    const auto *pgf0_407 = buffer.data(pgf0 + 407);
    const auto *pgf0_409 = buffer.data(pgf0 + 409);
    const auto *pgf0_410 = buffer.data(pgf0 + 410);
    const auto *pgf0_412 = buffer.data(pgf0 + 412);
    const auto *pgf0_413 = buffer.data(pgf0 + 413);
    const auto *pgf0_415 = buffer.data(pgf0 + 415);
    const auto *pgf0_416 = buffer.data(pgf0 + 416);
    const auto *pgf0_417 = buffer.data(pgf0 + 417);
    const auto *pgf0_418 = buffer.data(pgf0 + 418);
    const auto *pgf0_419 = buffer.data(pgf0 + 419);
    const auto *pgf0_420 = buffer.data(pgf0 + 420);
    const auto *pgf0_422 = buffer.data(pgf0 + 422);
    const auto *pgf0_423 = buffer.data(pgf0 + 423);
    const auto *pgf0_425 = buffer.data(pgf0 + 425);
    const auto *pgf0_426 = buffer.data(pgf0 + 426);
    const auto *pgf0_427 = buffer.data(pgf0 + 427);
    const auto *pgf0_428 = buffer.data(pgf0 + 428);
    const auto *pgf0_429 = buffer.data(pgf0 + 429);
    const auto *pgf0_433 = buffer.data(pgf0 + 433);
    const auto *pgf0_436 = buffer.data(pgf0 + 436);
    const auto *pgf0_437 = buffer.data(pgf0 + 437);
    const auto *pgf0_440 = buffer.data(pgf0 + 440);
    const auto *pgf0_442 = buffer.data(pgf0 + 442);
    const auto *pgf0_443 = buffer.data(pgf0 + 443);
    const auto *pgf0_445 = buffer.data(pgf0 + 445);
    const auto *pgf0_446 = buffer.data(pgf0 + 446);
    const auto *pgf0_447 = buffer.data(pgf0 + 447);
    const auto *pgf0_448 = buffer.data(pgf0 + 448);
    const auto *pgf0_449 = buffer.data(pgf0 + 449);

    const auto *pgf1_402 = buffer.data(pgf1 + 402);
    const auto *pgf1_405 = buffer.data(pgf1 + 405);
    const auto *pgf1_407 = buffer.data(pgf1 + 407);
    const auto *pgf1_409 = buffer.data(pgf1 + 409);
    const auto *pgf1_410 = buffer.data(pgf1 + 410);
    const auto *pgf1_412 = buffer.data(pgf1 + 412);
    const auto *pgf1_413 = buffer.data(pgf1 + 413);
    const auto *pgf1_415 = buffer.data(pgf1 + 415);
    const auto *pgf1_416 = buffer.data(pgf1 + 416);
    const auto *pgf1_417 = buffer.data(pgf1 + 417);
    const auto *pgf1_418 = buffer.data(pgf1 + 418);
    const auto *pgf1_419 = buffer.data(pgf1 + 419);
    const auto *pgf1_420 = buffer.data(pgf1 + 420);
    const auto *pgf1_422 = buffer.data(pgf1 + 422);
    const auto *pgf1_423 = buffer.data(pgf1 + 423);
    const auto *pgf1_425 = buffer.data(pgf1 + 425);
    const auto *pgf1_426 = buffer.data(pgf1 + 426);
    const auto *pgf1_427 = buffer.data(pgf1 + 427);
    const auto *pgf1_428 = buffer.data(pgf1 + 428);
    const auto *pgf1_429 = buffer.data(pgf1 + 429);
    const auto *pgf1_433 = buffer.data(pgf1 + 433);
    const auto *pgf1_436 = buffer.data(pgf1 + 436);
    const auto *pgf1_437 = buffer.data(pgf1 + 437);
    const auto *pgf1_440 = buffer.data(pgf1 + 440);
    const auto *pgf1_442 = buffer.data(pgf1 + 442);
    const auto *pgf1_443 = buffer.data(pgf1 + 443);
    const auto *pgf1_445 = buffer.data(pgf1 + 445);
    const auto *pgf1_446 = buffer.data(pgf1 + 446);
    const auto *pgf1_447 = buffer.data(pgf1 + 447);
    const auto *pgf1_448 = buffer.data(pgf1 + 448);
    const auto *pgf1_449 = buffer.data(pgf1 + 449);

    const auto *pgg_599 = buffer.data(pgg + 599);
    const auto *pgg_600 = buffer.data(pgg + 600);
    const auto *pgg_602 = buffer.data(pgg + 602);
    const auto *pgg_605 = buffer.data(pgg + 605);
    const auto *pgg_607 = buffer.data(pgg + 607);
    const auto *pgg_609 = buffer.data(pgg + 609);
    const auto *pgg_610 = buffer.data(pgg + 610);
    const auto *pgg_611 = buffer.data(pgg + 611);
    const auto *pgg_612 = buffer.data(pgg + 612);
    const auto *pgg_613 = buffer.data(pgg + 613);
    const auto *pgg_614 = buffer.data(pgg + 614);
    const auto *pgg_615 = buffer.data(pgg + 615);
    const auto *pgg_617 = buffer.data(pgg + 617);
    const auto *pgg_618 = buffer.data(pgg + 618);
    const auto *pgg_620 = buffer.data(pgg + 620);
    const auto *pgg_621 = buffer.data(pgg + 621);
    const auto *pgg_622 = buffer.data(pgg + 622);
    const auto *pgg_624 = buffer.data(pgg + 624);
    const auto *pgg_625 = buffer.data(pgg + 625);
    const auto *pgg_626 = buffer.data(pgg + 626);
    const auto *pgg_627 = buffer.data(pgg + 627);
    const auto *pgg_628 = buffer.data(pgg + 628);
    const auto *pgg_629 = buffer.data(pgg + 629);
    const auto *pgg_630 = buffer.data(pgg + 630);
    const auto *pgg_632 = buffer.data(pgg + 632);
    const auto *pgg_633 = buffer.data(pgg + 633);
    const auto *pgg_635 = buffer.data(pgg + 635);
    const auto *pgg_636 = buffer.data(pgg + 636);
    const auto *pgg_637 = buffer.data(pgg + 637);
    const auto *pgg_639 = buffer.data(pgg + 639);
    const auto *pgg_640 = buffer.data(pgg + 640);
    const auto *pgg_641 = buffer.data(pgg + 641);
    const auto *pgg_642 = buffer.data(pgg + 642);
    const auto *pgg_643 = buffer.data(pgg + 643);
    const auto *pgg_644 = buffer.data(pgg + 644);
    const auto *pgg_645 = buffer.data(pgg + 645);
    const auto *pgg_647 = buffer.data(pgg + 647);
    const auto *pgg_648 = buffer.data(pgg + 648);
    const auto *pgg_650 = buffer.data(pgg + 650);
    const auto *pgg_651 = buffer.data(pgg + 651);
    const auto *pgg_652 = buffer.data(pgg + 652);
    const auto *pgg_655 = buffer.data(pgg + 655);
    const auto *pgg_656 = buffer.data(pgg + 656);
    const auto *pgg_657 = buffer.data(pgg + 657);
    const auto *pgg_658 = buffer.data(pgg + 658);
    const auto *pgg_659 = buffer.data(pgg + 659);
    const auto *pgg_660 = buffer.data(pgg + 660);
    const auto *pgg_662 = buffer.data(pgg + 662);
    const auto *pgg_663 = buffer.data(pgg + 663);
    const auto *pgg_665 = buffer.data(pgg + 665);
    const auto *pgg_666 = buffer.data(pgg + 666);
    const auto *pgg_667 = buffer.data(pgg + 667);
    const auto *pgg_669 = buffer.data(pgg + 669);
    const auto *pgg_670 = buffer.data(pgg + 670);
    const auto *pgg_671 = buffer.data(pgg + 671);
    const auto *pgg_672 = buffer.data(pgg + 672);
    const auto *pgg_673 = buffer.data(pgg + 673);
    const auto *pgg_674 = buffer.data(pgg + 674);

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pb_x, pc_x, pc_y, pfh0_626, pfh0_627, \
                         pfh0_629, pfh1_626, pfh1_627, pfh1_629, \
                         pgg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = pb_x[k] * pfh0_626[k]
                   - f_9 * pc_x[k] * pfh1_626[k];

        t_837[k] = pb_x[k] * pfh0_627[k]
                   - f_9 * pc_x[k] * pfh1_627[k];

        t_838[k] = f_4 * pc_y[k] * pgg_599[k];

        t_839[k] = pb_x[k] * pfh0_629[k]
                   - f_9 * pc_x[k] * pfh1_629[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pa_z, pc_x, pc_y, pc_z, sgh0_210, sgh1_210, \
                         pfg_390, pgf0_402, pgf1_402, pgg_600, \
                         pgg_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = pa_z[k] * sgh0_210[k]
                   - f_9 * pc_z[k] * sgh1_210[k];

        t_841[k] = f_1 * pfg_390[k]
                   + f_4 * pc_y[k] * pgg_600[k];

        t_842[k] = f_15 * pgf0_402[k]
                   - f_16 * pgf1_402[k]
                   + f_4 * pc_x[k] * pgg_602[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pa_z, pc_x, pc_y, pc_z, sgh0_213, sgh1_213, \
                         pfg_392, pgf0_405, pgf1_405, pgg_602, \
                         pgg_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = pa_z[k] * sgh0_213[k]
                   - f_9 * pc_z[k] * sgh1_213[k];

        t_844[k] = f_1 * pfg_392[k]
                   + f_4 * pc_y[k] * pgg_602[k];

        t_845[k] = f_7 * pgf0_405[k]
                   - f_8 * pgf1_405[k]
                   + f_4 * pc_x[k] * pgg_605[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pa_z, pc_x, pc_y, pc_z, sgh0_216, sgh1_216, \
                         pfg_395, pgf0_407, pgf1_407, pgg_605, \
                         pgg_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = pa_z[k] * sgh0_216[k]
                   - f_9 * pc_z[k] * sgh1_216[k];

        t_847[k] = f_5 * pgf0_407[k]
                   - f_6 * pgf1_407[k]
                   + f_4 * pc_x[k] * pgg_607[k];

        t_848[k] = f_1 * pfg_395[k]
                   + f_4 * pc_y[k] * pgg_605[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, t_853, t_854, pc_x, pgf0_409, pgf1_409, \
                         pgg_609, pgg_610, pgg_611, pgg_612, pgg_613, \
                         pgg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_5 * pgf0_409[k]
                   - f_6 * pgf1_409[k]
                   + f_4 * pc_x[k] * pgg_609[k];

        t_850[k] = f_4 * pc_x[k] * pgg_610[k];

        t_851[k] = f_4 * pc_x[k] * pgg_611[k];

        t_852[k] = f_4 * pc_x[k] * pgg_612[k];

        t_853[k] = f_4 * pc_x[k] * pgg_613[k];

        t_854[k] = f_4 * pc_x[k] * pgg_614[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, pa_z, pc_z, sgh0_225, sgh0_226, sgh0_227, \
                         sgg_160, sgg_161, sgh1_225, sgh1_226, \
                         sgh1_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = pa_z[k] * sgh0_225[k]
                   - f_9 * pc_z[k] * sgh1_225[k];

        t_856[k] = pa_z[k] * sgh0_226[k]
                   + f_0 * sgg_160[k]
                   - f_9 * pc_z[k] * sgh1_226[k];

        t_857[k] = pa_z[k] * sgh0_227[k]
                   + f_10 * sgg_161[k]
                   - f_9 * pc_z[k] * sgh1_227[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, pa_z, pc_y, pc_z, sgh0_228, sgh0_230, sgg_162, \
                         sgg_164, sgh1_228, sgh1_230, pfg_404, \
                         pgg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = pa_z[k] * sgh0_228[k]
                   + f_11 * sgg_162[k]
                   - f_9 * pc_z[k] * sgh1_228[k];

        t_859[k] = f_1 * pfg_404[k]
                   + f_4 * pc_y[k] * pgg_614[k];

        t_860[k] = pa_z[k] * sgh0_230[k]
                   + f_14 * sgg_164[k]
                   - f_9 * pc_z[k] * sgh1_230[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, t_864, pc_x, pc_y, pfg_405, pgf0_410, pgf0_412, \
                         pgf0_413, pgf1_410, pgf1_412, pgf1_413, pgg_615, pgg_617, \
                         pgg_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = f_2 * pgf0_410[k]
                   - f_3 * pgf1_410[k]
                   + f_4 * pc_x[k] * pgg_615[k];

        t_862[k] = f_11 * pfg_405[k]
                   + f_4 * pc_y[k] * pgg_615[k];

        t_863[k] = f_15 * pgf0_412[k]
                   - f_16 * pgf1_412[k]
                   + f_4 * pc_x[k] * pgg_617[k];

        t_864[k] = f_7 * pgf0_413[k]
                   - f_8 * pgf1_413[k]
                   + f_4 * pc_x[k] * pgg_618[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, pc_x, pc_y, pfg_407, pgf0_415, pgf0_416, \
                         pgf1_415, pgf1_416, pgg_617, pgg_620, \
                         pgg_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_11 * pfg_407[k]
                   + f_4 * pc_y[k] * pgg_617[k];

        t_866[k] = f_7 * pgf0_415[k]
                   - f_8 * pgf1_415[k]
                   + f_4 * pc_x[k] * pgg_620[k];

        t_867[k] = f_5 * pgf0_416[k]
                   - f_6 * pgf1_416[k]
                   + f_4 * pc_x[k] * pgg_621[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, t_871, pc_x, pc_y, pfg_410, pgf0_417, pgf0_419, \
                         pgf1_417, pgf1_419, pgg_620, pgg_622, pgg_624, \
                         pgg_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = f_5 * pgf0_417[k]
                   - f_6 * pgf1_417[k]
                   + f_4 * pc_x[k] * pgg_622[k];

        t_869[k] = f_11 * pfg_410[k]
                   + f_4 * pc_y[k] * pgg_620[k];

        t_870[k] = f_5 * pgf0_419[k]
                   - f_6 * pgf1_419[k]
                   + f_4 * pc_x[k] * pgg_624[k];

        t_871[k] = f_4 * pc_x[k] * pgg_625[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, t_875, t_876, pc_x, pc_y, pfg_415, pgf0_416, \
                         pgf1_416, pgg_625, pgg_626, pgg_627, pgg_628, \
                         pgg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_4 * pc_x[k] * pgg_626[k];

        t_873[k] = f_4 * pc_x[k] * pgg_627[k];

        t_874[k] = f_4 * pc_x[k] * pgg_628[k];

        t_875[k] = f_4 * pc_x[k] * pgg_629[k];

        t_876[k] = f_11 * pfg_415[k]
                   + f_2 * pgf0_416[k]
                   - f_3 * pgf1_416[k]
                   + f_4 * pc_y[k] * pgg_625[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pc_y, pfg_416, pfg_417, pfg_418, pgf0_417, \
                         pgf0_418, pgf0_419, pgf1_417, pgf1_418, pgf1_419, pgg_626, pgg_627, \
                         pgg_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_11 * pfg_416[k]
                   + f_15 * pgf0_417[k]
                   - f_16 * pgf1_417[k]
                   + f_4 * pc_y[k] * pgg_626[k];

        t_878[k] = f_11 * pfg_417[k]
                   + f_7 * pgf0_418[k]
                   - f_8 * pgf1_418[k]
                   + f_4 * pc_y[k] * pgg_627[k];

        t_879[k] = f_11 * pfg_418[k]
                   + f_5 * pgf0_419[k]
                   - f_6 * pgf1_419[k]
                   + f_4 * pc_y[k] * pgg_628[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, pb_y, pc_x, pc_y, pdh0_356, pdh1_356, pfh0_587, \
                         pfg_419, pfh1_587, pgf0_420, pgf1_420, pgg_629, \
                         pgg_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_11 * pfg_419[k]
                   + f_4 * pc_y[k] * pgg_629[k];

        t_881[k] = f_17 * pdh0_356[k]
                   - f_18 * pdh1_356[k]
                   + pb_y[k] * pfh0_587[k]
                   - f_9 * pc_y[k] * pfh1_587[k];

        t_882[k] = f_2 * pgf0_420[k]
                   - f_3 * pgf1_420[k]
                   + f_4 * pc_x[k] * pgg_630[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, pc_x, pc_y, pfg_420, pfg_422, pgf0_422, \
                         pgf0_423, pgf1_422, pgf1_423, pgg_630, pgg_632, \
                         pgg_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_10 * pfg_420[k]
                   + f_4 * pc_y[k] * pgg_630[k];

        t_884[k] = f_15 * pgf0_422[k]
                   - f_16 * pgf1_422[k]
                   + f_4 * pc_x[k] * pgg_632[k];

        t_885[k] = f_7 * pgf0_423[k]
                   - f_8 * pgf1_423[k]
                   + f_4 * pc_x[k] * pgg_633[k];

        t_886[k] = f_10 * pfg_422[k]
                   + f_4 * pc_y[k] * pgg_632[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, t_890, pc_x, pc_y, pfg_425, pgf0_425, pgf0_426, \
                         pgf0_427, pgf1_425, pgf1_426, pgf1_427, pgg_635, pgg_636, \
                         pgg_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_7 * pgf0_425[k]
                   - f_8 * pgf1_425[k]
                   + f_4 * pc_x[k] * pgg_635[k];

        t_888[k] = f_5 * pgf0_426[k]
                   - f_6 * pgf1_426[k]
                   + f_4 * pc_x[k] * pgg_636[k];

        t_889[k] = f_5 * pgf0_427[k]
                   - f_6 * pgf1_427[k]
                   + f_4 * pc_x[k] * pgg_637[k];

        t_890[k] = f_10 * pfg_425[k]
                   + f_4 * pc_y[k] * pgg_635[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, t_895, t_896, pc_x, pgf0_429, pgf1_429, \
                         pgg_639, pgg_640, pgg_641, pgg_642, pgg_643, \
                         pgg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = f_5 * pgf0_429[k]
                   - f_6 * pgf1_429[k]
                   + f_4 * pc_x[k] * pgg_639[k];

        t_892[k] = f_4 * pc_x[k] * pgg_640[k];

        t_893[k] = f_4 * pc_x[k] * pgg_641[k];

        t_894[k] = f_4 * pc_x[k] * pgg_642[k];

        t_895[k] = f_4 * pc_x[k] * pgg_643[k];

        t_896[k] = f_4 * pc_x[k] * pgg_644[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, pc_y, pfg_430, pfg_431, pfg_432, pgf0_426, \
                         pgf0_427, pgf0_428, pgf1_426, pgf1_427, pgf1_428, pgg_640, pgg_641, \
                         pgg_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = f_10 * pfg_430[k]
                   + f_2 * pgf0_426[k]
                   - f_3 * pgf1_426[k]
                   + f_4 * pc_y[k] * pgg_640[k];

        t_898[k] = f_10 * pfg_431[k]
                   + f_15 * pgf0_427[k]
                   - f_16 * pgf1_427[k]
                   + f_4 * pc_y[k] * pgg_641[k];

        t_899[k] = f_10 * pfg_432[k]
                   + f_7 * pgf0_428[k]
                   - f_8 * pgf1_428[k]
                   + f_4 * pc_y[k] * pgg_642[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, pb_y, pc_y, pdh0_377, pdh1_377, pfh0_608, \
                         pfg_433, pfg_434, pfh1_608, pgf0_429, pgf1_429, pgg_643, \
                         pgg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_10 * pfg_433[k]
                   + f_5 * pgf0_429[k]
                   - f_6 * pgf1_429[k]
                   + f_4 * pc_y[k] * pgg_643[k];

        t_901[k] = f_10 * pfg_434[k]
                   + f_4 * pc_y[k] * pgg_644[k];

        t_902[k] = f_12 * pdh0_377[k]
                   - f_13 * pdh1_377[k]
                   + pb_y[k] * pfh0_608[k]
                   - f_9 * pc_y[k] * pfh1_608[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, t_906, pb_y, pc_x, pc_y, pfh0_609, pfh0_611, \
                         pfg_435, pfh1_609, pfh1_611, pgf0_433, pgf1_433, pgg_645, \
                         pgg_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = pb_y[k] * pfh0_609[k]
                   - f_9 * pc_y[k] * pfh1_609[k];

        t_904[k] = f_0 * pfg_435[k]
                   + f_4 * pc_y[k] * pgg_645[k];

        t_905[k] = pb_y[k] * pfh0_611[k]
                   - f_9 * pc_y[k] * pfh1_611[k];

        t_906[k] = f_7 * pgf0_433[k]
                   - f_8 * pgf1_433[k]
                   + f_4 * pc_x[k] * pgg_648[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, pb_y, pc_x, pc_y, pfh0_614, pfg_437, pfh1_614, \
                         pgf0_436, pgf1_436, pgg_647, pgg_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = f_0 * pfg_437[k]
                   + f_4 * pc_y[k] * pgg_647[k];

        t_908[k] = pb_y[k] * pfh0_614[k]
                   - f_9 * pc_y[k] * pfh1_614[k];

        t_909[k] = f_5 * pgf0_436[k]
                   - f_6 * pgf1_436[k]
                   + f_4 * pc_x[k] * pgg_651[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, pb_y, pc_x, pc_y, pfh0_618, pfg_440, \
                         pfh1_618, pgf0_437, pgf1_437, pgg_650, pgg_652, \
                         pgg_655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = f_5 * pgf0_437[k]
                   - f_6 * pgf1_437[k]
                   + f_4 * pc_x[k] * pgg_652[k];

        t_911[k] = f_0 * pfg_440[k]
                   + f_4 * pc_y[k] * pgg_650[k];

        t_912[k] = pb_y[k] * pfh0_618[k]
                   - f_9 * pc_y[k] * pfh1_618[k];

        t_913[k] = f_4 * pc_x[k] * pgg_655[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, pb_y, pc_x, pc_y, pfh0_624, \
                         pfg_445, pfh1_624, pgg_656, pgg_657, pgg_658, \
                         pgg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_4 * pc_x[k] * pgg_656[k];

        t_915[k] = f_4 * pc_x[k] * pgg_657[k];

        t_916[k] = f_4 * pc_x[k] * pgg_658[k];

        t_917[k] = f_4 * pc_x[k] * pgg_659[k];

        t_918[k] = pb_y[k] * pfh0_624[k]
                   + f_14 * pfg_445[k]
                   - f_9 * pc_y[k] * pfh1_624[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pb_y, pc_y, pfh0_625, pfh0_626, pfh0_627, \
                         pfg_446, pfg_447, pfg_448, pfh1_625, pfh1_626, \
                         pfh1_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = pb_y[k] * pfh0_625[k]
                   + f_1 * pfg_446[k]
                   - f_9 * pc_y[k] * pfh1_625[k];

        t_920[k] = pb_y[k] * pfh0_626[k]
                   + f_11 * pfg_447[k]
                   - f_9 * pc_y[k] * pfh1_626[k];

        t_921[k] = pb_y[k] * pfh0_627[k]
                   + f_10 * pfg_448[k]
                   - f_9 * pc_y[k] * pfh1_627[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, pb_y, pc_x, pc_y, pfh0_629, pfg_449, \
                         pfh1_629, pgf0_440, pgf1_440, pgg_659, \
                         pgg_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_0 * pfg_449[k]
                   + f_4 * pc_y[k] * pgg_659[k];

        t_923[k] = pb_y[k] * pfh0_629[k]
                   - f_9 * pc_y[k] * pfh1_629[k];

        t_924[k] = f_2 * pgf0_440[k]
                   - f_3 * pgf1_440[k]
                   + f_4 * pc_x[k] * pgg_660[k];

        t_925[k] = f_4 * pc_y[k] * pgg_660[k];
    }

#pragma omp simd aligned(t_926, t_927, t_928, t_929, pc_x, pc_y, pgf0_442, pgf0_443, pgf0_445, \
                         pgf1_442, pgf1_443, pgf1_445, pgg_662, pgg_663, \
                         pgg_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_926[k] = f_15 * pgf0_442[k]
                   - f_16 * pgf1_442[k]
                   + f_4 * pc_x[k] * pgg_662[k];

        t_927[k] = f_7 * pgf0_443[k]
                   - f_8 * pgf1_443[k]
                   + f_4 * pc_x[k] * pgg_663[k];

        t_928[k] = f_4 * pc_y[k] * pgg_662[k];

        t_929[k] = f_7 * pgf0_445[k]
                   - f_8 * pgf1_445[k]
                   + f_4 * pc_x[k] * pgg_665[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, pc_x, pc_y, pgf0_446, pgf0_447, pgf0_449, \
                         pgf1_446, pgf1_447, pgf1_449, pgg_665, pgg_666, pgg_667, \
                         pgg_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_5 * pgf0_446[k]
                   - f_6 * pgf1_446[k]
                   + f_4 * pc_x[k] * pgg_666[k];

        t_931[k] = f_5 * pgf0_447[k]
                   - f_6 * pgf1_447[k]
                   + f_4 * pc_x[k] * pgg_667[k];

        t_932[k] = f_4 * pc_y[k] * pgg_665[k];

        t_933[k] = f_5 * pgf0_449[k]
                   - f_6 * pgf1_449[k]
                   + f_4 * pc_x[k] * pgg_669[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, t_938, t_939, pc_x, pc_y, pgf0_446, \
                         pgf1_446, pgg_670, pgg_671, pgg_672, pgg_673, \
                         pgg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_4 * pc_x[k] * pgg_670[k];

        t_935[k] = f_4 * pc_x[k] * pgg_671[k];

        t_936[k] = f_4 * pc_x[k] * pgg_672[k];

        t_937[k] = f_4 * pc_x[k] * pgg_673[k];

        t_938[k] = f_4 * pc_x[k] * pgg_674[k];

        t_939[k] = f_2 * pgf0_446[k]
                   - f_3 * pgf1_446[k]
                   + f_4 * pc_y[k] * pgg_670[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, pc_y, pgf0_447, pgf0_448, pgf0_449, \
                         pgf1_447, pgf1_448, pgf1_449, pgg_671, pgg_672, pgg_673, \
                         pgg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_15 * pgf0_447[k]
                   - f_16 * pgf1_447[k]
                   + f_4 * pc_y[k] * pgg_671[k];

        t_941[k] = f_7 * pgf0_448[k]
                   - f_8 * pgf1_448[k]
                   + f_4 * pc_y[k] * pgg_672[k];

        t_942[k] = f_5 * pgf0_449[k]
                   - f_6 * pgf1_449[k]
                   + f_4 * pc_y[k] * pgg_673[k];

        t_943[k] = f_4 * pc_y[k] * pgg_674[k];
    }

#pragma omp simd aligned(t_944, pc_z, sgg_224, pfg_449, pgf0_449, pgf1_449, \
                         pgg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_0 * sgg_224[k]
                   + f_1 * pfg_449[k]
                   + f_2 * pgf0_449[k]
                   - f_3 * pgf1_449[k]
                   + f_4 * pc_z[k] * pgg_674[k];
    }
}

auto
compute_prim_pgh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sgh0,
                                                   const size_t sgg, const size_t sgh1,
                                                   const size_t pdh0, const size_t pdh1,
                                                   const size_t pfh0, const size_t pfg,
                                                   const size_t pfh1, const size_t pgf0,
                                                   const size_t pgf1, const size_t pgg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pgh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sgg, pdh0,
                                                              pdh1, pfh0, pfg, pfh1, pgf0, pgf1,
                                                              pgg, ncols, gamma, p, q);

    compute_prim_pgh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sgh0,
                                                              sgg, sgh1, pfh0, pfg, pfh1, pgf0,
                                                              pgf1, pgg, ncols, gamma, p, q);

    compute_prim_pgh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sgh0,
                                                              sgg, sgh1, pdh0, pdh1, pfh0, pfg,
                                                              pfh1, pgf0, pgf1, pgg, ncols,
                                                              gamma, p, q);

    compute_prim_pgh_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, sgh0,
                                                              sgg, sgh1, pdh0, pdh1, pfh0, pfg,
                                                              pfh1, pgf0, pgf1, pgg, ncols,
                                                              gamma, p, q);

    compute_prim_pgh_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, sgh0,
                                                              sgg, sgh1, pdh0, pdh1, pfh0, pfg,
                                                              pfh1, pgf0, pgf1, pgg, ncols,
                                                              gamma, p, q);

    compute_prim_pgh_three_center_electron_repulsion_0_piece5(buffer, target, pa, pb, pc, sgh0,
                                                              sgg, sgh1, pdh0, pdh1, pfh0, pfg,
                                                              pfh1, pgf0, pgf1, pgg, ncols,
                                                              gamma, p, q);

    compute_prim_pgh_three_center_electron_repulsion_0_piece6(buffer, target, pa, pb, pc, sgh0,
                                                              sgg, sgh1, pdh0, pdh1, pfh0, pfg,
                                                              pfh1, pgf0, pgf1, pgg, ncols,
                                                              gamma, p, q);

    compute_prim_pgh_three_center_electron_repulsion_0_piece7(buffer, target, pa, pb, pc, sgh0,
                                                              sgg, sgh1, pdh0, pdh1, pfh0, pfg,
                                                              pfh1, pgf0, pgf1, pgg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
