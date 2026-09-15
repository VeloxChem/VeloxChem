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


#include "SimdThreeCenterElectronRepulsionVrrRecPGG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pgg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgf,
                                                          const size_t pdg0, const size_t pdg1,
                                                          const size_t pfg0, const size_t pff,
                                                          const size_t pfg1, const size_t pgd0,
                                                          const size_t pgd1, const size_t pgf,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgf_0 = buffer.data(sgf + 0);
    const auto *sgf_6 = buffer.data(sgf + 6);
    const auto *sgf_9 = buffer.data(sgf + 9);
    const auto *sgf_16 = buffer.data(sgf + 16);
    const auto *sgf_29 = buffer.data(sgf + 29);
    const auto *sgf_36 = buffer.data(sgf + 36);
    const auto *sgf_39 = buffer.data(sgf + 39);
    const auto *sgf_56 = buffer.data(sgf + 56);
    const auto *sgf_59 = buffer.data(sgf + 59);
    const auto *sgf_60 = buffer.data(sgf + 60);
    const auto *sgf_66 = buffer.data(sgf + 66);
    const auto *sgf_69 = buffer.data(sgf + 69);
    const auto *sgf_79 = buffer.data(sgf + 79);
    const auto *sgf_86 = buffer.data(sgf + 86);

    const auto *pdg0_0 = buffer.data(pdg0 + 0);

    const auto *pdg1_0 = buffer.data(pdg1 + 0);

    const auto *pfg0_0 = buffer.data(pfg0 + 0);
    const auto *pfg0_3 = buffer.data(pfg0 + 3);
    const auto *pfg0_5 = buffer.data(pfg0 + 5);
    const auto *pfg0_6 = buffer.data(pfg0 + 6);
    const auto *pfg0_9 = buffer.data(pfg0 + 9);
    const auto *pfg0_10 = buffer.data(pfg0 + 10);
    const auto *pfg0_14 = buffer.data(pfg0 + 14);
    const auto *pfg0_15 = buffer.data(pfg0 + 15);
    const auto *pfg0_18 = buffer.data(pfg0 + 18);
    const auto *pfg0_21 = buffer.data(pfg0 + 21);
    const auto *pfg0_25 = buffer.data(pfg0 + 25);
    const auto *pfg0_30 = buffer.data(pfg0 + 30);
    const auto *pfg0_35 = buffer.data(pfg0 + 35);
    const auto *pfg0_39 = buffer.data(pfg0 + 39);
    const auto *pfg0_44 = buffer.data(pfg0 + 44);
    const auto *pfg0_45 = buffer.data(pfg0 + 45);
    const auto *pfg0_48 = buffer.data(pfg0 + 48);
    const auto *pfg0_51 = buffer.data(pfg0 + 51);
    const auto *pfg0_75 = buffer.data(pfg0 + 75);
    const auto *pfg0_80 = buffer.data(pfg0 + 80);

    const auto *pff_0 = buffer.data(pff + 0);
    const auto *pff_1 = buffer.data(pff + 1);
    const auto *pff_2 = buffer.data(pff + 2);
    const auto *pff_3 = buffer.data(pff + 3);
    const auto *pff_5 = buffer.data(pff + 5);
    const auto *pff_6 = buffer.data(pff + 6);
    const auto *pff_8 = buffer.data(pff + 8);
    const auto *pff_9 = buffer.data(pff + 9);
    const auto *pff_10 = buffer.data(pff + 10);
    const auto *pff_11 = buffer.data(pff + 11);
    const auto *pff_12 = buffer.data(pff + 12);
    const auto *pff_13 = buffer.data(pff + 13);
    const auto *pff_15 = buffer.data(pff + 15);
    const auto *pff_16 = buffer.data(pff + 16);
    const auto *pff_18 = buffer.data(pff + 18);
    const auto *pff_19 = buffer.data(pff + 19);
    const auto *pff_20 = buffer.data(pff + 20);
    const auto *pff_22 = buffer.data(pff + 22);
    const auto *pff_23 = buffer.data(pff + 23);
    const auto *pff_25 = buffer.data(pff + 25);
    const auto *pff_26 = buffer.data(pff + 26);
    const auto *pff_28 = buffer.data(pff + 28);
    const auto *pff_29 = buffer.data(pff + 29);
    const auto *pff_30 = buffer.data(pff + 30);
    const auto *pff_31 = buffer.data(pff + 31);
    const auto *pff_32 = buffer.data(pff + 32);
    const auto *pff_33 = buffer.data(pff + 33);
    const auto *pff_35 = buffer.data(pff + 35);
    const auto *pff_36 = buffer.data(pff + 36);
    const auto *pff_38 = buffer.data(pff + 38);
    const auto *pff_39 = buffer.data(pff + 39);
    const auto *pff_40 = buffer.data(pff + 40);
    const auto *pff_42 = buffer.data(pff + 42);
    const auto *pff_45 = buffer.data(pff + 45);
    const auto *pff_46 = buffer.data(pff + 46);
    const auto *pff_48 = buffer.data(pff + 48);
    const auto *pff_49 = buffer.data(pff + 49);
    const auto *pff_50 = buffer.data(pff + 50);
    const auto *pff_51 = buffer.data(pff + 51);
    const auto *pff_52 = buffer.data(pff + 52);
    const auto *pff_56 = buffer.data(pff + 56);
    const auto *pff_59 = buffer.data(pff + 59);
    const auto *pff_60 = buffer.data(pff + 60);
    const auto *pff_66 = buffer.data(pff + 66);
    const auto *pff_69 = buffer.data(pff + 69);
    const auto *pff_79 = buffer.data(pff + 79);
    const auto *pff_86 = buffer.data(pff + 86);

    const auto *pfg1_0 = buffer.data(pfg1 + 0);
    const auto *pfg1_3 = buffer.data(pfg1 + 3);
    const auto *pfg1_5 = buffer.data(pfg1 + 5);
    const auto *pfg1_6 = buffer.data(pfg1 + 6);
    const auto *pfg1_9 = buffer.data(pfg1 + 9);
    const auto *pfg1_10 = buffer.data(pfg1 + 10);
    const auto *pfg1_14 = buffer.data(pfg1 + 14);
    const auto *pfg1_15 = buffer.data(pfg1 + 15);
    const auto *pfg1_18 = buffer.data(pfg1 + 18);
    const auto *pfg1_21 = buffer.data(pfg1 + 21);
    const auto *pfg1_25 = buffer.data(pfg1 + 25);
    const auto *pfg1_30 = buffer.data(pfg1 + 30);
    const auto *pfg1_35 = buffer.data(pfg1 + 35);
    const auto *pfg1_39 = buffer.data(pfg1 + 39);
    const auto *pfg1_44 = buffer.data(pfg1 + 44);
    const auto *pfg1_45 = buffer.data(pfg1 + 45);
    const auto *pfg1_48 = buffer.data(pfg1 + 48);
    const auto *pfg1_51 = buffer.data(pfg1 + 51);
    const auto *pfg1_75 = buffer.data(pfg1 + 75);
    const auto *pfg1_80 = buffer.data(pfg1 + 80);

    const auto *pgd0_0 = buffer.data(pgd0 + 0);
    const auto *pgd0_3 = buffer.data(pgd0 + 3);
    const auto *pgd0_5 = buffer.data(pgd0 + 5);
    const auto *pgd0_9 = buffer.data(pgd0 + 9);
    const auto *pgd0_11 = buffer.data(pgd0 + 11);
    const auto *pgd0_17 = buffer.data(pgd0 + 17);
    const auto *pgd0_18 = buffer.data(pgd0 + 18);
    const auto *pgd0_21 = buffer.data(pgd0 + 21);
    const auto *pgd0_23 = buffer.data(pgd0 + 23);
    const auto *pgd0_29 = buffer.data(pgd0 + 29);
    const auto *pgd0_30 = buffer.data(pgd0 + 30);
    const auto *pgd0_33 = buffer.data(pgd0 + 33);
    const auto *pgd0_35 = buffer.data(pgd0 + 35);
    const auto *pgd0_36 = buffer.data(pgd0 + 36);
    const auto *pgd0_39 = buffer.data(pgd0 + 39);
    const auto *pgd0_41 = buffer.data(pgd0 + 41);
    const auto *pgd0_42 = buffer.data(pgd0 + 42);
    const auto *pgd0_45 = buffer.data(pgd0 + 45);
    const auto *pgd0_47 = buffer.data(pgd0 + 47);
    const auto *pgd0_48 = buffer.data(pgd0 + 48);

    const auto *pgd1_0 = buffer.data(pgd1 + 0);
    const auto *pgd1_3 = buffer.data(pgd1 + 3);
    const auto *pgd1_5 = buffer.data(pgd1 + 5);
    const auto *pgd1_9 = buffer.data(pgd1 + 9);
    const auto *pgd1_11 = buffer.data(pgd1 + 11);
    const auto *pgd1_17 = buffer.data(pgd1 + 17);
    const auto *pgd1_18 = buffer.data(pgd1 + 18);
    const auto *pgd1_21 = buffer.data(pgd1 + 21);
    const auto *pgd1_23 = buffer.data(pgd1 + 23);
    const auto *pgd1_29 = buffer.data(pgd1 + 29);
    const auto *pgd1_30 = buffer.data(pgd1 + 30);
    const auto *pgd1_33 = buffer.data(pgd1 + 33);
    const auto *pgd1_35 = buffer.data(pgd1 + 35);
    const auto *pgd1_36 = buffer.data(pgd1 + 36);
    const auto *pgd1_39 = buffer.data(pgd1 + 39);
    const auto *pgd1_41 = buffer.data(pgd1 + 41);
    const auto *pgd1_42 = buffer.data(pgd1 + 42);
    const auto *pgd1_45 = buffer.data(pgd1 + 45);
    const auto *pgd1_47 = buffer.data(pgd1 + 47);
    const auto *pgd1_48 = buffer.data(pgd1 + 48);

    const auto *pgf_0 = buffer.data(pgf + 0);
    const auto *pgf_1 = buffer.data(pgf + 1);
    const auto *pgf_2 = buffer.data(pgf + 2);
    const auto *pgf_3 = buffer.data(pgf + 3);
    const auto *pgf_5 = buffer.data(pgf + 5);
    const auto *pgf_6 = buffer.data(pgf + 6);
    const auto *pgf_8 = buffer.data(pgf + 8);
    const auto *pgf_9 = buffer.data(pgf + 9);
    const auto *pgf_10 = buffer.data(pgf + 10);
    const auto *pgf_12 = buffer.data(pgf + 12);
    const auto *pgf_13 = buffer.data(pgf + 13);
    const auto *pgf_15 = buffer.data(pgf + 15);
    const auto *pgf_16 = buffer.data(pgf + 16);
    const auto *pgf_18 = buffer.data(pgf + 18);
    const auto *pgf_19 = buffer.data(pgf + 19);
    const auto *pgf_20 = buffer.data(pgf + 20);
    const auto *pgf_22 = buffer.data(pgf + 22);
    const auto *pgf_23 = buffer.data(pgf + 23);
    const auto *pgf_25 = buffer.data(pgf + 25);
    const auto *pgf_26 = buffer.data(pgf + 26);
    const auto *pgf_28 = buffer.data(pgf + 28);
    const auto *pgf_29 = buffer.data(pgf + 29);
    const auto *pgf_30 = buffer.data(pgf + 30);
    const auto *pgf_31 = buffer.data(pgf + 31);
    const auto *pgf_32 = buffer.data(pgf + 32);
    const auto *pgf_33 = buffer.data(pgf + 33);
    const auto *pgf_35 = buffer.data(pgf + 35);
    const auto *pgf_36 = buffer.data(pgf + 36);
    const auto *pgf_38 = buffer.data(pgf + 38);
    const auto *pgf_39 = buffer.data(pgf + 39);
    const auto *pgf_40 = buffer.data(pgf + 40);
    const auto *pgf_42 = buffer.data(pgf + 42);
    const auto *pgf_43 = buffer.data(pgf + 43);
    const auto *pgf_45 = buffer.data(pgf + 45);
    const auto *pgf_46 = buffer.data(pgf + 46);
    const auto *pgf_48 = buffer.data(pgf + 48);
    const auto *pgf_49 = buffer.data(pgf + 49);
    const auto *pgf_50 = buffer.data(pgf + 50);
    const auto *pgf_51 = buffer.data(pgf + 51);
    const auto *pgf_52 = buffer.data(pgf + 52);
    const auto *pgf_53 = buffer.data(pgf + 53);
    const auto *pgf_55 = buffer.data(pgf + 55);
    const auto *pgf_56 = buffer.data(pgf + 56);
    const auto *pgf_58 = buffer.data(pgf + 58);
    const auto *pgf_59 = buffer.data(pgf + 59);
    const auto *pgf_60 = buffer.data(pgf + 60);
    const auto *pgf_61 = buffer.data(pgf + 61);
    const auto *pgf_62 = buffer.data(pgf + 62);
    const auto *pgf_63 = buffer.data(pgf + 63);
    const auto *pgf_65 = buffer.data(pgf + 65);
    const auto *pgf_66 = buffer.data(pgf + 66);
    const auto *pgf_68 = buffer.data(pgf + 68);
    const auto *pgf_69 = buffer.data(pgf + 69);
    const auto *pgf_70 = buffer.data(pgf + 70);
    const auto *pgf_72 = buffer.data(pgf + 72);
    const auto *pgf_73 = buffer.data(pgf + 73);
    const auto *pgf_75 = buffer.data(pgf + 75);
    const auto *pgf_76 = buffer.data(pgf + 76);
    const auto *pgf_78 = buffer.data(pgf + 78);
    const auto *pgf_79 = buffer.data(pgf + 79);
    const auto *pgf_80 = buffer.data(pgf + 80);
    const auto *pgf_81 = buffer.data(pgf + 81);
    const auto *pgf_82 = buffer.data(pgf + 82);
    const auto *pgf_86 = buffer.data(pgf + 86);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sgf_0, pff_0, pgd0_0, \
                         pgd1_0, pgf_0, pgf_1, pgf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgf_0[k]
                 + f_1 * pff_0[k]
                 + f_2 * pgd0_0[k]
                 - f_3 * pgd1_0[k]
                 + f_4 * pc_x[k] * pgf_0[k];

        t_1[k] = f_4 * pc_y[k] * pgf_0[k];

        t_2[k] = f_4 * pc_z[k] * pgf_0[k];

        t_3[k] = f_5 * pgd0_0[k]
                 - f_6 * pgd1_0[k]
                 + f_4 * pc_y[k] * pgf_1[k];

        t_4[k] = f_4 * pc_y[k] * pgf_2[k];

        t_5[k] = f_5 * pgd0_0[k]
                 - f_6 * pgd1_0[k]
                 + f_4 * pc_z[k] * pgf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, sgf_6, sgf_9, pff_6, pff_9, \
                         pgf_3, pgf_5, pgf_6, pgf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * sgf_6[k]
                 + f_1 * pff_6[k]
                 + f_4 * pc_x[k] * pgf_6[k];

        t_7[k] = f_4 * pc_z[k] * pgf_3[k];

        t_8[k] = f_4 * pc_y[k] * pgf_5[k];

        t_9[k] = f_0 * sgf_9[k]
                 + f_1 * pff_9[k]
                 + f_4 * pc_x[k] * pgf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pc_y, pc_z, pgd0_3, pgd0_5, pgd1_3, \
                         pgd1_5, pgf_6, pgf_8, pgf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * pgd0_3[k]
                  - f_3 * pgd1_3[k]
                  + f_4 * pc_y[k] * pgf_6[k];

        t_11[k] = f_4 * pc_z[k] * pgf_6[k];

        t_12[k] = f_5 * pgd0_5[k]
                  - f_6 * pgd1_5[k]
                  + f_4 * pc_y[k] * pgf_8[k];

        t_13[k] = f_4 * pc_y[k] * pgf_9[k];

        t_14[k] = f_2 * pgd0_5[k]
                  - f_3 * pgd1_5[k]
                  + f_4 * pc_z[k] * pgf_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pc_y, pc_z, pfg0_0, pfg0_3, pff_0, \
                         pff_1, pfg1_0, pfg1_3, pgf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * pfg0_0[k]
                  - f_7 * pc_y[k] * pfg1_0[k];

        t_16[k] = f_0 * pff_0[k]
                  + f_4 * pc_y[k] * pgf_10[k];

        t_17[k] = f_4 * pc_z[k] * pgf_10[k];

        t_18[k] = pb_y[k] * pfg0_3[k]
                  + f_8 * pff_1[k]
                  - f_7 * pc_y[k] * pfg1_3[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_y, pc_x, pc_y, pc_z, sgf_16, pfg0_5, \
                         pff_2, pff_16, pfg1_5, pgf_12, pgf_13, \
                         pgf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * pff_2[k]
                  + f_4 * pc_y[k] * pgf_12[k];

        t_20[k] = pb_y[k] * pfg0_5[k]
                  - f_7 * pc_y[k] * pfg1_5[k];

        t_21[k] = f_0 * sgf_16[k]
                  + f_9 * pff_16[k]
                  + f_4 * pc_x[k] * pgf_16[k];

        t_22[k] = f_4 * pc_z[k] * pgf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_y, pc_z, pfg0_9, pff_5, pff_6, \
                         pfg1_9, pgd0_9, pgd1_9, pgf_15, pgf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * pff_5[k]
                  + f_4 * pc_y[k] * pgf_15[k];

        t_24[k] = pb_y[k] * pfg0_9[k]
                  - f_7 * pc_y[k] * pfg1_9[k];

        t_25[k] = f_0 * pff_6[k]
                  + f_2 * pgd0_9[k]
                  - f_3 * pgd1_9[k]
                  + f_4 * pc_y[k] * pgf_16[k];

        t_26[k] = f_4 * pc_z[k] * pgf_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pc_y, pfg0_14, pff_8, pff_9, pfg1_14, \
                         pgd0_11, pgd1_11, pgf_18, pgf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * pff_8[k]
                  + f_5 * pgd0_11[k]
                  - f_6 * pgd1_11[k]
                  + f_4 * pc_y[k] * pgf_18[k];

        t_28[k] = f_0 * pff_9[k]
                  + f_4 * pc_y[k] * pgf_19[k];

        t_29[k] = pb_y[k] * pfg0_14[k]
                  - f_7 * pc_y[k] * pfg1_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_z, pc_y, pc_z, pfg0_0, pfg0_3, \
                         pff_0, pfg1_0, pfg1_3, pgf_20, pgf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * pfg0_0[k]
                  - f_7 * pc_z[k] * pfg1_0[k];

        t_31[k] = f_4 * pc_y[k] * pgf_20[k];

        t_32[k] = f_0 * pff_0[k]
                  + f_4 * pc_z[k] * pgf_20[k];

        t_33[k] = pb_z[k] * pfg0_3[k]
                  - f_7 * pc_z[k] * pfg1_3[k];

        t_34[k] = f_4 * pc_y[k] * pgf_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_z, pc_y, pc_z, pfg0_5, pfg0_6, pff_2, \
                         pff_3, pfg1_5, pfg1_6, pgf_23, pgf_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_z[k] * pfg0_5[k]
                  + f_8 * pff_2[k]
                  - f_7 * pc_z[k] * pfg1_5[k];

        t_36[k] = pb_z[k] * pfg0_6[k]
                  - f_7 * pc_z[k] * pfg1_6[k];

        t_37[k] = f_0 * pff_3[k]
                  + f_4 * pc_z[k] * pgf_23[k];

        t_38[k] = f_4 * pc_y[k] * pgf_25[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_z, pc_x, pc_z, sgf_29, pfg0_10, pff_6, pff_29, \
                         pfg1_10, pgf_26, pgf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * sgf_29[k]
                  + f_9 * pff_29[k]
                  + f_4 * pc_x[k] * pgf_29[k];

        t_40[k] = pb_z[k] * pfg0_10[k]
                  - f_7 * pc_z[k] * pfg1_10[k];

        t_41[k] = f_0 * pff_6[k]
                  + f_4 * pc_z[k] * pgf_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_y, pc_z, pdg0_0, pdg1_0, pfg0_15, \
                         pff_9, pfg1_15, pgd0_17, pgd1_17, pgf_28, \
                         pgf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * pgd0_17[k]
                  - f_6 * pgd1_17[k]
                  + f_4 * pc_y[k] * pgf_28[k];

        t_43[k] = f_4 * pc_y[k] * pgf_29[k];

        t_44[k] = f_0 * pff_9[k]
                  + f_2 * pgd0_17[k]
                  - f_3 * pgd1_17[k]
                  + f_4 * pc_z[k] * pgf_29[k];

        t_45[k] = f_10 * pdg0_0[k]
                  - f_11 * pdg1_0[k]
                  + pb_y[k] * pfg0_15[k]
                  - f_7 * pc_y[k] * pfg1_15[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pc_y, pc_z, pff_10, pff_11, pff_12, \
                         pgd0_18, pgd1_18, pgf_30, pgf_31, pgf_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_8 * pff_10[k]
                  + f_4 * pc_y[k] * pgf_30[k];

        t_47[k] = f_4 * pc_z[k] * pgf_30[k];

        t_48[k] = f_8 * pff_11[k]
                  + f_5 * pgd0_18[k]
                  - f_6 * pgd1_18[k]
                  + f_4 * pc_y[k] * pgf_31[k];

        t_49[k] = f_8 * pff_12[k]
                  + f_4 * pc_y[k] * pgf_32[k];

        t_50[k] = f_5 * pgd0_18[k]
                  - f_6 * pgd1_18[k]
                  + f_4 * pc_z[k] * pgf_32[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pc_x, pc_y, pc_z, sgf_36, sgf_39, pff_15, \
                         pff_36, pff_39, pgf_33, pgf_35, pgf_36, \
                         pgf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * sgf_36[k]
                  + f_8 * pff_36[k]
                  + f_4 * pc_x[k] * pgf_36[k];

        t_52[k] = f_4 * pc_z[k] * pgf_33[k];

        t_53[k] = f_8 * pff_15[k]
                  + f_4 * pc_y[k] * pgf_35[k];

        t_54[k] = f_0 * sgf_39[k]
                  + f_8 * pff_39[k]
                  + f_4 * pc_x[k] * pgf_39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pc_y, pc_z, pff_16, pff_18, pff_19, pgd0_21, \
                         pgd0_23, pgd1_21, pgd1_23, pgf_36, pgf_38, \
                         pgf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_8 * pff_16[k]
                  + f_2 * pgd0_21[k]
                  - f_3 * pgd1_21[k]
                  + f_4 * pc_y[k] * pgf_36[k];

        t_56[k] = f_4 * pc_z[k] * pgf_36[k];

        t_57[k] = f_8 * pff_18[k]
                  + f_5 * pgd0_23[k]
                  - f_6 * pgd1_23[k]
                  + f_4 * pc_y[k] * pgf_38[k];

        t_58[k] = f_8 * pff_19[k]
                  + f_4 * pc_y[k] * pgf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pb_y, pc_y, pc_z, pfg0_30, pff_10, pff_20, \
                         pfg1_30, pgd0_23, pgd1_23, pgf_39, pgf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_2 * pgd0_23[k]
                  - f_3 * pgd1_23[k]
                  + f_4 * pc_z[k] * pgf_39[k];

        t_60[k] = pb_y[k] * pfg0_30[k]
                  - f_7 * pc_y[k] * pfg1_30[k];

        t_61[k] = f_0 * pff_20[k]
                  + f_4 * pc_y[k] * pgf_40[k];

        t_62[k] = f_0 * pff_10[k]
                  + f_4 * pc_z[k] * pgf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pb_y, pb_z, pc_y, pc_z, pfg0_18, pfg0_21, \
                         pfg0_35, pff_22, pfg1_18, pfg1_21, pfg1_35, \
                         pgf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_z[k] * pfg0_18[k]
                  - f_7 * pc_z[k] * pfg1_18[k];

        t_64[k] = f_0 * pff_22[k]
                  + f_4 * pc_y[k] * pgf_42[k];

        t_65[k] = pb_y[k] * pfg0_35[k]
                  - f_7 * pc_y[k] * pfg1_35[k];

        t_66[k] = pb_z[k] * pfg0_21[k]
                  - f_7 * pc_z[k] * pfg1_21[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_y, pb_z, pc_y, pc_z, pfg0_25, pfg0_39, \
                         pff_13, pff_25, pfg1_25, pfg1_39, pgf_43, \
                         pgf_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_0 * pff_13[k]
                  + f_4 * pc_z[k] * pgf_43[k];

        t_68[k] = f_0 * pff_25[k]
                  + f_4 * pc_y[k] * pgf_45[k];

        t_69[k] = pb_y[k] * pfg0_39[k]
                  - f_7 * pc_y[k] * pfg1_39[k];

        t_70[k] = pb_z[k] * pfg0_25[k]
                  - f_7 * pc_z[k] * pfg1_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pc_y, pc_z, pff_16, pff_28, pff_29, pgd0_29, \
                         pgd1_29, pgf_46, pgf_48, pgf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_0 * pff_16[k]
                  + f_4 * pc_z[k] * pgf_46[k];

        t_72[k] = f_0 * pff_28[k]
                  + f_5 * pgd0_29[k]
                  - f_6 * pgd1_29[k]
                  + f_4 * pc_y[k] * pgf_48[k];

        t_73[k] = f_0 * pff_29[k]
                  + f_4 * pc_y[k] * pgf_49[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_y, pb_z, pc_y, pc_z, pdg0_0, pdg1_0, \
                         pfg0_30, pfg0_44, pff_20, pfg1_30, pfg1_44, \
                         pgf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pb_y[k] * pfg0_44[k]
                  - f_7 * pc_y[k] * pfg1_44[k];

        t_75[k] = f_10 * pdg0_0[k]
                  - f_11 * pdg1_0[k]
                  + pb_z[k] * pfg0_30[k]
                  - f_7 * pc_z[k] * pfg1_30[k];

        t_76[k] = f_4 * pc_y[k] * pgf_50[k];

        t_77[k] = f_8 * pff_20[k]
                  + f_4 * pc_z[k] * pgf_50[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pc_x, pc_y, pc_z, sgf_56, pff_22, pff_56, \
                         pgd0_30, pgd1_30, pgf_51, pgf_52, pgf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_5 * pgd0_30[k]
                  - f_6 * pgd1_30[k]
                  + f_4 * pc_y[k] * pgf_51[k];

        t_79[k] = f_4 * pc_y[k] * pgf_52[k];

        t_80[k] = f_8 * pff_22[k]
                  + f_5 * pgd0_30[k]
                  - f_6 * pgd1_30[k]
                  + f_4 * pc_z[k] * pgf_52[k];

        t_81[k] = f_0 * sgf_56[k]
                  + f_8 * pff_56[k]
                  + f_4 * pc_x[k] * pgf_56[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pc_x, pc_y, pc_z, sgf_59, pff_23, pff_59, \
                         pgd0_33, pgd1_33, pgf_53, pgf_55, pgf_56, \
                         pgf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_8 * pff_23[k]
                  + f_4 * pc_z[k] * pgf_53[k];

        t_83[k] = f_4 * pc_y[k] * pgf_55[k];

        t_84[k] = f_0 * sgf_59[k]
                  + f_8 * pff_59[k]
                  + f_4 * pc_x[k] * pgf_59[k];

        t_85[k] = f_2 * pgd0_33[k]
                  - f_3 * pgd1_33[k]
                  + f_4 * pc_y[k] * pgf_56[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_y, pc_z, pff_26, pff_29, pgd0_35, pgd1_35, \
                         pgf_56, pgf_58, pgf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_8 * pff_26[k]
                  + f_4 * pc_z[k] * pgf_56[k];

        t_87[k] = f_5 * pgd0_35[k]
                  - f_6 * pgd1_35[k]
                  + f_4 * pc_y[k] * pgf_58[k];

        t_88[k] = f_4 * pc_y[k] * pgf_59[k];

        t_89[k] = f_8 * pff_29[k]
                  + f_2 * pgd0_35[k]
                  - f_3 * pgd1_35[k]
                  + f_4 * pc_z[k] * pgf_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, pc_z, sgf_60, pff_30, pff_31, \
                         pff_60, pgd0_36, pgd1_36, pgf_60, pgf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * sgf_60[k]
                  + f_0 * pff_60[k]
                  + f_2 * pgd0_36[k]
                  - f_3 * pgd1_36[k]
                  + f_4 * pc_x[k] * pgf_60[k];

        t_91[k] = f_9 * pff_30[k]
                  + f_4 * pc_y[k] * pgf_60[k];

        t_92[k] = f_4 * pc_z[k] * pgf_60[k];

        t_93[k] = f_9 * pff_31[k]
                  + f_5 * pgd0_36[k]
                  - f_6 * pgd1_36[k]
                  + f_4 * pc_y[k] * pgf_61[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pc_x, pc_y, pc_z, sgf_66, pff_32, pff_66, \
                         pgd0_36, pgd1_36, pgf_62, pgf_63, pgf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * pff_32[k]
                  + f_4 * pc_y[k] * pgf_62[k];

        t_95[k] = f_5 * pgd0_36[k]
                  - f_6 * pgd1_36[k]
                  + f_4 * pc_z[k] * pgf_62[k];

        t_96[k] = f_0 * sgf_66[k]
                  + f_0 * pff_66[k]
                  + f_4 * pc_x[k] * pgf_66[k];

        t_97[k] = f_4 * pc_z[k] * pgf_63[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pc_x, pc_y, pc_z, sgf_69, pff_35, pff_36, \
                         pff_69, pgd0_39, pgd1_39, pgf_65, pgf_66, \
                         pgf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_9 * pff_35[k]
                  + f_4 * pc_y[k] * pgf_65[k];

        t_99[k] = f_0 * sgf_69[k]
                  + f_0 * pff_69[k]
                  + f_4 * pc_x[k] * pgf_69[k];

        t_100[k] = f_9 * pff_36[k]
                   + f_2 * pgd0_39[k]
                   - f_3 * pgd1_39[k]
                   + f_4 * pc_y[k] * pgf_66[k];

        t_101[k] = f_4 * pc_z[k] * pgf_66[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pb_z, pc_y, pc_z, pfg0_45, pff_38, \
                         pff_39, pfg1_45, pgd0_41, pgd1_41, pgf_68, \
                         pgf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * pff_38[k]
                   + f_5 * pgd0_41[k]
                   - f_6 * pgd1_41[k]
                   + f_4 * pc_y[k] * pgf_68[k];

        t_103[k] = f_9 * pff_39[k]
                   + f_4 * pc_y[k] * pgf_69[k];

        t_104[k] = f_2 * pgd0_41[k]
                   - f_3 * pgd1_41[k]
                   + f_4 * pc_z[k] * pgf_69[k];

        t_105[k] = pb_z[k] * pfg0_45[k]
                   - f_7 * pc_z[k] * pfg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, pfg0_48, pff_30, \
                         pff_40, pff_42, pfg1_48, pgf_70, pgf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * pff_40[k]
                   + f_4 * pc_y[k] * pgf_70[k];

        t_107[k] = f_0 * pff_30[k]
                   + f_4 * pc_z[k] * pgf_70[k];

        t_108[k] = pb_z[k] * pfg0_48[k]
                   - f_7 * pc_z[k] * pfg1_48[k];

        t_109[k] = f_8 * pff_42[k]
                   + f_4 * pc_y[k] * pgf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_z, pc_z, pfg0_51, pff_32, pff_33, pfg1_51, \
                         pgd0_42, pgd1_42, pgf_72, pgf_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_0 * pff_32[k]
                   + f_5 * pgd0_42[k]
                   - f_6 * pgd1_42[k]
                   + f_4 * pc_z[k] * pgf_72[k];

        t_111[k] = pb_z[k] * pfg0_51[k]
                   - f_7 * pc_z[k] * pfg1_51[k];

        t_112[k] = f_0 * pff_33[k]
                   + f_4 * pc_z[k] * pgf_73[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pc_x, pc_y, sgf_79, pff_45, pff_46, pff_79, \
                         pgd0_45, pgd1_45, pgf_75, pgf_76, pgf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * pff_45[k]
                   + f_4 * pc_y[k] * pgf_75[k];

        t_114[k] = f_0 * sgf_79[k]
                   + f_0 * pff_79[k]
                   + f_4 * pc_x[k] * pgf_79[k];

        t_115[k] = f_8 * pff_46[k]
                   + f_2 * pgd0_45[k]
                   - f_3 * pgd1_45[k]
                   + f_4 * pc_y[k] * pgf_76[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pc_y, pc_z, pff_36, pff_39, pff_48, \
                         pff_49, pgd0_47, pgd1_47, pgf_76, pgf_78, \
                         pgf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_0 * pff_36[k]
                   + f_4 * pc_z[k] * pgf_76[k];

        t_117[k] = f_8 * pff_48[k]
                   + f_5 * pgd0_47[k]
                   - f_6 * pgd1_47[k]
                   + f_4 * pc_y[k] * pgf_78[k];

        t_118[k] = f_8 * pff_49[k]
                   + f_4 * pc_y[k] * pgf_79[k];

        t_119[k] = f_0 * pff_39[k]
                   + f_2 * pgd0_47[k]
                   - f_3 * pgd1_47[k]
                   + f_4 * pc_z[k] * pgf_79[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pc_y, pc_z, pfg0_75, pff_40, \
                         pff_50, pff_51, pfg1_75, pgd0_48, pgd1_48, pgf_80, \
                         pgf_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pb_y[k] * pfg0_75[k]
                   - f_7 * pc_y[k] * pfg1_75[k];

        t_121[k] = f_0 * pff_50[k]
                   + f_4 * pc_y[k] * pgf_80[k];

        t_122[k] = f_8 * pff_40[k]
                   + f_4 * pc_z[k] * pgf_80[k];

        t_123[k] = f_0 * pff_51[k]
                   + f_5 * pgd0_48[k]
                   - f_6 * pgd1_48[k]
                   + f_4 * pc_y[k] * pgf_81[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pb_y, pc_x, pc_y, sgf_86, pfg0_80, pff_52, \
                         pff_86, pfg1_80, pgf_82, pgf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_0 * pff_52[k]
                   + f_4 * pc_y[k] * pgf_82[k];

        t_125[k] = pb_y[k] * pfg0_80[k]
                   - f_7 * pc_y[k] * pfg1_80[k];

        t_126[k] = f_0 * sgf_86[k]
                   + f_0 * pff_86[k]
                   + f_4 * pc_x[k] * pgf_86[k];
    }
}

static auto
compute_prim_pgg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgg0, const size_t sgf,
                                                          const size_t sgg1, const size_t pfg0,
                                                          const size_t pff, const size_t pfg1,
                                                          const size_t pgd0, const size_t pgd1,
                                                          const size_t pgf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgg0_0 = buffer.data(sgg0 + 0);
    const auto *sgg0_1 = buffer.data(sgg0 + 1);
    const auto *sgg0_3 = buffer.data(sgg0 + 3);
    const auto *sgg0_5 = buffer.data(sgg0 + 5);
    const auto *sgg0_10 = buffer.data(sgg0 + 10);
    const auto *sgg0_14 = buffer.data(sgg0 + 14);
    const auto *sgg0_150 = buffer.data(sgg0 + 150);
    const auto *sgg0_153 = buffer.data(sgg0 + 153);
    const auto *sgg0_160 = buffer.data(sgg0 + 160);
    const auto *sgg0_162 = buffer.data(sgg0 + 162);
    const auto *sgg0_164 = buffer.data(sgg0 + 164);
    const auto *sgg0_170 = buffer.data(sgg0 + 170);
    const auto *sgg0_175 = buffer.data(sgg0 + 175);
    const auto *sgg0_177 = buffer.data(sgg0 + 177);
    const auto *sgg0_179 = buffer.data(sgg0 + 179);
    const auto *sgg0_180 = buffer.data(sgg0 + 180);
    const auto *sgg0_183 = buffer.data(sgg0 + 183);
    const auto *sgg0_185 = buffer.data(sgg0 + 185);
    const auto *sgg0_190 = buffer.data(sgg0 + 190);
    const auto *sgg0_192 = buffer.data(sgg0 + 192);
    const auto *sgg0_194 = buffer.data(sgg0 + 194);
    const auto *sgg0_198 = buffer.data(sgg0 + 198);
    const auto *sgg0_205 = buffer.data(sgg0 + 205);
    const auto *sgg0_207 = buffer.data(sgg0 + 207);
    const auto *sgg0_209 = buffer.data(sgg0 + 209);
    const auto *sgg0_210 = buffer.data(sgg0 + 210);
    const auto *sgg0_215 = buffer.data(sgg0 + 215);
    const auto *sgg0_220 = buffer.data(sgg0 + 220);
    const auto *sgg0_222 = buffer.data(sgg0 + 222);
    const auto *sgg0_224 = buffer.data(sgg0 + 224);

    const auto *sgf_0 = buffer.data(sgf + 0);
    const auto *sgf_1 = buffer.data(sgf + 1);
    const auto *sgf_6 = buffer.data(sgf + 6);
    const auto *sgf_9 = buffer.data(sgf + 9);
    const auto *sgf_90 = buffer.data(sgf + 90);
    const auto *sgf_96 = buffer.data(sgf + 96);
    const auto *sgf_99 = buffer.data(sgf + 99);
    const auto *sgf_100 = buffer.data(sgf + 100);
    const auto *sgf_103 = buffer.data(sgf + 103);
    const auto *sgf_106 = buffer.data(sgf + 106);
    const auto *sgf_109 = buffer.data(sgf + 109);
    const auto *sgf_115 = buffer.data(sgf + 115);
    const auto *sgf_116 = buffer.data(sgf + 116);
    const auto *sgf_119 = buffer.data(sgf + 119);
    const auto *sgf_120 = buffer.data(sgf + 120);
    const auto *sgf_123 = buffer.data(sgf + 123);
    const auto *sgf_125 = buffer.data(sgf + 125);
    const auto *sgf_126 = buffer.data(sgf + 126);
    const auto *sgf_129 = buffer.data(sgf + 129);
    const auto *sgf_133 = buffer.data(sgf + 133);
    const auto *sgf_136 = buffer.data(sgf + 136);
    const auto *sgf_139 = buffer.data(sgf + 139);
    const auto *sgf_140 = buffer.data(sgf + 140);
    const auto *sgf_145 = buffer.data(sgf + 145);
    const auto *sgf_146 = buffer.data(sgf + 146);
    const auto *sgf_149 = buffer.data(sgf + 149);

    const auto *sgg1_0 = buffer.data(sgg1 + 0);
    const auto *sgg1_1 = buffer.data(sgg1 + 1);
    const auto *sgg1_3 = buffer.data(sgg1 + 3);
    const auto *sgg1_5 = buffer.data(sgg1 + 5);
    const auto *sgg1_10 = buffer.data(sgg1 + 10);
    const auto *sgg1_14 = buffer.data(sgg1 + 14);
    const auto *sgg1_150 = buffer.data(sgg1 + 150);
    const auto *sgg1_153 = buffer.data(sgg1 + 153);
    const auto *sgg1_160 = buffer.data(sgg1 + 160);
    const auto *sgg1_162 = buffer.data(sgg1 + 162);
    const auto *sgg1_164 = buffer.data(sgg1 + 164);
    const auto *sgg1_170 = buffer.data(sgg1 + 170);
    const auto *sgg1_175 = buffer.data(sgg1 + 175);
    const auto *sgg1_177 = buffer.data(sgg1 + 177);
    const auto *sgg1_179 = buffer.data(sgg1 + 179);
    const auto *sgg1_180 = buffer.data(sgg1 + 180);
    const auto *sgg1_183 = buffer.data(sgg1 + 183);
    const auto *sgg1_185 = buffer.data(sgg1 + 185);
    const auto *sgg1_190 = buffer.data(sgg1 + 190);
    const auto *sgg1_192 = buffer.data(sgg1 + 192);
    const auto *sgg1_194 = buffer.data(sgg1 + 194);
    const auto *sgg1_198 = buffer.data(sgg1 + 198);
    const auto *sgg1_205 = buffer.data(sgg1 + 205);
    const auto *sgg1_207 = buffer.data(sgg1 + 207);
    const auto *sgg1_209 = buffer.data(sgg1 + 209);
    const auto *sgg1_210 = buffer.data(sgg1 + 210);
    const auto *sgg1_215 = buffer.data(sgg1 + 215);
    const auto *sgg1_220 = buffer.data(sgg1 + 220);
    const auto *sgg1_222 = buffer.data(sgg1 + 222);
    const auto *sgg1_224 = buffer.data(sgg1 + 224);

    const auto *pfg0_84 = buffer.data(pfg0 + 84);
    const auto *pfg0_90 = buffer.data(pfg0 + 90);
    const auto *pfg0_93 = buffer.data(pfg0 + 93);
    const auto *pfg0_135 = buffer.data(pfg0 + 135);
    const auto *pfg0_140 = buffer.data(pfg0 + 140);

    const auto *pff_43 = buffer.data(pff + 43);
    const auto *pff_46 = buffer.data(pff + 46);
    const auto *pff_49 = buffer.data(pff + 49);
    const auto *pff_50 = buffer.data(pff + 50);
    const auto *pff_52 = buffer.data(pff + 52);
    const auto *pff_53 = buffer.data(pff + 53);
    const auto *pff_55 = buffer.data(pff + 55);
    const auto *pff_56 = buffer.data(pff + 56);
    const auto *pff_58 = buffer.data(pff + 58);
    const auto *pff_59 = buffer.data(pff + 59);
    const auto *pff_60 = buffer.data(pff + 60);
    const auto *pff_62 = buffer.data(pff + 62);
    const auto *pff_63 = buffer.data(pff + 63);
    const auto *pff_65 = buffer.data(pff + 65);
    const auto *pff_66 = buffer.data(pff + 66);
    const auto *pff_69 = buffer.data(pff + 69);
    const auto *pff_70 = buffer.data(pff + 70);
    const auto *pff_72 = buffer.data(pff + 72);
    const auto *pff_73 = buffer.data(pff + 73);
    const auto *pff_75 = buffer.data(pff + 75);
    const auto *pff_76 = buffer.data(pff + 76);
    const auto *pff_79 = buffer.data(pff + 79);
    const auto *pff_80 = buffer.data(pff + 80);
    const auto *pff_82 = buffer.data(pff + 82);
    const auto *pff_83 = buffer.data(pff + 83);
    const auto *pff_85 = buffer.data(pff + 85);
    const auto *pff_86 = buffer.data(pff + 86);
    const auto *pff_89 = buffer.data(pff + 89);
    const auto *pff_90 = buffer.data(pff + 90);
    const auto *pff_92 = buffer.data(pff + 92);
    const auto *pff_93 = buffer.data(pff + 93);
    const auto *pff_95 = buffer.data(pff + 95);
    const auto *pff_96 = buffer.data(pff + 96);
    const auto *pff_99 = buffer.data(pff + 99);
    const auto *pff_106 = buffer.data(pff + 106);
    const auto *pff_107 = buffer.data(pff + 107);
    const auto *pff_108 = buffer.data(pff + 108);
    const auto *pff_109 = buffer.data(pff + 109);
    const auto *pff_110 = buffer.data(pff + 110);
    const auto *pff_111 = buffer.data(pff + 111);
    const auto *pff_113 = buffer.data(pff + 113);
    const auto *pff_115 = buffer.data(pff + 115);
    const auto *pff_116 = buffer.data(pff + 116);
    const auto *pff_117 = buffer.data(pff + 117);
    const auto *pff_118 = buffer.data(pff + 118);
    const auto *pff_119 = buffer.data(pff + 119);

    const auto *pfg1_84 = buffer.data(pfg1 + 84);
    const auto *pfg1_90 = buffer.data(pfg1 + 90);
    const auto *pfg1_93 = buffer.data(pfg1 + 93);
    const auto *pfg1_135 = buffer.data(pfg1 + 135);
    const auto *pfg1_140 = buffer.data(pfg1 + 140);

    const auto *pgd0_51 = buffer.data(pgd0 + 51);
    const auto *pgd0_53 = buffer.data(pgd0 + 53);
    const auto *pgd0_54 = buffer.data(pgd0 + 54);
    const auto *pgd0_57 = buffer.data(pgd0 + 57);
    const auto *pgd0_59 = buffer.data(pgd0 + 59);
    const auto *pgd0_60 = buffer.data(pgd0 + 60);
    const auto *pgd0_84 = buffer.data(pgd0 + 84);
    const auto *pgd0_93 = buffer.data(pgd0 + 93);
    const auto *pgd0_96 = buffer.data(pgd0 + 96);
    const auto *pgd0_97 = buffer.data(pgd0 + 97);
    const auto *pgd0_99 = buffer.data(pgd0 + 99);
    const auto *pgd0_101 = buffer.data(pgd0 + 101);

    const auto *pgd1_51 = buffer.data(pgd1 + 51);
    const auto *pgd1_53 = buffer.data(pgd1 + 53);
    const auto *pgd1_54 = buffer.data(pgd1 + 54);
    const auto *pgd1_57 = buffer.data(pgd1 + 57);
    const auto *pgd1_59 = buffer.data(pgd1 + 59);
    const auto *pgd1_60 = buffer.data(pgd1 + 60);
    const auto *pgd1_84 = buffer.data(pgd1 + 84);
    const auto *pgd1_93 = buffer.data(pgd1 + 93);
    const auto *pgd1_96 = buffer.data(pgd1 + 96);
    const auto *pgd1_97 = buffer.data(pgd1 + 97);
    const auto *pgd1_99 = buffer.data(pgd1 + 99);
    const auto *pgd1_101 = buffer.data(pgd1 + 101);

    const auto *pgf_83 = buffer.data(pgf + 83);
    const auto *pgf_85 = buffer.data(pgf + 85);
    const auto *pgf_86 = buffer.data(pgf + 86);
    const auto *pgf_88 = buffer.data(pgf + 88);
    const auto *pgf_89 = buffer.data(pgf + 89);
    const auto *pgf_90 = buffer.data(pgf + 90);
    const auto *pgf_91 = buffer.data(pgf + 91);
    const auto *pgf_92 = buffer.data(pgf + 92);
    const auto *pgf_93 = buffer.data(pgf + 93);
    const auto *pgf_95 = buffer.data(pgf + 95);
    const auto *pgf_96 = buffer.data(pgf + 96);
    const auto *pgf_98 = buffer.data(pgf + 98);
    const auto *pgf_99 = buffer.data(pgf + 99);
    const auto *pgf_100 = buffer.data(pgf + 100);
    const auto *pgf_102 = buffer.data(pgf + 102);
    const auto *pgf_103 = buffer.data(pgf + 103);
    const auto *pgf_105 = buffer.data(pgf + 105);
    const auto *pgf_106 = buffer.data(pgf + 106);
    const auto *pgf_109 = buffer.data(pgf + 109);
    const auto *pgf_110 = buffer.data(pgf + 110);
    const auto *pgf_112 = buffer.data(pgf + 112);
    const auto *pgf_113 = buffer.data(pgf + 113);
    const auto *pgf_115 = buffer.data(pgf + 115);
    const auto *pgf_116 = buffer.data(pgf + 116);
    const auto *pgf_119 = buffer.data(pgf + 119);
    const auto *pgf_120 = buffer.data(pgf + 120);
    const auto *pgf_122 = buffer.data(pgf + 122);
    const auto *pgf_123 = buffer.data(pgf + 123);
    const auto *pgf_125 = buffer.data(pgf + 125);
    const auto *pgf_126 = buffer.data(pgf + 126);
    const auto *pgf_129 = buffer.data(pgf + 129);
    const auto *pgf_130 = buffer.data(pgf + 130);
    const auto *pgf_132 = buffer.data(pgf + 132);
    const auto *pgf_133 = buffer.data(pgf + 133);
    const auto *pgf_135 = buffer.data(pgf + 135);
    const auto *pgf_136 = buffer.data(pgf + 136);
    const auto *pgf_139 = buffer.data(pgf + 139);
    const auto *pgf_140 = buffer.data(pgf + 140);
    const auto *pgf_141 = buffer.data(pgf + 141);
    const auto *pgf_142 = buffer.data(pgf + 142);
    const auto *pgf_143 = buffer.data(pgf + 143);
    const auto *pgf_145 = buffer.data(pgf + 145);
    const auto *pgf_146 = buffer.data(pgf + 146);
    const auto *pgf_149 = buffer.data(pgf + 149);
    const auto *pgf_150 = buffer.data(pgf + 150);
    const auto *pgf_151 = buffer.data(pgf + 151);
    const auto *pgf_156 = buffer.data(pgf + 156);
    const auto *pgf_157 = buffer.data(pgf + 157);
    const auto *pgf_158 = buffer.data(pgf + 158);
    const auto *pgf_159 = buffer.data(pgf + 159);
    const auto *pgf_160 = buffer.data(pgf + 160);
    const auto *pgf_161 = buffer.data(pgf + 161);
    const auto *pgf_163 = buffer.data(pgf + 163);
    const auto *pgf_165 = buffer.data(pgf + 165);
    const auto *pgf_166 = buffer.data(pgf + 166);
    const auto *pgf_167 = buffer.data(pgf + 167);
    const auto *pgf_168 = buffer.data(pgf + 168);
    const auto *pgf_169 = buffer.data(pgf + 169);

#pragma omp simd aligned(t_127, t_128, t_129, pb_y, pc_y, pc_z, pfg0_84, pff_43, pff_55, \
                         pfg1_84, pgf_83, pgf_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_8 * pff_43[k]
                   + f_4 * pc_z[k] * pgf_83[k];

        t_128[k] = f_0 * pff_55[k]
                   + f_4 * pc_y[k] * pgf_85[k];

        t_129[k] = pb_y[k] * pfg0_84[k]
                   - f_7 * pc_y[k] * pfg1_84[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pc_y, pc_z, pff_46, pff_56, pff_58, pgd0_51, \
                         pgd0_53, pgd1_51, pgd1_53, pgf_86, pgf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * pff_56[k]
                   + f_2 * pgd0_51[k]
                   - f_3 * pgd1_51[k]
                   + f_4 * pc_y[k] * pgf_86[k];

        t_131[k] = f_8 * pff_46[k]
                   + f_4 * pc_z[k] * pgf_86[k];

        t_132[k] = f_0 * pff_58[k]
                   + f_5 * pgd0_53[k]
                   - f_6 * pgd1_53[k]
                   + f_4 * pc_y[k] * pgf_88[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pc_x, pc_y, pc_z, sgf_90, pff_49, pff_59, \
                         pff_90, pgd0_53, pgd0_54, pgd1_53, pgd1_54, pgf_89, \
                         pgf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_0 * pff_59[k]
                   + f_4 * pc_y[k] * pgf_89[k];

        t_134[k] = f_8 * pff_49[k]
                   + f_2 * pgd0_53[k]
                   - f_3 * pgd1_53[k]
                   + f_4 * pc_z[k] * pgf_89[k];

        t_135[k] = f_0 * sgf_90[k]
                   + f_0 * pff_90[k]
                   + f_2 * pgd0_54[k]
                   - f_3 * pgd1_54[k]
                   + f_4 * pc_x[k] * pgf_90[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pc_y, pc_z, pff_50, pff_52, \
                         pgd0_54, pgd1_54, pgf_90, pgf_91, pgf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_4 * pc_y[k] * pgf_90[k];

        t_137[k] = f_9 * pff_50[k]
                   + f_4 * pc_z[k] * pgf_90[k];

        t_138[k] = f_5 * pgd0_54[k]
                   - f_6 * pgd1_54[k]
                   + f_4 * pc_y[k] * pgf_91[k];

        t_139[k] = f_4 * pc_y[k] * pgf_92[k];

        t_140[k] = f_9 * pff_52[k]
                   + f_5 * pgd0_54[k]
                   - f_6 * pgd1_54[k]
                   + f_4 * pc_z[k] * pgf_92[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, sgf_96, sgf_99, pff_53, \
                         pff_96, pff_99, pgf_93, pgf_95, pgf_96, \
                         pgf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * sgf_96[k]
                   + f_0 * pff_96[k]
                   + f_4 * pc_x[k] * pgf_96[k];

        t_142[k] = f_9 * pff_53[k]
                   + f_4 * pc_z[k] * pgf_93[k];

        t_143[k] = f_4 * pc_y[k] * pgf_95[k];

        t_144[k] = f_0 * sgf_99[k]
                   + f_0 * pff_99[k]
                   + f_4 * pc_x[k] * pgf_99[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pc_y, pc_z, pff_56, pff_59, \
                         pgd0_57, pgd0_59, pgd1_57, pgd1_59, pgf_96, pgf_98, \
                         pgf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_2 * pgd0_57[k]
                   - f_3 * pgd1_57[k]
                   + f_4 * pc_y[k] * pgf_96[k];

        t_146[k] = f_9 * pff_56[k]
                   + f_4 * pc_z[k] * pgf_96[k];

        t_147[k] = f_5 * pgd0_59[k]
                   - f_6 * pgd1_59[k]
                   + f_4 * pc_y[k] * pgf_98[k];

        t_148[k] = f_4 * pc_y[k] * pgf_99[k];

        t_149[k] = f_9 * pff_59[k]
                   + f_2 * pgd0_59[k]
                   - f_3 * pgd1_59[k]
                   + f_4 * pc_z[k] * pgf_99[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_x, pc_x, pc_y, pc_z, sgg0_150, \
                         sgg0_153, sgf_100, sgf_103, sgg1_150, sgg1_153, pff_60, \
                         pgf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_x[k] * sgg0_150[k]
                   + f_1 * sgf_100[k]
                   - f_7 * pc_x[k] * sgg1_150[k];

        t_151[k] = f_1 * pff_60[k]
                   + f_4 * pc_y[k] * pgf_100[k];

        t_152[k] = f_4 * pc_z[k] * pgf_100[k];

        t_153[k] = pa_x[k] * sgg0_153[k]
                   + f_8 * sgf_103[k]
                   - f_7 * pc_x[k] * sgg1_153[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, sgf_106, pff_62, \
                         pgd0_60, pgd1_60, pgf_102, pgf_103, pgf_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_1 * pff_62[k]
                   + f_4 * pc_y[k] * pgf_102[k];

        t_155[k] = f_5 * pgd0_60[k]
                   - f_6 * pgd1_60[k]
                   + f_4 * pc_z[k] * pgf_102[k];

        t_156[k] = f_0 * sgf_106[k]
                   + f_4 * pc_x[k] * pgf_106[k];

        t_157[k] = f_4 * pc_z[k] * pgf_103[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_x, pc_x, pc_y, pc_z, sgg0_160, \
                         sgf_109, sgg1_160, pff_65, pgf_105, pgf_106, \
                         pgf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_1 * pff_65[k]
                   + f_4 * pc_y[k] * pgf_105[k];

        t_159[k] = f_0 * sgf_109[k]
                   + f_4 * pc_x[k] * pgf_109[k];

        t_160[k] = pa_x[k] * sgg0_160[k]
                   - f_7 * pc_x[k] * sgg1_160[k];

        t_161[k] = f_4 * pc_z[k] * pgf_106[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pa_x, pc_x, pc_y, sgg0_162, sgg0_164, sgg1_162, \
                         sgg1_164, pff_69, pgf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_x[k] * sgg0_162[k]
                   - f_7 * pc_x[k] * sgg1_162[k];

        t_163[k] = f_1 * pff_69[k]
                   + f_4 * pc_y[k] * pgf_109[k];

        t_164[k] = pa_x[k] * sgg0_164[k]
                   - f_7 * pc_x[k] * sgg1_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pb_z, pc_y, pc_z, pfg0_90, pfg0_93, \
                         pff_60, pff_70, pfg1_90, pfg1_93, pgf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pb_z[k] * pfg0_90[k]
                   - f_7 * pc_z[k] * pfg1_90[k];

        t_166[k] = f_9 * pff_70[k]
                   + f_4 * pc_y[k] * pgf_110[k];

        t_167[k] = f_0 * pff_60[k]
                   + f_4 * pc_z[k] * pgf_110[k];

        t_168[k] = pb_z[k] * pfg0_93[k]
                   - f_7 * pc_z[k] * pfg1_93[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pa_x, pc_x, pc_y, sgg0_170, sgf_115, sgf_116, \
                         sgg1_170, pff_72, pgf_112, pgf_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_9 * pff_72[k]
                   + f_4 * pc_y[k] * pgf_112[k];

        t_170[k] = pa_x[k] * sgg0_170[k]
                   + f_8 * sgf_115[k]
                   - f_7 * pc_x[k] * sgg1_170[k];

        t_171[k] = f_0 * sgf_116[k]
                   + f_4 * pc_x[k] * pgf_116[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_x, pc_x, pc_y, pc_z, sgg0_175, \
                         sgf_119, sgg1_175, pff_63, pff_75, pgf_113, pgf_115, \
                         pgf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_0 * pff_63[k]
                   + f_4 * pc_z[k] * pgf_113[k];

        t_173[k] = f_9 * pff_75[k]
                   + f_4 * pc_y[k] * pgf_115[k];

        t_174[k] = f_0 * sgf_119[k]
                   + f_4 * pc_x[k] * pgf_119[k];

        t_175[k] = pa_x[k] * sgg0_175[k]
                   - f_7 * pc_x[k] * sgg1_175[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_x, pc_x, pc_y, pc_z, sgg0_177, \
                         sgg0_179, sgg1_177, sgg1_179, pff_66, pff_79, pgf_116, \
                         pgf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_0 * pff_66[k]
                   + f_4 * pc_z[k] * pgf_116[k];

        t_177[k] = pa_x[k] * sgg0_177[k]
                   - f_7 * pc_x[k] * sgg1_177[k];

        t_178[k] = f_9 * pff_79[k]
                   + f_4 * pc_y[k] * pgf_119[k];

        t_179[k] = pa_x[k] * sgg0_179[k]
                   - f_7 * pc_x[k] * sgg1_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pa_x, pc_x, pc_y, pc_z, sgg0_180, sgf_120, \
                         sgg1_180, pff_70, pff_80, pgf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_x[k] * sgg0_180[k]
                   + f_1 * sgf_120[k]
                   - f_7 * pc_x[k] * sgg1_180[k];

        t_181[k] = f_8 * pff_80[k]
                   + f_4 * pc_y[k] * pgf_120[k];

        t_182[k] = f_8 * pff_70[k]
                   + f_4 * pc_z[k] * pgf_120[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pa_x, pc_x, pc_y, sgg0_183, sgg0_185, sgf_123, \
                         sgf_125, sgg1_183, sgg1_185, pff_82, pgf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pa_x[k] * sgg0_183[k]
                   + f_8 * sgf_123[k]
                   - f_7 * pc_x[k] * sgg1_183[k];

        t_184[k] = f_8 * pff_82[k]
                   + f_4 * pc_y[k] * pgf_122[k];

        t_185[k] = pa_x[k] * sgg0_185[k]
                   + f_8 * sgf_125[k]
                   - f_7 * pc_x[k] * sgg1_185[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, pc_y, pc_z, sgf_126, sgf_129, \
                         pff_73, pff_85, pgf_123, pgf_125, pgf_126, \
                         pgf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_0 * sgf_126[k]
                   + f_4 * pc_x[k] * pgf_126[k];

        t_187[k] = f_8 * pff_73[k]
                   + f_4 * pc_z[k] * pgf_123[k];

        t_188[k] = f_8 * pff_85[k]
                   + f_4 * pc_y[k] * pgf_125[k];

        t_189[k] = f_0 * sgf_129[k]
                   + f_4 * pc_x[k] * pgf_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pa_x, pc_x, pc_y, pc_z, sgg0_190, \
                         sgg0_192, sgg1_190, sgg1_192, pff_76, pff_89, pgf_126, \
                         pgf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pa_x[k] * sgg0_190[k]
                   - f_7 * pc_x[k] * sgg1_190[k];

        t_191[k] = f_8 * pff_76[k]
                   + f_4 * pc_z[k] * pgf_126[k];

        t_192[k] = pa_x[k] * sgg0_192[k]
                   - f_7 * pc_x[k] * sgg1_192[k];

        t_193[k] = f_8 * pff_89[k]
                   + f_4 * pc_y[k] * pgf_129[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pa_x, pb_y, pc_x, pc_y, pc_z, sgg0_194, \
                         sgg1_194, pfg0_135, pff_80, pff_90, pfg1_135, \
                         pgf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_x[k] * sgg0_194[k]
                   - f_7 * pc_x[k] * sgg1_194[k];

        t_195[k] = pb_y[k] * pfg0_135[k]
                   - f_7 * pc_y[k] * pfg1_135[k];

        t_196[k] = f_0 * pff_90[k]
                   + f_4 * pc_y[k] * pgf_130[k];

        t_197[k] = f_9 * pff_80[k]
                   + f_4 * pc_z[k] * pgf_130[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_x, pb_y, pc_x, pc_y, sgg0_198, sgf_133, \
                         sgg1_198, pfg0_140, pff_92, pfg1_140, \
                         pgf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_x[k] * sgg0_198[k]
                   + f_8 * sgf_133[k]
                   - f_7 * pc_x[k] * sgg1_198[k];

        t_199[k] = f_0 * pff_92[k]
                   + f_4 * pc_y[k] * pgf_132[k];

        t_200[k] = pb_y[k] * pfg0_140[k]
                   - f_7 * pc_y[k] * pfg1_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pc_x, pc_y, pc_z, sgf_136, sgf_139, \
                         pff_83, pff_95, pgf_133, pgf_135, pgf_136, \
                         pgf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_0 * sgf_136[k]
                   + f_4 * pc_x[k] * pgf_136[k];

        t_202[k] = f_9 * pff_83[k]
                   + f_4 * pc_z[k] * pgf_133[k];

        t_203[k] = f_0 * pff_95[k]
                   + f_4 * pc_y[k] * pgf_135[k];

        t_204[k] = f_0 * sgf_139[k]
                   + f_4 * pc_x[k] * pgf_139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pa_x, pc_x, pc_y, pc_z, sgg0_205, \
                         sgg0_207, sgg1_205, sgg1_207, pff_86, pff_99, pgf_136, \
                         pgf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = pa_x[k] * sgg0_205[k]
                   - f_7 * pc_x[k] * sgg1_205[k];

        t_206[k] = f_9 * pff_86[k]
                   + f_4 * pc_z[k] * pgf_136[k];

        t_207[k] = pa_x[k] * sgg0_207[k]
                   - f_7 * pc_x[k] * sgg1_207[k];

        t_208[k] = f_0 * pff_99[k]
                   + f_4 * pc_y[k] * pgf_139[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pa_x, pc_x, pc_y, pc_z, sgg0_209, \
                         sgg0_210, sgf_140, sgg1_209, sgg1_210, pff_90, \
                         pgf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pa_x[k] * sgg0_209[k]
                   - f_7 * pc_x[k] * sgg1_209[k];

        t_210[k] = pa_x[k] * sgg0_210[k]
                   + f_1 * sgf_140[k]
                   - f_7 * pc_x[k] * sgg1_210[k];

        t_211[k] = f_4 * pc_y[k] * pgf_140[k];

        t_212[k] = f_1 * pff_90[k]
                   + f_4 * pc_z[k] * pgf_140[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_x, pc_x, pc_y, sgg0_215, sgf_145, \
                         sgf_146, sgg1_215, pgd0_84, pgd1_84, pgf_141, pgf_142, \
                         pgf_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_5 * pgd0_84[k]
                   - f_6 * pgd1_84[k]
                   + f_4 * pc_y[k] * pgf_141[k];

        t_214[k] = f_4 * pc_y[k] * pgf_142[k];

        t_215[k] = pa_x[k] * sgg0_215[k]
                   + f_8 * sgf_145[k]
                   - f_7 * pc_x[k] * sgg1_215[k];

        t_216[k] = f_0 * sgf_146[k]
                   + f_4 * pc_x[k] * pgf_146[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_x, pc_x, pc_y, pc_z, sgg0_220, \
                         sgf_149, sgg1_220, pff_93, pgf_143, pgf_145, \
                         pgf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_1 * pff_93[k]
                   + f_4 * pc_z[k] * pgf_143[k];

        t_218[k] = f_4 * pc_y[k] * pgf_145[k];

        t_219[k] = f_0 * sgf_149[k]
                   + f_4 * pc_x[k] * pgf_149[k];

        t_220[k] = pa_x[k] * sgg0_220[k]
                   - f_7 * pc_x[k] * sgg1_220[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_x, pc_x, pc_y, pc_z, sgg0_222, \
                         sgg0_224, sgg1_222, sgg1_224, pff_96, pgf_146, \
                         pgf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_1 * pff_96[k]
                   + f_4 * pc_z[k] * pgf_146[k];

        t_222[k] = pa_x[k] * sgg0_222[k]
                   - f_7 * pc_x[k] * sgg1_222[k];

        t_223[k] = f_4 * pc_y[k] * pgf_149[k];

        t_224[k] = pa_x[k] * sgg0_224[k]
                   - f_7 * pc_x[k] * sgg1_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_y, pc_y, pc_z, sgg0_0, sgg0_1, sgg0_3, \
                         sgf_0, sgf_1, sgg1_0, sgg1_1, sgg1_3, \
                         pgf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = pa_y[k] * sgg0_0[k]
                   - f_7 * pc_y[k] * sgg1_0[k];

        t_226[k] = pa_y[k] * sgg0_1[k]
                   + f_0 * sgf_0[k]
                   - f_7 * pc_y[k] * sgg1_1[k];

        t_227[k] = f_4 * pc_z[k] * pgf_150[k];

        t_228[k] = pa_y[k] * sgg0_3[k]
                   + f_8 * sgf_1[k]
                   - f_7 * pc_y[k] * sgg1_3[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_x, pc_y, pc_z, sgg0_5, sgg1_5, \
                         pff_106, pff_107, pgf_151, pgf_156, pgf_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_4 * pc_z[k] * pgf_151[k];

        t_230[k] = pa_y[k] * sgg0_5[k]
                   - f_7 * pc_y[k] * sgg1_5[k];

        t_231[k] = f_1 * pff_106[k]
                   + f_4 * pc_x[k] * pgf_156[k];

        t_232[k] = f_1 * pff_107[k]
                   + f_4 * pc_x[k] * pgf_157[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_y, pc_x, pc_y, pc_z, sgg0_10, sgf_6, \
                         sgg1_10, pff_108, pff_109, pgf_156, pgf_158, \
                         pgf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_1 * pff_108[k]
                   + f_4 * pc_x[k] * pgf_158[k];

        t_234[k] = f_1 * pff_109[k]
                   + f_4 * pc_x[k] * pgf_159[k];

        t_235[k] = pa_y[k] * sgg0_10[k]
                   + f_1 * sgf_6[k]
                   - f_7 * pc_y[k] * sgg1_10[k];

        t_236[k] = f_4 * pc_z[k] * pgf_156[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_y, pc_y, pc_z, sgg0_14, sgf_9, sgg1_14, \
                         pgd0_93, pgd1_93, pgf_157, pgf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_5 * pgd0_93[k]
                   - f_6 * pgd1_93[k]
                   + f_4 * pc_z[k] * pgf_157[k];

        t_238[k] = f_0 * sgf_9[k]
                   + f_4 * pc_y[k] * pgf_159[k];

        t_239[k] = pa_y[k] * sgg0_14[k]
                   - f_7 * pc_y[k] * sgg1_14[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pc_x, pc_z, pff_110, pff_111, pgd0_96, pgd0_97, \
                         pgd1_96, pgd1_97, pgf_160, pgf_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_9 * pff_110[k]
                   + f_2 * pgd0_96[k]
                   - f_3 * pgd1_96[k]
                   + f_4 * pc_x[k] * pgf_160[k];

        t_241[k] = f_9 * pff_111[k]
                   + f_12 * pgd0_97[k]
                   - f_13 * pgd1_97[k]
                   + f_4 * pc_x[k] * pgf_161[k];

        t_242[k] = f_4 * pc_z[k] * pgf_160[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pc_x, pc_z, pff_113, pff_115, pgd0_99, pgd0_101, \
                         pgd1_99, pgd1_101, pgf_161, pgf_163, pgf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_9 * pff_113[k]
                   + f_5 * pgd0_99[k]
                   - f_6 * pgd1_99[k]
                   + f_4 * pc_x[k] * pgf_163[k];

        t_244[k] = f_4 * pc_z[k] * pgf_161[k];

        t_245[k] = f_9 * pff_115[k]
                   + f_5 * pgd0_101[k]
                   - f_6 * pgd1_101[k]
                   + f_4 * pc_x[k] * pgf_165[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, pff_116, pff_117, pff_118, pff_119, \
                         pgf_166, pgf_167, pgf_168, pgf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_9 * pff_116[k]
                   + f_4 * pc_x[k] * pgf_166[k];

        t_247[k] = f_9 * pff_117[k]
                   + f_4 * pc_x[k] * pgf_167[k];

        t_248[k] = f_9 * pff_118[k]
                   + f_4 * pc_x[k] * pgf_168[k];

        t_249[k] = f_9 * pff_119[k]
                   + f_4 * pc_x[k] * pgf_169[k];
    }
}

static auto
compute_prim_pgg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgg0, const size_t sgf,
                                                          const size_t sgg1, const size_t pdg0,
                                                          const size_t pdg1, const size_t pfg0,
                                                          const size_t pff, const size_t pfg1,
                                                          const size_t pgd0, const size_t pgd1,
                                                          const size_t pgf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);
    const auto f_14 = 1.0 / p;
    const auto f_15 = gamma / (p * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgg0_30 = buffer.data(sgg0 + 30);
    const auto *sgg0_35 = buffer.data(sgg0 + 35);
    const auto *sgg0_42 = buffer.data(sgg0 + 42);
    const auto *sgg0_44 = buffer.data(sgg0 + 44);
    const auto *sgg0_75 = buffer.data(sgg0 + 75);
    const auto *sgg0_76 = buffer.data(sgg0 + 76);
    const auto *sgg0_78 = buffer.data(sgg0 + 78);
    const auto *sgg0_80 = buffer.data(sgg0 + 80);
    const auto *sgg0_85 = buffer.data(sgg0 + 85);
    const auto *sgg0_87 = buffer.data(sgg0 + 87);
    const auto *sgg0_89 = buffer.data(sgg0 + 89);
    const auto *sgg0_135 = buffer.data(sgg0 + 135);
    const auto *sgg0_136 = buffer.data(sgg0 + 136);
    const auto *sgg0_138 = buffer.data(sgg0 + 138);
    const auto *sgg0_140 = buffer.data(sgg0 + 140);

    const auto *sgf_19 = buffer.data(sgf + 19);
    const auto *sgf_28 = buffer.data(sgf + 28);
    const auto *sgf_29 = buffer.data(sgf + 29);
    const auto *sgf_39 = buffer.data(sgf + 39);
    const auto *sgf_49 = buffer.data(sgf + 49);
    const auto *sgf_50 = buffer.data(sgf + 50);
    const auto *sgf_51 = buffer.data(sgf + 51);
    const auto *sgf_56 = buffer.data(sgf + 56);
    const auto *sgf_58 = buffer.data(sgf + 58);
    const auto *sgf_59 = buffer.data(sgf + 59);
    const auto *sgf_90 = buffer.data(sgf + 90);
    const auto *sgf_91 = buffer.data(sgf + 91);

    const auto *sgg1_30 = buffer.data(sgg1 + 30);
    const auto *sgg1_35 = buffer.data(sgg1 + 35);
    const auto *sgg1_42 = buffer.data(sgg1 + 42);
    const auto *sgg1_44 = buffer.data(sgg1 + 44);
    const auto *sgg1_75 = buffer.data(sgg1 + 75);
    const auto *sgg1_76 = buffer.data(sgg1 + 76);
    const auto *sgg1_78 = buffer.data(sgg1 + 78);
    const auto *sgg1_80 = buffer.data(sgg1 + 80);
    const auto *sgg1_85 = buffer.data(sgg1 + 85);
    const auto *sgg1_87 = buffer.data(sgg1 + 87);
    const auto *sgg1_89 = buffer.data(sgg1 + 89);
    const auto *sgg1_135 = buffer.data(sgg1 + 135);
    const auto *sgg1_136 = buffer.data(sgg1 + 136);
    const auto *sgg1_138 = buffer.data(sgg1 + 138);
    const auto *sgg1_140 = buffer.data(sgg1 + 140);

    const auto *pdg0_115 = buffer.data(pdg0 + 115);
    const auto *pdg0_145 = buffer.data(pdg0 + 145);

    const auto *pdg1_115 = buffer.data(pdg1 + 115);
    const auto *pdg1_145 = buffer.data(pdg1 + 145);

    const auto *pfg0_151 = buffer.data(pfg0 + 151);
    const auto *pfg0_153 = buffer.data(pfg0 + 153);
    const auto *pfg0_160 = buffer.data(pfg0 + 160);
    const auto *pfg0_166 = buffer.data(pfg0 + 166);
    const auto *pfg0_168 = buffer.data(pfg0 + 168);
    const auto *pfg0_175 = buffer.data(pfg0 + 175);
    const auto *pfg0_195 = buffer.data(pfg0 + 195);
    const auto *pfg0_196 = buffer.data(pfg0 + 196);
    const auto *pfg0_198 = buffer.data(pfg0 + 198);
    const auto *pfg0_205 = buffer.data(pfg0 + 205);
    const auto *pfg0_240 = buffer.data(pfg0 + 240);
    const auto *pfg0_241 = buffer.data(pfg0 + 241);
    const auto *pfg0_243 = buffer.data(pfg0 + 243);
    const auto *pfg0_245 = buffer.data(pfg0 + 245);
    const auto *pfg0_250 = buffer.data(pfg0 + 250);
    const auto *pfg0_252 = buffer.data(pfg0 + 252);
    const auto *pfg0_253 = buffer.data(pfg0 + 253);
    const auto *pfg0_254 = buffer.data(pfg0 + 254);
    const auto *pfg0_260 = buffer.data(pfg0 + 260);
    const auto *pfg0_265 = buffer.data(pfg0 + 265);
    const auto *pfg0_267 = buffer.data(pfg0 + 267);
    const auto *pfg0_268 = buffer.data(pfg0 + 268);
    const auto *pfg0_269 = buffer.data(pfg0 + 269);
    const auto *pfg0_271 = buffer.data(pfg0 + 271);
    const auto *pfg0_273 = buffer.data(pfg0 + 273);
    const auto *pfg0_280 = buffer.data(pfg0 + 280);
    const auto *pfg0_282 = buffer.data(pfg0 + 282);
    const auto *pfg0_283 = buffer.data(pfg0 + 283);
    const auto *pfg0_284 = buffer.data(pfg0 + 284);

    const auto *pff_100 = buffer.data(pff + 100);
    const auto *pff_101 = buffer.data(pff + 101);
    const auto *pff_106 = buffer.data(pff + 106);
    const auto *pff_109 = buffer.data(pff + 109);
    const auto *pff_110 = buffer.data(pff + 110);
    const auto *pff_111 = buffer.data(pff + 111);
    const auto *pff_116 = buffer.data(pff + 116);
    const auto *pff_117 = buffer.data(pff + 117);
    const auto *pff_119 = buffer.data(pff + 119);
    const auto *pff_120 = buffer.data(pff + 120);
    const auto *pff_121 = buffer.data(pff + 121);
    const auto *pff_126 = buffer.data(pff + 126);
    const auto *pff_127 = buffer.data(pff + 127);
    const auto *pff_128 = buffer.data(pff + 128);
    const auto *pff_129 = buffer.data(pff + 129);
    const auto *pff_130 = buffer.data(pff + 130);
    const auto *pff_131 = buffer.data(pff + 131);
    const auto *pff_133 = buffer.data(pff + 133);
    const auto *pff_135 = buffer.data(pff + 135);
    const auto *pff_136 = buffer.data(pff + 136);
    const auto *pff_137 = buffer.data(pff + 137);
    const auto *pff_138 = buffer.data(pff + 138);
    const auto *pff_139 = buffer.data(pff + 139);
    const auto *pff_140 = buffer.data(pff + 140);
    const auto *pff_141 = buffer.data(pff + 141);
    const auto *pff_145 = buffer.data(pff + 145);
    const auto *pff_146 = buffer.data(pff + 146);
    const auto *pff_147 = buffer.data(pff + 147);
    const auto *pff_148 = buffer.data(pff + 148);
    const auto *pff_149 = buffer.data(pff + 149);
    const auto *pff_150 = buffer.data(pff + 150);
    const auto *pff_151 = buffer.data(pff + 151);
    const auto *pff_156 = buffer.data(pff + 156);
    const auto *pff_157 = buffer.data(pff + 157);
    const auto *pff_158 = buffer.data(pff + 158);
    const auto *pff_159 = buffer.data(pff + 159);
    const auto *pff_160 = buffer.data(pff + 160);
    const auto *pff_161 = buffer.data(pff + 161);
    const auto *pff_163 = buffer.data(pff + 163);
    const auto *pff_165 = buffer.data(pff + 165);
    const auto *pff_166 = buffer.data(pff + 166);
    const auto *pff_167 = buffer.data(pff + 167);
    const auto *pff_168 = buffer.data(pff + 168);
    const auto *pff_169 = buffer.data(pff + 169);
    const auto *pff_175 = buffer.data(pff + 175);
    const auto *pff_176 = buffer.data(pff + 176);
    const auto *pff_177 = buffer.data(pff + 177);
    const auto *pff_178 = buffer.data(pff + 178);
    const auto *pff_179 = buffer.data(pff + 179);
    const auto *pff_180 = buffer.data(pff + 180);
    const auto *pff_181 = buffer.data(pff + 181);
    const auto *pff_183 = buffer.data(pff + 183);
    const auto *pff_185 = buffer.data(pff + 185);
    const auto *pff_186 = buffer.data(pff + 186);
    const auto *pff_187 = buffer.data(pff + 187);
    const auto *pff_188 = buffer.data(pff + 188);
    const auto *pff_189 = buffer.data(pff + 189);
    const auto *pff_196 = buffer.data(pff + 196);
    const auto *pff_197 = buffer.data(pff + 197);

    const auto *pfg1_151 = buffer.data(pfg1 + 151);
    const auto *pfg1_153 = buffer.data(pfg1 + 153);
    const auto *pfg1_160 = buffer.data(pfg1 + 160);
    const auto *pfg1_166 = buffer.data(pfg1 + 166);
    const auto *pfg1_168 = buffer.data(pfg1 + 168);
    const auto *pfg1_175 = buffer.data(pfg1 + 175);
    const auto *pfg1_195 = buffer.data(pfg1 + 195);
    const auto *pfg1_196 = buffer.data(pfg1 + 196);
    const auto *pfg1_198 = buffer.data(pfg1 + 198);
    const auto *pfg1_205 = buffer.data(pfg1 + 205);
    const auto *pfg1_240 = buffer.data(pfg1 + 240);
    const auto *pfg1_241 = buffer.data(pfg1 + 241);
    const auto *pfg1_243 = buffer.data(pfg1 + 243);
    const auto *pfg1_245 = buffer.data(pfg1 + 245);
    const auto *pfg1_250 = buffer.data(pfg1 + 250);
    const auto *pfg1_252 = buffer.data(pfg1 + 252);
    const auto *pfg1_253 = buffer.data(pfg1 + 253);
    const auto *pfg1_254 = buffer.data(pfg1 + 254);
    const auto *pfg1_260 = buffer.data(pfg1 + 260);
    const auto *pfg1_265 = buffer.data(pfg1 + 265);
    const auto *pfg1_267 = buffer.data(pfg1 + 267);
    const auto *pfg1_268 = buffer.data(pfg1 + 268);
    const auto *pfg1_269 = buffer.data(pfg1 + 269);
    const auto *pfg1_271 = buffer.data(pfg1 + 271);
    const auto *pfg1_273 = buffer.data(pfg1 + 273);
    const auto *pfg1_280 = buffer.data(pfg1 + 280);
    const auto *pfg1_282 = buffer.data(pfg1 + 282);
    const auto *pfg1_283 = buffer.data(pfg1 + 283);
    const auto *pfg1_284 = buffer.data(pfg1 + 284);

    const auto *pgd0_99 = buffer.data(pgd0 + 99);
    const auto *pgd0_101 = buffer.data(pgd0 + 101);
    const auto *pgd0_108 = buffer.data(pgd0 + 108);
    const auto *pgd0_109 = buffer.data(pgd0 + 109);
    const auto *pgd0_111 = buffer.data(pgd0 + 111);
    const auto *pgd0_113 = buffer.data(pgd0 + 113);
    const auto *pgd0_114 = buffer.data(pgd0 + 114);
    const auto *pgd0_117 = buffer.data(pgd0 + 117);
    const auto *pgd0_119 = buffer.data(pgd0 + 119);
    const auto *pgd0_138 = buffer.data(pgd0 + 138);
    const auto *pgd0_143 = buffer.data(pgd0 + 143);

    const auto *pgd1_99 = buffer.data(pgd1 + 99);
    const auto *pgd1_101 = buffer.data(pgd1 + 101);
    const auto *pgd1_108 = buffer.data(pgd1 + 108);
    const auto *pgd1_109 = buffer.data(pgd1 + 109);
    const auto *pgd1_111 = buffer.data(pgd1 + 111);
    const auto *pgd1_113 = buffer.data(pgd1 + 113);
    const auto *pgd1_114 = buffer.data(pgd1 + 114);
    const auto *pgd1_117 = buffer.data(pgd1 + 117);
    const auto *pgd1_119 = buffer.data(pgd1 + 119);
    const auto *pgd1_138 = buffer.data(pgd1 + 138);
    const auto *pgd1_143 = buffer.data(pgd1 + 143);

    const auto *pgf_166 = buffer.data(pgf + 166);
    const auto *pgf_167 = buffer.data(pgf + 167);
    const auto *pgf_169 = buffer.data(pgf + 169);
    const auto *pgf_170 = buffer.data(pgf + 170);
    const auto *pgf_171 = buffer.data(pgf + 171);
    const auto *pgf_176 = buffer.data(pgf + 176);
    const auto *pgf_177 = buffer.data(pgf + 177);
    const auto *pgf_178 = buffer.data(pgf + 178);
    const auto *pgf_179 = buffer.data(pgf + 179);
    const auto *pgf_180 = buffer.data(pgf + 180);
    const auto *pgf_181 = buffer.data(pgf + 181);
    const auto *pgf_183 = buffer.data(pgf + 183);
    const auto *pgf_185 = buffer.data(pgf + 185);
    const auto *pgf_186 = buffer.data(pgf + 186);
    const auto *pgf_187 = buffer.data(pgf + 187);
    const auto *pgf_188 = buffer.data(pgf + 188);
    const auto *pgf_189 = buffer.data(pgf + 189);
    const auto *pgf_190 = buffer.data(pgf + 190);
    const auto *pgf_191 = buffer.data(pgf + 191);
    const auto *pgf_195 = buffer.data(pgf + 195);
    const auto *pgf_196 = buffer.data(pgf + 196);
    const auto *pgf_197 = buffer.data(pgf + 197);
    const auto *pgf_198 = buffer.data(pgf + 198);
    const auto *pgf_199 = buffer.data(pgf + 199);
    const auto *pgf_200 = buffer.data(pgf + 200);
    const auto *pgf_201 = buffer.data(pgf + 201);
    const auto *pgf_206 = buffer.data(pgf + 206);
    const auto *pgf_207 = buffer.data(pgf + 207);
    const auto *pgf_208 = buffer.data(pgf + 208);
    const auto *pgf_209 = buffer.data(pgf + 209);
    const auto *pgf_210 = buffer.data(pgf + 210);
    const auto *pgf_211 = buffer.data(pgf + 211);
    const auto *pgf_216 = buffer.data(pgf + 216);
    const auto *pgf_217 = buffer.data(pgf + 217);
    const auto *pgf_218 = buffer.data(pgf + 218);
    const auto *pgf_219 = buffer.data(pgf + 219);
    const auto *pgf_220 = buffer.data(pgf + 220);
    const auto *pgf_221 = buffer.data(pgf + 221);
    const auto *pgf_226 = buffer.data(pgf + 226);
    const auto *pgf_227 = buffer.data(pgf + 227);
    const auto *pgf_228 = buffer.data(pgf + 228);
    const auto *pgf_229 = buffer.data(pgf + 229);
    const auto *pgf_230 = buffer.data(pgf + 230);
    const auto *pgf_231 = buffer.data(pgf + 231);
    const auto *pgf_235 = buffer.data(pgf + 235);
    const auto *pgf_236 = buffer.data(pgf + 236);
    const auto *pgf_237 = buffer.data(pgf + 237);
    const auto *pgf_238 = buffer.data(pgf + 238);
    const auto *pgf_239 = buffer.data(pgf + 239);
    const auto *pgf_240 = buffer.data(pgf + 240);
    const auto *pgf_241 = buffer.data(pgf + 241);
    const auto *pgf_246 = buffer.data(pgf + 246);
    const auto *pgf_247 = buffer.data(pgf + 247);

#pragma omp simd aligned(t_250, t_251, t_252, pb_x, pc_x, pc_z, pdg0_115, pdg1_115, pfg0_175, \
                         pfg1_175, pgd0_99, pgd1_99, pgf_166, pgf_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_14 * pdg0_115[k]
                   - f_15 * pdg1_115[k]
                   + pb_x[k] * pfg0_175[k]
                   - f_7 * pc_x[k] * pfg1_175[k];

        t_251[k] = f_4 * pc_z[k] * pgf_166[k];

        t_252[k] = f_5 * pgd0_99[k]
                   - f_6 * pgd1_99[k]
                   + f_4 * pc_z[k] * pgf_167[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, pa_y, pc_y, pc_z, sgg0_30, sgf_19, sgg1_30, \
                         pff_109, pgd0_101, pgd1_101, pgf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_0 * sgf_19[k]
                   + f_0 * pff_109[k]
                   + f_4 * pc_y[k] * pgf_169[k];

        t_254[k] = f_2 * pgd0_101[k]
                   - f_3 * pgd1_101[k]
                   + f_4 * pc_z[k] * pgf_169[k];

        t_255[k] = pa_y[k] * sgg0_30[k]
                   - f_7 * pc_y[k] * sgg1_30[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pb_z, pc_z, pfg0_151, pfg0_153, pff_100, \
                         pff_101, pfg1_151, pfg1_153, pgf_170, \
                         pgf_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = pb_z[k] * pfg0_151[k]
                   - f_7 * pc_z[k] * pfg1_151[k];

        t_257[k] = f_0 * pff_100[k]
                   + f_4 * pc_z[k] * pgf_170[k];

        t_258[k] = pb_z[k] * pfg0_153[k]
                   - f_7 * pc_z[k] * pfg1_153[k];

        t_259[k] = f_0 * pff_101[k]
                   + f_4 * pc_z[k] * pgf_171[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, sgg0_35, sgg1_35, \
                         pff_126, pff_127, pff_128, pgf_176, pgf_177, \
                         pgf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * sgg0_35[k]
                   - f_7 * pc_y[k] * sgg1_35[k];

        t_261[k] = f_9 * pff_126[k]
                   + f_4 * pc_x[k] * pgf_176[k];

        t_262[k] = f_9 * pff_127[k]
                   + f_4 * pc_x[k] * pgf_177[k];

        t_263[k] = f_9 * pff_128[k]
                   + f_4 * pc_x[k] * pgf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pb_z, pc_x, pc_z, pfg0_160, pff_106, pff_129, \
                         pfg1_160, pgf_176, pgf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_9 * pff_129[k]
                   + f_4 * pc_x[k] * pgf_179[k];

        t_265[k] = pb_z[k] * pfg0_160[k]
                   - f_7 * pc_z[k] * pfg1_160[k];

        t_266[k] = f_0 * pff_106[k]
                   + f_4 * pc_z[k] * pgf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_y, pc_y, sgg0_42, sgg0_44, sgf_28, sgf_29, \
                         sgg1_42, sgg1_44, pgf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pa_y[k] * sgg0_42[k]
                   + f_8 * sgf_28[k]
                   - f_7 * pc_y[k] * sgg1_42[k];

        t_268[k] = f_0 * sgf_29[k]
                   + f_4 * pc_y[k] * pgf_179[k];

        t_269[k] = pa_y[k] * sgg0_44[k]
                   - f_7 * pc_y[k] * sgg1_44[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pc_x, pc_z, pff_130, pff_131, pgd0_108, \
                         pgd0_109, pgd1_108, pgd1_109, pgf_180, \
                         pgf_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_8 * pff_130[k]
                   + f_2 * pgd0_108[k]
                   - f_3 * pgd1_108[k]
                   + f_4 * pc_x[k] * pgf_180[k];

        t_271[k] = f_8 * pff_131[k]
                   + f_12 * pgd0_109[k]
                   - f_13 * pgd1_109[k]
                   + f_4 * pc_x[k] * pgf_181[k];

        t_272[k] = f_4 * pc_z[k] * pgf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_x, pc_z, pff_133, pff_135, pgd0_111, \
                         pgd0_113, pgd1_111, pgd1_113, pgf_181, pgf_183, \
                         pgf_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_8 * pff_133[k]
                   + f_5 * pgd0_111[k]
                   - f_6 * pgd1_111[k]
                   + f_4 * pc_x[k] * pgf_183[k];

        t_274[k] = f_4 * pc_z[k] * pgf_181[k];

        t_275[k] = f_8 * pff_135[k]
                   + f_5 * pgd0_113[k]
                   - f_6 * pgd1_113[k]
                   + f_4 * pc_x[k] * pgf_185[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, pff_136, pff_137, pff_138, pff_139, \
                         pgf_186, pgf_187, pgf_188, pgf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_8 * pff_136[k]
                   + f_4 * pc_x[k] * pgf_186[k];

        t_277[k] = f_8 * pff_137[k]
                   + f_4 * pc_x[k] * pgf_187[k];

        t_278[k] = f_8 * pff_138[k]
                   + f_4 * pc_x[k] * pgf_188[k];

        t_279[k] = f_8 * pff_139[k]
                   + f_4 * pc_x[k] * pgf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_x, pc_x, pc_z, pdg0_145, pdg1_145, pfg0_205, \
                         pfg1_205, pgd0_111, pgd1_111, pgf_186, \
                         pgf_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_10 * pdg0_145[k]
                   - f_11 * pdg1_145[k]
                   + pb_x[k] * pfg0_205[k]
                   - f_7 * pc_x[k] * pfg1_205[k];

        t_281[k] = f_4 * pc_z[k] * pgf_186[k];

        t_282[k] = f_5 * pgd0_111[k]
                   - f_6 * pgd1_111[k]
                   + f_4 * pc_z[k] * pgf_187[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_x, pc_y, pc_z, sgf_39, pff_119, pff_140, \
                         pgd0_113, pgd0_114, pgd1_113, pgd1_114, pgf_189, \
                         pgf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_0 * sgf_39[k]
                   + f_8 * pff_119[k]
                   + f_4 * pc_y[k] * pgf_189[k];

        t_284[k] = f_2 * pgd0_113[k]
                   - f_3 * pgd1_113[k]
                   + f_4 * pc_z[k] * pgf_189[k];

        t_285[k] = f_8 * pff_140[k]
                   + f_2 * pgd0_114[k]
                   - f_3 * pgd1_114[k]
                   + f_4 * pc_x[k] * pgf_190[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_z, pc_z, pfg0_166, pfg0_168, pff_110, \
                         pff_111, pfg1_166, pfg1_168, pgf_190, \
                         pgf_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = pb_z[k] * pfg0_166[k]
                   - f_7 * pc_z[k] * pfg1_166[k];

        t_287[k] = f_0 * pff_110[k]
                   + f_4 * pc_z[k] * pgf_190[k];

        t_288[k] = pb_z[k] * pfg0_168[k]
                   - f_7 * pc_z[k] * pfg1_168[k];

        t_289[k] = f_0 * pff_111[k]
                   + f_4 * pc_z[k] * pgf_191[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pc_x, pff_145, pff_146, pff_147, pff_148, \
                         pgd0_119, pgd1_119, pgf_195, pgf_196, pgf_197, \
                         pgf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_8 * pff_145[k]
                   + f_5 * pgd0_119[k]
                   - f_6 * pgd1_119[k]
                   + f_4 * pc_x[k] * pgf_195[k];

        t_291[k] = f_8 * pff_146[k]
                   + f_4 * pc_x[k] * pgf_196[k];

        t_292[k] = f_8 * pff_147[k]
                   + f_4 * pc_x[k] * pgf_197[k];

        t_293[k] = f_8 * pff_148[k]
                   + f_4 * pc_x[k] * pgf_198[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pb_z, pc_x, pc_z, pfg0_175, pff_116, pff_149, \
                         pfg1_175, pgf_196, pgf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_8 * pff_149[k]
                   + f_4 * pc_x[k] * pgf_199[k];

        t_295[k] = pb_z[k] * pfg0_175[k]
                   - f_7 * pc_z[k] * pfg1_175[k];

        t_296[k] = f_0 * pff_116[k]
                   + f_4 * pc_z[k] * pgf_196[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pc_y, pc_z, sgf_49, pff_117, pff_119, pff_129, \
                         pgd0_117, pgd0_119, pgd1_117, pgd1_119, pgf_197, \
                         pgf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_0 * pff_117[k]
                   + f_5 * pgd0_117[k]
                   - f_6 * pgd1_117[k]
                   + f_4 * pc_z[k] * pgf_197[k];

        t_298[k] = f_0 * sgf_49[k]
                   + f_0 * pff_129[k]
                   + f_4 * pc_y[k] * pgf_199[k];

        t_299[k] = f_0 * pff_119[k]
                   + f_2 * pgd0_119[k]
                   - f_3 * pgd1_119[k]
                   + f_4 * pc_z[k] * pgf_199[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pa_y, pc_y, pc_z, sgg0_75, sgg0_76, sgf_50, \
                         sgg1_75, sgg1_76, pff_120, pgf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = pa_y[k] * sgg0_75[k]
                   - f_7 * pc_y[k] * sgg1_75[k];

        t_301[k] = pa_y[k] * sgg0_76[k]
                   + f_0 * sgf_50[k]
                   - f_7 * pc_y[k] * sgg1_76[k];

        t_302[k] = f_8 * pff_120[k]
                   + f_4 * pc_z[k] * pgf_200[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, pa_y, pc_y, pc_z, sgg0_78, sgg0_80, sgf_51, \
                         sgg1_78, sgg1_80, pff_121, pgf_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pa_y[k] * sgg0_78[k]
                   + f_8 * sgf_51[k]
                   - f_7 * pc_y[k] * sgg1_78[k];

        t_304[k] = f_8 * pff_121[k]
                   + f_4 * pc_z[k] * pgf_201[k];

        t_305[k] = pa_y[k] * sgg0_80[k]
                   - f_7 * pc_y[k] * sgg1_80[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pc_x, pff_156, pff_157, pff_158, pff_159, \
                         pgf_206, pgf_207, pgf_208, pgf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_8 * pff_156[k]
                   + f_4 * pc_x[k] * pgf_206[k];

        t_307[k] = f_8 * pff_157[k]
                   + f_4 * pc_x[k] * pgf_207[k];

        t_308[k] = f_8 * pff_158[k]
                   + f_4 * pc_x[k] * pgf_208[k];

        t_309[k] = f_8 * pff_159[k]
                   + f_4 * pc_x[k] * pgf_209[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, pa_y, pc_y, pc_z, sgg0_85, sgg0_87, sgf_56, \
                         sgf_58, sgg1_85, sgg1_87, pff_126, pgf_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = pa_y[k] * sgg0_85[k]
                   + f_1 * sgf_56[k]
                   - f_7 * pc_y[k] * sgg1_85[k];

        t_311[k] = f_8 * pff_126[k]
                   + f_4 * pc_z[k] * pgf_206[k];

        t_312[k] = pa_y[k] * sgg0_87[k]
                   + f_8 * sgf_58[k]
                   - f_7 * pc_y[k] * sgg1_87[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_y, pb_x, pc_x, pc_y, sgg0_89, sgf_59, \
                         sgg1_89, pfg0_240, pff_160, pfg1_240, \
                         pgf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_0 * sgf_59[k]
                   + f_4 * pc_y[k] * pgf_209[k];

        t_314[k] = pa_y[k] * sgg0_89[k]
                   - f_7 * pc_y[k] * sgg1_89[k];

        t_315[k] = pb_x[k] * pfg0_240[k]
                   + f_1 * pff_160[k]
                   - f_7 * pc_x[k] * pfg1_240[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pb_x, pc_x, pc_z, pfg0_241, pfg0_243, \
                         pff_161, pff_163, pfg1_241, pfg1_243, pgf_210, \
                         pgf_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pb_x[k] * pfg0_241[k]
                   + f_9 * pff_161[k]
                   - f_7 * pc_x[k] * pfg1_241[k];

        t_317[k] = f_4 * pc_z[k] * pgf_210[k];

        t_318[k] = pb_x[k] * pfg0_243[k]
                   + f_8 * pff_163[k]
                   - f_7 * pc_x[k] * pfg1_243[k];

        t_319[k] = f_4 * pc_z[k] * pgf_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pb_x, pc_x, pfg0_245, pff_165, pff_166, \
                         pff_167, pff_168, pfg1_245, pgf_216, pgf_217, \
                         pgf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = pb_x[k] * pfg0_245[k]
                   + f_8 * pff_165[k]
                   - f_7 * pc_x[k] * pfg1_245[k];

        t_321[k] = f_0 * pff_166[k]
                   + f_4 * pc_x[k] * pgf_216[k];

        t_322[k] = f_0 * pff_167[k]
                   + f_4 * pc_x[k] * pgf_217[k];

        t_323[k] = f_0 * pff_168[k]
                   + f_4 * pc_x[k] * pgf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pb_x, pc_x, pc_z, pfg0_250, pfg0_252, \
                         pff_169, pfg1_250, pfg1_252, pgf_216, \
                         pgf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_0 * pff_169[k]
                   + f_4 * pc_x[k] * pgf_219[k];

        t_325[k] = pb_x[k] * pfg0_250[k]
                   - f_7 * pc_x[k] * pfg1_250[k];

        t_326[k] = f_4 * pc_z[k] * pgf_216[k];

        t_327[k] = pb_x[k] * pfg0_252[k]
                   - f_7 * pc_x[k] * pfg1_252[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pb_x, pb_z, pc_x, pc_z, pfg0_195, \
                         pfg0_196, pfg0_253, pfg0_254, pfg1_195, pfg1_196, pfg1_253, \
                         pfg1_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = pb_x[k] * pfg0_253[k]
                   - f_7 * pc_x[k] * pfg1_253[k];

        t_329[k] = pb_x[k] * pfg0_254[k]
                   - f_7 * pc_x[k] * pfg1_254[k];

        t_330[k] = pb_z[k] * pfg0_195[k]
                   - f_7 * pc_z[k] * pfg1_195[k];

        t_331[k] = pb_z[k] * pfg0_196[k]
                   - f_7 * pc_z[k] * pfg1_196[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pb_z, pc_z, pfg0_198, pff_130, pff_131, \
                         pfg1_198, pgf_220, pgf_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_0 * pff_130[k]
                   + f_4 * pc_z[k] * pgf_220[k];

        t_333[k] = pb_z[k] * pfg0_198[k]
                   - f_7 * pc_z[k] * pfg1_198[k];

        t_334[k] = f_0 * pff_131[k]
                   + f_4 * pc_z[k] * pgf_221[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pb_x, pc_x, pfg0_260, pff_175, pff_176, \
                         pff_177, pff_178, pfg1_260, pgf_226, pgf_227, \
                         pgf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pb_x[k] * pfg0_260[k]
                   + f_8 * pff_175[k]
                   - f_7 * pc_x[k] * pfg1_260[k];

        t_336[k] = f_0 * pff_176[k]
                   + f_4 * pc_x[k] * pgf_226[k];

        t_337[k] = f_0 * pff_177[k]
                   + f_4 * pc_x[k] * pgf_227[k];

        t_338[k] = f_0 * pff_178[k]
                   + f_4 * pc_x[k] * pgf_228[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, pb_x, pc_x, pc_z, pfg0_265, pfg0_267, \
                         pff_136, pff_179, pfg1_265, pfg1_267, pgf_226, \
                         pgf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_0 * pff_179[k]
                   + f_4 * pc_x[k] * pgf_229[k];

        t_340[k] = pb_x[k] * pfg0_265[k]
                   - f_7 * pc_x[k] * pfg1_265[k];

        t_341[k] = f_0 * pff_136[k]
                   + f_4 * pc_z[k] * pgf_226[k];

        t_342[k] = pb_x[k] * pfg0_267[k]
                   - f_7 * pc_x[k] * pfg1_267[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pb_x, pc_x, pfg0_268, pfg0_269, pff_180, \
                         pfg1_268, pfg1_269, pgd0_138, pgd1_138, \
                         pgf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = pb_x[k] * pfg0_268[k]
                   - f_7 * pc_x[k] * pfg1_268[k];

        t_344[k] = pb_x[k] * pfg0_269[k]
                   - f_7 * pc_x[k] * pfg1_269[k];

        t_345[k] = f_0 * pff_180[k]
                   + f_2 * pgd0_138[k]
                   - f_3 * pgd1_138[k]
                   + f_4 * pc_x[k] * pgf_230[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, pb_x, pc_x, pc_z, pfg0_271, pfg0_273, pff_140, \
                         pff_181, pff_183, pfg1_271, pfg1_273, \
                         pgf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = pb_x[k] * pfg0_271[k]
                   + f_9 * pff_181[k]
                   - f_7 * pc_x[k] * pfg1_271[k];

        t_347[k] = f_8 * pff_140[k]
                   + f_4 * pc_z[k] * pgf_230[k];

        t_348[k] = pb_x[k] * pfg0_273[k]
                   + f_8 * pff_183[k]
                   - f_7 * pc_x[k] * pfg1_273[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pc_x, pc_z, pff_141, pff_185, pff_186, \
                         pff_187, pgd0_143, pgd1_143, pgf_231, pgf_235, pgf_236, \
                         pgf_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_8 * pff_141[k]
                   + f_4 * pc_z[k] * pgf_231[k];

        t_350[k] = f_0 * pff_185[k]
                   + f_5 * pgd0_143[k]
                   - f_6 * pgd1_143[k]
                   + f_4 * pc_x[k] * pgf_235[k];

        t_351[k] = f_0 * pff_186[k]
                   + f_4 * pc_x[k] * pgf_236[k];

        t_352[k] = f_0 * pff_187[k]
                   + f_4 * pc_x[k] * pgf_237[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, pb_x, pc_x, pc_z, pfg0_280, pff_146, \
                         pff_188, pff_189, pfg1_280, pgf_236, pgf_238, \
                         pgf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_0 * pff_188[k]
                   + f_4 * pc_x[k] * pgf_238[k];

        t_354[k] = f_0 * pff_189[k]
                   + f_4 * pc_x[k] * pgf_239[k];

        t_355[k] = pb_x[k] * pfg0_280[k]
                   - f_7 * pc_x[k] * pfg1_280[k];

        t_356[k] = f_8 * pff_146[k]
                   + f_4 * pc_z[k] * pgf_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, pa_y, pb_x, pc_x, pc_y, sgg0_135, \
                         sgg1_135, pfg0_282, pfg0_283, pfg0_284, pfg1_282, pfg1_283, \
                         pfg1_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pb_x[k] * pfg0_282[k]
                   - f_7 * pc_x[k] * pfg1_282[k];

        t_358[k] = pb_x[k] * pfg0_283[k]
                   - f_7 * pc_x[k] * pfg1_283[k];

        t_359[k] = pb_x[k] * pfg0_284[k]
                   - f_7 * pc_x[k] * pfg1_284[k];

        t_360[k] = pa_y[k] * sgg0_135[k]
                   - f_7 * pc_y[k] * sgg1_135[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, pa_y, pc_y, pc_z, sgg0_136, sgg0_138, sgf_90, \
                         sgf_91, sgg1_136, sgg1_138, pff_150, pgf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = pa_y[k] * sgg0_136[k]
                   + f_0 * sgf_90[k]
                   - f_7 * pc_y[k] * sgg1_136[k];

        t_362[k] = f_9 * pff_150[k]
                   + f_4 * pc_z[k] * pgf_240[k];

        t_363[k] = pa_y[k] * sgg0_138[k]
                   + f_8 * sgf_91[k]
                   - f_7 * pc_y[k] * sgg1_138[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, pa_y, pc_x, pc_y, pc_z, sgg0_140, \
                         sgg1_140, pff_151, pff_196, pff_197, pgf_241, pgf_246, \
                         pgf_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_9 * pff_151[k]
                   + f_4 * pc_z[k] * pgf_241[k];

        t_365[k] = pa_y[k] * sgg0_140[k]
                   - f_7 * pc_y[k] * sgg1_140[k];

        t_366[k] = f_0 * pff_196[k]
                   + f_4 * pc_x[k] * pgf_246[k];

        t_367[k] = f_0 * pff_197[k]
                   + f_4 * pc_x[k] * pgf_247[k];
    }
}

static auto
compute_prim_pgg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgg0, const size_t sgf,
                                                          const size_t sgg1, const size_t pdg0,
                                                          const size_t pdg1, const size_t pfg0,
                                                          const size_t pff, const size_t pfg1,
                                                          const size_t pgd0, const size_t pgd1,
                                                          const size_t pgf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);
    const auto f_14 = 1.0 / p;
    const auto f_15 = gamma / (p * q);

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
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgg0_0 = buffer.data(sgg0 + 0);
    const auto *sgg0_2 = buffer.data(sgg0 + 2);
    const auto *sgg0_3 = buffer.data(sgg0 + 3);
    const auto *sgg0_5 = buffer.data(sgg0 + 5);
    const auto *sgg0_10 = buffer.data(sgg0 + 10);
    const auto *sgg0_14 = buffer.data(sgg0 + 14);
    const auto *sgg0_15 = buffer.data(sgg0 + 15);
    const auto *sgg0_18 = buffer.data(sgg0 + 18);
    const auto *sgg0_25 = buffer.data(sgg0 + 25);
    const auto *sgg0_26 = buffer.data(sgg0 + 26);
    const auto *sgg0_27 = buffer.data(sgg0 + 27);
    const auto *sgg0_149 = buffer.data(sgg0 + 149);
    const auto *sgg0_210 = buffer.data(sgg0 + 210);
    const auto *sgg0_215 = buffer.data(sgg0 + 215);
    const auto *sgg0_220 = buffer.data(sgg0 + 220);
    const auto *sgg0_222 = buffer.data(sgg0 + 222);
    const auto *sgg0_224 = buffer.data(sgg0 + 224);

    const auto *sgf_0 = buffer.data(sgf + 0);
    const auto *sgf_2 = buffer.data(sgf + 2);
    const auto *sgf_9 = buffer.data(sgf + 9);
    const auto *sgf_16 = buffer.data(sgf + 16);
    const auto *sgf_17 = buffer.data(sgf + 17);
    const auto *sgf_99 = buffer.data(sgf + 99);
    const auto *sgf_106 = buffer.data(sgf + 106);
    const auto *sgf_109 = buffer.data(sgf + 109);
    const auto *sgf_119 = buffer.data(sgf + 119);
    const auto *sgf_129 = buffer.data(sgf + 129);
    const auto *sgf_136 = buffer.data(sgf + 136);
    const auto *sgf_139 = buffer.data(sgf + 139);
    const auto *sgf_146 = buffer.data(sgf + 146);
    const auto *sgf_148 = buffer.data(sgf + 148);
    const auto *sgf_149 = buffer.data(sgf + 149);

    const auto *sgg1_0 = buffer.data(sgg1 + 0);
    const auto *sgg1_2 = buffer.data(sgg1 + 2);
    const auto *sgg1_3 = buffer.data(sgg1 + 3);
    const auto *sgg1_5 = buffer.data(sgg1 + 5);
    const auto *sgg1_10 = buffer.data(sgg1 + 10);
    const auto *sgg1_14 = buffer.data(sgg1 + 14);
    const auto *sgg1_15 = buffer.data(sgg1 + 15);
    const auto *sgg1_18 = buffer.data(sgg1 + 18);
    const auto *sgg1_25 = buffer.data(sgg1 + 25);
    const auto *sgg1_26 = buffer.data(sgg1 + 26);
    const auto *sgg1_27 = buffer.data(sgg1 + 27);
    const auto *sgg1_149 = buffer.data(sgg1 + 149);
    const auto *sgg1_210 = buffer.data(sgg1 + 210);
    const auto *sgg1_215 = buffer.data(sgg1 + 215);
    const auto *sgg1_220 = buffer.data(sgg1 + 220);
    const auto *sgg1_222 = buffer.data(sgg1 + 222);
    const auto *sgg1_224 = buffer.data(sgg1 + 224);

    const auto *pdg0_145 = buffer.data(pdg0 + 145);
    const auto *pdg0_224 = buffer.data(pdg0 + 224);

    const auto *pdg1_145 = buffer.data(pdg1 + 145);
    const auto *pdg1_224 = buffer.data(pdg1 + 224);

    const auto *pfg0_240 = buffer.data(pfg0 + 240);
    const auto *pfg0_241 = buffer.data(pfg0 + 241);
    const auto *pfg0_243 = buffer.data(pfg0 + 243);
    const auto *pfg0_250 = buffer.data(pfg0 + 250);
    const auto *pfg0_252 = buffer.data(pfg0 + 252);
    const auto *pfg0_265 = buffer.data(pfg0 + 265);
    const auto *pfg0_295 = buffer.data(pfg0 + 295);
    const auto *pfg0_297 = buffer.data(pfg0 + 297);
    const auto *pfg0_302 = buffer.data(pfg0 + 302);
    const auto *pfg0_305 = buffer.data(pfg0 + 305);
    const auto *pfg0_314 = buffer.data(pfg0 + 314);
    const auto *pfg0_344 = buffer.data(pfg0 + 344);

    const auto *pff_156 = buffer.data(pff + 156);
    const auto *pff_160 = buffer.data(pff + 160);
    const auto *pff_161 = buffer.data(pff + 161);
    const auto *pff_166 = buffer.data(pff + 166);
    const auto *pff_167 = buffer.data(pff + 167);
    const auto *pff_169 = buffer.data(pff + 169);
    const auto *pff_170 = buffer.data(pff + 170);
    const auto *pff_171 = buffer.data(pff + 171);
    const auto *pff_176 = buffer.data(pff + 176);
    const auto *pff_177 = buffer.data(pff + 177);
    const auto *pff_179 = buffer.data(pff + 179);
    const auto *pff_180 = buffer.data(pff + 180);
    const auto *pff_181 = buffer.data(pff + 181);
    const auto *pff_186 = buffer.data(pff + 186);
    const auto *pff_187 = buffer.data(pff + 187);
    const auto *pff_189 = buffer.data(pff + 189);
    const auto *pff_190 = buffer.data(pff + 190);
    const auto *pff_191 = buffer.data(pff + 191);
    const auto *pff_196 = buffer.data(pff + 196);
    const auto *pff_198 = buffer.data(pff + 198);
    const auto *pff_199 = buffer.data(pff + 199);
    const auto *pff_200 = buffer.data(pff + 200);
    const auto *pff_202 = buffer.data(pff + 202);
    const auto *pff_206 = buffer.data(pff + 206);
    const auto *pff_207 = buffer.data(pff + 207);
    const auto *pff_208 = buffer.data(pff + 208);
    const auto *pff_209 = buffer.data(pff + 209);
    const auto *pff_216 = buffer.data(pff + 216);
    const auto *pff_217 = buffer.data(pff + 217);
    const auto *pff_218 = buffer.data(pff + 218);
    const auto *pff_219 = buffer.data(pff + 219);
    const auto *pff_220 = buffer.data(pff + 220);
    const auto *pff_222 = buffer.data(pff + 222);
    const auto *pff_223 = buffer.data(pff + 223);
    const auto *pff_225 = buffer.data(pff + 225);
    const auto *pff_226 = buffer.data(pff + 226);
    const auto *pff_227 = buffer.data(pff + 227);
    const auto *pff_228 = buffer.data(pff + 228);
    const auto *pff_229 = buffer.data(pff + 229);

    const auto *pfg1_240 = buffer.data(pfg1 + 240);
    const auto *pfg1_241 = buffer.data(pfg1 + 241);
    const auto *pfg1_243 = buffer.data(pfg1 + 243);
    const auto *pfg1_250 = buffer.data(pfg1 + 250);
    const auto *pfg1_252 = buffer.data(pfg1 + 252);
    const auto *pfg1_265 = buffer.data(pfg1 + 265);
    const auto *pfg1_295 = buffer.data(pfg1 + 295);
    const auto *pfg1_297 = buffer.data(pfg1 + 297);
    const auto *pfg1_302 = buffer.data(pfg1 + 302);
    const auto *pfg1_305 = buffer.data(pfg1 + 305);
    const auto *pfg1_314 = buffer.data(pfg1 + 314);
    const auto *pfg1_344 = buffer.data(pfg1 + 344);

    const auto *pgd0_150 = buffer.data(pgd0 + 150);
    const auto *pgd0_151 = buffer.data(pgd0 + 151);
    const auto *pgd0_153 = buffer.data(pgd0 + 153);
    const auto *pgd0_155 = buffer.data(pgd0 + 155);
    const auto *pgd0_161 = buffer.data(pgd0 + 161);
    const auto *pgd0_162 = buffer.data(pgd0 + 162);
    const auto *pgd0_163 = buffer.data(pgd0 + 163);
    const auto *pgd0_165 = buffer.data(pgd0 + 165);
    const auto *pgd0_167 = buffer.data(pgd0 + 167);
    const auto *pgd0_168 = buffer.data(pgd0 + 168);
    const auto *pgd0_169 = buffer.data(pgd0 + 169);
    const auto *pgd0_171 = buffer.data(pgd0 + 171);
    const auto *pgd0_173 = buffer.data(pgd0 + 173);
    const auto *pgd0_175 = buffer.data(pgd0 + 175);
    const auto *pgd0_177 = buffer.data(pgd0 + 177);
    const auto *pgd0_184 = buffer.data(pgd0 + 184);
    const auto *pgd0_185 = buffer.data(pgd0 + 185);
    const auto *pgd0_192 = buffer.data(pgd0 + 192);
    const auto *pgd0_194 = buffer.data(pgd0 + 194);
    const auto *pgd0_195 = buffer.data(pgd0 + 195);
    const auto *pgd0_196 = buffer.data(pgd0 + 196);
    const auto *pgd0_197 = buffer.data(pgd0 + 197);

    const auto *pgd1_150 = buffer.data(pgd1 + 150);
    const auto *pgd1_151 = buffer.data(pgd1 + 151);
    const auto *pgd1_153 = buffer.data(pgd1 + 153);
    const auto *pgd1_155 = buffer.data(pgd1 + 155);
    const auto *pgd1_161 = buffer.data(pgd1 + 161);
    const auto *pgd1_162 = buffer.data(pgd1 + 162);
    const auto *pgd1_163 = buffer.data(pgd1 + 163);
    const auto *pgd1_165 = buffer.data(pgd1 + 165);
    const auto *pgd1_167 = buffer.data(pgd1 + 167);
    const auto *pgd1_168 = buffer.data(pgd1 + 168);
    const auto *pgd1_169 = buffer.data(pgd1 + 169);
    const auto *pgd1_171 = buffer.data(pgd1 + 171);
    const auto *pgd1_173 = buffer.data(pgd1 + 173);
    const auto *pgd1_175 = buffer.data(pgd1 + 175);
    const auto *pgd1_177 = buffer.data(pgd1 + 177);
    const auto *pgd1_184 = buffer.data(pgd1 + 184);
    const auto *pgd1_185 = buffer.data(pgd1 + 185);
    const auto *pgd1_192 = buffer.data(pgd1 + 192);
    const auto *pgd1_194 = buffer.data(pgd1 + 194);
    const auto *pgd1_195 = buffer.data(pgd1 + 195);
    const auto *pgd1_196 = buffer.data(pgd1 + 196);
    const auto *pgd1_197 = buffer.data(pgd1 + 197);

    const auto *pgf_246 = buffer.data(pgf + 246);
    const auto *pgf_248 = buffer.data(pgf + 248);
    const auto *pgf_249 = buffer.data(pgf + 249);
    const auto *pgf_250 = buffer.data(pgf + 250);
    const auto *pgf_251 = buffer.data(pgf + 251);
    const auto *pgf_253 = buffer.data(pgf + 253);
    const auto *pgf_255 = buffer.data(pgf + 255);
    const auto *pgf_256 = buffer.data(pgf + 256);
    const auto *pgf_257 = buffer.data(pgf + 257);
    const auto *pgf_258 = buffer.data(pgf + 258);
    const auto *pgf_259 = buffer.data(pgf + 259);
    const auto *pgf_260 = buffer.data(pgf + 260);
    const auto *pgf_261 = buffer.data(pgf + 261);
    const auto *pgf_265 = buffer.data(pgf + 265);
    const auto *pgf_266 = buffer.data(pgf + 266);
    const auto *pgf_267 = buffer.data(pgf + 267);
    const auto *pgf_268 = buffer.data(pgf + 268);
    const auto *pgf_269 = buffer.data(pgf + 269);
    const auto *pgf_270 = buffer.data(pgf + 270);
    const auto *pgf_271 = buffer.data(pgf + 271);
    const auto *pgf_273 = buffer.data(pgf + 273);
    const auto *pgf_275 = buffer.data(pgf + 275);
    const auto *pgf_276 = buffer.data(pgf + 276);
    const auto *pgf_277 = buffer.data(pgf + 277);
    const auto *pgf_278 = buffer.data(pgf + 278);
    const auto *pgf_279 = buffer.data(pgf + 279);
    const auto *pgf_280 = buffer.data(pgf + 280);
    const auto *pgf_281 = buffer.data(pgf + 281);
    const auto *pgf_283 = buffer.data(pgf + 283);
    const auto *pgf_285 = buffer.data(pgf + 285);
    const auto *pgf_286 = buffer.data(pgf + 286);
    const auto *pgf_287 = buffer.data(pgf + 287);
    const auto *pgf_288 = buffer.data(pgf + 288);
    const auto *pgf_289 = buffer.data(pgf + 289);
    const auto *pgf_290 = buffer.data(pgf + 290);
    const auto *pgf_291 = buffer.data(pgf + 291);
    const auto *pgf_293 = buffer.data(pgf + 293);
    const auto *pgf_296 = buffer.data(pgf + 296);
    const auto *pgf_297 = buffer.data(pgf + 297);
    const auto *pgf_298 = buffer.data(pgf + 298);
    const auto *pgf_299 = buffer.data(pgf + 299);
    const auto *pgf_300 = buffer.data(pgf + 300);
    const auto *pgf_302 = buffer.data(pgf + 302);
    const auto *pgf_306 = buffer.data(pgf + 306);
    const auto *pgf_307 = buffer.data(pgf + 307);
    const auto *pgf_308 = buffer.data(pgf + 308);
    const auto *pgf_309 = buffer.data(pgf + 309);
    const auto *pgf_310 = buffer.data(pgf + 310);
    const auto *pgf_312 = buffer.data(pgf + 312);
    const auto *pgf_316 = buffer.data(pgf + 316);
    const auto *pgf_317 = buffer.data(pgf + 317);
    const auto *pgf_318 = buffer.data(pgf + 318);
    const auto *pgf_319 = buffer.data(pgf + 319);
    const auto *pgf_320 = buffer.data(pgf + 320);
    const auto *pgf_322 = buffer.data(pgf + 322);
    const auto *pgf_323 = buffer.data(pgf + 323);
    const auto *pgf_325 = buffer.data(pgf + 325);
    const auto *pgf_326 = buffer.data(pgf + 326);
    const auto *pgf_327 = buffer.data(pgf + 327);
    const auto *pgf_328 = buffer.data(pgf + 328);
    const auto *pgf_329 = buffer.data(pgf + 329);

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pb_x, pc_x, pc_z, pfg0_295, pff_156, \
                         pff_198, pff_199, pfg1_295, pgf_246, pgf_248, \
                         pgf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_0 * pff_198[k]
                   + f_4 * pc_x[k] * pgf_248[k];

        t_369[k] = f_0 * pff_199[k]
                   + f_4 * pc_x[k] * pgf_249[k];

        t_370[k] = pb_x[k] * pfg0_295[k]
                   - f_7 * pc_x[k] * pfg1_295[k];

        t_371[k] = f_9 * pff_156[k]
                   + f_4 * pc_z[k] * pgf_246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pa_y, pb_x, pc_x, pc_y, sgg0_149, sgf_99, \
                         sgg1_149, pfg0_297, pfg1_297, pgf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = pb_x[k] * pfg0_297[k]
                   - f_7 * pc_x[k] * pfg1_297[k];

        t_373[k] = f_0 * sgf_99[k]
                   + f_4 * pc_y[k] * pgf_249[k];

        t_374[k] = pa_y[k] * sgg0_149[k]
                   - f_7 * pc_y[k] * sgg1_149[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, pc_x, pc_z, pgd0_150, pgd0_151, \
                         pgd0_153, pgd1_150, pgd1_151, pgd1_153, pgf_250, pgf_251, \
                         pgf_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_2 * pgd0_150[k]
                   - f_3 * pgd1_150[k]
                   + f_4 * pc_x[k] * pgf_250[k];

        t_376[k] = f_12 * pgd0_151[k]
                   - f_13 * pgd1_151[k]
                   + f_4 * pc_x[k] * pgf_251[k];

        t_377[k] = f_4 * pc_z[k] * pgf_250[k];

        t_378[k] = f_5 * pgd0_153[k]
                   - f_6 * pgd1_153[k]
                   + f_4 * pc_x[k] * pgf_253[k];

        t_379[k] = f_4 * pc_z[k] * pgf_251[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, pc_x, pgd0_155, pgd1_155, pgf_255, \
                         pgf_256, pgf_257, pgf_258, pgf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_5 * pgd0_155[k]
                   - f_6 * pgd1_155[k]
                   + f_4 * pc_x[k] * pgf_255[k];

        t_381[k] = f_4 * pc_x[k] * pgf_256[k];

        t_382[k] = f_4 * pc_x[k] * pgf_257[k];

        t_383[k] = f_4 * pc_x[k] * pgf_258[k];

        t_384[k] = f_4 * pc_x[k] * pgf_259[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pc_y, pc_z, sgf_106, sgf_109, pff_166, \
                         pff_169, pgd0_153, pgd1_153, pgf_256, pgf_257, \
                         pgf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_0 * sgf_106[k]
                   + f_1 * pff_166[k]
                   + f_2 * pgd0_153[k]
                   - f_3 * pgd1_153[k]
                   + f_4 * pc_y[k] * pgf_256[k];

        t_386[k] = f_4 * pc_z[k] * pgf_256[k];

        t_387[k] = f_5 * pgd0_153[k]
                   - f_6 * pgd1_153[k]
                   + f_4 * pc_z[k] * pgf_257[k];

        t_388[k] = f_0 * sgf_109[k]
                   + f_1 * pff_169[k]
                   + f_4 * pc_y[k] * pgf_259[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, pb_z, pc_z, pfg0_240, pfg0_241, pff_160, \
                         pfg1_240, pfg1_241, pgd0_155, pgd1_155, pgf_259, \
                         pgf_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_2 * pgd0_155[k]
                   - f_3 * pgd1_155[k]
                   + f_4 * pc_z[k] * pgf_259[k];

        t_390[k] = pb_z[k] * pfg0_240[k]
                   - f_7 * pc_z[k] * pfg1_240[k];

        t_391[k] = pb_z[k] * pfg0_241[k]
                   - f_7 * pc_z[k] * pfg1_241[k];

        t_392[k] = f_0 * pff_160[k]
                   + f_4 * pc_z[k] * pgf_260[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pb_z, pc_x, pc_z, pfg0_243, pff_161, \
                         pfg1_243, pgd0_161, pgd1_161, pgf_261, pgf_265, \
                         pgf_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = pb_z[k] * pfg0_243[k]
                   - f_7 * pc_z[k] * pfg1_243[k];

        t_394[k] = f_0 * pff_161[k]
                   + f_4 * pc_z[k] * pgf_261[k];

        t_395[k] = f_5 * pgd0_161[k]
                   - f_6 * pgd1_161[k]
                   + f_4 * pc_x[k] * pgf_265[k];

        t_396[k] = f_4 * pc_x[k] * pgf_266[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, pb_z, pc_x, pc_z, pfg0_250, \
                         pff_166, pfg1_250, pgf_266, pgf_267, pgf_268, \
                         pgf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_4 * pc_x[k] * pgf_267[k];

        t_398[k] = f_4 * pc_x[k] * pgf_268[k];

        t_399[k] = f_4 * pc_x[k] * pgf_269[k];

        t_400[k] = pb_z[k] * pfg0_250[k]
                   - f_7 * pc_z[k] * pfg1_250[k];

        t_401[k] = f_0 * pff_166[k]
                   + f_4 * pc_z[k] * pgf_266[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pb_z, pc_y, pc_z, sgf_119, pfg0_252, pff_167, \
                         pff_169, pff_179, pfg1_252, pgd0_161, pgd1_161, \
                         pgf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pb_z[k] * pfg0_252[k]
                   + f_8 * pff_167[k]
                   - f_7 * pc_z[k] * pfg1_252[k];

        t_403[k] = f_0 * sgf_119[k]
                   + f_9 * pff_179[k]
                   + f_4 * pc_y[k] * pgf_269[k];

        t_404[k] = f_0 * pff_169[k]
                   + f_2 * pgd0_161[k]
                   - f_3 * pgd1_161[k]
                   + f_4 * pc_z[k] * pgf_269[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pc_x, pc_z, pff_170, pgd0_162, pgd0_163, \
                         pgd0_165, pgd1_162, pgd1_163, pgd1_165, pgf_270, pgf_271, \
                         pgf_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_2 * pgd0_162[k]
                   - f_3 * pgd1_162[k]
                   + f_4 * pc_x[k] * pgf_270[k];

        t_406[k] = f_12 * pgd0_163[k]
                   - f_13 * pgd1_163[k]
                   + f_4 * pc_x[k] * pgf_271[k];

        t_407[k] = f_8 * pff_170[k]
                   + f_4 * pc_z[k] * pgf_270[k];

        t_408[k] = f_5 * pgd0_165[k]
                   - f_6 * pgd1_165[k]
                   + f_4 * pc_x[k] * pgf_273[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, pc_z, pff_171, pgd0_167, \
                         pgd1_167, pgf_271, pgf_275, pgf_276, pgf_277, \
                         pgf_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_8 * pff_171[k]
                   + f_4 * pc_z[k] * pgf_271[k];

        t_410[k] = f_5 * pgd0_167[k]
                   - f_6 * pgd1_167[k]
                   + f_4 * pc_x[k] * pgf_275[k];

        t_411[k] = f_4 * pc_x[k] * pgf_276[k];

        t_412[k] = f_4 * pc_x[k] * pgf_277[k];

        t_413[k] = f_4 * pc_x[k] * pgf_278[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pb_z, pc_x, pc_z, pdg0_145, pdg1_145, pfg0_265, \
                         pff_176, pfg1_265, pgf_276, pgf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_4 * pc_x[k] * pgf_279[k];

        t_415[k] = f_10 * pdg0_145[k]
                   - f_11 * pdg1_145[k]
                   + pb_z[k] * pfg0_265[k]
                   - f_7 * pc_z[k] * pfg1_265[k];

        t_416[k] = f_8 * pff_176[k]
                   + f_4 * pc_z[k] * pgf_276[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_y, pc_z, sgf_129, pff_177, pff_179, pff_189, \
                         pgd0_165, pgd0_167, pgd1_165, pgd1_167, pgf_277, \
                         pgf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_8 * pff_177[k]
                   + f_5 * pgd0_165[k]
                   - f_6 * pgd1_165[k]
                   + f_4 * pc_z[k] * pgf_277[k];

        t_418[k] = f_0 * sgf_129[k]
                   + f_8 * pff_189[k]
                   + f_4 * pc_y[k] * pgf_279[k];

        t_419[k] = f_8 * pff_179[k]
                   + f_2 * pgd0_167[k]
                   - f_3 * pgd1_167[k]
                   + f_4 * pc_z[k] * pgf_279[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_z, pff_180, pgd0_168, pgd0_169, \
                         pgd0_171, pgd1_168, pgd1_169, pgd1_171, pgf_280, pgf_281, \
                         pgf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_2 * pgd0_168[k]
                   - f_3 * pgd1_168[k]
                   + f_4 * pc_x[k] * pgf_280[k];

        t_421[k] = f_12 * pgd0_169[k]
                   - f_13 * pgd1_169[k]
                   + f_4 * pc_x[k] * pgf_281[k];

        t_422[k] = f_9 * pff_180[k]
                   + f_4 * pc_z[k] * pgf_280[k];

        t_423[k] = f_5 * pgd0_171[k]
                   - f_6 * pgd1_171[k]
                   + f_4 * pc_x[k] * pgf_283[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, pc_x, pc_z, pff_181, pgd0_173, \
                         pgd1_173, pgf_281, pgf_285, pgf_286, pgf_287, \
                         pgf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_9 * pff_181[k]
                   + f_4 * pc_z[k] * pgf_281[k];

        t_425[k] = f_5 * pgd0_173[k]
                   - f_6 * pgd1_173[k]
                   + f_4 * pc_x[k] * pgf_285[k];

        t_426[k] = f_4 * pc_x[k] * pgf_286[k];

        t_427[k] = f_4 * pc_x[k] * pgf_287[k];

        t_428[k] = f_4 * pc_x[k] * pgf_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pc_x, pc_y, pc_z, sgf_136, pff_186, \
                         pff_187, pff_196, pgd0_171, pgd1_171, pgf_286, pgf_287, \
                         pgf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_4 * pc_x[k] * pgf_289[k];

        t_430[k] = f_0 * sgf_136[k]
                   + f_0 * pff_196[k]
                   + f_2 * pgd0_171[k]
                   - f_3 * pgd1_171[k]
                   + f_4 * pc_y[k] * pgf_286[k];

        t_431[k] = f_9 * pff_186[k]
                   + f_4 * pc_z[k] * pgf_286[k];

        t_432[k] = f_9 * pff_187[k]
                   + f_5 * pgd0_171[k]
                   - f_6 * pgd1_171[k]
                   + f_4 * pc_z[k] * pgf_287[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, pa_y, pc_y, pc_z, sgg0_210, sgf_139, sgg1_210, \
                         pff_189, pff_199, pgd0_173, pgd1_173, \
                         pgf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_0 * sgf_139[k]
                   + f_0 * pff_199[k]
                   + f_4 * pc_y[k] * pgf_289[k];

        t_434[k] = f_9 * pff_189[k]
                   + f_2 * pgd0_173[k]
                   - f_3 * pgd1_173[k]
                   + f_4 * pc_z[k] * pgf_289[k];

        t_435[k] = pa_y[k] * sgg0_210[k]
                   - f_7 * pc_y[k] * sgg1_210[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, pc_x, pc_z, pff_190, pff_191, pgd0_175, \
                         pgd0_177, pgd1_175, pgd1_177, pgf_290, pgf_291, \
                         pgf_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_12 * pgd0_175[k]
                   - f_13 * pgd1_175[k]
                   + f_4 * pc_x[k] * pgf_291[k];

        t_437[k] = f_1 * pff_190[k]
                   + f_4 * pc_z[k] * pgf_290[k];

        t_438[k] = f_5 * pgd0_177[k]
                   - f_6 * pgd1_177[k]
                   + f_4 * pc_x[k] * pgf_293[k];

        t_439[k] = f_1 * pff_191[k]
                   + f_4 * pc_z[k] * pgf_291[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pa_y, pc_x, pc_y, sgg0_215, \
                         sgg1_215, pgf_296, pgf_297, pgf_298, pgf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = pa_y[k] * sgg0_215[k]
                   - f_7 * pc_y[k] * sgg1_215[k];

        t_441[k] = f_4 * pc_x[k] * pgf_296[k];

        t_442[k] = f_4 * pc_x[k] * pgf_297[k];

        t_443[k] = f_4 * pc_x[k] * pgf_298[k];

        t_444[k] = f_4 * pc_x[k] * pgf_299[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, pa_y, pc_y, pc_z, sgg0_220, sgg0_222, sgf_146, \
                         sgf_148, sgg1_220, sgg1_222, pff_196, \
                         pgf_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pa_y[k] * sgg0_220[k]
                   + f_1 * sgf_146[k]
                   - f_7 * pc_y[k] * sgg1_220[k];

        t_446[k] = f_1 * pff_196[k]
                   + f_4 * pc_z[k] * pgf_296[k];

        t_447[k] = pa_y[k] * sgg0_222[k]
                   + f_8 * sgf_148[k]
                   - f_7 * pc_y[k] * sgg1_222[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pa_y, pa_z, pc_y, pc_z, sgg0_0, sgg0_224, \
                         sgf_149, sgg1_0, sgg1_224, pgf_299, pgf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_0 * sgf_149[k]
                   + f_4 * pc_y[k] * pgf_299[k];

        t_449[k] = pa_y[k] * sgg0_224[k]
                   - f_7 * pc_y[k] * sgg1_224[k];

        t_450[k] = pa_z[k] * sgg0_0[k]
                   - f_7 * pc_z[k] * sgg1_0[k];

        t_451[k] = f_4 * pc_y[k] * pgf_300[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pa_z, pc_y, pc_z, sgg0_2, sgg0_3, sgg0_5, \
                         sgf_0, sgf_2, sgg1_2, sgg1_3, sgg1_5, \
                         pgf_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = pa_z[k] * sgg0_2[k]
                   + f_0 * sgf_0[k]
                   - f_7 * pc_z[k] * sgg1_2[k];

        t_453[k] = pa_z[k] * sgg0_3[k]
                   - f_7 * pc_z[k] * sgg1_3[k];

        t_454[k] = f_4 * pc_y[k] * pgf_302[k];

        t_455[k] = pa_z[k] * sgg0_5[k]
                   + f_8 * sgf_2[k]
                   - f_7 * pc_z[k] * sgg1_5[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, pff_206, pff_207, pff_208, pff_209, \
                         pgf_306, pgf_307, pgf_308, pgf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_1 * pff_206[k]
                   + f_4 * pc_x[k] * pgf_306[k];

        t_457[k] = f_1 * pff_207[k]
                   + f_4 * pc_x[k] * pgf_307[k];

        t_458[k] = f_1 * pff_208[k]
                   + f_4 * pc_x[k] * pgf_308[k];

        t_459[k] = f_1 * pff_209[k]
                   + f_4 * pc_x[k] * pgf_309[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pa_z, pc_y, pc_z, sgg0_10, sgg1_10, \
                         pgd0_184, pgd0_185, pgd1_184, pgd1_185, pgf_307, pgf_308, \
                         pgf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = pa_z[k] * sgg0_10[k]
                   - f_7 * pc_z[k] * sgg1_10[k];

        t_461[k] = f_12 * pgd0_184[k]
                   - f_13 * pgd1_184[k]
                   + f_4 * pc_y[k] * pgf_307[k];

        t_462[k] = f_5 * pgd0_185[k]
                   - f_6 * pgd1_185[k]
                   + f_4 * pc_y[k] * pgf_308[k];

        t_463[k] = f_4 * pc_y[k] * pgf_309[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pa_z, pc_y, pc_z, sgg0_14, sgg0_15, sgf_9, \
                         sgg1_14, sgg1_15, pff_200, pgf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = pa_z[k] * sgg0_14[k]
                   + f_1 * sgf_9[k]
                   - f_7 * pc_z[k] * sgg1_14[k];

        t_465[k] = pa_z[k] * sgg0_15[k]
                   - f_7 * pc_z[k] * sgg1_15[k];

        t_466[k] = f_0 * pff_200[k]
                   + f_4 * pc_y[k] * pgf_310[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, pa_z, pb_y, pc_y, pc_z, sgg0_18, sgg1_18, \
                         pfg0_302, pfg0_305, pff_202, pfg1_302, pfg1_305, \
                         pgf_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = pb_y[k] * pfg0_302[k]
                   - f_7 * pc_y[k] * pfg1_302[k];

        t_468[k] = pa_z[k] * sgg0_18[k]
                   - f_7 * pc_z[k] * sgg1_18[k];

        t_469[k] = f_0 * pff_202[k]
                   + f_4 * pc_y[k] * pgf_312[k];

        t_470[k] = pb_y[k] * pfg0_305[k]
                   - f_7 * pc_y[k] * pfg1_305[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, pff_216, pff_217, pff_218, pff_219, \
                         pgf_316, pgf_317, pgf_318, pgf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_9 * pff_216[k]
                   + f_4 * pc_x[k] * pgf_316[k];

        t_472[k] = f_9 * pff_217[k]
                   + f_4 * pc_x[k] * pgf_317[k];

        t_473[k] = f_9 * pff_218[k]
                   + f_4 * pc_x[k] * pgf_318[k];

        t_474[k] = f_9 * pff_219[k]
                   + f_4 * pc_x[k] * pgf_319[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, pa_z, pc_z, sgg0_25, sgg0_26, sgg0_27, sgf_16, \
                         sgf_17, sgg1_25, sgg1_26, sgg1_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = pa_z[k] * sgg0_25[k]
                   - f_7 * pc_z[k] * sgg1_25[k];

        t_476[k] = pa_z[k] * sgg0_26[k]
                   + f_0 * sgf_16[k]
                   - f_7 * pc_z[k] * sgg1_26[k];

        t_477[k] = pa_z[k] * sgg0_27[k]
                   + f_8 * sgf_17[k]
                   - f_7 * pc_z[k] * sgg1_27[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, pb_y, pc_x, pc_y, pfg0_314, pff_209, \
                         pff_220, pfg1_314, pgd0_192, pgd1_192, pgf_319, \
                         pgf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_0 * pff_209[k]
                   + f_4 * pc_y[k] * pgf_319[k];

        t_479[k] = pb_y[k] * pfg0_314[k]
                   - f_7 * pc_y[k] * pfg1_314[k];

        t_480[k] = f_9 * pff_220[k]
                   + f_2 * pgd0_192[k]
                   - f_3 * pgd1_192[k]
                   + f_4 * pc_x[k] * pgf_320[k];

        t_481[k] = f_4 * pc_y[k] * pgf_320[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pff_222, pff_223, pgd0_194, \
                         pgd0_195, pgd1_194, pgd1_195, pgf_322, \
                         pgf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_9 * pff_222[k]
                   + f_12 * pgd0_194[k]
                   - f_13 * pgd1_194[k]
                   + f_4 * pc_x[k] * pgf_322[k];

        t_483[k] = f_9 * pff_223[k]
                   + f_5 * pgd0_195[k]
                   - f_6 * pgd1_195[k]
                   + f_4 * pc_x[k] * pgf_323[k];

        t_484[k] = f_4 * pc_y[k] * pgf_322[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pc_x, pff_225, pff_226, pff_227, pff_228, \
                         pgd0_197, pgd1_197, pgf_325, pgf_326, pgf_327, \
                         pgf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_9 * pff_225[k]
                   + f_5 * pgd0_197[k]
                   - f_6 * pgd1_197[k]
                   + f_4 * pc_x[k] * pgf_325[k];

        t_486[k] = f_9 * pff_226[k]
                   + f_4 * pc_x[k] * pgf_326[k];

        t_487[k] = f_9 * pff_227[k]
                   + f_4 * pc_x[k] * pgf_327[k];

        t_488[k] = f_9 * pff_228[k]
                   + f_4 * pc_x[k] * pgf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_x, pc_y, pff_229, pgd0_195, pgd0_196, \
                         pgd1_195, pgd1_196, pgf_326, pgf_327, \
                         pgf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_9 * pff_229[k]
                   + f_4 * pc_x[k] * pgf_329[k];

        t_490[k] = f_2 * pgd0_195[k]
                   - f_3 * pgd1_195[k]
                   + f_4 * pc_y[k] * pgf_326[k];

        t_491[k] = f_12 * pgd0_196[k]
                   - f_13 * pgd1_196[k]
                   + f_4 * pc_y[k] * pgf_327[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, pb_x, pc_x, pc_y, pdg0_224, pdg1_224, pfg0_344, \
                         pfg1_344, pgd0_197, pgd1_197, pgf_328, \
                         pgf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_5 * pgd0_197[k]
                   - f_6 * pgd1_197[k]
                   + f_4 * pc_y[k] * pgf_328[k];

        t_493[k] = f_4 * pc_y[k] * pgf_329[k];

        t_494[k] = f_14 * pdg0_224[k]
                   - f_15 * pdg1_224[k]
                   + pb_x[k] * pfg0_344[k]
                   - f_7 * pc_x[k] * pfg1_344[k];
    }
}

static auto
compute_prim_pgg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgg0, const size_t sgf,
                                                          const size_t sgg1, const size_t pdg0,
                                                          const size_t pdg1, const size_t pfg0,
                                                          const size_t pff, const size_t pfg1,
                                                          const size_t pgd0, const size_t pgd1,
                                                          const size_t pgf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgg0_45 = buffer.data(sgg0 + 45);
    const auto *sgg0_47 = buffer.data(sgg0 + 47);
    const auto *sgg0_48 = buffer.data(sgg0 + 48);
    const auto *sgg0_50 = buffer.data(sgg0 + 50);
    const auto *sgg0_55 = buffer.data(sgg0 + 55);
    const auto *sgg0_56 = buffer.data(sgg0 + 56);
    const auto *sgg0_57 = buffer.data(sgg0 + 57);
    const auto *sgg0_59 = buffer.data(sgg0 + 59);
    const auto *sgg0_90 = buffer.data(sgg0 + 90);
    const auto *sgg0_92 = buffer.data(sgg0 + 92);
    const auto *sgg0_93 = buffer.data(sgg0 + 93);
    const auto *sgg0_95 = buffer.data(sgg0 + 95);
    const auto *sgg0_100 = buffer.data(sgg0 + 100);
    const auto *sgg0_150 = buffer.data(sgg0 + 150);
    const auto *sgg0_153 = buffer.data(sgg0 + 153);
    const auto *sgg0_160 = buffer.data(sgg0 + 160);
    const auto *sgg0_161 = buffer.data(sgg0 + 161);
    const auto *sgg0_162 = buffer.data(sgg0 + 162);
    const auto *sgg0_164 = buffer.data(sgg0 + 164);

    const auto *sgf_30 = buffer.data(sgf + 30);
    const auto *sgf_32 = buffer.data(sgf + 32);
    const auto *sgf_36 = buffer.data(sgf + 36);
    const auto *sgf_37 = buffer.data(sgf + 37);
    const auto *sgf_39 = buffer.data(sgf + 39);
    const auto *sgf_60 = buffer.data(sgf + 60);
    const auto *sgf_62 = buffer.data(sgf + 62);
    const auto *sgf_106 = buffer.data(sgf + 106);
    const auto *sgf_107 = buffer.data(sgf + 107);
    const auto *sgf_109 = buffer.data(sgf + 109);

    const auto *sgg1_45 = buffer.data(sgg1 + 45);
    const auto *sgg1_47 = buffer.data(sgg1 + 47);
    const auto *sgg1_48 = buffer.data(sgg1 + 48);
    const auto *sgg1_50 = buffer.data(sgg1 + 50);
    const auto *sgg1_55 = buffer.data(sgg1 + 55);
    const auto *sgg1_56 = buffer.data(sgg1 + 56);
    const auto *sgg1_57 = buffer.data(sgg1 + 57);
    const auto *sgg1_59 = buffer.data(sgg1 + 59);
    const auto *sgg1_90 = buffer.data(sgg1 + 90);
    const auto *sgg1_92 = buffer.data(sgg1 + 92);
    const auto *sgg1_93 = buffer.data(sgg1 + 93);
    const auto *sgg1_95 = buffer.data(sgg1 + 95);
    const auto *sgg1_100 = buffer.data(sgg1 + 100);
    const auto *sgg1_150 = buffer.data(sgg1 + 150);
    const auto *sgg1_153 = buffer.data(sgg1 + 153);
    const auto *sgg1_160 = buffer.data(sgg1 + 160);
    const auto *sgg1_161 = buffer.data(sgg1 + 161);
    const auto *sgg1_162 = buffer.data(sgg1 + 162);
    const auto *sgg1_164 = buffer.data(sgg1 + 164);

    const auto *pdg0_269 = buffer.data(pdg0 + 269);

    const auto *pdg1_269 = buffer.data(pdg1 + 269);

    const auto *pfg0_330 = buffer.data(pfg0 + 330);
    const auto *pfg0_332 = buffer.data(pfg0 + 332);
    const auto *pfg0_335 = buffer.data(pfg0 + 335);
    const auto *pfg0_344 = buffer.data(pfg0 + 344);
    const auto *pfg0_375 = buffer.data(pfg0 + 375);
    const auto *pfg0_377 = buffer.data(pfg0 + 377);
    const auto *pfg0_380 = buffer.data(pfg0 + 380);
    const auto *pfg0_389 = buffer.data(pfg0 + 389);
    const auto *pfg0_401 = buffer.data(pfg0 + 401);
    const auto *pfg0_402 = buffer.data(pfg0 + 402);
    const auto *pfg0_404 = buffer.data(pfg0 + 404);
    const auto *pfg0_407 = buffer.data(pfg0 + 407);
    const auto *pfg0_410 = buffer.data(pfg0 + 410);
    const auto *pfg0_415 = buffer.data(pfg0 + 415);
    const auto *pfg0_416 = buffer.data(pfg0 + 416);
    const auto *pfg0_417 = buffer.data(pfg0 + 417);
    const auto *pfg0_419 = buffer.data(pfg0 + 419);
    const auto *pfg0_423 = buffer.data(pfg0 + 423);
    const auto *pfg0_430 = buffer.data(pfg0 + 430);
    const auto *pfg0_431 = buffer.data(pfg0 + 431);
    const auto *pfg0_432 = buffer.data(pfg0 + 432);
    const auto *pfg0_434 = buffer.data(pfg0 + 434);
    const auto *pfg0_435 = buffer.data(pfg0 + 435);
    const auto *pfg0_437 = buffer.data(pfg0 + 437);
    const auto *pfg0_438 = buffer.data(pfg0 + 438);
    const auto *pfg0_440 = buffer.data(pfg0 + 440);
    const auto *pfg0_445 = buffer.data(pfg0 + 445);
    const auto *pfg0_446 = buffer.data(pfg0 + 446);
    const auto *pfg0_447 = buffer.data(pfg0 + 447);
    const auto *pfg0_449 = buffer.data(pfg0 + 449);

    const auto *pff_210 = buffer.data(pff + 210);
    const auto *pff_212 = buffer.data(pff + 212);
    const auto *pff_219 = buffer.data(pff + 219);
    const auto *pff_220 = buffer.data(pff + 220);
    const auto *pff_222 = buffer.data(pff + 222);
    const auto *pff_226 = buffer.data(pff + 226);
    const auto *pff_227 = buffer.data(pff + 227);
    const auto *pff_228 = buffer.data(pff + 228);
    const auto *pff_229 = buffer.data(pff + 229);
    const auto *pff_230 = buffer.data(pff + 230);
    const auto *pff_232 = buffer.data(pff + 232);
    const auto *pff_236 = buffer.data(pff + 236);
    const auto *pff_237 = buffer.data(pff + 237);
    const auto *pff_238 = buffer.data(pff + 238);
    const auto *pff_239 = buffer.data(pff + 239);
    const auto *pff_240 = buffer.data(pff + 240);
    const auto *pff_242 = buffer.data(pff + 242);
    const auto *pff_243 = buffer.data(pff + 243);
    const auto *pff_246 = buffer.data(pff + 246);
    const auto *pff_247 = buffer.data(pff + 247);
    const auto *pff_248 = buffer.data(pff + 248);
    const auto *pff_249 = buffer.data(pff + 249);
    const auto *pff_250 = buffer.data(pff + 250);
    const auto *pff_252 = buffer.data(pff + 252);
    const auto *pff_253 = buffer.data(pff + 253);
    const auto *pff_255 = buffer.data(pff + 255);
    const auto *pff_256 = buffer.data(pff + 256);
    const auto *pff_257 = buffer.data(pff + 257);
    const auto *pff_258 = buffer.data(pff + 258);
    const auto *pff_259 = buffer.data(pff + 259);
    const auto *pff_260 = buffer.data(pff + 260);
    const auto *pff_262 = buffer.data(pff + 262);
    const auto *pff_266 = buffer.data(pff + 266);
    const auto *pff_267 = buffer.data(pff + 267);
    const auto *pff_268 = buffer.data(pff + 268);
    const auto *pff_269 = buffer.data(pff + 269);
    const auto *pff_270 = buffer.data(pff + 270);
    const auto *pff_272 = buffer.data(pff + 272);
    const auto *pff_273 = buffer.data(pff + 273);
    const auto *pff_275 = buffer.data(pff + 275);
    const auto *pff_276 = buffer.data(pff + 276);
    const auto *pff_277 = buffer.data(pff + 277);
    const auto *pff_278 = buffer.data(pff + 278);
    const auto *pff_279 = buffer.data(pff + 279);
    const auto *pff_283 = buffer.data(pff + 283);
    const auto *pff_286 = buffer.data(pff + 286);
    const auto *pff_287 = buffer.data(pff + 287);
    const auto *pff_288 = buffer.data(pff + 288);
    const auto *pff_289 = buffer.data(pff + 289);
    const auto *pff_290 = buffer.data(pff + 290);
    const auto *pff_292 = buffer.data(pff + 292);
    const auto *pff_293 = buffer.data(pff + 293);
    const auto *pff_295 = buffer.data(pff + 295);
    const auto *pff_296 = buffer.data(pff + 296);
    const auto *pff_297 = buffer.data(pff + 297);
    const auto *pff_298 = buffer.data(pff + 298);
    const auto *pff_299 = buffer.data(pff + 299);

    const auto *pfg1_330 = buffer.data(pfg1 + 330);
    const auto *pfg1_332 = buffer.data(pfg1 + 332);
    const auto *pfg1_335 = buffer.data(pfg1 + 335);
    const auto *pfg1_344 = buffer.data(pfg1 + 344);
    const auto *pfg1_375 = buffer.data(pfg1 + 375);
    const auto *pfg1_377 = buffer.data(pfg1 + 377);
    const auto *pfg1_380 = buffer.data(pfg1 + 380);
    const auto *pfg1_389 = buffer.data(pfg1 + 389);
    const auto *pfg1_401 = buffer.data(pfg1 + 401);
    const auto *pfg1_402 = buffer.data(pfg1 + 402);
    const auto *pfg1_404 = buffer.data(pfg1 + 404);
    const auto *pfg1_407 = buffer.data(pfg1 + 407);
    const auto *pfg1_410 = buffer.data(pfg1 + 410);
    const auto *pfg1_415 = buffer.data(pfg1 + 415);
    const auto *pfg1_416 = buffer.data(pfg1 + 416);
    const auto *pfg1_417 = buffer.data(pfg1 + 417);
    const auto *pfg1_419 = buffer.data(pfg1 + 419);
    const auto *pfg1_423 = buffer.data(pfg1 + 423);
    const auto *pfg1_430 = buffer.data(pfg1 + 430);
    const auto *pfg1_431 = buffer.data(pfg1 + 431);
    const auto *pfg1_432 = buffer.data(pfg1 + 432);
    const auto *pfg1_434 = buffer.data(pfg1 + 434);
    const auto *pfg1_435 = buffer.data(pfg1 + 435);
    const auto *pfg1_437 = buffer.data(pfg1 + 437);
    const auto *pfg1_438 = buffer.data(pfg1 + 438);
    const auto *pfg1_440 = buffer.data(pfg1 + 440);
    const auto *pfg1_445 = buffer.data(pfg1 + 445);
    const auto *pfg1_446 = buffer.data(pfg1 + 446);
    const auto *pfg1_447 = buffer.data(pfg1 + 447);
    const auto *pfg1_449 = buffer.data(pfg1 + 449);

    const auto *pgd0_207 = buffer.data(pgd0 + 207);
    const auto *pgd0_208 = buffer.data(pgd0 + 208);
    const auto *pgd0_209 = buffer.data(pgd0 + 209);
    const auto *pgd0_210 = buffer.data(pgd0 + 210);
    const auto *pgd0_212 = buffer.data(pgd0 + 212);
    const auto *pgd0_213 = buffer.data(pgd0 + 213);
    const auto *pgd0_214 = buffer.data(pgd0 + 214);
    const auto *pgd0_215 = buffer.data(pgd0 + 215);
    const auto *pgd0_222 = buffer.data(pgd0 + 222);
    const auto *pgd0_225 = buffer.data(pgd0 + 225);
    const auto *pgd0_242 = buffer.data(pgd0 + 242);
    const auto *pgd0_245 = buffer.data(pgd0 + 245);
    const auto *pgd0_246 = buffer.data(pgd0 + 246);

    const auto *pgd1_207 = buffer.data(pgd1 + 207);
    const auto *pgd1_208 = buffer.data(pgd1 + 208);
    const auto *pgd1_209 = buffer.data(pgd1 + 209);
    const auto *pgd1_210 = buffer.data(pgd1 + 210);
    const auto *pgd1_212 = buffer.data(pgd1 + 212);
    const auto *pgd1_213 = buffer.data(pgd1 + 213);
    const auto *pgd1_214 = buffer.data(pgd1 + 214);
    const auto *pgd1_215 = buffer.data(pgd1 + 215);
    const auto *pgd1_222 = buffer.data(pgd1 + 222);
    const auto *pgd1_225 = buffer.data(pgd1 + 225);
    const auto *pgd1_242 = buffer.data(pgd1 + 242);
    const auto *pgd1_245 = buffer.data(pgd1 + 245);
    const auto *pgd1_246 = buffer.data(pgd1 + 246);

    const auto *pgf_330 = buffer.data(pgf + 330);
    const auto *pgf_332 = buffer.data(pgf + 332);
    const auto *pgf_336 = buffer.data(pgf + 336);
    const auto *pgf_337 = buffer.data(pgf + 337);
    const auto *pgf_338 = buffer.data(pgf + 338);
    const auto *pgf_339 = buffer.data(pgf + 339);
    const auto *pgf_340 = buffer.data(pgf + 340);
    const auto *pgf_342 = buffer.data(pgf + 342);
    const auto *pgf_343 = buffer.data(pgf + 343);
    const auto *pgf_346 = buffer.data(pgf + 346);
    const auto *pgf_347 = buffer.data(pgf + 347);
    const auto *pgf_348 = buffer.data(pgf + 348);
    const auto *pgf_349 = buffer.data(pgf + 349);
    const auto *pgf_350 = buffer.data(pgf + 350);
    const auto *pgf_352 = buffer.data(pgf + 352);
    const auto *pgf_353 = buffer.data(pgf + 353);
    const auto *pgf_355 = buffer.data(pgf + 355);
    const auto *pgf_356 = buffer.data(pgf + 356);
    const auto *pgf_357 = buffer.data(pgf + 357);
    const auto *pgf_358 = buffer.data(pgf + 358);
    const auto *pgf_359 = buffer.data(pgf + 359);
    const auto *pgf_360 = buffer.data(pgf + 360);
    const auto *pgf_362 = buffer.data(pgf + 362);
    const auto *pgf_366 = buffer.data(pgf + 366);
    const auto *pgf_367 = buffer.data(pgf + 367);
    const auto *pgf_368 = buffer.data(pgf + 368);
    const auto *pgf_369 = buffer.data(pgf + 369);
    const auto *pgf_370 = buffer.data(pgf + 370);
    const auto *pgf_372 = buffer.data(pgf + 372);
    const auto *pgf_373 = buffer.data(pgf + 373);
    const auto *pgf_376 = buffer.data(pgf + 376);
    const auto *pgf_377 = buffer.data(pgf + 377);
    const auto *pgf_378 = buffer.data(pgf + 378);
    const auto *pgf_379 = buffer.data(pgf + 379);
    const auto *pgf_380 = buffer.data(pgf + 380);
    const auto *pgf_382 = buffer.data(pgf + 382);
    const auto *pgf_386 = buffer.data(pgf + 386);
    const auto *pgf_387 = buffer.data(pgf + 387);
    const auto *pgf_388 = buffer.data(pgf + 388);
    const auto *pgf_389 = buffer.data(pgf + 389);
    const auto *pgf_390 = buffer.data(pgf + 390);
    const auto *pgf_392 = buffer.data(pgf + 392);
    const auto *pgf_396 = buffer.data(pgf + 396);
    const auto *pgf_397 = buffer.data(pgf + 397);
    const auto *pgf_398 = buffer.data(pgf + 398);
    const auto *pgf_399 = buffer.data(pgf + 399);
    const auto *pgf_400 = buffer.data(pgf + 400);
    const auto *pgf_402 = buffer.data(pgf + 402);
    const auto *pgf_405 = buffer.data(pgf + 405);
    const auto *pgf_406 = buffer.data(pgf + 406);
    const auto *pgf_407 = buffer.data(pgf + 407);
    const auto *pgf_408 = buffer.data(pgf + 408);
    const auto *pgf_409 = buffer.data(pgf + 409);
    const auto *pgf_410 = buffer.data(pgf + 410);

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pa_z, pc_y, pc_z, sgg0_45, sgg0_47, \
                         sgg0_48, sgf_30, sgg1_45, sgg1_47, sgg1_48, pff_210, \
                         pgf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * sgg0_45[k]
                   - f_7 * pc_z[k] * sgg1_45[k];

        t_496[k] = f_8 * pff_210[k]
                   + f_4 * pc_y[k] * pgf_330[k];

        t_497[k] = pa_z[k] * sgg0_47[k]
                   + f_0 * sgf_30[k]
                   - f_7 * pc_z[k] * sgg1_47[k];

        t_498[k] = pa_z[k] * sgg0_48[k]
                   - f_7 * pc_z[k] * sgg1_48[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pa_z, pc_x, pc_y, pc_z, sgg0_50, sgf_32, \
                         sgg1_50, pff_212, pff_236, pgf_332, pgf_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_8 * pff_212[k]
                   + f_4 * pc_y[k] * pgf_332[k];

        t_500[k] = pa_z[k] * sgg0_50[k]
                   + f_8 * sgf_32[k]
                   - f_7 * pc_z[k] * sgg1_50[k];

        t_501[k] = f_8 * pff_236[k]
                   + f_4 * pc_x[k] * pgf_336[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, pa_z, pc_x, pc_z, sgg0_55, sgg1_55, \
                         pff_237, pff_238, pff_239, pgf_337, pgf_338, \
                         pgf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_8 * pff_237[k]
                   + f_4 * pc_x[k] * pgf_337[k];

        t_503[k] = f_8 * pff_238[k]
                   + f_4 * pc_x[k] * pgf_338[k];

        t_504[k] = f_8 * pff_239[k]
                   + f_4 * pc_x[k] * pgf_339[k];

        t_505[k] = pa_z[k] * sgg0_55[k]
                   - f_7 * pc_z[k] * sgg1_55[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pa_z, pc_y, pc_z, sgg0_56, sgg0_57, sgf_36, \
                         sgf_37, sgg1_56, sgg1_57, pff_219, pgf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = pa_z[k] * sgg0_56[k]
                   + f_0 * sgf_36[k]
                   - f_7 * pc_z[k] * sgg1_56[k];

        t_507[k] = pa_z[k] * sgg0_57[k]
                   + f_8 * sgf_37[k]
                   - f_7 * pc_z[k] * sgg1_57[k];

        t_508[k] = f_8 * pff_219[k]
                   + f_4 * pc_y[k] * pgf_339[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_z, pb_y, pc_y, pc_z, sgg0_59, sgf_39, \
                         sgg1_59, pfg0_330, pff_220, pfg1_330, \
                         pgf_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = pa_z[k] * sgg0_59[k]
                   + f_1 * sgf_39[k]
                   - f_7 * pc_z[k] * sgg1_59[k];

        t_510[k] = pb_y[k] * pfg0_330[k]
                   - f_7 * pc_y[k] * pfg1_330[k];

        t_511[k] = f_0 * pff_220[k]
                   + f_4 * pc_y[k] * pgf_340[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pb_y, pc_x, pc_y, pfg0_332, pff_222, pff_243, \
                         pfg1_332, pgd0_207, pgd1_207, pgf_342, \
                         pgf_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pb_y[k] * pfg0_332[k]
                   - f_7 * pc_y[k] * pfg1_332[k];

        t_513[k] = f_8 * pff_243[k]
                   + f_5 * pgd0_207[k]
                   - f_6 * pgd1_207[k]
                   + f_4 * pc_x[k] * pgf_343[k];

        t_514[k] = f_0 * pff_222[k]
                   + f_4 * pc_y[k] * pgf_342[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pb_y, pc_x, pc_y, pfg0_335, pff_246, \
                         pff_247, pff_248, pfg1_335, pgf_346, pgf_347, \
                         pgf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = pb_y[k] * pfg0_335[k]
                   - f_7 * pc_y[k] * pfg1_335[k];

        t_516[k] = f_8 * pff_246[k]
                   + f_4 * pc_x[k] * pgf_346[k];

        t_517[k] = f_8 * pff_247[k]
                   + f_4 * pc_x[k] * pgf_347[k];

        t_518[k] = f_8 * pff_248[k]
                   + f_4 * pc_x[k] * pgf_348[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_x, pc_y, pff_226, pff_227, pff_249, pgd0_207, \
                         pgd0_208, pgd1_207, pgd1_208, pgf_346, pgf_347, \
                         pgf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_8 * pff_249[k]
                   + f_4 * pc_x[k] * pgf_349[k];

        t_520[k] = f_0 * pff_226[k]
                   + f_2 * pgd0_207[k]
                   - f_3 * pgd1_207[k]
                   + f_4 * pc_y[k] * pgf_346[k];

        t_521[k] = f_0 * pff_227[k]
                   + f_12 * pgd0_208[k]
                   - f_13 * pgd1_208[k]
                   + f_4 * pc_y[k] * pgf_347[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, pb_y, pc_y, pfg0_344, pff_228, pff_229, \
                         pfg1_344, pgd0_209, pgd1_209, pgf_348, \
                         pgf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_0 * pff_228[k]
                   + f_5 * pgd0_209[k]
                   - f_6 * pgd1_209[k]
                   + f_4 * pc_y[k] * pgf_348[k];

        t_523[k] = f_0 * pff_229[k]
                   + f_4 * pc_y[k] * pgf_349[k];

        t_524[k] = pb_y[k] * pfg0_344[k]
                   - f_7 * pc_y[k] * pfg1_344[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, pc_x, pc_y, pff_250, pff_252, pgd0_210, \
                         pgd0_212, pgd1_210, pgd1_212, pgf_350, \
                         pgf_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_8 * pff_250[k]
                   + f_2 * pgd0_210[k]
                   - f_3 * pgd1_210[k]
                   + f_4 * pc_x[k] * pgf_350[k];

        t_526[k] = f_4 * pc_y[k] * pgf_350[k];

        t_527[k] = f_8 * pff_252[k]
                   + f_12 * pgd0_212[k]
                   - f_13 * pgd1_212[k]
                   + f_4 * pc_x[k] * pgf_352[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, pc_x, pc_y, pff_253, pff_255, pgd0_213, \
                         pgd0_215, pgd1_213, pgd1_215, pgf_352, pgf_353, \
                         pgf_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_8 * pff_253[k]
                   + f_5 * pgd0_213[k]
                   - f_6 * pgd1_213[k]
                   + f_4 * pc_x[k] * pgf_353[k];

        t_529[k] = f_4 * pc_y[k] * pgf_352[k];

        t_530[k] = f_8 * pff_255[k]
                   + f_5 * pgd0_215[k]
                   - f_6 * pgd1_215[k]
                   + f_4 * pc_x[k] * pgf_355[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, pc_x, pff_256, pff_257, pff_258, pff_259, \
                         pgf_356, pgf_357, pgf_358, pgf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_8 * pff_256[k]
                   + f_4 * pc_x[k] * pgf_356[k];

        t_532[k] = f_8 * pff_257[k]
                   + f_4 * pc_x[k] * pgf_357[k];

        t_533[k] = f_8 * pff_258[k]
                   + f_4 * pc_x[k] * pgf_358[k];

        t_534[k] = f_8 * pff_259[k]
                   + f_4 * pc_x[k] * pgf_359[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, pc_y, pgd0_213, pgd0_214, pgd0_215, \
                         pgd1_213, pgd1_214, pgd1_215, pgf_356, pgf_357, pgf_358, \
                         pgf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = f_2 * pgd0_213[k]
                   - f_3 * pgd1_213[k]
                   + f_4 * pc_y[k] * pgf_356[k];

        t_536[k] = f_12 * pgd0_214[k]
                   - f_13 * pgd1_214[k]
                   + f_4 * pc_y[k] * pgf_357[k];

        t_537[k] = f_5 * pgd0_215[k]
                   - f_6 * pgd1_215[k]
                   + f_4 * pc_y[k] * pgf_358[k];

        t_538[k] = f_4 * pc_y[k] * pgf_359[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pa_z, pb_x, pc_x, pc_y, pc_z, sgg0_90, sgg1_90, \
                         pdg0_269, pdg1_269, pfg0_389, pff_230, pfg1_389, \
                         pgf_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_10 * pdg0_269[k]
                   - f_11 * pdg1_269[k]
                   + pb_x[k] * pfg0_389[k]
                   - f_7 * pc_x[k] * pfg1_389[k];

        t_540[k] = pa_z[k] * sgg0_90[k]
                   - f_7 * pc_z[k] * sgg1_90[k];

        t_541[k] = f_9 * pff_230[k]
                   + f_4 * pc_y[k] * pgf_360[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pa_z, pc_y, pc_z, sgg0_92, sgg0_93, sgf_60, \
                         sgg1_92, sgg1_93, pff_232, pgf_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = pa_z[k] * sgg0_92[k]
                   + f_0 * sgf_60[k]
                   - f_7 * pc_z[k] * sgg1_92[k];

        t_543[k] = pa_z[k] * sgg0_93[k]
                   - f_7 * pc_z[k] * sgg1_93[k];

        t_544[k] = f_9 * pff_232[k]
                   + f_4 * pc_y[k] * pgf_362[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pa_z, pc_x, pc_z, sgg0_95, sgf_62, \
                         sgg1_95, pff_266, pff_267, pff_268, pgf_366, pgf_367, \
                         pgf_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = pa_z[k] * sgg0_95[k]
                   + f_8 * sgf_62[k]
                   - f_7 * pc_z[k] * sgg1_95[k];

        t_546[k] = f_0 * pff_266[k]
                   + f_4 * pc_x[k] * pgf_366[k];

        t_547[k] = f_0 * pff_267[k]
                   + f_4 * pc_x[k] * pgf_367[k];

        t_548[k] = f_0 * pff_268[k]
                   + f_4 * pc_x[k] * pgf_368[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, pa_z, pb_x, pc_x, pc_z, sgg0_100, \
                         sgg1_100, pfg0_401, pfg0_402, pff_269, pfg1_401, pfg1_402, \
                         pgf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_0 * pff_269[k]
                   + f_4 * pc_x[k] * pgf_369[k];

        t_550[k] = pa_z[k] * sgg0_100[k]
                   - f_7 * pc_z[k] * sgg1_100[k];

        t_551[k] = pb_x[k] * pfg0_401[k]
                   - f_7 * pc_x[k] * pfg1_401[k];

        t_552[k] = pb_x[k] * pfg0_402[k]
                   - f_7 * pc_x[k] * pfg1_402[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pb_x, pc_x, pc_y, pfg0_404, pff_239, \
                         pff_240, pff_270, pfg1_404, pgd0_222, pgd1_222, pgf_369, \
                         pgf_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_9 * pff_239[k]
                   + f_4 * pc_y[k] * pgf_369[k];

        t_554[k] = pb_x[k] * pfg0_404[k]
                   - f_7 * pc_x[k] * pfg1_404[k];

        t_555[k] = f_0 * pff_270[k]
                   + f_2 * pgd0_222[k]
                   - f_3 * pgd1_222[k]
                   + f_4 * pc_x[k] * pgf_370[k];

        t_556[k] = f_8 * pff_240[k]
                   + f_4 * pc_y[k] * pgf_370[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pb_x, pc_x, pc_y, pfg0_407, pff_242, pff_272, \
                         pff_273, pfg1_407, pgd0_225, pgd1_225, pgf_372, \
                         pgf_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = pb_x[k] * pfg0_407[k]
                   + f_9 * pff_272[k]
                   - f_7 * pc_x[k] * pfg1_407[k];

        t_558[k] = f_0 * pff_273[k]
                   + f_5 * pgd0_225[k]
                   - f_6 * pgd1_225[k]
                   + f_4 * pc_x[k] * pgf_373[k];

        t_559[k] = f_8 * pff_242[k]
                   + f_4 * pc_y[k] * pgf_372[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, pb_x, pc_x, pfg0_410, pff_275, pff_276, \
                         pff_277, pff_278, pfg1_410, pgf_376, pgf_377, \
                         pgf_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = pb_x[k] * pfg0_410[k]
                   + f_8 * pff_275[k]
                   - f_7 * pc_x[k] * pfg1_410[k];

        t_561[k] = f_0 * pff_276[k]
                   + f_4 * pc_x[k] * pgf_376[k];

        t_562[k] = f_0 * pff_277[k]
                   + f_4 * pc_x[k] * pgf_377[k];

        t_563[k] = f_0 * pff_278[k]
                   + f_4 * pc_x[k] * pgf_378[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, pb_x, pc_x, pfg0_415, pfg0_416, pfg0_417, \
                         pff_279, pfg1_415, pfg1_416, pfg1_417, \
                         pgf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_0 * pff_279[k]
                   + f_4 * pc_x[k] * pgf_379[k];

        t_565[k] = pb_x[k] * pfg0_415[k]
                   - f_7 * pc_x[k] * pfg1_415[k];

        t_566[k] = pb_x[k] * pfg0_416[k]
                   - f_7 * pc_x[k] * pfg1_416[k];

        t_567[k] = pb_x[k] * pfg0_417[k]
                   - f_7 * pc_x[k] * pfg1_417[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, pb_x, pb_y, pc_x, pc_y, pfg0_375, \
                         pfg0_419, pff_249, pff_250, pfg1_375, pfg1_419, pgf_379, \
                         pgf_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_8 * pff_249[k]
                   + f_4 * pc_y[k] * pgf_379[k];

        t_569[k] = pb_x[k] * pfg0_419[k]
                   - f_7 * pc_x[k] * pfg1_419[k];

        t_570[k] = pb_y[k] * pfg0_375[k]
                   - f_7 * pc_y[k] * pfg1_375[k];

        t_571[k] = f_0 * pff_250[k]
                   + f_4 * pc_y[k] * pgf_380[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pb_x, pb_y, pc_x, pc_y, pfg0_377, pfg0_423, \
                         pff_252, pff_283, pfg1_377, pfg1_423, \
                         pgf_382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = pb_y[k] * pfg0_377[k]
                   - f_7 * pc_y[k] * pfg1_377[k];

        t_573[k] = pb_x[k] * pfg0_423[k]
                   + f_8 * pff_283[k]
                   - f_7 * pc_x[k] * pfg1_423[k];

        t_574[k] = f_0 * pff_252[k]
                   + f_4 * pc_y[k] * pgf_382[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, pb_y, pc_x, pc_y, pfg0_380, pff_286, \
                         pff_287, pff_288, pfg1_380, pgf_386, pgf_387, \
                         pgf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = pb_y[k] * pfg0_380[k]
                   - f_7 * pc_y[k] * pfg1_380[k];

        t_576[k] = f_0 * pff_286[k]
                   + f_4 * pc_x[k] * pgf_386[k];

        t_577[k] = f_0 * pff_287[k]
                   + f_4 * pc_x[k] * pgf_387[k];

        t_578[k] = f_0 * pff_288[k]
                   + f_4 * pc_x[k] * pgf_388[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, t_582, pb_x, pc_x, pfg0_430, pfg0_431, pfg0_432, \
                         pff_289, pfg1_430, pfg1_431, pfg1_432, \
                         pgf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_0 * pff_289[k]
                   + f_4 * pc_x[k] * pgf_389[k];

        t_580[k] = pb_x[k] * pfg0_430[k]
                   - f_7 * pc_x[k] * pfg1_430[k];

        t_581[k] = pb_x[k] * pfg0_431[k]
                   - f_7 * pc_x[k] * pfg1_431[k];

        t_582[k] = pb_x[k] * pfg0_432[k]
                   - f_7 * pc_x[k] * pfg1_432[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, t_586, pb_x, pc_x, pc_y, pfg0_434, pfg0_435, \
                         pff_259, pff_290, pfg1_434, pfg1_435, pgf_389, \
                         pgf_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_0 * pff_259[k]
                   + f_4 * pc_y[k] * pgf_389[k];

        t_584[k] = pb_x[k] * pfg0_434[k]
                   - f_7 * pc_x[k] * pfg1_434[k];

        t_585[k] = pb_x[k] * pfg0_435[k]
                   + f_1 * pff_290[k]
                   - f_7 * pc_x[k] * pfg1_435[k];

        t_586[k] = f_4 * pc_y[k] * pgf_390[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pb_x, pc_x, pc_y, pfg0_437, pfg0_438, pff_292, \
                         pff_293, pfg1_437, pfg1_438, pgf_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = pb_x[k] * pfg0_437[k]
                   + f_9 * pff_292[k]
                   - f_7 * pc_x[k] * pfg1_437[k];

        t_588[k] = pb_x[k] * pfg0_438[k]
                   + f_8 * pff_293[k]
                   - f_7 * pc_x[k] * pfg1_438[k];

        t_589[k] = f_4 * pc_y[k] * pgf_392[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pb_x, pc_x, pfg0_440, pff_295, pff_296, \
                         pff_297, pff_298, pfg1_440, pgf_396, pgf_397, \
                         pgf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = pb_x[k] * pfg0_440[k]
                   + f_8 * pff_295[k]
                   - f_7 * pc_x[k] * pfg1_440[k];

        t_591[k] = f_0 * pff_296[k]
                   + f_4 * pc_x[k] * pgf_396[k];

        t_592[k] = f_0 * pff_297[k]
                   + f_4 * pc_x[k] * pgf_397[k];

        t_593[k] = f_0 * pff_298[k]
                   + f_4 * pc_x[k] * pgf_398[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, t_598, pb_x, pc_x, pc_y, pfg0_445, \
                         pfg0_446, pfg0_447, pff_299, pfg1_445, pfg1_446, pfg1_447, \
                         pgf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_0 * pff_299[k]
                   + f_4 * pc_x[k] * pgf_399[k];

        t_595[k] = pb_x[k] * pfg0_445[k]
                   - f_7 * pc_x[k] * pfg1_445[k];

        t_596[k] = pb_x[k] * pfg0_446[k]
                   - f_7 * pc_x[k] * pfg1_446[k];

        t_597[k] = pb_x[k] * pfg0_447[k]
                   - f_7 * pc_x[k] * pfg1_447[k];

        t_598[k] = f_4 * pc_y[k] * pgf_399[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, pa_z, pb_x, pc_x, pc_y, pc_z, sgg0_150, \
                         sgg1_150, pfg0_449, pff_260, pfg1_449, \
                         pgf_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = pb_x[k] * pfg0_449[k]
                   - f_7 * pc_x[k] * pfg1_449[k];

        t_600[k] = pa_z[k] * sgg0_150[k]
                   - f_7 * pc_z[k] * sgg1_150[k];

        t_601[k] = f_1 * pff_260[k]
                   + f_4 * pc_y[k] * pgf_400[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, pa_z, pc_x, pc_y, pc_z, sgg0_153, sgg1_153, \
                         pff_262, pgd0_242, pgd1_242, pgf_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_12 * pgd0_242[k]
                   - f_13 * pgd1_242[k]
                   + f_4 * pc_x[k] * pgf_402[k];

        t_603[k] = pa_z[k] * sgg0_153[k]
                   - f_7 * pc_z[k] * sgg1_153[k];

        t_604[k] = f_1 * pff_262[k]
                   + f_4 * pc_y[k] * pgf_402[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, pc_x, pgd0_245, pgd1_245, pgf_405, \
                         pgf_406, pgf_407, pgf_408, pgf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_5 * pgd0_245[k]
                   - f_6 * pgd1_245[k]
                   + f_4 * pc_x[k] * pgf_405[k];

        t_606[k] = f_4 * pc_x[k] * pgf_406[k];

        t_607[k] = f_4 * pc_x[k] * pgf_407[k];

        t_608[k] = f_4 * pc_x[k] * pgf_408[k];

        t_609[k] = f_4 * pc_x[k] * pgf_409[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pa_z, pc_z, sgg0_160, sgg0_161, sgg0_162, \
                         sgf_106, sgf_107, sgg1_160, sgg1_161, \
                         sgg1_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pa_z[k] * sgg0_160[k]
                   - f_7 * pc_z[k] * sgg1_160[k];

        t_611[k] = pa_z[k] * sgg0_161[k]
                   + f_0 * sgf_106[k]
                   - f_7 * pc_z[k] * sgg1_161[k];

        t_612[k] = pa_z[k] * sgg0_162[k]
                   + f_8 * sgf_107[k]
                   - f_7 * pc_z[k] * sgg1_162[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, pa_z, pc_x, pc_y, pc_z, sgg0_164, sgf_109, \
                         sgg1_164, pff_269, pgd0_246, pgd1_246, pgf_409, \
                         pgf_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_1 * pff_269[k]
                   + f_4 * pc_y[k] * pgf_409[k];

        t_614[k] = pa_z[k] * sgg0_164[k]
                   + f_1 * sgf_109[k]
                   - f_7 * pc_z[k] * sgg1_164[k];

        t_615[k] = f_2 * pgd0_246[k]
                   - f_3 * pgd1_246[k]
                   + f_4 * pc_x[k] * pgf_410[k];
    }
}

static auto
compute_prim_pgg_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgf,
                                                          const size_t pdg0, const size_t pdg1,
                                                          const size_t pfg0, const size_t pff,
                                                          const size_t pfg1, const size_t pgd0,
                                                          const size_t pgd1, const size_t pgf,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 0.5 / p;
    const auto f_11 = 0.5 * gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);
    const auto f_14 = 1.0 / p;
    const auto f_15 = gamma / (p * q);

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgf_149 = buffer.data(sgf + 149);

    const auto *pdg0_254 = buffer.data(pdg0 + 254);
    const auto *pdg0_269 = buffer.data(pdg0 + 269);

    const auto *pdg1_254 = buffer.data(pdg1 + 254);
    const auto *pdg1_269 = buffer.data(pdg1 + 269);

    const auto *pfg0_419 = buffer.data(pfg0 + 419);
    const auto *pfg0_434 = buffer.data(pfg0 + 434);
    const auto *pfg0_435 = buffer.data(pfg0 + 435);
    const auto *pfg0_437 = buffer.data(pfg0 + 437);
    const auto *pfg0_440 = buffer.data(pfg0 + 440);
    const auto *pfg0_445 = buffer.data(pfg0 + 445);
    const auto *pfg0_446 = buffer.data(pfg0 + 446);
    const auto *pfg0_447 = buffer.data(pfg0 + 447);
    const auto *pfg0_449 = buffer.data(pfg0 + 449);

    const auto *pff_270 = buffer.data(pff + 270);
    const auto *pff_272 = buffer.data(pff + 272);
    const auto *pff_276 = buffer.data(pff + 276);
    const auto *pff_277 = buffer.data(pff + 277);
    const auto *pff_278 = buffer.data(pff + 278);
    const auto *pff_279 = buffer.data(pff + 279);
    const auto *pff_280 = buffer.data(pff + 280);
    const auto *pff_282 = buffer.data(pff + 282);
    const auto *pff_286 = buffer.data(pff + 286);
    const auto *pff_287 = buffer.data(pff + 287);
    const auto *pff_288 = buffer.data(pff + 288);
    const auto *pff_289 = buffer.data(pff + 289);
    const auto *pff_290 = buffer.data(pff + 290);
    const auto *pff_292 = buffer.data(pff + 292);
    const auto *pff_296 = buffer.data(pff + 296);
    const auto *pff_297 = buffer.data(pff + 297);
    const auto *pff_298 = buffer.data(pff + 298);
    const auto *pff_299 = buffer.data(pff + 299);

    const auto *pfg1_419 = buffer.data(pfg1 + 419);
    const auto *pfg1_434 = buffer.data(pfg1 + 434);
    const auto *pfg1_435 = buffer.data(pfg1 + 435);
    const auto *pfg1_437 = buffer.data(pfg1 + 437);
    const auto *pfg1_440 = buffer.data(pfg1 + 440);
    const auto *pfg1_445 = buffer.data(pfg1 + 445);
    const auto *pfg1_446 = buffer.data(pfg1 + 446);
    const auto *pfg1_447 = buffer.data(pfg1 + 447);
    const auto *pfg1_449 = buffer.data(pfg1 + 449);

    const auto *pgd0_248 = buffer.data(pgd0 + 248);
    const auto *pgd0_249 = buffer.data(pgd0 + 249);
    const auto *pgd0_250 = buffer.data(pgd0 + 250);
    const auto *pgd0_251 = buffer.data(pgd0 + 251);
    const auto *pgd0_252 = buffer.data(pgd0 + 252);
    const auto *pgd0_254 = buffer.data(pgd0 + 254);
    const auto *pgd0_255 = buffer.data(pgd0 + 255);
    const auto *pgd0_256 = buffer.data(pgd0 + 256);
    const auto *pgd0_257 = buffer.data(pgd0 + 257);
    const auto *pgd0_261 = buffer.data(pgd0 + 261);
    const auto *pgd0_264 = buffer.data(pgd0 + 264);
    const auto *pgd0_266 = buffer.data(pgd0 + 266);
    const auto *pgd0_267 = buffer.data(pgd0 + 267);
    const auto *pgd0_268 = buffer.data(pgd0 + 268);
    const auto *pgd0_269 = buffer.data(pgd0 + 269);

    const auto *pgd1_248 = buffer.data(pgd1 + 248);
    const auto *pgd1_249 = buffer.data(pgd1 + 249);
    const auto *pgd1_250 = buffer.data(pgd1 + 250);
    const auto *pgd1_251 = buffer.data(pgd1 + 251);
    const auto *pgd1_252 = buffer.data(pgd1 + 252);
    const auto *pgd1_254 = buffer.data(pgd1 + 254);
    const auto *pgd1_255 = buffer.data(pgd1 + 255);
    const auto *pgd1_256 = buffer.data(pgd1 + 256);
    const auto *pgd1_257 = buffer.data(pgd1 + 257);
    const auto *pgd1_261 = buffer.data(pgd1 + 261);
    const auto *pgd1_264 = buffer.data(pgd1 + 264);
    const auto *pgd1_266 = buffer.data(pgd1 + 266);
    const auto *pgd1_267 = buffer.data(pgd1 + 267);
    const auto *pgd1_268 = buffer.data(pgd1 + 268);
    const auto *pgd1_269 = buffer.data(pgd1 + 269);

    const auto *pgf_410 = buffer.data(pgf + 410);
    const auto *pgf_412 = buffer.data(pgf + 412);
    const auto *pgf_413 = buffer.data(pgf + 413);
    const auto *pgf_415 = buffer.data(pgf + 415);
    const auto *pgf_416 = buffer.data(pgf + 416);
    const auto *pgf_417 = buffer.data(pgf + 417);
    const auto *pgf_418 = buffer.data(pgf + 418);
    const auto *pgf_419 = buffer.data(pgf + 419);
    const auto *pgf_420 = buffer.data(pgf + 420);
    const auto *pgf_422 = buffer.data(pgf + 422);
    const auto *pgf_423 = buffer.data(pgf + 423);
    const auto *pgf_425 = buffer.data(pgf + 425);
    const auto *pgf_426 = buffer.data(pgf + 426);
    const auto *pgf_427 = buffer.data(pgf + 427);
    const auto *pgf_428 = buffer.data(pgf + 428);
    const auto *pgf_429 = buffer.data(pgf + 429);
    const auto *pgf_430 = buffer.data(pgf + 430);
    const auto *pgf_432 = buffer.data(pgf + 432);
    const auto *pgf_433 = buffer.data(pgf + 433);
    const auto *pgf_436 = buffer.data(pgf + 436);
    const auto *pgf_437 = buffer.data(pgf + 437);
    const auto *pgf_438 = buffer.data(pgf + 438);
    const auto *pgf_439 = buffer.data(pgf + 439);
    const auto *pgf_440 = buffer.data(pgf + 440);
    const auto *pgf_442 = buffer.data(pgf + 442);
    const auto *pgf_443 = buffer.data(pgf + 443);
    const auto *pgf_445 = buffer.data(pgf + 445);
    const auto *pgf_446 = buffer.data(pgf + 446);
    const auto *pgf_447 = buffer.data(pgf + 447);
    const auto *pgf_448 = buffer.data(pgf + 448);
    const auto *pgf_449 = buffer.data(pgf + 449);

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pc_x, pc_y, pff_270, pff_272, pgd0_248, \
                         pgd0_249, pgd1_248, pgd1_249, pgf_410, pgf_412, \
                         pgf_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = f_9 * pff_270[k]
                   + f_4 * pc_y[k] * pgf_410[k];

        t_617[k] = f_12 * pgd0_248[k]
                   - f_13 * pgd1_248[k]
                   + f_4 * pc_x[k] * pgf_412[k];

        t_618[k] = f_5 * pgd0_249[k]
                   - f_6 * pgd1_249[k]
                   + f_4 * pc_x[k] * pgf_413[k];

        t_619[k] = f_9 * pff_272[k]
                   + f_4 * pc_y[k] * pgf_412[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, t_624, pc_x, pgd0_251, pgd1_251, pgf_415, \
                         pgf_416, pgf_417, pgf_418, pgf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_5 * pgd0_251[k]
                   - f_6 * pgd1_251[k]
                   + f_4 * pc_x[k] * pgf_415[k];

        t_621[k] = f_4 * pc_x[k] * pgf_416[k];

        t_622[k] = f_4 * pc_x[k] * pgf_417[k];

        t_623[k] = f_4 * pc_x[k] * pgf_418[k];

        t_624[k] = f_4 * pc_x[k] * pgf_419[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, pc_y, pff_276, pff_277, pff_278, pgd0_249, \
                         pgd0_250, pgd0_251, pgd1_249, pgd1_250, pgd1_251, pgf_416, pgf_417, \
                         pgf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_9 * pff_276[k]
                   + f_2 * pgd0_249[k]
                   - f_3 * pgd1_249[k]
                   + f_4 * pc_y[k] * pgf_416[k];

        t_626[k] = f_9 * pff_277[k]
                   + f_12 * pgd0_250[k]
                   - f_13 * pgd1_250[k]
                   + f_4 * pc_y[k] * pgf_417[k];

        t_627[k] = f_9 * pff_278[k]
                   + f_5 * pgd0_251[k]
                   - f_6 * pgd1_251[k]
                   + f_4 * pc_y[k] * pgf_418[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, pb_y, pc_x, pc_y, pdg0_254, pdg1_254, pfg0_419, \
                         pff_279, pfg1_419, pgd0_252, pgd1_252, pgf_419, \
                         pgf_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_9 * pff_279[k]
                   + f_4 * pc_y[k] * pgf_419[k];

        t_629[k] = f_14 * pdg0_254[k]
                   - f_15 * pdg1_254[k]
                   + pb_y[k] * pfg0_419[k]
                   - f_7 * pc_y[k] * pfg1_419[k];

        t_630[k] = f_2 * pgd0_252[k]
                   - f_3 * pgd1_252[k]
                   + f_4 * pc_x[k] * pgf_420[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, pc_x, pc_y, pff_280, pff_282, pgd0_254, \
                         pgd0_255, pgd1_254, pgd1_255, pgf_420, pgf_422, \
                         pgf_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_8 * pff_280[k]
                   + f_4 * pc_y[k] * pgf_420[k];

        t_632[k] = f_12 * pgd0_254[k]
                   - f_13 * pgd1_254[k]
                   + f_4 * pc_x[k] * pgf_422[k];

        t_633[k] = f_5 * pgd0_255[k]
                   - f_6 * pgd1_255[k]
                   + f_4 * pc_x[k] * pgf_423[k];

        t_634[k] = f_8 * pff_282[k]
                   + f_4 * pc_y[k] * pgf_422[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, t_639, pc_x, pgd0_257, pgd1_257, pgf_425, \
                         pgf_426, pgf_427, pgf_428, pgf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_5 * pgd0_257[k]
                   - f_6 * pgd1_257[k]
                   + f_4 * pc_x[k] * pgf_425[k];

        t_636[k] = f_4 * pc_x[k] * pgf_426[k];

        t_637[k] = f_4 * pc_x[k] * pgf_427[k];

        t_638[k] = f_4 * pc_x[k] * pgf_428[k];

        t_639[k] = f_4 * pc_x[k] * pgf_429[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pff_286, pff_287, pff_288, pgd0_255, \
                         pgd0_256, pgd0_257, pgd1_255, pgd1_256, pgd1_257, pgf_426, pgf_427, \
                         pgf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_8 * pff_286[k]
                   + f_2 * pgd0_255[k]
                   - f_3 * pgd1_255[k]
                   + f_4 * pc_y[k] * pgf_426[k];

        t_641[k] = f_8 * pff_287[k]
                   + f_12 * pgd0_256[k]
                   - f_13 * pgd1_256[k]
                   + f_4 * pc_y[k] * pgf_427[k];

        t_642[k] = f_8 * pff_288[k]
                   + f_5 * pgd0_257[k]
                   - f_6 * pgd1_257[k]
                   + f_4 * pc_y[k] * pgf_428[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pb_y, pc_y, pdg0_269, pdg1_269, pfg0_434, \
                         pfg0_435, pff_289, pff_290, pfg1_434, pfg1_435, pgf_429, \
                         pgf_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_8 * pff_289[k]
                   + f_4 * pc_y[k] * pgf_429[k];

        t_644[k] = f_10 * pdg0_269[k]
                   - f_11 * pdg1_269[k]
                   + pb_y[k] * pfg0_434[k]
                   - f_7 * pc_y[k] * pfg1_434[k];

        t_645[k] = pb_y[k] * pfg0_435[k]
                   - f_7 * pc_y[k] * pfg1_435[k];

        t_646[k] = f_0 * pff_290[k]
                   + f_4 * pc_y[k] * pgf_430[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, pb_y, pc_x, pc_y, pfg0_437, pfg0_440, \
                         pff_292, pfg1_437, pfg1_440, pgd0_261, pgd1_261, pgf_432, \
                         pgf_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = pb_y[k] * pfg0_437[k]
                   - f_7 * pc_y[k] * pfg1_437[k];

        t_648[k] = f_5 * pgd0_261[k]
                   - f_6 * pgd1_261[k]
                   + f_4 * pc_x[k] * pgf_433[k];

        t_649[k] = f_0 * pff_292[k]
                   + f_4 * pc_y[k] * pgf_432[k];

        t_650[k] = pb_y[k] * pfg0_440[k]
                   - f_7 * pc_y[k] * pfg1_440[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, t_654, t_655, pb_y, pc_x, pc_y, pfg0_445, \
                         pff_296, pfg1_445, pgf_436, pgf_437, pgf_438, \
                         pgf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_4 * pc_x[k] * pgf_436[k];

        t_652[k] = f_4 * pc_x[k] * pgf_437[k];

        t_653[k] = f_4 * pc_x[k] * pgf_438[k];

        t_654[k] = f_4 * pc_x[k] * pgf_439[k];

        t_655[k] = pb_y[k] * pfg0_445[k]
                   + f_1 * pff_296[k]
                   - f_7 * pc_y[k] * pfg1_445[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, t_659, pb_y, pc_y, pfg0_446, pfg0_447, pfg0_449, \
                         pff_297, pff_298, pff_299, pfg1_446, pfg1_447, pfg1_449, \
                         pgf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = pb_y[k] * pfg0_446[k]
                   + f_9 * pff_297[k]
                   - f_7 * pc_y[k] * pfg1_446[k];

        t_657[k] = pb_y[k] * pfg0_447[k]
                   + f_8 * pff_298[k]
                   - f_7 * pc_y[k] * pfg1_447[k];

        t_658[k] = f_0 * pff_299[k]
                   + f_4 * pc_y[k] * pgf_439[k];

        t_659[k] = pb_y[k] * pfg0_449[k]
                   - f_7 * pc_y[k] * pfg1_449[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, pc_x, pc_y, pgd0_264, pgd0_266, \
                         pgd0_267, pgd1_264, pgd1_266, pgd1_267, pgf_440, pgf_442, \
                         pgf_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_2 * pgd0_264[k]
                   - f_3 * pgd1_264[k]
                   + f_4 * pc_x[k] * pgf_440[k];

        t_661[k] = f_4 * pc_y[k] * pgf_440[k];

        t_662[k] = f_12 * pgd0_266[k]
                   - f_13 * pgd1_266[k]
                   + f_4 * pc_x[k] * pgf_442[k];

        t_663[k] = f_5 * pgd0_267[k]
                   - f_6 * pgd1_267[k]
                   + f_4 * pc_x[k] * pgf_443[k];

        t_664[k] = f_4 * pc_y[k] * pgf_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, pc_x, pgd0_269, pgd1_269, pgf_445, \
                         pgf_446, pgf_447, pgf_448, pgf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_5 * pgd0_269[k]
                   - f_6 * pgd1_269[k]
                   + f_4 * pc_x[k] * pgf_445[k];

        t_666[k] = f_4 * pc_x[k] * pgf_446[k];

        t_667[k] = f_4 * pc_x[k] * pgf_447[k];

        t_668[k] = f_4 * pc_x[k] * pgf_448[k];

        t_669[k] = f_4 * pc_x[k] * pgf_449[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pc_y, pgd0_267, pgd0_268, pgd0_269, \
                         pgd1_267, pgd1_268, pgd1_269, pgf_446, pgf_447, pgf_448, \
                         pgf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_2 * pgd0_267[k]
                   - f_3 * pgd1_267[k]
                   + f_4 * pc_y[k] * pgf_446[k];

        t_671[k] = f_12 * pgd0_268[k]
                   - f_13 * pgd1_268[k]
                   + f_4 * pc_y[k] * pgf_447[k];

        t_672[k] = f_5 * pgd0_269[k]
                   - f_6 * pgd1_269[k]
                   + f_4 * pc_y[k] * pgf_448[k];

        t_673[k] = f_4 * pc_y[k] * pgf_449[k];
    }

#pragma omp simd aligned(t_674, pc_z, sgf_149, pff_299, pgd0_269, pgd1_269, \
                         pgf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_0 * sgf_149[k]
                   + f_1 * pff_299[k]
                   + f_2 * pgd0_269[k]
                   - f_3 * pgd1_269[k]
                   + f_4 * pc_z[k] * pgf_449[k];
    }
}

auto
compute_prim_pgg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sgg0,
                                                   const size_t sgf, const size_t sgg1,
                                                   const size_t pdg0, const size_t pdg1,
                                                   const size_t pfg0, const size_t pff,
                                                   const size_t pfg1, const size_t pgd0,
                                                   const size_t pgd1, const size_t pgf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pgg_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sgf, pdg0,
                                                              pdg1, pfg0, pff, pfg1, pgd0, pgd1,
                                                              pgf, ncols, gamma, p, q);

    compute_prim_pgg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sgg0,
                                                              sgf, sgg1, pfg0, pff, pfg1, pgd0,
                                                              pgd1, pgf, ncols, gamma, p, q);

    compute_prim_pgg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sgg0,
                                                              sgf, sgg1, pdg0, pdg1, pfg0, pff,
                                                              pfg1, pgd0, pgd1, pgf, ncols,
                                                              gamma, p, q);

    compute_prim_pgg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, sgg0,
                                                              sgf, sgg1, pdg0, pdg1, pfg0, pff,
                                                              pfg1, pgd0, pgd1, pgf, ncols,
                                                              gamma, p, q);

    compute_prim_pgg_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, sgg0,
                                                              sgf, sgg1, pdg0, pdg1, pfg0, pff,
                                                              pfg1, pgd0, pgd1, pgf, ncols,
                                                              gamma, p, q);

    compute_prim_pgg_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sgf, pdg0,
                                                              pdg1, pfg0, pff, pfg1, pgd0, pgd1,
                                                              pgf, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
