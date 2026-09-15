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


#include "SimdThreeCenterElectronRepulsionVrrRecPGI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgh,
                                                          const size_t pdi0, const size_t pdi1,
                                                          const size_t pfi0, const size_t pfh,
                                                          const size_t pfi1, const size_t pgg0,
                                                          const size_t pgg1, const size_t pgh,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);

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

    const auto *sgh_0 = buffer.data(sgh + 0);
    const auto *sgh_15 = buffer.data(sgh + 15);
    const auto *sgh_17 = buffer.data(sgh + 17);
    const auto *sgh_18 = buffer.data(sgh + 18);
    const auto *sgh_20 = buffer.data(sgh + 20);
    const auto *sgh_36 = buffer.data(sgh + 36);
    const auto *sgh_38 = buffer.data(sgh + 38);
    const auto *sgh_39 = buffer.data(sgh + 39);
    const auto *sgh_59 = buffer.data(sgh + 59);
    const auto *sgh_60 = buffer.data(sgh + 60);
    const auto *sgh_62 = buffer.data(sgh + 62);
    const auto *sgh_78 = buffer.data(sgh + 78);
    const auto *sgh_80 = buffer.data(sgh + 80);
    const auto *sgh_81 = buffer.data(sgh + 81);
    const auto *sgh_83 = buffer.data(sgh + 83);

    const auto *pdi0_0 = buffer.data(pdi0 + 0);

    const auto *pdi1_0 = buffer.data(pdi1 + 0);

    const auto *pfi0_0 = buffer.data(pfi0 + 0);
    const auto *pfi0_3 = buffer.data(pfi0 + 3);
    const auto *pfi0_5 = buffer.data(pfi0 + 5);
    const auto *pfi0_6 = buffer.data(pfi0 + 6);
    const auto *pfi0_9 = buffer.data(pfi0 + 9);
    const auto *pfi0_10 = buffer.data(pfi0 + 10);
    const auto *pfi0_12 = buffer.data(pfi0 + 12);
    const auto *pfi0_14 = buffer.data(pfi0 + 14);
    const auto *pfi0_15 = buffer.data(pfi0 + 15);
    const auto *pfi0_20 = buffer.data(pfi0 + 20);
    const auto *pfi0_21 = buffer.data(pfi0 + 21);
    const auto *pfi0_27 = buffer.data(pfi0 + 27);
    const auto *pfi0_28 = buffer.data(pfi0 + 28);
    const auto *pfi0_31 = buffer.data(pfi0 + 31);
    const auto *pfi0_34 = buffer.data(pfi0 + 34);
    const auto *pfi0_38 = buffer.data(pfi0 + 38);
    const auto *pfi0_56 = buffer.data(pfi0 + 56);
    const auto *pfi0_61 = buffer.data(pfi0 + 61);
    const auto *pfi0_65 = buffer.data(pfi0 + 65);
    const auto *pfi0_68 = buffer.data(pfi0 + 68);

    const auto *pfh_0 = buffer.data(pfh + 0);
    const auto *pfh_1 = buffer.data(pfh + 1);
    const auto *pfh_2 = buffer.data(pfh + 2);
    const auto *pfh_3 = buffer.data(pfh + 3);
    const auto *pfh_5 = buffer.data(pfh + 5);
    const auto *pfh_6 = buffer.data(pfh + 6);
    const auto *pfh_8 = buffer.data(pfh + 8);
    const auto *pfh_9 = buffer.data(pfh + 9);
    const auto *pfh_10 = buffer.data(pfh + 10);
    const auto *pfh_14 = buffer.data(pfh + 14);
    const auto *pfh_15 = buffer.data(pfh + 15);
    const auto *pfh_17 = buffer.data(pfh + 17);
    const auto *pfh_18 = buffer.data(pfh + 18);
    const auto *pfh_19 = buffer.data(pfh + 19);
    const auto *pfh_20 = buffer.data(pfh + 20);
    const auto *pfh_21 = buffer.data(pfh + 21);
    const auto *pfh_22 = buffer.data(pfh + 22);
    const auto *pfh_23 = buffer.data(pfh + 23);
    const auto *pfh_24 = buffer.data(pfh + 24);
    const auto *pfh_26 = buffer.data(pfh + 26);
    const auto *pfh_27 = buffer.data(pfh + 27);
    const auto *pfh_29 = buffer.data(pfh + 29);
    const auto *pfh_30 = buffer.data(pfh + 30);
    const auto *pfh_35 = buffer.data(pfh + 35);
    const auto *pfh_36 = buffer.data(pfh + 36);
    const auto *pfh_38 = buffer.data(pfh + 38);
    const auto *pfh_39 = buffer.data(pfh + 39);
    const auto *pfh_40 = buffer.data(pfh + 40);
    const auto *pfh_41 = buffer.data(pfh + 41);
    const auto *pfh_42 = buffer.data(pfh + 42);
    const auto *pfh_44 = buffer.data(pfh + 44);
    const auto *pfh_47 = buffer.data(pfh + 47);
    const auto *pfh_50 = buffer.data(pfh + 50);
    const auto *pfh_59 = buffer.data(pfh + 59);
    const auto *pfh_60 = buffer.data(pfh + 60);
    const auto *pfh_62 = buffer.data(pfh + 62);
    const auto *pfh_78 = buffer.data(pfh + 78);
    const auto *pfh_80 = buffer.data(pfh + 80);
    const auto *pfh_81 = buffer.data(pfh + 81);
    const auto *pfh_83 = buffer.data(pfh + 83);

    const auto *pfi1_0 = buffer.data(pfi1 + 0);
    const auto *pfi1_3 = buffer.data(pfi1 + 3);
    const auto *pfi1_5 = buffer.data(pfi1 + 5);
    const auto *pfi1_6 = buffer.data(pfi1 + 6);
    const auto *pfi1_9 = buffer.data(pfi1 + 9);
    const auto *pfi1_10 = buffer.data(pfi1 + 10);
    const auto *pfi1_12 = buffer.data(pfi1 + 12);
    const auto *pfi1_14 = buffer.data(pfi1 + 14);
    const auto *pfi1_15 = buffer.data(pfi1 + 15);
    const auto *pfi1_20 = buffer.data(pfi1 + 20);
    const auto *pfi1_21 = buffer.data(pfi1 + 21);
    const auto *pfi1_27 = buffer.data(pfi1 + 27);
    const auto *pfi1_28 = buffer.data(pfi1 + 28);
    const auto *pfi1_31 = buffer.data(pfi1 + 31);
    const auto *pfi1_34 = buffer.data(pfi1 + 34);
    const auto *pfi1_38 = buffer.data(pfi1 + 38);
    const auto *pfi1_56 = buffer.data(pfi1 + 56);
    const auto *pfi1_61 = buffer.data(pfi1 + 61);
    const auto *pfi1_65 = buffer.data(pfi1 + 65);
    const auto *pfi1_68 = buffer.data(pfi1 + 68);

    const auto *pgg0_0 = buffer.data(pgg0 + 0);
    const auto *pgg0_1 = buffer.data(pgg0 + 1);
    const auto *pgg0_2 = buffer.data(pgg0 + 2);
    const auto *pgg0_3 = buffer.data(pgg0 + 3);
    const auto *pgg0_5 = buffer.data(pgg0 + 5);
    const auto *pgg0_10 = buffer.data(pgg0 + 10);
    const auto *pgg0_12 = buffer.data(pgg0 + 12);
    const auto *pgg0_13 = buffer.data(pgg0 + 13);
    const auto *pgg0_14 = buffer.data(pgg0 + 14);
    const auto *pgg0_25 = buffer.data(pgg0 + 25);
    const auto *pgg0_27 = buffer.data(pgg0 + 27);
    const auto *pgg0_28 = buffer.data(pgg0 + 28);
    const auto *pgg0_29 = buffer.data(pgg0 + 29);
    const auto *pgg0_35 = buffer.data(pgg0 + 35);
    const auto *pgg0_42 = buffer.data(pgg0 + 42);
    const auto *pgg0_43 = buffer.data(pgg0 + 43);
    const auto *pgg0_44 = buffer.data(pgg0 + 44);
    const auto *pgg0_45 = buffer.data(pgg0 + 45);
    const auto *pgg0_46 = buffer.data(pgg0 + 46);
    const auto *pgg0_47 = buffer.data(pgg0 + 47);
    const auto *pgg0_48 = buffer.data(pgg0 + 48);
    const auto *pgg0_50 = buffer.data(pgg0 + 50);
    const auto *pgg0_55 = buffer.data(pgg0 + 55);
    const auto *pgg0_57 = buffer.data(pgg0 + 57);
    const auto *pgg0_58 = buffer.data(pgg0 + 58);
    const auto *pgg0_59 = buffer.data(pgg0 + 59);

    const auto *pgg1_0 = buffer.data(pgg1 + 0);
    const auto *pgg1_1 = buffer.data(pgg1 + 1);
    const auto *pgg1_2 = buffer.data(pgg1 + 2);
    const auto *pgg1_3 = buffer.data(pgg1 + 3);
    const auto *pgg1_5 = buffer.data(pgg1 + 5);
    const auto *pgg1_10 = buffer.data(pgg1 + 10);
    const auto *pgg1_12 = buffer.data(pgg1 + 12);
    const auto *pgg1_13 = buffer.data(pgg1 + 13);
    const auto *pgg1_14 = buffer.data(pgg1 + 14);
    const auto *pgg1_25 = buffer.data(pgg1 + 25);
    const auto *pgg1_27 = buffer.data(pgg1 + 27);
    const auto *pgg1_28 = buffer.data(pgg1 + 28);
    const auto *pgg1_29 = buffer.data(pgg1 + 29);
    const auto *pgg1_35 = buffer.data(pgg1 + 35);
    const auto *pgg1_42 = buffer.data(pgg1 + 42);
    const auto *pgg1_43 = buffer.data(pgg1 + 43);
    const auto *pgg1_44 = buffer.data(pgg1 + 44);
    const auto *pgg1_45 = buffer.data(pgg1 + 45);
    const auto *pgg1_46 = buffer.data(pgg1 + 46);
    const auto *pgg1_47 = buffer.data(pgg1 + 47);
    const auto *pgg1_48 = buffer.data(pgg1 + 48);
    const auto *pgg1_50 = buffer.data(pgg1 + 50);
    const auto *pgg1_55 = buffer.data(pgg1 + 55);
    const auto *pgg1_57 = buffer.data(pgg1 + 57);
    const auto *pgg1_58 = buffer.data(pgg1 + 58);
    const auto *pgg1_59 = buffer.data(pgg1 + 59);

    const auto *pgh_0 = buffer.data(pgh + 0);
    const auto *pgh_1 = buffer.data(pgh + 1);
    const auto *pgh_2 = buffer.data(pgh + 2);
    const auto *pgh_3 = buffer.data(pgh + 3);
    const auto *pgh_5 = buffer.data(pgh + 5);
    const auto *pgh_6 = buffer.data(pgh + 6);
    const auto *pgh_8 = buffer.data(pgh + 8);
    const auto *pgh_9 = buffer.data(pgh + 9);
    const auto *pgh_10 = buffer.data(pgh + 10);
    const auto *pgh_14 = buffer.data(pgh + 14);
    const auto *pgh_15 = buffer.data(pgh + 15);
    const auto *pgh_17 = buffer.data(pgh + 17);
    const auto *pgh_18 = buffer.data(pgh + 18);
    const auto *pgh_19 = buffer.data(pgh + 19);
    const auto *pgh_20 = buffer.data(pgh + 20);
    const auto *pgh_21 = buffer.data(pgh + 21);
    const auto *pgh_23 = buffer.data(pgh + 23);
    const auto *pgh_24 = buffer.data(pgh + 24);
    const auto *pgh_26 = buffer.data(pgh + 26);
    const auto *pgh_27 = buffer.data(pgh + 27);
    const auto *pgh_30 = buffer.data(pgh + 30);
    const auto *pgh_31 = buffer.data(pgh + 31);
    const auto *pgh_35 = buffer.data(pgh + 35);
    const auto *pgh_36 = buffer.data(pgh + 36);
    const auto *pgh_38 = buffer.data(pgh + 38);
    const auto *pgh_39 = buffer.data(pgh + 39);
    const auto *pgh_40 = buffer.data(pgh + 40);
    const auto *pgh_41 = buffer.data(pgh + 41);
    const auto *pgh_42 = buffer.data(pgh + 42);
    const auto *pgh_44 = buffer.data(pgh + 44);
    const auto *pgh_45 = buffer.data(pgh + 45);
    const auto *pgh_47 = buffer.data(pgh + 47);
    const auto *pgh_48 = buffer.data(pgh + 48);
    const auto *pgh_50 = buffer.data(pgh + 50);
    const auto *pgh_51 = buffer.data(pgh + 51);
    const auto *pgh_52 = buffer.data(pgh + 52);
    const auto *pgh_56 = buffer.data(pgh + 56);
    const auto *pgh_57 = buffer.data(pgh + 57);
    const auto *pgh_59 = buffer.data(pgh + 59);
    const auto *pgh_60 = buffer.data(pgh + 60);
    const auto *pgh_61 = buffer.data(pgh + 61);
    const auto *pgh_62 = buffer.data(pgh + 62);
    const auto *pgh_63 = buffer.data(pgh + 63);
    const auto *pgh_64 = buffer.data(pgh + 64);
    const auto *pgh_65 = buffer.data(pgh + 65);
    const auto *pgh_66 = buffer.data(pgh + 66);
    const auto *pgh_68 = buffer.data(pgh + 68);
    const auto *pgh_69 = buffer.data(pgh + 69);
    const auto *pgh_71 = buffer.data(pgh + 71);
    const auto *pgh_72 = buffer.data(pgh + 72);
    const auto *pgh_73 = buffer.data(pgh + 73);
    const auto *pgh_77 = buffer.data(pgh + 77);
    const auto *pgh_78 = buffer.data(pgh + 78);
    const auto *pgh_80 = buffer.data(pgh + 80);
    const auto *pgh_81 = buffer.data(pgh + 81);
    const auto *pgh_82 = buffer.data(pgh + 82);
    const auto *pgh_83 = buffer.data(pgh + 83);
    const auto *pgh_84 = buffer.data(pgh + 84);
    const auto *pgh_86 = buffer.data(pgh + 86);
    const auto *pgh_87 = buffer.data(pgh + 87);
    const auto *pgh_89 = buffer.data(pgh + 89);
    const auto *pgh_90 = buffer.data(pgh + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sgh_0, pfh_0, pgg0_0, \
                         pgg1_0, pgh_0, pgh_1, pgh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgh_0[k]
                 + f_1 * pfh_0[k]
                 + f_2 * pgg0_0[k]
                 - f_3 * pgg1_0[k]
                 + f_4 * pc_x[k] * pgh_0[k];

        t_1[k] = f_4 * pc_y[k] * pgh_0[k];

        t_2[k] = f_4 * pc_z[k] * pgh_0[k];

        t_3[k] = f_5 * pgg0_0[k]
                 - f_6 * pgg1_0[k]
                 + f_4 * pc_y[k] * pgh_1[k];

        t_4[k] = f_4 * pc_y[k] * pgh_2[k];

        t_5[k] = f_5 * pgg0_0[k]
                 - f_6 * pgg1_0[k]
                 + f_4 * pc_z[k] * pgh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, pgg0_1, pgg0_2, pgg0_3, pgg1_1, \
                         pgg1_2, pgg1_3, pgh_3, pgh_5, pgh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * pgg0_1[k]
                 - f_8 * pgg1_1[k]
                 + f_4 * pc_y[k] * pgh_3[k];

        t_7[k] = f_4 * pc_z[k] * pgh_3[k];

        t_8[k] = f_4 * pc_y[k] * pgh_5[k];

        t_9[k] = f_7 * pgg0_2[k]
                 - f_8 * pgg1_2[k]
                 + f_4 * pc_z[k] * pgh_5[k];

        t_10[k] = f_9 * pgg0_3[k]
                  - f_10 * pgg1_3[k]
                  + f_4 * pc_y[k] * pgh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, sgh_15, pfh_15, \
                         pgg0_5, pgg1_5, pgh_6, pgh_8, pgh_9, pgh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_4 * pc_z[k] * pgh_6[k];

        t_12[k] = f_5 * pgg0_5[k]
                  - f_6 * pgg1_5[k]
                  + f_4 * pc_y[k] * pgh_8[k];

        t_13[k] = f_4 * pc_y[k] * pgh_9[k];

        t_14[k] = f_9 * pgg0_5[k]
                  - f_10 * pgg1_5[k]
                  + f_4 * pc_z[k] * pgh_9[k];

        t_15[k] = f_0 * sgh_15[k]
                  + f_1 * pfh_15[k]
                  + f_4 * pc_x[k] * pgh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pc_x, pc_y, pc_z, sgh_17, sgh_18, pfh_17, \
                         pfh_18, pgh_10, pgh_14, pgh_17, pgh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * pc_z[k] * pgh_10[k];

        t_17[k] = f_0 * sgh_17[k]
                  + f_1 * pfh_17[k]
                  + f_4 * pc_x[k] * pgh_17[k];

        t_18[k] = f_0 * sgh_18[k]
                  + f_1 * pfh_18[k]
                  + f_4 * pc_x[k] * pgh_18[k];

        t_19[k] = f_4 * pc_y[k] * pgh_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, sgh_20, pfh_20, pgg0_10, \
                         pgg0_12, pgg1_10, pgg1_12, pgh_15, pgh_17, \
                         pgh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * sgh_20[k]
                  + f_1 * pfh_20[k]
                  + f_4 * pc_x[k] * pgh_20[k];

        t_21[k] = f_2 * pgg0_10[k]
                  - f_3 * pgg1_10[k]
                  + f_4 * pc_y[k] * pgh_15[k];

        t_22[k] = f_4 * pc_z[k] * pgh_15[k];

        t_23[k] = f_9 * pgg0_12[k]
                  - f_10 * pgg1_12[k]
                  + f_4 * pc_y[k] * pgh_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pc_y, pc_z, pgg0_13, pgg0_14, pgg1_13, \
                         pgg1_14, pgh_18, pgh_19, pgh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * pgg0_13[k]
                  - f_8 * pgg1_13[k]
                  + f_4 * pc_y[k] * pgh_18[k];

        t_25[k] = f_5 * pgg0_14[k]
                  - f_6 * pgg1_14[k]
                  + f_4 * pc_y[k] * pgh_19[k];

        t_26[k] = f_4 * pc_y[k] * pgh_20[k];

        t_27[k] = f_2 * pgg0_14[k]
                  - f_3 * pgg1_14[k]
                  + f_4 * pc_z[k] * pgh_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_y, pc_z, pfi0_0, pfi0_3, pfh_0, \
                         pfh_1, pfi1_0, pfi1_3, pgh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * pfi0_0[k]
                  - f_11 * pc_y[k] * pfi1_0[k];

        t_29[k] = f_0 * pfh_0[k]
                  + f_4 * pc_y[k] * pgh_21[k];

        t_30[k] = f_4 * pc_z[k] * pgh_21[k];

        t_31[k] = pb_y[k] * pfi0_3[k]
                  + f_12 * pfh_1[k]
                  - f_11 * pc_y[k] * pfi1_3[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_y, pc_y, pc_z, pfi0_5, pfi0_6, pfh_2, \
                         pfh_3, pfi1_5, pfi1_6, pgh_23, pgh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pfh_2[k]
                  + f_4 * pc_y[k] * pgh_23[k];

        t_33[k] = pb_y[k] * pfi0_5[k]
                  - f_11 * pc_y[k] * pfi1_5[k];

        t_34[k] = pb_y[k] * pfi0_6[k]
                  + f_13 * pfh_3[k]
                  - f_11 * pc_y[k] * pfi1_6[k];

        t_35[k] = f_4 * pc_z[k] * pgh_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, pfi0_9, pfi0_10, pfh_5, \
                         pfh_6, pfi1_9, pfi1_10, pgh_26, pgh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * pfh_5[k]
                  + f_4 * pc_y[k] * pgh_26[k];

        t_37[k] = pb_y[k] * pfi0_9[k]
                  - f_11 * pc_y[k] * pfi1_9[k];

        t_38[k] = pb_y[k] * pfi0_10[k]
                  + f_1 * pfh_6[k]
                  - f_11 * pc_y[k] * pfi1_10[k];

        t_39[k] = f_4 * pc_z[k] * pgh_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_y, pc_y, pfi0_12, pfi0_14, pfh_8, pfh_9, \
                         pfi1_12, pfi1_14, pgh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * pfi0_12[k]
                  + f_12 * pfh_8[k]
                  - f_11 * pc_y[k] * pfi1_12[k];

        t_41[k] = f_0 * pfh_9[k]
                  + f_4 * pc_y[k] * pgh_30[k];

        t_42[k] = pb_y[k] * pfi0_14[k]
                  - f_11 * pc_y[k] * pfi1_14[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pc_x, pc_z, sgh_36, sgh_38, sgh_39, pfh_36, \
                         pfh_38, pfh_39, pgh_31, pgh_36, pgh_38, \
                         pgh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * sgh_36[k]
                  + f_13 * pfh_36[k]
                  + f_4 * pc_x[k] * pgh_36[k];

        t_44[k] = f_4 * pc_z[k] * pgh_31[k];

        t_45[k] = f_0 * sgh_38[k]
                  + f_13 * pfh_38[k]
                  + f_4 * pc_x[k] * pgh_38[k];

        t_46[k] = f_0 * sgh_39[k]
                  + f_13 * pfh_39[k]
                  + f_4 * pc_x[k] * pgh_39[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_y, pc_y, pc_z, pfi0_20, pfh_14, pfh_15, \
                         pfi1_20, pgg0_25, pgg1_25, pgh_35, pgh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * pfh_14[k]
                  + f_4 * pc_y[k] * pgh_35[k];

        t_48[k] = pb_y[k] * pfi0_20[k]
                  - f_11 * pc_y[k] * pfi1_20[k];

        t_49[k] = f_0 * pfh_15[k]
                  + f_2 * pgg0_25[k]
                  - f_3 * pgg1_25[k]
                  + f_4 * pc_y[k] * pgh_36[k];

        t_50[k] = f_4 * pc_z[k] * pgh_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pc_y, pfh_17, pfh_18, pfh_19, pgg0_27, pgg0_28, \
                         pgg0_29, pgg1_27, pgg1_28, pgg1_29, pgh_38, pgh_39, \
                         pgh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * pfh_17[k]
                  + f_9 * pgg0_27[k]
                  - f_10 * pgg1_27[k]
                  + f_4 * pc_y[k] * pgh_38[k];

        t_52[k] = f_0 * pfh_18[k]
                  + f_7 * pgg0_28[k]
                  - f_8 * pgg1_28[k]
                  + f_4 * pc_y[k] * pgh_39[k];

        t_53[k] = f_0 * pfh_19[k]
                  + f_5 * pgg0_29[k]
                  - f_6 * pgg1_29[k]
                  + f_4 * pc_y[k] * pgh_40[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_y, pb_z, pc_y, pc_z, pfi0_0, pfi0_27, \
                         pfh_20, pfi1_0, pfi1_27, pgh_41, pgh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_0 * pfh_20[k]
                  + f_4 * pc_y[k] * pgh_41[k];

        t_55[k] = pb_y[k] * pfi0_27[k]
                  - f_11 * pc_y[k] * pfi1_27[k];

        t_56[k] = pb_z[k] * pfi0_0[k]
                  - f_11 * pc_z[k] * pfi1_0[k];

        t_57[k] = f_4 * pc_y[k] * pgh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_z, pc_y, pc_z, pfi0_3, pfi0_5, pfh_0, \
                         pfh_2, pfi1_3, pfi1_5, pgh_42, pgh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_0 * pfh_0[k]
                  + f_4 * pc_z[k] * pgh_42[k];

        t_59[k] = pb_z[k] * pfi0_3[k]
                  - f_11 * pc_z[k] * pfi1_3[k];

        t_60[k] = f_4 * pc_y[k] * pgh_44[k];

        t_61[k] = pb_z[k] * pfi0_5[k]
                  + f_12 * pfh_2[k]
                  - f_11 * pc_z[k] * pfi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_z, pc_y, pc_z, pfi0_6, pfi0_9, pfh_3, \
                         pfh_5, pfi1_6, pfi1_9, pgh_45, pgh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_z[k] * pfi0_6[k]
                  - f_11 * pc_z[k] * pfi1_6[k];

        t_63[k] = f_0 * pfh_3[k]
                  + f_4 * pc_z[k] * pgh_45[k];

        t_64[k] = f_4 * pc_y[k] * pgh_47[k];

        t_65[k] = pb_z[k] * pfi0_9[k]
                  + f_13 * pfh_5[k]
                  - f_11 * pc_z[k] * pfi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_z, pc_y, pc_z, pfi0_10, pfh_6, pfi1_10, \
                         pgg0_35, pgg1_35, pgh_48, pgh_50, pgh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_z[k] * pfi0_10[k]
                  - f_11 * pc_z[k] * pfi1_10[k];

        t_67[k] = f_0 * pfh_6[k]
                  + f_4 * pc_z[k] * pgh_48[k];

        t_68[k] = f_5 * pgg0_35[k]
                  - f_6 * pgg1_35[k]
                  + f_4 * pc_y[k] * pgh_50[k];

        t_69[k] = f_4 * pc_y[k] * pgh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_z, pc_z, pfi0_14, pfi0_15, pfh_9, pfh_10, \
                         pfi1_14, pfi1_15, pgh_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_z[k] * pfi0_14[k]
                  + f_1 * pfh_9[k]
                  - f_11 * pc_z[k] * pfi1_14[k];

        t_71[k] = pb_z[k] * pfi0_15[k]
                  - f_11 * pc_z[k] * pfi1_15[k];

        t_72[k] = f_0 * pfh_10[k]
                  + f_4 * pc_z[k] * pgh_52[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, pc_y, sgh_59, sgh_60, sgh_62, pfh_59, \
                         pfh_60, pfh_62, pgh_56, pgh_59, pgh_60, \
                         pgh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_0 * sgh_59[k]
                  + f_13 * pfh_59[k]
                  + f_4 * pc_x[k] * pgh_59[k];

        t_74[k] = f_0 * sgh_60[k]
                  + f_13 * pfh_60[k]
                  + f_4 * pc_x[k] * pgh_60[k];

        t_75[k] = f_4 * pc_y[k] * pgh_56[k];

        t_76[k] = f_0 * sgh_62[k]
                  + f_13 * pfh_62[k]
                  + f_4 * pc_x[k] * pgh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, pfi0_21, pfh_15, pfi1_21, \
                         pgg0_42, pgg1_42, pgh_57, pgh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * pfi0_21[k]
                  - f_11 * pc_z[k] * pfi1_21[k];

        t_78[k] = f_0 * pfh_15[k]
                  + f_4 * pc_z[k] * pgh_57[k];

        t_79[k] = f_9 * pgg0_42[k]
                  - f_10 * pgg1_42[k]
                  + f_4 * pc_y[k] * pgh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, pfh_20, pgg0_43, pgg0_44, \
                         pgg1_43, pgg1_44, pgh_60, pgh_61, pgh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_7 * pgg0_43[k]
                  - f_8 * pgg1_43[k]
                  + f_4 * pc_y[k] * pgh_60[k];

        t_81[k] = f_5 * pgg0_44[k]
                  - f_6 * pgg1_44[k]
                  + f_4 * pc_y[k] * pgh_61[k];

        t_82[k] = f_4 * pc_y[k] * pgh_62[k];

        t_83[k] = f_0 * pfh_20[k]
                  + f_2 * pgg0_44[k]
                  - f_3 * pgg1_44[k]
                  + f_4 * pc_z[k] * pgh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_y, pc_y, pc_z, pdi0_0, pdi1_0, pfi0_28, pfh_21, \
                         pfi1_28, pgh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_14 * pdi0_0[k]
                  - f_15 * pdi1_0[k]
                  + pb_y[k] * pfi0_28[k]
                  - f_11 * pc_y[k] * pfi1_28[k];

        t_85[k] = f_12 * pfh_21[k]
                  + f_4 * pc_y[k] * pgh_63[k];

        t_86[k] = f_4 * pc_z[k] * pgh_63[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_y, pc_z, pfh_22, pfh_23, pfh_24, pgg0_45, \
                         pgg0_46, pgg1_45, pgg1_46, pgh_64, pgh_65, \
                         pgh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_12 * pfh_22[k]
                  + f_5 * pgg0_45[k]
                  - f_6 * pgg1_45[k]
                  + f_4 * pc_y[k] * pgh_64[k];

        t_88[k] = f_12 * pfh_23[k]
                  + f_4 * pc_y[k] * pgh_65[k];

        t_89[k] = f_5 * pgg0_45[k]
                  - f_6 * pgg1_45[k]
                  + f_4 * pc_z[k] * pgh_65[k];

        t_90[k] = f_12 * pfh_24[k]
                  + f_7 * pgg0_46[k]
                  - f_8 * pgg1_46[k]
                  + f_4 * pc_y[k] * pgh_66[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pc_y, pc_z, pfh_26, pfh_27, pgg0_47, \
                         pgg0_48, pgg1_47, pgg1_48, pgh_66, pgh_68, \
                         pgh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_4 * pc_z[k] * pgh_66[k];

        t_92[k] = f_12 * pfh_26[k]
                  + f_4 * pc_y[k] * pgh_68[k];

        t_93[k] = f_7 * pgg0_47[k]
                  - f_8 * pgg1_47[k]
                  + f_4 * pc_z[k] * pgh_68[k];

        t_94[k] = f_12 * pfh_27[k]
                  + f_9 * pgg0_48[k]
                  - f_10 * pgg1_48[k]
                  + f_4 * pc_y[k] * pgh_69[k];

        t_95[k] = f_4 * pc_z[k] * pgh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, sgh_78, pfh_29, pfh_30, \
                         pfh_78, pgg0_50, pgg1_50, pgh_71, pgh_72, \
                         pgh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_12 * pfh_29[k]
                  + f_5 * pgg0_50[k]
                  - f_6 * pgg1_50[k]
                  + f_4 * pc_y[k] * pgh_71[k];

        t_97[k] = f_12 * pfh_30[k]
                  + f_4 * pc_y[k] * pgh_72[k];

        t_98[k] = f_9 * pgg0_50[k]
                  - f_10 * pgg1_50[k]
                  + f_4 * pc_z[k] * pgh_72[k];

        t_99[k] = f_0 * sgh_78[k]
                  + f_12 * pfh_78[k]
                  + f_4 * pc_x[k] * pgh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, sgh_80, sgh_81, pfh_35, \
                         pfh_80, pfh_81, pgh_73, pgh_77, pgh_80, \
                         pgh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_4 * pc_z[k] * pgh_73[k];

        t_101[k] = f_0 * sgh_80[k]
                   + f_12 * pfh_80[k]
                   + f_4 * pc_x[k] * pgh_80[k];

        t_102[k] = f_0 * sgh_81[k]
                   + f_12 * pfh_81[k]
                   + f_4 * pc_x[k] * pgh_81[k];

        t_103[k] = f_12 * pfh_35[k]
                   + f_4 * pc_y[k] * pgh_77[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pc_x, pc_y, pc_z, sgh_83, pfh_36, pfh_83, \
                         pgg0_55, pgg1_55, pgh_78, pgh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_0 * sgh_83[k]
                   + f_12 * pfh_83[k]
                   + f_4 * pc_x[k] * pgh_83[k];

        t_105[k] = f_12 * pfh_36[k]
                   + f_2 * pgg0_55[k]
                   - f_3 * pgg1_55[k]
                   + f_4 * pc_y[k] * pgh_78[k];

        t_106[k] = f_4 * pc_z[k] * pgh_78[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, pfh_38, pfh_39, pfh_40, pgg0_57, pgg0_58, \
                         pgg0_59, pgg1_57, pgg1_58, pgg1_59, pgh_80, pgh_81, \
                         pgh_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_12 * pfh_38[k]
                   + f_9 * pgg0_57[k]
                   - f_10 * pgg1_57[k]
                   + f_4 * pc_y[k] * pgh_80[k];

        t_108[k] = f_12 * pfh_39[k]
                   + f_7 * pgg0_58[k]
                   - f_8 * pgg1_58[k]
                   + f_4 * pc_y[k] * pgh_81[k];

        t_109[k] = f_12 * pfh_40[k]
                   + f_5 * pgg0_59[k]
                   - f_6 * pgg1_59[k]
                   + f_4 * pc_y[k] * pgh_82[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_y, pc_y, pc_z, pfi0_56, pfh_41, \
                         pfh_42, pfi1_56, pgg0_59, pgg1_59, pgh_83, \
                         pgh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_12 * pfh_41[k]
                   + f_4 * pc_y[k] * pgh_83[k];

        t_111[k] = f_2 * pgg0_59[k]
                   - f_3 * pgg1_59[k]
                   + f_4 * pc_z[k] * pgh_83[k];

        t_112[k] = pb_y[k] * pfi0_56[k]
                   - f_11 * pc_y[k] * pfi1_56[k];

        t_113[k] = f_0 * pfh_42[k]
                   + f_4 * pc_y[k] * pgh_84[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_y, pb_z, pc_y, pc_z, pfi0_31, pfi0_61, \
                         pfh_21, pfh_44, pfi1_31, pfi1_61, pgh_84, \
                         pgh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * pfh_21[k]
                   + f_4 * pc_z[k] * pgh_84[k];

        t_115[k] = pb_z[k] * pfi0_31[k]
                   - f_11 * pc_z[k] * pfi1_31[k];

        t_116[k] = f_0 * pfh_44[k]
                   + f_4 * pc_y[k] * pgh_86[k];

        t_117[k] = pb_y[k] * pfi0_61[k]
                   - f_11 * pc_y[k] * pfi1_61[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pb_y, pb_z, pc_y, pc_z, pfi0_34, pfi0_65, \
                         pfh_24, pfh_47, pfi1_34, pfi1_65, pgh_87, \
                         pgh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_z[k] * pfi0_34[k]
                   - f_11 * pc_z[k] * pfi1_34[k];

        t_119[k] = f_0 * pfh_24[k]
                   + f_4 * pc_z[k] * pgh_87[k];

        t_120[k] = f_0 * pfh_47[k]
                   + f_4 * pc_y[k] * pgh_89[k];

        t_121[k] = pb_y[k] * pfi0_65[k]
                   - f_11 * pc_y[k] * pfi1_65[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pb_y, pb_z, pc_y, pc_z, pfi0_38, pfi0_68, \
                         pfh_27, pfh_50, pfi1_38, pfi1_68, pgh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = pb_z[k] * pfi0_38[k]
                   - f_11 * pc_z[k] * pfi1_38[k];

        t_123[k] = f_0 * pfh_27[k]
                   + f_4 * pc_z[k] * pgh_90[k];

        t_124[k] = pb_y[k] * pfi0_68[k]
                   + f_12 * pfh_50[k]
                   - f_11 * pc_y[k] * pfi1_68[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgh,
                                                          const size_t pdi0, const size_t pdi1,
                                                          const size_t pfi0, const size_t pfh,
                                                          const size_t pfi1, const size_t pgg0,
                                                          const size_t pgg1, const size_t pgh,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh_101 = buffer.data(sgh + 101);
    const auto *sgh_102 = buffer.data(sgh + 102);
    const auto *sgh_120 = buffer.data(sgh + 120);
    const auto *sgh_122 = buffer.data(sgh + 122);
    const auto *sgh_123 = buffer.data(sgh + 123);
    const auto *sgh_125 = buffer.data(sgh + 125);
    const auto *sgh_126 = buffer.data(sgh + 126);
    const auto *sgh_141 = buffer.data(sgh + 141);
    const auto *sgh_143 = buffer.data(sgh + 143);
    const auto *sgh_144 = buffer.data(sgh + 144);
    const auto *sgh_146 = buffer.data(sgh + 146);
    const auto *sgh_164 = buffer.data(sgh + 164);
    const auto *sgh_165 = buffer.data(sgh + 165);
    const auto *sgh_167 = buffer.data(sgh + 167);

    const auto *pdi0_0 = buffer.data(pdi0 + 0);

    const auto *pdi1_0 = buffer.data(pdi1 + 0);

    const auto *pfi0_43 = buffer.data(pfi0 + 43);
    const auto *pfi0_49 = buffer.data(pfi0 + 49);
    const auto *pfi0_56 = buffer.data(pfi0 + 56);
    const auto *pfi0_70 = buffer.data(pfi0 + 70);
    const auto *pfi0_76 = buffer.data(pfi0 + 76);
    const auto *pfi0_83 = buffer.data(pfi0 + 83);
    const auto *pfi0_84 = buffer.data(pfi0 + 84);
    const auto *pfi0_87 = buffer.data(pfi0 + 87);
    const auto *pfi0_90 = buffer.data(pfi0 + 90);
    const auto *pfi0_94 = buffer.data(pfi0 + 94);
    const auto *pfi0_99 = buffer.data(pfi0 + 99);
    const auto *pfi0_140 = buffer.data(pfi0 + 140);
    const auto *pfi0_145 = buffer.data(pfi0 + 145);
    const auto *pfi0_149 = buffer.data(pfi0 + 149);
    const auto *pfi0_154 = buffer.data(pfi0 + 154);

    const auto *pfh_31 = buffer.data(pfh + 31);
    const auto *pfh_36 = buffer.data(pfh + 36);
    const auto *pfh_42 = buffer.data(pfh + 42);
    const auto *pfh_44 = buffer.data(pfh + 44);
    const auto *pfh_45 = buffer.data(pfh + 45);
    const auto *pfh_47 = buffer.data(pfh + 47);
    const auto *pfh_48 = buffer.data(pfh + 48);
    const auto *pfh_51 = buffer.data(pfh + 51);
    const auto *pfh_52 = buffer.data(pfh + 52);
    const auto *pfh_56 = buffer.data(pfh + 56);
    const auto *pfh_57 = buffer.data(pfh + 57);
    const auto *pfh_59 = buffer.data(pfh + 59);
    const auto *pfh_60 = buffer.data(pfh + 60);
    const auto *pfh_61 = buffer.data(pfh + 61);
    const auto *pfh_62 = buffer.data(pfh + 62);
    const auto *pfh_63 = buffer.data(pfh + 63);
    const auto *pfh_64 = buffer.data(pfh + 64);
    const auto *pfh_65 = buffer.data(pfh + 65);
    const auto *pfh_66 = buffer.data(pfh + 66);
    const auto *pfh_68 = buffer.data(pfh + 68);
    const auto *pfh_69 = buffer.data(pfh + 69);
    const auto *pfh_71 = buffer.data(pfh + 71);
    const auto *pfh_72 = buffer.data(pfh + 72);
    const auto *pfh_73 = buffer.data(pfh + 73);
    const auto *pfh_77 = buffer.data(pfh + 77);
    const auto *pfh_78 = buffer.data(pfh + 78);
    const auto *pfh_80 = buffer.data(pfh + 80);
    const auto *pfh_81 = buffer.data(pfh + 81);
    const auto *pfh_82 = buffer.data(pfh + 82);
    const auto *pfh_83 = buffer.data(pfh + 83);
    const auto *pfh_84 = buffer.data(pfh + 84);
    const auto *pfh_86 = buffer.data(pfh + 86);
    const auto *pfh_87 = buffer.data(pfh + 87);
    const auto *pfh_89 = buffer.data(pfh + 89);
    const auto *pfh_90 = buffer.data(pfh + 90);
    const auto *pfh_92 = buffer.data(pfh + 92);
    const auto *pfh_93 = buffer.data(pfh + 93);
    const auto *pfh_98 = buffer.data(pfh + 98);
    const auto *pfh_99 = buffer.data(pfh + 99);
    const auto *pfh_101 = buffer.data(pfh + 101);
    const auto *pfh_102 = buffer.data(pfh + 102);
    const auto *pfh_103 = buffer.data(pfh + 103);
    const auto *pfh_104 = buffer.data(pfh + 104);
    const auto *pfh_105 = buffer.data(pfh + 105);
    const auto *pfh_106 = buffer.data(pfh + 106);
    const auto *pfh_107 = buffer.data(pfh + 107);
    const auto *pfh_108 = buffer.data(pfh + 108);
    const auto *pfh_110 = buffer.data(pfh + 110);
    const auto *pfh_111 = buffer.data(pfh + 111);
    const auto *pfh_113 = buffer.data(pfh + 113);
    const auto *pfh_114 = buffer.data(pfh + 114);
    const auto *pfh_120 = buffer.data(pfh + 120);
    const auto *pfh_122 = buffer.data(pfh + 122);
    const auto *pfh_123 = buffer.data(pfh + 123);
    const auto *pfh_125 = buffer.data(pfh + 125);
    const auto *pfh_126 = buffer.data(pfh + 126);
    const auto *pfh_141 = buffer.data(pfh + 141);
    const auto *pfh_143 = buffer.data(pfh + 143);
    const auto *pfh_144 = buffer.data(pfh + 144);
    const auto *pfh_146 = buffer.data(pfh + 146);
    const auto *pfh_164 = buffer.data(pfh + 164);
    const auto *pfh_165 = buffer.data(pfh + 165);
    const auto *pfh_167 = buffer.data(pfh + 167);

    const auto *pfi1_43 = buffer.data(pfi1 + 43);
    const auto *pfi1_49 = buffer.data(pfi1 + 49);
    const auto *pfi1_56 = buffer.data(pfi1 + 56);
    const auto *pfi1_70 = buffer.data(pfi1 + 70);
    const auto *pfi1_76 = buffer.data(pfi1 + 76);
    const auto *pfi1_83 = buffer.data(pfi1 + 83);
    const auto *pfi1_84 = buffer.data(pfi1 + 84);
    const auto *pfi1_87 = buffer.data(pfi1 + 87);
    const auto *pfi1_90 = buffer.data(pfi1 + 90);
    const auto *pfi1_94 = buffer.data(pfi1 + 94);
    const auto *pfi1_99 = buffer.data(pfi1 + 99);
    const auto *pfi1_140 = buffer.data(pfi1 + 140);
    const auto *pfi1_145 = buffer.data(pfi1 + 145);
    const auto *pfi1_149 = buffer.data(pfi1 + 149);
    const auto *pfi1_154 = buffer.data(pfi1 + 154);

    const auto *pgg0_72 = buffer.data(pgg0 + 72);
    const auto *pgg0_73 = buffer.data(pgg0 + 73);
    const auto *pgg0_74 = buffer.data(pgg0 + 74);
    const auto *pgg0_75 = buffer.data(pgg0 + 75);
    const auto *pgg0_76 = buffer.data(pgg0 + 76);
    const auto *pgg0_77 = buffer.data(pgg0 + 77);
    const auto *pgg0_78 = buffer.data(pgg0 + 78);
    const auto *pgg0_80 = buffer.data(pgg0 + 80);
    const auto *pgg0_85 = buffer.data(pgg0 + 85);
    const auto *pgg0_87 = buffer.data(pgg0 + 87);
    const auto *pgg0_88 = buffer.data(pgg0 + 88);
    const auto *pgg0_89 = buffer.data(pgg0 + 89);
    const auto *pgg0_90 = buffer.data(pgg0 + 90);
    const auto *pgg0_91 = buffer.data(pgg0 + 91);
    const auto *pgg0_92 = buffer.data(pgg0 + 92);
    const auto *pgg0_93 = buffer.data(pgg0 + 93);
    const auto *pgg0_95 = buffer.data(pgg0 + 95);
    const auto *pgg0_100 = buffer.data(pgg0 + 100);
    const auto *pgg0_102 = buffer.data(pgg0 + 102);
    const auto *pgg0_103 = buffer.data(pgg0 + 103);
    const auto *pgg0_104 = buffer.data(pgg0 + 104);
    const auto *pgg0_105 = buffer.data(pgg0 + 105);
    const auto *pgg0_107 = buffer.data(pgg0 + 107);
    const auto *pgg0_110 = buffer.data(pgg0 + 110);
    const auto *pgg0_115 = buffer.data(pgg0 + 115);
    const auto *pgg0_117 = buffer.data(pgg0 + 117);
    const auto *pgg0_118 = buffer.data(pgg0 + 118);
    const auto *pgg0_119 = buffer.data(pgg0 + 119);
    const auto *pgg0_120 = buffer.data(pgg0 + 120);
    const auto *pgg0_121 = buffer.data(pgg0 + 121);
    const auto *pgg0_123 = buffer.data(pgg0 + 123);
    const auto *pgg0_125 = buffer.data(pgg0 + 125);

    const auto *pgg1_72 = buffer.data(pgg1 + 72);
    const auto *pgg1_73 = buffer.data(pgg1 + 73);
    const auto *pgg1_74 = buffer.data(pgg1 + 74);
    const auto *pgg1_75 = buffer.data(pgg1 + 75);
    const auto *pgg1_76 = buffer.data(pgg1 + 76);
    const auto *pgg1_77 = buffer.data(pgg1 + 77);
    const auto *pgg1_78 = buffer.data(pgg1 + 78);
    const auto *pgg1_80 = buffer.data(pgg1 + 80);
    const auto *pgg1_85 = buffer.data(pgg1 + 85);
    const auto *pgg1_87 = buffer.data(pgg1 + 87);
    const auto *pgg1_88 = buffer.data(pgg1 + 88);
    const auto *pgg1_89 = buffer.data(pgg1 + 89);
    const auto *pgg1_90 = buffer.data(pgg1 + 90);
    const auto *pgg1_91 = buffer.data(pgg1 + 91);
    const auto *pgg1_92 = buffer.data(pgg1 + 92);
    const auto *pgg1_93 = buffer.data(pgg1 + 93);
    const auto *pgg1_95 = buffer.data(pgg1 + 95);
    const auto *pgg1_100 = buffer.data(pgg1 + 100);
    const auto *pgg1_102 = buffer.data(pgg1 + 102);
    const auto *pgg1_103 = buffer.data(pgg1 + 103);
    const auto *pgg1_104 = buffer.data(pgg1 + 104);
    const auto *pgg1_105 = buffer.data(pgg1 + 105);
    const auto *pgg1_107 = buffer.data(pgg1 + 107);
    const auto *pgg1_110 = buffer.data(pgg1 + 110);
    const auto *pgg1_115 = buffer.data(pgg1 + 115);
    const auto *pgg1_117 = buffer.data(pgg1 + 117);
    const auto *pgg1_118 = buffer.data(pgg1 + 118);
    const auto *pgg1_119 = buffer.data(pgg1 + 119);
    const auto *pgg1_120 = buffer.data(pgg1 + 120);
    const auto *pgg1_121 = buffer.data(pgg1 + 121);
    const auto *pgg1_123 = buffer.data(pgg1 + 123);
    const auto *pgg1_125 = buffer.data(pgg1 + 125);

    const auto *pgh_93 = buffer.data(pgh + 93);
    const auto *pgh_94 = buffer.data(pgh + 94);
    const auto *pgh_98 = buffer.data(pgh + 98);
    const auto *pgh_99 = buffer.data(pgh + 99);
    const auto *pgh_101 = buffer.data(pgh + 101);
    const auto *pgh_102 = buffer.data(pgh + 102);
    const auto *pgh_103 = buffer.data(pgh + 103);
    const auto *pgh_104 = buffer.data(pgh + 104);
    const auto *pgh_105 = buffer.data(pgh + 105);
    const auto *pgh_106 = buffer.data(pgh + 106);
    const auto *pgh_107 = buffer.data(pgh + 107);
    const auto *pgh_108 = buffer.data(pgh + 108);
    const auto *pgh_110 = buffer.data(pgh + 110);
    const auto *pgh_111 = buffer.data(pgh + 111);
    const auto *pgh_113 = buffer.data(pgh + 113);
    const auto *pgh_114 = buffer.data(pgh + 114);
    const auto *pgh_115 = buffer.data(pgh + 115);
    const auto *pgh_119 = buffer.data(pgh + 119);
    const auto *pgh_120 = buffer.data(pgh + 120);
    const auto *pgh_122 = buffer.data(pgh + 122);
    const auto *pgh_123 = buffer.data(pgh + 123);
    const auto *pgh_124 = buffer.data(pgh + 124);
    const auto *pgh_125 = buffer.data(pgh + 125);
    const auto *pgh_126 = buffer.data(pgh + 126);
    const auto *pgh_127 = buffer.data(pgh + 127);
    const auto *pgh_128 = buffer.data(pgh + 128);
    const auto *pgh_129 = buffer.data(pgh + 129);
    const auto *pgh_131 = buffer.data(pgh + 131);
    const auto *pgh_132 = buffer.data(pgh + 132);
    const auto *pgh_134 = buffer.data(pgh + 134);
    const auto *pgh_135 = buffer.data(pgh + 135);
    const auto *pgh_136 = buffer.data(pgh + 136);
    const auto *pgh_140 = buffer.data(pgh + 140);
    const auto *pgh_141 = buffer.data(pgh + 141);
    const auto *pgh_143 = buffer.data(pgh + 143);
    const auto *pgh_144 = buffer.data(pgh + 144);
    const auto *pgh_145 = buffer.data(pgh + 145);
    const auto *pgh_146 = buffer.data(pgh + 146);
    const auto *pgh_147 = buffer.data(pgh + 147);
    const auto *pgh_149 = buffer.data(pgh + 149);
    const auto *pgh_150 = buffer.data(pgh + 150);
    const auto *pgh_152 = buffer.data(pgh + 152);
    const auto *pgh_153 = buffer.data(pgh + 153);
    const auto *pgh_155 = buffer.data(pgh + 155);
    const auto *pgh_156 = buffer.data(pgh + 156);
    const auto *pgh_157 = buffer.data(pgh + 157);
    const auto *pgh_161 = buffer.data(pgh + 161);
    const auto *pgh_162 = buffer.data(pgh + 162);
    const auto *pgh_164 = buffer.data(pgh + 164);
    const auto *pgh_165 = buffer.data(pgh + 165);
    const auto *pgh_166 = buffer.data(pgh + 166);
    const auto *pgh_167 = buffer.data(pgh + 167);
    const auto *pgh_168 = buffer.data(pgh + 168);
    const auto *pgh_169 = buffer.data(pgh + 169);
    const auto *pgh_170 = buffer.data(pgh + 170);
    const auto *pgh_171 = buffer.data(pgh + 171);
    const auto *pgh_173 = buffer.data(pgh + 173);
    const auto *pgh_174 = buffer.data(pgh + 174);
    const auto *pgh_176 = buffer.data(pgh + 176);
    const auto *pgh_177 = buffer.data(pgh + 177);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pb_z, pc_y, pc_z, pfi0_43, pfi0_70, \
                         pfh_31, pfh_51, pfi1_43, pfi1_70, pgh_93, \
                         pgh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * pfh_51[k]
                   + f_4 * pc_y[k] * pgh_93[k];

        t_126[k] = pb_y[k] * pfi0_70[k]
                   - f_11 * pc_y[k] * pfi1_70[k];

        t_127[k] = pb_z[k] * pfi0_43[k]
                   - f_11 * pc_z[k] * pfi1_43[k];

        t_128[k] = f_0 * pfh_31[k]
                   + f_4 * pc_z[k] * pgh_94[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, sgh_101, sgh_102, pfh_56, pfh_101, \
                         pfh_102, pgh_98, pgh_101, pgh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_0 * sgh_101[k]
                   + f_12 * pfh_101[k]
                   + f_4 * pc_x[k] * pgh_101[k];

        t_130[k] = f_0 * sgh_102[k]
                   + f_12 * pfh_102[k]
                   + f_4 * pc_x[k] * pgh_102[k];

        t_131[k] = f_0 * pfh_56[k]
                   + f_4 * pc_y[k] * pgh_98[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pb_z, pc_y, pc_z, pfi0_49, pfi0_76, \
                         pfh_36, pfi1_49, pfi1_76, pgh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_y[k] * pfi0_76[k]
                   - f_11 * pc_y[k] * pfi1_76[k];

        t_133[k] = pb_z[k] * pfi0_49[k]
                   - f_11 * pc_z[k] * pfi1_49[k];

        t_134[k] = f_0 * pfh_36[k]
                   + f_4 * pc_z[k] * pgh_99[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pc_y, pfh_59, pfh_60, pfh_61, pgg0_72, pgg0_73, \
                         pgg0_74, pgg1_72, pgg1_73, pgg1_74, pgh_101, pgh_102, \
                         pgh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * pfh_59[k]
                   + f_9 * pgg0_72[k]
                   - f_10 * pgg1_72[k]
                   + f_4 * pc_y[k] * pgh_101[k];

        t_136[k] = f_0 * pfh_60[k]
                   + f_7 * pgg0_73[k]
                   - f_8 * pgg1_73[k]
                   + f_4 * pc_y[k] * pgh_102[k];

        t_137[k] = f_0 * pfh_61[k]
                   + f_5 * pgg0_74[k]
                   - f_6 * pgg1_74[k]
                   + f_4 * pc_y[k] * pgh_103[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_y, pb_z, pc_y, pc_z, pdi0_0, pdi1_0, pfi0_56, \
                         pfi0_83, pfh_62, pfi1_56, pfi1_83, pgh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * pfh_62[k]
                   + f_4 * pc_y[k] * pgh_104[k];

        t_139[k] = pb_y[k] * pfi0_83[k]
                   - f_11 * pc_y[k] * pfi1_83[k];

        t_140[k] = f_14 * pdi0_0[k]
                   - f_15 * pdi1_0[k]
                   + pb_z[k] * pfi0_56[k]
                   - f_11 * pc_z[k] * pfi1_56[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, pc_y, pc_z, pfh_42, pfh_44, \
                         pgg0_75, pgg1_75, pgh_105, pgh_106, pgh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_4 * pc_y[k] * pgh_105[k];

        t_142[k] = f_12 * pfh_42[k]
                   + f_4 * pc_z[k] * pgh_105[k];

        t_143[k] = f_5 * pgg0_75[k]
                   - f_6 * pgg1_75[k]
                   + f_4 * pc_y[k] * pgh_106[k];

        t_144[k] = f_4 * pc_y[k] * pgh_107[k];

        t_145[k] = f_12 * pfh_44[k]
                   + f_5 * pgg0_75[k]
                   - f_6 * pgg1_75[k]
                   + f_4 * pc_z[k] * pgh_107[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_y, pc_z, pfh_45, pfh_47, pgg0_76, \
                         pgg0_77, pgg1_76, pgg1_77, pgh_108, pgh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_7 * pgg0_76[k]
                   - f_8 * pgg1_76[k]
                   + f_4 * pc_y[k] * pgh_108[k];

        t_147[k] = f_12 * pfh_45[k]
                   + f_4 * pc_z[k] * pgh_108[k];

        t_148[k] = f_4 * pc_y[k] * pgh_110[k];

        t_149[k] = f_12 * pfh_47[k]
                   + f_7 * pgg0_77[k]
                   - f_8 * pgg1_77[k]
                   + f_4 * pc_z[k] * pgh_110[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pc_y, pc_z, pfh_48, pfh_51, \
                         pgg0_78, pgg0_80, pgg1_78, pgg1_80, pgh_111, pgh_113, \
                         pgh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_9 * pgg0_78[k]
                   - f_10 * pgg1_78[k]
                   + f_4 * pc_y[k] * pgh_111[k];

        t_151[k] = f_12 * pfh_48[k]
                   + f_4 * pc_z[k] * pgh_111[k];

        t_152[k] = f_5 * pgg0_80[k]
                   - f_6 * pgg1_80[k]
                   + f_4 * pc_y[k] * pgh_113[k];

        t_153[k] = f_4 * pc_y[k] * pgh_114[k];

        t_154[k] = f_12 * pfh_51[k]
                   + f_9 * pgg0_80[k]
                   - f_10 * pgg1_80[k]
                   + f_4 * pc_z[k] * pgh_114[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pc_x, pc_z, sgh_120, sgh_122, pfh_52, pfh_120, \
                         pfh_122, pgh_115, pgh_120, pgh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_0 * sgh_120[k]
                   + f_12 * pfh_120[k]
                   + f_4 * pc_x[k] * pgh_120[k];

        t_156[k] = f_12 * pfh_52[k]
                   + f_4 * pc_z[k] * pgh_115[k];

        t_157[k] = f_0 * sgh_122[k]
                   + f_12 * pfh_122[k]
                   + f_4 * pc_x[k] * pgh_122[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pc_x, pc_y, sgh_123, sgh_125, pfh_123, \
                         pfh_125, pgg0_85, pgg1_85, pgh_119, pgh_120, pgh_123, \
                         pgh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_0 * sgh_123[k]
                   + f_12 * pfh_123[k]
                   + f_4 * pc_x[k] * pgh_123[k];

        t_159[k] = f_4 * pc_y[k] * pgh_119[k];

        t_160[k] = f_0 * sgh_125[k]
                   + f_12 * pfh_125[k]
                   + f_4 * pc_x[k] * pgh_125[k];

        t_161[k] = f_2 * pgg0_85[k]
                   - f_3 * pgg1_85[k]
                   + f_4 * pc_y[k] * pgh_120[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pc_y, pc_z, pfh_57, pgg0_87, pgg0_88, pgg1_87, \
                         pgg1_88, pgh_120, pgh_122, pgh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_12 * pfh_57[k]
                   + f_4 * pc_z[k] * pgh_120[k];

        t_163[k] = f_9 * pgg0_87[k]
                   - f_10 * pgg1_87[k]
                   + f_4 * pc_y[k] * pgh_122[k];

        t_164[k] = f_7 * pgg0_88[k]
                   - f_8 * pgg1_88[k]
                   + f_4 * pc_y[k] * pgh_123[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pc_y, pc_z, pfh_62, pgg0_89, pgg1_89, pgh_124, \
                         pgh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_5 * pgg0_89[k]
                   - f_6 * pgg1_89[k]
                   + f_4 * pc_y[k] * pgh_124[k];

        t_166[k] = f_4 * pc_y[k] * pgh_125[k];

        t_167[k] = f_12 * pfh_62[k]
                   + f_2 * pgg0_89[k]
                   - f_3 * pgg1_89[k]
                   + f_4 * pc_z[k] * pgh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, sgh_126, pfh_63, \
                         pfh_64, pfh_126, pgg0_90, pgg1_90, pgh_126, \
                         pgh_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_0 * sgh_126[k]
                   + f_0 * pfh_126[k]
                   + f_2 * pgg0_90[k]
                   - f_3 * pgg1_90[k]
                   + f_4 * pc_x[k] * pgh_126[k];

        t_169[k] = f_13 * pfh_63[k]
                   + f_4 * pc_y[k] * pgh_126[k];

        t_170[k] = f_4 * pc_z[k] * pgh_126[k];

        t_171[k] = f_13 * pfh_64[k]
                   + f_5 * pgg0_90[k]
                   - f_6 * pgg1_90[k]
                   + f_4 * pc_y[k] * pgh_127[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_y, pc_z, pfh_65, pfh_66, pgg0_90, \
                         pgg0_91, pgg1_90, pgg1_91, pgh_128, pgh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_13 * pfh_65[k]
                   + f_4 * pc_y[k] * pgh_128[k];

        t_173[k] = f_5 * pgg0_90[k]
                   - f_6 * pgg1_90[k]
                   + f_4 * pc_z[k] * pgh_128[k];

        t_174[k] = f_13 * pfh_66[k]
                   + f_7 * pgg0_91[k]
                   - f_8 * pgg1_91[k]
                   + f_4 * pc_y[k] * pgh_129[k];

        t_175[k] = f_4 * pc_z[k] * pgh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_y, pc_z, pfh_68, pfh_69, pgg0_92, \
                         pgg0_93, pgg1_92, pgg1_93, pgh_131, pgh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * pfh_68[k]
                   + f_4 * pc_y[k] * pgh_131[k];

        t_177[k] = f_7 * pgg0_92[k]
                   - f_8 * pgg1_92[k]
                   + f_4 * pc_z[k] * pgh_131[k];

        t_178[k] = f_13 * pfh_69[k]
                   + f_9 * pgg0_93[k]
                   - f_10 * pgg1_93[k]
                   + f_4 * pc_y[k] * pgh_132[k];

        t_179[k] = f_4 * pc_z[k] * pgh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, sgh_141, pfh_71, \
                         pfh_72, pfh_141, pgg0_95, pgg1_95, pgh_134, pgh_135, \
                         pgh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_13 * pfh_71[k]
                   + f_5 * pgg0_95[k]
                   - f_6 * pgg1_95[k]
                   + f_4 * pc_y[k] * pgh_134[k];

        t_181[k] = f_13 * pfh_72[k]
                   + f_4 * pc_y[k] * pgh_135[k];

        t_182[k] = f_9 * pgg0_95[k]
                   - f_10 * pgg1_95[k]
                   + f_4 * pc_z[k] * pgh_135[k];

        t_183[k] = f_0 * sgh_141[k]
                   + f_0 * pfh_141[k]
                   + f_4 * pc_x[k] * pgh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pc_x, pc_y, pc_z, sgh_143, sgh_144, \
                         pfh_77, pfh_143, pfh_144, pgh_136, pgh_140, pgh_143, \
                         pgh_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_4 * pc_z[k] * pgh_136[k];

        t_185[k] = f_0 * sgh_143[k]
                   + f_0 * pfh_143[k]
                   + f_4 * pc_x[k] * pgh_143[k];

        t_186[k] = f_0 * sgh_144[k]
                   + f_0 * pfh_144[k]
                   + f_4 * pc_x[k] * pgh_144[k];

        t_187[k] = f_13 * pfh_77[k]
                   + f_4 * pc_y[k] * pgh_140[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pc_x, pc_y, pc_z, sgh_146, pfh_78, pfh_146, \
                         pgg0_100, pgg1_100, pgh_141, pgh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_0 * sgh_146[k]
                   + f_0 * pfh_146[k]
                   + f_4 * pc_x[k] * pgh_146[k];

        t_189[k] = f_13 * pfh_78[k]
                   + f_2 * pgg0_100[k]
                   - f_3 * pgg1_100[k]
                   + f_4 * pc_y[k] * pgh_141[k];

        t_190[k] = f_4 * pc_z[k] * pgh_141[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pc_y, pfh_80, pfh_81, pfh_82, pgg0_102, \
                         pgg0_103, pgg0_104, pgg1_102, pgg1_103, pgg1_104, pgh_143, pgh_144, \
                         pgh_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_13 * pfh_80[k]
                   + f_9 * pgg0_102[k]
                   - f_10 * pgg1_102[k]
                   + f_4 * pc_y[k] * pgh_143[k];

        t_192[k] = f_13 * pfh_81[k]
                   + f_7 * pgg0_103[k]
                   - f_8 * pgg1_103[k]
                   + f_4 * pc_y[k] * pgh_144[k];

        t_193[k] = f_13 * pfh_82[k]
                   + f_5 * pgg0_104[k]
                   - f_6 * pgg1_104[k]
                   + f_4 * pc_y[k] * pgh_145[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pb_z, pc_y, pc_z, pfi0_84, pfh_83, \
                         pfh_84, pfi1_84, pgg0_104, pgg1_104, pgh_146, \
                         pgh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_13 * pfh_83[k]
                   + f_4 * pc_y[k] * pgh_146[k];

        t_195[k] = f_2 * pgg0_104[k]
                   - f_3 * pgg1_104[k]
                   + f_4 * pc_z[k] * pgh_146[k];

        t_196[k] = pb_z[k] * pfi0_84[k]
                   - f_11 * pc_z[k] * pfi1_84[k];

        t_197[k] = f_12 * pfh_84[k]
                   + f_4 * pc_y[k] * pgh_147[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pb_z, pc_y, pc_z, pfi0_87, pfh_63, \
                         pfh_65, pfh_86, pfi1_87, pgg0_105, pgg1_105, pgh_147, \
                         pgh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_0 * pfh_63[k]
                   + f_4 * pc_z[k] * pgh_147[k];

        t_199[k] = pb_z[k] * pfi0_87[k]
                   - f_11 * pc_z[k] * pfi1_87[k];

        t_200[k] = f_12 * pfh_86[k]
                   + f_4 * pc_y[k] * pgh_149[k];

        t_201[k] = f_0 * pfh_65[k]
                   + f_5 * pgg0_105[k]
                   - f_6 * pgg1_105[k]
                   + f_4 * pc_z[k] * pgh_149[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_z, pc_y, pc_z, pfi0_90, pfh_66, \
                         pfh_68, pfh_89, pfi1_90, pgg0_107, pgg1_107, pgh_150, \
                         pgh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = pb_z[k] * pfi0_90[k]
                   - f_11 * pc_z[k] * pfi1_90[k];

        t_203[k] = f_0 * pfh_66[k]
                   + f_4 * pc_z[k] * pgh_150[k];

        t_204[k] = f_12 * pfh_89[k]
                   + f_4 * pc_y[k] * pgh_152[k];

        t_205[k] = f_0 * pfh_68[k]
                   + f_7 * pgg0_107[k]
                   - f_8 * pgg1_107[k]
                   + f_4 * pc_z[k] * pgh_152[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pb_z, pc_y, pc_z, pfi0_94, pfh_69, pfh_92, \
                         pfi1_94, pgg0_110, pgg1_110, pgh_153, \
                         pgh_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_z[k] * pfi0_94[k]
                   - f_11 * pc_z[k] * pfi1_94[k];

        t_207[k] = f_0 * pfh_69[k]
                   + f_4 * pc_z[k] * pgh_153[k];

        t_208[k] = f_12 * pfh_92[k]
                   + f_5 * pgg0_110[k]
                   - f_6 * pgg1_110[k]
                   + f_4 * pc_y[k] * pgh_155[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pb_z, pc_y, pc_z, pfi0_99, pfh_72, \
                         pfh_73, pfh_93, pfi1_99, pgg0_110, pgg1_110, pgh_156, \
                         pgh_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_12 * pfh_93[k]
                   + f_4 * pc_y[k] * pgh_156[k];

        t_210[k] = f_0 * pfh_72[k]
                   + f_9 * pgg0_110[k]
                   - f_10 * pgg1_110[k]
                   + f_4 * pc_z[k] * pgh_156[k];

        t_211[k] = pb_z[k] * pfi0_99[k]
                   - f_11 * pc_z[k] * pfi1_99[k];

        t_212[k] = f_0 * pfh_73[k]
                   + f_4 * pc_z[k] * pgh_157[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, pc_x, pc_y, sgh_164, sgh_165, pfh_98, pfh_164, \
                         pfh_165, pgh_161, pgh_164, pgh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_0 * sgh_164[k]
                   + f_0 * pfh_164[k]
                   + f_4 * pc_x[k] * pgh_164[k];

        t_214[k] = f_0 * sgh_165[k]
                   + f_0 * pfh_165[k]
                   + f_4 * pc_x[k] * pgh_165[k];

        t_215[k] = f_12 * pfh_98[k]
                   + f_4 * pc_y[k] * pgh_161[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pc_x, pc_y, pc_z, sgh_167, pfh_78, pfh_99, \
                         pfh_167, pgg0_115, pgg1_115, pgh_162, \
                         pgh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_0 * sgh_167[k]
                   + f_0 * pfh_167[k]
                   + f_4 * pc_x[k] * pgh_167[k];

        t_217[k] = f_12 * pfh_99[k]
                   + f_2 * pgg0_115[k]
                   - f_3 * pgg1_115[k]
                   + f_4 * pc_y[k] * pgh_162[k];

        t_218[k] = f_0 * pfh_78[k]
                   + f_4 * pc_z[k] * pgh_162[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pc_y, pfh_101, pfh_102, pfh_103, pgg0_117, \
                         pgg0_118, pgg0_119, pgg1_117, pgg1_118, pgg1_119, pgh_164, pgh_165, \
                         pgh_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_12 * pfh_101[k]
                   + f_9 * pgg0_117[k]
                   - f_10 * pgg1_117[k]
                   + f_4 * pc_y[k] * pgh_164[k];

        t_220[k] = f_12 * pfh_102[k]
                   + f_7 * pgg0_118[k]
                   - f_8 * pgg1_118[k]
                   + f_4 * pc_y[k] * pgh_165[k];

        t_221[k] = f_12 * pfh_103[k]
                   + f_5 * pgg0_119[k]
                   - f_6 * pgg1_119[k]
                   + f_4 * pc_y[k] * pgh_166[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pb_y, pc_y, pc_z, pfi0_140, pfh_83, \
                         pfh_104, pfh_105, pfi1_140, pgg0_119, pgg1_119, pgh_167, \
                         pgh_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_12 * pfh_104[k]
                   + f_4 * pc_y[k] * pgh_167[k];

        t_223[k] = f_0 * pfh_83[k]
                   + f_2 * pgg0_119[k]
                   - f_3 * pgg1_119[k]
                   + f_4 * pc_z[k] * pgh_167[k];

        t_224[k] = pb_y[k] * pfi0_140[k]
                   - f_11 * pc_y[k] * pfi1_140[k];

        t_225[k] = f_0 * pfh_105[k]
                   + f_4 * pc_y[k] * pgh_168[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_y, pc_z, pfh_84, pfh_106, pfh_107, pgg0_120, \
                         pgg1_120, pgh_168, pgh_169, pgh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_12 * pfh_84[k]
                   + f_4 * pc_z[k] * pgh_168[k];

        t_227[k] = f_0 * pfh_106[k]
                   + f_5 * pgg0_120[k]
                   - f_6 * pgg1_120[k]
                   + f_4 * pc_y[k] * pgh_169[k];

        t_228[k] = f_0 * pfh_107[k]
                   + f_4 * pc_y[k] * pgh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_y, pc_y, pc_z, pfi0_145, pfh_87, \
                         pfh_108, pfh_110, pfi1_145, pgg0_121, pgg1_121, pgh_171, \
                         pgh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * pfi0_145[k]
                   - f_11 * pc_y[k] * pfi1_145[k];

        t_230[k] = f_0 * pfh_108[k]
                   + f_7 * pgg0_121[k]
                   - f_8 * pgg1_121[k]
                   + f_4 * pc_y[k] * pgh_171[k];

        t_231[k] = f_12 * pfh_87[k]
                   + f_4 * pc_z[k] * pgh_171[k];

        t_232[k] = f_0 * pfh_110[k]
                   + f_4 * pc_y[k] * pgh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_y, pc_y, pc_z, pfi0_149, pfh_90, pfh_111, \
                         pfi1_149, pgg0_123, pgg1_123, pgh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * pfi0_149[k]
                   - f_11 * pc_y[k] * pfi1_149[k];

        t_234[k] = f_0 * pfh_111[k]
                   + f_9 * pgg0_123[k]
                   - f_10 * pgg1_123[k]
                   + f_4 * pc_y[k] * pgh_174[k];

        t_235[k] = f_12 * pfh_90[k]
                   + f_4 * pc_z[k] * pgh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pb_y, pc_y, pfi0_154, pfh_113, pfh_114, \
                         pfi1_154, pgg0_125, pgg1_125, pgh_176, \
                         pgh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_0 * pfh_113[k]
                   + f_5 * pgg0_125[k]
                   - f_6 * pgg1_125[k]
                   + f_4 * pc_y[k] * pgh_176[k];

        t_237[k] = f_0 * pfh_114[k]
                   + f_4 * pc_y[k] * pgh_177[k];

        t_238[k] = pb_y[k] * pfi0_154[k]
                   - f_11 * pc_y[k] * pfi1_154[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgi0, const size_t sgh,
                                                          const size_t sgi1, const size_t pfi0,
                                                          const size_t pfh, const size_t pfi1,
                                                          const size_t pgg0, const size_t pgg1,
                                                          const size_t pgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_16 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_280 = buffer.data(sgi0 + 280);
    const auto *sgi0_283 = buffer.data(sgi0 + 283);
    const auto *sgi0_286 = buffer.data(sgi0 + 286);
    const auto *sgi0_290 = buffer.data(sgi0 + 290);
    const auto *sgi0_292 = buffer.data(sgi0 + 292);
    const auto *sgi0_301 = buffer.data(sgi0 + 301);
    const auto *sgi0_303 = buffer.data(sgi0 + 303);
    const auto *sgi0_304 = buffer.data(sgi0 + 304);
    const auto *sgi0_305 = buffer.data(sgi0 + 305);
    const auto *sgi0_307 = buffer.data(sgi0 + 307);
    const auto *sgi0_313 = buffer.data(sgi0 + 313);
    const auto *sgi0_317 = buffer.data(sgi0 + 317);
    const auto *sgi0_320 = buffer.data(sgi0 + 320);
    const auto *sgi0_322 = buffer.data(sgi0 + 322);
    const auto *sgi0_329 = buffer.data(sgi0 + 329);
    const auto *sgi0_331 = buffer.data(sgi0 + 331);
    const auto *sgi0_332 = buffer.data(sgi0 + 332);
    const auto *sgi0_333 = buffer.data(sgi0 + 333);
    const auto *sgi0_335 = buffer.data(sgi0 + 335);
    const auto *sgi0_336 = buffer.data(sgi0 + 336);
    const auto *sgi0_339 = buffer.data(sgi0 + 339);
    const auto *sgi0_341 = buffer.data(sgi0 + 341);
    const auto *sgi0_342 = buffer.data(sgi0 + 342);
    const auto *sgi0_345 = buffer.data(sgi0 + 345);
    const auto *sgi0_346 = buffer.data(sgi0 + 346);
    const auto *sgi0_348 = buffer.data(sgi0 + 348);
    const auto *sgi0_350 = buffer.data(sgi0 + 350);
    const auto *sgi0_357 = buffer.data(sgi0 + 357);

    const auto *sgh_183 = buffer.data(sgh + 183);
    const auto *sgh_185 = buffer.data(sgh + 185);
    const auto *sgh_186 = buffer.data(sgh + 186);
    const auto *sgh_189 = buffer.data(sgh + 189);
    const auto *sgh_204 = buffer.data(sgh + 204);
    const auto *sgh_206 = buffer.data(sgh + 206);
    const auto *sgh_207 = buffer.data(sgh + 207);
    const auto *sgh_209 = buffer.data(sgh + 209);
    const auto *sgh_210 = buffer.data(sgh + 210);
    const auto *sgh_213 = buffer.data(sgh + 213);
    const auto *sgh_216 = buffer.data(sgh + 216);
    const auto *sgh_220 = buffer.data(sgh + 220);
    const auto *sgh_222 = buffer.data(sgh + 222);
    const auto *sgh_225 = buffer.data(sgh + 225);
    const auto *sgh_227 = buffer.data(sgh + 227);
    const auto *sgh_228 = buffer.data(sgh + 228);
    const auto *sgh_230 = buffer.data(sgh + 230);
    const auto *sgh_236 = buffer.data(sgh + 236);
    const auto *sgh_240 = buffer.data(sgh + 240);
    const auto *sgh_243 = buffer.data(sgh + 243);
    const auto *sgh_245 = buffer.data(sgh + 245);
    const auto *sgh_246 = buffer.data(sgh + 246);
    const auto *sgh_248 = buffer.data(sgh + 248);
    const auto *sgh_249 = buffer.data(sgh + 249);
    const auto *sgh_251 = buffer.data(sgh + 251);
    const auto *sgh_252 = buffer.data(sgh + 252);
    const auto *sgh_255 = buffer.data(sgh + 255);
    const auto *sgh_257 = buffer.data(sgh + 257);
    const auto *sgh_258 = buffer.data(sgh + 258);
    const auto *sgh_261 = buffer.data(sgh + 261);
    const auto *sgh_262 = buffer.data(sgh + 262);
    const auto *sgh_264 = buffer.data(sgh + 264);
    const auto *sgh_266 = buffer.data(sgh + 266);
    const auto *sgh_267 = buffer.data(sgh + 267);
    const auto *sgh_269 = buffer.data(sgh + 269);
    const auto *sgh_270 = buffer.data(sgh + 270);
    const auto *sgh_272 = buffer.data(sgh + 272);

    const auto *sgi1_280 = buffer.data(sgi1 + 280);
    const auto *sgi1_283 = buffer.data(sgi1 + 283);
    const auto *sgi1_286 = buffer.data(sgi1 + 286);
    const auto *sgi1_290 = buffer.data(sgi1 + 290);
    const auto *sgi1_292 = buffer.data(sgi1 + 292);
    const auto *sgi1_301 = buffer.data(sgi1 + 301);
    const auto *sgi1_303 = buffer.data(sgi1 + 303);
    const auto *sgi1_304 = buffer.data(sgi1 + 304);
    const auto *sgi1_305 = buffer.data(sgi1 + 305);
    const auto *sgi1_307 = buffer.data(sgi1 + 307);
    const auto *sgi1_313 = buffer.data(sgi1 + 313);
    const auto *sgi1_317 = buffer.data(sgi1 + 317);
    const auto *sgi1_320 = buffer.data(sgi1 + 320);
    const auto *sgi1_322 = buffer.data(sgi1 + 322);
    const auto *sgi1_329 = buffer.data(sgi1 + 329);
    const auto *sgi1_331 = buffer.data(sgi1 + 331);
    const auto *sgi1_332 = buffer.data(sgi1 + 332);
    const auto *sgi1_333 = buffer.data(sgi1 + 333);
    const auto *sgi1_335 = buffer.data(sgi1 + 335);
    const auto *sgi1_336 = buffer.data(sgi1 + 336);
    const auto *sgi1_339 = buffer.data(sgi1 + 339);
    const auto *sgi1_341 = buffer.data(sgi1 + 341);
    const auto *sgi1_342 = buffer.data(sgi1 + 342);
    const auto *sgi1_345 = buffer.data(sgi1 + 345);
    const auto *sgi1_346 = buffer.data(sgi1 + 346);
    const auto *sgi1_348 = buffer.data(sgi1 + 348);
    const auto *sgi1_350 = buffer.data(sgi1 + 350);
    const auto *sgi1_357 = buffer.data(sgi1 + 357);

    const auto *pfi0_160 = buffer.data(pfi0 + 160);
    const auto *pfi0_168 = buffer.data(pfi0 + 168);
    const auto *pfi0_171 = buffer.data(pfi0 + 171);
    const auto *pfi0_174 = buffer.data(pfi0 + 174);
    const auto *pfi0_178 = buffer.data(pfi0 + 178);

    const auto *pfh_94 = buffer.data(pfh + 94);
    const auto *pfh_99 = buffer.data(pfh + 99);
    const auto *pfh_104 = buffer.data(pfh + 104);
    const auto *pfh_105 = buffer.data(pfh + 105);
    const auto *pfh_107 = buffer.data(pfh + 107);
    const auto *pfh_108 = buffer.data(pfh + 108);
    const auto *pfh_110 = buffer.data(pfh + 110);
    const auto *pfh_111 = buffer.data(pfh + 111);
    const auto *pfh_114 = buffer.data(pfh + 114);
    const auto *pfh_115 = buffer.data(pfh + 115);
    const auto *pfh_119 = buffer.data(pfh + 119);
    const auto *pfh_120 = buffer.data(pfh + 120);
    const auto *pfh_122 = buffer.data(pfh + 122);
    const auto *pfh_123 = buffer.data(pfh + 123);
    const auto *pfh_124 = buffer.data(pfh + 124);
    const auto *pfh_125 = buffer.data(pfh + 125);
    const auto *pfh_126 = buffer.data(pfh + 126);
    const auto *pfh_128 = buffer.data(pfh + 128);
    const auto *pfh_129 = buffer.data(pfh + 129);
    const auto *pfh_131 = buffer.data(pfh + 131);
    const auto *pfh_132 = buffer.data(pfh + 132);
    const auto *pfh_135 = buffer.data(pfh + 135);
    const auto *pfh_136 = buffer.data(pfh + 136);
    const auto *pfh_140 = buffer.data(pfh + 140);
    const auto *pfh_141 = buffer.data(pfh + 141);
    const auto *pfh_146 = buffer.data(pfh + 146);
    const auto *pfh_147 = buffer.data(pfh + 147);
    const auto *pfh_149 = buffer.data(pfh + 149);
    const auto *pfh_150 = buffer.data(pfh + 150);
    const auto *pfh_152 = buffer.data(pfh + 152);
    const auto *pfh_153 = buffer.data(pfh + 153);
    const auto *pfh_156 = buffer.data(pfh + 156);
    const auto *pfh_157 = buffer.data(pfh + 157);
    const auto *pfh_161 = buffer.data(pfh + 161);
    const auto *pfh_162 = buffer.data(pfh + 162);
    const auto *pfh_167 = buffer.data(pfh + 167);
    const auto *pfh_168 = buffer.data(pfh + 168);
    const auto *pfh_170 = buffer.data(pfh + 170);
    const auto *pfh_173 = buffer.data(pfh + 173);
    const auto *pfh_177 = buffer.data(pfh + 177);
    const auto *pfh_182 = buffer.data(pfh + 182);
    const auto *pfh_183 = buffer.data(pfh + 183);
    const auto *pfh_185 = buffer.data(pfh + 185);
    const auto *pfh_186 = buffer.data(pfh + 186);
    const auto *pfh_189 = buffer.data(pfh + 189);
    const auto *pfh_204 = buffer.data(pfh + 204);
    const auto *pfh_206 = buffer.data(pfh + 206);
    const auto *pfh_207 = buffer.data(pfh + 207);
    const auto *pfh_209 = buffer.data(pfh + 209);

    const auto *pfi1_160 = buffer.data(pfi1 + 160);
    const auto *pfi1_168 = buffer.data(pfi1 + 168);
    const auto *pfi1_171 = buffer.data(pfi1 + 171);
    const auto *pfi1_174 = buffer.data(pfi1 + 174);
    const auto *pfi1_178 = buffer.data(pfi1 + 178);

    const auto *pgg0_130 = buffer.data(pgg0 + 130);
    const auto *pgg0_132 = buffer.data(pgg0 + 132);
    const auto *pgg0_133 = buffer.data(pgg0 + 133);
    const auto *pgg0_134 = buffer.data(pgg0 + 134);
    const auto *pgg0_135 = buffer.data(pgg0 + 135);
    const auto *pgg0_136 = buffer.data(pgg0 + 136);
    const auto *pgg0_137 = buffer.data(pgg0 + 137);
    const auto *pgg0_138 = buffer.data(pgg0 + 138);
    const auto *pgg0_140 = buffer.data(pgg0 + 140);
    const auto *pgg0_145 = buffer.data(pgg0 + 145);
    const auto *pgg0_147 = buffer.data(pgg0 + 147);
    const auto *pgg0_148 = buffer.data(pgg0 + 148);
    const auto *pgg0_149 = buffer.data(pgg0 + 149);
    const auto *pgg0_150 = buffer.data(pgg0 + 150);
    const auto *pgg0_152 = buffer.data(pgg0 + 152);
    const auto *pgg0_155 = buffer.data(pgg0 + 155);

    const auto *pgg1_130 = buffer.data(pgg1 + 130);
    const auto *pgg1_132 = buffer.data(pgg1 + 132);
    const auto *pgg1_133 = buffer.data(pgg1 + 133);
    const auto *pgg1_134 = buffer.data(pgg1 + 134);
    const auto *pgg1_135 = buffer.data(pgg1 + 135);
    const auto *pgg1_136 = buffer.data(pgg1 + 136);
    const auto *pgg1_137 = buffer.data(pgg1 + 137);
    const auto *pgg1_138 = buffer.data(pgg1 + 138);
    const auto *pgg1_140 = buffer.data(pgg1 + 140);
    const auto *pgg1_145 = buffer.data(pgg1 + 145);
    const auto *pgg1_147 = buffer.data(pgg1 + 147);
    const auto *pgg1_148 = buffer.data(pgg1 + 148);
    const auto *pgg1_149 = buffer.data(pgg1 + 149);
    const auto *pgg1_150 = buffer.data(pgg1 + 150);
    const auto *pgg1_152 = buffer.data(pgg1 + 152);
    const auto *pgg1_155 = buffer.data(pgg1 + 155);

    const auto *pgh_178 = buffer.data(pgh + 178);
    const auto *pgh_182 = buffer.data(pgh + 182);
    const auto *pgh_183 = buffer.data(pgh + 183);
    const auto *pgh_185 = buffer.data(pgh + 185);
    const auto *pgh_186 = buffer.data(pgh + 186);
    const auto *pgh_187 = buffer.data(pgh + 187);
    const auto *pgh_188 = buffer.data(pgh + 188);
    const auto *pgh_189 = buffer.data(pgh + 189);
    const auto *pgh_190 = buffer.data(pgh + 190);
    const auto *pgh_191 = buffer.data(pgh + 191);
    const auto *pgh_192 = buffer.data(pgh + 192);
    const auto *pgh_194 = buffer.data(pgh + 194);
    const auto *pgh_195 = buffer.data(pgh + 195);
    const auto *pgh_197 = buffer.data(pgh + 197);
    const auto *pgh_198 = buffer.data(pgh + 198);
    const auto *pgh_199 = buffer.data(pgh + 199);
    const auto *pgh_203 = buffer.data(pgh + 203);
    const auto *pgh_204 = buffer.data(pgh + 204);
    const auto *pgh_206 = buffer.data(pgh + 206);
    const auto *pgh_207 = buffer.data(pgh + 207);
    const auto *pgh_208 = buffer.data(pgh + 208);
    const auto *pgh_209 = buffer.data(pgh + 209);
    const auto *pgh_210 = buffer.data(pgh + 210);
    const auto *pgh_212 = buffer.data(pgh + 212);
    const auto *pgh_213 = buffer.data(pgh + 213);
    const auto *pgh_215 = buffer.data(pgh + 215);
    const auto *pgh_216 = buffer.data(pgh + 216);
    const auto *pgh_219 = buffer.data(pgh + 219);
    const auto *pgh_220 = buffer.data(pgh + 220);
    const auto *pgh_224 = buffer.data(pgh + 224);
    const auto *pgh_225 = buffer.data(pgh + 225);
    const auto *pgh_227 = buffer.data(pgh + 227);
    const auto *pgh_228 = buffer.data(pgh + 228);
    const auto *pgh_230 = buffer.data(pgh + 230);
    const auto *pgh_231 = buffer.data(pgh + 231);
    const auto *pgh_233 = buffer.data(pgh + 233);
    const auto *pgh_234 = buffer.data(pgh + 234);
    const auto *pgh_236 = buffer.data(pgh + 236);
    const auto *pgh_237 = buffer.data(pgh + 237);
    const auto *pgh_240 = buffer.data(pgh + 240);
    const auto *pgh_241 = buffer.data(pgh + 241);
    const auto *pgh_245 = buffer.data(pgh + 245);
    const auto *pgh_246 = buffer.data(pgh + 246);
    const auto *pgh_248 = buffer.data(pgh + 248);
    const auto *pgh_249 = buffer.data(pgh + 249);
    const auto *pgh_251 = buffer.data(pgh + 251);
    const auto *pgh_252 = buffer.data(pgh + 252);
    const auto *pgh_254 = buffer.data(pgh + 254);
    const auto *pgh_255 = buffer.data(pgh + 255);
    const auto *pgh_257 = buffer.data(pgh + 257);
    const auto *pgh_258 = buffer.data(pgh + 258);
    const auto *pgh_261 = buffer.data(pgh + 261);
    const auto *pgh_262 = buffer.data(pgh + 262);
    const auto *pgh_266 = buffer.data(pgh + 266);
    const auto *pgh_267 = buffer.data(pgh + 267);
    const auto *pgh_269 = buffer.data(pgh + 269);
    const auto *pgh_270 = buffer.data(pgh + 270);
    const auto *pgh_272 = buffer.data(pgh + 272);

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, sgh_183, sgh_185, pfh_94, pfh_183, \
                         pfh_185, pgh_178, pgh_183, pgh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_0 * sgh_183[k]
                   + f_0 * pfh_183[k]
                   + f_4 * pc_x[k] * pgh_183[k];

        t_240[k] = f_12 * pfh_94[k]
                   + f_4 * pc_z[k] * pgh_178[k];

        t_241[k] = f_0 * sgh_185[k]
                   + f_0 * pfh_185[k]
                   + f_4 * pc_x[k] * pgh_185[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pb_y, pc_x, pc_y, sgh_186, pfi0_160, pfh_119, \
                         pfh_186, pfi1_160, pgh_182, pgh_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_0 * sgh_186[k]
                   + f_0 * pfh_186[k]
                   + f_4 * pc_x[k] * pgh_186[k];

        t_243[k] = f_0 * pfh_119[k]
                   + f_4 * pc_y[k] * pgh_182[k];

        t_244[k] = pb_y[k] * pfi0_160[k]
                   - f_11 * pc_y[k] * pfi1_160[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, pfh_99, pfh_120, pfh_122, pgg0_130, \
                         pgg0_132, pgg1_130, pgg1_132, pgh_183, \
                         pgh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_0 * pfh_120[k]
                   + f_2 * pgg0_130[k]
                   - f_3 * pgg1_130[k]
                   + f_4 * pc_y[k] * pgh_183[k];

        t_246[k] = f_12 * pfh_99[k]
                   + f_4 * pc_z[k] * pgh_183[k];

        t_247[k] = f_0 * pfh_122[k]
                   + f_9 * pgg0_132[k]
                   - f_10 * pgg1_132[k]
                   + f_4 * pc_y[k] * pgh_185[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, pfh_123, pfh_124, pfh_125, pgg0_133, \
                         pgg0_134, pgg1_133, pgg1_134, pgh_186, pgh_187, \
                         pgh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_0 * pfh_123[k]
                   + f_7 * pgg0_133[k]
                   - f_8 * pgg1_133[k]
                   + f_4 * pc_y[k] * pgh_186[k];

        t_249[k] = f_0 * pfh_124[k]
                   + f_5 * pgg0_134[k]
                   - f_6 * pgg1_134[k]
                   + f_4 * pc_y[k] * pgh_187[k];

        t_250[k] = f_0 * pfh_125[k]
                   + f_4 * pc_y[k] * pgh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, pc_x, pc_y, pc_z, sgh_189, pfh_104, pfh_189, \
                         pgg0_134, pgg0_135, pgg1_134, pgg1_135, pgh_188, \
                         pgh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_12 * pfh_104[k]
                   + f_2 * pgg0_134[k]
                   - f_3 * pgg1_134[k]
                   + f_4 * pc_z[k] * pgh_188[k];

        t_252[k] = f_0 * sgh_189[k]
                   + f_0 * pfh_189[k]
                   + f_2 * pgg0_135[k]
                   - f_3 * pgg1_135[k]
                   + f_4 * pc_x[k] * pgh_189[k];

        t_253[k] = f_4 * pc_y[k] * pgh_189[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pc_y, pc_z, pfh_105, pfh_107, pgg0_135, \
                         pgg1_135, pgh_189, pgh_190, pgh_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * pfh_105[k]
                   + f_4 * pc_z[k] * pgh_189[k];

        t_255[k] = f_5 * pgg0_135[k]
                   - f_6 * pgg1_135[k]
                   + f_4 * pc_y[k] * pgh_190[k];

        t_256[k] = f_4 * pc_y[k] * pgh_191[k];

        t_257[k] = f_13 * pfh_107[k]
                   + f_5 * pgg0_135[k]
                   - f_6 * pgg1_135[k]
                   + f_4 * pc_z[k] * pgh_191[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, pc_y, pc_z, pfh_108, pfh_110, pgg0_136, \
                         pgg0_137, pgg1_136, pgg1_137, pgh_192, \
                         pgh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_7 * pgg0_136[k]
                   - f_8 * pgg1_136[k]
                   + f_4 * pc_y[k] * pgh_192[k];

        t_259[k] = f_13 * pfh_108[k]
                   + f_4 * pc_z[k] * pgh_192[k];

        t_260[k] = f_4 * pc_y[k] * pgh_194[k];

        t_261[k] = f_13 * pfh_110[k]
                   + f_7 * pgg0_137[k]
                   - f_8 * pgg1_137[k]
                   + f_4 * pc_z[k] * pgh_194[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, pc_y, pc_z, pfh_111, pfh_114, \
                         pgg0_138, pgg0_140, pgg1_138, pgg1_140, pgh_195, pgh_197, \
                         pgh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_9 * pgg0_138[k]
                   - f_10 * pgg1_138[k]
                   + f_4 * pc_y[k] * pgh_195[k];

        t_263[k] = f_13 * pfh_111[k]
                   + f_4 * pc_z[k] * pgh_195[k];

        t_264[k] = f_5 * pgg0_140[k]
                   - f_6 * pgg1_140[k]
                   + f_4 * pc_y[k] * pgh_197[k];

        t_265[k] = f_4 * pc_y[k] * pgh_198[k];

        t_266[k] = f_13 * pfh_114[k]
                   + f_9 * pgg0_140[k]
                   - f_10 * pgg1_140[k]
                   + f_4 * pc_z[k] * pgh_198[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_x, pc_z, sgh_204, sgh_206, pfh_115, pfh_204, \
                         pfh_206, pgh_199, pgh_204, pgh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_0 * sgh_204[k]
                   + f_0 * pfh_204[k]
                   + f_4 * pc_x[k] * pgh_204[k];

        t_268[k] = f_13 * pfh_115[k]
                   + f_4 * pc_z[k] * pgh_199[k];

        t_269[k] = f_0 * sgh_206[k]
                   + f_0 * pfh_206[k]
                   + f_4 * pc_x[k] * pgh_206[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pc_x, pc_y, sgh_207, sgh_209, pfh_207, \
                         pfh_209, pgg0_145, pgg1_145, pgh_203, pgh_204, pgh_207, \
                         pgh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * sgh_207[k]
                   + f_0 * pfh_207[k]
                   + f_4 * pc_x[k] * pgh_207[k];

        t_271[k] = f_4 * pc_y[k] * pgh_203[k];

        t_272[k] = f_0 * sgh_209[k]
                   + f_0 * pfh_209[k]
                   + f_4 * pc_x[k] * pgh_209[k];

        t_273[k] = f_2 * pgg0_145[k]
                   - f_3 * pgg1_145[k]
                   + f_4 * pc_y[k] * pgh_204[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, pc_y, pc_z, pfh_120, pgg0_147, pgg0_148, \
                         pgg1_147, pgg1_148, pgh_204, pgh_206, \
                         pgh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_13 * pfh_120[k]
                   + f_4 * pc_z[k] * pgh_204[k];

        t_275[k] = f_9 * pgg0_147[k]
                   - f_10 * pgg1_147[k]
                   + f_4 * pc_y[k] * pgh_206[k];

        t_276[k] = f_7 * pgg0_148[k]
                   - f_8 * pgg1_148[k]
                   + f_4 * pc_y[k] * pgh_207[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, pa_x, pc_x, pc_y, pc_z, sgi0_280, \
                         sgh_210, sgi1_280, pfh_125, pgg0_149, pgg1_149, pgh_208, \
                         pgh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_5 * pgg0_149[k]
                   - f_6 * pgg1_149[k]
                   + f_4 * pc_y[k] * pgh_208[k];

        t_278[k] = f_4 * pc_y[k] * pgh_209[k];

        t_279[k] = f_13 * pfh_125[k]
                   + f_2 * pgg0_149[k]
                   - f_3 * pgg1_149[k]
                   + f_4 * pc_z[k] * pgh_209[k];

        t_280[k] = pa_x[k] * sgi0_280[k]
                   + f_16 * sgh_210[k]
                   - f_11 * pc_x[k] * sgi1_280[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pa_x, pc_x, pc_y, pc_z, sgi0_283, \
                         sgh_213, sgi1_283, pfh_126, pfh_128, pgh_210, \
                         pgh_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_1 * pfh_126[k]
                   + f_4 * pc_y[k] * pgh_210[k];

        t_282[k] = f_4 * pc_z[k] * pgh_210[k];

        t_283[k] = pa_x[k] * sgi0_283[k]
                   + f_1 * sgh_213[k]
                   - f_11 * pc_x[k] * sgi1_283[k];

        t_284[k] = f_1 * pfh_128[k]
                   + f_4 * pc_y[k] * pgh_212[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pa_x, pc_x, pc_z, sgi0_286, sgh_216, sgi1_286, \
                         pgg0_150, pgg1_150, pgh_212, pgh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_5 * pgg0_150[k]
                   - f_6 * pgg1_150[k]
                   + f_4 * pc_z[k] * pgh_212[k];

        t_286[k] = pa_x[k] * sgi0_286[k]
                   + f_13 * sgh_216[k]
                   - f_11 * pc_x[k] * sgi1_286[k];

        t_287[k] = f_4 * pc_z[k] * pgh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pa_x, pc_x, pc_y, pc_z, sgi0_290, \
                         sgh_220, sgi1_290, pfh_131, pgg0_152, pgg1_152, pgh_215, \
                         pgh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_1 * pfh_131[k]
                   + f_4 * pc_y[k] * pgh_215[k];

        t_289[k] = f_7 * pgg0_152[k]
                   - f_8 * pgg1_152[k]
                   + f_4 * pc_z[k] * pgh_215[k];

        t_290[k] = pa_x[k] * sgi0_290[k]
                   + f_12 * sgh_220[k]
                   - f_11 * pc_x[k] * sgi1_290[k];

        t_291[k] = f_4 * pc_z[k] * pgh_216[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, pa_x, pc_x, pc_y, pc_z, sgi0_292, sgh_222, \
                         sgi1_292, pfh_135, pgg0_155, pgg1_155, \
                         pgh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = pa_x[k] * sgi0_292[k]
                   + f_12 * sgh_222[k]
                   - f_11 * pc_x[k] * sgi1_292[k];

        t_293[k] = f_1 * pfh_135[k]
                   + f_4 * pc_y[k] * pgh_219[k];

        t_294[k] = f_9 * pgg0_155[k]
                   - f_10 * pgg1_155[k]
                   + f_4 * pc_z[k] * pgh_219[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pc_x, pc_z, sgh_225, sgh_227, sgh_228, \
                         pgh_220, pgh_225, pgh_227, pgh_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_0 * sgh_225[k]
                   + f_4 * pc_x[k] * pgh_225[k];

        t_296[k] = f_4 * pc_z[k] * pgh_220[k];

        t_297[k] = f_0 * sgh_227[k]
                   + f_4 * pc_x[k] * pgh_227[k];

        t_298[k] = f_0 * sgh_228[k]
                   + f_4 * pc_x[k] * pgh_228[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pa_x, pc_x, pc_y, pc_z, sgi0_301, \
                         sgh_230, sgi1_301, pfh_140, pgh_224, pgh_225, \
                         pgh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_1 * pfh_140[k]
                   + f_4 * pc_y[k] * pgh_224[k];

        t_300[k] = f_0 * sgh_230[k]
                   + f_4 * pc_x[k] * pgh_230[k];

        t_301[k] = pa_x[k] * sgi0_301[k]
                   - f_11 * pc_x[k] * sgi1_301[k];

        t_302[k] = f_4 * pc_z[k] * pgh_225[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pa_x, pc_x, pc_y, sgi0_303, sgi0_304, \
                         sgi0_305, sgi1_303, sgi1_304, sgi1_305, pfh_146, \
                         pgh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pa_x[k] * sgi0_303[k]
                   - f_11 * pc_x[k] * sgi1_303[k];

        t_304[k] = pa_x[k] * sgi0_304[k]
                   - f_11 * pc_x[k] * sgi1_304[k];

        t_305[k] = pa_x[k] * sgi0_305[k]
                   - f_11 * pc_x[k] * sgi1_305[k];

        t_306[k] = f_1 * pfh_146[k]
                   + f_4 * pc_y[k] * pgh_230[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pa_x, pb_z, pc_x, pc_y, pc_z, sgi0_307, \
                         sgi1_307, pfi0_168, pfh_126, pfh_147, pfi1_168, \
                         pgh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = pa_x[k] * sgi0_307[k]
                   - f_11 * pc_x[k] * sgi1_307[k];

        t_308[k] = pb_z[k] * pfi0_168[k]
                   - f_11 * pc_z[k] * pfi1_168[k];

        t_309[k] = f_13 * pfh_147[k]
                   + f_4 * pc_y[k] * pgh_231[k];

        t_310[k] = f_0 * pfh_126[k]
                   + f_4 * pc_z[k] * pgh_231[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, pa_x, pb_z, pc_x, pc_y, pc_z, sgi0_313, sgh_236, \
                         sgi1_313, pfi0_171, pfh_149, pfi1_171, \
                         pgh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = pb_z[k] * pfi0_171[k]
                   - f_11 * pc_z[k] * pfi1_171[k];

        t_312[k] = f_13 * pfh_149[k]
                   + f_4 * pc_y[k] * pgh_233[k];

        t_313[k] = pa_x[k] * sgi0_313[k]
                   + f_1 * sgh_236[k]
                   - f_11 * pc_x[k] * sgi1_313[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, pb_z, pc_y, pc_z, pfi0_174, pfh_129, pfh_152, \
                         pfi1_174, pgh_234, pgh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pb_z[k] * pfi0_174[k]
                   - f_11 * pc_z[k] * pfi1_174[k];

        t_315[k] = f_0 * pfh_129[k]
                   + f_4 * pc_z[k] * pgh_234[k];

        t_316[k] = f_13 * pfh_152[k]
                   + f_4 * pc_y[k] * pgh_236[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pa_x, pb_z, pc_x, pc_z, sgi0_317, sgh_240, \
                         sgi1_317, pfi0_178, pfh_132, pfi1_178, \
                         pgh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = pa_x[k] * sgi0_317[k]
                   + f_13 * sgh_240[k]
                   - f_11 * pc_x[k] * sgi1_317[k];

        t_318[k] = pb_z[k] * pfi0_178[k]
                   - f_11 * pc_z[k] * pfi1_178[k];

        t_319[k] = f_0 * pfh_132[k]
                   + f_4 * pc_z[k] * pgh_237[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pa_x, pc_x, pc_y, sgi0_320, sgi0_322, sgh_243, \
                         sgh_245, sgi1_320, sgi1_322, pfh_156, \
                         pgh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = pa_x[k] * sgi0_320[k]
                   + f_12 * sgh_243[k]
                   - f_11 * pc_x[k] * sgi1_320[k];

        t_321[k] = f_13 * pfh_156[k]
                   + f_4 * pc_y[k] * pgh_240[k];

        t_322[k] = pa_x[k] * sgi0_322[k]
                   + f_12 * sgh_245[k]
                   - f_11 * pc_x[k] * sgi1_322[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_z, sgh_246, sgh_248, sgh_249, \
                         pfh_136, pgh_241, pgh_246, pgh_248, pgh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_0 * sgh_246[k]
                   + f_4 * pc_x[k] * pgh_246[k];

        t_324[k] = f_0 * pfh_136[k]
                   + f_4 * pc_z[k] * pgh_241[k];

        t_325[k] = f_0 * sgh_248[k]
                   + f_4 * pc_x[k] * pgh_248[k];

        t_326[k] = f_0 * sgh_249[k]
                   + f_4 * pc_x[k] * pgh_249[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pa_x, pc_x, pc_y, pc_z, sgi0_329, \
                         sgh_251, sgi1_329, pfh_141, pfh_161, pgh_245, pgh_246, \
                         pgh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_13 * pfh_161[k]
                   + f_4 * pc_y[k] * pgh_245[k];

        t_328[k] = f_0 * sgh_251[k]
                   + f_4 * pc_x[k] * pgh_251[k];

        t_329[k] = pa_x[k] * sgi0_329[k]
                   - f_11 * pc_x[k] * sgi1_329[k];

        t_330[k] = f_0 * pfh_141[k]
                   + f_4 * pc_z[k] * pgh_246[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pa_x, pc_x, pc_y, sgi0_331, sgi0_332, \
                         sgi0_333, sgi1_331, sgi1_332, sgi1_333, pfh_167, \
                         pgh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_x[k] * sgi0_331[k]
                   - f_11 * pc_x[k] * sgi1_331[k];

        t_332[k] = pa_x[k] * sgi0_332[k]
                   - f_11 * pc_x[k] * sgi1_332[k];

        t_333[k] = pa_x[k] * sgi0_333[k]
                   - f_11 * pc_x[k] * sgi1_333[k];

        t_334[k] = f_13 * pfh_167[k]
                   + f_4 * pc_y[k] * pgh_251[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_x, pc_x, pc_y, pc_z, sgi0_335, \
                         sgi0_336, sgh_252, sgi1_335, sgi1_336, pfh_147, pfh_168, \
                         pgh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pa_x[k] * sgi0_335[k]
                   - f_11 * pc_x[k] * sgi1_335[k];

        t_336[k] = pa_x[k] * sgi0_336[k]
                   + f_16 * sgh_252[k]
                   - f_11 * pc_x[k] * sgi1_336[k];

        t_337[k] = f_12 * pfh_168[k]
                   + f_4 * pc_y[k] * pgh_252[k];

        t_338[k] = f_12 * pfh_147[k]
                   + f_4 * pc_z[k] * pgh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_x, pc_x, pc_y, sgi0_339, sgi0_341, sgh_255, \
                         sgh_257, sgi1_339, sgi1_341, pfh_170, \
                         pgh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_x[k] * sgi0_339[k]
                   + f_1 * sgh_255[k]
                   - f_11 * pc_x[k] * sgi1_339[k];

        t_340[k] = f_12 * pfh_170[k]
                   + f_4 * pc_y[k] * pgh_254[k];

        t_341[k] = pa_x[k] * sgi0_341[k]
                   + f_1 * sgh_257[k]
                   - f_11 * pc_x[k] * sgi1_341[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_x, pc_x, pc_y, pc_z, sgi0_342, sgh_258, \
                         sgi1_342, pfh_150, pfh_173, pgh_255, pgh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_x[k] * sgi0_342[k]
                   + f_13 * sgh_258[k]
                   - f_11 * pc_x[k] * sgi1_342[k];

        t_343[k] = f_12 * pfh_150[k]
                   + f_4 * pc_z[k] * pgh_255[k];

        t_344[k] = f_12 * pfh_173[k]
                   + f_4 * pc_y[k] * pgh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pa_x, pc_x, pc_z, sgi0_345, sgi0_346, sgh_261, \
                         sgh_262, sgi1_345, sgi1_346, pfh_153, \
                         pgh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = pa_x[k] * sgi0_345[k]
                   + f_13 * sgh_261[k]
                   - f_11 * pc_x[k] * sgi1_345[k];

        t_346[k] = pa_x[k] * sgi0_346[k]
                   + f_12 * sgh_262[k]
                   - f_11 * pc_x[k] * sgi1_346[k];

        t_347[k] = f_12 * pfh_153[k]
                   + f_4 * pc_z[k] * pgh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pa_x, pc_x, pc_y, sgi0_348, sgi0_350, sgh_264, \
                         sgh_266, sgi1_348, sgi1_350, pfh_177, \
                         pgh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = pa_x[k] * sgi0_348[k]
                   + f_12 * sgh_264[k]
                   - f_11 * pc_x[k] * sgi1_348[k];

        t_349[k] = f_12 * pfh_177[k]
                   + f_4 * pc_y[k] * pgh_261[k];

        t_350[k] = pa_x[k] * sgi0_350[k]
                   + f_12 * sgh_266[k]
                   - f_11 * pc_x[k] * sgi1_350[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, pc_x, pc_z, sgh_267, sgh_269, sgh_270, \
                         pfh_157, pgh_262, pgh_267, pgh_269, pgh_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_0 * sgh_267[k]
                   + f_4 * pc_x[k] * pgh_267[k];

        t_352[k] = f_12 * pfh_157[k]
                   + f_4 * pc_z[k] * pgh_262[k];

        t_353[k] = f_0 * sgh_269[k]
                   + f_4 * pc_x[k] * pgh_269[k];

        t_354[k] = f_0 * sgh_270[k]
                   + f_4 * pc_x[k] * pgh_270[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pa_x, pc_x, pc_y, pc_z, sgi0_357, \
                         sgh_272, sgi1_357, pfh_162, pfh_182, pgh_266, pgh_267, \
                         pgh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_12 * pfh_182[k]
                   + f_4 * pc_y[k] * pgh_266[k];

        t_356[k] = f_0 * sgh_272[k]
                   + f_4 * pc_x[k] * pgh_272[k];

        t_357[k] = pa_x[k] * sgi0_357[k]
                   - f_11 * pc_x[k] * sgi1_357[k];

        t_358[k] = f_12 * pfh_162[k]
                   + f_4 * pc_z[k] * pgh_267[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgi0, const size_t sgh,
                                                          const size_t sgi1, const size_t pdi0,
                                                          const size_t pdi1, const size_t pfi0,
                                                          const size_t pfh, const size_t pfi1,
                                                          const size_t pgg0, const size_t pgg1,
                                                          const size_t pgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 1.0 / p;
    const auto f_20 = gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_0 = buffer.data(sgi0 + 0);
    const auto *sgi0_1 = buffer.data(sgi0 + 1);
    const auto *sgi0_3 = buffer.data(sgi0 + 3);
    const auto *sgi0_5 = buffer.data(sgi0 + 5);
    const auto *sgi0_6 = buffer.data(sgi0 + 6);
    const auto *sgi0_8 = buffer.data(sgi0 + 8);
    const auto *sgi0_9 = buffer.data(sgi0 + 9);
    const auto *sgi0_10 = buffer.data(sgi0 + 10);
    const auto *sgi0_12 = buffer.data(sgi0 + 12);
    const auto *sgi0_13 = buffer.data(sgi0 + 13);
    const auto *sgi0_14 = buffer.data(sgi0 + 14);
    const auto *sgi0_21 = buffer.data(sgi0 + 21);
    const auto *sgi0_27 = buffer.data(sgi0 + 27);
    const auto *sgi0_56 = buffer.data(sgi0 + 56);
    const auto *sgi0_359 = buffer.data(sgi0 + 359);
    const auto *sgi0_360 = buffer.data(sgi0 + 360);
    const auto *sgi0_361 = buffer.data(sgi0 + 361);
    const auto *sgi0_363 = buffer.data(sgi0 + 363);
    const auto *sgi0_367 = buffer.data(sgi0 + 367);
    const auto *sgi0_370 = buffer.data(sgi0 + 370);
    const auto *sgi0_374 = buffer.data(sgi0 + 374);
    const auto *sgi0_376 = buffer.data(sgi0 + 376);
    const auto *sgi0_385 = buffer.data(sgi0 + 385);
    const auto *sgi0_387 = buffer.data(sgi0 + 387);
    const auto *sgi0_388 = buffer.data(sgi0 + 388);
    const auto *sgi0_389 = buffer.data(sgi0 + 389);
    const auto *sgi0_391 = buffer.data(sgi0 + 391);
    const auto *sgi0_392 = buffer.data(sgi0 + 392);
    const auto *sgi0_397 = buffer.data(sgi0 + 397);
    const auto *sgi0_401 = buffer.data(sgi0 + 401);
    const auto *sgi0_406 = buffer.data(sgi0 + 406);
    const auto *sgi0_413 = buffer.data(sgi0 + 413);
    const auto *sgi0_415 = buffer.data(sgi0 + 415);
    const auto *sgi0_416 = buffer.data(sgi0 + 416);
    const auto *sgi0_417 = buffer.data(sgi0 + 417);
    const auto *sgi0_419 = buffer.data(sgi0 + 419);

    const auto *sgh_0 = buffer.data(sgh + 0);
    const auto *sgh_1 = buffer.data(sgh + 1);
    const auto *sgh_3 = buffer.data(sgh + 3);
    const auto *sgh_5 = buffer.data(sgh + 5);
    const auto *sgh_6 = buffer.data(sgh + 6);
    const auto *sgh_8 = buffer.data(sgh + 8);
    const auto *sgh_9 = buffer.data(sgh + 9);
    const auto *sgh_15 = buffer.data(sgh + 15);
    const auto *sgh_20 = buffer.data(sgh + 20);
    const auto *sgh_41 = buffer.data(sgh + 41);
    const auto *sgh_276 = buffer.data(sgh + 276);
    const auto *sgh_279 = buffer.data(sgh + 279);
    const auto *sgh_283 = buffer.data(sgh + 283);
    const auto *sgh_285 = buffer.data(sgh + 285);
    const auto *sgh_288 = buffer.data(sgh + 288);
    const auto *sgh_290 = buffer.data(sgh + 290);
    const auto *sgh_291 = buffer.data(sgh + 291);
    const auto *sgh_293 = buffer.data(sgh + 293);
    const auto *sgh_294 = buffer.data(sgh + 294);
    const auto *sgh_299 = buffer.data(sgh + 299);
    const auto *sgh_303 = buffer.data(sgh + 303);
    const auto *sgh_308 = buffer.data(sgh + 308);
    const auto *sgh_309 = buffer.data(sgh + 309);
    const auto *sgh_311 = buffer.data(sgh + 311);
    const auto *sgh_312 = buffer.data(sgh + 312);
    const auto *sgh_314 = buffer.data(sgh + 314);

    const auto *sgi1_0 = buffer.data(sgi1 + 0);
    const auto *sgi1_1 = buffer.data(sgi1 + 1);
    const auto *sgi1_3 = buffer.data(sgi1 + 3);
    const auto *sgi1_5 = buffer.data(sgi1 + 5);
    const auto *sgi1_6 = buffer.data(sgi1 + 6);
    const auto *sgi1_8 = buffer.data(sgi1 + 8);
    const auto *sgi1_9 = buffer.data(sgi1 + 9);
    const auto *sgi1_10 = buffer.data(sgi1 + 10);
    const auto *sgi1_12 = buffer.data(sgi1 + 12);
    const auto *sgi1_13 = buffer.data(sgi1 + 13);
    const auto *sgi1_14 = buffer.data(sgi1 + 14);
    const auto *sgi1_21 = buffer.data(sgi1 + 21);
    const auto *sgi1_27 = buffer.data(sgi1 + 27);
    const auto *sgi1_56 = buffer.data(sgi1 + 56);
    const auto *sgi1_359 = buffer.data(sgi1 + 359);
    const auto *sgi1_360 = buffer.data(sgi1 + 360);
    const auto *sgi1_361 = buffer.data(sgi1 + 361);
    const auto *sgi1_363 = buffer.data(sgi1 + 363);
    const auto *sgi1_367 = buffer.data(sgi1 + 367);
    const auto *sgi1_370 = buffer.data(sgi1 + 370);
    const auto *sgi1_374 = buffer.data(sgi1 + 374);
    const auto *sgi1_376 = buffer.data(sgi1 + 376);
    const auto *sgi1_385 = buffer.data(sgi1 + 385);
    const auto *sgi1_387 = buffer.data(sgi1 + 387);
    const auto *sgi1_388 = buffer.data(sgi1 + 388);
    const auto *sgi1_389 = buffer.data(sgi1 + 389);
    const auto *sgi1_391 = buffer.data(sgi1 + 391);
    const auto *sgi1_392 = buffer.data(sgi1 + 392);
    const auto *sgi1_397 = buffer.data(sgi1 + 397);
    const auto *sgi1_401 = buffer.data(sgi1 + 401);
    const auto *sgi1_406 = buffer.data(sgi1 + 406);
    const auto *sgi1_413 = buffer.data(sgi1 + 413);
    const auto *sgi1_415 = buffer.data(sgi1 + 415);
    const auto *sgi1_416 = buffer.data(sgi1 + 416);
    const auto *sgi1_417 = buffer.data(sgi1 + 417);
    const auto *sgi1_419 = buffer.data(sgi1 + 419);

    const auto *pdi0_217 = buffer.data(pdi0 + 217);

    const auto *pdi1_217 = buffer.data(pdi1 + 217);

    const auto *pfi0_252 = buffer.data(pfi0 + 252);
    const auto *pfi0_257 = buffer.data(pfi0 + 257);
    const auto *pfi0_261 = buffer.data(pfi0 + 261);
    const auto *pfi0_266 = buffer.data(pfi0 + 266);
    const auto *pfi0_281 = buffer.data(pfi0 + 281);
    const auto *pfi0_283 = buffer.data(pfi0 + 283);
    const auto *pfi0_329 = buffer.data(pfi0 + 329);

    const auto *pfh_168 = buffer.data(pfh + 168);
    const auto *pfh_171 = buffer.data(pfh + 171);
    const auto *pfh_174 = buffer.data(pfh + 174);
    const auto *pfh_178 = buffer.data(pfh + 178);
    const auto *pfh_183 = buffer.data(pfh + 183);
    const auto *pfh_188 = buffer.data(pfh + 188);
    const auto *pfh_189 = buffer.data(pfh + 189);
    const auto *pfh_191 = buffer.data(pfh + 191);
    const auto *pfh_192 = buffer.data(pfh + 192);
    const auto *pfh_194 = buffer.data(pfh + 194);
    const auto *pfh_195 = buffer.data(pfh + 195);
    const auto *pfh_198 = buffer.data(pfh + 198);
    const auto *pfh_199 = buffer.data(pfh + 199);
    const auto *pfh_203 = buffer.data(pfh + 203);
    const auto *pfh_204 = buffer.data(pfh + 204);
    const auto *pfh_209 = buffer.data(pfh + 209);
    const auto *pfh_210 = buffer.data(pfh + 210);
    const auto *pfh_211 = buffer.data(pfh + 211);
    const auto *pfh_225 = buffer.data(pfh + 225);
    const auto *pfh_226 = buffer.data(pfh + 226);
    const auto *pfh_227 = buffer.data(pfh + 227);
    const auto *pfh_228 = buffer.data(pfh + 228);
    const auto *pfh_229 = buffer.data(pfh + 229);
    const auto *pfh_230 = buffer.data(pfh + 230);
    const auto *pfh_231 = buffer.data(pfh + 231);
    const auto *pfh_232 = buffer.data(pfh + 232);
    const auto *pfh_234 = buffer.data(pfh + 234);
    const auto *pfh_236 = buffer.data(pfh + 236);
    const auto *pfh_237 = buffer.data(pfh + 237);
    const auto *pfh_239 = buffer.data(pfh + 239);
    const auto *pfh_240 = buffer.data(pfh + 240);
    const auto *pfh_241 = buffer.data(pfh + 241);
    const auto *pfh_243 = buffer.data(pfh + 243);
    const auto *pfh_244 = buffer.data(pfh + 244);
    const auto *pfh_245 = buffer.data(pfh + 245);
    const auto *pfh_246 = buffer.data(pfh + 246);
    const auto *pfh_247 = buffer.data(pfh + 247);
    const auto *pfh_248 = buffer.data(pfh + 248);
    const auto *pfh_249 = buffer.data(pfh + 249);
    const auto *pfh_250 = buffer.data(pfh + 250);
    const auto *pfh_251 = buffer.data(pfh + 251);

    const auto *pfi1_252 = buffer.data(pfi1 + 252);
    const auto *pfi1_257 = buffer.data(pfi1 + 257);
    const auto *pfi1_261 = buffer.data(pfi1 + 261);
    const auto *pfi1_266 = buffer.data(pfi1 + 266);
    const auto *pfi1_281 = buffer.data(pfi1 + 281);
    const auto *pfi1_283 = buffer.data(pfi1 + 283);
    const auto *pfi1_329 = buffer.data(pfi1 + 329);

    const auto *pgg0_210 = buffer.data(pgg0 + 210);
    const auto *pgg0_211 = buffer.data(pgg0 + 211);
    const auto *pgg0_213 = buffer.data(pgg0 + 213);
    const auto *pgg0_215 = buffer.data(pgg0 + 215);
    const auto *pgg0_235 = buffer.data(pgg0 + 235);
    const auto *pgg0_236 = buffer.data(pgg0 + 236);
    const auto *pgg0_237 = buffer.data(pgg0 + 237);
    const auto *pgg0_240 = buffer.data(pgg0 + 240);
    const auto *pgg0_241 = buffer.data(pgg0 + 241);
    const auto *pgg0_243 = buffer.data(pgg0 + 243);
    const auto *pgg0_245 = buffer.data(pgg0 + 245);
    const auto *pgg0_246 = buffer.data(pgg0 + 246);
    const auto *pgg0_248 = buffer.data(pgg0 + 248);
    const auto *pgg0_249 = buffer.data(pgg0 + 249);
    const auto *pgg0_250 = buffer.data(pgg0 + 250);
    const auto *pgg0_251 = buffer.data(pgg0 + 251);
    const auto *pgg0_252 = buffer.data(pgg0 + 252);
    const auto *pgg0_253 = buffer.data(pgg0 + 253);
    const auto *pgg0_254 = buffer.data(pgg0 + 254);

    const auto *pgg1_210 = buffer.data(pgg1 + 210);
    const auto *pgg1_211 = buffer.data(pgg1 + 211);
    const auto *pgg1_213 = buffer.data(pgg1 + 213);
    const auto *pgg1_215 = buffer.data(pgg1 + 215);
    const auto *pgg1_235 = buffer.data(pgg1 + 235);
    const auto *pgg1_236 = buffer.data(pgg1 + 236);
    const auto *pgg1_237 = buffer.data(pgg1 + 237);
    const auto *pgg1_240 = buffer.data(pgg1 + 240);
    const auto *pgg1_241 = buffer.data(pgg1 + 241);
    const auto *pgg1_243 = buffer.data(pgg1 + 243);
    const auto *pgg1_245 = buffer.data(pgg1 + 245);
    const auto *pgg1_246 = buffer.data(pgg1 + 246);
    const auto *pgg1_248 = buffer.data(pgg1 + 248);
    const auto *pgg1_249 = buffer.data(pgg1 + 249);
    const auto *pgg1_250 = buffer.data(pgg1 + 250);
    const auto *pgg1_251 = buffer.data(pgg1 + 251);
    const auto *pgg1_252 = buffer.data(pgg1 + 252);
    const auto *pgg1_253 = buffer.data(pgg1 + 253);
    const auto *pgg1_254 = buffer.data(pgg1 + 254);

    const auto *pgh_272 = buffer.data(pgh + 272);
    const auto *pgh_273 = buffer.data(pgh + 273);
    const auto *pgh_275 = buffer.data(pgh + 275);
    const auto *pgh_276 = buffer.data(pgh + 276);
    const auto *pgh_278 = buffer.data(pgh + 278);
    const auto *pgh_279 = buffer.data(pgh + 279);
    const auto *pgh_282 = buffer.data(pgh + 282);
    const auto *pgh_283 = buffer.data(pgh + 283);
    const auto *pgh_287 = buffer.data(pgh + 287);
    const auto *pgh_288 = buffer.data(pgh + 288);
    const auto *pgh_290 = buffer.data(pgh + 290);
    const auto *pgh_291 = buffer.data(pgh + 291);
    const auto *pgh_293 = buffer.data(pgh + 293);
    const auto *pgh_294 = buffer.data(pgh + 294);
    const auto *pgh_295 = buffer.data(pgh + 295);
    const auto *pgh_296 = buffer.data(pgh + 296);
    const auto *pgh_297 = buffer.data(pgh + 297);
    const auto *pgh_299 = buffer.data(pgh + 299);
    const auto *pgh_300 = buffer.data(pgh + 300);
    const auto *pgh_302 = buffer.data(pgh + 302);
    const auto *pgh_303 = buffer.data(pgh + 303);
    const auto *pgh_304 = buffer.data(pgh + 304);
    const auto *pgh_308 = buffer.data(pgh + 308);
    const auto *pgh_309 = buffer.data(pgh + 309);
    const auto *pgh_311 = buffer.data(pgh + 311);
    const auto *pgh_312 = buffer.data(pgh + 312);
    const auto *pgh_314 = buffer.data(pgh + 314);
    const auto *pgh_315 = buffer.data(pgh + 315);
    const auto *pgh_316 = buffer.data(pgh + 316);
    const auto *pgh_318 = buffer.data(pgh + 318);
    const auto *pgh_321 = buffer.data(pgh + 321);
    const auto *pgh_330 = buffer.data(pgh + 330);
    const auto *pgh_331 = buffer.data(pgh + 331);
    const auto *pgh_332 = buffer.data(pgh + 332);
    const auto *pgh_333 = buffer.data(pgh + 333);
    const auto *pgh_334 = buffer.data(pgh + 334);
    const auto *pgh_335 = buffer.data(pgh + 335);
    const auto *pgh_336 = buffer.data(pgh + 336);
    const auto *pgh_337 = buffer.data(pgh + 337);
    const auto *pgh_339 = buffer.data(pgh + 339);
    const auto *pgh_341 = buffer.data(pgh + 341);
    const auto *pgh_342 = buffer.data(pgh + 342);
    const auto *pgh_344 = buffer.data(pgh + 344);
    const auto *pgh_345 = buffer.data(pgh + 345);
    const auto *pgh_346 = buffer.data(pgh + 346);
    const auto *pgh_348 = buffer.data(pgh + 348);
    const auto *pgh_349 = buffer.data(pgh + 349);
    const auto *pgh_350 = buffer.data(pgh + 350);
    const auto *pgh_351 = buffer.data(pgh + 351);
    const auto *pgh_352 = buffer.data(pgh + 352);
    const auto *pgh_353 = buffer.data(pgh + 353);
    const auto *pgh_354 = buffer.data(pgh + 354);
    const auto *pgh_355 = buffer.data(pgh + 355);
    const auto *pgh_356 = buffer.data(pgh + 356);
    const auto *pgh_357 = buffer.data(pgh + 357);
    const auto *pgh_358 = buffer.data(pgh + 358);

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pa_x, pc_x, pc_y, sgi0_359, sgi0_360, \
                         sgi0_361, sgi1_359, sgi1_360, sgi1_361, pfh_188, \
                         pgh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = pa_x[k] * sgi0_359[k]
                   - f_11 * pc_x[k] * sgi1_359[k];

        t_360[k] = pa_x[k] * sgi0_360[k]
                   - f_11 * pc_x[k] * sgi1_360[k];

        t_361[k] = pa_x[k] * sgi0_361[k]
                   - f_11 * pc_x[k] * sgi1_361[k];

        t_362[k] = f_12 * pfh_188[k]
                   + f_4 * pc_y[k] * pgh_272[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pa_x, pb_y, pc_x, pc_y, pc_z, sgi0_363, \
                         sgi1_363, pfi0_252, pfh_168, pfh_189, pfi1_252, \
                         pgh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pa_x[k] * sgi0_363[k]
                   - f_11 * pc_x[k] * sgi1_363[k];

        t_364[k] = pb_y[k] * pfi0_252[k]
                   - f_11 * pc_y[k] * pfi1_252[k];

        t_365[k] = f_0 * pfh_189[k]
                   + f_4 * pc_y[k] * pgh_273[k];

        t_366[k] = f_13 * pfh_168[k]
                   + f_4 * pc_z[k] * pgh_273[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pa_x, pb_y, pc_x, pc_y, sgi0_367, sgh_276, \
                         sgi1_367, pfi0_257, pfh_191, pfi1_257, \
                         pgh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = pa_x[k] * sgi0_367[k]
                   + f_1 * sgh_276[k]
                   - f_11 * pc_x[k] * sgi1_367[k];

        t_368[k] = f_0 * pfh_191[k]
                   + f_4 * pc_y[k] * pgh_275[k];

        t_369[k] = pb_y[k] * pfi0_257[k]
                   - f_11 * pc_y[k] * pfi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pa_x, pc_x, pc_y, pc_z, sgi0_370, sgh_279, \
                         sgi1_370, pfh_171, pfh_194, pgh_276, pgh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_x[k] * sgi0_370[k]
                   + f_13 * sgh_279[k]
                   - f_11 * pc_x[k] * sgi1_370[k];

        t_371[k] = f_13 * pfh_171[k]
                   + f_4 * pc_z[k] * pgh_276[k];

        t_372[k] = f_0 * pfh_194[k]
                   + f_4 * pc_y[k] * pgh_278[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pa_x, pb_y, pc_x, pc_y, pc_z, sgi0_374, sgh_283, \
                         sgi1_374, pfi0_261, pfh_174, pfi1_261, \
                         pgh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = pb_y[k] * pfi0_261[k]
                   - f_11 * pc_y[k] * pfi1_261[k];

        t_374[k] = pa_x[k] * sgi0_374[k]
                   + f_12 * sgh_283[k]
                   - f_11 * pc_x[k] * sgi1_374[k];

        t_375[k] = f_13 * pfh_174[k]
                   + f_4 * pc_z[k] * pgh_279[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pa_x, pb_y, pc_x, pc_y, sgi0_376, sgh_285, \
                         sgi1_376, pfi0_266, pfh_198, pfi1_266, \
                         pgh_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pa_x[k] * sgi0_376[k]
                   + f_12 * sgh_285[k]
                   - f_11 * pc_x[k] * sgi1_376[k];

        t_377[k] = f_0 * pfh_198[k]
                   + f_4 * pc_y[k] * pgh_282[k];

        t_378[k] = pb_y[k] * pfi0_266[k]
                   - f_11 * pc_y[k] * pfi1_266[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_z, sgh_288, sgh_290, sgh_291, \
                         pfh_178, pgh_283, pgh_288, pgh_290, pgh_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_0 * sgh_288[k]
                   + f_4 * pc_x[k] * pgh_288[k];

        t_380[k] = f_13 * pfh_178[k]
                   + f_4 * pc_z[k] * pgh_283[k];

        t_381[k] = f_0 * sgh_290[k]
                   + f_4 * pc_x[k] * pgh_290[k];

        t_382[k] = f_0 * sgh_291[k]
                   + f_4 * pc_x[k] * pgh_291[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, pa_x, pc_x, pc_y, pc_z, sgi0_385, \
                         sgh_293, sgi1_385, pfh_183, pfh_203, pgh_287, pgh_288, \
                         pgh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_0 * pfh_203[k]
                   + f_4 * pc_y[k] * pgh_287[k];

        t_384[k] = f_0 * sgh_293[k]
                   + f_4 * pc_x[k] * pgh_293[k];

        t_385[k] = pa_x[k] * sgi0_385[k]
                   - f_11 * pc_x[k] * sgi1_385[k];

        t_386[k] = f_13 * pfh_183[k]
                   + f_4 * pc_z[k] * pgh_288[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_x, pc_x, pc_y, sgi0_387, sgi0_388, \
                         sgi0_389, sgi1_387, sgi1_388, sgi1_389, pfh_209, \
                         pgh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = pa_x[k] * sgi0_387[k]
                   - f_11 * pc_x[k] * sgi1_387[k];

        t_388[k] = pa_x[k] * sgi0_388[k]
                   - f_11 * pc_x[k] * sgi1_388[k];

        t_389[k] = pa_x[k] * sgi0_389[k]
                   - f_11 * pc_x[k] * sgi1_389[k];

        t_390[k] = f_0 * pfh_209[k]
                   + f_4 * pc_y[k] * pgh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_x, pc_x, pc_y, pc_z, sgi0_391, \
                         sgi0_392, sgh_294, sgi1_391, sgi1_392, pfh_189, \
                         pgh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_x[k] * sgi0_391[k]
                   - f_11 * pc_x[k] * sgi1_391[k];

        t_392[k] = pa_x[k] * sgi0_392[k]
                   + f_16 * sgh_294[k]
                   - f_11 * pc_x[k] * sgi1_392[k];

        t_393[k] = f_4 * pc_y[k] * pgh_294[k];

        t_394[k] = f_1 * pfh_189[k]
                   + f_4 * pc_z[k] * pgh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pa_x, pc_x, pc_y, sgi0_397, sgh_299, sgi1_397, \
                         pgg0_210, pgg1_210, pgh_295, pgh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_5 * pgg0_210[k]
                   - f_6 * pgg1_210[k]
                   + f_4 * pc_y[k] * pgh_295[k];

        t_396[k] = f_4 * pc_y[k] * pgh_296[k];

        t_397[k] = pa_x[k] * sgi0_397[k]
                   + f_1 * sgh_299[k]
                   - f_11 * pc_x[k] * sgi1_397[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, pa_x, pc_x, pc_y, pc_z, sgi0_401, \
                         sgh_303, sgi1_401, pfh_192, pgg0_211, pgg1_211, pgh_297, \
                         pgh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_7 * pgg0_211[k]
                   - f_8 * pgg1_211[k]
                   + f_4 * pc_y[k] * pgh_297[k];

        t_399[k] = f_1 * pfh_192[k]
                   + f_4 * pc_z[k] * pgh_297[k];

        t_400[k] = f_4 * pc_y[k] * pgh_299[k];

        t_401[k] = pa_x[k] * sgi0_401[k]
                   + f_13 * sgh_303[k]
                   - f_11 * pc_x[k] * sgi1_401[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, pc_y, pc_z, pfh_195, pgg0_213, pgg0_215, \
                         pgg1_213, pgg1_215, pgh_300, pgh_302, \
                         pgh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_9 * pgg0_213[k]
                   - f_10 * pgg1_213[k]
                   + f_4 * pc_y[k] * pgh_300[k];

        t_403[k] = f_1 * pfh_195[k]
                   + f_4 * pc_z[k] * pgh_300[k];

        t_404[k] = f_5 * pgg0_215[k]
                   - f_6 * pgg1_215[k]
                   + f_4 * pc_y[k] * pgh_302[k];

        t_405[k] = f_4 * pc_y[k] * pgh_303[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pa_x, pc_x, pc_z, sgi0_406, sgh_308, \
                         sgh_309, sgh_311, sgi1_406, pfh_199, pgh_304, pgh_309, \
                         pgh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pa_x[k] * sgi0_406[k]
                   + f_12 * sgh_308[k]
                   - f_11 * pc_x[k] * sgi1_406[k];

        t_407[k] = f_0 * sgh_309[k]
                   + f_4 * pc_x[k] * pgh_309[k];

        t_408[k] = f_1 * pfh_199[k]
                   + f_4 * pc_z[k] * pgh_304[k];

        t_409[k] = f_0 * sgh_311[k]
                   + f_4 * pc_x[k] * pgh_311[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pa_x, pc_x, pc_y, sgi0_413, sgh_312, \
                         sgh_314, sgi1_413, pgh_308, pgh_312, pgh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_0 * sgh_312[k]
                   + f_4 * pc_x[k] * pgh_312[k];

        t_411[k] = f_4 * pc_y[k] * pgh_308[k];

        t_412[k] = f_0 * sgh_314[k]
                   + f_4 * pc_x[k] * pgh_314[k];

        t_413[k] = pa_x[k] * sgi0_413[k]
                   - f_11 * pc_x[k] * sgi1_413[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pa_x, pc_x, pc_z, sgi0_415, sgi0_416, \
                         sgi0_417, sgi1_415, sgi1_416, sgi1_417, pfh_204, \
                         pgh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_1 * pfh_204[k]
                   + f_4 * pc_z[k] * pgh_309[k];

        t_415[k] = pa_x[k] * sgi0_415[k]
                   - f_11 * pc_x[k] * sgi1_415[k];

        t_416[k] = pa_x[k] * sgi0_416[k]
                   - f_11 * pc_x[k] * sgi1_416[k];

        t_417[k] = pa_x[k] * sgi0_417[k]
                   - f_11 * pc_x[k] * sgi1_417[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pa_x, pa_y, pc_x, pc_y, sgi0_0, sgi0_1, \
                         sgi0_419, sgh_0, sgi1_0, sgi1_1, sgi1_419, \
                         pgh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_4 * pc_y[k] * pgh_314[k];

        t_419[k] = pa_x[k] * sgi0_419[k]
                   - f_11 * pc_x[k] * sgi1_419[k];

        t_420[k] = pa_y[k] * sgi0_0[k]
                   - f_11 * pc_y[k] * sgi1_0[k];

        t_421[k] = pa_y[k] * sgi0_1[k]
                   + f_0 * sgh_0[k]
                   - f_11 * pc_y[k] * sgi1_1[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, pa_y, pc_y, pc_z, sgi0_3, sgi0_5, sgh_1, \
                         sgi1_3, sgi1_5, pgh_315, pgh_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_4 * pc_z[k] * pgh_315[k];

        t_423[k] = pa_y[k] * sgi0_3[k]
                   + f_12 * sgh_1[k]
                   - f_11 * pc_y[k] * sgi1_3[k];

        t_424[k] = f_4 * pc_z[k] * pgh_316[k];

        t_425[k] = pa_y[k] * sgi0_5[k]
                   - f_11 * pc_y[k] * sgi1_5[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pa_y, pc_y, pc_z, sgi0_6, sgi0_8, sgi0_9, \
                         sgh_3, sgh_5, sgi1_6, sgi1_8, sgi1_9, \
                         pgh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = pa_y[k] * sgi0_6[k]
                   + f_13 * sgh_3[k]
                   - f_11 * pc_y[k] * sgi1_6[k];

        t_427[k] = f_4 * pc_z[k] * pgh_318[k];

        t_428[k] = pa_y[k] * sgi0_8[k]
                   + f_0 * sgh_5[k]
                   - f_11 * pc_y[k] * sgi1_8[k];

        t_429[k] = pa_y[k] * sgi0_9[k]
                   - f_11 * pc_y[k] * sgi1_9[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pa_y, pc_y, pc_z, sgi0_10, sgi0_12, sgh_6, \
                         sgh_8, sgi1_10, sgi1_12, pgh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = pa_y[k] * sgi0_10[k]
                   + f_1 * sgh_6[k]
                   - f_11 * pc_y[k] * sgi1_10[k];

        t_431[k] = f_4 * pc_z[k] * pgh_321[k];

        t_432[k] = pa_y[k] * sgi0_12[k]
                   + f_12 * sgh_8[k]
                   - f_11 * pc_y[k] * sgi1_12[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pa_y, pc_x, pc_y, sgi0_13, sgi0_14, \
                         sgh_9, sgi1_13, sgi1_14, pfh_225, pfh_226, pgh_330, \
                         pgh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = pa_y[k] * sgi0_13[k]
                   + f_0 * sgh_9[k]
                   - f_11 * pc_y[k] * sgi1_13[k];

        t_434[k] = pa_y[k] * sgi0_14[k]
                   - f_11 * pc_y[k] * sgi1_14[k];

        t_435[k] = f_1 * pfh_225[k]
                   + f_4 * pc_x[k] * pgh_330[k];

        t_436[k] = f_1 * pfh_226[k]
                   + f_4 * pc_x[k] * pgh_331[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_x, pfh_227, pfh_228, pfh_229, pfh_230, \
                         pgh_332, pgh_333, pgh_334, pgh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_1 * pfh_227[k]
                   + f_4 * pc_x[k] * pgh_332[k];

        t_438[k] = f_1 * pfh_228[k]
                   + f_4 * pc_x[k] * pgh_333[k];

        t_439[k] = f_1 * pfh_229[k]
                   + f_4 * pc_x[k] * pgh_334[k];

        t_440[k] = f_1 * pfh_230[k]
                   + f_4 * pc_x[k] * pgh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pa_y, pc_y, pc_z, sgi0_21, sgh_15, sgi1_21, \
                         pgg0_235, pgg1_235, pgh_330, pgh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = pa_y[k] * sgi0_21[k]
                   + f_16 * sgh_15[k]
                   - f_11 * pc_y[k] * sgi1_21[k];

        t_442[k] = f_4 * pc_z[k] * pgh_330[k];

        t_443[k] = f_5 * pgg0_235[k]
                   - f_6 * pgg1_235[k]
                   + f_4 * pc_z[k] * pgh_331[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pc_y, pc_z, sgh_20, pgg0_236, pgg0_237, \
                         pgg1_236, pgg1_237, pgh_332, pgh_333, \
                         pgh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_7 * pgg0_236[k]
                   - f_8 * pgg1_236[k]
                   + f_4 * pc_z[k] * pgh_332[k];

        t_445[k] = f_9 * pgg0_237[k]
                   - f_10 * pgg1_237[k]
                   + f_4 * pc_z[k] * pgh_333[k];

        t_446[k] = f_0 * sgh_20[k]
                   + f_4 * pc_y[k] * pgh_335[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pa_y, pc_x, pc_y, sgi0_27, sgi1_27, pfh_231, \
                         pfh_232, pgg0_240, pgg0_241, pgg1_240, pgg1_241, pgh_336, \
                         pgh_337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = pa_y[k] * sgi0_27[k]
                   - f_11 * pc_y[k] * sgi1_27[k];

        t_448[k] = f_13 * pfh_231[k]
                   + f_2 * pgg0_240[k]
                   - f_3 * pgg1_240[k]
                   + f_4 * pc_x[k] * pgh_336[k];

        t_449[k] = f_13 * pfh_232[k]
                   + f_17 * pgg0_241[k]
                   - f_18 * pgg1_241[k]
                   + f_4 * pc_x[k] * pgh_337[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, pc_x, pc_z, pfh_234, pfh_236, pgg0_243, \
                         pgg0_245, pgg1_243, pgg1_245, pgh_336, pgh_337, pgh_339, \
                         pgh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_4 * pc_z[k] * pgh_336[k];

        t_451[k] = f_13 * pfh_234[k]
                   + f_9 * pgg0_243[k]
                   - f_10 * pgg1_243[k]
                   + f_4 * pc_x[k] * pgh_339[k];

        t_452[k] = f_4 * pc_z[k] * pgh_337[k];

        t_453[k] = f_13 * pfh_236[k]
                   + f_9 * pgg0_245[k]
                   - f_10 * pgg1_245[k]
                   + f_4 * pc_x[k] * pgh_341[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, pc_x, pc_z, pfh_237, pfh_239, pgg0_246, \
                         pgg0_248, pgg1_246, pgg1_248, pgh_339, pgh_342, \
                         pgh_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_13 * pfh_237[k]
                   + f_7 * pgg0_246[k]
                   - f_8 * pgg1_246[k]
                   + f_4 * pc_x[k] * pgh_342[k];

        t_455[k] = f_4 * pc_z[k] * pgh_339[k];

        t_456[k] = f_13 * pfh_239[k]
                   + f_7 * pgg0_248[k]
                   - f_8 * pgg1_248[k]
                   + f_4 * pc_x[k] * pgh_344[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, pc_x, pc_z, pfh_240, pfh_241, pgg0_249, \
                         pgg0_250, pgg1_249, pgg1_250, pgh_342, pgh_345, \
                         pgh_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_13 * pfh_240[k]
                   + f_7 * pgg0_249[k]
                   - f_8 * pgg1_249[k]
                   + f_4 * pc_x[k] * pgh_345[k];

        t_458[k] = f_13 * pfh_241[k]
                   + f_5 * pgg0_250[k]
                   - f_6 * pgg1_250[k]
                   + f_4 * pc_x[k] * pgh_346[k];

        t_459[k] = f_4 * pc_z[k] * pgh_342[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_x, pfh_243, pfh_244, pfh_245, pgg0_252, \
                         pgg0_253, pgg0_254, pgg1_252, pgg1_253, pgg1_254, pgh_348, pgh_349, \
                         pgh_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_13 * pfh_243[k]
                   + f_5 * pgg0_252[k]
                   - f_6 * pgg1_252[k]
                   + f_4 * pc_x[k] * pgh_348[k];

        t_461[k] = f_13 * pfh_244[k]
                   + f_5 * pgg0_253[k]
                   - f_6 * pgg1_253[k]
                   + f_4 * pc_x[k] * pgh_349[k];

        t_462[k] = f_13 * pfh_245[k]
                   + f_5 * pgg0_254[k]
                   - f_6 * pgg1_254[k]
                   + f_4 * pc_x[k] * pgh_350[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, t_467, pc_x, pfh_246, pfh_247, pfh_248, \
                         pfh_249, pfh_250, pgh_351, pgh_352, pgh_353, pgh_354, \
                         pgh_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_13 * pfh_246[k]
                   + f_4 * pc_x[k] * pgh_351[k];

        t_464[k] = f_13 * pfh_247[k]
                   + f_4 * pc_x[k] * pgh_352[k];

        t_465[k] = f_13 * pfh_248[k]
                   + f_4 * pc_x[k] * pgh_353[k];

        t_466[k] = f_13 * pfh_249[k]
                   + f_4 * pc_x[k] * pgh_354[k];

        t_467[k] = f_13 * pfh_250[k]
                   + f_4 * pc_x[k] * pgh_355[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pb_x, pc_x, pc_z, pdi0_217, pdi1_217, pfi0_329, \
                         pfh_251, pfi1_329, pgh_351, pgh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_13 * pfh_251[k]
                   + f_4 * pc_x[k] * pgh_356[k];

        t_469[k] = f_19 * pdi0_217[k]
                   - f_20 * pdi1_217[k]
                   + pb_x[k] * pfi0_329[k]
                   - f_11 * pc_x[k] * pfi1_329[k];

        t_470[k] = f_4 * pc_z[k] * pgh_351[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_z, pgg0_250, pgg0_251, pgg0_252, pgg1_250, \
                         pgg1_251, pgg1_252, pgh_352, pgh_353, \
                         pgh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_5 * pgg0_250[k]
                   - f_6 * pgg1_250[k]
                   + f_4 * pc_z[k] * pgh_352[k];

        t_472[k] = f_7 * pgg0_251[k]
                   - f_8 * pgg1_251[k]
                   + f_4 * pc_z[k] * pgh_353[k];

        t_473[k] = f_9 * pgg0_252[k]
                   - f_10 * pgg1_252[k]
                   + f_4 * pc_z[k] * pgh_354[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pa_y, pc_y, pc_z, sgi0_56, sgh_41, sgi1_56, \
                         pfh_230, pgg0_254, pgg1_254, pgh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_0 * sgh_41[k]
                   + f_0 * pfh_230[k]
                   + f_4 * pc_y[k] * pgh_356[k];

        t_475[k] = f_2 * pgg0_254[k]
                   - f_3 * pgg1_254[k]
                   + f_4 * pc_z[k] * pgh_356[k];

        t_476[k] = pa_y[k] * sgi0_56[k]
                   - f_11 * pc_y[k] * sgi1_56[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pb_z, pc_z, pfi0_281, pfi0_283, pfh_210, \
                         pfh_211, pfi1_281, pfi1_283, pgh_357, \
                         pgh_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = pb_z[k] * pfi0_281[k]
                   - f_11 * pc_z[k] * pfi1_281[k];

        t_478[k] = f_0 * pfh_210[k]
                   + f_4 * pc_z[k] * pgh_357[k];

        t_479[k] = pb_z[k] * pfi0_283[k]
                   - f_11 * pc_z[k] * pfi1_283[k];

        t_480[k] = f_0 * pfh_211[k]
                   + f_4 * pc_z[k] * pgh_358[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgi0, const size_t sgh,
                                                          const size_t sgi1, const size_t pdi0,
                                                          const size_t pdi1, const size_t pfi0,
                                                          const size_t pfh, const size_t pfi1,
                                                          const size_t pgg0, const size_t pgg1,
                                                          const size_t pgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_21 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_61 = buffer.data(sgi0 + 61);
    const auto *sgi0_64 = buffer.data(sgi0 + 64);
    const auto *sgi0_65 = buffer.data(sgi0 + 65);
    const auto *sgi0_68 = buffer.data(sgi0 + 68);
    const auto *sgi0_69 = buffer.data(sgi0 + 69);
    const auto *sgi0_70 = buffer.data(sgi0 + 70);
    const auto *sgi0_79 = buffer.data(sgi0 + 79);
    const auto *sgi0_80 = buffer.data(sgi0 + 80);
    const auto *sgi0_81 = buffer.data(sgi0 + 81);
    const auto *sgi0_83 = buffer.data(sgi0 + 83);
    const auto *sgi0_140 = buffer.data(sgi0 + 140);
    const auto *sgi0_141 = buffer.data(sgi0 + 141);
    const auto *sgi0_143 = buffer.data(sgi0 + 143);
    const auto *sgi0_145 = buffer.data(sgi0 + 145);
    const auto *sgi0_146 = buffer.data(sgi0 + 146);
    const auto *sgi0_148 = buffer.data(sgi0 + 148);
    const auto *sgi0_149 = buffer.data(sgi0 + 149);
    const auto *sgi0_150 = buffer.data(sgi0 + 150);
    const auto *sgi0_152 = buffer.data(sgi0 + 152);
    const auto *sgi0_153 = buffer.data(sgi0 + 153);
    const auto *sgi0_154 = buffer.data(sgi0 + 154);
    const auto *sgi0_161 = buffer.data(sgi0 + 161);
    const auto *sgi0_163 = buffer.data(sgi0 + 163);
    const auto *sgi0_164 = buffer.data(sgi0 + 164);
    const auto *sgi0_165 = buffer.data(sgi0 + 165);
    const auto *sgi0_167 = buffer.data(sgi0 + 167);

    const auto *sgh_47 = buffer.data(sgh + 47);
    const auto *sgh_50 = buffer.data(sgh + 50);
    const auto *sgh_51 = buffer.data(sgh + 51);
    const auto *sgh_59 = buffer.data(sgh + 59);
    const auto *sgh_60 = buffer.data(sgh + 60);
    const auto *sgh_61 = buffer.data(sgh + 61);
    const auto *sgh_62 = buffer.data(sgh + 62);
    const auto *sgh_83 = buffer.data(sgh + 83);
    const auto *sgh_104 = buffer.data(sgh + 104);
    const auto *sgh_105 = buffer.data(sgh + 105);
    const auto *sgh_106 = buffer.data(sgh + 106);
    const auto *sgh_108 = buffer.data(sgh + 108);
    const auto *sgh_110 = buffer.data(sgh + 110);
    const auto *sgh_111 = buffer.data(sgh + 111);
    const auto *sgh_113 = buffer.data(sgh + 113);
    const auto *sgh_114 = buffer.data(sgh + 114);
    const auto *sgh_120 = buffer.data(sgh + 120);
    const auto *sgh_122 = buffer.data(sgh + 122);
    const auto *sgh_123 = buffer.data(sgh + 123);
    const auto *sgh_124 = buffer.data(sgh + 124);
    const auto *sgh_125 = buffer.data(sgh + 125);

    const auto *sgi1_61 = buffer.data(sgi1 + 61);
    const auto *sgi1_64 = buffer.data(sgi1 + 64);
    const auto *sgi1_65 = buffer.data(sgi1 + 65);
    const auto *sgi1_68 = buffer.data(sgi1 + 68);
    const auto *sgi1_69 = buffer.data(sgi1 + 69);
    const auto *sgi1_70 = buffer.data(sgi1 + 70);
    const auto *sgi1_79 = buffer.data(sgi1 + 79);
    const auto *sgi1_80 = buffer.data(sgi1 + 80);
    const auto *sgi1_81 = buffer.data(sgi1 + 81);
    const auto *sgi1_83 = buffer.data(sgi1 + 83);
    const auto *sgi1_140 = buffer.data(sgi1 + 140);
    const auto *sgi1_141 = buffer.data(sgi1 + 141);
    const auto *sgi1_143 = buffer.data(sgi1 + 143);
    const auto *sgi1_145 = buffer.data(sgi1 + 145);
    const auto *sgi1_146 = buffer.data(sgi1 + 146);
    const auto *sgi1_148 = buffer.data(sgi1 + 148);
    const auto *sgi1_149 = buffer.data(sgi1 + 149);
    const auto *sgi1_150 = buffer.data(sgi1 + 150);
    const auto *sgi1_152 = buffer.data(sgi1 + 152);
    const auto *sgi1_153 = buffer.data(sgi1 + 153);
    const auto *sgi1_154 = buffer.data(sgi1 + 154);
    const auto *sgi1_161 = buffer.data(sgi1 + 161);
    const auto *sgi1_163 = buffer.data(sgi1 + 163);
    const auto *sgi1_164 = buffer.data(sgi1 + 164);
    const auto *sgi1_165 = buffer.data(sgi1 + 165);
    const auto *sgi1_167 = buffer.data(sgi1 + 167);

    const auto *pdi0_273 = buffer.data(pdi0 + 273);

    const auto *pdi1_273 = buffer.data(pdi1 + 273);

    const auto *pfi0_286 = buffer.data(pfi0 + 286);
    const auto *pfi0_290 = buffer.data(pfi0 + 290);
    const auto *pfi0_301 = buffer.data(pfi0 + 301);
    const auto *pfi0_309 = buffer.data(pfi0 + 309);
    const auto *pfi0_311 = buffer.data(pfi0 + 311);
    const auto *pfi0_314 = buffer.data(pfi0 + 314);
    const auto *pfi0_318 = buffer.data(pfi0 + 318);
    const auto *pfi0_329 = buffer.data(pfi0 + 329);
    const auto *pfi0_385 = buffer.data(pfi0 + 385);
    const auto *pfi0_448 = buffer.data(pfi0 + 448);
    const auto *pfi0_449 = buffer.data(pfi0 + 449);
    const auto *pfi0_451 = buffer.data(pfi0 + 451);

    const auto *pfh_213 = buffer.data(pfh + 213);
    const auto *pfh_216 = buffer.data(pfh + 216);
    const auto *pfh_225 = buffer.data(pfh + 225);
    const auto *pfh_231 = buffer.data(pfh + 231);
    const auto *pfh_232 = buffer.data(pfh + 232);
    const auto *pfh_234 = buffer.data(pfh + 234);
    const auto *pfh_237 = buffer.data(pfh + 237);
    const auto *pfh_246 = buffer.data(pfh + 246);
    const auto *pfh_247 = buffer.data(pfh + 247);
    const auto *pfh_248 = buffer.data(pfh + 248);
    const auto *pfh_249 = buffer.data(pfh + 249);
    const auto *pfh_251 = buffer.data(pfh + 251);
    const auto *pfh_252 = buffer.data(pfh + 252);
    const auto *pfh_253 = buffer.data(pfh + 253);
    const auto *pfh_255 = buffer.data(pfh + 255);
    const auto *pfh_258 = buffer.data(pfh + 258);
    const auto *pfh_267 = buffer.data(pfh + 267);
    const auto *pfh_268 = buffer.data(pfh + 268);
    const auto *pfh_269 = buffer.data(pfh + 269);
    const auto *pfh_270 = buffer.data(pfh + 270);
    const auto *pfh_271 = buffer.data(pfh + 271);
    const auto *pfh_272 = buffer.data(pfh + 272);
    const auto *pfh_273 = buffer.data(pfh + 273);
    const auto *pfh_274 = buffer.data(pfh + 274);
    const auto *pfh_276 = buffer.data(pfh + 276);
    const auto *pfh_278 = buffer.data(pfh + 278);
    const auto *pfh_279 = buffer.data(pfh + 279);
    const auto *pfh_281 = buffer.data(pfh + 281);
    const auto *pfh_282 = buffer.data(pfh + 282);
    const auto *pfh_283 = buffer.data(pfh + 283);
    const auto *pfh_285 = buffer.data(pfh + 285);
    const auto *pfh_286 = buffer.data(pfh + 286);
    const auto *pfh_287 = buffer.data(pfh + 287);
    const auto *pfh_288 = buffer.data(pfh + 288);
    const auto *pfh_289 = buffer.data(pfh + 289);
    const auto *pfh_290 = buffer.data(pfh + 290);
    const auto *pfh_291 = buffer.data(pfh + 291);
    const auto *pfh_292 = buffer.data(pfh + 292);
    const auto *pfh_293 = buffer.data(pfh + 293);
    const auto *pfh_294 = buffer.data(pfh + 294);
    const auto *pfh_299 = buffer.data(pfh + 299);
    const auto *pfh_302 = buffer.data(pfh + 302);
    const auto *pfh_303 = buffer.data(pfh + 303);
    const auto *pfh_306 = buffer.data(pfh + 306);
    const auto *pfh_307 = buffer.data(pfh + 307);
    const auto *pfh_308 = buffer.data(pfh + 308);
    const auto *pfh_309 = buffer.data(pfh + 309);
    const auto *pfh_310 = buffer.data(pfh + 310);
    const auto *pfh_311 = buffer.data(pfh + 311);
    const auto *pfh_312 = buffer.data(pfh + 312);
    const auto *pfh_313 = buffer.data(pfh + 313);
    const auto *pfh_314 = buffer.data(pfh + 314);
    const auto *pfh_330 = buffer.data(pfh + 330);
    const auto *pfh_331 = buffer.data(pfh + 331);
    const auto *pfh_332 = buffer.data(pfh + 332);
    const auto *pfh_333 = buffer.data(pfh + 333);
    const auto *pfh_334 = buffer.data(pfh + 334);
    const auto *pfh_335 = buffer.data(pfh + 335);
    const auto *pfh_336 = buffer.data(pfh + 336);
    const auto *pfh_337 = buffer.data(pfh + 337);
    const auto *pfh_339 = buffer.data(pfh + 339);

    const auto *pfi1_286 = buffer.data(pfi1 + 286);
    const auto *pfi1_290 = buffer.data(pfi1 + 290);
    const auto *pfi1_301 = buffer.data(pfi1 + 301);
    const auto *pfi1_309 = buffer.data(pfi1 + 309);
    const auto *pfi1_311 = buffer.data(pfi1 + 311);
    const auto *pfi1_314 = buffer.data(pfi1 + 314);
    const auto *pfi1_318 = buffer.data(pfi1 + 318);
    const auto *pfi1_329 = buffer.data(pfi1 + 329);
    const auto *pfi1_385 = buffer.data(pfi1 + 385);
    const auto *pfi1_448 = buffer.data(pfi1 + 448);
    const auto *pfi1_449 = buffer.data(pfi1 + 449);
    const auto *pfi1_451 = buffer.data(pfi1 + 451);

    const auto *pgg0_270 = buffer.data(pgg0 + 270);
    const auto *pgg0_271 = buffer.data(pgg0 + 271);
    const auto *pgg0_273 = buffer.data(pgg0 + 273);
    const auto *pgg0_275 = buffer.data(pgg0 + 275);
    const auto *pgg0_276 = buffer.data(pgg0 + 276);
    const auto *pgg0_278 = buffer.data(pgg0 + 278);
    const auto *pgg0_279 = buffer.data(pgg0 + 279);
    const auto *pgg0_280 = buffer.data(pgg0 + 280);
    const auto *pgg0_281 = buffer.data(pgg0 + 281);
    const auto *pgg0_282 = buffer.data(pgg0 + 282);
    const auto *pgg0_283 = buffer.data(pgg0 + 283);
    const auto *pgg0_284 = buffer.data(pgg0 + 284);
    const auto *pgg0_285 = buffer.data(pgg0 + 285);
    const auto *pgg0_290 = buffer.data(pgg0 + 290);
    const auto *pgg0_293 = buffer.data(pgg0 + 293);
    const auto *pgg0_294 = buffer.data(pgg0 + 294);
    const auto *pgg0_295 = buffer.data(pgg0 + 295);
    const auto *pgg0_296 = buffer.data(pgg0 + 296);
    const auto *pgg0_297 = buffer.data(pgg0 + 297);
    const auto *pgg0_298 = buffer.data(pgg0 + 298);
    const auto *pgg0_299 = buffer.data(pgg0 + 299);

    const auto *pgg1_270 = buffer.data(pgg1 + 270);
    const auto *pgg1_271 = buffer.data(pgg1 + 271);
    const auto *pgg1_273 = buffer.data(pgg1 + 273);
    const auto *pgg1_275 = buffer.data(pgg1 + 275);
    const auto *pgg1_276 = buffer.data(pgg1 + 276);
    const auto *pgg1_278 = buffer.data(pgg1 + 278);
    const auto *pgg1_279 = buffer.data(pgg1 + 279);
    const auto *pgg1_280 = buffer.data(pgg1 + 280);
    const auto *pgg1_281 = buffer.data(pgg1 + 281);
    const auto *pgg1_282 = buffer.data(pgg1 + 282);
    const auto *pgg1_283 = buffer.data(pgg1 + 283);
    const auto *pgg1_284 = buffer.data(pgg1 + 284);
    const auto *pgg1_285 = buffer.data(pgg1 + 285);
    const auto *pgg1_290 = buffer.data(pgg1 + 290);
    const auto *pgg1_293 = buffer.data(pgg1 + 293);
    const auto *pgg1_294 = buffer.data(pgg1 + 294);
    const auto *pgg1_295 = buffer.data(pgg1 + 295);
    const auto *pgg1_296 = buffer.data(pgg1 + 296);
    const auto *pgg1_297 = buffer.data(pgg1 + 297);
    const auto *pgg1_298 = buffer.data(pgg1 + 298);
    const auto *pgg1_299 = buffer.data(pgg1 + 299);

    const auto *pgh_360 = buffer.data(pgh + 360);
    const auto *pgh_363 = buffer.data(pgh + 363);
    const auto *pgh_372 = buffer.data(pgh + 372);
    const auto *pgh_373 = buffer.data(pgh + 373);
    const auto *pgh_374 = buffer.data(pgh + 374);
    const auto *pgh_375 = buffer.data(pgh + 375);
    const auto *pgh_376 = buffer.data(pgh + 376);
    const auto *pgh_377 = buffer.data(pgh + 377);
    const auto *pgh_378 = buffer.data(pgh + 378);
    const auto *pgh_379 = buffer.data(pgh + 379);
    const auto *pgh_381 = buffer.data(pgh + 381);
    const auto *pgh_383 = buffer.data(pgh + 383);
    const auto *pgh_384 = buffer.data(pgh + 384);
    const auto *pgh_386 = buffer.data(pgh + 386);
    const auto *pgh_387 = buffer.data(pgh + 387);
    const auto *pgh_388 = buffer.data(pgh + 388);
    const auto *pgh_390 = buffer.data(pgh + 390);
    const auto *pgh_391 = buffer.data(pgh + 391);
    const auto *pgh_392 = buffer.data(pgh + 392);
    const auto *pgh_393 = buffer.data(pgh + 393);
    const auto *pgh_394 = buffer.data(pgh + 394);
    const auto *pgh_395 = buffer.data(pgh + 395);
    const auto *pgh_396 = buffer.data(pgh + 396);
    const auto *pgh_397 = buffer.data(pgh + 397);
    const auto *pgh_398 = buffer.data(pgh + 398);
    const auto *pgh_399 = buffer.data(pgh + 399);
    const auto *pgh_400 = buffer.data(pgh + 400);
    const auto *pgh_402 = buffer.data(pgh + 402);
    const auto *pgh_404 = buffer.data(pgh + 404);
    const auto *pgh_405 = buffer.data(pgh + 405);
    const auto *pgh_407 = buffer.data(pgh + 407);
    const auto *pgh_408 = buffer.data(pgh + 408);
    const auto *pgh_411 = buffer.data(pgh + 411);
    const auto *pgh_412 = buffer.data(pgh + 412);
    const auto *pgh_413 = buffer.data(pgh + 413);
    const auto *pgh_414 = buffer.data(pgh + 414);
    const auto *pgh_415 = buffer.data(pgh + 415);
    const auto *pgh_416 = buffer.data(pgh + 416);
    const auto *pgh_417 = buffer.data(pgh + 417);
    const auto *pgh_418 = buffer.data(pgh + 418);
    const auto *pgh_419 = buffer.data(pgh + 419);
    const auto *pgh_420 = buffer.data(pgh + 420);
    const auto *pgh_421 = buffer.data(pgh + 421);
    const auto *pgh_423 = buffer.data(pgh + 423);
    const auto *pgh_426 = buffer.data(pgh + 426);
    const auto *pgh_435 = buffer.data(pgh + 435);
    const auto *pgh_436 = buffer.data(pgh + 436);
    const auto *pgh_437 = buffer.data(pgh + 437);
    const auto *pgh_438 = buffer.data(pgh + 438);
    const auto *pgh_439 = buffer.data(pgh + 439);
    const auto *pgh_440 = buffer.data(pgh + 440);
    const auto *pgh_441 = buffer.data(pgh + 441);
    const auto *pgh_442 = buffer.data(pgh + 442);

#pragma omp simd aligned(t_481, t_482, t_483, pa_y, pb_z, pc_y, pc_z, sgi0_61, sgi1_61, \
                         pfi0_286, pfh_213, pfi1_286, pgh_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = pa_y[k] * sgi0_61[k]
                   - f_11 * pc_y[k] * sgi1_61[k];

        t_482[k] = pb_z[k] * pfi0_286[k]
                   - f_11 * pc_z[k] * pfi1_286[k];

        t_483[k] = f_0 * pfh_213[k]
                   + f_4 * pc_z[k] * pgh_360[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, pa_y, pb_z, pc_y, pc_z, sgi0_64, sgi0_65, \
                         sgh_47, sgi1_64, sgi1_65, pfi0_290, pfi1_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = pa_y[k] * sgi0_64[k]
                   + f_0 * sgh_47[k]
                   - f_11 * pc_y[k] * sgi1_64[k];

        t_485[k] = pa_y[k] * sgi0_65[k]
                   - f_11 * pc_y[k] * sgi1_65[k];

        t_486[k] = pb_z[k] * pfi0_290[k]
                   - f_11 * pc_z[k] * pfi1_290[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, pa_y, pc_y, pc_z, sgi0_68, sgi0_69, sgh_50, \
                         sgh_51, sgi1_68, sgi1_69, pfh_216, pgh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_0 * pfh_216[k]
                   + f_4 * pc_z[k] * pgh_363[k];

        t_488[k] = pa_y[k] * sgi0_68[k]
                   + f_12 * sgh_50[k]
                   - f_11 * pc_y[k] * sgi1_68[k];

        t_489[k] = pa_y[k] * sgi0_69[k]
                   + f_0 * sgh_51[k]
                   - f_11 * pc_y[k] * sgi1_69[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, pa_y, pc_x, pc_y, sgi0_70, sgi1_70, \
                         pfh_267, pfh_268, pfh_269, pgh_372, pgh_373, \
                         pgh_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = pa_y[k] * sgi0_70[k]
                   - f_11 * pc_y[k] * sgi1_70[k];

        t_491[k] = f_13 * pfh_267[k]
                   + f_4 * pc_x[k] * pgh_372[k];

        t_492[k] = f_13 * pfh_268[k]
                   + f_4 * pc_x[k] * pgh_373[k];

        t_493[k] = f_13 * pfh_269[k]
                   + f_4 * pc_x[k] * pgh_374[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pb_z, pc_x, pc_z, pfi0_301, pfh_270, \
                         pfh_271, pfh_272, pfi1_301, pgh_375, pgh_376, \
                         pgh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_13 * pfh_270[k]
                   + f_4 * pc_x[k] * pgh_375[k];

        t_495[k] = f_13 * pfh_271[k]
                   + f_4 * pc_x[k] * pgh_376[k];

        t_496[k] = f_13 * pfh_272[k]
                   + f_4 * pc_x[k] * pgh_377[k];

        t_497[k] = pb_z[k] * pfi0_301[k]
                   - f_11 * pc_z[k] * pfi1_301[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pa_y, pc_y, pc_z, sgi0_79, sgi0_80, sgh_59, \
                         sgh_60, sgi1_79, sgi1_80, pfh_225, pgh_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_0 * pfh_225[k]
                   + f_4 * pc_z[k] * pgh_372[k];

        t_499[k] = pa_y[k] * sgi0_79[k]
                   + f_1 * sgh_59[k]
                   - f_11 * pc_y[k] * sgi1_79[k];

        t_500[k] = pa_y[k] * sgi0_80[k]
                   + f_13 * sgh_60[k]
                   - f_11 * pc_y[k] * sgi1_80[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pa_y, pc_y, sgi0_81, sgi0_83, sgh_61, sgh_62, \
                         sgi1_81, sgi1_83, pgh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = pa_y[k] * sgi0_81[k]
                   + f_12 * sgh_61[k]
                   - f_11 * pc_y[k] * sgi1_81[k];

        t_502[k] = f_0 * sgh_62[k]
                   + f_4 * pc_y[k] * pgh_377[k];

        t_503[k] = pa_y[k] * sgi0_83[k]
                   - f_11 * pc_y[k] * sgi1_83[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pc_x, pc_z, pfh_273, pfh_274, pgg0_270, \
                         pgg0_271, pgg1_270, pgg1_271, pgh_378, \
                         pgh_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_12 * pfh_273[k]
                   + f_2 * pgg0_270[k]
                   - f_3 * pgg1_270[k]
                   + f_4 * pc_x[k] * pgh_378[k];

        t_505[k] = f_12 * pfh_274[k]
                   + f_17 * pgg0_271[k]
                   - f_18 * pgg1_271[k]
                   + f_4 * pc_x[k] * pgh_379[k];

        t_506[k] = f_4 * pc_z[k] * pgh_378[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, pc_x, pc_z, pfh_276, pfh_278, pgg0_273, \
                         pgg0_275, pgg1_273, pgg1_275, pgh_379, pgh_381, \
                         pgh_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_12 * pfh_276[k]
                   + f_9 * pgg0_273[k]
                   - f_10 * pgg1_273[k]
                   + f_4 * pc_x[k] * pgh_381[k];

        t_508[k] = f_4 * pc_z[k] * pgh_379[k];

        t_509[k] = f_12 * pfh_278[k]
                   + f_9 * pgg0_275[k]
                   - f_10 * pgg1_275[k]
                   + f_4 * pc_x[k] * pgh_383[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, pc_x, pc_z, pfh_279, pfh_281, pgg0_276, \
                         pgg0_278, pgg1_276, pgg1_278, pgh_381, pgh_384, \
                         pgh_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = f_12 * pfh_279[k]
                   + f_7 * pgg0_276[k]
                   - f_8 * pgg1_276[k]
                   + f_4 * pc_x[k] * pgh_384[k];

        t_511[k] = f_4 * pc_z[k] * pgh_381[k];

        t_512[k] = f_12 * pfh_281[k]
                   + f_7 * pgg0_278[k]
                   - f_8 * pgg1_278[k]
                   + f_4 * pc_x[k] * pgh_386[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_z, pfh_282, pfh_283, pgg0_279, \
                         pgg0_280, pgg1_279, pgg1_280, pgh_384, pgh_387, \
                         pgh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_12 * pfh_282[k]
                   + f_7 * pgg0_279[k]
                   - f_8 * pgg1_279[k]
                   + f_4 * pc_x[k] * pgh_387[k];

        t_514[k] = f_12 * pfh_283[k]
                   + f_5 * pgg0_280[k]
                   - f_6 * pgg1_280[k]
                   + f_4 * pc_x[k] * pgh_388[k];

        t_515[k] = f_4 * pc_z[k] * pgh_384[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pc_x, pfh_285, pfh_286, pfh_287, pgg0_282, \
                         pgg0_283, pgg0_284, pgg1_282, pgg1_283, pgg1_284, pgh_390, pgh_391, \
                         pgh_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_12 * pfh_285[k]
                   + f_5 * pgg0_282[k]
                   - f_6 * pgg1_282[k]
                   + f_4 * pc_x[k] * pgh_390[k];

        t_517[k] = f_12 * pfh_286[k]
                   + f_5 * pgg0_283[k]
                   - f_6 * pgg1_283[k]
                   + f_4 * pc_x[k] * pgh_391[k];

        t_518[k] = f_12 * pfh_287[k]
                   + f_5 * pgg0_284[k]
                   - f_6 * pgg1_284[k]
                   + f_4 * pc_x[k] * pgh_392[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, t_523, pc_x, pfh_288, pfh_289, pfh_290, \
                         pfh_291, pfh_292, pgh_393, pgh_394, pgh_395, pgh_396, \
                         pgh_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_12 * pfh_288[k]
                   + f_4 * pc_x[k] * pgh_393[k];

        t_520[k] = f_12 * pfh_289[k]
                   + f_4 * pc_x[k] * pgh_394[k];

        t_521[k] = f_12 * pfh_290[k]
                   + f_4 * pc_x[k] * pgh_395[k];

        t_522[k] = f_12 * pfh_291[k]
                   + f_4 * pc_x[k] * pgh_396[k];

        t_523[k] = f_12 * pfh_292[k]
                   + f_4 * pc_x[k] * pgh_397[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, pb_x, pc_x, pc_z, pdi0_273, pdi1_273, pfi0_385, \
                         pfh_293, pfi1_385, pgh_393, pgh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_12 * pfh_293[k]
                   + f_4 * pc_x[k] * pgh_398[k];

        t_525[k] = f_14 * pdi0_273[k]
                   - f_15 * pdi1_273[k]
                   + pb_x[k] * pfi0_385[k]
                   - f_11 * pc_x[k] * pfi1_385[k];

        t_526[k] = f_4 * pc_z[k] * pgh_393[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pc_z, pgg0_280, pgg0_281, pgg0_282, pgg1_280, \
                         pgg1_281, pgg1_282, pgh_394, pgh_395, \
                         pgh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_5 * pgg0_280[k]
                   - f_6 * pgg1_280[k]
                   + f_4 * pc_z[k] * pgh_394[k];

        t_528[k] = f_7 * pgg0_281[k]
                   - f_8 * pgg1_281[k]
                   + f_4 * pc_z[k] * pgh_395[k];

        t_529[k] = f_9 * pgg0_282[k]
                   - f_10 * pgg1_282[k]
                   + f_4 * pc_z[k] * pgh_396[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_y, pc_z, sgh_83, pfh_251, pfh_294, \
                         pgg0_284, pgg0_285, pgg1_284, pgg1_285, pgh_398, \
                         pgh_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_0 * sgh_83[k]
                   + f_12 * pfh_251[k]
                   + f_4 * pc_y[k] * pgh_398[k];

        t_531[k] = f_2 * pgg0_284[k]
                   - f_3 * pgg1_284[k]
                   + f_4 * pc_z[k] * pgh_398[k];

        t_532[k] = f_12 * pfh_294[k]
                   + f_2 * pgg0_285[k]
                   - f_3 * pgg1_285[k]
                   + f_4 * pc_x[k] * pgh_399[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pb_z, pc_z, pfi0_309, pfi0_311, pfh_231, \
                         pfh_232, pfi1_309, pfi1_311, pgh_399, \
                         pgh_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = pb_z[k] * pfi0_309[k]
                   - f_11 * pc_z[k] * pfi1_309[k];

        t_534[k] = f_0 * pfh_231[k]
                   + f_4 * pc_z[k] * pgh_399[k];

        t_535[k] = pb_z[k] * pfi0_311[k]
                   - f_11 * pc_z[k] * pfi1_311[k];

        t_536[k] = f_0 * pfh_232[k]
                   + f_4 * pc_z[k] * pgh_400[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pb_z, pc_x, pc_z, pfi0_314, pfh_234, pfh_299, \
                         pfi1_314, pgg0_290, pgg1_290, pgh_402, \
                         pgh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_12 * pfh_299[k]
                   + f_9 * pgg0_290[k]
                   - f_10 * pgg1_290[k]
                   + f_4 * pc_x[k] * pgh_404[k];

        t_538[k] = pb_z[k] * pfi0_314[k]
                   - f_11 * pc_z[k] * pfi1_314[k];

        t_539[k] = f_0 * pfh_234[k]
                   + f_4 * pc_z[k] * pgh_402[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pb_z, pc_x, pc_z, pfi0_318, pfh_302, pfh_303, \
                         pfi1_318, pgg0_293, pgg0_294, pgg1_293, pgg1_294, pgh_407, \
                         pgh_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_12 * pfh_302[k]
                   + f_7 * pgg0_293[k]
                   - f_8 * pgg1_293[k]
                   + f_4 * pc_x[k] * pgh_407[k];

        t_541[k] = f_12 * pfh_303[k]
                   + f_7 * pgg0_294[k]
                   - f_8 * pgg1_294[k]
                   + f_4 * pc_x[k] * pgh_408[k];

        t_542[k] = pb_z[k] * pfi0_318[k]
                   - f_11 * pc_z[k] * pfi1_318[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pc_x, pc_z, pfh_237, pfh_306, pfh_307, pgg0_297, \
                         pgg0_298, pgg1_297, pgg1_298, pgh_405, pgh_411, \
                         pgh_412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_0 * pfh_237[k]
                   + f_4 * pc_z[k] * pgh_405[k];

        t_544[k] = f_12 * pfh_306[k]
                   + f_5 * pgg0_297[k]
                   - f_6 * pgg1_297[k]
                   + f_4 * pc_x[k] * pgh_411[k];

        t_545[k] = f_12 * pfh_307[k]
                   + f_5 * pgg0_298[k]
                   - f_6 * pgg1_298[k]
                   + f_4 * pc_x[k] * pgh_412[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, pc_x, pfh_308, pfh_309, pfh_310, pfh_311, \
                         pgg0_299, pgg1_299, pgh_413, pgh_414, pgh_415, \
                         pgh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = f_12 * pfh_308[k]
                   + f_5 * pgg0_299[k]
                   - f_6 * pgg1_299[k]
                   + f_4 * pc_x[k] * pgh_413[k];

        t_547[k] = f_12 * pfh_309[k]
                   + f_4 * pc_x[k] * pgh_414[k];

        t_548[k] = f_12 * pfh_310[k]
                   + f_4 * pc_x[k] * pgh_415[k];

        t_549[k] = f_12 * pfh_311[k]
                   + f_4 * pc_x[k] * pgh_416[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, pb_z, pc_x, pc_z, pfi0_329, pfh_312, \
                         pfh_313, pfh_314, pfi1_329, pgh_417, pgh_418, \
                         pgh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_12 * pfh_312[k]
                   + f_4 * pc_x[k] * pgh_417[k];

        t_551[k] = f_12 * pfh_313[k]
                   + f_4 * pc_x[k] * pgh_418[k];

        t_552[k] = f_12 * pfh_314[k]
                   + f_4 * pc_x[k] * pgh_419[k];

        t_553[k] = pb_z[k] * pfi0_329[k]
                   - f_11 * pc_z[k] * pfi1_329[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pc_z, pfh_246, pfh_247, pfh_248, pgg0_295, \
                         pgg0_296, pgg1_295, pgg1_296, pgh_414, pgh_415, \
                         pgh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_0 * pfh_246[k]
                   + f_4 * pc_z[k] * pgh_414[k];

        t_555[k] = f_0 * pfh_247[k]
                   + f_5 * pgg0_295[k]
                   - f_6 * pgg1_295[k]
                   + f_4 * pc_z[k] * pgh_415[k];

        t_556[k] = f_0 * pfh_248[k]
                   + f_7 * pgg0_296[k]
                   - f_8 * pgg1_296[k]
                   + f_4 * pc_z[k] * pgh_416[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_y, pc_z, sgh_104, pfh_249, pfh_251, pfh_272, \
                         pgg0_297, pgg0_299, pgg1_297, pgg1_299, pgh_417, \
                         pgh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_0 * pfh_249[k]
                   + f_9 * pgg0_297[k]
                   - f_10 * pgg1_297[k]
                   + f_4 * pc_z[k] * pgh_417[k];

        t_558[k] = f_0 * sgh_104[k]
                   + f_0 * pfh_272[k]
                   + f_4 * pc_y[k] * pgh_419[k];

        t_559[k] = f_0 * pfh_251[k]
                   + f_2 * pgg0_299[k]
                   - f_3 * pgg1_299[k]
                   + f_4 * pc_z[k] * pgh_419[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pa_y, pc_y, pc_z, sgi0_140, sgi0_141, sgh_105, \
                         sgi1_140, sgi1_141, pfh_252, pgh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = pa_y[k] * sgi0_140[k]
                   - f_11 * pc_y[k] * sgi1_140[k];

        t_561[k] = pa_y[k] * sgi0_141[k]
                   + f_0 * sgh_105[k]
                   - f_11 * pc_y[k] * sgi1_141[k];

        t_562[k] = f_12 * pfh_252[k]
                   + f_4 * pc_z[k] * pgh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pa_y, pc_y, pc_z, sgi0_143, sgi0_145, sgh_106, \
                         sgi1_143, sgi1_145, pfh_253, pgh_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = pa_y[k] * sgi0_143[k]
                   + f_12 * sgh_106[k]
                   - f_11 * pc_y[k] * sgi1_143[k];

        t_564[k] = f_12 * pfh_253[k]
                   + f_4 * pc_z[k] * pgh_421[k];

        t_565[k] = pa_y[k] * sgi0_145[k]
                   - f_11 * pc_y[k] * sgi1_145[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pa_y, pc_y, pc_z, sgi0_146, sgi0_148, sgh_108, \
                         sgh_110, sgi1_146, sgi1_148, pfh_255, \
                         pgh_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = pa_y[k] * sgi0_146[k]
                   + f_13 * sgh_108[k]
                   - f_11 * pc_y[k] * sgi1_146[k];

        t_567[k] = f_12 * pfh_255[k]
                   + f_4 * pc_z[k] * pgh_423[k];

        t_568[k] = pa_y[k] * sgi0_148[k]
                   + f_0 * sgh_110[k]
                   - f_11 * pc_y[k] * sgi1_148[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pa_y, pc_y, pc_z, sgi0_149, sgi0_150, sgh_111, \
                         sgi1_149, sgi1_150, pfh_258, pgh_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = pa_y[k] * sgi0_149[k]
                   - f_11 * pc_y[k] * sgi1_149[k];

        t_570[k] = pa_y[k] * sgi0_150[k]
                   + f_1 * sgh_111[k]
                   - f_11 * pc_y[k] * sgi1_150[k];

        t_571[k] = f_12 * pfh_258[k]
                   + f_4 * pc_z[k] * pgh_426[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pa_y, pc_y, sgi0_152, sgi0_153, sgi0_154, \
                         sgh_113, sgh_114, sgi1_152, sgi1_153, \
                         sgi1_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = pa_y[k] * sgi0_152[k]
                   + f_12 * sgh_113[k]
                   - f_11 * pc_y[k] * sgi1_152[k];

        t_573[k] = pa_y[k] * sgi0_153[k]
                   + f_0 * sgh_114[k]
                   - f_11 * pc_y[k] * sgi1_153[k];

        t_574[k] = pa_y[k] * sgi0_154[k]
                   - f_11 * pc_y[k] * sgi1_154[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, pc_x, pfh_330, pfh_331, pfh_332, \
                         pfh_333, pfh_334, pgh_435, pgh_436, pgh_437, pgh_438, \
                         pgh_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_12 * pfh_330[k]
                   + f_4 * pc_x[k] * pgh_435[k];

        t_576[k] = f_12 * pfh_331[k]
                   + f_4 * pc_x[k] * pgh_436[k];

        t_577[k] = f_12 * pfh_332[k]
                   + f_4 * pc_x[k] * pgh_437[k];

        t_578[k] = f_12 * pfh_333[k]
                   + f_4 * pc_x[k] * pgh_438[k];

        t_579[k] = f_12 * pfh_334[k]
                   + f_4 * pc_x[k] * pgh_439[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, pa_y, pc_x, pc_y, pc_z, sgi0_161, sgh_120, \
                         sgi1_161, pfh_267, pfh_335, pgh_435, pgh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_12 * pfh_335[k]
                   + f_4 * pc_x[k] * pgh_440[k];

        t_581[k] = pa_y[k] * sgi0_161[k]
                   + f_16 * sgh_120[k]
                   - f_11 * pc_y[k] * sgi1_161[k];

        t_582[k] = f_12 * pfh_267[k]
                   + f_4 * pc_z[k] * pgh_435[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, pa_y, pc_y, sgi0_163, sgi0_164, sgi0_165, \
                         sgh_122, sgh_123, sgh_124, sgi1_163, sgi1_164, \
                         sgi1_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = pa_y[k] * sgi0_163[k]
                   + f_1 * sgh_122[k]
                   - f_11 * pc_y[k] * sgi1_163[k];

        t_584[k] = pa_y[k] * sgi0_164[k]
                   + f_13 * sgh_123[k]
                   - f_11 * pc_y[k] * sgi1_164[k];

        t_585[k] = pa_y[k] * sgi0_165[k]
                   + f_12 * sgh_124[k]
                   - f_11 * pc_y[k] * sgi1_165[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, pa_y, pb_x, pc_x, pc_y, sgi0_167, sgh_125, \
                         sgi1_167, pfi0_448, pfh_336, pfi1_448, \
                         pgh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_0 * sgh_125[k]
                   + f_4 * pc_y[k] * pgh_440[k];

        t_587[k] = pa_y[k] * sgi0_167[k]
                   - f_11 * pc_y[k] * sgi1_167[k];

        t_588[k] = pb_x[k] * pfi0_448[k]
                   + f_16 * pfh_336[k]
                   - f_11 * pc_x[k] * pfi1_448[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, pb_x, pc_x, pc_z, pfi0_449, pfi0_451, \
                         pfh_337, pfh_339, pfi1_449, pfi1_451, pgh_441, \
                         pgh_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = pb_x[k] * pfi0_449[k]
                   + f_21 * pfh_337[k]
                   - f_11 * pc_x[k] * pfi1_449[k];

        t_590[k] = f_4 * pc_z[k] * pgh_441[k];

        t_591[k] = pb_x[k] * pfi0_451[k]
                   + f_1 * pfh_339[k]
                   - f_11 * pc_x[k] * pfi1_451[k];

        t_592[k] = f_4 * pc_z[k] * pgh_442[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgi0, const size_t sgh,
                                                          const size_t sgi1, const size_t pfi0,
                                                          const size_t pfh, const size_t pfi1,
                                                          const size_t pgg0, const size_t pgg1,
                                                          const size_t pgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_21 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_252 = buffer.data(sgi0 + 252);
    const auto *sgi0_253 = buffer.data(sgi0 + 253);
    const auto *sgi0_255 = buffer.data(sgi0 + 255);
    const auto *sgi0_257 = buffer.data(sgi0 + 257);
    const auto *sgi0_258 = buffer.data(sgi0 + 258);
    const auto *sgi0_260 = buffer.data(sgi0 + 260);
    const auto *sgi0_261 = buffer.data(sgi0 + 261);
    const auto *sgi0_262 = buffer.data(sgi0 + 262);
    const auto *sgi0_264 = buffer.data(sgi0 + 264);
    const auto *sgi0_265 = buffer.data(sgi0 + 265);
    const auto *sgi0_266 = buffer.data(sgi0 + 266);
    const auto *sgi0_279 = buffer.data(sgi0 + 279);

    const auto *sgh_189 = buffer.data(sgh + 189);
    const auto *sgh_190 = buffer.data(sgh + 190);
    const auto *sgh_192 = buffer.data(sgh + 192);
    const auto *sgh_194 = buffer.data(sgh + 194);
    const auto *sgh_195 = buffer.data(sgh + 195);
    const auto *sgh_197 = buffer.data(sgh + 197);
    const auto *sgh_198 = buffer.data(sgh + 198);
    const auto *sgh_209 = buffer.data(sgh + 209);

    const auto *sgi1_252 = buffer.data(sgi1 + 252);
    const auto *sgi1_253 = buffer.data(sgi1 + 253);
    const auto *sgi1_255 = buffer.data(sgi1 + 255);
    const auto *sgi1_257 = buffer.data(sgi1 + 257);
    const auto *sgi1_258 = buffer.data(sgi1 + 258);
    const auto *sgi1_260 = buffer.data(sgi1 + 260);
    const auto *sgi1_261 = buffer.data(sgi1 + 261);
    const auto *sgi1_262 = buffer.data(sgi1 + 262);
    const auto *sgi1_264 = buffer.data(sgi1 + 264);
    const auto *sgi1_265 = buffer.data(sgi1 + 265);
    const auto *sgi1_266 = buffer.data(sgi1 + 266);
    const auto *sgi1_279 = buffer.data(sgi1 + 279);

    const auto *pfi0_364 = buffer.data(pfi0 + 364);
    const auto *pfi0_365 = buffer.data(pfi0 + 365);
    const auto *pfi0_367 = buffer.data(pfi0 + 367);
    const auto *pfi0_370 = buffer.data(pfi0 + 370);
    const auto *pfi0_374 = buffer.data(pfi0 + 374);
    const auto *pfi0_453 = buffer.data(pfi0 + 453);
    const auto *pfi0_454 = buffer.data(pfi0 + 454);
    const auto *pfi0_456 = buffer.data(pfi0 + 456);
    const auto *pfi0_457 = buffer.data(pfi0 + 457);
    const auto *pfi0_458 = buffer.data(pfi0 + 458);
    const auto *pfi0_460 = buffer.data(pfi0 + 460);
    const auto *pfi0_461 = buffer.data(pfi0 + 461);
    const auto *pfi0_462 = buffer.data(pfi0 + 462);
    const auto *pfi0_469 = buffer.data(pfi0 + 469);
    const auto *pfi0_471 = buffer.data(pfi0 + 471);
    const auto *pfi0_472 = buffer.data(pfi0 + 472);
    const auto *pfi0_473 = buffer.data(pfi0 + 473);
    const auto *pfi0_474 = buffer.data(pfi0 + 474);
    const auto *pfi0_475 = buffer.data(pfi0 + 475);
    const auto *pfi0_481 = buffer.data(pfi0 + 481);
    const auto *pfi0_484 = buffer.data(pfi0 + 484);
    const auto *pfi0_485 = buffer.data(pfi0 + 485);
    const auto *pfi0_488 = buffer.data(pfi0 + 488);
    const auto *pfi0_489 = buffer.data(pfi0 + 489);
    const auto *pfi0_490 = buffer.data(pfi0 + 490);
    const auto *pfi0_497 = buffer.data(pfi0 + 497);
    const auto *pfi0_499 = buffer.data(pfi0 + 499);
    const auto *pfi0_500 = buffer.data(pfi0 + 500);
    const auto *pfi0_501 = buffer.data(pfi0 + 501);
    const auto *pfi0_502 = buffer.data(pfi0 + 502);
    const auto *pfi0_503 = buffer.data(pfi0 + 503);
    const auto *pfi0_505 = buffer.data(pfi0 + 505);
    const auto *pfi0_507 = buffer.data(pfi0 + 507);
    const auto *pfi0_510 = buffer.data(pfi0 + 510);
    const auto *pfi0_512 = buffer.data(pfi0 + 512);
    const auto *pfi0_514 = buffer.data(pfi0 + 514);
    const auto *pfi0_516 = buffer.data(pfi0 + 516);
    const auto *pfi0_517 = buffer.data(pfi0 + 517);
    const auto *pfi0_525 = buffer.data(pfi0 + 525);
    const auto *pfi0_527 = buffer.data(pfi0 + 527);
    const auto *pfi0_528 = buffer.data(pfi0 + 528);
    const auto *pfi0_529 = buffer.data(pfi0 + 529);
    const auto *pfi0_530 = buffer.data(pfi0 + 530);
    const auto *pfi0_531 = buffer.data(pfi0 + 531);
    const auto *pfi0_553 = buffer.data(pfi0 + 553);
    const auto *pfi0_555 = buffer.data(pfi0 + 555);
    const auto *pfi0_556 = buffer.data(pfi0 + 556);
    const auto *pfi0_557 = buffer.data(pfi0 + 557);

    const auto *pfh_273 = buffer.data(pfh + 273);
    const auto *pfh_274 = buffer.data(pfh + 274);
    const auto *pfh_276 = buffer.data(pfh + 276);
    const auto *pfh_279 = buffer.data(pfh + 279);
    const auto *pfh_288 = buffer.data(pfh + 288);
    const auto *pfh_294 = buffer.data(pfh + 294);
    const auto *pfh_295 = buffer.data(pfh + 295);
    const auto *pfh_297 = buffer.data(pfh + 297);
    const auto *pfh_300 = buffer.data(pfh + 300);
    const auto *pfh_309 = buffer.data(pfh + 309);
    const auto *pfh_315 = buffer.data(pfh + 315);
    const auto *pfh_316 = buffer.data(pfh + 316);
    const auto *pfh_318 = buffer.data(pfh + 318);
    const auto *pfh_321 = buffer.data(pfh + 321);
    const auto *pfh_330 = buffer.data(pfh + 330);
    const auto *pfh_341 = buffer.data(pfh + 341);
    const auto *pfh_342 = buffer.data(pfh + 342);
    const auto *pfh_344 = buffer.data(pfh + 344);
    const auto *pfh_345 = buffer.data(pfh + 345);
    const auto *pfh_346 = buffer.data(pfh + 346);
    const auto *pfh_348 = buffer.data(pfh + 348);
    const auto *pfh_349 = buffer.data(pfh + 349);
    const auto *pfh_350 = buffer.data(pfh + 350);
    const auto *pfh_351 = buffer.data(pfh + 351);
    const auto *pfh_352 = buffer.data(pfh + 352);
    const auto *pfh_353 = buffer.data(pfh + 353);
    const auto *pfh_354 = buffer.data(pfh + 354);
    const auto *pfh_355 = buffer.data(pfh + 355);
    const auto *pfh_356 = buffer.data(pfh + 356);
    const auto *pfh_362 = buffer.data(pfh + 362);
    const auto *pfh_365 = buffer.data(pfh + 365);
    const auto *pfh_366 = buffer.data(pfh + 366);
    const auto *pfh_369 = buffer.data(pfh + 369);
    const auto *pfh_370 = buffer.data(pfh + 370);
    const auto *pfh_371 = buffer.data(pfh + 371);
    const auto *pfh_372 = buffer.data(pfh + 372);
    const auto *pfh_373 = buffer.data(pfh + 373);
    const auto *pfh_374 = buffer.data(pfh + 374);
    const auto *pfh_375 = buffer.data(pfh + 375);
    const auto *pfh_376 = buffer.data(pfh + 376);
    const auto *pfh_377 = buffer.data(pfh + 377);
    const auto *pfh_378 = buffer.data(pfh + 378);
    const auto *pfh_379 = buffer.data(pfh + 379);
    const auto *pfh_381 = buffer.data(pfh + 381);
    const auto *pfh_383 = buffer.data(pfh + 383);
    const auto *pfh_384 = buffer.data(pfh + 384);
    const auto *pfh_386 = buffer.data(pfh + 386);
    const auto *pfh_387 = buffer.data(pfh + 387);
    const auto *pfh_388 = buffer.data(pfh + 388);
    const auto *pfh_390 = buffer.data(pfh + 390);
    const auto *pfh_391 = buffer.data(pfh + 391);
    const auto *pfh_392 = buffer.data(pfh + 392);
    const auto *pfh_393 = buffer.data(pfh + 393);
    const auto *pfh_394 = buffer.data(pfh + 394);
    const auto *pfh_395 = buffer.data(pfh + 395);
    const auto *pfh_396 = buffer.data(pfh + 396);
    const auto *pfh_397 = buffer.data(pfh + 397);
    const auto *pfh_398 = buffer.data(pfh + 398);
    const auto *pfh_414 = buffer.data(pfh + 414);
    const auto *pfh_415 = buffer.data(pfh + 415);
    const auto *pfh_416 = buffer.data(pfh + 416);
    const auto *pfh_417 = buffer.data(pfh + 417);
    const auto *pfh_418 = buffer.data(pfh + 418);
    const auto *pfh_419 = buffer.data(pfh + 419);

    const auto *pfi1_364 = buffer.data(pfi1 + 364);
    const auto *pfi1_365 = buffer.data(pfi1 + 365);
    const auto *pfi1_367 = buffer.data(pfi1 + 367);
    const auto *pfi1_370 = buffer.data(pfi1 + 370);
    const auto *pfi1_374 = buffer.data(pfi1 + 374);
    const auto *pfi1_453 = buffer.data(pfi1 + 453);
    const auto *pfi1_454 = buffer.data(pfi1 + 454);
    const auto *pfi1_456 = buffer.data(pfi1 + 456);
    const auto *pfi1_457 = buffer.data(pfi1 + 457);
    const auto *pfi1_458 = buffer.data(pfi1 + 458);
    const auto *pfi1_460 = buffer.data(pfi1 + 460);
    const auto *pfi1_461 = buffer.data(pfi1 + 461);
    const auto *pfi1_462 = buffer.data(pfi1 + 462);
    const auto *pfi1_469 = buffer.data(pfi1 + 469);
    const auto *pfi1_471 = buffer.data(pfi1 + 471);
    const auto *pfi1_472 = buffer.data(pfi1 + 472);
    const auto *pfi1_473 = buffer.data(pfi1 + 473);
    const auto *pfi1_474 = buffer.data(pfi1 + 474);
    const auto *pfi1_475 = buffer.data(pfi1 + 475);
    const auto *pfi1_481 = buffer.data(pfi1 + 481);
    const auto *pfi1_484 = buffer.data(pfi1 + 484);
    const auto *pfi1_485 = buffer.data(pfi1 + 485);
    const auto *pfi1_488 = buffer.data(pfi1 + 488);
    const auto *pfi1_489 = buffer.data(pfi1 + 489);
    const auto *pfi1_490 = buffer.data(pfi1 + 490);
    const auto *pfi1_497 = buffer.data(pfi1 + 497);
    const auto *pfi1_499 = buffer.data(pfi1 + 499);
    const auto *pfi1_500 = buffer.data(pfi1 + 500);
    const auto *pfi1_501 = buffer.data(pfi1 + 501);
    const auto *pfi1_502 = buffer.data(pfi1 + 502);
    const auto *pfi1_503 = buffer.data(pfi1 + 503);
    const auto *pfi1_505 = buffer.data(pfi1 + 505);
    const auto *pfi1_507 = buffer.data(pfi1 + 507);
    const auto *pfi1_510 = buffer.data(pfi1 + 510);
    const auto *pfi1_512 = buffer.data(pfi1 + 512);
    const auto *pfi1_514 = buffer.data(pfi1 + 514);
    const auto *pfi1_516 = buffer.data(pfi1 + 516);
    const auto *pfi1_517 = buffer.data(pfi1 + 517);
    const auto *pfi1_525 = buffer.data(pfi1 + 525);
    const auto *pfi1_527 = buffer.data(pfi1 + 527);
    const auto *pfi1_528 = buffer.data(pfi1 + 528);
    const auto *pfi1_529 = buffer.data(pfi1 + 529);
    const auto *pfi1_530 = buffer.data(pfi1 + 530);
    const auto *pfi1_531 = buffer.data(pfi1 + 531);
    const auto *pfi1_553 = buffer.data(pfi1 + 553);
    const auto *pfi1_555 = buffer.data(pfi1 + 555);
    const auto *pfi1_556 = buffer.data(pfi1 + 556);
    const auto *pfi1_557 = buffer.data(pfi1 + 557);

    const auto *pgg0_345 = buffer.data(pgg0 + 345);
    const auto *pgg0_350 = buffer.data(pgg0 + 350);
    const auto *pgg0_354 = buffer.data(pgg0 + 354);
    const auto *pgg0_359 = buffer.data(pgg0 + 359);
    const auto *pgg0_375 = buffer.data(pgg0 + 375);
    const auto *pgg0_376 = buffer.data(pgg0 + 376);
    const auto *pgg0_378 = buffer.data(pgg0 + 378);
    const auto *pgg0_380 = buffer.data(pgg0 + 380);
    const auto *pgg0_381 = buffer.data(pgg0 + 381);
    const auto *pgg0_383 = buffer.data(pgg0 + 383);
    const auto *pgg0_384 = buffer.data(pgg0 + 384);
    const auto *pgg0_385 = buffer.data(pgg0 + 385);
    const auto *pgg0_387 = buffer.data(pgg0 + 387);

    const auto *pgg1_345 = buffer.data(pgg1 + 345);
    const auto *pgg1_350 = buffer.data(pgg1 + 350);
    const auto *pgg1_354 = buffer.data(pgg1 + 354);
    const auto *pgg1_359 = buffer.data(pgg1 + 359);
    const auto *pgg1_375 = buffer.data(pgg1 + 375);
    const auto *pgg1_376 = buffer.data(pgg1 + 376);
    const auto *pgg1_378 = buffer.data(pgg1 + 378);
    const auto *pgg1_380 = buffer.data(pgg1 + 380);
    const auto *pgg1_381 = buffer.data(pgg1 + 381);
    const auto *pgg1_383 = buffer.data(pgg1 + 383);
    const auto *pgg1_384 = buffer.data(pgg1 + 384);
    const auto *pgg1_385 = buffer.data(pgg1 + 385);
    const auto *pgg1_387 = buffer.data(pgg1 + 387);

    const auto *pgh_444 = buffer.data(pgh + 444);
    const auto *pgh_447 = buffer.data(pgh + 447);
    const auto *pgh_456 = buffer.data(pgh + 456);
    const auto *pgh_457 = buffer.data(pgh + 457);
    const auto *pgh_458 = buffer.data(pgh + 458);
    const auto *pgh_459 = buffer.data(pgh + 459);
    const auto *pgh_460 = buffer.data(pgh + 460);
    const auto *pgh_461 = buffer.data(pgh + 461);
    const auto *pgh_462 = buffer.data(pgh + 462);
    const auto *pgh_463 = buffer.data(pgh + 463);
    const auto *pgh_465 = buffer.data(pgh + 465);
    const auto *pgh_468 = buffer.data(pgh + 468);
    const auto *pgh_477 = buffer.data(pgh + 477);
    const auto *pgh_478 = buffer.data(pgh + 478);
    const auto *pgh_479 = buffer.data(pgh + 479);
    const auto *pgh_480 = buffer.data(pgh + 480);
    const auto *pgh_481 = buffer.data(pgh + 481);
    const auto *pgh_482 = buffer.data(pgh + 482);
    const auto *pgh_483 = buffer.data(pgh + 483);
    const auto *pgh_484 = buffer.data(pgh + 484);
    const auto *pgh_486 = buffer.data(pgh + 486);
    const auto *pgh_488 = buffer.data(pgh + 488);
    const auto *pgh_489 = buffer.data(pgh + 489);
    const auto *pgh_492 = buffer.data(pgh + 492);
    const auto *pgh_497 = buffer.data(pgh + 497);
    const auto *pgh_498 = buffer.data(pgh + 498);
    const auto *pgh_499 = buffer.data(pgh + 499);
    const auto *pgh_500 = buffer.data(pgh + 500);
    const auto *pgh_501 = buffer.data(pgh + 501);
    const auto *pgh_502 = buffer.data(pgh + 502);
    const auto *pgh_503 = buffer.data(pgh + 503);
    const auto *pgh_504 = buffer.data(pgh + 504);
    const auto *pgh_505 = buffer.data(pgh + 505);
    const auto *pgh_507 = buffer.data(pgh + 507);
    const auto *pgh_510 = buffer.data(pgh + 510);
    const auto *pgh_519 = buffer.data(pgh + 519);
    const auto *pgh_520 = buffer.data(pgh + 520);
    const auto *pgh_521 = buffer.data(pgh + 521);
    const auto *pgh_522 = buffer.data(pgh + 522);
    const auto *pgh_523 = buffer.data(pgh + 523);
    const auto *pgh_524 = buffer.data(pgh + 524);
    const auto *pgh_525 = buffer.data(pgh + 525);
    const auto *pgh_526 = buffer.data(pgh + 526);
    const auto *pgh_528 = buffer.data(pgh + 528);
    const auto *pgh_530 = buffer.data(pgh + 530);
    const auto *pgh_531 = buffer.data(pgh + 531);
    const auto *pgh_533 = buffer.data(pgh + 533);
    const auto *pgh_534 = buffer.data(pgh + 534);
    const auto *pgh_535 = buffer.data(pgh + 535);
    const auto *pgh_537 = buffer.data(pgh + 537);

#pragma omp simd aligned(t_593, t_594, t_595, pb_x, pc_x, pc_z, pfi0_453, pfi0_454, pfh_341, \
                         pfh_342, pfi1_453, pfi1_454, pgh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pb_x[k] * pfi0_453[k]
                   + f_1 * pfh_341[k]
                   - f_11 * pc_x[k] * pfi1_453[k];

        t_594[k] = pb_x[k] * pfi0_454[k]
                   + f_13 * pfh_342[k]
                   - f_11 * pc_x[k] * pfi1_454[k];

        t_595[k] = f_4 * pc_z[k] * pgh_444[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, pb_x, pc_x, pfi0_456, pfi0_457, pfi0_458, \
                         pfh_344, pfh_345, pfh_346, pfi1_456, pfi1_457, \
                         pfi1_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = pb_x[k] * pfi0_456[k]
                   + f_13 * pfh_344[k]
                   - f_11 * pc_x[k] * pfi1_456[k];

        t_597[k] = pb_x[k] * pfi0_457[k]
                   + f_13 * pfh_345[k]
                   - f_11 * pc_x[k] * pfi1_457[k];

        t_598[k] = pb_x[k] * pfi0_458[k]
                   + f_12 * pfh_346[k]
                   - f_11 * pc_x[k] * pfi1_458[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, pb_x, pc_x, pc_z, pfi0_460, pfi0_461, pfh_348, \
                         pfh_349, pfi1_460, pfi1_461, pgh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = f_4 * pc_z[k] * pgh_447[k];

        t_600[k] = pb_x[k] * pfi0_460[k]
                   + f_12 * pfh_348[k]
                   - f_11 * pc_x[k] * pfi1_460[k];

        t_601[k] = pb_x[k] * pfi0_461[k]
                   + f_12 * pfh_349[k]
                   - f_11 * pc_x[k] * pfi1_461[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, pb_x, pc_x, pfi0_462, pfh_350, pfh_351, \
                         pfh_352, pfh_353, pfi1_462, pgh_456, pgh_457, \
                         pgh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = pb_x[k] * pfi0_462[k]
                   + f_12 * pfh_350[k]
                   - f_11 * pc_x[k] * pfi1_462[k];

        t_603[k] = f_0 * pfh_351[k]
                   + f_4 * pc_x[k] * pgh_456[k];

        t_604[k] = f_0 * pfh_352[k]
                   + f_4 * pc_x[k] * pgh_457[k];

        t_605[k] = f_0 * pfh_353[k]
                   + f_4 * pc_x[k] * pgh_458[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, t_609, pb_x, pc_x, pfi0_469, pfh_354, pfh_355, \
                         pfh_356, pfi1_469, pgh_459, pgh_460, pgh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = f_0 * pfh_354[k]
                   + f_4 * pc_x[k] * pgh_459[k];

        t_607[k] = f_0 * pfh_355[k]
                   + f_4 * pc_x[k] * pgh_460[k];

        t_608[k] = f_0 * pfh_356[k]
                   + f_4 * pc_x[k] * pgh_461[k];

        t_609[k] = pb_x[k] * pfi0_469[k]
                   - f_11 * pc_x[k] * pfi1_469[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pb_x, pc_x, pc_z, pfi0_471, pfi0_472, \
                         pfi0_473, pfi1_471, pfi1_472, pfi1_473, \
                         pgh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_4 * pc_z[k] * pgh_456[k];

        t_611[k] = pb_x[k] * pfi0_471[k]
                   - f_11 * pc_x[k] * pfi1_471[k];

        t_612[k] = pb_x[k] * pfi0_472[k]
                   - f_11 * pc_x[k] * pfi1_472[k];

        t_613[k] = pb_x[k] * pfi0_473[k]
                   - f_11 * pc_x[k] * pfi1_473[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, pb_x, pb_z, pc_x, pc_z, pfi0_364, \
                         pfi0_365, pfi0_474, pfi0_475, pfi1_364, pfi1_365, pfi1_474, \
                         pfi1_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = pb_x[k] * pfi0_474[k]
                   - f_11 * pc_x[k] * pfi1_474[k];

        t_615[k] = pb_x[k] * pfi0_475[k]
                   - f_11 * pc_x[k] * pfi1_475[k];

        t_616[k] = pb_z[k] * pfi0_364[k]
                   - f_11 * pc_z[k] * pfi1_364[k];

        t_617[k] = pb_z[k] * pfi0_365[k]
                   - f_11 * pc_z[k] * pfi1_365[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, pb_z, pc_z, pfi0_367, pfh_273, pfh_274, \
                         pfi1_367, pgh_462, pgh_463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_0 * pfh_273[k]
                   + f_4 * pc_z[k] * pgh_462[k];

        t_619[k] = pb_z[k] * pfi0_367[k]
                   - f_11 * pc_z[k] * pfi1_367[k];

        t_620[k] = f_0 * pfh_274[k]
                   + f_4 * pc_z[k] * pgh_463[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pb_x, pb_z, pc_x, pc_z, pfi0_370, pfi0_481, \
                         pfh_276, pfh_362, pfi1_370, pfi1_481, \
                         pgh_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = pb_x[k] * pfi0_481[k]
                   + f_1 * pfh_362[k]
                   - f_11 * pc_x[k] * pfi1_481[k];

        t_622[k] = pb_z[k] * pfi0_370[k]
                   - f_11 * pc_z[k] * pfi1_370[k];

        t_623[k] = f_0 * pfh_276[k]
                   + f_4 * pc_z[k] * pgh_465[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pb_x, pb_z, pc_x, pc_z, pfi0_374, pfi0_484, \
                         pfi0_485, pfh_365, pfh_366, pfi1_374, pfi1_484, \
                         pfi1_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = pb_x[k] * pfi0_484[k]
                   + f_13 * pfh_365[k]
                   - f_11 * pc_x[k] * pfi1_484[k];

        t_625[k] = pb_x[k] * pfi0_485[k]
                   + f_13 * pfh_366[k]
                   - f_11 * pc_x[k] * pfi1_485[k];

        t_626[k] = pb_z[k] * pfi0_374[k]
                   - f_11 * pc_z[k] * pfi1_374[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pb_x, pc_x, pc_z, pfi0_488, pfi0_489, pfh_279, \
                         pfh_369, pfh_370, pfi1_488, pfi1_489, \
                         pgh_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_0 * pfh_279[k]
                   + f_4 * pc_z[k] * pgh_468[k];

        t_628[k] = pb_x[k] * pfi0_488[k]
                   + f_12 * pfh_369[k]
                   - f_11 * pc_x[k] * pfi1_488[k];

        t_629[k] = pb_x[k] * pfi0_489[k]
                   + f_12 * pfh_370[k]
                   - f_11 * pc_x[k] * pfi1_489[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pb_x, pc_x, pfi0_490, pfh_371, pfh_372, \
                         pfh_373, pfh_374, pfi1_490, pgh_477, pgh_478, \
                         pgh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = pb_x[k] * pfi0_490[k]
                   + f_12 * pfh_371[k]
                   - f_11 * pc_x[k] * pfi1_490[k];

        t_631[k] = f_0 * pfh_372[k]
                   + f_4 * pc_x[k] * pgh_477[k];

        t_632[k] = f_0 * pfh_373[k]
                   + f_4 * pc_x[k] * pgh_478[k];

        t_633[k] = f_0 * pfh_374[k]
                   + f_4 * pc_x[k] * pgh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pb_x, pc_x, pfi0_497, pfh_375, pfh_376, \
                         pfh_377, pfi1_497, pgh_480, pgh_481, pgh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_0 * pfh_375[k]
                   + f_4 * pc_x[k] * pgh_480[k];

        t_635[k] = f_0 * pfh_376[k]
                   + f_4 * pc_x[k] * pgh_481[k];

        t_636[k] = f_0 * pfh_377[k]
                   + f_4 * pc_x[k] * pgh_482[k];

        t_637[k] = pb_x[k] * pfi0_497[k]
                   - f_11 * pc_x[k] * pfi1_497[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pb_x, pc_x, pc_z, pfi0_499, pfi0_500, \
                         pfi0_501, pfh_288, pfi1_499, pfi1_500, pfi1_501, \
                         pgh_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_0 * pfh_288[k]
                   + f_4 * pc_z[k] * pgh_477[k];

        t_639[k] = pb_x[k] * pfi0_499[k]
                   - f_11 * pc_x[k] * pfi1_499[k];

        t_640[k] = pb_x[k] * pfi0_500[k]
                   - f_11 * pc_x[k] * pfi1_500[k];

        t_641[k] = pb_x[k] * pfi0_501[k]
                   - f_11 * pc_x[k] * pfi1_501[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pb_x, pc_x, pfi0_502, pfi0_503, pfh_378, \
                         pfi1_502, pfi1_503, pgg0_345, pgg1_345, \
                         pgh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = pb_x[k] * pfi0_502[k]
                   - f_11 * pc_x[k] * pfi1_502[k];

        t_643[k] = pb_x[k] * pfi0_503[k]
                   - f_11 * pc_x[k] * pfi1_503[k];

        t_644[k] = f_0 * pfh_378[k]
                   + f_2 * pgg0_345[k]
                   - f_3 * pgg1_345[k]
                   + f_4 * pc_x[k] * pgh_483[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pb_x, pc_x, pc_z, pfi0_505, pfi0_507, pfh_294, \
                         pfh_379, pfh_381, pfi1_505, pfi1_507, \
                         pgh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = pb_x[k] * pfi0_505[k]
                   + f_21 * pfh_379[k]
                   - f_11 * pc_x[k] * pfi1_505[k];

        t_646[k] = f_12 * pfh_294[k]
                   + f_4 * pc_z[k] * pgh_483[k];

        t_647[k] = pb_x[k] * pfi0_507[k]
                   + f_1 * pfh_381[k]
                   - f_11 * pc_x[k] * pfi1_507[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pb_x, pc_x, pc_z, pfi0_510, pfh_295, pfh_383, \
                         pfh_384, pfi1_510, pgg0_350, pgg1_350, pgh_484, \
                         pgh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_12 * pfh_295[k]
                   + f_4 * pc_z[k] * pgh_484[k];

        t_649[k] = f_0 * pfh_383[k]
                   + f_9 * pgg0_350[k]
                   - f_10 * pgg1_350[k]
                   + f_4 * pc_x[k] * pgh_488[k];

        t_650[k] = pb_x[k] * pfi0_510[k]
                   + f_13 * pfh_384[k]
                   - f_11 * pc_x[k] * pfi1_510[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pb_x, pc_x, pc_z, pfi0_512, pfh_297, pfh_386, \
                         pfh_387, pfi1_512, pgg0_354, pgg1_354, pgh_486, \
                         pgh_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_12 * pfh_297[k]
                   + f_4 * pc_z[k] * pgh_486[k];

        t_652[k] = pb_x[k] * pfi0_512[k]
                   + f_13 * pfh_386[k]
                   - f_11 * pc_x[k] * pfi1_512[k];

        t_653[k] = f_0 * pfh_387[k]
                   + f_7 * pgg0_354[k]
                   - f_8 * pgg1_354[k]
                   + f_4 * pc_x[k] * pgh_492[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pb_x, pc_x, pc_z, pfi0_514, pfi0_516, pfh_300, \
                         pfh_388, pfh_390, pfi1_514, pfi1_516, \
                         pgh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = pb_x[k] * pfi0_514[k]
                   + f_12 * pfh_388[k]
                   - f_11 * pc_x[k] * pfi1_514[k];

        t_655[k] = f_12 * pfh_300[k]
                   + f_4 * pc_z[k] * pgh_489[k];

        t_656[k] = pb_x[k] * pfi0_516[k]
                   + f_12 * pfh_390[k]
                   - f_11 * pc_x[k] * pfi1_516[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pb_x, pc_x, pfi0_517, pfh_391, pfh_392, pfh_393, \
                         pfi1_517, pgg0_359, pgg1_359, pgh_497, \
                         pgh_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = pb_x[k] * pfi0_517[k]
                   + f_12 * pfh_391[k]
                   - f_11 * pc_x[k] * pfi1_517[k];

        t_658[k] = f_0 * pfh_392[k]
                   + f_5 * pgg0_359[k]
                   - f_6 * pgg1_359[k]
                   + f_4 * pc_x[k] * pgh_497[k];

        t_659[k] = f_0 * pfh_393[k]
                   + f_4 * pc_x[k] * pgh_498[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, pc_x, pfh_394, pfh_395, pfh_396, \
                         pfh_397, pfh_398, pgh_499, pgh_500, pgh_501, pgh_502, \
                         pgh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_0 * pfh_394[k]
                   + f_4 * pc_x[k] * pgh_499[k];

        t_661[k] = f_0 * pfh_395[k]
                   + f_4 * pc_x[k] * pgh_500[k];

        t_662[k] = f_0 * pfh_396[k]
                   + f_4 * pc_x[k] * pgh_501[k];

        t_663[k] = f_0 * pfh_397[k]
                   + f_4 * pc_x[k] * pgh_502[k];

        t_664[k] = f_0 * pfh_398[k]
                   + f_4 * pc_x[k] * pgh_503[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pb_x, pc_x, pc_z, pfi0_525, pfi0_527, \
                         pfi0_528, pfh_309, pfi1_525, pfi1_527, pfi1_528, \
                         pgh_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = pb_x[k] * pfi0_525[k]
                   - f_11 * pc_x[k] * pfi1_525[k];

        t_666[k] = f_12 * pfh_309[k]
                   + f_4 * pc_z[k] * pgh_498[k];

        t_667[k] = pb_x[k] * pfi0_527[k]
                   - f_11 * pc_x[k] * pfi1_527[k];

        t_668[k] = pb_x[k] * pfi0_528[k]
                   - f_11 * pc_x[k] * pfi1_528[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, pa_y, pb_x, pc_x, pc_y, sgi0_252, \
                         sgi1_252, pfi0_529, pfi0_530, pfi0_531, pfi1_529, pfi1_530, \
                         pfi1_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = pb_x[k] * pfi0_529[k]
                   - f_11 * pc_x[k] * pfi1_529[k];

        t_670[k] = pb_x[k] * pfi0_530[k]
                   - f_11 * pc_x[k] * pfi1_530[k];

        t_671[k] = pb_x[k] * pfi0_531[k]
                   - f_11 * pc_x[k] * pfi1_531[k];

        t_672[k] = pa_y[k] * sgi0_252[k]
                   - f_11 * pc_y[k] * sgi1_252[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, pa_y, pc_y, pc_z, sgi0_253, sgi0_255, sgh_189, \
                         sgh_190, sgi1_253, sgi1_255, pfh_315, \
                         pgh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = pa_y[k] * sgi0_253[k]
                   + f_0 * sgh_189[k]
                   - f_11 * pc_y[k] * sgi1_253[k];

        t_674[k] = f_13 * pfh_315[k]
                   + f_4 * pc_z[k] * pgh_504[k];

        t_675[k] = pa_y[k] * sgi0_255[k]
                   + f_12 * sgh_190[k]
                   - f_11 * pc_y[k] * sgi1_255[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pa_y, pc_y, pc_z, sgi0_257, sgi0_258, \
                         sgh_192, sgi1_257, sgi1_258, pfh_316, pfh_318, pgh_505, \
                         pgh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_13 * pfh_316[k]
                   + f_4 * pc_z[k] * pgh_505[k];

        t_677[k] = pa_y[k] * sgi0_257[k]
                   - f_11 * pc_y[k] * sgi1_257[k];

        t_678[k] = pa_y[k] * sgi0_258[k]
                   + f_13 * sgh_192[k]
                   - f_11 * pc_y[k] * sgi1_258[k];

        t_679[k] = f_13 * pfh_318[k]
                   + f_4 * pc_z[k] * pgh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pa_y, pc_y, sgi0_260, sgi0_261, sgi0_262, \
                         sgh_194, sgh_195, sgi1_260, sgi1_261, \
                         sgi1_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = pa_y[k] * sgi0_260[k]
                   + f_0 * sgh_194[k]
                   - f_11 * pc_y[k] * sgi1_260[k];

        t_681[k] = pa_y[k] * sgi0_261[k]
                   - f_11 * pc_y[k] * sgi1_261[k];

        t_682[k] = pa_y[k] * sgi0_262[k]
                   + f_1 * sgh_195[k]
                   - f_11 * pc_y[k] * sgi1_262[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pa_y, pc_y, pc_z, sgi0_264, sgi0_265, sgh_197, \
                         sgh_198, sgi1_264, sgi1_265, pfh_321, \
                         pgh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * pfh_321[k]
                   + f_4 * pc_z[k] * pgh_510[k];

        t_684[k] = pa_y[k] * sgi0_264[k]
                   + f_12 * sgh_197[k]
                   - f_11 * pc_y[k] * sgi1_264[k];

        t_685[k] = pa_y[k] * sgi0_265[k]
                   + f_0 * sgh_198[k]
                   - f_11 * pc_y[k] * sgi1_265[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pa_y, pc_x, pc_y, sgi0_266, sgi1_266, \
                         pfh_414, pfh_415, pfh_416, pgh_519, pgh_520, \
                         pgh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = pa_y[k] * sgi0_266[k]
                   - f_11 * pc_y[k] * sgi1_266[k];

        t_687[k] = f_0 * pfh_414[k]
                   + f_4 * pc_x[k] * pgh_519[k];

        t_688[k] = f_0 * pfh_415[k]
                   + f_4 * pc_x[k] * pgh_520[k];

        t_689[k] = f_0 * pfh_416[k]
                   + f_4 * pc_x[k] * pgh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pb_x, pc_x, pfi0_553, pfh_417, pfh_418, \
                         pfh_419, pfi1_553, pgh_522, pgh_523, pgh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_0 * pfh_417[k]
                   + f_4 * pc_x[k] * pgh_522[k];

        t_691[k] = f_0 * pfh_418[k]
                   + f_4 * pc_x[k] * pgh_523[k];

        t_692[k] = f_0 * pfh_419[k]
                   + f_4 * pc_x[k] * pgh_524[k];

        t_693[k] = pb_x[k] * pfi0_553[k]
                   - f_11 * pc_x[k] * pfi1_553[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pb_x, pc_x, pc_z, pfi0_555, pfi0_556, \
                         pfi0_557, pfh_330, pfi1_555, pfi1_556, pfi1_557, \
                         pgh_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * pfh_330[k]
                   + f_4 * pc_z[k] * pgh_519[k];

        t_695[k] = pb_x[k] * pfi0_555[k]
                   - f_11 * pc_x[k] * pfi1_555[k];

        t_696[k] = pb_x[k] * pfi0_556[k]
                   - f_11 * pc_x[k] * pfi1_556[k];

        t_697[k] = pb_x[k] * pfi0_557[k]
                   - f_11 * pc_x[k] * pfi1_557[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pa_y, pc_x, pc_y, sgi0_279, sgh_209, sgi1_279, \
                         pgg0_375, pgg1_375, pgh_524, pgh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_0 * sgh_209[k]
                   + f_4 * pc_y[k] * pgh_524[k];

        t_699[k] = pa_y[k] * sgi0_279[k]
                   - f_11 * pc_y[k] * sgi1_279[k];

        t_700[k] = f_2 * pgg0_375[k]
                   - f_3 * pgg1_375[k]
                   + f_4 * pc_x[k] * pgh_525[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pc_x, pc_z, pgg0_376, pgg0_378, pgg1_376, \
                         pgg1_378, pgh_525, pgh_526, pgh_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_17 * pgg0_376[k]
                   - f_18 * pgg1_376[k]
                   + f_4 * pc_x[k] * pgh_526[k];

        t_702[k] = f_4 * pc_z[k] * pgh_525[k];

        t_703[k] = f_9 * pgg0_378[k]
                   - f_10 * pgg1_378[k]
                   + f_4 * pc_x[k] * pgh_528[k];

        t_704[k] = f_4 * pc_z[k] * pgh_526[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pc_x, pc_z, pgg0_380, pgg0_381, pgg0_383, \
                         pgg1_380, pgg1_381, pgg1_383, pgh_528, pgh_530, pgh_531, \
                         pgh_533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_9 * pgg0_380[k]
                   - f_10 * pgg1_380[k]
                   + f_4 * pc_x[k] * pgh_530[k];

        t_706[k] = f_7 * pgg0_381[k]
                   - f_8 * pgg1_381[k]
                   + f_4 * pc_x[k] * pgh_531[k];

        t_707[k] = f_4 * pc_z[k] * pgh_528[k];

        t_708[k] = f_7 * pgg0_383[k]
                   - f_8 * pgg1_383[k]
                   + f_4 * pc_x[k] * pgh_533[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, pc_x, pc_z, pgg0_384, pgg0_385, pgg0_387, \
                         pgg1_384, pgg1_385, pgg1_387, pgh_531, pgh_534, pgh_535, \
                         pgh_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_7 * pgg0_384[k]
                   - f_8 * pgg1_384[k]
                   + f_4 * pc_x[k] * pgh_534[k];

        t_710[k] = f_5 * pgg0_385[k]
                   - f_6 * pgg1_385[k]
                   + f_4 * pc_x[k] * pgh_535[k];

        t_711[k] = f_4 * pc_z[k] * pgh_531[k];

        t_712[k] = f_5 * pgg0_387[k]
                   - f_6 * pgg1_387[k]
                   + f_4 * pc_x[k] * pgh_537[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgi0, const size_t sgh,
                                                          const size_t sgi1, const size_t pdi0,
                                                          const size_t pdi1, const size_t pfi0,
                                                          const size_t pfh, const size_t pfi1,
                                                          const size_t pgg0, const size_t pgg1,
                                                          const size_t pgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);

    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_392 = buffer.data(sgi0 + 392);
    const auto *sgi0_397 = buffer.data(sgi0 + 397);
    const auto *sgi0_401 = buffer.data(sgi0 + 401);
    const auto *sgi0_406 = buffer.data(sgi0 + 406);
    const auto *sgi0_413 = buffer.data(sgi0 + 413);

    const auto *sgh_225 = buffer.data(sgh + 225);
    const auto *sgh_230 = buffer.data(sgh + 230);
    const auto *sgh_251 = buffer.data(sgh + 251);
    const auto *sgh_272 = buffer.data(sgh + 272);
    const auto *sgh_288 = buffer.data(sgh + 288);
    const auto *sgh_293 = buffer.data(sgh + 293);
    const auto *sgh_309 = buffer.data(sgh + 309);

    const auto *sgi1_392 = buffer.data(sgi1 + 392);
    const auto *sgi1_397 = buffer.data(sgi1 + 397);
    const auto *sgi1_401 = buffer.data(sgi1 + 401);
    const auto *sgi1_406 = buffer.data(sgi1 + 406);
    const auto *sgi1_413 = buffer.data(sgi1 + 413);

    const auto *pdi0_273 = buffer.data(pdi0 + 273);

    const auto *pdi1_273 = buffer.data(pdi1 + 273);

    const auto *pfi0_448 = buffer.data(pfi0 + 448);
    const auto *pfi0_449 = buffer.data(pfi0 + 449);
    const auto *pfi0_451 = buffer.data(pfi0 + 451);
    const auto *pfi0_454 = buffer.data(pfi0 + 454);
    const auto *pfi0_458 = buffer.data(pfi0 + 458);
    const auto *pfi0_469 = buffer.data(pfi0 + 469);
    const auto *pfi0_471 = buffer.data(pfi0 + 471);
    const auto *pfi0_472 = buffer.data(pfi0 + 472);
    const auto *pfi0_473 = buffer.data(pfi0 + 473);
    const auto *pfi0_497 = buffer.data(pfi0 + 497);

    const auto *pfh_336 = buffer.data(pfh + 336);
    const auto *pfh_337 = buffer.data(pfh + 337);
    const auto *pfh_339 = buffer.data(pfh + 339);
    const auto *pfh_342 = buffer.data(pfh + 342);
    const auto *pfh_351 = buffer.data(pfh + 351);
    const auto *pfh_352 = buffer.data(pfh + 352);
    const auto *pfh_353 = buffer.data(pfh + 353);
    const auto *pfh_354 = buffer.data(pfh + 354);
    const auto *pfh_356 = buffer.data(pfh + 356);
    const auto *pfh_357 = buffer.data(pfh + 357);
    const auto *pfh_358 = buffer.data(pfh + 358);
    const auto *pfh_360 = buffer.data(pfh + 360);
    const auto *pfh_363 = buffer.data(pfh + 363);
    const auto *pfh_372 = buffer.data(pfh + 372);
    const auto *pfh_373 = buffer.data(pfh + 373);
    const auto *pfh_374 = buffer.data(pfh + 374);
    const auto *pfh_375 = buffer.data(pfh + 375);
    const auto *pfh_377 = buffer.data(pfh + 377);
    const auto *pfh_378 = buffer.data(pfh + 378);
    const auto *pfh_379 = buffer.data(pfh + 379);
    const auto *pfh_381 = buffer.data(pfh + 381);
    const auto *pfh_384 = buffer.data(pfh + 384);
    const auto *pfh_393 = buffer.data(pfh + 393);
    const auto *pfh_394 = buffer.data(pfh + 394);
    const auto *pfh_395 = buffer.data(pfh + 395);
    const auto *pfh_396 = buffer.data(pfh + 396);
    const auto *pfh_398 = buffer.data(pfh + 398);
    const auto *pfh_399 = buffer.data(pfh + 399);
    const auto *pfh_400 = buffer.data(pfh + 400);
    const auto *pfh_402 = buffer.data(pfh + 402);
    const auto *pfh_405 = buffer.data(pfh + 405);
    const auto *pfh_414 = buffer.data(pfh + 414);
    const auto *pfh_419 = buffer.data(pfh + 419);

    const auto *pfi1_448 = buffer.data(pfi1 + 448);
    const auto *pfi1_449 = buffer.data(pfi1 + 449);
    const auto *pfi1_451 = buffer.data(pfi1 + 451);
    const auto *pfi1_454 = buffer.data(pfi1 + 454);
    const auto *pfi1_458 = buffer.data(pfi1 + 458);
    const auto *pfi1_469 = buffer.data(pfi1 + 469);
    const auto *pfi1_471 = buffer.data(pfi1 + 471);
    const auto *pfi1_472 = buffer.data(pfi1 + 472);
    const auto *pfi1_473 = buffer.data(pfi1 + 473);
    const auto *pfi1_497 = buffer.data(pfi1 + 497);

    const auto *pgg0_385 = buffer.data(pgg0 + 385);
    const auto *pgg0_386 = buffer.data(pgg0 + 386);
    const auto *pgg0_387 = buffer.data(pgg0 + 387);
    const auto *pgg0_388 = buffer.data(pgg0 + 388);
    const auto *pgg0_389 = buffer.data(pgg0 + 389);
    const auto *pgg0_395 = buffer.data(pgg0 + 395);
    const auto *pgg0_398 = buffer.data(pgg0 + 398);
    const auto *pgg0_399 = buffer.data(pgg0 + 399);
    const auto *pgg0_402 = buffer.data(pgg0 + 402);
    const auto *pgg0_403 = buffer.data(pgg0 + 403);
    const auto *pgg0_404 = buffer.data(pgg0 + 404);
    const auto *pgg0_405 = buffer.data(pgg0 + 405);
    const auto *pgg0_406 = buffer.data(pgg0 + 406);
    const auto *pgg0_408 = buffer.data(pgg0 + 408);
    const auto *pgg0_410 = buffer.data(pgg0 + 410);
    const auto *pgg0_411 = buffer.data(pgg0 + 411);
    const auto *pgg0_413 = buffer.data(pgg0 + 413);
    const auto *pgg0_414 = buffer.data(pgg0 + 414);
    const auto *pgg0_415 = buffer.data(pgg0 + 415);
    const auto *pgg0_416 = buffer.data(pgg0 + 416);
    const auto *pgg0_417 = buffer.data(pgg0 + 417);
    const auto *pgg0_418 = buffer.data(pgg0 + 418);
    const auto *pgg0_419 = buffer.data(pgg0 + 419);
    const auto *pgg0_420 = buffer.data(pgg0 + 420);
    const auto *pgg0_421 = buffer.data(pgg0 + 421);
    const auto *pgg0_423 = buffer.data(pgg0 + 423);
    const auto *pgg0_425 = buffer.data(pgg0 + 425);
    const auto *pgg0_426 = buffer.data(pgg0 + 426);
    const auto *pgg0_428 = buffer.data(pgg0 + 428);
    const auto *pgg0_429 = buffer.data(pgg0 + 429);
    const auto *pgg0_430 = buffer.data(pgg0 + 430);
    const auto *pgg0_431 = buffer.data(pgg0 + 431);
    const auto *pgg0_432 = buffer.data(pgg0 + 432);
    const auto *pgg0_433 = buffer.data(pgg0 + 433);
    const auto *pgg0_434 = buffer.data(pgg0 + 434);
    const auto *pgg0_436 = buffer.data(pgg0 + 436);
    const auto *pgg0_438 = buffer.data(pgg0 + 438);
    const auto *pgg0_441 = buffer.data(pgg0 + 441);
    const auto *pgg0_443 = buffer.data(pgg0 + 443);
    const auto *pgg0_445 = buffer.data(pgg0 + 445);
    const auto *pgg0_447 = buffer.data(pgg0 + 447);
    const auto *pgg0_448 = buffer.data(pgg0 + 448);

    const auto *pgg1_385 = buffer.data(pgg1 + 385);
    const auto *pgg1_386 = buffer.data(pgg1 + 386);
    const auto *pgg1_387 = buffer.data(pgg1 + 387);
    const auto *pgg1_388 = buffer.data(pgg1 + 388);
    const auto *pgg1_389 = buffer.data(pgg1 + 389);
    const auto *pgg1_395 = buffer.data(pgg1 + 395);
    const auto *pgg1_398 = buffer.data(pgg1 + 398);
    const auto *pgg1_399 = buffer.data(pgg1 + 399);
    const auto *pgg1_402 = buffer.data(pgg1 + 402);
    const auto *pgg1_403 = buffer.data(pgg1 + 403);
    const auto *pgg1_404 = buffer.data(pgg1 + 404);
    const auto *pgg1_405 = buffer.data(pgg1 + 405);
    const auto *pgg1_406 = buffer.data(pgg1 + 406);
    const auto *pgg1_408 = buffer.data(pgg1 + 408);
    const auto *pgg1_410 = buffer.data(pgg1 + 410);
    const auto *pgg1_411 = buffer.data(pgg1 + 411);
    const auto *pgg1_413 = buffer.data(pgg1 + 413);
    const auto *pgg1_414 = buffer.data(pgg1 + 414);
    const auto *pgg1_415 = buffer.data(pgg1 + 415);
    const auto *pgg1_416 = buffer.data(pgg1 + 416);
    const auto *pgg1_417 = buffer.data(pgg1 + 417);
    const auto *pgg1_418 = buffer.data(pgg1 + 418);
    const auto *pgg1_419 = buffer.data(pgg1 + 419);
    const auto *pgg1_420 = buffer.data(pgg1 + 420);
    const auto *pgg1_421 = buffer.data(pgg1 + 421);
    const auto *pgg1_423 = buffer.data(pgg1 + 423);
    const auto *pgg1_425 = buffer.data(pgg1 + 425);
    const auto *pgg1_426 = buffer.data(pgg1 + 426);
    const auto *pgg1_428 = buffer.data(pgg1 + 428);
    const auto *pgg1_429 = buffer.data(pgg1 + 429);
    const auto *pgg1_430 = buffer.data(pgg1 + 430);
    const auto *pgg1_431 = buffer.data(pgg1 + 431);
    const auto *pgg1_432 = buffer.data(pgg1 + 432);
    const auto *pgg1_433 = buffer.data(pgg1 + 433);
    const auto *pgg1_434 = buffer.data(pgg1 + 434);
    const auto *pgg1_436 = buffer.data(pgg1 + 436);
    const auto *pgg1_438 = buffer.data(pgg1 + 438);
    const auto *pgg1_441 = buffer.data(pgg1 + 441);
    const auto *pgg1_443 = buffer.data(pgg1 + 443);
    const auto *pgg1_445 = buffer.data(pgg1 + 445);
    const auto *pgg1_447 = buffer.data(pgg1 + 447);
    const auto *pgg1_448 = buffer.data(pgg1 + 448);

    const auto *pgh_538 = buffer.data(pgh + 538);
    const auto *pgh_539 = buffer.data(pgh + 539);
    const auto *pgh_540 = buffer.data(pgh + 540);
    const auto *pgh_541 = buffer.data(pgh + 541);
    const auto *pgh_542 = buffer.data(pgh + 542);
    const auto *pgh_543 = buffer.data(pgh + 543);
    const auto *pgh_544 = buffer.data(pgh + 544);
    const auto *pgh_545 = buffer.data(pgh + 545);
    const auto *pgh_546 = buffer.data(pgh + 546);
    const auto *pgh_547 = buffer.data(pgh + 547);
    const auto *pgh_549 = buffer.data(pgh + 549);
    const auto *pgh_551 = buffer.data(pgh + 551);
    const auto *pgh_552 = buffer.data(pgh + 552);
    const auto *pgh_554 = buffer.data(pgh + 554);
    const auto *pgh_555 = buffer.data(pgh + 555);
    const auto *pgh_558 = buffer.data(pgh + 558);
    const auto *pgh_559 = buffer.data(pgh + 559);
    const auto *pgh_560 = buffer.data(pgh + 560);
    const auto *pgh_561 = buffer.data(pgh + 561);
    const auto *pgh_562 = buffer.data(pgh + 562);
    const auto *pgh_563 = buffer.data(pgh + 563);
    const auto *pgh_564 = buffer.data(pgh + 564);
    const auto *pgh_565 = buffer.data(pgh + 565);
    const auto *pgh_566 = buffer.data(pgh + 566);
    const auto *pgh_567 = buffer.data(pgh + 567);
    const auto *pgh_568 = buffer.data(pgh + 568);
    const auto *pgh_570 = buffer.data(pgh + 570);
    const auto *pgh_572 = buffer.data(pgh + 572);
    const auto *pgh_573 = buffer.data(pgh + 573);
    const auto *pgh_575 = buffer.data(pgh + 575);
    const auto *pgh_576 = buffer.data(pgh + 576);
    const auto *pgh_577 = buffer.data(pgh + 577);
    const auto *pgh_579 = buffer.data(pgh + 579);
    const auto *pgh_580 = buffer.data(pgh + 580);
    const auto *pgh_581 = buffer.data(pgh + 581);
    const auto *pgh_582 = buffer.data(pgh + 582);
    const auto *pgh_583 = buffer.data(pgh + 583);
    const auto *pgh_584 = buffer.data(pgh + 584);
    const auto *pgh_585 = buffer.data(pgh + 585);
    const auto *pgh_586 = buffer.data(pgh + 586);
    const auto *pgh_587 = buffer.data(pgh + 587);
    const auto *pgh_588 = buffer.data(pgh + 588);
    const auto *pgh_589 = buffer.data(pgh + 589);
    const auto *pgh_591 = buffer.data(pgh + 591);
    const auto *pgh_593 = buffer.data(pgh + 593);
    const auto *pgh_594 = buffer.data(pgh + 594);
    const auto *pgh_596 = buffer.data(pgh + 596);
    const auto *pgh_597 = buffer.data(pgh + 597);
    const auto *pgh_598 = buffer.data(pgh + 598);
    const auto *pgh_600 = buffer.data(pgh + 600);
    const auto *pgh_601 = buffer.data(pgh + 601);
    const auto *pgh_602 = buffer.data(pgh + 602);
    const auto *pgh_603 = buffer.data(pgh + 603);
    const auto *pgh_604 = buffer.data(pgh + 604);
    const auto *pgh_605 = buffer.data(pgh + 605);
    const auto *pgh_606 = buffer.data(pgh + 606);
    const auto *pgh_607 = buffer.data(pgh + 607);
    const auto *pgh_608 = buffer.data(pgh + 608);
    const auto *pgh_609 = buffer.data(pgh + 609);
    const auto *pgh_610 = buffer.data(pgh + 610);
    const auto *pgh_612 = buffer.data(pgh + 612);
    const auto *pgh_615 = buffer.data(pgh + 615);
    const auto *pgh_617 = buffer.data(pgh + 617);
    const auto *pgh_619 = buffer.data(pgh + 619);
    const auto *pgh_621 = buffer.data(pgh + 621);
    const auto *pgh_622 = buffer.data(pgh + 622);
    const auto *pgh_624 = buffer.data(pgh + 624);
    const auto *pgh_625 = buffer.data(pgh + 625);
    const auto *pgh_626 = buffer.data(pgh + 626);
    const auto *pgh_627 = buffer.data(pgh + 627);
    const auto *pgh_628 = buffer.data(pgh + 628);
    const auto *pgh_629 = buffer.data(pgh + 629);

#pragma omp simd aligned(t_713, t_714, t_715, t_716, t_717, pc_x, pgg0_388, pgg0_389, \
                         pgg1_388, pgg1_389, pgh_538, pgh_539, pgh_540, pgh_541, \
                         pgh_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_713[k] = f_5 * pgg0_388[k]
                   - f_6 * pgg1_388[k]
                   + f_4 * pc_x[k] * pgh_538[k];

        t_714[k] = f_5 * pgg0_389[k]
                   - f_6 * pgg1_389[k]
                   + f_4 * pc_x[k] * pgh_539[k];

        t_715[k] = f_4 * pc_x[k] * pgh_540[k];

        t_716[k] = f_4 * pc_x[k] * pgh_541[k];

        t_717[k] = f_4 * pc_x[k] * pgh_542[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, t_722, pc_x, pc_y, pc_z, sgh_225, \
                         pfh_351, pgg0_385, pgg1_385, pgh_540, pgh_543, pgh_544, \
                         pgh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_4 * pc_x[k] * pgh_543[k];

        t_719[k] = f_4 * pc_x[k] * pgh_544[k];

        t_720[k] = f_4 * pc_x[k] * pgh_545[k];

        t_721[k] = f_0 * sgh_225[k]
                   + f_1 * pfh_351[k]
                   + f_2 * pgg0_385[k]
                   - f_3 * pgg1_385[k]
                   + f_4 * pc_y[k] * pgh_540[k];

        t_722[k] = f_4 * pc_z[k] * pgh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_z, pgg0_385, pgg0_386, pgg0_387, pgg1_385, \
                         pgg1_386, pgg1_387, pgh_541, pgh_542, \
                         pgh_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_5 * pgg0_385[k]
                   - f_6 * pgg1_385[k]
                   + f_4 * pc_z[k] * pgh_541[k];

        t_724[k] = f_7 * pgg0_386[k]
                   - f_8 * pgg1_386[k]
                   + f_4 * pc_z[k] * pgh_542[k];

        t_725[k] = f_9 * pgg0_387[k]
                   - f_10 * pgg1_387[k]
                   + f_4 * pc_z[k] * pgh_543[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pb_z, pc_y, pc_z, sgh_230, pfi0_448, \
                         pfi0_449, pfh_356, pfi1_448, pfi1_449, pgg0_389, pgg1_389, \
                         pgh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_0 * sgh_230[k]
                   + f_1 * pfh_356[k]
                   + f_4 * pc_y[k] * pgh_545[k];

        t_727[k] = f_2 * pgg0_389[k]
                   - f_3 * pgg1_389[k]
                   + f_4 * pc_z[k] * pgh_545[k];

        t_728[k] = pb_z[k] * pfi0_448[k]
                   - f_11 * pc_z[k] * pfi1_448[k];

        t_729[k] = pb_z[k] * pfi0_449[k]
                   - f_11 * pc_z[k] * pfi1_449[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pb_z, pc_x, pc_z, pfi0_451, pfh_336, \
                         pfh_337, pfi1_451, pgg0_395, pgg1_395, pgh_546, pgh_547, \
                         pgh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_0 * pfh_336[k]
                   + f_4 * pc_z[k] * pgh_546[k];

        t_731[k] = pb_z[k] * pfi0_451[k]
                   - f_11 * pc_z[k] * pfi1_451[k];

        t_732[k] = f_0 * pfh_337[k]
                   + f_4 * pc_z[k] * pgh_547[k];

        t_733[k] = f_9 * pgg0_395[k]
                   - f_10 * pgg1_395[k]
                   + f_4 * pc_x[k] * pgh_551[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pb_z, pc_x, pc_z, pfi0_454, pfh_339, pfi1_454, \
                         pgg0_398, pgg1_398, pgh_549, pgh_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pb_z[k] * pfi0_454[k]
                   - f_11 * pc_z[k] * pfi1_454[k];

        t_735[k] = f_0 * pfh_339[k]
                   + f_4 * pc_z[k] * pgh_549[k];

        t_736[k] = f_7 * pgg0_398[k]
                   - f_8 * pgg1_398[k]
                   + f_4 * pc_x[k] * pgh_554[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pb_z, pc_x, pc_z, pfi0_458, pfh_342, pfi1_458, \
                         pgg0_399, pgg1_399, pgh_552, pgh_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_7 * pgg0_399[k]
                   - f_8 * pgg1_399[k]
                   + f_4 * pc_x[k] * pgh_555[k];

        t_738[k] = pb_z[k] * pfi0_458[k]
                   - f_11 * pc_z[k] * pfi1_458[k];

        t_739[k] = f_0 * pfh_342[k]
                   + f_4 * pc_z[k] * pgh_552[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, pgg0_402, pgg0_403, pgg0_404, \
                         pgg1_402, pgg1_403, pgg1_404, pgh_558, pgh_559, pgh_560, \
                         pgh_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_5 * pgg0_402[k]
                   - f_6 * pgg1_402[k]
                   + f_4 * pc_x[k] * pgh_558[k];

        t_741[k] = f_5 * pgg0_403[k]
                   - f_6 * pgg1_403[k]
                   + f_4 * pc_x[k] * pgh_559[k];

        t_742[k] = f_5 * pgg0_404[k]
                   - f_6 * pgg1_404[k]
                   + f_4 * pc_x[k] * pgh_560[k];

        t_743[k] = f_4 * pc_x[k] * pgh_561[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, t_749, pb_z, pc_x, pc_z, pfi0_469, \
                         pfi1_469, pgh_562, pgh_563, pgh_564, pgh_565, \
                         pgh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_4 * pc_x[k] * pgh_562[k];

        t_745[k] = f_4 * pc_x[k] * pgh_563[k];

        t_746[k] = f_4 * pc_x[k] * pgh_564[k];

        t_747[k] = f_4 * pc_x[k] * pgh_565[k];

        t_748[k] = f_4 * pc_x[k] * pgh_566[k];

        t_749[k] = pb_z[k] * pfi0_469[k]
                   - f_11 * pc_z[k] * pfi1_469[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pb_z, pc_z, pfi0_471, pfi0_472, pfh_351, \
                         pfh_352, pfh_353, pfi1_471, pfi1_472, \
                         pgh_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_0 * pfh_351[k]
                   + f_4 * pc_z[k] * pgh_561[k];

        t_751[k] = pb_z[k] * pfi0_471[k]
                   + f_12 * pfh_352[k]
                   - f_11 * pc_z[k] * pfi1_471[k];

        t_752[k] = pb_z[k] * pfi0_472[k]
                   + f_13 * pfh_353[k]
                   - f_11 * pc_z[k] * pfi1_472[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, pb_z, pc_y, pc_z, sgh_251, pfi0_473, pfh_354, \
                         pfh_356, pfh_377, pfi1_473, pgg0_404, pgg1_404, \
                         pgh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = pb_z[k] * pfi0_473[k]
                   + f_1 * pfh_354[k]
                   - f_11 * pc_z[k] * pfi1_473[k];

        t_754[k] = f_0 * sgh_251[k]
                   + f_13 * pfh_377[k]
                   + f_4 * pc_y[k] * pgh_566[k];

        t_755[k] = f_0 * pfh_356[k]
                   + f_2 * pgg0_404[k]
                   - f_3 * pgg1_404[k]
                   + f_4 * pc_z[k] * pgh_566[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, pc_z, pfh_357, pgg0_405, pgg0_406, \
                         pgg0_408, pgg1_405, pgg1_406, pgg1_408, pgh_567, pgh_568, \
                         pgh_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_2 * pgg0_405[k]
                   - f_3 * pgg1_405[k]
                   + f_4 * pc_x[k] * pgh_567[k];

        t_757[k] = f_17 * pgg0_406[k]
                   - f_18 * pgg1_406[k]
                   + f_4 * pc_x[k] * pgh_568[k];

        t_758[k] = f_12 * pfh_357[k]
                   + f_4 * pc_z[k] * pgh_567[k];

        t_759[k] = f_9 * pgg0_408[k]
                   - f_10 * pgg1_408[k]
                   + f_4 * pc_x[k] * pgh_570[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, pc_x, pc_z, pfh_358, pfh_360, pgg0_410, \
                         pgg0_411, pgg1_410, pgg1_411, pgh_568, pgh_570, pgh_572, \
                         pgh_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_12 * pfh_358[k]
                   + f_4 * pc_z[k] * pgh_568[k];

        t_761[k] = f_9 * pgg0_410[k]
                   - f_10 * pgg1_410[k]
                   + f_4 * pc_x[k] * pgh_572[k];

        t_762[k] = f_7 * pgg0_411[k]
                   - f_8 * pgg1_411[k]
                   + f_4 * pc_x[k] * pgh_573[k];

        t_763[k] = f_12 * pfh_360[k]
                   + f_4 * pc_z[k] * pgh_570[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pc_x, pgg0_413, pgg0_414, pgg0_415, pgg1_413, \
                         pgg1_414, pgg1_415, pgh_575, pgh_576, \
                         pgh_577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_7 * pgg0_413[k]
                   - f_8 * pgg1_413[k]
                   + f_4 * pc_x[k] * pgh_575[k];

        t_765[k] = f_7 * pgg0_414[k]
                   - f_8 * pgg1_414[k]
                   + f_4 * pc_x[k] * pgh_576[k];

        t_766[k] = f_5 * pgg0_415[k]
                   - f_6 * pgg1_415[k]
                   + f_4 * pc_x[k] * pgh_577[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pc_x, pc_z, pfh_363, pgg0_417, pgg0_418, \
                         pgg1_417, pgg1_418, pgh_573, pgh_579, \
                         pgh_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_12 * pfh_363[k]
                   + f_4 * pc_z[k] * pgh_573[k];

        t_768[k] = f_5 * pgg0_417[k]
                   - f_6 * pgg1_417[k]
                   + f_4 * pc_x[k] * pgh_579[k];

        t_769[k] = f_5 * pgg0_418[k]
                   - f_6 * pgg1_418[k]
                   + f_4 * pc_x[k] * pgh_580[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, t_775, pc_x, pgg0_419, pgg1_419, \
                         pgh_581, pgh_582, pgh_583, pgh_584, pgh_585, \
                         pgh_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_5 * pgg0_419[k]
                   - f_6 * pgg1_419[k]
                   + f_4 * pc_x[k] * pgh_581[k];

        t_771[k] = f_4 * pc_x[k] * pgh_582[k];

        t_772[k] = f_4 * pc_x[k] * pgh_583[k];

        t_773[k] = f_4 * pc_x[k] * pgh_584[k];

        t_774[k] = f_4 * pc_x[k] * pgh_585[k];

        t_775[k] = f_4 * pc_x[k] * pgh_586[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pb_z, pc_x, pc_z, pdi0_273, pdi1_273, pfi0_497, \
                         pfh_372, pfi1_497, pgh_582, pgh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_4 * pc_x[k] * pgh_587[k];

        t_777[k] = f_14 * pdi0_273[k]
                   - f_15 * pdi1_273[k]
                   + pb_z[k] * pfi0_497[k]
                   - f_11 * pc_z[k] * pfi1_497[k];

        t_778[k] = f_12 * pfh_372[k]
                   + f_4 * pc_z[k] * pgh_582[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pc_z, pfh_373, pfh_374, pfh_375, pgg0_415, \
                         pgg0_416, pgg0_417, pgg1_415, pgg1_416, pgg1_417, pgh_583, pgh_584, \
                         pgh_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_12 * pfh_373[k]
                   + f_5 * pgg0_415[k]
                   - f_6 * pgg1_415[k]
                   + f_4 * pc_z[k] * pgh_583[k];

        t_780[k] = f_12 * pfh_374[k]
                   + f_7 * pgg0_416[k]
                   - f_8 * pgg1_416[k]
                   + f_4 * pc_z[k] * pgh_584[k];

        t_781[k] = f_12 * pfh_375[k]
                   + f_9 * pgg0_417[k]
                   - f_10 * pgg1_417[k]
                   + f_4 * pc_z[k] * pgh_585[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pc_x, pc_y, pc_z, sgh_272, pfh_377, pfh_398, \
                         pgg0_419, pgg0_420, pgg1_419, pgg1_420, pgh_587, \
                         pgh_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_0 * sgh_272[k]
                   + f_12 * pfh_398[k]
                   + f_4 * pc_y[k] * pgh_587[k];

        t_783[k] = f_12 * pfh_377[k]
                   + f_2 * pgg0_419[k]
                   - f_3 * pgg1_419[k]
                   + f_4 * pc_z[k] * pgh_587[k];

        t_784[k] = f_2 * pgg0_420[k]
                   - f_3 * pgg1_420[k]
                   + f_4 * pc_x[k] * pgh_588[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, pc_x, pc_z, pfh_378, pfh_379, pgg0_421, \
                         pgg0_423, pgg1_421, pgg1_423, pgh_588, pgh_589, \
                         pgh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_17 * pgg0_421[k]
                   - f_18 * pgg1_421[k]
                   + f_4 * pc_x[k] * pgh_589[k];

        t_786[k] = f_13 * pfh_378[k]
                   + f_4 * pc_z[k] * pgh_588[k];

        t_787[k] = f_9 * pgg0_423[k]
                   - f_10 * pgg1_423[k]
                   + f_4 * pc_x[k] * pgh_591[k];

        t_788[k] = f_13 * pfh_379[k]
                   + f_4 * pc_z[k] * pgh_589[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, pc_x, pc_z, pfh_381, pgg0_425, pgg0_426, \
                         pgg1_425, pgg1_426, pgh_591, pgh_593, \
                         pgh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_9 * pgg0_425[k]
                   - f_10 * pgg1_425[k]
                   + f_4 * pc_x[k] * pgh_593[k];

        t_790[k] = f_7 * pgg0_426[k]
                   - f_8 * pgg1_426[k]
                   + f_4 * pc_x[k] * pgh_594[k];

        t_791[k] = f_13 * pfh_381[k]
                   + f_4 * pc_z[k] * pgh_591[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, pc_x, pgg0_428, pgg0_429, pgg0_430, pgg1_428, \
                         pgg1_429, pgg1_430, pgh_596, pgh_597, \
                         pgh_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_7 * pgg0_428[k]
                   - f_8 * pgg1_428[k]
                   + f_4 * pc_x[k] * pgh_596[k];

        t_793[k] = f_7 * pgg0_429[k]
                   - f_8 * pgg1_429[k]
                   + f_4 * pc_x[k] * pgh_597[k];

        t_794[k] = f_5 * pgg0_430[k]
                   - f_6 * pgg1_430[k]
                   + f_4 * pc_x[k] * pgh_598[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, pc_x, pc_z, pfh_384, pgg0_432, pgg0_433, \
                         pgg1_432, pgg1_433, pgh_594, pgh_600, \
                         pgh_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_13 * pfh_384[k]
                   + f_4 * pc_z[k] * pgh_594[k];

        t_796[k] = f_5 * pgg0_432[k]
                   - f_6 * pgg1_432[k]
                   + f_4 * pc_x[k] * pgh_600[k];

        t_797[k] = f_5 * pgg0_433[k]
                   - f_6 * pgg1_433[k]
                   + f_4 * pc_x[k] * pgh_601[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, t_801, t_802, t_803, pc_x, pgg0_434, pgg1_434, \
                         pgh_602, pgh_603, pgh_604, pgh_605, pgh_606, \
                         pgh_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = f_5 * pgg0_434[k]
                   - f_6 * pgg1_434[k]
                   + f_4 * pc_x[k] * pgh_602[k];

        t_799[k] = f_4 * pc_x[k] * pgh_603[k];

        t_800[k] = f_4 * pc_x[k] * pgh_604[k];

        t_801[k] = f_4 * pc_x[k] * pgh_605[k];

        t_802[k] = f_4 * pc_x[k] * pgh_606[k];

        t_803[k] = f_4 * pc_x[k] * pgh_607[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, t_807, pc_x, pc_y, pc_z, sgh_288, pfh_393, \
                         pfh_394, pfh_414, pgg0_430, pgg1_430, pgh_603, pgh_604, \
                         pgh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_4 * pc_x[k] * pgh_608[k];

        t_805[k] = f_0 * sgh_288[k]
                   + f_0 * pfh_414[k]
                   + f_2 * pgg0_430[k]
                   - f_3 * pgg1_430[k]
                   + f_4 * pc_y[k] * pgh_603[k];

        t_806[k] = f_13 * pfh_393[k]
                   + f_4 * pc_z[k] * pgh_603[k];

        t_807[k] = f_13 * pfh_394[k]
                   + f_5 * pgg0_430[k]
                   - f_6 * pgg1_430[k]
                   + f_4 * pc_z[k] * pgh_604[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, pc_y, pc_z, sgh_293, pfh_395, pfh_396, pfh_419, \
                         pgg0_431, pgg0_432, pgg1_431, pgg1_432, pgh_605, pgh_606, \
                         pgh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_13 * pfh_395[k]
                   + f_7 * pgg0_431[k]
                   - f_8 * pgg1_431[k]
                   + f_4 * pc_z[k] * pgh_605[k];

        t_809[k] = f_13 * pfh_396[k]
                   + f_9 * pgg0_432[k]
                   - f_10 * pgg1_432[k]
                   + f_4 * pc_z[k] * pgh_606[k];

        t_810[k] = f_0 * sgh_293[k]
                   + f_0 * pfh_419[k]
                   + f_4 * pc_y[k] * pgh_608[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, pa_y, pc_x, pc_y, pc_z, sgi0_392, sgi1_392, \
                         pfh_398, pgg0_434, pgg0_436, pgg1_434, pgg1_436, pgh_608, \
                         pgh_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_13 * pfh_398[k]
                   + f_2 * pgg0_434[k]
                   - f_3 * pgg1_434[k]
                   + f_4 * pc_z[k] * pgh_608[k];

        t_812[k] = pa_y[k] * sgi0_392[k]
                   - f_11 * pc_y[k] * sgi1_392[k];

        t_813[k] = f_17 * pgg0_436[k]
                   - f_18 * pgg1_436[k]
                   + f_4 * pc_x[k] * pgh_610[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, pc_x, pc_z, pfh_399, pfh_400, pgg0_438, \
                         pgg1_438, pgh_609, pgh_610, pgh_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_1 * pfh_399[k]
                   + f_4 * pc_z[k] * pgh_609[k];

        t_815[k] = f_9 * pgg0_438[k]
                   - f_10 * pgg1_438[k]
                   + f_4 * pc_x[k] * pgh_612[k];

        t_816[k] = f_1 * pfh_400[k]
                   + f_4 * pc_z[k] * pgh_610[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pa_y, pc_x, pc_y, pc_z, sgi0_397, sgi1_397, \
                         pfh_402, pgg0_441, pgg1_441, pgh_612, \
                         pgh_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = pa_y[k] * sgi0_397[k]
                   - f_11 * pc_y[k] * sgi1_397[k];

        t_818[k] = f_7 * pgg0_441[k]
                   - f_8 * pgg1_441[k]
                   + f_4 * pc_x[k] * pgh_615[k];

        t_819[k] = f_1 * pfh_402[k]
                   + f_4 * pc_z[k] * pgh_612[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, pa_y, pc_x, pc_y, sgi0_401, sgi1_401, pgg0_443, \
                         pgg0_445, pgg1_443, pgg1_445, pgh_617, \
                         pgh_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_7 * pgg0_443[k]
                   - f_8 * pgg1_443[k]
                   + f_4 * pc_x[k] * pgh_617[k];

        t_821[k] = pa_y[k] * sgi0_401[k]
                   - f_11 * pc_y[k] * sgi1_401[k];

        t_822[k] = f_5 * pgg0_445[k]
                   - f_6 * pgg1_445[k]
                   + f_4 * pc_x[k] * pgh_619[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pc_x, pc_z, pfh_405, pgg0_447, pgg0_448, \
                         pgg1_447, pgg1_448, pgh_615, pgh_621, \
                         pgh_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_1 * pfh_405[k]
                   + f_4 * pc_z[k] * pgh_615[k];

        t_824[k] = f_5 * pgg0_447[k]
                   - f_6 * pgg1_447[k]
                   + f_4 * pc_x[k] * pgh_621[k];

        t_825[k] = f_5 * pgg0_448[k]
                   - f_6 * pgg1_448[k]
                   + f_4 * pc_x[k] * pgh_622[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, t_830, t_831, pa_y, pc_x, pc_y, sgi0_406, \
                         sgi1_406, pgh_624, pgh_625, pgh_626, pgh_627, \
                         pgh_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = pa_y[k] * sgi0_406[k]
                   - f_11 * pc_y[k] * sgi1_406[k];

        t_827[k] = f_4 * pc_x[k] * pgh_624[k];

        t_828[k] = f_4 * pc_x[k] * pgh_625[k];

        t_829[k] = f_4 * pc_x[k] * pgh_626[k];

        t_830[k] = f_4 * pc_x[k] * pgh_627[k];

        t_831[k] = f_4 * pc_x[k] * pgh_628[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, pa_y, pc_x, pc_y, pc_z, sgi0_413, sgh_309, \
                         sgi1_413, pfh_414, pgh_624, pgh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_4 * pc_x[k] * pgh_629[k];

        t_833[k] = pa_y[k] * sgi0_413[k]
                   + f_16 * sgh_309[k]
                   - f_11 * pc_y[k] * sgi1_413[k];

        t_834[k] = f_1 * pfh_414[k]
                   + f_4 * pc_z[k] * pgh_624[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgi0, const size_t sgh,
                                                          const size_t sgi1, const size_t pdi0,
                                                          const size_t pdi1, const size_t pfi0,
                                                          const size_t pfh, const size_t pfi1,
                                                          const size_t pgg0, const size_t pgg1,
                                                          const size_t pgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 1.0 / p;
    const auto f_20 = gamma / (p * q);

    auto *t_835 = buffer.data(target + 835);
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
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_0 = buffer.data(sgi0 + 0);
    const auto *sgi0_2 = buffer.data(sgi0 + 2);
    const auto *sgi0_3 = buffer.data(sgi0 + 3);
    const auto *sgi0_5 = buffer.data(sgi0 + 5);
    const auto *sgi0_6 = buffer.data(sgi0 + 6);
    const auto *sgi0_7 = buffer.data(sgi0 + 7);
    const auto *sgi0_9 = buffer.data(sgi0 + 9);
    const auto *sgi0_10 = buffer.data(sgi0 + 10);
    const auto *sgi0_11 = buffer.data(sgi0 + 11);
    const auto *sgi0_12 = buffer.data(sgi0 + 12);
    const auto *sgi0_14 = buffer.data(sgi0 + 14);
    const auto *sgi0_21 = buffer.data(sgi0 + 21);
    const auto *sgi0_27 = buffer.data(sgi0 + 27);
    const auto *sgi0_28 = buffer.data(sgi0 + 28);
    const auto *sgi0_31 = buffer.data(sgi0 + 31);
    const auto *sgi0_34 = buffer.data(sgi0 + 34);
    const auto *sgi0_35 = buffer.data(sgi0 + 35);
    const auto *sgi0_38 = buffer.data(sgi0 + 38);
    const auto *sgi0_39 = buffer.data(sgi0 + 39);
    const auto *sgi0_40 = buffer.data(sgi0 + 40);
    const auto *sgi0_49 = buffer.data(sgi0 + 49);
    const auto *sgi0_50 = buffer.data(sgi0 + 50);
    const auto *sgi0_51 = buffer.data(sgi0 + 51);
    const auto *sgi0_52 = buffer.data(sgi0 + 52);
    const auto *sgi0_53 = buffer.data(sgi0 + 53);
    const auto *sgi0_84 = buffer.data(sgi0 + 84);
    const auto *sgi0_86 = buffer.data(sgi0 + 86);
    const auto *sgi0_87 = buffer.data(sgi0 + 87);
    const auto *sgi0_89 = buffer.data(sgi0 + 89);
    const auto *sgi0_90 = buffer.data(sgi0 + 90);
    const auto *sgi0_91 = buffer.data(sgi0 + 91);
    const auto *sgi0_93 = buffer.data(sgi0 + 93);
    const auto *sgi0_94 = buffer.data(sgi0 + 94);
    const auto *sgi0_95 = buffer.data(sgi0 + 95);
    const auto *sgi0_96 = buffer.data(sgi0 + 96);
    const auto *sgi0_98 = buffer.data(sgi0 + 98);
    const auto *sgi0_105 = buffer.data(sgi0 + 105);
    const auto *sgi0_106 = buffer.data(sgi0 + 106);
    const auto *sgi0_107 = buffer.data(sgi0 + 107);
    const auto *sgi0_108 = buffer.data(sgi0 + 108);
    const auto *sgi0_415 = buffer.data(sgi0 + 415);
    const auto *sgi0_416 = buffer.data(sgi0 + 416);
    const auto *sgi0_417 = buffer.data(sgi0 + 417);
    const auto *sgi0_419 = buffer.data(sgi0 + 419);

    const auto *sgh_0 = buffer.data(sgh + 0);
    const auto *sgh_2 = buffer.data(sgh + 2);
    const auto *sgh_3 = buffer.data(sgh + 3);
    const auto *sgh_5 = buffer.data(sgh + 5);
    const auto *sgh_6 = buffer.data(sgh + 6);
    const auto *sgh_7 = buffer.data(sgh + 7);
    const auto *sgh_9 = buffer.data(sgh + 9);
    const auto *sgh_20 = buffer.data(sgh + 20);
    const auto *sgh_24 = buffer.data(sgh + 24);
    const auto *sgh_27 = buffer.data(sgh + 27);
    const auto *sgh_28 = buffer.data(sgh + 28);
    const auto *sgh_36 = buffer.data(sgh + 36);
    const auto *sgh_37 = buffer.data(sgh + 37);
    const auto *sgh_38 = buffer.data(sgh + 38);
    const auto *sgh_39 = buffer.data(sgh + 39);
    const auto *sgh_63 = buffer.data(sgh + 63);
    const auto *sgh_65 = buffer.data(sgh + 65);
    const auto *sgh_66 = buffer.data(sgh + 66);
    const auto *sgh_68 = buffer.data(sgh + 68);
    const auto *sgh_69 = buffer.data(sgh + 69);
    const auto *sgh_70 = buffer.data(sgh + 70);
    const auto *sgh_72 = buffer.data(sgh + 72);
    const auto *sgh_78 = buffer.data(sgh + 78);
    const auto *sgh_79 = buffer.data(sgh + 79);
    const auto *sgh_80 = buffer.data(sgh + 80);
    const auto *sgh_311 = buffer.data(sgh + 311);
    const auto *sgh_312 = buffer.data(sgh + 312);
    const auto *sgh_313 = buffer.data(sgh + 313);
    const auto *sgh_314 = buffer.data(sgh + 314);

    const auto *sgi1_0 = buffer.data(sgi1 + 0);
    const auto *sgi1_2 = buffer.data(sgi1 + 2);
    const auto *sgi1_3 = buffer.data(sgi1 + 3);
    const auto *sgi1_5 = buffer.data(sgi1 + 5);
    const auto *sgi1_6 = buffer.data(sgi1 + 6);
    const auto *sgi1_7 = buffer.data(sgi1 + 7);
    const auto *sgi1_9 = buffer.data(sgi1 + 9);
    const auto *sgi1_10 = buffer.data(sgi1 + 10);
    const auto *sgi1_11 = buffer.data(sgi1 + 11);
    const auto *sgi1_12 = buffer.data(sgi1 + 12);
    const auto *sgi1_14 = buffer.data(sgi1 + 14);
    const auto *sgi1_21 = buffer.data(sgi1 + 21);
    const auto *sgi1_27 = buffer.data(sgi1 + 27);
    const auto *sgi1_28 = buffer.data(sgi1 + 28);
    const auto *sgi1_31 = buffer.data(sgi1 + 31);
    const auto *sgi1_34 = buffer.data(sgi1 + 34);
    const auto *sgi1_35 = buffer.data(sgi1 + 35);
    const auto *sgi1_38 = buffer.data(sgi1 + 38);
    const auto *sgi1_39 = buffer.data(sgi1 + 39);
    const auto *sgi1_40 = buffer.data(sgi1 + 40);
    const auto *sgi1_49 = buffer.data(sgi1 + 49);
    const auto *sgi1_50 = buffer.data(sgi1 + 50);
    const auto *sgi1_51 = buffer.data(sgi1 + 51);
    const auto *sgi1_52 = buffer.data(sgi1 + 52);
    const auto *sgi1_53 = buffer.data(sgi1 + 53);
    const auto *sgi1_84 = buffer.data(sgi1 + 84);
    const auto *sgi1_86 = buffer.data(sgi1 + 86);
    const auto *sgi1_87 = buffer.data(sgi1 + 87);
    const auto *sgi1_89 = buffer.data(sgi1 + 89);
    const auto *sgi1_90 = buffer.data(sgi1 + 90);
    const auto *sgi1_91 = buffer.data(sgi1 + 91);
    const auto *sgi1_93 = buffer.data(sgi1 + 93);
    const auto *sgi1_94 = buffer.data(sgi1 + 94);
    const auto *sgi1_95 = buffer.data(sgi1 + 95);
    const auto *sgi1_96 = buffer.data(sgi1 + 96);
    const auto *sgi1_98 = buffer.data(sgi1 + 98);
    const auto *sgi1_105 = buffer.data(sgi1 + 105);
    const auto *sgi1_106 = buffer.data(sgi1 + 106);
    const auto *sgi1_107 = buffer.data(sgi1 + 107);
    const auto *sgi1_108 = buffer.data(sgi1 + 108);
    const auto *sgi1_415 = buffer.data(sgi1 + 415);
    const auto *sgi1_416 = buffer.data(sgi1 + 416);
    const auto *sgi1_417 = buffer.data(sgi1 + 417);
    const auto *sgi1_419 = buffer.data(sgi1 + 419);

    const auto *pdi0_419 = buffer.data(pdi0 + 419);

    const auto *pdi1_419 = buffer.data(pdi1 + 419);

    const auto *pfi0_562 = buffer.data(pfi0 + 562);
    const auto *pfi0_565 = buffer.data(pfi0 + 565);
    const auto *pfi0_569 = buffer.data(pfi0 + 569);
    const auto *pfi0_574 = buffer.data(pfi0 + 574);
    const auto *pfi0_587 = buffer.data(pfi0 + 587);
    const auto *pfi0_643 = buffer.data(pfi0 + 643);

    const auto *pfh_420 = buffer.data(pfh + 420);
    const auto *pfh_422 = buffer.data(pfh + 422);
    const auto *pfh_425 = buffer.data(pfh + 425);
    const auto *pfh_429 = buffer.data(pfh + 429);
    const auto *pfh_435 = buffer.data(pfh + 435);
    const auto *pfh_436 = buffer.data(pfh + 436);
    const auto *pfh_437 = buffer.data(pfh + 437);
    const auto *pfh_438 = buffer.data(pfh + 438);
    const auto *pfh_439 = buffer.data(pfh + 439);
    const auto *pfh_440 = buffer.data(pfh + 440);
    const auto *pfh_441 = buffer.data(pfh + 441);
    const auto *pfh_443 = buffer.data(pfh + 443);
    const auto *pfh_446 = buffer.data(pfh + 446);
    const auto *pfh_450 = buffer.data(pfh + 450);
    const auto *pfh_456 = buffer.data(pfh + 456);
    const auto *pfh_457 = buffer.data(pfh + 457);
    const auto *pfh_458 = buffer.data(pfh + 458);
    const auto *pfh_459 = buffer.data(pfh + 459);
    const auto *pfh_460 = buffer.data(pfh + 460);
    const auto *pfh_461 = buffer.data(pfh + 461);
    const auto *pfh_462 = buffer.data(pfh + 462);
    const auto *pfh_464 = buffer.data(pfh + 464);
    const auto *pfh_465 = buffer.data(pfh + 465);
    const auto *pfh_467 = buffer.data(pfh + 467);
    const auto *pfh_468 = buffer.data(pfh + 468);
    const auto *pfh_469 = buffer.data(pfh + 469);
    const auto *pfh_471 = buffer.data(pfh + 471);
    const auto *pfh_472 = buffer.data(pfh + 472);
    const auto *pfh_473 = buffer.data(pfh + 473);
    const auto *pfh_474 = buffer.data(pfh + 474);
    const auto *pfh_476 = buffer.data(pfh + 476);
    const auto *pfh_477 = buffer.data(pfh + 477);
    const auto *pfh_478 = buffer.data(pfh + 478);
    const auto *pfh_479 = buffer.data(pfh + 479);
    const auto *pfh_480 = buffer.data(pfh + 480);
    const auto *pfh_481 = buffer.data(pfh + 481);
    const auto *pfh_482 = buffer.data(pfh + 482);
    const auto *pfh_498 = buffer.data(pfh + 498);
    const auto *pfh_499 = buffer.data(pfh + 499);
    const auto *pfh_500 = buffer.data(pfh + 500);
    const auto *pfh_501 = buffer.data(pfh + 501);
    const auto *pfh_502 = buffer.data(pfh + 502);
    const auto *pfh_503 = buffer.data(pfh + 503);

    const auto *pfi1_562 = buffer.data(pfi1 + 562);
    const auto *pfi1_565 = buffer.data(pfi1 + 565);
    const auto *pfi1_569 = buffer.data(pfi1 + 569);
    const auto *pfi1_574 = buffer.data(pfi1 + 574);
    const auto *pfi1_587 = buffer.data(pfi1 + 587);
    const auto *pfi1_643 = buffer.data(pfi1 + 643);

    const auto *pgg0_461 = buffer.data(pgg0 + 461);
    const auto *pgg0_462 = buffer.data(pgg0 + 462);
    const auto *pgg0_463 = buffer.data(pgg0 + 463);
    const auto *pgg0_464 = buffer.data(pgg0 + 464);
    const auto *pgg0_480 = buffer.data(pgg0 + 480);
    const auto *pgg0_482 = buffer.data(pgg0 + 482);
    const auto *pgg0_483 = buffer.data(pgg0 + 483);
    const auto *pgg0_485 = buffer.data(pgg0 + 485);
    const auto *pgg0_486 = buffer.data(pgg0 + 486);
    const auto *pgg0_487 = buffer.data(pgg0 + 487);
    const auto *pgg0_489 = buffer.data(pgg0 + 489);
    const auto *pgg0_490 = buffer.data(pgg0 + 490);
    const auto *pgg0_491 = buffer.data(pgg0 + 491);
    const auto *pgg0_492 = buffer.data(pgg0 + 492);
    const auto *pgg0_493 = buffer.data(pgg0 + 493);
    const auto *pgg0_494 = buffer.data(pgg0 + 494);

    const auto *pgg1_461 = buffer.data(pgg1 + 461);
    const auto *pgg1_462 = buffer.data(pgg1 + 462);
    const auto *pgg1_463 = buffer.data(pgg1 + 463);
    const auto *pgg1_464 = buffer.data(pgg1 + 464);
    const auto *pgg1_480 = buffer.data(pgg1 + 480);
    const auto *pgg1_482 = buffer.data(pgg1 + 482);
    const auto *pgg1_483 = buffer.data(pgg1 + 483);
    const auto *pgg1_485 = buffer.data(pgg1 + 485);
    const auto *pgg1_486 = buffer.data(pgg1 + 486);
    const auto *pgg1_487 = buffer.data(pgg1 + 487);
    const auto *pgg1_489 = buffer.data(pgg1 + 489);
    const auto *pgg1_490 = buffer.data(pgg1 + 490);
    const auto *pgg1_491 = buffer.data(pgg1 + 491);
    const auto *pgg1_492 = buffer.data(pgg1 + 492);
    const auto *pgg1_493 = buffer.data(pgg1 + 493);
    const auto *pgg1_494 = buffer.data(pgg1 + 494);

    const auto *pgh_629 = buffer.data(pgh + 629);
    const auto *pgh_630 = buffer.data(pgh + 630);
    const auto *pgh_632 = buffer.data(pgh + 632);
    const auto *pgh_635 = buffer.data(pgh + 635);
    const auto *pgh_639 = buffer.data(pgh + 639);
    const auto *pgh_645 = buffer.data(pgh + 645);
    const auto *pgh_646 = buffer.data(pgh + 646);
    const auto *pgh_647 = buffer.data(pgh + 647);
    const auto *pgh_648 = buffer.data(pgh + 648);
    const auto *pgh_649 = buffer.data(pgh + 649);
    const auto *pgh_650 = buffer.data(pgh + 650);
    const auto *pgh_651 = buffer.data(pgh + 651);
    const auto *pgh_653 = buffer.data(pgh + 653);
    const auto *pgh_656 = buffer.data(pgh + 656);
    const auto *pgh_660 = buffer.data(pgh + 660);
    const auto *pgh_666 = buffer.data(pgh + 666);
    const auto *pgh_667 = buffer.data(pgh + 667);
    const auto *pgh_668 = buffer.data(pgh + 668);
    const auto *pgh_669 = buffer.data(pgh + 669);
    const auto *pgh_670 = buffer.data(pgh + 670);
    const auto *pgh_671 = buffer.data(pgh + 671);
    const auto *pgh_672 = buffer.data(pgh + 672);
    const auto *pgh_674 = buffer.data(pgh + 674);
    const auto *pgh_675 = buffer.data(pgh + 675);
    const auto *pgh_677 = buffer.data(pgh + 677);
    const auto *pgh_678 = buffer.data(pgh + 678);
    const auto *pgh_679 = buffer.data(pgh + 679);
    const auto *pgh_681 = buffer.data(pgh + 681);
    const auto *pgh_682 = buffer.data(pgh + 682);
    const auto *pgh_683 = buffer.data(pgh + 683);
    const auto *pgh_684 = buffer.data(pgh + 684);
    const auto *pgh_686 = buffer.data(pgh + 686);
    const auto *pgh_687 = buffer.data(pgh + 687);
    const auto *pgh_688 = buffer.data(pgh + 688);
    const auto *pgh_689 = buffer.data(pgh + 689);
    const auto *pgh_690 = buffer.data(pgh + 690);
    const auto *pgh_691 = buffer.data(pgh + 691);
    const auto *pgh_692 = buffer.data(pgh + 692);
    const auto *pgh_693 = buffer.data(pgh + 693);
    const auto *pgh_695 = buffer.data(pgh + 695);
    const auto *pgh_698 = buffer.data(pgh + 698);
    const auto *pgh_702 = buffer.data(pgh + 702);
    const auto *pgh_708 = buffer.data(pgh + 708);
    const auto *pgh_709 = buffer.data(pgh + 709);
    const auto *pgh_710 = buffer.data(pgh + 710);
    const auto *pgh_711 = buffer.data(pgh + 711);
    const auto *pgh_712 = buffer.data(pgh + 712);
    const auto *pgh_713 = buffer.data(pgh + 713);

#pragma omp simd aligned(t_835, t_836, t_837, pa_y, pc_y, sgi0_415, sgi0_416, sgi0_417, \
                         sgh_311, sgh_312, sgh_313, sgi1_415, sgi1_416, \
                         sgi1_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = pa_y[k] * sgi0_415[k]
                   + f_1 * sgh_311[k]
                   - f_11 * pc_y[k] * sgi1_415[k];

        t_836[k] = pa_y[k] * sgi0_416[k]
                   + f_13 * sgh_312[k]
                   - f_11 * pc_y[k] * sgi1_416[k];

        t_837[k] = pa_y[k] * sgi0_417[k]
                   + f_12 * sgh_313[k]
                   - f_11 * pc_y[k] * sgi1_417[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, t_841, pa_y, pa_z, pc_y, pc_z, sgi0_0, sgi0_419, \
                         sgh_314, sgi1_0, sgi1_419, pgh_629, pgh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_0 * sgh_314[k]
                   + f_4 * pc_y[k] * pgh_629[k];

        t_839[k] = pa_y[k] * sgi0_419[k]
                   - f_11 * pc_y[k] * sgi1_419[k];

        t_840[k] = pa_z[k] * sgi0_0[k]
                   - f_11 * pc_z[k] * sgi1_0[k];

        t_841[k] = f_4 * pc_y[k] * pgh_630[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, t_845, pa_z, pc_y, pc_z, sgi0_2, sgi0_3, sgi0_5, \
                         sgh_0, sgh_2, sgi1_2, sgi1_3, sgi1_5, \
                         pgh_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = pa_z[k] * sgi0_2[k]
                   + f_0 * sgh_0[k]
                   - f_11 * pc_z[k] * sgi1_2[k];

        t_843[k] = pa_z[k] * sgi0_3[k]
                   - f_11 * pc_z[k] * sgi1_3[k];

        t_844[k] = f_4 * pc_y[k] * pgh_632[k];

        t_845[k] = pa_z[k] * sgi0_5[k]
                   + f_12 * sgh_2[k]
                   - f_11 * pc_z[k] * sgi1_5[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, t_849, pa_z, pc_y, pc_z, sgi0_6, sgi0_7, sgi0_9, \
                         sgh_3, sgh_5, sgi1_6, sgi1_7, sgi1_9, \
                         pgh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = pa_z[k] * sgi0_6[k]
                   - f_11 * pc_z[k] * sgi1_6[k];

        t_847[k] = pa_z[k] * sgi0_7[k]
                   + f_0 * sgh_3[k]
                   - f_11 * pc_z[k] * sgi1_7[k];

        t_848[k] = f_4 * pc_y[k] * pgh_635[k];

        t_849[k] = pa_z[k] * sgi0_9[k]
                   + f_13 * sgh_5[k]
                   - f_11 * pc_z[k] * sgi1_9[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, pa_z, pc_y, pc_z, sgi0_10, sgi0_11, \
                         sgi0_12, sgh_6, sgh_7, sgi1_10, sgi1_11, sgi1_12, \
                         pgh_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = pa_z[k] * sgi0_10[k]
                   - f_11 * pc_z[k] * sgi1_10[k];

        t_851[k] = pa_z[k] * sgi0_11[k]
                   + f_0 * sgh_6[k]
                   - f_11 * pc_z[k] * sgi1_11[k];

        t_852[k] = pa_z[k] * sgi0_12[k]
                   + f_12 * sgh_7[k]
                   - f_11 * pc_z[k] * sgi1_12[k];

        t_853[k] = f_4 * pc_y[k] * pgh_639[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pa_z, pc_x, pc_z, sgi0_14, sgh_9, \
                         sgi1_14, pfh_435, pfh_436, pfh_437, pgh_645, pgh_646, \
                         pgh_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = pa_z[k] * sgi0_14[k]
                   + f_1 * sgh_9[k]
                   - f_11 * pc_z[k] * sgi1_14[k];

        t_855[k] = f_1 * pfh_435[k]
                   + f_4 * pc_x[k] * pgh_645[k];

        t_856[k] = f_1 * pfh_436[k]
                   + f_4 * pc_x[k] * pgh_646[k];

        t_857[k] = f_1 * pfh_437[k]
                   + f_4 * pc_x[k] * pgh_647[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pa_z, pc_x, pc_z, sgi0_21, sgi1_21, \
                         pfh_438, pfh_439, pfh_440, pgh_648, pgh_649, \
                         pgh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_1 * pfh_438[k]
                   + f_4 * pc_x[k] * pgh_648[k];

        t_859[k] = f_1 * pfh_439[k]
                   + f_4 * pc_x[k] * pgh_649[k];

        t_860[k] = f_1 * pfh_440[k]
                   + f_4 * pc_x[k] * pgh_650[k];

        t_861[k] = pa_z[k] * sgi0_21[k]
                   - f_11 * pc_z[k] * sgi1_21[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, pc_y, pgg0_461, pgg0_462, pgg0_463, pgg1_461, \
                         pgg1_462, pgg1_463, pgh_646, pgh_647, \
                         pgh_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_17 * pgg0_461[k]
                   - f_18 * pgg1_461[k]
                   + f_4 * pc_y[k] * pgh_646[k];

        t_863[k] = f_9 * pgg0_462[k]
                   - f_10 * pgg1_462[k]
                   + f_4 * pc_y[k] * pgh_647[k];

        t_864[k] = f_7 * pgg0_463[k]
                   - f_8 * pgg1_463[k]
                   + f_4 * pc_y[k] * pgh_648[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pa_z, pc_y, pc_z, sgi0_27, sgi0_28, \
                         sgh_20, sgi1_27, sgi1_28, pgg0_464, pgg1_464, pgh_649, \
                         pgh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_5 * pgg0_464[k]
                   - f_6 * pgg1_464[k]
                   + f_4 * pc_y[k] * pgh_649[k];

        t_866[k] = f_4 * pc_y[k] * pgh_650[k];

        t_867[k] = pa_z[k] * sgi0_27[k]
                   + f_16 * sgh_20[k]
                   - f_11 * pc_z[k] * sgi1_27[k];

        t_868[k] = pa_z[k] * sgi0_28[k]
                   - f_11 * pc_z[k] * sgi1_28[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pa_z, pb_y, pc_y, pc_z, sgi0_31, sgi1_31, \
                         pfi0_562, pfh_420, pfh_422, pfi1_562, pgh_651, \
                         pgh_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_0 * pfh_420[k]
                   + f_4 * pc_y[k] * pgh_651[k];

        t_870[k] = pb_y[k] * pfi0_562[k]
                   - f_11 * pc_y[k] * pfi1_562[k];

        t_871[k] = pa_z[k] * sgi0_31[k]
                   - f_11 * pc_z[k] * sgi1_31[k];

        t_872[k] = f_0 * pfh_422[k]
                   + f_4 * pc_y[k] * pgh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pa_z, pb_y, pc_y, pc_z, sgi0_34, sgi0_35, \
                         sgh_24, sgi1_34, sgi1_35, pfi0_565, pfi1_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pb_y[k] * pfi0_565[k]
                   - f_11 * pc_y[k] * pfi1_565[k];

        t_874[k] = pa_z[k] * sgi0_34[k]
                   - f_11 * pc_z[k] * sgi1_34[k];

        t_875[k] = pa_z[k] * sgi0_35[k]
                   + f_0 * sgh_24[k]
                   - f_11 * pc_z[k] * sgi1_35[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pa_z, pb_y, pc_y, pc_z, sgi0_38, sgi1_38, \
                         pfi0_569, pfh_425, pfi1_569, pgh_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_0 * pfh_425[k]
                   + f_4 * pc_y[k] * pgh_656[k];

        t_877[k] = pb_y[k] * pfi0_569[k]
                   - f_11 * pc_y[k] * pfi1_569[k];

        t_878[k] = pa_z[k] * sgi0_38[k]
                   - f_11 * pc_z[k] * sgi1_38[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pa_z, pc_y, pc_z, sgi0_39, sgi0_40, sgh_27, \
                         sgh_28, sgi1_39, sgi1_40, pfh_429, pgh_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = pa_z[k] * sgi0_39[k]
                   + f_0 * sgh_27[k]
                   - f_11 * pc_z[k] * sgi1_39[k];

        t_880[k] = pa_z[k] * sgi0_40[k]
                   + f_12 * sgh_28[k]
                   - f_11 * pc_z[k] * sgi1_40[k];

        t_881[k] = f_0 * pfh_429[k]
                   + f_4 * pc_y[k] * pgh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pb_y, pc_x, pc_y, pfi0_574, pfh_456, \
                         pfh_457, pfh_458, pfi1_574, pgh_666, pgh_667, \
                         pgh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = pb_y[k] * pfi0_574[k]
                   - f_11 * pc_y[k] * pfi1_574[k];

        t_883[k] = f_13 * pfh_456[k]
                   + f_4 * pc_x[k] * pgh_666[k];

        t_884[k] = f_13 * pfh_457[k]
                   + f_4 * pc_x[k] * pgh_667[k];

        t_885[k] = f_13 * pfh_458[k]
                   + f_4 * pc_x[k] * pgh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pa_z, pc_x, pc_z, sgi0_49, sgi1_49, \
                         pfh_459, pfh_460, pfh_461, pgh_669, pgh_670, \
                         pgh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_13 * pfh_459[k]
                   + f_4 * pc_x[k] * pgh_669[k];

        t_887[k] = f_13 * pfh_460[k]
                   + f_4 * pc_x[k] * pgh_670[k];

        t_888[k] = f_13 * pfh_461[k]
                   + f_4 * pc_x[k] * pgh_671[k];

        t_889[k] = pa_z[k] * sgi0_49[k]
                   - f_11 * pc_z[k] * sgi1_49[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, pa_z, pc_z, sgi0_50, sgi0_51, sgi0_52, sgh_36, \
                         sgh_37, sgh_38, sgi1_50, sgi1_51, sgi1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = pa_z[k] * sgi0_50[k]
                   + f_0 * sgh_36[k]
                   - f_11 * pc_z[k] * sgi1_50[k];

        t_891[k] = pa_z[k] * sgi0_51[k]
                   + f_12 * sgh_37[k]
                   - f_11 * pc_z[k] * sgi1_51[k];

        t_892[k] = pa_z[k] * sgi0_52[k]
                   + f_13 * sgh_38[k]
                   - f_11 * pc_z[k] * sgi1_52[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pa_z, pb_y, pc_y, pc_z, sgi0_53, sgh_39, \
                         sgi1_53, pfi0_587, pfh_440, pfi1_587, \
                         pgh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = pa_z[k] * sgi0_53[k]
                   + f_1 * sgh_39[k]
                   - f_11 * pc_z[k] * sgi1_53[k];

        t_894[k] = f_0 * pfh_440[k]
                   + f_4 * pc_y[k] * pgh_671[k];

        t_895[k] = pb_y[k] * pfi0_587[k]
                   - f_11 * pc_y[k] * pfi1_587[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_x, pc_y, pfh_462, pfh_464, pgg0_480, \
                         pgg0_482, pgg1_480, pgg1_482, pgh_672, \
                         pgh_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_13 * pfh_462[k]
                   + f_2 * pgg0_480[k]
                   - f_3 * pgg1_480[k]
                   + f_4 * pc_x[k] * pgh_672[k];

        t_897[k] = f_4 * pc_y[k] * pgh_672[k];

        t_898[k] = f_13 * pfh_464[k]
                   + f_17 * pgg0_482[k]
                   - f_18 * pgg1_482[k]
                   + f_4 * pc_x[k] * pgh_674[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_y, pfh_465, pfh_467, pgg0_483, \
                         pgg0_485, pgg1_483, pgg1_485, pgh_674, pgh_675, \
                         pgh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_13 * pfh_465[k]
                   + f_9 * pgg0_483[k]
                   - f_10 * pgg1_483[k]
                   + f_4 * pc_x[k] * pgh_675[k];

        t_900[k] = f_4 * pc_y[k] * pgh_674[k];

        t_901[k] = f_13 * pfh_467[k]
                   + f_9 * pgg0_485[k]
                   - f_10 * pgg1_485[k]
                   + f_4 * pc_x[k] * pgh_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pfh_468, pfh_469, pgg0_486, \
                         pgg0_487, pgg1_486, pgg1_487, pgh_677, pgh_678, \
                         pgh_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_13 * pfh_468[k]
                   + f_7 * pgg0_486[k]
                   - f_8 * pgg1_486[k]
                   + f_4 * pc_x[k] * pgh_678[k];

        t_903[k] = f_13 * pfh_469[k]
                   + f_7 * pgg0_487[k]
                   - f_8 * pgg1_487[k]
                   + f_4 * pc_x[k] * pgh_679[k];

        t_904[k] = f_4 * pc_y[k] * pgh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pfh_471, pfh_472, pfh_473, pgg0_489, \
                         pgg0_490, pgg0_491, pgg1_489, pgg1_490, pgg1_491, pgh_681, pgh_682, \
                         pgh_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_13 * pfh_471[k]
                   + f_7 * pgg0_489[k]
                   - f_8 * pgg1_489[k]
                   + f_4 * pc_x[k] * pgh_681[k];

        t_906[k] = f_13 * pfh_472[k]
                   + f_5 * pgg0_490[k]
                   - f_6 * pgg1_490[k]
                   + f_4 * pc_x[k] * pgh_682[k];

        t_907[k] = f_13 * pfh_473[k]
                   + f_5 * pgg0_491[k]
                   - f_6 * pgg1_491[k]
                   + f_4 * pc_x[k] * pgh_683[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, pfh_474, pfh_476, pgg0_492, \
                         pgg0_494, pgg1_492, pgg1_494, pgh_681, pgh_684, \
                         pgh_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_13 * pfh_474[k]
                   + f_5 * pgg0_492[k]
                   - f_6 * pgg1_492[k]
                   + f_4 * pc_x[k] * pgh_684[k];

        t_909[k] = f_4 * pc_y[k] * pgh_681[k];

        t_910[k] = f_13 * pfh_476[k]
                   + f_5 * pgg0_494[k]
                   - f_6 * pgg1_494[k]
                   + f_4 * pc_x[k] * pgh_686[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pc_x, pfh_477, pfh_478, pfh_479, \
                         pfh_480, pfh_481, pgh_687, pgh_688, pgh_689, pgh_690, \
                         pgh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_13 * pfh_477[k]
                   + f_4 * pc_x[k] * pgh_687[k];

        t_912[k] = f_13 * pfh_478[k]
                   + f_4 * pc_x[k] * pgh_688[k];

        t_913[k] = f_13 * pfh_479[k]
                   + f_4 * pc_x[k] * pgh_689[k];

        t_914[k] = f_13 * pfh_480[k]
                   + f_4 * pc_x[k] * pgh_690[k];

        t_915[k] = f_13 * pfh_481[k]
                   + f_4 * pc_x[k] * pgh_691[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, pc_x, pc_y, pfh_482, pgg0_490, pgg0_491, \
                         pgg1_490, pgg1_491, pgh_687, pgh_688, \
                         pgh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_13 * pfh_482[k]
                   + f_4 * pc_x[k] * pgh_692[k];

        t_917[k] = f_2 * pgg0_490[k]
                   - f_3 * pgg1_490[k]
                   + f_4 * pc_y[k] * pgh_687[k];

        t_918[k] = f_17 * pgg0_491[k]
                   - f_18 * pgg1_491[k]
                   + f_4 * pc_y[k] * pgh_688[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, t_922, pc_y, pgg0_492, pgg0_493, pgg0_494, \
                         pgg1_492, pgg1_493, pgg1_494, pgh_689, pgh_690, pgh_691, \
                         pgh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_9 * pgg0_492[k]
                   - f_10 * pgg1_492[k]
                   + f_4 * pc_y[k] * pgh_689[k];

        t_920[k] = f_7 * pgg0_493[k]
                   - f_8 * pgg1_493[k]
                   + f_4 * pc_y[k] * pgh_690[k];

        t_921[k] = f_5 * pgg0_494[k]
                   - f_6 * pgg1_494[k]
                   + f_4 * pc_y[k] * pgh_691[k];

        t_922[k] = f_4 * pc_y[k] * pgh_692[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, pa_z, pb_x, pc_x, pc_y, pc_z, sgi0_84, sgi1_84, \
                         pdi0_419, pdi1_419, pfi0_643, pfh_441, pfi1_643, \
                         pgh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = f_19 * pdi0_419[k]
                   - f_20 * pdi1_419[k]
                   + pb_x[k] * pfi0_643[k]
                   - f_11 * pc_x[k] * pfi1_643[k];

        t_924[k] = pa_z[k] * sgi0_84[k]
                   - f_11 * pc_z[k] * sgi1_84[k];

        t_925[k] = f_12 * pfh_441[k]
                   + f_4 * pc_y[k] * pgh_693[k];
    }

#pragma omp simd aligned(t_926, t_927, t_928, pa_z, pc_y, pc_z, sgi0_86, sgi0_87, sgh_63, \
                         sgi1_86, sgi1_87, pfh_443, pgh_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_926[k] = pa_z[k] * sgi0_86[k]
                   + f_0 * sgh_63[k]
                   - f_11 * pc_z[k] * sgi1_86[k];

        t_927[k] = pa_z[k] * sgi0_87[k]
                   - f_11 * pc_z[k] * sgi1_87[k];

        t_928[k] = f_12 * pfh_443[k]
                   + f_4 * pc_y[k] * pgh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pa_z, pc_z, sgi0_89, sgi0_90, sgi0_91, sgh_65, \
                         sgh_66, sgi1_89, sgi1_90, sgi1_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = pa_z[k] * sgi0_89[k]
                   + f_12 * sgh_65[k]
                   - f_11 * pc_z[k] * sgi1_89[k];

        t_930[k] = pa_z[k] * sgi0_90[k]
                   - f_11 * pc_z[k] * sgi1_90[k];

        t_931[k] = pa_z[k] * sgi0_91[k]
                   + f_0 * sgh_66[k]
                   - f_11 * pc_z[k] * sgi1_91[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pa_z, pc_y, pc_z, sgi0_93, sgi0_94, sgh_68, \
                         sgi1_93, sgi1_94, pfh_446, pgh_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * pfh_446[k]
                   + f_4 * pc_y[k] * pgh_698[k];

        t_933[k] = pa_z[k] * sgi0_93[k]
                   + f_13 * sgh_68[k]
                   - f_11 * pc_z[k] * sgi1_93[k];

        t_934[k] = pa_z[k] * sgi0_94[k]
                   - f_11 * pc_z[k] * sgi1_94[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pa_z, pc_y, pc_z, sgi0_95, sgi0_96, sgh_69, \
                         sgh_70, sgi1_95, sgi1_96, pfh_450, pgh_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = pa_z[k] * sgi0_95[k]
                   + f_0 * sgh_69[k]
                   - f_11 * pc_z[k] * sgi1_95[k];

        t_936[k] = pa_z[k] * sgi0_96[k]
                   + f_12 * sgh_70[k]
                   - f_11 * pc_z[k] * sgi1_96[k];

        t_937[k] = f_12 * pfh_450[k]
                   + f_4 * pc_y[k] * pgh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pa_z, pc_x, pc_z, sgi0_98, sgh_72, \
                         sgi1_98, pfh_498, pfh_499, pfh_500, pgh_708, pgh_709, \
                         pgh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = pa_z[k] * sgi0_98[k]
                   + f_1 * sgh_72[k]
                   - f_11 * pc_z[k] * sgi1_98[k];

        t_939[k] = f_12 * pfh_498[k]
                   + f_4 * pc_x[k] * pgh_708[k];

        t_940[k] = f_12 * pfh_499[k]
                   + f_4 * pc_x[k] * pgh_709[k];

        t_941[k] = f_12 * pfh_500[k]
                   + f_4 * pc_x[k] * pgh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pa_z, pc_x, pc_z, sgi0_105, sgi1_105, \
                         pfh_501, pfh_502, pfh_503, pgh_711, pgh_712, \
                         pgh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_12 * pfh_501[k]
                   + f_4 * pc_x[k] * pgh_711[k];

        t_943[k] = f_12 * pfh_502[k]
                   + f_4 * pc_x[k] * pgh_712[k];

        t_944[k] = f_12 * pfh_503[k]
                   + f_4 * pc_x[k] * pgh_713[k];

        t_945[k] = pa_z[k] * sgi0_105[k]
                   - f_11 * pc_z[k] * sgi1_105[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pa_z, pc_z, sgi0_106, sgi0_107, sgi0_108, \
                         sgh_78, sgh_79, sgh_80, sgi1_106, sgi1_107, \
                         sgi1_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = pa_z[k] * sgi0_106[k]
                   + f_0 * sgh_78[k]
                   - f_11 * pc_z[k] * sgi1_106[k];

        t_947[k] = pa_z[k] * sgi0_107[k]
                   + f_12 * sgh_79[k]
                   - f_11 * pc_z[k] * sgi1_107[k];

        t_948[k] = pa_z[k] * sgi0_108[k]
                   + f_13 * sgh_80[k]
                   - f_11 * pc_z[k] * sgi1_108[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgi0, const size_t sgh,
                                                          const size_t sgi1, const size_t pdi0,
                                                          const size_t pdi1, const size_t pfi0,
                                                          const size_t pfh, const size_t pfi1,
                                                          const size_t pgg0, const size_t pgg1,
                                                          const size_t pgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_21 = 2.5 / q;

    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_109 = buffer.data(sgi0 + 109);
    const auto *sgi0_111 = buffer.data(sgi0 + 111);
    const auto *sgi0_168 = buffer.data(sgi0 + 168);
    const auto *sgi0_170 = buffer.data(sgi0 + 170);
    const auto *sgi0_171 = buffer.data(sgi0 + 171);
    const auto *sgi0_173 = buffer.data(sgi0 + 173);
    const auto *sgi0_174 = buffer.data(sgi0 + 174);
    const auto *sgi0_175 = buffer.data(sgi0 + 175);
    const auto *sgi0_177 = buffer.data(sgi0 + 177);
    const auto *sgi0_178 = buffer.data(sgi0 + 178);
    const auto *sgi0_179 = buffer.data(sgi0 + 179);
    const auto *sgi0_180 = buffer.data(sgi0 + 180);
    const auto *sgi0_182 = buffer.data(sgi0 + 182);
    const auto *sgi0_189 = buffer.data(sgi0 + 189);

    const auto *sgh_81 = buffer.data(sgh + 81);
    const auto *sgh_83 = buffer.data(sgh + 83);
    const auto *sgh_126 = buffer.data(sgh + 126);
    const auto *sgh_128 = buffer.data(sgh + 128);
    const auto *sgh_129 = buffer.data(sgh + 129);
    const auto *sgh_131 = buffer.data(sgh + 131);
    const auto *sgh_132 = buffer.data(sgh + 132);
    const auto *sgh_133 = buffer.data(sgh + 133);
    const auto *sgh_135 = buffer.data(sgh + 135);

    const auto *sgi1_109 = buffer.data(sgi1 + 109);
    const auto *sgi1_111 = buffer.data(sgi1 + 111);
    const auto *sgi1_168 = buffer.data(sgi1 + 168);
    const auto *sgi1_170 = buffer.data(sgi1 + 170);
    const auto *sgi1_171 = buffer.data(sgi1 + 171);
    const auto *sgi1_173 = buffer.data(sgi1 + 173);
    const auto *sgi1_174 = buffer.data(sgi1 + 174);
    const auto *sgi1_175 = buffer.data(sgi1 + 175);
    const auto *sgi1_177 = buffer.data(sgi1 + 177);
    const auto *sgi1_178 = buffer.data(sgi1 + 178);
    const auto *sgi1_179 = buffer.data(sgi1 + 179);
    const auto *sgi1_180 = buffer.data(sgi1 + 180);
    const auto *sgi1_182 = buffer.data(sgi1 + 182);
    const auto *sgi1_189 = buffer.data(sgi1 + 189);

    const auto *pdi0_503 = buffer.data(pdi0 + 503);

    const auto *pdi1_503 = buffer.data(pdi1 + 503);

    const auto *pfi0_616 = buffer.data(pfi0 + 616);
    const auto *pfi0_618 = buffer.data(pfi0 + 618);
    const auto *pfi0_621 = buffer.data(pfi0 + 621);
    const auto *pfi0_623 = buffer.data(pfi0 + 623);
    const auto *pfi0_625 = buffer.data(pfi0 + 625);
    const auto *pfi0_627 = buffer.data(pfi0 + 627);
    const auto *pfi0_628 = buffer.data(pfi0 + 628);
    const auto *pfi0_630 = buffer.data(pfi0 + 630);
    const auto *pfi0_643 = buffer.data(pfi0 + 643);
    const auto *pfi0_727 = buffer.data(pfi0 + 727);
    const auto *pfi0_750 = buffer.data(pfi0 + 750);
    const auto *pfi0_751 = buffer.data(pfi0 + 751);
    const auto *pfi0_752 = buffer.data(pfi0 + 752);
    const auto *pfi0_753 = buffer.data(pfi0 + 753);
    const auto *pfi0_755 = buffer.data(pfi0 + 755);
    const auto *pfi0_758 = buffer.data(pfi0 + 758);
    const auto *pfi0_761 = buffer.data(pfi0 + 761);
    const auto *pfi0_763 = buffer.data(pfi0 + 763);
    const auto *pfi0_765 = buffer.data(pfi0 + 765);
    const auto *pfi0_767 = buffer.data(pfi0 + 767);
    const auto *pfi0_768 = buffer.data(pfi0 + 768);
    const auto *pfi0_770 = buffer.data(pfi0 + 770);
    const auto *pfi0_777 = buffer.data(pfi0 + 777);
    const auto *pfi0_778 = buffer.data(pfi0 + 778);
    const auto *pfi0_779 = buffer.data(pfi0 + 779);

    const auto *pfh_461 = buffer.data(pfh + 461);
    const auto *pfh_462 = buffer.data(pfh + 462);
    const auto *pfh_464 = buffer.data(pfh + 464);
    const auto *pfh_466 = buffer.data(pfh + 466);
    const auto *pfh_467 = buffer.data(pfh + 467);
    const auto *pfh_469 = buffer.data(pfh + 469);
    const auto *pfh_470 = buffer.data(pfh + 470);
    const auto *pfh_471 = buffer.data(pfh + 471);
    const auto *pfh_477 = buffer.data(pfh + 477);
    const auto *pfh_478 = buffer.data(pfh + 478);
    const auto *pfh_479 = buffer.data(pfh + 479);
    const auto *pfh_480 = buffer.data(pfh + 480);
    const auto *pfh_481 = buffer.data(pfh + 481);
    const auto *pfh_482 = buffer.data(pfh + 482);
    const auto *pfh_483 = buffer.data(pfh + 483);
    const auto *pfh_485 = buffer.data(pfh + 485);
    const auto *pfh_488 = buffer.data(pfh + 488);
    const auto *pfh_492 = buffer.data(pfh + 492);
    const auto *pfh_503 = buffer.data(pfh + 503);
    const auto *pfh_504 = buffer.data(pfh + 504);
    const auto *pfh_506 = buffer.data(pfh + 506);
    const auto *pfh_507 = buffer.data(pfh + 507);
    const auto *pfh_509 = buffer.data(pfh + 509);
    const auto *pfh_510 = buffer.data(pfh + 510);
    const auto *pfh_513 = buffer.data(pfh + 513);
    const auto *pfh_514 = buffer.data(pfh + 514);
    const auto *pfh_519 = buffer.data(pfh + 519);
    const auto *pfh_520 = buffer.data(pfh + 520);
    const auto *pfh_521 = buffer.data(pfh + 521);
    const auto *pfh_522 = buffer.data(pfh + 522);
    const auto *pfh_523 = buffer.data(pfh + 523);
    const auto *pfh_524 = buffer.data(pfh + 524);
    const auto *pfh_525 = buffer.data(pfh + 525);
    const auto *pfh_527 = buffer.data(pfh + 527);
    const auto *pfh_528 = buffer.data(pfh + 528);
    const auto *pfh_530 = buffer.data(pfh + 530);
    const auto *pfh_531 = buffer.data(pfh + 531);
    const auto *pfh_532 = buffer.data(pfh + 532);
    const auto *pfh_534 = buffer.data(pfh + 534);
    const auto *pfh_535 = buffer.data(pfh + 535);
    const auto *pfh_536 = buffer.data(pfh + 536);
    const auto *pfh_537 = buffer.data(pfh + 537);
    const auto *pfh_539 = buffer.data(pfh + 539);
    const auto *pfh_540 = buffer.data(pfh + 540);
    const auto *pfh_541 = buffer.data(pfh + 541);
    const auto *pfh_542 = buffer.data(pfh + 542);
    const auto *pfh_543 = buffer.data(pfh + 543);
    const auto *pfh_544 = buffer.data(pfh + 544);
    const auto *pfh_545 = buffer.data(pfh + 545);
    const auto *pfh_561 = buffer.data(pfh + 561);
    const auto *pfh_562 = buffer.data(pfh + 562);
    const auto *pfh_563 = buffer.data(pfh + 563);
    const auto *pfh_564 = buffer.data(pfh + 564);
    const auto *pfh_565 = buffer.data(pfh + 565);
    const auto *pfh_566 = buffer.data(pfh + 566);
    const auto *pfh_567 = buffer.data(pfh + 567);
    const auto *pfh_569 = buffer.data(pfh + 569);
    const auto *pfh_570 = buffer.data(pfh + 570);
    const auto *pfh_572 = buffer.data(pfh + 572);
    const auto *pfh_573 = buffer.data(pfh + 573);
    const auto *pfh_574 = buffer.data(pfh + 574);
    const auto *pfh_576 = buffer.data(pfh + 576);
    const auto *pfh_577 = buffer.data(pfh + 577);
    const auto *pfh_578 = buffer.data(pfh + 578);
    const auto *pfh_579 = buffer.data(pfh + 579);
    const auto *pfh_581 = buffer.data(pfh + 581);
    const auto *pfh_582 = buffer.data(pfh + 582);
    const auto *pfh_583 = buffer.data(pfh + 583);
    const auto *pfh_584 = buffer.data(pfh + 584);
    const auto *pfh_585 = buffer.data(pfh + 585);
    const auto *pfh_586 = buffer.data(pfh + 586);
    const auto *pfh_587 = buffer.data(pfh + 587);

    const auto *pfi1_616 = buffer.data(pfi1 + 616);
    const auto *pfi1_618 = buffer.data(pfi1 + 618);
    const auto *pfi1_621 = buffer.data(pfi1 + 621);
    const auto *pfi1_623 = buffer.data(pfi1 + 623);
    const auto *pfi1_625 = buffer.data(pfi1 + 625);
    const auto *pfi1_627 = buffer.data(pfi1 + 627);
    const auto *pfi1_628 = buffer.data(pfi1 + 628);
    const auto *pfi1_630 = buffer.data(pfi1 + 630);
    const auto *pfi1_643 = buffer.data(pfi1 + 643);
    const auto *pfi1_727 = buffer.data(pfi1 + 727);
    const auto *pfi1_750 = buffer.data(pfi1 + 750);
    const auto *pfi1_751 = buffer.data(pfi1 + 751);
    const auto *pfi1_752 = buffer.data(pfi1 + 752);
    const auto *pfi1_753 = buffer.data(pfi1 + 753);
    const auto *pfi1_755 = buffer.data(pfi1 + 755);
    const auto *pfi1_758 = buffer.data(pfi1 + 758);
    const auto *pfi1_761 = buffer.data(pfi1 + 761);
    const auto *pfi1_763 = buffer.data(pfi1 + 763);
    const auto *pfi1_765 = buffer.data(pfi1 + 765);
    const auto *pfi1_767 = buffer.data(pfi1 + 767);
    const auto *pfi1_768 = buffer.data(pfi1 + 768);
    const auto *pfi1_770 = buffer.data(pfi1 + 770);
    const auto *pfi1_777 = buffer.data(pfi1 + 777);
    const auto *pfi1_778 = buffer.data(pfi1 + 778);
    const auto *pfi1_779 = buffer.data(pfi1 + 779);

    const auto *pgg0_513 = buffer.data(pgg0 + 513);
    const auto *pgg0_516 = buffer.data(pgg0 + 516);
    const auto *pgg0_520 = buffer.data(pgg0 + 520);
    const auto *pgg0_521 = buffer.data(pgg0 + 521);
    const auto *pgg0_522 = buffer.data(pgg0 + 522);
    const auto *pgg0_523 = buffer.data(pgg0 + 523);
    const auto *pgg0_524 = buffer.data(pgg0 + 524);
    const auto *pgg0_525 = buffer.data(pgg0 + 525);
    const auto *pgg0_527 = buffer.data(pgg0 + 527);
    const auto *pgg0_528 = buffer.data(pgg0 + 528);
    const auto *pgg0_530 = buffer.data(pgg0 + 530);
    const auto *pgg0_531 = buffer.data(pgg0 + 531);
    const auto *pgg0_532 = buffer.data(pgg0 + 532);
    const auto *pgg0_534 = buffer.data(pgg0 + 534);
    const auto *pgg0_535 = buffer.data(pgg0 + 535);
    const auto *pgg0_536 = buffer.data(pgg0 + 536);
    const auto *pgg0_537 = buffer.data(pgg0 + 537);
    const auto *pgg0_538 = buffer.data(pgg0 + 538);
    const auto *pgg0_539 = buffer.data(pgg0 + 539);
    const auto *pgg0_555 = buffer.data(pgg0 + 555);
    const auto *pgg0_558 = buffer.data(pgg0 + 558);
    const auto *pgg0_561 = buffer.data(pgg0 + 561);
    const auto *pgg0_565 = buffer.data(pgg0 + 565);

    const auto *pgg1_513 = buffer.data(pgg1 + 513);
    const auto *pgg1_516 = buffer.data(pgg1 + 516);
    const auto *pgg1_520 = buffer.data(pgg1 + 520);
    const auto *pgg1_521 = buffer.data(pgg1 + 521);
    const auto *pgg1_522 = buffer.data(pgg1 + 522);
    const auto *pgg1_523 = buffer.data(pgg1 + 523);
    const auto *pgg1_524 = buffer.data(pgg1 + 524);
    const auto *pgg1_525 = buffer.data(pgg1 + 525);
    const auto *pgg1_527 = buffer.data(pgg1 + 527);
    const auto *pgg1_528 = buffer.data(pgg1 + 528);
    const auto *pgg1_530 = buffer.data(pgg1 + 530);
    const auto *pgg1_531 = buffer.data(pgg1 + 531);
    const auto *pgg1_532 = buffer.data(pgg1 + 532);
    const auto *pgg1_534 = buffer.data(pgg1 + 534);
    const auto *pgg1_535 = buffer.data(pgg1 + 535);
    const auto *pgg1_536 = buffer.data(pgg1 + 536);
    const auto *pgg1_537 = buffer.data(pgg1 + 537);
    const auto *pgg1_538 = buffer.data(pgg1 + 538);
    const auto *pgg1_539 = buffer.data(pgg1 + 539);
    const auto *pgg1_555 = buffer.data(pgg1 + 555);
    const auto *pgg1_558 = buffer.data(pgg1 + 558);
    const auto *pgg1_561 = buffer.data(pgg1 + 561);
    const auto *pgg1_565 = buffer.data(pgg1 + 565);

    const auto *pgh_713 = buffer.data(pgh + 713);
    const auto *pgh_714 = buffer.data(pgh + 714);
    const auto *pgh_716 = buffer.data(pgh + 716);
    const auto *pgh_717 = buffer.data(pgh + 717);
    const auto *pgh_719 = buffer.data(pgh + 719);
    const auto *pgh_720 = buffer.data(pgh + 720);
    const auto *pgh_723 = buffer.data(pgh + 723);
    const auto *pgh_724 = buffer.data(pgh + 724);
    const auto *pgh_729 = buffer.data(pgh + 729);
    const auto *pgh_730 = buffer.data(pgh + 730);
    const auto *pgh_731 = buffer.data(pgh + 731);
    const auto *pgh_732 = buffer.data(pgh + 732);
    const auto *pgh_733 = buffer.data(pgh + 733);
    const auto *pgh_734 = buffer.data(pgh + 734);
    const auto *pgh_735 = buffer.data(pgh + 735);
    const auto *pgh_737 = buffer.data(pgh + 737);
    const auto *pgh_738 = buffer.data(pgh + 738);
    const auto *pgh_740 = buffer.data(pgh + 740);
    const auto *pgh_741 = buffer.data(pgh + 741);
    const auto *pgh_742 = buffer.data(pgh + 742);
    const auto *pgh_744 = buffer.data(pgh + 744);
    const auto *pgh_745 = buffer.data(pgh + 745);
    const auto *pgh_746 = buffer.data(pgh + 746);
    const auto *pgh_747 = buffer.data(pgh + 747);
    const auto *pgh_749 = buffer.data(pgh + 749);
    const auto *pgh_750 = buffer.data(pgh + 750);
    const auto *pgh_751 = buffer.data(pgh + 751);
    const auto *pgh_752 = buffer.data(pgh + 752);
    const auto *pgh_753 = buffer.data(pgh + 753);
    const auto *pgh_754 = buffer.data(pgh + 754);
    const auto *pgh_755 = buffer.data(pgh + 755);
    const auto *pgh_756 = buffer.data(pgh + 756);
    const auto *pgh_758 = buffer.data(pgh + 758);
    const auto *pgh_761 = buffer.data(pgh + 761);
    const auto *pgh_765 = buffer.data(pgh + 765);
    const auto *pgh_771 = buffer.data(pgh + 771);
    const auto *pgh_772 = buffer.data(pgh + 772);
    const auto *pgh_773 = buffer.data(pgh + 773);
    const auto *pgh_774 = buffer.data(pgh + 774);
    const auto *pgh_775 = buffer.data(pgh + 775);
    const auto *pgh_776 = buffer.data(pgh + 776);
    const auto *pgh_777 = buffer.data(pgh + 777);
    const auto *pgh_779 = buffer.data(pgh + 779);
    const auto *pgh_780 = buffer.data(pgh + 780);
    const auto *pgh_782 = buffer.data(pgh + 782);
    const auto *pgh_783 = buffer.data(pgh + 783);
    const auto *pgh_786 = buffer.data(pgh + 786);
    const auto *pgh_787 = buffer.data(pgh + 787);
    const auto *pgh_792 = buffer.data(pgh + 792);
    const auto *pgh_793 = buffer.data(pgh + 793);
    const auto *pgh_794 = buffer.data(pgh + 794);
    const auto *pgh_795 = buffer.data(pgh + 795);
    const auto *pgh_796 = buffer.data(pgh + 796);
    const auto *pgh_797 = buffer.data(pgh + 797);

#pragma omp simd aligned(t_949, t_950, t_951, pa_z, pc_y, pc_z, sgi0_109, sgi0_111, sgh_81, \
                         sgh_83, sgi1_109, sgi1_111, pfh_461, pgh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = pa_z[k] * sgi0_109[k]
                   + f_1 * sgh_81[k]
                   - f_11 * pc_z[k] * sgi1_109[k];

        t_950[k] = f_12 * pfh_461[k]
                   + f_4 * pc_y[k] * pgh_713[k];

        t_951[k] = pa_z[k] * sgi0_111[k]
                   + f_16 * sgh_83[k]
                   - f_11 * pc_z[k] * sgi1_111[k];
    }

#pragma omp simd aligned(t_952, t_953, t_954, pb_y, pc_y, pfi0_616, pfi0_618, pfh_462, \
                         pfi1_616, pfi1_618, pgh_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_952[k] = pb_y[k] * pfi0_616[k]
                   - f_11 * pc_y[k] * pfi1_616[k];

        t_953[k] = f_0 * pfh_462[k]
                   + f_4 * pc_y[k] * pgh_714[k];

        t_954[k] = pb_y[k] * pfi0_618[k]
                   - f_11 * pc_y[k] * pfi1_618[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, pb_y, pc_x, pc_y, pfi0_621, pfh_464, pfh_507, \
                         pfi1_621, pgg0_513, pgg1_513, pgh_716, \
                         pgh_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = f_12 * pfh_507[k]
                   + f_9 * pgg0_513[k]
                   - f_10 * pgg1_513[k]
                   + f_4 * pc_x[k] * pgh_717[k];

        t_956[k] = f_0 * pfh_464[k]
                   + f_4 * pc_y[k] * pgh_716[k];

        t_957[k] = pb_y[k] * pfi0_621[k]
                   - f_11 * pc_y[k] * pfi1_621[k];
    }

#pragma omp simd aligned(t_958, t_959, t_960, pb_y, pc_x, pc_y, pfi0_623, pfh_466, pfh_467, \
                         pfh_510, pfi1_623, pgg0_516, pgg1_516, pgh_719, \
                         pgh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_958[k] = f_12 * pfh_510[k]
                   + f_7 * pgg0_516[k]
                   - f_8 * pgg1_516[k]
                   + f_4 * pc_x[k] * pgh_720[k];

        t_959[k] = pb_y[k] * pfi0_623[k]
                   + f_12 * pfh_466[k]
                   - f_11 * pc_y[k] * pfi1_623[k];

        t_960[k] = f_0 * pfh_467[k]
                   + f_4 * pc_y[k] * pgh_719[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pb_y, pc_x, pc_y, pfi0_625, pfi0_627, pfh_469, \
                         pfh_514, pfi1_625, pfi1_627, pgg0_520, pgg1_520, \
                         pgh_724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = pb_y[k] * pfi0_625[k]
                   - f_11 * pc_y[k] * pfi1_625[k];

        t_962[k] = f_12 * pfh_514[k]
                   + f_5 * pgg0_520[k]
                   - f_6 * pgg1_520[k]
                   + f_4 * pc_x[k] * pgh_724[k];

        t_963[k] = pb_y[k] * pfi0_627[k]
                   + f_13 * pfh_469[k]
                   - f_11 * pc_y[k] * pfi1_627[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pb_y, pc_x, pc_y, pfi0_628, pfi0_630, \
                         pfh_470, pfh_471, pfh_519, pfi1_628, pfi1_630, pgh_723, \
                         pgh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pb_y[k] * pfi0_628[k]
                   + f_12 * pfh_470[k]
                   - f_11 * pc_y[k] * pfi1_628[k];

        t_965[k] = f_0 * pfh_471[k]
                   + f_4 * pc_y[k] * pgh_723[k];

        t_966[k] = pb_y[k] * pfi0_630[k]
                   - f_11 * pc_y[k] * pfi1_630[k];

        t_967[k] = f_12 * pfh_519[k]
                   + f_4 * pc_x[k] * pgh_729[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, pc_x, pfh_520, pfh_521, pfh_522, \
                         pfh_523, pfh_524, pgh_730, pgh_731, pgh_732, pgh_733, \
                         pgh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_12 * pfh_520[k]
                   + f_4 * pc_x[k] * pgh_730[k];

        t_969[k] = f_12 * pfh_521[k]
                   + f_4 * pc_x[k] * pgh_731[k];

        t_970[k] = f_12 * pfh_522[k]
                   + f_4 * pc_x[k] * pgh_732[k];

        t_971[k] = f_12 * pfh_523[k]
                   + f_4 * pc_x[k] * pgh_733[k];

        t_972[k] = f_12 * pfh_524[k]
                   + f_4 * pc_x[k] * pgh_734[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, pc_y, pfh_477, pfh_478, pfh_479, pgg0_520, \
                         pgg0_521, pgg0_522, pgg1_520, pgg1_521, pgg1_522, pgh_729, pgh_730, \
                         pgh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_0 * pfh_477[k]
                   + f_2 * pgg0_520[k]
                   - f_3 * pgg1_520[k]
                   + f_4 * pc_y[k] * pgh_729[k];

        t_974[k] = f_0 * pfh_478[k]
                   + f_17 * pgg0_521[k]
                   - f_18 * pgg1_521[k]
                   + f_4 * pc_y[k] * pgh_730[k];

        t_975[k] = f_0 * pfh_479[k]
                   + f_9 * pgg0_522[k]
                   - f_10 * pgg1_522[k]
                   + f_4 * pc_y[k] * pgh_731[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_y, pfh_480, pfh_481, pfh_482, pgg0_523, \
                         pgg0_524, pgg1_523, pgg1_524, pgh_732, pgh_733, \
                         pgh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_0 * pfh_480[k]
                   + f_7 * pgg0_523[k]
                   - f_8 * pgg1_523[k]
                   + f_4 * pc_y[k] * pgh_732[k];

        t_977[k] = f_0 * pfh_481[k]
                   + f_5 * pgg0_524[k]
                   - f_6 * pgg1_524[k]
                   + f_4 * pc_y[k] * pgh_733[k];

        t_978[k] = f_0 * pfh_482[k]
                   + f_4 * pc_y[k] * pgh_734[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, pb_y, pc_x, pc_y, pfi0_643, pfh_525, pfi1_643, \
                         pgg0_525, pgg1_525, pgh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = pb_y[k] * pfi0_643[k]
                   - f_11 * pc_y[k] * pfi1_643[k];

        t_980[k] = f_12 * pfh_525[k]
                   + f_2 * pgg0_525[k]
                   - f_3 * pgg1_525[k]
                   + f_4 * pc_x[k] * pgh_735[k];

        t_981[k] = f_4 * pc_y[k] * pgh_735[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, pc_x, pc_y, pfh_527, pfh_528, pgg0_527, \
                         pgg0_528, pgg1_527, pgg1_528, pgh_737, \
                         pgh_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_12 * pfh_527[k]
                   + f_17 * pgg0_527[k]
                   - f_18 * pgg1_527[k]
                   + f_4 * pc_x[k] * pgh_737[k];

        t_983[k] = f_12 * pfh_528[k]
                   + f_9 * pgg0_528[k]
                   - f_10 * pgg1_528[k]
                   + f_4 * pc_x[k] * pgh_738[k];

        t_984[k] = f_4 * pc_y[k] * pgh_737[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, pc_x, pfh_530, pfh_531, pfh_532, pgg0_530, \
                         pgg0_531, pgg0_532, pgg1_530, pgg1_531, pgg1_532, pgh_740, pgh_741, \
                         pgh_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_12 * pfh_530[k]
                   + f_9 * pgg0_530[k]
                   - f_10 * pgg1_530[k]
                   + f_4 * pc_x[k] * pgh_740[k];

        t_986[k] = f_12 * pfh_531[k]
                   + f_7 * pgg0_531[k]
                   - f_8 * pgg1_531[k]
                   + f_4 * pc_x[k] * pgh_741[k];

        t_987[k] = f_12 * pfh_532[k]
                   + f_7 * pgg0_532[k]
                   - f_8 * pgg1_532[k]
                   + f_4 * pc_x[k] * pgh_742[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, pc_x, pc_y, pfh_534, pfh_535, pgg0_534, \
                         pgg0_535, pgg1_534, pgg1_535, pgh_740, pgh_744, \
                         pgh_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_4 * pc_y[k] * pgh_740[k];

        t_989[k] = f_12 * pfh_534[k]
                   + f_7 * pgg0_534[k]
                   - f_8 * pgg1_534[k]
                   + f_4 * pc_x[k] * pgh_744[k];

        t_990[k] = f_12 * pfh_535[k]
                   + f_5 * pgg0_535[k]
                   - f_6 * pgg1_535[k]
                   + f_4 * pc_x[k] * pgh_745[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, pc_x, pc_y, pfh_536, pfh_537, pgg0_536, \
                         pgg0_537, pgg1_536, pgg1_537, pgh_744, pgh_746, \
                         pgh_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_12 * pfh_536[k]
                   + f_5 * pgg0_536[k]
                   - f_6 * pgg1_536[k]
                   + f_4 * pc_x[k] * pgh_746[k];

        t_992[k] = f_12 * pfh_537[k]
                   + f_5 * pgg0_537[k]
                   - f_6 * pgg1_537[k]
                   + f_4 * pc_x[k] * pgh_747[k];

        t_993[k] = f_4 * pc_y[k] * pgh_744[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, t_997, pc_x, pfh_539, pfh_540, pfh_541, pfh_542, \
                         pgg0_539, pgg1_539, pgh_749, pgh_750, pgh_751, \
                         pgh_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = f_12 * pfh_539[k]
                   + f_5 * pgg0_539[k]
                   - f_6 * pgg1_539[k]
                   + f_4 * pc_x[k] * pgh_749[k];

        t_995[k] = f_12 * pfh_540[k]
                   + f_4 * pc_x[k] * pgh_750[k];

        t_996[k] = f_12 * pfh_541[k]
                   + f_4 * pc_x[k] * pgh_751[k];

        t_997[k] = f_12 * pfh_542[k]
                   + f_4 * pc_x[k] * pgh_752[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, pc_x, pc_y, pfh_543, pfh_544, pfh_545, \
                         pgg0_535, pgg1_535, pgh_750, pgh_753, pgh_754, \
                         pgh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_12 * pfh_543[k]
                   + f_4 * pc_x[k] * pgh_753[k];

        t_999[k] = f_12 * pfh_544[k]
                   + f_4 * pc_x[k] * pgh_754[k];

        t_1000[k] = f_12 * pfh_545[k]
                    + f_4 * pc_x[k] * pgh_755[k];

        t_1001[k] = f_2 * pgg0_535[k]
                    - f_3 * pgg1_535[k]
                    + f_4 * pc_y[k] * pgh_750[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, pgg0_536, pgg0_537, pgg0_538, pgg1_536, \
                         pgg1_537, pgg1_538, pgh_751, pgh_752, \
                         pgh_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_17 * pgg0_536[k]
                    - f_18 * pgg1_536[k]
                    + f_4 * pc_y[k] * pgh_751[k];

        t_1003[k] = f_9 * pgg0_537[k]
                    - f_10 * pgg1_537[k]
                    + f_4 * pc_y[k] * pgh_752[k];

        t_1004[k] = f_7 * pgg0_538[k]
                    - f_8 * pgg1_538[k]
                    + f_4 * pc_y[k] * pgh_753[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, pb_x, pc_x, pc_y, pdi0_503, pdi1_503, \
                         pfi0_727, pfi1_727, pgg0_539, pgg1_539, pgh_754, \
                         pgh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_5 * pgg0_539[k]
                    - f_6 * pgg1_539[k]
                    + f_4 * pc_y[k] * pgh_754[k];

        t_1006[k] = f_4 * pc_y[k] * pgh_755[k];

        t_1007[k] = f_14 * pdi0_503[k]
                    - f_15 * pdi1_503[k]
                    + pb_x[k] * pfi0_727[k]
                    - f_11 * pc_x[k] * pfi1_727[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pa_z, pc_y, pc_z, sgi0_168, sgi0_170, \
                         sgi0_171, sgh_126, sgi1_168, sgi1_170, sgi1_171, pfh_483, \
                         pgh_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = pa_z[k] * sgi0_168[k]
                    - f_11 * pc_z[k] * sgi1_168[k];

        t_1009[k] = f_13 * pfh_483[k]
                    + f_4 * pc_y[k] * pgh_756[k];

        t_1010[k] = pa_z[k] * sgi0_170[k]
                    + f_0 * sgh_126[k]
                    - f_11 * pc_z[k] * sgi1_170[k];

        t_1011[k] = pa_z[k] * sgi0_171[k]
                    - f_11 * pc_z[k] * sgi1_171[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pa_z, pc_y, pc_z, sgi0_173, sgi0_174, \
                         sgh_128, sgi1_173, sgi1_174, pfh_485, \
                         pgh_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_13 * pfh_485[k]
                    + f_4 * pc_y[k] * pgh_758[k];

        t_1013[k] = pa_z[k] * sgi0_173[k]
                    + f_12 * sgh_128[k]
                    - f_11 * pc_z[k] * sgi1_173[k];

        t_1014[k] = pa_z[k] * sgi0_174[k]
                    - f_11 * pc_z[k] * sgi1_174[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, pa_z, pc_y, pc_z, sgi0_175, sgi0_177, \
                         sgh_129, sgh_131, sgi1_175, sgi1_177, pfh_488, \
                         pgh_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = pa_z[k] * sgi0_175[k]
                    + f_0 * sgh_129[k]
                    - f_11 * pc_z[k] * sgi1_175[k];

        t_1016[k] = f_13 * pfh_488[k]
                    + f_4 * pc_y[k] * pgh_761[k];

        t_1017[k] = pa_z[k] * sgi0_177[k]
                    + f_13 * sgh_131[k]
                    - f_11 * pc_z[k] * sgi1_177[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, pa_z, pc_z, sgi0_178, sgi0_179, sgi0_180, \
                         sgh_132, sgh_133, sgi1_178, sgi1_179, \
                         sgi1_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = pa_z[k] * sgi0_178[k]
                    - f_11 * pc_z[k] * sgi1_178[k];

        t_1019[k] = pa_z[k] * sgi0_179[k]
                    + f_0 * sgh_132[k]
                    - f_11 * pc_z[k] * sgi1_179[k];

        t_1020[k] = pa_z[k] * sgi0_180[k]
                    + f_12 * sgh_133[k]
                    - f_11 * pc_z[k] * sgi1_180[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, pa_z, pc_x, pc_y, pc_z, sgi0_182, sgh_135, \
                         sgi1_182, pfh_492, pfh_561, pgh_765, pgh_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_13 * pfh_492[k]
                    + f_4 * pc_y[k] * pgh_765[k];

        t_1022[k] = pa_z[k] * sgi0_182[k]
                    + f_1 * sgh_135[k]
                    - f_11 * pc_z[k] * sgi1_182[k];

        t_1023[k] = f_0 * pfh_561[k]
                    + f_4 * pc_x[k] * pgh_771[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, t_1027, t_1028, pc_x, pfh_562, pfh_563, \
                         pfh_564, pfh_565, pfh_566, pgh_772, pgh_773, pgh_774, pgh_775, \
                         pgh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_0 * pfh_562[k]
                    + f_4 * pc_x[k] * pgh_772[k];

        t_1025[k] = f_0 * pfh_563[k]
                    + f_4 * pc_x[k] * pgh_773[k];

        t_1026[k] = f_0 * pfh_564[k]
                    + f_4 * pc_x[k] * pgh_774[k];

        t_1027[k] = f_0 * pfh_565[k]
                    + f_4 * pc_x[k] * pgh_775[k];

        t_1028[k] = f_0 * pfh_566[k]
                    + f_4 * pc_x[k] * pgh_776[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, pa_z, pb_x, pc_x, pc_z, sgi0_189, \
                         sgi1_189, pfi0_750, pfi0_751, pfi0_752, pfi1_750, pfi1_751, \
                         pfi1_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = pa_z[k] * sgi0_189[k]
                    - f_11 * pc_z[k] * sgi1_189[k];

        t_1030[k] = pb_x[k] * pfi0_750[k]
                    - f_11 * pc_x[k] * pfi1_750[k];

        t_1031[k] = pb_x[k] * pfi0_751[k]
                    - f_11 * pc_x[k] * pfi1_751[k];

        t_1032[k] = pb_x[k] * pfi0_752[k]
                    - f_11 * pc_x[k] * pfi1_752[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, pb_x, pc_x, pc_y, pfi0_753, pfi0_755, \
                         pfh_503, pfi1_753, pfi1_755, pgh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = pb_x[k] * pfi0_753[k]
                    - f_11 * pc_x[k] * pfi1_753[k];

        t_1034[k] = f_13 * pfh_503[k]
                    + f_4 * pc_y[k] * pgh_776[k];

        t_1035[k] = pb_x[k] * pfi0_755[k]
                    - f_11 * pc_x[k] * pfi1_755[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, pb_x, pc_x, pc_y, pfi0_758, pfh_504, pfh_567, \
                         pfh_569, pfi1_758, pgg0_555, pgg1_555, \
                         pgh_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = f_0 * pfh_567[k]
                    + f_2 * pgg0_555[k]
                    - f_3 * pgg1_555[k]
                    + f_4 * pc_x[k] * pgh_777[k];

        t_1037[k] = f_12 * pfh_504[k]
                    + f_4 * pc_y[k] * pgh_777[k];

        t_1038[k] = pb_x[k] * pfi0_758[k]
                    + f_21 * pfh_569[k]
                    - f_11 * pc_x[k] * pfi1_758[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pb_x, pc_x, pc_y, pfi0_761, pfh_506, pfh_570, \
                         pfh_572, pfi1_761, pgg0_558, pgg1_558, pgh_779, \
                         pgh_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_0 * pfh_570[k]
                    + f_9 * pgg0_558[k]
                    - f_10 * pgg1_558[k]
                    + f_4 * pc_x[k] * pgh_780[k];

        t_1040[k] = f_12 * pfh_506[k]
                    + f_4 * pc_y[k] * pgh_779[k];

        t_1041[k] = pb_x[k] * pfi0_761[k]
                    + f_1 * pfh_572[k]
                    - f_11 * pc_x[k] * pfi1_761[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, pb_x, pc_x, pc_y, pfi0_763, pfh_509, pfh_573, \
                         pfh_574, pfi1_763, pgg0_561, pgg1_561, pgh_782, \
                         pgh_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_0 * pfh_573[k]
                    + f_7 * pgg0_561[k]
                    - f_8 * pgg1_561[k]
                    + f_4 * pc_x[k] * pgh_783[k];

        t_1043[k] = pb_x[k] * pfi0_763[k]
                    + f_13 * pfh_574[k]
                    - f_11 * pc_x[k] * pfi1_763[k];

        t_1044[k] = f_12 * pfh_509[k]
                    + f_4 * pc_y[k] * pgh_782[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, pb_x, pc_x, pfi0_765, pfi0_767, pfh_576, \
                         pfh_577, pfh_578, pfi1_765, pfi1_767, pgg0_565, pgg1_565, \
                         pgh_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = pb_x[k] * pfi0_765[k]
                    + f_13 * pfh_576[k]
                    - f_11 * pc_x[k] * pfi1_765[k];

        t_1046[k] = f_0 * pfh_577[k]
                    + f_5 * pgg0_565[k]
                    - f_6 * pgg1_565[k]
                    + f_4 * pc_x[k] * pgh_787[k];

        t_1047[k] = pb_x[k] * pfi0_767[k]
                    + f_12 * pfh_578[k]
                    - f_11 * pc_x[k] * pfi1_767[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, pb_x, pc_x, pc_y, pfi0_768, pfi0_770, \
                         pfh_513, pfh_579, pfh_581, pfi1_768, pfi1_770, \
                         pgh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = pb_x[k] * pfi0_768[k]
                    + f_12 * pfh_579[k]
                    - f_11 * pc_x[k] * pfi1_768[k];

        t_1049[k] = f_12 * pfh_513[k]
                    + f_4 * pc_y[k] * pgh_786[k];

        t_1050[k] = pb_x[k] * pfi0_770[k]
                    + f_12 * pfh_581[k]
                    - f_11 * pc_x[k] * pfi1_770[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, t_1054, t_1055, pc_x, pfh_582, pfh_583, \
                         pfh_584, pfh_585, pfh_586, pgh_792, pgh_793, pgh_794, pgh_795, \
                         pgh_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_0 * pfh_582[k]
                    + f_4 * pc_x[k] * pgh_792[k];

        t_1052[k] = f_0 * pfh_583[k]
                    + f_4 * pc_x[k] * pgh_793[k];

        t_1053[k] = f_0 * pfh_584[k]
                    + f_4 * pc_x[k] * pgh_794[k];

        t_1054[k] = f_0 * pfh_585[k]
                    + f_4 * pc_x[k] * pgh_795[k];

        t_1055[k] = f_0 * pfh_586[k]
                    + f_4 * pc_x[k] * pgh_796[k];
    }

#pragma omp simd aligned(t_1056, t_1057, t_1058, t_1059, pb_x, pc_x, pfi0_777, pfi0_778, \
                         pfi0_779, pfh_587, pfi1_777, pfi1_778, pfi1_779, \
                         pgh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1056[k] = f_0 * pfh_587[k]
                    + f_4 * pc_x[k] * pgh_797[k];

        t_1057[k] = pb_x[k] * pfi0_777[k]
                    - f_11 * pc_x[k] * pfi1_777[k];

        t_1058[k] = pb_x[k] * pfi0_778[k]
                    - f_11 * pc_x[k] * pfi1_778[k];

        t_1059[k] = pb_x[k] * pfi0_779[k]
                    - f_11 * pc_x[k] * pfi1_779[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sgi0, const size_t sgh,
                                                          const size_t sgi1, const size_t pdi0,
                                                          const size_t pdi1, const size_t pfi0,
                                                          const size_t pfh, const size_t pfi1,
                                                          const size_t pgg0, const size_t pgg1,
                                                          const size_t pgh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 1.0 / p;
    const auto f_20 = gamma / (p * q);
    const auto f_21 = 2.5 / q;

    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);
    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgi0_280 = buffer.data(sgi0 + 280);
    const auto *sgi0_283 = buffer.data(sgi0 + 283);
    const auto *sgi0_286 = buffer.data(sgi0 + 286);
    const auto *sgi0_290 = buffer.data(sgi0 + 290);
    const auto *sgi0_301 = buffer.data(sgi0 + 301);
    const auto *sgi0_302 = buffer.data(sgi0 + 302);
    const auto *sgi0_303 = buffer.data(sgi0 + 303);
    const auto *sgi0_304 = buffer.data(sgi0 + 304);
    const auto *sgi0_305 = buffer.data(sgi0 + 305);
    const auto *sgi0_307 = buffer.data(sgi0 + 307);

    const auto *sgh_225 = buffer.data(sgh + 225);
    const auto *sgh_226 = buffer.data(sgh + 226);
    const auto *sgh_227 = buffer.data(sgh + 227);
    const auto *sgh_228 = buffer.data(sgh + 228);
    const auto *sgh_230 = buffer.data(sgh + 230);

    const auto *sgi1_280 = buffer.data(sgi1 + 280);
    const auto *sgi1_283 = buffer.data(sgi1 + 283);
    const auto *sgi1_286 = buffer.data(sgi1 + 286);
    const auto *sgi1_290 = buffer.data(sgi1 + 290);
    const auto *sgi1_301 = buffer.data(sgi1 + 301);
    const auto *sgi1_302 = buffer.data(sgi1 + 302);
    const auto *sgi1_303 = buffer.data(sgi1 + 303);
    const auto *sgi1_304 = buffer.data(sgi1 + 304);
    const auto *sgi1_305 = buffer.data(sgi1 + 305);
    const auto *sgi1_307 = buffer.data(sgi1 + 307);

    const auto *pdi0_475 = buffer.data(pdi0 + 475);

    const auto *pdi1_475 = buffer.data(pdi1 + 475);

    const auto *pfi0_700 = buffer.data(pfi0 + 700);
    const auto *pfi0_702 = buffer.data(pfi0 + 702);
    const auto *pfi0_705 = buffer.data(pfi0 + 705);
    const auto *pfi0_709 = buffer.data(pfi0 + 709);
    const auto *pfi0_714 = buffer.data(pfi0 + 714);
    const auto *pfi0_780 = buffer.data(pfi0 + 780);
    const auto *pfi0_781 = buffer.data(pfi0 + 781);
    const auto *pfi0_783 = buffer.data(pfi0 + 783);
    const auto *pfi0_787 = buffer.data(pfi0 + 787);
    const auto *pfi0_790 = buffer.data(pfi0 + 790);
    const auto *pfi0_791 = buffer.data(pfi0 + 791);
    const auto *pfi0_794 = buffer.data(pfi0 + 794);
    const auto *pfi0_795 = buffer.data(pfi0 + 795);
    const auto *pfi0_796 = buffer.data(pfi0 + 796);
    const auto *pfi0_805 = buffer.data(pfi0 + 805);
    const auto *pfi0_806 = buffer.data(pfi0 + 806);
    const auto *pfi0_807 = buffer.data(pfi0 + 807);
    const auto *pfi0_808 = buffer.data(pfi0 + 808);
    const auto *pfi0_809 = buffer.data(pfi0 + 809);
    const auto *pfi0_811 = buffer.data(pfi0 + 811);
    const auto *pfi0_812 = buffer.data(pfi0 + 812);
    const auto *pfi0_814 = buffer.data(pfi0 + 814);
    const auto *pfi0_815 = buffer.data(pfi0 + 815);
    const auto *pfi0_817 = buffer.data(pfi0 + 817);
    const auto *pfi0_818 = buffer.data(pfi0 + 818);
    const auto *pfi0_819 = buffer.data(pfi0 + 819);
    const auto *pfi0_821 = buffer.data(pfi0 + 821);
    const auto *pfi0_822 = buffer.data(pfi0 + 822);
    const auto *pfi0_823 = buffer.data(pfi0 + 823);
    const auto *pfi0_824 = buffer.data(pfi0 + 824);
    const auto *pfi0_826 = buffer.data(pfi0 + 826);
    const auto *pfi0_833 = buffer.data(pfi0 + 833);
    const auto *pfi0_834 = buffer.data(pfi0 + 834);
    const auto *pfi0_835 = buffer.data(pfi0 + 835);
    const auto *pfi0_836 = buffer.data(pfi0 + 836);
    const auto *pfi0_837 = buffer.data(pfi0 + 837);
    const auto *pfi0_839 = buffer.data(pfi0 + 839);

    const auto *pfh_524 = buffer.data(pfh + 524);
    const auto *pfh_525 = buffer.data(pfh + 525);
    const auto *pfh_527 = buffer.data(pfh + 527);
    const auto *pfh_530 = buffer.data(pfh + 530);
    const auto *pfh_534 = buffer.data(pfh + 534);
    const auto *pfh_545 = buffer.data(pfh + 545);
    const auto *pfh_546 = buffer.data(pfh + 546);
    const auto *pfh_548 = buffer.data(pfh + 548);
    const auto *pfh_551 = buffer.data(pfh + 551);
    const auto *pfh_555 = buffer.data(pfh + 555);
    const auto *pfh_566 = buffer.data(pfh + 566);
    const auto *pfh_567 = buffer.data(pfh + 567);
    const auto *pfh_569 = buffer.data(pfh + 569);
    const auto *pfh_572 = buffer.data(pfh + 572);
    const auto *pfh_576 = buffer.data(pfh + 576);
    const auto *pfh_582 = buffer.data(pfh + 582);
    const auto *pfh_583 = buffer.data(pfh + 583);
    const auto *pfh_584 = buffer.data(pfh + 584);
    const auto *pfh_585 = buffer.data(pfh + 585);
    const auto *pfh_586 = buffer.data(pfh + 586);
    const auto *pfh_587 = buffer.data(pfh + 587);
    const auto *pfh_588 = buffer.data(pfh + 588);
    const auto *pfh_590 = buffer.data(pfh + 590);
    const auto *pfh_591 = buffer.data(pfh + 591);
    const auto *pfh_594 = buffer.data(pfh + 594);
    const auto *pfh_595 = buffer.data(pfh + 595);
    const auto *pfh_598 = buffer.data(pfh + 598);
    const auto *pfh_599 = buffer.data(pfh + 599);
    const auto *pfh_600 = buffer.data(pfh + 600);
    const auto *pfh_603 = buffer.data(pfh + 603);
    const auto *pfh_604 = buffer.data(pfh + 604);
    const auto *pfh_605 = buffer.data(pfh + 605);
    const auto *pfh_606 = buffer.data(pfh + 606);
    const auto *pfh_607 = buffer.data(pfh + 607);
    const auto *pfh_608 = buffer.data(pfh + 608);
    const auto *pfh_609 = buffer.data(pfh + 609);
    const auto *pfh_611 = buffer.data(pfh + 611);
    const auto *pfh_612 = buffer.data(pfh + 612);
    const auto *pfh_614 = buffer.data(pfh + 614);
    const auto *pfh_615 = buffer.data(pfh + 615);
    const auto *pfh_616 = buffer.data(pfh + 616);
    const auto *pfh_618 = buffer.data(pfh + 618);
    const auto *pfh_619 = buffer.data(pfh + 619);
    const auto *pfh_620 = buffer.data(pfh + 620);
    const auto *pfh_621 = buffer.data(pfh + 621);
    const auto *pfh_623 = buffer.data(pfh + 623);
    const auto *pfh_624 = buffer.data(pfh + 624);
    const auto *pfh_625 = buffer.data(pfh + 625);
    const auto *pfh_626 = buffer.data(pfh + 626);
    const auto *pfh_627 = buffer.data(pfh + 627);
    const auto *pfh_628 = buffer.data(pfh + 628);
    const auto *pfh_629 = buffer.data(pfh + 629);

    const auto *pfi1_700 = buffer.data(pfi1 + 700);
    const auto *pfi1_702 = buffer.data(pfi1 + 702);
    const auto *pfi1_705 = buffer.data(pfi1 + 705);
    const auto *pfi1_709 = buffer.data(pfi1 + 709);
    const auto *pfi1_714 = buffer.data(pfi1 + 714);
    const auto *pfi1_780 = buffer.data(pfi1 + 780);
    const auto *pfi1_781 = buffer.data(pfi1 + 781);
    const auto *pfi1_783 = buffer.data(pfi1 + 783);
    const auto *pfi1_787 = buffer.data(pfi1 + 787);
    const auto *pfi1_790 = buffer.data(pfi1 + 790);
    const auto *pfi1_791 = buffer.data(pfi1 + 791);
    const auto *pfi1_794 = buffer.data(pfi1 + 794);
    const auto *pfi1_795 = buffer.data(pfi1 + 795);
    const auto *pfi1_796 = buffer.data(pfi1 + 796);
    const auto *pfi1_805 = buffer.data(pfi1 + 805);
    const auto *pfi1_806 = buffer.data(pfi1 + 806);
    const auto *pfi1_807 = buffer.data(pfi1 + 807);
    const auto *pfi1_808 = buffer.data(pfi1 + 808);
    const auto *pfi1_809 = buffer.data(pfi1 + 809);
    const auto *pfi1_811 = buffer.data(pfi1 + 811);
    const auto *pfi1_812 = buffer.data(pfi1 + 812);
    const auto *pfi1_814 = buffer.data(pfi1 + 814);
    const auto *pfi1_815 = buffer.data(pfi1 + 815);
    const auto *pfi1_817 = buffer.data(pfi1 + 817);
    const auto *pfi1_818 = buffer.data(pfi1 + 818);
    const auto *pfi1_819 = buffer.data(pfi1 + 819);
    const auto *pfi1_821 = buffer.data(pfi1 + 821);
    const auto *pfi1_822 = buffer.data(pfi1 + 822);
    const auto *pfi1_823 = buffer.data(pfi1 + 823);
    const auto *pfi1_824 = buffer.data(pfi1 + 824);
    const auto *pfi1_826 = buffer.data(pfi1 + 826);
    const auto *pfi1_833 = buffer.data(pfi1 + 833);
    const auto *pfi1_834 = buffer.data(pfi1 + 834);
    const auto *pfi1_835 = buffer.data(pfi1 + 835);
    const auto *pfi1_836 = buffer.data(pfi1 + 836);
    const auto *pfi1_837 = buffer.data(pfi1 + 837);
    const auto *pfi1_839 = buffer.data(pfi1 + 839);

    const auto *pgg0_602 = buffer.data(pgg0 + 602);
    const auto *pgg0_605 = buffer.data(pgg0 + 605);
    const auto *pgg0_607 = buffer.data(pgg0 + 607);
    const auto *pgg0_609 = buffer.data(pgg0 + 609);
    const auto *pgg0_611 = buffer.data(pgg0 + 611);
    const auto *pgg0_612 = buffer.data(pgg0 + 612);
    const auto *pgg0_614 = buffer.data(pgg0 + 614);
    const auto *pgg0_615 = buffer.data(pgg0 + 615);
    const auto *pgg0_617 = buffer.data(pgg0 + 617);
    const auto *pgg0_618 = buffer.data(pgg0 + 618);
    const auto *pgg0_620 = buffer.data(pgg0 + 620);
    const auto *pgg0_621 = buffer.data(pgg0 + 621);
    const auto *pgg0_622 = buffer.data(pgg0 + 622);
    const auto *pgg0_624 = buffer.data(pgg0 + 624);
    const auto *pgg0_625 = buffer.data(pgg0 + 625);
    const auto *pgg0_626 = buffer.data(pgg0 + 626);
    const auto *pgg0_627 = buffer.data(pgg0 + 627);
    const auto *pgg0_628 = buffer.data(pgg0 + 628);
    const auto *pgg0_629 = buffer.data(pgg0 + 629);
    const auto *pgg0_630 = buffer.data(pgg0 + 630);
    const auto *pgg0_632 = buffer.data(pgg0 + 632);
    const auto *pgg0_633 = buffer.data(pgg0 + 633);
    const auto *pgg0_635 = buffer.data(pgg0 + 635);

    const auto *pgg1_602 = buffer.data(pgg1 + 602);
    const auto *pgg1_605 = buffer.data(pgg1 + 605);
    const auto *pgg1_607 = buffer.data(pgg1 + 607);
    const auto *pgg1_609 = buffer.data(pgg1 + 609);
    const auto *pgg1_611 = buffer.data(pgg1 + 611);
    const auto *pgg1_612 = buffer.data(pgg1 + 612);
    const auto *pgg1_614 = buffer.data(pgg1 + 614);
    const auto *pgg1_615 = buffer.data(pgg1 + 615);
    const auto *pgg1_617 = buffer.data(pgg1 + 617);
    const auto *pgg1_618 = buffer.data(pgg1 + 618);
    const auto *pgg1_620 = buffer.data(pgg1 + 620);
    const auto *pgg1_621 = buffer.data(pgg1 + 621);
    const auto *pgg1_622 = buffer.data(pgg1 + 622);
    const auto *pgg1_624 = buffer.data(pgg1 + 624);
    const auto *pgg1_625 = buffer.data(pgg1 + 625);
    const auto *pgg1_626 = buffer.data(pgg1 + 626);
    const auto *pgg1_627 = buffer.data(pgg1 + 627);
    const auto *pgg1_628 = buffer.data(pgg1 + 628);
    const auto *pgg1_629 = buffer.data(pgg1 + 629);
    const auto *pgg1_630 = buffer.data(pgg1 + 630);
    const auto *pgg1_632 = buffer.data(pgg1 + 632);
    const auto *pgg1_633 = buffer.data(pgg1 + 633);
    const auto *pgg1_635 = buffer.data(pgg1 + 635);

    const auto *pgh_797 = buffer.data(pgh + 797);
    const auto *pgh_798 = buffer.data(pgh + 798);
    const auto *pgh_800 = buffer.data(pgh + 800);
    const auto *pgh_803 = buffer.data(pgh + 803);
    const auto *pgh_807 = buffer.data(pgh + 807);
    const auto *pgh_813 = buffer.data(pgh + 813);
    const auto *pgh_814 = buffer.data(pgh + 814);
    const auto *pgh_815 = buffer.data(pgh + 815);
    const auto *pgh_816 = buffer.data(pgh + 816);
    const auto *pgh_817 = buffer.data(pgh + 817);
    const auto *pgh_818 = buffer.data(pgh + 818);
    const auto *pgh_819 = buffer.data(pgh + 819);
    const auto *pgh_821 = buffer.data(pgh + 821);
    const auto *pgh_824 = buffer.data(pgh + 824);
    const auto *pgh_828 = buffer.data(pgh + 828);
    const auto *pgh_834 = buffer.data(pgh + 834);
    const auto *pgh_835 = buffer.data(pgh + 835);
    const auto *pgh_836 = buffer.data(pgh + 836);
    const auto *pgh_837 = buffer.data(pgh + 837);
    const auto *pgh_838 = buffer.data(pgh + 838);
    const auto *pgh_839 = buffer.data(pgh + 839);
    const auto *pgh_840 = buffer.data(pgh + 840);
    const auto *pgh_842 = buffer.data(pgh + 842);
    const auto *pgh_845 = buffer.data(pgh + 845);
    const auto *pgh_847 = buffer.data(pgh + 847);
    const auto *pgh_849 = buffer.data(pgh + 849);
    const auto *pgh_851 = buffer.data(pgh + 851);
    const auto *pgh_852 = buffer.data(pgh + 852);
    const auto *pgh_854 = buffer.data(pgh + 854);
    const auto *pgh_855 = buffer.data(pgh + 855);
    const auto *pgh_856 = buffer.data(pgh + 856);
    const auto *pgh_857 = buffer.data(pgh + 857);
    const auto *pgh_858 = buffer.data(pgh + 858);
    const auto *pgh_859 = buffer.data(pgh + 859);
    const auto *pgh_860 = buffer.data(pgh + 860);
    const auto *pgh_861 = buffer.data(pgh + 861);
    const auto *pgh_863 = buffer.data(pgh + 863);
    const auto *pgh_864 = buffer.data(pgh + 864);
    const auto *pgh_866 = buffer.data(pgh + 866);
    const auto *pgh_867 = buffer.data(pgh + 867);
    const auto *pgh_868 = buffer.data(pgh + 868);
    const auto *pgh_870 = buffer.data(pgh + 870);
    const auto *pgh_871 = buffer.data(pgh + 871);
    const auto *pgh_872 = buffer.data(pgh + 872);
    const auto *pgh_873 = buffer.data(pgh + 873);
    const auto *pgh_875 = buffer.data(pgh + 875);
    const auto *pgh_876 = buffer.data(pgh + 876);
    const auto *pgh_877 = buffer.data(pgh + 877);
    const auto *pgh_878 = buffer.data(pgh + 878);
    const auto *pgh_879 = buffer.data(pgh + 879);
    const auto *pgh_880 = buffer.data(pgh + 880);
    const auto *pgh_881 = buffer.data(pgh + 881);
    const auto *pgh_882 = buffer.data(pgh + 882);
    const auto *pgh_884 = buffer.data(pgh + 884);
    const auto *pgh_885 = buffer.data(pgh + 885);
    const auto *pgh_887 = buffer.data(pgh + 887);

#pragma omp simd aligned(t_1060, t_1061, t_1062, t_1063, pb_x, pc_x, pc_y, pfi0_780, pfi0_781, \
                         pfi0_783, pfh_524, pfi1_780, pfi1_781, pfi1_783, \
                         pgh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1060[k] = pb_x[k] * pfi0_780[k]
                    - f_11 * pc_x[k] * pfi1_780[k];

        t_1061[k] = pb_x[k] * pfi0_781[k]
                    - f_11 * pc_x[k] * pfi1_781[k];

        t_1062[k] = f_12 * pfh_524[k]
                    + f_4 * pc_y[k] * pgh_797[k];

        t_1063[k] = pb_x[k] * pfi0_783[k]
                    - f_11 * pc_x[k] * pfi1_783[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, pb_y, pc_y, pfi0_700, pfi0_702, pfh_525, \
                         pfi1_700, pfi1_702, pgh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = pb_y[k] * pfi0_700[k]
                    - f_11 * pc_y[k] * pfi1_700[k];

        t_1065[k] = f_0 * pfh_525[k]
                    + f_4 * pc_y[k] * pgh_798[k];

        t_1066[k] = pb_y[k] * pfi0_702[k]
                    - f_11 * pc_y[k] * pfi1_702[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, pb_x, pb_y, pc_x, pc_y, pfi0_705, pfi0_787, \
                         pfh_527, pfh_591, pfi1_705, pfi1_787, \
                         pgh_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = pb_x[k] * pfi0_787[k]
                    + f_1 * pfh_591[k]
                    - f_11 * pc_x[k] * pfi1_787[k];

        t_1068[k] = f_0 * pfh_527[k]
                    + f_4 * pc_y[k] * pgh_800[k];

        t_1069[k] = pb_y[k] * pfi0_705[k]
                    - f_11 * pc_y[k] * pfi1_705[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pb_x, pc_x, pc_y, pfi0_790, pfi0_791, \
                         pfh_530, pfh_594, pfh_595, pfi1_790, pfi1_791, \
                         pgh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = pb_x[k] * pfi0_790[k]
                    + f_13 * pfh_594[k]
                    - f_11 * pc_x[k] * pfi1_790[k];

        t_1071[k] = pb_x[k] * pfi0_791[k]
                    + f_13 * pfh_595[k]
                    - f_11 * pc_x[k] * pfi1_791[k];

        t_1072[k] = f_0 * pfh_530[k]
                    + f_4 * pc_y[k] * pgh_803[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pb_x, pb_y, pc_x, pc_y, pfi0_709, pfi0_794, \
                         pfi0_795, pfh_598, pfh_599, pfi1_709, pfi1_794, \
                         pfi1_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = pb_y[k] * pfi0_709[k]
                    - f_11 * pc_y[k] * pfi1_709[k];

        t_1074[k] = pb_x[k] * pfi0_794[k]
                    + f_12 * pfh_598[k]
                    - f_11 * pc_x[k] * pfi1_794[k];

        t_1075[k] = pb_x[k] * pfi0_795[k]
                    + f_12 * pfh_599[k]
                    - f_11 * pc_x[k] * pfi1_795[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pb_x, pb_y, pc_x, pc_y, pfi0_714, pfi0_796, \
                         pfh_534, pfh_600, pfi1_714, pfi1_796, \
                         pgh_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = pb_x[k] * pfi0_796[k]
                    + f_12 * pfh_600[k]
                    - f_11 * pc_x[k] * pfi1_796[k];

        t_1077[k] = f_0 * pfh_534[k]
                    + f_4 * pc_y[k] * pgh_807[k];

        t_1078[k] = pb_y[k] * pfi0_714[k]
                    - f_11 * pc_y[k] * pfi1_714[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, t_1083, pc_x, pfh_603, pfh_604, \
                         pfh_605, pfh_606, pfh_607, pgh_813, pgh_814, pgh_815, pgh_816, \
                         pgh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_0 * pfh_603[k]
                    + f_4 * pc_x[k] * pgh_813[k];

        t_1080[k] = f_0 * pfh_604[k]
                    + f_4 * pc_x[k] * pgh_814[k];

        t_1081[k] = f_0 * pfh_605[k]
                    + f_4 * pc_x[k] * pgh_815[k];

        t_1082[k] = f_0 * pfh_606[k]
                    + f_4 * pc_x[k] * pgh_816[k];

        t_1083[k] = f_0 * pfh_607[k]
                    + f_4 * pc_x[k] * pgh_817[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, t_1087, pb_x, pc_x, pfi0_805, pfi0_806, \
                         pfi0_807, pfh_608, pfi1_805, pfi1_806, pfi1_807, \
                         pgh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_0 * pfh_608[k]
                    + f_4 * pc_x[k] * pgh_818[k];

        t_1085[k] = pb_x[k] * pfi0_805[k]
                    - f_11 * pc_x[k] * pfi1_805[k];

        t_1086[k] = pb_x[k] * pfi0_806[k]
                    - f_11 * pc_x[k] * pfi1_806[k];

        t_1087[k] = pb_x[k] * pfi0_807[k]
                    - f_11 * pc_x[k] * pfi1_807[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, t_1091, pb_x, pc_x, pc_y, pfi0_808, pfi0_809, \
                         pfi0_811, pfh_545, pfi1_808, pfi1_809, pfi1_811, \
                         pgh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = pb_x[k] * pfi0_808[k]
                    - f_11 * pc_x[k] * pfi1_808[k];

        t_1089[k] = pb_x[k] * pfi0_809[k]
                    - f_11 * pc_x[k] * pfi1_809[k];

        t_1090[k] = f_0 * pfh_545[k]
                    + f_4 * pc_y[k] * pgh_818[k];

        t_1091[k] = pb_x[k] * pfi0_811[k]
                    - f_11 * pc_x[k] * pfi1_811[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pb_x, pc_x, pc_y, pfi0_812, pfi0_814, \
                         pfh_609, pfh_611, pfi1_812, pfi1_814, \
                         pgh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = pb_x[k] * pfi0_812[k]
                    + f_16 * pfh_609[k]
                    - f_11 * pc_x[k] * pfi1_812[k];

        t_1093[k] = f_4 * pc_y[k] * pgh_819[k];

        t_1094[k] = pb_x[k] * pfi0_814[k]
                    + f_21 * pfh_611[k]
                    - f_11 * pc_x[k] * pfi1_814[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pb_x, pc_x, pc_y, pfi0_815, pfi0_817, \
                         pfh_612, pfh_614, pfi1_815, pfi1_817, \
                         pgh_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = pb_x[k] * pfi0_815[k]
                    + f_1 * pfh_612[k]
                    - f_11 * pc_x[k] * pfi1_815[k];

        t_1096[k] = f_4 * pc_y[k] * pgh_821[k];

        t_1097[k] = pb_x[k] * pfi0_817[k]
                    + f_1 * pfh_614[k]
                    - f_11 * pc_x[k] * pfi1_817[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, pb_x, pc_x, pc_y, pfi0_818, pfi0_819, \
                         pfh_615, pfh_616, pfi1_818, pfi1_819, \
                         pgh_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = pb_x[k] * pfi0_818[k]
                    + f_13 * pfh_615[k]
                    - f_11 * pc_x[k] * pfi1_818[k];

        t_1099[k] = pb_x[k] * pfi0_819[k]
                    + f_13 * pfh_616[k]
                    - f_11 * pc_x[k] * pfi1_819[k];

        t_1100[k] = f_4 * pc_y[k] * pgh_824[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, pb_x, pc_x, pfi0_821, pfi0_822, pfi0_823, \
                         pfh_618, pfh_619, pfh_620, pfi1_821, pfi1_822, \
                         pfi1_823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = pb_x[k] * pfi0_821[k]
                    + f_13 * pfh_618[k]
                    - f_11 * pc_x[k] * pfi1_821[k];

        t_1102[k] = pb_x[k] * pfi0_822[k]
                    + f_12 * pfh_619[k]
                    - f_11 * pc_x[k] * pfi1_822[k];

        t_1103[k] = pb_x[k] * pfi0_823[k]
                    + f_12 * pfh_620[k]
                    - f_11 * pc_x[k] * pfi1_823[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pb_x, pc_x, pc_y, pfi0_824, pfi0_826, \
                         pfh_621, pfh_623, pfh_624, pfi1_824, pfi1_826, pgh_828, \
                         pgh_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = pb_x[k] * pfi0_824[k]
                    + f_12 * pfh_621[k]
                    - f_11 * pc_x[k] * pfi1_824[k];

        t_1105[k] = f_4 * pc_y[k] * pgh_828[k];

        t_1106[k] = pb_x[k] * pfi0_826[k]
                    + f_12 * pfh_623[k]
                    - f_11 * pc_x[k] * pfi1_826[k];

        t_1107[k] = f_0 * pfh_624[k]
                    + f_4 * pc_x[k] * pgh_834[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, t_1111, t_1112, pc_x, pfh_625, pfh_626, \
                         pfh_627, pfh_628, pfh_629, pgh_835, pgh_836, pgh_837, pgh_838, \
                         pgh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_0 * pfh_625[k]
                    + f_4 * pc_x[k] * pgh_835[k];

        t_1109[k] = f_0 * pfh_626[k]
                    + f_4 * pc_x[k] * pgh_836[k];

        t_1110[k] = f_0 * pfh_627[k]
                    + f_4 * pc_x[k] * pgh_837[k];

        t_1111[k] = f_0 * pfh_628[k]
                    + f_4 * pc_x[k] * pgh_838[k];

        t_1112[k] = f_0 * pfh_629[k]
                    + f_4 * pc_x[k] * pgh_839[k];
    }

#pragma omp simd aligned(t_1113, t_1114, t_1115, t_1116, pb_x, pc_x, pfi0_833, pfi0_834, \
                         pfi0_835, pfi0_836, pfi1_833, pfi1_834, pfi1_835, \
                         pfi1_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1113[k] = pb_x[k] * pfi0_833[k]
                    - f_11 * pc_x[k] * pfi1_833[k];

        t_1114[k] = pb_x[k] * pfi0_834[k]
                    - f_11 * pc_x[k] * pfi1_834[k];

        t_1115[k] = pb_x[k] * pfi0_835[k]
                    - f_11 * pc_x[k] * pfi1_835[k];

        t_1116[k] = pb_x[k] * pfi0_836[k]
                    - f_11 * pc_x[k] * pfi1_836[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, t_1120, pa_z, pb_x, pc_x, pc_y, pc_z, \
                         sgi0_280, sgi1_280, pfi0_837, pfi0_839, pfi1_837, pfi1_839, \
                         pgh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = pb_x[k] * pfi0_837[k]
                    - f_11 * pc_x[k] * pfi1_837[k];

        t_1118[k] = f_4 * pc_y[k] * pgh_839[k];

        t_1119[k] = pb_x[k] * pfi0_839[k]
                    - f_11 * pc_x[k] * pfi1_839[k];

        t_1120[k] = pa_z[k] * sgi0_280[k]
                    - f_11 * pc_z[k] * sgi1_280[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, t_1124, pa_z, pc_x, pc_y, pc_z, sgi0_283, \
                         sgi1_283, pfh_546, pfh_548, pgg0_602, pgg1_602, pgh_840, \
                         pgh_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_1 * pfh_546[k]
                    + f_4 * pc_y[k] * pgh_840[k];

        t_1122[k] = f_17 * pgg0_602[k]
                    - f_18 * pgg1_602[k]
                    + f_4 * pc_x[k] * pgh_842[k];

        t_1123[k] = pa_z[k] * sgi0_283[k]
                    - f_11 * pc_z[k] * sgi1_283[k];

        t_1124[k] = f_1 * pfh_548[k]
                    + f_4 * pc_y[k] * pgh_842[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, pa_z, pc_x, pc_z, sgi0_286, sgi1_286, \
                         pgg0_605, pgg0_607, pgg1_605, pgg1_607, pgh_845, \
                         pgh_847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = f_9 * pgg0_605[k]
                    - f_10 * pgg1_605[k]
                    + f_4 * pc_x[k] * pgh_845[k];

        t_1126[k] = pa_z[k] * sgi0_286[k]
                    - f_11 * pc_z[k] * sgi1_286[k];

        t_1127[k] = f_7 * pgg0_607[k]
                    - f_8 * pgg1_607[k]
                    + f_4 * pc_x[k] * pgh_847[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pa_z, pc_x, pc_y, pc_z, sgi0_290, sgi1_290, \
                         pfh_551, pgg0_609, pgg1_609, pgh_845, \
                         pgh_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_1 * pfh_551[k]
                    + f_4 * pc_y[k] * pgh_845[k];

        t_1129[k] = f_7 * pgg0_609[k]
                    - f_8 * pgg1_609[k]
                    + f_4 * pc_x[k] * pgh_849[k];

        t_1130[k] = pa_z[k] * sgi0_290[k]
                    - f_11 * pc_z[k] * sgi1_290[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pc_x, pc_y, pfh_555, pgg0_611, pgg0_612, \
                         pgg1_611, pgg1_612, pgh_849, pgh_851, \
                         pgh_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_5 * pgg0_611[k]
                    - f_6 * pgg1_611[k]
                    + f_4 * pc_x[k] * pgh_851[k];

        t_1132[k] = f_5 * pgg0_612[k]
                    - f_6 * pgg1_612[k]
                    + f_4 * pc_x[k] * pgh_852[k];

        t_1133[k] = f_1 * pfh_555[k]
                    + f_4 * pc_y[k] * pgh_849[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, t_1138, t_1139, pc_x, pgg0_614, \
                         pgg1_614, pgh_854, pgh_855, pgh_856, pgh_857, pgh_858, \
                         pgh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_5 * pgg0_614[k]
                    - f_6 * pgg1_614[k]
                    + f_4 * pc_x[k] * pgh_854[k];

        t_1135[k] = f_4 * pc_x[k] * pgh_855[k];

        t_1136[k] = f_4 * pc_x[k] * pgh_856[k];

        t_1137[k] = f_4 * pc_x[k] * pgh_857[k];

        t_1138[k] = f_4 * pc_x[k] * pgh_858[k];

        t_1139[k] = f_4 * pc_x[k] * pgh_859[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, pa_z, pc_x, pc_z, sgi0_301, sgi0_302, \
                         sgi0_303, sgh_225, sgh_226, sgi1_301, sgi1_302, sgi1_303, \
                         pgh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_4 * pc_x[k] * pgh_860[k];

        t_1141[k] = pa_z[k] * sgi0_301[k]
                    - f_11 * pc_z[k] * sgi1_301[k];

        t_1142[k] = pa_z[k] * sgi0_302[k]
                    + f_0 * sgh_225[k]
                    - f_11 * pc_z[k] * sgi1_302[k];

        t_1143[k] = pa_z[k] * sgi0_303[k]
                    + f_12 * sgh_226[k]
                    - f_11 * pc_z[k] * sgi1_303[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pa_z, pc_y, pc_z, sgi0_304, sgi0_305, \
                         sgh_227, sgh_228, sgi1_304, sgi1_305, pfh_566, \
                         pgh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = pa_z[k] * sgi0_304[k]
                    + f_13 * sgh_227[k]
                    - f_11 * pc_z[k] * sgi1_304[k];

        t_1145[k] = pa_z[k] * sgi0_305[k]
                    + f_1 * sgh_228[k]
                    - f_11 * pc_z[k] * sgi1_305[k];

        t_1146[k] = f_1 * pfh_566[k]
                    + f_4 * pc_y[k] * pgh_860[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, pa_z, pc_x, pc_y, pc_z, sgi0_307, sgh_230, \
                         sgi1_307, pfh_567, pgg0_615, pgg1_615, \
                         pgh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = pa_z[k] * sgi0_307[k]
                    + f_16 * sgh_230[k]
                    - f_11 * pc_z[k] * sgi1_307[k];

        t_1148[k] = f_2 * pgg0_615[k]
                    - f_3 * pgg1_615[k]
                    + f_4 * pc_x[k] * pgh_861[k];

        t_1149[k] = f_13 * pfh_567[k]
                    + f_4 * pc_y[k] * pgh_861[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, t_1153, pc_x, pc_y, pfh_569, pgg0_617, \
                         pgg0_618, pgg0_620, pgg1_617, pgg1_618, pgg1_620, pgh_863, pgh_864, \
                         pgh_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_17 * pgg0_617[k]
                    - f_18 * pgg1_617[k]
                    + f_4 * pc_x[k] * pgh_863[k];

        t_1151[k] = f_9 * pgg0_618[k]
                    - f_10 * pgg1_618[k]
                    + f_4 * pc_x[k] * pgh_864[k];

        t_1152[k] = f_13 * pfh_569[k]
                    + f_4 * pc_y[k] * pgh_863[k];

        t_1153[k] = f_9 * pgg0_620[k]
                    - f_10 * pgg1_620[k]
                    + f_4 * pc_x[k] * pgh_866[k];
    }

#pragma omp simd aligned(t_1154, t_1155, t_1156, pc_x, pc_y, pfh_572, pgg0_621, pgg0_622, \
                         pgg1_621, pgg1_622, pgh_866, pgh_867, \
                         pgh_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1154[k] = f_7 * pgg0_621[k]
                    - f_8 * pgg1_621[k]
                    + f_4 * pc_x[k] * pgh_867[k];

        t_1155[k] = f_7 * pgg0_622[k]
                    - f_8 * pgg1_622[k]
                    + f_4 * pc_x[k] * pgh_868[k];

        t_1156[k] = f_13 * pfh_572[k]
                    + f_4 * pc_y[k] * pgh_866[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, pc_x, pgg0_624, pgg0_625, pgg0_626, pgg1_624, \
                         pgg1_625, pgg1_626, pgh_870, pgh_871, \
                         pgh_872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_7 * pgg0_624[k]
                    - f_8 * pgg1_624[k]
                    + f_4 * pc_x[k] * pgh_870[k];

        t_1158[k] = f_5 * pgg0_625[k]
                    - f_6 * pgg1_625[k]
                    + f_4 * pc_x[k] * pgh_871[k];

        t_1159[k] = f_5 * pgg0_626[k]
                    - f_6 * pgg1_626[k]
                    + f_4 * pc_x[k] * pgh_872[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, t_1163, pc_x, pc_y, pfh_576, pgg0_627, \
                         pgg0_629, pgg1_627, pgg1_629, pgh_870, pgh_873, pgh_875, \
                         pgh_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_5 * pgg0_627[k]
                    - f_6 * pgg1_627[k]
                    + f_4 * pc_x[k] * pgh_873[k];

        t_1161[k] = f_13 * pfh_576[k]
                    + f_4 * pc_y[k] * pgh_870[k];

        t_1162[k] = f_5 * pgg0_629[k]
                    - f_6 * pgg1_629[k]
                    + f_4 * pc_x[k] * pgh_875[k];

        t_1163[k] = f_4 * pc_x[k] * pgh_876[k];
    }

#pragma omp simd aligned(t_1164, t_1165, t_1166, t_1167, t_1168, pc_x, pgh_877, pgh_878, \
                         pgh_879, pgh_880, pgh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1164[k] = f_4 * pc_x[k] * pgh_877[k];

        t_1165[k] = f_4 * pc_x[k] * pgh_878[k];

        t_1166[k] = f_4 * pc_x[k] * pgh_879[k];

        t_1167[k] = f_4 * pc_x[k] * pgh_880[k];

        t_1168[k] = f_4 * pc_x[k] * pgh_881[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, pc_y, pfh_582, pfh_583, pfh_584, pgg0_625, \
                         pgg0_626, pgg0_627, pgg1_625, pgg1_626, pgg1_627, pgh_876, pgh_877, \
                         pgh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_13 * pfh_582[k]
                    + f_2 * pgg0_625[k]
                    - f_3 * pgg1_625[k]
                    + f_4 * pc_y[k] * pgh_876[k];

        t_1170[k] = f_13 * pfh_583[k]
                    + f_17 * pgg0_626[k]
                    - f_18 * pgg1_626[k]
                    + f_4 * pc_y[k] * pgh_877[k];

        t_1171[k] = f_13 * pfh_584[k]
                    + f_9 * pgg0_627[k]
                    - f_10 * pgg1_627[k]
                    + f_4 * pc_y[k] * pgh_878[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, pc_y, pfh_585, pfh_586, pfh_587, pgg0_628, \
                         pgg0_629, pgg1_628, pgg1_629, pgh_879, pgh_880, \
                         pgh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = f_13 * pfh_585[k]
                    + f_7 * pgg0_628[k]
                    - f_8 * pgg1_628[k]
                    + f_4 * pc_y[k] * pgh_879[k];

        t_1173[k] = f_13 * pfh_586[k]
                    + f_5 * pgg0_629[k]
                    - f_6 * pgg1_629[k]
                    + f_4 * pc_y[k] * pgh_880[k];

        t_1174[k] = f_13 * pfh_587[k]
                    + f_4 * pc_y[k] * pgh_881[k];
    }

#pragma omp simd aligned(t_1175, t_1176, t_1177, pb_y, pc_x, pc_y, pdi0_475, pdi1_475, \
                         pfi0_783, pfh_588, pfi1_783, pgg0_630, pgg1_630, \
                         pgh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1175[k] = f_19 * pdi0_475[k]
                    - f_20 * pdi1_475[k]
                    + pb_y[k] * pfi0_783[k]
                    - f_11 * pc_y[k] * pfi1_783[k];

        t_1176[k] = f_2 * pgg0_630[k]
                    - f_3 * pgg1_630[k]
                    + f_4 * pc_x[k] * pgh_882[k];

        t_1177[k] = f_12 * pfh_588[k]
                    + f_4 * pc_y[k] * pgh_882[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, t_1181, pc_x, pc_y, pfh_590, pgg0_632, \
                         pgg0_633, pgg0_635, pgg1_632, pgg1_633, pgg1_635, pgh_884, pgh_885, \
                         pgh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = f_17 * pgg0_632[k]
                    - f_18 * pgg1_632[k]
                    + f_4 * pc_x[k] * pgh_884[k];

        t_1179[k] = f_9 * pgg0_633[k]
                    - f_10 * pgg1_633[k]
                    + f_4 * pc_x[k] * pgh_885[k];

        t_1180[k] = f_12 * pfh_590[k]
                    + f_4 * pc_y[k] * pgh_884[k];

        t_1181[k] = f_9 * pgg0_635[k]
                    - f_10 * pgg1_635[k]
                    + f_4 * pc_x[k] * pgh_887[k];
    }
}

static auto
compute_prim_pgi_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sgh,
                                                           const size_t pdi0, const size_t pdi1,
                                                           const size_t pfi0, const size_t pfh,
                                                           const size_t pfi1, const size_t pgg0,
                                                           const size_t pgg1, const size_t pgh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_21 = 2.5 / q;

    auto *t_1182 = buffer.data(target + 1182);
    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);
    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh_314 = buffer.data(sgh + 314);

    const auto *pdi0_503 = buffer.data(pdi0 + 503);

    const auto *pdi1_503 = buffer.data(pdi1 + 503);

    const auto *pfi0_811 = buffer.data(pfi0 + 811);
    const auto *pfi0_812 = buffer.data(pfi0 + 812);
    const auto *pfi0_814 = buffer.data(pfi0 + 814);
    const auto *pfi0_817 = buffer.data(pfi0 + 817);
    const auto *pfi0_821 = buffer.data(pfi0 + 821);
    const auto *pfi0_826 = buffer.data(pfi0 + 826);
    const auto *pfi0_833 = buffer.data(pfi0 + 833);
    const auto *pfi0_834 = buffer.data(pfi0 + 834);
    const auto *pfi0_835 = buffer.data(pfi0 + 835);
    const auto *pfi0_836 = buffer.data(pfi0 + 836);
    const auto *pfi0_837 = buffer.data(pfi0 + 837);
    const auto *pfi0_839 = buffer.data(pfi0 + 839);

    const auto *pfh_593 = buffer.data(pfh + 593);
    const auto *pfh_597 = buffer.data(pfh + 597);
    const auto *pfh_603 = buffer.data(pfh + 603);
    const auto *pfh_604 = buffer.data(pfh + 604);
    const auto *pfh_605 = buffer.data(pfh + 605);
    const auto *pfh_606 = buffer.data(pfh + 606);
    const auto *pfh_607 = buffer.data(pfh + 607);
    const auto *pfh_608 = buffer.data(pfh + 608);
    const auto *pfh_609 = buffer.data(pfh + 609);
    const auto *pfh_611 = buffer.data(pfh + 611);
    const auto *pfh_614 = buffer.data(pfh + 614);
    const auto *pfh_618 = buffer.data(pfh + 618);
    const auto *pfh_624 = buffer.data(pfh + 624);
    const auto *pfh_625 = buffer.data(pfh + 625);
    const auto *pfh_626 = buffer.data(pfh + 626);
    const auto *pfh_627 = buffer.data(pfh + 627);
    const auto *pfh_628 = buffer.data(pfh + 628);
    const auto *pfh_629 = buffer.data(pfh + 629);

    const auto *pfi1_811 = buffer.data(pfi1 + 811);
    const auto *pfi1_812 = buffer.data(pfi1 + 812);
    const auto *pfi1_814 = buffer.data(pfi1 + 814);
    const auto *pfi1_817 = buffer.data(pfi1 + 817);
    const auto *pfi1_821 = buffer.data(pfi1 + 821);
    const auto *pfi1_826 = buffer.data(pfi1 + 826);
    const auto *pfi1_833 = buffer.data(pfi1 + 833);
    const auto *pfi1_834 = buffer.data(pfi1 + 834);
    const auto *pfi1_835 = buffer.data(pfi1 + 835);
    const auto *pfi1_836 = buffer.data(pfi1 + 836);
    const auto *pfi1_837 = buffer.data(pfi1 + 837);
    const auto *pfi1_839 = buffer.data(pfi1 + 839);

    const auto *pgg0_636 = buffer.data(pgg0 + 636);
    const auto *pgg0_637 = buffer.data(pgg0 + 637);
    const auto *pgg0_639 = buffer.data(pgg0 + 639);
    const auto *pgg0_640 = buffer.data(pgg0 + 640);
    const auto *pgg0_641 = buffer.data(pgg0 + 641);
    const auto *pgg0_642 = buffer.data(pgg0 + 642);
    const auto *pgg0_643 = buffer.data(pgg0 + 643);
    const auto *pgg0_644 = buffer.data(pgg0 + 644);
    const auto *pgg0_648 = buffer.data(pgg0 + 648);
    const auto *pgg0_651 = buffer.data(pgg0 + 651);
    const auto *pgg0_652 = buffer.data(pgg0 + 652);
    const auto *pgg0_655 = buffer.data(pgg0 + 655);
    const auto *pgg0_656 = buffer.data(pgg0 + 656);
    const auto *pgg0_657 = buffer.data(pgg0 + 657);
    const auto *pgg0_660 = buffer.data(pgg0 + 660);
    const auto *pgg0_662 = buffer.data(pgg0 + 662);
    const auto *pgg0_663 = buffer.data(pgg0 + 663);
    const auto *pgg0_665 = buffer.data(pgg0 + 665);
    const auto *pgg0_666 = buffer.data(pgg0 + 666);
    const auto *pgg0_667 = buffer.data(pgg0 + 667);
    const auto *pgg0_669 = buffer.data(pgg0 + 669);
    const auto *pgg0_670 = buffer.data(pgg0 + 670);
    const auto *pgg0_671 = buffer.data(pgg0 + 671);
    const auto *pgg0_672 = buffer.data(pgg0 + 672);
    const auto *pgg0_673 = buffer.data(pgg0 + 673);
    const auto *pgg0_674 = buffer.data(pgg0 + 674);

    const auto *pgg1_636 = buffer.data(pgg1 + 636);
    const auto *pgg1_637 = buffer.data(pgg1 + 637);
    const auto *pgg1_639 = buffer.data(pgg1 + 639);
    const auto *pgg1_640 = buffer.data(pgg1 + 640);
    const auto *pgg1_641 = buffer.data(pgg1 + 641);
    const auto *pgg1_642 = buffer.data(pgg1 + 642);
    const auto *pgg1_643 = buffer.data(pgg1 + 643);
    const auto *pgg1_644 = buffer.data(pgg1 + 644);
    const auto *pgg1_648 = buffer.data(pgg1 + 648);
    const auto *pgg1_651 = buffer.data(pgg1 + 651);
    const auto *pgg1_652 = buffer.data(pgg1 + 652);
    const auto *pgg1_655 = buffer.data(pgg1 + 655);
    const auto *pgg1_656 = buffer.data(pgg1 + 656);
    const auto *pgg1_657 = buffer.data(pgg1 + 657);
    const auto *pgg1_660 = buffer.data(pgg1 + 660);
    const auto *pgg1_662 = buffer.data(pgg1 + 662);
    const auto *pgg1_663 = buffer.data(pgg1 + 663);
    const auto *pgg1_665 = buffer.data(pgg1 + 665);
    const auto *pgg1_666 = buffer.data(pgg1 + 666);
    const auto *pgg1_667 = buffer.data(pgg1 + 667);
    const auto *pgg1_669 = buffer.data(pgg1 + 669);
    const auto *pgg1_670 = buffer.data(pgg1 + 670);
    const auto *pgg1_671 = buffer.data(pgg1 + 671);
    const auto *pgg1_672 = buffer.data(pgg1 + 672);
    const auto *pgg1_673 = buffer.data(pgg1 + 673);
    const auto *pgg1_674 = buffer.data(pgg1 + 674);

    const auto *pgh_887 = buffer.data(pgh + 887);
    const auto *pgh_888 = buffer.data(pgh + 888);
    const auto *pgh_889 = buffer.data(pgh + 889);
    const auto *pgh_891 = buffer.data(pgh + 891);
    const auto *pgh_892 = buffer.data(pgh + 892);
    const auto *pgh_893 = buffer.data(pgh + 893);
    const auto *pgh_894 = buffer.data(pgh + 894);
    const auto *pgh_896 = buffer.data(pgh + 896);
    const auto *pgh_897 = buffer.data(pgh + 897);
    const auto *pgh_898 = buffer.data(pgh + 898);
    const auto *pgh_899 = buffer.data(pgh + 899);
    const auto *pgh_900 = buffer.data(pgh + 900);
    const auto *pgh_901 = buffer.data(pgh + 901);
    const auto *pgh_902 = buffer.data(pgh + 902);
    const auto *pgh_903 = buffer.data(pgh + 903);
    const auto *pgh_905 = buffer.data(pgh + 905);
    const auto *pgh_906 = buffer.data(pgh + 906);
    const auto *pgh_908 = buffer.data(pgh + 908);
    const auto *pgh_909 = buffer.data(pgh + 909);
    const auto *pgh_910 = buffer.data(pgh + 910);
    const auto *pgh_912 = buffer.data(pgh + 912);
    const auto *pgh_913 = buffer.data(pgh + 913);
    const auto *pgh_914 = buffer.data(pgh + 914);
    const auto *pgh_915 = buffer.data(pgh + 915);
    const auto *pgh_918 = buffer.data(pgh + 918);
    const auto *pgh_919 = buffer.data(pgh + 919);
    const auto *pgh_920 = buffer.data(pgh + 920);
    const auto *pgh_921 = buffer.data(pgh + 921);
    const auto *pgh_922 = buffer.data(pgh + 922);
    const auto *pgh_923 = buffer.data(pgh + 923);
    const auto *pgh_924 = buffer.data(pgh + 924);
    const auto *pgh_926 = buffer.data(pgh + 926);
    const auto *pgh_927 = buffer.data(pgh + 927);
    const auto *pgh_929 = buffer.data(pgh + 929);
    const auto *pgh_930 = buffer.data(pgh + 930);
    const auto *pgh_931 = buffer.data(pgh + 931);
    const auto *pgh_933 = buffer.data(pgh + 933);
    const auto *pgh_934 = buffer.data(pgh + 934);
    const auto *pgh_935 = buffer.data(pgh + 935);
    const auto *pgh_936 = buffer.data(pgh + 936);
    const auto *pgh_938 = buffer.data(pgh + 938);
    const auto *pgh_939 = buffer.data(pgh + 939);
    const auto *pgh_940 = buffer.data(pgh + 940);
    const auto *pgh_941 = buffer.data(pgh + 941);
    const auto *pgh_942 = buffer.data(pgh + 942);
    const auto *pgh_943 = buffer.data(pgh + 943);
    const auto *pgh_944 = buffer.data(pgh + 944);

#pragma omp simd aligned(t_1182, t_1183, t_1184, pc_x, pc_y, pfh_593, pgg0_636, pgg0_637, \
                         pgg1_636, pgg1_637, pgh_887, pgh_888, \
                         pgh_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = f_7 * pgg0_636[k]
                    - f_8 * pgg1_636[k]
                    + f_4 * pc_x[k] * pgh_888[k];

        t_1183[k] = f_7 * pgg0_637[k]
                    - f_8 * pgg1_637[k]
                    + f_4 * pc_x[k] * pgh_889[k];

        t_1184[k] = f_12 * pfh_593[k]
                    + f_4 * pc_y[k] * pgh_887[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pc_x, pgg0_639, pgg0_640, pgg0_641, pgg1_639, \
                         pgg1_640, pgg1_641, pgh_891, pgh_892, \
                         pgh_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_7 * pgg0_639[k]
                    - f_8 * pgg1_639[k]
                    + f_4 * pc_x[k] * pgh_891[k];

        t_1186[k] = f_5 * pgg0_640[k]
                    - f_6 * pgg1_640[k]
                    + f_4 * pc_x[k] * pgh_892[k];

        t_1187[k] = f_5 * pgg0_641[k]
                    - f_6 * pgg1_641[k]
                    + f_4 * pc_x[k] * pgh_893[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, t_1191, pc_x, pc_y, pfh_597, pgg0_642, \
                         pgg0_644, pgg1_642, pgg1_644, pgh_891, pgh_894, pgh_896, \
                         pgh_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = f_5 * pgg0_642[k]
                    - f_6 * pgg1_642[k]
                    + f_4 * pc_x[k] * pgh_894[k];

        t_1189[k] = f_12 * pfh_597[k]
                    + f_4 * pc_y[k] * pgh_891[k];

        t_1190[k] = f_5 * pgg0_644[k]
                    - f_6 * pgg1_644[k]
                    + f_4 * pc_x[k] * pgh_896[k];

        t_1191[k] = f_4 * pc_x[k] * pgh_897[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, t_1195, t_1196, pc_x, pgh_898, pgh_899, \
                         pgh_900, pgh_901, pgh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = f_4 * pc_x[k] * pgh_898[k];

        t_1193[k] = f_4 * pc_x[k] * pgh_899[k];

        t_1194[k] = f_4 * pc_x[k] * pgh_900[k];

        t_1195[k] = f_4 * pc_x[k] * pgh_901[k];

        t_1196[k] = f_4 * pc_x[k] * pgh_902[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, pc_y, pfh_603, pfh_604, pfh_605, pgg0_640, \
                         pgg0_641, pgg0_642, pgg1_640, pgg1_641, pgg1_642, pgh_897, pgh_898, \
                         pgh_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = f_12 * pfh_603[k]
                    + f_2 * pgg0_640[k]
                    - f_3 * pgg1_640[k]
                    + f_4 * pc_y[k] * pgh_897[k];

        t_1198[k] = f_12 * pfh_604[k]
                    + f_17 * pgg0_641[k]
                    - f_18 * pgg1_641[k]
                    + f_4 * pc_y[k] * pgh_898[k];

        t_1199[k] = f_12 * pfh_605[k]
                    + f_9 * pgg0_642[k]
                    - f_10 * pgg1_642[k]
                    + f_4 * pc_y[k] * pgh_899[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, pc_y, pfh_606, pfh_607, pfh_608, pgg0_643, \
                         pgg0_644, pgg1_643, pgg1_644, pgh_900, pgh_901, \
                         pgh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = f_12 * pfh_606[k]
                    + f_7 * pgg0_643[k]
                    - f_8 * pgg1_643[k]
                    + f_4 * pc_y[k] * pgh_900[k];

        t_1201[k] = f_12 * pfh_607[k]
                    + f_5 * pgg0_644[k]
                    - f_6 * pgg1_644[k]
                    + f_4 * pc_y[k] * pgh_901[k];

        t_1202[k] = f_12 * pfh_608[k]
                    + f_4 * pc_y[k] * pgh_902[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, t_1206, pb_y, pc_y, pdi0_503, pdi1_503, \
                         pfi0_811, pfi0_812, pfi0_814, pfh_609, pfi1_811, pfi1_812, pfi1_814, \
                         pgh_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = f_14 * pdi0_503[k]
                    - f_15 * pdi1_503[k]
                    + pb_y[k] * pfi0_811[k]
                    - f_11 * pc_y[k] * pfi1_811[k];

        t_1204[k] = pb_y[k] * pfi0_812[k]
                    - f_11 * pc_y[k] * pfi1_812[k];

        t_1205[k] = f_0 * pfh_609[k]
                    + f_4 * pc_y[k] * pgh_903[k];

        t_1206[k] = pb_y[k] * pfi0_814[k]
                    - f_11 * pc_y[k] * pfi1_814[k];
    }

#pragma omp simd aligned(t_1207, t_1208, t_1209, pb_y, pc_x, pc_y, pfi0_817, pfh_611, \
                         pfi1_817, pgg0_648, pgg1_648, pgh_905, \
                         pgh_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1207[k] = f_9 * pgg0_648[k]
                    - f_10 * pgg1_648[k]
                    + f_4 * pc_x[k] * pgh_906[k];

        t_1208[k] = f_0 * pfh_611[k]
                    + f_4 * pc_y[k] * pgh_905[k];

        t_1209[k] = pb_y[k] * pfi0_817[k]
                    - f_11 * pc_y[k] * pfi1_817[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, pc_x, pc_y, pfh_614, pgg0_651, pgg0_652, \
                         pgg1_651, pgg1_652, pgh_908, pgh_909, \
                         pgh_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = f_7 * pgg0_651[k]
                    - f_8 * pgg1_651[k]
                    + f_4 * pc_x[k] * pgh_909[k];

        t_1211[k] = f_7 * pgg0_652[k]
                    - f_8 * pgg1_652[k]
                    + f_4 * pc_x[k] * pgh_910[k];

        t_1212[k] = f_0 * pfh_614[k]
                    + f_4 * pc_y[k] * pgh_908[k];
    }

#pragma omp simd aligned(t_1213, t_1214, t_1215, pb_y, pc_x, pc_y, pfi0_821, pfi1_821, \
                         pgg0_655, pgg0_656, pgg1_655, pgg1_656, pgh_913, \
                         pgh_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1213[k] = pb_y[k] * pfi0_821[k]
                    - f_11 * pc_y[k] * pfi1_821[k];

        t_1214[k] = f_5 * pgg0_655[k]
                    - f_6 * pgg1_655[k]
                    + f_4 * pc_x[k] * pgh_913[k];

        t_1215[k] = f_5 * pgg0_656[k]
                    - f_6 * pgg1_656[k]
                    + f_4 * pc_x[k] * pgh_914[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, t_1219, pb_y, pc_x, pc_y, pfi0_826, pfh_618, \
                         pfi1_826, pgg0_657, pgg1_657, pgh_912, pgh_915, \
                         pgh_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_5 * pgg0_657[k]
                    - f_6 * pgg1_657[k]
                    + f_4 * pc_x[k] * pgh_915[k];

        t_1217[k] = f_0 * pfh_618[k]
                    + f_4 * pc_y[k] * pgh_912[k];

        t_1218[k] = pb_y[k] * pfi0_826[k]
                    - f_11 * pc_y[k] * pfi1_826[k];

        t_1219[k] = f_4 * pc_x[k] * pgh_918[k];
    }

#pragma omp simd aligned(t_1220, t_1221, t_1222, t_1223, t_1224, pc_x, pgh_919, pgh_920, \
                         pgh_921, pgh_922, pgh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1220[k] = f_4 * pc_x[k] * pgh_919[k];

        t_1221[k] = f_4 * pc_x[k] * pgh_920[k];

        t_1222[k] = f_4 * pc_x[k] * pgh_921[k];

        t_1223[k] = f_4 * pc_x[k] * pgh_922[k];

        t_1224[k] = f_4 * pc_x[k] * pgh_923[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, pb_y, pc_y, pfi0_833, pfi0_834, pfi0_835, \
                         pfh_624, pfh_625, pfh_626, pfi1_833, pfi1_834, \
                         pfi1_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = pb_y[k] * pfi0_833[k]
                    + f_16 * pfh_624[k]
                    - f_11 * pc_y[k] * pfi1_833[k];

        t_1226[k] = pb_y[k] * pfi0_834[k]
                    + f_21 * pfh_625[k]
                    - f_11 * pc_y[k] * pfi1_834[k];

        t_1227[k] = pb_y[k] * pfi0_835[k]
                    + f_1 * pfh_626[k]
                    - f_11 * pc_y[k] * pfi1_835[k];
    }

#pragma omp simd aligned(t_1228, t_1229, t_1230, t_1231, pb_y, pc_y, pfi0_836, pfi0_837, \
                         pfi0_839, pfh_627, pfh_628, pfh_629, pfi1_836, pfi1_837, pfi1_839, \
                         pgh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1228[k] = pb_y[k] * pfi0_836[k]
                    + f_13 * pfh_627[k]
                    - f_11 * pc_y[k] * pfi1_836[k];

        t_1229[k] = pb_y[k] * pfi0_837[k]
                    + f_12 * pfh_628[k]
                    - f_11 * pc_y[k] * pfi1_837[k];

        t_1230[k] = f_0 * pfh_629[k]
                    + f_4 * pc_y[k] * pgh_923[k];

        t_1231[k] = pb_y[k] * pfi0_839[k]
                    - f_11 * pc_y[k] * pfi1_839[k];
    }

#pragma omp simd aligned(t_1232, t_1233, t_1234, t_1235, t_1236, pc_x, pc_y, pgg0_660, \
                         pgg0_662, pgg0_663, pgg1_660, pgg1_662, pgg1_663, pgh_924, pgh_926, \
                         pgh_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1232[k] = f_2 * pgg0_660[k]
                    - f_3 * pgg1_660[k]
                    + f_4 * pc_x[k] * pgh_924[k];

        t_1233[k] = f_4 * pc_y[k] * pgh_924[k];

        t_1234[k] = f_17 * pgg0_662[k]
                    - f_18 * pgg1_662[k]
                    + f_4 * pc_x[k] * pgh_926[k];

        t_1235[k] = f_9 * pgg0_663[k]
                    - f_10 * pgg1_663[k]
                    + f_4 * pc_x[k] * pgh_927[k];

        t_1236[k] = f_4 * pc_y[k] * pgh_926[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, pc_x, pc_y, pgg0_665, pgg0_666, \
                         pgg0_667, pgg1_665, pgg1_666, pgg1_667, pgh_929, pgh_930, \
                         pgh_931 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_9 * pgg0_665[k]
                    - f_10 * pgg1_665[k]
                    + f_4 * pc_x[k] * pgh_929[k];

        t_1238[k] = f_7 * pgg0_666[k]
                    - f_8 * pgg1_666[k]
                    + f_4 * pc_x[k] * pgh_930[k];

        t_1239[k] = f_7 * pgg0_667[k]
                    - f_8 * pgg1_667[k]
                    + f_4 * pc_x[k] * pgh_931[k];

        t_1240[k] = f_4 * pc_y[k] * pgh_929[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, pgg0_669, pgg0_670, pgg0_671, pgg1_669, \
                         pgg1_670, pgg1_671, pgh_933, pgh_934, \
                         pgh_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_7 * pgg0_669[k]
                    - f_8 * pgg1_669[k]
                    + f_4 * pc_x[k] * pgh_933[k];

        t_1242[k] = f_5 * pgg0_670[k]
                    - f_6 * pgg1_670[k]
                    + f_4 * pc_x[k] * pgh_934[k];

        t_1243[k] = f_5 * pgg0_671[k]
                    - f_6 * pgg1_671[k]
                    + f_4 * pc_x[k] * pgh_935[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, t_1247, t_1248, pc_x, pc_y, pgg0_672, \
                         pgg0_674, pgg1_672, pgg1_674, pgh_933, pgh_936, pgh_938, pgh_939, \
                         pgh_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_5 * pgg0_672[k]
                    - f_6 * pgg1_672[k]
                    + f_4 * pc_x[k] * pgh_936[k];

        t_1245[k] = f_4 * pc_y[k] * pgh_933[k];

        t_1246[k] = f_5 * pgg0_674[k]
                    - f_6 * pgg1_674[k]
                    + f_4 * pc_x[k] * pgh_938[k];

        t_1247[k] = f_4 * pc_x[k] * pgh_939[k];

        t_1248[k] = f_4 * pc_x[k] * pgh_940[k];
    }

#pragma omp simd aligned(t_1249, t_1250, t_1251, t_1252, t_1253, pc_x, pc_y, pgg0_670, \
                         pgg1_670, pgh_939, pgh_941, pgh_942, pgh_943, \
                         pgh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1249[k] = f_4 * pc_x[k] * pgh_941[k];

        t_1250[k] = f_4 * pc_x[k] * pgh_942[k];

        t_1251[k] = f_4 * pc_x[k] * pgh_943[k];

        t_1252[k] = f_4 * pc_x[k] * pgh_944[k];

        t_1253[k] = f_2 * pgg0_670[k]
                    - f_3 * pgg1_670[k]
                    + f_4 * pc_y[k] * pgh_939[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, pgg0_671, pgg0_672, pgg0_673, pgg1_671, \
                         pgg1_672, pgg1_673, pgh_940, pgh_941, \
                         pgh_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_17 * pgg0_671[k]
                    - f_18 * pgg1_671[k]
                    + f_4 * pc_y[k] * pgh_940[k];

        t_1255[k] = f_9 * pgg0_672[k]
                    - f_10 * pgg1_672[k]
                    + f_4 * pc_y[k] * pgh_941[k];

        t_1256[k] = f_7 * pgg0_673[k]
                    - f_8 * pgg1_673[k]
                    + f_4 * pc_y[k] * pgh_942[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pc_y, pc_z, sgh_314, pfh_629, pgg0_674, \
                         pgg1_674, pgh_943, pgh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_5 * pgg0_674[k]
                    - f_6 * pgg1_674[k]
                    + f_4 * pc_y[k] * pgh_943[k];

        t_1258[k] = f_4 * pc_y[k] * pgh_944[k];

        t_1259[k] = f_0 * sgh_314[k]
                    + f_1 * pfh_629[k]
                    + f_2 * pgg0_674[k]
                    - f_3 * pgg1_674[k]
                    + f_4 * pc_z[k] * pgh_944[k];
    }
}

auto
compute_prim_pgi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sgi0,
                                                   const size_t sgh, const size_t sgi1,
                                                   const size_t pdi0, const size_t pdi1,
                                                   const size_t pfi0, const size_t pfh,
                                                   const size_t pfi1, const size_t pgg0,
                                                   const size_t pgg1, const size_t pgh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pgi_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sgh, pdi0,
                                                              pdi1, pfi0, pfh, pfi1, pgg0, pgg1,
                                                              pgh, ncols, gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sgh, pdi0,
                                                              pdi1, pfi0, pfh, pfi1, pgg0, pgg1,
                                                              pgh, ncols, gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sgi0,
                                                              sgh, sgi1, pfi0, pfh, pfi1, pgg0,
                                                              pgg1, pgh, ncols, gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, sgi0,
                                                              sgh, sgi1, pdi0, pdi1, pfi0, pfh,
                                                              pfi1, pgg0, pgg1, pgh, ncols,
                                                              gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, sgi0,
                                                              sgh, sgi1, pdi0, pdi1, pfi0, pfh,
                                                              pfi1, pgg0, pgg1, pgh, ncols,
                                                              gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pb, pc, sgi0,
                                                              sgh, sgi1, pfi0, pfh, pfi1, pgg0,
                                                              pgg1, pgh, ncols, gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pb, pc, sgi0,
                                                              sgh, sgi1, pdi0, pdi1, pfi0, pfh,
                                                              pfi1, pgg0, pgg1, pgh, ncols,
                                                              gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece7(buffer, target, pa, pb, pc, sgi0,
                                                              sgh, sgi1, pdi0, pdi1, pfi0, pfh,
                                                              pfi1, pgg0, pgg1, pgh, ncols,
                                                              gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece8(buffer, target, pa, pb, pc, sgi0,
                                                              sgh, sgi1, pdi0, pdi1, pfi0, pfh,
                                                              pfi1, pgg0, pgg1, pgh, ncols,
                                                              gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece9(buffer, target, pa, pb, pc, sgi0,
                                                              sgh, sgi1, pdi0, pdi1, pfi0, pfh,
                                                              pfi1, pgg0, pgg1, pgh, ncols,
                                                              gamma, p, q);

    compute_prim_pgi_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, sgh,
                                                               pdi0, pdi1, pfi0, pfh, pfi1,
                                                               pgg0, pgg1, pgh, ncols, gamma, p,
                                                               q);
}

}  // namespace simdt3ceri
