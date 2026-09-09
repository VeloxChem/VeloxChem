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


#include "SimdThreeCenterElectronRepulsionVrrRecPSL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_psl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ssl0,
                                                          const size_t ssk, const size_t ssl1,
                                                          const size_t psi0, const size_t psi1,
                                                          const size_t psk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_11 = 2.5 / gamma;
    const auto f_12 = 2.5 * p / (gamma * q);
    const auto f_13 = 0.5 / q;
    const auto f_14 = 3.0 / gamma;
    const auto f_15 = 3.0 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssl0_0 = buffer.data(ssl0 + 0);
    const auto *ssl0_3 = buffer.data(ssl0 + 3);
    const auto *ssl0_5 = buffer.data(ssl0 + 5);
    const auto *ssl0_6 = buffer.data(ssl0 + 6);
    const auto *ssl0_9 = buffer.data(ssl0 + 9);
    const auto *ssl0_10 = buffer.data(ssl0 + 10);
    const auto *ssl0_14 = buffer.data(ssl0 + 14);
    const auto *ssl0_15 = buffer.data(ssl0 + 15);
    const auto *ssl0_20 = buffer.data(ssl0 + 20);
    const auto *ssl0_21 = buffer.data(ssl0 + 21);
    const auto *ssl0_27 = buffer.data(ssl0 + 27);
    const auto *ssl0_36 = buffer.data(ssl0 + 36);
    const auto *ssl0_38 = buffer.data(ssl0 + 38);
    const auto *ssl0_39 = buffer.data(ssl0 + 39);
    const auto *ssl0_40 = buffer.data(ssl0 + 40);
    const auto *ssl0_41 = buffer.data(ssl0 + 41);
    const auto *ssl0_42 = buffer.data(ssl0 + 42);
    const auto *ssl0_44 = buffer.data(ssl0 + 44);

    const auto *ssk_0 = buffer.data(ssk + 0);
    const auto *ssk_28 = buffer.data(ssk + 28);
    const auto *ssk_30 = buffer.data(ssk + 30);
    const auto *ssk_31 = buffer.data(ssk + 31);
    const auto *ssk_32 = buffer.data(ssk + 32);
    const auto *ssk_33 = buffer.data(ssk + 33);
    const auto *ssk_35 = buffer.data(ssk + 35);

    const auto *ssl1_0 = buffer.data(ssl1 + 0);
    const auto *ssl1_3 = buffer.data(ssl1 + 3);
    const auto *ssl1_5 = buffer.data(ssl1 + 5);
    const auto *ssl1_6 = buffer.data(ssl1 + 6);
    const auto *ssl1_9 = buffer.data(ssl1 + 9);
    const auto *ssl1_10 = buffer.data(ssl1 + 10);
    const auto *ssl1_14 = buffer.data(ssl1 + 14);
    const auto *ssl1_15 = buffer.data(ssl1 + 15);
    const auto *ssl1_20 = buffer.data(ssl1 + 20);
    const auto *ssl1_21 = buffer.data(ssl1 + 21);
    const auto *ssl1_27 = buffer.data(ssl1 + 27);
    const auto *ssl1_36 = buffer.data(ssl1 + 36);
    const auto *ssl1_38 = buffer.data(ssl1 + 38);
    const auto *ssl1_39 = buffer.data(ssl1 + 39);
    const auto *ssl1_40 = buffer.data(ssl1 + 40);
    const auto *ssl1_41 = buffer.data(ssl1 + 41);
    const auto *ssl1_42 = buffer.data(ssl1 + 42);
    const auto *ssl1_44 = buffer.data(ssl1 + 44);

    const auto *psi0_0 = buffer.data(psi0 + 0);
    const auto *psi0_1 = buffer.data(psi0 + 1);
    const auto *psi0_2 = buffer.data(psi0 + 2);
    const auto *psi0_3 = buffer.data(psi0 + 3);
    const auto *psi0_5 = buffer.data(psi0 + 5);
    const auto *psi0_6 = buffer.data(psi0 + 6);
    const auto *psi0_8 = buffer.data(psi0 + 8);
    const auto *psi0_9 = buffer.data(psi0 + 9);
    const auto *psi0_10 = buffer.data(psi0 + 10);
    const auto *psi0_12 = buffer.data(psi0 + 12);
    const auto *psi0_13 = buffer.data(psi0 + 13);
    const auto *psi0_14 = buffer.data(psi0 + 14);
    const auto *psi0_29 = buffer.data(psi0 + 29);
    const auto *psi0_31 = buffer.data(psi0 + 31);
    const auto *psi0_34 = buffer.data(psi0 + 34);
    const auto *psi0_36 = buffer.data(psi0 + 36);
    const auto *psi0_38 = buffer.data(psi0 + 38);
    const auto *psi0_40 = buffer.data(psi0 + 40);
    const auto *psi0_41 = buffer.data(psi0 + 41);
    const auto *psi0_43 = buffer.data(psi0 + 43);
    const auto *psi0_45 = buffer.data(psi0 + 45);
    const auto *psi0_46 = buffer.data(psi0 + 46);
    const auto *psi0_47 = buffer.data(psi0 + 47);
    const auto *psi0_49 = buffer.data(psi0 + 49);
    const auto *psi0_50 = buffer.data(psi0 + 50);
    const auto *psi0_51 = buffer.data(psi0 + 51);
    const auto *psi0_52 = buffer.data(psi0 + 52);
    const auto *psi0_53 = buffer.data(psi0 + 53);
    const auto *psi0_54 = buffer.data(psi0 + 54);
    const auto *psi0_58 = buffer.data(psi0 + 58);
    const auto *psi0_61 = buffer.data(psi0 + 61);
    const auto *psi0_63 = buffer.data(psi0 + 63);
    const auto *psi0_65 = buffer.data(psi0 + 65);
    const auto *psi0_67 = buffer.data(psi0 + 67);
    const auto *psi0_68 = buffer.data(psi0 + 68);
    const auto *psi0_70 = buffer.data(psi0 + 70);
    const auto *psi0_72 = buffer.data(psi0 + 72);
    const auto *psi0_73 = buffer.data(psi0 + 73);
    const auto *psi0_74 = buffer.data(psi0 + 74);
    const auto *psi0_76 = buffer.data(psi0 + 76);
    const auto *psi0_78 = buffer.data(psi0 + 78);
    const auto *psi0_79 = buffer.data(psi0 + 79);
    const auto *psi0_80 = buffer.data(psi0 + 80);
    const auto *psi0_81 = buffer.data(psi0 + 81);
    const auto *psi0_83 = buffer.data(psi0 + 83);

    const auto *psi1_0 = buffer.data(psi1 + 0);
    const auto *psi1_1 = buffer.data(psi1 + 1);
    const auto *psi1_2 = buffer.data(psi1 + 2);
    const auto *psi1_3 = buffer.data(psi1 + 3);
    const auto *psi1_5 = buffer.data(psi1 + 5);
    const auto *psi1_6 = buffer.data(psi1 + 6);
    const auto *psi1_8 = buffer.data(psi1 + 8);
    const auto *psi1_9 = buffer.data(psi1 + 9);
    const auto *psi1_10 = buffer.data(psi1 + 10);
    const auto *psi1_12 = buffer.data(psi1 + 12);
    const auto *psi1_13 = buffer.data(psi1 + 13);
    const auto *psi1_14 = buffer.data(psi1 + 14);
    const auto *psi1_29 = buffer.data(psi1 + 29);
    const auto *psi1_31 = buffer.data(psi1 + 31);
    const auto *psi1_34 = buffer.data(psi1 + 34);
    const auto *psi1_36 = buffer.data(psi1 + 36);
    const auto *psi1_38 = buffer.data(psi1 + 38);
    const auto *psi1_40 = buffer.data(psi1 + 40);
    const auto *psi1_41 = buffer.data(psi1 + 41);
    const auto *psi1_43 = buffer.data(psi1 + 43);
    const auto *psi1_45 = buffer.data(psi1 + 45);
    const auto *psi1_46 = buffer.data(psi1 + 46);
    const auto *psi1_47 = buffer.data(psi1 + 47);
    const auto *psi1_49 = buffer.data(psi1 + 49);
    const auto *psi1_50 = buffer.data(psi1 + 50);
    const auto *psi1_51 = buffer.data(psi1 + 51);
    const auto *psi1_52 = buffer.data(psi1 + 52);
    const auto *psi1_53 = buffer.data(psi1 + 53);
    const auto *psi1_54 = buffer.data(psi1 + 54);
    const auto *psi1_58 = buffer.data(psi1 + 58);
    const auto *psi1_61 = buffer.data(psi1 + 61);
    const auto *psi1_63 = buffer.data(psi1 + 63);
    const auto *psi1_65 = buffer.data(psi1 + 65);
    const auto *psi1_67 = buffer.data(psi1 + 67);
    const auto *psi1_68 = buffer.data(psi1 + 68);
    const auto *psi1_70 = buffer.data(psi1 + 70);
    const auto *psi1_72 = buffer.data(psi1 + 72);
    const auto *psi1_73 = buffer.data(psi1 + 73);
    const auto *psi1_74 = buffer.data(psi1 + 74);
    const auto *psi1_76 = buffer.data(psi1 + 76);
    const auto *psi1_78 = buffer.data(psi1 + 78);
    const auto *psi1_79 = buffer.data(psi1 + 79);
    const auto *psi1_80 = buffer.data(psi1 + 80);
    const auto *psi1_81 = buffer.data(psi1 + 81);
    const auto *psi1_83 = buffer.data(psi1 + 83);

    const auto *psk_0 = buffer.data(psk + 0);
    const auto *psk_1 = buffer.data(psk + 1);
    const auto *psk_2 = buffer.data(psk + 2);
    const auto *psk_3 = buffer.data(psk + 3);
    const auto *psk_5 = buffer.data(psk + 5);
    const auto *psk_6 = buffer.data(psk + 6);
    const auto *psk_8 = buffer.data(psk + 8);
    const auto *psk_9 = buffer.data(psk + 9);
    const auto *psk_10 = buffer.data(psk + 10);
    const auto *psk_12 = buffer.data(psk + 12);
    const auto *psk_13 = buffer.data(psk + 13);
    const auto *psk_14 = buffer.data(psk + 14);
    const auto *psk_15 = buffer.data(psk + 15);
    const auto *psk_17 = buffer.data(psk + 17);
    const auto *psk_18 = buffer.data(psk + 18);
    const auto *psk_19 = buffer.data(psk + 19);
    const auto *psk_20 = buffer.data(psk + 20);
    const auto *psk_21 = buffer.data(psk + 21);
    const auto *psk_27 = buffer.data(psk + 27);
    const auto *psk_28 = buffer.data(psk + 28);
    const auto *psk_30 = buffer.data(psk + 30);
    const auto *psk_31 = buffer.data(psk + 31);
    const auto *psk_32 = buffer.data(psk + 32);
    const auto *psk_33 = buffer.data(psk + 33);
    const auto *psk_35 = buffer.data(psk + 35);
    const auto *psk_36 = buffer.data(psk + 36);
    const auto *psk_37 = buffer.data(psk + 37);
    const auto *psk_39 = buffer.data(psk + 39);
    const auto *psk_42 = buffer.data(psk + 42);
    const auto *psk_44 = buffer.data(psk + 44);
    const auto *psk_46 = buffer.data(psk + 46);
    const auto *psk_48 = buffer.data(psk + 48);
    const auto *psk_49 = buffer.data(psk + 49);
    const auto *psk_51 = buffer.data(psk + 51);
    const auto *psk_53 = buffer.data(psk + 53);
    const auto *psk_54 = buffer.data(psk + 54);
    const auto *psk_55 = buffer.data(psk + 55);
    const auto *psk_57 = buffer.data(psk + 57);
    const auto *psk_59 = buffer.data(psk + 59);
    const auto *psk_60 = buffer.data(psk + 60);
    const auto *psk_61 = buffer.data(psk + 61);
    const auto *psk_62 = buffer.data(psk + 62);
    const auto *psk_64 = buffer.data(psk + 64);
    const auto *psk_65 = buffer.data(psk + 65);
    const auto *psk_66 = buffer.data(psk + 66);
    const auto *psk_67 = buffer.data(psk + 67);
    const auto *psk_68 = buffer.data(psk + 68);
    const auto *psk_69 = buffer.data(psk + 69);
    const auto *psk_70 = buffer.data(psk + 70);
    const auto *psk_71 = buffer.data(psk + 71);
    const auto *psk_72 = buffer.data(psk + 72);
    const auto *psk_74 = buffer.data(psk + 74);
    const auto *psk_77 = buffer.data(psk + 77);
    const auto *psk_79 = buffer.data(psk + 79);
    const auto *psk_81 = buffer.data(psk + 81);
    const auto *psk_83 = buffer.data(psk + 83);
    const auto *psk_84 = buffer.data(psk + 84);
    const auto *psk_86 = buffer.data(psk + 86);
    const auto *psk_88 = buffer.data(psk + 88);
    const auto *psk_89 = buffer.data(psk + 89);
    const auto *psk_90 = buffer.data(psk + 90);
    const auto *psk_92 = buffer.data(psk + 92);
    const auto *psk_94 = buffer.data(psk + 94);
    const auto *psk_95 = buffer.data(psk + 95);
    const auto *psk_96 = buffer.data(psk + 96);
    const auto *psk_97 = buffer.data(psk + 97);
    const auto *psk_99 = buffer.data(psk + 99);
    const auto *psk_100 = buffer.data(psk + 100);
    const auto *psk_101 = buffer.data(psk + 101);
    const auto *psk_102 = buffer.data(psk + 102);
    const auto *psk_103 = buffer.data(psk + 103);
    const auto *psk_104 = buffer.data(psk + 104);
    const auto *psk_105 = buffer.data(psk + 105);
    const auto *psk_106 = buffer.data(psk + 106);
    const auto *psk_107 = buffer.data(psk + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pc_x, pc_y, pc_z, ssl0_0, ssk_0, ssl1_0, \
                         psi0_0, psi1_0, psk_0, psk_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pa_x[k] * ssl0_0[k]
                 + f_0 * ssk_0[k]
                 - f_1 * pc_x[k] * ssl1_0[k];

        t_1[k] = f_2 * pc_y[k] * psk_0[k];

        t_2[k] = f_2 * pc_z[k] * psk_0[k];

        t_3[k] = f_3 * psi0_0[k]
                 - f_4 * psi1_0[k]
                 + f_2 * pc_y[k] * psk_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_y, pc_z, psi0_0, psi0_1, psi1_0, psi1_1, \
                         psk_2, psk_3, psk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * psk_2[k];

        t_5[k] = f_3 * psi0_0[k]
                 - f_4 * psi1_0[k]
                 + f_2 * pc_z[k] * psk_2[k];

        t_6[k] = f_5 * psi0_1[k]
                 - f_6 * psi1_1[k]
                 + f_2 * pc_y[k] * psk_3[k];

        t_7[k] = f_2 * pc_z[k] * psk_3[k];

        t_8[k] = f_2 * pc_y[k] * psk_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pc_y, pc_z, psi0_2, psi0_3, psi0_5, psi1_2, \
                         psi1_3, psi1_5, psk_5, psk_6, psk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * psi0_2[k]
                 - f_6 * psi1_2[k]
                 + f_2 * pc_z[k] * psk_5[k];

        t_10[k] = f_7 * psi0_3[k]
                  - f_8 * psi1_3[k]
                  + f_2 * pc_y[k] * psk_6[k];

        t_11[k] = f_2 * pc_z[k] * psk_6[k];

        t_12[k] = f_3 * psi0_5[k]
                  - f_4 * psi1_5[k]
                  + f_2 * pc_y[k] * psk_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pc_y, pc_z, psi0_5, psi0_6, psi0_8, \
                         psi1_5, psi1_6, psi1_8, psk_9, psk_10, \
                         psk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * pc_y[k] * psk_9[k];

        t_14[k] = f_7 * psi0_5[k]
                  - f_8 * psi1_5[k]
                  + f_2 * pc_z[k] * psk_9[k];

        t_15[k] = f_9 * psi0_6[k]
                  - f_10 * psi1_6[k]
                  + f_2 * pc_y[k] * psk_10[k];

        t_16[k] = f_2 * pc_z[k] * psk_10[k];

        t_17[k] = f_5 * psi0_8[k]
                  - f_6 * psi1_8[k]
                  + f_2 * pc_y[k] * psk_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_y, pc_z, psi0_9, psi0_10, psi1_9, \
                         psi1_10, psk_13, psk_14, psk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * psi0_9[k]
                  - f_4 * psi1_9[k]
                  + f_2 * pc_y[k] * psk_13[k];

        t_19[k] = f_2 * pc_y[k] * psk_14[k];

        t_20[k] = f_9 * psi0_9[k]
                  - f_10 * psi1_9[k]
                  + f_2 * pc_z[k] * psk_14[k];

        t_21[k] = f_11 * psi0_10[k]
                  - f_12 * psi1_10[k]
                  + f_2 * pc_y[k] * psk_15[k];

        t_22[k] = f_2 * pc_z[k] * psk_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pc_y, psi0_12, psi0_13, psi0_14, psi1_12, \
                         psi1_13, psi1_14, psk_17, psk_18, psk_19, \
                         psk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_7 * psi0_12[k]
                  - f_8 * psi1_12[k]
                  + f_2 * pc_y[k] * psk_17[k];

        t_24[k] = f_5 * psi0_13[k]
                  - f_6 * psi1_13[k]
                  + f_2 * pc_y[k] * psk_18[k];

        t_25[k] = f_3 * psi0_14[k]
                  - f_4 * psi1_14[k]
                  + f_2 * pc_y[k] * psk_19[k];

        t_26[k] = f_2 * pc_y[k] * psk_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_z, ssk_28, ssk_30, psi0_14, psi1_14, \
                         psk_20, psk_21, psk_28, psk_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_11 * psi0_14[k]
                  - f_12 * psi1_14[k]
                  + f_2 * pc_z[k] * psk_20[k];

        t_28[k] = f_13 * ssk_28[k]
                  + f_2 * pc_x[k] * psk_28[k];

        t_29[k] = f_2 * pc_z[k] * psk_21[k];

        t_30[k] = f_13 * ssk_30[k]
                  + f_2 * pc_x[k] * psk_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, ssk_31, ssk_32, ssk_33, \
                         ssk_35, psk_27, psk_31, psk_32, psk_33, \
                         psk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_13 * ssk_31[k]
                  + f_2 * pc_x[k] * psk_31[k];

        t_32[k] = f_13 * ssk_32[k]
                  + f_2 * pc_x[k] * psk_32[k];

        t_33[k] = f_13 * ssk_33[k]
                  + f_2 * pc_x[k] * psk_33[k];

        t_34[k] = f_2 * pc_y[k] * psk_27[k];

        t_35[k] = f_13 * ssk_35[k]
                  + f_2 * pc_x[k] * psk_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pc_x, pc_z, ssl0_36, ssl0_38, ssl0_39, \
                         ssl1_36, ssl1_38, ssl1_39, psk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * ssl0_36[k]
                  - f_1 * pc_x[k] * ssl1_36[k];

        t_37[k] = f_2 * pc_z[k] * psk_28[k];

        t_38[k] = pa_x[k] * ssl0_38[k]
                  - f_1 * pc_x[k] * ssl1_38[k];

        t_39[k] = pa_x[k] * ssl0_39[k]
                  - f_1 * pc_x[k] * ssl1_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pc_x, pc_y, ssl0_40, ssl0_41, ssl0_42, \
                         ssl1_40, ssl1_41, ssl1_42, psk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_x[k] * ssl0_40[k]
                  - f_1 * pc_x[k] * ssl1_40[k];

        t_41[k] = pa_x[k] * ssl0_41[k]
                  - f_1 * pc_x[k] * ssl1_41[k];

        t_42[k] = pa_x[k] * ssl0_42[k]
                  - f_1 * pc_x[k] * ssl1_42[k];

        t_43[k] = f_2 * pc_y[k] * psk_35[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_x, pa_y, pc_x, pc_y, ssl0_0, ssl0_44, ssl1_0, \
                         ssl1_44, psi0_29, psi1_29, psk_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_x[k] * ssl0_44[k]
                  - f_1 * pc_x[k] * ssl1_44[k];

        t_45[k] = pa_y[k] * ssl0_0[k]
                  - f_1 * pc_y[k] * ssl1_0[k];

        t_46[k] = f_14 * psi0_29[k]
                  - f_15 * psi1_29[k]
                  + f_2 * pc_x[k] * psk_37[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_x, pc_y, pc_z, ssl0_5, ssl1_5, \
                         psi0_31, psi1_31, psk_36, psk_37, psk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_2 * pc_z[k] * psk_36[k];

        t_48[k] = f_11 * psi0_31[k]
                  - f_12 * psi1_31[k]
                  + f_2 * pc_x[k] * psk_39[k];

        t_49[k] = f_2 * pc_z[k] * psk_37[k];

        t_50[k] = pa_y[k] * ssl0_5[k]
                  - f_1 * pc_y[k] * ssl1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pc_x, pc_z, psi0_34, psi0_36, psi1_34, psi1_36, \
                         psk_39, psk_42, psk_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_9 * psi0_34[k]
                  - f_10 * psi1_34[k]
                  + f_2 * pc_x[k] * psk_42[k];

        t_52[k] = f_2 * pc_z[k] * psk_39[k];

        t_53[k] = f_9 * psi0_36[k]
                  - f_10 * psi1_36[k]
                  + f_2 * pc_x[k] * psk_44[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pa_y, pc_x, pc_y, pc_z, ssl0_9, ssl1_9, psi0_38, \
                         psi1_38, psk_42, psk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_y[k] * ssl0_9[k]
                  - f_1 * pc_y[k] * ssl1_9[k];

        t_55[k] = f_7 * psi0_38[k]
                  - f_8 * psi1_38[k]
                  + f_2 * pc_x[k] * psk_46[k];

        t_56[k] = f_2 * pc_z[k] * psk_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_y, pc_x, pc_y, ssl0_14, ssl1_14, psi0_40, \
                         psi0_41, psi1_40, psi1_41, psk_48, psk_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_7 * psi0_40[k]
                  - f_8 * psi1_40[k]
                  + f_2 * pc_x[k] * psk_48[k];

        t_58[k] = f_7 * psi0_41[k]
                  - f_8 * psi1_41[k]
                  + f_2 * pc_x[k] * psk_49[k];

        t_59[k] = pa_y[k] * ssl0_14[k]
                  - f_1 * pc_y[k] * ssl1_14[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_z, psi0_43, psi0_45, psi0_46, \
                         psi1_43, psi1_45, psi1_46, psk_46, psk_51, psk_53, \
                         psk_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_5 * psi0_43[k]
                  - f_6 * psi1_43[k]
                  + f_2 * pc_x[k] * psk_51[k];

        t_61[k] = f_2 * pc_z[k] * psk_46[k];

        t_62[k] = f_5 * psi0_45[k]
                  - f_6 * psi1_45[k]
                  + f_2 * pc_x[k] * psk_53[k];

        t_63[k] = f_5 * psi0_46[k]
                  - f_6 * psi1_46[k]
                  + f_2 * pc_x[k] * psk_54[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_y, pc_x, pc_y, ssl0_20, ssl1_20, psi0_47, \
                         psi0_49, psi1_47, psi1_49, psk_55, psk_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * psi0_47[k]
                  - f_6 * psi1_47[k]
                  + f_2 * pc_x[k] * psk_55[k];

        t_65[k] = pa_y[k] * ssl0_20[k]
                  - f_1 * pc_y[k] * ssl1_20[k];

        t_66[k] = f_3 * psi0_49[k]
                  - f_4 * psi1_49[k]
                  + f_2 * pc_x[k] * psk_57[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pc_x, pc_z, psi0_51, psi0_52, psi0_53, \
                         psi1_51, psi1_52, psi1_53, psk_51, psk_59, psk_60, \
                         psk_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_2 * pc_z[k] * psk_51[k];

        t_68[k] = f_3 * psi0_51[k]
                  - f_4 * psi1_51[k]
                  + f_2 * pc_x[k] * psk_59[k];

        t_69[k] = f_3 * psi0_52[k]
                  - f_4 * psi1_52[k]
                  + f_2 * pc_x[k] * psk_60[k];

        t_70[k] = f_3 * psi0_53[k]
                  - f_4 * psi1_53[k]
                  + f_2 * pc_x[k] * psk_61[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, pa_y, pc_x, pc_y, ssl0_27, ssl1_27, \
                         psi0_54, psi1_54, psk_62, psk_64, psk_65, \
                         psk_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_3 * psi0_54[k]
                  - f_4 * psi1_54[k]
                  + f_2 * pc_x[k] * psk_62[k];

        t_72[k] = pa_y[k] * ssl0_27[k]
                  - f_1 * pc_y[k] * ssl1_27[k];

        t_73[k] = f_2 * pc_x[k] * psk_64[k];

        t_74[k] = f_2 * pc_x[k] * psk_65[k];

        t_75[k] = f_2 * pc_x[k] * psk_66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, pc_x, psk_67, psk_68, psk_69, psk_70, \
                         psk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * pc_x[k] * psk_67[k];

        t_77[k] = f_2 * pc_x[k] * psk_68[k];

        t_78[k] = f_2 * pc_x[k] * psk_69[k];

        t_79[k] = f_2 * pc_x[k] * psk_70[k];

        t_80[k] = f_2 * pc_x[k] * psk_71[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_y, pc_y, pc_z, ssl0_36, ssk_28, ssl1_36, \
                         psi0_49, psi1_49, psk_64, psk_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_y[k] * ssl0_36[k]
                  + f_0 * ssk_28[k]
                  - f_1 * pc_y[k] * ssl1_36[k];

        t_82[k] = f_2 * pc_z[k] * psk_64[k];

        t_83[k] = f_3 * psi0_49[k]
                  - f_4 * psi1_49[k]
                  + f_2 * pc_z[k] * psk_65[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_z, psi0_50, psi0_51, psi0_52, psi1_50, psi1_51, \
                         psi1_52, psk_66, psk_67, psk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_5 * psi0_50[k]
                  - f_6 * psi1_50[k]
                  + f_2 * pc_z[k] * psk_66[k];

        t_85[k] = f_7 * psi0_51[k]
                  - f_8 * psi1_51[k]
                  + f_2 * pc_z[k] * psk_67[k];

        t_86[k] = f_9 * psi0_52[k]
                  - f_10 * psi1_52[k]
                  + f_2 * pc_z[k] * psk_68[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_y, pc_y, pc_z, ssl0_44, ssk_35, ssl1_44, \
                         psi0_53, psi1_53, psk_69, psk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_11 * psi0_53[k]
                  - f_12 * psi1_53[k]
                  + f_2 * pc_z[k] * psk_69[k];

        t_88[k] = f_13 * ssk_35[k]
                  + f_2 * pc_y[k] * psk_71[k];

        t_89[k] = pa_y[k] * ssl0_44[k]
                  - f_1 * pc_y[k] * ssl1_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_z, pc_x, pc_y, pc_z, ssl0_0, ssl0_3, \
                         ssl1_0, ssl1_3, psi0_58, psi1_58, psk_72, \
                         psk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_z[k] * ssl0_0[k]
                  - f_1 * pc_z[k] * ssl1_0[k];

        t_91[k] = f_2 * pc_y[k] * psk_72[k];

        t_92[k] = f_14 * psi0_58[k]
                  - f_15 * psi1_58[k]
                  + f_2 * pc_x[k] * psk_74[k];

        t_93[k] = pa_z[k] * ssl0_3[k]
                  - f_1 * pc_z[k] * ssl1_3[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_z, pc_x, pc_y, pc_z, ssl0_6, ssl1_6, psi0_61, \
                         psi1_61, psk_74, psk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_2 * pc_y[k] * psk_74[k];

        t_95[k] = f_11 * psi0_61[k]
                  - f_12 * psi1_61[k]
                  + f_2 * pc_x[k] * psk_77[k];

        t_96[k] = pa_z[k] * ssl0_6[k]
                  - f_1 * pc_z[k] * ssl1_6[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pc_x, pc_y, psi0_63, psi0_65, psi1_63, psi1_65, \
                         psk_77, psk_79, psk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_9 * psi0_63[k]
                  - f_10 * psi1_63[k]
                  + f_2 * pc_x[k] * psk_79[k];

        t_98[k] = f_2 * pc_y[k] * psk_77[k];

        t_99[k] = f_9 * psi0_65[k]
                  - f_10 * psi1_65[k]
                  + f_2 * pc_x[k] * psk_81[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pc_x, pc_z, ssl0_10, ssl1_10, psi0_67, \
                         psi0_68, psi1_67, psi1_68, psk_83, psk_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * ssl0_10[k]
                   - f_1 * pc_z[k] * ssl1_10[k];

        t_101[k] = f_7 * psi0_67[k]
                   - f_8 * psi1_67[k]
                   + f_2 * pc_x[k] * psk_83[k];

        t_102[k] = f_7 * psi0_68[k]
                   - f_8 * psi1_68[k]
                   + f_2 * pc_x[k] * psk_84[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pa_z, pc_x, pc_y, pc_z, ssl0_15, ssl1_15, \
                         psi0_70, psi1_70, psk_81, psk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_2 * pc_y[k] * psk_81[k];

        t_104[k] = f_7 * psi0_70[k]
                   - f_8 * psi1_70[k]
                   + f_2 * pc_x[k] * psk_86[k];

        t_105[k] = pa_z[k] * ssl0_15[k]
                   - f_1 * pc_z[k] * ssl1_15[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_x, pc_y, psi0_72, psi0_73, psi0_74, \
                         psi1_72, psi1_73, psi1_74, psk_86, psk_88, psk_89, \
                         psk_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_5 * psi0_72[k]
                   - f_6 * psi1_72[k]
                   + f_2 * pc_x[k] * psk_88[k];

        t_107[k] = f_5 * psi0_73[k]
                   - f_6 * psi1_73[k]
                   + f_2 * pc_x[k] * psk_89[k];

        t_108[k] = f_5 * psi0_74[k]
                   - f_6 * psi1_74[k]
                   + f_2 * pc_x[k] * psk_90[k];

        t_109[k] = f_2 * pc_y[k] * psk_86[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pc_x, pc_z, ssl0_21, ssl1_21, psi0_76, \
                         psi0_78, psi1_76, psi1_78, psk_92, psk_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_5 * psi0_76[k]
                   - f_6 * psi1_76[k]
                   + f_2 * pc_x[k] * psk_92[k];

        t_111[k] = pa_z[k] * ssl0_21[k]
                   - f_1 * pc_z[k] * ssl1_21[k];

        t_112[k] = f_3 * psi0_78[k]
                   - f_4 * psi1_78[k]
                   + f_2 * pc_x[k] * psk_94[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pc_x, pc_y, psi0_79, psi0_80, psi0_81, \
                         psi1_79, psi1_80, psi1_81, psk_92, psk_95, psk_96, \
                         psk_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_3 * psi0_79[k]
                   - f_4 * psi1_79[k]
                   + f_2 * pc_x[k] * psk_95[k];

        t_114[k] = f_3 * psi0_80[k]
                   - f_4 * psi1_80[k]
                   + f_2 * pc_x[k] * psk_96[k];

        t_115[k] = f_3 * psi0_81[k]
                   - f_4 * psi1_81[k]
                   + f_2 * pc_x[k] * psk_97[k];

        t_116[k] = f_2 * pc_y[k] * psk_92[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, t_122, pc_x, psi0_83, psi1_83, \
                         psk_99, psk_100, psk_101, psk_102, psk_103, \
                         psk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_3 * psi0_83[k]
                   - f_4 * psi1_83[k]
                   + f_2 * pc_x[k] * psk_99[k];

        t_118[k] = f_2 * pc_x[k] * psk_100[k];

        t_119[k] = f_2 * pc_x[k] * psk_101[k];

        t_120[k] = f_2 * pc_x[k] * psk_102[k];

        t_121[k] = f_2 * pc_x[k] * psk_103[k];

        t_122[k] = f_2 * pc_x[k] * psk_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_z, pc_x, pc_z, ssl0_36, ssl1_36, \
                         psk_105, psk_106, psk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_2 * pc_x[k] * psk_105[k];

        t_124[k] = f_2 * pc_x[k] * psk_106[k];

        t_125[k] = f_2 * pc_x[k] * psk_107[k];

        t_126[k] = pa_z[k] * ssl0_36[k]
                   - f_1 * pc_z[k] * ssl1_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pc_y, psi0_78, psi0_79, psi0_80, psi1_78, \
                         psi1_79, psi1_80, psk_101, psk_102, psk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_14 * psi0_78[k]
                   - f_15 * psi1_78[k]
                   + f_2 * pc_y[k] * psk_101[k];

        t_128[k] = f_11 * psi0_79[k]
                   - f_12 * psi1_79[k]
                   + f_2 * pc_y[k] * psk_102[k];

        t_129[k] = f_9 * psi0_80[k]
                   - f_10 * psi1_80[k]
                   + f_2 * pc_y[k] * psk_103[k];
    }
}

static auto
compute_prim_psl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ssl0,
                                                          const size_t ssk, const size_t ssl1,
                                                          const size_t psi0, const size_t psi1,
                                                          const size_t psk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 0.5 / gamma;
    const auto f_4 = 0.5 * p / (gamma * q);
    const auto f_5 = 1.0 / gamma;
    const auto f_6 = p / (gamma * q);
    const auto f_7 = 1.5 / gamma;
    const auto f_8 = 1.5 * p / (gamma * q);

    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssl0_44 = buffer.data(ssl0 + 44);

    const auto *ssk_35 = buffer.data(ssk + 35);

    const auto *ssl1_44 = buffer.data(ssl1 + 44);

    const auto *psi0_81 = buffer.data(psi0 + 81);
    const auto *psi0_82 = buffer.data(psi0 + 82);
    const auto *psi0_83 = buffer.data(psi0 + 83);

    const auto *psi1_81 = buffer.data(psi1 + 81);
    const auto *psi1_82 = buffer.data(psi1 + 82);
    const auto *psi1_83 = buffer.data(psi1 + 83);

    const auto *psk_104 = buffer.data(psk + 104);
    const auto *psk_105 = buffer.data(psk + 105);
    const auto *psk_106 = buffer.data(psk + 106);
    const auto *psk_107 = buffer.data(psk + 107);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_y, psi0_81, psi0_82, psi0_83, psi1_81, \
                         psi1_82, psi1_83, psk_104, psk_105, psk_106, \
                         psk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * psi0_81[k]
                   - f_8 * psi1_81[k]
                   + f_2 * pc_y[k] * psk_104[k];

        t_131[k] = f_5 * psi0_82[k]
                   - f_6 * psi1_82[k]
                   + f_2 * pc_y[k] * psk_105[k];

        t_132[k] = f_3 * psi0_83[k]
                   - f_4 * psi1_83[k]
                   + f_2 * pc_y[k] * psk_106[k];

        t_133[k] = f_2 * pc_y[k] * psk_107[k];
    }

#pragma omp simd aligned(t_134, pa_z, pc_z, ssl0_44, ssk_35, ssl1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_z[k] * ssl0_44[k]
                   + f_0 * ssk_35[k]
                   - f_1 * pc_z[k] * ssl1_44[k];
    }
}

auto
compute_prim_psl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ssl0, const size_t ssk,
                                                   const size_t ssl1, const size_t psi0,
                                                   const size_t psi1, const size_t psk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_psl_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, ssl0, ssk,
                                                              ssl1, psi0, psi1, psk, ncols,
                                                              gamma, p, q);

    compute_prim_psl_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, ssl0, ssk,
                                                              ssl1, psi0, psi1, psk, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
