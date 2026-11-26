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


#include "SimdThreeCenterElectronRepulsionVrrRecGPI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t fpi0, const size_t fph,
                                                          const size_t fpi1, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t gpg0, const size_t gpg1,
                                                          const size_t gph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_14 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fpi0_0 = buffer.data(fpi0 + 0);
    const auto *fpi0_3 = buffer.data(fpi0 + 3);
    const auto *fpi0_5 = buffer.data(fpi0 + 5);
    const auto *fpi0_6 = buffer.data(fpi0 + 6);
    const auto *fpi0_9 = buffer.data(fpi0 + 9);
    const auto *fpi0_10 = buffer.data(fpi0 + 10);
    const auto *fpi0_14 = buffer.data(fpi0 + 14);
    const auto *fpi0_20 = buffer.data(fpi0 + 20);
    const auto *fpi0_27 = buffer.data(fpi0 + 27);

    const auto *fph_0 = buffer.data(fph + 0);
    const auto *fph_1 = buffer.data(fph + 1);
    const auto *fph_3 = buffer.data(fph + 3);
    const auto *fph_5 = buffer.data(fph + 5);
    const auto *fph_6 = buffer.data(fph + 6);
    const auto *fph_9 = buffer.data(fph + 9);
    const auto *fph_15 = buffer.data(fph + 15);
    const auto *fph_17 = buffer.data(fph + 17);
    const auto *fph_18 = buffer.data(fph + 18);
    const auto *fph_20 = buffer.data(fph + 20);
    const auto *fph_21 = buffer.data(fph + 21);
    const auto *fph_26 = buffer.data(fph + 26);
    const auto *fph_30 = buffer.data(fph + 30);
    const auto *fph_36 = buffer.data(fph + 36);
    const auto *fph_38 = buffer.data(fph + 38);
    const auto *fph_39 = buffer.data(fph + 39);
    const auto *fph_41 = buffer.data(fph + 41);
    const auto *fph_57 = buffer.data(fph + 57);
    const auto *fph_59 = buffer.data(fph + 59);
    const auto *fph_60 = buffer.data(fph + 60);
    const auto *fph_62 = buffer.data(fph + 62);
    const auto *fph_78 = buffer.data(fph + 78);
    const auto *fph_80 = buffer.data(fph + 80);
    const auto *fph_81 = buffer.data(fph + 81);
    const auto *fph_82 = buffer.data(fph + 82);
    const auto *fph_84 = buffer.data(fph + 84);
    const auto *fph_87 = buffer.data(fph + 87);
    const auto *fph_90 = buffer.data(fph + 90);
    const auto *fph_94 = buffer.data(fph + 94);

    const auto *fpi1_0 = buffer.data(fpi1 + 0);
    const auto *fpi1_3 = buffer.data(fpi1 + 3);
    const auto *fpi1_5 = buffer.data(fpi1 + 5);
    const auto *fpi1_6 = buffer.data(fpi1 + 6);
    const auto *fpi1_9 = buffer.data(fpi1 + 9);
    const auto *fpi1_10 = buffer.data(fpi1 + 10);
    const auto *fpi1_14 = buffer.data(fpi1 + 14);
    const auto *fpi1_20 = buffer.data(fpi1 + 20);
    const auto *fpi1_27 = buffer.data(fpi1 + 27);

    const auto *gsi0_0 = buffer.data(gsi0 + 0);
    const auto *gsi0_3 = buffer.data(gsi0 + 3);
    const auto *gsi0_5 = buffer.data(gsi0 + 5);
    const auto *gsi0_6 = buffer.data(gsi0 + 6);
    const auto *gsi0_9 = buffer.data(gsi0 + 9);
    const auto *gsi0_10 = buffer.data(gsi0 + 10);
    const auto *gsi0_12 = buffer.data(gsi0 + 12);
    const auto *gsi0_14 = buffer.data(gsi0 + 14);
    const auto *gsi0_21 = buffer.data(gsi0 + 21);
    const auto *gsi0_23 = buffer.data(gsi0 + 23);
    const auto *gsi0_24 = buffer.data(gsi0 + 24);
    const auto *gsi0_25 = buffer.data(gsi0 + 25);
    const auto *gsi0_27 = buffer.data(gsi0 + 27);

    const auto *gsh_0 = buffer.data(gsh + 0);
    const auto *gsh_1 = buffer.data(gsh + 1);
    const auto *gsh_2 = buffer.data(gsh + 2);
    const auto *gsh_3 = buffer.data(gsh + 3);
    const auto *gsh_5 = buffer.data(gsh + 5);
    const auto *gsh_6 = buffer.data(gsh + 6);
    const auto *gsh_8 = buffer.data(gsh + 8);
    const auto *gsh_9 = buffer.data(gsh + 9);
    const auto *gsh_10 = buffer.data(gsh + 10);
    const auto *gsh_14 = buffer.data(gsh + 14);
    const auto *gsh_15 = buffer.data(gsh + 15);
    const auto *gsh_17 = buffer.data(gsh + 17);
    const auto *gsh_18 = buffer.data(gsh + 18);
    const auto *gsh_19 = buffer.data(gsh + 19);
    const auto *gsh_20 = buffer.data(gsh + 20);
    const auto *gsh_21 = buffer.data(gsh + 21);
    const auto *gsh_26 = buffer.data(gsh + 26);
    const auto *gsh_30 = buffer.data(gsh + 30);
    const auto *gsh_36 = buffer.data(gsh + 36);
    const auto *gsh_38 = buffer.data(gsh + 38);
    const auto *gsh_39 = buffer.data(gsh + 39);
    const auto *gsh_40 = buffer.data(gsh + 40);

    const auto *gsi1_0 = buffer.data(gsi1 + 0);
    const auto *gsi1_3 = buffer.data(gsi1 + 3);
    const auto *gsi1_5 = buffer.data(gsi1 + 5);
    const auto *gsi1_6 = buffer.data(gsi1 + 6);
    const auto *gsi1_9 = buffer.data(gsi1 + 9);
    const auto *gsi1_10 = buffer.data(gsi1 + 10);
    const auto *gsi1_12 = buffer.data(gsi1 + 12);
    const auto *gsi1_14 = buffer.data(gsi1 + 14);
    const auto *gsi1_21 = buffer.data(gsi1 + 21);
    const auto *gsi1_23 = buffer.data(gsi1 + 23);
    const auto *gsi1_24 = buffer.data(gsi1 + 24);
    const auto *gsi1_25 = buffer.data(gsi1 + 25);
    const auto *gsi1_27 = buffer.data(gsi1 + 27);

    const auto *gpg0_0 = buffer.data(gpg0 + 0);
    const auto *gpg0_1 = buffer.data(gpg0 + 1);
    const auto *gpg0_2 = buffer.data(gpg0 + 2);
    const auto *gpg0_3 = buffer.data(gpg0 + 3);
    const auto *gpg0_5 = buffer.data(gpg0 + 5);
    const auto *gpg0_10 = buffer.data(gpg0 + 10);
    const auto *gpg0_12 = buffer.data(gpg0 + 12);
    const auto *gpg0_13 = buffer.data(gpg0 + 13);
    const auto *gpg0_14 = buffer.data(gpg0 + 14);
    const auto *gpg0_35 = buffer.data(gpg0 + 35);
    const auto *gpg0_42 = buffer.data(gpg0 + 42);
    const auto *gpg0_43 = buffer.data(gpg0 + 43);
    const auto *gpg0_44 = buffer.data(gpg0 + 44);
    const auto *gpg0_48 = buffer.data(gpg0 + 48);
    const auto *gpg0_55 = buffer.data(gpg0 + 55);
    const auto *gpg0_56 = buffer.data(gpg0 + 56);
    const auto *gpg0_57 = buffer.data(gpg0 + 57);
    const auto *gpg0_60 = buffer.data(gpg0 + 60);
    const auto *gpg0_62 = buffer.data(gpg0 + 62);
    const auto *gpg0_63 = buffer.data(gpg0 + 63);
    const auto *gpg0_65 = buffer.data(gpg0 + 65);
    const auto *gpg0_66 = buffer.data(gpg0 + 66);
    const auto *gpg0_70 = buffer.data(gpg0 + 70);

    const auto *gpg1_0 = buffer.data(gpg1 + 0);
    const auto *gpg1_1 = buffer.data(gpg1 + 1);
    const auto *gpg1_2 = buffer.data(gpg1 + 2);
    const auto *gpg1_3 = buffer.data(gpg1 + 3);
    const auto *gpg1_5 = buffer.data(gpg1 + 5);
    const auto *gpg1_10 = buffer.data(gpg1 + 10);
    const auto *gpg1_12 = buffer.data(gpg1 + 12);
    const auto *gpg1_13 = buffer.data(gpg1 + 13);
    const auto *gpg1_14 = buffer.data(gpg1 + 14);
    const auto *gpg1_35 = buffer.data(gpg1 + 35);
    const auto *gpg1_42 = buffer.data(gpg1 + 42);
    const auto *gpg1_43 = buffer.data(gpg1 + 43);
    const auto *gpg1_44 = buffer.data(gpg1 + 44);
    const auto *gpg1_48 = buffer.data(gpg1 + 48);
    const auto *gpg1_55 = buffer.data(gpg1 + 55);
    const auto *gpg1_56 = buffer.data(gpg1 + 56);
    const auto *gpg1_57 = buffer.data(gpg1 + 57);
    const auto *gpg1_60 = buffer.data(gpg1 + 60);
    const auto *gpg1_62 = buffer.data(gpg1 + 62);
    const auto *gpg1_63 = buffer.data(gpg1 + 63);
    const auto *gpg1_65 = buffer.data(gpg1 + 65);
    const auto *gpg1_66 = buffer.data(gpg1 + 66);
    const auto *gpg1_70 = buffer.data(gpg1 + 70);

    const auto *gph_0 = buffer.data(gph + 0);
    const auto *gph_1 = buffer.data(gph + 1);
    const auto *gph_2 = buffer.data(gph + 2);
    const auto *gph_3 = buffer.data(gph + 3);
    const auto *gph_5 = buffer.data(gph + 5);
    const auto *gph_6 = buffer.data(gph + 6);
    const auto *gph_8 = buffer.data(gph + 8);
    const auto *gph_9 = buffer.data(gph + 9);
    const auto *gph_10 = buffer.data(gph + 10);
    const auto *gph_14 = buffer.data(gph + 14);
    const auto *gph_15 = buffer.data(gph + 15);
    const auto *gph_17 = buffer.data(gph + 17);
    const auto *gph_18 = buffer.data(gph + 18);
    const auto *gph_19 = buffer.data(gph + 19);
    const auto *gph_20 = buffer.data(gph + 20);
    const auto *gph_21 = buffer.data(gph + 21);
    const auto *gph_23 = buffer.data(gph + 23);
    const auto *gph_24 = buffer.data(gph + 24);
    const auto *gph_26 = buffer.data(gph + 26);
    const auto *gph_27 = buffer.data(gph + 27);
    const auto *gph_30 = buffer.data(gph + 30);
    const auto *gph_31 = buffer.data(gph + 31);
    const auto *gph_35 = buffer.data(gph + 35);
    const auto *gph_36 = buffer.data(gph + 36);
    const auto *gph_38 = buffer.data(gph + 38);
    const auto *gph_39 = buffer.data(gph + 39);
    const auto *gph_41 = buffer.data(gph + 41);
    const auto *gph_42 = buffer.data(gph + 42);
    const auto *gph_44 = buffer.data(gph + 44);
    const auto *gph_45 = buffer.data(gph + 45);
    const auto *gph_47 = buffer.data(gph + 47);
    const auto *gph_48 = buffer.data(gph + 48);
    const auto *gph_50 = buffer.data(gph + 50);
    const auto *gph_51 = buffer.data(gph + 51);
    const auto *gph_52 = buffer.data(gph + 52);
    const auto *gph_56 = buffer.data(gph + 56);
    const auto *gph_57 = buffer.data(gph + 57);
    const auto *gph_59 = buffer.data(gph + 59);
    const auto *gph_60 = buffer.data(gph + 60);
    const auto *gph_61 = buffer.data(gph + 61);
    const auto *gph_62 = buffer.data(gph + 62);
    const auto *gph_63 = buffer.data(gph + 63);
    const auto *gph_64 = buffer.data(gph + 64);
    const auto *gph_66 = buffer.data(gph + 66);
    const auto *gph_68 = buffer.data(gph + 68);
    const auto *gph_69 = buffer.data(gph + 69);
    const auto *gph_70 = buffer.data(gph + 70);
    const auto *gph_72 = buffer.data(gph + 72);
    const auto *gph_73 = buffer.data(gph + 73);
    const auto *gph_78 = buffer.data(gph + 78);
    const auto *gph_79 = buffer.data(gph + 79);
    const auto *gph_80 = buffer.data(gph + 80);
    const auto *gph_81 = buffer.data(gph + 81);
    const auto *gph_82 = buffer.data(gph + 82);
    const auto *gph_83 = buffer.data(gph + 83);
    const auto *gph_84 = buffer.data(gph + 84);
    const auto *gph_85 = buffer.data(gph + 85);
    const auto *gph_86 = buffer.data(gph + 86);
    const auto *gph_87 = buffer.data(gph + 87);
    const auto *gph_89 = buffer.data(gph + 89);
    const auto *gph_90 = buffer.data(gph + 90);
    const auto *gph_91 = buffer.data(gph + 91);
    const auto *gph_93 = buffer.data(gph + 93);
    const auto *gph_94 = buffer.data(gph + 94);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fph_0, gsh_0, gpg0_0, \
                         gpg1_0, gph_0, gph_1, gph_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fph_0[k]
                 + f_1 * gsh_0[k]
                 + f_2 * gpg0_0[k]
                 - f_3 * gpg1_0[k]
                 + f_4 * pc_x[k] * gph_0[k];

        t_1[k] = f_4 * pc_y[k] * gph_0[k];

        t_2[k] = f_4 * pc_z[k] * gph_0[k];

        t_3[k] = f_5 * gpg0_0[k]
                 - f_6 * gpg1_0[k]
                 + f_4 * pc_y[k] * gph_1[k];

        t_4[k] = f_4 * pc_y[k] * gph_2[k];

        t_5[k] = f_5 * gpg0_0[k]
                 - f_6 * gpg1_0[k]
                 + f_4 * pc_z[k] * gph_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, gpg0_1, gpg0_2, gpg0_3, gpg1_1, \
                         gpg1_2, gpg1_3, gph_3, gph_5, gph_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * gpg0_1[k]
                 - f_8 * gpg1_1[k]
                 + f_4 * pc_y[k] * gph_3[k];

        t_7[k] = f_4 * pc_z[k] * gph_3[k];

        t_8[k] = f_4 * pc_y[k] * gph_5[k];

        t_9[k] = f_7 * gpg0_2[k]
                 - f_8 * gpg1_2[k]
                 + f_4 * pc_z[k] * gph_5[k];

        t_10[k] = f_9 * gpg0_3[k]
                  - f_10 * gpg1_3[k]
                  + f_4 * pc_y[k] * gph_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, fph_15, gsh_15, \
                         gpg0_5, gpg1_5, gph_6, gph_8, gph_9, gph_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_4 * pc_z[k] * gph_6[k];

        t_12[k] = f_5 * gpg0_5[k]
                  - f_6 * gpg1_5[k]
                  + f_4 * pc_y[k] * gph_8[k];

        t_13[k] = f_4 * pc_y[k] * gph_9[k];

        t_14[k] = f_9 * gpg0_5[k]
                  - f_10 * gpg1_5[k]
                  + f_4 * pc_z[k] * gph_9[k];

        t_15[k] = f_0 * fph_15[k]
                  + f_1 * gsh_15[k]
                  + f_4 * pc_x[k] * gph_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pc_x, pc_y, pc_z, fph_17, fph_18, gsh_17, \
                         gsh_18, gph_10, gph_14, gph_17, gph_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * pc_z[k] * gph_10[k];

        t_17[k] = f_0 * fph_17[k]
                  + f_1 * gsh_17[k]
                  + f_4 * pc_x[k] * gph_17[k];

        t_18[k] = f_0 * fph_18[k]
                  + f_1 * gsh_18[k]
                  + f_4 * pc_x[k] * gph_18[k];

        t_19[k] = f_4 * pc_y[k] * gph_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, fph_20, gsh_20, gpg0_10, \
                         gpg0_12, gpg1_10, gpg1_12, gph_15, gph_17, \
                         gph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * fph_20[k]
                  + f_1 * gsh_20[k]
                  + f_4 * pc_x[k] * gph_20[k];

        t_21[k] = f_2 * gpg0_10[k]
                  - f_3 * gpg1_10[k]
                  + f_4 * pc_y[k] * gph_15[k];

        t_22[k] = f_4 * pc_z[k] * gph_15[k];

        t_23[k] = f_9 * gpg0_12[k]
                  - f_10 * gpg1_12[k]
                  + f_4 * pc_y[k] * gph_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pc_y, pc_z, gpg0_13, gpg0_14, gpg1_13, \
                         gpg1_14, gph_18, gph_19, gph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * gpg0_13[k]
                  - f_8 * gpg1_13[k]
                  + f_4 * pc_y[k] * gph_18[k];

        t_25[k] = f_5 * gpg0_14[k]
                  - f_6 * gpg1_14[k]
                  + f_4 * pc_y[k] * gph_19[k];

        t_26[k] = f_4 * pc_y[k] * gph_20[k];

        t_27[k] = f_2 * gpg0_14[k]
                  - f_3 * gpg1_14[k]
                  + f_4 * pc_z[k] * gph_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_y, pc_z, gsi0_0, gsi0_3, gsh_0, \
                         gsh_1, gsi1_0, gsi1_3, gph_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * gsi0_0[k]
                  - f_11 * pc_y[k] * gsi1_0[k];

        t_29[k] = f_1 * gsh_0[k]
                  + f_4 * pc_y[k] * gph_21[k];

        t_30[k] = f_4 * pc_z[k] * gph_21[k];

        t_31[k] = pb_y[k] * gsi0_3[k]
                  + f_12 * gsh_1[k]
                  - f_11 * pc_y[k] * gsi1_3[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_y, pc_y, pc_z, gsi0_5, gsi0_6, gsh_2, \
                         gsh_3, gsi1_5, gsi1_6, gph_23, gph_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * gsh_2[k]
                  + f_4 * pc_y[k] * gph_23[k];

        t_33[k] = pb_y[k] * gsi0_5[k]
                  - f_11 * pc_y[k] * gsi1_5[k];

        t_34[k] = pb_y[k] * gsi0_6[k]
                  + f_13 * gsh_3[k]
                  - f_11 * pc_y[k] * gsi1_6[k];

        t_35[k] = f_4 * pc_z[k] * gph_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, gsi0_9, gsi0_10, gsh_5, \
                         gsh_6, gsi1_9, gsi1_10, gph_26, gph_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * gsh_5[k]
                  + f_4 * pc_y[k] * gph_26[k];

        t_37[k] = pb_y[k] * gsi0_9[k]
                  - f_11 * pc_y[k] * gsi1_9[k];

        t_38[k] = pb_y[k] * gsi0_10[k]
                  + f_0 * gsh_6[k]
                  - f_11 * pc_y[k] * gsi1_10[k];

        t_39[k] = f_4 * pc_z[k] * gph_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_x, pc_y, fph_36, gsi0_12, gsi0_14, \
                         gsh_8, gsh_9, gsi1_12, gsi1_14, gph_30, \
                         gph_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * gsi0_12[k]
                  + f_12 * gsh_8[k]
                  - f_11 * pc_y[k] * gsi1_12[k];

        t_41[k] = f_1 * gsh_9[k]
                  + f_4 * pc_y[k] * gph_30[k];

        t_42[k] = pb_y[k] * gsi0_14[k]
                  - f_11 * pc_y[k] * gsi1_14[k];

        t_43[k] = f_0 * fph_36[k]
                  + f_4 * pc_x[k] * gph_36[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pc_x, pc_y, pc_z, fph_38, fph_39, gsh_14, \
                         gph_31, gph_35, gph_38, gph_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_4 * pc_z[k] * gph_31[k];

        t_45[k] = f_0 * fph_38[k]
                  + f_4 * pc_x[k] * gph_38[k];

        t_46[k] = f_0 * fph_39[k]
                  + f_4 * pc_x[k] * gph_39[k];

        t_47[k] = f_1 * gsh_14[k]
                  + f_4 * pc_y[k] * gph_35[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_y, pc_x, pc_y, pc_z, fph_41, gsi0_21, gsh_15, \
                         gsi1_21, gph_36, gph_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * fph_41[k]
                  + f_4 * pc_x[k] * gph_41[k];

        t_49[k] = pb_y[k] * gsi0_21[k]
                  + f_14 * gsh_15[k]
                  - f_11 * pc_y[k] * gsi1_21[k];

        t_50[k] = f_4 * pc_z[k] * gph_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pc_y, gsi0_23, gsi0_24, gsi0_25, gsh_17, \
                         gsh_18, gsh_19, gsi1_23, gsi1_24, gsi1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_y[k] * gsi0_23[k]
                  + f_0 * gsh_17[k]
                  - f_11 * pc_y[k] * gsi1_23[k];

        t_52[k] = pb_y[k] * gsi0_24[k]
                  + f_13 * gsh_18[k]
                  - f_11 * pc_y[k] * gsi1_24[k];

        t_53[k] = pb_y[k] * gsi0_25[k]
                  + f_12 * gsh_19[k]
                  - f_11 * pc_y[k] * gsi1_25[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_y, pb_z, pc_y, pc_z, gsi0_0, gsi0_27, \
                         gsh_20, gsi1_0, gsi1_27, gph_41, gph_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * gsh_20[k]
                  + f_4 * pc_y[k] * gph_41[k];

        t_55[k] = pb_y[k] * gsi0_27[k]
                  - f_11 * pc_y[k] * gsi1_27[k];

        t_56[k] = pb_z[k] * gsi0_0[k]
                  - f_11 * pc_z[k] * gsi1_0[k];

        t_57[k] = f_4 * pc_y[k] * gph_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_z, pc_y, pc_z, gsi0_3, gsi0_5, gsh_0, \
                         gsh_2, gsi1_3, gsi1_5, gph_42, gph_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * gsh_0[k]
                  + f_4 * pc_z[k] * gph_42[k];

        t_59[k] = pb_z[k] * gsi0_3[k]
                  - f_11 * pc_z[k] * gsi1_3[k];

        t_60[k] = f_4 * pc_y[k] * gph_44[k];

        t_61[k] = pb_z[k] * gsi0_5[k]
                  + f_12 * gsh_2[k]
                  - f_11 * pc_z[k] * gsi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_z, pc_y, pc_z, gsi0_6, gsi0_9, gsh_3, \
                         gsh_5, gsi1_6, gsi1_9, gph_45, gph_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_z[k] * gsi0_6[k]
                  - f_11 * pc_z[k] * gsi1_6[k];

        t_63[k] = f_1 * gsh_3[k]
                  + f_4 * pc_z[k] * gph_45[k];

        t_64[k] = f_4 * pc_y[k] * gph_47[k];

        t_65[k] = pb_z[k] * gsi0_9[k]
                  + f_13 * gsh_5[k]
                  - f_11 * pc_z[k] * gsi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_z, pc_y, pc_z, gsi0_10, gsh_6, gsi1_10, \
                         gpg0_35, gpg1_35, gph_48, gph_50, gph_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_z[k] * gsi0_10[k]
                  - f_11 * pc_z[k] * gsi1_10[k];

        t_67[k] = f_1 * gsh_6[k]
                  + f_4 * pc_z[k] * gph_48[k];

        t_68[k] = f_5 * gpg0_35[k]
                  - f_6 * gpg1_35[k]
                  + f_4 * pc_y[k] * gph_50[k];

        t_69[k] = f_4 * pc_y[k] * gph_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_z, pc_x, pc_z, fph_57, fph_59, gsi0_14, \
                         gsh_9, gsh_10, gsi1_14, gph_52, gph_57, \
                         gph_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_z[k] * gsi0_14[k]
                  + f_0 * gsh_9[k]
                  - f_11 * pc_z[k] * gsi1_14[k];

        t_71[k] = f_0 * fph_57[k]
                  + f_4 * pc_x[k] * gph_57[k];

        t_72[k] = f_1 * gsh_10[k]
                  + f_4 * pc_z[k] * gph_52[k];

        t_73[k] = f_0 * fph_59[k]
                  + f_4 * pc_x[k] * gph_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_x, pc_y, pc_z, fph_60, fph_62, \
                         gsi0_21, gsi1_21, gph_56, gph_60, gph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * fph_60[k]
                  + f_4 * pc_x[k] * gph_60[k];

        t_75[k] = f_4 * pc_y[k] * gph_56[k];

        t_76[k] = f_0 * fph_62[k]
                  + f_4 * pc_x[k] * gph_62[k];

        t_77[k] = pb_z[k] * gsi0_21[k]
                  - f_11 * pc_z[k] * gsi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, gsh_15, gpg0_42, gpg0_43, gpg1_42, \
                         gpg1_43, gph_57, gph_59, gph_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * gsh_15[k]
                  + f_4 * pc_z[k] * gph_57[k];

        t_79[k] = f_9 * gpg0_42[k]
                  - f_10 * gpg1_42[k]
                  + f_4 * pc_y[k] * gph_59[k];

        t_80[k] = f_7 * gpg0_43[k]
                  - f_8 * gpg1_43[k]
                  + f_4 * pc_y[k] * gph_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_z, pc_y, pc_z, gsi0_27, gsh_20, gsi1_27, \
                         gpg0_44, gpg1_44, gph_61, gph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * gpg0_44[k]
                  - f_6 * gpg1_44[k]
                  + f_4 * pc_y[k] * gph_61[k];

        t_82[k] = f_4 * pc_y[k] * gph_62[k];

        t_83[k] = pb_z[k] * gsi0_27[k]
                  + f_14 * gsh_20[k]
                  - f_11 * pc_z[k] * gsi1_27[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_y, pc_y, pc_z, fpi0_0, fpi0_3, \
                         fph_0, fph_1, fpi1_0, fpi1_3, gph_63, gph_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pa_y[k] * fpi0_0[k]
                  - f_11 * pc_y[k] * fpi1_0[k];

        t_85[k] = f_1 * fph_0[k]
                  + f_4 * pc_y[k] * gph_63[k];

        t_86[k] = f_4 * pc_z[k] * gph_63[k];

        t_87[k] = pa_y[k] * fpi0_3[k]
                  + f_12 * fph_1[k]
                  - f_11 * pc_y[k] * fpi1_3[k];

        t_88[k] = f_4 * pc_z[k] * gph_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_y, pc_z, fpi0_5, fpi0_6, fph_3, \
                         fph_5, fpi1_5, fpi1_6, gph_66, gph_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * fpi0_5[k]
                  - f_11 * pc_y[k] * fpi1_5[k];

        t_90[k] = pa_y[k] * fpi0_6[k]
                  + f_13 * fph_3[k]
                  - f_11 * pc_y[k] * fpi1_6[k];

        t_91[k] = f_4 * pc_z[k] * gph_66[k];

        t_92[k] = f_1 * fph_5[k]
                  + f_4 * pc_y[k] * gph_68[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pc_y, pc_z, fpi0_9, fpi0_10, fph_6, \
                         fpi1_9, fpi1_10, gpg0_48, gpg1_48, gph_69, \
                         gph_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_y[k] * fpi0_9[k]
                  - f_11 * pc_y[k] * fpi1_9[k];

        t_94[k] = pa_y[k] * fpi0_10[k]
                  + f_0 * fph_6[k]
                  - f_11 * pc_y[k] * fpi1_10[k];

        t_95[k] = f_4 * pc_z[k] * gph_69[k];

        t_96[k] = f_5 * gpg0_48[k]
                  - f_6 * gpg1_48[k]
                  + f_4 * pc_z[k] * gph_70[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_y, pc_x, pc_y, pc_z, fpi0_14, fph_9, \
                         fph_78, fpi1_14, gsh_36, gph_72, gph_73, \
                         gph_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * fph_9[k]
                  + f_4 * pc_y[k] * gph_72[k];

        t_98[k] = pa_y[k] * fpi0_14[k]
                  - f_11 * pc_y[k] * fpi1_14[k];

        t_99[k] = f_13 * fph_78[k]
                  + f_1 * gsh_36[k]
                  + f_4 * pc_x[k] * gph_78[k];

        t_100[k] = f_4 * pc_z[k] * gph_73[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_x, fph_80, fph_81, fph_82, gsh_38, gsh_39, \
                         gsh_40, gph_80, gph_81, gph_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_13 * fph_80[k]
                   + f_1 * gsh_38[k]
                   + f_4 * pc_x[k] * gph_80[k];

        t_102[k] = f_13 * fph_81[k]
                   + f_1 * gsh_39[k]
                   + f_4 * pc_x[k] * gph_81[k];

        t_103[k] = f_13 * fph_82[k]
                   + f_1 * gsh_40[k]
                   + f_4 * pc_x[k] * gph_82[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pc_y, pc_z, fpi0_20, fph_15, \
                         fpi1_20, gpg0_55, gpg1_55, gph_78, gph_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * fpi0_20[k]
                   - f_11 * pc_y[k] * fpi1_20[k];

        t_105[k] = f_1 * fph_15[k]
                   + f_2 * gpg0_55[k]
                   - f_3 * gpg1_55[k]
                   + f_4 * pc_y[k] * gph_78[k];

        t_106[k] = f_4 * pc_z[k] * gph_78[k];

        t_107[k] = f_5 * gpg0_55[k]
                   - f_6 * gpg1_55[k]
                   + f_4 * pc_z[k] * gph_79[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_y, pc_z, fph_20, gpg0_56, gpg0_57, gpg1_56, \
                         gpg1_57, gph_80, gph_81, gph_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_7 * gpg0_56[k]
                   - f_8 * gpg1_56[k]
                   + f_4 * pc_z[k] * gph_80[k];

        t_109[k] = f_9 * gpg0_57[k]
                   - f_10 * gpg1_57[k]
                   + f_4 * pc_z[k] * gph_81[k];

        t_110[k] = f_1 * fph_20[k]
                   + f_4 * pc_y[k] * gph_83[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_y, pc_x, pc_y, pc_z, fpi0_27, fph_21, \
                         fph_84, fpi1_27, gsh_21, gpg0_60, gpg1_60, \
                         gph_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pa_y[k] * fpi0_27[k]
                   - f_11 * pc_y[k] * fpi1_27[k];

        t_112[k] = f_13 * fph_84[k]
                   + f_2 * gpg0_60[k]
                   - f_3 * gpg1_60[k]
                   + f_4 * pc_x[k] * gph_84[k];

        t_113[k] = f_1 * fph_21[k]
                   + f_1 * gsh_21[k]
                   + f_4 * pc_y[k] * gph_84[k];

        t_114[k] = f_4 * pc_z[k] * gph_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pc_x, pc_z, fph_87, gpg0_60, gpg0_63, gpg1_60, \
                         gpg1_63, gph_85, gph_86, gph_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_13 * fph_87[k]
                   + f_9 * gpg0_63[k]
                   - f_10 * gpg1_63[k]
                   + f_4 * pc_x[k] * gph_87[k];

        t_116[k] = f_4 * pc_z[k] * gph_85[k];

        t_117[k] = f_5 * gpg0_60[k]
                   - f_6 * gpg1_60[k]
                   + f_4 * pc_z[k] * gph_86[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_y, pc_z, fph_26, fph_90, gsh_26, \
                         gpg0_66, gpg1_66, gph_87, gph_89, gph_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_13 * fph_90[k]
                   + f_7 * gpg0_66[k]
                   - f_8 * gpg1_66[k]
                   + f_4 * pc_x[k] * gph_90[k];

        t_119[k] = f_4 * pc_z[k] * gph_87[k];

        t_120[k] = f_1 * fph_26[k]
                   + f_1 * gsh_26[k]
                   + f_4 * pc_y[k] * gph_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pc_x, pc_z, fph_94, gpg0_62, gpg0_70, gpg1_62, \
                         gpg1_70, gph_89, gph_90, gph_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_7 * gpg0_62[k]
                   - f_8 * gpg1_62[k]
                   + f_4 * pc_z[k] * gph_89[k];

        t_122[k] = f_13 * fph_94[k]
                   + f_5 * gpg0_70[k]
                   - f_6 * gpg1_70[k]
                   + f_4 * pc_x[k] * gph_94[k];

        t_123[k] = f_4 * pc_z[k] * gph_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pc_y, pc_z, fph_30, gsh_30, gpg0_63, gpg0_65, \
                         gpg1_63, gpg1_65, gph_91, gph_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * gpg0_63[k]
                   - f_6 * gpg1_63[k]
                   + f_4 * pc_z[k] * gph_91[k];

        t_125[k] = f_1 * fph_30[k]
                   + f_1 * gsh_30[k]
                   + f_4 * pc_y[k] * gph_93[k];

        t_126[k] = f_9 * gpg0_65[k]
                   - f_10 * gpg1_65[k]
                   + f_4 * pc_z[k] * gph_93[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dpi1,
                                                          const size_t fpi0, const size_t fph,
                                                          const size_t fpi1, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t gpg0, const size_t gpg1,
                                                          const size_t gph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_15 = 1.0 / p;
    const auto f_16 = gamma / (p * q);
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpi0_133 = buffer.data(dpi0 + 133);

    const auto *dpi1_133 = buffer.data(dpi1 + 133);

    const auto *fpi0_0 = buffer.data(fpi0 + 0);
    const auto *fpi0_3 = buffer.data(fpi0 + 3);
    const auto *fpi0_5 = buffer.data(fpi0 + 5);
    const auto *fpi0_6 = buffer.data(fpi0 + 6);
    const auto *fpi0_9 = buffer.data(fpi0 + 9);
    const auto *fpi0_10 = buffer.data(fpi0 + 10);
    const auto *fpi0_14 = buffer.data(fpi0 + 14);
    const auto *fpi0_15 = buffer.data(fpi0 + 15);
    const auto *fpi0_21 = buffer.data(fpi0 + 21);
    const auto *fpi0_28 = buffer.data(fpi0 + 28);
    const auto *fpi0_31 = buffer.data(fpi0 + 31);
    const auto *fpi0_34 = buffer.data(fpi0 + 34);
    const auto *fpi0_35 = buffer.data(fpi0 + 35);
    const auto *fpi0_38 = buffer.data(fpi0 + 38);
    const auto *fpi0_39 = buffer.data(fpi0 + 39);
    const auto *fpi0_40 = buffer.data(fpi0 + 40);
    const auto *fpi0_49 = buffer.data(fpi0 + 49);
    const auto *fpi0_56 = buffer.data(fpi0 + 56);
    const auto *fpi0_61 = buffer.data(fpi0 + 61);
    const auto *fpi0_65 = buffer.data(fpi0 + 65);
    const auto *fpi0_68 = buffer.data(fpi0 + 68);
    const auto *fpi0_70 = buffer.data(fpi0 + 70);
    const auto *fpi0_83 = buffer.data(fpi0 + 83);
    const auto *fpi0_133 = buffer.data(fpi0 + 133);

    const auto *fph_0 = buffer.data(fph + 0);
    const auto *fph_2 = buffer.data(fph + 2);
    const auto *fph_5 = buffer.data(fph + 5);
    const auto *fph_9 = buffer.data(fph + 9);
    const auto *fph_20 = buffer.data(fph + 20);
    const auto *fph_21 = buffer.data(fph + 21);
    const auto *fph_24 = buffer.data(fph + 24);
    const auto *fph_27 = buffer.data(fph + 27);
    const auto *fph_28 = buffer.data(fph + 28);
    const auto *fph_41 = buffer.data(fph + 41);
    const auto *fph_42 = buffer.data(fph + 42);
    const auto *fph_47 = buffer.data(fph + 47);
    const auto *fph_50 = buffer.data(fph + 50);
    const auto *fph_51 = buffer.data(fph + 51);
    const auto *fph_62 = buffer.data(fph + 62);
    const auto *fph_99 = buffer.data(fph + 99);
    const auto *fph_101 = buffer.data(fph + 101);
    const auto *fph_102 = buffer.data(fph + 102);
    const auto *fph_103 = buffer.data(fph + 103);
    const auto *fph_104 = buffer.data(fph + 104);
    const auto *fph_120 = buffer.data(fph + 120);
    const auto *fph_122 = buffer.data(fph + 122);
    const auto *fph_123 = buffer.data(fph + 123);
    const auto *fph_124 = buffer.data(fph + 124);
    const auto *fph_125 = buffer.data(fph + 125);
    const auto *fph_142 = buffer.data(fph + 142);
    const auto *fph_143 = buffer.data(fph + 143);
    const auto *fph_144 = buffer.data(fph + 144);
    const auto *fph_146 = buffer.data(fph + 146);
    const auto *fph_162 = buffer.data(fph + 162);
    const auto *fph_163 = buffer.data(fph + 163);
    const auto *fph_164 = buffer.data(fph + 164);
    const auto *fph_165 = buffer.data(fph + 165);
    const auto *fph_167 = buffer.data(fph + 167);
    const auto *fph_168 = buffer.data(fph + 168);
    const auto *fph_173 = buffer.data(fph + 173);
    const auto *fph_177 = buffer.data(fph + 177);
    const auto *fph_182 = buffer.data(fph + 182);
    const auto *fph_183 = buffer.data(fph + 183);
    const auto *fph_184 = buffer.data(fph + 184);
    const auto *fph_185 = buffer.data(fph + 185);
    const auto *fph_186 = buffer.data(fph + 186);
    const auto *fph_188 = buffer.data(fph + 188);

    const auto *fpi1_0 = buffer.data(fpi1 + 0);
    const auto *fpi1_3 = buffer.data(fpi1 + 3);
    const auto *fpi1_5 = buffer.data(fpi1 + 5);
    const auto *fpi1_6 = buffer.data(fpi1 + 6);
    const auto *fpi1_9 = buffer.data(fpi1 + 9);
    const auto *fpi1_10 = buffer.data(fpi1 + 10);
    const auto *fpi1_14 = buffer.data(fpi1 + 14);
    const auto *fpi1_15 = buffer.data(fpi1 + 15);
    const auto *fpi1_21 = buffer.data(fpi1 + 21);
    const auto *fpi1_28 = buffer.data(fpi1 + 28);
    const auto *fpi1_31 = buffer.data(fpi1 + 31);
    const auto *fpi1_34 = buffer.data(fpi1 + 34);
    const auto *fpi1_35 = buffer.data(fpi1 + 35);
    const auto *fpi1_38 = buffer.data(fpi1 + 38);
    const auto *fpi1_39 = buffer.data(fpi1 + 39);
    const auto *fpi1_40 = buffer.data(fpi1 + 40);
    const auto *fpi1_49 = buffer.data(fpi1 + 49);
    const auto *fpi1_56 = buffer.data(fpi1 + 56);
    const auto *fpi1_61 = buffer.data(fpi1 + 61);
    const auto *fpi1_65 = buffer.data(fpi1 + 65);
    const auto *fpi1_68 = buffer.data(fpi1 + 68);
    const auto *fpi1_70 = buffer.data(fpi1 + 70);
    const auto *fpi1_83 = buffer.data(fpi1 + 83);
    const auto *fpi1_133 = buffer.data(fpi1 + 133);

    const auto *gsi0_31 = buffer.data(gsi0 + 31);
    const auto *gsi0_34 = buffer.data(gsi0 + 34);
    const auto *gsi0_38 = buffer.data(gsi0 + 38);
    const auto *gsi0_49 = buffer.data(gsi0 + 49);
    const auto *gsi0_51 = buffer.data(gsi0 + 51);
    const auto *gsi0_52 = buffer.data(gsi0 + 52);
    const auto *gsi0_53 = buffer.data(gsi0 + 53);
    const auto *gsi0_61 = buffer.data(gsi0 + 61);
    const auto *gsi0_65 = buffer.data(gsi0 + 65);
    const auto *gsi0_70 = buffer.data(gsi0 + 70);
    const auto *gsi0_78 = buffer.data(gsi0 + 78);
    const auto *gsi0_79 = buffer.data(gsi0 + 79);
    const auto *gsi0_80 = buffer.data(gsi0 + 80);
    const auto *gsi0_81 = buffer.data(gsi0 + 81);
    const auto *gsi0_83 = buffer.data(gsi0 + 83);

    const auto *gsh_21 = buffer.data(gsh + 21);
    const auto *gsh_22 = buffer.data(gsh + 22);
    const auto *gsh_24 = buffer.data(gsh + 24);
    const auto *gsh_27 = buffer.data(gsh + 27);
    const auto *gsh_31 = buffer.data(gsh + 31);
    const auto *gsh_36 = buffer.data(gsh + 36);
    const auto *gsh_37 = buffer.data(gsh + 37);
    const auto *gsh_38 = buffer.data(gsh + 38);
    const auto *gsh_39 = buffer.data(gsh + 39);
    const auto *gsh_41 = buffer.data(gsh + 41);
    const auto *gsh_42 = buffer.data(gsh + 42);
    const auto *gsh_44 = buffer.data(gsh + 44);
    const auto *gsh_47 = buffer.data(gsh + 47);
    const auto *gsh_51 = buffer.data(gsh + 51);
    const auto *gsh_56 = buffer.data(gsh + 56);
    const auto *gsh_58 = buffer.data(gsh + 58);
    const auto *gsh_59 = buffer.data(gsh + 59);
    const auto *gsh_60 = buffer.data(gsh + 60);
    const auto *gsh_61 = buffer.data(gsh + 61);
    const auto *gsh_62 = buffer.data(gsh + 62);

    const auto *gsi1_31 = buffer.data(gsi1 + 31);
    const auto *gsi1_34 = buffer.data(gsi1 + 34);
    const auto *gsi1_38 = buffer.data(gsi1 + 38);
    const auto *gsi1_49 = buffer.data(gsi1 + 49);
    const auto *gsi1_51 = buffer.data(gsi1 + 51);
    const auto *gsi1_52 = buffer.data(gsi1 + 52);
    const auto *gsi1_53 = buffer.data(gsi1 + 53);
    const auto *gsi1_61 = buffer.data(gsi1 + 61);
    const auto *gsi1_65 = buffer.data(gsi1 + 65);
    const auto *gsi1_70 = buffer.data(gsi1 + 70);
    const auto *gsi1_78 = buffer.data(gsi1 + 78);
    const auto *gsi1_79 = buffer.data(gsi1 + 79);
    const auto *gsi1_80 = buffer.data(gsi1 + 80);
    const auto *gsi1_81 = buffer.data(gsi1 + 81);
    const auto *gsi1_83 = buffer.data(gsi1 + 83);

    const auto *gpg0_70 = buffer.data(gpg0 + 70);
    const auto *gpg0_71 = buffer.data(gpg0 + 71);
    const auto *gpg0_72 = buffer.data(gpg0 + 72);
    const auto *gpg0_74 = buffer.data(gpg0 + 74);
    const auto *gpg0_92 = buffer.data(gpg0 + 92);
    const auto *gpg0_94 = buffer.data(gpg0 + 94);
    const auto *gpg0_95 = buffer.data(gpg0 + 95);
    const auto *gpg0_101 = buffer.data(gpg0 + 101);
    const auto *gpg0_102 = buffer.data(gpg0 + 102);
    const auto *gpg0_103 = buffer.data(gpg0 + 103);
    const auto *gpg0_104 = buffer.data(gpg0 + 104);
    const auto *gpg0_120 = buffer.data(gpg0 + 120);
    const auto *gpg0_121 = buffer.data(gpg0 + 121);
    const auto *gpg0_122 = buffer.data(gpg0 + 122);
    const auto *gpg0_123 = buffer.data(gpg0 + 123);
    const auto *gpg0_124 = buffer.data(gpg0 + 124);
    const auto *gpg0_125 = buffer.data(gpg0 + 125);
    const auto *gpg0_129 = buffer.data(gpg0 + 129);
    const auto *gpg0_130 = buffer.data(gpg0 + 130);
    const auto *gpg0_131 = buffer.data(gpg0 + 131);
    const auto *gpg0_132 = buffer.data(gpg0 + 132);
    const auto *gpg0_134 = buffer.data(gpg0 + 134);

    const auto *gpg1_70 = buffer.data(gpg1 + 70);
    const auto *gpg1_71 = buffer.data(gpg1 + 71);
    const auto *gpg1_72 = buffer.data(gpg1 + 72);
    const auto *gpg1_74 = buffer.data(gpg1 + 74);
    const auto *gpg1_92 = buffer.data(gpg1 + 92);
    const auto *gpg1_94 = buffer.data(gpg1 + 94);
    const auto *gpg1_95 = buffer.data(gpg1 + 95);
    const auto *gpg1_101 = buffer.data(gpg1 + 101);
    const auto *gpg1_102 = buffer.data(gpg1 + 102);
    const auto *gpg1_103 = buffer.data(gpg1 + 103);
    const auto *gpg1_104 = buffer.data(gpg1 + 104);
    const auto *gpg1_120 = buffer.data(gpg1 + 120);
    const auto *gpg1_121 = buffer.data(gpg1 + 121);
    const auto *gpg1_122 = buffer.data(gpg1 + 122);
    const auto *gpg1_123 = buffer.data(gpg1 + 123);
    const auto *gpg1_124 = buffer.data(gpg1 + 124);
    const auto *gpg1_125 = buffer.data(gpg1 + 125);
    const auto *gpg1_129 = buffer.data(gpg1 + 129);
    const auto *gpg1_130 = buffer.data(gpg1 + 130);
    const auto *gpg1_131 = buffer.data(gpg1 + 131);
    const auto *gpg1_132 = buffer.data(gpg1 + 132);
    const auto *gpg1_134 = buffer.data(gpg1 + 134);

    const auto *gph_94 = buffer.data(gph + 94);
    const auto *gph_99 = buffer.data(gph + 99);
    const auto *gph_100 = buffer.data(gph + 100);
    const auto *gph_101 = buffer.data(gph + 101);
    const auto *gph_102 = buffer.data(gph + 102);
    const auto *gph_103 = buffer.data(gph + 103);
    const auto *gph_104 = buffer.data(gph + 104);
    const auto *gph_105 = buffer.data(gph + 105);
    const auto *gph_106 = buffer.data(gph + 106);
    const auto *gph_108 = buffer.data(gph + 108);
    const auto *gph_110 = buffer.data(gph + 110);
    const auto *gph_111 = buffer.data(gph + 111);
    const auto *gph_114 = buffer.data(gph + 114);
    const auto *gph_115 = buffer.data(gph + 115);
    const auto *gph_120 = buffer.data(gph + 120);
    const auto *gph_122 = buffer.data(gph + 122);
    const auto *gph_123 = buffer.data(gph + 123);
    const auto *gph_124 = buffer.data(gph + 124);
    const auto *gph_125 = buffer.data(gph + 125);
    const auto *gph_126 = buffer.data(gph + 126);
    const auto *gph_128 = buffer.data(gph + 128);
    const auto *gph_130 = buffer.data(gph + 130);
    const auto *gph_131 = buffer.data(gph + 131);
    const auto *gph_133 = buffer.data(gph + 133);
    const auto *gph_134 = buffer.data(gph + 134);
    const auto *gph_135 = buffer.data(gph + 135);
    const auto *gph_140 = buffer.data(gph + 140);
    const auto *gph_142 = buffer.data(gph + 142);
    const auto *gph_143 = buffer.data(gph + 143);
    const auto *gph_144 = buffer.data(gph + 144);
    const auto *gph_145 = buffer.data(gph + 145);
    const auto *gph_146 = buffer.data(gph + 146);
    const auto *gph_147 = buffer.data(gph + 147);
    const auto *gph_149 = buffer.data(gph + 149);
    const auto *gph_152 = buffer.data(gph + 152);
    const auto *gph_156 = buffer.data(gph + 156);
    const auto *gph_161 = buffer.data(gph + 161);
    const auto *gph_162 = buffer.data(gph + 162);
    const auto *gph_163 = buffer.data(gph + 163);
    const auto *gph_164 = buffer.data(gph + 164);
    const auto *gph_165 = buffer.data(gph + 165);
    const auto *gph_167 = buffer.data(gph + 167);
    const auto *gph_168 = buffer.data(gph + 168);
    const auto *gph_169 = buffer.data(gph + 169);
    const auto *gph_170 = buffer.data(gph + 170);
    const auto *gph_171 = buffer.data(gph + 171);
    const auto *gph_172 = buffer.data(gph + 172);
    const auto *gph_173 = buffer.data(gph + 173);
    const auto *gph_174 = buffer.data(gph + 174);
    const auto *gph_175 = buffer.data(gph + 175);
    const auto *gph_176 = buffer.data(gph + 176);
    const auto *gph_177 = buffer.data(gph + 177);
    const auto *gph_182 = buffer.data(gph + 182);
    const auto *gph_183 = buffer.data(gph + 183);
    const auto *gph_184 = buffer.data(gph + 184);
    const auto *gph_185 = buffer.data(gph + 185);
    const auto *gph_186 = buffer.data(gph + 186);
    const auto *gph_188 = buffer.data(gph + 188);

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pc_x, pc_z, fph_99, fph_101, \
                         fph_102, fph_103, gph_94, gph_99, gph_101, gph_102, \
                         gph_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_13 * fph_99[k]
                   + f_4 * pc_x[k] * gph_99[k];

        t_128[k] = f_4 * pc_z[k] * gph_94[k];

        t_129[k] = f_13 * fph_101[k]
                   + f_4 * pc_x[k] * gph_101[k];

        t_130[k] = f_13 * fph_102[k]
                   + f_4 * pc_x[k] * gph_102[k];

        t_131[k] = f_13 * fph_103[k]
                   + f_4 * pc_x[k] * gph_103[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pa_x, pc_x, pc_z, dpi0_133, dpi1_133, fpi0_133, \
                         fph_104, fpi1_133, gph_99, gph_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_13 * fph_104[k]
                   + f_4 * pc_x[k] * gph_104[k];

        t_133[k] = f_15 * dpi0_133[k]
                   - f_16 * dpi1_133[k]
                   + pa_x[k] * fpi0_133[k]
                   - f_11 * pc_x[k] * fpi1_133[k];

        t_134[k] = f_4 * pc_z[k] * gph_99[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pc_z, gpg0_70, gpg0_71, gpg0_72, gpg1_70, \
                         gpg1_71, gpg1_72, gph_100, gph_101, gph_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_5 * gpg0_70[k]
                   - f_6 * gpg1_70[k]
                   + f_4 * pc_z[k] * gph_100[k];

        t_136[k] = f_7 * gpg0_71[k]
                   - f_8 * gpg1_71[k]
                   + f_4 * pc_z[k] * gph_101[k];

        t_137[k] = f_9 * gpg0_72[k]
                   - f_10 * gpg1_72[k]
                   + f_4 * pc_z[k] * gph_102[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_y, pc_y, pc_z, fpi0_56, fph_41, \
                         fph_42, fpi1_56, gsh_41, gpg0_74, gpg1_74, gph_104, \
                         gph_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_1 * fph_41[k]
                   + f_1 * gsh_41[k]
                   + f_4 * pc_y[k] * gph_104[k];

        t_139[k] = f_2 * gpg0_74[k]
                   - f_3 * gpg1_74[k]
                   + f_4 * pc_z[k] * gph_104[k];

        t_140[k] = pa_y[k] * fpi0_56[k]
                   - f_11 * pc_y[k] * fpi1_56[k];

        t_141[k] = f_1 * fph_42[k]
                   + f_4 * pc_y[k] * gph_105[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pa_y, pb_z, pc_y, pc_z, fpi0_61, fpi1_61, \
                         gsi0_31, gsh_21, gsh_22, gsi1_31, gph_105, \
                         gph_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_1 * gsh_21[k]
                   + f_4 * pc_z[k] * gph_105[k];

        t_143[k] = pb_z[k] * gsi0_31[k]
                   - f_11 * pc_z[k] * gsi1_31[k];

        t_144[k] = f_1 * gsh_22[k]
                   + f_4 * pc_z[k] * gph_106[k];

        t_145[k] = pa_y[k] * fpi0_61[k]
                   - f_11 * pc_y[k] * fpi1_61[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_y, pb_z, pc_y, pc_z, fpi0_65, fph_47, \
                         fpi1_65, gsi0_34, gsh_24, gsi1_34, gph_108, \
                         gph_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pb_z[k] * gsi0_34[k]
                   - f_11 * pc_z[k] * gsi1_34[k];

        t_147[k] = f_1 * gsh_24[k]
                   + f_4 * pc_z[k] * gph_108[k];

        t_148[k] = f_1 * fph_47[k]
                   + f_4 * pc_y[k] * gph_110[k];

        t_149[k] = pa_y[k] * fpi0_65[k]
                   - f_11 * pc_y[k] * fpi1_65[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_y, pb_z, pc_y, pc_z, fpi0_68, fph_50, \
                         fpi1_68, gsi0_38, gsh_27, gsi1_38, gph_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_z[k] * gsi0_38[k]
                   - f_11 * pc_z[k] * gsi1_38[k];

        t_151[k] = f_1 * gsh_27[k]
                   + f_4 * pc_z[k] * gph_111[k];

        t_152[k] = pa_y[k] * fpi0_68[k]
                   + f_12 * fph_50[k]
                   - f_11 * pc_y[k] * fpi1_68[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_y, pc_x, pc_y, pc_z, fpi0_70, fph_51, \
                         fph_120, fpi1_70, gsh_31, gph_114, gph_115, \
                         gph_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_1 * fph_51[k]
                   + f_4 * pc_y[k] * gph_114[k];

        t_154[k] = pa_y[k] * fpi0_70[k]
                   - f_11 * pc_y[k] * fpi1_70[k];

        t_155[k] = f_13 * fph_120[k]
                   + f_4 * pc_x[k] * gph_120[k];

        t_156[k] = f_1 * gsh_31[k]
                   + f_4 * pc_z[k] * gph_115[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, fph_122, fph_123, fph_124, fph_125, \
                         gph_122, gph_123, gph_124, gph_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_13 * fph_122[k]
                   + f_4 * pc_x[k] * gph_122[k];

        t_158[k] = f_13 * fph_123[k]
                   + f_4 * pc_x[k] * gph_123[k];

        t_159[k] = f_13 * fph_124[k]
                   + f_4 * pc_x[k] * gph_124[k];

        t_160[k] = f_13 * fph_125[k]
                   + f_4 * pc_x[k] * gph_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_z, pc_z, gsi0_49, gsi0_51, gsi0_52, \
                         gsh_36, gsh_37, gsh_38, gsi1_49, gsi1_51, gsi1_52, \
                         gph_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_z[k] * gsi0_49[k]
                   - f_11 * pc_z[k] * gsi1_49[k];

        t_162[k] = f_1 * gsh_36[k]
                   + f_4 * pc_z[k] * gph_120[k];

        t_163[k] = pb_z[k] * gsi0_51[k]
                   + f_12 * gsh_37[k]
                   - f_11 * pc_z[k] * gsi1_51[k];

        t_164[k] = pb_z[k] * gsi0_52[k]
                   + f_13 * gsh_38[k]
                   - f_11 * pc_z[k] * gsi1_52[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_y, pb_z, pc_y, pc_z, fpi0_83, fph_62, \
                         fpi1_83, gsi0_53, gsh_39, gsi1_53, gph_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pb_z[k] * gsi0_53[k]
                   + f_0 * gsh_39[k]
                   - f_11 * pc_z[k] * gsi1_53[k];

        t_166[k] = f_1 * fph_62[k]
                   + f_4 * pc_y[k] * gph_125[k];

        t_167[k] = pa_y[k] * fpi0_83[k]
                   - f_11 * pc_y[k] * fpi1_83[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pa_z, pc_y, pc_z, fpi0_0, fpi0_3, \
                         fph_0, fpi1_0, fpi1_3, gph_126, gph_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * fpi0_0[k]
                   - f_11 * pc_z[k] * fpi1_0[k];

        t_169[k] = f_4 * pc_y[k] * gph_126[k];

        t_170[k] = f_1 * fph_0[k]
                   + f_4 * pc_z[k] * gph_126[k];

        t_171[k] = pa_z[k] * fpi0_3[k]
                   - f_11 * pc_z[k] * fpi1_3[k];

        t_172[k] = f_4 * pc_y[k] * gph_128[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_z, pc_y, pc_z, fpi0_5, fpi0_6, fph_2, \
                         fpi1_5, fpi1_6, gpg0_92, gpg1_92, gph_130, \
                         gph_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = pa_z[k] * fpi0_5[k]
                   + f_12 * fph_2[k]
                   - f_11 * pc_z[k] * fpi1_5[k];

        t_174[k] = pa_z[k] * fpi0_6[k]
                   - f_11 * pc_z[k] * fpi1_6[k];

        t_175[k] = f_5 * gpg0_92[k]
                   - f_6 * gpg1_92[k]
                   + f_4 * pc_y[k] * gph_130[k];

        t_176[k] = f_4 * pc_y[k] * gph_131[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_z, pc_y, pc_z, fpi0_9, fpi0_10, fph_5, \
                         fpi1_9, fpi1_10, gpg0_94, gpg1_94, gph_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pa_z[k] * fpi0_9[k]
                   + f_13 * fph_5[k]
                   - f_11 * pc_z[k] * fpi1_9[k];

        t_178[k] = pa_z[k] * fpi0_10[k]
                   - f_11 * pc_z[k] * fpi1_10[k];

        t_179[k] = f_7 * gpg0_94[k]
                   - f_8 * gpg1_94[k]
                   + f_4 * pc_y[k] * gph_133[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pc_y, pc_z, fpi0_14, fpi0_15, \
                         fph_9, fpi1_14, fpi1_15, gpg0_95, gpg1_95, gph_134, \
                         gph_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_5 * gpg0_95[k]
                   - f_6 * gpg1_95[k]
                   + f_4 * pc_y[k] * gph_134[k];

        t_181[k] = f_4 * pc_y[k] * gph_135[k];

        t_182[k] = pa_z[k] * fpi0_14[k]
                   + f_0 * fph_9[k]
                   - f_11 * pc_z[k] * fpi1_14[k];

        t_183[k] = pa_z[k] * fpi0_15[k]
                   - f_11 * pc_z[k] * fpi1_15[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pc_x, pc_y, fph_142, fph_143, fph_144, \
                         gsh_58, gsh_59, gsh_60, gph_140, gph_142, gph_143, \
                         gph_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_13 * fph_142[k]
                   + f_1 * gsh_58[k]
                   + f_4 * pc_x[k] * gph_142[k];

        t_185[k] = f_13 * fph_143[k]
                   + f_1 * gsh_59[k]
                   + f_4 * pc_x[k] * gph_143[k];

        t_186[k] = f_13 * fph_144[k]
                   + f_1 * gsh_60[k]
                   + f_4 * pc_x[k] * gph_144[k];

        t_187[k] = f_4 * pc_y[k] * gph_140[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_z, pc_x, pc_y, pc_z, fpi0_21, fph_146, \
                         fpi1_21, gsh_62, gpg0_101, gpg1_101, gph_142, \
                         gph_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_13 * fph_146[k]
                   + f_1 * gsh_62[k]
                   + f_4 * pc_x[k] * gph_146[k];

        t_189[k] = pa_z[k] * fpi0_21[k]
                   - f_11 * pc_z[k] * fpi1_21[k];

        t_190[k] = f_17 * gpg0_101[k]
                   - f_18 * gpg1_101[k]
                   + f_4 * pc_y[k] * gph_142[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pc_y, gpg0_102, gpg0_103, gpg0_104, \
                         gpg1_102, gpg1_103, gpg1_104, gph_143, gph_144, gph_145, \
                         gph_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_9 * gpg0_102[k]
                   - f_10 * gpg1_102[k]
                   + f_4 * pc_y[k] * gph_143[k];

        t_192[k] = f_7 * gpg0_103[k]
                   - f_8 * gpg1_103[k]
                   + f_4 * pc_y[k] * gph_144[k];

        t_193[k] = f_5 * gpg0_104[k]
                   - f_6 * gpg1_104[k]
                   + f_4 * pc_y[k] * gph_145[k];

        t_194[k] = f_4 * pc_y[k] * gph_146[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_z, pc_y, pc_z, fpi0_28, fph_20, \
                         fph_21, fpi1_28, gsh_42, gpg0_104, gpg1_104, gph_146, \
                         gph_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_1 * fph_20[k]
                   + f_2 * gpg0_104[k]
                   - f_3 * gpg1_104[k]
                   + f_4 * pc_z[k] * gph_146[k];

        t_196[k] = pa_z[k] * fpi0_28[k]
                   - f_11 * pc_z[k] * fpi1_28[k];

        t_197[k] = f_1 * gsh_42[k]
                   + f_4 * pc_y[k] * gph_147[k];

        t_198[k] = f_1 * fph_21[k]
                   + f_4 * pc_z[k] * gph_147[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_z, pb_y, pc_y, pc_z, fpi0_31, fpi0_34, \
                         fpi1_31, fpi1_34, gsi0_61, gsh_44, gsi1_61, \
                         gph_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pa_z[k] * fpi0_31[k]
                   - f_11 * pc_z[k] * fpi1_31[k];

        t_200[k] = f_1 * gsh_44[k]
                   + f_4 * pc_y[k] * gph_149[k];

        t_201[k] = pb_y[k] * gsi0_61[k]
                   - f_11 * pc_y[k] * gsi1_61[k];

        t_202[k] = pa_z[k] * fpi0_34[k]
                   - f_11 * pc_z[k] * fpi1_34[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pa_z, pb_y, pc_y, pc_z, fpi0_35, fph_24, \
                         fpi1_35, gsi0_65, gsh_47, gsi1_65, gph_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = pa_z[k] * fpi0_35[k]
                   + f_1 * fph_24[k]
                   - f_11 * pc_z[k] * fpi1_35[k];

        t_204[k] = f_1 * gsh_47[k]
                   + f_4 * pc_y[k] * gph_152[k];

        t_205[k] = pb_y[k] * gsi0_65[k]
                   - f_11 * pc_y[k] * gsi1_65[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pa_z, pc_z, fpi0_38, fpi0_39, fpi0_40, fph_27, \
                         fph_28, fpi1_38, fpi1_39, fpi1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pa_z[k] * fpi0_38[k]
                   - f_11 * pc_z[k] * fpi1_38[k];

        t_207[k] = pa_z[k] * fpi0_39[k]
                   + f_1 * fph_27[k]
                   - f_11 * pc_z[k] * fpi1_39[k];

        t_208[k] = pa_z[k] * fpi0_40[k]
                   + f_12 * fph_28[k]
                   - f_11 * pc_z[k] * fpi1_40[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pb_y, pc_x, pc_y, fph_162, fph_163, \
                         gsi0_70, gsh_51, gsi1_70, gph_156, gph_162, \
                         gph_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_1 * gsh_51[k]
                   + f_4 * pc_y[k] * gph_156[k];

        t_210[k] = pb_y[k] * gsi0_70[k]
                   - f_11 * pc_y[k] * gsi1_70[k];

        t_211[k] = f_13 * fph_162[k]
                   + f_4 * pc_x[k] * gph_162[k];

        t_212[k] = f_13 * fph_163[k]
                   + f_4 * pc_x[k] * gph_163[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, fph_164, fph_165, fph_167, \
                         gsh_56, gph_161, gph_164, gph_165, gph_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_13 * fph_164[k]
                   + f_4 * pc_x[k] * gph_164[k];

        t_214[k] = f_13 * fph_165[k]
                   + f_4 * pc_x[k] * gph_165[k];

        t_215[k] = f_1 * gsh_56[k]
                   + f_4 * pc_y[k] * gph_161[k];

        t_216[k] = f_13 * fph_167[k]
                   + f_4 * pc_x[k] * gph_167[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pa_z, pb_y, pc_y, pc_z, fpi0_49, fpi1_49, \
                         gsi0_78, gsi0_79, gsh_58, gsh_59, gsi1_78, \
                         gsi1_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pa_z[k] * fpi0_49[k]
                   - f_11 * pc_z[k] * fpi1_49[k];

        t_218[k] = pb_y[k] * gsi0_78[k]
                   + f_19 * gsh_58[k]
                   - f_11 * pc_y[k] * gsi1_78[k];

        t_219[k] = pb_y[k] * gsi0_79[k]
                   + f_0 * gsh_59[k]
                   - f_11 * pc_y[k] * gsi1_79[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pb_y, pc_y, gsi0_80, gsi0_81, gsi0_83, \
                         gsh_60, gsh_61, gsh_62, gsi1_80, gsi1_81, gsi1_83, \
                         gph_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pb_y[k] * gsi0_80[k]
                   + f_13 * gsh_60[k]
                   - f_11 * pc_y[k] * gsi1_80[k];

        t_221[k] = pb_y[k] * gsi0_81[k]
                   + f_12 * gsh_61[k]
                   - f_11 * pc_y[k] * gsi1_81[k];

        t_222[k] = f_1 * gsh_62[k]
                   + f_4 * pc_y[k] * gph_167[k];

        t_223[k] = pb_y[k] * gsi0_83[k]
                   - f_11 * pc_y[k] * gsi1_83[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pc_x, pc_y, pc_z, fph_42, fph_168, \
                         gsh_42, gpg0_120, gpg1_120, gph_168, gph_169, \
                         gph_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_13 * fph_168[k]
                   + f_2 * gpg0_120[k]
                   - f_3 * gpg1_120[k]
                   + f_4 * pc_x[k] * gph_168[k];

        t_225[k] = f_4 * pc_y[k] * gph_168[k];

        t_226[k] = f_1 * fph_42[k]
                   + f_1 * gsh_42[k]
                   + f_4 * pc_z[k] * gph_168[k];

        t_227[k] = f_5 * gpg0_120[k]
                   - f_6 * gpg1_120[k]
                   + f_4 * pc_y[k] * gph_169[k];

        t_228[k] = f_4 * pc_y[k] * gph_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pc_x, pc_y, fph_173, gpg0_121, gpg0_122, \
                         gpg0_125, gpg1_121, gpg1_122, gpg1_125, gph_171, gph_172, \
                         gph_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_13 * fph_173[k]
                   + f_9 * gpg0_125[k]
                   - f_10 * gpg1_125[k]
                   + f_4 * pc_x[k] * gph_173[k];

        t_230[k] = f_7 * gpg0_121[k]
                   - f_8 * gpg1_121[k]
                   + f_4 * pc_y[k] * gph_171[k];

        t_231[k] = f_5 * gpg0_122[k]
                   - f_6 * gpg1_122[k]
                   + f_4 * pc_y[k] * gph_172[k];

        t_232[k] = f_4 * pc_y[k] * gph_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, pc_y, fph_177, gpg0_123, gpg0_124, \
                         gpg0_129, gpg1_123, gpg1_124, gpg1_129, gph_174, gph_175, \
                         gph_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_13 * fph_177[k]
                   + f_7 * gpg0_129[k]
                   - f_8 * gpg1_129[k]
                   + f_4 * pc_x[k] * gph_177[k];

        t_234[k] = f_9 * gpg0_123[k]
                   - f_10 * gpg1_123[k]
                   + f_4 * pc_y[k] * gph_174[k];

        t_235[k] = f_7 * gpg0_124[k]
                   - f_8 * gpg1_124[k]
                   + f_4 * pc_y[k] * gph_175[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pc_x, pc_y, fph_182, fph_183, gpg0_125, \
                         gpg0_134, gpg1_125, gpg1_134, gph_176, gph_177, gph_182, \
                         gph_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_5 * gpg0_125[k]
                   - f_6 * gpg1_125[k]
                   + f_4 * pc_y[k] * gph_176[k];

        t_237[k] = f_4 * pc_y[k] * gph_177[k];

        t_238[k] = f_13 * fph_182[k]
                   + f_5 * gpg0_134[k]
                   - f_6 * gpg1_134[k]
                   + f_4 * pc_x[k] * gph_182[k];

        t_239[k] = f_13 * fph_183[k]
                   + f_4 * pc_x[k] * gph_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, pc_y, fph_184, fph_185, \
                         fph_186, fph_188, gph_182, gph_184, gph_185, gph_186, \
                         gph_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_13 * fph_184[k]
                   + f_4 * pc_x[k] * gph_184[k];

        t_241[k] = f_13 * fph_185[k]
                   + f_4 * pc_x[k] * gph_185[k];

        t_242[k] = f_13 * fph_186[k]
                   + f_4 * pc_x[k] * gph_186[k];

        t_243[k] = f_4 * pc_y[k] * gph_182[k];

        t_244[k] = f_13 * fph_188[k]
                   + f_4 * pc_x[k] * gph_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, gpg0_130, gpg0_131, gpg0_132, gpg1_130, \
                         gpg1_131, gpg1_132, gph_183, gph_184, \
                         gph_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_2 * gpg0_130[k]
                   - f_3 * gpg1_130[k]
                   + f_4 * pc_y[k] * gph_183[k];

        t_246[k] = f_17 * gpg0_131[k]
                   - f_18 * gpg1_131[k]
                   + f_4 * pc_y[k] * gph_184[k];

        t_247[k] = f_9 * gpg0_132[k]
                   - f_10 * gpg1_132[k]
                   + f_4 * pc_y[k] * gph_185[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dpi1,
                                                          const size_t fpi0, const size_t fph,
                                                          const size_t fpi1, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t gpg0, const size_t gpg1,
                                                          const size_t gph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_14 = 3.0 / q;
    const auto f_15 = 1.0 / p;
    const auto f_16 = gamma / (p * q);
    const auto f_20 = 0.5 / p;
    const auto f_21 = 0.5 * gamma / (p * q);

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
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpi0_0 = buffer.data(dpi0 + 0);
    const auto *dpi0_251 = buffer.data(dpi0 + 251);
    const auto *dpi0_301 = buffer.data(dpi0 + 301);

    const auto *dpi1_0 = buffer.data(dpi1 + 0);
    const auto *dpi1_251 = buffer.data(dpi1 + 251);
    const auto *dpi1_301 = buffer.data(dpi1 + 301);

    const auto *fpi0_84 = buffer.data(fpi0 + 84);
    const auto *fpi0_87 = buffer.data(fpi0 + 87);
    const auto *fpi0_90 = buffer.data(fpi0 + 90);
    const auto *fpi0_94 = buffer.data(fpi0 + 94);
    const auto *fpi0_99 = buffer.data(fpi0 + 99);
    const auto *fpi0_105 = buffer.data(fpi0 + 105);
    const auto *fpi0_168 = buffer.data(fpi0 + 168);
    const auto *fpi0_173 = buffer.data(fpi0 + 173);
    const auto *fpi0_177 = buffer.data(fpi0 + 177);
    const auto *fpi0_180 = buffer.data(fpi0 + 180);
    const auto *fpi0_182 = buffer.data(fpi0 + 182);
    const auto *fpi0_188 = buffer.data(fpi0 + 188);
    const auto *fpi0_251 = buffer.data(fpi0 + 251);
    const auto *fpi0_301 = buffer.data(fpi0 + 301);

    const auto *fph_63 = buffer.data(fph + 63);
    const auto *fph_66 = buffer.data(fph + 66);
    const auto *fph_68 = buffer.data(fph + 68);
    const auto *fph_69 = buffer.data(fph + 69);
    const auto *fph_72 = buffer.data(fph + 72);
    const auto *fph_78 = buffer.data(fph + 78);
    const auto *fph_83 = buffer.data(fph + 83);
    const auto *fph_84 = buffer.data(fph + 84);
    const auto *fph_89 = buffer.data(fph + 89);
    const auto *fph_93 = buffer.data(fph + 93);
    const auto *fph_104 = buffer.data(fph + 104);
    const auto *fph_105 = buffer.data(fph + 105);
    const auto *fph_110 = buffer.data(fph + 110);
    const auto *fph_114 = buffer.data(fph + 114);
    const auto *fph_125 = buffer.data(fph + 125);
    const auto *fph_126 = buffer.data(fph + 126);
    const auto *fph_128 = buffer.data(fph + 128);
    const auto *fph_131 = buffer.data(fph + 131);
    const auto *fph_134 = buffer.data(fph + 134);
    const auto *fph_135 = buffer.data(fph + 135);
    const auto *fph_143 = buffer.data(fph + 143);
    const auto *fph_144 = buffer.data(fph + 144);
    const auto *fph_145 = buffer.data(fph + 145);
    const auto *fph_146 = buffer.data(fph + 146);
    const auto *fph_192 = buffer.data(fph + 192);
    const auto *fph_195 = buffer.data(fph + 195);
    const auto *fph_199 = buffer.data(fph + 199);
    const auto *fph_204 = buffer.data(fph + 204);
    const auto *fph_206 = buffer.data(fph + 206);
    const auto *fph_207 = buffer.data(fph + 207);
    const auto *fph_208 = buffer.data(fph + 208);
    const auto *fph_209 = buffer.data(fph + 209);
    const auto *fph_210 = buffer.data(fph + 210);
    const auto *fph_213 = buffer.data(fph + 213);
    const auto *fph_216 = buffer.data(fph + 216);
    const auto *fph_220 = buffer.data(fph + 220);
    const auto *fph_225 = buffer.data(fph + 225);
    const auto *fph_227 = buffer.data(fph + 227);
    const auto *fph_228 = buffer.data(fph + 228);
    const auto *fph_229 = buffer.data(fph + 229);
    const auto *fph_230 = buffer.data(fph + 230);
    const auto *fph_246 = buffer.data(fph + 246);
    const auto *fph_248 = buffer.data(fph + 248);
    const auto *fph_249 = buffer.data(fph + 249);
    const auto *fph_250 = buffer.data(fph + 250);
    const auto *fph_251 = buffer.data(fph + 251);
    const auto *fph_268 = buffer.data(fph + 268);
    const auto *fph_269 = buffer.data(fph + 269);
    const auto *fph_270 = buffer.data(fph + 270);
    const auto *fph_271 = buffer.data(fph + 271);

    const auto *fpi1_84 = buffer.data(fpi1 + 84);
    const auto *fpi1_87 = buffer.data(fpi1 + 87);
    const auto *fpi1_90 = buffer.data(fpi1 + 90);
    const auto *fpi1_94 = buffer.data(fpi1 + 94);
    const auto *fpi1_99 = buffer.data(fpi1 + 99);
    const auto *fpi1_105 = buffer.data(fpi1 + 105);
    const auto *fpi1_168 = buffer.data(fpi1 + 168);
    const auto *fpi1_173 = buffer.data(fpi1 + 173);
    const auto *fpi1_177 = buffer.data(fpi1 + 177);
    const auto *fpi1_180 = buffer.data(fpi1 + 180);
    const auto *fpi1_182 = buffer.data(fpi1 + 182);
    const auto *fpi1_188 = buffer.data(fpi1 + 188);
    const auto *fpi1_251 = buffer.data(fpi1 + 251);
    const auto *fpi1_301 = buffer.data(fpi1 + 301);

    const auto *gsi0_84 = buffer.data(gsi0 + 84);
    const auto *gsi0_87 = buffer.data(gsi0 + 87);
    const auto *gsi0_89 = buffer.data(gsi0 + 89);
    const auto *gsi0_90 = buffer.data(gsi0 + 90);
    const auto *gsi0_93 = buffer.data(gsi0 + 93);
    const auto *gsi0_94 = buffer.data(gsi0 + 94);
    const auto *gsi0_96 = buffer.data(gsi0 + 96);
    const auto *gsi0_98 = buffer.data(gsi0 + 98);
    const auto *gsi0_105 = buffer.data(gsi0 + 105);
    const auto *gsi0_107 = buffer.data(gsi0 + 107);
    const auto *gsi0_108 = buffer.data(gsi0 + 108);
    const auto *gsi0_109 = buffer.data(gsi0 + 109);
    const auto *gsi0_111 = buffer.data(gsi0 + 111);

    const auto *gsh_63 = buffer.data(gsh + 63);
    const auto *gsh_64 = buffer.data(gsh + 64);
    const auto *gsh_65 = buffer.data(gsh + 65);
    const auto *gsh_66 = buffer.data(gsh + 66);
    const auto *gsh_68 = buffer.data(gsh + 68);
    const auto *gsh_69 = buffer.data(gsh + 69);
    const auto *gsh_70 = buffer.data(gsh + 70);
    const auto *gsh_72 = buffer.data(gsh + 72);
    const auto *gsh_73 = buffer.data(gsh + 73);
    const auto *gsh_78 = buffer.data(gsh + 78);
    const auto *gsh_79 = buffer.data(gsh + 79);
    const auto *gsh_80 = buffer.data(gsh + 80);
    const auto *gsh_81 = buffer.data(gsh + 81);
    const auto *gsh_82 = buffer.data(gsh + 82);
    const auto *gsh_83 = buffer.data(gsh + 83);
    const auto *gsh_100 = buffer.data(gsh + 100);
    const auto *gsh_101 = buffer.data(gsh + 101);
    const auto *gsh_102 = buffer.data(gsh + 102);
    const auto *gsh_103 = buffer.data(gsh + 103);

    const auto *gsi1_84 = buffer.data(gsi1 + 84);
    const auto *gsi1_87 = buffer.data(gsi1 + 87);
    const auto *gsi1_89 = buffer.data(gsi1 + 89);
    const auto *gsi1_90 = buffer.data(gsi1 + 90);
    const auto *gsi1_93 = buffer.data(gsi1 + 93);
    const auto *gsi1_94 = buffer.data(gsi1 + 94);
    const auto *gsi1_96 = buffer.data(gsi1 + 96);
    const auto *gsi1_98 = buffer.data(gsi1 + 98);
    const auto *gsi1_105 = buffer.data(gsi1 + 105);
    const auto *gsi1_107 = buffer.data(gsi1 + 107);
    const auto *gsi1_108 = buffer.data(gsi1 + 108);
    const auto *gsi1_109 = buffer.data(gsi1 + 109);
    const auto *gsi1_111 = buffer.data(gsi1 + 111);

    const auto *gpg0_133 = buffer.data(gpg0 + 133);
    const auto *gpg0_134 = buffer.data(gpg0 + 134);
    const auto *gpg0_135 = buffer.data(gpg0 + 135);
    const auto *gpg0_137 = buffer.data(gpg0 + 137);
    const auto *gpg0_138 = buffer.data(gpg0 + 138);
    const auto *gpg0_140 = buffer.data(gpg0 + 140);
    const auto *gpg0_141 = buffer.data(gpg0 + 141);
    const auto *gpg0_145 = buffer.data(gpg0 + 145);
    const auto *gpg0_146 = buffer.data(gpg0 + 146);
    const auto *gpg0_147 = buffer.data(gpg0 + 147);
    const auto *gpg0_149 = buffer.data(gpg0 + 149);
    const auto *gpg0_150 = buffer.data(gpg0 + 150);
    const auto *gpg0_152 = buffer.data(gpg0 + 152);
    const auto *gpg0_153 = buffer.data(gpg0 + 153);
    const auto *gpg0_155 = buffer.data(gpg0 + 155);
    const auto *gpg0_156 = buffer.data(gpg0 + 156);
    const auto *gpg0_160 = buffer.data(gpg0 + 160);
    const auto *gpg0_161 = buffer.data(gpg0 + 161);
    const auto *gpg0_162 = buffer.data(gpg0 + 162);
    const auto *gpg0_164 = buffer.data(gpg0 + 164);
    const auto *gpg0_192 = buffer.data(gpg0 + 192);
    const auto *gpg0_193 = buffer.data(gpg0 + 193);
    const auto *gpg0_194 = buffer.data(gpg0 + 194);

    const auto *gpg1_133 = buffer.data(gpg1 + 133);
    const auto *gpg1_134 = buffer.data(gpg1 + 134);
    const auto *gpg1_135 = buffer.data(gpg1 + 135);
    const auto *gpg1_137 = buffer.data(gpg1 + 137);
    const auto *gpg1_138 = buffer.data(gpg1 + 138);
    const auto *gpg1_140 = buffer.data(gpg1 + 140);
    const auto *gpg1_141 = buffer.data(gpg1 + 141);
    const auto *gpg1_145 = buffer.data(gpg1 + 145);
    const auto *gpg1_146 = buffer.data(gpg1 + 146);
    const auto *gpg1_147 = buffer.data(gpg1 + 147);
    const auto *gpg1_149 = buffer.data(gpg1 + 149);
    const auto *gpg1_150 = buffer.data(gpg1 + 150);
    const auto *gpg1_152 = buffer.data(gpg1 + 152);
    const auto *gpg1_153 = buffer.data(gpg1 + 153);
    const auto *gpg1_155 = buffer.data(gpg1 + 155);
    const auto *gpg1_156 = buffer.data(gpg1 + 156);
    const auto *gpg1_160 = buffer.data(gpg1 + 160);
    const auto *gpg1_161 = buffer.data(gpg1 + 161);
    const auto *gpg1_162 = buffer.data(gpg1 + 162);
    const auto *gpg1_164 = buffer.data(gpg1 + 164);
    const auto *gpg1_192 = buffer.data(gpg1 + 192);
    const auto *gpg1_193 = buffer.data(gpg1 + 193);
    const auto *gpg1_194 = buffer.data(gpg1 + 194);

    const auto *gph_186 = buffer.data(gph + 186);
    const auto *gph_187 = buffer.data(gph + 187);
    const auto *gph_188 = buffer.data(gph + 188);
    const auto *gph_189 = buffer.data(gph + 189);
    const auto *gph_190 = buffer.data(gph + 190);
    const auto *gph_191 = buffer.data(gph + 191);
    const auto *gph_192 = buffer.data(gph + 192);
    const auto *gph_194 = buffer.data(gph + 194);
    const auto *gph_195 = buffer.data(gph + 195);
    const auto *gph_196 = buffer.data(gph + 196);
    const auto *gph_198 = buffer.data(gph + 198);
    const auto *gph_199 = buffer.data(gph + 199);
    const auto *gph_204 = buffer.data(gph + 204);
    const auto *gph_205 = buffer.data(gph + 205);
    const auto *gph_206 = buffer.data(gph + 206);
    const auto *gph_207 = buffer.data(gph + 207);
    const auto *gph_208 = buffer.data(gph + 208);
    const auto *gph_209 = buffer.data(gph + 209);
    const auto *gph_210 = buffer.data(gph + 210);
    const auto *gph_211 = buffer.data(gph + 211);
    const auto *gph_212 = buffer.data(gph + 212);
    const auto *gph_213 = buffer.data(gph + 213);
    const auto *gph_215 = buffer.data(gph + 215);
    const auto *gph_216 = buffer.data(gph + 216);
    const auto *gph_217 = buffer.data(gph + 217);
    const auto *gph_219 = buffer.data(gph + 219);
    const auto *gph_220 = buffer.data(gph + 220);
    const auto *gph_225 = buffer.data(gph + 225);
    const auto *gph_226 = buffer.data(gph + 226);
    const auto *gph_227 = buffer.data(gph + 227);
    const auto *gph_228 = buffer.data(gph + 228);
    const auto *gph_229 = buffer.data(gph + 229);
    const auto *gph_230 = buffer.data(gph + 230);
    const auto *gph_231 = buffer.data(gph + 231);
    const auto *gph_232 = buffer.data(gph + 232);
    const auto *gph_234 = buffer.data(gph + 234);
    const auto *gph_236 = buffer.data(gph + 236);
    const auto *gph_237 = buffer.data(gph + 237);
    const auto *gph_240 = buffer.data(gph + 240);
    const auto *gph_241 = buffer.data(gph + 241);
    const auto *gph_246 = buffer.data(gph + 246);
    const auto *gph_248 = buffer.data(gph + 248);
    const auto *gph_249 = buffer.data(gph + 249);
    const auto *gph_250 = buffer.data(gph + 250);
    const auto *gph_251 = buffer.data(gph + 251);
    const auto *gph_252 = buffer.data(gph + 252);
    const auto *gph_254 = buffer.data(gph + 254);
    const auto *gph_255 = buffer.data(gph + 255);
    const auto *gph_257 = buffer.data(gph + 257);
    const auto *gph_258 = buffer.data(gph + 258);
    const auto *gph_261 = buffer.data(gph + 261);
    const auto *gph_267 = buffer.data(gph + 267);
    const auto *gph_268 = buffer.data(gph + 268);
    const auto *gph_269 = buffer.data(gph + 269);
    const auto *gph_270 = buffer.data(gph + 270);
    const auto *gph_271 = buffer.data(gph + 271);
    const auto *gph_272 = buffer.data(gph + 272);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, gpg0_133, gpg0_134, gpg1_133, gpg1_134, \
                         gph_186, gph_187, gph_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_7 * gpg0_133[k]
                   - f_8 * gpg1_133[k]
                   + f_4 * pc_y[k] * gph_186[k];

        t_249[k] = f_5 * gpg0_134[k]
                   - f_6 * gpg1_134[k]
                   + f_4 * pc_y[k] * gph_187[k];

        t_250[k] = f_4 * pc_y[k] * gph_188[k];
    }

#pragma omp simd aligned(t_251, t_252, pa_x, pa_y, pc_x, pc_y, dpi0_0, dpi0_251, dpi1_0, \
                         dpi1_251, fpi0_84, fpi0_251, fpi1_84, \
                         fpi1_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_15 * dpi0_251[k]
                   - f_16 * dpi1_251[k]
                   + pa_x[k] * fpi0_251[k]
                   - f_11 * pc_x[k] * fpi1_251[k];

        t_252[k] = f_20 * dpi0_0[k]
                   - f_21 * dpi1_0[k]
                   + pa_y[k] * fpi0_84[k]
                   - f_11 * pc_y[k] * fpi1_84[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, pc_z, fph_63, fph_192, \
                         gsh_66, gpg0_138, gpg1_138, gph_189, gph_190, \
                         gph_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_12 * fph_63[k]
                   + f_4 * pc_y[k] * gph_189[k];

        t_254[k] = f_4 * pc_z[k] * gph_189[k];

        t_255[k] = f_12 * fph_192[k]
                   + f_1 * gsh_66[k]
                   + f_9 * gpg0_138[k]
                   - f_10 * gpg1_138[k]
                   + f_4 * pc_x[k] * gph_192[k];

        t_256[k] = f_4 * pc_z[k] * gph_190[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, fph_195, gsh_69, gpg0_135, gpg0_141, \
                         gpg1_135, gpg1_141, gph_191, gph_192, \
                         gph_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_5 * gpg0_135[k]
                   - f_6 * gpg1_135[k]
                   + f_4 * pc_z[k] * gph_191[k];

        t_258[k] = f_12 * fph_195[k]
                   + f_1 * gsh_69[k]
                   + f_7 * gpg0_141[k]
                   - f_8 * gpg1_141[k]
                   + f_4 * pc_x[k] * gph_195[k];

        t_259[k] = f_4 * pc_z[k] * gph_192[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pc_x, pc_y, pc_z, fph_68, fph_199, gsh_73, \
                         gpg0_137, gpg0_145, gpg1_137, gpg1_145, gph_194, \
                         gph_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_12 * fph_68[k]
                   + f_4 * pc_y[k] * gph_194[k];

        t_261[k] = f_7 * gpg0_137[k]
                   - f_8 * gpg1_137[k]
                   + f_4 * pc_z[k] * gph_194[k];

        t_262[k] = f_12 * fph_199[k]
                   + f_1 * gsh_73[k]
                   + f_5 * gpg0_145[k]
                   - f_6 * gpg1_145[k]
                   + f_4 * pc_x[k] * gph_199[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pc_y, pc_z, fph_72, gpg0_138, gpg0_140, \
                         gpg1_138, gpg1_140, gph_195, gph_196, \
                         gph_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_4 * pc_z[k] * gph_195[k];

        t_264[k] = f_5 * gpg0_138[k]
                   - f_6 * gpg1_138[k]
                   + f_4 * pc_z[k] * gph_196[k];

        t_265[k] = f_12 * fph_72[k]
                   + f_4 * pc_y[k] * gph_198[k];

        t_266[k] = f_9 * gpg0_140[k]
                   - f_10 * gpg1_140[k]
                   + f_4 * pc_z[k] * gph_198[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_z, fph_204, fph_206, fph_207, \
                         gsh_78, gsh_80, gsh_81, gph_199, gph_204, gph_206, \
                         gph_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_12 * fph_204[k]
                   + f_1 * gsh_78[k]
                   + f_4 * pc_x[k] * gph_204[k];

        t_268[k] = f_4 * pc_z[k] * gph_199[k];

        t_269[k] = f_12 * fph_206[k]
                   + f_1 * gsh_80[k]
                   + f_4 * pc_x[k] * gph_206[k];

        t_270[k] = f_12 * fph_207[k]
                   + f_1 * gsh_81[k]
                   + f_4 * pc_x[k] * gph_207[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, pc_x, pc_y, fph_78, fph_208, fph_209, gsh_82, \
                         gsh_83, gpg0_145, gpg1_145, gph_204, gph_208, \
                         gph_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_12 * fph_208[k]
                   + f_1 * gsh_82[k]
                   + f_4 * pc_x[k] * gph_208[k];

        t_272[k] = f_12 * fph_209[k]
                   + f_1 * gsh_83[k]
                   + f_4 * pc_x[k] * gph_209[k];

        t_273[k] = f_12 * fph_78[k]
                   + f_2 * gpg0_145[k]
                   - f_3 * gpg1_145[k]
                   + f_4 * pc_y[k] * gph_204[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, pc_z, gpg0_145, gpg0_146, gpg0_147, \
                         gpg1_145, gpg1_146, gpg1_147, gph_204, gph_205, gph_206, \
                         gph_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_4 * pc_z[k] * gph_204[k];

        t_275[k] = f_5 * gpg0_145[k]
                   - f_6 * gpg1_145[k]
                   + f_4 * pc_z[k] * gph_205[k];

        t_276[k] = f_7 * gpg0_146[k]
                   - f_8 * gpg1_146[k]
                   + f_4 * pc_z[k] * gph_206[k];

        t_277[k] = f_9 * gpg0_147[k]
                   - f_10 * gpg1_147[k]
                   + f_4 * pc_z[k] * gph_207[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pc_x, pc_y, pc_z, fph_83, fph_210, gpg0_149, \
                         gpg0_150, gpg1_149, gpg1_150, gph_209, \
                         gph_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_12 * fph_83[k]
                   + f_4 * pc_y[k] * gph_209[k];

        t_279[k] = f_2 * gpg0_149[k]
                   - f_3 * gpg1_149[k]
                   + f_4 * pc_z[k] * gph_209[k];

        t_280[k] = f_12 * fph_210[k]
                   + f_2 * gpg0_150[k]
                   - f_3 * gpg1_150[k]
                   + f_4 * pc_x[k] * gph_210[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pc_x, pc_y, pc_z, fph_84, fph_213, \
                         gsh_63, gpg0_153, gpg1_153, gph_210, gph_211, \
                         gph_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_12 * fph_84[k]
                   + f_1 * gsh_63[k]
                   + f_4 * pc_y[k] * gph_210[k];

        t_282[k] = f_4 * pc_z[k] * gph_210[k];

        t_283[k] = f_12 * fph_213[k]
                   + f_9 * gpg0_153[k]
                   - f_10 * gpg1_153[k]
                   + f_4 * pc_x[k] * gph_213[k];

        t_284[k] = f_4 * pc_z[k] * gph_211[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pc_x, pc_z, fph_216, gpg0_150, gpg0_156, \
                         gpg1_150, gpg1_156, gph_212, gph_213, \
                         gph_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_5 * gpg0_150[k]
                   - f_6 * gpg1_150[k]
                   + f_4 * pc_z[k] * gph_212[k];

        t_286[k] = f_12 * fph_216[k]
                   + f_7 * gpg0_156[k]
                   - f_8 * gpg1_156[k]
                   + f_4 * pc_x[k] * gph_216[k];

        t_287[k] = f_4 * pc_z[k] * gph_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_x, pc_y, pc_z, fph_89, fph_220, gsh_68, \
                         gpg0_152, gpg0_160, gpg1_152, gpg1_160, gph_215, \
                         gph_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_12 * fph_89[k]
                   + f_1 * gsh_68[k]
                   + f_4 * pc_y[k] * gph_215[k];

        t_289[k] = f_7 * gpg0_152[k]
                   - f_8 * gpg1_152[k]
                   + f_4 * pc_z[k] * gph_215[k];

        t_290[k] = f_12 * fph_220[k]
                   + f_5 * gpg0_160[k]
                   - f_6 * gpg1_160[k]
                   + f_4 * pc_x[k] * gph_220[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_y, pc_z, fph_93, gsh_72, gpg0_153, \
                         gpg0_155, gpg1_153, gpg1_155, gph_216, gph_217, \
                         gph_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_4 * pc_z[k] * gph_216[k];

        t_292[k] = f_5 * gpg0_153[k]
                   - f_6 * gpg1_153[k]
                   + f_4 * pc_z[k] * gph_217[k];

        t_293[k] = f_12 * fph_93[k]
                   + f_1 * gsh_72[k]
                   + f_4 * pc_y[k] * gph_219[k];

        t_294[k] = f_9 * gpg0_155[k]
                   - f_10 * gpg1_155[k]
                   + f_4 * pc_z[k] * gph_219[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, pc_x, pc_z, fph_225, fph_227, \
                         fph_228, fph_229, gph_220, gph_225, gph_227, gph_228, \
                         gph_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_12 * fph_225[k]
                   + f_4 * pc_x[k] * gph_225[k];

        t_296[k] = f_4 * pc_z[k] * gph_220[k];

        t_297[k] = f_12 * fph_227[k]
                   + f_4 * pc_x[k] * gph_227[k];

        t_298[k] = f_12 * fph_228[k]
                   + f_4 * pc_x[k] * gph_228[k];

        t_299[k] = f_12 * fph_229[k]
                   + f_4 * pc_x[k] * gph_229[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pa_x, pc_x, pc_z, dpi0_301, dpi1_301, fpi0_301, \
                         fph_230, fpi1_301, gph_225, gph_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_12 * fph_230[k]
                   + f_4 * pc_x[k] * gph_230[k];

        t_301[k] = f_20 * dpi0_301[k]
                   - f_21 * dpi1_301[k]
                   + pa_x[k] * fpi0_301[k]
                   - f_11 * pc_x[k] * fpi1_301[k];

        t_302[k] = f_4 * pc_z[k] * gph_225[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, pc_z, gpg0_160, gpg0_161, gpg0_162, gpg1_160, \
                         gpg1_161, gpg1_162, gph_226, gph_227, \
                         gph_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_5 * gpg0_160[k]
                   - f_6 * gpg1_160[k]
                   + f_4 * pc_z[k] * gph_226[k];

        t_304[k] = f_7 * gpg0_161[k]
                   - f_8 * gpg1_161[k]
                   + f_4 * pc_z[k] * gph_227[k];

        t_305[k] = f_9 * gpg0_162[k]
                   - f_10 * gpg1_162[k]
                   + f_4 * pc_z[k] * gph_228[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pb_z, pc_y, pc_z, fph_104, fph_105, \
                         gsi0_84, gsh_83, gsi1_84, gpg0_164, gpg1_164, gph_230, \
                         gph_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_12 * fph_104[k]
                   + f_1 * gsh_83[k]
                   + f_4 * pc_y[k] * gph_230[k];

        t_307[k] = f_2 * gpg0_164[k]
                   - f_3 * gpg1_164[k]
                   + f_4 * pc_z[k] * gph_230[k];

        t_308[k] = pb_z[k] * gsi0_84[k]
                   - f_11 * pc_z[k] * gsi1_84[k];

        t_309[k] = f_12 * fph_105[k]
                   + f_4 * pc_y[k] * gph_231[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pb_z, pc_z, gsi0_87, gsi0_89, gsh_63, \
                         gsh_64, gsh_65, gsi1_87, gsi1_89, gph_231, \
                         gph_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_1 * gsh_63[k]
                   + f_4 * pc_z[k] * gph_231[k];

        t_311[k] = pb_z[k] * gsi0_87[k]
                   - f_11 * pc_z[k] * gsi1_87[k];

        t_312[k] = f_1 * gsh_64[k]
                   + f_4 * pc_z[k] * gph_232[k];

        t_313[k] = pb_z[k] * gsi0_89[k]
                   + f_12 * gsh_65[k]
                   - f_11 * pc_z[k] * gsi1_89[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pb_z, pc_y, pc_z, fph_110, gsi0_90, \
                         gsi0_93, gsh_66, gsh_68, gsi1_90, gsi1_93, gph_234, \
                         gph_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pb_z[k] * gsi0_90[k]
                   - f_11 * pc_z[k] * gsi1_90[k];

        t_315[k] = f_1 * gsh_66[k]
                   + f_4 * pc_z[k] * gph_234[k];

        t_316[k] = f_12 * fph_110[k]
                   + f_4 * pc_y[k] * gph_236[k];

        t_317[k] = pb_z[k] * gsi0_93[k]
                   + f_13 * gsh_68[k]
                   - f_11 * pc_z[k] * gsi1_93[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pb_z, pc_y, pc_z, fph_114, gsi0_94, \
                         gsi0_96, gsh_69, gsh_70, gsi1_94, gsi1_96, gph_237, \
                         gph_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * gsi0_94[k]
                   - f_11 * pc_z[k] * gsi1_94[k];

        t_319[k] = f_1 * gsh_69[k]
                   + f_4 * pc_z[k] * gph_237[k];

        t_320[k] = pb_z[k] * gsi0_96[k]
                   + f_12 * gsh_70[k]
                   - f_11 * pc_z[k] * gsi1_96[k];

        t_321[k] = f_12 * fph_114[k]
                   + f_4 * pc_y[k] * gph_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pb_z, pc_x, pc_z, fph_246, fph_248, \
                         gsi0_98, gsh_72, gsh_73, gsi1_98, gph_241, gph_246, \
                         gph_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = pb_z[k] * gsi0_98[k]
                   + f_0 * gsh_72[k]
                   - f_11 * pc_z[k] * gsi1_98[k];

        t_323[k] = f_12 * fph_246[k]
                   + f_4 * pc_x[k] * gph_246[k];

        t_324[k] = f_1 * gsh_73[k]
                   + f_4 * pc_z[k] * gph_241[k];

        t_325[k] = f_12 * fph_248[k]
                   + f_4 * pc_x[k] * gph_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pb_z, pc_x, pc_z, fph_249, fph_250, \
                         fph_251, gsi0_105, gsi1_105, gph_249, gph_250, \
                         gph_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_12 * fph_249[k]
                   + f_4 * pc_x[k] * gph_249[k];

        t_327[k] = f_12 * fph_250[k]
                   + f_4 * pc_x[k] * gph_250[k];

        t_328[k] = f_12 * fph_251[k]
                   + f_4 * pc_x[k] * gph_251[k];

        t_329[k] = pb_z[k] * gsi0_105[k]
                   - f_11 * pc_z[k] * gsi1_105[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pb_z, pc_z, gsi0_107, gsi0_108, gsh_78, gsh_79, \
                         gsh_80, gsi1_107, gsi1_108, gph_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_1 * gsh_78[k]
                   + f_4 * pc_z[k] * gph_246[k];

        t_331[k] = pb_z[k] * gsi0_107[k]
                   + f_12 * gsh_79[k]
                   - f_11 * pc_z[k] * gsi1_107[k];

        t_332[k] = pb_z[k] * gsi0_108[k]
                   + f_13 * gsh_80[k]
                   - f_11 * pc_z[k] * gsi1_108[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pb_z, pc_y, pc_z, fph_125, gsi0_109, gsi0_111, \
                         gsh_81, gsh_83, gsi1_109, gsi1_111, gph_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pb_z[k] * gsi0_109[k]
                   + f_0 * gsh_81[k]
                   - f_11 * pc_z[k] * gsi1_109[k];

        t_334[k] = f_12 * fph_125[k]
                   + f_4 * pc_y[k] * gph_251[k];

        t_335[k] = pb_z[k] * gsi0_111[k]
                   + f_14 * gsh_83[k]
                   - f_11 * pc_z[k] * gsi1_111[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pa_y, pa_z, pc_y, pc_z, fpi0_87, \
                         fpi0_168, fph_63, fph_126, fpi1_87, fpi1_168, \
                         gph_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_y[k] * fpi0_168[k]
                   - f_11 * pc_y[k] * fpi1_168[k];

        t_337[k] = f_1 * fph_126[k]
                   + f_4 * pc_y[k] * gph_252[k];

        t_338[k] = f_1 * fph_63[k]
                   + f_4 * pc_z[k] * gph_252[k];

        t_339[k] = pa_z[k] * fpi0_87[k]
                   - f_11 * pc_z[k] * fpi1_87[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pa_y, pa_z, pc_y, pc_z, fpi0_90, \
                         fpi0_173, fph_66, fph_128, fpi1_90, fpi1_173, gph_254, \
                         gph_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_1 * fph_128[k]
                   + f_4 * pc_y[k] * gph_254[k];

        t_341[k] = pa_y[k] * fpi0_173[k]
                   - f_11 * pc_y[k] * fpi1_173[k];

        t_342[k] = pa_z[k] * fpi0_90[k]
                   - f_11 * pc_z[k] * fpi1_90[k];

        t_343[k] = f_1 * fph_66[k]
                   + f_4 * pc_z[k] * gph_255[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pa_y, pa_z, pc_y, pc_z, fpi0_94, \
                         fpi0_177, fph_69, fph_131, fpi1_94, fpi1_177, gph_257, \
                         gph_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_1 * fph_131[k]
                   + f_4 * pc_y[k] * gph_257[k];

        t_345[k] = pa_y[k] * fpi0_177[k]
                   - f_11 * pc_y[k] * fpi1_177[k];

        t_346[k] = pa_z[k] * fpi0_94[k]
                   - f_11 * pc_z[k] * fpi1_94[k];

        t_347[k] = f_1 * fph_69[k]
                   + f_4 * pc_z[k] * gph_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pa_y, pc_y, fpi0_180, fpi0_182, fph_134, \
                         fph_135, fpi1_180, fpi1_182, gph_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = pa_y[k] * fpi0_180[k]
                   + f_12 * fph_134[k]
                   - f_11 * pc_y[k] * fpi1_180[k];

        t_349[k] = f_1 * fph_135[k]
                   + f_4 * pc_y[k] * gph_261[k];

        t_350[k] = pa_y[k] * fpi0_182[k]
                   - f_11 * pc_y[k] * fpi1_182[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pa_z, pc_x, pc_z, fpi0_99, fph_268, fph_269, \
                         fpi1_99, gsh_100, gsh_101, gph_268, gph_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = pa_z[k] * fpi0_99[k]
                   - f_11 * pc_z[k] * fpi1_99[k];

        t_352[k] = f_12 * fph_268[k]
                   + f_1 * gsh_100[k]
                   + f_4 * pc_x[k] * gph_268[k];

        t_353[k] = f_12 * fph_269[k]
                   + f_1 * gsh_101[k]
                   + f_4 * pc_x[k] * gph_269[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pa_y, pc_x, pc_y, fpi0_188, fph_270, fph_271, \
                         fpi1_188, gsh_102, gsh_103, gph_270, gph_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_12 * fph_270[k]
                   + f_1 * gsh_102[k]
                   + f_4 * pc_x[k] * gph_270[k];

        t_355[k] = f_12 * fph_271[k]
                   + f_1 * gsh_103[k]
                   + f_4 * pc_x[k] * gph_271[k];

        t_356[k] = pa_y[k] * fpi0_188[k]
                   - f_11 * pc_y[k] * fpi1_188[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pa_z, pc_y, pc_z, fpi0_105, fph_78, fph_143, \
                         fpi1_105, gpg0_192, gpg1_192, gph_267, \
                         gph_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pa_z[k] * fpi0_105[k]
                   - f_11 * pc_z[k] * fpi1_105[k];

        t_358[k] = f_1 * fph_78[k]
                   + f_4 * pc_z[k] * gph_267[k];

        t_359[k] = f_1 * fph_143[k]
                   + f_9 * gpg0_192[k]
                   - f_10 * gpg1_192[k]
                   + f_4 * pc_y[k] * gph_269[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pc_y, fph_144, fph_145, fph_146, gpg0_193, \
                         gpg0_194, gpg1_193, gpg1_194, gph_270, gph_271, \
                         gph_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_1 * fph_144[k]
                   + f_7 * gpg0_193[k]
                   - f_8 * gpg1_193[k]
                   + f_4 * pc_y[k] * gph_270[k];

        t_361[k] = f_1 * fph_145[k]
                   + f_5 * gpg0_194[k]
                   - f_6 * gpg1_194[k]
                   + f_4 * pc_y[k] * gph_271[k];

        t_362[k] = f_1 * fph_146[k]
                   + f_4 * pc_y[k] * gph_272[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dpi1,
                                                          const size_t fpi0, const size_t fph,
                                                          const size_t fpi1, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t gpg0, const size_t gpg1,
                                                          const size_t gph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_14 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;
    const auto f_20 = 0.5 / p;
    const auto f_21 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpi0_0 = buffer.data(dpi0 + 0);
    const auto *dpi0_387 = buffer.data(dpi0 + 387);
    const auto *dpi0_388 = buffer.data(dpi0 + 388);
    const auto *dpi0_389 = buffer.data(dpi0 + 389);

    const auto *dpi1_0 = buffer.data(dpi1 + 0);
    const auto *dpi1_387 = buffer.data(dpi1 + 387);
    const auto *dpi1_388 = buffer.data(dpi1 + 388);
    const auto *dpi1_389 = buffer.data(dpi1 + 389);

    const auto *fpi0_113 = buffer.data(fpi0 + 113);
    const auto *fpi0_115 = buffer.data(fpi0 + 115);
    const auto *fpi0_118 = buffer.data(fpi0 + 118);
    const auto *fpi0_122 = buffer.data(fpi0 + 122);
    const auto *fpi0_133 = buffer.data(fpi0 + 133);
    const auto *fpi0_168 = buffer.data(fpi0 + 168);
    const auto *fpi0_195 = buffer.data(fpi0 + 195);
    const auto *fpi0_224 = buffer.data(fpi0 + 224);
    const auto *fpi0_226 = buffer.data(fpi0 + 226);
    const auto *fpi0_229 = buffer.data(fpi0 + 229);
    const auto *fpi0_233 = buffer.data(fpi0 + 233);
    const auto *fpi0_236 = buffer.data(fpi0 + 236);
    const auto *fpi0_238 = buffer.data(fpi0 + 238);
    const auto *fpi0_251 = buffer.data(fpi0 + 251);
    const auto *fpi0_387 = buffer.data(fpi0 + 387);
    const auto *fpi0_388 = buffer.data(fpi0 + 388);
    const auto *fpi0_389 = buffer.data(fpi0 + 389);

    const auto *fph_84 = buffer.data(fph + 84);
    const auto *fph_87 = buffer.data(fph + 87);
    const auto *fph_90 = buffer.data(fph + 90);
    const auto *fph_99 = buffer.data(fph + 99);
    const auto *fph_104 = buffer.data(fph + 104);
    const auto *fph_108 = buffer.data(fph + 108);
    const auto *fph_111 = buffer.data(fph + 111);
    const auto *fph_120 = buffer.data(fph + 120);
    const auto *fph_126 = buffer.data(fph + 126);
    const auto *fph_146 = buffer.data(fph + 146);
    const auto *fph_147 = buffer.data(fph + 147);
    const auto *fph_149 = buffer.data(fph + 149);
    const auto *fph_152 = buffer.data(fph + 152);
    const auto *fph_156 = buffer.data(fph + 156);
    const auto *fph_167 = buffer.data(fph + 167);
    const auto *fph_168 = buffer.data(fph + 168);
    const auto *fph_170 = buffer.data(fph + 170);
    const auto *fph_173 = buffer.data(fph + 173);
    const auto *fph_176 = buffer.data(fph + 176);
    const auto *fph_177 = buffer.data(fph + 177);
    const auto *fph_183 = buffer.data(fph + 183);
    const auto *fph_185 = buffer.data(fph + 185);
    const auto *fph_186 = buffer.data(fph + 186);
    const auto *fph_187 = buffer.data(fph + 187);
    const auto *fph_188 = buffer.data(fph + 188);
    const auto *fph_273 = buffer.data(fph + 273);
    const auto *fph_278 = buffer.data(fph + 278);
    const auto *fph_282 = buffer.data(fph + 282);
    const auto *fph_285 = buffer.data(fph + 285);
    const auto *fph_287 = buffer.data(fph + 287);
    const auto *fph_288 = buffer.data(fph + 288);
    const auto *fph_289 = buffer.data(fph + 289);
    const auto *fph_290 = buffer.data(fph + 290);
    const auto *fph_291 = buffer.data(fph + 291);
    const auto *fph_292 = buffer.data(fph + 292);
    const auto *fph_293 = buffer.data(fph + 293);
    const auto *fph_297 = buffer.data(fph + 297);
    const auto *fph_300 = buffer.data(fph + 300);
    const auto *fph_304 = buffer.data(fph + 304);
    const auto *fph_309 = buffer.data(fph + 309);
    const auto *fph_310 = buffer.data(fph + 310);
    const auto *fph_311 = buffer.data(fph + 311);
    const auto *fph_312 = buffer.data(fph + 312);
    const auto *fph_313 = buffer.data(fph + 313);
    const auto *fph_314 = buffer.data(fph + 314);
    const auto *fph_320 = buffer.data(fph + 320);
    const auto *fph_324 = buffer.data(fph + 324);
    const auto *fph_329 = buffer.data(fph + 329);
    const auto *fph_330 = buffer.data(fph + 330);
    const auto *fph_331 = buffer.data(fph + 331);
    const auto *fph_332 = buffer.data(fph + 332);
    const auto *fph_333 = buffer.data(fph + 333);
    const auto *fph_335 = buffer.data(fph + 335);
    const auto *fph_351 = buffer.data(fph + 351);
    const auto *fph_352 = buffer.data(fph + 352);
    const auto *fph_353 = buffer.data(fph + 353);
    const auto *fph_354 = buffer.data(fph + 354);
    const auto *fph_356 = buffer.data(fph + 356);

    const auto *fpi1_113 = buffer.data(fpi1 + 113);
    const auto *fpi1_115 = buffer.data(fpi1 + 115);
    const auto *fpi1_118 = buffer.data(fpi1 + 118);
    const auto *fpi1_122 = buffer.data(fpi1 + 122);
    const auto *fpi1_133 = buffer.data(fpi1 + 133);
    const auto *fpi1_168 = buffer.data(fpi1 + 168);
    const auto *fpi1_195 = buffer.data(fpi1 + 195);
    const auto *fpi1_224 = buffer.data(fpi1 + 224);
    const auto *fpi1_226 = buffer.data(fpi1 + 226);
    const auto *fpi1_229 = buffer.data(fpi1 + 229);
    const auto *fpi1_233 = buffer.data(fpi1 + 233);
    const auto *fpi1_236 = buffer.data(fpi1 + 236);
    const auto *fpi1_238 = buffer.data(fpi1 + 238);
    const auto *fpi1_251 = buffer.data(fpi1 + 251);
    const auto *fpi1_387 = buffer.data(fpi1 + 387);
    const auto *fpi1_388 = buffer.data(fpi1 + 388);
    const auto *fpi1_389 = buffer.data(fpi1 + 389);

    const auto *gsi0_140 = buffer.data(gsi0 + 140);
    const auto *gsi0_143 = buffer.data(gsi0 + 143);
    const auto *gsi0_145 = buffer.data(gsi0 + 145);
    const auto *gsi0_146 = buffer.data(gsi0 + 146);
    const auto *gsi0_147 = buffer.data(gsi0 + 147);
    const auto *gsi0_149 = buffer.data(gsi0 + 149);
    const auto *gsi0_150 = buffer.data(gsi0 + 150);
    const auto *gsi0_151 = buffer.data(gsi0 + 151);
    const auto *gsi0_152 = buffer.data(gsi0 + 152);
    const auto *gsi0_154 = buffer.data(gsi0 + 154);
    const auto *gsi0_161 = buffer.data(gsi0 + 161);
    const auto *gsi0_162 = buffer.data(gsi0 + 162);
    const auto *gsi0_163 = buffer.data(gsi0 + 163);
    const auto *gsi0_164 = buffer.data(gsi0 + 164);
    const auto *gsi0_165 = buffer.data(gsi0 + 165);

    const auto *gsh_86 = buffer.data(gsh + 86);
    const auto *gsh_87 = buffer.data(gsh + 87);
    const auto *gsh_89 = buffer.data(gsh + 89);
    const auto *gsh_90 = buffer.data(gsh + 90);
    const auto *gsh_93 = buffer.data(gsh + 93);
    const auto *gsh_99 = buffer.data(gsh + 99);
    const auto *gsh_104 = buffer.data(gsh + 104);
    const auto *gsh_105 = buffer.data(gsh + 105);
    const auto *gsh_106 = buffer.data(gsh + 106);
    const auto *gsh_107 = buffer.data(gsh + 107);
    const auto *gsh_108 = buffer.data(gsh + 108);
    const auto *gsh_109 = buffer.data(gsh + 109);
    const auto *gsh_110 = buffer.data(gsh + 110);
    const auto *gsh_111 = buffer.data(gsh + 111);
    const auto *gsh_112 = buffer.data(gsh + 112);
    const auto *gsh_113 = buffer.data(gsh + 113);
    const auto *gsh_114 = buffer.data(gsh + 114);
    const auto *gsh_119 = buffer.data(gsh + 119);
    const auto *gsh_120 = buffer.data(gsh + 120);
    const auto *gsh_121 = buffer.data(gsh + 121);
    const auto *gsh_122 = buffer.data(gsh + 122);
    const auto *gsh_123 = buffer.data(gsh + 123);
    const auto *gsh_124 = buffer.data(gsh + 124);
    const auto *gsh_125 = buffer.data(gsh + 125);

    const auto *gsi1_140 = buffer.data(gsi1 + 140);
    const auto *gsi1_143 = buffer.data(gsi1 + 143);
    const auto *gsi1_145 = buffer.data(gsi1 + 145);
    const auto *gsi1_146 = buffer.data(gsi1 + 146);
    const auto *gsi1_147 = buffer.data(gsi1 + 147);
    const auto *gsi1_149 = buffer.data(gsi1 + 149);
    const auto *gsi1_150 = buffer.data(gsi1 + 150);
    const auto *gsi1_151 = buffer.data(gsi1 + 151);
    const auto *gsi1_152 = buffer.data(gsi1 + 152);
    const auto *gsi1_154 = buffer.data(gsi1 + 154);
    const auto *gsi1_161 = buffer.data(gsi1 + 161);
    const auto *gsi1_162 = buffer.data(gsi1 + 162);
    const auto *gsi1_163 = buffer.data(gsi1 + 163);
    const auto *gsi1_164 = buffer.data(gsi1 + 164);
    const auto *gsi1_165 = buffer.data(gsi1 + 165);

    const auto *gpg0_195 = buffer.data(gpg0 + 195);
    const auto *gpg0_200 = buffer.data(gpg0 + 200);
    const auto *gpg0_204 = buffer.data(gpg0 + 204);
    const auto *gpg0_207 = buffer.data(gpg0 + 207);
    const auto *gpg0_209 = buffer.data(gpg0 + 209);
    const auto *gpg0_213 = buffer.data(gpg0 + 213);
    const auto *gpg0_216 = buffer.data(gpg0 + 216);
    const auto *gpg0_220 = buffer.data(gpg0 + 220);
    const auto *gpg0_222 = buffer.data(gpg0 + 222);
    const auto *gpg0_223 = buffer.data(gpg0 + 223);
    const auto *gpg0_224 = buffer.data(gpg0 + 224);
    const auto *gpg0_225 = buffer.data(gpg0 + 225);
    const auto *gpg0_226 = buffer.data(gpg0 + 226);
    const auto *gpg0_227 = buffer.data(gpg0 + 227);
    const auto *gpg0_228 = buffer.data(gpg0 + 228);
    const auto *gpg0_229 = buffer.data(gpg0 + 229);
    const auto *gpg0_230 = buffer.data(gpg0 + 230);
    const auto *gpg0_234 = buffer.data(gpg0 + 234);
    const auto *gpg0_235 = buffer.data(gpg0 + 235);
    const auto *gpg0_236 = buffer.data(gpg0 + 236);
    const auto *gpg0_237 = buffer.data(gpg0 + 237);
    const auto *gpg0_238 = buffer.data(gpg0 + 238);
    const auto *gpg0_239 = buffer.data(gpg0 + 239);

    const auto *gpg1_195 = buffer.data(gpg1 + 195);
    const auto *gpg1_200 = buffer.data(gpg1 + 200);
    const auto *gpg1_204 = buffer.data(gpg1 + 204);
    const auto *gpg1_207 = buffer.data(gpg1 + 207);
    const auto *gpg1_209 = buffer.data(gpg1 + 209);
    const auto *gpg1_213 = buffer.data(gpg1 + 213);
    const auto *gpg1_216 = buffer.data(gpg1 + 216);
    const auto *gpg1_220 = buffer.data(gpg1 + 220);
    const auto *gpg1_222 = buffer.data(gpg1 + 222);
    const auto *gpg1_223 = buffer.data(gpg1 + 223);
    const auto *gpg1_224 = buffer.data(gpg1 + 224);
    const auto *gpg1_225 = buffer.data(gpg1 + 225);
    const auto *gpg1_226 = buffer.data(gpg1 + 226);
    const auto *gpg1_227 = buffer.data(gpg1 + 227);
    const auto *gpg1_228 = buffer.data(gpg1 + 228);
    const auto *gpg1_229 = buffer.data(gpg1 + 229);
    const auto *gpg1_230 = buffer.data(gpg1 + 230);
    const auto *gpg1_234 = buffer.data(gpg1 + 234);
    const auto *gpg1_235 = buffer.data(gpg1 + 235);
    const auto *gpg1_236 = buffer.data(gpg1 + 236);
    const auto *gpg1_237 = buffer.data(gpg1 + 237);
    const auto *gpg1_238 = buffer.data(gpg1 + 238);
    const auto *gpg1_239 = buffer.data(gpg1 + 239);

    const auto *gph_273 = buffer.data(gph + 273);
    const auto *gph_275 = buffer.data(gph + 275);
    const auto *gph_276 = buffer.data(gph + 276);
    const auto *gph_278 = buffer.data(gph + 278);
    const auto *gph_279 = buffer.data(gph + 279);
    const auto *gph_282 = buffer.data(gph + 282);
    const auto *gph_285 = buffer.data(gph + 285);
    const auto *gph_287 = buffer.data(gph + 287);
    const auto *gph_288 = buffer.data(gph + 288);
    const auto *gph_289 = buffer.data(gph + 289);
    const auto *gph_290 = buffer.data(gph + 290);
    const auto *gph_291 = buffer.data(gph + 291);
    const auto *gph_292 = buffer.data(gph + 292);
    const auto *gph_293 = buffer.data(gph + 293);
    const auto *gph_294 = buffer.data(gph + 294);
    const auto *gph_296 = buffer.data(gph + 296);
    const auto *gph_297 = buffer.data(gph + 297);
    const auto *gph_299 = buffer.data(gph + 299);
    const auto *gph_300 = buffer.data(gph + 300);
    const auto *gph_303 = buffer.data(gph + 303);
    const auto *gph_304 = buffer.data(gph + 304);
    const auto *gph_309 = buffer.data(gph + 309);
    const auto *gph_310 = buffer.data(gph + 310);
    const auto *gph_311 = buffer.data(gph + 311);
    const auto *gph_312 = buffer.data(gph + 312);
    const auto *gph_313 = buffer.data(gph + 313);
    const auto *gph_314 = buffer.data(gph + 314);
    const auto *gph_315 = buffer.data(gph + 315);
    const auto *gph_316 = buffer.data(gph + 316);
    const auto *gph_317 = buffer.data(gph + 317);
    const auto *gph_318 = buffer.data(gph + 318);
    const auto *gph_319 = buffer.data(gph + 319);
    const auto *gph_320 = buffer.data(gph + 320);
    const auto *gph_321 = buffer.data(gph + 321);
    const auto *gph_322 = buffer.data(gph + 322);
    const auto *gph_323 = buffer.data(gph + 323);
    const auto *gph_324 = buffer.data(gph + 324);
    const auto *gph_329 = buffer.data(gph + 329);
    const auto *gph_330 = buffer.data(gph + 330);
    const auto *gph_331 = buffer.data(gph + 331);
    const auto *gph_332 = buffer.data(gph + 332);
    const auto *gph_333 = buffer.data(gph + 333);
    const auto *gph_334 = buffer.data(gph + 334);
    const auto *gph_335 = buffer.data(gph + 335);
    const auto *gph_336 = buffer.data(gph + 336);
    const auto *gph_338 = buffer.data(gph + 338);
    const auto *gph_341 = buffer.data(gph + 341);
    const auto *gph_345 = buffer.data(gph + 345);
    const auto *gph_350 = buffer.data(gph + 350);
    const auto *gph_351 = buffer.data(gph + 351);
    const auto *gph_352 = buffer.data(gph + 352);
    const auto *gph_353 = buffer.data(gph + 353);
    const auto *gph_354 = buffer.data(gph + 354);
    const auto *gph_356 = buffer.data(gph + 356);

#pragma omp simd aligned(t_363, t_364, t_365, pa_y, pa_z, pc_x, pc_y, pc_z, fpi0_113, \
                         fpi0_195, fph_273, fpi1_113, fpi1_195, gpg0_195, gpg1_195, \
                         gph_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pa_y[k] * fpi0_195[k]
                   - f_11 * pc_y[k] * fpi1_195[k];

        t_364[k] = f_12 * fph_273[k]
                   + f_2 * gpg0_195[k]
                   - f_3 * gpg1_195[k]
                   + f_4 * pc_x[k] * gph_273[k];

        t_365[k] = pa_z[k] * fpi0_113[k]
                   - f_11 * pc_z[k] * fpi1_113[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, pa_z, pc_y, pc_z, fpi0_115, fph_84, fph_149, \
                         fpi1_115, gsh_86, gph_273, gph_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_1 * fph_84[k]
                   + f_4 * pc_z[k] * gph_273[k];

        t_367[k] = pa_z[k] * fpi0_115[k]
                   - f_11 * pc_z[k] * fpi1_115[k];

        t_368[k] = f_1 * fph_149[k]
                   + f_1 * gsh_86[k]
                   + f_4 * pc_y[k] * gph_275[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pa_z, pc_x, pc_z, fpi0_118, fph_87, fph_278, \
                         fpi1_118, gpg0_200, gpg1_200, gph_276, \
                         gph_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_12 * fph_278[k]
                   + f_9 * gpg0_200[k]
                   - f_10 * gpg1_200[k]
                   + f_4 * pc_x[k] * gph_278[k];

        t_370[k] = pa_z[k] * fpi0_118[k]
                   - f_11 * pc_z[k] * fpi1_118[k];

        t_371[k] = f_1 * fph_87[k]
                   + f_4 * pc_z[k] * gph_276[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pa_z, pc_x, pc_y, pc_z, fpi0_122, fph_152, \
                         fph_282, fpi1_122, gsh_89, gpg0_204, gpg1_204, gph_278, \
                         gph_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_1 * fph_152[k]
                   + f_1 * gsh_89[k]
                   + f_4 * pc_y[k] * gph_278[k];

        t_373[k] = f_12 * fph_282[k]
                   + f_7 * gpg0_204[k]
                   - f_8 * gpg1_204[k]
                   + f_4 * pc_x[k] * gph_282[k];

        t_374[k] = pa_z[k] * fpi0_122[k]
                   - f_11 * pc_z[k] * fpi1_122[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_y, pc_z, fph_90, fph_156, fph_285, \
                         gsh_93, gpg0_207, gpg1_207, gph_279, gph_282, \
                         gph_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_1 * fph_90[k]
                   + f_4 * pc_z[k] * gph_279[k];

        t_376[k] = f_12 * fph_285[k]
                   + f_5 * gpg0_207[k]
                   - f_6 * gpg1_207[k]
                   + f_4 * pc_x[k] * gph_285[k];

        t_377[k] = f_1 * fph_156[k]
                   + f_1 * gsh_93[k]
                   + f_4 * pc_y[k] * gph_282[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, fph_287, fph_288, fph_289, fph_290, \
                         gpg0_209, gpg1_209, gph_287, gph_288, gph_289, \
                         gph_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_12 * fph_287[k]
                   + f_5 * gpg0_209[k]
                   - f_6 * gpg1_209[k]
                   + f_4 * pc_x[k] * gph_287[k];

        t_379[k] = f_12 * fph_288[k]
                   + f_4 * pc_x[k] * gph_288[k];

        t_380[k] = f_12 * fph_289[k]
                   + f_4 * pc_x[k] * gph_289[k];

        t_381[k] = f_12 * fph_290[k]
                   + f_4 * pc_x[k] * gph_290[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pa_z, pc_x, pc_z, fpi0_133, fph_291, \
                         fph_292, fph_293, fpi1_133, gph_291, gph_292, \
                         gph_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_12 * fph_291[k]
                   + f_4 * pc_x[k] * gph_291[k];

        t_383[k] = f_12 * fph_292[k]
                   + f_4 * pc_x[k] * gph_292[k];

        t_384[k] = f_12 * fph_293[k]
                   + f_4 * pc_x[k] * gph_293[k];

        t_385[k] = pa_z[k] * fpi0_133[k]
                   - f_11 * pc_z[k] * fpi1_133[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, pa_x, pc_x, pc_z, dpi0_387, dpi0_388, dpi1_387, \
                         dpi1_388, fpi0_387, fpi0_388, fph_99, fpi1_387, fpi1_388, \
                         gph_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_1 * fph_99[k]
                   + f_4 * pc_z[k] * gph_288[k];

        t_387[k] = f_20 * dpi0_387[k]
                   - f_21 * dpi1_387[k]
                   + pa_x[k] * fpi0_387[k]
                   - f_11 * pc_x[k] * fpi1_387[k];

        t_388[k] = f_20 * dpi0_388[k]
                   - f_21 * dpi1_388[k]
                   + pa_x[k] * fpi0_388[k]
                   - f_11 * pc_x[k] * fpi1_388[k];
    }

#pragma omp simd aligned(t_389, t_390, pa_x, pc_x, pc_y, dpi0_389, dpi1_389, fpi0_389, \
                         fph_167, fpi1_389, gsh_104, gph_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_20 * dpi0_389[k]
                   - f_21 * dpi1_389[k]
                   + pa_x[k] * fpi0_389[k]
                   - f_11 * pc_x[k] * fpi1_389[k];

        t_390[k] = f_1 * fph_167[k]
                   + f_1 * gsh_104[k]
                   + f_4 * pc_y[k] * gph_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pa_y, pc_y, pc_z, fpi0_224, fph_104, fph_168, \
                         fpi1_224, gpg0_209, gpg1_209, gph_293, \
                         gph_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_1 * fph_104[k]
                   + f_2 * gpg0_209[k]
                   - f_3 * gpg1_209[k]
                   + f_4 * pc_z[k] * gph_293[k];

        t_392[k] = pa_y[k] * fpi0_224[k]
                   - f_11 * pc_y[k] * fpi1_224[k];

        t_393[k] = f_1 * fph_168[k]
                   + f_4 * pc_y[k] * gph_294[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pa_y, pc_x, pc_y, fpi0_226, fph_170, fph_297, \
                         fpi1_226, gpg0_213, gpg1_213, gph_296, \
                         gph_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = pa_y[k] * fpi0_226[k]
                   - f_11 * pc_y[k] * fpi1_226[k];

        t_395[k] = f_12 * fph_297[k]
                   + f_9 * gpg0_213[k]
                   - f_10 * gpg1_213[k]
                   + f_4 * pc_x[k] * gph_297[k];

        t_396[k] = f_1 * fph_170[k]
                   + f_4 * pc_y[k] * gph_296[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, pa_y, pc_x, pc_y, pc_z, fpi0_229, fph_108, \
                         fph_300, fpi1_229, gsh_87, gpg0_216, gpg1_216, gph_297, \
                         gph_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = pa_y[k] * fpi0_229[k]
                   - f_11 * pc_y[k] * fpi1_229[k];

        t_398[k] = f_12 * fph_300[k]
                   + f_7 * gpg0_216[k]
                   - f_8 * gpg1_216[k]
                   + f_4 * pc_x[k] * gph_300[k];

        t_399[k] = f_1 * fph_108[k]
                   + f_1 * gsh_87[k]
                   + f_4 * pc_z[k] * gph_297[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pa_y, pc_x, pc_y, fpi0_233, fph_173, fph_304, \
                         fpi1_233, gpg0_220, gpg1_220, gph_299, \
                         gph_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_1 * fph_173[k]
                   + f_4 * pc_y[k] * gph_299[k];

        t_401[k] = pa_y[k] * fpi0_233[k]
                   - f_11 * pc_y[k] * fpi1_233[k];

        t_402[k] = f_12 * fph_304[k]
                   + f_5 * gpg0_220[k]
                   - f_6 * gpg1_220[k]
                   + f_4 * pc_x[k] * gph_304[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, pa_y, pc_y, pc_z, fpi0_236, fph_111, fph_176, \
                         fph_177, fpi1_236, gsh_90, gph_300, gph_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_1 * fph_111[k]
                   + f_1 * gsh_90[k]
                   + f_4 * pc_z[k] * gph_300[k];

        t_404[k] = pa_y[k] * fpi0_236[k]
                   + f_12 * fph_176[k]
                   - f_11 * pc_y[k] * fpi1_236[k];

        t_405[k] = f_1 * fph_177[k]
                   + f_4 * pc_y[k] * gph_303[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pa_y, pc_x, pc_y, fpi0_238, fph_309, \
                         fph_310, fph_311, fpi1_238, gph_309, gph_310, \
                         gph_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pa_y[k] * fpi0_238[k]
                   - f_11 * pc_y[k] * fpi1_238[k];

        t_407[k] = f_12 * fph_309[k]
                   + f_4 * pc_x[k] * gph_309[k];

        t_408[k] = f_12 * fph_310[k]
                   + f_4 * pc_x[k] * gph_310[k];

        t_409[k] = f_12 * fph_311[k]
                   + f_4 * pc_x[k] * gph_311[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, fph_183, fph_312, fph_313, \
                         fph_314, gpg0_220, gpg1_220, gph_309, gph_312, gph_313, \
                         gph_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_12 * fph_312[k]
                   + f_4 * pc_x[k] * gph_312[k];

        t_411[k] = f_12 * fph_313[k]
                   + f_4 * pc_x[k] * gph_313[k];

        t_412[k] = f_12 * fph_314[k]
                   + f_4 * pc_x[k] * gph_314[k];

        t_413[k] = f_1 * fph_183[k]
                   + f_2 * gpg0_220[k]
                   - f_3 * gpg1_220[k]
                   + f_4 * pc_y[k] * gph_309[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, fph_120, fph_185, fph_186, gsh_99, \
                         gpg0_222, gpg0_223, gpg1_222, gpg1_223, gph_309, gph_311, \
                         gph_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_1 * fph_120[k]
                   + f_1 * gsh_99[k]
                   + f_4 * pc_z[k] * gph_309[k];

        t_415[k] = f_1 * fph_185[k]
                   + f_9 * gpg0_222[k]
                   - f_10 * gpg1_222[k]
                   + f_4 * pc_y[k] * gph_311[k];

        t_416[k] = f_1 * fph_186[k]
                   + f_7 * gpg0_223[k]
                   - f_8 * gpg1_223[k]
                   + f_4 * pc_y[k] * gph_312[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_y, pc_y, fpi0_251, fph_187, fph_188, \
                         fpi1_251, gpg0_224, gpg1_224, gph_313, \
                         gph_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_1 * fph_187[k]
                   + f_5 * gpg0_224[k]
                   - f_6 * gpg1_224[k]
                   + f_4 * pc_y[k] * gph_313[k];

        t_418[k] = f_1 * fph_188[k]
                   + f_4 * pc_y[k] * gph_314[k];

        t_419[k] = pa_y[k] * fpi0_251[k]
                   - f_11 * pc_y[k] * fpi1_251[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pa_z, pc_y, pc_z, dpi0_0, dpi1_0, \
                         fpi0_168, fph_126, fpi1_168, gpg0_225, gpg1_225, gph_315, \
                         gph_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_20 * dpi0_0[k]
                   - f_21 * dpi1_0[k]
                   + pa_z[k] * fpi0_168[k]
                   - f_11 * pc_z[k] * fpi1_168[k];

        t_421[k] = f_4 * pc_y[k] * gph_315[k];

        t_422[k] = f_12 * fph_126[k]
                   + f_4 * pc_z[k] * gph_315[k];

        t_423[k] = f_5 * gpg0_225[k]
                   - f_6 * gpg1_225[k]
                   + f_4 * pc_y[k] * gph_316[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, fph_320, gsh_110, gpg0_226, \
                         gpg0_230, gpg1_226, gpg1_230, gph_317, gph_318, \
                         gph_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_4 * pc_y[k] * gph_317[k];

        t_425[k] = f_12 * fph_320[k]
                   + f_1 * gsh_110[k]
                   + f_9 * gpg0_230[k]
                   - f_10 * gpg1_230[k]
                   + f_4 * pc_x[k] * gph_320[k];

        t_426[k] = f_7 * gpg0_226[k]
                   - f_8 * gpg1_226[k]
                   + f_4 * pc_y[k] * gph_318[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_x, pc_y, fph_324, gsh_114, gpg0_227, \
                         gpg0_234, gpg1_227, gpg1_234, gph_319, gph_320, \
                         gph_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_5 * gpg0_227[k]
                   - f_6 * gpg1_227[k]
                   + f_4 * pc_y[k] * gph_319[k];

        t_428[k] = f_4 * pc_y[k] * gph_320[k];

        t_429[k] = f_12 * fph_324[k]
                   + f_1 * gsh_114[k]
                   + f_7 * gpg0_234[k]
                   - f_8 * gpg1_234[k]
                   + f_4 * pc_x[k] * gph_324[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, pc_y, gpg0_228, gpg0_229, gpg0_230, \
                         gpg1_228, gpg1_229, gpg1_230, gph_321, gph_322, gph_323, \
                         gph_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_9 * gpg0_228[k]
                   - f_10 * gpg1_228[k]
                   + f_4 * pc_y[k] * gph_321[k];

        t_431[k] = f_7 * gpg0_229[k]
                   - f_8 * gpg1_229[k]
                   + f_4 * pc_y[k] * gph_322[k];

        t_432[k] = f_5 * gpg0_230[k]
                   - f_6 * gpg1_230[k]
                   + f_4 * pc_y[k] * gph_323[k];

        t_433[k] = f_4 * pc_y[k] * gph_324[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, pc_x, fph_329, fph_330, fph_331, gsh_119, \
                         gsh_120, gsh_121, gpg0_239, gpg1_239, gph_329, gph_330, \
                         gph_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_12 * fph_329[k]
                   + f_1 * gsh_119[k]
                   + f_5 * gpg0_239[k]
                   - f_6 * gpg1_239[k]
                   + f_4 * pc_x[k] * gph_329[k];

        t_435[k] = f_12 * fph_330[k]
                   + f_1 * gsh_120[k]
                   + f_4 * pc_x[k] * gph_330[k];

        t_436[k] = f_12 * fph_331[k]
                   + f_1 * gsh_121[k]
                   + f_4 * pc_x[k] * gph_331[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_x, pc_y, fph_332, fph_333, fph_335, \
                         gsh_122, gsh_123, gsh_125, gph_329, gph_332, gph_333, \
                         gph_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_12 * fph_332[k]
                   + f_1 * gsh_122[k]
                   + f_4 * pc_x[k] * gph_332[k];

        t_438[k] = f_12 * fph_333[k]
                   + f_1 * gsh_123[k]
                   + f_4 * pc_x[k] * gph_333[k];

        t_439[k] = f_4 * pc_y[k] * gph_329[k];

        t_440[k] = f_12 * fph_335[k]
                   + f_1 * gsh_125[k]
                   + f_4 * pc_x[k] * gph_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, gpg0_235, gpg0_236, gpg0_237, gpg1_235, \
                         gpg1_236, gpg1_237, gph_330, gph_331, \
                         gph_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_2 * gpg0_235[k]
                   - f_3 * gpg1_235[k]
                   + f_4 * pc_y[k] * gph_330[k];

        t_442[k] = f_17 * gpg0_236[k]
                   - f_18 * gpg1_236[k]
                   + f_4 * pc_y[k] * gph_331[k];

        t_443[k] = f_9 * gpg0_237[k]
                   - f_10 * gpg1_237[k]
                   + f_4 * pc_y[k] * gph_332[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pc_y, pc_z, fph_146, gpg0_238, gpg0_239, \
                         gpg1_238, gpg1_239, gph_333, gph_334, \
                         gph_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_7 * gpg0_238[k]
                   - f_8 * gpg1_238[k]
                   + f_4 * pc_y[k] * gph_333[k];

        t_445[k] = f_5 * gpg0_239[k]
                   - f_6 * gpg1_239[k]
                   + f_4 * pc_y[k] * gph_334[k];

        t_446[k] = f_4 * pc_y[k] * gph_335[k];

        t_447[k] = f_12 * fph_146[k]
                   + f_2 * gpg0_239[k]
                   - f_3 * gpg1_239[k]
                   + f_4 * pc_z[k] * gph_335[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pb_y, pc_y, pc_z, fph_147, gsi0_140, \
                         gsi0_143, gsh_105, gsh_106, gsi1_140, gsi1_143, \
                         gph_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pb_y[k] * gsi0_140[k]
                   - f_11 * pc_y[k] * gsi1_140[k];

        t_449[k] = f_1 * gsh_105[k]
                   + f_4 * pc_y[k] * gph_336[k];

        t_450[k] = f_12 * fph_147[k]
                   + f_4 * pc_z[k] * gph_336[k];

        t_451[k] = pb_y[k] * gsi0_143[k]
                   + f_12 * gsh_106[k]
                   - f_11 * pc_y[k] * gsi1_143[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pb_y, pc_y, gsi0_145, gsi0_146, gsi0_147, \
                         gsh_107, gsh_108, gsh_109, gsi1_145, gsi1_146, gsi1_147, \
                         gph_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_1 * gsh_107[k]
                   + f_4 * pc_y[k] * gph_338[k];

        t_453[k] = pb_y[k] * gsi0_145[k]
                   - f_11 * pc_y[k] * gsi1_145[k];

        t_454[k] = pb_y[k] * gsi0_146[k]
                   + f_13 * gsh_108[k]
                   - f_11 * pc_y[k] * gsi1_146[k];

        t_455[k] = pb_y[k] * gsi0_147[k]
                   + f_12 * gsh_109[k]
                   - f_11 * pc_y[k] * gsi1_147[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pb_y, pc_y, gsi0_149, gsi0_150, gsi0_151, \
                         gsh_110, gsh_111, gsh_112, gsi1_149, gsi1_150, gsi1_151, \
                         gph_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_1 * gsh_110[k]
                   + f_4 * pc_y[k] * gph_341[k];

        t_457[k] = pb_y[k] * gsi0_149[k]
                   - f_11 * pc_y[k] * gsi1_149[k];

        t_458[k] = pb_y[k] * gsi0_150[k]
                   + f_0 * gsh_111[k]
                   - f_11 * pc_y[k] * gsi1_150[k];

        t_459[k] = pb_y[k] * gsi0_151[k]
                   + f_13 * gsh_112[k]
                   - f_11 * pc_y[k] * gsi1_151[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pb_y, pc_x, pc_y, fph_351, gsi0_152, \
                         gsi0_154, gsh_113, gsh_114, gsi1_152, gsi1_154, gph_345, \
                         gph_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = pb_y[k] * gsi0_152[k]
                   + f_12 * gsh_113[k]
                   - f_11 * pc_y[k] * gsi1_152[k];

        t_461[k] = f_1 * gsh_114[k]
                   + f_4 * pc_y[k] * gph_345[k];

        t_462[k] = pb_y[k] * gsi0_154[k]
                   - f_11 * pc_y[k] * gsi1_154[k];

        t_463[k] = f_12 * fph_351[k]
                   + f_4 * pc_x[k] * gph_351[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, pc_x, pc_y, fph_352, fph_353, fph_354, \
                         gsh_119, gph_350, gph_352, gph_353, gph_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_12 * fph_352[k]
                   + f_4 * pc_x[k] * gph_352[k];

        t_465[k] = f_12 * fph_353[k]
                   + f_4 * pc_x[k] * gph_353[k];

        t_466[k] = f_12 * fph_354[k]
                   + f_4 * pc_x[k] * gph_354[k];

        t_467[k] = f_1 * gsh_119[k]
                   + f_4 * pc_y[k] * gph_350[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pb_y, pc_x, pc_y, fph_356, gsi0_161, gsi0_162, \
                         gsh_120, gsh_121, gsi1_161, gsi1_162, \
                         gph_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_12 * fph_356[k]
                   + f_4 * pc_x[k] * gph_356[k];

        t_469[k] = pb_y[k] * gsi0_161[k]
                   + f_14 * gsh_120[k]
                   - f_11 * pc_y[k] * gsi1_161[k];

        t_470[k] = pb_y[k] * gsi0_162[k]
                   + f_19 * gsh_121[k]
                   - f_11 * pc_y[k] * gsi1_162[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pb_y, pc_y, gsi0_163, gsi0_164, gsi0_165, \
                         gsh_122, gsh_123, gsh_124, gsi1_163, gsi1_164, \
                         gsi1_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pb_y[k] * gsi0_163[k]
                   + f_0 * gsh_122[k]
                   - f_11 * pc_y[k] * gsi1_163[k];

        t_472[k] = pb_y[k] * gsi0_164[k]
                   + f_13 * gsh_123[k]
                   - f_11 * pc_y[k] * gsi1_164[k];

        t_473[k] = pb_y[k] * gsi0_165[k]
                   + f_12 * gsh_124[k]
                   - f_11 * pc_y[k] * gsi1_165[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dpi1,
                                                          const size_t fpi0, const size_t fph,
                                                          const size_t fpi1, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t gpg0, const size_t gpg1,
                                                          const size_t gph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_14 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_20 = 0.5 / p;
    const auto f_21 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpi0_503 = buffer.data(dpi0 + 503);

    const auto *dpi1_503 = buffer.data(dpi1 + 503);

    const auto *fpi0_252 = buffer.data(fpi0 + 252);
    const auto *fpi0_255 = buffer.data(fpi0 + 255);
    const auto *fpi0_503 = buffer.data(fpi0 + 503);
    const auto *fpi0_532 = buffer.data(fpi0 + 532);
    const auto *fpi0_535 = buffer.data(fpi0 + 535);
    const auto *fpi0_538 = buffer.data(fpi0 + 538);
    const auto *fpi0_542 = buffer.data(fpi0 + 542);
    const auto *fpi0_553 = buffer.data(fpi0 + 553);
    const auto *fpi0_555 = buffer.data(fpi0 + 555);
    const auto *fpi0_556 = buffer.data(fpi0 + 556);
    const auto *fpi0_557 = buffer.data(fpi0 + 557);
    const auto *fpi0_558 = buffer.data(fpi0 + 558);
    const auto *fpi0_559 = buffer.data(fpi0 + 559);
    const auto *fpi0_565 = buffer.data(fpi0 + 565);
    const auto *fpi0_569 = buffer.data(fpi0 + 569);
    const auto *fpi0_572 = buffer.data(fpi0 + 572);
    const auto *fpi0_574 = buffer.data(fpi0 + 574);
    const auto *fpi0_581 = buffer.data(fpi0 + 581);
    const auto *fpi0_583 = buffer.data(fpi0 + 583);
    const auto *fpi0_584 = buffer.data(fpi0 + 584);
    const auto *fpi0_585 = buffer.data(fpi0 + 585);
    const auto *fpi0_587 = buffer.data(fpi0 + 587);

    const auto *fph_168 = buffer.data(fph + 168);
    const auto *fph_189 = buffer.data(fph + 189);
    const auto *fph_194 = buffer.data(fph + 194);
    const auto *fph_198 = buffer.data(fph + 198);
    const auto *fph_204 = buffer.data(fph + 204);
    const auto *fph_209 = buffer.data(fph + 209);
    const auto *fph_210 = buffer.data(fph + 210);
    const auto *fph_215 = buffer.data(fph + 215);
    const auto *fph_219 = buffer.data(fph + 219);
    const auto *fph_231 = buffer.data(fph + 231);
    const auto *fph_236 = buffer.data(fph + 236);
    const auto *fph_240 = buffer.data(fph + 240);
    const auto *fph_251 = buffer.data(fph + 251);
    const auto *fph_252 = buffer.data(fph + 252);
    const auto *fph_357 = buffer.data(fph + 357);
    const auto *fph_362 = buffer.data(fph + 362);
    const auto *fph_366 = buffer.data(fph + 366);
    const auto *fph_371 = buffer.data(fph + 371);
    const auto *fph_372 = buffer.data(fph + 372);
    const auto *fph_373 = buffer.data(fph + 373);
    const auto *fph_374 = buffer.data(fph + 374);
    const auto *fph_375 = buffer.data(fph + 375);
    const auto *fph_377 = buffer.data(fph + 377);
    const auto *fph_378 = buffer.data(fph + 378);
    const auto *fph_381 = buffer.data(fph + 381);
    const auto *fph_384 = buffer.data(fph + 384);
    const auto *fph_388 = buffer.data(fph + 388);
    const auto *fph_393 = buffer.data(fph + 393);
    const auto *fph_395 = buffer.data(fph + 395);
    const auto *fph_396 = buffer.data(fph + 396);
    const auto *fph_397 = buffer.data(fph + 397);
    const auto *fph_398 = buffer.data(fph + 398);
    const auto *fph_399 = buffer.data(fph + 399);
    const auto *fph_402 = buffer.data(fph + 402);
    const auto *fph_405 = buffer.data(fph + 405);
    const auto *fph_409 = buffer.data(fph + 409);
    const auto *fph_414 = buffer.data(fph + 414);
    const auto *fph_416 = buffer.data(fph + 416);
    const auto *fph_417 = buffer.data(fph + 417);
    const auto *fph_418 = buffer.data(fph + 418);
    const auto *fph_419 = buffer.data(fph + 419);
    const auto *fph_425 = buffer.data(fph + 425);
    const auto *fph_429 = buffer.data(fph + 429);
    const auto *fph_432 = buffer.data(fph + 432);
    const auto *fph_434 = buffer.data(fph + 434);
    const auto *fph_435 = buffer.data(fph + 435);
    const auto *fph_437 = buffer.data(fph + 437);
    const auto *fph_438 = buffer.data(fph + 438);
    const auto *fph_439 = buffer.data(fph + 439);
    const auto *fph_440 = buffer.data(fph + 440);

    const auto *fpi1_252 = buffer.data(fpi1 + 252);
    const auto *fpi1_255 = buffer.data(fpi1 + 255);
    const auto *fpi1_503 = buffer.data(fpi1 + 503);
    const auto *fpi1_532 = buffer.data(fpi1 + 532);
    const auto *fpi1_535 = buffer.data(fpi1 + 535);
    const auto *fpi1_538 = buffer.data(fpi1 + 538);
    const auto *fpi1_542 = buffer.data(fpi1 + 542);
    const auto *fpi1_553 = buffer.data(fpi1 + 553);
    const auto *fpi1_555 = buffer.data(fpi1 + 555);
    const auto *fpi1_556 = buffer.data(fpi1 + 556);
    const auto *fpi1_557 = buffer.data(fpi1 + 557);
    const auto *fpi1_558 = buffer.data(fpi1 + 558);
    const auto *fpi1_559 = buffer.data(fpi1 + 559);
    const auto *fpi1_565 = buffer.data(fpi1 + 565);
    const auto *fpi1_569 = buffer.data(fpi1 + 569);
    const auto *fpi1_572 = buffer.data(fpi1 + 572);
    const auto *fpi1_574 = buffer.data(fpi1 + 574);
    const auto *fpi1_581 = buffer.data(fpi1 + 581);
    const auto *fpi1_583 = buffer.data(fpi1 + 583);
    const auto *fpi1_584 = buffer.data(fpi1 + 584);
    const auto *fpi1_585 = buffer.data(fpi1 + 585);
    const auto *fpi1_587 = buffer.data(fpi1 + 587);

    const auto *gsi0_167 = buffer.data(gsi0 + 167);
    const auto *gsi0_168 = buffer.data(gsi0 + 168);
    const auto *gsi0_171 = buffer.data(gsi0 + 171);
    const auto *gsi0_174 = buffer.data(gsi0 + 174);
    const auto *gsi0_178 = buffer.data(gsi0 + 178);

    const auto *gsh_105 = buffer.data(gsh + 105);
    const auto *gsh_125 = buffer.data(gsh + 125);
    const auto *gsh_126 = buffer.data(gsh + 126);
    const auto *gsh_127 = buffer.data(gsh + 127);
    const auto *gsh_129 = buffer.data(gsh + 129);
    const auto *gsh_131 = buffer.data(gsh + 131);
    const auto *gsh_132 = buffer.data(gsh + 132);
    const auto *gsh_135 = buffer.data(gsh + 135);
    const auto *gsh_136 = buffer.data(gsh + 136);
    const auto *gsh_141 = buffer.data(gsh + 141);
    const auto *gsh_143 = buffer.data(gsh + 143);
    const auto *gsh_144 = buffer.data(gsh + 144);
    const auto *gsh_145 = buffer.data(gsh + 145);
    const auto *gsh_146 = buffer.data(gsh + 146);

    const auto *gsi1_167 = buffer.data(gsi1 + 167);
    const auto *gsi1_168 = buffer.data(gsi1 + 168);
    const auto *gsi1_171 = buffer.data(gsi1 + 171);
    const auto *gsi1_174 = buffer.data(gsi1 + 174);
    const auto *gsi1_178 = buffer.data(gsi1 + 178);

    const auto *gpg0_255 = buffer.data(gpg0 + 255);
    const auto *gpg0_256 = buffer.data(gpg0 + 256);
    const auto *gpg0_257 = buffer.data(gpg0 + 257);
    const auto *gpg0_258 = buffer.data(gpg0 + 258);
    const auto *gpg0_259 = buffer.data(gpg0 + 259);
    const auto *gpg0_260 = buffer.data(gpg0 + 260);
    const auto *gpg0_264 = buffer.data(gpg0 + 264);
    const auto *gpg0_265 = buffer.data(gpg0 + 265);
    const auto *gpg0_266 = buffer.data(gpg0 + 266);
    const auto *gpg0_267 = buffer.data(gpg0 + 267);
    const auto *gpg0_268 = buffer.data(gpg0 + 268);
    const auto *gpg0_269 = buffer.data(gpg0 + 269);
    const auto *gpg0_270 = buffer.data(gpg0 + 270);
    const auto *gpg0_272 = buffer.data(gpg0 + 272);
    const auto *gpg0_273 = buffer.data(gpg0 + 273);
    const auto *gpg0_275 = buffer.data(gpg0 + 275);
    const auto *gpg0_276 = buffer.data(gpg0 + 276);
    const auto *gpg0_280 = buffer.data(gpg0 + 280);
    const auto *gpg0_281 = buffer.data(gpg0 + 281);
    const auto *gpg0_282 = buffer.data(gpg0 + 282);
    const auto *gpg0_284 = buffer.data(gpg0 + 284);
    const auto *gpg0_285 = buffer.data(gpg0 + 285);
    const auto *gpg0_287 = buffer.data(gpg0 + 287);
    const auto *gpg0_288 = buffer.data(gpg0 + 288);
    const auto *gpg0_290 = buffer.data(gpg0 + 290);

    const auto *gpg1_255 = buffer.data(gpg1 + 255);
    const auto *gpg1_256 = buffer.data(gpg1 + 256);
    const auto *gpg1_257 = buffer.data(gpg1 + 257);
    const auto *gpg1_258 = buffer.data(gpg1 + 258);
    const auto *gpg1_259 = buffer.data(gpg1 + 259);
    const auto *gpg1_260 = buffer.data(gpg1 + 260);
    const auto *gpg1_264 = buffer.data(gpg1 + 264);
    const auto *gpg1_265 = buffer.data(gpg1 + 265);
    const auto *gpg1_266 = buffer.data(gpg1 + 266);
    const auto *gpg1_267 = buffer.data(gpg1 + 267);
    const auto *gpg1_268 = buffer.data(gpg1 + 268);
    const auto *gpg1_269 = buffer.data(gpg1 + 269);
    const auto *gpg1_270 = buffer.data(gpg1 + 270);
    const auto *gpg1_272 = buffer.data(gpg1 + 272);
    const auto *gpg1_273 = buffer.data(gpg1 + 273);
    const auto *gpg1_275 = buffer.data(gpg1 + 275);
    const auto *gpg1_276 = buffer.data(gpg1 + 276);
    const auto *gpg1_280 = buffer.data(gpg1 + 280);
    const auto *gpg1_281 = buffer.data(gpg1 + 281);
    const auto *gpg1_282 = buffer.data(gpg1 + 282);
    const auto *gpg1_284 = buffer.data(gpg1 + 284);
    const auto *gpg1_285 = buffer.data(gpg1 + 285);
    const auto *gpg1_287 = buffer.data(gpg1 + 287);
    const auto *gpg1_288 = buffer.data(gpg1 + 288);
    const auto *gpg1_290 = buffer.data(gpg1 + 290);

    const auto *gph_356 = buffer.data(gph + 356);
    const auto *gph_357 = buffer.data(gph + 357);
    const auto *gph_358 = buffer.data(gph + 358);
    const auto *gph_359 = buffer.data(gph + 359);
    const auto *gph_360 = buffer.data(gph + 360);
    const auto *gph_361 = buffer.data(gph + 361);
    const auto *gph_362 = buffer.data(gph + 362);
    const auto *gph_363 = buffer.data(gph + 363);
    const auto *gph_364 = buffer.data(gph + 364);
    const auto *gph_365 = buffer.data(gph + 365);
    const auto *gph_366 = buffer.data(gph + 366);
    const auto *gph_371 = buffer.data(gph + 371);
    const auto *gph_372 = buffer.data(gph + 372);
    const auto *gph_373 = buffer.data(gph + 373);
    const auto *gph_374 = buffer.data(gph + 374);
    const auto *gph_375 = buffer.data(gph + 375);
    const auto *gph_376 = buffer.data(gph + 376);
    const auto *gph_377 = buffer.data(gph + 377);
    const auto *gph_378 = buffer.data(gph + 378);
    const auto *gph_379 = buffer.data(gph + 379);
    const auto *gph_380 = buffer.data(gph + 380);
    const auto *gph_381 = buffer.data(gph + 381);
    const auto *gph_383 = buffer.data(gph + 383);
    const auto *gph_384 = buffer.data(gph + 384);
    const auto *gph_385 = buffer.data(gph + 385);
    const auto *gph_387 = buffer.data(gph + 387);
    const auto *gph_388 = buffer.data(gph + 388);
    const auto *gph_393 = buffer.data(gph + 393);
    const auto *gph_394 = buffer.data(gph + 394);
    const auto *gph_395 = buffer.data(gph + 395);
    const auto *gph_396 = buffer.data(gph + 396);
    const auto *gph_397 = buffer.data(gph + 397);
    const auto *gph_398 = buffer.data(gph + 398);
    const auto *gph_399 = buffer.data(gph + 399);
    const auto *gph_400 = buffer.data(gph + 400);
    const auto *gph_401 = buffer.data(gph + 401);
    const auto *gph_402 = buffer.data(gph + 402);
    const auto *gph_404 = buffer.data(gph + 404);
    const auto *gph_405 = buffer.data(gph + 405);
    const auto *gph_406 = buffer.data(gph + 406);
    const auto *gph_408 = buffer.data(gph + 408);
    const auto *gph_409 = buffer.data(gph + 409);
    const auto *gph_414 = buffer.data(gph + 414);
    const auto *gph_416 = buffer.data(gph + 416);
    const auto *gph_417 = buffer.data(gph + 417);
    const auto *gph_418 = buffer.data(gph + 418);
    const auto *gph_419 = buffer.data(gph + 419);
    const auto *gph_420 = buffer.data(gph + 420);
    const auto *gph_421 = buffer.data(gph + 421);
    const auto *gph_423 = buffer.data(gph + 423);
    const auto *gph_425 = buffer.data(gph + 425);
    const auto *gph_426 = buffer.data(gph + 426);
    const auto *gph_429 = buffer.data(gph + 429);
    const auto *gph_430 = buffer.data(gph + 430);
    const auto *gph_435 = buffer.data(gph + 435);
    const auto *gph_437 = buffer.data(gph + 437);
    const auto *gph_438 = buffer.data(gph + 438);
    const auto *gph_439 = buffer.data(gph + 439);
    const auto *gph_440 = buffer.data(gph + 440);
    const auto *gph_441 = buffer.data(gph + 441);

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_y, pc_x, pc_y, fph_357, gsi0_167, \
                         gsh_125, gsi1_167, gpg0_255, gpg1_255, gph_356, \
                         gph_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_1 * gsh_125[k]
                   + f_4 * pc_y[k] * gph_356[k];

        t_475[k] = pb_y[k] * gsi0_167[k]
                   - f_11 * pc_y[k] * gsi1_167[k];

        t_476[k] = f_12 * fph_357[k]
                   + f_2 * gpg0_255[k]
                   - f_3 * gpg1_255[k]
                   + f_4 * pc_x[k] * gph_357[k];

        t_477[k] = f_4 * pc_y[k] * gph_357[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_y, pc_z, fph_168, gsh_105, gpg0_255, \
                         gpg1_255, gph_357, gph_358, gph_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_12 * fph_168[k]
                   + f_1 * gsh_105[k]
                   + f_4 * pc_z[k] * gph_357[k];

        t_479[k] = f_5 * gpg0_255[k]
                   - f_6 * gpg1_255[k]
                   + f_4 * pc_y[k] * gph_358[k];

        t_480[k] = f_4 * pc_y[k] * gph_359[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, fph_362, gpg0_256, gpg0_257, \
                         gpg0_260, gpg1_256, gpg1_257, gpg1_260, gph_360, gph_361, \
                         gph_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_12 * fph_362[k]
                   + f_9 * gpg0_260[k]
                   - f_10 * gpg1_260[k]
                   + f_4 * pc_x[k] * gph_362[k];

        t_482[k] = f_7 * gpg0_256[k]
                   - f_8 * gpg1_256[k]
                   + f_4 * pc_y[k] * gph_360[k];

        t_483[k] = f_5 * gpg0_257[k]
                   - f_6 * gpg1_257[k]
                   + f_4 * pc_y[k] * gph_361[k];

        t_484[k] = f_4 * pc_y[k] * gph_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_y, fph_366, gpg0_258, gpg0_259, \
                         gpg0_264, gpg1_258, gpg1_259, gpg1_264, gph_363, gph_364, \
                         gph_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_12 * fph_366[k]
                   + f_7 * gpg0_264[k]
                   - f_8 * gpg1_264[k]
                   + f_4 * pc_x[k] * gph_366[k];

        t_486[k] = f_9 * gpg0_258[k]
                   - f_10 * gpg1_258[k]
                   + f_4 * pc_y[k] * gph_363[k];

        t_487[k] = f_7 * gpg0_259[k]
                   - f_8 * gpg1_259[k]
                   + f_4 * pc_y[k] * gph_364[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pc_x, pc_y, fph_371, fph_372, gpg0_260, \
                         gpg0_269, gpg1_260, gpg1_269, gph_365, gph_366, gph_371, \
                         gph_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_5 * gpg0_260[k]
                   - f_6 * gpg1_260[k]
                   + f_4 * pc_y[k] * gph_365[k];

        t_489[k] = f_4 * pc_y[k] * gph_366[k];

        t_490[k] = f_12 * fph_371[k]
                   + f_5 * gpg0_269[k]
                   - f_6 * gpg1_269[k]
                   + f_4 * pc_x[k] * gph_371[k];

        t_491[k] = f_12 * fph_372[k]
                   + f_4 * pc_x[k] * gph_372[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, pc_x, pc_y, fph_373, fph_374, \
                         fph_375, fph_377, gph_371, gph_373, gph_374, gph_375, \
                         gph_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_12 * fph_373[k]
                   + f_4 * pc_x[k] * gph_373[k];

        t_493[k] = f_12 * fph_374[k]
                   + f_4 * pc_x[k] * gph_374[k];

        t_494[k] = f_12 * fph_375[k]
                   + f_4 * pc_x[k] * gph_375[k];

        t_495[k] = f_4 * pc_y[k] * gph_371[k];

        t_496[k] = f_12 * fph_377[k]
                   + f_4 * pc_x[k] * gph_377[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, pc_y, gpg0_265, gpg0_266, gpg0_267, gpg1_265, \
                         gpg1_266, gpg1_267, gph_372, gph_373, \
                         gph_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = f_2 * gpg0_265[k]
                   - f_3 * gpg1_265[k]
                   + f_4 * pc_y[k] * gph_372[k];

        t_498[k] = f_17 * gpg0_266[k]
                   - f_18 * gpg1_266[k]
                   + f_4 * pc_y[k] * gph_373[k];

        t_499[k] = f_9 * gpg0_267[k]
                   - f_10 * gpg1_267[k]
                   + f_4 * pc_y[k] * gph_374[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pc_y, gpg0_268, gpg0_269, gpg1_268, gpg1_269, \
                         gph_375, gph_376, gph_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_7 * gpg0_268[k]
                   - f_8 * gpg1_268[k]
                   + f_4 * pc_y[k] * gph_375[k];

        t_501[k] = f_5 * gpg0_269[k]
                   - f_6 * gpg1_269[k]
                   + f_4 * pc_y[k] * gph_376[k];

        t_502[k] = f_4 * pc_y[k] * gph_377[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pa_x, pc_x, pc_y, dpi0_503, dpi1_503, fpi0_503, \
                         fph_189, fph_378, fpi1_503, gsh_126, gpg0_270, gpg1_270, \
                         gph_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_20 * dpi0_503[k]
                   - f_21 * dpi1_503[k]
                   + pa_x[k] * fpi0_503[k]
                   - f_11 * pc_x[k] * fpi1_503[k];

        t_504[k] = f_1 * fph_378[k]
                   + f_1 * gsh_126[k]
                   + f_2 * gpg0_270[k]
                   - f_3 * gpg1_270[k]
                   + f_4 * pc_x[k] * gph_378[k];

        t_505[k] = f_13 * fph_189[k]
                   + f_4 * pc_y[k] * gph_378[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, pc_x, pc_z, fph_381, gsh_129, gpg0_270, \
                         gpg0_273, gpg1_270, gpg1_273, gph_378, gph_379, gph_380, \
                         gph_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_4 * pc_z[k] * gph_378[k];

        t_507[k] = f_1 * fph_381[k]
                   + f_1 * gsh_129[k]
                   + f_9 * gpg0_273[k]
                   - f_10 * gpg1_273[k]
                   + f_4 * pc_x[k] * gph_381[k];

        t_508[k] = f_4 * pc_z[k] * gph_379[k];

        t_509[k] = f_5 * gpg0_270[k]
                   - f_6 * gpg1_270[k]
                   + f_4 * pc_z[k] * gph_380[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, pc_x, pc_y, pc_z, fph_194, fph_384, gsh_132, \
                         gpg0_276, gpg1_276, gph_381, gph_383, \
                         gph_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = f_1 * fph_384[k]
                   + f_1 * gsh_132[k]
                   + f_7 * gpg0_276[k]
                   - f_8 * gpg1_276[k]
                   + f_4 * pc_x[k] * gph_384[k];

        t_511[k] = f_4 * pc_z[k] * gph_381[k];

        t_512[k] = f_13 * fph_194[k]
                   + f_4 * pc_y[k] * gph_383[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pc_x, pc_z, fph_388, gsh_136, gpg0_272, \
                         gpg0_280, gpg1_272, gpg1_280, gph_383, gph_384, \
                         gph_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_7 * gpg0_272[k]
                   - f_8 * gpg1_272[k]
                   + f_4 * pc_z[k] * gph_383[k];

        t_514[k] = f_1 * fph_388[k]
                   + f_1 * gsh_136[k]
                   + f_5 * gpg0_280[k]
                   - f_6 * gpg1_280[k]
                   + f_4 * pc_x[k] * gph_388[k];

        t_515[k] = f_4 * pc_z[k] * gph_384[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pc_y, pc_z, fph_198, gpg0_273, gpg0_275, \
                         gpg1_273, gpg1_275, gph_385, gph_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_5 * gpg0_273[k]
                   - f_6 * gpg1_273[k]
                   + f_4 * pc_z[k] * gph_385[k];

        t_517[k] = f_13 * fph_198[k]
                   + f_4 * pc_y[k] * gph_387[k];

        t_518[k] = f_9 * gpg0_275[k]
                   - f_10 * gpg1_275[k]
                   + f_4 * pc_z[k] * gph_387[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, pc_x, pc_z, fph_393, fph_395, fph_396, \
                         gsh_141, gsh_143, gsh_144, gph_388, gph_393, gph_395, \
                         gph_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_1 * fph_393[k]
                   + f_1 * gsh_141[k]
                   + f_4 * pc_x[k] * gph_393[k];

        t_520[k] = f_4 * pc_z[k] * gph_388[k];

        t_521[k] = f_1 * fph_395[k]
                   + f_1 * gsh_143[k]
                   + f_4 * pc_x[k] * gph_395[k];

        t_522[k] = f_1 * fph_396[k]
                   + f_1 * gsh_144[k]
                   + f_4 * pc_x[k] * gph_396[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, fph_204, fph_397, fph_398, gsh_145, \
                         gsh_146, gpg0_280, gpg1_280, gph_393, gph_397, \
                         gph_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_1 * fph_397[k]
                   + f_1 * gsh_145[k]
                   + f_4 * pc_x[k] * gph_397[k];

        t_524[k] = f_1 * fph_398[k]
                   + f_1 * gsh_146[k]
                   + f_4 * pc_x[k] * gph_398[k];

        t_525[k] = f_13 * fph_204[k]
                   + f_2 * gpg0_280[k]
                   - f_3 * gpg1_280[k]
                   + f_4 * pc_y[k] * gph_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_z, gpg0_280, gpg0_281, gpg0_282, \
                         gpg1_280, gpg1_281, gpg1_282, gph_393, gph_394, gph_395, \
                         gph_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_4 * pc_z[k] * gph_393[k];

        t_527[k] = f_5 * gpg0_280[k]
                   - f_6 * gpg1_280[k]
                   + f_4 * pc_z[k] * gph_394[k];

        t_528[k] = f_7 * gpg0_281[k]
                   - f_8 * gpg1_281[k]
                   + f_4 * pc_z[k] * gph_395[k];

        t_529[k] = f_9 * gpg0_282[k]
                   - f_10 * gpg1_282[k]
                   + f_4 * pc_z[k] * gph_396[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pa_x, pc_x, pc_y, pc_z, fpi0_532, fph_209, \
                         fph_399, fpi1_532, gpg0_284, gpg1_284, \
                         gph_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_13 * fph_209[k]
                   + f_4 * pc_y[k] * gph_398[k];

        t_531[k] = f_2 * gpg0_284[k]
                   - f_3 * gpg1_284[k]
                   + f_4 * pc_z[k] * gph_398[k];

        t_532[k] = pa_x[k] * fpi0_532[k]
                   + f_14 * fph_399[k]
                   - f_11 * pc_x[k] * fpi1_532[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pa_x, pc_x, pc_y, pc_z, fpi0_535, \
                         fph_210, fph_402, fpi1_535, gsh_126, gph_399, \
                         gph_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_13 * fph_210[k]
                   + f_1 * gsh_126[k]
                   + f_4 * pc_y[k] * gph_399[k];

        t_534[k] = f_4 * pc_z[k] * gph_399[k];

        t_535[k] = pa_x[k] * fpi0_535[k]
                   + f_0 * fph_402[k]
                   - f_11 * pc_x[k] * fpi1_535[k];

        t_536[k] = f_4 * pc_z[k] * gph_400[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pa_x, pc_x, pc_z, fpi0_538, fph_405, fpi1_538, \
                         gpg0_285, gpg1_285, gph_401, gph_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_5 * gpg0_285[k]
                   - f_6 * gpg1_285[k]
                   + f_4 * pc_z[k] * gph_401[k];

        t_538[k] = pa_x[k] * fpi0_538[k]
                   + f_13 * fph_405[k]
                   - f_11 * pc_x[k] * fpi1_538[k];

        t_539[k] = f_4 * pc_z[k] * gph_402[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pa_x, pc_x, pc_y, pc_z, fpi0_542, fph_215, \
                         fph_409, fpi1_542, gsh_131, gpg0_287, gpg1_287, \
                         gph_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_13 * fph_215[k]
                   + f_1 * gsh_131[k]
                   + f_4 * pc_y[k] * gph_404[k];

        t_541[k] = f_7 * gpg0_287[k]
                   - f_8 * gpg1_287[k]
                   + f_4 * pc_z[k] * gph_404[k];

        t_542[k] = pa_x[k] * fpi0_542[k]
                   + f_12 * fph_409[k]
                   - f_11 * pc_x[k] * fpi1_542[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, pc_y, pc_z, fph_219, gsh_135, gpg0_288, \
                         gpg0_290, gpg1_288, gpg1_290, gph_405, gph_406, \
                         gph_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_4 * pc_z[k] * gph_405[k];

        t_544[k] = f_5 * gpg0_288[k]
                   - f_6 * gpg1_288[k]
                   + f_4 * pc_z[k] * gph_406[k];

        t_545[k] = f_13 * fph_219[k]
                   + f_1 * gsh_135[k]
                   + f_4 * pc_y[k] * gph_408[k];

        t_546[k] = f_9 * gpg0_290[k]
                   - f_10 * gpg1_290[k]
                   + f_4 * pc_z[k] * gph_408[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, pc_x, pc_z, fph_414, fph_416, \
                         fph_417, fph_418, gph_409, gph_414, gph_416, gph_417, \
                         gph_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_1 * fph_414[k]
                   + f_4 * pc_x[k] * gph_414[k];

        t_548[k] = f_4 * pc_z[k] * gph_409[k];

        t_549[k] = f_1 * fph_416[k]
                   + f_4 * pc_x[k] * gph_416[k];

        t_550[k] = f_1 * fph_417[k]
                   + f_4 * pc_x[k] * gph_417[k];

        t_551[k] = f_1 * fph_418[k]
                   + f_4 * pc_x[k] * gph_418[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pa_x, pc_x, pc_z, fpi0_553, fpi0_555, \
                         fph_419, fpi1_553, fpi1_555, gph_414, \
                         gph_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_1 * fph_419[k]
                   + f_4 * pc_x[k] * gph_419[k];

        t_553[k] = pa_x[k] * fpi0_553[k]
                   - f_11 * pc_x[k] * fpi1_553[k];

        t_554[k] = f_4 * pc_z[k] * gph_414[k];

        t_555[k] = pa_x[k] * fpi0_555[k]
                   - f_11 * pc_x[k] * fpi1_555[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, pa_x, pc_x, fpi0_556, fpi0_557, fpi0_558, \
                         fpi0_559, fpi1_556, fpi1_557, fpi1_558, \
                         fpi1_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = pa_x[k] * fpi0_556[k]
                   - f_11 * pc_x[k] * fpi1_556[k];

        t_557[k] = pa_x[k] * fpi0_557[k]
                   - f_11 * pc_x[k] * fpi1_557[k];

        t_558[k] = pa_x[k] * fpi0_558[k]
                   - f_11 * pc_x[k] * fpi1_558[k];

        t_559[k] = pa_x[k] * fpi0_559[k]
                   - f_11 * pc_x[k] * fpi1_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, pb_z, pc_y, pc_z, fph_231, gsi0_168, \
                         gsi0_171, gsh_126, gsi1_168, gsi1_171, \
                         gph_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = pb_z[k] * gsi0_168[k]
                   - f_11 * pc_z[k] * gsi1_168[k];

        t_561[k] = f_13 * fph_231[k]
                   + f_4 * pc_y[k] * gph_420[k];

        t_562[k] = f_1 * gsh_126[k]
                   + f_4 * pc_z[k] * gph_420[k];

        t_563[k] = pb_z[k] * gsi0_171[k]
                   - f_11 * pc_z[k] * gsi1_171[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_x, pb_z, pc_x, pc_z, fpi0_565, fph_425, \
                         fpi1_565, gsi0_174, gsh_127, gsi1_174, \
                         gph_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_1 * gsh_127[k]
                   + f_4 * pc_z[k] * gph_421[k];

        t_565[k] = pa_x[k] * fpi0_565[k]
                   + f_0 * fph_425[k]
                   - f_11 * pc_x[k] * fpi1_565[k];

        t_566[k] = pb_z[k] * gsi0_174[k]
                   - f_11 * pc_z[k] * gsi1_174[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pa_x, pc_x, pc_y, pc_z, fpi0_569, fph_236, \
                         fph_429, fpi1_569, gsh_129, gph_423, gph_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_1 * gsh_129[k]
                   + f_4 * pc_z[k] * gph_423[k];

        t_568[k] = f_13 * fph_236[k]
                   + f_4 * pc_y[k] * gph_425[k];

        t_569[k] = pa_x[k] * fpi0_569[k]
                   + f_13 * fph_429[k]
                   - f_11 * pc_x[k] * fpi1_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pa_x, pb_z, pc_x, pc_z, fpi0_572, fph_432, \
                         fpi1_572, gsi0_178, gsh_132, gsi1_178, \
                         gph_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = pb_z[k] * gsi0_178[k]
                   - f_11 * pc_z[k] * gsi1_178[k];

        t_571[k] = f_1 * gsh_132[k]
                   + f_4 * pc_z[k] * gph_426[k];

        t_572[k] = pa_x[k] * fpi0_572[k]
                   + f_12 * fph_432[k]
                   - f_11 * pc_x[k] * fpi1_572[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, pa_x, pc_x, pc_y, fpi0_574, fph_240, fph_434, \
                         fph_435, fpi1_574, gph_429, gph_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_13 * fph_240[k]
                   + f_4 * pc_y[k] * gph_429[k];

        t_574[k] = pa_x[k] * fpi0_574[k]
                   + f_12 * fph_434[k]
                   - f_11 * pc_x[k] * fpi1_574[k];

        t_575[k] = f_1 * fph_435[k]
                   + f_4 * pc_x[k] * gph_435[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pc_x, pc_z, fph_437, fph_438, fph_439, \
                         gsh_136, gph_430, gph_437, gph_438, gph_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_1 * gsh_136[k]
                   + f_4 * pc_z[k] * gph_430[k];

        t_577[k] = f_1 * fph_437[k]
                   + f_4 * pc_x[k] * gph_437[k];

        t_578[k] = f_1 * fph_438[k]
                   + f_4 * pc_x[k] * gph_438[k];

        t_579[k] = f_1 * fph_439[k]
                   + f_4 * pc_x[k] * gph_439[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pa_x, pc_x, pc_z, fpi0_581, fpi0_583, \
                         fph_440, fpi1_581, fpi1_583, gsh_141, gph_435, \
                         gph_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_1 * fph_440[k]
                   + f_4 * pc_x[k] * gph_440[k];

        t_581[k] = pa_x[k] * fpi0_581[k]
                   - f_11 * pc_x[k] * fpi1_581[k];

        t_582[k] = f_1 * gsh_141[k]
                   + f_4 * pc_z[k] * gph_435[k];

        t_583[k] = pa_x[k] * fpi0_583[k]
                   - f_11 * pc_x[k] * fpi1_583[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pa_x, pc_x, pc_y, fpi0_584, fpi0_585, \
                         fpi0_587, fph_251, fpi1_584, fpi1_585, fpi1_587, \
                         gph_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = pa_x[k] * fpi0_584[k]
                   - f_11 * pc_x[k] * fpi1_584[k];

        t_585[k] = pa_x[k] * fpi0_585[k]
                   - f_11 * pc_x[k] * fpi1_585[k];

        t_586[k] = f_13 * fph_251[k]
                   + f_4 * pc_y[k] * gph_440[k];

        t_587[k] = pa_x[k] * fpi0_587[k]
                   - f_11 * pc_x[k] * fpi1_587[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_z, pc_y, pc_z, fpi0_252, fpi0_255, \
                         fph_189, fph_252, fpi1_252, fpi1_255, \
                         gph_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_z[k] * fpi0_252[k]
                   - f_11 * pc_z[k] * fpi1_252[k];

        t_589[k] = f_12 * fph_252[k]
                   + f_4 * pc_y[k] * gph_441[k];

        t_590[k] = f_1 * fph_189[k]
                   + f_4 * pc_z[k] * gph_441[k];

        t_591[k] = pa_z[k] * fpi0_255[k]
                   - f_11 * pc_z[k] * fpi1_255[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fpi0,
                                                          const size_t fph, const size_t fpi1,
                                                          const size_t gsh, const size_t gpg0,
                                                          const size_t gpg1, const size_t gph,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fpi0_258 = buffer.data(fpi0 + 258);
    const auto *fpi0_262 = buffer.data(fpi0 + 262);
    const auto *fpi0_267 = buffer.data(fpi0 + 267);
    const auto *fpi0_280 = buffer.data(fpi0 + 280);
    const auto *fpi0_281 = buffer.data(fpi0 + 281);
    const auto *fpi0_283 = buffer.data(fpi0 + 283);
    const auto *fpi0_286 = buffer.data(fpi0 + 286);
    const auto *fpi0_290 = buffer.data(fpi0 + 290);
    const auto *fpi0_420 = buffer.data(fpi0 + 420);
    const auto *fpi0_425 = buffer.data(fpi0 + 425);
    const auto *fpi0_429 = buffer.data(fpi0 + 429);
    const auto *fpi0_434 = buffer.data(fpi0 + 434);
    const auto *fpi0_440 = buffer.data(fpi0 + 440);
    const auto *fpi0_621 = buffer.data(fpi0 + 621);
    const auto *fpi0_625 = buffer.data(fpi0 + 625);
    const auto *fpi0_628 = buffer.data(fpi0 + 628);
    const auto *fpi0_630 = buffer.data(fpi0 + 630);
    const auto *fpi0_637 = buffer.data(fpi0 + 637);
    const auto *fpi0_639 = buffer.data(fpi0 + 639);
    const auto *fpi0_640 = buffer.data(fpi0 + 640);
    const auto *fpi0_641 = buffer.data(fpi0 + 641);
    const auto *fpi0_642 = buffer.data(fpi0 + 642);
    const auto *fpi0_643 = buffer.data(fpi0 + 643);
    const auto *fpi0_649 = buffer.data(fpi0 + 649);
    const auto *fpi0_653 = buffer.data(fpi0 + 653);
    const auto *fpi0_656 = buffer.data(fpi0 + 656);
    const auto *fpi0_658 = buffer.data(fpi0 + 658);
    const auto *fpi0_665 = buffer.data(fpi0 + 665);
    const auto *fpi0_666 = buffer.data(fpi0 + 666);
    const auto *fpi0_667 = buffer.data(fpi0 + 667);
    const auto *fpi0_668 = buffer.data(fpi0 + 668);
    const auto *fpi0_669 = buffer.data(fpi0 + 669);
    const auto *fpi0_671 = buffer.data(fpi0 + 671);

    const auto *fph_192 = buffer.data(fph + 192);
    const auto *fph_195 = buffer.data(fph + 195);
    const auto *fph_204 = buffer.data(fph + 204);
    const auto *fph_209 = buffer.data(fph + 209);
    const auto *fph_210 = buffer.data(fph + 210);
    const auto *fph_213 = buffer.data(fph + 213);
    const auto *fph_216 = buffer.data(fph + 216);
    const auto *fph_225 = buffer.data(fph + 225);
    const auto *fph_231 = buffer.data(fph + 231);
    const auto *fph_234 = buffer.data(fph + 234);
    const auto *fph_237 = buffer.data(fph + 237);
    const auto *fph_252 = buffer.data(fph + 252);
    const auto *fph_254 = buffer.data(fph + 254);
    const auto *fph_255 = buffer.data(fph + 255);
    const auto *fph_257 = buffer.data(fph + 257);
    const auto *fph_258 = buffer.data(fph + 258);
    const auto *fph_261 = buffer.data(fph + 261);
    const auto *fph_267 = buffer.data(fph + 267);
    const auto *fph_269 = buffer.data(fph + 269);
    const auto *fph_270 = buffer.data(fph + 270);
    const auto *fph_271 = buffer.data(fph + 271);
    const auto *fph_272 = buffer.data(fph + 272);
    const auto *fph_275 = buffer.data(fph + 275);
    const auto *fph_278 = buffer.data(fph + 278);
    const auto *fph_282 = buffer.data(fph + 282);
    const auto *fph_294 = buffer.data(fph + 294);
    const auto *fph_296 = buffer.data(fph + 296);
    const auto *fph_299 = buffer.data(fph + 299);
    const auto *fph_303 = buffer.data(fph + 303);
    const auto *fph_314 = buffer.data(fph + 314);
    const auto *fph_315 = buffer.data(fph + 315);
    const auto *fph_317 = buffer.data(fph + 317);
    const auto *fph_320 = buffer.data(fph + 320);
    const auto *fph_324 = buffer.data(fph + 324);
    const auto *fph_330 = buffer.data(fph + 330);
    const auto *fph_332 = buffer.data(fph + 332);
    const auto *fph_333 = buffer.data(fph + 333);
    const auto *fph_334 = buffer.data(fph + 334);
    const auto *fph_335 = buffer.data(fph + 335);
    const auto *fph_446 = buffer.data(fph + 446);
    const auto *fph_450 = buffer.data(fph + 450);
    const auto *fph_453 = buffer.data(fph + 453);
    const auto *fph_455 = buffer.data(fph + 455);
    const auto *fph_457 = buffer.data(fph + 457);
    const auto *fph_458 = buffer.data(fph + 458);
    const auto *fph_459 = buffer.data(fph + 459);
    const auto *fph_460 = buffer.data(fph + 460);
    const auto *fph_461 = buffer.data(fph + 461);
    const auto *fph_467 = buffer.data(fph + 467);
    const auto *fph_471 = buffer.data(fph + 471);
    const auto *fph_474 = buffer.data(fph + 474);
    const auto *fph_476 = buffer.data(fph + 476);
    const auto *fph_477 = buffer.data(fph + 477);
    const auto *fph_478 = buffer.data(fph + 478);
    const auto *fph_479 = buffer.data(fph + 479);
    const auto *fph_480 = buffer.data(fph + 480);
    const auto *fph_481 = buffer.data(fph + 481);
    const auto *fph_482 = buffer.data(fph + 482);
    const auto *fph_483 = buffer.data(fph + 483);
    const auto *fph_486 = buffer.data(fph + 486);
    const auto *fph_488 = buffer.data(fph + 488);
    const auto *fph_489 = buffer.data(fph + 489);
    const auto *fph_492 = buffer.data(fph + 492);
    const auto *fph_493 = buffer.data(fph + 493);
    const auto *fph_495 = buffer.data(fph + 495);
    const auto *fph_497 = buffer.data(fph + 497);
    const auto *fph_498 = buffer.data(fph + 498);
    const auto *fph_499 = buffer.data(fph + 499);
    const auto *fph_500 = buffer.data(fph + 500);
    const auto *fph_501 = buffer.data(fph + 501);
    const auto *fph_502 = buffer.data(fph + 502);
    const auto *fph_503 = buffer.data(fph + 503);
    const auto *fph_507 = buffer.data(fph + 507);
    const auto *fph_510 = buffer.data(fph + 510);
    const auto *fph_514 = buffer.data(fph + 514);
    const auto *fph_516 = buffer.data(fph + 516);
    const auto *fph_519 = buffer.data(fph + 519);
    const auto *fph_520 = buffer.data(fph + 520);
    const auto *fph_521 = buffer.data(fph + 521);
    const auto *fph_522 = buffer.data(fph + 522);
    const auto *fph_523 = buffer.data(fph + 523);
    const auto *fph_525 = buffer.data(fph + 525);

    const auto *fpi1_258 = buffer.data(fpi1 + 258);
    const auto *fpi1_262 = buffer.data(fpi1 + 262);
    const auto *fpi1_267 = buffer.data(fpi1 + 267);
    const auto *fpi1_280 = buffer.data(fpi1 + 280);
    const auto *fpi1_281 = buffer.data(fpi1 + 281);
    const auto *fpi1_283 = buffer.data(fpi1 + 283);
    const auto *fpi1_286 = buffer.data(fpi1 + 286);
    const auto *fpi1_290 = buffer.data(fpi1 + 290);
    const auto *fpi1_420 = buffer.data(fpi1 + 420);
    const auto *fpi1_425 = buffer.data(fpi1 + 425);
    const auto *fpi1_429 = buffer.data(fpi1 + 429);
    const auto *fpi1_434 = buffer.data(fpi1 + 434);
    const auto *fpi1_440 = buffer.data(fpi1 + 440);
    const auto *fpi1_621 = buffer.data(fpi1 + 621);
    const auto *fpi1_625 = buffer.data(fpi1 + 625);
    const auto *fpi1_628 = buffer.data(fpi1 + 628);
    const auto *fpi1_630 = buffer.data(fpi1 + 630);
    const auto *fpi1_637 = buffer.data(fpi1 + 637);
    const auto *fpi1_639 = buffer.data(fpi1 + 639);
    const auto *fpi1_640 = buffer.data(fpi1 + 640);
    const auto *fpi1_641 = buffer.data(fpi1 + 641);
    const auto *fpi1_642 = buffer.data(fpi1 + 642);
    const auto *fpi1_643 = buffer.data(fpi1 + 643);
    const auto *fpi1_649 = buffer.data(fpi1 + 649);
    const auto *fpi1_653 = buffer.data(fpi1 + 653);
    const auto *fpi1_656 = buffer.data(fpi1 + 656);
    const auto *fpi1_658 = buffer.data(fpi1 + 658);
    const auto *fpi1_665 = buffer.data(fpi1 + 665);
    const auto *fpi1_666 = buffer.data(fpi1 + 666);
    const auto *fpi1_667 = buffer.data(fpi1 + 667);
    const auto *fpi1_668 = buffer.data(fpi1 + 668);
    const auto *fpi1_669 = buffer.data(fpi1 + 669);
    const auto *fpi1_671 = buffer.data(fpi1 + 671);

    const auto *gsh_147 = buffer.data(gsh + 147);
    const auto *gsh_149 = buffer.data(gsh + 149);
    const auto *gsh_150 = buffer.data(gsh + 150);
    const auto *gsh_152 = buffer.data(gsh + 152);
    const auto *gsh_153 = buffer.data(gsh + 153);
    const auto *gsh_156 = buffer.data(gsh + 156);
    const auto *gsh_159 = buffer.data(gsh + 159);
    const auto *gsh_161 = buffer.data(gsh + 161);
    const auto *gsh_163 = buffer.data(gsh + 163);
    const auto *gsh_164 = buffer.data(gsh + 164);
    const auto *gsh_165 = buffer.data(gsh + 165);
    const auto *gsh_166 = buffer.data(gsh + 166);
    const auto *gsh_167 = buffer.data(gsh + 167);
    const auto *gsh_171 = buffer.data(gsh + 171);
    const auto *gsh_174 = buffer.data(gsh + 174);
    const auto *gsh_178 = buffer.data(gsh + 178);
    const auto *gsh_180 = buffer.data(gsh + 180);
    const auto *gsh_183 = buffer.data(gsh + 183);
    const auto *gsh_184 = buffer.data(gsh + 184);
    const auto *gsh_185 = buffer.data(gsh + 185);
    const auto *gsh_186 = buffer.data(gsh + 186);
    const auto *gsh_187 = buffer.data(gsh + 187);

    const auto *gpg0_320 = buffer.data(gpg0 + 320);
    const auto *gpg0_324 = buffer.data(gpg0 + 324);
    const auto *gpg0_325 = buffer.data(gpg0 + 325);
    const auto *gpg0_327 = buffer.data(gpg0 + 327);
    const auto *gpg0_328 = buffer.data(gpg0 + 328);
    const auto *gpg0_329 = buffer.data(gpg0 + 329);
    const auto *gpg0_345 = buffer.data(gpg0 + 345);
    const auto *gpg0_348 = buffer.data(gpg0 + 348);
    const auto *gpg0_351 = buffer.data(gpg0 + 351);
    const auto *gpg0_355 = buffer.data(gpg0 + 355);
    const auto *gpg0_363 = buffer.data(gpg0 + 363);
    const auto *gpg0_366 = buffer.data(gpg0 + 366);
    const auto *gpg0_370 = buffer.data(gpg0 + 370);
    const auto *gpg0_372 = buffer.data(gpg0 + 372);
    const auto *gpg0_373 = buffer.data(gpg0 + 373);
    const auto *gpg0_374 = buffer.data(gpg0 + 374);
    const auto *gpg0_375 = buffer.data(gpg0 + 375);

    const auto *gpg1_320 = buffer.data(gpg1 + 320);
    const auto *gpg1_324 = buffer.data(gpg1 + 324);
    const auto *gpg1_325 = buffer.data(gpg1 + 325);
    const auto *gpg1_327 = buffer.data(gpg1 + 327);
    const auto *gpg1_328 = buffer.data(gpg1 + 328);
    const auto *gpg1_329 = buffer.data(gpg1 + 329);
    const auto *gpg1_345 = buffer.data(gpg1 + 345);
    const auto *gpg1_348 = buffer.data(gpg1 + 348);
    const auto *gpg1_351 = buffer.data(gpg1 + 351);
    const auto *gpg1_355 = buffer.data(gpg1 + 355);
    const auto *gpg1_363 = buffer.data(gpg1 + 363);
    const auto *gpg1_366 = buffer.data(gpg1 + 366);
    const auto *gpg1_370 = buffer.data(gpg1 + 370);
    const auto *gpg1_372 = buffer.data(gpg1 + 372);
    const auto *gpg1_373 = buffer.data(gpg1 + 373);
    const auto *gpg1_374 = buffer.data(gpg1 + 374);
    const auto *gpg1_375 = buffer.data(gpg1 + 375);

    const auto *gph_443 = buffer.data(gph + 443);
    const auto *gph_444 = buffer.data(gph + 444);
    const auto *gph_446 = buffer.data(gph + 446);
    const auto *gph_447 = buffer.data(gph + 447);
    const auto *gph_450 = buffer.data(gph + 450);
    const auto *gph_453 = buffer.data(gph + 453);
    const auto *gph_455 = buffer.data(gph + 455);
    const auto *gph_456 = buffer.data(gph + 456);
    const auto *gph_457 = buffer.data(gph + 457);
    const auto *gph_458 = buffer.data(gph + 458);
    const auto *gph_459 = buffer.data(gph + 459);
    const auto *gph_460 = buffer.data(gph + 460);
    const auto *gph_461 = buffer.data(gph + 461);
    const auto *gph_462 = buffer.data(gph + 462);
    const auto *gph_464 = buffer.data(gph + 464);
    const auto *gph_465 = buffer.data(gph + 465);
    const auto *gph_467 = buffer.data(gph + 467);
    const auto *gph_468 = buffer.data(gph + 468);
    const auto *gph_471 = buffer.data(gph + 471);
    const auto *gph_477 = buffer.data(gph + 477);
    const auto *gph_478 = buffer.data(gph + 478);
    const auto *gph_479 = buffer.data(gph + 479);
    const auto *gph_480 = buffer.data(gph + 480);
    const auto *gph_481 = buffer.data(gph + 481);
    const auto *gph_482 = buffer.data(gph + 482);
    const auto *gph_483 = buffer.data(gph + 483);
    const auto *gph_485 = buffer.data(gph + 485);
    const auto *gph_486 = buffer.data(gph + 486);
    const auto *gph_488 = buffer.data(gph + 488);
    const auto *gph_489 = buffer.data(gph + 489);
    const auto *gph_492 = buffer.data(gph + 492);
    const auto *gph_493 = buffer.data(gph + 493);
    const auto *gph_498 = buffer.data(gph + 498);
    const auto *gph_499 = buffer.data(gph + 499);
    const auto *gph_500 = buffer.data(gph + 500);
    const auto *gph_501 = buffer.data(gph + 501);
    const auto *gph_502 = buffer.data(gph + 502);
    const auto *gph_503 = buffer.data(gph + 503);
    const auto *gph_504 = buffer.data(gph + 504);
    const auto *gph_506 = buffer.data(gph + 506);
    const auto *gph_507 = buffer.data(gph + 507);
    const auto *gph_509 = buffer.data(gph + 509);
    const auto *gph_510 = buffer.data(gph + 510);
    const auto *gph_513 = buffer.data(gph + 513);
    const auto *gph_514 = buffer.data(gph + 514);
    const auto *gph_516 = buffer.data(gph + 516);
    const auto *gph_519 = buffer.data(gph + 519);
    const auto *gph_520 = buffer.data(gph + 520);
    const auto *gph_521 = buffer.data(gph + 521);
    const auto *gph_522 = buffer.data(gph + 522);
    const auto *gph_523 = buffer.data(gph + 523);
    const auto *gph_524 = buffer.data(gph + 524);
    const auto *gph_525 = buffer.data(gph + 525);

#pragma omp simd aligned(t_592, t_593, t_594, pa_z, pc_x, pc_y, pc_z, fpi0_258, fph_254, \
                         fph_446, fpi1_258, gsh_152, gpg0_320, gpg1_320, gph_443, \
                         gph_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_12 * fph_254[k]
                   + f_4 * pc_y[k] * gph_443[k];

        t_593[k] = f_1 * fph_446[k]
                   + f_1 * gsh_152[k]
                   + f_9 * gpg0_320[k]
                   - f_10 * gpg1_320[k]
                   + f_4 * pc_x[k] * gph_446[k];

        t_594[k] = pa_z[k] * fpi0_258[k]
                   - f_11 * pc_z[k] * fpi1_258[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pc_x, pc_y, pc_z, fph_192, fph_257, fph_450, \
                         gsh_156, gpg0_324, gpg1_324, gph_444, gph_446, \
                         gph_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_1 * fph_192[k]
                   + f_4 * pc_z[k] * gph_444[k];

        t_596[k] = f_12 * fph_257[k]
                   + f_4 * pc_y[k] * gph_446[k];

        t_597[k] = f_1 * fph_450[k]
                   + f_1 * gsh_156[k]
                   + f_7 * gpg0_324[k]
                   - f_8 * gpg1_324[k]
                   + f_4 * pc_x[k] * gph_450[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pa_z, pc_x, pc_z, fpi0_262, fph_195, fph_453, \
                         fpi1_262, gsh_159, gpg0_327, gpg1_327, gph_447, \
                         gph_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = pa_z[k] * fpi0_262[k]
                   - f_11 * pc_z[k] * fpi1_262[k];

        t_599[k] = f_1 * fph_195[k]
                   + f_4 * pc_z[k] * gph_447[k];

        t_600[k] = f_1 * fph_453[k]
                   + f_1 * gsh_159[k]
                   + f_5 * gpg0_327[k]
                   - f_6 * gpg1_327[k]
                   + f_4 * pc_x[k] * gph_453[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, pa_z, pc_x, pc_y, pc_z, fpi0_267, fph_261, \
                         fph_455, fpi1_267, gsh_161, gpg0_329, gpg1_329, gph_450, \
                         gph_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_12 * fph_261[k]
                   + f_4 * pc_y[k] * gph_450[k];

        t_602[k] = f_1 * fph_455[k]
                   + f_1 * gsh_161[k]
                   + f_5 * gpg0_329[k]
                   - f_6 * gpg1_329[k]
                   + f_4 * pc_x[k] * gph_455[k];

        t_603[k] = pa_z[k] * fpi0_267[k]
                   - f_11 * pc_z[k] * fpi1_267[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pc_x, fph_457, fph_458, fph_459, gsh_163, \
                         gsh_164, gsh_165, gph_457, gph_458, gph_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_1 * fph_457[k]
                   + f_1 * gsh_163[k]
                   + f_4 * pc_x[k] * gph_457[k];

        t_605[k] = f_1 * fph_458[k]
                   + f_1 * gsh_164[k]
                   + f_4 * pc_x[k] * gph_458[k];

        t_606[k] = f_1 * fph_459[k]
                   + f_1 * gsh_165[k]
                   + f_4 * pc_x[k] * gph_459[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_x, pc_y, fph_267, fph_460, fph_461, gsh_166, \
                         gsh_167, gpg0_325, gpg1_325, gph_456, gph_460, \
                         gph_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_1 * fph_460[k]
                   + f_1 * gsh_166[k]
                   + f_4 * pc_x[k] * gph_460[k];

        t_608[k] = f_1 * fph_461[k]
                   + f_1 * gsh_167[k]
                   + f_4 * pc_x[k] * gph_461[k];

        t_609[k] = f_12 * fph_267[k]
                   + f_2 * gpg0_325[k]
                   - f_3 * gpg1_325[k]
                   + f_4 * pc_y[k] * gph_456[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_y, pc_z, fph_204, fph_269, fph_270, gpg0_327, \
                         gpg0_328, gpg1_327, gpg1_328, gph_456, gph_458, \
                         gph_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_1 * fph_204[k]
                   + f_4 * pc_z[k] * gph_456[k];

        t_611[k] = f_12 * fph_269[k]
                   + f_9 * gpg0_327[k]
                   - f_10 * gpg1_327[k]
                   + f_4 * pc_y[k] * gph_458[k];

        t_612[k] = f_12 * fph_270[k]
                   + f_7 * gpg0_328[k]
                   - f_8 * gpg1_328[k]
                   + f_4 * pc_y[k] * gph_459[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_z, pc_y, pc_z, fpi0_280, fph_209, \
                         fph_271, fph_272, fpi1_280, gpg0_329, gpg1_329, gph_460, \
                         gph_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_12 * fph_271[k]
                   + f_5 * gpg0_329[k]
                   - f_6 * gpg1_329[k]
                   + f_4 * pc_y[k] * gph_460[k];

        t_614[k] = f_12 * fph_272[k]
                   + f_4 * pc_y[k] * gph_461[k];

        t_615[k] = f_1 * fph_209[k]
                   + f_2 * gpg0_329[k]
                   - f_3 * gpg1_329[k]
                   + f_4 * pc_z[k] * gph_461[k];

        t_616[k] = pa_z[k] * fpi0_280[k]
                   - f_11 * pc_z[k] * fpi1_280[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pa_z, pc_y, pc_z, fpi0_281, fpi0_283, \
                         fph_210, fph_275, fpi1_281, fpi1_283, gsh_149, gph_462, \
                         gph_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = pa_z[k] * fpi0_281[k]
                   - f_11 * pc_z[k] * fpi1_281[k];

        t_618[k] = f_1 * fph_210[k]
                   + f_4 * pc_z[k] * gph_462[k];

        t_619[k] = pa_z[k] * fpi0_283[k]
                   - f_11 * pc_z[k] * fpi1_283[k];

        t_620[k] = f_12 * fph_275[k]
                   + f_1 * gsh_149[k]
                   + f_4 * pc_y[k] * gph_464[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pa_x, pa_z, pc_x, pc_z, fpi0_286, fpi0_621, \
                         fph_213, fph_467, fpi1_286, fpi1_621, \
                         gph_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = pa_x[k] * fpi0_621[k]
                   + f_0 * fph_467[k]
                   - f_11 * pc_x[k] * fpi1_621[k];

        t_622[k] = pa_z[k] * fpi0_286[k]
                   - f_11 * pc_z[k] * fpi1_286[k];

        t_623[k] = f_1 * fph_213[k]
                   + f_4 * pc_z[k] * gph_465[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_x, pa_z, pc_x, pc_y, pc_z, fpi0_290, \
                         fpi0_625, fph_278, fph_471, fpi1_290, fpi1_625, gsh_152, \
                         gph_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_12 * fph_278[k]
                   + f_1 * gsh_152[k]
                   + f_4 * pc_y[k] * gph_467[k];

        t_625[k] = pa_x[k] * fpi0_625[k]
                   + f_13 * fph_471[k]
                   - f_11 * pc_x[k] * fpi1_625[k];

        t_626[k] = pa_z[k] * fpi0_290[k]
                   - f_11 * pc_z[k] * fpi1_290[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_x, pc_x, pc_y, pc_z, fpi0_628, fph_216, \
                         fph_282, fph_474, fpi1_628, gsh_156, gph_468, \
                         gph_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_1 * fph_216[k]
                   + f_4 * pc_z[k] * gph_468[k];

        t_628[k] = pa_x[k] * fpi0_628[k]
                   + f_12 * fph_474[k]
                   - f_11 * pc_x[k] * fpi1_628[k];

        t_629[k] = f_12 * fph_282[k]
                   + f_1 * gsh_156[k]
                   + f_4 * pc_y[k] * gph_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pa_x, pc_x, fpi0_630, fph_476, fph_477, \
                         fph_478, fph_479, fpi1_630, gph_477, gph_478, \
                         gph_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = pa_x[k] * fpi0_630[k]
                   + f_12 * fph_476[k]
                   - f_11 * pc_x[k] * fpi1_630[k];

        t_631[k] = f_1 * fph_477[k]
                   + f_4 * pc_x[k] * gph_477[k];

        t_632[k] = f_1 * fph_478[k]
                   + f_4 * pc_x[k] * gph_478[k];

        t_633[k] = f_1 * fph_479[k]
                   + f_4 * pc_x[k] * gph_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pa_x, pc_x, fpi0_637, fph_480, fph_481, \
                         fph_482, fpi1_637, gph_480, gph_481, gph_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_1 * fph_480[k]
                   + f_4 * pc_x[k] * gph_480[k];

        t_635[k] = f_1 * fph_481[k]
                   + f_4 * pc_x[k] * gph_481[k];

        t_636[k] = f_1 * fph_482[k]
                   + f_4 * pc_x[k] * gph_482[k];

        t_637[k] = pa_x[k] * fpi0_637[k]
                   - f_11 * pc_x[k] * fpi1_637[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pa_x, pc_x, pc_z, fpi0_639, fpi0_640, \
                         fpi0_641, fph_225, fpi1_639, fpi1_640, fpi1_641, \
                         gph_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_1 * fph_225[k]
                   + f_4 * pc_z[k] * gph_477[k];

        t_639[k] = pa_x[k] * fpi0_639[k]
                   - f_11 * pc_x[k] * fpi1_639[k];

        t_640[k] = pa_x[k] * fpi0_640[k]
                   - f_11 * pc_x[k] * fpi1_640[k];

        t_641[k] = pa_x[k] * fpi0_641[k]
                   - f_11 * pc_x[k] * fpi1_641[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pa_x, pc_x, pc_y, fpi0_642, fpi0_643, \
                         fph_294, fph_483, fpi1_642, fpi1_643, gpg0_345, gpg1_345, \
                         gph_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = pa_x[k] * fpi0_642[k]
                   - f_11 * pc_x[k] * fpi1_642[k];

        t_643[k] = pa_x[k] * fpi0_643[k]
                   - f_11 * pc_x[k] * fpi1_643[k];

        t_644[k] = f_1 * fph_483[k]
                   + f_2 * gpg0_345[k]
                   - f_3 * gpg1_345[k]
                   + f_4 * pc_x[k] * gph_483[k];

        t_645[k] = f_12 * fph_294[k]
                   + f_4 * pc_y[k] * gph_483[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, fph_231, fph_296, fph_486, \
                         gsh_147, gpg0_348, gpg1_348, gph_483, gph_485, \
                         gph_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_1 * fph_231[k]
                   + f_1 * gsh_147[k]
                   + f_4 * pc_z[k] * gph_483[k];

        t_647[k] = f_1 * fph_486[k]
                   + f_9 * gpg0_348[k]
                   - f_10 * gpg1_348[k]
                   + f_4 * pc_x[k] * gph_486[k];

        t_648[k] = f_12 * fph_296[k]
                   + f_4 * pc_y[k] * gph_485[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pa_x, pc_x, pc_z, fpi0_649, fph_234, fph_488, \
                         fph_489, fpi1_649, gsh_150, gpg0_351, gpg1_351, gph_486, \
                         gph_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = pa_x[k] * fpi0_649[k]
                   + f_0 * fph_488[k]
                   - f_11 * pc_x[k] * fpi1_649[k];

        t_650[k] = f_1 * fph_489[k]
                   + f_7 * gpg0_351[k]
                   - f_8 * gpg1_351[k]
                   + f_4 * pc_x[k] * gph_489[k];

        t_651[k] = f_1 * fph_234[k]
                   + f_1 * gsh_150[k]
                   + f_4 * pc_z[k] * gph_486[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, pa_x, pc_x, pc_y, fpi0_653, fph_299, fph_492, \
                         fph_493, fpi1_653, gpg0_355, gpg1_355, gph_488, \
                         gph_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_12 * fph_299[k]
                   + f_4 * pc_y[k] * gph_488[k];

        t_653[k] = pa_x[k] * fpi0_653[k]
                   + f_13 * fph_492[k]
                   - f_11 * pc_x[k] * fpi1_653[k];

        t_654[k] = f_1 * fph_493[k]
                   + f_5 * gpg0_355[k]
                   - f_6 * gpg1_355[k]
                   + f_4 * pc_x[k] * gph_493[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, pa_x, pc_x, pc_y, pc_z, fpi0_656, fph_237, \
                         fph_303, fph_495, fpi1_656, gsh_153, gph_489, \
                         gph_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_1 * fph_237[k]
                   + f_1 * gsh_153[k]
                   + f_4 * pc_z[k] * gph_489[k];

        t_656[k] = pa_x[k] * fpi0_656[k]
                   + f_12 * fph_495[k]
                   - f_11 * pc_x[k] * fpi1_656[k];

        t_657[k] = f_12 * fph_303[k]
                   + f_4 * pc_y[k] * gph_492[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pa_x, pc_x, fpi0_658, fph_497, fph_498, \
                         fph_499, fph_500, fpi1_658, gph_498, gph_499, \
                         gph_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = pa_x[k] * fpi0_658[k]
                   + f_12 * fph_497[k]
                   - f_11 * pc_x[k] * fpi1_658[k];

        t_659[k] = f_1 * fph_498[k]
                   + f_4 * pc_x[k] * gph_498[k];

        t_660[k] = f_1 * fph_499[k]
                   + f_4 * pc_x[k] * gph_499[k];

        t_661[k] = f_1 * fph_500[k]
                   + f_4 * pc_x[k] * gph_500[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, pa_x, pc_x, fpi0_665, fph_501, fph_502, \
                         fph_503, fpi1_665, gph_501, gph_502, gph_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_1 * fph_501[k]
                   + f_4 * pc_x[k] * gph_501[k];

        t_663[k] = f_1 * fph_502[k]
                   + f_4 * pc_x[k] * gph_502[k];

        t_664[k] = f_1 * fph_503[k]
                   + f_4 * pc_x[k] * gph_503[k];

        t_665[k] = pa_x[k] * fpi0_665[k]
                   - f_11 * pc_x[k] * fpi1_665[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pa_x, pc_x, fpi0_666, fpi0_667, fpi0_668, \
                         fpi0_669, fpi1_666, fpi1_667, fpi1_668, \
                         fpi1_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = pa_x[k] * fpi0_666[k]
                   - f_11 * pc_x[k] * fpi1_666[k];

        t_667[k] = pa_x[k] * fpi0_667[k]
                   - f_11 * pc_x[k] * fpi1_667[k];

        t_668[k] = pa_x[k] * fpi0_668[k]
                   - f_11 * pc_x[k] * fpi1_668[k];

        t_669[k] = pa_x[k] * fpi0_669[k]
                   - f_11 * pc_x[k] * fpi1_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pa_x, pa_y, pc_x, pc_y, fpi0_420, \
                         fpi0_671, fph_314, fph_315, fpi1_420, fpi1_671, gph_503, \
                         gph_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_12 * fph_314[k]
                   + f_4 * pc_y[k] * gph_503[k];

        t_671[k] = pa_x[k] * fpi0_671[k]
                   - f_11 * pc_x[k] * fpi1_671[k];

        t_672[k] = pa_y[k] * fpi0_420[k]
                   - f_11 * pc_y[k] * fpi1_420[k];

        t_673[k] = f_1 * fph_315[k]
                   + f_4 * pc_y[k] * gph_504[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, pc_x, pc_y, pc_z, fph_252, fph_317, fph_507, \
                         gsh_171, gpg0_363, gpg1_363, gph_504, gph_506, \
                         gph_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_12 * fph_252[k]
                   + f_4 * pc_z[k] * gph_504[k];

        t_675[k] = f_1 * fph_507[k]
                   + f_1 * gsh_171[k]
                   + f_9 * gpg0_363[k]
                   - f_10 * gpg1_363[k]
                   + f_4 * pc_x[k] * gph_507[k];

        t_676[k] = f_1 * fph_317[k]
                   + f_4 * pc_y[k] * gph_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pa_y, pc_x, pc_y, pc_z, fpi0_425, fph_255, \
                         fph_510, fpi1_425, gsh_174, gpg0_366, gpg1_366, gph_507, \
                         gph_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = pa_y[k] * fpi0_425[k]
                   - f_11 * pc_y[k] * fpi1_425[k];

        t_678[k] = f_1 * fph_510[k]
                   + f_1 * gsh_174[k]
                   + f_7 * gpg0_366[k]
                   - f_8 * gpg1_366[k]
                   + f_4 * pc_x[k] * gph_510[k];

        t_679[k] = f_12 * fph_255[k]
                   + f_4 * pc_z[k] * gph_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pa_y, pc_x, pc_y, fpi0_429, fph_320, fph_514, \
                         fpi1_429, gsh_178, gpg0_370, gpg1_370, gph_509, \
                         gph_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_1 * fph_320[k]
                   + f_4 * pc_y[k] * gph_509[k];

        t_681[k] = pa_y[k] * fpi0_429[k]
                   - f_11 * pc_y[k] * fpi1_429[k];

        t_682[k] = f_1 * fph_514[k]
                   + f_1 * gsh_178[k]
                   + f_5 * gpg0_370[k]
                   - f_6 * gpg1_370[k]
                   + f_4 * pc_x[k] * gph_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, fph_258, fph_324, fph_516, \
                         gsh_180, gpg0_372, gpg1_372, gph_510, gph_513, \
                         gph_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_12 * fph_258[k]
                   + f_4 * pc_z[k] * gph_510[k];

        t_684[k] = f_1 * fph_516[k]
                   + f_1 * gsh_180[k]
                   + f_5 * gpg0_372[k]
                   - f_6 * gpg1_372[k]
                   + f_4 * pc_x[k] * gph_516[k];

        t_685[k] = f_1 * fph_324[k]
                   + f_4 * pc_y[k] * gph_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, pa_y, pc_x, pc_y, fpi0_434, fph_519, fph_520, \
                         fpi1_434, gsh_183, gsh_184, gph_519, gph_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = pa_y[k] * fpi0_434[k]
                   - f_11 * pc_y[k] * fpi1_434[k];

        t_687[k] = f_1 * fph_519[k]
                   + f_1 * gsh_183[k]
                   + f_4 * pc_x[k] * gph_519[k];

        t_688[k] = f_1 * fph_520[k]
                   + f_1 * gsh_184[k]
                   + f_4 * pc_x[k] * gph_520[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, pc_x, fph_521, fph_522, fph_523, gsh_185, \
                         gsh_186, gsh_187, gph_521, gph_522, gph_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = f_1 * fph_521[k]
                   + f_1 * gsh_185[k]
                   + f_4 * pc_x[k] * gph_521[k];

        t_690[k] = f_1 * fph_522[k]
                   + f_1 * gsh_186[k]
                   + f_4 * pc_x[k] * gph_522[k];

        t_691[k] = f_1 * fph_523[k]
                   + f_1 * gsh_187[k]
                   + f_4 * pc_x[k] * gph_523[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, pa_y, pc_y, pc_z, fpi0_440, fph_267, fph_330, \
                         fpi1_440, gpg0_370, gpg1_370, gph_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = pa_y[k] * fpi0_440[k]
                   - f_11 * pc_y[k] * fpi1_440[k];

        t_693[k] = f_1 * fph_330[k]
                   + f_2 * gpg0_370[k]
                   - f_3 * gpg1_370[k]
                   + f_4 * pc_y[k] * gph_519[k];

        t_694[k] = f_12 * fph_267[k]
                   + f_4 * pc_z[k] * gph_519[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, pc_y, fph_332, fph_333, fph_334, gpg0_372, \
                         gpg0_373, gpg0_374, gpg1_372, gpg1_373, gpg1_374, gph_521, gph_522, \
                         gph_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_1 * fph_332[k]
                   + f_9 * gpg0_372[k]
                   - f_10 * gpg1_372[k]
                   + f_4 * pc_y[k] * gph_521[k];

        t_696[k] = f_1 * fph_333[k]
                   + f_7 * gpg0_373[k]
                   - f_8 * gpg1_373[k]
                   + f_4 * pc_y[k] * gph_522[k];

        t_697[k] = f_1 * fph_334[k]
                   + f_5 * gpg0_374[k]
                   - f_6 * gpg1_374[k]
                   + f_4 * pc_y[k] * gph_523[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, pc_y, pc_z, fph_272, fph_335, fph_525, \
                         gpg0_374, gpg0_375, gpg1_374, gpg1_375, gph_524, \
                         gph_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_1 * fph_335[k]
                   + f_4 * pc_y[k] * gph_524[k];

        t_699[k] = f_12 * fph_272[k]
                   + f_2 * gpg0_374[k]
                   - f_3 * gpg1_374[k]
                   + f_4 * pc_z[k] * gph_524[k];

        t_700[k] = f_1 * fph_525[k]
                   + f_2 * gpg0_375[k]
                   - f_3 * gpg1_375[k]
                   + f_4 * pc_x[k] * gph_525[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t fpi0, const size_t fph,
                                                          const size_t fpi1, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t gpg0, const size_t gpg1,
                                                          const size_t gph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_14 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fpi0_476 = buffer.data(fpi0 + 476);
    const auto *fpi0_478 = buffer.data(fpi0 + 478);
    const auto *fpi0_481 = buffer.data(fpi0 + 481);
    const auto *fpi0_485 = buffer.data(fpi0 + 485);
    const auto *fpi0_490 = buffer.data(fpi0 + 490);
    const auto *fpi0_703 = buffer.data(fpi0 + 703);
    const auto *fpi0_706 = buffer.data(fpi0 + 706);
    const auto *fpi0_710 = buffer.data(fpi0 + 710);
    const auto *fpi0_712 = buffer.data(fpi0 + 712);
    const auto *fpi0_721 = buffer.data(fpi0 + 721);
    const auto *fpi0_723 = buffer.data(fpi0 + 723);
    const auto *fpi0_724 = buffer.data(fpi0 + 724);
    const auto *fpi0_725 = buffer.data(fpi0 + 725);
    const auto *fpi0_726 = buffer.data(fpi0 + 726);
    const auto *fpi0_727 = buffer.data(fpi0 + 727);
    const auto *fpi0_731 = buffer.data(fpi0 + 731);
    const auto *fpi0_734 = buffer.data(fpi0 + 734);
    const auto *fpi0_738 = buffer.data(fpi0 + 738);
    const auto *fpi0_740 = buffer.data(fpi0 + 740);
    const auto *fpi0_749 = buffer.data(fpi0 + 749);
    const auto *fpi0_750 = buffer.data(fpi0 + 750);
    const auto *fpi0_751 = buffer.data(fpi0 + 751);
    const auto *fpi0_752 = buffer.data(fpi0 + 752);
    const auto *fpi0_753 = buffer.data(fpi0 + 753);
    const auto *fpi0_755 = buffer.data(fpi0 + 755);
    const auto *fpi0_787 = buffer.data(fpi0 + 787);
    const auto *fpi0_790 = buffer.data(fpi0 + 790);
    const auto *fpi0_791 = buffer.data(fpi0 + 791);
    const auto *fpi0_794 = buffer.data(fpi0 + 794);
    const auto *fpi0_795 = buffer.data(fpi0 + 795);
    const auto *fpi0_796 = buffer.data(fpi0 + 796);
    const auto *fpi0_805 = buffer.data(fpi0 + 805);
    const auto *fpi0_806 = buffer.data(fpi0 + 806);
    const auto *fpi0_807 = buffer.data(fpi0 + 807);
    const auto *fpi0_808 = buffer.data(fpi0 + 808);
    const auto *fpi0_809 = buffer.data(fpi0 + 809);
    const auto *fpi0_811 = buffer.data(fpi0 + 811);
    const auto *fpi0_812 = buffer.data(fpi0 + 812);
    const auto *fpi0_817 = buffer.data(fpi0 + 817);

    const auto *fph_273 = buffer.data(fph + 273);
    const auto *fph_276 = buffer.data(fph + 276);
    const auto *fph_279 = buffer.data(fph + 279);
    const auto *fph_288 = buffer.data(fph + 288);
    const auto *fph_297 = buffer.data(fph + 297);
    const auto *fph_300 = buffer.data(fph + 300);
    const auto *fph_315 = buffer.data(fph + 315);
    const auto *fph_335 = buffer.data(fph + 335);
    const auto *fph_336 = buffer.data(fph + 336);
    const auto *fph_338 = buffer.data(fph + 338);
    const auto *fph_341 = buffer.data(fph + 341);
    const auto *fph_345 = buffer.data(fph + 345);
    const auto *fph_357 = buffer.data(fph + 357);
    const auto *fph_359 = buffer.data(fph + 359);
    const auto *fph_362 = buffer.data(fph + 362);
    const auto *fph_366 = buffer.data(fph + 366);
    const auto *fph_377 = buffer.data(fph + 377);
    const auto *fph_528 = buffer.data(fph + 528);
    const auto *fph_530 = buffer.data(fph + 530);
    const auto *fph_531 = buffer.data(fph + 531);
    const auto *fph_534 = buffer.data(fph + 534);
    const auto *fph_535 = buffer.data(fph + 535);
    const auto *fph_537 = buffer.data(fph + 537);
    const auto *fph_539 = buffer.data(fph + 539);
    const auto *fph_540 = buffer.data(fph + 540);
    const auto *fph_541 = buffer.data(fph + 541);
    const auto *fph_542 = buffer.data(fph + 542);
    const auto *fph_543 = buffer.data(fph + 543);
    const auto *fph_544 = buffer.data(fph + 544);
    const auto *fph_545 = buffer.data(fph + 545);
    const auto *fph_549 = buffer.data(fph + 549);
    const auto *fph_552 = buffer.data(fph + 552);
    const auto *fph_556 = buffer.data(fph + 556);
    const auto *fph_558 = buffer.data(fph + 558);
    const auto *fph_561 = buffer.data(fph + 561);
    const auto *fph_562 = buffer.data(fph + 562);
    const auto *fph_563 = buffer.data(fph + 563);
    const auto *fph_564 = buffer.data(fph + 564);
    const auto *fph_565 = buffer.data(fph + 565);
    const auto *fph_566 = buffer.data(fph + 566);
    const auto *fph_567 = buffer.data(fph + 567);
    const auto *fph_572 = buffer.data(fph + 572);
    const auto *fph_576 = buffer.data(fph + 576);
    const auto *fph_581 = buffer.data(fph + 581);
    const auto *fph_582 = buffer.data(fph + 582);
    const auto *fph_583 = buffer.data(fph + 583);
    const auto *fph_584 = buffer.data(fph + 584);
    const auto *fph_585 = buffer.data(fph + 585);
    const auto *fph_587 = buffer.data(fph + 587);
    const auto *fph_591 = buffer.data(fph + 591);
    const auto *fph_594 = buffer.data(fph + 594);
    const auto *fph_595 = buffer.data(fph + 595);
    const auto *fph_598 = buffer.data(fph + 598);
    const auto *fph_599 = buffer.data(fph + 599);
    const auto *fph_600 = buffer.data(fph + 600);
    const auto *fph_603 = buffer.data(fph + 603);
    const auto *fph_604 = buffer.data(fph + 604);
    const auto *fph_605 = buffer.data(fph + 605);
    const auto *fph_606 = buffer.data(fph + 606);
    const auto *fph_608 = buffer.data(fph + 608);
    const auto *fph_609 = buffer.data(fph + 609);
    const auto *fph_614 = buffer.data(fph + 614);

    const auto *fpi1_476 = buffer.data(fpi1 + 476);
    const auto *fpi1_478 = buffer.data(fpi1 + 478);
    const auto *fpi1_481 = buffer.data(fpi1 + 481);
    const auto *fpi1_485 = buffer.data(fpi1 + 485);
    const auto *fpi1_490 = buffer.data(fpi1 + 490);
    const auto *fpi1_703 = buffer.data(fpi1 + 703);
    const auto *fpi1_706 = buffer.data(fpi1 + 706);
    const auto *fpi1_710 = buffer.data(fpi1 + 710);
    const auto *fpi1_712 = buffer.data(fpi1 + 712);
    const auto *fpi1_721 = buffer.data(fpi1 + 721);
    const auto *fpi1_723 = buffer.data(fpi1 + 723);
    const auto *fpi1_724 = buffer.data(fpi1 + 724);
    const auto *fpi1_725 = buffer.data(fpi1 + 725);
    const auto *fpi1_726 = buffer.data(fpi1 + 726);
    const auto *fpi1_727 = buffer.data(fpi1 + 727);
    const auto *fpi1_731 = buffer.data(fpi1 + 731);
    const auto *fpi1_734 = buffer.data(fpi1 + 734);
    const auto *fpi1_738 = buffer.data(fpi1 + 738);
    const auto *fpi1_740 = buffer.data(fpi1 + 740);
    const auto *fpi1_749 = buffer.data(fpi1 + 749);
    const auto *fpi1_750 = buffer.data(fpi1 + 750);
    const auto *fpi1_751 = buffer.data(fpi1 + 751);
    const auto *fpi1_752 = buffer.data(fpi1 + 752);
    const auto *fpi1_753 = buffer.data(fpi1 + 753);
    const auto *fpi1_755 = buffer.data(fpi1 + 755);
    const auto *fpi1_787 = buffer.data(fpi1 + 787);
    const auto *fpi1_790 = buffer.data(fpi1 + 790);
    const auto *fpi1_791 = buffer.data(fpi1 + 791);
    const auto *fpi1_794 = buffer.data(fpi1 + 794);
    const auto *fpi1_795 = buffer.data(fpi1 + 795);
    const auto *fpi1_796 = buffer.data(fpi1 + 796);
    const auto *fpi1_805 = buffer.data(fpi1 + 805);
    const auto *fpi1_806 = buffer.data(fpi1 + 806);
    const auto *fpi1_807 = buffer.data(fpi1 + 807);
    const auto *fpi1_808 = buffer.data(fpi1 + 808);
    const auto *fpi1_809 = buffer.data(fpi1 + 809);
    const auto *fpi1_811 = buffer.data(fpi1 + 811);
    const auto *fpi1_812 = buffer.data(fpi1 + 812);
    const auto *fpi1_817 = buffer.data(fpi1 + 817);

    const auto *gsi0_252 = buffer.data(gsi0 + 252);
    const auto *gsi0_257 = buffer.data(gsi0 + 257);
    const auto *gsi0_261 = buffer.data(gsi0 + 261);
    const auto *gsi0_266 = buffer.data(gsi0 + 266);

    const auto *gsh_168 = buffer.data(gsh + 168);
    const auto *gsh_170 = buffer.data(gsh + 170);
    const auto *gsh_171 = buffer.data(gsh + 171);
    const auto *gsh_173 = buffer.data(gsh + 173);
    const auto *gsh_174 = buffer.data(gsh + 174);
    const auto *gsh_177 = buffer.data(gsh + 177);
    const auto *gsh_189 = buffer.data(gsh + 189);
    const auto *gsh_191 = buffer.data(gsh + 191);
    const auto *gsh_194 = buffer.data(gsh + 194);
    const auto *gsh_198 = buffer.data(gsh + 198);
    const auto *gsh_203 = buffer.data(gsh + 203);
    const auto *gsh_204 = buffer.data(gsh + 204);
    const auto *gsh_205 = buffer.data(gsh + 205);
    const auto *gsh_206 = buffer.data(gsh + 206);
    const auto *gsh_207 = buffer.data(gsh + 207);
    const auto *gsh_209 = buffer.data(gsh + 209);

    const auto *gsi1_252 = buffer.data(gsi1 + 252);
    const auto *gsi1_257 = buffer.data(gsi1 + 257);
    const auto *gsi1_261 = buffer.data(gsi1 + 261);
    const auto *gsi1_266 = buffer.data(gsi1 + 266);

    const auto *gpg0_380 = buffer.data(gpg0 + 380);
    const auto *gpg0_384 = buffer.data(gpg0 + 384);
    const auto *gpg0_389 = buffer.data(gpg0 + 389);
    const auto *gpg0_405 = buffer.data(gpg0 + 405);
    const auto *gpg0_406 = buffer.data(gpg0 + 406);
    const auto *gpg0_407 = buffer.data(gpg0 + 407);
    const auto *gpg0_408 = buffer.data(gpg0 + 408);
    const auto *gpg0_409 = buffer.data(gpg0 + 409);
    const auto *gpg0_410 = buffer.data(gpg0 + 410);
    const auto *gpg0_414 = buffer.data(gpg0 + 414);
    const auto *gpg0_415 = buffer.data(gpg0 + 415);
    const auto *gpg0_416 = buffer.data(gpg0 + 416);
    const auto *gpg0_417 = buffer.data(gpg0 + 417);
    const auto *gpg0_418 = buffer.data(gpg0 + 418);
    const auto *gpg0_419 = buffer.data(gpg0 + 419);
    const auto *gpg0_435 = buffer.data(gpg0 + 435);

    const auto *gpg1_380 = buffer.data(gpg1 + 380);
    const auto *gpg1_384 = buffer.data(gpg1 + 384);
    const auto *gpg1_389 = buffer.data(gpg1 + 389);
    const auto *gpg1_405 = buffer.data(gpg1 + 405);
    const auto *gpg1_406 = buffer.data(gpg1 + 406);
    const auto *gpg1_407 = buffer.data(gpg1 + 407);
    const auto *gpg1_408 = buffer.data(gpg1 + 408);
    const auto *gpg1_409 = buffer.data(gpg1 + 409);
    const auto *gpg1_410 = buffer.data(gpg1 + 410);
    const auto *gpg1_414 = buffer.data(gpg1 + 414);
    const auto *gpg1_415 = buffer.data(gpg1 + 415);
    const auto *gpg1_416 = buffer.data(gpg1 + 416);
    const auto *gpg1_417 = buffer.data(gpg1 + 417);
    const auto *gpg1_418 = buffer.data(gpg1 + 418);
    const auto *gpg1_419 = buffer.data(gpg1 + 419);
    const auto *gpg1_435 = buffer.data(gpg1 + 435);

    const auto *gph_525 = buffer.data(gph + 525);
    const auto *gph_527 = buffer.data(gph + 527);
    const auto *gph_528 = buffer.data(gph + 528);
    const auto *gph_530 = buffer.data(gph + 530);
    const auto *gph_531 = buffer.data(gph + 531);
    const auto *gph_534 = buffer.data(gph + 534);
    const auto *gph_539 = buffer.data(gph + 539);
    const auto *gph_540 = buffer.data(gph + 540);
    const auto *gph_541 = buffer.data(gph + 541);
    const auto *gph_542 = buffer.data(gph + 542);
    const auto *gph_543 = buffer.data(gph + 543);
    const auto *gph_544 = buffer.data(gph + 544);
    const auto *gph_545 = buffer.data(gph + 545);
    const auto *gph_546 = buffer.data(gph + 546);
    const auto *gph_548 = buffer.data(gph + 548);
    const auto *gph_549 = buffer.data(gph + 549);
    const auto *gph_551 = buffer.data(gph + 551);
    const auto *gph_552 = buffer.data(gph + 552);
    const auto *gph_555 = buffer.data(gph + 555);
    const auto *gph_561 = buffer.data(gph + 561);
    const auto *gph_562 = buffer.data(gph + 562);
    const auto *gph_563 = buffer.data(gph + 563);
    const auto *gph_564 = buffer.data(gph + 564);
    const auto *gph_565 = buffer.data(gph + 565);
    const auto *gph_566 = buffer.data(gph + 566);
    const auto *gph_567 = buffer.data(gph + 567);
    const auto *gph_568 = buffer.data(gph + 568);
    const auto *gph_569 = buffer.data(gph + 569);
    const auto *gph_570 = buffer.data(gph + 570);
    const auto *gph_571 = buffer.data(gph + 571);
    const auto *gph_572 = buffer.data(gph + 572);
    const auto *gph_573 = buffer.data(gph + 573);
    const auto *gph_574 = buffer.data(gph + 574);
    const auto *gph_575 = buffer.data(gph + 575);
    const auto *gph_576 = buffer.data(gph + 576);
    const auto *gph_581 = buffer.data(gph + 581);
    const auto *gph_582 = buffer.data(gph + 582);
    const auto *gph_583 = buffer.data(gph + 583);
    const auto *gph_584 = buffer.data(gph + 584);
    const auto *gph_585 = buffer.data(gph + 585);
    const auto *gph_586 = buffer.data(gph + 586);
    const auto *gph_587 = buffer.data(gph + 587);
    const auto *gph_588 = buffer.data(gph + 588);
    const auto *gph_590 = buffer.data(gph + 590);
    const auto *gph_593 = buffer.data(gph + 593);
    const auto *gph_597 = buffer.data(gph + 597);
    const auto *gph_602 = buffer.data(gph + 602);
    const auto *gph_603 = buffer.data(gph + 603);
    const auto *gph_604 = buffer.data(gph + 604);
    const auto *gph_605 = buffer.data(gph + 605);
    const auto *gph_606 = buffer.data(gph + 606);
    const auto *gph_608 = buffer.data(gph + 608);
    const auto *gph_609 = buffer.data(gph + 609);
    const auto *gph_610 = buffer.data(gph + 610);
    const auto *gph_611 = buffer.data(gph + 611);

#pragma omp simd aligned(t_701, t_702, t_703, pa_x, pc_x, pc_y, pc_z, fpi0_703, fph_273, \
                         fph_336, fph_528, fpi1_703, gsh_168, gph_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_1 * fph_336[k]
                   + f_1 * gsh_168[k]
                   + f_4 * pc_y[k] * gph_525[k];

        t_702[k] = f_12 * fph_273[k]
                   + f_4 * pc_z[k] * gph_525[k];

        t_703[k] = pa_x[k] * fpi0_703[k]
                   + f_0 * fph_528[k]
                   - f_11 * pc_x[k] * fpi1_703[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, pa_x, pc_x, pc_y, fpi0_706, fph_338, fph_530, \
                         fph_531, fpi1_706, gsh_170, gpg0_380, gpg1_380, gph_527, \
                         gph_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = f_1 * fph_338[k]
                   + f_1 * gsh_170[k]
                   + f_4 * pc_y[k] * gph_527[k];

        t_705[k] = f_1 * fph_530[k]
                   + f_9 * gpg0_380[k]
                   - f_10 * gpg1_380[k]
                   + f_4 * pc_x[k] * gph_530[k];

        t_706[k] = pa_x[k] * fpi0_706[k]
                   + f_13 * fph_531[k]
                   - f_11 * pc_x[k] * fpi1_706[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, pc_x, pc_y, pc_z, fph_276, fph_341, fph_534, \
                         gsh_173, gpg0_384, gpg1_384, gph_528, gph_530, \
                         gph_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = f_12 * fph_276[k]
                   + f_4 * pc_z[k] * gph_528[k];

        t_708[k] = f_1 * fph_341[k]
                   + f_1 * gsh_173[k]
                   + f_4 * pc_y[k] * gph_530[k];

        t_709[k] = f_1 * fph_534[k]
                   + f_7 * gpg0_384[k]
                   - f_8 * gpg1_384[k]
                   + f_4 * pc_x[k] * gph_534[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, pa_x, pc_x, pc_z, fpi0_710, fpi0_712, fph_279, \
                         fph_535, fph_537, fpi1_710, fpi1_712, \
                         gph_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = pa_x[k] * fpi0_710[k]
                   + f_12 * fph_535[k]
                   - f_11 * pc_x[k] * fpi1_710[k];

        t_711[k] = f_12 * fph_279[k]
                   + f_4 * pc_z[k] * gph_531[k];

        t_712[k] = pa_x[k] * fpi0_712[k]
                   + f_12 * fph_537[k]
                   - f_11 * pc_x[k] * fpi1_712[k];
    }

#pragma omp simd aligned(t_713, t_714, t_715, pc_x, pc_y, fph_345, fph_539, fph_540, gsh_177, \
                         gpg0_389, gpg1_389, gph_534, gph_539, \
                         gph_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_713[k] = f_1 * fph_345[k]
                   + f_1 * gsh_177[k]
                   + f_4 * pc_y[k] * gph_534[k];

        t_714[k] = f_1 * fph_539[k]
                   + f_5 * gpg0_389[k]
                   - f_6 * gpg1_389[k]
                   + f_4 * pc_x[k] * gph_539[k];

        t_715[k] = f_1 * fph_540[k]
                   + f_4 * pc_x[k] * gph_540[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, t_720, pc_x, fph_541, fph_542, fph_543, \
                         fph_544, fph_545, gph_541, gph_542, gph_543, gph_544, \
                         gph_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_1 * fph_541[k]
                   + f_4 * pc_x[k] * gph_541[k];

        t_717[k] = f_1 * fph_542[k]
                   + f_4 * pc_x[k] * gph_542[k];

        t_718[k] = f_1 * fph_543[k]
                   + f_4 * pc_x[k] * gph_543[k];

        t_719[k] = f_1 * fph_544[k]
                   + f_4 * pc_x[k] * gph_544[k];

        t_720[k] = f_1 * fph_545[k]
                   + f_4 * pc_x[k] * gph_545[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pa_x, pc_x, pc_z, fpi0_721, fpi0_723, \
                         fpi0_724, fph_288, fpi1_721, fpi1_723, fpi1_724, \
                         gph_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = pa_x[k] * fpi0_721[k]
                   - f_11 * pc_x[k] * fpi1_721[k];

        t_722[k] = f_12 * fph_288[k]
                   + f_4 * pc_z[k] * gph_540[k];

        t_723[k] = pa_x[k] * fpi0_723[k]
                   - f_11 * pc_x[k] * fpi1_723[k];

        t_724[k] = pa_x[k] * fpi0_724[k]
                   - f_11 * pc_x[k] * fpi1_724[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pa_x, pa_y, pc_x, pc_y, fpi0_476, \
                         fpi0_725, fpi0_726, fpi0_727, fpi1_476, fpi1_725, fpi1_726, \
                         fpi1_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = pa_x[k] * fpi0_725[k]
                   - f_11 * pc_x[k] * fpi1_725[k];

        t_726[k] = pa_x[k] * fpi0_726[k]
                   - f_11 * pc_x[k] * fpi1_726[k];

        t_727[k] = pa_x[k] * fpi0_727[k]
                   - f_11 * pc_x[k] * fpi1_727[k];

        t_728[k] = pa_y[k] * fpi0_476[k]
                   - f_11 * pc_y[k] * fpi1_476[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pa_x, pa_y, pc_x, pc_y, fpi0_478, fpi0_731, \
                         fph_357, fph_549, fpi1_478, fpi1_731, \
                         gph_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_1 * fph_357[k]
                   + f_4 * pc_y[k] * gph_546[k];

        t_730[k] = pa_y[k] * fpi0_478[k]
                   - f_11 * pc_y[k] * fpi1_478[k];

        t_731[k] = pa_x[k] * fpi0_731[k]
                   + f_0 * fph_549[k]
                   - f_11 * pc_x[k] * fpi1_731[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pa_x, pa_y, pc_x, pc_y, fpi0_481, fpi0_734, \
                         fph_359, fph_552, fpi1_481, fpi1_734, \
                         gph_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_1 * fph_359[k]
                   + f_4 * pc_y[k] * gph_548[k];

        t_733[k] = pa_y[k] * fpi0_481[k]
                   - f_11 * pc_y[k] * fpi1_481[k];

        t_734[k] = pa_x[k] * fpi0_734[k]
                   + f_13 * fph_552[k]
                   - f_11 * pc_x[k] * fpi1_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pa_y, pc_y, pc_z, fpi0_485, fph_297, fph_362, \
                         fpi1_485, gsh_171, gph_549, gph_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_12 * fph_297[k]
                   + f_1 * gsh_171[k]
                   + f_4 * pc_z[k] * gph_549[k];

        t_736[k] = f_1 * fph_362[k]
                   + f_4 * pc_y[k] * gph_551[k];

        t_737[k] = pa_y[k] * fpi0_485[k]
                   - f_11 * pc_y[k] * fpi1_485[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pa_x, pc_x, pc_z, fpi0_738, fpi0_740, fph_300, \
                         fph_556, fph_558, fpi1_738, fpi1_740, gsh_174, \
                         gph_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pa_x[k] * fpi0_738[k]
                   + f_12 * fph_556[k]
                   - f_11 * pc_x[k] * fpi1_738[k];

        t_739[k] = f_12 * fph_300[k]
                   + f_1 * gsh_174[k]
                   + f_4 * pc_z[k] * gph_552[k];

        t_740[k] = pa_x[k] * fpi0_740[k]
                   + f_12 * fph_558[k]
                   - f_11 * pc_x[k] * fpi1_740[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pa_y, pc_x, pc_y, fpi0_490, fph_366, \
                         fph_561, fph_562, fpi1_490, gph_555, gph_561, \
                         gph_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_1 * fph_366[k]
                   + f_4 * pc_y[k] * gph_555[k];

        t_742[k] = pa_y[k] * fpi0_490[k]
                   - f_11 * pc_y[k] * fpi1_490[k];

        t_743[k] = f_1 * fph_561[k]
                   + f_4 * pc_x[k] * gph_561[k];

        t_744[k] = f_1 * fph_562[k]
                   + f_4 * pc_x[k] * gph_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, fph_563, fph_564, fph_565, fph_566, \
                         gph_563, gph_564, gph_565, gph_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_1 * fph_563[k]
                   + f_4 * pc_x[k] * gph_563[k];

        t_746[k] = f_1 * fph_564[k]
                   + f_4 * pc_x[k] * gph_564[k];

        t_747[k] = f_1 * fph_565[k]
                   + f_4 * pc_x[k] * gph_565[k];

        t_748[k] = f_1 * fph_566[k]
                   + f_4 * pc_x[k] * gph_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, t_752, pa_x, pc_x, fpi0_749, fpi0_750, fpi0_751, \
                         fpi0_752, fpi1_749, fpi1_750, fpi1_751, \
                         fpi1_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = pa_x[k] * fpi0_749[k]
                   - f_11 * pc_x[k] * fpi1_749[k];

        t_750[k] = pa_x[k] * fpi0_750[k]
                   - f_11 * pc_x[k] * fpi1_750[k];

        t_751[k] = pa_x[k] * fpi0_751[k]
                   - f_11 * pc_x[k] * fpi1_751[k];

        t_752[k] = pa_x[k] * fpi0_752[k]
                   - f_11 * pc_x[k] * fpi1_752[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, pa_x, pc_x, pc_y, fpi0_753, fpi0_755, fph_377, \
                         fpi1_753, fpi1_755, gph_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = pa_x[k] * fpi0_753[k]
                   - f_11 * pc_x[k] * fpi1_753[k];

        t_754[k] = f_1 * fph_377[k]
                   + f_4 * pc_y[k] * gph_566[k];

        t_755[k] = pa_x[k] * fpi0_755[k]
                   - f_11 * pc_x[k] * fpi1_755[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, t_760, pc_x, pc_y, pc_z, fph_315, \
                         fph_567, gsh_189, gpg0_405, gpg1_405, gph_567, gph_568, \
                         gph_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_1 * fph_567[k]
                   + f_1 * gsh_189[k]
                   + f_2 * gpg0_405[k]
                   - f_3 * gpg1_405[k]
                   + f_4 * pc_x[k] * gph_567[k];

        t_757[k] = f_4 * pc_y[k] * gph_567[k];

        t_758[k] = f_13 * fph_315[k]
                   + f_4 * pc_z[k] * gph_567[k];

        t_759[k] = f_5 * gpg0_405[k]
                   - f_6 * gpg1_405[k]
                   + f_4 * pc_y[k] * gph_568[k];

        t_760[k] = f_4 * pc_y[k] * gph_569[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_x, pc_y, fph_572, gsh_194, gpg0_406, \
                         gpg0_407, gpg0_410, gpg1_406, gpg1_407, gpg1_410, gph_570, gph_571, \
                         gph_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_1 * fph_572[k]
                   + f_1 * gsh_194[k]
                   + f_9 * gpg0_410[k]
                   - f_10 * gpg1_410[k]
                   + f_4 * pc_x[k] * gph_572[k];

        t_762[k] = f_7 * gpg0_406[k]
                   - f_8 * gpg1_406[k]
                   + f_4 * pc_y[k] * gph_570[k];

        t_763[k] = f_5 * gpg0_407[k]
                   - f_6 * gpg1_407[k]
                   + f_4 * pc_y[k] * gph_571[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pc_x, pc_y, fph_576, gsh_198, gpg0_408, \
                         gpg0_414, gpg1_408, gpg1_414, gph_572, gph_573, \
                         gph_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_4 * pc_y[k] * gph_572[k];

        t_765[k] = f_1 * fph_576[k]
                   + f_1 * gsh_198[k]
                   + f_7 * gpg0_414[k]
                   - f_8 * gpg1_414[k]
                   + f_4 * pc_x[k] * gph_576[k];

        t_766[k] = f_9 * gpg0_408[k]
                   - f_10 * gpg1_408[k]
                   + f_4 * pc_y[k] * gph_573[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pc_y, gpg0_409, gpg0_410, gpg1_409, gpg1_410, \
                         gph_574, gph_575, gph_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_7 * gpg0_409[k]
                   - f_8 * gpg1_409[k]
                   + f_4 * pc_y[k] * gph_574[k];

        t_768[k] = f_5 * gpg0_410[k]
                   - f_6 * gpg1_410[k]
                   + f_4 * pc_y[k] * gph_575[k];

        t_769[k] = f_4 * pc_y[k] * gph_576[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, pc_x, fph_581, fph_582, fph_583, gsh_203, \
                         gsh_204, gsh_205, gpg0_419, gpg1_419, gph_581, gph_582, \
                         gph_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_1 * fph_581[k]
                   + f_1 * gsh_203[k]
                   + f_5 * gpg0_419[k]
                   - f_6 * gpg1_419[k]
                   + f_4 * pc_x[k] * gph_581[k];

        t_771[k] = f_1 * fph_582[k]
                   + f_1 * gsh_204[k]
                   + f_4 * pc_x[k] * gph_582[k];

        t_772[k] = f_1 * fph_583[k]
                   + f_1 * gsh_205[k]
                   + f_4 * pc_x[k] * gph_583[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, t_776, pc_x, pc_y, fph_584, fph_585, fph_587, \
                         gsh_206, gsh_207, gsh_209, gph_581, gph_584, gph_585, \
                         gph_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_1 * fph_584[k]
                   + f_1 * gsh_206[k]
                   + f_4 * pc_x[k] * gph_584[k];

        t_774[k] = f_1 * fph_585[k]
                   + f_1 * gsh_207[k]
                   + f_4 * pc_x[k] * gph_585[k];

        t_775[k] = f_4 * pc_y[k] * gph_581[k];

        t_776[k] = f_1 * fph_587[k]
                   + f_1 * gsh_209[k]
                   + f_4 * pc_x[k] * gph_587[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, gpg0_415, gpg0_416, gpg0_417, gpg1_415, \
                         gpg1_416, gpg1_417, gph_582, gph_583, \
                         gph_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_2 * gpg0_415[k]
                   - f_3 * gpg1_415[k]
                   + f_4 * pc_y[k] * gph_582[k];

        t_778[k] = f_17 * gpg0_416[k]
                   - f_18 * gpg1_416[k]
                   + f_4 * pc_y[k] * gph_583[k];

        t_779[k] = f_9 * gpg0_417[k]
                   - f_10 * gpg1_417[k]
                   + f_4 * pc_y[k] * gph_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, fph_335, gpg0_418, gpg0_419, \
                         gpg1_418, gpg1_419, gph_585, gph_586, \
                         gph_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_7 * gpg0_418[k]
                   - f_8 * gpg1_418[k]
                   + f_4 * pc_y[k] * gph_585[k];

        t_781[k] = f_5 * gpg0_419[k]
                   - f_6 * gpg1_419[k]
                   + f_4 * pc_y[k] * gph_586[k];

        t_782[k] = f_4 * pc_y[k] * gph_587[k];

        t_783[k] = f_13 * fph_335[k]
                   + f_2 * gpg0_419[k]
                   - f_3 * gpg1_419[k]
                   + f_4 * pc_z[k] * gph_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, pb_y, pc_y, pc_z, fph_336, gsi0_252, gsh_189, \
                         gsi1_252, gph_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = pb_y[k] * gsi0_252[k]
                   - f_11 * pc_y[k] * gsi1_252[k];

        t_785[k] = f_1 * gsh_189[k]
                   + f_4 * pc_y[k] * gph_588[k];

        t_786[k] = f_13 * fph_336[k]
                   + f_4 * pc_z[k] * gph_588[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pa_x, pb_y, pc_x, pc_y, fpi0_787, fph_591, \
                         fpi1_787, gsi0_257, gsh_191, gsi1_257, \
                         gph_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = pa_x[k] * fpi0_787[k]
                   + f_0 * fph_591[k]
                   - f_11 * pc_x[k] * fpi1_787[k];

        t_788[k] = f_1 * gsh_191[k]
                   + f_4 * pc_y[k] * gph_590[k];

        t_789[k] = pb_y[k] * gsi0_257[k]
                   - f_11 * pc_y[k] * gsi1_257[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, pa_x, pc_x, pc_y, fpi0_790, fpi0_791, fph_594, \
                         fph_595, fpi1_790, fpi1_791, gsh_194, \
                         gph_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = pa_x[k] * fpi0_790[k]
                   + f_13 * fph_594[k]
                   - f_11 * pc_x[k] * fpi1_790[k];

        t_791[k] = pa_x[k] * fpi0_791[k]
                   + f_13 * fph_595[k]
                   - f_11 * pc_x[k] * fpi1_791[k];

        t_792[k] = f_1 * gsh_194[k]
                   + f_4 * pc_y[k] * gph_593[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, pa_x, pb_y, pc_x, pc_y, fpi0_794, fpi0_795, \
                         fph_598, fph_599, fpi1_794, fpi1_795, gsi0_261, \
                         gsi1_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = pb_y[k] * gsi0_261[k]
                   - f_11 * pc_y[k] * gsi1_261[k];

        t_794[k] = pa_x[k] * fpi0_794[k]
                   + f_12 * fph_598[k]
                   - f_11 * pc_x[k] * fpi1_794[k];

        t_795[k] = pa_x[k] * fpi0_795[k]
                   + f_12 * fph_599[k]
                   - f_11 * pc_x[k] * fpi1_795[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, pa_x, pb_y, pc_x, pc_y, fpi0_796, fph_600, \
                         fpi1_796, gsi0_266, gsh_198, gsi1_266, \
                         gph_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = pa_x[k] * fpi0_796[k]
                   + f_12 * fph_600[k]
                   - f_11 * pc_x[k] * fpi1_796[k];

        t_797[k] = f_1 * gsh_198[k]
                   + f_4 * pc_y[k] * gph_597[k];

        t_798[k] = pb_y[k] * gsi0_266[k]
                   - f_11 * pc_y[k] * gsi1_266[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, t_802, pc_x, fph_603, fph_604, fph_605, fph_606, \
                         gph_603, gph_604, gph_605, gph_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = f_1 * fph_603[k]
                   + f_4 * pc_x[k] * gph_603[k];

        t_800[k] = f_1 * fph_604[k]
                   + f_4 * pc_x[k] * gph_604[k];

        t_801[k] = f_1 * fph_605[k]
                   + f_4 * pc_x[k] * gph_605[k];

        t_802[k] = f_1 * fph_606[k]
                   + f_4 * pc_x[k] * gph_606[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, t_806, pa_x, pc_x, pc_y, fpi0_805, fpi0_806, \
                         fph_608, fpi1_805, fpi1_806, gsh_203, gph_602, \
                         gph_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_1 * gsh_203[k]
                   + f_4 * pc_y[k] * gph_602[k];

        t_804[k] = f_1 * fph_608[k]
                   + f_4 * pc_x[k] * gph_608[k];

        t_805[k] = pa_x[k] * fpi0_805[k]
                   - f_11 * pc_x[k] * fpi1_805[k];

        t_806[k] = pa_x[k] * fpi0_806[k]
                   - f_11 * pc_x[k] * fpi1_806[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, pa_x, pc_x, pc_y, fpi0_807, fpi0_808, \
                         fpi0_809, fpi1_807, fpi1_808, fpi1_809, gsh_209, \
                         gph_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = pa_x[k] * fpi0_807[k]
                   - f_11 * pc_x[k] * fpi1_807[k];

        t_808[k] = pa_x[k] * fpi0_808[k]
                   - f_11 * pc_x[k] * fpi1_808[k];

        t_809[k] = pa_x[k] * fpi0_809[k]
                   - f_11 * pc_x[k] * fpi1_809[k];

        t_810[k] = f_1 * gsh_209[k]
                   + f_4 * pc_y[k] * gph_608[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, pa_x, pc_x, pc_y, pc_z, fpi0_811, \
                         fpi0_812, fph_357, fph_609, fpi1_811, fpi1_812, gsh_189, \
                         gph_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = pa_x[k] * fpi0_811[k]
                   - f_11 * pc_x[k] * fpi1_811[k];

        t_812[k] = pa_x[k] * fpi0_812[k]
                   + f_14 * fph_609[k]
                   - f_11 * pc_x[k] * fpi1_812[k];

        t_813[k] = f_4 * pc_y[k] * gph_609[k];

        t_814[k] = f_13 * fph_357[k]
                   + f_1 * gsh_189[k]
                   + f_4 * pc_z[k] * gph_609[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pa_x, pc_x, pc_y, fpi0_817, fph_614, fpi1_817, \
                         gpg0_435, gpg1_435, gph_610, gph_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_5 * gpg0_435[k]
                   - f_6 * gpg1_435[k]
                   + f_4 * pc_y[k] * gph_610[k];

        t_816[k] = f_4 * pc_y[k] * gph_611[k];

        t_817[k] = pa_x[k] * fpi0_817[k]
                   + f_0 * fph_614[k]
                   - f_11 * pc_x[k] * fpi1_817[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t fpi0, const size_t fph,
                                                          const size_t fpi1, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t gpg0, const size_t gpg1,
                                                          const size_t gph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_14 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fpi0_504 = buffer.data(fpi0 + 504);
    const auto *fpi0_505 = buffer.data(fpi0 + 505);
    const auto *fpi0_507 = buffer.data(fpi0 + 507);
    const auto *fpi0_508 = buffer.data(fpi0 + 508);
    const auto *fpi0_510 = buffer.data(fpi0 + 510);
    const auto *fpi0_511 = buffer.data(fpi0 + 511);
    const auto *fpi0_512 = buffer.data(fpi0 + 512);
    const auto *fpi0_514 = buffer.data(fpi0 + 514);
    const auto *fpi0_515 = buffer.data(fpi0 + 515);
    const auto *fpi0_516 = buffer.data(fpi0 + 516);
    const auto *fpi0_517 = buffer.data(fpi0 + 517);
    const auto *fpi0_821 = buffer.data(fpi0 + 821);
    const auto *fpi0_826 = buffer.data(fpi0 + 826);
    const auto *fpi0_833 = buffer.data(fpi0 + 833);
    const auto *fpi0_834 = buffer.data(fpi0 + 834);
    const auto *fpi0_835 = buffer.data(fpi0 + 835);
    const auto *fpi0_836 = buffer.data(fpi0 + 836);
    const auto *fpi0_837 = buffer.data(fpi0 + 837);
    const auto *fpi0_839 = buffer.data(fpi0 + 839);

    const auto *fph_379 = buffer.data(fph + 379);
    const auto *fph_381 = buffer.data(fph + 381);
    const auto *fph_382 = buffer.data(fph + 382);
    const auto *fph_384 = buffer.data(fph + 384);
    const auto *fph_385 = buffer.data(fph + 385);
    const auto *fph_386 = buffer.data(fph + 386);
    const auto *fph_398 = buffer.data(fph + 398);
    const auto *fph_414 = buffer.data(fph + 414);
    const auto *fph_419 = buffer.data(fph + 419);
    const auto *fph_440 = buffer.data(fph + 440);
    const auto *fph_618 = buffer.data(fph + 618);
    const auto *fph_623 = buffer.data(fph + 623);
    const auto *fph_624 = buffer.data(fph + 624);
    const auto *fph_625 = buffer.data(fph + 625);
    const auto *fph_626 = buffer.data(fph + 626);
    const auto *fph_627 = buffer.data(fph + 627);
    const auto *fph_629 = buffer.data(fph + 629);

    const auto *fpi1_504 = buffer.data(fpi1 + 504);
    const auto *fpi1_505 = buffer.data(fpi1 + 505);
    const auto *fpi1_507 = buffer.data(fpi1 + 507);
    const auto *fpi1_508 = buffer.data(fpi1 + 508);
    const auto *fpi1_510 = buffer.data(fpi1 + 510);
    const auto *fpi1_511 = buffer.data(fpi1 + 511);
    const auto *fpi1_512 = buffer.data(fpi1 + 512);
    const auto *fpi1_514 = buffer.data(fpi1 + 514);
    const auto *fpi1_515 = buffer.data(fpi1 + 515);
    const auto *fpi1_516 = buffer.data(fpi1 + 516);
    const auto *fpi1_517 = buffer.data(fpi1 + 517);
    const auto *fpi1_821 = buffer.data(fpi1 + 821);
    const auto *fpi1_826 = buffer.data(fpi1 + 826);
    const auto *fpi1_833 = buffer.data(fpi1 + 833);
    const auto *fpi1_834 = buffer.data(fpi1 + 834);
    const auto *fpi1_835 = buffer.data(fpi1 + 835);
    const auto *fpi1_836 = buffer.data(fpi1 + 836);
    const auto *fpi1_837 = buffer.data(fpi1 + 837);
    const auto *fpi1_839 = buffer.data(fpi1 + 839);

    const auto *gsi0_280 = buffer.data(gsi0 + 280);
    const auto *gsi0_281 = buffer.data(gsi0 + 281);
    const auto *gsi0_283 = buffer.data(gsi0 + 283);
    const auto *gsi0_285 = buffer.data(gsi0 + 285);
    const auto *gsi0_286 = buffer.data(gsi0 + 286);
    const auto *gsi0_288 = buffer.data(gsi0 + 288);
    const auto *gsi0_289 = buffer.data(gsi0 + 289);
    const auto *gsi0_290 = buffer.data(gsi0 + 290);
    const auto *gsi0_292 = buffer.data(gsi0 + 292);
    const auto *gsi0_293 = buffer.data(gsi0 + 293);
    const auto *gsi0_294 = buffer.data(gsi0 + 294);
    const auto *gsi0_301 = buffer.data(gsi0 + 301);
    const auto *gsi0_303 = buffer.data(gsi0 + 303);
    const auto *gsi0_304 = buffer.data(gsi0 + 304);
    const auto *gsi0_305 = buffer.data(gsi0 + 305);
    const auto *gsi0_307 = buffer.data(gsi0 + 307);
    const auto *gsi0_310 = buffer.data(gsi0 + 310);
    const auto *gsi0_313 = buffer.data(gsi0 + 313);
    const auto *gsi0_317 = buffer.data(gsi0 + 317);
    const auto *gsi0_322 = buffer.data(gsi0 + 322);

    const auto *gsh_210 = buffer.data(gsh + 210);
    const auto *gsh_211 = buffer.data(gsh + 211);
    const auto *gsh_213 = buffer.data(gsh + 213);
    const auto *gsh_215 = buffer.data(gsh + 215);
    const auto *gsh_216 = buffer.data(gsh + 216);
    const auto *gsh_218 = buffer.data(gsh + 218);
    const auto *gsh_219 = buffer.data(gsh + 219);
    const auto *gsh_220 = buffer.data(gsh + 220);
    const auto *gsh_222 = buffer.data(gsh + 222);
    const auto *gsh_223 = buffer.data(gsh + 223);
    const auto *gsh_224 = buffer.data(gsh + 224);
    const auto *gsh_225 = buffer.data(gsh + 225);
    const auto *gsh_226 = buffer.data(gsh + 226);
    const auto *gsh_227 = buffer.data(gsh + 227);
    const auto *gsh_228 = buffer.data(gsh + 228);
    const auto *gsh_229 = buffer.data(gsh + 229);
    const auto *gsh_230 = buffer.data(gsh + 230);
    const auto *gsh_233 = buffer.data(gsh + 233);
    const auto *gsh_236 = buffer.data(gsh + 236);
    const auto *gsh_240 = buffer.data(gsh + 240);
    const auto *gsh_245 = buffer.data(gsh + 245);

    const auto *gsi1_280 = buffer.data(gsi1 + 280);
    const auto *gsi1_281 = buffer.data(gsi1 + 281);
    const auto *gsi1_283 = buffer.data(gsi1 + 283);
    const auto *gsi1_285 = buffer.data(gsi1 + 285);
    const auto *gsi1_286 = buffer.data(gsi1 + 286);
    const auto *gsi1_288 = buffer.data(gsi1 + 288);
    const auto *gsi1_289 = buffer.data(gsi1 + 289);
    const auto *gsi1_290 = buffer.data(gsi1 + 290);
    const auto *gsi1_292 = buffer.data(gsi1 + 292);
    const auto *gsi1_293 = buffer.data(gsi1 + 293);
    const auto *gsi1_294 = buffer.data(gsi1 + 294);
    const auto *gsi1_301 = buffer.data(gsi1 + 301);
    const auto *gsi1_303 = buffer.data(gsi1 + 303);
    const auto *gsi1_304 = buffer.data(gsi1 + 304);
    const auto *gsi1_305 = buffer.data(gsi1 + 305);
    const auto *gsi1_307 = buffer.data(gsi1 + 307);
    const auto *gsi1_310 = buffer.data(gsi1 + 310);
    const auto *gsi1_313 = buffer.data(gsi1 + 313);
    const auto *gsi1_317 = buffer.data(gsi1 + 317);
    const auto *gsi1_322 = buffer.data(gsi1 + 322);

    const auto *gpg0_436 = buffer.data(gpg0 + 436);
    const auto *gpg0_437 = buffer.data(gpg0 + 437);
    const auto *gpg0_438 = buffer.data(gpg0 + 438);
    const auto *gpg0_439 = buffer.data(gpg0 + 439);
    const auto *gpg0_440 = buffer.data(gpg0 + 440);
    const auto *gpg0_465 = buffer.data(gpg0 + 465);
    const auto *gpg0_466 = buffer.data(gpg0 + 466);
    const auto *gpg0_468 = buffer.data(gpg0 + 468);
    const auto *gpg0_470 = buffer.data(gpg0 + 470);
    const auto *gpg0_471 = buffer.data(gpg0 + 471);
    const auto *gpg0_473 = buffer.data(gpg0 + 473);
    const auto *gpg0_474 = buffer.data(gpg0 + 474);
    const auto *gpg0_475 = buffer.data(gpg0 + 475);
    const auto *gpg0_476 = buffer.data(gpg0 + 476);
    const auto *gpg0_477 = buffer.data(gpg0 + 477);
    const auto *gpg0_478 = buffer.data(gpg0 + 478);
    const auto *gpg0_479 = buffer.data(gpg0 + 479);
    const auto *gpg0_485 = buffer.data(gpg0 + 485);
    const auto *gpg0_488 = buffer.data(gpg0 + 488);
    const auto *gpg0_489 = buffer.data(gpg0 + 489);
    const auto *gpg0_492 = buffer.data(gpg0 + 492);
    const auto *gpg0_493 = buffer.data(gpg0 + 493);
    const auto *gpg0_494 = buffer.data(gpg0 + 494);

    const auto *gpg1_436 = buffer.data(gpg1 + 436);
    const auto *gpg1_437 = buffer.data(gpg1 + 437);
    const auto *gpg1_438 = buffer.data(gpg1 + 438);
    const auto *gpg1_439 = buffer.data(gpg1 + 439);
    const auto *gpg1_440 = buffer.data(gpg1 + 440);
    const auto *gpg1_465 = buffer.data(gpg1 + 465);
    const auto *gpg1_466 = buffer.data(gpg1 + 466);
    const auto *gpg1_468 = buffer.data(gpg1 + 468);
    const auto *gpg1_470 = buffer.data(gpg1 + 470);
    const auto *gpg1_471 = buffer.data(gpg1 + 471);
    const auto *gpg1_473 = buffer.data(gpg1 + 473);
    const auto *gpg1_474 = buffer.data(gpg1 + 474);
    const auto *gpg1_475 = buffer.data(gpg1 + 475);
    const auto *gpg1_476 = buffer.data(gpg1 + 476);
    const auto *gpg1_477 = buffer.data(gpg1 + 477);
    const auto *gpg1_478 = buffer.data(gpg1 + 478);
    const auto *gpg1_479 = buffer.data(gpg1 + 479);
    const auto *gpg1_485 = buffer.data(gpg1 + 485);
    const auto *gpg1_488 = buffer.data(gpg1 + 488);
    const auto *gpg1_489 = buffer.data(gpg1 + 489);
    const auto *gpg1_492 = buffer.data(gpg1 + 492);
    const auto *gpg1_493 = buffer.data(gpg1 + 493);
    const auto *gpg1_494 = buffer.data(gpg1 + 494);

    const auto *gph_612 = buffer.data(gph + 612);
    const auto *gph_613 = buffer.data(gph + 613);
    const auto *gph_614 = buffer.data(gph + 614);
    const auto *gph_615 = buffer.data(gph + 615);
    const auto *gph_616 = buffer.data(gph + 616);
    const auto *gph_617 = buffer.data(gph + 617);
    const auto *gph_618 = buffer.data(gph + 618);
    const auto *gph_623 = buffer.data(gph + 623);
    const auto *gph_624 = buffer.data(gph + 624);
    const auto *gph_625 = buffer.data(gph + 625);
    const auto *gph_626 = buffer.data(gph + 626);
    const auto *gph_627 = buffer.data(gph + 627);
    const auto *gph_629 = buffer.data(gph + 629);
    const auto *gph_630 = buffer.data(gph + 630);
    const auto *gph_631 = buffer.data(gph + 631);
    const auto *gph_633 = buffer.data(gph + 633);
    const auto *gph_636 = buffer.data(gph + 636);
    const auto *gph_645 = buffer.data(gph + 645);
    const auto *gph_646 = buffer.data(gph + 646);
    const auto *gph_647 = buffer.data(gph + 647);
    const auto *gph_648 = buffer.data(gph + 648);
    const auto *gph_649 = buffer.data(gph + 649);
    const auto *gph_650 = buffer.data(gph + 650);
    const auto *gph_651 = buffer.data(gph + 651);
    const auto *gph_652 = buffer.data(gph + 652);
    const auto *gph_654 = buffer.data(gph + 654);
    const auto *gph_656 = buffer.data(gph + 656);
    const auto *gph_657 = buffer.data(gph + 657);
    const auto *gph_659 = buffer.data(gph + 659);
    const auto *gph_660 = buffer.data(gph + 660);
    const auto *gph_661 = buffer.data(gph + 661);
    const auto *gph_663 = buffer.data(gph + 663);
    const auto *gph_664 = buffer.data(gph + 664);
    const auto *gph_665 = buffer.data(gph + 665);
    const auto *gph_666 = buffer.data(gph + 666);
    const auto *gph_667 = buffer.data(gph + 667);
    const auto *gph_668 = buffer.data(gph + 668);
    const auto *gph_669 = buffer.data(gph + 669);
    const auto *gph_670 = buffer.data(gph + 670);
    const auto *gph_671 = buffer.data(gph + 671);
    const auto *gph_672 = buffer.data(gph + 672);
    const auto *gph_673 = buffer.data(gph + 673);
    const auto *gph_675 = buffer.data(gph + 675);
    const auto *gph_677 = buffer.data(gph + 677);
    const auto *gph_678 = buffer.data(gph + 678);
    const auto *gph_680 = buffer.data(gph + 680);
    const auto *gph_681 = buffer.data(gph + 681);
    const auto *gph_684 = buffer.data(gph + 684);
    const auto *gph_685 = buffer.data(gph + 685);
    const auto *gph_686 = buffer.data(gph + 686);
    const auto *gph_687 = buffer.data(gph + 687);
    const auto *gph_688 = buffer.data(gph + 688);
    const auto *gph_689 = buffer.data(gph + 689);
    const auto *gph_690 = buffer.data(gph + 690);
    const auto *gph_691 = buffer.data(gph + 691);
    const auto *gph_692 = buffer.data(gph + 692);

#pragma omp simd aligned(t_818, t_819, t_820, pc_y, gpg0_436, gpg0_437, gpg1_436, gpg1_437, \
                         gph_612, gph_613, gph_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_7 * gpg0_436[k]
                   - f_8 * gpg1_436[k]
                   + f_4 * pc_y[k] * gph_612[k];

        t_819[k] = f_5 * gpg0_437[k]
                   - f_6 * gpg1_437[k]
                   + f_4 * pc_y[k] * gph_613[k];

        t_820[k] = f_4 * pc_y[k] * gph_614[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pa_x, pc_x, pc_y, fpi0_821, fph_618, fpi1_821, \
                         gpg0_438, gpg0_439, gpg1_438, gpg1_439, gph_615, \
                         gph_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = pa_x[k] * fpi0_821[k]
                   + f_13 * fph_618[k]
                   - f_11 * pc_x[k] * fpi1_821[k];

        t_822[k] = f_9 * gpg0_438[k]
                   - f_10 * gpg1_438[k]
                   + f_4 * pc_y[k] * gph_615[k];

        t_823[k] = f_7 * gpg0_439[k]
                   - f_8 * gpg1_439[k]
                   + f_4 * pc_y[k] * gph_616[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, t_827, pa_x, pc_x, pc_y, fpi0_826, fph_623, \
                         fph_624, fpi1_826, gpg0_440, gpg1_440, gph_617, gph_618, \
                         gph_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_5 * gpg0_440[k]
                   - f_6 * gpg1_440[k]
                   + f_4 * pc_y[k] * gph_617[k];

        t_825[k] = f_4 * pc_y[k] * gph_618[k];

        t_826[k] = pa_x[k] * fpi0_826[k]
                   + f_12 * fph_623[k]
                   - f_11 * pc_x[k] * fpi1_826[k];

        t_827[k] = f_1 * fph_624[k]
                   + f_4 * pc_x[k] * gph_624[k];
    }

#pragma omp simd aligned(t_828, t_829, t_830, t_831, t_832, pc_x, pc_y, fph_625, fph_626, \
                         fph_627, fph_629, gph_623, gph_625, gph_626, gph_627, \
                         gph_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_828[k] = f_1 * fph_625[k]
                   + f_4 * pc_x[k] * gph_625[k];

        t_829[k] = f_1 * fph_626[k]
                   + f_4 * pc_x[k] * gph_626[k];

        t_830[k] = f_1 * fph_627[k]
                   + f_4 * pc_x[k] * gph_627[k];

        t_831[k] = f_4 * pc_y[k] * gph_623[k];

        t_832[k] = f_1 * fph_629[k]
                   + f_4 * pc_x[k] * gph_629[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, t_836, pa_x, pc_x, fpi0_833, fpi0_834, fpi0_835, \
                         fpi0_836, fpi1_833, fpi1_834, fpi1_835, \
                         fpi1_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = pa_x[k] * fpi0_833[k]
                   - f_11 * pc_x[k] * fpi1_833[k];

        t_834[k] = pa_x[k] * fpi0_834[k]
                   - f_11 * pc_x[k] * fpi1_834[k];

        t_835[k] = pa_x[k] * fpi0_835[k]
                   - f_11 * pc_x[k] * fpi1_835[k];

        t_836[k] = pa_x[k] * fpi0_836[k]
                   - f_11 * pc_x[k] * fpi1_836[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, t_840, pa_x, pb_x, pc_x, pc_y, fpi0_837, \
                         fpi0_839, fpi1_837, fpi1_839, gsi0_280, gsh_210, gsi1_280, \
                         gph_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = pa_x[k] * fpi0_837[k]
                   - f_11 * pc_x[k] * fpi1_837[k];

        t_838[k] = f_4 * pc_y[k] * gph_629[k];

        t_839[k] = pa_x[k] * fpi0_839[k]
                   - f_11 * pc_x[k] * fpi1_839[k];

        t_840[k] = pb_x[k] * gsi0_280[k]
                   + f_14 * gsh_210[k]
                   - f_11 * pc_x[k] * gsi1_280[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pb_x, pc_x, pc_z, gsi0_281, gsi0_283, \
                         gsh_211, gsh_213, gsi1_281, gsi1_283, gph_630, \
                         gph_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = pb_x[k] * gsi0_281[k]
                   + f_19 * gsh_211[k]
                   - f_11 * pc_x[k] * gsi1_281[k];

        t_842[k] = f_4 * pc_z[k] * gph_630[k];

        t_843[k] = pb_x[k] * gsi0_283[k]
                   + f_0 * gsh_213[k]
                   - f_11 * pc_x[k] * gsi1_283[k];

        t_844[k] = f_4 * pc_z[k] * gph_631[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pb_x, pc_x, pc_z, gsi0_285, gsi0_286, gsh_215, \
                         gsh_216, gsi1_285, gsi1_286, gph_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = pb_x[k] * gsi0_285[k]
                   + f_0 * gsh_215[k]
                   - f_11 * pc_x[k] * gsi1_285[k];

        t_846[k] = pb_x[k] * gsi0_286[k]
                   + f_13 * gsh_216[k]
                   - f_11 * pc_x[k] * gsi1_286[k];

        t_847[k] = f_4 * pc_z[k] * gph_633[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pb_x, pc_x, gsi0_288, gsi0_289, gsi0_290, \
                         gsh_218, gsh_219, gsh_220, gsi1_288, gsi1_289, \
                         gsi1_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = pb_x[k] * gsi0_288[k]
                   + f_13 * gsh_218[k]
                   - f_11 * pc_x[k] * gsi1_288[k];

        t_849[k] = pb_x[k] * gsi0_289[k]
                   + f_13 * gsh_219[k]
                   - f_11 * pc_x[k] * gsi1_289[k];

        t_850[k] = pb_x[k] * gsi0_290[k]
                   + f_12 * gsh_220[k]
                   - f_11 * pc_x[k] * gsi1_290[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pb_x, pc_x, pc_z, gsi0_292, gsi0_293, gsh_222, \
                         gsh_223, gsi1_292, gsi1_293, gph_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_4 * pc_z[k] * gph_636[k];

        t_852[k] = pb_x[k] * gsi0_292[k]
                   + f_12 * gsh_222[k]
                   - f_11 * pc_x[k] * gsi1_292[k];

        t_853[k] = pb_x[k] * gsi0_293[k]
                   + f_12 * gsh_223[k]
                   - f_11 * pc_x[k] * gsi1_293[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pb_x, pc_x, gsi0_294, gsh_224, gsh_225, \
                         gsh_226, gsh_227, gsi1_294, gph_645, gph_646, \
                         gph_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = pb_x[k] * gsi0_294[k]
                   + f_12 * gsh_224[k]
                   - f_11 * pc_x[k] * gsi1_294[k];

        t_855[k] = f_1 * gsh_225[k]
                   + f_4 * pc_x[k] * gph_645[k];

        t_856[k] = f_1 * gsh_226[k]
                   + f_4 * pc_x[k] * gph_646[k];

        t_857[k] = f_1 * gsh_227[k]
                   + f_4 * pc_x[k] * gph_647[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pb_x, pc_x, gsi0_301, gsh_228, gsh_229, \
                         gsh_230, gsi1_301, gph_648, gph_649, gph_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_1 * gsh_228[k]
                   + f_4 * pc_x[k] * gph_648[k];

        t_859[k] = f_1 * gsh_229[k]
                   + f_4 * pc_x[k] * gph_649[k];

        t_860[k] = f_1 * gsh_230[k]
                   + f_4 * pc_x[k] * gph_650[k];

        t_861[k] = pb_x[k] * gsi0_301[k]
                   - f_11 * pc_x[k] * gsi1_301[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pb_x, pc_x, pc_z, gsi0_303, gsi0_304, \
                         gsi0_305, gsi1_303, gsi1_304, gsi1_305, \
                         gph_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_4 * pc_z[k] * gph_645[k];

        t_863[k] = pb_x[k] * gsi0_303[k]
                   - f_11 * pc_x[k] * gsi1_303[k];

        t_864[k] = pb_x[k] * gsi0_304[k]
                   - f_11 * pc_x[k] * gsi1_304[k];

        t_865[k] = pb_x[k] * gsi0_305[k]
                   - f_11 * pc_x[k] * gsi1_305[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pb_x, pc_x, pc_y, fph_398, gsi0_307, gsi1_307, \
                         gpg0_465, gpg1_465, gph_650, gph_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_0 * fph_398[k]
                   + f_4 * pc_y[k] * gph_650[k];

        t_867[k] = pb_x[k] * gsi0_307[k]
                   - f_11 * pc_x[k] * gsi1_307[k];

        t_868[k] = f_2 * gpg0_465[k]
                   - f_3 * gpg1_465[k]
                   + f_4 * pc_x[k] * gph_651[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_z, gpg0_466, gpg0_468, gpg1_466, \
                         gpg1_468, gph_651, gph_652, gph_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_17 * gpg0_466[k]
                   - f_18 * gpg1_466[k]
                   + f_4 * pc_x[k] * gph_652[k];

        t_870[k] = f_4 * pc_z[k] * gph_651[k];

        t_871[k] = f_9 * gpg0_468[k]
                   - f_10 * gpg1_468[k]
                   + f_4 * pc_x[k] * gph_654[k];

        t_872[k] = f_4 * pc_z[k] * gph_652[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pc_x, pc_z, gpg0_470, gpg0_471, gpg0_473, \
                         gpg1_470, gpg1_471, gpg1_473, gph_654, gph_656, gph_657, \
                         gph_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_9 * gpg0_470[k]
                   - f_10 * gpg1_470[k]
                   + f_4 * pc_x[k] * gph_656[k];

        t_874[k] = f_7 * gpg0_471[k]
                   - f_8 * gpg1_471[k]
                   + f_4 * pc_x[k] * gph_657[k];

        t_875[k] = f_4 * pc_z[k] * gph_654[k];

        t_876[k] = f_7 * gpg0_473[k]
                   - f_8 * gpg1_473[k]
                   + f_4 * pc_x[k] * gph_659[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, t_880, pc_x, pc_z, gpg0_474, gpg0_475, gpg0_477, \
                         gpg1_474, gpg1_475, gpg1_477, gph_657, gph_660, gph_661, \
                         gph_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_7 * gpg0_474[k]
                   - f_8 * gpg1_474[k]
                   + f_4 * pc_x[k] * gph_660[k];

        t_878[k] = f_5 * gpg0_475[k]
                   - f_6 * gpg1_475[k]
                   + f_4 * pc_x[k] * gph_661[k];

        t_879[k] = f_4 * pc_z[k] * gph_657[k];

        t_880[k] = f_5 * gpg0_477[k]
                   - f_6 * gpg1_477[k]
                   + f_4 * pc_x[k] * gph_663[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, t_884, t_885, pc_x, gpg0_478, gpg0_479, \
                         gpg1_478, gpg1_479, gph_664, gph_665, gph_666, gph_667, \
                         gph_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_5 * gpg0_478[k]
                   - f_6 * gpg1_478[k]
                   + f_4 * pc_x[k] * gph_664[k];

        t_882[k] = f_5 * gpg0_479[k]
                   - f_6 * gpg1_479[k]
                   + f_4 * pc_x[k] * gph_665[k];

        t_883[k] = f_4 * pc_x[k] * gph_666[k];

        t_884[k] = f_4 * pc_x[k] * gph_667[k];

        t_885[k] = f_4 * pc_x[k] * gph_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, t_890, pc_x, pc_y, pc_z, fph_414, \
                         gsh_225, gpg0_475, gpg1_475, gph_666, gph_669, gph_670, \
                         gph_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_4 * pc_x[k] * gph_669[k];

        t_887[k] = f_4 * pc_x[k] * gph_670[k];

        t_888[k] = f_4 * pc_x[k] * gph_671[k];

        t_889[k] = f_0 * fph_414[k]
                   + f_1 * gsh_225[k]
                   + f_2 * gpg0_475[k]
                   - f_3 * gpg1_475[k]
                   + f_4 * pc_y[k] * gph_666[k];

        t_890[k] = f_4 * pc_z[k] * gph_666[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, pc_z, gpg0_475, gpg0_476, gpg0_477, gpg1_475, \
                         gpg1_476, gpg1_477, gph_667, gph_668, \
                         gph_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = f_5 * gpg0_475[k]
                   - f_6 * gpg1_475[k]
                   + f_4 * pc_z[k] * gph_667[k];

        t_892[k] = f_7 * gpg0_476[k]
                   - f_8 * gpg1_476[k]
                   + f_4 * pc_z[k] * gph_668[k];

        t_893[k] = f_9 * gpg0_477[k]
                   - f_10 * gpg1_477[k]
                   + f_4 * pc_z[k] * gph_669[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pb_z, pc_y, pc_z, fph_419, gsi0_280, \
                         gsi0_281, gsh_230, gsi1_280, gsi1_281, gpg0_479, gpg1_479, \
                         gph_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_0 * fph_419[k]
                   + f_1 * gsh_230[k]
                   + f_4 * pc_y[k] * gph_671[k];

        t_895[k] = f_2 * gpg0_479[k]
                   - f_3 * gpg1_479[k]
                   + f_4 * pc_z[k] * gph_671[k];

        t_896[k] = pb_z[k] * gsi0_280[k]
                   - f_11 * pc_z[k] * gsi1_280[k];

        t_897[k] = pb_z[k] * gsi0_281[k]
                   - f_11 * pc_z[k] * gsi1_281[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, pb_z, pc_x, pc_z, gsi0_283, gsh_210, \
                         gsh_211, gsi1_283, gpg0_485, gpg1_485, gph_672, gph_673, \
                         gph_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_1 * gsh_210[k]
                   + f_4 * pc_z[k] * gph_672[k];

        t_899[k] = pb_z[k] * gsi0_283[k]
                   - f_11 * pc_z[k] * gsi1_283[k];

        t_900[k] = f_1 * gsh_211[k]
                   + f_4 * pc_z[k] * gph_673[k];

        t_901[k] = f_9 * gpg0_485[k]
                   - f_10 * gpg1_485[k]
                   + f_4 * pc_x[k] * gph_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pb_z, pc_x, pc_z, gsi0_286, gsh_213, gsi1_286, \
                         gpg0_488, gpg1_488, gph_675, gph_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = pb_z[k] * gsi0_286[k]
                   - f_11 * pc_z[k] * gsi1_286[k];

        t_903[k] = f_1 * gsh_213[k]
                   + f_4 * pc_z[k] * gph_675[k];

        t_904[k] = f_7 * gpg0_488[k]
                   - f_8 * gpg1_488[k]
                   + f_4 * pc_x[k] * gph_680[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pb_z, pc_x, pc_z, gsi0_290, gsh_216, gsi1_290, \
                         gpg0_489, gpg1_489, gph_678, gph_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_7 * gpg0_489[k]
                   - f_8 * gpg1_489[k]
                   + f_4 * pc_x[k] * gph_681[k];

        t_906[k] = pb_z[k] * gsi0_290[k]
                   - f_11 * pc_z[k] * gsi1_290[k];

        t_907[k] = f_1 * gsh_216[k]
                   + f_4 * pc_z[k] * gph_678[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, t_911, pc_x, gpg0_492, gpg0_493, gpg0_494, \
                         gpg1_492, gpg1_493, gpg1_494, gph_684, gph_685, gph_686, \
                         gph_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_5 * gpg0_492[k]
                   - f_6 * gpg1_492[k]
                   + f_4 * pc_x[k] * gph_684[k];

        t_909[k] = f_5 * gpg0_493[k]
                   - f_6 * gpg1_493[k]
                   + f_4 * pc_x[k] * gph_685[k];

        t_910[k] = f_5 * gpg0_494[k]
                   - f_6 * gpg1_494[k]
                   + f_4 * pc_x[k] * gph_686[k];

        t_911[k] = f_4 * pc_x[k] * gph_687[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, t_915, t_916, t_917, pb_z, pc_x, pc_z, gsi0_301, \
                         gsi1_301, gph_688, gph_689, gph_690, gph_691, \
                         gph_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_4 * pc_x[k] * gph_688[k];

        t_913[k] = f_4 * pc_x[k] * gph_689[k];

        t_914[k] = f_4 * pc_x[k] * gph_690[k];

        t_915[k] = f_4 * pc_x[k] * gph_691[k];

        t_916[k] = f_4 * pc_x[k] * gph_692[k];

        t_917[k] = pb_z[k] * gsi0_301[k]
                   - f_11 * pc_z[k] * gsi1_301[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pb_z, pc_z, gsi0_303, gsi0_304, gsh_225, \
                         gsh_226, gsh_227, gsi1_303, gsi1_304, \
                         gph_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_1 * gsh_225[k]
                   + f_4 * pc_z[k] * gph_687[k];

        t_919[k] = pb_z[k] * gsi0_303[k]
                   + f_12 * gsh_226[k]
                   - f_11 * pc_z[k] * gsi1_303[k];

        t_920[k] = pb_z[k] * gsi0_304[k]
                   + f_13 * gsh_227[k]
                   - f_11 * pc_z[k] * gsi1_304[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, pb_z, pc_y, pc_z, fph_440, gsi0_305, gsi0_307, \
                         gsh_228, gsh_230, gsi1_305, gsi1_307, \
                         gph_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = pb_z[k] * gsi0_305[k]
                   + f_0 * gsh_228[k]
                   - f_11 * pc_z[k] * gsi1_305[k];

        t_922[k] = f_0 * fph_440[k]
                   + f_4 * pc_y[k] * gph_692[k];

        t_923[k] = pb_z[k] * gsi0_307[k]
                   + f_14 * gsh_230[k]
                   - f_11 * pc_z[k] * gsi1_307[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pa_z, pb_x, pc_x, pc_z, fpi0_504, fpi0_505, \
                         fpi1_504, fpi1_505, gsi0_310, gsh_233, \
                         gsi1_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = pa_z[k] * fpi0_504[k]
                   - f_11 * pc_z[k] * fpi1_504[k];

        t_925[k] = pa_z[k] * fpi0_505[k]
                   - f_11 * pc_z[k] * fpi1_505[k];

        t_926[k] = pb_x[k] * gsi0_310[k]
                   + f_19 * gsh_233[k]
                   - f_11 * pc_x[k] * gsi1_310[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, pa_z, pb_x, pc_x, pc_z, fpi0_507, fpi0_508, \
                         fph_379, fpi1_507, fpi1_508, gsi0_313, gsh_236, \
                         gsi1_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = pa_z[k] * fpi0_507[k]
                   - f_11 * pc_z[k] * fpi1_507[k];

        t_928[k] = pa_z[k] * fpi0_508[k]
                   + f_1 * fph_379[k]
                   - f_11 * pc_z[k] * fpi1_508[k];

        t_929[k] = pb_x[k] * gsi0_313[k]
                   + f_0 * gsh_236[k]
                   - f_11 * pc_x[k] * gsi1_313[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, pa_z, pc_z, fpi0_510, fpi0_511, fpi0_512, \
                         fph_381, fph_382, fpi1_510, fpi1_511, \
                         fpi1_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = pa_z[k] * fpi0_510[k]
                   - f_11 * pc_z[k] * fpi1_510[k];

        t_931[k] = pa_z[k] * fpi0_511[k]
                   + f_1 * fph_381[k]
                   - f_11 * pc_z[k] * fpi1_511[k];

        t_932[k] = pa_z[k] * fpi0_512[k]
                   + f_12 * fph_382[k]
                   - f_11 * pc_z[k] * fpi1_512[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, pa_z, pb_x, pc_x, pc_z, fpi0_514, fpi0_515, \
                         fph_384, fpi1_514, fpi1_515, gsi0_317, gsh_240, \
                         gsi1_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = pb_x[k] * gsi0_317[k]
                   + f_13 * gsh_240[k]
                   - f_11 * pc_x[k] * gsi1_317[k];

        t_934[k] = pa_z[k] * fpi0_514[k]
                   - f_11 * pc_z[k] * fpi1_514[k];

        t_935[k] = pa_z[k] * fpi0_515[k]
                   + f_1 * fph_384[k]
                   - f_11 * pc_z[k] * fpi1_515[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, pa_z, pb_x, pc_x, pc_z, fpi0_516, fpi0_517, \
                         fph_385, fph_386, fpi1_516, fpi1_517, gsi0_322, gsh_245, \
                         gsi1_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = pa_z[k] * fpi0_516[k]
                   + f_12 * fph_385[k]
                   - f_11 * pc_z[k] * fpi1_516[k];

        t_937[k] = pa_z[k] * fpi0_517[k]
                   + f_13 * fph_386[k]
                   - f_11 * pc_z[k] * fpi1_517[k];

        t_938[k] = pb_x[k] * gsi0_322[k]
                   + f_12 * gsh_245[k]
                   - f_11 * pc_x[k] * gsi1_322[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dpi1,
                                                          const size_t fpi0, const size_t fph,
                                                          const size_t fpi1, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t gpg0, const size_t gpg1,
                                                          const size_t gph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_14 = 3.0 / q;
    const auto f_15 = 1.0 / p;
    const auto f_16 = gamma / (p * q);
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;
    const auto f_20 = 0.5 / p;
    const auto f_21 = 0.5 * gamma / (p * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpi0_301 = buffer.data(dpi0 + 301);
    const auto *dpi0_419 = buffer.data(dpi0 + 419);

    const auto *dpi1_301 = buffer.data(dpi1 + 301);
    const auto *dpi1_419 = buffer.data(dpi1 + 419);

    const auto *fpi0_525 = buffer.data(fpi0 + 525);
    const auto *fpi0_532 = buffer.data(fpi0 + 532);
    const auto *fpi0_533 = buffer.data(fpi0 + 533);
    const auto *fpi0_535 = buffer.data(fpi0 + 535);
    const auto *fpi0_538 = buffer.data(fpi0 + 538);
    const auto *fpi0_542 = buffer.data(fpi0 + 542);
    const auto *fpi0_553 = buffer.data(fpi0 + 553);
    const auto *fpi0_555 = buffer.data(fpi0 + 555);
    const auto *fpi0_556 = buffer.data(fpi0 + 556);
    const auto *fpi0_557 = buffer.data(fpi0 + 557);
    const auto *fpi0_637 = buffer.data(fpi0 + 637);
    const auto *fpi0_671 = buffer.data(fpi0 + 671);

    const auto *fph_393 = buffer.data(fph + 393);
    const auto *fph_414 = buffer.data(fph + 414);
    const auto *fph_415 = buffer.data(fph + 415);
    const auto *fph_416 = buffer.data(fph + 416);
    const auto *fph_417 = buffer.data(fph + 417);
    const auto *fph_419 = buffer.data(fph + 419);
    const auto *fph_435 = buffer.data(fph + 435);
    const auto *fph_456 = buffer.data(fph + 456);
    const auto *fph_461 = buffer.data(fph + 461);
    const auto *fph_477 = buffer.data(fph + 477);
    const auto *fph_482 = buffer.data(fph + 482);
    const auto *fph_498 = buffer.data(fph + 498);
    const auto *fph_500 = buffer.data(fph + 500);
    const auto *fph_501 = buffer.data(fph + 501);
    const auto *fph_502 = buffer.data(fph + 502);
    const auto *fph_503 = buffer.data(fph + 503);
    const auto *fph_524 = buffer.data(fph + 524);

    const auto *fpi1_525 = buffer.data(fpi1 + 525);
    const auto *fpi1_532 = buffer.data(fpi1 + 532);
    const auto *fpi1_533 = buffer.data(fpi1 + 533);
    const auto *fpi1_535 = buffer.data(fpi1 + 535);
    const auto *fpi1_538 = buffer.data(fpi1 + 538);
    const auto *fpi1_542 = buffer.data(fpi1 + 542);
    const auto *fpi1_553 = buffer.data(fpi1 + 553);
    const auto *fpi1_555 = buffer.data(fpi1 + 555);
    const auto *fpi1_556 = buffer.data(fpi1 + 556);
    const auto *fpi1_557 = buffer.data(fpi1 + 557);
    const auto *fpi1_637 = buffer.data(fpi1 + 637);
    const auto *fpi1_671 = buffer.data(fpi1 + 671);

    const auto *gsi0_331 = buffer.data(gsi0 + 331);
    const auto *gsi0_332 = buffer.data(gsi0 + 332);
    const auto *gsi0_333 = buffer.data(gsi0 + 333);
    const auto *gsi0_335 = buffer.data(gsi0 + 335);
    const auto *gsi0_336 = buffer.data(gsi0 + 336);
    const auto *gsi0_337 = buffer.data(gsi0 + 337);
    const auto *gsi0_338 = buffer.data(gsi0 + 338);
    const auto *gsi0_339 = buffer.data(gsi0 + 339);
    const auto *gsi0_340 = buffer.data(gsi0 + 340);
    const auto *gsi0_341 = buffer.data(gsi0 + 341);
    const auto *gsi0_342 = buffer.data(gsi0 + 342);
    const auto *gsi0_343 = buffer.data(gsi0 + 343);
    const auto *gsi0_344 = buffer.data(gsi0 + 344);
    const auto *gsi0_345 = buffer.data(gsi0 + 345);
    const auto *gsi0_346 = buffer.data(gsi0 + 346);
    const auto *gsi0_347 = buffer.data(gsi0 + 347);
    const auto *gsi0_348 = buffer.data(gsi0 + 348);
    const auto *gsi0_349 = buffer.data(gsi0 + 349);
    const auto *gsi0_350 = buffer.data(gsi0 + 350);
    const auto *gsi0_357 = buffer.data(gsi0 + 357);
    const auto *gsi0_359 = buffer.data(gsi0 + 359);
    const auto *gsi0_360 = buffer.data(gsi0 + 360);
    const auto *gsi0_361 = buffer.data(gsi0 + 361);
    const auto *gsi0_363 = buffer.data(gsi0 + 363);

    const auto *gsh_246 = buffer.data(gsh + 246);
    const auto *gsh_247 = buffer.data(gsh + 247);
    const auto *gsh_248 = buffer.data(gsh + 248);
    const auto *gsh_249 = buffer.data(gsh + 249);
    const auto *gsh_250 = buffer.data(gsh + 250);
    const auto *gsh_251 = buffer.data(gsh + 251);
    const auto *gsh_252 = buffer.data(gsh + 252);
    const auto *gsh_253 = buffer.data(gsh + 253);
    const auto *gsh_254 = buffer.data(gsh + 254);
    const auto *gsh_255 = buffer.data(gsh + 255);
    const auto *gsh_256 = buffer.data(gsh + 256);
    const auto *gsh_257 = buffer.data(gsh + 257);
    const auto *gsh_258 = buffer.data(gsh + 258);
    const auto *gsh_259 = buffer.data(gsh + 259);
    const auto *gsh_260 = buffer.data(gsh + 260);
    const auto *gsh_261 = buffer.data(gsh + 261);
    const auto *gsh_262 = buffer.data(gsh + 262);
    const auto *gsh_263 = buffer.data(gsh + 263);
    const auto *gsh_264 = buffer.data(gsh + 264);
    const auto *gsh_265 = buffer.data(gsh + 265);
    const auto *gsh_266 = buffer.data(gsh + 266);
    const auto *gsh_267 = buffer.data(gsh + 267);
    const auto *gsh_268 = buffer.data(gsh + 268);
    const auto *gsh_269 = buffer.data(gsh + 269);
    const auto *gsh_270 = buffer.data(gsh + 270);
    const auto *gsh_271 = buffer.data(gsh + 271);
    const auto *gsh_272 = buffer.data(gsh + 272);

    const auto *gsi1_331 = buffer.data(gsi1 + 331);
    const auto *gsi1_332 = buffer.data(gsi1 + 332);
    const auto *gsi1_333 = buffer.data(gsi1 + 333);
    const auto *gsi1_335 = buffer.data(gsi1 + 335);
    const auto *gsi1_336 = buffer.data(gsi1 + 336);
    const auto *gsi1_337 = buffer.data(gsi1 + 337);
    const auto *gsi1_338 = buffer.data(gsi1 + 338);
    const auto *gsi1_339 = buffer.data(gsi1 + 339);
    const auto *gsi1_340 = buffer.data(gsi1 + 340);
    const auto *gsi1_341 = buffer.data(gsi1 + 341);
    const auto *gsi1_342 = buffer.data(gsi1 + 342);
    const auto *gsi1_343 = buffer.data(gsi1 + 343);
    const auto *gsi1_344 = buffer.data(gsi1 + 344);
    const auto *gsi1_345 = buffer.data(gsi1 + 345);
    const auto *gsi1_346 = buffer.data(gsi1 + 346);
    const auto *gsi1_347 = buffer.data(gsi1 + 347);
    const auto *gsi1_348 = buffer.data(gsi1 + 348);
    const auto *gsi1_349 = buffer.data(gsi1 + 349);
    const auto *gsi1_350 = buffer.data(gsi1 + 350);
    const auto *gsi1_357 = buffer.data(gsi1 + 357);
    const auto *gsi1_359 = buffer.data(gsi1 + 359);
    const auto *gsi1_360 = buffer.data(gsi1 + 360);
    const auto *gsi1_361 = buffer.data(gsi1 + 361);
    const auto *gsi1_363 = buffer.data(gsi1 + 363);

    const auto *gpg0_512 = buffer.data(gpg0 + 512);
    const auto *gpg0_514 = buffer.data(gpg0 + 514);
    const auto *gpg0_515 = buffer.data(gpg0 + 515);
    const auto *gpg0_517 = buffer.data(gpg0 + 517);
    const auto *gpg0_518 = buffer.data(gpg0 + 518);
    const auto *gpg0_519 = buffer.data(gpg0 + 519);
    const auto *gpg0_521 = buffer.data(gpg0 + 521);
    const auto *gpg0_522 = buffer.data(gpg0 + 522);
    const auto *gpg0_523 = buffer.data(gpg0 + 523);
    const auto *gpg0_524 = buffer.data(gpg0 + 524);
    const auto *gpg0_525 = buffer.data(gpg0 + 525);
    const auto *gpg0_526 = buffer.data(gpg0 + 526);
    const auto *gpg0_527 = buffer.data(gpg0 + 527);
    const auto *gpg0_528 = buffer.data(gpg0 + 528);
    const auto *gpg0_529 = buffer.data(gpg0 + 529);
    const auto *gpg0_530 = buffer.data(gpg0 + 530);
    const auto *gpg0_531 = buffer.data(gpg0 + 531);
    const auto *gpg0_532 = buffer.data(gpg0 + 532);
    const auto *gpg0_533 = buffer.data(gpg0 + 533);
    const auto *gpg0_534 = buffer.data(gpg0 + 534);
    const auto *gpg0_535 = buffer.data(gpg0 + 535);
    const auto *gpg0_536 = buffer.data(gpg0 + 536);
    const auto *gpg0_537 = buffer.data(gpg0 + 537);
    const auto *gpg0_538 = buffer.data(gpg0 + 538);
    const auto *gpg0_539 = buffer.data(gpg0 + 539);
    const auto *gpg0_555 = buffer.data(gpg0 + 555);
    const auto *gpg0_556 = buffer.data(gpg0 + 556);
    const auto *gpg0_557 = buffer.data(gpg0 + 557);
    const auto *gpg0_558 = buffer.data(gpg0 + 558);
    const auto *gpg0_559 = buffer.data(gpg0 + 559);
    const auto *gpg0_560 = buffer.data(gpg0 + 560);
    const auto *gpg0_561 = buffer.data(gpg0 + 561);
    const auto *gpg0_562 = buffer.data(gpg0 + 562);
    const auto *gpg0_563 = buffer.data(gpg0 + 563);
    const auto *gpg0_564 = buffer.data(gpg0 + 564);
    const auto *gpg0_565 = buffer.data(gpg0 + 565);
    const auto *gpg0_566 = buffer.data(gpg0 + 566);
    const auto *gpg0_567 = buffer.data(gpg0 + 567);
    const auto *gpg0_568 = buffer.data(gpg0 + 568);
    const auto *gpg0_569 = buffer.data(gpg0 + 569);

    const auto *gpg1_512 = buffer.data(gpg1 + 512);
    const auto *gpg1_514 = buffer.data(gpg1 + 514);
    const auto *gpg1_515 = buffer.data(gpg1 + 515);
    const auto *gpg1_517 = buffer.data(gpg1 + 517);
    const auto *gpg1_518 = buffer.data(gpg1 + 518);
    const auto *gpg1_519 = buffer.data(gpg1 + 519);
    const auto *gpg1_521 = buffer.data(gpg1 + 521);
    const auto *gpg1_522 = buffer.data(gpg1 + 522);
    const auto *gpg1_523 = buffer.data(gpg1 + 523);
    const auto *gpg1_524 = buffer.data(gpg1 + 524);
    const auto *gpg1_525 = buffer.data(gpg1 + 525);
    const auto *gpg1_526 = buffer.data(gpg1 + 526);
    const auto *gpg1_527 = buffer.data(gpg1 + 527);
    const auto *gpg1_528 = buffer.data(gpg1 + 528);
    const auto *gpg1_529 = buffer.data(gpg1 + 529);
    const auto *gpg1_530 = buffer.data(gpg1 + 530);
    const auto *gpg1_531 = buffer.data(gpg1 + 531);
    const auto *gpg1_532 = buffer.data(gpg1 + 532);
    const auto *gpg1_533 = buffer.data(gpg1 + 533);
    const auto *gpg1_534 = buffer.data(gpg1 + 534);
    const auto *gpg1_535 = buffer.data(gpg1 + 535);
    const auto *gpg1_536 = buffer.data(gpg1 + 536);
    const auto *gpg1_537 = buffer.data(gpg1 + 537);
    const auto *gpg1_538 = buffer.data(gpg1 + 538);
    const auto *gpg1_539 = buffer.data(gpg1 + 539);
    const auto *gpg1_555 = buffer.data(gpg1 + 555);
    const auto *gpg1_556 = buffer.data(gpg1 + 556);
    const auto *gpg1_557 = buffer.data(gpg1 + 557);
    const auto *gpg1_558 = buffer.data(gpg1 + 558);
    const auto *gpg1_559 = buffer.data(gpg1 + 559);
    const auto *gpg1_560 = buffer.data(gpg1 + 560);
    const auto *gpg1_561 = buffer.data(gpg1 + 561);
    const auto *gpg1_562 = buffer.data(gpg1 + 562);
    const auto *gpg1_563 = buffer.data(gpg1 + 563);
    const auto *gpg1_564 = buffer.data(gpg1 + 564);
    const auto *gpg1_565 = buffer.data(gpg1 + 565);
    const auto *gpg1_566 = buffer.data(gpg1 + 566);
    const auto *gpg1_567 = buffer.data(gpg1 + 567);
    const auto *gpg1_568 = buffer.data(gpg1 + 568);
    const auto *gpg1_569 = buffer.data(gpg1 + 569);

    const auto *gph_708 = buffer.data(gph + 708);
    const auto *gph_709 = buffer.data(gph + 709);
    const auto *gph_710 = buffer.data(gph + 710);
    const auto *gph_711 = buffer.data(gph + 711);
    const auto *gph_712 = buffer.data(gph + 712);
    const auto *gph_713 = buffer.data(gph + 713);
    const auto *gph_716 = buffer.data(gph + 716);
    const auto *gph_718 = buffer.data(gph + 718);
    const auto *gph_719 = buffer.data(gph + 719);
    const auto *gph_721 = buffer.data(gph + 721);
    const auto *gph_722 = buffer.data(gph + 722);
    const auto *gph_723 = buffer.data(gph + 723);
    const auto *gph_725 = buffer.data(gph + 725);
    const auto *gph_726 = buffer.data(gph + 726);
    const auto *gph_727 = buffer.data(gph + 727);
    const auto *gph_728 = buffer.data(gph + 728);
    const auto *gph_729 = buffer.data(gph + 729);
    const auto *gph_730 = buffer.data(gph + 730);
    const auto *gph_731 = buffer.data(gph + 731);
    const auto *gph_732 = buffer.data(gph + 732);
    const auto *gph_733 = buffer.data(gph + 733);
    const auto *gph_734 = buffer.data(gph + 734);
    const auto *gph_735 = buffer.data(gph + 735);
    const auto *gph_736 = buffer.data(gph + 736);
    const auto *gph_737 = buffer.data(gph + 737);
    const auto *gph_738 = buffer.data(gph + 738);
    const auto *gph_739 = buffer.data(gph + 739);
    const auto *gph_740 = buffer.data(gph + 740);
    const auto *gph_741 = buffer.data(gph + 741);
    const auto *gph_742 = buffer.data(gph + 742);
    const auto *gph_743 = buffer.data(gph + 743);
    const auto *gph_744 = buffer.data(gph + 744);
    const auto *gph_745 = buffer.data(gph + 745);
    const auto *gph_746 = buffer.data(gph + 746);
    const auto *gph_747 = buffer.data(gph + 747);
    const auto *gph_748 = buffer.data(gph + 748);
    const auto *gph_749 = buffer.data(gph + 749);
    const auto *gph_750 = buffer.data(gph + 750);
    const auto *gph_751 = buffer.data(gph + 751);
    const auto *gph_752 = buffer.data(gph + 752);
    const auto *gph_753 = buffer.data(gph + 753);
    const auto *gph_754 = buffer.data(gph + 754);
    const auto *gph_755 = buffer.data(gph + 755);
    const auto *gph_771 = buffer.data(gph + 771);
    const auto *gph_772 = buffer.data(gph + 772);
    const auto *gph_773 = buffer.data(gph + 773);
    const auto *gph_774 = buffer.data(gph + 774);
    const auto *gph_775 = buffer.data(gph + 775);
    const auto *gph_776 = buffer.data(gph + 776);
    const auto *gph_777 = buffer.data(gph + 777);
    const auto *gph_778 = buffer.data(gph + 778);
    const auto *gph_779 = buffer.data(gph + 779);
    const auto *gph_780 = buffer.data(gph + 780);
    const auto *gph_781 = buffer.data(gph + 781);
    const auto *gph_782 = buffer.data(gph + 782);
    const auto *gph_783 = buffer.data(gph + 783);
    const auto *gph_784 = buffer.data(gph + 784);
    const auto *gph_785 = buffer.data(gph + 785);
    const auto *gph_786 = buffer.data(gph + 786);
    const auto *gph_787 = buffer.data(gph + 787);
    const auto *gph_788 = buffer.data(gph + 788);
    const auto *gph_789 = buffer.data(gph + 789);
    const auto *gph_790 = buffer.data(gph + 790);
    const auto *gph_791 = buffer.data(gph + 791);
    const auto *gph_792 = buffer.data(gph + 792);
    const auto *gph_793 = buffer.data(gph + 793);
    const auto *gph_794 = buffer.data(gph + 794);
    const auto *gph_795 = buffer.data(gph + 795);
    const auto *gph_796 = buffer.data(gph + 796);
    const auto *gph_797 = buffer.data(gph + 797);

#pragma omp simd aligned(t_939, t_940, t_941, t_942, t_943, pc_x, gsh_246, gsh_247, gsh_248, \
                         gsh_249, gsh_250, gph_708, gph_709, gph_710, gph_711, \
                         gph_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_1 * gsh_246[k]
                   + f_4 * pc_x[k] * gph_708[k];

        t_940[k] = f_1 * gsh_247[k]
                   + f_4 * pc_x[k] * gph_709[k];

        t_941[k] = f_1 * gsh_248[k]
                   + f_4 * pc_x[k] * gph_710[k];

        t_942[k] = f_1 * gsh_249[k]
                   + f_4 * pc_x[k] * gph_711[k];

        t_943[k] = f_1 * gsh_250[k]
                   + f_4 * pc_x[k] * gph_712[k];
    }

#pragma omp simd aligned(t_944, t_945, t_946, t_947, pa_z, pb_x, pc_x, pc_z, fpi0_525, \
                         fph_393, fpi1_525, gsi0_331, gsh_251, gsi1_331, gph_708, \
                         gph_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_1 * gsh_251[k]
                   + f_4 * pc_x[k] * gph_713[k];

        t_945[k] = pa_z[k] * fpi0_525[k]
                   - f_11 * pc_z[k] * fpi1_525[k];

        t_946[k] = f_1 * fph_393[k]
                   + f_4 * pc_z[k] * gph_708[k];

        t_947[k] = pb_x[k] * gsi0_331[k]
                   - f_11 * pc_x[k] * gsi1_331[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, t_951, pb_x, pc_x, pc_y, fph_461, gsi0_332, \
                         gsi0_333, gsi0_335, gsi1_332, gsi1_333, gsi1_335, \
                         gph_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = pb_x[k] * gsi0_332[k]
                   - f_11 * pc_x[k] * gsi1_332[k];

        t_949[k] = pb_x[k] * gsi0_333[k]
                   - f_11 * pc_x[k] * gsi1_333[k];

        t_950[k] = f_13 * fph_461[k]
                   + f_4 * pc_y[k] * gph_713[k];

        t_951[k] = pb_x[k] * gsi0_335[k]
                   - f_11 * pc_x[k] * gsi1_335[k];
    }

#pragma omp simd aligned(t_952, t_953, t_954, t_955, pa_z, pc_x, pc_z, fpi0_532, fpi0_533, \
                         fpi0_535, fpi1_532, fpi1_533, fpi1_535, gpg0_512, gpg1_512, \
                         gph_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_952[k] = pa_z[k] * fpi0_532[k]
                   - f_11 * pc_z[k] * fpi1_532[k];

        t_953[k] = pa_z[k] * fpi0_533[k]
                   - f_11 * pc_z[k] * fpi1_533[k];

        t_954[k] = f_17 * gpg0_512[k]
                   - f_18 * gpg1_512[k]
                   + f_4 * pc_x[k] * gph_716[k];

        t_955[k] = pa_z[k] * fpi0_535[k]
                   - f_11 * pc_z[k] * fpi1_535[k];
    }

#pragma omp simd aligned(t_956, t_957, t_958, pa_z, pc_x, pc_z, fpi0_538, fpi1_538, gpg0_514, \
                         gpg0_515, gpg1_514, gpg1_515, gph_718, \
                         gph_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_956[k] = f_9 * gpg0_514[k]
                   - f_10 * gpg1_514[k]
                   + f_4 * pc_x[k] * gph_718[k];

        t_957[k] = f_9 * gpg0_515[k]
                   - f_10 * gpg1_515[k]
                   + f_4 * pc_x[k] * gph_719[k];

        t_958[k] = pa_z[k] * fpi0_538[k]
                   - f_11 * pc_z[k] * fpi1_538[k];
    }

#pragma omp simd aligned(t_959, t_960, t_961, pc_x, gpg0_517, gpg0_518, gpg0_519, gpg1_517, \
                         gpg1_518, gpg1_519, gph_721, gph_722, \
                         gph_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_959[k] = f_7 * gpg0_517[k]
                   - f_8 * gpg1_517[k]
                   + f_4 * pc_x[k] * gph_721[k];

        t_960[k] = f_7 * gpg0_518[k]
                   - f_8 * gpg1_518[k]
                   + f_4 * pc_x[k] * gph_722[k];

        t_961[k] = f_7 * gpg0_519[k]
                   - f_8 * gpg1_519[k]
                   + f_4 * pc_x[k] * gph_723[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, pa_z, pc_x, pc_z, fpi0_542, fpi1_542, gpg0_521, \
                         gpg0_522, gpg1_521, gpg1_522, gph_725, \
                         gph_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = pa_z[k] * fpi0_542[k]
                   - f_11 * pc_z[k] * fpi1_542[k];

        t_963[k] = f_5 * gpg0_521[k]
                   - f_6 * gpg1_521[k]
                   + f_4 * pc_x[k] * gph_725[k];

        t_964[k] = f_5 * gpg0_522[k]
                   - f_6 * gpg1_522[k]
                   + f_4 * pc_x[k] * gph_726[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, pc_x, gpg0_523, gpg0_524, \
                         gpg1_523, gpg1_524, gph_727, gph_728, gph_729, gph_730, \
                         gph_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = f_5 * gpg0_523[k]
                   - f_6 * gpg1_523[k]
                   + f_4 * pc_x[k] * gph_727[k];

        t_966[k] = f_5 * gpg0_524[k]
                   - f_6 * gpg1_524[k]
                   + f_4 * pc_x[k] * gph_728[k];

        t_967[k] = f_4 * pc_x[k] * gph_729[k];

        t_968[k] = f_4 * pc_x[k] * gph_730[k];

        t_969[k] = f_4 * pc_x[k] * gph_731[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, t_974, pa_z, pc_x, pc_z, fpi0_553, \
                         fph_414, fpi1_553, gph_729, gph_732, gph_733, \
                         gph_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = f_4 * pc_x[k] * gph_732[k];

        t_971[k] = f_4 * pc_x[k] * gph_733[k];

        t_972[k] = f_4 * pc_x[k] * gph_734[k];

        t_973[k] = pa_z[k] * fpi0_553[k]
                   - f_11 * pc_z[k] * fpi1_553[k];

        t_974[k] = f_1 * fph_414[k]
                   + f_4 * pc_z[k] * gph_729[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, pa_z, pc_z, fpi0_555, fpi0_556, fpi0_557, \
                         fph_415, fph_416, fph_417, fpi1_555, fpi1_556, \
                         fpi1_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = pa_z[k] * fpi0_555[k]
                   + f_12 * fph_415[k]
                   - f_11 * pc_z[k] * fpi1_555[k];

        t_976[k] = pa_z[k] * fpi0_556[k]
                   + f_13 * fph_416[k]
                   - f_11 * pc_z[k] * fpi1_556[k];

        t_977[k] = pa_z[k] * fpi0_557[k]
                   + f_0 * fph_417[k]
                   - f_11 * pc_z[k] * fpi1_557[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, pc_x, pc_y, pc_z, fph_419, fph_482, gsh_251, \
                         gpg0_524, gpg0_525, gpg1_524, gpg1_525, gph_734, \
                         gph_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_13 * fph_482[k]
                   + f_1 * gsh_251[k]
                   + f_4 * pc_y[k] * gph_734[k];

        t_979[k] = f_1 * fph_419[k]
                   + f_2 * gpg0_524[k]
                   - f_3 * gpg1_524[k]
                   + f_4 * pc_z[k] * gph_734[k];

        t_980[k] = f_2 * gpg0_525[k]
                   - f_3 * gpg1_525[k]
                   + f_4 * pc_x[k] * gph_735[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, pc_x, gpg0_526, gpg0_527, gpg0_528, gpg1_526, \
                         gpg1_527, gpg1_528, gph_736, gph_737, \
                         gph_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_17 * gpg0_526[k]
                   - f_18 * gpg1_526[k]
                   + f_4 * pc_x[k] * gph_736[k];

        t_982[k] = f_17 * gpg0_527[k]
                   - f_18 * gpg1_527[k]
                   + f_4 * pc_x[k] * gph_737[k];

        t_983[k] = f_9 * gpg0_528[k]
                   - f_10 * gpg1_528[k]
                   + f_4 * pc_x[k] * gph_738[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, pc_x, gpg0_529, gpg0_530, gpg0_531, gpg1_529, \
                         gpg1_530, gpg1_531, gph_739, gph_740, \
                         gph_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_9 * gpg0_529[k]
                   - f_10 * gpg1_529[k]
                   + f_4 * pc_x[k] * gph_739[k];

        t_985[k] = f_9 * gpg0_530[k]
                   - f_10 * gpg1_530[k]
                   + f_4 * pc_x[k] * gph_740[k];

        t_986[k] = f_7 * gpg0_531[k]
                   - f_8 * gpg1_531[k]
                   + f_4 * pc_x[k] * gph_741[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, pc_x, gpg0_532, gpg0_533, gpg0_534, gpg1_532, \
                         gpg1_533, gpg1_534, gph_742, gph_743, \
                         gph_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_7 * gpg0_532[k]
                   - f_8 * gpg1_532[k]
                   + f_4 * pc_x[k] * gph_742[k];

        t_988[k] = f_7 * gpg0_533[k]
                   - f_8 * gpg1_533[k]
                   + f_4 * pc_x[k] * gph_743[k];

        t_989[k] = f_7 * gpg0_534[k]
                   - f_8 * gpg1_534[k]
                   + f_4 * pc_x[k] * gph_744[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, pc_x, gpg0_535, gpg0_536, gpg0_537, gpg1_535, \
                         gpg1_536, gpg1_537, gph_745, gph_746, \
                         gph_747 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = f_5 * gpg0_535[k]
                   - f_6 * gpg1_535[k]
                   + f_4 * pc_x[k] * gph_745[k];

        t_991[k] = f_5 * gpg0_536[k]
                   - f_6 * gpg1_536[k]
                   + f_4 * pc_x[k] * gph_746[k];

        t_992[k] = f_5 * gpg0_537[k]
                   - f_6 * gpg1_537[k]
                   + f_4 * pc_x[k] * gph_747[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, t_996, t_997, pc_x, gpg0_538, gpg0_539, \
                         gpg1_538, gpg1_539, gph_748, gph_749, gph_750, gph_751, \
                         gph_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = f_5 * gpg0_538[k]
                   - f_6 * gpg1_538[k]
                   + f_4 * pc_x[k] * gph_748[k];

        t_994[k] = f_5 * gpg0_539[k]
                   - f_6 * gpg1_539[k]
                   + f_4 * pc_x[k] * gph_749[k];

        t_995[k] = f_4 * pc_x[k] * gph_750[k];

        t_996[k] = f_4 * pc_x[k] * gph_751[k];

        t_997[k] = f_4 * pc_x[k] * gph_752[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, pc_x, pc_y, fph_498, gpg0_535, \
                         gpg1_535, gph_750, gph_753, gph_754, gph_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_4 * pc_x[k] * gph_753[k];

        t_999[k] = f_4 * pc_x[k] * gph_754[k];

        t_1000[k] = f_4 * pc_x[k] * gph_755[k];

        t_1001[k] = f_13 * fph_498[k]
                    + f_2 * gpg0_535[k]
                    - f_3 * gpg1_535[k]
                    + f_4 * pc_y[k] * gph_750[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, pc_z, fph_435, fph_500, fph_501, \
                         gsh_246, gpg0_537, gpg0_538, gpg1_537, gpg1_538, gph_750, gph_752, \
                         gph_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_1 * fph_435[k]
                    + f_1 * gsh_246[k]
                    + f_4 * pc_z[k] * gph_750[k];

        t_1003[k] = f_13 * fph_500[k]
                    + f_9 * gpg0_537[k]
                    - f_10 * gpg1_537[k]
                    + f_4 * pc_y[k] * gph_752[k];

        t_1004[k] = f_13 * fph_501[k]
                    + f_7 * gpg0_538[k]
                    - f_8 * gpg1_538[k]
                    + f_4 * pc_y[k] * gph_753[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, pa_y, pc_y, dpi0_419, dpi1_419, fpi0_671, \
                         fph_502, fph_503, fpi1_671, gpg0_539, gpg1_539, gph_754, \
                         gph_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_13 * fph_502[k]
                    + f_5 * gpg0_539[k]
                    - f_6 * gpg1_539[k]
                    + f_4 * pc_y[k] * gph_754[k];

        t_1006[k] = f_13 * fph_503[k]
                    + f_4 * pc_y[k] * gph_755[k];

        t_1007[k] = f_15 * dpi0_419[k]
                    - f_16 * dpi1_419[k]
                    + pa_y[k] * fpi0_671[k]
                    - f_11 * pc_y[k] * fpi1_671[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, pb_x, pc_x, gsi0_336, gsi0_337, gsi0_338, \
                         gsh_252, gsh_253, gsh_254, gsi1_336, gsi1_337, \
                         gsi1_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = pb_x[k] * gsi0_336[k]
                    + f_14 * gsh_252[k]
                    - f_11 * pc_x[k] * gsi1_336[k];

        t_1009[k] = pb_x[k] * gsi0_337[k]
                    + f_19 * gsh_253[k]
                    - f_11 * pc_x[k] * gsi1_337[k];

        t_1010[k] = pb_x[k] * gsi0_338[k]
                    + f_19 * gsh_254[k]
                    - f_11 * pc_x[k] * gsi1_338[k];
    }

#pragma omp simd aligned(t_1011, t_1012, t_1013, pb_x, pc_x, gsi0_339, gsi0_340, gsi0_341, \
                         gsh_255, gsh_256, gsh_257, gsi1_339, gsi1_340, \
                         gsi1_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1011[k] = pb_x[k] * gsi0_339[k]
                    + f_0 * gsh_255[k]
                    - f_11 * pc_x[k] * gsi1_339[k];

        t_1012[k] = pb_x[k] * gsi0_340[k]
                    + f_0 * gsh_256[k]
                    - f_11 * pc_x[k] * gsi1_340[k];

        t_1013[k] = pb_x[k] * gsi0_341[k]
                    + f_0 * gsh_257[k]
                    - f_11 * pc_x[k] * gsi1_341[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, pb_x, pc_x, gsi0_342, gsi0_343, gsi0_344, \
                         gsh_258, gsh_259, gsh_260, gsi1_342, gsi1_343, \
                         gsi1_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = pb_x[k] * gsi0_342[k]
                    + f_13 * gsh_258[k]
                    - f_11 * pc_x[k] * gsi1_342[k];

        t_1015[k] = pb_x[k] * gsi0_343[k]
                    + f_13 * gsh_259[k]
                    - f_11 * pc_x[k] * gsi1_343[k];

        t_1016[k] = pb_x[k] * gsi0_344[k]
                    + f_13 * gsh_260[k]
                    - f_11 * pc_x[k] * gsi1_344[k];
    }

#pragma omp simd aligned(t_1017, t_1018, t_1019, pb_x, pc_x, gsi0_345, gsi0_346, gsi0_347, \
                         gsh_261, gsh_262, gsh_263, gsi1_345, gsi1_346, \
                         gsi1_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1017[k] = pb_x[k] * gsi0_345[k]
                    + f_13 * gsh_261[k]
                    - f_11 * pc_x[k] * gsi1_345[k];

        t_1018[k] = pb_x[k] * gsi0_346[k]
                    + f_12 * gsh_262[k]
                    - f_11 * pc_x[k] * gsi1_346[k];

        t_1019[k] = pb_x[k] * gsi0_347[k]
                    + f_12 * gsh_263[k]
                    - f_11 * pc_x[k] * gsi1_347[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, pb_x, pc_x, gsi0_348, gsi0_349, gsi0_350, \
                         gsh_264, gsh_265, gsh_266, gsi1_348, gsi1_349, \
                         gsi1_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = pb_x[k] * gsi0_348[k]
                    + f_12 * gsh_264[k]
                    - f_11 * pc_x[k] * gsi1_348[k];

        t_1021[k] = pb_x[k] * gsi0_349[k]
                    + f_12 * gsh_265[k]
                    - f_11 * pc_x[k] * gsi1_349[k];

        t_1022[k] = pb_x[k] * gsi0_350[k]
                    + f_12 * gsh_266[k]
                    - f_11 * pc_x[k] * gsi1_350[k];
    }

#pragma omp simd aligned(t_1023, t_1024, t_1025, t_1026, t_1027, pc_x, gsh_267, gsh_268, \
                         gsh_269, gsh_270, gsh_271, gph_771, gph_772, gph_773, gph_774, \
                         gph_775 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1023[k] = f_1 * gsh_267[k]
                    + f_4 * pc_x[k] * gph_771[k];

        t_1024[k] = f_1 * gsh_268[k]
                    + f_4 * pc_x[k] * gph_772[k];

        t_1025[k] = f_1 * gsh_269[k]
                    + f_4 * pc_x[k] * gph_773[k];

        t_1026[k] = f_1 * gsh_270[k]
                    + f_4 * pc_x[k] * gph_774[k];

        t_1027[k] = f_1 * gsh_271[k]
                    + f_4 * pc_x[k] * gph_775[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, t_1031, pb_x, pc_x, pc_z, fph_456, gsi0_357, \
                         gsi0_359, gsh_272, gsi1_357, gsi1_359, gph_771, \
                         gph_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_1 * gsh_272[k]
                    + f_4 * pc_x[k] * gph_776[k];

        t_1029[k] = pb_x[k] * gsi0_357[k]
                    - f_11 * pc_x[k] * gsi1_357[k];

        t_1030[k] = f_12 * fph_456[k]
                    + f_4 * pc_z[k] * gph_771[k];

        t_1031[k] = pb_x[k] * gsi0_359[k]
                    - f_11 * pc_x[k] * gsi1_359[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, t_1035, pb_x, pc_x, pc_y, fph_524, gsi0_360, \
                         gsi0_361, gsi0_363, gsi1_360, gsi1_361, gsi1_363, \
                         gph_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = pb_x[k] * gsi0_360[k]
                    - f_11 * pc_x[k] * gsi1_360[k];

        t_1033[k] = pb_x[k] * gsi0_361[k]
                    - f_11 * pc_x[k] * gsi1_361[k];

        t_1034[k] = f_12 * fph_524[k]
                    + f_4 * pc_y[k] * gph_776[k];

        t_1035[k] = pb_x[k] * gsi0_363[k]
                    - f_11 * pc_x[k] * gsi1_363[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, pc_x, gpg0_555, gpg0_556, gpg0_557, gpg1_555, \
                         gpg1_556, gpg1_557, gph_777, gph_778, \
                         gph_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = f_2 * gpg0_555[k]
                    - f_3 * gpg1_555[k]
                    + f_4 * pc_x[k] * gph_777[k];

        t_1037[k] = f_17 * gpg0_556[k]
                    - f_18 * gpg1_556[k]
                    + f_4 * pc_x[k] * gph_778[k];

        t_1038[k] = f_17 * gpg0_557[k]
                    - f_18 * gpg1_557[k]
                    + f_4 * pc_x[k] * gph_779[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pc_x, gpg0_558, gpg0_559, gpg0_560, gpg1_558, \
                         gpg1_559, gpg1_560, gph_780, gph_781, \
                         gph_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_9 * gpg0_558[k]
                    - f_10 * gpg1_558[k]
                    + f_4 * pc_x[k] * gph_780[k];

        t_1040[k] = f_9 * gpg0_559[k]
                    - f_10 * gpg1_559[k]
                    + f_4 * pc_x[k] * gph_781[k];

        t_1041[k] = f_9 * gpg0_560[k]
                    - f_10 * gpg1_560[k]
                    + f_4 * pc_x[k] * gph_782[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, pc_x, gpg0_561, gpg0_562, gpg0_563, gpg1_561, \
                         gpg1_562, gpg1_563, gph_783, gph_784, \
                         gph_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_7 * gpg0_561[k]
                    - f_8 * gpg1_561[k]
                    + f_4 * pc_x[k] * gph_783[k];

        t_1043[k] = f_7 * gpg0_562[k]
                    - f_8 * gpg1_562[k]
                    + f_4 * pc_x[k] * gph_784[k];

        t_1044[k] = f_7 * gpg0_563[k]
                    - f_8 * gpg1_563[k]
                    + f_4 * pc_x[k] * gph_785[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, pc_x, gpg0_564, gpg0_565, gpg0_566, gpg1_564, \
                         gpg1_565, gpg1_566, gph_786, gph_787, \
                         gph_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = f_7 * gpg0_564[k]
                    - f_8 * gpg1_564[k]
                    + f_4 * pc_x[k] * gph_786[k];

        t_1046[k] = f_5 * gpg0_565[k]
                    - f_6 * gpg1_565[k]
                    + f_4 * pc_x[k] * gph_787[k];

        t_1047[k] = f_5 * gpg0_566[k]
                    - f_6 * gpg1_566[k]
                    + f_4 * pc_x[k] * gph_788[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, t_1051, pc_x, gpg0_567, gpg0_568, gpg0_569, \
                         gpg1_567, gpg1_568, gpg1_569, gph_789, gph_790, gph_791, \
                         gph_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = f_5 * gpg0_567[k]
                    - f_6 * gpg1_567[k]
                    + f_4 * pc_x[k] * gph_789[k];

        t_1049[k] = f_5 * gpg0_568[k]
                    - f_6 * gpg1_568[k]
                    + f_4 * pc_x[k] * gph_790[k];

        t_1050[k] = f_5 * gpg0_569[k]
                    - f_6 * gpg1_569[k]
                    + f_4 * pc_x[k] * gph_791[k];

        t_1051[k] = f_4 * pc_x[k] * gph_792[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, t_1055, t_1056, pc_x, gph_793, gph_794, \
                         gph_795, gph_796, gph_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_4 * pc_x[k] * gph_793[k];

        t_1053[k] = f_4 * pc_x[k] * gph_794[k];

        t_1054[k] = f_4 * pc_x[k] * gph_795[k];

        t_1055[k] = f_4 * pc_x[k] * gph_796[k];

        t_1056[k] = f_4 * pc_x[k] * gph_797[k];
    }

#pragma omp simd aligned(t_1057, t_1058, pa_z, pc_z, dpi0_301, dpi1_301, fpi0_637, fph_477, \
                         fpi1_637, gph_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = f_20 * dpi0_301[k]
                    - f_21 * dpi1_301[k]
                    + pa_z[k] * fpi0_637[k]
                    - f_11 * pc_z[k] * fpi1_637[k];

        t_1058[k] = f_12 * fph_477[k]
                    + f_4 * pc_z[k] * gph_792[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dpi1,
                                                          const size_t fpi0, const size_t fph,
                                                          const size_t fpi1, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t gpg0, const size_t gpg1,
                                                          const size_t gph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_14 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_20 = 0.5 / p;
    const auto f_21 = 0.5 * gamma / (p * q);

    auto *t_1059 = buffer.data(target + 1059);
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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpi0_503 = buffer.data(dpi0 + 503);

    const auto *dpi1_503 = buffer.data(dpi1 + 503);

    const auto *fpi0_755 = buffer.data(fpi0 + 755);
    const auto *fpi0_756 = buffer.data(fpi0 + 756);
    const auto *fpi0_757 = buffer.data(fpi0 + 757);
    const auto *fpi0_758 = buffer.data(fpi0 + 758);
    const auto *fpi0_759 = buffer.data(fpi0 + 759);
    const auto *fpi0_760 = buffer.data(fpi0 + 760);
    const auto *fpi0_761 = buffer.data(fpi0 + 761);
    const auto *fpi0_762 = buffer.data(fpi0 + 762);
    const auto *fpi0_763 = buffer.data(fpi0 + 763);
    const auto *fpi0_764 = buffer.data(fpi0 + 764);
    const auto *fpi0_765 = buffer.data(fpi0 + 765);
    const auto *fpi0_766 = buffer.data(fpi0 + 766);
    const auto *fpi0_767 = buffer.data(fpi0 + 767);
    const auto *fpi0_768 = buffer.data(fpi0 + 768);
    const auto *fpi0_769 = buffer.data(fpi0 + 769);
    const auto *fpi0_770 = buffer.data(fpi0 + 770);
    const auto *fpi0_783 = buffer.data(fpi0 + 783);
    const auto *fpi0_812 = buffer.data(fpi0 + 812);
    const auto *fpi0_814 = buffer.data(fpi0 + 814);
    const auto *fpi0_817 = buffer.data(fpi0 + 817);
    const auto *fpi0_821 = buffer.data(fpi0 + 821);
    const auto *fpi0_826 = buffer.data(fpi0 + 826);
    const auto *fpi0_833 = buffer.data(fpi0 + 833);
    const auto *fpi0_835 = buffer.data(fpi0 + 835);
    const auto *fpi0_836 = buffer.data(fpi0 + 836);

    const auto *fph_482 = buffer.data(fph + 482);
    const auto *fph_498 = buffer.data(fph + 498);
    const auto *fph_519 = buffer.data(fph + 519);
    const auto *fph_540 = buffer.data(fph + 540);
    const auto *fph_542 = buffer.data(fph + 542);
    const auto *fph_543 = buffer.data(fph + 543);
    const auto *fph_544 = buffer.data(fph + 544);
    const auto *fph_545 = buffer.data(fph + 545);
    const auto *fph_561 = buffer.data(fph + 561);
    const auto *fph_563 = buffer.data(fph + 563);
    const auto *fph_564 = buffer.data(fph + 564);
    const auto *fph_565 = buffer.data(fph + 565);
    const auto *fph_566 = buffer.data(fph + 566);
    const auto *fph_567 = buffer.data(fph + 567);
    const auto *fph_568 = buffer.data(fph + 568);
    const auto *fph_569 = buffer.data(fph + 569);
    const auto *fph_570 = buffer.data(fph + 570);
    const auto *fph_571 = buffer.data(fph + 571);
    const auto *fph_572 = buffer.data(fph + 572);
    const auto *fph_573 = buffer.data(fph + 573);
    const auto *fph_574 = buffer.data(fph + 574);
    const auto *fph_575 = buffer.data(fph + 575);
    const auto *fph_576 = buffer.data(fph + 576);
    const auto *fph_587 = buffer.data(fph + 587);
    const auto *fph_603 = buffer.data(fph + 603);
    const auto *fph_605 = buffer.data(fph + 605);
    const auto *fph_606 = buffer.data(fph + 606);
    const auto *fph_607 = buffer.data(fph + 607);
    const auto *fph_608 = buffer.data(fph + 608);
    const auto *fph_624 = buffer.data(fph + 624);
    const auto *fph_626 = buffer.data(fph + 626);
    const auto *fph_627 = buffer.data(fph + 627);

    const auto *fpi1_755 = buffer.data(fpi1 + 755);
    const auto *fpi1_756 = buffer.data(fpi1 + 756);
    const auto *fpi1_757 = buffer.data(fpi1 + 757);
    const auto *fpi1_758 = buffer.data(fpi1 + 758);
    const auto *fpi1_759 = buffer.data(fpi1 + 759);
    const auto *fpi1_760 = buffer.data(fpi1 + 760);
    const auto *fpi1_761 = buffer.data(fpi1 + 761);
    const auto *fpi1_762 = buffer.data(fpi1 + 762);
    const auto *fpi1_763 = buffer.data(fpi1 + 763);
    const auto *fpi1_764 = buffer.data(fpi1 + 764);
    const auto *fpi1_765 = buffer.data(fpi1 + 765);
    const auto *fpi1_766 = buffer.data(fpi1 + 766);
    const auto *fpi1_767 = buffer.data(fpi1 + 767);
    const auto *fpi1_768 = buffer.data(fpi1 + 768);
    const auto *fpi1_769 = buffer.data(fpi1 + 769);
    const auto *fpi1_770 = buffer.data(fpi1 + 770);
    const auto *fpi1_783 = buffer.data(fpi1 + 783);
    const auto *fpi1_812 = buffer.data(fpi1 + 812);
    const auto *fpi1_814 = buffer.data(fpi1 + 814);
    const auto *fpi1_817 = buffer.data(fpi1 + 817);
    const auto *fpi1_821 = buffer.data(fpi1 + 821);
    const auto *fpi1_826 = buffer.data(fpi1 + 826);
    const auto *fpi1_833 = buffer.data(fpi1 + 833);
    const auto *fpi1_835 = buffer.data(fpi1 + 835);
    const auto *fpi1_836 = buffer.data(fpi1 + 836);

    const auto *gsi0_385 = buffer.data(gsi0 + 385);
    const auto *gsi0_387 = buffer.data(gsi0 + 387);
    const auto *gsi0_388 = buffer.data(gsi0 + 388);
    const auto *gsi0_389 = buffer.data(gsi0 + 389);

    const auto *gsh_267 = buffer.data(gsh + 267);
    const auto *gsh_269 = buffer.data(gsh + 269);
    const auto *gsh_270 = buffer.data(gsh + 270);
    const auto *gsh_271 = buffer.data(gsh + 271);
    const auto *gsh_272 = buffer.data(gsh + 272);
    const auto *gsh_288 = buffer.data(gsh + 288);
    const auto *gsh_289 = buffer.data(gsh + 289);
    const auto *gsh_290 = buffer.data(gsh + 290);
    const auto *gsh_291 = buffer.data(gsh + 291);
    const auto *gsh_292 = buffer.data(gsh + 292);
    const auto *gsh_293 = buffer.data(gsh + 293);

    const auto *gsi1_385 = buffer.data(gsi1 + 385);
    const auto *gsi1_387 = buffer.data(gsi1 + 387);
    const auto *gsi1_388 = buffer.data(gsi1 + 388);
    const auto *gsi1_389 = buffer.data(gsi1 + 389);

    const auto *gpg0_567 = buffer.data(gpg0 + 567);
    const auto *gpg0_568 = buffer.data(gpg0 + 568);
    const auto *gpg0_569 = buffer.data(gpg0 + 569);
    const auto *gpg0_570 = buffer.data(gpg0 + 570);
    const auto *gpg0_571 = buffer.data(gpg0 + 571);
    const auto *gpg0_572 = buffer.data(gpg0 + 572);
    const auto *gpg0_573 = buffer.data(gpg0 + 573);
    const auto *gpg0_574 = buffer.data(gpg0 + 574);
    const auto *gpg0_575 = buffer.data(gpg0 + 575);
    const auto *gpg0_576 = buffer.data(gpg0 + 576);
    const auto *gpg0_577 = buffer.data(gpg0 + 577);
    const auto *gpg0_578 = buffer.data(gpg0 + 578);
    const auto *gpg0_579 = buffer.data(gpg0 + 579);
    const auto *gpg0_580 = buffer.data(gpg0 + 580);
    const auto *gpg0_581 = buffer.data(gpg0 + 581);
    const auto *gpg0_582 = buffer.data(gpg0 + 582);
    const auto *gpg0_583 = buffer.data(gpg0 + 583);
    const auto *gpg0_584 = buffer.data(gpg0 + 584);
    const auto *gpg0_600 = buffer.data(gpg0 + 600);
    const auto *gpg0_601 = buffer.data(gpg0 + 601);
    const auto *gpg0_602 = buffer.data(gpg0 + 602);
    const auto *gpg0_603 = buffer.data(gpg0 + 603);
    const auto *gpg0_604 = buffer.data(gpg0 + 604);
    const auto *gpg0_605 = buffer.data(gpg0 + 605);
    const auto *gpg0_606 = buffer.data(gpg0 + 606);
    const auto *gpg0_607 = buffer.data(gpg0 + 607);
    const auto *gpg0_608 = buffer.data(gpg0 + 608);
    const auto *gpg0_609 = buffer.data(gpg0 + 609);
    const auto *gpg0_610 = buffer.data(gpg0 + 610);
    const auto *gpg0_611 = buffer.data(gpg0 + 611);
    const auto *gpg0_612 = buffer.data(gpg0 + 612);
    const auto *gpg0_613 = buffer.data(gpg0 + 613);
    const auto *gpg0_614 = buffer.data(gpg0 + 614);
    const auto *gpg0_616 = buffer.data(gpg0 + 616);
    const auto *gpg0_618 = buffer.data(gpg0 + 618);
    const auto *gpg0_619 = buffer.data(gpg0 + 619);
    const auto *gpg0_621 = buffer.data(gpg0 + 621);
    const auto *gpg0_622 = buffer.data(gpg0 + 622);
    const auto *gpg0_623 = buffer.data(gpg0 + 623);
    const auto *gpg0_625 = buffer.data(gpg0 + 625);
    const auto *gpg0_626 = buffer.data(gpg0 + 626);
    const auto *gpg0_627 = buffer.data(gpg0 + 627);
    const auto *gpg0_628 = buffer.data(gpg0 + 628);

    const auto *gpg1_567 = buffer.data(gpg1 + 567);
    const auto *gpg1_568 = buffer.data(gpg1 + 568);
    const auto *gpg1_569 = buffer.data(gpg1 + 569);
    const auto *gpg1_570 = buffer.data(gpg1 + 570);
    const auto *gpg1_571 = buffer.data(gpg1 + 571);
    const auto *gpg1_572 = buffer.data(gpg1 + 572);
    const auto *gpg1_573 = buffer.data(gpg1 + 573);
    const auto *gpg1_574 = buffer.data(gpg1 + 574);
    const auto *gpg1_575 = buffer.data(gpg1 + 575);
    const auto *gpg1_576 = buffer.data(gpg1 + 576);
    const auto *gpg1_577 = buffer.data(gpg1 + 577);
    const auto *gpg1_578 = buffer.data(gpg1 + 578);
    const auto *gpg1_579 = buffer.data(gpg1 + 579);
    const auto *gpg1_580 = buffer.data(gpg1 + 580);
    const auto *gpg1_581 = buffer.data(gpg1 + 581);
    const auto *gpg1_582 = buffer.data(gpg1 + 582);
    const auto *gpg1_583 = buffer.data(gpg1 + 583);
    const auto *gpg1_584 = buffer.data(gpg1 + 584);
    const auto *gpg1_600 = buffer.data(gpg1 + 600);
    const auto *gpg1_601 = buffer.data(gpg1 + 601);
    const auto *gpg1_602 = buffer.data(gpg1 + 602);
    const auto *gpg1_603 = buffer.data(gpg1 + 603);
    const auto *gpg1_604 = buffer.data(gpg1 + 604);
    const auto *gpg1_605 = buffer.data(gpg1 + 605);
    const auto *gpg1_606 = buffer.data(gpg1 + 606);
    const auto *gpg1_607 = buffer.data(gpg1 + 607);
    const auto *gpg1_608 = buffer.data(gpg1 + 608);
    const auto *gpg1_609 = buffer.data(gpg1 + 609);
    const auto *gpg1_610 = buffer.data(gpg1 + 610);
    const auto *gpg1_611 = buffer.data(gpg1 + 611);
    const auto *gpg1_612 = buffer.data(gpg1 + 612);
    const auto *gpg1_613 = buffer.data(gpg1 + 613);
    const auto *gpg1_614 = buffer.data(gpg1 + 614);
    const auto *gpg1_616 = buffer.data(gpg1 + 616);
    const auto *gpg1_618 = buffer.data(gpg1 + 618);
    const auto *gpg1_619 = buffer.data(gpg1 + 619);
    const auto *gpg1_621 = buffer.data(gpg1 + 621);
    const auto *gpg1_622 = buffer.data(gpg1 + 622);
    const auto *gpg1_623 = buffer.data(gpg1 + 623);
    const auto *gpg1_625 = buffer.data(gpg1 + 625);
    const auto *gpg1_626 = buffer.data(gpg1 + 626);
    const auto *gpg1_627 = buffer.data(gpg1 + 627);
    const auto *gpg1_628 = buffer.data(gpg1 + 628);

    const auto *gph_794 = buffer.data(gph + 794);
    const auto *gph_795 = buffer.data(gph + 795);
    const auto *gph_796 = buffer.data(gph + 796);
    const auto *gph_797 = buffer.data(gph + 797);
    const auto *gph_798 = buffer.data(gph + 798);
    const auto *gph_799 = buffer.data(gph + 799);
    const auto *gph_800 = buffer.data(gph + 800);
    const auto *gph_801 = buffer.data(gph + 801);
    const auto *gph_802 = buffer.data(gph + 802);
    const auto *gph_803 = buffer.data(gph + 803);
    const auto *gph_804 = buffer.data(gph + 804);
    const auto *gph_805 = buffer.data(gph + 805);
    const auto *gph_806 = buffer.data(gph + 806);
    const auto *gph_807 = buffer.data(gph + 807);
    const auto *gph_808 = buffer.data(gph + 808);
    const auto *gph_809 = buffer.data(gph + 809);
    const auto *gph_810 = buffer.data(gph + 810);
    const auto *gph_811 = buffer.data(gph + 811);
    const auto *gph_812 = buffer.data(gph + 812);
    const auto *gph_813 = buffer.data(gph + 813);
    const auto *gph_814 = buffer.data(gph + 814);
    const auto *gph_815 = buffer.data(gph + 815);
    const auto *gph_816 = buffer.data(gph + 816);
    const auto *gph_817 = buffer.data(gph + 817);
    const auto *gph_818 = buffer.data(gph + 818);
    const auto *gph_834 = buffer.data(gph + 834);
    const auto *gph_835 = buffer.data(gph + 835);
    const auto *gph_836 = buffer.data(gph + 836);
    const auto *gph_837 = buffer.data(gph + 837);
    const auto *gph_838 = buffer.data(gph + 838);
    const auto *gph_839 = buffer.data(gph + 839);
    const auto *gph_840 = buffer.data(gph + 840);
    const auto *gph_841 = buffer.data(gph + 841);
    const auto *gph_842 = buffer.data(gph + 842);
    const auto *gph_843 = buffer.data(gph + 843);
    const auto *gph_844 = buffer.data(gph + 844);
    const auto *gph_845 = buffer.data(gph + 845);
    const auto *gph_846 = buffer.data(gph + 846);
    const auto *gph_847 = buffer.data(gph + 847);
    const auto *gph_848 = buffer.data(gph + 848);
    const auto *gph_849 = buffer.data(gph + 849);
    const auto *gph_850 = buffer.data(gph + 850);
    const auto *gph_851 = buffer.data(gph + 851);
    const auto *gph_852 = buffer.data(gph + 852);
    const auto *gph_853 = buffer.data(gph + 853);
    const auto *gph_854 = buffer.data(gph + 854);
    const auto *gph_855 = buffer.data(gph + 855);
    const auto *gph_856 = buffer.data(gph + 856);
    const auto *gph_857 = buffer.data(gph + 857);
    const auto *gph_858 = buffer.data(gph + 858);
    const auto *gph_859 = buffer.data(gph + 859);
    const auto *gph_860 = buffer.data(gph + 860);
    const auto *gph_862 = buffer.data(gph + 862);
    const auto *gph_864 = buffer.data(gph + 864);
    const auto *gph_865 = buffer.data(gph + 865);
    const auto *gph_867 = buffer.data(gph + 867);
    const auto *gph_868 = buffer.data(gph + 868);
    const auto *gph_869 = buffer.data(gph + 869);
    const auto *gph_871 = buffer.data(gph + 871);
    const auto *gph_872 = buffer.data(gph + 872);
    const auto *gph_873 = buffer.data(gph + 873);
    const auto *gph_874 = buffer.data(gph + 874);
    const auto *gph_876 = buffer.data(gph + 876);
    const auto *gph_877 = buffer.data(gph + 877);
    const auto *gph_878 = buffer.data(gph + 878);
    const auto *gph_879 = buffer.data(gph + 879);
    const auto *gph_880 = buffer.data(gph + 880);
    const auto *gph_881 = buffer.data(gph + 881);

#pragma omp simd aligned(t_1059, t_1060, pc_y, fph_542, fph_543, gsh_269, gsh_270, gpg0_567, \
                         gpg0_568, gpg1_567, gpg1_568, gph_794, \
                         gph_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = f_12 * fph_542[k]
                    + f_1 * gsh_269[k]
                    + f_9 * gpg0_567[k]
                    - f_10 * gpg1_567[k]
                    + f_4 * pc_y[k] * gph_794[k];

        t_1060[k] = f_12 * fph_543[k]
                    + f_1 * gsh_270[k]
                    + f_7 * gpg0_568[k]
                    - f_8 * gpg1_568[k]
                    + f_4 * pc_y[k] * gph_795[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pc_y, pc_z, fph_482, fph_544, fph_545, \
                         gsh_271, gsh_272, gpg0_569, gpg1_569, gph_796, \
                         gph_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_12 * fph_544[k]
                    + f_1 * gsh_271[k]
                    + f_5 * gpg0_569[k]
                    - f_6 * gpg1_569[k]
                    + f_4 * pc_y[k] * gph_796[k];

        t_1062[k] = f_12 * fph_545[k]
                    + f_1 * gsh_272[k]
                    + f_4 * pc_y[k] * gph_797[k];

        t_1063[k] = f_12 * fph_482[k]
                    + f_2 * gpg0_569[k]
                    - f_3 * gpg1_569[k]
                    + f_4 * pc_z[k] * gph_797[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, pc_x, gpg0_570, gpg0_571, gpg0_572, gpg1_570, \
                         gpg1_571, gpg1_572, gph_798, gph_799, \
                         gph_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_2 * gpg0_570[k]
                    - f_3 * gpg1_570[k]
                    + f_4 * pc_x[k] * gph_798[k];

        t_1065[k] = f_17 * gpg0_571[k]
                    - f_18 * gpg1_571[k]
                    + f_4 * pc_x[k] * gph_799[k];

        t_1066[k] = f_17 * gpg0_572[k]
                    - f_18 * gpg1_572[k]
                    + f_4 * pc_x[k] * gph_800[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, pc_x, gpg0_573, gpg0_574, gpg0_575, gpg1_573, \
                         gpg1_574, gpg1_575, gph_801, gph_802, \
                         gph_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_9 * gpg0_573[k]
                    - f_10 * gpg1_573[k]
                    + f_4 * pc_x[k] * gph_801[k];

        t_1068[k] = f_9 * gpg0_574[k]
                    - f_10 * gpg1_574[k]
                    + f_4 * pc_x[k] * gph_802[k];

        t_1069[k] = f_9 * gpg0_575[k]
                    - f_10 * gpg1_575[k]
                    + f_4 * pc_x[k] * gph_803[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, gpg0_576, gpg0_577, gpg0_578, gpg1_576, \
                         gpg1_577, gpg1_578, gph_804, gph_805, \
                         gph_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_7 * gpg0_576[k]
                    - f_8 * gpg1_576[k]
                    + f_4 * pc_x[k] * gph_804[k];

        t_1071[k] = f_7 * gpg0_577[k]
                    - f_8 * gpg1_577[k]
                    + f_4 * pc_x[k] * gph_805[k];

        t_1072[k] = f_7 * gpg0_578[k]
                    - f_8 * gpg1_578[k]
                    + f_4 * pc_x[k] * gph_806[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_x, gpg0_579, gpg0_580, gpg0_581, gpg1_579, \
                         gpg1_580, gpg1_581, gph_807, gph_808, \
                         gph_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_7 * gpg0_579[k]
                    - f_8 * gpg1_579[k]
                    + f_4 * pc_x[k] * gph_807[k];

        t_1074[k] = f_5 * gpg0_580[k]
                    - f_6 * gpg1_580[k]
                    + f_4 * pc_x[k] * gph_808[k];

        t_1075[k] = f_5 * gpg0_581[k]
                    - f_6 * gpg1_581[k]
                    + f_4 * pc_x[k] * gph_809[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, t_1079, pc_x, gpg0_582, gpg0_583, gpg0_584, \
                         gpg1_582, gpg1_583, gpg1_584, gph_810, gph_811, gph_812, \
                         gph_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_5 * gpg0_582[k]
                    - f_6 * gpg1_582[k]
                    + f_4 * pc_x[k] * gph_810[k];

        t_1077[k] = f_5 * gpg0_583[k]
                    - f_6 * gpg1_583[k]
                    + f_4 * pc_x[k] * gph_811[k];

        t_1078[k] = f_5 * gpg0_584[k]
                    - f_6 * gpg1_584[k]
                    + f_4 * pc_x[k] * gph_812[k];

        t_1079[k] = f_4 * pc_x[k] * gph_813[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, t_1084, pc_x, gph_814, gph_815, \
                         gph_816, gph_817, gph_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = f_4 * pc_x[k] * gph_814[k];

        t_1081[k] = f_4 * pc_x[k] * gph_815[k];

        t_1082[k] = f_4 * pc_x[k] * gph_816[k];

        t_1083[k] = f_4 * pc_x[k] * gph_817[k];

        t_1084[k] = f_4 * pc_x[k] * gph_818[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_y, pc_z, fph_498, fph_561, fph_563, \
                         gsh_267, gpg0_580, gpg0_582, gpg1_580, gpg1_582, gph_813, \
                         gph_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_12 * fph_561[k]
                    + f_2 * gpg0_580[k]
                    - f_3 * gpg1_580[k]
                    + f_4 * pc_y[k] * gph_813[k];

        t_1086[k] = f_12 * fph_498[k]
                    + f_1 * gsh_267[k]
                    + f_4 * pc_z[k] * gph_813[k];

        t_1087[k] = f_12 * fph_563[k]
                    + f_9 * gpg0_582[k]
                    - f_10 * gpg1_582[k]
                    + f_4 * pc_y[k] * gph_815[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_y, fph_564, fph_565, fph_566, gpg0_583, \
                         gpg0_584, gpg1_583, gpg1_584, gph_816, gph_817, \
                         gph_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_12 * fph_564[k]
                    + f_7 * gpg0_583[k]
                    - f_8 * gpg1_583[k]
                    + f_4 * pc_y[k] * gph_816[k];

        t_1089[k] = f_12 * fph_565[k]
                    + f_5 * gpg0_584[k]
                    - f_6 * gpg1_584[k]
                    + f_4 * pc_y[k] * gph_817[k];

        t_1090[k] = f_12 * fph_566[k]
                    + f_4 * pc_y[k] * gph_818[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pa_y, pc_y, dpi0_503, dpi1_503, fpi0_755, \
                         fpi0_756, fpi0_757, fph_567, fpi1_755, fpi1_756, \
                         fpi1_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_20 * dpi0_503[k]
                    - f_21 * dpi1_503[k]
                    + pa_y[k] * fpi0_755[k]
                    - f_11 * pc_y[k] * fpi1_755[k];

        t_1092[k] = pa_y[k] * fpi0_756[k]
                    - f_11 * pc_y[k] * fpi1_756[k];

        t_1093[k] = pa_y[k] * fpi0_757[k]
                    + f_1 * fph_567[k]
                    - f_11 * pc_y[k] * fpi1_757[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, t_1097, pa_y, pc_y, fpi0_758, fpi0_759, \
                         fpi0_760, fpi0_761, fph_568, fph_569, fpi1_758, fpi1_759, fpi1_760, \
                         fpi1_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = pa_y[k] * fpi0_758[k]
                    - f_11 * pc_y[k] * fpi1_758[k];

        t_1095[k] = pa_y[k] * fpi0_759[k]
                    + f_12 * fph_568[k]
                    - f_11 * pc_y[k] * fpi1_759[k];

        t_1096[k] = pa_y[k] * fpi0_760[k]
                    + f_1 * fph_569[k]
                    - f_11 * pc_y[k] * fpi1_760[k];

        t_1097[k] = pa_y[k] * fpi0_761[k]
                    - f_11 * pc_y[k] * fpi1_761[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, pa_y, pc_y, fpi0_762, fpi0_763, fpi0_764, \
                         fph_570, fph_571, fph_572, fpi1_762, fpi1_763, \
                         fpi1_764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = pa_y[k] * fpi0_762[k]
                    + f_13 * fph_570[k]
                    - f_11 * pc_y[k] * fpi1_762[k];

        t_1099[k] = pa_y[k] * fpi0_763[k]
                    + f_12 * fph_571[k]
                    - f_11 * pc_y[k] * fpi1_763[k];

        t_1100[k] = pa_y[k] * fpi0_764[k]
                    + f_1 * fph_572[k]
                    - f_11 * pc_y[k] * fpi1_764[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, pa_y, pc_y, fpi0_765, fpi0_766, fpi0_767, \
                         fph_573, fph_574, fpi1_765, fpi1_766, \
                         fpi1_767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = pa_y[k] * fpi0_765[k]
                    - f_11 * pc_y[k] * fpi1_765[k];

        t_1102[k] = pa_y[k] * fpi0_766[k]
                    + f_0 * fph_573[k]
                    - f_11 * pc_y[k] * fpi1_766[k];

        t_1103[k] = pa_y[k] * fpi0_767[k]
                    + f_13 * fph_574[k]
                    - f_11 * pc_y[k] * fpi1_767[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, pa_y, pc_y, fpi0_768, fpi0_769, fpi0_770, \
                         fph_575, fph_576, fpi1_768, fpi1_769, \
                         fpi1_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = pa_y[k] * fpi0_768[k]
                    + f_12 * fph_575[k]
                    - f_11 * pc_y[k] * fpi1_768[k];

        t_1105[k] = pa_y[k] * fpi0_769[k]
                    + f_1 * fph_576[k]
                    - f_11 * pc_y[k] * fpi1_769[k];

        t_1106[k] = pa_y[k] * fpi0_770[k]
                    - f_11 * pc_y[k] * fpi1_770[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, t_1110, t_1111, pc_x, gsh_288, gsh_289, \
                         gsh_290, gsh_291, gsh_292, gph_834, gph_835, gph_836, gph_837, \
                         gph_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = f_1 * gsh_288[k]
                    + f_4 * pc_x[k] * gph_834[k];

        t_1108[k] = f_1 * gsh_289[k]
                    + f_4 * pc_x[k] * gph_835[k];

        t_1109[k] = f_1 * gsh_290[k]
                    + f_4 * pc_x[k] * gph_836[k];

        t_1110[k] = f_1 * gsh_291[k]
                    + f_4 * pc_x[k] * gph_837[k];

        t_1111[k] = f_1 * gsh_292[k]
                    + f_4 * pc_x[k] * gph_838[k];
    }

#pragma omp simd aligned(t_1112, t_1113, t_1114, t_1115, pb_x, pc_x, pc_z, fph_519, gsi0_385, \
                         gsi0_387, gsh_293, gsi1_385, gsi1_387, gph_834, \
                         gph_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1112[k] = f_1 * gsh_293[k]
                    + f_4 * pc_x[k] * gph_839[k];

        t_1113[k] = pb_x[k] * gsi0_385[k]
                    - f_11 * pc_x[k] * gsi1_385[k];

        t_1114[k] = f_13 * fph_519[k]
                    + f_4 * pc_z[k] * gph_834[k];

        t_1115[k] = pb_x[k] * gsi0_387[k]
                    - f_11 * pc_x[k] * gsi1_387[k];
    }

#pragma omp simd aligned(t_1116, t_1117, t_1118, t_1119, pa_y, pb_x, pc_x, pc_y, fpi0_783, \
                         fph_587, fpi1_783, gsi0_388, gsi0_389, gsi1_388, gsi1_389, \
                         gph_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1116[k] = pb_x[k] * gsi0_388[k]
                    - f_11 * pc_x[k] * gsi1_388[k];

        t_1117[k] = pb_x[k] * gsi0_389[k]
                    - f_11 * pc_x[k] * gsi1_389[k];

        t_1118[k] = f_1 * fph_587[k]
                    + f_4 * pc_y[k] * gph_839[k];

        t_1119[k] = pa_y[k] * fpi0_783[k]
                    - f_11 * pc_y[k] * fpi1_783[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, pc_x, gpg0_600, gpg0_601, gpg0_602, gpg1_600, \
                         gpg1_601, gpg1_602, gph_840, gph_841, \
                         gph_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = f_2 * gpg0_600[k]
                    - f_3 * gpg1_600[k]
                    + f_4 * pc_x[k] * gph_840[k];

        t_1121[k] = f_17 * gpg0_601[k]
                    - f_18 * gpg1_601[k]
                    + f_4 * pc_x[k] * gph_841[k];

        t_1122[k] = f_17 * gpg0_602[k]
                    - f_18 * gpg1_602[k]
                    + f_4 * pc_x[k] * gph_842[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, pc_x, gpg0_603, gpg0_604, gpg0_605, gpg1_603, \
                         gpg1_604, gpg1_605, gph_843, gph_844, \
                         gph_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_9 * gpg0_603[k]
                    - f_10 * gpg1_603[k]
                    + f_4 * pc_x[k] * gph_843[k];

        t_1124[k] = f_9 * gpg0_604[k]
                    - f_10 * gpg1_604[k]
                    + f_4 * pc_x[k] * gph_844[k];

        t_1125[k] = f_9 * gpg0_605[k]
                    - f_10 * gpg1_605[k]
                    + f_4 * pc_x[k] * gph_845[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, pc_x, gpg0_606, gpg0_607, gpg0_608, gpg1_606, \
                         gpg1_607, gpg1_608, gph_846, gph_847, \
                         gph_848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_7 * gpg0_606[k]
                    - f_8 * gpg1_606[k]
                    + f_4 * pc_x[k] * gph_846[k];

        t_1127[k] = f_7 * gpg0_607[k]
                    - f_8 * gpg1_607[k]
                    + f_4 * pc_x[k] * gph_847[k];

        t_1128[k] = f_7 * gpg0_608[k]
                    - f_8 * gpg1_608[k]
                    + f_4 * pc_x[k] * gph_848[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, pc_x, gpg0_609, gpg0_610, gpg0_611, gpg1_609, \
                         gpg1_610, gpg1_611, gph_849, gph_850, \
                         gph_851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_7 * gpg0_609[k]
                    - f_8 * gpg1_609[k]
                    + f_4 * pc_x[k] * gph_849[k];

        t_1130[k] = f_5 * gpg0_610[k]
                    - f_6 * gpg1_610[k]
                    + f_4 * pc_x[k] * gph_850[k];

        t_1131[k] = f_5 * gpg0_611[k]
                    - f_6 * gpg1_611[k]
                    + f_4 * pc_x[k] * gph_851[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, t_1135, pc_x, gpg0_612, gpg0_613, gpg0_614, \
                         gpg1_612, gpg1_613, gpg1_614, gph_852, gph_853, gph_854, \
                         gph_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_5 * gpg0_612[k]
                    - f_6 * gpg1_612[k]
                    + f_4 * pc_x[k] * gph_852[k];

        t_1133[k] = f_5 * gpg0_613[k]
                    - f_6 * gpg1_613[k]
                    + f_4 * pc_x[k] * gph_853[k];

        t_1134[k] = f_5 * gpg0_614[k]
                    - f_6 * gpg1_614[k]
                    + f_4 * pc_x[k] * gph_854[k];

        t_1135[k] = f_4 * pc_x[k] * gph_855[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, t_1139, t_1140, pc_x, gph_856, gph_857, \
                         gph_858, gph_859, gph_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_4 * pc_x[k] * gph_856[k];

        t_1137[k] = f_4 * pc_x[k] * gph_857[k];

        t_1138[k] = f_4 * pc_x[k] * gph_858[k];

        t_1139[k] = f_4 * pc_x[k] * gph_859[k];

        t_1140[k] = f_4 * pc_x[k] * gph_860[k];
    }

#pragma omp simd aligned(t_1141, t_1142, t_1143, pc_y, pc_z, fph_540, fph_603, fph_605, \
                         gsh_288, gsh_290, gpg0_610, gpg0_612, gpg1_610, gpg1_612, gph_855, \
                         gph_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1141[k] = f_1 * fph_603[k]
                    + f_1 * gsh_288[k]
                    + f_2 * gpg0_610[k]
                    - f_3 * gpg1_610[k]
                    + f_4 * pc_y[k] * gph_855[k];

        t_1142[k] = f_13 * fph_540[k]
                    + f_4 * pc_z[k] * gph_855[k];

        t_1143[k] = f_1 * fph_605[k]
                    + f_1 * gsh_290[k]
                    + f_9 * gpg0_612[k]
                    - f_10 * gpg1_612[k]
                    + f_4 * pc_y[k] * gph_857[k];
    }

#pragma omp simd aligned(t_1144, t_1145, pc_y, fph_606, fph_607, gsh_291, gsh_292, gpg0_613, \
                         gpg0_614, gpg1_613, gpg1_614, gph_858, \
                         gph_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_1 * fph_606[k]
                    + f_1 * gsh_291[k]
                    + f_7 * gpg0_613[k]
                    - f_8 * gpg1_613[k]
                    + f_4 * pc_y[k] * gph_858[k];

        t_1145[k] = f_1 * fph_607[k]
                    + f_1 * gsh_292[k]
                    + f_5 * gpg0_614[k]
                    - f_6 * gpg1_614[k]
                    + f_4 * pc_y[k] * gph_859[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, pa_y, pc_y, pc_z, fpi0_812, fph_545, fph_608, \
                         fpi1_812, gsh_293, gpg0_614, gpg1_614, \
                         gph_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_1 * fph_608[k]
                    + f_1 * gsh_293[k]
                    + f_4 * pc_y[k] * gph_860[k];

        t_1147[k] = f_13 * fph_545[k]
                    + f_2 * gpg0_614[k]
                    - f_3 * gpg1_614[k]
                    + f_4 * pc_z[k] * gph_860[k];

        t_1148[k] = pa_y[k] * fpi0_812[k]
                    - f_11 * pc_y[k] * fpi1_812[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, pa_y, pc_x, pc_y, fpi0_814, fpi1_814, \
                         gpg0_616, gpg0_618, gpg1_616, gpg1_618, gph_862, \
                         gph_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_17 * gpg0_616[k]
                    - f_18 * gpg1_616[k]
                    + f_4 * pc_x[k] * gph_862[k];

        t_1150[k] = pa_y[k] * fpi0_814[k]
                    - f_11 * pc_y[k] * fpi1_814[k];

        t_1151[k] = f_9 * gpg0_618[k]
                    - f_10 * gpg1_618[k]
                    + f_4 * pc_x[k] * gph_864[k];
    }

#pragma omp simd aligned(t_1152, t_1153, t_1154, pa_y, pc_x, pc_y, fpi0_817, fpi1_817, \
                         gpg0_619, gpg0_621, gpg1_619, gpg1_621, gph_865, \
                         gph_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1152[k] = f_9 * gpg0_619[k]
                    - f_10 * gpg1_619[k]
                    + f_4 * pc_x[k] * gph_865[k];

        t_1153[k] = pa_y[k] * fpi0_817[k]
                    - f_11 * pc_y[k] * fpi1_817[k];

        t_1154[k] = f_7 * gpg0_621[k]
                    - f_8 * gpg1_621[k]
                    + f_4 * pc_x[k] * gph_867[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, pa_y, pc_x, pc_y, fpi0_821, fpi1_821, \
                         gpg0_622, gpg0_623, gpg1_622, gpg1_623, gph_868, \
                         gph_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_7 * gpg0_622[k]
                    - f_8 * gpg1_622[k]
                    + f_4 * pc_x[k] * gph_868[k];

        t_1156[k] = f_7 * gpg0_623[k]
                    - f_8 * gpg1_623[k]
                    + f_4 * pc_x[k] * gph_869[k];

        t_1157[k] = pa_y[k] * fpi0_821[k]
                    - f_11 * pc_y[k] * fpi1_821[k];
    }

#pragma omp simd aligned(t_1158, t_1159, t_1160, pc_x, gpg0_625, gpg0_626, gpg0_627, gpg1_625, \
                         gpg1_626, gpg1_627, gph_871, gph_872, \
                         gph_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1158[k] = f_5 * gpg0_625[k]
                    - f_6 * gpg1_625[k]
                    + f_4 * pc_x[k] * gph_871[k];

        t_1159[k] = f_5 * gpg0_626[k]
                    - f_6 * gpg1_626[k]
                    + f_4 * pc_x[k] * gph_872[k];

        t_1160[k] = f_5 * gpg0_627[k]
                    - f_6 * gpg1_627[k]
                    + f_4 * pc_x[k] * gph_873[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, t_1164, t_1165, pa_y, pc_x, pc_y, fpi0_826, \
                         fpi1_826, gpg0_628, gpg1_628, gph_874, gph_876, gph_877, \
                         gph_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_5 * gpg0_628[k]
                    - f_6 * gpg1_628[k]
                    + f_4 * pc_x[k] * gph_874[k];

        t_1162[k] = pa_y[k] * fpi0_826[k]
                    - f_11 * pc_y[k] * fpi1_826[k];

        t_1163[k] = f_4 * pc_x[k] * gph_876[k];

        t_1164[k] = f_4 * pc_x[k] * gph_877[k];

        t_1165[k] = f_4 * pc_x[k] * gph_878[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, pa_y, pc_x, pc_y, fpi0_833, fph_624, \
                         fpi1_833, gph_879, gph_880, gph_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_4 * pc_x[k] * gph_879[k];

        t_1167[k] = f_4 * pc_x[k] * gph_880[k];

        t_1168[k] = f_4 * pc_x[k] * gph_881[k];

        t_1169[k] = pa_y[k] * fpi0_833[k]
                    + f_14 * fph_624[k]
                    - f_11 * pc_y[k] * fpi1_833[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, pa_y, pc_y, pc_z, fpi0_835, fpi0_836, \
                         fph_561, fph_626, fph_627, fpi1_835, fpi1_836, gsh_288, \
                         gph_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_13 * fph_561[k]
                    + f_1 * gsh_288[k]
                    + f_4 * pc_z[k] * gph_876[k];

        t_1171[k] = pa_y[k] * fpi0_835[k]
                    + f_0 * fph_626[k]
                    - f_11 * pc_y[k] * fpi1_835[k];

        t_1172[k] = pa_y[k] * fpi0_836[k]
                    + f_13 * fph_627[k]
                    - f_11 * pc_y[k] * fpi1_836[k];
    }
}

static auto
compute_prim_gpi_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pb, const size_t pc,
                                                           const size_t fpi0, const size_t fph,
                                                           const size_t fpi1, const size_t gsi0,
                                                           const size_t gsh, const size_t gsi1,
                                                           const size_t gpg0, const size_t gpg1,
                                                           const size_t gph, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
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
    const auto f_14 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;

    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fpi0_837 = buffer.data(fpi0 + 837);
    const auto *fpi0_839 = buffer.data(fpi0 + 839);

    const auto *fph_628 = buffer.data(fph + 628);
    const auto *fph_629 = buffer.data(fph + 629);

    const auto *fpi1_837 = buffer.data(fpi1 + 837);
    const auto *fpi1_839 = buffer.data(fpi1 + 839);

    const auto *gsi0_392 = buffer.data(gsi0 + 392);
    const auto *gsi0_394 = buffer.data(gsi0 + 394);
    const auto *gsi0_395 = buffer.data(gsi0 + 395);
    const auto *gsi0_397 = buffer.data(gsi0 + 397);
    const auto *gsi0_398 = buffer.data(gsi0 + 398);
    const auto *gsi0_399 = buffer.data(gsi0 + 399);
    const auto *gsi0_401 = buffer.data(gsi0 + 401);
    const auto *gsi0_402 = buffer.data(gsi0 + 402);
    const auto *gsi0_403 = buffer.data(gsi0 + 403);
    const auto *gsi0_404 = buffer.data(gsi0 + 404);
    const auto *gsi0_406 = buffer.data(gsi0 + 406);
    const auto *gsi0_413 = buffer.data(gsi0 + 413);
    const auto *gsi0_414 = buffer.data(gsi0 + 414);
    const auto *gsi0_415 = buffer.data(gsi0 + 415);
    const auto *gsi0_416 = buffer.data(gsi0 + 416);
    const auto *gsi0_417 = buffer.data(gsi0 + 417);
    const auto *gsi0_419 = buffer.data(gsi0 + 419);

    const auto *gsh_294 = buffer.data(gsh + 294);
    const auto *gsh_296 = buffer.data(gsh + 296);
    const auto *gsh_297 = buffer.data(gsh + 297);
    const auto *gsh_299 = buffer.data(gsh + 299);
    const auto *gsh_300 = buffer.data(gsh + 300);
    const auto *gsh_301 = buffer.data(gsh + 301);
    const auto *gsh_303 = buffer.data(gsh + 303);
    const auto *gsh_304 = buffer.data(gsh + 304);
    const auto *gsh_305 = buffer.data(gsh + 305);
    const auto *gsh_306 = buffer.data(gsh + 306);
    const auto *gsh_308 = buffer.data(gsh + 308);
    const auto *gsh_309 = buffer.data(gsh + 309);
    const auto *gsh_310 = buffer.data(gsh + 310);
    const auto *gsh_311 = buffer.data(gsh + 311);
    const auto *gsh_312 = buffer.data(gsh + 312);
    const auto *gsh_313 = buffer.data(gsh + 313);
    const auto *gsh_314 = buffer.data(gsh + 314);

    const auto *gsi1_392 = buffer.data(gsi1 + 392);
    const auto *gsi1_394 = buffer.data(gsi1 + 394);
    const auto *gsi1_395 = buffer.data(gsi1 + 395);
    const auto *gsi1_397 = buffer.data(gsi1 + 397);
    const auto *gsi1_398 = buffer.data(gsi1 + 398);
    const auto *gsi1_399 = buffer.data(gsi1 + 399);
    const auto *gsi1_401 = buffer.data(gsi1 + 401);
    const auto *gsi1_402 = buffer.data(gsi1 + 402);
    const auto *gsi1_403 = buffer.data(gsi1 + 403);
    const auto *gsi1_404 = buffer.data(gsi1 + 404);
    const auto *gsi1_406 = buffer.data(gsi1 + 406);
    const auto *gsi1_413 = buffer.data(gsi1 + 413);
    const auto *gsi1_414 = buffer.data(gsi1 + 414);
    const auto *gsi1_415 = buffer.data(gsi1 + 415);
    const auto *gsi1_416 = buffer.data(gsi1 + 416);
    const auto *gsi1_417 = buffer.data(gsi1 + 417);
    const auto *gsi1_419 = buffer.data(gsi1 + 419);

    const auto *gpg0_648 = buffer.data(gpg0 + 648);
    const auto *gpg0_651 = buffer.data(gpg0 + 651);
    const auto *gpg0_652 = buffer.data(gpg0 + 652);
    const auto *gpg0_655 = buffer.data(gpg0 + 655);
    const auto *gpg0_656 = buffer.data(gpg0 + 656);
    const auto *gpg0_657 = buffer.data(gpg0 + 657);
    const auto *gpg0_660 = buffer.data(gpg0 + 660);
    const auto *gpg0_662 = buffer.data(gpg0 + 662);
    const auto *gpg0_663 = buffer.data(gpg0 + 663);
    const auto *gpg0_665 = buffer.data(gpg0 + 665);
    const auto *gpg0_666 = buffer.data(gpg0 + 666);
    const auto *gpg0_667 = buffer.data(gpg0 + 667);
    const auto *gpg0_669 = buffer.data(gpg0 + 669);
    const auto *gpg0_670 = buffer.data(gpg0 + 670);
    const auto *gpg0_671 = buffer.data(gpg0 + 671);
    const auto *gpg0_672 = buffer.data(gpg0 + 672);
    const auto *gpg0_673 = buffer.data(gpg0 + 673);
    const auto *gpg0_674 = buffer.data(gpg0 + 674);

    const auto *gpg1_648 = buffer.data(gpg1 + 648);
    const auto *gpg1_651 = buffer.data(gpg1 + 651);
    const auto *gpg1_652 = buffer.data(gpg1 + 652);
    const auto *gpg1_655 = buffer.data(gpg1 + 655);
    const auto *gpg1_656 = buffer.data(gpg1 + 656);
    const auto *gpg1_657 = buffer.data(gpg1 + 657);
    const auto *gpg1_660 = buffer.data(gpg1 + 660);
    const auto *gpg1_662 = buffer.data(gpg1 + 662);
    const auto *gpg1_663 = buffer.data(gpg1 + 663);
    const auto *gpg1_665 = buffer.data(gpg1 + 665);
    const auto *gpg1_666 = buffer.data(gpg1 + 666);
    const auto *gpg1_667 = buffer.data(gpg1 + 667);
    const auto *gpg1_669 = buffer.data(gpg1 + 669);
    const auto *gpg1_670 = buffer.data(gpg1 + 670);
    const auto *gpg1_671 = buffer.data(gpg1 + 671);
    const auto *gpg1_672 = buffer.data(gpg1 + 672);
    const auto *gpg1_673 = buffer.data(gpg1 + 673);
    const auto *gpg1_674 = buffer.data(gpg1 + 674);

    const auto *gph_881 = buffer.data(gph + 881);
    const auto *gph_882 = buffer.data(gph + 882);
    const auto *gph_884 = buffer.data(gph + 884);
    const auto *gph_887 = buffer.data(gph + 887);
    const auto *gph_891 = buffer.data(gph + 891);
    const auto *gph_897 = buffer.data(gph + 897);
    const auto *gph_898 = buffer.data(gph + 898);
    const auto *gph_899 = buffer.data(gph + 899);
    const auto *gph_900 = buffer.data(gph + 900);
    const auto *gph_901 = buffer.data(gph + 901);
    const auto *gph_902 = buffer.data(gph + 902);
    const auto *gph_903 = buffer.data(gph + 903);
    const auto *gph_905 = buffer.data(gph + 905);
    const auto *gph_906 = buffer.data(gph + 906);
    const auto *gph_908 = buffer.data(gph + 908);
    const auto *gph_909 = buffer.data(gph + 909);
    const auto *gph_910 = buffer.data(gph + 910);
    const auto *gph_912 = buffer.data(gph + 912);
    const auto *gph_913 = buffer.data(gph + 913);
    const auto *gph_914 = buffer.data(gph + 914);
    const auto *gph_915 = buffer.data(gph + 915);
    const auto *gph_918 = buffer.data(gph + 918);
    const auto *gph_919 = buffer.data(gph + 919);
    const auto *gph_920 = buffer.data(gph + 920);
    const auto *gph_921 = buffer.data(gph + 921);
    const auto *gph_922 = buffer.data(gph + 922);
    const auto *gph_923 = buffer.data(gph + 923);
    const auto *gph_924 = buffer.data(gph + 924);
    const auto *gph_926 = buffer.data(gph + 926);
    const auto *gph_927 = buffer.data(gph + 927);
    const auto *gph_929 = buffer.data(gph + 929);
    const auto *gph_930 = buffer.data(gph + 930);
    const auto *gph_931 = buffer.data(gph + 931);
    const auto *gph_933 = buffer.data(gph + 933);
    const auto *gph_934 = buffer.data(gph + 934);
    const auto *gph_935 = buffer.data(gph + 935);
    const auto *gph_936 = buffer.data(gph + 936);
    const auto *gph_938 = buffer.data(gph + 938);
    const auto *gph_939 = buffer.data(gph + 939);
    const auto *gph_940 = buffer.data(gph + 940);
    const auto *gph_941 = buffer.data(gph + 941);
    const auto *gph_942 = buffer.data(gph + 942);
    const auto *gph_943 = buffer.data(gph + 943);
    const auto *gph_944 = buffer.data(gph + 944);

#pragma omp simd aligned(t_1173, t_1174, t_1175, pa_y, pc_y, fpi0_837, fpi0_839, fph_628, \
                         fph_629, fpi1_837, fpi1_839, gph_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = pa_y[k] * fpi0_837[k]
                    + f_12 * fph_628[k]
                    - f_11 * pc_y[k] * fpi1_837[k];

        t_1174[k] = f_1 * fph_629[k]
                    + f_4 * pc_y[k] * gph_881[k];

        t_1175[k] = pa_y[k] * fpi0_839[k]
                    - f_11 * pc_y[k] * fpi1_839[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, pb_x, pc_x, pc_y, gsi0_392, gsi0_394, \
                         gsh_294, gsh_296, gsi1_392, gsi1_394, \
                         gph_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = pb_x[k] * gsi0_392[k]
                    + f_14 * gsh_294[k]
                    - f_11 * pc_x[k] * gsi1_392[k];

        t_1177[k] = f_4 * pc_y[k] * gph_882[k];

        t_1178[k] = pb_x[k] * gsi0_394[k]
                    + f_19 * gsh_296[k]
                    - f_11 * pc_x[k] * gsi1_394[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pb_x, pc_x, pc_y, gsi0_395, gsi0_397, \
                         gsh_297, gsh_299, gsi1_395, gsi1_397, \
                         gph_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = pb_x[k] * gsi0_395[k]
                    + f_0 * gsh_297[k]
                    - f_11 * pc_x[k] * gsi1_395[k];

        t_1180[k] = f_4 * pc_y[k] * gph_884[k];

        t_1181[k] = pb_x[k] * gsi0_397[k]
                    + f_0 * gsh_299[k]
                    - f_11 * pc_x[k] * gsi1_397[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pb_x, pc_x, pc_y, gsi0_398, gsi0_399, \
                         gsh_300, gsh_301, gsi1_398, gsi1_399, \
                         gph_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = pb_x[k] * gsi0_398[k]
                    + f_13 * gsh_300[k]
                    - f_11 * pc_x[k] * gsi1_398[k];

        t_1183[k] = pb_x[k] * gsi0_399[k]
                    + f_13 * gsh_301[k]
                    - f_11 * pc_x[k] * gsi1_399[k];

        t_1184[k] = f_4 * pc_y[k] * gph_887[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pb_x, pc_x, gsi0_401, gsi0_402, gsi0_403, \
                         gsh_303, gsh_304, gsh_305, gsi1_401, gsi1_402, \
                         gsi1_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = pb_x[k] * gsi0_401[k]
                    + f_13 * gsh_303[k]
                    - f_11 * pc_x[k] * gsi1_401[k];

        t_1186[k] = pb_x[k] * gsi0_402[k]
                    + f_12 * gsh_304[k]
                    - f_11 * pc_x[k] * gsi1_402[k];

        t_1187[k] = pb_x[k] * gsi0_403[k]
                    + f_12 * gsh_305[k]
                    - f_11 * pc_x[k] * gsi1_403[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, t_1191, pb_x, pc_x, pc_y, gsi0_404, gsi0_406, \
                         gsh_306, gsh_308, gsh_309, gsi1_404, gsi1_406, gph_891, \
                         gph_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = pb_x[k] * gsi0_404[k]
                    + f_12 * gsh_306[k]
                    - f_11 * pc_x[k] * gsi1_404[k];

        t_1189[k] = f_4 * pc_y[k] * gph_891[k];

        t_1190[k] = pb_x[k] * gsi0_406[k]
                    + f_12 * gsh_308[k]
                    - f_11 * pc_x[k] * gsi1_406[k];

        t_1191[k] = f_1 * gsh_309[k]
                    + f_4 * pc_x[k] * gph_897[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, t_1195, t_1196, pc_x, gsh_310, gsh_311, \
                         gsh_312, gsh_313, gsh_314, gph_898, gph_899, gph_900, gph_901, \
                         gph_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = f_1 * gsh_310[k]
                    + f_4 * pc_x[k] * gph_898[k];

        t_1193[k] = f_1 * gsh_311[k]
                    + f_4 * pc_x[k] * gph_899[k];

        t_1194[k] = f_1 * gsh_312[k]
                    + f_4 * pc_x[k] * gph_900[k];

        t_1195[k] = f_1 * gsh_313[k]
                    + f_4 * pc_x[k] * gph_901[k];

        t_1196[k] = f_1 * gsh_314[k]
                    + f_4 * pc_x[k] * gph_902[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, t_1200, pb_x, pc_x, gsi0_413, gsi0_414, \
                         gsi0_415, gsi0_416, gsi1_413, gsi1_414, gsi1_415, \
                         gsi1_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = pb_x[k] * gsi0_413[k]
                    - f_11 * pc_x[k] * gsi1_413[k];

        t_1198[k] = pb_x[k] * gsi0_414[k]
                    - f_11 * pc_x[k] * gsi1_414[k];

        t_1199[k] = pb_x[k] * gsi0_415[k]
                    - f_11 * pc_x[k] * gsi1_415[k];

        t_1200[k] = pb_x[k] * gsi0_416[k]
                    - f_11 * pc_x[k] * gsi1_416[k];
    }

#pragma omp simd aligned(t_1201, t_1202, t_1203, t_1204, pb_x, pb_y, pc_x, pc_y, gsi0_392, \
                         gsi0_417, gsi0_419, gsi1_392, gsi1_417, gsi1_419, \
                         gph_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1201[k] = pb_x[k] * gsi0_417[k]
                    - f_11 * pc_x[k] * gsi1_417[k];

        t_1202[k] = f_4 * pc_y[k] * gph_902[k];

        t_1203[k] = pb_x[k] * gsi0_419[k]
                    - f_11 * pc_x[k] * gsi1_419[k];

        t_1204[k] = pb_y[k] * gsi0_392[k]
                    - f_11 * pc_y[k] * gsi1_392[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, t_1208, pb_y, pc_x, pc_y, gsi0_394, gsh_294, \
                         gsh_296, gsi1_394, gpg0_648, gpg1_648, gph_903, gph_905, \
                         gph_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_1 * gsh_294[k]
                    + f_4 * pc_y[k] * gph_903[k];

        t_1206[k] = pb_y[k] * gsi0_394[k]
                    - f_11 * pc_y[k] * gsi1_394[k];

        t_1207[k] = f_9 * gpg0_648[k]
                    - f_10 * gpg1_648[k]
                    + f_4 * pc_x[k] * gph_906[k];

        t_1208[k] = f_1 * gsh_296[k]
                    + f_4 * pc_y[k] * gph_905[k];
    }

#pragma omp simd aligned(t_1209, t_1210, t_1211, pb_y, pc_x, pc_y, gsi0_397, gsi1_397, \
                         gpg0_651, gpg0_652, gpg1_651, gpg1_652, gph_909, \
                         gph_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1209[k] = pb_y[k] * gsi0_397[k]
                    - f_11 * pc_y[k] * gsi1_397[k];

        t_1210[k] = f_7 * gpg0_651[k]
                    - f_8 * gpg1_651[k]
                    + f_4 * pc_x[k] * gph_909[k];

        t_1211[k] = f_7 * gpg0_652[k]
                    - f_8 * gpg1_652[k]
                    + f_4 * pc_x[k] * gph_910[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, pb_y, pc_x, pc_y, gsi0_401, gsh_299, \
                         gsi1_401, gpg0_655, gpg1_655, gph_908, \
                         gph_913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_1 * gsh_299[k]
                    + f_4 * pc_y[k] * gph_908[k];

        t_1213[k] = pb_y[k] * gsi0_401[k]
                    - f_11 * pc_y[k] * gsi1_401[k];

        t_1214[k] = f_5 * gpg0_655[k]
                    - f_6 * gpg1_655[k]
                    + f_4 * pc_x[k] * gph_913[k];
    }

#pragma omp simd aligned(t_1215, t_1216, t_1217, pc_x, pc_y, gsh_303, gpg0_656, gpg0_657, \
                         gpg1_656, gpg1_657, gph_912, gph_914, \
                         gph_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1215[k] = f_5 * gpg0_656[k]
                    - f_6 * gpg1_656[k]
                    + f_4 * pc_x[k] * gph_914[k];

        t_1216[k] = f_5 * gpg0_657[k]
                    - f_6 * gpg1_657[k]
                    + f_4 * pc_x[k] * gph_915[k];

        t_1217[k] = f_1 * gsh_303[k]
                    + f_4 * pc_y[k] * gph_912[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, t_1221, t_1222, t_1223, pb_y, pc_x, pc_y, \
                         gsi0_406, gsi1_406, gph_918, gph_919, gph_920, gph_921, \
                         gph_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = pb_y[k] * gsi0_406[k]
                    - f_11 * pc_y[k] * gsi1_406[k];

        t_1219[k] = f_4 * pc_x[k] * gph_918[k];

        t_1220[k] = f_4 * pc_x[k] * gph_919[k];

        t_1221[k] = f_4 * pc_x[k] * gph_920[k];

        t_1222[k] = f_4 * pc_x[k] * gph_921[k];

        t_1223[k] = f_4 * pc_x[k] * gph_922[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pb_y, pc_x, pc_y, gsi0_413, gsi0_414, \
                         gsh_309, gsh_310, gsi1_413, gsi1_414, \
                         gph_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = f_4 * pc_x[k] * gph_923[k];

        t_1225[k] = pb_y[k] * gsi0_413[k]
                    + f_14 * gsh_309[k]
                    - f_11 * pc_y[k] * gsi1_413[k];

        t_1226[k] = pb_y[k] * gsi0_414[k]
                    + f_19 * gsh_310[k]
                    - f_11 * pc_y[k] * gsi1_414[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, pb_y, pc_y, gsi0_415, gsi0_416, gsi0_417, \
                         gsh_311, gsh_312, gsh_313, gsi1_415, gsi1_416, \
                         gsi1_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = pb_y[k] * gsi0_415[k]
                    + f_0 * gsh_311[k]
                    - f_11 * pc_y[k] * gsi1_415[k];

        t_1228[k] = pb_y[k] * gsi0_416[k]
                    + f_13 * gsh_312[k]
                    - f_11 * pc_y[k] * gsi1_416[k];

        t_1229[k] = pb_y[k] * gsi0_417[k]
                    + f_12 * gsh_313[k]
                    - f_11 * pc_y[k] * gsi1_417[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, pb_y, pc_x, pc_y, gsi0_419, gsh_314, \
                         gsi1_419, gpg0_660, gpg1_660, gph_923, \
                         gph_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_1 * gsh_314[k]
                    + f_4 * pc_y[k] * gph_923[k];

        t_1231[k] = pb_y[k] * gsi0_419[k]
                    - f_11 * pc_y[k] * gsi1_419[k];

        t_1232[k] = f_2 * gpg0_660[k]
                    - f_3 * gpg1_660[k]
                    + f_4 * pc_x[k] * gph_924[k];

        t_1233[k] = f_4 * pc_y[k] * gph_924[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, t_1237, pc_x, pc_y, gpg0_662, gpg0_663, \
                         gpg0_665, gpg1_662, gpg1_663, gpg1_665, gph_926, gph_927, \
                         gph_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = f_17 * gpg0_662[k]
                    - f_18 * gpg1_662[k]
                    + f_4 * pc_x[k] * gph_926[k];

        t_1235[k] = f_9 * gpg0_663[k]
                    - f_10 * gpg1_663[k]
                    + f_4 * pc_x[k] * gph_927[k];

        t_1236[k] = f_4 * pc_y[k] * gph_926[k];

        t_1237[k] = f_9 * gpg0_665[k]
                    - f_10 * gpg1_665[k]
                    + f_4 * pc_x[k] * gph_929[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, t_1241, pc_x, pc_y, gpg0_666, gpg0_667, \
                         gpg0_669, gpg1_666, gpg1_667, gpg1_669, gph_929, gph_930, gph_931, \
                         gph_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = f_7 * gpg0_666[k]
                    - f_8 * gpg1_666[k]
                    + f_4 * pc_x[k] * gph_930[k];

        t_1239[k] = f_7 * gpg0_667[k]
                    - f_8 * gpg1_667[k]
                    + f_4 * pc_x[k] * gph_931[k];

        t_1240[k] = f_4 * pc_y[k] * gph_929[k];

        t_1241[k] = f_7 * gpg0_669[k]
                    - f_8 * gpg1_669[k]
                    + f_4 * pc_x[k] * gph_933[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, pc_x, pc_y, gpg0_670, gpg0_671, \
                         gpg0_672, gpg1_670, gpg1_671, gpg1_672, gph_933, gph_934, gph_935, \
                         gph_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = f_5 * gpg0_670[k]
                    - f_6 * gpg1_670[k]
                    + f_4 * pc_x[k] * gph_934[k];

        t_1243[k] = f_5 * gpg0_671[k]
                    - f_6 * gpg1_671[k]
                    + f_4 * pc_x[k] * gph_935[k];

        t_1244[k] = f_5 * gpg0_672[k]
                    - f_6 * gpg1_672[k]
                    + f_4 * pc_x[k] * gph_936[k];

        t_1245[k] = f_4 * pc_y[k] * gph_933[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, t_1250, t_1251, pc_x, gpg0_674, \
                         gpg1_674, gph_938, gph_939, gph_940, gph_941, gph_942, \
                         gph_943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = f_5 * gpg0_674[k]
                    - f_6 * gpg1_674[k]
                    + f_4 * pc_x[k] * gph_938[k];

        t_1247[k] = f_4 * pc_x[k] * gph_939[k];

        t_1248[k] = f_4 * pc_x[k] * gph_940[k];

        t_1249[k] = f_4 * pc_x[k] * gph_941[k];

        t_1250[k] = f_4 * pc_x[k] * gph_942[k];

        t_1251[k] = f_4 * pc_x[k] * gph_943[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, pc_x, pc_y, gpg0_670, gpg0_671, \
                         gpg0_672, gpg1_670, gpg1_671, gpg1_672, gph_939, gph_940, gph_941, \
                         gph_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = f_4 * pc_x[k] * gph_944[k];

        t_1253[k] = f_2 * gpg0_670[k]
                    - f_3 * gpg1_670[k]
                    + f_4 * pc_y[k] * gph_939[k];

        t_1254[k] = f_17 * gpg0_671[k]
                    - f_18 * gpg1_671[k]
                    + f_4 * pc_y[k] * gph_940[k];

        t_1255[k] = f_9 * gpg0_672[k]
                    - f_10 * gpg1_672[k]
                    + f_4 * pc_y[k] * gph_941[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pc_y, pc_z, fph_629, gsh_314, \
                         gpg0_673, gpg0_674, gpg1_673, gpg1_674, gph_942, gph_943, \
                         gph_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = f_7 * gpg0_673[k]
                    - f_8 * gpg1_673[k]
                    + f_4 * pc_y[k] * gph_942[k];

        t_1257[k] = f_5 * gpg0_674[k]
                    - f_6 * gpg1_674[k]
                    + f_4 * pc_y[k] * gph_943[k];

        t_1258[k] = f_4 * pc_y[k] * gph_944[k];

        t_1259[k] = f_0 * fph_629[k]
                    + f_1 * gsh_314[k]
                    + f_2 * gpg0_674[k]
                    - f_3 * gpg1_674[k]
                    + f_4 * pc_z[k] * gph_944[k];
    }
}

auto
compute_prim_gpi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t dpi0,
                                                   const size_t dpi1, const size_t fpi0,
                                                   const size_t fph, const size_t fpi1,
                                                   const size_t gsi0, const size_t gsh,
                                                   const size_t gsi1, const size_t gpg0,
                                                   const size_t gpg1, const size_t gph,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gpi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, fpi0,
                                                              fph, fpi1, gsi0, gsh, gsi1, gpg0,
                                                              gpg1, gph, ncols, gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, dpi0,
                                                              dpi1, fpi0, fph, fpi1, gsi0, gsh,
                                                              gsi1, gpg0, gpg1, gph, ncols,
                                                              gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, dpi0,
                                                              dpi1, fpi0, fph, fpi1, gsi0, gsh,
                                                              gsi1, gpg0, gpg1, gph, ncols,
                                                              gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, dpi0,
                                                              dpi1, fpi0, fph, fpi1, gsi0, gsh,
                                                              gsi1, gpg0, gpg1, gph, ncols,
                                                              gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, dpi0,
                                                              dpi1, fpi0, fph, fpi1, gsi0, gsh,
                                                              gsi1, gpg0, gpg1, gph, ncols,
                                                              gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, fpi0, fph,
                                                              fpi1, gsh, gpg0, gpg1, gph, ncols,
                                                              gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pb, pc, fpi0,
                                                              fph, fpi1, gsi0, gsh, gsi1, gpg0,
                                                              gpg1, gph, ncols, gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece7(buffer, target, pa, pb, pc, fpi0,
                                                              fph, fpi1, gsi0, gsh, gsi1, gpg0,
                                                              gpg1, gph, ncols, gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece8(buffer, target, pa, pb, pc, dpi0,
                                                              dpi1, fpi0, fph, fpi1, gsi0, gsh,
                                                              gsi1, gpg0, gpg1, gph, ncols,
                                                              gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece9(buffer, target, pa, pb, pc, dpi0,
                                                              dpi1, fpi0, fph, fpi1, gsi0, gsh,
                                                              gsi1, gpg0, gpg1, gph, ncols,
                                                              gamma, p, q);

    compute_prim_gpi_three_center_electron_repulsion_0_piece10(buffer, target, pa, pb, pc, fpi0,
                                                               fph, fpi1, gsi0, gsh, gsi1, gpg0,
                                                               gpg1, gph, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
